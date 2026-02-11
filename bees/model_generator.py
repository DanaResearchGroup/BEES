#!/usr/bin/env python3

"""
Model Generator Module
----------------------
Coordinates reaction template generation with kinetic database queries
to produce complete biochemical reaction models.

This module orchestrates:
1. Loading kinetic parameters from database
2. Generating reaction templates from EC numbers
3. Inferring products when not in database
4. Creating complete reaction objects with kinetics
"""

from typing import List, Dict, Optional
from dataclasses import dataclass
import os
from types import SimpleNamespace
import json
import time

from bees.reaction_template import (
    create_reaction_from_database,
    ReactionTemplate,
)
from db.reaction_database import ReactionDatabase, KineticData
from bees.common import (
    BEES_PATH,
    GENERAL_COFACTORS,
    get_ontology_equivalents,
    get_coenzyme_like_flags,
)
from bees.reaction_utils import (
    get_ec_aliases,
    validate_reaction_reactants
)
from bees.kinetics_estimator import build_estimator


@dataclass
class GeneratedReaction:
    """
    Complete reaction with template and kinetics.
    
    Attributes:
        enzyme_label (str): Enzyme name
        substrate_label (str): Substrate name
        template (ReactionTemplate): Reaction template with EC classification
        kinetics (KineticData): Kinetic parameters from database (or None)
        reactant_labels (List[str]): List of reactant names
        product_labels (List[str]): List of product names
        stoichiometry (Dict[str, int]): Stoichiometric coefficients
        rate_law (str or None): Rate law type if available from database, None otherwise
    """
    enzyme_label: str
    substrate_label: str
    ec_number: Optional[str]
    template: ReactionTemplate
    kinetics: Optional[KineticData]
    reactant_labels: List[str]
    product_labels: List[str]
    stoichiometry: Dict[str, int]
    rate_law: Optional[str] = None  # Only set if kinetics available from database
    
    def __repr__(self):
        return f"Reaction: {' + '.join(self.reactant_labels)} → {' + '.join(self.product_labels)}"


class ModelGenerator:
    """
    Generates biochemical reaction models from species and enzymes.
    
    Workflow:
    1. Load kinetic database
    2. For each enzyme-substrate pair:
       a. Query database for kinetic parameters
       b. Generate reaction template from EC number
       c. Infer products if not in database
       d. Create complete reaction
    3. Export model to output files
    
    Attributes:
        bees_object: Validated input object containing species, enzymes, etc.
        logger: Logger instance for tracking operations
        output_directory (str): Directory for output files
        kinetic_db (ReactionDatabase): Loaded reaction database
        reactions (List[GeneratedReaction]): Generated reactions
    """
    
    def __init__(self, bees_object, logger, output_directory):
        """
        Initialize model generator.
        
        Args:
            bees_object: Validated BEES input object
            logger: Logger instance
            output_directory (str): Path to output directory
        """
        self.bees_object = bees_object
        self.logger = logger
        self.output_directory = output_directory
        self.kinetic_db = None
        self.reactions: List[GeneratedReaction] = []
        self.kinetics_estimator = None
    
    def load_kinetic_database(self, db_path: str, ontology: Optional[Dict[str, List[str]]] = None) -> int:
        """
        Load kinetic database from CSV file.
        
        Args:
            db_path (str): Path to kinetic database CSV
            ontology (dict, optional): Chemical ontology for alias matching.
            
        Returns:
            int: Number of reactions loaded
            
        Raises:
            FileNotFoundError: If database file not found
            ValueError: If CSV format is invalid
        """
        self.logger.info(f"Loading kinetic database from {db_path}")
        self.kinetic_db = ReactionDatabase(logger=self.logger, ontology=ontology)
        num_reactions = self.kinetic_db.load_from_csv(db_path)
        
        # Log summary statistics
        summary = self.kinetic_db.summary()
        self.logger.info(f"Database loaded: {summary['total_reactions']} reactions, "
                        f"{summary['unique_enzymes']} unique enzymes, "
                        f"{summary['unique_substrates']} unique substrates")
        self.logger.debug(f"  - Reactions with Km: {summary['reactions_with_km']}")
        self.logger.debug(f"  - Reactions with kcat: {summary['reactions_with_kcat']}")
        self.logger.debug(f"  - Reactions with ΔG: {summary['reactions_with_delta_g']}")
        
        return num_reactions
    
    def generate_reactions(self) -> List[GeneratedReaction]:
        """
        Generate all reactions from enzyme-substrate combinations using iterative discovery.
        
        Iteratively discovers reactions where products from one reaction become substrates
        for another reaction, building the complete reaction network.
        
        Returns:
            List[GeneratedReaction]: List of generated reactions
        """
        self.logger.info("=" * 60)
        self.logger.info("REACTION GENERATION")
        self.logger.info("=" * 60)
        
        def reaction_signature(reaction: GeneratedReaction) -> tuple:
            stoich_tuple = tuple(sorted(reaction.stoichiometry.items()))
            reactants_tuple = tuple(sorted(reaction.reactant_labels))
            products_tuple = tuple(sorted(reaction.product_labels))
            return (reactants_tuple, products_tuple, stoich_tuple)

        # Get candidate substrates:
        # - reactive species
        # - exclude solvents
        # - exclude general cofactors as "triggers" (they can still be used as reactants via DB stoichiometry)
        initial_substrates = [
            s
            for s in self.bees_object.species
            if s.reactive and not s.solvent and s.label.lower().strip() not in GENERAL_COFACTORS
        ]
        reactive_species_labels = [
            s.label for s in self.bees_object.species if s.reactive and not s.solvent
        ]
        excluded_general_cofactors = [
            s.label
            for s in self.bees_object.species
            if s.reactive and not s.solvent and s.label.lower().strip() in GENERAL_COFACTORS
        ]

        # Collect reactive enzymes for discovery
        enzymes = [e for e in self.bees_object.enzymes if e.reactive]
        
        substrate_names = [s.label for s in initial_substrates]
        enzyme_names = [e.label for e in enzymes]
        self.logger.info(
            f"Reactive substrates ({len(initial_substrates)}): {', '.join(substrate_names)}"
        )
        self.logger.info(
            f"Reactive enzymes ({len(enzymes)}): {', '.join(enzyme_names)}"
        )
        if excluded_general_cofactors:
            self.logger.info(
                f"General cofactors (available but not triggers): {', '.join(excluded_general_cofactors)}"
            )
        self.logger.info("")
        
        # Build kinetics estimator once (if enabled/configured)
        try:
            if hasattr(self.bees_object, "settings") and getattr(self.bees_object.settings, "estimate_kinetics", False):
                include_sd = getattr(self.bees_object.settings, "kinetics_include_sd", False)
                self.kinetics_estimator = build_estimator(
                    getattr(self.bees_object.settings, "kinetics_estimator", None),
                    include_sd=include_sd,
                )
                if self.kinetics_estimator:
                    self.logger.info(f"Kinetics estimation enabled: {self.kinetics_estimator.name}")
                else:
                    self.logger.info("Kinetics estimation enabled, but no estimator selected.")
        except Exception as e:
            self.logger.warning(f"Failed to initialize kinetics estimator: {e}")
            self.kinetics_estimator = None

        # Initialize available species from input (including ontology equivalents)
        available_species_labels_lc = set()
        provided_species_labels_lc = set()  # Only species explicitly provided or generated
        for s in getattr(self.bees_object, "species", []):
            if getattr(s, "label", None):
                label_lc = s.label.lower().strip()
                provided_species_labels_lc.add(label_lc)
                available_species_labels_lc.update(get_ontology_equivalents(s.label))
        
        # Iterative reaction discovery: products become substrates for next iteration
        existing_signatures = {reaction_signature(r) for r in self.reactions}
        iteration = 0
        max_iterations = 10  # Safety limit to prevent infinite loops
        reaction_count = len(self.reactions)
        
        while iteration < max_iterations:
            iteration += 1
            new_reactions_this_iteration = 0
            substrates_this_iteration = initial_substrates.copy()
            new_products_added = 0
            
            # Add products from previous reactions as potential substrates
            for reaction in self.reactions:
                for product_label in reaction.product_labels:
                    # Check if this product is already a substrate we're tracking
                    product_lc = product_label.lower().strip()
                    if product_lc not in provided_species_labels_lc:
                        # Create a virtual substrate object for products that aren't in input
                        # We'll use it for reaction generation but won't add to actual species list
                        virtual_substrate = SimpleNamespace(
                            label=product_label,
                            reactive=True,
                            solvent=False,
                            smiles=None
                        )
                        # Only add if it's not a general cofactor or coenzyme-like carrier
                        added_to_substrates = False
                        coenzyme_flags = get_coenzyme_like_flags(product_label)
                        is_coenzyme_like = coenzyme_flags["is_coenzyme_like"]
                        if product_lc not in GENERAL_COFACTORS and not is_coenzyme_like:
                            substrates_this_iteration.append(virtual_substrate)
                            new_products_added += 1
                            added_to_substrates = True
                        
                        provided_species_labels_lc.add(product_lc)
                        available_species_labels_lc.update(get_ontology_equivalents(product_label))
                        
                        if coenzyme_flags["is_debug_target"]:
                            pass  # Debug target handling
            
            self.logger.info(f"Iteration {iteration}: Checking {len(substrates_this_iteration)} potential substrate(s) x {len(enzymes)} enzyme(s)...")
            
            for enzyme in enzymes:
                for substrate in substrates_this_iteration:
                    new_reactions = self._generate_reactions(
                        enzyme, 
                        substrate, 
                        available_species_labels_lc,
                        provided_species_labels_lc
                    )
                    for reaction in new_reactions:
                        # Deduplicate against existing reactions
                        signature = reaction_signature(reaction)
                        if signature in existing_signatures:
                            continue
                        self.reactions.append(reaction)
                        existing_signatures.add(signature)
                        reaction_count += 1
                        new_reactions_this_iteration += 1
                        self.logger.info(f"  [{reaction_count}] {reaction}")
                        if reaction.kinetics:
                            kin = reaction.kinetics
                            substrate_label_for_display = reaction.substrate_label
                            enzyme_label_for_display = reaction.enzyme_label
                            # Build kinetics line with substrate-specific labels
                            kin_parts = []
                            km_per = getattr(kin, "km_per_substrate", None)
                            if km_per:
                                for rname, val in km_per.items():
                                    kin_parts.append(f"Km({rname})={val} mM")
                            elif kin.km is not None:
                                kin_parts.append(f"Km({substrate_label_for_display})={kin.km} mM")
                            if kin.kcat is not None:
                                kin_parts.append(f"kcat({enzyme_label_for_display})={kin.kcat} 1/s")
                            if kin.vmax is not None:
                                kin_parts.append(f"Vmax={kin.vmax} mM/s")
                            if kin.delta_g is not None:
                                kin_parts.append(f"dG={kin.delta_g} kJ/mol")
                            # SD info
                            sd_parts = []
                            km_sd_per = getattr(kin, "km_sd_per_substrate", None)
                            if km_sd_per:
                                for rname, val in km_sd_per.items():
                                    sd_parts.append(f"Km_sd({rname})={val:.4g} mM")
                            elif getattr(kin, "km_sd", None) is not None:
                                sd_parts.append(f"Km_sd={kin.km_sd:.4g} mM")
                            if getattr(kin, "kcat_sd", None) is not None:
                                sd_parts.append(f"kcat_sd={kin.kcat_sd:.4g} 1/s")
                            sd_part = f" [SD: {', '.join(sd_parts)}]" if sd_parts else ""
                            self.logger.info(
                                f"      ├─ Kinetics: {', '.join(kin_parts)} "
                                f"(source={kin.source}){sd_part}"
                            )
                        else:
                            self.logger.info(f"      ├─ No database match (will need parameter estimation)")
                        self.logger.info(f"      └─ Type: {reaction.template.template_type} "
                                       f"({reaction.template.ec_class.name})")
            
            if new_reactions_this_iteration == 0:
                self.logger.info(f"No new reactions found in iteration {iteration}. Stopping discovery.")
                break
        
        self.logger.info("")
        self.logger.info(f"Generated {len(self.reactions)} total reaction(s) in {iteration} iteration(s)")
        return self.reactions
    
    def _resolve_smiles_for_compound(
        self,
        compound_label: str,
        kinetic_data,
        substrate_label: str,
    ) -> Optional[str]:
        """Resolve SMILES for a compound from species, compound_smiles, or ontology."""
        # 1. Species label match
        compound_lc = compound_label.lower().strip()
        for species in getattr(self.bees_object, "species", []):
            if hasattr(species, "label") and species.label.lower().strip() == compound_lc:
                if hasattr(species, "smiles") and species.smiles:
                    return species.smiles
                break
        # 2. compound_smiles from kinetic_data (with ontology/stoichiometry)
        compound_smiles = getattr(kinetic_data, "compound_smiles", None)
        if isinstance(compound_smiles, dict):
            s = compound_smiles.get(compound_label)
            if s:
                return s
            equivalents_lc = {str(e).lower().strip() for e in get_ontology_equivalents(compound_label)}
            equivalents_lc.add(compound_lc)
            for name, smiles in compound_smiles.items():
                if not smiles:
                    continue
                if str(name).lower().strip() in equivalents_lc:
                    return smiles
        return None
    
    def _generate_reactions(
        self, 
        enzyme, 
        substrate,
        available_species_labels_lc: Optional[set] = None,
        provided_species_labels_lc: Optional[set] = None
    ) -> List[GeneratedReaction]:
        """
        Generate reactions from enzyme and substrate.
        
        Steps:
        1. Query database for *all* potential reaction matches
        2. For each match:
           a. Generate reaction template from EC number
           b. Merge database products with template inference
           c. Create complete reaction object
        
        Args:
            enzyme: Enzyme object from input
            substrate: Substrate species object from input
            available_species_labels_lc: Set of available species labels (lowercase) including
                                         input species and products from previous reactions
            provided_species_labels_lc: Set of species labels (lowercase) that were actually
                                        provided in input or generated as products. Used for
                                        prioritizing database matches.
            
        Returns:
            List[GeneratedReaction]: List of generated reactions
        """
        enzyme_label = enzyme.label
        substrate_label = substrate.label
        ec_number = enzyme.ecnumber
        
        self.logger.debug(f"Generating reactions for: {substrate_label} + {enzyme_label} ({ec_number})")
        
        # Step 1: Query database for *reaction existence* and kinetic parameters.
        # Support EC number aliasing: try alternative EC numbers if primary fails
        ec_numbers_to_try = get_ec_aliases(ec_number)
        
        all_kinetic_data = []
        if self.kinetic_db:
            # Prepare temperature and pH ranges for filtering
            temp_range = None
            if isinstance(self.bees_object.environment.temperature, tuple):
                temp_range = self.bees_object.environment.temperature
            elif self.bees_object.environment.temperature is not None:
                temp_val = self.bees_object.environment.temperature
                temp_range = (temp_val * 0.95, temp_val * 1.05)
            
            ph_range = None
            if isinstance(self.bees_object.environment.pH, tuple):
                ph_range = self.bees_object.environment.pH
            elif self.bees_object.environment.pH is not None:
                ph_val = self.bees_object.environment.pH
                ph_range = (max(0, ph_val - 0.5), min(14, ph_val + 0.5))
            
            # Try primary EC number and aliases, returning all matches
            for ec_to_try in ec_numbers_to_try:
                matches = self.kinetic_db.query_by_enzyme_substrate(
                    ec_number=ec_to_try,
                    substrate_label=substrate_label,
                    strict=False,
                    temperature_range=temp_range,
                    ph_range=ph_range,
                    available_species_labels_lc=provided_species_labels_lc,
                    return_all=True
                )
                if matches:
                    all_kinetic_data.extend(matches)
        
        if not all_kinetic_data:
            return []

        generated_reactions = []
        for kinetic_data in all_kinetic_data:
            # Step 1b: Optional kinetics estimation (CatPred: kcat=concatenated reactants, Km=per substrate)
            if self.kinetics_estimator is not None:
                try:
                    enzyme_seq = getattr(enzyme, "amino_acid_sequence", None)
                    stoich = kinetic_data.stoichiometry or {}
                    reactants = [k for k, v in stoich.items() if v < 0]
                    reactant_smiles: Dict[str, str] = {}
                    smiles_mode = getattr(
                        getattr(self.bees_object, "settings", None), "smiles_mode", "auto"
                    )
                    for rname in reactants:
                        smi = self._resolve_smiles_for_compound(
                            compound_label=rname,
                            kinetic_data=kinetic_data,
                            substrate_label=substrate_label,
                        )
                        if not smi and smiles_mode == "interactive":
                            print(f"\n[SMILES needed] Compound: \"{rname}\" "
                                  f"(enzyme: {enzyme_label}, EC {ec_number})")
                            user_smiles = input("  Enter SMILES (or press Enter to skip): ").strip()
                            if user_smiles:
                                smi = user_smiles
                                self.logger.info(f"  User provided SMILES for '{rname}': {user_smiles}")
                        if smi:
                            reactant_smiles[rname] = smi

                    if enzyme_seq and reactant_smiles:
                        est = self.kinetics_estimator.estimate(
                            enzyme_sequence=enzyme_seq,
                            reactant_smiles=reactant_smiles,
                            inhibitor_smiles=None,
                        )
                        if est.kcat is not None or est.km is not None or getattr(est, "km_per_substrate", None):
                            kinetic_data.kcat = est.kcat if est.kcat is not None else kinetic_data.kcat
                            kinetic_data.source = est.source
                            if getattr(est, "kcat_sd", None) is not None:
                                kinetic_data.kcat_sd = est.kcat_sd
                            if getattr(est, "ki_sd", None) is not None:
                                kinetic_data.ki_sd = est.ki_sd
                            if getattr(est, "km_per_substrate", None):
                                kinetic_data.km_per_substrate = est.km_per_substrate
                                if getattr(est, "km_sd_per_substrate", None):
                                    kinetic_data.km_sd_per_substrate = est.km_sd_per_substrate
                            if est.km is not None:
                                kinetic_data.km = est.km
                                if getattr(est, "km_sd", None) is not None:
                                    kinetic_data.km_sd = est.km_sd
                        else:
                            self.logger.debug(f"  Kinetics estimation returned no values (source: {est.source})")
                    else:
                        missing = []
                        if not enzyme_seq:
                            missing.append("sequence")
                        if not reactant_smiles:
                            missing.append("SMILES for reactants")
                        self.logger.debug(f"  Skipping kinetics estimation: missing {', '.join(missing)}")
                except Exception as e:
                    self.logger.warning(f"  Kinetics estimation failed: {e}")

            try:
                # Build a template
                template = create_reaction_from_database(
                    substrate=substrate_label,
                    enzyme_label=enzyme_label,
                    ec_number=ec_number,
                    cofactor=None,
                    database_products=None,
                )
                
                # Override topology from database stoichiometry
                if kinetic_data.stoichiometry:
                    template.stoichiometry = dict(kinetic_data.stoichiometry)
                    template.reactants = [s for s, coeff in template.stoichiometry.items() if coeff < 0]
                    template.products = [s for s, coeff in template.stoichiometry.items() if coeff > 0]
                    template.cofactors = [
                        s for s in template.reactants if s.lower().strip() in GENERAL_COFACTORS
                    ]

                # Admission rule: only include the reaction if *all* reactants are available.
                all_available, missing_reactants = validate_reaction_reactants(
                    reactants=template.reactants,
                    available_species_labels_lc=available_species_labels_lc,
                    enzyme_label=enzyme_label,
                    ec_number=ec_number
                )
                
                if not all_available:
                    continue
                
                # Build complete reaction
                rate_law = self._determine_rate_law(template, kinetic_data)
                
                reaction = GeneratedReaction(
                    enzyme_label=enzyme_label,
                    substrate_label=substrate_label,
                    ec_number=ec_number,
                    template=template,
                    kinetics=kinetic_data,
                    reactant_labels=template.reactants,
                    product_labels=template.products,
                    stoichiometry=template.stoichiometry,
                    rate_law=rate_law
                )
                generated_reactions.append(reaction)
                
            except Exception as e:
                self.logger.warning(f"Failed to create reaction for {substrate_label}: {e}")
                continue
        
        return generated_reactions
    
    def _determine_rate_law(
        self, 
        template: ReactionTemplate, 
        kinetics: Optional[KineticData]
    ) -> Optional[str]:
        """
        Determine appropriate rate law for reaction.
        
        Only returns a rate law if kinetic data is available from the database.
        If no database match, returns None.
        
        Args:
            template (ReactionTemplate): Reaction template
            kinetics (KineticData): Kinetic parameters (or None)
            
        Returns:
            str or None: Rate law type if kinetics available, None otherwise
        """
        # Assign rate law if we have kinetics (DB or estimated)
        has_km = kinetics.km is not None or (getattr(kinetics, "km_per_substrate", None) and kinetics.km_per_substrate)
        if kinetics and (has_km or kinetics.kcat is not None):
            if template.reversible:
                return "Reversible-MM"
            return "Michaelis-Menten"
        
        # No rate law if not in database
        return None
    
    def export_reactions_summary(self, filename="reactions_summary.txt") -> str:
        """
        Export human-readable summary of generated reactions.
        
        Args:
            filename (str): Output filename
            
        Returns:
            str: Path to output file
        """
        output_path = os.path.join(self.output_directory, filename)
        
        with open(output_path, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("BEES GENERATED REACTIONS SUMMARY\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"Total Reactions: {len(self.reactions)}\n")
            f.write(f"Project: {self.bees_object.project}\n")
            f.write(f"Database: {self.bees_object.database.name}\n")
            f.write(f"Temperature: {self.bees_object.environment.temperature} K\n")
            f.write(f"pH: {self.bees_object.environment.pH}\n")
            f.write("\n" + "=" * 80 + "\n\n")

            # Coenzyme participation note
            f.write("COENZYME PARTICIPATION NOTE\n")
            f.write("-" * 80 + "\n")
            f.write(
                "Energy carriers and redox coenzymes (e.g., ATP, NADH) are treated as\n"
                "secondary reactants. They can participate in reaction stoichiometry, but\n"
                "they do not trigger new reaction discovery on their own.\n"
            )
            f.write(
                "Common coenzymes/cofactors considered in this run: "
                f"{', '.join(sorted(GENERAL_COFACTORS))}\n"
            )
            f.write("\n" + "=" * 80 + "\n\n")
            
            # Add network schema section
            f.write("REACTION NETWORK SCHEMA\n")
            f.write("-" * 80 + "\n")
            f.write("This section shows the connectivity of species in the reaction network.\n")
            f.write("Species are connected through reactions where products of one reaction\n")
            f.write("may be reactants in another.\n\n")
            
            # Collect all unique species from reactions
            all_species = set()
            species_as_reactant = {}  # species -> list of reaction indices
            species_as_product = {}   # species -> list of reaction indices
            
            for i, rxn in enumerate(self.reactions, 1):
                for reactant in rxn.reactant_labels:
                    all_species.add(reactant)
                    if reactant not in species_as_reactant:
                        species_as_reactant[reactant] = []
                    species_as_reactant[reactant].append(i)
                
                for product in rxn.product_labels:
                    all_species.add(product)
                    if product not in species_as_product:
                        species_as_product[product] = []
                    species_as_product[product].append(i)
            
            # Sort species for consistent output
            sorted_species = sorted(all_species)
            
            f.write(f"Total Unique Species in Network: {len(sorted_species)}\n\n")
            
            # Show species connectivity
            f.write("Species Connectivity:\n")
            for species in sorted_species:
                f.write(f"  {species}:\n")
                if species in species_as_reactant:
                    rxns_consuming = species_as_reactant[species]
                    f.write(f"    Consumed in reaction(s): {', '.join(f'R{i}' for i in rxns_consuming)}\n")
                if species in species_as_product:
                    rxns_producing = species_as_product[species]
                    f.write(f"    Produced in reaction(s): {', '.join(f'R{i}' for i in rxns_producing)}\n")
                if species not in species_as_reactant and species not in species_as_product:
                    f.write(f"    (Only in stoichiometry, not in main reaction equation)\n")
            
            f.write("\n" + "-" * 80 + "\n\n")
            
            for i, rxn in enumerate(self.reactions, 1):
                f.write(f"Reaction {i}:\n")
                f.write(f"  Enzyme: {rxn.enzyme_label}\n")
                f.write(f"  EC Number: {rxn.ec_number or 'N/A'}\n")
                f.write(f"  EC Class: {rxn.template.ec_class.name}\n")
                f.write(f"  Type: {rxn.template.template_type}\n")
                f.write(f"  Equation: {' + '.join(rxn.reactant_labels)} → {' + '.join(rxn.product_labels)}\n")
                f.write(f"  Stoichiometry: {rxn.stoichiometry}\n")
                
                # Only show rate law if kinetics data is available from database
                if rxn.kinetics:
                    f.write(f"  Rate Law: {rxn.rate_law}\n")
                
                if rxn.kinetics:
                    substrate = rxn.substrate_label
                    enzyme = rxn.enzyme_label
                    kin = rxn.kinetics
                    f.write(f"  Kinetic Parameters:\n")
                    km_per = getattr(kin, "km_per_substrate", None)
                    if km_per:
                        for rname, val in km_per.items():
                            f.write(f"    Km({rname}) = {val} mM\n")
                    elif kin.km is not None:
                        f.write(f"    Km({substrate}) = {kin.km} mM\n")
                    if rxn.kinetics.kcat is not None:
                        f.write(f"    kcat({enzyme}) = {rxn.kinetics.kcat} 1/s\n")
                    if rxn.kinetics.vmax is not None:
                        f.write(f"    Vmax = {rxn.kinetics.vmax} mM/s\n")
                    if rxn.kinetics.delta_g is not None:
                        f.write(f"    dG = {rxn.kinetics.delta_g} kJ/mol\n")
                    if rxn.kinetics.temperature is not None:
                        f.write(f"    Temperature = {rxn.kinetics.temperature} K\n")
                    if rxn.kinetics.ph is not None:
                        f.write(f"    pH = {rxn.kinetics.ph}\n")
                    f.write(f"    Source: {kin.source}\n")
                    km_sd_per = getattr(kin, "km_sd_per_substrate", None)
                    if km_sd_per:
                        for rname, val in km_sd_per.items():
                            f.write(f"    Km({rname}) SD = {val:.4g} mM\n")
                    elif getattr(kin, "km_sd", None) is not None:
                        f.write(f"    Km({substrate}) SD = {kin.km_sd:.4g} mM\n")
                    if getattr(kin, "kcat_sd", None) is not None:
                        f.write(f"    kcat({enzyme}) SD = {kin.kcat_sd:.4g} 1/s\n")
                    if getattr(kin, "ki_sd", None) is not None:
                        f.write(f"    Ki SD = {kin.ki_sd:.4g} mM\n")
                else:
                    f.write(f"  Kinetic Parameters: NOT FOUND IN DATABASE\n")
                    f.write(f"    (Parameters will need to be estimated)\n")
                
                if rxn.template.cofactors:
                    f.write(f"  Cofactors: {', '.join(rxn.template.cofactors)}\n")
                
                f.write("\n")
        
        self.logger.info(f"Exported reactions summary to {output_path}")
        return output_path
    
    def get_reactions(self) -> List[GeneratedReaction]:
        """
        Get list of generated reactions.
        
        Returns:
            List[GeneratedReaction]: Generated reactions
        """
        return self.reactions

