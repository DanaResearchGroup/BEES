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

from typing import List, Dict, Optional, Any
from dataclasses import dataclass
import os

from bees.reaction_template import (
    create_reaction_from_database,
    ReactionTemplate,
    determine_template_from_ec
)
from db.kinetic_database import KineticDatabase, KineticData


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
        rate_law (str): Rate law type (e.g., "Michaelis-Menten")
    """
    enzyme_label: str
    substrate_label: str
    template: ReactionTemplate
    kinetics: Optional[KineticData]
    reactant_labels: List[str]
    product_labels: List[str]
    stoichiometry: Dict[str, int]
    rate_law: str = "Michaelis-Menten"  # Default
    
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
        kinetic_db (KineticDatabase): Loaded kinetic parameter database
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
    
    def load_kinetic_database(self, db_path: str) -> int:
        """
        Load kinetic database from CSV file.
        
        Args:
            db_path (str): Path to kinetic database CSV
            
        Returns:
            int: Number of reactions loaded
            
        Raises:
            FileNotFoundError: If database file not found
            ValueError: If CSV format is invalid
        """
        self.logger.info(f"Loading kinetic database from {db_path}")
        self.kinetic_db = KineticDatabase(logger=self.logger)
        num_reactions = self.kinetic_db.load_from_csv(db_path)
        
        # Log summary statistics
        summary = self.kinetic_db.summary()
        self.logger.info(f"Database loaded: {summary['total_reactions']} reactions, "
                        f"{summary['unique_enzymes']} unique enzymes, "
                        f"{summary['unique_substrates']} unique substrates")
        self.logger.info(f"  - Reactions with Km: {summary['reactions_with_km']}")
        self.logger.info(f"  - Reactions with kcat: {summary['reactions_with_kcat']}")
        self.logger.info(f"  - Reactions with ΔG: {summary['reactions_with_delta_g']}")
        
        return num_reactions
    
    def generate_reactions(self) -> List[GeneratedReaction]:
        """
        Generate all reactions from enzyme-substrate combinations.
        
        Creates a reaction for each reactive enzyme × reactive substrate pair
        by querying the database and generating templates.
        
        Returns:
            List[GeneratedReaction]: List of generated reactions
        """
        self.logger.info("=" * 60)
        self.logger.info("REACTION GENERATION")
        self.logger.info("=" * 60)
        
        # Get reactive species (substrates) - exclude solvents
        substrates = [s for s in self.bees_object.species if s.reactive and not s.solvent]
        enzymes = [e for e in self.bees_object.enzymes if e.reactive]
        
        self.logger.info(f"Found {len(substrates)} reactive substrate(s) and {len(enzymes)} reactive enzyme(s)")
        self.logger.info(f"Generating reactions for {len(enzymes) * len(substrates)} enzyme-substrate pair(s)...")
        self.logger.info("")
        
        # Generate reactions for each enzyme-substrate pair
        reaction_count = 0
        for enzyme in enzymes:
            for substrate in substrates:
                reaction = self._generate_single_reaction(enzyme, substrate)
                if reaction:
                    self.reactions.append(reaction)
                    reaction_count += 1
                    self.logger.info(f"  [{reaction_count}] {reaction}")
                    if reaction.kinetics:
                        self.logger.info(f"      ├─ Database match: Km={reaction.kinetics.km} mM, "
                                       f"kcat={reaction.kinetics.kcat} 1/s")
                    else:
                        self.logger.info(f"      ├─ No database match (will need parameter estimation)")
                    self.logger.info(f"      └─ Type: {reaction.template.template_type} "
                                   f"({reaction.template.ec_class.name})")
        
        self.logger.info("")
        self.logger.info(f"Generated {len(self.reactions)} total reaction(s)")
        return self.reactions
    
    def _generate_single_reaction(
        self, 
        enzyme, 
        substrate
    ) -> Optional[GeneratedReaction]:
        """
        Generate a single reaction from enzyme and substrate.
        
        Steps:
        1. Query database for kinetic parameters
        2. Generate reaction template from EC number
        3. Merge database products with template inference
        4. Create complete reaction object
        
        Args:
            enzyme: Enzyme object from input
            substrate: Substrate species object from input
            
        Returns:
            GeneratedReaction or None: Generated reaction if successful
        """
        enzyme_label = enzyme.label
        substrate_label = substrate.label
        ec_number = enzyme.ecnumber
        
        self.logger.debug(f"Generating reaction: {substrate_label} + {enzyme_label} ({ec_number})")
        
        # Step 1: Query database for kinetic parameters
        kinetic_data = None
        if self.kinetic_db:
            # Prepare temperature and pH ranges for filtering
            temp_range = None
            if isinstance(self.bees_object.environment.temperature, tuple):
                temp_range = self.bees_object.environment.temperature
            elif self.bees_object.environment.temperature is not None:
                # Single value - create a small range around it (5% tolerance)
                temp_val = self.bees_object.environment.temperature
                temp_range = (temp_val * 0.95, temp_val * 1.05)
            
            ph_range = None
            if isinstance(self.bees_object.environment.pH, tuple):
                ph_range = self.bees_object.environment.pH
            elif self.bees_object.environment.pH is not None:
                # Single value - create a small range around it (0.5 pH units)
                ph_val = self.bees_object.environment.pH
                ph_range = (max(0, ph_val - 0.5), min(14, ph_val + 0.5))
            
            kinetic_data = self.kinetic_db.query_by_enzyme_substrate(
                ec_number=ec_number,
                substrate_label=substrate_label,
                strict=False,  # Don't require exact cofactor match
                temperature_range=temp_range,
                ph_range=ph_range
            )
            
            if kinetic_data:
                self.logger.debug(f"  Found kinetic data in database: {kinetic_data}")
        
        # Step 2: Generate reaction template from EC number
        cofactor = kinetic_data.cofactors if kinetic_data else None
        database_products = kinetic_data.products if kinetic_data else None
        
        try:
            template = create_reaction_from_database(
                substrate=substrate_label,
                enzyme_label=enzyme_label,
                ec_number=ec_number,
                cofactor=cofactor,
                database_products=database_products
            )
        except Exception as e:
            self.logger.warning(f"Failed to create reaction template for {substrate_label} + "
                              f"{enzyme_label}: {e}")
            return None
        
        # Step 3: Build complete reaction
        reactants = template.reactants
        products = template.products
        stoichiometry = template.stoichiometry
        
        # Determine rate law based on EC class or available kinetic data
        rate_law = self._determine_rate_law(template, kinetic_data)
        
        reaction = GeneratedReaction(
            enzyme_label=enzyme_label,
            substrate_label=substrate_label,
            template=template,
            kinetics=kinetic_data,
            reactant_labels=reactants,
            product_labels=products,
            stoichiometry=stoichiometry,
            rate_law=rate_law
        )
        
        return reaction
    
    def _determine_rate_law(
        self, 
        template: ReactionTemplate, 
        kinetics: Optional[KineticData]
    ) -> str:
        """
        Determine appropriate rate law for reaction.
        
        Args:
            template (ReactionTemplate): Reaction template
            kinetics (KineticData): Kinetic parameters (or None)
            
        Returns:
            str: Rate law type
        """
        # If we have Km/Vmax or kcat, use Michaelis-Menten
        if kinetics and (kinetics.km is not None or kinetics.kcat is not None):
            if template.reversible:
                return "Reversible-MM"
            return "Michaelis-Menten"
        
        # For reversible reactions without kinetics, use reversible mass action
        if template.reversible:
            return "Reversible-MassAction"
        
        # Default to Michaelis-Menten (parameters will need estimation)
        return "Michaelis-Menten"
    
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
                f.write(f"  EC Class: {rxn.template.ec_class.name}\n")
                f.write(f"  Type: {rxn.template.template_type}\n")
                f.write(f"  Description: {rxn.template.description}\n")
                f.write(f"  Equation: {' + '.join(rxn.reactant_labels)} → {' + '.join(rxn.product_labels)}\n")
                f.write(f"  Stoichiometry: {rxn.stoichiometry}\n")
                f.write(f"  Rate Law: {rxn.rate_law}\n")
                
                if rxn.kinetics:
                    f.write(f"  Kinetic Parameters (from database):\n")
                    if rxn.kinetics.km is not None:
                        f.write(f"    Km = {rxn.kinetics.km} mM\n")
                    if rxn.kinetics.kcat is not None:
                        f.write(f"    kcat = {rxn.kinetics.kcat} 1/s\n")
                    if rxn.kinetics.vmax is not None:
                        f.write(f"    Vmax = {rxn.kinetics.vmax} mM/s\n")
                    if rxn.kinetics.delta_g is not None:
                        f.write(f"    ΔG = {rxn.kinetics.delta_g} kJ/mol\n")
                    if rxn.kinetics.temperature is not None:
                        f.write(f"    Temperature = {rxn.kinetics.temperature} K\n")
                    if rxn.kinetics.ph is not None:
                        f.write(f"    pH = {rxn.kinetics.ph}\n")
                    f.write(f"    Source: {rxn.kinetics.source}\n")
                else:
                    f.write(f"  Kinetic Parameters: NOT FOUND IN DATABASE\n")
                    f.write(f"    (Parameters will need to be estimated)\n")
                
                f.write(f"  Reversible: {rxn.template.reversible}\n")
                
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

