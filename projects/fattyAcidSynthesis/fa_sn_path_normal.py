#!/usr/bin/env python3

"""
Fatty Acid Synthesis (FAS) Normal Pathway Module
------------------------------------------------
Defines the "Golden Path" for de novo lipogenesis in cells.

This module loads the pathway definition from a YAML file (fa_sn_path_normal.yaml),


Pathway Overview:
1. Citrate → Acetyl-CoA + Oxaloacetate (ATP citrate lyase)
2. Acetyl-CoA → Malonyl-CoA (Acetyl-CoA carboxylase)
3. Acetyl-CoA + Malonyl-CoA → Butyryl-CoA (C4) - Initial condensation
4-9. Iterative elongation: C4 → C6 → C8 → C10 → C12 → C14 → C16 (Palmitate)

The pathway includes all elongation steps from C2 to C16, providing a complete
model of the fatty acid synthesis process.

The pathway represents the core fatty acid synthesis route that will be used
as the baseline for perturbations (inhibitors, mutations, etc.).

Data Source: Hepatokin1 SBML Model (Berndt et al., 2019)
"""

import logging
import os
from typing import List, Optional, Dict, Any
from dataclasses import dataclass, field
import yaml

# Logger
logger = logging.getLogger('BEES')


@dataclass
class PathwayReaction:
    """
    Represents a single reaction in the fatty acid synthesis pathway.
    
    Attributes:
        step_number (int): Sequential step number in the pathway
        enzyme_name (str): Name of the enzyme catalyzing this reaction
        ec_number (str): EC number (e.g., "EC 2.3.3.8")
        reactants (List[str]): List of reactant metabolite names (consumed in reaction)
        products (List[str]): List of product metabolite names (produced in reaction)
        stoichiometry (Dict[str, int]): Stoichiometric coefficients (key=metabolite name, value=coefficient)
                                        Negative for reactants, positive for products
        cofactors (List[str]): Required cofactors/coenzymes (e.g., ["Mg2+", "Biotin"])
                               Note: Cofactors are required but may not be consumed (vs reactants)
        reversible (bool): Whether the reaction is reversible
        description (str): Human-readable description of the reaction
        default_km (Optional[float]): Default Km value in mM (primary substrate, if known)
        km_values (Dict[str, float]): Km values for multiple substrates in mM (key=substrate name, value=Km)
        default_kcat (Optional[float]): Default kcat value in 1/s (if known)
        default_ki (Optional[float]): Default Ki value in mM for inhibitors (if known)
        parameter_source (Optional[str]): Source/reference for kinetic parameters (e.g., "literature", "BRENDA", "experimental")
        references (List[str]): List of reference citations for kinetic parameters
        metabolite_smiles (Dict[str, str]): SMILES strings for metabolites (key=metabolite name, value=SMILES)
        rate_law_type (Optional[str]): Type of rate law (e.g., "Michaelis-Menten", "Hill", "Custom")
        rate_law_equation (Optional[str]): Human-readable rate law equation as string
        rate_law_mathml (Optional[str]): Original MathML representation of rate law from SBML
    """
    step_number: int
    enzyme_name: str
    ec_number: str
    reactants: List[str] = field(default_factory=list)
    products: List[str] = field(default_factory=list)
    stoichiometry: Dict[str, int] = field(default_factory=dict)  # Metabolite -> coefficient
    cofactors: List[str] = field(default_factory=list)
    reversible: bool = False
    description: str = ""
    default_km: Optional[float] = None  # mM (primary substrate Km)
    km_values: Dict[str, float] = field(default_factory=dict)  # Substrate name -> Km in mM (for multiple substrates)
    default_kcat: Optional[float] = None  # 1/s
    default_ki: Optional[float] = None  # mM
    parameter_source: Optional[str] = None  # Source/reference for parameters
    references: List[str] = field(default_factory=list)  # List of reference citations
    metabolite_smiles: Dict[str, str] = field(default_factory=dict)  # Metabolite name -> SMILES
    rate_law_type: Optional[str] = None  # Type of rate law (e.g., "Michaelis-Menten", "Hill", "Custom")
    rate_law_equation: Optional[str] = None  # Human-readable rate law equation
    rate_law_mathml: Optional[str] = None  # Original MathML from SBML
    
    def __repr__(self):
        reactants_str = " + ".join(self.reactants)
        products_str = " + ".join(self.products)
        return f"Step {self.step_number}: {reactants_str} → {products_str} ({self.enzyme_name})"


@dataclass
class FASPathway:
    """
    Container for the complete fatty acid synthesis pathway.
    
    This represents the "Golden Path" or normal (unperturbed) state of
    de novo lipogenesis in cancer cells.
    
    Attributes:
        name (str): Pathway name
        reactions (List[PathwayReaction]): Ordered list of pathway reactions
        entry_metabolite (str): Starting metabolite (Citrate)
        exit_metabolite (str): Final product (Palmitate)
        intermediate_metabolites (List[str]): Intermediate metabolites in the pathway
    """
    name: str
    reactions: List[PathwayReaction] = field(default_factory=list)
    entry_metabolite: str = "Citrate"
    exit_metabolite: str = "Palmitate"
    intermediate_metabolites: List[str] = field(default_factory=list)
    
    def get_reaction_by_step(self, step_number: int) -> Optional[PathwayReaction]:
        """Get reaction by step number."""
        for rxn in self.reactions:
            if rxn.step_number == step_number:
                return rxn
        return None
    
    def get_reaction_by_enzyme(self, enzyme_name: str) -> Optional[PathwayReaction]:
        """Get reaction by enzyme name."""
        for rxn in self.reactions:
            if rxn.enzyme_name.lower() == enzyme_name.lower():
                return rxn
        return None
    
    def get_reaction_by_ec(self, ec_number: str) -> Optional[PathwayReaction]:
        """Get reaction by EC number."""
        for rxn in self.reactions:
            if rxn.ec_number == ec_number:
                return rxn
        return None
    
    def get_reaction_by_metabolite(self, metabolite: str) -> List[PathwayReaction]:
        """Get all reactions involving a specific metabolite (as reactant or product)."""
        matching = []
        for rxn in self.reactions:
            if metabolite in rxn.reactants or metabolite in rxn.products:
                matching.append(rxn)
        return matching
    
    def validate_pathway(self) -> bool:
        """
        Validate that the pathway is connected and makes sense.
        
        Returns:
            bool: True if pathway is valid
        """
        if not self.reactions:
            logger.warning("Pathway has no reactions")
            return False
        
        # Check that entry metabolite appears as reactant in first reaction
        first_rxn = self.reactions[0]
        if self.entry_metabolite not in first_rxn.reactants:
            logger.warning(f"Entry metabolite {self.entry_metabolite} not found in first reaction")
            return False
        
        # Check that exit metabolite appears as product in last reaction
        last_rxn = self.reactions[-1]
        if self.exit_metabolite not in last_rxn.products:
            logger.warning(f"Exit metabolite {self.exit_metabolite} not found in last reaction")
            return False
        
        # Check connectivity (each product should be reactant in next step or be exit metabolite)
        for i in range(len(self.reactions) - 1):
            current_rxn = self.reactions[i]
            next_rxn = self.reactions[i + 1]
            
            # At least one product of current should be reactant of next
            common_metabolites = set(current_rxn.products) & set(next_rxn.reactants)
            if not common_metabolites and self.exit_metabolite not in current_rxn.products:
                logger.warning(
                    f"Pathway discontinuity: Step {current_rxn.step_number} products "
                    f"({current_rxn.products}) don't connect to Step {next_rxn.step_number} "
                    f"reactants ({next_rxn.reactants})"
                )
                return False
        
        return True
    
    def __repr__(self):
        return f"FASPathway(name={self.name}, {len(self.reactions)} reactions, {self.entry_metabolite} → {self.exit_metabolite})"


def _load_pathway_from_yaml(yaml_path: Optional[str] = None) -> Dict[str, Any]:
    """
    Load pathway definition from YAML file.
    
    Args:
        yaml_path: Path to YAML file. If None, uses default location.
        
    Returns:
        dict: Parsed YAML data
        
    Raises:
        FileNotFoundError: If YAML file not found
        yaml.YAMLError: If YAML parsing fails
    """
    if yaml_path is None:
        # Default to YAML file in same directory as this module
        module_dir = os.path.dirname(os.path.abspath(__file__))
        yaml_path = os.path.join(module_dir, "fa_sn_path_normal.yaml")
    
    if not os.path.exists(yaml_path):
        raise FileNotFoundError(
            f"Pathway YAML file not found: {yaml_path}. "
            f"Please ensure the file exists or provide a valid path."
        )
    
    with open(yaml_path, 'r') as f:
        data = yaml.safe_load(f)
    
    logger.info(f"Loaded pathway definition from {yaml_path}")
    return data


def _create_reaction_from_yaml(rxn_data: Dict[str, Any]) -> PathwayReaction:
    """
    Create a PathwayReaction from YAML data.
    
    Args:
        rxn_data: Dictionary containing reaction data from YAML
        
    Returns:
        PathwayReaction: Created reaction object
    """
    # Extract kinetic parameters
    kin_params = rxn_data.get('kinetic_parameters', {})
    km_values = kin_params.get('km_values', {})
    default_km = kin_params.get('default_km')
    default_kcat = kin_params.get('default_kcat')
    default_ki = kin_params.get('default_ki')
    
    # Extract rate law information
    rate_law_type = rxn_data.get('rate_law_type')
    rate_law_equation = rxn_data.get('rate_law_equation')
    rate_law_mathml = rxn_data.get('rate_law_mathml')
    
    # Extract other fields
    metabolite_smiles = rxn_data.get('metabolite_smiles', {})
    
    return PathwayReaction(
        step_number=rxn_data['step_number'],
        enzyme_name=rxn_data['enzyme_name'],
        ec_number=rxn_data['ec_number'],
        reactants=rxn_data.get('reactants', []),
        products=rxn_data.get('products', []),
        stoichiometry=rxn_data.get('stoichiometry', {}),
        cofactors=rxn_data.get('cofactors', []),
        reversible=rxn_data.get('reversible', False),
        description=rxn_data.get('description', ''),
        default_km=default_km,
        km_values=km_values,
        default_kcat=default_kcat,
        default_ki=default_ki,
        parameter_source=rxn_data.get('parameter_source'),
        references=rxn_data.get('references', []),
        metabolite_smiles=metabolite_smiles,
        rate_law_type=rate_law_type,
        rate_law_equation=rate_law_equation,
        rate_law_mathml=rate_law_mathml
    )


def create_fa_sn_path_normal(yaml_path: Optional[str] = None) -> FASPathway:
    """
    Create the normal (unperturbed) fatty acid synthesis pathway.
    
    This function loads the pathway definition from a YAML file, which contains
    the complete FAS pathway including all elongation steps from C2 to C16.
    
    The pathway represents the "Golden Path" for de novo lipogenesis:
    Citrate → Acetyl-CoA → Malonyl-CoA → [C4 → C6 → C8 → C10 → C12 → C14 → C16] → Palmitate
    
    Args:
        yaml_path: Optional path to YAML file. If None, uses default location.
        
    Returns:
        FASPathway: Complete pathway definition loaded from YAML
        
    Raises:
        FileNotFoundError: If YAML file not found
        ValueError: If YAML data is invalid
    """
    # Load YAML data
    data = _load_pathway_from_yaml(yaml_path)
    
    pathway_info = data.get('pathway', {})
    reactions_data = data.get('reactions', [])
    intermediate_metabolites = data.get('intermediate_metabolites', [])
    
    if not reactions_data:
        raise ValueError("No reactions found in YAML file")
    
    # Create PathwayReaction objects from YAML data
    reactions = []
    for rxn_data in reactions_data:
        reaction = _create_reaction_from_yaml(rxn_data)
        reactions.append(reaction)
    
    # Create pathway object
    pathway = FASPathway(
        name=pathway_info.get('name', 'Fatty Acid Synthesis Normal Path'),
        reactions=reactions,
        entry_metabolite=pathway_info.get('entry_metabolite', 'Citrate'),
        exit_metabolite=pathway_info.get('exit_metabolite', 'Palmitate'),
        intermediate_metabolites=intermediate_metabolites
    )
    
    # Validate pathway
    if not pathway.validate_pathway():
        logger.warning("Pathway validation failed, but continuing anyway")
    
    logger.info(f"Created FAS pathway from YAML: {pathway}")
    return pathway


def get_fa_sn_path_normal() -> FASPathway:
    """
    Get the normal fatty acid synthesis pathway.
    
    Convenience function that returns the standard pathway definition.
    This is the main entry point for accessing the Golden Path.
    
    Returns:
        FASPathway: The normal (unperturbed) pathway
    """
    return create_fa_sn_path_normal()


# Module-level pathway instance (created on import)
_FA_SN_PATH_NORMAL: Optional[FASPathway] = None


def get_pathway() -> FASPathway:
    """
    Get the cached pathway instance.
    
    Returns:
        FASPathway: The normal pathway (cached after first call)
    """
    global _FA_SN_PATH_NORMAL
    if _FA_SN_PATH_NORMAL is None:
        _FA_SN_PATH_NORMAL = create_fa_sn_path_normal()
    return _FA_SN_PATH_NORMAL
