#!/usr/bin/env python3

"""
Reaction Template Module
------------------------
Determines reaction types and templates based on enzyme EC numbers.
Provides simple product inference for enzymatic reactions.

EC Number Classification:
- EC 1.x.x.x: Oxidoreductases (electron transfer)
- EC 2.x.x.x: Transferases (group transfer)
- EC 3.x.x.x: Hydrolases (hydrolysis)
- EC 4.x.x.x: Lyases (addition/removal without hydrolysis)
- EC 5.x.x.x: Isomerases (intramolecular rearrangement)
- EC 6.x.x.x: Ligases (bond formation with ATP)
- EC 7.x.x.x: Translocases (movement across membranes)


#TODO Future possible enhancements:
- SMARTS-based pattern matching for substrate structure analysis
- More detailed sub-class templates
- Cofactor-specific product determination
"""

import logging
import re
from typing import List, Optional, Dict, Any, Tuple
from dataclasses import dataclass, field
from enum import Enum

# Logger
logger = logging.getLogger('BEES')


class ECClass(Enum):
    """EC Number main classes."""
    OXIDOREDUCTASE = 1  # Electron transfer reactions
    TRANSFERASE = 2      # Group transfer reactions
    HYDROLASE = 3        # Hydrolysis reactions
    LYASE = 4            # Addition/elimination reactions
    ISOMERASE = 5        # Isomerization reactions
    LIGASE = 6           # Bond formation with ATP cleavage
    TRANSLOCASE = 7      # Translocation reactions
    UNKNOWN = 0


@dataclass
class ReactionTemplate:
    """
    Reaction template describing the type and stoichiometry of a biochemical reaction.
    
    Attributes:
        template_type (str): Type of reaction (e.g., "phosphorylation", "hydrolysis")
        ec_class (ECClass): Main EC class
        reactants (List[str]): List of reactant labels
        products (List[str]): List of product labels (may be inferred)
        stoichiometry (Dict[str, int]): Stoichiometric coefficients (negative=consumed, positive=produced)
        cofactors (List[str]): Required cofactors
        description (str): Human-readable description
        reversible (bool): Whether reaction is reversible
    """
    template_type: str
    ec_class: ECClass
    reactants: List[str] = field(default_factory=list)
    products: List[str] = field(default_factory=list)
    stoichiometry: Dict[str, int] = field(default_factory=dict)
    cofactors: List[str] = field(default_factory=list)
    description: str = ""
    reversible: bool = False
    
    def __repr__(self):
        return (f"ReactionTemplate(type={self.template_type}, class={self.ec_class.name}, "
                f"reactants={self.reactants} → products={self.products})")


def parse_ec_number(ec_number: str) -> Tuple[int, int, int, int]:
    """
    Parse EC number string into components.
    
    Args:
        ec_number (str): EC number in format "EC 1.2.3.4" or "1.2.3.4"
        
    Returns:
        Tuple[int, int, int, int]: (class, subclass, sub-subclass, serial)
        
    Raises:
        ValueError: If EC number format is invalid
    """
    # Remove "EC " prefix if present
    ec_str = ec_number.strip().upper()
    if ec_str.startswith("EC "):
        ec_str = ec_str[3:].strip()
    
    # Parse components
    match = re.match(r'^(\d+)\.(\d+)\.(\d+)\.(\d+)$', ec_str)
    if not match:
        raise ValueError(f"Invalid EC number format: {ec_number}. Expected 'EC X.X.X.X' or 'X.X.X.X'")
    
    return tuple(int(x) for x in match.groups())


def get_ec_class(ec_number: str) -> ECClass:
    """
    Get the main EC class from an EC number.
    
    Args:
        ec_number (str): EC number
        
    Returns:
        ECClass: Main enzyme class
    """
    try:
        main_class, _, _, _ = parse_ec_number(ec_number)
        return ECClass(main_class)
    except (ValueError, KeyError):
        logger.warning(f"Could not determine EC class for {ec_number}")
        return ECClass.UNKNOWN


def determine_template_from_ec(ec_number: str, substrate: Optional[str] = None) -> ReactionTemplate:
    """
    Determine reaction template based on EC number.
    
    Uses EC classification to infer reaction type. More specific templates
    can be determined by analyzing EC sub-classes.
    
    Args:
        ec_number (str): Enzyme EC number
        substrate (str, optional): Substrate name for more specific inference
        
    Returns:
        ReactionTemplate: Inferred reaction template
    """
    try:
        main_class, sub_class, _, _ = parse_ec_number(ec_number)
        ec_class = ECClass(main_class)
    except (ValueError, KeyError):
        logger.warning(f"Invalid EC number {ec_number}, using generic template")
        return ReactionTemplate(
            template_type="generic",
            ec_class=ECClass.UNKNOWN,
            description="Generic enzymatic reaction"
        )
    
    # Determine template based on EC class
    if ec_class == ECClass.OXIDOREDUCTASE:
        # EC 1: Oxidoreductases
        if sub_class == 1:
            template_type = "oxidation_CH-OH"
            description = "Oxidation of CH-OH group"
        elif sub_class == 2:
            template_type = "oxidation_C=O"
            description = "Oxidation of aldehyde or ketone"
        elif sub_class == 3:
            template_type = "oxidation_CH-CH"
            description = "Oxidation of CH-CH group"
        else:
            template_type = "oxidoreduction"
            description = "Oxidoreduction reaction"
        
        return ReactionTemplate(
            template_type=template_type,
            ec_class=ec_class,
            description=description,
            reversible=True
        )
    
    elif ec_class == ECClass.TRANSFERASE:
        # EC 2: Transferases
        if sub_class == 7:
            # Phosphotransferases (kinases)
            template_type = "phosphorylation"
            description = "Phosphate group transfer from ATP"
            cofactors = ["ATP", "Mg2+"]
        elif sub_class == 3:
            template_type = "acyl_transfer"
            description = "Acyl group transfer"
            cofactors = []
        elif sub_class == 4:
            template_type = "glycosyl_transfer"
            description = "Glycosyl group transfer"
            cofactors = []
        else:
            template_type = "group_transfer"
            description = "Generic group transfer"
            cofactors = []
        
        return ReactionTemplate(
            template_type=template_type,
            ec_class=ec_class,
            cofactors=cofactors,
            description=description,
            reversible=False
        )
    
    elif ec_class == ECClass.HYDROLASE:
        # EC 3: Hydrolases
        if sub_class == 1:
            template_type = "ester_hydrolysis"
            description = "Hydrolysis of ester bonds"
        elif sub_class == 2:
            template_type = "glycoside_hydrolysis"
            description = "Hydrolysis of glycosidic bonds"
        elif sub_class == 4:
            template_type = "peptide_hydrolysis"
            description = "Hydrolysis of peptide bonds"
        else:
            template_type = "hydrolysis"
            description = "Hydrolysis reaction"
        
        return ReactionTemplate(
            template_type=template_type,
            ec_class=ec_class,
            description=description,
            reversible=False
        )
    
    elif ec_class == ECClass.LYASE:
        # EC 4: Lyases
        template_type = "elimination"
        description = "Addition or elimination reaction"
        return ReactionTemplate(
            template_type=template_type,
            ec_class=ec_class,
            description=description,
            reversible=True
        )
    
    elif ec_class == ECClass.ISOMERASE:
        # EC 5: Isomerases
        template_type = "isomerization"
        description = "Intramolecular rearrangement"
        return ReactionTemplate(
            template_type=template_type,
            ec_class=ec_class,
            description=description,
            reversible=True
        )
    
    elif ec_class == ECClass.LIGASE:
        # EC 6: Ligases
        template_type = "ligation"
        description = "Bond formation coupled to ATP hydrolysis"
        return ReactionTemplate(
            template_type=template_type,
            ec_class=ec_class,
            cofactors=["ATP", "Mg2+"],
            description=description,
            reversible=False
        )
    
    elif ec_class == ECClass.TRANSLOCASE:
        # EC 7: Translocases
        template_type = "translocation"
        description = "Movement across membrane"
        return ReactionTemplate(
            template_type=template_type,
            ec_class=ec_class,
            description=description,
            reversible=False
        )
    
    else:
        # Unknown
        return ReactionTemplate(
            template_type="generic",
            ec_class=ECClass.UNKNOWN,
            description="Generic enzymatic reaction"
        )


def infer_products(
    substrate: str,
    enzyme_label: str,
    template: ReactionTemplate,
    cofactor: Optional[str] = None,
    database_products: Optional[str] = None
) -> List[str]:
    """
    Infer reaction products based on template and substrate.
    
    This is a simple rule-based approach. Future versions will use:
    - SMARTS pattern matching on substrate structure
    - Chemical transformation rules
    - Product prediction ML models
    
    Args:
        substrate (str): Substrate name
        enzyme_label (str): Enzyme name
        template (ReactionTemplate): Reaction template
        cofactor (str, optional): Cofactor if present
        database_products (str, optional): Products from database (if available)
        
    Returns:
        List[str]: Predicted product names
    """
    # If database has products, use them
    if database_products:
        # Split by common delimiters
        products = [p.strip() for p in re.split(r'[+,;]', database_products)]
        return [p for p in products if p]  # Remove empty strings
    
    # Otherwise, infer based on template
    products = []
    
    if template.template_type == "phosphorylation":
        # Substrate + ATP → Substrate-phosphate + ADP
        products.append(f"{substrate}-phosphate")
        products.append("ADP")
    
    elif template.template_type == "hydrolysis":
        # Substrate + H2O → Product fragments
        if "glycos" in template.template_type.lower() or "sacchar" in substrate.lower():
            # Glycoside hydrolysis typically produces monosaccharides
            products.append(f"{substrate} fragments")
        else:
            products.append(f"{substrate} (hydrolyzed)")
    
    elif template.template_type == "ester_hydrolysis":
        # Ester + H2O → Alcohol + Acid
        if "ATP" in substrate or "ADP" in substrate:
            products.append("ADP" if "ATP" in substrate else "AMP")
            products.append("Pi")
        else:
            products.append(f"{substrate} (hydrolyzed)")
    
    elif template.template_type.startswith("oxidation"):
        # Substrate + NAD+ → Substrate-oxidized + NADH
        if cofactor and "NAD" in cofactor.upper():
            products.append(f"{substrate} (oxidized)")
            products.append("NADH" if "NAD+" in cofactor else "NADPH")
        else:
            products.append(f"{substrate} (oxidized)")
    
    elif template.template_type == "isomerization":
        # Substrate → Substrate-isomer
        products.append(f"{substrate} (isomer)")
    
    elif template.template_type == "ligation":
        # Substrate1 + Substrate2 + ATP → Product + ADP + Pi
        products.append(f"{substrate} (ligated)")
        products.append("ADP")
        products.append("Pi")
    
    elif template.template_type == "group_transfer":
        # Substrate + Acceptor → Substrate-group + Acceptor+group
        products.append(f"{substrate} (modified)")
    
    else:
        # Generic/unknown
        products.append(f"{substrate} (product)")
    
    logger.debug(f"Inferred products for {substrate} + {enzyme_label}: {products}")
    return products


def create_reaction_from_database(
    substrate: str,
    enzyme_label: str,
    ec_number: str,
    cofactor: Optional[str] = None,
    database_products: Optional[str] = None
) -> ReactionTemplate:
    """
    Create a complete reaction template with products from database information.
    
    Args:
        substrate (str): Substrate name
        enzyme_label (str): Enzyme name
        ec_number (str): EC number
        cofactor (str, optional): Cofactor name
        database_products (str, optional): Products from database
        
    Returns:
        ReactionTemplate: Complete reaction template with inferred/known products
    """
    # Determine template from EC number
    template = determine_template_from_ec(ec_number, substrate)
    
    # Set reactants
    template.reactants = [substrate]
    if cofactor and cofactor.lower() not in ['none', 'null', '']:
        template.reactants.append(cofactor)
    
    # Infer products
    template.products = infer_products(
        substrate=substrate,
        enzyme_label=enzyme_label,
        template=template,
        cofactor=cofactor,
        database_products=database_products
    )
    
    # Create stoichiometry
    for reactant in template.reactants:
        template.stoichiometry[reactant] = -1
    for product in template.products:
        template.stoichiometry[product] = 1
    
    logger.info(f"Created reaction template: {template}")
    return template




