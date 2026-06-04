#!/usr/bin/env python3

from typing import List, Set, Optional
from bees.common import (
    COFACTORS_ALWAYS_AVAILABLE,
    get_ontology_equivalents,
    get_coenzyme_like_flags,
    EC_ALIASES,
    ENZYME_DOMAIN_COFACTORS,
)


def get_ec_aliases(ec_number: Optional[str]) -> List[str]:
    """
    Get all EC number aliases for a given EC number.
    
    Args:
        ec_number (str): Primary EC number (e.g., "EC 2.3.1.85")
        
    Returns:
        List[str]: List of EC numbers to try, including the primary EC number first
    """
    if not ec_number:
        return []

    ec_numbers_to_try = [ec_number]
    if ec_number in EC_ALIASES:
        ec_numbers_to_try.extend(EC_ALIASES[ec_number])
    return ec_numbers_to_try


def get_enzyme_domain_cofactors(enzyme_label: str) -> List[str]:
    """
    Get domain cofactors for an enzyme (cofactors that are part of the enzyme structure).
    
    Args:
        enzyme_label (str): Enzyme name/label
        
    Returns:
        List[str]: List of cofactor patterns that are part of this enzyme's structure
    """
    enzyme_lc = enzyme_label.lower()
    domain_cofactors = []
    for enzyme_pattern, cofactor_patterns in ENZYME_DOMAIN_COFACTORS.items():
        if enzyme_pattern in enzyme_lc:
            domain_cofactors.extend(cofactor_patterns)
    return domain_cofactors


def check_reactant_availability(
    reactant: str,
    available_species_labels_lc: Set[str],
    enzyme_label: Optional[str] = None,
) -> tuple[bool, Optional[str]]:
    """
    Check whether a reactant is "available" for reaction generation.

    Checks are applied in this order:
    1) Always-available cofactors (e.g., H2O, H+, Pi; see COFACTORS_ALWAYS_AVAILABLE)
    2) Direct match in `available_species_labels_lc`
    3) Enzyme domain cofactors (if `enzyme_label` given; skipped for acyl-ACP reactants)
    4) Ontology equivalents
    
    Args:
        reactant (str): Reactant name to check
        available_species_labels_lc (Set[str]): Set of available species (lowercase)
        enzyme_label (str, optional): Enzyme label for domain cofactor checking
        
    Returns:
        tuple[bool, Optional[str]]: (is_available, reason)
            - is_available: True if reactant is available
            - reason: One of:
              "always_available_cofactor", "direct_match", "domain_cofactor",
              "ontology", or None if unavailable.
    """
    r_lc = str(reactant).lower().strip()
    coenzyme_flags = get_coenzyme_like_flags(reactant)

    # Check if it's an always-available cofactor (implicitly available).
    # NOTE: This uses COFACTORS_ALWAYS_AVAILABLE, which is a restricted subset
    # of GENERAL_COFACTORS (H2O, H+, Pi, inorganic ions, etc.), so that high‑
    # energy carriers like ATP / NAD(H)/NADP(H) still need to be provided in
    # the input species list or produced in the network.
    if r_lc in COFACTORS_ALWAYS_AVAILABLE:
        return True, "always_available_cofactor"
    # Check if it's directly available
    if r_lc in available_species_labels_lc:
        return True, "direct_match"

    # Check if it's a domain cofactor (part of enzyme structure)
    # Skip for acyl-ACP: acetyl-ACP, hexanoyl-ACP etc. must come from prior reactions, not the domain
    if enzyme_label and not coenzyme_flags.get("is_acyl_acp", False):
        domain_cofactors = get_enzyme_domain_cofactors(enzyme_label)
        # Ensure patterns are lowercase for consistent matching
        domain_cofactors_lc = [pattern.lower() for pattern in domain_cofactors]
        if any(pattern in r_lc for pattern in domain_cofactors_lc):
            return True, "domain_cofactor"

    # Check ontology equivalents
    equivalents = get_ontology_equivalents(r_lc)
    if any(eq in available_species_labels_lc for eq in equivalents):
        return True, "ontology"

    return False, None


def validate_reaction_reactants(
    reactants: List[str],
    available_species_labels_lc: Set[str],
    enzyme_label: Optional[str] = None,
) -> tuple[bool, List[str]]:
    """
    Validate that all reactants for a reaction are available.
    
    Args:
        reactants (List[str]): List of reactant names
        available_species_labels_lc (Set[str]): Set of available species (lowercase)
        enzyme_label (str, optional): Enzyme label for domain cofactor checking
        
    Returns:
        tuple[bool, List[str]]: (all_available, missing_reactants)
            - all_available: True if all reactants are available
            - missing_reactants: List of reactants that are not available
    """
    missing_reactants = []
    
    for reactant in reactants:
        is_available, reason = check_reactant_availability(
            reactant,
            available_species_labels_lc,
            enzyme_label=enzyme_label,
        )
        
        if not is_available:
            missing_reactants.append(reactant)
    
    return len(missing_reactants) == 0, missing_reactants
