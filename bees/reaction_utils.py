#!/usr/bin/env python3

"""
Reaction Utilities Module
-------------------------
General utility functions for reaction validation, reactant checking, and cofactor handling.
These functions are pathway-agnostic and can be used across different biochemical pathways.
"""

from typing import List, Set, Optional, Dict
import json
import time
from bees.common import (
    GENERAL_COFACTORS,
    get_ontology_equivalents,
    get_coenzyme_like_flags,
    EC_ALIASES,
    COFACTOR_SUBSTITUTIONS,
    ENZYME_DOMAIN_COFACTORS
)


def get_ec_aliases(ec_number: str) -> List[str]:
    """
    Get all EC number aliases for a given EC number.
    
    Args:
        ec_number (str): Primary EC number (e.g., "EC 2.3.1.85")
        
    Returns:
        List[str]: List of EC numbers to try, including the primary EC number first
    """
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
    ec_number: Optional[str] = None
) -> tuple[bool, Optional[str]]:
    """
    Check if a reactant is available, considering:
    - General cofactors (implicitly available)
    - Enzyme domain cofactors (part of enzyme structure)
    - Cofactor substitutions
    - Ontology equivalents
    
    Args:
        reactant (str): Reactant name to check
        available_species_labels_lc (Set[str]): Set of available species (lowercase)
        enzyme_label (str, optional): Enzyme label for domain cofactor checking
        ec_number (str, optional): EC number for context-specific substitutions
        
    Returns:
        tuple[bool, Optional[str]]: (is_available, reason)
            - is_available: True if reactant is available
            - reason: Optional explanation (e.g., "general_cofactor", "domain_cofactor", "substitution", "ontology")
    """
    r_lc = str(reactant).lower().strip()
    coenzyme_flags = get_coenzyme_like_flags(reactant)
    matched_domain_pattern = None

    # Check if it's a general cofactor (implicitly available)
    if r_lc in GENERAL_COFACTORS:
        reason = "general_cofactor"
        is_available = True
    # Check if it's directly available
    elif r_lc in available_species_labels_lc:
        reason = "direct_match"
        is_available = True
    # Check if it's a domain cofactor (part of enzyme structure)
    # Skip for acyl-ACP: acetyl-ACP, hexanoyl-ACP etc. must come from prior reactions, not the domain
    elif enzyme_label and not coenzyme_flags.get("is_acyl_acp", False):
        domain_cofactors = get_enzyme_domain_cofactors(enzyme_label)
        # Ensure patterns are lowercase for consistent matching
        domain_cofactors_lc = [pattern.lower() for pattern in domain_cofactors]
        for pattern in domain_cofactors_lc:
            if pattern in r_lc:
                matched_domain_pattern = pattern
                break
        if matched_domain_pattern:
            reason = "domain_cofactor"
            is_available = True
        else:
            reason = None
            is_available = False
    # Check for cofactor substitutions
    elif r_lc in COFACTOR_SUBSTITUTIONS:
        for substitute in COFACTOR_SUBSTITUTIONS[r_lc]:
            if substitute in available_species_labels_lc:
                reason = "substitution"
                is_available = True
                break
        else:
            reason = None
            is_available = False
    # Check ontology equivalents
    else:
        equivalents = get_ontology_equivalents(r_lc)
        if any(eq in available_species_labels_lc for eq in equivalents):
            reason = "ontology"
            is_available = True
        else:
            reason = None
            is_available = False

    return is_available, reason


def validate_reaction_reactants(
    reactants: List[str],
    available_species_labels_lc: Set[str],
    enzyme_label: Optional[str] = None,
    ec_number: Optional[str] = None
) -> tuple[bool, List[str]]:
    """
    Validate that all reactants for a reaction are available.
    
    Args:
        reactants (List[str]): List of reactant names
        available_species_labels_lc (Set[str]): Set of available species (lowercase)
        enzyme_label (str, optional): Enzyme label for domain cofactor checking
        ec_number (str, optional): EC number for context-specific substitutions
        
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
            ec_number=ec_number
        )
        
        if not is_available:
            missing_reactants.append(reactant)
    
    return len(missing_reactants) == 0, missing_reactants
