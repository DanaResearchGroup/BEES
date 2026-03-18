#!/usr/bin/env python3

"""
Flux Calculator Module
----------------------
Calculates species production/consumption rates for the rate-based
model enlargement algorithm.

Computes the characteristic rate R_char of the core model and
identifies edge species whose flux exceeds the user-specified
tolerance, marking them for promotion to the core.


"""

import math
from dataclasses import dataclass
from typing import Dict, List, Set

from bees.model_generator import GeneratedReaction


@dataclass
class SpeciesFlux:
    """Flux information for a single species at a given time point."""
    label: str
    rate: float  # mM/s  (dC/dt)
    normalized_rate: float = 0.0  # |rate| / R_char
    time: float = 0.0


def compute_mm_rate(
    reaction: GeneratedReaction,
    concentrations: Dict[str, float],
    enzyme_concentrations: Dict[str, float],
) -> float:
    """
    Compute the reaction rate for a single reaction using its assigned
    rate law (Michaelis-Menten only for now)

    For Michaelis-Menten:
        v = kcat * [E] * prod_i( [S_i] / (Km_i + [S_i]) )
        (Vmax = kcat * [E])


    Returns 0.0 when required kinetic parameters are missing.
    """
    kin = reaction.kinetics
    if kin is None or reaction.rate_law is None:
        return 0.0

    # Enzyme concentration (mM)
    enzyme_key = reaction.enzyme_label.lower().strip()
    e_conc = enzyme_concentrations.get(enzyme_key, 0.0)

    # Determine Vmax-equivalent
    kcat = kin.kcat  # 1/s
    vmax = kin.vmax  # mM/s

    if kcat is not None and e_conc > 0:
        v_max_eff = kcat * e_conc
    elif vmax is not None:
        v_max_eff = vmax
    else:
        return 0.0

    # Substrate saturation terms  prod_i( [S_i] / (Km_i + [S_i]) )
    km_per = getattr(kin, "km_per_substrate", None) or {}
    km_single = kin.km  # fallback single Km

    saturation = 1.0
    for reactant in reaction.reactant_labels:
        reactant_lc = reactant.lower().strip()
        s_conc = concentrations.get(reactant_lc, 0.0)
        s_conc = max(s_conc, 0.0)  # guard against negative from ODE

        # Find the appropriate Km (case-insensitive fallback for DB key mismatch)
        km_val = km_per.get(reactant)
        if km_val is None:
            km_val = next(
                (v for k, v in km_per.items() if k.lower().strip() == reactant_lc),
                None,
            )
        if km_val is None:
            km_val = km_single
        if km_val is None or km_val <= 0:
            # Without Km, assume saturated (saturation factor = 1)
            continue

        saturation *= s_conc / (km_val + s_conc)
        if saturation == 0.0:
            return 0.0

    return v_max_eff * saturation


def calculate_species_rates(
    reactions: List[GeneratedReaction],
    concentrations: Dict[str, float],
    enzyme_concentrations: Dict[str, float],
) -> Dict[str, float]:
    """
    Calculate net dC/dt for every species that appears in *reactions*.

    Args:
        reactions: List of reactions to evaluate.
        concentrations: Current species concentrations (label_lc -> mM).
        enzyme_concentrations: Enzyme concentrations (label_lc -> mM).

    Returns:
        Dictionary mapping species label (lowercase) to its net rate (mM/s).
    """
    rates: Dict[str, float] = {}

    for rxn in reactions:
        v = compute_mm_rate(rxn, concentrations, enzyme_concentrations)
        if v == 0.0:
            continue

        stoich = rxn.stoichiometry
        for species_label, coeff in stoich.items():
            lc = species_label.lower().strip()
            rates[lc] = rates.get(lc, 0.0) + coeff * v

    return rates


def calculate_characteristic_rate(core_species_rates: Dict[str, float]) -> float:
    """
    Compute the characteristic rate R_char of the core model.

    R_char = sqrt( sum_j  R_j^2 )   for species j in core.

    Args:
        core_species_rates: label_lc -> dC/dt for core species.

    Returns:
        R_char in mM/s.
    """
    if not core_species_rates:
        return 0.0
    return math.sqrt(sum(r * r for r in core_species_rates.values()))


def identify_significant_species(
    edge_species_rates: Dict[str, float],
    r_char: float,
    tol_move_to_core: float,
) -> List[SpeciesFlux]:
    """
    Return edge species whose absolute rate exceeds epsilon * R_char.

    Args:
        edge_species_rates: label_lc -> dC/dt for edge species.
        r_char: Characteristic rate of the core model (mM/s).
        tol_move_to_core: Tolerance epsilon.

    Returns:
        Sorted list (descending by |rate|) of SpeciesFlux objects that
        exceed the threshold.
    """
    if r_char <= 0.0 or not edge_species_rates:
        return []

    threshold = tol_move_to_core * r_char
    significant: List[SpeciesFlux] = []

    for label_lc, rate in edge_species_rates.items():
        abs_rate = abs(rate)
        if abs_rate >= threshold:
            significant.append(SpeciesFlux(
                label=label_lc,
                rate=rate,
                normalized_rate=abs_rate / r_char,
            ))

    significant.sort(key=lambda sf: abs(sf.rate), reverse=True)
    return significant


def identify_insignificant_species(
    edge_species_rates: Dict[str, float],
    r_char: float,
    tol_keep_in_edge: float,
) -> Set[str]:
    """
    Return edge species whose absolute rate is below the keep-in-edge
    tolerance, suitable for pruning.

    Args:
        edge_species_rates: label_lc -> dC/dt for edge species.
        r_char: Characteristic rate of the core model (mM/s).
        tol_keep_in_edge: Tolerance below which species are pruned.

    Returns:
        Set of lowercase labels to prune.
    """
    if r_char <= 0.0:
        return set()

    threshold = tol_keep_in_edge * r_char
    to_remove: Set[str] = set()

    for label_lc, rate in edge_species_rates.items():
        if abs(rate) < threshold:
            to_remove.add(label_lc)

    return to_remove
