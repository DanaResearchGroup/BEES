#!/usr/bin/env python3

"""
this module is used to calculate the flux of a species in a reaction network and the characteristic rate of the core mode
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Set

from bees.reaction_generator import GeneratedReaction


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

        # Find the appropriate Km.
        #
        # IMPORTANT: when `km_per_substrate` exists, only apply saturation terms
        # for substrates that are explicitly present in it. Do NOT fall back to a
        # global/single Km for other reactants (e.g., H+), because that can
        # incorrectly suppress flux by treating buffered/auxiliary reactants as
        # kinetic substrates.
        km_val = None
        if km_per:
            km_val = km_per.get(reactant)
            if km_val is None:
                km_val = next(
                    (v for k, v in km_per.items() if k.lower().strip() == reactant_lc),
                    None,
                )
            if km_val is None:
                # No per-substrate Km for this reactant -> assume saturated (factor=1).
                continue
        else:
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


def identify_significant_species_at_interrupt(
    edge_rates: Dict[str, float],
    char_rate: float,
    tol_move_to_core: float,
    max_objects: int = 10,
    abs_flux_floor: float = 1e-12,
) -> List[SpeciesFlux]:
    """
     At the exact moment the solver is interrupted (``t_interrupt``), compute
    ``rr_i = |R_i| / R_char`` for each edge species *i*.  Species whose
    ``rr_i >= toleranceMoveToCore`` are candidates.  The list is sorted by
    ``rr_i`` descending and truncated to ``max_objects``.

    Flat-core safeguard (``char_rate <= 0`` but edge flux exists): ratio-based
    promotion is meaningless because any tiny ``|R_i|`` would produce an
    infinite ratio.  Instead, promote edge species whose ``|R_i|`` exceeds the
    absolute flux floor ``abs_flux_floor``, sorted by ``|R_i|`` descending and
    capped to ``max_objects``.  This avoids spurious promotions from numerical
    noise while still catching species with real flux.

    Args:
        edge_rates: label_lc -> instantaneous dC/dt (mM/s) at interrupt time.
        char_rate: Instantaneous R_char at interrupt time.
        tol_move_to_core: Tolerance epsilon (toleranceMoveToCore).
        max_objects: Maximum number of species to return per interrupt.
        abs_flux_floor: Absolute |rate| threshold used when R_char is zero.

    Returns:
        Sorted list (descending by rr_i or |rate|) of SpeciesFlux candidates,
        truncated to *max_objects*.
    """
    if not edge_rates:
        return []

    if char_rate <= 0.0:
        candidates: List[SpeciesFlux] = []
        for label, rate in edge_rates.items():
            if abs(rate) > abs_flux_floor:
                candidates.append(
                    SpeciesFlux(
                        label=label,
                        rate=rate,
                        normalized_rate=float("inf"),
                    )
                )
        candidates.sort(key=lambda sf: abs(sf.rate), reverse=True)
        return candidates[:max_objects]

    candidates = []
    for label, rate in edge_rates.items():
        rr = abs(rate) / char_rate
        if rr >= tol_move_to_core:
            candidates.append(
                SpeciesFlux(label=label, rate=rate, normalized_rate=rr)
            )

    candidates.sort(key=lambda sf: sf.normalized_rate, reverse=True)
    return candidates[:max_objects]


def identify_insignificant_species_from_peak_ratios(
    max_edge_rate_ratio: Dict[str, float],
    max_char_rate: float,
    tol_keep_in_edge: float,
    ineligible_for_prune: Optional[Set[str]] = None,
) -> Set[str]:
    """
    Prune edge species whose *peak* |R_edge|/R_char falls below tol_keep_in_edge.

    prunes using aggregated maximum rate ratios over the run, not a
    single end-time snapshot. Species listed in ineligible_for_prune are
    skipped (e.g. too young to prune).
    """
    if max_char_rate <= 0.0 or tol_keep_in_edge <= 0.0:
        return set()

    ineligible = ineligible_for_prune or set()
    to_remove: Set[str] = set()

    for label_lc, rr in max_edge_rate_ratio.items():
        if label_lc in ineligible:
            continue
        if rr < tol_keep_in_edge:
            to_remove.add(label_lc)

    return to_remove
