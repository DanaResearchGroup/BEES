#!/usr/bin/env python3

"""Flux and characteristic-rate calculations for the core/edge reaction network."""

import math
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Set, Tuple

from bees.cofactors import is_rate_law_exempt_cofactor
from bees.reaction_generator import GeneratedReaction


@dataclass
class SpeciesFlux:
    label: str
    rate: float  # mM/s  (dC/dt)
    normalized_rate: float = 0.0  # |rate| / R_char


def explicit_product_km(reaction: GeneratedReaction, label: str) -> Optional[float]:
    """
    Strict product-Km lookup: km_per_substrate ONLY, no fallback to kin.km.

    kin.km is the substrate Km; using it for a product would produce fake
    product inhibition. Returns None when no explicit entry is found.
    """
    kin = reaction.kinetics
    if kin is None:
        return None
    km_per = getattr(kin, "km_per_substrate", None) or {}
    if not km_per:
        return None
    km = km_per.get(label)
    if km is None:
        label_lc = label.lower().strip()
        km = next(
            (v for k, v in km_per.items() if k.lower().strip() == label_lc),
            None,
        )
    return km if (km is not None and km > 0) else None


def _product_labels_requiring_km(reaction: GeneratedReaction) -> List[str]:
    """Product labels that must have explicit Kms for CM / reversible MM.

    Buffered cofactors (H2O, H+, CO2, …) are omitted — they stay in stoichiometry
    and Q/Keq but must not gate the rate-law form.
    """
    products = getattr(reaction, "product_labels", None) or []
    return [p for p in products if not is_rate_law_exempt_cofactor(p)]


def has_complete_explicit_product_kms(reaction: GeneratedReaction) -> bool:
    """True iff every non-exempt product has a positive explicit Km in km_per_substrate.

    Exempt = buffered always-available cofactors (see is_rate_law_exempt_cofactor).
    Returns False when the reaction has no products at all.
    """
    products = getattr(reaction, "product_labels", None) or []
    if not products:
        return False
    required = _product_labels_requiring_km(reaction)
    # Only buffered products (e.g. pure CO2 release with nothing else) — treat
    # as complete so we do not silently fall back to legacy for lack of a
    # meaningless H2O/CO2 Km.
    if not required:
        return True
    return all(explicit_product_km(reaction, p) is not None for p in required)


def _compute_mm_rate_legacy(
    reaction: GeneratedReaction,
    concentrations: Dict[str, float],
    v_max_eff: float,
) -> float:
    """
    Legacy irreversible MM: v = Vmax * prod_i( S_i/(Km_i+S_i) )^νi.
    Products are ignored entirely. Used as the fallback when product Kms
    are incomplete.
    """
    kin = reaction.kinetics
    km_per = getattr(kin, "km_per_substrate", None) or {}
    km_single = kin.km
    stoich = reaction.stoichiometry

    saturation = 1.0
    for reactant in reaction.reactant_labels:
        reactant_lc = reactant.lower().strip()
        s_conc = max(concentrations.get(reactant_lc, 0.0), 0.0)

        km_val = None
        if km_per:
            km_val = km_per.get(reactant)
            if km_val is None:
                km_val = next(
                    (v for k, v in km_per.items() if k.lower().strip() == reactant_lc),
                    None,
                )
            if km_val is None:
                continue  # cofactor-skip: assume saturated
        else:
            km_val = km_single
        if km_val is None or km_val <= 0:
            continue

        nu = abs(stoich.get(reactant, 1))
        ratio = s_conc / (km_val + s_conc)
        saturation *= ratio ** nu
        if saturation == 0.0:
            return 0.0

    return v_max_eff * saturation


def compute_mm_rate(
    reaction: GeneratedReaction,
    concentrations: Dict[str, float],
    enzyme_concentrations: Dict[str, float],
) -> float:
    """
    Compute the irreversible reaction rate.

    When all product Kms are explicitly available in km_per_substrate, uses
    the common-modular (CM) symmetric denominator with forward-only flux:

        numerator   = Vmax * prod_i (S_i/Km_s,i)^νi
        denominator = prod_i (1+S_i/Km_s,i)^νi + prod_j (1+P_j/Km_p,j)^νj - 1
        v = numerator / denominator

    Falls back to legacy v = Vmax * prod_i (S_i/(Km_i+S_i))^νi when any
    product Km is missing or the reaction has no products.

    Returns 0.0 when required kinetic parameters are missing.
    """
    kin = reaction.kinetics
    if kin is None or reaction.rate_law is None:
        return 0.0

    enzyme_key = reaction.enzyme_label.lower().strip()
    e_conc = enzyme_concentrations.get(enzyme_key, 0.0)

    kcat = kin.kcat
    vmax = kin.vmax
    if kcat is not None and e_conc > 0:
        v_max_eff = kcat * e_conc
    elif vmax is not None:
        v_max_eff = vmax
    else:
        return 0.0

    if not has_complete_explicit_product_kms(reaction):
        return _compute_mm_rate_legacy(reaction, concentrations, v_max_eff)

    km_per = getattr(kin, "km_per_substrate", None) or {}
    km_single = kin.km
    stoich = reaction.stoichiometry

    # ---- Substrate terms ------------------------------------------------
    # Numerator:   prod_i (S_i/Km_s,i)^νi
    # Denominator: prod_i (1 + S_i/Km_s,i)^νi
    # Substrate Km semantics preserved: cofactor-skip when not in km_per.
    sub_sat_num = 1.0
    sub_sat_den = 1.0
    for reactant in reaction.reactant_labels:
        reactant_lc = reactant.lower().strip()
        s_conc = max(concentrations.get(reactant_lc, 0.0), 0.0)

        km_val = None
        if km_per:
            km_val = km_per.get(reactant)
            if km_val is None:
                km_val = next(
                    (v for k, v in km_per.items() if k.lower().strip() == reactant_lc),
                    None,
                )
            if km_val is None:
                continue  # cofactor-skip: assume saturated (factor 1 in num & den)
        else:
            km_val = km_single
        if km_val is None or km_val <= 0:
            continue

        nu = abs(stoich.get(reactant, 1))
        ratio = s_conc / km_val
        sub_sat_num *= ratio ** nu
        sub_sat_den *= (1.0 + ratio) ** nu
        if sub_sat_num == 0.0:
            return 0.0

    # ---- Product denominator terms --------------------------------------
    # prod_j (1 + P_j/Km_p,j)^νj  — skip buffered cofactors (H2O/H+/CO2/…)
    prod_sat_den = 1.0
    for product in _product_labels_requiring_km(reaction):
        km_p = explicit_product_km(reaction, product)
        p_conc = max(concentrations.get(product.lower().strip(), 0.0), 0.0)
        nu = abs(stoich.get(product, 1))
        prod_sat_den *= (1.0 + p_conc / km_p) ** nu

    denom = sub_sat_den + prod_sat_den - 1.0
    if denom <= 0.0 or not math.isfinite(denom):
        return 0.0

    return v_max_eff * sub_sat_num / denom


# Below this magnitude, |1 − Q/Keq| is treated as exact zero — protects the
# rate-law evaluation from catastrophic cancellation near equilibrium.
_EQUILIBRIUM_CLAMP = 1e-12


def _km_for(reaction: GeneratedReaction, label: str) -> Optional[float]:
    """Find Km for a single label using the same lookup rules as compute_mm_rate."""
    kin = reaction.kinetics
    if kin is None:
        return None
    km_per = getattr(kin, "km_per_substrate", None) or {}
    if km_per:
        km = km_per.get(label)
        if km is None:
            label_lc = label.lower().strip()
            km = next(
                (v for k, v in km_per.items() if k.lower().strip() == label_lc),
                None,
            )
        return km
    return kin.km


def compute_reversible_mm_rate(
    reaction: GeneratedReaction,
    concentrations: Dict[str, float],
    enzyme_concentrations: Dict[str, float],
) -> float:
    """
    Reversible Michaelis-Menten rate — common-modular (CM) kinetics
    (Liebermeister, Uhlendorf & Klipp 2010), with stoichiometric exponents.
    When all ν = 1 this coincides with convenience kinetics (2006).

    Requires reaction.thermo with irreversible=False and a finite Keq.
    Falls back to compute_mm_rate() when thermo is absent/irreversible, or
    when any *non-exempt* product Km is missing. Buffered cofactors (H2O,
    H+, CO2, …) are omitted from saturation terms but still enter Q via
    stoichiometry.

    Full general form with stoichiometric exponents νi:

        numerator   = ∏_i ([S_i]/Km_s,i)^νi  ·  (1 − Q/Keq)
        denominator = ∏_i (1 + [S_i]/Km_s,i)^νi
                    + ∏_j (1 + [P_j]/Km_p,j)^νj  − 1
        v = Vmax_f · numerator / denominator

    Q is computed in log-space from stoichiometry to avoid cancellation near
    equilibrium. Negative v means net reverse flux (Q > Keq).
    """
    kin = reaction.kinetics
    thermo = reaction.thermo
    if (
        kin is None
        or reaction.rate_law is None
        or thermo is None
        or thermo.irreversible
        or not math.isfinite(thermo.keq)
        or thermo.keq <= 0.0
    ):
        return compute_mm_rate(reaction, concentrations, enzyme_concentrations)

    enzyme_key = reaction.enzyme_label.lower().strip()
    e_conc = enzyme_concentrations.get(enzyme_key, 0.0)
    kcat_fwd = kin.kcat
    if kcat_fwd is None or e_conc <= 0:
        return compute_mm_rate(reaction, concentrations, enzyme_concentrations)

    vmax_f = kcat_fwd * e_conc

    # Stoichiometric coefficients per label (abs values; sign tracked by
    # reactant_labels / product_labels membership).
    stoich = reaction.stoichiometry  # label → signed int

    # ---- Substrate saturation (numerator and denominator terms) ------------
    # Buffered cofactors without Km are skipped (assume saturated), matching
    # legacy MM. Missing Km on a non-exempt substrate → fall back.
    sub_sat_num = 1.0   # ∏ ([S_i]/Km_s,i)^νi
    sub_sat_den = 1.0   # ∏ (1 + [S_i]/Km_s,i)^νi
    for reactant in reaction.reactant_labels:
        km = _km_for(reaction, reactant)
        if km is None or km <= 0:
            if is_rate_law_exempt_cofactor(reactant):
                continue
            return compute_mm_rate(reaction, concentrations, enzyme_concentrations)
        c = max(concentrations.get(reactant.lower().strip(), 0.0), 0.0)
        nu = abs(stoich.get(reactant, 1))        # stoichiometric exponent
        ratio = c / km
        sub_sat_num *= ratio ** nu
        sub_sat_den *= (1.0 + ratio) ** nu

    # ---- Product saturation (denominator term only) ------------------------
    # Skip buffered cofactors; require Km for every regulatory / organic product.
    prod_sat_den = 1.0   # ∏ (1 + [P_j]/Km_p,j)^νj
    for product in _product_labels_requiring_km(reaction):
        km_p = _km_for(reaction, product)
        if km_p is None or km_p <= 0:
            # No product Km → Haldane was not applied; fall back to forward-only.
            return compute_mm_rate(reaction, concentrations, enzyme_concentrations)
        c = max(concentrations.get(product.lower().strip(), 0.0), 0.0)
        nu = abs(stoich.get(product, 1))
        prod_sat_den *= (1.0 + c / km_p) ** nu

    # ---- Disequilibrium ratio Q/Keq (log-space to avoid cancellation) ------
    # Q = ∏ [P_j]^νp,j / ∏ [S_i]^νs,i  (stoich coefficients are signed).
    # Buffered always-available species (H2O, H+, CO2, …) are omitted: their
    # activities are already baked into the biochemical K'eq from equilibrator.
    # Including [H2O]=55 (mM pool) would push every dehydration far reverse.
    log_q = 0.0
    finite = True
    for label, coeff in stoich.items():
        if is_rate_law_exempt_cofactor(label):
            continue
        c = max(concentrations.get(label.lower().strip(), 0.0), 0.0)
        if c <= 0.0:
            if coeff > 0:
                # A product is absent → Q = 0 → full forward driving force
                log_q = float("-inf")
            else:
                # A substrate is absent → v = 0 (numerator already 0)
                return 0.0
            finite = False
            break
        log_q += coeff * math.log(c)

    if finite:
        log_q_over_keq = log_q - math.log(thermo.keq)
        q_over_keq = math.exp(log_q_over_keq) if log_q_over_keq < 700 else float("inf")
    else:
        q_over_keq = 0.0  # Q = 0 → disequilibrium = 1

    disequilibrium = 1.0 - q_over_keq
    if abs(disequilibrium) < _EQUILIBRIUM_CLAMP:
        return 0.0

    denom = sub_sat_den + prod_sat_den - 1.0
    if denom <= 0.0 or not math.isfinite(denom):
        return 0.0

    return vmax_f * sub_sat_num * disequilibrium / denom


@dataclass
class RateLawSpec:
    """Symbolic form of the rate law the simulator evaluates (feedback factor excluded).

    form: "reversible" (compute_reversible_mm_rate), "product_inhibited" (CM branch of
    compute_mm_rate), "mm" (legacy irreversible MM) or "zero". Exporters build their
    formulas from this so written models evaluate exactly like the simulator.
    """

    form: str
    kcat: Optional[float] = None
    vmax: Optional[float] = None
    substrates: List[Tuple[str, float, float]] = field(default_factory=list)  # (label, Km, nu)
    products: List[Tuple[str, float, float]] = field(default_factory=list)  # (label, Km, nu)
    q_substrates: List[Tuple[str, float]] = field(default_factory=list)  # (label, nu)
    q_products: List[Tuple[str, float]] = field(default_factory=list)  # (label, nu)
    keq: Optional[float] = None


def _vmax_source(reaction: GeneratedReaction, e_conc: float) -> Tuple[Optional[float], Optional[float]]:
    kin = reaction.kinetics
    if kin.kcat is not None and e_conc > 0:
        return kin.kcat, None
    if kin.vmax is not None:
        return None, kin.vmax
    return None, None


def _substrate_terms(reaction: GeneratedReaction) -> List[Tuple[str, float, float]]:
    """Substrates with a usable Km (others are treated as saturated), as in compute_mm_rate."""
    stoich = reaction.stoichiometry
    terms = []
    for reactant in reaction.reactant_labels:
        km = _km_for(reaction, reactant)
        if km is None or km <= 0:
            continue
        terms.append((reactant, float(km), float(abs(stoich.get(reactant, 1)))))
    return terms


def _irreversible_spec(reaction: GeneratedReaction, e_conc: float) -> RateLawSpec:
    kcat, vmax = _vmax_source(reaction, e_conc)
    if kcat is None and vmax is None:
        return RateLawSpec("zero")
    substrates = _substrate_terms(reaction)
    if not has_complete_explicit_product_kms(reaction):
        return RateLawSpec("mm", kcat=kcat, vmax=vmax, substrates=substrates)
    stoich = reaction.stoichiometry
    products = [
        (p, float(explicit_product_km(reaction, p)), float(abs(stoich.get(p, 1))))
        for p in _product_labels_requiring_km(reaction)
    ]
    return RateLawSpec("product_inhibited", kcat=kcat, vmax=vmax,
                       substrates=substrates, products=products)


def _reversible_spec(reaction: GeneratedReaction, e_conc: float) -> Optional[RateLawSpec]:
    """Mirror of compute_reversible_mm_rate; None where it falls back to compute_mm_rate."""
    kin = reaction.kinetics
    thermo = reaction.thermo
    keq = getattr(thermo, "keq", None)
    if not isinstance(keq, (int, float)) or not math.isfinite(keq) or keq <= 0.0:
        return None
    if kin.kcat is None or e_conc <= 0:
        return None
    stoich = reaction.stoichiometry
    substrates = []
    for reactant in reaction.reactant_labels:
        km = _km_for(reaction, reactant)
        if km is None or km <= 0:
            if is_rate_law_exempt_cofactor(reactant):
                continue
            return None
        substrates.append((reactant, float(km), float(abs(stoich.get(reactant, 1)))))
    products = []
    for product in _product_labels_requiring_km(reaction):
        km_p = _km_for(reaction, product)
        if km_p is None or km_p <= 0:
            return None
        products.append((product, float(km_p), float(abs(stoich.get(product, 1)))))
    q_s = [(lab, float(-c)) for lab, c in stoich.items() if c < 0 and not is_rate_law_exempt_cofactor(lab)]
    q_p = [(lab, float(c)) for lab, c in stoich.items() if c > 0 and not is_rate_law_exempt_cofactor(lab)]
    return RateLawSpec("reversible", kcat=kin.kcat, substrates=substrates, products=products,
                       q_substrates=q_s, q_products=q_p, keq=float(keq))


def rate_law_spec(reaction: GeneratedReaction, enzyme_concentrations: Dict[str, float]) -> RateLawSpec:
    """Which rate law the simulator uses for this reaction, with its parameters (see RateLawSpec)."""
    kin = reaction.kinetics
    if kin is None or reaction.rate_law is None:
        return RateLawSpec("zero")
    thermo = getattr(reaction, "thermo", None)
    e_conc = enzyme_concentrations.get(reaction.enzyme_label.lower().strip(), 0.0)
    is_reversible = (
        getattr(reaction.template, "reversible", False)
        and thermo is not None
        and not thermo.irreversible
    )
    if is_reversible:
        spec = _reversible_spec(reaction, e_conc)
        if spec is not None:
            return spec
    return _irreversible_spec(reaction, e_conc)


def format_rate_law(
    spec: RateLawSpec,
    conc: Callable[[str], Optional[str]],
    km: Callable[[str, str], str],
    vf: str,
    keq: str = "Keq",
) -> str:
    """Infix formula for ``spec`` (``^`` = power). ``conc`` returns None for species absent
    from the written model (concentration 0); ``km(label, "S"|"P")`` names Km parameters;
    ``vf`` is the Vmax expression (kcat*[E] or Vmax) and ``keq`` the Keq token.

    Reversible form: vf * prod_i (S_i/Km_i)^nu_i * (1 - Q/Keq) / (den_S + den_P - 1),
    written without dividing by concentrations so it stays finite at S = 0.
    """
    def pw(base: str, nu: float) -> str:
        return base if nu == 1 else f"({base})^{nu:g}"

    def prod(parts: List[str]) -> str:
        return " * ".join(parts) if parts else "1"

    if spec.form == "zero":
        return "0"

    if spec.form == "mm":
        parts = [vf]
        for lab, _, nu in spec.substrates:
            c = conc(lab)
            if c is None:
                return "0"
            parts.append(pw(f"{c} / ({km(lab, 'S')} + {c})", nu))
        return " * ".join(parts)

    for lab, _, _ in spec.substrates:
        if conc(lab) is None:
            return "0"
    sub_den = prod([pw(f"(1 + {conc(lab)} / {km(lab, 'S')})", nu) for lab, _, nu in spec.substrates])
    prod_den = prod([
        pw(f"(1 + {conc(lab)} / {km(lab, 'P')})", nu)
        for lab, _, nu in spec.products if conc(lab) is not None
    ])
    den = f"({sub_den} + {prod_den} - 1)"

    if spec.form == "product_inhibited":
        sub_num = prod([pw(f"({conc(lab)} / {km(lab, 'S')})", nu) for lab, _, nu in spec.substrates])
        return f"{vf} * {sub_num} / {den}"

    q_s = dict(spec.q_substrates)
    sat_labels = {lab for lab, _, _ in spec.substrates}
    q_s_only = [(lab, nu) for lab, nu in spec.q_substrates if lab not in sat_labels]
    if any(conc(lab) is None for lab, _ in q_s_only):
        return "0"
    both = [(lab, nu) for lab, _, nu in spec.substrates if lab in q_s]
    fwd = prod([pw(conc(lab), q_s[lab]) for lab, _ in both])
    if any(conc(lab) is None for lab, _ in spec.q_products):
        driving = fwd
    else:
        rev_den = " * ".join([keq] + [pw(conc(lab), nu) for lab, nu in q_s_only])
        rev = prod([pw(conc(lab), nu) for lab, nu in spec.q_products])
        driving = f"{fwd} - {rev} / ({rev_den})"
    factors = [vf] + [pw(f"({conc(lab)} / {km(lab, 'S')})", nu)
                      for lab, _, nu in spec.substrates if lab not in q_s]
    formula = " * ".join(factors + [f"({driving})"])
    if both:
        formula += f" / ({prod([pw(km(lab, 'S'), nu) for lab, nu in both])})"
    return f"{formula} / {den}"


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


def identify_significant_species_from_peak_ratios(
    max_edge_rate_ratio: Dict[str, float],
    tol_move_to_core: float,
    max_objects: int = 10,
) -> List[SpeciesFlux]:
    """
    Promote edge species whose *peak* |R_edge|/R_char over a full (non-interrupted)
    run reached tol_move_to_core.

    RMG adds any edge species whose rate ratio exceeds toleranceMoveToCore at any
    time during the simulation, not only at interrupt instants
    (rmgpy/solver/base.pyx, ``invalid_objects``). Without this, a species with
    tol_move_to_core <= rr < toleranceInterruptSimulation would never be promoted.

    Returns candidates sorted by peak ratio descending, truncated to *max_objects*.
    ``rate`` is not available from peak ratios and is set to 0.0.
    """
    candidates = [
        SpeciesFlux(label=label, rate=0.0, normalized_rate=rr)
        for label, rr in max_edge_rate_ratio.items()
        if rr >= tol_move_to_core
    ]
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

    Uses aggregated peak rate ratios (not a single end-time snapshot); skips species in ineligible_for_prune.
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
