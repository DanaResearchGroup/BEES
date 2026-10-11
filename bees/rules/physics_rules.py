"""BEES physics rules — always-on physical-chemistry corrections.

Opt-in FAS calibrations: bees.rules.calibrations.fas.
"""

from __future__ import annotations

import math
from typing import Iterable, Optional

from bees.cofactors import COFACTORS_ALWAYS_AVAILABLE
from bees.rules.base import RULES, Reference, Rule
from bees.thermodynamics import STANDARD_CONCENTRATION_MM
from bees.rules.helpers import (  # noqa: F401  — re-export detect_acyl_chain_length
    detect_acyl_chain_length,
    _CO2_LABELS,
    _CO2_SMILES,
)

DGR_IRREVERSIBLE_KJMOL = 30.0
_DG_PER_CH2_KJMOL = -1.0
_HYDROPHOBIC_REFERENCE_CHAIN = 2
_HYDROPHOBIC_TEMPERATURE_K = 298.15
_R_KJ = 8.314462618e-3

def compute_haldane_reverse_kcat(
    kcat_fwd: float,
    keq: float,
    km_substrates: Iterable[float],
    km_products: Iterable[float],
    nu_substrates: Optional[Iterable[float]] = None,
    nu_products: Optional[Iterable[float]] = None,
) -> Optional[float]:
    """kcat_rev = kcat_fwd · ∏(Km_p/c°)^νp / (Keq · ∏(Km_s/c°)^νs).

    Kms in mM; Keq dimensionless against the 1 M standard state (c° = 1000 mM).
    ν defaults to 1 per Km. Matches the saturated reverse limit of the CM rate law.
    """
    if not (
        math.isfinite(keq)
        and keq > 0.0
        and math.isfinite(kcat_fwd)
        and kcat_fwd > 0.0
    ):
        return None

    def _log_prod(kms, nus) -> Optional[float]:
        kms = list(kms)
        nus = [1.0] * len(kms) if nus is None else list(nus)
        if len(nus) != len(kms):
            return None
        total = 0.0
        for km, nu in zip(kms, nus):
            if km is None or km <= 0 or not math.isfinite(km):
                return None
            total += nu * math.log(km / STANDARD_CONCENTRATION_MM)
        return total

    log_p = _log_prod(km_products, nu_products)
    log_s = _log_prod(km_substrates, nu_substrates)
    if log_p is None or log_s is None:
        return None
    log_val = math.log(kcat_fwd) + log_p - math.log(keq) - log_s
    if log_val > 700.0:
        return None
    val = math.exp(log_val)
    if not (math.isfinite(val) and val > 0.0):
        return None
    return val

def _has_co2_product(reaction) -> bool:
    """True if CO2 is a product (label or SMILES)."""
    stoich = getattr(reaction, "stoichiometry", None) or {}
    product_labels = getattr(reaction, "product_labels", None) or [
        lab for lab, c in stoich.items() if c > 0
    ]
    kin = getattr(reaction, "kinetics", None)
    smiles_map = (getattr(kin, "compound_smiles", None) or {}) if kin is not None else {}

    co2_canon = None
    for lab in product_labels:
        if stoich.get(lab, 0) <= 0 and lab in stoich:
            continue
        if lab.lower().strip() in _CO2_LABELS:
            return True
        smi = smiles_map.get(lab) or smiles_map.get(lab.lower().strip())
        if smi:
            from bees.common import canonical_smiles
            if co2_canon is None:
                co2_canon = canonical_smiles(_CO2_SMILES)
            if canonical_smiles(smi) == co2_canon:
                return True
    return False

class DgrIrreversibility(Rule):
    """Flag forward-irreversible when ΔG°' < −cutoff or Keq invalid.

    Strongly endergonic steps (ΔG°' > +cutoff) are left reversible so they
    are not forced downhill the wrong way. Invalid/non-positive Keq still
    cannot support a reversible rate law.
    """

    def applies_to(self, reaction) -> bool:
        thermo = getattr(reaction, "thermo", None)
        if thermo is None:
            return False
        return isinstance(getattr(thermo, "dgr_prime_kJmol", None), (int, float)) and \
            isinstance(getattr(thermo, "keq", None), (int, float))

    def apply(self, reaction) -> None:
        thermo = reaction.thermo
        dgr = thermo.dgr_prime_kJmol
        keq = thermo.keq
        cutoff = float(self.params["dgr_kjmol_cutoff"])
        if (
            (math.isfinite(dgr) and dgr < -cutoff)
            or not math.isfinite(keq)
            or keq <= 0.0
        ):
            thermo.irreversible = True
            thermo.kcat_rev = None

class DecarboxylationIrreversible(Rule):
    """Flag irreversible when CO2 is a product."""

    def applies_to(self, reaction) -> bool:
        thermo = getattr(reaction, "thermo", None)
        if thermo is None or getattr(thermo, "irreversible", True):
            return False
        return _has_co2_product(reaction)

    def apply(self, reaction) -> None:
        reaction.thermo.irreversible = True
        reaction.thermo.kcat_rev = None

class HydrophobicChainLengthKm(Rule):
    """Scale acyl-substrate Km by chain length."""

    def applies_to(self, reaction) -> bool:
        kin = getattr(reaction, "kinetics", None)
        if kin is None:
            return False
        km_per = getattr(kin, "km_per_substrate", None)
        if not km_per:
            return False
        n_ref = int(self.params["reference_chain_length"])
        stoich = getattr(reaction, "stoichiometry", None) or {}
        smiles_map = getattr(kin, "compound_smiles", None) or {}
        for lab in km_per:
            if stoich.get(lab, 0) >= 0:
                continue
            if lab.lower().strip() in COFACTORS_ALWAYS_AVAILABLE:
                continue
            n = detect_acyl_chain_length(lab, smiles_map.get(lab))
            if n is not None and n >= n_ref:
                return True
        return False

    def apply(self, reaction) -> None:
        kin = reaction.kinetics
        km_per = kin.km_per_substrate
        stoich = getattr(reaction, "stoichiometry", None) or {}
        smiles_map = getattr(kin, "compound_smiles", None) or {}
        dG = float(self.params["dG_per_CH2_kJmol"])
        n_ref = int(self.params["reference_chain_length"])
        rt = _R_KJ * float(self.params["temperature_K"])
        for lab, km in list(km_per.items()):
            if km is None or km <= 0:
                continue
            if stoich.get(lab, 0) >= 0:
                continue
            if lab.lower().strip() in COFACTORS_ALWAYS_AVAILABLE:
                continue
            n = detect_acyl_chain_length(lab, smiles_map.get(lab))
            if n is None or n < n_ref:
                continue
            km_per[lab] = km * math.exp(dG * (n - n_ref) / rt)

class HaldaneReverseKcat(Rule):
    """Recompute kcat_rev via Haldane from the Kms the rate law uses.

    Runs after the Km rules, so kcat_rev is the saturated reverse turnover of the
    simulated law. Only this reaction's participants count; buffered cofactors have
    no Km (their activity enters as c/c°).
    """

    def applies_to(self, reaction) -> bool:
        thermo = getattr(reaction, "thermo", None)
        if thermo is None or getattr(thermo, "irreversible", True):
            return False
        kin = getattr(reaction, "kinetics", None)
        if kin is None:
            return False
        kcat = getattr(kin, "kcat", None)
        if not isinstance(kcat, (int, float)):
            return False
        keq = getattr(thermo, "keq", None)
        if not isinstance(keq, (int, float)):
            return False
        return kcat > 0

    def apply(self, reaction) -> None:
        thermo = reaction.thermo
        kin = reaction.kinetics
        stoich = getattr(reaction, "stoichiometry", None) or {}
        km_by_lc = {
            str(lab).lower().strip(): km
            for lab, km in (getattr(kin, "km_per_substrate", None) or {}).items()
        }
        substrates, products = [], []
        for lab, coeff in stoich.items():
            lab_lc = lab.lower().strip()
            km = km_by_lc.get(lab_lc)
            if coeff == 0 or lab_lc in COFACTORS_ALWAYS_AVAILABLE or not km or km <= 0:
                continue
            (substrates if coeff < 0 else products).append((km, abs(coeff)))
        kcat_rev = compute_haldane_reverse_kcat(
            kcat_fwd=kin.kcat,
            keq=thermo.keq,
            km_substrates=[km for km, _ in substrates],
            km_products=[km for km, _ in products],
            nu_substrates=[nu for _, nu in substrates],
            nu_products=[nu for _, nu in products],
        )
        if kcat_rev is None:
            thermo.kcat_rev = None
            thermo.irreversible = True
        else:
            thermo.kcat_rev = kcat_rev

_noor_2014 = Reference(
    authors=(
        "Noor, E.",
        "Bar-Even, A.",
        "Flamholz, A.",
        "Reznik, E.",
        "Liebermeister, W.",
        "Milo, R.",
    ),
    title=(
        "Pathway thermodynamics highlights kinetic obstacles in central metabolism"
    ),
    year="2014",
    journal="PLoS Computational Biology",
    volume="10",
    pages="e1003483",
    doi="10.1371/journal.pcbi.1003483",
)
_haldane_1930 = Reference(
    authors=("Haldane, J. B. S.",),
    title="Enzymes",
    year="1930",
    journal="Longmans, Green and Co.",
)
_ruppe_fox_2018 = Reference(
    authors=("Ruppe, A.", "Fox, J. M."),
    title=(
        "Analysis of interdependent kinetic controls of fatty acid synthases"
    ),
    year="2018",
    journal="ACS Catalysis",
    volume="8",
    pages="11722-11734",
    doi="10.1021/acscatal.8b03171",
)
_tanford_1980 = Reference(
    authors=("Tanford, C.",),
    title=(
        "The Hydrophobic Effect: Formation of Micelles and Biological Membranes"
    ),
    year="1980",
    journal="Wiley",
)
RULES.register(
    DgrIrreversibility(
        name="dgr_irreversibility",
        description=(
            "Flag forward-irreversible when ΔG°' < −30 kJ/mol or Keq invalid."
        ),
        reference=_noor_2014,
        reference_type="theoretical",
        params={"dgr_kjmol_cutoff": DGR_IRREVERSIBLE_KJMOL},
    )
)

RULES.register(
    DecarboxylationIrreversible(
        name="decarboxylation_irreversible",
        description="Flag irreversible when CO₂ is a product (gas escape).",
        reference=_ruppe_fox_2018,
        reference_type="experimental",
        params={},
    )
)

RULES.register(
    HydrophobicChainLengthKm(
        name="hydrophobic_chain_length_km",
        description="Scale acyl-substrate Km by exp(ΔG_CH2·(n−n_ref)/RT); ΔG_CH2≈−1.0 kJ/mol.",
        reference=_tanford_1980,
        reference_type="textbook",
        params={
            "dG_per_CH2_kJmol": _DG_PER_CH2_KJMOL,
            "reference_chain_length": _HYDROPHOBIC_REFERENCE_CHAIN,
            "temperature_K": _HYDROPHOBIC_TEMPERATURE_K,
        },
    )
)

RULES.register(
    HaldaneReverseKcat(
        name="haldane_reverse_kcat",
        description="Compute kcat_rev via Haldane; flag irreversible if unevaluable.",
        reference=_haldane_1930,
        reference_type="textbook",
        params={},
    )
)
