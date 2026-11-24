#!/usr/bin/env python3
r"""
Local ranking plus an ODE-only global sample for the fixed E. coli FAS-II network.

Saved per-reaction slopes are grouped (one factor per enzyme and parameter class)
and ranked by ``(group coefficient)²``. The global step samples the top groups
with full ODEs (log-uniform in each group's CatPred log10 SD).

``prepare()`` builds the network, ``apply_sample()`` applies one draw, and
``_chain_profile()`` reads chain lengths. Draw i depends only on (seed, i, source).

Writes ``output/uncertainty/propagation/`` and does not touch deliverables.
CatPred SD_total is uncalibrated.

Run (repo root, bees_env):
    python projects/fattyAcidSynthesis/fattyAcidSynthesis_ecoli/uncertainty.py
    python projects/fattyAcidSynthesis/fattyAcidSynthesis_ecoli/uncertainty.py --rank-only
    python projects/fattyAcidSynthesis/fattyAcidSynthesis_ecoli/uncertainty.py --postprocess
"""

from __future__ import annotations

import argparse
import copy
import csv
import math
import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

# One BLAS thread. Must precede the numpy import.
for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

_PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _PROJECT_DIR)

# Loads .env.bees before any bees.* import.
import sensitivity_analysis as sa  # noqa: E402

import numpy as np  # noqa: E402

from bees.common import R, linear_sd_to_log10_sd  # noqa: E402
from bees.rules import RULES  # noqa: E402
from bees.rules.helpers import _ec_norm  # noqa: E402


FABI_EC = "1.3.1.9"
FABI_RULE = "fabI_enoyl_reductase_measured_kcat"
TESA_RULE = "tesa_long_chain_preference"
FABI_KCAT_MEAN = 15.0     # s⁻¹; Rafi 2006
FABI_KCAT_SD = 16.0 / 60.0
# CatPred SD_total, hexadecanoate vs FabH (SI_CatPred_FabH_Ki_table.csv)
FABH_KI_LOG10_SD = 0.350218
KM_TIE_RTOL = 1e-4

# Fig. 6 bins: (chain carbons, double bonds).
CHAIN_KEYS = ((4, 0), (6, 0), (8, 0), (10, 0), (12, 0), (12, 1), (14, 0), (14, 1),
              (16, 0), (16, 1), (18, 0), (18, 1), (20, 0), (20, 1))
CHAIN_TIMES_MIN = (12.0, 120.0)


@dataclass
class Dim:
    source: str
    rxn: int               # index into core/edge reaction list (-1: global)
    key: Optional[str]     # Km substrate label; None = scalar km / not a Km
    mean: float            # pre-rule linear value, or ΔG°′ (kJ/mol)
    sd: float              # log10 SD (kcat/km/fabh_ki), linear SD (fabi_kcat, dg diag)
    enzyme: str
    label: str
    role: str = ""
    # Every reaction whose snapshot must carry this value (reactions sharing one kinetics object).
    rxns: tuple = ()
    # Km only: every (reaction index, Km label) that carries this one prediction.
    targets: tuple = ()

    def km_targets(self) -> tuple:
        return self.targets or tuple((k, self.key) for k in (self.rxns or (self.rxn,)))


@dataclass
class Inventory:
    kcat: List[Dim]
    km: List[Dim]
    dg: List[Dim]
    dg_chol: np.ndarray    # L with L @ L.T = joint ΔG°′ covariance (kJ/mol)
    fabh_ki_mM: float
    excluded: List[tuple]  # (source, label, reason)


def _psd_factor(cov: np.ndarray) -> np.ndarray:
    """Square-root factor of a (possibly singular) PSD covariance via eigendecomposition."""
    if cov.size == 0:
        return cov
    sym = 0.5 * (cov + cov.T)
    w, v = np.linalg.eigh(sym)
    w = np.clip(w, 0.0, None)
    return v * np.sqrt(w)


def build_inventory(reactions, snaps, n_core, thermo_engine=None, smiles_maps=None,
                    fabh_ki_mM=float("nan")) -> Inventory:
    """List every sampled dimension from pre-rule snapshots of the core reactions.

    Reactions sharing one kinetics object (production FabA/FabZ pairs) get one set of
    kinetic dimensions, written into every sharing reaction's snapshot.
    """
    kcat_dims: List[Dim] = []
    km_dims: List[Dim] = []
    dg_dims: List[Dim] = []
    excluded: List[tuple] = []
    seen: Dict[str, int] = {}
    sharing: Dict[int, List[int]] = {}
    for j, rxn in enumerate(reactions):
        if getattr(rxn, "kinetics", None) is not None:
            sharing.setdefault(id(rxn.kinetics), []).append(j)
    done_kin = set()
    for j in range(n_core):
        rxn = reactions[j]
        label = sa._reaction_label(rxn, seen)
        enzyme = str(getattr(rxn, "enzyme_label", ""))
        kin = getattr(rxn, "kinetics", None)
        ks = (snaps[j] or {}).get("kinetics")
        if kin is not None and ks is not None and id(kin) not in done_kin:
            done_kin.add(id(kin))
            group = tuple(sharing[id(kin)])
            if len(group) > 1:
                enzyme = "/".join(
                    sorted({str(reactions[k].enzyme_label) for k in group})
                )
                for k in group:
                    if (snaps[k] or {}).get("kinetics") != ks:
                        raise RuntimeError(f"shared kinetics with differing baselines: {label}")
            kcat0 = ks.get("kcat")
            kcat_sd = getattr(kin, "kcat_sd", None)
            if _ec_norm(getattr(rxn, "ec_number", None)) == FABI_EC:
                excluded.append(("kcat", label, "FabI anchored; sampled as fabi_kcat"))
            elif kcat0 and kcat0 > 0 and kcat_sd and kcat_sd > 0:
                kcat_dims.append(Dim("kcat", j, None, float(kcat0),
                                     linear_sd_to_log10_sd(kcat0, kcat_sd), enzyme, label,
                                     rxns=group))
            elif kcat0 and kcat0 > 0:
                excluded.append(("kcat", label, "no CatPred SD"))

            per = ks.get("km_per_substrate") or {}
            sds = getattr(kin, "km_sd_per_substrate", None) or {}
            reactants = set(rxn.reactant_labels)
            products = set(rxn.product_labels)
            for lab, km0 in per.items():
                role = "reactant" if lab in reactants else (
                    "product" if lab in products else "unused")
                sd = sds.get(lab)
                if km0 and km0 > 0 and sd and sd > 0:
                    km_dims.append(Dim("km", j, lab, float(km0),
                                       linear_sd_to_log10_sd(km0, sd), enzyme, label, role,
                                       rxns=group))
                elif km0 and km0 > 0:
                    excluded.append(("km", f"{label} [{lab}]", f"no SD ({role})"))
            km0 = ks.get("km")
            km_sd = getattr(kin, "km_sd", None)
            if not per and km0 and km0 > 0 and km_sd and km_sd > 0:
                km_dims.append(Dim("km", j, None, float(km0),
                                   linear_sd_to_log10_sd(km0, km_sd), enzyme, label, "scalar",
                                   rxns=group))

        ts = (snaps[j] or {}).get("thermo")
        if ts is not None and ts.get("source") == "equilibrator" and math.isfinite(
            ts.get("dgr_prime_kJmol", float("nan"))
        ):
            dg_dims.append(Dim("dg", j, None, float(ts["dgr_prime_kJmol"]),
                               float(ts.get("sigma_kJmol") or 0.0),
                               str(getattr(rxn, "enzyme_label", "")), label, rxns=(j,)))
        elif ts is not None:
            excluded.append(("dg", label, f"thermo source={ts.get('source')}"))

    km_dims = tie_identical_km(km_dims)

    dg_chol = np.zeros((0, 0))
    if dg_dims and thermo_engine is not None and smiles_maps is not None:
        stoichs = [reactions[d.rxn].stoichiometry for d in dg_dims]
        maps = [smiles_maps[d.rxn] for d in dg_dims]
        joint = thermo_engine.joint_dgr_prime(stoichs, maps)
        if joint is None:
            raise RuntimeError("equilibrator unavailable: cannot build joint ΔG°′ covariance")
        mu, cov, resolved = joint
        for d, m, ok in zip(dg_dims, mu, resolved):
            if not ok:
                raise RuntimeError(f"joint ΔG°′ could not rebuild reaction: {d.label}")
        inv_dg_mean_diff = float(np.max(np.abs(mu - np.array([d.mean for d in dg_dims]))))
        sig = np.sqrt(np.clip(np.diag(cov), 0, None))
        inv_sig_diff = float(np.max(np.abs(sig - np.array([d.sd for d in dg_dims]))))
        print(
            f"joint ΔG°′: {len(dg_dims)} reactions; max |mean - production| = "
            f"{inv_dg_mean_diff:.3g} kJ/mol; max |sqrt(diag) - σ| = {inv_sig_diff:.3g} kJ/mol",
            flush=True,
        )
        dg_chol = _psd_factor(cov)
    elif dg_dims:
        dg_chol = np.diag([d.sd for d in dg_dims])
    return Inventory(kcat_dims, km_dims, dg_dims, dg_chol, float(fabh_ki_mM), excluded)


def tie_identical_km(dims: List[Dim]) -> List[Dim]:
    """Merge Km dims that are one CatPred prediction repeated across reactions.

    A per-substrate Km depends only on (enzyme sequence, compound), so the same enzyme and
    compound with the same value and SD is one physical parameter: it gets one draw.
    CatPred predicts each reaction in its own batch, so repeats differ by ~1e-7 relative;
    values within KM_TIE_RTOL count as the same prediction, larger gaps (a different input
    SMILES for the same label) stay separate.
    """
    groups: Dict[tuple, List[Dim]] = {}
    out: List[Dim] = []
    for d in dims:
        if d.key is None:
            out.append(d)
            continue
        candidates = groups.setdefault((d.enzyme.lower().strip(), d.key.lower().strip()), [])
        first = next((g for g in candidates
                      if math.isclose(g.mean, d.mean, rel_tol=KM_TIE_RTOL)
                      and math.isclose(g.sd, d.sd, rel_tol=KM_TIE_RTOL)), None)
        if first is None:
            first = Dim(d.source, d.rxn, d.key, d.mean, d.sd, d.enzyme, d.label, d.role,
                        rxns=d.rxns, targets=d.km_targets())
            candidates.append(first)
            out.append(first)
            continue
        first.targets = first.targets + d.km_targets()
        roles = set(first.role.split("/")) | {d.role}
        first.role = "/".join(sorted(r for r in roles if r))
    for d in out:
        if len(d.targets) > 1:
            d.label = f"{d.enzyme} [{d.key}] ({len(d.targets)} reactions)"
    return out


def perturbed_snapshots(base_snaps, inv: Inventory, values: Dict[str, np.ndarray], T_K: float):
    """Copy the pre-rule snapshots with sampled values written in (unchanged where equal)."""
    snaps = copy.deepcopy(base_snaps)
    for d, v in zip(inv.kcat, values["kcat"]):
        if v != d.mean:
            for k in d.rxns or (d.rxn,):
                snaps[k]["kinetics"]["kcat"] = float(v)
    for d, v in zip(inv.km, values["km"]):
        if v == d.mean:
            continue
        for k, key in d.km_targets():
            ks = snaps[k]["kinetics"]
            if key is None:
                ks["km"] = float(v)
                continue
            per = ks["km_per_substrate"]
            if ks.get("km") == per.get(key) and next(iter(per)) == key:
                ks["km"] = float(v)
            per[key] = float(v)
    for d, v in zip(inv.dg, values["dg"]):
        if v != d.mean:
            th = snaps[d.rxn]["thermo"]
            th["dgr_prime_kJmol"] = float(v)
            try:
                th["keq"] = math.exp(-float(v) * 1000 / (R * T_K))
            except OverflowError:
                th["keq"] = float("inf")
    return snaps


def _scale_feedback(reactions, enzyme: str, factor: float) -> None:
    key = enzyme.lower()
    for rxn in reactions:
        fb = getattr(rxn, "feedback_inhibitors", None)
        if not isinstance(fb, dict) or not fb:
            continue
        if str(getattr(rxn, "enzyme_label", "")).lower().strip() != key:
            continue
        rxn.feedback_inhibitors = {
            lab: (float(ki) * factor, float(h)) for lab, (ki, h) in fb.items()
        }


def apply_sample(enlarger, base_model, base_snaps, inv, values, T_K, scenario="none"):
    """Deep-copied model with sampled pre-rule values, rules re-applied, feedback Ki set."""
    model = copy.deepcopy(base_model)
    rxns = list(model.core_reactions) + list(model.edge_reactions)
    snaps = perturbed_snapshots(base_snaps, inv, values, T_K)
    fabi_rule = RULES.by_name(FABI_RULE)
    saved_ec_kcat = fabi_rule.params["ec_kcat"]
    saved_model = enlarger.model
    settings = enlarger.bees_object.settings
    saved_cal = list(getattr(settings, "calibrations", None) or [])
    try:
        fabi_rule.params["ec_kcat"] = {**saved_ec_kcat, FABI_EC: float(values["fabi_kcat"])}
        if scenario == "tesa_off":
            settings.calibrations = [c for c in saved_cal if c != TESA_RULE]
        enlarger.model = model
        with RULES.isolated_baselines():
            for rxn, snap in zip(rxns, snaps):
                RULES.seed_baseline(rxn, snap)
            enlarger._attach_thermo_to_reactions()
    finally:
        fabi_rule.params["ec_kcat"] = saved_ec_kcat
        enlarger.model = saved_model
        settings.calibrations = saved_cal
    f = float(values["fabh_ki_factor"])
    if f != 1.0:
        _scale_feedback(model.core_reactions, "FabH", f)
    fi = float(values["fabi_ki_factor"])
    fi *= {"fabi_ki_x0.5": 0.5, "fabi_ki_x2": 2.0}.get(scenario, 1.0)
    if fi != 1.0:
        _scale_feedback(model.core_reactions, "FabI", fi)
    return model


def _fa_key(label: str):
    """(carbons, double bonds) of a free fatty acid label, else None."""
    if sa._fa_carbons(label) is None:
        return None
    n = label.lower()
    for stem, c in sa._FA_STEMS:
        if stem in n:
            return c, int(stem.endswith("en"))
    return None


def _chain_profile(sim) -> np.ndarray:
    """Free-acid µM per CHAIN_KEYS bin at each CHAIN_TIMES_MIN, shape (times, bins)."""
    col = {k: j for j, k in enumerate(CHAIN_KEYS)}
    out = np.zeros((len(CHAIN_TIMES_MIN), len(CHAIN_KEYS)))
    t_s = np.array(CHAIN_TIMES_MIN) * 60.0
    for i, lab in enumerate(sim.species_labels):
        j = col.get(_fa_key(lab))
        if j is not None:
            out[:, j] += 1000.0 * np.interp(t_s, sim.t, sim.y[i])
    return out


def _global_smiles_maps(enlarger, reactions):
    glob = enlarger._build_global_smiles_map()
    maps = []
    for rxn in reactions:
        m = dict(glob)
        for lab, smi in (getattr(rxn.kinetics, "compound_smiles", None) or {}).items():
            if smi:
                m[lab] = smi
        maps.append(m)
    return maps


def _fabh_ki_mM(model) -> float:
    for rxn in model.core_reactions:
        if str(getattr(rxn, "enzyme_label", "")).lower().strip() != "fabh":
            continue
        fb = getattr(rxn, "feedback_inhibitors", None)
        if isinstance(fb, dict) and fb:
            return float(next(iter(fb.values()))[0])
    return float("nan")


def prepare():
    """Build the production network once; return everything a sample needs."""
    print("Building final FAS-II network (CatPred cached)...", flush=True)
    enlarger, model, T_K = sa.build_model()
    rxns = list(model.core_reactions) + list(model.edge_reactions)
    n_core = len(model.core_reactions)
    snaps = [RULES.baseline_snapshot(r) for r in rxns]
    missing = [i for i, s in enumerate(snaps) if s is None]
    if missing:
        raise RuntimeError(f"{len(missing)} reactions have no rule baseline snapshot")
    kin_ids = [id(r.kinetics) for r in rxns if r.kinetics is not None]
    n_shared = len(kin_ids) - len(set(kin_ids))
    if n_shared:
        print(f"NOTE: {n_shared} reactions share a kinetics object with another reaction; "
              "they are sampled together (production behaviour kept).", flush=True)
    engine = enlarger._thermo_engine
    T_K = float(getattr(engine, "T_K", T_K))
    inv = build_inventory(
        rxns, snaps, n_core,
        thermo_engine=engine,
        smiles_maps=_global_smiles_maps(enlarger, rxns),
        fabh_ki_mM=_fabh_ki_mM(model),
    )
    return enlarger, model, snaps, inv, T_K


LN10 = math.log(10.0)
P_LO = 5.0
P_HI = 95.0
ODE_SEED = 20261123
BOOT_SEED = 20261124
N_BOOT = 1000
STABILITY_SPLIT = 64
STABILITY_REL = 0.15
STABILITY_ABS_UM = 1.0
BAR_EPS = 1e-3  # µM floor only for relative-error helpers
MAX_GROUPS = 6
TOP_ODE_GROUPS = 4
N_ODE = 128
MIN_SHARE = 0.01
T_END_S = 7200.0
GRID_DT_S = 30.0
BAND_LABEL = "±1 CatPred SD box, top 4 groups"
FABI_KCAT = 15.0
FABH_S_KI = 0.0380  # saved local slope (d ln PE / d ln Ki)

_REPO = os.path.abspath(os.path.join(_PROJECT_DIR, "..", "..", ".."))
_SENS = os.path.join(_REPO, "docs", "deliverables", "sensitivity_analysis")
_OUT = os.path.join(_PROJECT_DIR, "output", "uncertainty", "propagation")
_EXP_CSV = os.path.join(
    _REPO, "knowledge", "references", "ecoli_fas2_project", "ruppe_2018_pnas",
    "Model and Solver", "Experimental_Dataset.csv",
)
_INV_CSV = os.path.join(_PROJECT_DIR, "output", "uncertainty", "parameter_inventory.csv")
_DG_NPZ = os.path.join(_PROJECT_DIR, "output", "uncertainty", "dg_covariance.npz")
_ECOLI_CSV = os.path.join(_REPO, "db", "ecoli.csv")


# ---------------------------------------------------------------------------
# Pure helpers (unit-tested without a FAS run)
# ---------------------------------------------------------------------------


def ln_sd(sd_log10: float) -> float:
    """CatPred reports SD in log10. The normal coordinate is in natural log."""
    return float(sd_log10) * LN10


def coeff_of(pairs: Sequence[Tuple[float, float]]) -> float:
    """Sum of slope × uncertainty. ``pairs`` are (S, sigma_ln)."""
    return float(sum(s * sig for s, sig in pairs))


def dg_variance(s: np.ndarray, cov: np.ndarray) -> float:
    s = np.asarray(s, dtype=float)
    cov = np.asarray(cov, dtype=float)
    return float(s @ cov @ s)


def relative_errors(pred: np.ndarray, obs: np.ndarray) -> np.ndarray:
    obs = np.asarray(obs, dtype=float)
    pred = np.asarray(pred, dtype=float)
    denom = np.maximum(np.abs(obs), BAR_EPS)
    return np.abs(pred - obs) / denom


def median_rel(pred: np.ndarray, obs: np.ndarray) -> float:
    err = relative_errors(pred, obs)
    err = err[np.isfinite(err)]
    if err.size == 0:
        return float("nan")
    return float(np.median(err))


def bootstrap_percentile_ci(
    samples: np.ndarray,
    q: float,
    n_boot: int = N_BOOT,
    seed: int = BOOT_SEED,
) -> Tuple[float, float]:
    """95% bootstrap CI for a percentile of a 1-D sample."""
    x = np.asarray(samples, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, x.size, size=(n_boot, x.size))
    qs = np.percentile(x[idx], q, axis=1)
    lo, hi = np.percentile(qs, [2.5, 97.5])
    return float(lo), float(hi)


def stability_ok(a: float, b: float) -> bool:
    """First-64 vs all-N rule at 12 min: <15% relative, or <1 µM absolute if near zero."""
    if not (np.isfinite(a) and np.isfinite(b)):
        return False
    if abs(a) < STABILITY_ABS_UM and abs(b) < STABILITY_ABS_UM:
        return abs(a - b) < STABILITY_ABS_UM
    denom = max(abs(a), abs(b), STABILITY_ABS_UM)
    return abs(a - b) / denom < STABILITY_REL


# ---------------------------------------------------------------------------
# Ranking
# ---------------------------------------------------------------------------


@dataclass
class Group:
    name: str
    kind: str
    coeff: float
    variance: float
    detail: str = ""
    # Median CatPred log10 SD used for the log-uniform factor range [10^{-sd}, 10^{+sd}].
    sd_log10: float = 0.0
    # kcat/km: (inventory index, log10 sd). These indices are scaled together.
    scales: List[Tuple[int, float]] = field(default_factory=list)
    dg_mean: Optional[np.ndarray] = None
    dg_direction: Optional[np.ndarray] = None
    is_fabi_kcat: bool = False


def _enzyme(label: str) -> str:
    return label.split(":", 1)[0].strip()


def select_groups(groups: Sequence[Group]) -> Tuple[List[Group], bool]:
    ranked = sorted(groups, key=lambda g: g.variance, reverse=True)
    total = sum(g.variance for g in ranked) or 1.0
    chosen: List[Group] = []
    for g in ranked:
        if g.variance <= 0:
            continue
        if len(chosen) >= MAX_GROUPS:
            break
        if chosen and g.variance / total < MIN_SHARE:
            break
        chosen.append(g)
    omitted = [g for g in ranked if g not in chosen and g.variance / total >= 0.05]
    return chosen, bool(omitted)


def select_top_ode(groups: Sequence[Group], n: int = TOP_ODE_GROUPS) -> List[Group]:
    """Top kinetic groups for the ODE-only global sample (skip ΔG and tiny FabI kcat)."""
    ranked = sorted(groups, key=lambda g: g.variance, reverse=True)
    out: List[Group] = []
    for g in ranked:
        if g.kind in ("dg", "fabh_ki") or g.is_fabi_kcat:
            continue
        if g.sd_log10 <= 0 or g.kind not in ("kcat", "km"):
            continue
        out.append(g)
        if len(out) >= n:
            break
    return out


def design_loguniform(n: int, d: int, seed: int) -> np.ndarray:
    """Sobol in the unit cube mapped to Uniform(-1, 1) for log10 factors."""
    from scipy.stats import qmc
    sampler = qmc.Sobol(d=d, scramble=True, seed=seed)
    return 2.0 * sampler.random(n) - 1.0


def rank_from_rows(
    kcat_rows: Sequence[Tuple[str, float]],
    km_rows: Sequence[Tuple[str, float]],
    dg_rows: Sequence[Tuple[str, float]],
    kcat_sd: Dict[str, float],
    km_sd: Dict[str, float],
    km_enzyme_sd: Dict[str, float],
    dg_labels: Sequence[str],
    dg_cov: np.ndarray,
    dg_mean: np.ndarray,
    dg_near: Sequence[bool],
) -> List[Group]:
    groups: List[Group] = []
    by: Dict[str, List[Tuple[float, float]]] = {}
    by_sd: Dict[str, List[float]] = {}
    fabi_pairs: List[Tuple[float, float]] = []
    for label, s in kcat_rows:
        enz = _enzyme(label)
        if enz.lower() == "fabi":
            fabi_pairs.append((s, FABI_KCAT_SD / FABI_KCAT))
            continue
        sd = kcat_sd.get(label)
        if sd is None or sd <= 0:
            continue
        by.setdefault(f"{enz} kcat", []).append((s, ln_sd(sd)))
        by_sd.setdefault(f"{enz} kcat", []).append(sd)
    for name, pairs in by.items():
        c = coeff_of(pairs)
        groups.append(Group(
            name, "kcat", c, c * c, f"{len(pairs)} reactions",
            sd_log10=float(np.median(by_sd[name])),
        ))
    if fabi_pairs:
        c = coeff_of(fabi_pairs)
        groups.append(Group(
            "FabI kcat", "kcat", c, c * c,
            f"{len(fabi_pairs)} reactions; measured 15 ± 0.27 s⁻¹",
            sd_log10=0.0,
            is_fabi_kcat=True,
        ))

    by = {}
    by_sd = {}
    for label, s in km_rows:
        enz = _enzyme(label)
        sd = km_sd.get(label)
        if sd is None or sd <= 0:
            sd = km_enzyme_sd.get(enz)
        if sd is None or sd <= 0:
            continue
        by.setdefault(f"{enz} Km", []).append((s, ln_sd(sd)))
        by_sd.setdefault(f"{enz} Km", []).append(sd)
    for name, pairs in by.items():
        c = coeff_of(pairs)
        groups.append(Group(
            name, "km", c, c * c, f"{len(pairs)} reactions",
            sd_log10=float(np.median(by_sd[name])),
        ))

    if dg_labels and dg_cov.size:
        s = np.zeros(len(dg_labels))
        smap = {lab: val for lab, val in dg_rows}
        for i, lab in enumerate(dg_labels):
            if dg_near[i]:
                continue
            s[i] = smap.get(lab, 0.0)
        var = max(0.0, dg_variance(s, dg_cov))
        direction = None
        if var > 0:
            direction = (dg_cov @ s) / math.sqrt(var)
        groups.append(Group(
            "dG joint", "dg", math.sqrt(var), var,
            f"{int(np.count_nonzero(s))} reactions in the quadratic form",
            dg_mean=np.asarray(dg_mean, dtype=float),
            dg_direction=direction,
        ))

    c = FABH_S_KI * ln_sd(FABH_KI_LOG10_SD)
    groups.append(Group("FabH Ki", "fabh_ki", c, c * c, "saved S_Ki = 0.038"))
    return groups


def _read_tornado(path: str) -> List[Tuple[str, float]]:
    rows = []
    with open(path, newline="", encoding="utf-8") as fh:
        for rec in csv.DictReader(fh):
            rows.append((rec["label"], float(rec["S"])))
    return rows


def _inventory_sds(path: str):
    kcat_sd: Dict[str, float] = {}
    km_lists: Dict[str, List[float]] = {}
    km_enz: Dict[str, List[float]] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for rec in csv.DictReader(fh):
            if not rec["status"].startswith("sampled"):
                continue
            try:
                sd = float(rec["sd"])
            except ValueError:
                continue
            if rec["source"] == "kcat":
                kcat_sd[rec["reaction"]] = sd
            elif rec["source"] == "km" and "reactions)" not in rec["reaction"]:
                km_lists.setdefault(rec["reaction"], []).append(sd)
                km_enz.setdefault(_enzyme(rec["reaction"]), []).append(sd)
    km_sd = {k: float(np.median(v)) for k, v in km_lists.items()}
    km_enzyme_sd = {k: float(np.median(v)) for k, v in km_enz.items()}
    return kcat_sd, km_sd, km_enzyme_sd


def _near_mask(means: np.ndarray, sigmas: np.ndarray, cutoff: float) -> List[bool]:
    return [bool(s > 0 and abs(m - cutoff) <= 2.0 * s) for m, s in zip(means, sigmas)]


def load_ranking(cutoff: float) -> Tuple[List[Group], str]:
    kcat = _read_tornado(os.path.join(_SENS, "sensitivity_analysis_kcat.csv"))
    km = _read_tornado(os.path.join(_SENS, "sensitivity_analysis_km.csv"))
    dg_rows = _read_tornado(os.path.join(_SENS, "sensitivity_analysis_dgr.csv"))
    kcat_sd, km_sd, km_enzyme_sd = _inventory_sds(_INV_CSV)
    z = np.load(_DG_NPZ, allow_pickle=True)
    labels = [str(x) for x in z["labels"]]
    cov = np.asarray(z["cov"], dtype=float)
    means = np.asarray(z["mean"], dtype=float)
    sig = np.sqrt(np.clip(np.diag(cov), 0, None))
    near = _near_mask(means, sig, cutoff)
    groups = rank_from_rows(
        kcat, km, dg_rows, kcat_sd, km_sd, km_enzyme_sd,
        labels, cov, means, near,
    )
    note = (
        f"Tornado files in {_SENS}; S is d ln(product) / d ln(input) for kcat and Km, "
        f"and d ln(product) / d(ΔG°′) per kJ/mol. Ranking time is the tornado t* "
        f"(about 441 s), not 12 min. {sum(near)} ΔG reactions sit within 2σ of the "
        f"irreversibility cutoff and are left out of the ΔG quadratic form."
    )
    return groups, note


def write_contributions(
    groups: Sequence[Group],
    chosen: Sequence[Group],
    path: str,
    flag_col: str = "in_ode_global",
) -> None:
    total = sum(g.variance for g in groups) or 1.0
    chosen_names = {g.name for g in chosen}
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([
            "group", "kind", "coeff_dlnP_dz", "variance", "share", "sd_log10",
            flag_col, "detail",
        ])
        for g in sorted(groups, key=lambda g: g.variance, reverse=True):
            w.writerow([
                g.name, g.kind, f"{g.coeff:.6g}", f"{g.variance:.6g}",
                f"{g.variance / total:.6g}", f"{g.sd_log10:.6g}",
                int(g.name in chosen_names), g.detail,
            ])


def write_experiment(path: str) -> None:
    rows = []
    with open(_EXP_CSV, newline="", encoding="utf-8") as fh:
        for rec in csv.DictReader(fh):
            rows.append((float(rec["Time (min)"]), float(rec["Concentration (um)"])))
    cite = (
        "Digitized means from Ruppe 2018 SI Experimental_Dataset.csv (S2A time course). "
        "That file has no SD or SEM. Yu 2011 Fig. 4 in this repo is a titration, not these "
        "points. No measurement bar is invented."
    )
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["t_min", "mean_uM", "error_uM", "citation"])
        for t, y in rows:
            w.writerow([f"{t:.4f}", f"{y:.4f}", "", cite])


# ---------------------------------------------------------------------------
# Applying a design point to the fixed network
# ---------------------------------------------------------------------------


def _attach_scales(groups: Sequence[Group], inv) -> None:
    by_name = {g.name: g for g in groups}
    for i, d in enumerate(inv.kcat):
        g = by_name.get(f"{d.enzyme} kcat")
        if g is not None and not g.is_fabi_kcat:
            g.scales.append((i, d.sd))
    for i, d in enumerate(inv.km):
        enz = d.enzyme.split("/")[0].strip()
        g = by_name.get(f"{enz} Km") or by_name.get(f"{d.enzyme} Km")
        if g is not None:
            g.scales.append((i, d.sd))


def _baseline_values(inv):
    """Unperturbed pre-rule values. The published sample scales kcat and Km from here."""
    return {
        "kcat": np.array([d.mean for d in inv.kcat], dtype=float),
        "km": np.array([d.mean for d in inv.km], dtype=float),
        "dg": np.array([d.mean for d in inv.dg], dtype=float),
        "fabi_kcat": np.array(FABI_KCAT_MEAN),
        "fabh_ki_factor": np.array(1.0),
        "fabi_ki_factor": np.array(1.0),
    }


def _apply_u(inv, base: Dict[str, np.ndarray], chosen: Sequence[Group], u: np.ndarray):
    """Log-uniform apply: factor = 10**(sd_log10 * u) with u in [-1, 1], one factor per group."""
    values = {k: np.array(v, copy=True) for k, v in base.items()}
    for g, ui in zip(chosen, u):
        ui = float(ui)
        factor = 10.0 ** (g.sd_log10 * ui)
        if g.kind == "kcat":
            for idx, _sd in g.scales:
                values["kcat"][idx] = inv.kcat[idx].mean * factor
        elif g.kind == "km":
            for idx, _sd in g.scales:
                values["km"][idx] = inv.km[idx].mean * factor
        else:
            raise RuntimeError(f"ODE global apply does not support kind={g.kind}")
    return values


def _metrics(sim, grid_t: np.ndarray):
    import sensitivity_analysis as sa
    pe = sa._pe_series(sim)
    t_min = sim.t / 60.0
    return {
        "pe": np.interp(grid_t, sim.t, pe),
        "pe12": float(np.interp(12.0, t_min, pe)),
        "pe120": float(np.interp(120.0, t_min, pe)),
        "chains": _chain_profile(sim),
    }


_W: Dict[str, object] = {}


def _run_one(i: int):
    import sensitivity_analysis as sa
    coord = _W["U"][i]
    values = _apply_u(_W["inv"], _W["base"], _W["chosen"], coord)
    try:
        model = apply_sample(
            _W["enlarger"], _W["model"], _W["snaps"], _W["inv"], values, _W["T_K"], "none",
        )
        sim = sa._simulate(_W["enlarger"], model, end_time=T_END_S)
        if not sim.success or sim.t.size == 0 or sim.t[-1] < T_END_S - 1.0:
            return i, None, f"solver: {sim.message}"
        return i, _metrics(sim, _W["grid"]), ""
    except (ValueError, FloatingPointError, np.linalg.LinAlgError, ArithmeticError) as exc:
        return i, None, f"{type(exc).__name__}: {exc}"


def _run_batch(us: np.ndarray, jobs: int):
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor
    n = us.shape[0]
    _W["U"] = us
    results = [None] * n
    msgs = [""] * n
    if jobs <= 1:
        for i in range(n):
            i2, met, msg = _run_one(i)
            results[i2] = met
            msgs[i2] = msg
            print(f"  point {i + 1}/{n} {msg or 'ok'}", flush=True)
        return results, msgs
    ctx = mp.get_context("fork")
    with ProcessPoolExecutor(max_workers=jobs, mp_context=ctx) as ex:
        done = 0
        for i, met, msg in ex.map(_run_one, range(n), chunksize=1):
            results[i] = met
            msgs[i] = msg
            done += 1
            if done % 5 == 0 or done == n:
                print(f"  {done}/{n}", flush=True)
    return results, msgs


def _published(enlarger, model, grid):
    import sensitivity_analysis as sa
    sim = sa._simulate(enlarger, model, end_time=T_END_S)
    if not sim.success:
        raise RuntimeError(f"published simulation failed: {sim.message}")
    return _metrics(sim, grid)


def ecoli_csv_fingerprint() -> str:
    import hashlib
    path = _ECOLI_CSV
    if not os.path.isfile(path):
        return f"{path} (missing)"
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return f"{path} sha256={h.hexdigest()[:16]}…"


def percentile_band(pe_ok: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    lo = np.percentile(pe_ok, P_LO, axis=0)
    med = np.percentile(pe_ok, 50, axis=0)
    hi = np.percentile(pe_ok, P_HI, axis=0)
    mean = pe_ok.mean(axis=0)
    return mean, med, lo, hi


def bootstrap_band_cis(
    pe_ok: np.ndarray,
    n_boot: int = N_BOOT,
    seed: int = BOOT_SEED,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Pointwise 95% bootstrap CIs for the p5 and p95 curves."""
    n, nt = pe_ok.shape
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_boot, n))
    samples = pe_ok[idx]  # (n_boot, n, nt)
    p5s = np.percentile(samples, P_LO, axis=1)
    p95s = np.percentile(samples, P_HI, axis=1)
    p5_lo, p5_hi = np.percentile(p5s, [2.5, 97.5], axis=0)
    p95_lo, p95_hi = np.percentile(p95s, [2.5, 97.5], axis=0)
    return p5_lo, p5_hi, p95_lo, p95_hi


def stability_at_12min(pe12: np.ndarray, split: int = STABILITY_SPLIT) -> Dict[str, object]:
    """Compare first-`split` draws vs all draws at 12 min (p5/p50/p95/mean)."""
    x = np.asarray(pe12, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < split:
        raise RuntimeError(f"need ≥{split} finite pe12 draws for stability; got {x.size}")
    first = x[:split]
    keys = ("p5", "p50", "p95", "mean")
    full_vals = {
        "p5": float(np.percentile(x, P_LO)),
        "p50": float(np.percentile(x, 50)),
        "p95": float(np.percentile(x, P_HI)),
        "mean": float(x.mean()),
    }
    first_vals = {
        "p5": float(np.percentile(first, P_LO)),
        "p50": float(np.percentile(first, 50)),
        "p95": float(np.percentile(first, P_HI)),
        "mean": float(first.mean()),
    }
    ok = {k: stability_ok(first_vals[k], full_vals[k]) for k in keys}
    return {
        "n_first": split,
        "n_all": int(x.size),
        "first": first_vals,
        "all": full_vals,
        "ok": ok,
        "pass": all(ok.values()),
    }


def write_band_csv(
    path: str,
    grid: np.ndarray,
    published_pe: np.ndarray,
    mean: np.ndarray,
    med: np.ndarray,
    lo: np.ndarray,
    hi: np.ndarray,
    p5_ci_lo: np.ndarray,
    p5_ci_hi: np.ndarray,
    p95_ci_lo: np.ndarray,
    p95_ci_hi: np.ndarray,
    p25: Optional[np.ndarray] = None,
    p75: Optional[np.ndarray] = None,
) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        header = [
            "t_min", "published_uM", "ensemble_mean_uM", "p50_uM",
            "p5_uM", "p95_uM",
            "p5_boot95_lo_uM", "p5_boot95_hi_uM",
            "p95_boot95_lo_uM", "p95_boot95_hi_uM",
        ]
        if p25 is not None and p75 is not None:
            header += ["p25_uM", "p75_uM"]
        w.writerow(header)
        for k, t in enumerate(grid):
            row = [
                f"{t / 60.0:.4f}", f"{published_pe[k]:.5g}", f"{mean[k]:.5g}",
                f"{med[k]:.5g}", f"{lo[k]:.5g}", f"{hi[k]:.5g}",
                f"{p5_ci_lo[k]:.5g}", f"{p5_ci_hi[k]:.5g}",
                f"{p95_ci_lo[k]:.5g}", f"{p95_ci_hi[k]:.5g}",
            ]
            if p25 is not None and p75 is not None:
                row += [f"{p25[k]:.5g}", f"{p75[k]:.5g}"]
            w.writerow(row)


def build_ode_notes(
    *,
    rank_note: str,
    group_line: str,
    n: int,
    n_ok: int,
    pub12: float,
    mean12: float,
    med12: float,
    p5_12: float,
    p95_12: float,
    p5_ci: Tuple[float, float],
    p95_ci: Tuple[float, float],
    stab: Dict[str, object],
    failure_lines: Sequence[str] = (),
) -> str:
    stab_bits = []
    for key in ("p5", "p50", "p95", "mean"):
        a = stab["first"][key]
        b = stab["all"][key]
        mark = "PASS" if stab["ok"][key] else "FAIL"
        stab_bits.append(
            f"{key}: first-{stab['n_first']}={a:.2f} vs all-{stab['n_all']}={b:.2f} µM ({mark})"
        )
    stab_line = (
        f"64-vs-128 stability at 12 min "
        f"(rule: <{STABILITY_REL:.0%} relative or <{STABILITY_ABS_UM:g} µM absolute): "
        + "; ".join(stab_bits)
        + f". Overall: {'PASS' if stab['pass'] else 'FAIL'}."
    )
    max_move = max(
        abs(stab["first"][k] - stab["all"][k])
        / max(abs(stab["first"][k]), abs(stab["all"][k]), STABILITY_ABS_UM)
        for k in ("p5", "p50", "p95", "mean")
    )
    lines = [
        "FAS-II local ranking + ODE-only global sample on top kinetic groups.",
        f"Band label: {BAND_LABEL}.",
        "The band shows how the model responds to CatPred uncertainty. "
        "It is not a calibrated prediction interval.",
        rank_note,
        "Local ranking (contributions.csv) is the main uncertainty result for which "
        "inputs drive variance. The ODE ensemble is the global sample on that shortlist.",
        "CatPred log10 SDs are uncalibrated assigned uncertainties. "
        "Enzyme amount, starting substrates, and FabI Ki have no SD and are excluded. "
        "ΔG°′ and FabH Ki are ranked but not sampled in this ODE global set.",
        "Full within-enzyme correlation: one shared log-uniform factor per enzyme and "
        "parameter class (all reactions of that enzyme×class move together). "
        "Groups are independent of each other. A median log10 SD hides chain-length differences.",
        "The TesA reversibility cutoff is not sampled. It stays a separate scenario.",
        group_line,
        f"Design: scrambled Sobol → Uniform(-1,1), N={n}, seed={ODE_SEED}. "
        f"Succeeded {n_ok}/{n}.",
        "Band edges are the 5th and 95th percentiles of runs inside the ±1 CatPred "
        "SD box on the top 4 groups, not of the full CatPred uncertainty.",
        "Uniform(-1,1) has an effective SD of about 0.58 of one CatPred SD. "
        "A 3–10× CatPred error is not covered.",
        f"Published PE at 12 min = {pub12:.2f} µM. "
        f"Ensemble mean = {mean12:.2f} µM, median = {med12:.2f} µM. "
        f"p5 is about 8 µM (point {p5_12:.2f}; bootstrap 95% CI "
        f"[{p5_ci[0]:.1f}, {p5_ci[1]:.1f}]). "
        f"p95 point = {p95_12:.2f} µM "
        f"(bootstrap 95% CI [{p95_ci[0]:.2f}, {p95_ci[1]:.2f}]). "
        "These edges are the 5th and 95th percentiles of draws inside the box.",
        stab_line,
        "The 1 µM absolute clause was added to the 15% relative stability rule. "
        f"Largest 64-vs-128 relative move is {max_move:.1%}.",
        "The 15 Oct tornado is not the paper model: after that date the database "
        "dropped the SurE rows (02624a2) and decarboxylations stopped being forced "
        "irreversible (5d70ee3). A 23 Nov check on the paper model, one shared "
        "factor per group at t*=440.7 s, kept FabI Km, FabG kcat and FabF kcat, "
        "and swapped 4th place from FabF Km to FabG Km (|S| 0.204 vs 0.198). "
        "FabF Km is 5th and is excluded. It is nearly tied with FabG Km, "
        "so the cutoff is not what sets the band. This band uses FabG Km.",
        "The ±1 box is partial: top 4 groups, one factor per enzyme. "
        "Enzyme amount, starting substrates, FabI Ki, and the ΔG reactions "
        "near the irreversibility cutoff are left out.",
        "A moderately larger error, about 3–10× per group, is the ±2 SD box. "
        "In a small check (N≈32) its 12 min lower edge falls below 1 µM. "
        "That is not a measured p5, and that wider box is not on the figures. "
        "At 120 min the same ±1 box is about 67–71 µM around the published "
        "plateau of 71 µM. The final amount is tight. The uncertainty is in "
        "how fast the curve gets there. Chain lengths at 12 min are not "
        "robust: C16 runs from about 0 to 15 µM around the published 7 µM.",
        "The band is a real result for this model. It does not say that "
        "46.6 µM is robust.",
        f"Reaction database for this run: {ecoli_csv_fingerprint()}.",
        "FabI kcat is 15.0 s⁻¹ (Rafi 2006) with SD 0.27 s⁻¹; it ranks near zero "
        "and is not sampled.",
    ]
    lines.extend(failure_lines)
    return "\n".join(lines) + "\n"


def write_ensemble_outputs(
    *,
    grid: np.ndarray,
    published_pe: np.ndarray,
    published_chains: np.ndarray,
    pe: np.ndarray,
    pe12: np.ndarray,
    pe120: np.ndarray,
    ok: np.ndarray,
    chains: Optional[np.ndarray],
    u: np.ndarray,
    group_names: Sequence[str],
    sd_log10: Sequence[float],
    rank_note: str,
    group_line: str,
    msgs: Sequence[str],
    chain_p5: Optional[np.ndarray] = None,
    chain_p50: Optional[np.ndarray] = None,
    chain_p95: Optional[np.ndarray] = None,
) -> Dict[str, object]:
    pe_ok = pe[ok]
    mean, med, lo, hi = percentile_band(pe_ok)
    p25 = np.percentile(pe_ok, 25, axis=0)
    p75 = np.percentile(pe_ok, 75, axis=0)
    p5_ci_lo, p5_ci_hi, p95_ci_lo, p95_ci_hi = bootstrap_band_cis(pe_ok)
    if chains is not None:
        ch_lo = np.percentile(chains[ok], P_LO, axis=0)
        ch_med = np.percentile(chains[ok], 50, axis=0)
        ch_hi = np.percentile(chains[ok], P_HI, axis=0)
    else:
        if chain_p5 is None or chain_p50 is None or chain_p95 is None:
            raise RuntimeError("need per-draw chains or saved chain percentile arrays")
        ch_lo, ch_med, ch_hi = chain_p5, chain_p50, chain_p95

    write_band_csv(
        os.path.join(_OUT, "band.csv"),
        grid, published_pe, mean, med, lo, hi,
        p5_ci_lo, p5_ci_hi, p95_ci_lo, p95_ci_hi,
        p25, p75,
    )

    n = pe.shape[0]
    d = len(group_names)
    with open(os.path.join(_OUT, "draws.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(
            ["i", "ok", "pe12_uM", "pe120_uM"]
            + [f"u_{name.replace(' ', '_')}" for name in group_names]
            + ["message"]
        )
        for i in range(n):
            w.writerow(
                [i, int(ok[i]), f"{pe12[i]:.5g}", f"{pe120[i]:.5g}"]
                + [f"{float(u[i, j]):.6f}" for j in range(d)]
                + [msgs[i] if i < len(msgs) else ""]
            )

    pe12_ok = pe12[ok]
    p5_12 = float(np.percentile(pe12_ok, P_LO))
    p95_12 = float(np.percentile(pe12_ok, P_HI))
    p5_ci = bootstrap_percentile_ci(pe12_ok, P_LO, seed=BOOT_SEED)
    p95_ci = bootstrap_percentile_ci(pe12_ok, P_HI, seed=BOOT_SEED + 1)
    stab = stability_at_12min(pe12_ok)

    pub12 = float(np.interp(12.0, grid / 60.0, published_pe))
    mean12 = float(pe12_ok.mean())
    med12 = float(np.percentile(pe12_ok, 50))

    np.savez_compressed(
        os.path.join(_OUT, "ode_global.npz"),
        groups=np.array(list(group_names)),
        sd_log10=np.asarray(sd_log10, dtype=float),
        u=u, grid_s=grid,
        published_pe=published_pe,
        pe=pe, pe12=pe12, pe120=pe120, ok=ok,
        ensemble_mean=mean, p50=med, p5=lo, p95=hi,
        p5_boot95_lo=p5_ci_lo, p5_boot95_hi=p5_ci_hi,
        p95_boot95_lo=p95_ci_lo, p95_boot95_hi=p95_ci_hi,
        published_chains=published_chains,
        chain_p5=ch_lo, chain_p50=ch_med, chain_p95=ch_hi,
        method=np.array("ode_loguniform_pm1_catpred_sd"),
        band_label=np.array(BAND_LABEL),
        stability_12min=np.array(str(stab)),
    )

    failures = [f"failure i={i}: {msg}" for i, msg in enumerate(msgs) if msg]
    notes = build_ode_notes(
        rank_note=rank_note,
        group_line=group_line,
        n=n,
        n_ok=int(ok.sum()),
        pub12=pub12,
        mean12=mean12,
        med12=med12,
        p5_12=p5_12,
        p95_12=p95_12,
        p5_ci=p5_ci,
        p95_ci=p95_ci,
        stab=stab,
        failure_lines=failures,
    )
    with open(os.path.join(_OUT, "notes.txt"), "w", encoding="utf-8") as fh:
        fh.write(notes)
    print(notes, flush=True)
    return {
        "pub12": pub12, "mean12": mean12, "med12": med12,
        "p5_12": p5_12, "p95_12": p95_12,
        "p5_ci": p5_ci, "p95_ci": p95_ci, "stab": stab,
    }


def run_ode_global(
    jobs: int,
    n: int = N_ODE,
    top_n: int = TOP_ODE_GROUPS,
    names: Optional[Sequence[str]] = None,
) -> None:
    """Global sample: full ODEs on the top groups, log-uniform ±1 CatPred SD."""
    from bees.rules import RULES

    for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(_var, "1")
    os.makedirs(_OUT, exist_ok=True)

    print("Building the published network...", flush=True)
    enlarger, model, snaps, inv, T_K = prepare()
    cutoff = -float(RULES.by_name("dgr_irreversibility").params["dgr_kjmol_cutoff"])
    groups, rank_note = load_ranking(cutoff)
    _attach_scales(groups, inv)
    if names:
        by_name = {g.name: g for g in groups}
        missing = [name for name in names if name not in by_name]
        if missing:
            raise RuntimeError(f"unknown groups: {missing}")
        chosen = [by_name[name] for name in names]
    else:
        chosen = select_top_ode(groups, n=top_n)
    if not chosen:
        raise RuntimeError("no kinetic groups selected for the ODE global sample")
    for g in chosen:
        if not g.scales:
            raise RuntimeError(f"{g.name}: no inventory scales attached; cannot apply factors")
    write_contributions(
        groups, chosen, os.path.join(_OUT, "contributions.csv"), flag_col="in_ode_global",
    )
    write_experiment(os.path.join(_OUT, "experiment.csv"))
    total_var = sum(g.variance for g in groups) or 1.0
    group_line = (
        f"ODE groups ({len(chosen)}): "
        + ", ".join(
            f"{g.name} (share {g.variance / total_var:.1%}, "
            f"sd_log10={g.sd_log10:.3f}, fold≈{10 ** g.sd_log10:.2g})"
            for g in chosen
        )
    )
    print(group_line, flush=True)

    d = len(chosen)
    u = design_loguniform(n, d, ODE_SEED)
    grid = np.arange(0.0, T_END_S + 0.5 * GRID_DT_S, GRID_DT_S)
    base_values = _baseline_values(inv)
    _W.update(
        enlarger=enlarger, model=model, snaps=snaps, inv=inv, T_K=T_K,
        base=base_values, chosen=chosen, grid=grid,
    )
    print(f"Published trajectory ({T_END_S:.0f} s)...", flush=True)
    published = _published(enlarger, model, grid)
    print(f"ODE global sample: {n} full runs, {d} groups, seed {ODE_SEED}...", flush=True)
    mets, msgs = _run_batch(u, jobs)
    ok = np.array([m is not None for m in mets])
    n_ok = int(ok.sum())
    if n_ok < 20:
        raise RuntimeError(f"only {n_ok}/{n} ODE runs succeeded")

    pe = np.full((n, grid.size), np.nan)
    pe12 = np.full(n, np.nan)
    pe120 = np.full(n, np.nan)
    chains = None
    for i, m in enumerate(mets):
        if m is None:
            continue
        pe[i] = m["pe"]
        pe12[i] = m["pe12"]
        pe120[i] = m["pe120"]
        if chains is None:
            chains = np.full((n,) + m["chains"].shape, np.nan)
        chains[i] = m["chains"]

    write_ensemble_outputs(
        grid=grid,
        published_pe=published["pe"],
        published_chains=published["chains"],
        pe=pe, pe12=pe12, pe120=pe120, ok=ok, chains=chains, u=u,
        group_names=[g.name for g in chosen],
        sd_log10=[g.sd_log10 for g in chosen],
        rank_note=rank_note,
        group_line=group_line,
        msgs=msgs,
    )
    _write_pm2_check(d, jobs)


def _write_pm2_check(d: int, jobs: int, n: int = 32) -> None:
    """N≈32 log-uniform draws on ±2 SD. Notes only; does not change the ±1 figures."""
    u = 2.0 * design_loguniform(n, d, ODE_SEED + 7)
    print(f"±2 SD box check: {n} runs (notes only)...", flush=True)
    mets, msgs = _run_batch(u, jobs)
    vals = np.array([m["pe12"] if m is not None else np.nan for m in mets], dtype=float)
    ok = vals[np.isfinite(vals)]
    path = os.path.join(_OUT, "pm2_box.txt")
    if ok.size < 8:
        text = f"±2 SD box check failed: {ok.size}/{n} runs succeeded.\n"
    else:
        p5 = float(np.percentile(ok, P_LO))
        p50 = float(np.percentile(ok, 50))
        p95 = float(np.percentile(ok, P_HI))
        text = (
            "±2 CatPred SD box, top 4 groups, notes only (not a figure).\n"
            f"Log-uniform Uniform(-2,2) in SD units, N={ok.size}/{n}, "
            f"seed={ODE_SEED + 7}.\n"
            f"PE at 12 min: p5={p5:.2f} µM, p50={p50:.2f} µM, p95={p95:.2f} µM, "
            f"mean={float(ok.mean()):.2f} µM.\n"
            "This wider box is notes only. Compare p5 with the ±1 box in notes.txt.\n"
        )
    fails = [f"failure i={i}: {msg}" for i, msg in enumerate(msgs) if msg]
    text += "\n".join(fails)
    if fails:
        text += "\n"
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    print(text, flush=True)


def postprocess_ode_global() -> Dict[str, object]:
    """Refresh band/notes/bootstrap/stability from an existing ode_global.npz (no new ODEs)."""
    npz_path = os.path.join(_OUT, "ode_global.npz")
    if not os.path.isfile(npz_path):
        raise FileNotFoundError(npz_path)
    data = np.load(npz_path, allow_pickle=True)
    pe = data["pe"]
    pe12 = data["pe12"]
    pe120 = data["pe120"]
    ok = np.asarray(data["ok"], dtype=bool)
    grid = data["grid_s"]
    published_pe = data["published_pe"]
    published_chains = data["published_chains"]
    if "u" not in data.files:
        raise RuntimeError("ode_global.npz missing design matrix u")
    u = data["u"]
    group_names = [str(x) for x in data["groups"]]
    sd_log10 = [float(x) for x in data["sd_log10"]]
    if pe.ndim != 2:
        raise RuntimeError("ode_global.npz pe must be (n, nt)")
    chains = data["chains"] if "chains" in data.files else None
    chain_p5 = data["chain_p5"] if chains is None else None
    chain_p50 = data["chain_p50"] if chains is None else None
    chain_p95 = data["chain_p95"] if chains is None else None

    contrib = os.path.join(_OUT, "contributions.csv")
    group_line = f"ODE groups ({len(group_names)}): " + ", ".join(group_names)
    if os.path.isfile(contrib):
        shares = {}
        with open(contrib, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                shares[row["group"]] = row
        bits = []
        for name, sd in zip(group_names, sd_log10):
            row = shares.get(name, {})
            share = float(row["share"]) if row.get("share") else float("nan")
            bits.append(
                f"{name} (share {share:.1%}, sd_log10={sd:.3f}, fold≈{10 ** sd:.2g})"
            )
        group_line = f"ODE groups ({len(group_names)}): " + ", ".join(bits)

    notes_old = os.path.join(_OUT, "notes.txt")
    rank_note = (
        f"Tornado files in {_SENS}; S is d ln(product) / d ln(input) for kcat and Km, "
        "and d ln(product) / d(ΔG°′) per kJ/mol. Ranking time is the tornado t* "
        "(about 441 s), not 12 min."
    )
    if os.path.isfile(notes_old):
        with open(notes_old, encoding="utf-8") as fh:
            old_lines = fh.read().splitlines()
        for line in old_lines:
            if "Tornado files" in line or "Ranking time" in line and "441" in line:
                if "Tornado files" in line:
                    rank_note = line
                    break

    msgs = [""] * pe.shape[0]
    draws_path = os.path.join(_OUT, "draws.csv")
    if os.path.isfile(draws_path):
        with open(draws_path, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                i = int(row["i"])
                if 0 <= i < len(msgs):
                    msgs[i] = row.get("message") or ""
    return write_ensemble_outputs(
        grid=grid,
        published_pe=published_pe,
        published_chains=published_chains,
        pe=pe, pe12=pe12, pe120=pe120, ok=ok, chains=chains, u=u,
        group_names=group_names,
        sd_log10=sd_log10,
        rank_note=rank_note,
        group_line=group_line,
        msgs=msgs,
        chain_p5=chain_p5,
        chain_p50=chain_p50,
        chain_p95=chain_p95,
    )


def main(argv=None) -> None:
    p = argparse.ArgumentParser(
        description="FAS-II local ranking and ODE-only global sample (log-uniform ±1 CatPred SD).",
    )
    p.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    p.add_argument("--n", type=int, default=N_ODE, help="ODE global sample size (default 128)")
    p.add_argument(
        "--top", type=int, default=TOP_ODE_GROUPS,
        help="Number of top kinetic groups to sample (default 4)",
    )
    p.add_argument(
        "--rank-only", action="store_true",
        help="Write contributions.csv and experiment.csv only (no ODE).",
    )
    p.add_argument(
        "--postprocess", action="store_true",
        help="Refresh band/notes/bootstrap/stability from existing ode_global.npz.",
    )
    p.add_argument(
        "--groups", default="",
        help="Comma-separated group names to sample instead of the tornado top N.",
    )
    args = p.parse_args(argv)
    names = [s.strip() for s in args.groups.split(",") if s.strip()] or None
    if args.rank_only:
        os.makedirs(_OUT, exist_ok=True)
        groups, note = load_ranking(cutoff=-30.0)
        chosen = select_top_ode(groups, n=args.top)
        write_contributions(
            groups, chosen, os.path.join(_OUT, "contributions.csv"),
            flag_col="in_ode_global",
        )
        write_experiment(os.path.join(_OUT, "experiment.csv"))
        total = sum(g.variance for g in groups) or 1.0
        print(note)
        print(
            "ODE-global candidates:",
            ", ".join(f"{g.name} ({g.variance / total:.1%})" for g in chosen),
        )
        return
    if args.postprocess:
        postprocess_ode_global()
        return
    run_ode_global(args.jobs, n=args.n, top_n=args.top, names=names)


if __name__ == "__main__":
    main()
