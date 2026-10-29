#!/usr/bin/env python3
r"""
Monte Carlo propagation of parameter uncertainty through the fixed E. coli FAS-II network.

Each sample perturbs the PRE-RULE parameter values, then re-runs the rule layer, so
calibrations/laws (TesA logistic, hydrophobic Km, ΔG°′ irreversibility cutoff) act on the
sampled values exactly as in production. The network itself is fixed (CatPred is not re-run).

Sources (``--sources``):
    kcat       CatPred kcat, log10-normal with CatPred SD_total (FabI excluded: anchored)
    km         CatPred Km (forward + reverse/product), log10-normal with SD_total
    dg         equilibrator ΔG°′, multivariate normal with the joint covariance
    fabi_kcat  FabI kcat 15 ± 0.27 s⁻¹ (Rafi 2006 Table 2 fit error), normal
    fabh_ki    FabH feedback Ki, log10-normal SD 0.350218 (CatPred C16 SD_total)
    fabi_kcat_3x  FabI kcat log-uniform over 15/3..15×3 s⁻¹ (measured on C12-ACP only,
                  applied to every chain length)
    fabh_ki_10x   FabH feedback Ki log-uniform over ×0.1..×10 (Fig. 8 range)
    fabi_ki_10x   FabI feedback Ki log-uniform over ×0.1..×10 (no reported error)
fabi_kcat/fabi_kcat_3x and fabh_ki/fabh_ki_10x set the same parameter: use one of each pair.

Draw i is seeded from SeedSequence(seed, spawn_key=(i, source)), so a sample's draws do not
depend on N, worker count, or which other sources are on.

Every sample also records the free-acid chain-length distribution at 12 min and at the end
of the run (120 min), for the product-distribution band.

A Km that is one CatPred prediction (same enzyme, compound, value and SD) repeated across
reactions is a single parameter and gets a single draw (e.g. FabG-NADPH in all 14 FabG
reactions). Scheme "independent" draws every remaining parameter independently; "shared"
uses one factor per enzyme and parameter class, the fully correlated bound. CatPred SD_total
is an uncalibrated model uncertainty, so the band is a range under that uncertainty, not a
calibrated prediction interval.

Run (repo root, bees_env):
    python projects/fattyAcidSynthesis/fattyAcidSynthesis_ecoli/monte_carlo_uncertainty.py --zero-variance
    python projects/fattyAcidSynthesis/fattyAcidSynthesis_ecoli/monte_carlo_uncertainty.py --pilot
    python projects/fattyAcidSynthesis/fattyAcidSynthesis_ecoli/monte_carlo_uncertainty.py --n 500
    python projects/fattyAcidSynthesis/fattyAcidSynthesis_ecoli/monte_carlo_uncertainty.py --n 500 --sources kcat
    python projects/fattyAcidSynthesis/fattyAcidSynthesis_ecoli/monte_carlo_uncertainty.py --n 500 --scheme shared
    python projects/fattyAcidSynthesis/fattyAcidSynthesis_ecoli/monte_carlo_uncertainty.py --scenario fabi_ki_x2 --zero-variance
    python projects/fattyAcidSynthesis/fattyAcidSynthesis_ecoli/monte_carlo_uncertainty.py --suite review
"""

from __future__ import annotations

import argparse
import copy
import csv
import math
import multiprocessing as mp
import os
import signal
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

# One BLAS thread per worker: with one OpenBLAS pool per core in every forked worker the
# machine is oversubscribed and each draw runs ~10x slower. Must precede the numpy import.
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

OUT_DIR = os.path.join(_PROJECT_DIR, "output", "uncertainty")

# Append only: the index of a source seeds its draws.
SOURCES = ("kcat", "km", "dg", "fabi_kcat", "fabh_ki", "fabi_kcat_3x", "fabh_ki_10x",
           "fabi_ki_10x")
DEFAULT_SOURCES = ("kcat", "km", "dg", "fabi_kcat", "fabh_ki")
ESTIMATOR_SOURCES = ("kcat", "km", "dg")
WIDE_CORRECTED_SOURCES = ("fabi_kcat_3x", "fabh_ki_10x", "fabi_ki_10x")
EXCLUSIVE_SOURCES = (("fabi_kcat", "fabi_kcat_3x"), ("fabh_ki", "fabh_ki_10x"))
SCENARIOS = ("none", "fabi_ki_x0.5", "fabi_ki_x2", "tesa_off")

FABI_EC = "1.3.1.9"
FABI_RULE = "fabI_enoyl_reductase_measured_kcat"
TESA_RULE = "tesa_long_chain_preference"
FABI_KCAT_MEAN = 15.0     # s⁻¹; Rafi 2006
FABI_KCAT_SD = 16.0 / 60.0
FABH_KI_LOG10_SD = 0.350218  # CatPred SD_total, hexadecanoate vs FabH (SI_CatPred_FabH_Ki_table.csv)
WIDE_FABI_KCAT_FOLD = 3.0
WIDE_KI_FOLD = 10.0

# Fig. 6 bins: (chain carbons, double bonds).
CHAIN_KEYS = ((4, 0), (6, 0), (8, 0), (10, 0), (12, 0), (12, 1), (14, 0), (14, 1),
              (16, 0), (16, 1), (18, 0), (18, 1), (20, 0), (20, 1))
CHAIN_TIMES_MIN = (12.0, 120.0)

T_END_S = 7200.0          # same horizon as the production S2A script (exact baseline match)
GRID_DT_S = 6.0
T_STAR_S = 375.0          # sensitivity-analysis t* (first-order cross-check)
NEAR_CUTOFF_NSIGMA = 2.0
PERCENTILES = (5, 25, 50, 75, 95)


@dataclass
class Dim:
    source: str
    rxn: int               # index into core+edge reaction list (-1: global)
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

    def size(self, source: str) -> int:
        return {
            "kcat": len(self.kcat),
            "km": len(self.km),
            "dg": self.dg_chol.shape[1] if self.dg else 0,
            "fabi_kcat": 1,
            "fabh_ki": 1,
            "fabi_kcat_3x": 1,
            "fabh_ki_10x": 1,
            "fabi_ki_10x": 1,
        }[source]


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
    compound with the identical value and SD is one physical parameter: it gets one draw.
    """
    merged: Dict[tuple, Dim] = {}
    out: List[Dim] = []
    for d in dims:
        if d.key is None:
            out.append(d)
            continue
        # SD went through a log10 -> linear -> log10 round trip: compare at 9 significant digits.
        k = (d.enzyme.lower().strip(), d.key.lower().strip(), f"{d.mean:.9g}", f"{d.sd:.9g}")
        first = merged.get(k)
        if first is None:
            merged[k] = Dim(d.source, d.rxn, d.key, d.mean, d.sd, d.enzyme, d.label, d.role,
                            rxns=d.rxns, targets=d.km_targets())
            out.append(merged[k])
            continue
        first.targets = first.targets + d.km_targets()
        roles = set(first.role.split("/")) | {d.role}
        first.role = "/".join(sorted(r for r in roles if r))
    for d in out:
        if len(d.targets) > 1:
            d.label = f"{d.enzyme} [{d.key}] ({len(d.targets)} reactions)"
    return out


def _sample_rng(seed: int, i: int, source: str) -> np.random.Generator:
    return np.random.default_rng(
        np.random.SeedSequence(seed, spawn_key=(i, SOURCES.index(source)))
    )


def _group_keys(dims: Sequence[Dim]) -> List[str]:
    return [d.enzyme.lower().strip() for d in dims]


def draw_sample(inv: Inventory, seed: int, i: int, sources: Sequence[str],
                scheme: str = "independent") -> Dict[str, np.ndarray]:
    """Linear parameter values for sample ``i``. Inactive sources stay at their means."""
    out: Dict[str, np.ndarray] = {}
    for src, dims in (("kcat", inv.kcat), ("km", inv.km)):
        mean = np.array([d.mean for d in dims], dtype=float)
        if src in sources and dims:
            sd = np.array([d.sd for d in dims], dtype=float)
            rng = _sample_rng(seed, i, src)
            if scheme == "shared":
                keys = _group_keys(dims)
                uniq = sorted(set(keys))
                zg = dict(zip(uniq, rng.standard_normal(len(uniq))))
                z = np.array([zg[k] for k in keys])
            else:
                z = rng.standard_normal(len(dims))
            out[src] = mean * np.power(10.0, sd * z)
        else:
            out[src] = mean
    dg_mean = np.array([d.mean for d in inv.dg], dtype=float)
    if "dg" in sources and inv.dg:
        z = _sample_rng(seed, i, "dg").standard_normal(inv.dg_chol.shape[1])
        out["dg"] = dg_mean + inv.dg_chol @ z
    else:
        out["dg"] = dg_mean
    if "fabi_kcat" in sources:
        z = _sample_rng(seed, i, "fabi_kcat").standard_normal()
        out["fabi_kcat"] = np.array(max(1e-6, FABI_KCAT_MEAN + FABI_KCAT_SD * z))
    elif "fabi_kcat_3x" in sources:
        out["fabi_kcat"] = np.array(
            FABI_KCAT_MEAN * _log_uniform(seed, i, "fabi_kcat_3x", WIDE_FABI_KCAT_FOLD))
    else:
        out["fabi_kcat"] = np.array(FABI_KCAT_MEAN)
    if "fabh_ki" in sources:
        z = _sample_rng(seed, i, "fabh_ki").standard_normal()
        out["fabh_ki_factor"] = np.array(10.0 ** (FABH_KI_LOG10_SD * z))
    elif "fabh_ki_10x" in sources:
        out["fabh_ki_factor"] = np.array(_log_uniform(seed, i, "fabh_ki_10x", WIDE_KI_FOLD))
    else:
        out["fabh_ki_factor"] = np.array(1.0)
    if "fabi_ki_10x" in sources:
        out["fabi_ki_factor"] = np.array(_log_uniform(seed, i, "fabi_ki_10x", WIDE_KI_FOLD))
    else:
        out["fabi_ki_factor"] = np.array(1.0)
    return out


def _log_uniform(seed: int, i: int, source: str, fold: float) -> float:
    """Factor drawn log-uniformly from [1/fold, fold]."""
    u = _sample_rng(seed, i, source).uniform(-1.0, 1.0)
    return float(fold ** u)


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


def _metrics(sim, grid_t):
    pe = sa._pe_series(sim)
    t_min = sim.t / 60.0
    bees_at = np.interp(sa.S2A_EXP_T_MIN, t_min, pe)
    rmse = float(np.sqrt(np.mean((bees_at - sa.S2A_EXP_UM) ** 2)))
    return {
        "pe_grid": np.interp(grid_t, sim.t, pe),
        "pe12": float(np.interp(12.0, t_min, pe)),
        "pe_tstar": float(np.interp(T_STAR_S, sim.t, pe)),
        "rmse": rmse,
        "bees_at_exp": bees_at,
        "chains": _chain_profile(sim),
    }


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


# Fork workers inherit these from the parent.
_W: Dict[str, object] = {}

STATUS_CONVERGED = "converged"
STATUS_TIMEOUT = "timeout"
STATUS_FAILED = "solver_failure"


class SampleTimeout(Exception):
    pass


def _on_alarm(signum, frame):
    raise SampleTimeout()


def _run_sample(i: int):
    """(i, metrics or None, irreversible flags, status, message, wall seconds)."""
    enl = _W["enlarger"]
    values = {k: v[i] for k, v in _W["values"].items()}
    timeout = float(_W.get("timeout") or 0.0)
    irrev = None
    t0 = time.time()
    if timeout > 0:
        signal.signal(signal.SIGALRM, _on_alarm)
        signal.setitimer(signal.ITIMER_REAL, timeout)
    try:
        model = apply_sample(enl, _W["model"], _W["snaps"], _W["inv"], values,
                             _W["T_K"], _W["scenario"])
        irrev = np.array([bool(getattr(r.thermo, "irreversible", True))
                          if r.thermo is not None else True
                          for r in model.core_reactions])
        sim = sa._simulate(enl, model, end_time=T_END_S)
        if not sim.success or sim.t.size == 0 or sim.t[-1] < T_END_S - 1e-6:
            return (i, None, irrev, STATUS_FAILED, f"solver: {sim.message}",
                    time.time() - t0)
        return i, _metrics(sim, _W["grid_t"]), irrev, STATUS_CONVERGED, "", time.time() - t0
    except SampleTimeout:
        return (i, None, irrev, STATUS_TIMEOUT, f"exceeded {timeout:g} s wall-clock",
                time.time() - t0)
    except (ValueError, FloatingPointError, np.linalg.LinAlgError, ArithmeticError) as exc:
        return (i, None, irrev, STATUS_FAILED, f"{type(exc).__name__}: {exc}",
                time.time() - t0)
    finally:
        if timeout > 0:
            signal.setitimer(signal.ITIMER_REAL, 0)


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


def _tag(args) -> str:
    srcs = "+".join(args.sources) if args.sources else "none"
    return f"{args.scheme}_{srcs}_{args.scenario}_n{args.n}_s{args.seed}"


def write_inventory(inv: Inventory, model, out_dir: str) -> None:
    path = os.path.join(out_dir, "parameter_inventory.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["source", "reaction", "enzyme", "key", "role", "pre_rule_mean",
                    "sd", "sd_units", "status", "n_reactions"])
        for d in inv.kcat:
            w.writerow(["kcat", d.label, d.enzyme, "", "", f"{d.mean:.6g}", f"{d.sd:.6g}",
                        "log10", "sampled", len(d.rxns or (d.rxn,))])
        for d in inv.km:
            w.writerow(["km", d.label, d.enzyme, d.key or "", d.role, f"{d.mean:.6g}",
                        f"{d.sd:.6g}", "log10", "sampled", len(d.km_targets())])
        for d in inv.dg:
            w.writerow(["dg", d.label, d.enzyme, "", "", f"{d.mean:.6g}", f"{d.sd:.6g}",
                        "kJ/mol", "sampled (joint covariance)"])
        w.writerow(["fabi_kcat", "FabI (all reactions)", "FabI", "", "", FABI_KCAT_MEAN,
                    f"{FABI_KCAT_SD:.6g}", "1/s", "sampled (normal)"])
        w.writerow(["fabh_ki", "FabH feedback (all inhibitors)", "FabH", "", "",
                    f"{inv.fabh_ki_mM:.6g}", FABH_KI_LOG10_SD, "log10", "sampled"])
        w.writerow(["fabi_kcat_3x", "FabI (all reactions)", "FabI", "", "", FABI_KCAT_MEAN,
                    f"x1/{WIDE_FABI_KCAT_FOLD:g}..x{WIDE_FABI_KCAT_FOLD:g}", "fold",
                    "sampled (log-uniform)"])
        w.writerow(["fabh_ki_10x", "FabH feedback (all inhibitors)", "FabH", "", "",
                    f"{inv.fabh_ki_mM:.6g}", f"x1/{WIDE_KI_FOLD:g}..x{WIDE_KI_FOLD:g}", "fold",
                    "sampled (log-uniform)"])
        w.writerow(["fabi_ki_10x", "FabI feedback (all inhibitors)", "FabI", "", "", "",
                    f"x1/{WIDE_KI_FOLD:g}..x{WIDE_KI_FOLD:g}", "fold",
                    "sampled (log-uniform)"])
        for src, lab, why in inv.excluded:
            w.writerow([src, lab, "", "", "", "", "", "", f"fixed: {why}"])
    print(f"wrote {path}")
    if inv.dg:
        path = os.path.join(out_dir, "dg_covariance.npz")
        np.savez_compressed(path, labels=np.array([d.label for d in inv.dg]),
                            mean=np.array([d.mean for d in inv.dg]),
                            cov=inv.dg_chol @ inv.dg_chol.T)
        print(f"wrote {path}")

    rxns = list(model.core_reactions)
    path = os.path.join(out_dir, "near_cutoff_reactions.csv")
    cutoff = -float(RULES.by_name("dgr_irreversibility").params["dgr_kjmol_cutoff"])
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["reaction", "enzyme", "dG_kJmol", "sigma_kJmol", "z_from_cutoff",
                    "baseline_irreversible"])
        for d in inv.dg:
            if d.sd > 0 and abs(d.mean - cutoff) <= NEAR_CUTOFF_NSIGMA * d.sd:
                th = rxns[d.rxn].thermo
                w.writerow([d.label, d.enzyme, f"{d.mean:.3f}", f"{d.sd:.3f}",
                            f"{(d.mean - cutoff) / d.sd:+.2f}", bool(th.irreversible)])
    print(f"wrote {path}")


def setup():
    """One network build shared by every configuration run in this process."""
    os.makedirs(OUT_DIR, exist_ok=True)
    enlarger, model, snaps, inv, T_K = prepare()
    print(
        f"Inventory: kcat={len(inv.kcat)} km={len(inv.km)} "
        f"(covering {sum(len(d.km_targets()) for d in inv.km)} reaction Km entries) "
        f"dg={len(inv.dg)} "
        f"(rank {inv.dg_chol.shape[1] if inv.dg else 0}) fabi_kcat=1 fabh_ki=1; "
        f"fixed={len(inv.excluded)}; FabH Ki={inv.fabh_ki_mM * 1000:.3g} µM",
        flush=True,
    )
    write_inventory(inv, model, OUT_DIR)
    grid_t = np.arange(0.0, T_END_S + 0.5 * GRID_DT_S, GRID_DT_S)
    print("Baseline (production rules, no resampling)...", flush=True)
    base_sim = sa._simulate(enlarger, model, end_time=T_END_S)
    base = _metrics(base_sim, grid_t)
    base_irrev = np.array([bool(r.thermo.irreversible) if r.thermo is not None else True
                           for r in model.core_reactions])
    print(f"  baseline RMSE={base['rmse']:.3f} µM  PE@12={base['pe12']:.2f} µM", flush=True)
    return dict(enlarger=enlarger, model=model, snaps=snaps, inv=inv, T_K=T_K,
                grid_t=grid_t, base=base, base_irrev=base_irrev)


def run(args, ctx):
    enlarger, model, snaps, inv, T_K = (ctx[k] for k in
                                        ("enlarger", "model", "snaps", "inv", "T_K"))
    grid_t, base, base_irrev = ctx["grid_t"], ctx["base"], ctx["base_irrev"]
    n = args.n
    draws = [draw_sample(inv, args.seed, i, args.sources, args.scheme) for i in range(n)]
    values = {k: np.stack([d[k] for d in draws]) for k in draws[0]}

    _W.update(enlarger=enlarger, model=model, snaps=snaps, inv=inv, values=values,
              T_K=T_K, scenario=args.scenario, grid_t=grid_t, timeout=args.timeout)

    if not args.sources and args.scenario == "none":
        zm = apply_sample(enlarger, model, snaps, inv, {k: v[0] for k, v in values.items()},
                          T_K, args.scenario)
        _check_identical(model, zm)

    jobs = max(1, min(args.jobs, n))
    print(f"Running {n} samples on {jobs} workers (tag {_tag(args)})...", flush=True)
    t0 = time.time()
    results = [None] * n
    if jobs == 1:
        for i in range(n):
            results[i] = _run_sample(i)
    else:
        ctx = mp.get_context("fork")
        with ProcessPoolExecutor(max_workers=jobs, mp_context=ctx) as ex:
            for k, res in enumerate(ex.map(_run_sample, range(n), chunksize=1)):
                results[res[0]] = res
                if (k + 1) % max(1, n // 10) == 0:
                    done = [r for r in results if r is not None]
                    n_to = sum(r[3] == STATUS_TIMEOUT for r in done)
                    n_fail = sum(r[3] == STATUS_FAILED for r in done)
                    print(f"  {k + 1}/{n} done ({time.time() - t0:.0f} s; "
                          f"timeouts {n_to}, solver failures {n_fail})", flush=True)
    print(f"Samples finished in {time.time() - t0:.0f} s", flush=True)
    save_outputs(args, inv, model, values, results, base, base_irrev, grid_t)


def _check_identical(base_model, zm) -> None:
    """Zero-variance copy must match the production model parameter-for-parameter."""
    diffs = 0
    for a, b in zip(base_model.core_reactions, zm.core_reactions):
        ka, kb = a.kinetics, b.kinetics
        if ka is not None and (ka.kcat != kb.kcat or ka.km != kb.km
                               or (ka.km_per_substrate or {}) != (kb.km_per_substrate or {})):
            diffs += 1
        ta, tb = a.thermo, b.thermo
        if ta is not None and (ta.irreversible != tb.irreversible
                               or not (ta.keq == tb.keq or (math.isnan(ta.keq) and math.isnan(tb.keq)))):
            diffs += 1
        if a.feedback_inhibitors != b.feedback_inhibitors:
            diffs += 1
        if getattr(a.template, "reversible", None) != getattr(b.template, "reversible", None):
            diffs += 1
    print(f"Zero-variance parameter check: {diffs} differing reactions", flush=True)
    if diffs:
        raise RuntimeError("zero-variance model differs from production model")


def _bootstrap_ci(x, q, n_boot=2000, seed=0):
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(x), size=(n_boot, len(x)))
    qs = np.percentile(x[idx], q, axis=1)
    return np.percentile(qs, [2.5, 97.5])


def save_outputs(args, inv, model, values, results, base, base_irrev, grid_t):
    tag = _tag(args)
    n = len(results)
    ok = np.array([r[1] is not None for r in results])
    pe = np.full((n, grid_t.size), np.nan)
    pe12 = np.full(n, np.nan)
    pe_tstar = np.full(n, np.nan)
    rmse = np.full(n, np.nan)
    at_exp = np.full((n, sa.S2A_EXP_T_MIN.size), np.nan)
    chains = np.full((n, len(CHAIN_TIMES_MIN), len(CHAIN_KEYS)), np.nan)
    irrev = np.zeros((n, len(base_irrev)), dtype=bool)
    msgs = [""] * n
    status = np.array([r[3] for r in results])
    wall_s = np.array([r[5] for r in results], dtype=float)
    for i, m, ir, _st, msg, _dt in results:
        if ir is not None:
            irrev[i] = ir
        msgs[i] = msg
        if m is None:
            continue
        pe[i] = m["pe_grid"]
        pe12[i] = m["pe12"]
        pe_tstar[i] = m["pe_tstar"]
        rmse[i] = m["rmse"]
        at_exp[i] = m["bees_at_exp"]
        chains[i] = m["chains"]
    flips = (irrev != base_irrev[None, :]).sum(axis=1)

    npz = os.path.join(OUT_DIR, f"mc_{tag}.npz")
    np.savez_compressed(
        npz, t_s=grid_t, pe=pe, ok=ok, status=status, wall_s=wall_s,
        timeout_s=float(args.timeout), pe12=pe12, pe_tstar=pe_tstar, rmse=rmse,
        pe_at_exp=at_exp, irreversible=irrev, baseline_irreversible=base_irrev,
        baseline_pe=base["pe_grid"], baseline_pe12=base["pe12"], baseline_rmse=base["rmse"],
        exp_t_min=sa.S2A_EXP_T_MIN, exp_um=sa.S2A_EXP_UM,
        chains=chains, baseline_chains=base["chains"],
        chain_keys=np.array(CHAIN_KEYS), chain_times_min=np.array(CHAIN_TIMES_MIN),
        **{f"draw_{k}": v for k, v in values.items()},
    )
    print(f"wrote {npz}")

    good = pe[ok]
    path = os.path.join(OUT_DIR, f"mc_{tag}_percentiles.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["t_min"] + [f"p{q}" for q in PERCENTILES] + ["mean", "baseline"])
        if good.size:
            pct = np.percentile(good, PERCENTILES, axis=0)
            mean = good.mean(axis=0)
        else:
            pct = np.full((len(PERCENTILES), grid_t.size), np.nan)
            mean = np.full(grid_t.size, np.nan)
        for k, t in enumerate(grid_t):
            w.writerow([f"{t / 60:.4f}"] + [f"{pct[q, k]:.5g}" for q in range(len(PERCENTILES))]
                       + [f"{mean[k]:.5g}", f"{base['pe_grid'][k]:.5g}"])
    print(f"wrote {path}")

    path = os.path.join(OUT_DIR, f"mc_{tag}_draws.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["i", "ok", "status", "wall_s", "pe12_uM", "rmse_uM",
                    "n_reversibility_flips", "fabi_kcat",
                    "fabh_ki_factor", "fabi_ki_factor", "median_log10_kcat_ratio",
                    "median_log10_km_ratio", "message"])
        km_mean = np.array([d.mean for d in inv.km]) if inv.km else None
        kc_mean = np.array([d.mean for d in inv.kcat]) if inv.kcat else None
        for i in range(n):
            mk = (np.median(np.log10(values["kcat"][i] / kc_mean)) if kc_mean is not None else 0.0)
            mm = (np.median(np.log10(values["km"][i] / km_mean)) if km_mean is not None else 0.0)
            w.writerow([i, int(ok[i]), status[i], f"{wall_s[i]:.1f}", f"{pe12[i]:.4f}",
                        f"{rmse[i]:.4f}", int(flips[i]),
                        f"{float(values['fabi_kcat'][i]):.5f}",
                        f"{float(values['fabh_ki_factor'][i]):.5f}",
                        f"{float(values['fabi_ki_factor'][i]):.5f}", f"{mk:+.4f}",
                        f"{mm:+.4f}", msgs[i]])
    print(f"wrote {path}")

    lines = [
        f"FAS-II Monte Carlo: {tag}",
        f"sources={list(args.sources)} scheme={args.scheme} scenario={args.scenario} "
        f"seed={args.seed} N={n}",
        f"sampled dims: kcat={len(inv.kcat)} km={len(inv.km)} dg={len(inv.dg)} "
        f"fabi_kcat=1 fabh_ki=1 (inactive sources held at mean)",
        f"convergence: {int(ok.sum())}/{n} draws converged ({100.0 * ok.mean():.1f}%); "
        f"timeouts {int((status == STATUS_TIMEOUT).sum())} "
        f"(limit {args.timeout:g} s per draw), "
        f"solver failures {int((status == STATUS_FAILED).sum())}",
        "wall time per draw (s): " + "  ".join(
            f"p{q}={v:.1f}" for q, v in zip((50, 90, 99), np.percentile(wall_s, [50, 90, 99]))
        ) + f"  max={wall_s.max():.1f}",
        "all percentiles below use converged draws only",
        f"baseline: RMSE={base['rmse']:.3f} µM  PE@12={base['pe12']:.2f} µM",
    ]
    if ok.any():
        p = np.percentile(pe12[ok], PERCENTILES)
        lines += [
            "PE@12 min (µM): " + "  ".join(f"p{q}={v:.2f}" for q, v in zip(PERCENTILES, p)),
            f"PE@12 5-95% band width = {p[4] - p[0]:.2f} µM",
            f"PE@12 median - baseline = {p[2] - base['pe12']:+.2f} µM",
            f"PE@12 ratio p95/p5 = {p[4] / p[0]:.2f}" if p[0] > 0 else "PE@12 p5 <= 0",
            "RMSE (µM): " + "  ".join(
                f"p{q}={v:.2f}" for q, v in zip(PERCENTILES, np.percentile(rmse[ok], PERCENTILES))
            ),
            f"fraction of samples with RMSE <= baseline RMSE: "
            f"{np.mean(rmse[ok] <= base['rmse'] + 1e-12):.3f}",
        ]
        lo = np.percentile(at_exp[ok], 5, axis=0)
        hi = np.percentile(at_exp[ok], 95, axis=0)
        inside = (sa.S2A_EXP_UM >= lo) & (sa.S2A_EXP_UM <= hi)
        lines.append(f"experimental points inside 5-95% band: {int(inside.sum())}/{inside.size}")
        lines.append("  t_min   exp    p5     p50    p95   inside")
        med = np.percentile(at_exp[ok], 50, axis=0)
        for t, e, a, b, c, ins in zip(sa.S2A_EXP_T_MIN, sa.S2A_EXP_UM, lo, med, hi, inside):
            lines.append(f"  {t:5.2f} {e:6.2f} {a:6.2f} {b:6.2f} {c:6.2f}   {int(ins)}")
        lines.append(f"var ln PE(t*={T_STAR_S:g} s) = {np.var(np.log(pe_tstar[ok])):.4f}")
        lines.append(
            f"reversibility flips per sample: mean={flips[ok].mean():.2f} "
            f"max={flips.max()} samples with >=1 flip={int((flips > 0).sum())}"
        )
        flipped = np.where((irrev != base_irrev[None, :]).any(axis=0))[0]
        seen: Dict[str, int] = {}
        labels = [sa._reaction_label(r, seen) for r in model.core_reactions]
        for j in flipped:
            frac = float((irrev[:, j] != base_irrev[j]).mean())
            lines.append(f"  flip {frac:5.1%}  {labels[j]}")
        if ok.sum() >= 100:
            lines.append("convergence (PE@12 percentiles, 95% bootstrap CI):")
            idx_ok = np.where(ok)[0]
            for m in (100, 250, 500):
                sel = pe12[idx_ok[idx_ok < m]] if m <= n else None
                if sel is None or sel.size < 50:
                    continue
                p5, p95 = np.percentile(sel, [5, 95])
                c5 = _bootstrap_ci(sel, 5)
                c95 = _bootstrap_ci(sel, 95)
                lines.append(
                    f"  first {m:3d} draws (ok={sel.size}): p5={p5:.2f} [{c5[0]:.2f},{c5[1]:.2f}]"
                    f"  p95={p95:.2f} [{c95[0]:.2f},{c95[1]:.2f}]"
                )
    if (~ok).any():
        lines.append("failures:")
        for i in np.where(~ok)[0][:20]:
            lines.append(f"  i={i}: {msgs[i]}")
    text = "\n".join(lines) + "\n"
    path = os.path.join(OUT_DIR, f"mc_{tag}_stats.txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    print(text)
    print(f"wrote {path}")


def main(argv=None):
    p = argparse.ArgumentParser(description="FAS-II Monte Carlo uncertainty band.")
    p.add_argument("--n", type=int, default=500)
    p.add_argument("--seed", type=int, default=20261001)
    p.add_argument("--scheme", choices=("independent", "shared"), default="independent",
                   help="shared: one CatPred factor per enzyme and parameter class (SI only)")
    p.add_argument("--sources", nargs="*", default=list(DEFAULT_SOURCES), choices=SOURCES)
    p.add_argument("--scenario", choices=SCENARIOS, default="none")
    p.add_argument("--pilot", action="store_true", help="N=30")
    p.add_argument("--zero-variance", action="store_true",
                   help="N=1, no sources; must reproduce the production baseline")
    p.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    p.add_argument("--timeout", type=float, default=600.0,
                   help="wall-clock seconds per draw before it is recorded as a timeout "
                   "(0 = no limit). With 27 workers, sampled draws take ~90 s median and "
                   "~200 s at p99; a tighter limit drops valid stiff draws and biases the band")
    p.add_argument("--suite", choices=("pilot", "full", "review"), default=None,
                   help="pilot: zero-variance + N=30 CatPred-only + N=30 all sources. "
                   "full: zero-variance, main run, single-source runs, shared scheme, "
                   "and the three fixed scenarios (one network build). "
                   "review: zero-variance, estimators only (kcat, km, dg; corrected "
                   "parameters at nominal) independent and enzyme-correlated, estimators "
                   "+ wide corrected ranges (stress test), and each estimator source alone.")
    args = p.parse_args(argv)
    if args.pilot:
        args.n = 30
    if args.zero_variance:
        args.n = 1
        args.sources = []
    args.sources = [s for s in SOURCES if s in args.sources]
    for a, b in EXCLUSIVE_SOURCES:
        if a in args.sources and b in args.sources:
            p.error(f"--sources: {a} and {b} set the same parameter; pick one")

    def cfg(**kw):
        c = argparse.Namespace(**vars(args))
        for k, v in kw.items():
            setattr(c, k, v)
        return c

    zero = cfg(n=1, sources=[], scheme="independent", scenario="none")
    if args.suite == "pilot":
        configs = [zero,
                   cfg(n=30, sources=["kcat", "km"], scheme="independent", scenario="none"),
                   cfg(n=30, sources=list(DEFAULT_SOURCES), scheme="independent",
                       scenario="none")]
    elif args.suite == "full":
        configs = [zero, cfg(sources=list(DEFAULT_SOURCES), scheme="independent",
                             scenario="none")]
        configs += [cfg(sources=[s], scheme="independent", scenario="none")
                    for s in DEFAULT_SOURCES]
        configs.append(cfg(sources=list(DEFAULT_SOURCES), scheme="shared", scenario="none"))
        configs += [cfg(n=1, sources=[], scheme="independent", scenario=s)
                    for s in SCENARIOS if s != "none"]
    elif args.suite == "review":
        configs = [zero,
                   cfg(sources=list(ESTIMATOR_SOURCES), scheme="independent", scenario="none"),
                   cfg(sources=list(ESTIMATOR_SOURCES), scheme="shared", scenario="none"),
                   cfg(sources=list(ESTIMATOR_SOURCES + WIDE_CORRECTED_SOURCES),
                       scheme="independent", scenario="none")]
        configs += [cfg(sources=[s], scheme="independent", scenario="none")
                    for s in ESTIMATOR_SOURCES]
    else:
        configs = [args]
    ctx = setup()
    for c in configs:
        run(c, ctx)


if __name__ == "__main__":
    main()
