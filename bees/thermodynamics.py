"""ΔG°′ / Keq via equilibrator-api. Reverse kcat is filled in by the rules layer."""

from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import pickle
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from bees.common import canonical_smiles, smiles_to_inchi, smiles_to_inchikey,R

class CompoundSubstitutor:
    """Replace unresolved SMILES before equilibrator lookup. None = pass-through."""

    def substitute(self, label: str, smiles: Optional[str]) -> Optional[str]:
        """Return a replacement SMILES, or None to leave it unchanged."""
        return None

_DEFAULT_CACHE_DIR = Path(os.environ.get(
    "BEES_THERMO_CACHE",
    str(Path.home() / ".cache" / "bees" / "equilibrator"),
))

_DGR_IRREVERSIBLE_KJMOL = 30.0  # ΔG°' < −this ⇒ forward-irreversible (rule layer)
_SIGMA_WARN_KJMOL = 10.0        # σ above this ⇒ Keq uncertain by >50×

logger = logging.getLogger("BEES")

def _is_disabled() -> bool:
    return os.environ.get("BEES_DISABLE_THERMO", "").strip() not in ("", "0", "false", "False")

@dataclass
class ThermoData:
    """Per-reaction ΔG°′ / Keq / kcat_rev. Attached as reaction.thermo (runtime-only)."""
    dgr_prime_kJmol: float
    sigma_kJmol: float
    keq: float
    kcat_rev: Optional[float]
    irreversible: bool
    source: str

def _disabled_thermo() -> ThermoData:
    return ThermoData(
        dgr_prime_kJmol=float("nan"),
        sigma_kJmol=float("nan"),
        keq=float("nan"),
        kcat_rev=None,
        irreversible=True,
        source="disabled",
    )

def _fallback_thermo(reason: str) -> ThermoData:
    logger.info("Thermo fallback: %s", reason)
    return ThermoData(
        dgr_prime_kJmol=float("nan"),
        sigma_kJmol=float("nan"),
        keq=float("nan"),
        kcat_rev=None,
        irreversible=True,
        source="fallback",
    )

class _DiskCache:
    """Pickle-on-disk cache for equilibrator lookups. Caches compound hits and misses."""

    def __init__(self, root: Path = _DEFAULT_CACHE_DIR):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)

    def _path(self, namespace: str, key: str) -> Path:
        h = hashlib.sha256(key.encode("utf-8")).hexdigest()
        return self.root / namespace / f"{h}.pkl"

    def get(self, namespace: str, key: str, default=...):
        p = self._path(namespace, key)
        if not p.exists():
            return default
        try:
            with p.open("rb") as f:
                return pickle.load(f)
        except Exception:
            return default

    def set(self, namespace: str, key: str, value) -> None:
        p = self._path(namespace, key)
        p.parent.mkdir(parents=True, exist_ok=True)
        try:
            with p.open("wb") as f:
                pickle.dump(value, f)
        except Exception as e:
            logger.warning("Thermo cache write failed (%s): %s", p, e)

class ThermoEngine:
    """Compute ΔG°′ / Keq for reactions via equilibrator-api. One instance per run."""

    def __init__(
        self,
        pH: float = 7.0,
        ionic_strength_M: float = 0.25,
        pMg: float = 3.0,
        T_K: float = 310.15,
        irreversible_cutoff_kJmol: float = _DGR_IRREVERSIBLE_KJMOL,
        cache: Optional[_DiskCache] = None,
        substitutor: Optional["CompoundSubstitutor"] = None,
    ):
        self.pH = float(pH)
        self.ionic_strength_M = float(ionic_strength_M)
        self.pMg = float(pMg)
        self.T_K = float(T_K)
        self.irreversible_cutoff_kJmol = float(irreversible_cutoff_kJmol)
        self.cache = cache or _DiskCache()
        self._cc = None
        self._cc_load_failed = False
        # None = pass-through; enlarger injects a CompoundSubstitutor at init.
        self.substitutor: Optional[CompoundSubstitutor] = substitutor

    def _load_cc(self):
        """Load ComponentContribution and set environmental conditions."""
        if self._cc is not None or self._cc_load_failed:
            return self._cc
        try:
            from equilibrator_api import ComponentContribution, Q_  # type: ignore
        except ImportError as e:
            logger.warning(
                "equilibrator-api not installed (%s). "
                "Run: conda install equilibrator-api  "
                "or set BEES_DISABLE_THERMO=1 to suppress this warning.",
                e,
            )
            self._cc_load_failed = True
            return None
        try:
            logger.info(
                "Loading ComponentContribution (first run downloads ~1.3 GB from Zenodo; "
                "subsequent loads are fast)..."
            )
            cc = ComponentContribution()
            cc.p_h = Q_(self.pH)
            cc.p_mg = Q_(self.pMg)
            cc.ionic_strength = Q_(f"{self.ionic_strength_M} M")
            cc.temperature = Q_(f"{self.T_K} K")
            self._cc = cc
            logger.info("ComponentContribution ready (pH=%.1f, I=%.3f M, T=%.1f K).",
                        self.pH, self.ionic_strength_M, self.T_K)
        except Exception as e:
            logger.warning("ComponentContribution init failed (%s); thermo disabled.", e)
            self._cc_load_failed = True
            return None
        return self._cc

    def _resolve_compound(self, smiles: Optional[str]):
        """SMILES → equilibrator Compound. Returns None on miss."""
        cc = self._load_cc()
        if cc is None:
            return None

        canon = canonical_smiles(smiles)
        if not canon:
            return None

        inchi = smiles_to_inchi(canon)
        ikey = smiles_to_inchikey(canon)
        cache_key = inchi or ikey or canon
        if not cache_key:
            return None

        # Do not pickle Compound objects (SQLAlchemy DetachedInstanceError); cache misses only.
        cached = self.cache.get("compound", cache_key, default=...)
        if cached is None:
            return None  # cached miss

        cpd = None
        try:
            if inchi:
                cpd = cc.get_compound_by_inchi(inchi)

            # InChIKey 14-char prefix fallback (connectivity; ignores stereo/protonation).
            if cpd is None and ikey:
                ikey_prefix = ikey.split("-")[0]  # first 14-char block
                matches = cc.search_compound_by_inchi_key(ikey_prefix)
                if matches:
                    if len(matches) == 1:
                        cpd = matches[0]
                    else:
                        for m in matches:
                            if getattr(m, "inchi_key", None) == ikey:
                                cpd = m
                                break
                        if cpd is None:
                            cpd = matches[0]
                            logger.debug(
                                "Ambiguous InChIKey prefix %s → %d matches; "
                                "using first (%s).",
                                ikey_prefix, len(matches),
                                getattr(cpd, "inchi_key", "?"),
                            )

            if cpd is None:
                try:
                    cpd = cc.get_compound(canon)
                except Exception:
                    cpd = None
        except Exception as e:
            logger.debug("Compound lookup failed for smiles=%r: %s", smiles, e)
            cpd = None

        if cpd is None:
            self.cache.set("compound", cache_key, None)  # persist miss only
            logger.info(
                "Compound not found in eQuilibrator database: smiles=%r (inchi=%r)",
                smiles, inchi,
            )
        return cpd

    def _build_cc_reaction(
        self,
        stoichiometry: Dict[str, int],
        smiles_map: Dict[str, str],
    ):
        """Build the equilibrator Reaction for a stoichiometry.

        Returns ``("ok", Reaction)``, ``("null", None)`` when substitution cancels the
        reaction (ΔG°′ = 0 exactly), or ``(None, None)`` when it cannot be built.
        H+ and H2O must stay in the stoichiometry (needed for is_balanced(); equilibrator accounts for both).
        """
        try:
            from equilibrator_api import Reaction  # type: ignore
        except ImportError:
            return None, None

        effective_smi_map: Dict[str, Optional[str]] = {}
        for label in stoichiometry:
            smi = smiles_map.get(label)
            sub_smi = self.substitutor.substitute(label, smi) if self.substitutor else None
            effective_smi = sub_smi if sub_smi is not None else smi
            if sub_smi is not None and sub_smi != smi:
                logger.debug(
                    "CompoundSubstitutor: replaced SMILES for %r (was %r)",
                    label, smi,
                )
            effective_smi_map[label] = effective_smi

        # Carrier swaps collapse to a null reaction under substitution; short-circuit before equilibrator.
        smiles_coeffs: Dict[str, float] = {}
        all_have_smiles = True
        for label, coeff in stoichiometry.items():
            can = canonical_smiles(effective_smi_map.get(label))
            if not can:
                all_have_smiles = False
                break
            smiles_coeffs[can] = smiles_coeffs.get(can, 0) + coeff
        if all_have_smiles and not any(c != 0 for c in smiles_coeffs.values()):
            return "null", None

        compound_map: Dict[str, object] = {}
        for label in stoichiometry:
            effective_smi = effective_smi_map[label]
            cpd = self._resolve_compound(effective_smi)
            if cpd is None:
                logger.debug(
                    "Cannot compute ΔG°': compound not resolved for label=%r (smiles=%r)",
                    label, effective_smi,
                )
                return None, None
            compound_map[label] = cpd

        # Accumulate coeffs per CC compound (don't overwrite when labels resolve to the same id).
        accumulated: Dict[object, float] = {}
        for lab, coeff in stoichiometry.items():
            cpd = compound_map[lab]
            accumulated[cpd] = accumulated.get(cpd, 0) + coeff
        # Drop net-zero spectators from CoA substitution.
        accumulated = {cpd: c for cpd, c in accumulated.items() if c != 0}
        if not accumulated:
            return "null", None
        rxn = Reaction(accumulated)

        # Check balance before computing — unbalanced reactions return garbage ΔG°' with no warning.
        if not rxn.is_balanced():
            logger.warning(
                "Reaction is not balanced (atoms/charge); skipping ΔG°' computation. "
                "Stoichiometry: %s", dict(stoichiometry)
            )
            return None, None
        return "ok", rxn

    def _compute_dgr_prime(
        self,
        stoichiometry: Dict[str, int],
        smiles_map: Dict[str, str],
    ) -> Tuple[Optional[float], Optional[float]]:
        """Build equilibrator Reaction and return (ΔG°′ kJ/mol, σ kJ/mol), or (None, None)."""
        cc = self._cc
        if cc is None:
            return None, None

        try:
            status, rxn = self._build_cc_reaction(stoichiometry, smiles_map)
            if status == "null":
                return 0.0, 0.0
            if status is None:
                return None, None

            res = cc.standard_dg_prime(rxn)
            dgr = float(res.value.m_as("kJ/mol"))
            sigma = float(res.error.m_as("kJ/mol"))
            return dgr, sigma

        except Exception as e:
            logger.warning(
                "standard_dg_prime failed for stoichiometry %s: %s",
                dict(stoichiometry), e,
            )
            return None, None

    def joint_dgr_prime(
        self,
        stoichiometries: List[Dict[str, int]],
        smiles_maps: List[Dict[str, str]],
    ):
        """Joint ΔG°′ (kJ/mol) and covariance ((kJ/mol)²) for several reactions.

        Reactions share group-contribution parameters, so their ΔG°′ errors are correlated;
        the scalar σ from ``compute_keq`` drops that. Null (carrier-swap) reactions get
        ΔG°′ = 0 with zero variance; unbuildable ones get NaN, zero variance and
        ``resolved=False``. Returns None when equilibrator is unavailable or disabled.
        Not cached.
        """
        import numpy as np

        if _is_disabled():
            return None
        cc = self._load_cc()
        if cc is None:
            return None

        n = len(stoichiometries)
        dg = np.full(n, np.nan)
        cov = np.zeros((n, n))
        resolved = np.zeros(n, dtype=bool)
        cc_rxns = []
        cc_idx = []
        for i, (stoich, smap) in enumerate(zip(stoichiometries, smiles_maps)):
            try:
                status, rxn = self._build_cc_reaction(stoich, smap)
            except Exception as e:
                logger.warning("joint ΔG°': cannot build reaction %s: %s", dict(stoich), e)
                continue
            if status == "null":
                dg[i] = 0.0
                resolved[i] = True
            elif status == "ok":
                cc_rxns.append(rxn)
                cc_idx.append(i)

        if cc_rxns:
            values, covariance = cc.standard_dg_prime_multi(
                cc_rxns, uncertainty_representation="cov"
            )
            vals = np.asarray(values.m_as("kJ/mol"), dtype=float).ravel()
            cov_sub = np.asarray(covariance.m_as("kJ**2/mol**2"), dtype=float)
            idx = np.asarray(cc_idx)
            dg[idx] = vals
            cov[np.ix_(idx, idx)] = cov_sub
            resolved[idx] = True
        return dg, cov, resolved

    def compute_keq(
        self,
        stoichiometry: Dict[str, int],
        smiles_map: Dict[str, str],
    ) -> ThermoData:
        """Compute ΔG°′ and Keq. kcat_rev is always None; irreversibility is decided by the rule layer."""
        if _is_disabled():
            return _disabled_thermo()

        cc = self._load_cc()
        if cc is None:
            return _fallback_thermo("equilibrator unavailable")

        canonical_pairs = sorted(
            ((canonical_smiles(smiles_map.get(lab)) or lab, c)
             for lab, c in stoichiometry.items()),
            key=lambda p: p[0],
        )
        rxn_key = json.dumps(
            [canonical_pairs, self.pH, self.ionic_strength_M, self.pMg, self.T_K],
            default=str,
        )

        cached = self.cache.get("reaction", rxn_key, default=...)
        if cached is not ... and cached is not None and cached != (None, None):
            result = cached
        else:
            result = self._compute_dgr_prime(stoichiometry, smiles_map)
            # Cache reaction hits only (substitutor is not in the key; caching misses would trap later runs).
            if result is not None and result != (None, None):
                self.cache.set("reaction", rxn_key, result)

        if result is None or result == (None, None):
            stoich_str = " + ".join(
                f"{c} {lab}" for lab, c in stoichiometry.items() if c != 0
            )
            return _fallback_thermo(
                f"equilibrator returned no ΔG°' for: {stoich_str}"
            )

        dgr, sigma = result

        if sigma is not None and sigma > _SIGMA_WARN_KJMOL:
            logger.warning(
                "High ΔG°' uncertainty (σ=%.1f kJ/mol) for reaction %s — "
                "Keq may be unreliable.",
                sigma, canonical_pairs,
            )

        try:
            keq = math.exp(-dgr * 1000 / (R * self.T_K))
        except OverflowError:
            keq = float("inf")

        # Irreversibility is rule-layer; compute_keq returns irreversible=False.
        td = ThermoData(
            dgr_prime_kJmol=float(dgr),
            sigma_kJmol=float(sigma) if sigma is not None else float("nan"),
            keq=float(keq),
            kcat_rev=None,
            irreversible=False,
            source="equilibrator",
        )
        logger.debug(
            "thermo: ΔG°'=%.2f kJ/mol (σ=%.2f), Keq=%.3g, irreversible=%s",
            td.dgr_prime_kJmol, td.sigma_kJmol, td.keq, td.irreversible,
        )
        return td
