#!/usr/bin/env python3

"""
Enzyme Sharing Module
---------------------
Factors applied on top of each reaction's own rate law (RateLawOptions):

- Shared-enzyme competition: reactions of one enzyme draw on the same free enzyme.
  Rapid-equilibrium partition function built from the reactions' own CM terms:
  each side prod(1 + L/K)^nu expands into bound forms (complexes); the core
  reactions' complexes are merged, each distinct complex counted once, and

      Z_r = D_r + sum_{k not in own(r)} w_k,   c_r = D_r / Z_r,
      w_k = mult_k * prod_{L in k} ([L] / K_canon[e, L])^m_L.

  D_r is the denominator reaction r's law actually divides by; its own complexes
  are already inside D_r with its own Km, so they are not added again.
"""

import itertools
import logging
import math
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

from bees.flux_calculator import _km_for, rate_law_spec

logger = logging.getLogger("BEES")

# Warn once per process about the same enzyme/ligand Km disagreement.
_WARNED: Set[Tuple[str, str]] = set()

# Canonical K spread (max/min) above which a disagreement is logged.
_KM_SPREAD_WARN = 2.0

Complex = Tuple[Tuple[int, int], ...]  # sorted ((species_idx, copies), ...)


def _species_index(label: str, label_to_idx: Dict[str, int], alias: Dict[str, str]) -> Optional[int]:
    lc = label.lower().strip()
    return label_to_idx.get(alias.get(lc, lc), label_to_idx.get(lc))


def _side_complexes(side: Sequence[Tuple[int, int]]) -> List[Complex]:
    """Every bound form of one CM side prod(1 + L/K)^nu (free enzyme excluded)."""
    out = []
    for copies in itertools.product(*(range(nu + 1) for _, nu in side)):
        k = tuple(sorted((idx, m) for (idx, _), m in zip(side, copies) if m > 0))
        if k:
            out.append(k)
    return out


def _law_sides(rxn, enzyme_conc_map, label_to_idx, alias):
    """(substrate side, product side) of the law the simulator uses: [(idx, Km, nu)]."""
    spec = rate_law_spec(rxn, enzyme_conc_map)
    if spec.form == "zero":
        return None

    def side(terms):
        out = []
        for lab, km, nu in terms:
            idx = _species_index(lab, label_to_idx, alias)
            if idx is not None and km > 0:
                out.append((idx, float(km), int(round(nu))))
        return out

    return side(spec.substrates), side(spec.products)


class PartitionCompetition:
    """Per-enzyme partition function; see module docstring."""

    def __init__(
        self,
        reactions: list,
        n_core_reactions: int,
        species_labels: List[str],
        label_to_idx: Dict[str, int],
        alias: Dict[str, str],
        enzyme_conc_map: Dict[str, float],
        include_products: bool = True,
    ):
        n_sp = len(species_labels)
        self._labels = list(species_labels)
        self.include_products = include_products

        sides: Dict[int, tuple] = {}
        enzyme_of: Dict[int, str] = {}
        for j, rxn in enumerate(reactions):
            s = _law_sides(rxn, enzyme_conc_map, label_to_idx, alias)
            if s is None:
                continue
            sides[j] = s
            enzyme_of[j] = str(rxn.enzyme_label).lower().strip()

        # Canonical K and the largest nu per (enzyme, ligand) from core reactions only.
        sub_km: Dict[Tuple[str, int], List[float]] = {}
        prod_km: Dict[Tuple[str, int], List[float]] = {}
        nu_max: Dict[Tuple[str, int], int] = {}
        pool: Dict[str, Dict[Complex, bool]] = {}  # enzyme -> complex -> on a substrate side
        for j, (subs, prods) in sides.items():
            if j >= n_core_reactions:
                continue
            e = enzyme_of[j]
            for terms, store, is_sub in ((subs, sub_km, True), (prods, prod_km, False)):
                for idx, km, nu in terms:
                    store.setdefault((e, idx), []).append(km)
                    if nu_max.get((e, idx), nu) != nu:
                        logger.debug(f"Competition: {e}/{species_labels[idx]} nu disagrees; using the larger")
                    nu_max[(e, idx)] = max(nu_max.get((e, idx), 0), nu)
                for k in _side_complexes([(idx, nu) for idx, _, nu in terms]):
                    pool.setdefault(e, {})
                    pool[e][k] = pool[e].get(k, False) or is_sub

        k_canon: Dict[Tuple[str, int], float] = {}
        for key in set(sub_km) | set(prod_km):
            kms = sub_km.get(key) or prod_km[key]
            k_canon[key] = math.exp(sum(math.log(k) for k in kms) / len(kms))
            if max(kms) > _KM_SPREAD_WARN * min(kms) and key not in _WARNED:
                _WARNED.add(key)
                logger.warning(
                    f"Competition: {key[0]} Km for {species_labels[key[1]]} spans "
                    f"{min(kms):.3g}-{max(kms):.3g} mM; using geometric mean {k_canon[key]:.3g}"
                )

        self.groups = []
        for e, complexes in pool.items():
            keys = [k for k, on_sub in complexes.items() if on_sub or include_products]
            if not keys:
                continue
            members, notown = [], []
            for j, (subs, prods) in sides.items():
                if enzyme_of[j] != e:
                    continue
                own = set(_side_complexes([(i, nu) for i, _, nu in subs]))
                own |= set(_side_complexes([(i, nu) for i, _, nu in prods]))
                row = [k not in own for k in keys]
                if any(row):
                    members.append(j)
                    notown.append(row)
            if not members:
                continue
            width = max(len(k) for k in keys)
            lig = np.full((len(keys), width), n_sp, dtype=np.int64)
            kc = np.ones((len(keys), width), dtype=np.float64)
            m = np.zeros((len(keys), width), dtype=np.float64)
            mult = np.ones(len(keys), dtype=np.float64)
            for a, k in enumerate(keys):
                for b, (idx, copies) in enumerate(k):
                    lig[a, b] = idx
                    kc[a, b] = k_canon[(e, idx)]
                    m[a, b] = copies
                    mult[a] *= math.comb(nu_max[(e, idx)], copies)
            self.groups.append({
                "enzyme": e,
                "keys": keys,
                "rows": np.array(members, dtype=np.int64),
                "notown": np.array(notown, dtype=np.float64),
                "lig": lig, "kc": kc, "m": m, "mult": mult,
            })

    def _weights(self, g, y_ext: np.ndarray) -> np.ndarray:
        s = np.maximum(y_ext[g["lig"]], 0.0)
        if s.ndim == 3:
            return g["mult"][:, None] * np.prod((s / g["kc"][:, :, None]) ** g["m"][:, :, None], axis=1)
        return g["mult"] * np.prod((s / g["kc"]) ** g["m"], axis=1)

    def factor(self, y_ext: np.ndarray, den: np.ndarray) -> np.ndarray:
        """c_r for every reaction (1 where the enzyme has no other bound forms)."""
        c = np.ones_like(den)
        for g in self.groups:
            others = g["notown"] @ self._weights(g, y_ext)
            d = den[g["rows"]]
            c[g["rows"]] = d / (d + others)
        return c

    def complex_label(self, k: Complex) -> str:
        return " + ".join(
            self._labels[idx] if copies == 1 else f"{copies}x {self._labels[idx]}" for idx, copies in k
        )

    def other_terms(self, reaction_index: int) -> List[Tuple[float, List[Tuple[str, int, float]]]]:
        """Non-own complexes of one reaction: (multiplicity, [(label, copies, canonical K)])."""
        for g in self.groups:
            pos = np.flatnonzero(g["rows"] == reaction_index)
            if pos.size == 0:
                continue
            out = []
            for a, keep in enumerate(g["notown"][int(pos[0])]):
                if not keep:
                    continue
                ligands = [
                    (self._labels[idx], int(copies), float(g["kc"][a, b]))
                    for b, (idx, copies) in enumerate(g["keys"][a])
                ]
                out.append((float(g["mult"][a]), ligands))
            return out
        return []

    def contributions(self, y_ext: np.ndarray) -> Dict[str, Dict[str, float]]:
        """{enzyme: {complex: w_k}} over the enzyme's merged complexes (y_ext 1-D)."""
        out = {}
        for g in self.groups:
            w = self._weights(g, y_ext)
            out[g["enzyme"]] = {self.complex_label(k): float(x) for k, x in zip(g["keys"], w)}
        return out

    def breakdown(self, y_ext: np.ndarray, den: np.ndarray) -> Dict[int, dict]:
        """Per competing reaction j: enzyme, D_r, and {complex: w_k} of its non-own part."""
        out = {}
        for g in self.groups:
            w = self._weights(g, y_ext)
            for row, j in zip(g["notown"], g["rows"]):
                out[int(j)] = {
                    "enzyme": g["enzyme"],
                    "D": float(den[j]),
                    "others": {
                        self.complex_label(k): float(x)
                        for k, x, keep in zip(g["keys"], w, row) if keep
                    },
                }
        return out


def format_sharing_terms(terms, conc, kname) -> Optional[str]:
    """Sum of complex weights as an infix formula (``^`` = power).

    ``conc(label)`` and ``kname(label)`` name the concentration and canonical K;
    ``conc`` returning None (species absent from the written model) gives None.
    """
    parts = []
    for mult, ligands in terms:
        bits = []
        if mult != 1.0:
            bits.append(f"{mult:g}")
        for label, copies, _k in ligands:
            token = conc(label)
            if token is None:
                return None
            base = f"{token} / {kname(label)}"
            bits.append(base if copies == 1 else f"({base})^{copies}")
        parts.append(" * ".join(bits))
    return " + ".join(parts)


# The earlier experiment patch (acyl substrates only, FAS enzymes); verification only.
_LEGACY_ENZYMES = frozenset({"faba", "fabb", "fabf", "fabg", "fabi", "fabz", "tesa"})


def _legacy_acyl_substrate(rxn) -> Optional[str]:
    """Longest carrier-bound acyl substrate (malonyl donor excluded), or None."""
    from bees.rules.helpers import _CARRIER_SUFFIXES, detect_acyl_chain_length

    smiles_map = getattr(rxn.kinetics, "compound_smiles", None)
    if not isinstance(smiles_map, dict):
        smiles_map = {}
    best, best_n = None, -1
    for lab in rxn.reactant_labels:
        lc = lab.lower().strip()
        if "malonyl" in lc or not lc.endswith(_CARRIER_SUFFIXES):
            continue
        n = detect_acyl_chain_length(lab, smiles_map.get(lab))
        if n is not None and n > best_n:
            best, best_n = lab, n
    return best


class LegacyCompetition:
    """c_r = (1 + x_r) / (1 + x_r + sum of other core x_j), x = (S/Km)^nu on the acyl substrate."""

    def __init__(self, reactions: list, n_core_reactions: int, label_to_idx, alias):
        by_enzyme: Dict[str, list] = {}
        for j, rxn in enumerate(reactions):
            enzyme = str(rxn.enzyme_label).lower().strip()
            if enzyme not in _LEGACY_ENZYMES or rxn.kinetics is None or rxn.rate_law is None:
                continue
            sub = _legacy_acyl_substrate(rxn)
            if sub is None:
                continue
            km = _km_for(rxn, sub)
            if km is None or km <= 0:
                continue
            idx = _species_index(sub, label_to_idx, alias)
            if idx is None:
                continue
            nu = float(abs(rxn.stoichiometry.get(sub, 1)))
            by_enzyme.setdefault(enzyme, []).append((j, idx, float(km), nu))

        self.groups = []
        for members in by_enzyme.values():
            if len(members) < 2:
                continue
            sp_idx = np.array([m[1] for m in members], dtype=np.int64)
            first_core: Dict[int, int] = {}
            for k, (j, idx, _, _) in enumerate(members):
                if j < n_core_reactions and idx not in first_core:
                    first_core[idx] = k
            distinct = np.zeros(len(members), dtype=bool)
            distinct[list(first_core.values())] = True
            rep = np.array([first_core.get(idx, -1) for idx in sp_idx], dtype=np.int64)
            self.groups.append((
                np.array([m[0] for m in members], dtype=np.int64),
                sp_idx,
                np.array([m[2] for m in members], dtype=np.float64),
                np.array([m[3] for m in members], dtype=np.float64),
                distinct,
                rep,
            ))

    def factor(self, y_ext: np.ndarray, shape: tuple) -> np.ndarray:
        c = np.ones(shape)
        is_mat = len(shape) > 1
        for rxn_idx, sp_idx, km, nu, distinct, rep in self.groups:
            s = np.maximum(y_ext[sp_idx], 0.0)
            if is_mat:
                x = (s / km[:, np.newaxis]) ** nu[:, np.newaxis]
                has_rep = (rep >= 0)[:, np.newaxis]
            else:
                x = (s / km) ** nu
                has_rep = rep >= 0
            total = x[distinct].sum(axis=0)
            own = np.where(has_rep, x[np.maximum(rep, 0)], 0.0)
            others = np.maximum(total - own, 0.0)
            c[rxn_idx] = (1.0 + x) / (1.0 + x + others)
        return c
