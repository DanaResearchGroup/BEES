#!/usr/bin/env python3

"""Reversible enzyme inhibition: rate-law factors and a vectorised evaluator.

Concentrations mM. Each inhibitor I of a reaction binds the free enzyme (constant Ki)
and/or the enzyme-substrate complex (constant Ki'), so against a substrate S

    v = Vmax * S / (Km * a + S * b),   a = 1 + sum_i I_i/Ki_i,   b = 1 + sum_i I_i/Ki'_i

which is the Michaelis-Menten law with

    Vmax_app = Vmax / b,   Km_app = Km * a / b

    competitive     Ki' = inf          Km_app = Km (1 + I/Ki), Vmax unchanged
    uncompetitive   Ki  = inf          both scaled by 1 / (1 + I/Ki')
    noncompetitive  Ki' = Ki           Vmax_app = Vmax / (1 + I/Ki), Km unchanged
    mixed           Ki' = alpha * Ki   alpha -> inf is competitive, alpha = 1 noncompetitive
    tight_binding   Morrison: the inhibitor depletes itself against [E]; Ki is the apparent
                    constant and Vmax scales with the free-enzyme fraction

Inhibitors that share a site add in a and b; different substrates and tight binders multiply.
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

INHIBITION_TYPES = ("competitive", "uncompetitive", "noncompetitive", "mixed", "tight_binding")
# Types whose effect depends on which substrate they act against.
SUBSTRATE_SPECIFIC_TYPES = ("competitive", "uncompetitive", "mixed")


@dataclass(frozen=True)
class InhibitorSpec:
    """One resolved inhibitor of a reaction (attached by the enlarger as ``rxn.inhibitors``)."""

    inhibitor: str  # species label, lowercase
    type: str
    ki: float  # mM
    alpha: Optional[float] = None  # mixed only
    substrate: Optional[str] = None  # lowercase; None = the reaction's first substrate with a Km


def inverse_site_constants(type: str, ki: float, alpha: Optional[float] = None) -> Tuple[float, float]:
    """(1/Ki, 1/Ki'): affinity for the free enzyme and for the enzyme-substrate complex (0 = no binding)."""
    if ki <= 0 or not math.isfinite(ki):
        raise ValueError(f"Ki must be positive and finite, got {ki}")
    k = 1.0 / ki
    if type == "competitive":
        return k, 0.0
    if type == "uncompetitive":
        return 0.0, k
    if type == "noncompetitive":
        return k, k
    if type == "mixed":
        if alpha is None or alpha <= 0:
            raise ValueError("mixed inhibition needs alpha > 0")
        return k, k / alpha
    raise ValueError(f"No Michaelis-Menten site constants for inhibition type {type!r}")


def mm_scales(type: str, ki: float, inhibitor: float, alpha: Optional[float] = None) -> Tuple[float, float]:
    """(Vmax factor, Km factor) of one inhibitor at concentration ``inhibitor``."""
    inv_ka, inv_kb = inverse_site_constants(type, ki, alpha)
    i = max(inhibitor, 0.0)
    a = 1.0 + i * inv_ka
    b = 1.0 + i * inv_kb
    return 1.0 / b, a / b


def inhibited_mm_rate(
    vmax: float, km: float, s: float, inhibitor: float, type: str, ki: float,
    alpha: Optional[float] = None,
) -> float:
    """Single-substrate Michaelis-Menten rate with one reversible inhibitor."""
    f_v, f_km = mm_scales(type, ki, inhibitor, alpha)
    s = max(s, 0.0)
    return vmax * f_v * s / (km * f_km + s)


def morrison_free_fraction(inhibitor, enzyme, ki):
    """Fraction of free enzyme (v/v0) for a tight binder; ``inhibitor`` is the total concentration.

    1 - EI/E with EI = 2 E I / (E + I + Ki + sqrt((E + I + Ki)^2 - 4 E I)), written without
    the subtraction of the usual closed form, so it stays accurate for E << Ki. With E = 0 it
    is exactly the classic 1 / (1 + I/Ki).
    """
    i = np.maximum(inhibitor, 0.0)
    s = enzyme + i + ki
    disc = np.maximum(s * s - 4.0 * enzyme * i, 0.0)
    return 1.0 - 2.0 * i / (s + np.sqrt(disc))


class InhibitionTable:
    """Vectorised (Vmax factor, Km factors) of every inhibited reaction at a state ``y``.

    ``entries`` are (reaction index, inhibitor species index, slot, InhibitorSpec, enzyme conc);
    ``slot`` is the column of the substrate matrix the inhibitor acts on, or ``n_slots`` for
    inhibitors that act on the enzyme as a whole (noncompetitive, tight binding).
    """

    def __init__(self, n_rx: int, n_slots: int, entries: Sequence[tuple]):
        self.n_rx = n_rx
        self.n_slots = n_slots
        mm = [e for e in entries if e[3].type != "tight_binding"]
        tight = [e for e in entries if e[3].type == "tight_binding"]
        self.active = bool(entries)

        self._mm_rx = np.array([e[0] for e in mm], dtype=np.int64)
        self._mm_sp = np.array([e[1] for e in mm], dtype=np.int64)
        self._mm_slot = np.array([e[2] for e in mm], dtype=np.int64)
        consts = [inverse_site_constants(e[3].type, e[3].ki, e[3].alpha) for e in mm]
        self._inv_ka = np.array([c[0] for c in consts], dtype=np.float64)
        self._inv_kb = np.array([c[1] for c in consts], dtype=np.float64)

        self._t_rx = np.array([e[0] for e in tight], dtype=np.int64)
        self._t_sp = np.array([e[1] for e in tight], dtype=np.int64)
        self._t_ki = np.array([e[3].ki for e in tight], dtype=np.float64)
        self._t_e = np.array([e[4] for e in tight], dtype=np.float64)
        self.inhibited_reactions = sorted({e[0] for e in entries})

    def scales(self, y_ext: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """(v_scale (n_rx[, n_t]), km_scale (n_rx, n_slots[, n_t])) for y_ext of shape (n_sp+1[, n_t])."""
        extra = y_ext.shape[1:]
        v_scale = np.ones((self.n_rx,) + extra)
        km_scale = np.ones((self.n_rx, self.n_slots) + extra)
        if not self.active:
            return v_scale, km_scale

        if len(self._mm_rx):
            conc = np.maximum(y_ext[self._mm_sp], 0.0)
            inv_ka = self._inv_ka.reshape((-1,) + (1,) * len(extra))
            inv_kb = self._inv_kb.reshape((-1,) + (1,) * len(extra))
            a = np.zeros((self.n_rx, self.n_slots + 1) + extra)
            b = np.zeros_like(a)
            np.add.at(a, (self._mm_rx, self._mm_slot), conc * inv_ka)
            np.add.at(b, (self._mm_rx, self._mm_slot), conc * inv_kb)
            a += 1.0
            b += 1.0
            km_scale = a[:, :self.n_slots] / b[:, :self.n_slots]
            v_scale = np.prod(1.0 / b, axis=1)

        if len(self._t_rx):
            conc = y_ext[self._t_sp]
            e = self._t_e.reshape((-1,) + (1,) * len(extra))
            ki = self._t_ki.reshape((-1,) + (1,) * len(extra))
            np.multiply.at(v_scale, self._t_rx, morrison_free_fraction(conc, e, ki))
        return v_scale, km_scale


def resolve_slot(substrate: Optional[str], slot_labels: List[str], alias_to_model: Dict[str, str]) -> Optional[int]:
    """Column of a reaction's substrate list an inhibitor acts on (None when not a substrate of it)."""
    if not slot_labels:
        return None
    if substrate is None:
        return 0
    want = alias_to_model.get(substrate, substrate)
    for k, lab in enumerate(slot_labels):
        lc = lab.lower().strip()
        if lc == substrate or alias_to_model.get(lc, lc) == want:
            return k
    return None


def morrison_expression(inhibitor: str, enzyme: str, ki: str) -> str:
    """Formula text of morrison_free_fraction for exporters (arguments are tokens)."""
    s = f"({enzyme} + {inhibitor} + {ki})"
    return f"(1 - 2 * {inhibitor} / ({s} + sqrt({s}^2 - 4 * {enzyme} * {inhibitor})))"


def format_inhibition_terms(
    specs: Sequence[InhibitorSpec],
    conc,
    param,
    first_substrate: Optional[str],
) -> Tuple[str, Dict[str, str]]:
    """Formula text of the inhibition factors, for exporters.

    ``conc(label)`` returns the inhibitor's concentration token (None when it is not a model
    species), ``param(kind, spec, value)`` registers the constant ``kind`` ("Ki" or "Ki'") of
    ``spec`` and returns its token. Returns
    (Vmax factor expression, {lowercase substrate label: Km factor expression}) mirroring ``scales``;
    tight binders are not handled here (they need [E]).
    """
    sites: Dict[Optional[str], Tuple[List[str], List[str]]] = {}
    for sp in specs:
        if sp.type == "tight_binding":
            continue
        token = conc(sp.inhibitor)
        if token is None:
            continue
        inv_ka, inv_kb = inverse_site_constants(sp.type, sp.ki, sp.alpha)
        sub = None
        if sp.type in SUBSTRATE_SPECIFIC_TYPES:
            sub = sp.substrate if sp.substrate is not None else first_substrate
            if sub is None:
                continue
            sub = sub.lower().strip()
        a_terms, b_terms = sites.setdefault(sub, ([], []))
        if inv_ka:
            a_terms.append(token + " / " + param("Ki", sp, 1.0 / inv_ka))
        if inv_kb:
            b_terms.append(token + " / " + param("Ki'", sp, 1.0 / inv_kb))

    v_parts: List[str] = []
    km_factor: Dict[str, str] = {}
    for sub, (a_terms, b_terms) in sites.items():
        a = " + ".join(["1"] + a_terms)
        b = " + ".join(["1"] + b_terms)
        if b_terms:
            v_parts.append(f"(1 / ({b}))")
        if sub is not None and (a_terms or b_terms):
            km_factor[sub] = f"(({a}) / ({b}))"
    return (" * ".join(v_parts) if v_parts else "1"), km_factor
