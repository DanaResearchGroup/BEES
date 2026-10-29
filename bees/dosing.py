#!/usr/bin/env python3

"""Time-dependent dosing: release profiles of a carrier and first-order elimination.

Pure functions of time. Units: mM, s. A profile is the fraction F(tau) of the dose released
tau seconds after ``t_start`` (F -> 1), so the species receives the additive source
u(t) = dose * dF/dt, and optionally loses ``elimination_k * [species]``.

    bolus             all of the dose at t_start (an impulse, not a rate)
    constant          duration                     F = tau / duration, then 0
    first_order       k | half_life                F = 1 - exp(-k tau)
    biexponential     fraction_fast,               F = 1 - [phi exp(-k_f tau) + (1 - phi) exp(-k_s tau)]
                      k_fast | half_life_fast,
                      k_slow | half_life_slow
    weibull           tau, beta                    F = 1 - exp(-(tau / tau_c)^beta)
    korsmeyer_peppas  k, n                         F = min(k tau^n, 1)
    table             table [[t, rate], ...]       piecewise-linear rate, scaled to integrate to the dose

The Weibull (beta < 1) and Korsmeyer-Peppas (n < 1) rates are unbounded at tau = 0. Their
first ``SINGULAR_START`` seconds are delivered as an impulse of dose * F(SINGULAR_START), and
the rate is used from there on.
"""

import math
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

PROFILES = (
    "bolus", "constant", "first_order", "biexponential", "weibull", "korsmeyer_peppas", "table",
)

# Seconds delivered as an impulse for profiles whose rate diverges at tau = 0.
SINGULAR_START = 1e-3


def _rate_from(params: Mapping, key: str, what: str) -> float:
    """A first-order rate constant given as ``key`` (1/s) or ``half_life``-style alternative."""
    half = {"k": "half_life", "k_fast": "half_life_fast", "k_slow": "half_life_slow"}[key]
    if params.get(key) is not None:
        value = float(params[key])
    elif params.get(half) is not None:
        t_half = float(params[half])
        if t_half <= 0:
            raise ValueError(f"{what}: {half} must be positive")
        value = math.log(2.0) / t_half
    else:
        raise ValueError(f"{what}: give '{key}' (1/s) or '{half}' (s)")
    if value <= 0:
        raise ValueError(f"{what}: {key} must be positive")
    return value


def normalize_params(profile: str, params: Optional[Mapping]) -> Dict[str, object]:
    """Validate ``params`` and reduce them to the canonical keys the profile functions use."""
    params = dict(params or {})
    what = f"dosing profile '{profile}'"
    if profile not in PROFILES:
        raise ValueError(f"Unknown dosing profile '{profile}'. Choose from {list(PROFILES)}")
    if profile == "bolus":
        return {}
    if profile == "constant":
        duration = float(params.get("duration", 0.0))
        if duration <= 0:
            raise ValueError(f"{what}: 'duration' (s) must be positive")
        return {"duration": duration}
    if profile == "first_order":
        return {"k": _rate_from(params, "k", what)}
    if profile == "biexponential":
        phi = params.get("fraction_fast")
        if phi is None or not 0.0 <= float(phi) <= 1.0:
            raise ValueError(f"{what}: 'fraction_fast' must be in [0, 1]")
        return {
            "fraction_fast": float(phi),
            "k_fast": _rate_from(params, "k_fast", what),
            "k_slow": _rate_from(params, "k_slow", what),
        }
    if profile == "weibull":
        tau, beta = float(params.get("tau", 0.0)), float(params.get("beta", 0.0))
        if tau <= 0 or beta <= 0:
            raise ValueError(f"{what}: 'tau' (s) and 'beta' must be positive")
        return {"tau": tau, "beta": beta}
    if profile == "korsmeyer_peppas":
        k, n = float(params.get("k", 0.0)), float(params.get("n", 0.0))
        if k <= 0 or n <= 0:
            raise ValueError(f"{what}: 'k' and 'n' must be positive")
        return {"k": k, "n": n}
    rows = params.get("table")
    if not rows or len(rows) < 2:
        raise ValueError(f"{what}: 'table' needs at least two [t, rate] rows")
    arr = np.array([[float(x) for x in row] for row in rows], dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ValueError(f"{what}: 'table' rows must be [t, rate] pairs")
    if np.any(arr[:, 0] < 0) or np.any(np.diff(arr[:, 0]) <= 0):
        raise ValueError(f"{what}: table times must be >= 0 and strictly increasing")
    if np.any(arr[:, 1] < 0) or not np.sum(0.5 * (arr[1:, 1] + arr[:-1, 1]) * np.diff(arr[:, 0])) > 0:
        raise ValueError(f"{what}: table rates must be >= 0 and release something")
    return {"table": arr}


def _table_cumulative(tab: np.ndarray, tau: np.ndarray) -> np.ndarray:
    """Integral of the piecewise-linear rate from tab[0, 0] to tau (unnormalised)."""
    t, r = tab[:, 0], tab[:, 1]
    seg = np.concatenate([[0.0], np.cumsum(0.5 * (r[1:] + r[:-1]) * np.diff(t))])
    tau_c = np.clip(tau, t[0], t[-1])
    i = np.clip(np.searchsorted(t, tau_c, side="right") - 1, 0, len(t) - 2)
    dt = tau_c - t[i]
    slope = (r[i + 1] - r[i]) / (t[i + 1] - t[i])
    return seg[i] + r[i] * dt + 0.5 * slope * dt * dt


def cumulative_fraction(profile: str, params: Mapping, tau) -> np.ndarray:
    """Fraction of the dose released tau seconds after t_start (params as from normalize_params)."""
    tau = np.asarray(tau, dtype=float)
    pos = np.maximum(tau, 0.0)
    if profile == "bolus":
        return np.where(tau >= 0.0, 1.0, 0.0)
    if profile == "constant":
        return np.minimum(pos / params["duration"], 1.0)
    if profile == "first_order":
        return -np.expm1(-params["k"] * pos)
    if profile == "biexponential":
        phi = params["fraction_fast"]
        return 1.0 - (phi * np.exp(-params["k_fast"] * pos) + (1.0 - phi) * np.exp(-params["k_slow"] * pos))
    if profile == "weibull":
        return -np.expm1(-((pos / params["tau"]) ** params["beta"]))
    if profile == "korsmeyer_peppas":
        return np.minimum(params["k"] * pos ** params["n"], 1.0)
    tab = params["table"]
    return _table_cumulative(tab, tau) / _table_cumulative(tab, np.array(tab[-1, 0]))


def release_rate(profile: str, params: Mapping, tau) -> np.ndarray:
    """dF/dtau (1/s); 0 for tau <= 0 (a bolus has no rate, it is an impulse)."""
    tau = np.asarray(tau, dtype=float)
    pos = np.maximum(tau, 0.0)
    live = tau > 0.0
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        if profile == "bolus":
            rate = np.zeros_like(tau)
        elif profile == "constant":
            rate = np.where(tau < params["duration"], 1.0 / params["duration"], 0.0)
        elif profile == "first_order":
            rate = params["k"] * np.exp(-params["k"] * pos)
        elif profile == "biexponential":
            phi = params["fraction_fast"]
            rate = (phi * params["k_fast"] * np.exp(-params["k_fast"] * pos)
                    + (1.0 - phi) * params["k_slow"] * np.exp(-params["k_slow"] * pos))
        elif profile == "weibull":
            x = pos / params["tau"]
            rate = params["beta"] / params["tau"] * x ** (params["beta"] - 1.0) * np.exp(-(x ** params["beta"]))
        elif profile == "korsmeyer_peppas":
            n = params["n"]
            t_end = params["k"] ** (-1.0 / n)
            rate = np.where(tau < t_end, params["k"] * n * pos ** (n - 1.0), 0.0)
        else:
            tab = params["table"]
            area = _table_cumulative(tab, np.array(tab[-1, 0]))
            inside = (tau >= tab[0, 0]) & (tau <= tab[-1, 0])
            rate = np.where(inside, np.interp(tau, tab[:, 0], tab[:, 1]) / area, 0.0)
    return np.where(live, rate, 0.0)


def is_singular_at_start(profile: str, params: Mapping) -> bool:
    """True when the release rate diverges at tau = 0 (Weibull beta < 1, Korsmeyer-Peppas n < 1)."""
    if profile == "weibull":
        return params["beta"] < 1.0
    if profile == "korsmeyer_peppas":
        return params["n"] < 1.0
    return False


def release_end(profile: str, params: Mapping) -> Optional[float]:
    """tau after which nothing more is released (None for profiles with an infinite tail)."""
    if profile == "bolus":
        return 0.0
    if profile == "constant":
        return float(params["duration"])
    if profile == "korsmeyer_peppas":
        return float(params["k"] ** (-1.0 / params["n"]))
    if profile == "table":
        return float(params["table"][-1, 0])
    return None


@dataclass(frozen=True)
class DoseSchedule:
    """One dosing entry: ``dose`` mM of ``species`` released from ``t_start`` by ``profile``."""

    species: str
    profile: str
    dose: float
    params: Mapping = field(default_factory=dict)
    t_start: float = 0.0
    elimination_k: float = 0.0  # 1/s, first-order loss of the free species

    def __post_init__(self):
        if self.dose < 0:
            raise ValueError(f"Dose of '{self.species}' must be >= 0")
        object.__setattr__(self, "params", normalize_params(self.profile, self.params))

    def with_time_scale(self, factor: float) -> "DoseSchedule":
        """Copy whose release times are multiplied by ``factor`` (elimination and t_start unchanged)."""
        p = dict(self.params)
        for key in ("duration", "tau"):
            if key in p:
                p[key] *= factor
        for key in ("k", "k_fast", "k_slow"):
            if key in p and self.profile != "korsmeyer_peppas":
                p[key] /= factor
        if self.profile == "korsmeyer_peppas":
            p["k"] = p["k"] * factor ** (-p["n"])
        if "table" in p:
            tab = p["table"].copy()
            tab[:, 0] *= factor
            p["table"] = tab
        return DoseSchedule(self.species, self.profile, self.dose, p, self.t_start, self.elimination_k)

    @property
    def singular(self) -> bool:
        return is_singular_at_start(self.profile, self.params)

    def released(self, t) -> np.ndarray:
        """Cumulative amount (mM) released by absolute time ``t``."""
        return self.dose * cumulative_fraction(self.profile, self.params, np.asarray(t, dtype=float) - self.t_start)

    def rate(self, t) -> np.ndarray:
        """Release rate (mM/s) at absolute time ``t``; the impulse parts are not included."""
        tau = np.asarray(t, dtype=float) - self.t_start
        rate = self.dose * release_rate(self.profile, self.params, tau)
        if self.singular:
            rate = np.where(tau >= SINGULAR_START, rate, 0.0)
        return rate

    def impulses(self) -> List[Tuple[float, float]]:
        """(absolute time, amount in mM) of the parts delivered instantaneously."""
        if self.profile == "bolus":
            return [(self.t_start, self.dose)]
        if self.singular:
            return [(self.t_start, self.dose * float(cumulative_fraction(self.profile, self.params, SINGULAR_START)))]
        return []

    def breakpoints(self) -> List[float]:
        """Absolute times at which the rate is not smooth."""
        pts = [self.t_start]
        if self.singular:
            pts.append(self.t_start + SINGULAR_START)
        end = release_end(self.profile, self.params)
        if end is not None:
            pts.append(self.t_start + end)
        if self.profile == "table":
            pts.extend(self.t_start + self.params["table"][:, 0])
        return sorted({float(p) for p in pts})


def schedules_from_spec(entries: Sequence) -> List[DoseSchedule]:
    """DoseSchedules from validated ``dosing`` input entries (see schema.Dosing)."""
    out = []
    for e in entries or []:
        k = getattr(e, "elimination_k", None)
        if k is None and getattr(e, "elimination_half_life", None) is not None:
            k = math.log(2.0) / float(e.elimination_half_life)
        out.append(DoseSchedule(
            species=e.species, profile=e.profile, dose=float(e.dose),
            params=dict(getattr(e, "params", None) or {}),
            t_start=float(getattr(e, "t_start", 0.0) or 0.0),
            elimination_k=float(k or 0.0),
        ))
    return out


class DosePlan:
    """A set of schedules laid over one simulation window of the species vector.

    The window starts at absolute time ``t_offset`` (a resumed run), so run time t is
    absolute time t + t_offset. ``y0_add`` is what is already delivered at the start,
    ``events`` are impulses inside the window, ``breakpoints`` the run times at which the
    source is not smooth, and ``elimination`` the per-species first-order loss.
    """

    def __init__(self, schedules: Sequence[DoseSchedule], label_to_idx: Mapping[str, int],
                 n_species: int, t_offset: float = 0.0):
        self.t_offset = float(t_offset)
        self.n_species = n_species
        self.y0_add = np.zeros(n_species)
        self.elimination = np.zeros(n_species)
        self.events: List[Tuple[float, int, float]] = []
        self.missing: List[str] = []
        self._sources: List[Tuple[int, DoseSchedule]] = []
        points = set()
        for sch in schedules:
            idx = label_to_idx.get(sch.species.lower().strip())
            if idx is None:
                if sch.species not in self.missing:
                    self.missing.append(sch.species)
                continue
            if sch.elimination_k:
                self.elimination[idx] = sch.elimination_k
            for t_imp, amount in sch.impulses():
                if t_imp <= self.t_offset:
                    if self.t_offset == 0.0:
                        self.y0_add[idx] += amount
                else:
                    self.events.append((t_imp - self.t_offset, idx, amount))
            if sch.profile != "bolus" and sch.dose > 0:
                self._sources.append((idx, sch))
            points.update(t - self.t_offset for t in sch.breakpoints())
        self.events.sort()
        self.breakpoints = sorted(p for p in points if p > 0.0)
        self.active = bool(self._sources or self.events or self.y0_add.any() or self.elimination.any())

    def source(self, t: float) -> np.ndarray:
        """Additive source u(t) (mM/s) for every species at run time ``t``."""
        u = np.zeros(self.n_species)
        T = t + self.t_offset
        for idx, sch in self._sources:
            u[idx] += float(sch.rate(T))
        return u

    def dydt(self, t: float, y: np.ndarray) -> np.ndarray:
        """Source minus elimination: what the doses add to dy/dt."""
        out = self.source(t)
        if self.elimination.any():
            out -= self.elimination * y
        return out

    def next_break(self, t: float, end_time: float) -> float:
        """The next run time after ``t`` where the source or an impulse changes (capped at end_time)."""
        for p in self.breakpoints:
            if p > t * (1.0 + 1e-12) + 1e-15:
                return min(p, end_time)
        return end_time

    def impulses_at(self, t: float) -> List[Tuple[int, float]]:
        """(species index, amount) of the impulses due at run time ``t``."""
        tol = 1e-9 * max(1.0, abs(t))
        return [(i, a) for te, i, a in self.events if abs(te - t) <= tol]
