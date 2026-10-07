"""
Shared-enzyme competition against explicit mass-action mechanisms.

The mechanism's bound forms are generated from the test reactions' CM terms (every subset
of each side's ligands, merged across reactions, each counted once), not chosen by hand.
Binding is rapid (k_off = 1000 x kcat, Km = (k_off + kcat)/k_on); catalysis converts a
reaction's full substrate form into its full product form (or releases free enzyme when
the law has no product side). Rates are compared at quasi-steady state.

To run:  pytest -v tests/test_competition_mechanism.py
"""

import itertools
import math
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from bees.core_edge_model import RateLawOptions
from bees.simulator import _VectorizedRHS
from bees.thermodynamics import ThermoData

E_TOT = 1e-3
KCAT = 1.0
K_OFF = 1000.0 * KCAT


def _spec(subs, prods=None, kcat=KCAT, keq=None):
    return {"subs": dict(subs), "prods": dict(prods or {}), "kcat": kcat, "keq": keq}


def _subsets(names):
    return [frozenset(c) for n in range(1, len(names) + 1) for c in itertools.combinations(sorted(names), n)]


def _mechanism_fluxes(specs, conc, e_tot=E_TOT):
    """Steady-state net catalytic flux per reaction with ligand concentrations clamped."""
    forms = {frozenset()}
    k_of = {}
    for s in specs:
        for side in (s["subs"], s["prods"]):
            forms.update(_subsets(side))
            k_of.update(side)
    forms = sorted(forms, key=lambda f: (len(f), sorted(f)))
    idx = {f: i for i, f in enumerate(forms)}
    n = len(forms)
    m = np.zeros((n, n))  # dx/dt = m @ x

    def edge(src, dst, rate):
        m[idx[dst], idx[src]] += rate
        m[idx[src], idx[src]] -= rate

    for f in forms:
        for lig in k_of:
            g = f | {lig}
            if lig in f or g not in idx:
                continue
            k_on = (K_OFF + KCAT) / k_of[lig]
            edge(f, g, k_on * conc[lig])
            edge(g, f, K_OFF)

    steps = []
    for s in specs:
        full_sub = frozenset(s["subs"])
        full_prod = frozenset(s["prods"])
        edge(full_sub, full_prod, s["kcat"])
        kcat_r = 0.0
        if s["keq"] is not None:
            kcat_r = s["kcat"] * math.prod(s["prods"].values()) / (s["keq"] * math.prod(s["subs"].values()))
            edge(full_prod, full_sub, kcat_r)
        steps.append((idx[full_sub], idx[full_prod], s["kcat"], kcat_r))

    x0 = np.zeros(n)
    x0[idx[frozenset()]] = e_tot
    sol = solve_ivp(lambda t, x: m @ x, (0.0, 50.0 / KCAT), x0, method="BDF",
                    rtol=1e-10, atol=1e-16, jac=m)
    x = sol.y[:, -1]
    return np.array([kf * x[i] - kr * x[j] for i, j, kf, kr in steps])


def _code_rates(specs, conc, options, extra_products=None):
    """The simulator's rates for the same reactions (one enzyme)."""
    rxns = []
    for k, s in enumerate(specs):
        prods = list(s["prods"]) or [f"out{k}"]
        km = {**s["subs"], **s["prods"]}
        thermo = None
        if s["keq"] is not None:
            thermo = ThermoData(dgr_prime_kJmol=-1.0, sigma_kJmol=1.0, keq=s["keq"], kcat_rev=None,
                                irreversible=False, source="equilibrator")
        stoich = {lab: -1 for lab in s["subs"]}
        stoich.update({lab: 1 for lab in prods})
        rxns.append(SimpleNamespace(
            enzyme_label="Enz", reactant_labels=list(s["subs"]), product_labels=prods,
            stoichiometry=stoich, rate_law="Michaelis-Menten",
            kinetics=SimpleNamespace(km=None, kcat=s["kcat"], vmax=None, km_per_substrate=km),
            template=SimpleNamespace(reversible=s["keq"] is not None), thermo=thermo,
            feedback_inhibitors=None,
        ))
    species = sorted({lab for r in rxns for lab in r.stoichiometry})
    rhs = _VectorizedRHS(
        species_labels=species, reactions=rxns,
        alias_to_model_label={lab.lower(): lab.lower() for lab in species},
        enzyme_conc_map={"enz": E_TOT}, constant_mask=np.zeros(len(species), dtype=bool),
        options=options,
    )
    y = np.array([conc.get(lab, 0.0) for lab in species])
    return rhs.compute_v(y), rhs, y


ON = RateLawOptions(enzyme_competition=True)
OFF = RateLawOptions()


def _rel_err(got, ref):
    return np.max(np.abs(got - ref) / np.abs(ref))


def test_two_competing_substrates():
    specs = [_spec({"A": 1.0}), _spec({"B": 2.0})]
    conc = {"A": 2.0, "B": 3.0}
    ref = _mechanism_fluxes(specs, conc)
    assert _rel_err(_code_rates(specs, conc, ON)[0], ref) < 0.01
    assert _rel_err(_code_rates(specs, conc, OFF)[0], ref) > 0.10


def test_reversible_pure_branch():
    """S <-> P1 and S <-> P2 (s = 10, p1 = p2 = 0.5): c = 11.5/12 = 0.96."""
    specs = [_spec({"S": 1.0}, {"P1": 1.0}, keq=100.0), _spec({"S": 1.0}, {"P2": 1.0}, keq=100.0)]
    conc = {"S": 10.0, "P1": 0.5, "P2": 0.5}
    ref = _mechanism_fluxes(specs, conc)
    v, rhs, y = _code_rates(specs, conc, ON)
    assert _rel_err(v, ref) < 0.01
    assert rhs.competition_factor(y)[0] == pytest.approx(11.5 / 12.0, rel=1e-12)


def _fabg():
    """oxo + NADPH <-> hyd + NADP for two chain lengths; all K = 1 so concentrations are ratios."""
    return [
        _spec({"oxo": 1.0, "NADPH": 1.0}, {"hyd": 1.0, "NADP": 1.0}, keq=340.0),
        _spec({"oxo2": 1.0, "NADPH": 1.0}, {"hyd2": 1.0, "NADP": 1.0}, keq=340.0),
    ]


@pytest.mark.parametrize("case, nadp, k2, h2, d_hand, c_hand", [
    ("A", 0.0, 0.0, 0.5, 7.86, 0.94),
    ("B", 1.0, 0.25, 0.25, 8.91, 0.80),
])
def test_fabg_like(case, nadp, k2, h2, d_hand, c_hand):
    specs = _fabg()
    conc = {"oxo": 0.1, "NADPH": 6.1, "hyd": 0.05, "NADP": nadp, "oxo2": k2, "hyd2": h2}
    ref = _mechanism_fluxes(specs, conc)
    v, rhs, y = _code_rates(specs, conc, ON)
    assert abs(v[0] - ref[0]) / abs(ref[0]) < 0.01
    _, den = rhs.compute_v_and_den(y)
    assert den[0] == pytest.approx(d_hand, abs=0.005)
    assert rhs.competition_factor(y)[0] == pytest.approx(c_hand, abs=0.005)
    # Earlier patch form: acyl substrates only, c = (1 + k)/(1 + k + k2).
    v_off = _code_rates(specs, conc, OFF)[0]
    c_old = (1 + conc["oxo"]) / (1 + conc["oxo"] + conc["oxo2"])
    assert abs(v_off[0] * c_old - ref[0]) / abs(ref[0]) > 0.01


def test_holo_acp_partner_forms():
    """TesA-like: acyl_n -> FA_n + holo-ACP on two chain lengths (product-inhibited)."""
    specs = [
        _spec({"acyl1": 1.0}, {"FA1": 1.0, "ACP": 1.0}),
        _spec({"acyl2": 1.0}, {"FA2": 1.0, "ACP": 1.0}),
    ]
    conc = {"acyl1": 0.5, "acyl2": 1.0, "FA1": 2.0, "FA2": 3.0, "ACP": 1.5}
    ref = _mechanism_fluxes(specs, conc)
    v, rhs, y = _code_rates(specs, conc, ON)
    assert _rel_err(v, ref) < 0.01
    # "F = 1" variant: other chain's FA.ACP form ignored.
    v_off = _code_rates(specs, conc, OFF)[0]
    _, den = rhs.compute_v_and_den(y)
    d1 = den[0]
    c_f1 = d1 / (d1 + conc["acyl2"] + conc["FA2"])
    assert abs(v_off[0] * c_f1 - ref[0]) / ref[0] > 0.01


def test_mixed_group():
    """A+B -> P and A -> Q: the A-only reaction runs at its own MM times 1/(1+b)."""
    specs = [_spec({"A": 1.0, "B": 1.0}), _spec({"A": 1.0})]
    conc = {"A": 2.0, "B": 0.5}
    ref = _mechanism_fluxes(specs, conc)
    v, rhs, y = _code_rates(specs, conc, ON)
    assert _rel_err(v, ref) < 0.01
    assert rhs.competition_factor(y)[1] == pytest.approx(1.0 / 1.5, rel=1e-12)

