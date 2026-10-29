"""
Tests for time-dependent dosing (bees/dosing.py and its use in the simulator).

To run: pytest -v tests/test_dosing.py
"""

import math
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.integrate import quad

from bees.core_edge_model import CoreEdgeModel, SpeciesData
from bees.dosing import (
    SINGULAR_START,
    DosePlan,
    DoseSchedule,
    cumulative_fraction,
    normalize_params,
    rate_formula,
    schedules_from_spec,
)
from bees.simulator import ODESimulator

DOSE = 2.0
H = 3600.0

PROFILES = {
    "bolus": {},
    "constant": {"duration": 500.0},
    "first_order": {"half_life": 300.0},
    "biexponential": {"fraction_fast": 0.1, "half_life_fast": 2 * H, "half_life_slow": 45 * H},
    "weibull": {"tau": 24 * H, "beta": 0.6},
    "korsmeyer_peppas": {"k": 0.02, "n": 0.45},
    "table": {"table": [[0, 0], [100, 2], [300, 2], [400, 0.5], [900, 0]]},
}


def _model():
    """S -> P by E plus an inert drug species D that only enters by dosing."""
    model = CoreEdgeModel()
    for lab, c, enz, const in [("S", 10.0, False, False), ("P", 0.0, False, False),
                               ("E", 0.01, True, True), ("D", 0.0, False, False)]:
        model.add_core_species(SpeciesData(label=lab, concentration=c, initial_concentration=c,
                                           is_enzyme=enz, constant=const))
    kin = SimpleNamespace(km=1.0, kcat=10.0, vmax=None, km_per_substrate=None)
    model.add_reaction(SimpleNamespace(
        enzyme_label="E", substrate_label="S", reactant_labels=["S"], product_labels=["P"],
        stoichiometry={"S": -1, "P": 1}, rate_law="Michaelis-Menten", kinetics=kin,
        ec_number="EC 1.1.1.1", template=SimpleNamespace(reversible=False), thermo=None,
    ))
    return model


def _run(schedules, end_time=2000.0, stepwise=False, **kwargs):
    sim = ODESimulator(_model(), dosing=schedules)
    extra = {"interrupt_simulation_tol": 1e9} if stepwise else {"time_step": 5.0}
    res = sim.simulate(end_time=end_time, rtol=1e-10, atol=1e-13, **extra, **kwargs)
    assert res.success, res.message
    return res, res.y[res.species_labels.index("D")]


@pytest.mark.parametrize("profile", sorted(PROFILES))
def test_profile_releases_exactly_the_dose(profile):
    sch = DoseSchedule("D", profile, DOSE, PROFILES[profile])
    assert float(sch.released(1e12)) == pytest.approx(DOSE, rel=1e-6)
    # rate integral + impulses = dose (an infinite tail is cut after 40 slow half-lives)
    pts = sorted(set(sch.breakpoints()) | set(np.geomspace(1e-2, 40 * 45 * H, 60)))
    area = sum(quad(lambda t: float(sch.rate(t)), a, b, limit=400)[0] for a, b in zip(pts[:-1], pts[1:]))
    assert area + sum(a for _, a in sch.impulses()) == pytest.approx(DOSE, rel=2e-3)


@pytest.mark.parametrize("profile", sorted(set(PROFILES) - {"bolus"}))
def test_rate_is_the_derivative_of_the_release(profile):
    sch = DoseSchedule("D", profile, DOSE, PROFILES[profile])
    t = np.array([1.0, 50.0, 250.0, 350.0, 700.0])
    h = 1e-3
    numeric = (sch.released(t + h) - sch.released(t - h)) / (2 * h)
    assert sch.rate(t) == pytest.approx(numeric, rel=1e-4)


def test_profile_shapes():
    t = np.array([0.0, 100.0, 1000.0])
    assert cumulative_fraction("first_order", {"k": 0.01}, t) == pytest.approx(1 - np.exp(-0.01 * t))
    bi = {"fraction_fast": 0.1, "k_fast": 1e-3, "k_slow": 1e-5}
    assert cumulative_fraction("biexponential", bi, t) == pytest.approx(
        1 - (0.1 * np.exp(-1e-3 * t) + 0.9 * np.exp(-1e-5 * t)))
    assert cumulative_fraction("weibull", {"tau": 200.0, "beta": 0.6}, t) == pytest.approx(
        1 - np.exp(-(t / 200.0) ** 0.6))
    kp = {"k": 0.05, "n": 0.5}
    assert cumulative_fraction("korsmeyer_peppas", kp, np.array([100.0, 1e4])) == pytest.approx([0.5, 1.0])
    assert cumulative_fraction("constant", {"duration": 400.0}, np.array([100.0, 800.0])) == pytest.approx([0.25, 1.0])


def test_half_life_and_rate_constant_are_equivalent():
    a = normalize_params("first_order", {"half_life": 100.0})
    b = normalize_params("first_order", {"k": math.log(2) / 100.0})
    assert a["k"] == pytest.approx(b["k"])


@pytest.mark.parametrize("profile,params", [
    ("nope", {}), ("constant", {}), ("first_order", {}), ("first_order", {"k": -1}),
    ("biexponential", {"fraction_fast": 1.5, "k_fast": 1, "k_slow": 1}), ("weibull", {"tau": 1}),
    ("korsmeyer_peppas", {"k": 1}), ("table", {"table": [[0, 1]]}),
    ("table", {"table": [[0, 1], [0, 2]]}), ("table", {"table": [[0, 0], [1, 0]]}),
])
def test_invalid_parameters_are_rejected(profile, params):
    with pytest.raises(ValueError):
        normalize_params(profile, params)


def test_time_scale_stretches_the_release_only():
    for profile in set(PROFILES) - {"bolus"}:
        base = DoseSchedule("D", profile, DOSE, PROFILES[profile], t_start=10.0, elimination_k=1e-3)
        slow = base.with_time_scale(10.0)
        t = np.array([20.0, 500.0, 3000.0])
        # F_slow(10 + 10 tau) = F_base(10 + tau)
        assert slow.released(10.0 + 10.0 * (t - 10.0)) == pytest.approx(base.released(t), rel=1e-9)
        assert slow.elimination_k == base.elimination_k and slow.t_start == base.t_start


@pytest.mark.parametrize("profile", sorted(set(PROFILES) - {"bolus"}))
def test_rate_formula_matches_rate(profile):
    sch = DoseSchedule("D", profile, DOSE, PROFILES[profile], t_start=7.0)
    text = (rate_formula(sch).replace("piecewise", "_pw").replace("^", "**")
            .replace("&&", " and ").replace("exp", "_exp"))

    def _pw(*a):
        for val, cond in zip(a[0::2], a[1::2]):
            if cond:
                return val
        return a[-1]

    for t in [8.0, 120.0, 330.0, 450.0, 5e4]:
        got = eval(text, {"_pw": _pw, "_exp": math.exp, "time": t})
        assert got == pytest.approx(float(sch.rate(t)), rel=1e-9, abs=1e-18)


def test_first_order_input_with_elimination_matches_closed_form():
    k, k_el = 0.01, 0.002
    sch = DoseSchedule("D", "first_order", DOSE, {"k": k}, elimination_k=k_el)
    exact = lambda t: DOSE * k / (k - k_el) * (np.exp(-k_el * t) - np.exp(-k * t))
    for stepwise in (False, True):
        res, d = _run([sch], stepwise=stepwise)
        assert d == pytest.approx(exact(res.t), abs=1e-7)


def test_bolus_with_elimination_decays_exponentially():
    k_el = 0.003
    sch = DoseSchedule("D", "bolus", DOSE, elimination_k=k_el)
    for stepwise in (False, True):
        res, d = _run([sch], stepwise=stepwise)
        assert d[0] == pytest.approx(DOSE)
        assert d == pytest.approx(DOSE * np.exp(-k_el * res.t), rel=1e-6)


def test_late_bolus_is_applied_once_at_its_time():
    sch = DoseSchedule("D", "bolus", DOSE, t_start=500.0)
    res, d = _run([sch])
    assert d[res.t < 500.0] == pytest.approx(0.0, abs=1e-14)
    assert d[res.t >= 500.0] == pytest.approx(DOSE)
    assert np.all(np.diff(res.t) > 0)
    res, d = _run([sch], stepwise=True)
    assert d[res.t < 500.0].max() == 0.0 and d[-1] == pytest.approx(DOSE)


@pytest.mark.parametrize("stepwise", [False, True])
@pytest.mark.parametrize("profile", ["weibull", "korsmeyer_peppas", "constant", "table", "biexponential"])
def test_simulated_amount_is_the_released_amount(profile, stepwise):
    sch = DoseSchedule("D", profile, DOSE, PROFILES[profile], t_start=30.0)
    res, d = _run([sch], stepwise=stepwise)
    assert d[-1] == pytest.approx(float(sch.released(res.t[-1])), rel=2e-6, abs=1e-9)
    if not stepwise:
        # the grid point at t_start already holds the impulse that stands for the first millisecond
        after = res.t != 30.0
        assert d[after] == pytest.approx(sch.released(res.t)[after], rel=1e-4, abs=2e-6)


def test_two_entries_for_one_species_add():
    a = DoseSchedule("D", "bolus", 1.0)
    b = DoseSchedule("D", "first_order", 0.5, {"k": 0.01})
    res, d = _run([a, b])
    assert d[-1] == pytest.approx(1.0 + 0.5 * (1 - math.exp(-0.01 * res.t[-1])), rel=1e-6)


def test_dosing_does_not_disturb_the_chemistry():
    plain = ODESimulator(_model()).simulate(end_time=300.0, time_step=10.0, rtol=1e-9, atol=1e-12)
    res, _ = _run([DoseSchedule("D", "weibull", DOSE, PROFILES["weibull"])], end_time=300.0)
    i = res.species_labels.index("P")
    assert res.y[i, -1] == pytest.approx(plain.y[plain.species_labels.index("P"), -1], rel=1e-8)


def test_resumed_run_continues_the_schedule_from_t_offset():
    """The enlarger resumes an interrupted run from the model's concentrations at t_resume."""
    sch = DoseSchedule("D", "first_order", DOSE, {"k": 0.01})
    t_resume = 150.0
    model = _model()
    first = ODESimulator(model, dosing=[sch]).simulate(end_time=t_resume, time_step=5.0, rtol=1e-10, atol=1e-13)
    assert first.success
    rest = ODESimulator(model, dosing=[sch]).simulate(
        end_time=300.0 - t_resume, time_step=5.0, rtol=1e-10, atol=1e-13, t_offset=t_resume)
    d = rest.y[rest.species_labels.index("D")]
    assert d[-1] == pytest.approx(float(sch.released(300.0)), rel=1e-7)

    b = DoseSchedule("D", "bolus", DOSE)
    model = _model()
    ODESimulator(model, dosing=[b]).simulate(end_time=t_resume, time_step=5.0)
    rest = ODESimulator(model, dosing=[b]).simulate(end_time=100.0, time_step=5.0, t_offset=t_resume)
    assert rest.y[rest.species_labels.index("D")][-1] == pytest.approx(DOSE)  # not added twice


def test_plan_collects_start_amount_events_and_breakpoints():
    labels = {"d": 0, "x": 1}
    sch = [DoseSchedule("D", "bolus", 1.0), DoseSchedule("D", "bolus", 2.0, t_start=50.0),
           DoseSchedule("D", "constant", 3.0, {"duration": 100.0}, t_start=10.0),
           DoseSchedule("missing", "bolus", 1.0)]
    plan = DosePlan(sch, labels, 2)
    assert plan.y0_add == pytest.approx([1.0, 0.0])
    assert plan.events == [(50.0, 0, 2.0)]
    assert plan.breakpoints == [10.0, 50.0, 110.0]
    assert plan.missing == ["missing"]
    assert plan.next_break(0.0, 1000.0) == 10.0 and plan.next_break(110.0, 1000.0) == 1000.0
    late = DosePlan(sch, labels, 2, t_offset=60.0)
    assert late.y0_add.sum() == 0.0 and late.events == [] and late.breakpoints == [50.0]


def test_singular_profiles_deliver_their_first_instant_as_an_impulse():
    sch = DoseSchedule("D", "weibull", DOSE, {"tau": 1000.0, "beta": 0.5})
    (t_imp, amount), = sch.impulses()
    assert t_imp == 0.0 and amount == pytest.approx(DOSE * (1 - math.exp(-math.sqrt(SINGULAR_START / 1000.0))))
    assert float(sch.rate(SINGULAR_START / 2)) == 0.0
    smooth = DoseSchedule("D", "weibull", DOSE, {"tau": 1000.0, "beta": 1.5})
    assert smooth.impulses() == []


def test_schedules_from_spec_resolves_elimination():
    entries = [SimpleNamespace(species="D", profile="first_order", dose=1.0, params={"k": 0.1},
                               t_start=5.0, elimination_k=None, elimination_half_life=100.0)]
    (sch,) = schedules_from_spec(entries)
    assert sch.elimination_k == pytest.approx(math.log(2) / 100.0) and sch.t_start == 5.0
