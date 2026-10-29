"""
Tests for reversible enzyme inhibition (bees/inhibition.py and its use in the simulator).

To run: pytest -v tests/test_inhibition.py
"""

from types import SimpleNamespace

import numpy as np
import pytest

from bees.core_edge_model import CoreEdgeModel, SpeciesData
from bees.inhibition import (
    InhibitionTable,
    InhibitorSpec,
    format_inhibition_terms,
    inhibited_mm_rate,
    inverse_site_constants,
    mm_scales,
    morrison_free_fraction,
    morrison_expression,
)
from bees.simulator import ODESimulator

VMAX, KM, S, KI = 0.1, 1.0, 3.0, 0.5


def _analytic(type_, i, alpha=None):
    """Textbook rate of single-substrate Michaelis-Menten with one reversible inhibitor."""
    if type_ == "competitive":
        return VMAX * S / (KM * (1 + i / KI) + S)
    if type_ == "uncompetitive":
        return VMAX * S / (KM + S * (1 + i / KI))
    if type_ == "noncompetitive":
        return VMAX * S / ((KM + S) * (1 + i / KI))
    return VMAX * S / (KM * (1 + i / KI) + S * (1 + i / (alpha * KI)))


@pytest.mark.parametrize("type_,alpha", [
    ("competitive", None), ("uncompetitive", None), ("noncompetitive", None), ("mixed", 3.0), ("mixed", 0.4),
])
@pytest.mark.parametrize("i", [0.0, 0.2, 2.0, 50.0])
def test_factors_match_analytic_michaelis_menten(type_, alpha, i):
    assert inhibited_mm_rate(VMAX, KM, S, i, type_, KI, alpha) == pytest.approx(_analytic(type_, i, alpha), rel=1e-12)


def test_competitive_scales_km_only_and_noncompetitive_vmax_only():
    assert mm_scales("competitive", KI, 1.0) == pytest.approx((1.0, 1 + 1.0 / KI))
    assert mm_scales("noncompetitive", KI, 1.0) == pytest.approx((1 / (1 + 1.0 / KI), 1.0))
    f_v, f_km = mm_scales("uncompetitive", KI, 1.0)
    assert f_v == pytest.approx(f_km) == pytest.approx(1 / (1 + 1.0 / KI))


def test_mixed_limits_are_competitive_and_noncompetitive():
    assert inhibited_mm_rate(VMAX, KM, S, 2.0, "mixed", KI, 1e9) == pytest.approx(
        inhibited_mm_rate(VMAX, KM, S, 2.0, "competitive", KI), rel=1e-6)
    assert inhibited_mm_rate(VMAX, KM, S, 2.0, "mixed", KI, 1.0) == pytest.approx(
        inhibited_mm_rate(VMAX, KM, S, 2.0, "noncompetitive", KI), rel=1e-12)


def test_invalid_constants_rejected():
    with pytest.raises(ValueError):
        inverse_site_constants("competitive", 0.0)
    with pytest.raises(ValueError):
        inverse_site_constants("mixed", 1.0)
    with pytest.raises(ValueError):
        inverse_site_constants("tight_binding", 1.0)


class TestMorrison:
    def test_reduces_to_classic_ki_when_enzyme_is_far_below_ki(self):
        i = np.array([0.0, 1e-3, 0.1, 1.0, 20.0])
        classic = 1.0 / (1.0 + i / 0.5)
        assert morrison_free_fraction(i, 1e-9, 0.5) == pytest.approx(classic, rel=1e-6)
        assert morrison_free_fraction(i, 0.0, 0.5) == pytest.approx(classic, rel=1e-12)

    def test_stoichiometric_titration_when_enzyme_is_far_above_ki(self):
        """E >> Ki: each inhibitor molecule shuts one enzyme, so free fraction = 1 - I/E until I = E."""
        e, ki = 1.0, 1e-6
        assert morrison_free_fraction(0.25, e, ki) == pytest.approx(0.75, rel=1e-4)
        assert morrison_free_fraction(0.9, e, ki) == pytest.approx(0.1, rel=1e-4)
        assert morrison_free_fraction(2.0, e, ki) < 1e-5

    def test_matches_closed_form_quadratic(self):
        e, i, ki = 0.3, 0.2, 0.15
        bound = ((e + i + ki) - np.sqrt((e + i + ki) ** 2 - 4 * e * i)) / 2
        assert morrison_free_fraction(i, e, ki) == pytest.approx(1 - bound / e, rel=1e-12)

    def test_accurate_deep_in_the_tight_regime(self):
        """No cancellation error for E << Ki (the methotrexate regime: E ~ 1e-6 mM, Ki = 1e-7 mM)."""
        e, ki, i = 1e-9, 1e-7, 1e-4
        assert morrison_free_fraction(i, e, ki) == pytest.approx(ki / (ki + i), rel=1e-4)

    def test_expression_evaluates_like_the_function(self):
        e, i, ki = 0.3, 0.2, 0.15
        text = morrison_expression("i", "e", "ki").replace("sqrt", "np.sqrt").replace("^", "**")
        assert eval(text, {"np": np}, {"i": i, "e": e, "ki": ki}) == pytest.approx(morrison_free_fraction(i, e, ki))


def _model(inhibitors, i_conc, *, km=KM, kcat=10.0, e_conc=0.01, species=("S",), km_per=None, extra=(), stoichiometry=None):
    """S -> P by E (kcat 10 / s so that Vmax = 0.1 mM/s) with a constant inhibitor species I."""
    model = CoreEdgeModel()
    for lab, c, enz, const in [("S", S, False, False), ("P", 0.0, False, False),
                               ("E", e_conc, True, True), ("I", i_conc, False, True), *extra]:
        model.add_core_species(SpeciesData(label=lab, concentration=c, initial_concentration=c,
                                           is_enzyme=enz, constant=const))
    kin = SimpleNamespace(km=km, kcat=kcat, vmax=None, km_per_substrate=km_per)
    rxn = SimpleNamespace(
        enzyme_label="E", substrate_label="S", reactant_labels=list(species), product_labels=["P"],
        stoichiometry=stoichiometry or {"S": -1, "P": 1}, rate_law="Michaelis-Menten", kinetics=kin,
        ec_number="EC 1.1.1.1", template=SimpleNamespace(reversible=False), thermo=None,
        inhibitors=inhibitors,
    )
    model.add_reaction(rxn)
    return model


def _rate(model):
    rhs = ODESimulator(model)._build_rhs()
    return rhs.compute_v(np.array(model.get_all_concentration_vector()))[0]


@pytest.mark.parametrize("type_,alpha", [
    ("competitive", None), ("uncompetitive", None), ("noncompetitive", None), ("mixed", 3.0),
])
def test_simulator_rate_matches_analytic(type_, alpha):
    model = _model([InhibitorSpec("i", type_, KI, alpha, "s")], 2.0)
    assert _rate(model) == pytest.approx(_analytic(type_, 2.0, alpha), rel=1e-12)


def test_substrate_defaults_to_the_first_substrate_with_a_km():
    model = _model([InhibitorSpec("i", "competitive", KI)], 2.0)
    assert _rate(model) == pytest.approx(_analytic("competitive", 2.0), rel=1e-12)


def test_no_inhibitor_present_leaves_the_rate_alone():
    plain = _rate(_model(None, 0.0))
    assert _rate(_model([InhibitorSpec("i", "competitive", KI, None, "s")], 0.0)) == pytest.approx(plain, rel=1e-14)


def test_tight_binding_in_the_simulator_uses_enzyme_concentration():
    e = 0.01
    model = _model([InhibitorSpec("i", "tight_binding", KI)], 2.0, e_conc=e)
    v0 = _rate(_model(None, 0.0, e_conc=e))
    assert _rate(model) / v0 == pytest.approx(morrison_free_fraction(2.0, e, KI), rel=1e-12)


def test_formatted_factors_put_inhibitors_of_one_site_in_one_sum():
    table_site = format_inhibition_terms(
        [InhibitorSpec("i", "competitive", 0.5, None, "s"), InhibitorSpec("j", "competitive", 2.0, None, "s")],
        lambda lab: lab, lambda kind, sp, value: f"K{sp.inhibitor}", "s",
    )
    v_factor, km_factor = table_site
    assert v_factor == "1"
    assert km_factor["s"] == "((1 + i / Ki + j / Kj) / (1))"


def test_two_competitive_inhibitors_of_one_substrate_are_mutually_exclusive():
    model = _model([InhibitorSpec("i", "competitive", 0.5, None, "s"),
                    InhibitorSpec("j", "competitive", 2.0, None, "s")], 1.0,
                   extra=[("J", 4.0, False, True)])
    expected = VMAX * S / (KM * (1 + 1.0 / 0.5 + 4.0 / 2.0) + S)
    assert _rate(model) == pytest.approx(expected, rel=1e-12)


def test_competitive_inhibition_with_a_second_substrate():
    """Only the substrate named by the inhibitor sees the larger Km."""
    model = _model(
        [InhibitorSpec("i", "competitive", KI, None, "s2")], 2.0,
        species=("S", "S2"), km_per={"S": 1.0, "S2": 0.5}, extra=[("S2", 2.0, False, True)],
        stoichiometry={"S": -1, "S2": -1, "P": 1},
    )
    expected = 0.1 * (S / (1.0 + S)) * (2.0 / (0.5 * (1 + 2.0 / KI) + 2.0))
    assert _rate(model) == pytest.approx(expected, rel=1e-12)


def test_unresolved_inhibitors_are_reported():
    model = _model([InhibitorSpec("missing", "competitive", KI), InhibitorSpec("i", "competitive", KI, None, "nope")], 1.0)
    unresolved = ODESimulator(model)._build_rhs().unresolved_inhibitors
    assert [u[1] for u in unresolved] == ["missing", "i"]
    assert "not a model species" in unresolved[0][2]
    assert "not a substrate" in unresolved[1][2]


def test_state_matrix_agrees_with_state_vector():
    """The trajectory path (n_sp x n_t) and the solver path (n_sp) must give the same rates."""
    model = _model([InhibitorSpec("i", "mixed", KI, 2.0, "s"), InhibitorSpec("i", "tight_binding", 0.3)], 1.0)
    rhs = ODESimulator(model)._build_rhs()
    y = np.array(model.get_all_concentration_vector())
    cols = np.stack([y, y * np.array([0.5, 1.0, 1.0, 3.0])], axis=1)
    mat = rhs.compute_v(cols)
    assert mat[0, 0] == pytest.approx(rhs.compute_v(cols[:, 0])[0], rel=1e-12)
    assert mat[0, 1] == pytest.approx(rhs.compute_v(cols[:, 1])[0], rel=1e-12)


def test_table_scales_shapes():
    spec = InhibitorSpec("i", "competitive", 1.0)
    table = InhibitionTable(3, 2, [(1, 0, 0, spec, 0.0)])
    v, km = table.scales(np.array([2.0, 1.0]))
    assert v.shape == (3,) and km.shape == (3, 2)
    assert km[1, 0] == pytest.approx(3.0) and km[0, 0] == 1.0 and v[1] == 1.0


def test_reversible_reaction_uses_the_inhibited_km():
    """A reversible (CM) law reads Km off the reaction; competitive inhibition must reach it."""
    from bees.flux_calculator import compute_reversible_mm_rate

    def make(i_conc, inhibitors):
        model = _model(inhibitors, i_conc, km=None, km_per={"S": KM, "P": 2.0})
        rxn = model.core_reactions[0]
        rxn.template = SimpleNamespace(reversible=True)
        rxn.thermo = SimpleNamespace(irreversible=False, keq=50.0, source="test")
        model.set_all_concentrations([S, 0.3, 0.01, i_conc])
        return model, rxn

    model, rxn = make(2.0, [InhibitorSpec("i", "competitive", KI, None, "s")])
    got = _rate(model)
    plain_model, plain_rxn = make(2.0, None)
    plain_rxn.kinetics.km_per_substrate = {"S": KM * (1 + 2.0 / KI), "P": 2.0}
    conc = {"s": S, "p": 0.3}
    assert got == pytest.approx(compute_reversible_mm_rate(plain_rxn, conc, {"e": 0.01}), rel=1e-10)
    assert got < _rate(make(0.0, None)[0])
