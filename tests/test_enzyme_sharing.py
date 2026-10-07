"""
Tests for shared-enzyme competition (bees/enzyme_sharing.py) and its plumbing.

To run:  pytest -v tests/test_enzyme_sharing.py
"""

from types import SimpleNamespace

import numpy as np
import pytest

import bees.enzyme_sharing as sharing
from bees.core_edge_model import CoreEdgeModel, RateLawOptions, SpeciesData
from bees.simulator import ODESimulator, _VectorizedRHS, rate_law_options_of
from bees.thermodynamics import ThermoData

ON = RateLawOptions(enzyme_competition=True)
E = 0.01
KCAT = 10.0
VMAX = KCAT * E


def _rxn(enzyme, subs, prods, km, kcat=KCAT, keq=None, nu=None):
    """Reaction stub; Km for a product makes the law product-inhibited (or reversible with keq)."""
    nu = nu or {}
    stoich = {s: -nu.get(s, 1) for s in subs}
    stoich.update({p: nu.get(p, 1) for p in prods})
    thermo = None
    if keq is not None:
        thermo = ThermoData(dgr_prime_kJmol=-1.0, sigma_kJmol=1.0, keq=keq, kcat_rev=None,
                            irreversible=False, source="equilibrator")
    return SimpleNamespace(
        enzyme_label=enzyme, reactant_labels=list(subs), product_labels=list(prods),
        stoichiometry=stoich, rate_law="Michaelis-Menten",
        kinetics=SimpleNamespace(km=None, kcat=kcat, vmax=None, km_per_substrate=dict(km)),
        template=SimpleNamespace(reversible=keq is not None), thermo=thermo,
        feedback_inhibitors=None, ec_number=None,
    )


def _rhs(reactions, species, options=ON, n_core=None, enz=None):
    enz = enz or {r.enzyme_label.lower(): E for r in reactions}
    return _VectorizedRHS(
        species_labels=species,
        reactions=reactions,
        alias_to_model_label={s.lower(): s.lower() for s in species},
        enzyme_conc_map=enz,
        constant_mask=np.zeros(len(species), dtype=bool),
        n_core_reactions=len(reactions) if n_core is None else n_core,
        options=options,
    )


def _y(species, conc):
    return np.array([conc.get(s, 0.0) for s in species], dtype=float)


SP2 = ["A", "B", "P", "Q"]


def _two_mm(enzyme2="E1", ka=1.0, kb=2.0):
    return [_rxn("E1", ["A"], ["P"], {"A": ka}), _rxn(enzyme2, ["B"], ["Q"], {"B": kb})]


# ---------------------------------------------------------------------------
# Partition-function competition
# ---------------------------------------------------------------------------

class TestCompetition:
    def test_single_reaction_is_unchanged(self):
        rxns = [_rxn("E1", ["A"], ["P"], {"A": 1.0})]
        y = _y(SP2, {"A": 2.0})
        v_on = _rhs(rxns, SP2).compute_v(y)
        v_off = _rhs(rxns, SP2, options=None).compute_v(y)
        assert v_on[0] == v_off[0] == pytest.approx(VMAX * 2.0 / 3.0)

    def test_different_enzymes_do_not_compete(self):
        y = _y(SP2, {"A": 2.0, "B": 3.0})
        v_on = _rhs(_two_mm(enzyme2="E2"), SP2).compute_v(y)
        v_off = _rhs(_two_mm(enzyme2="E2"), SP2, options=None).compute_v(y)
        np.testing.assert_array_equal(v_on, v_off)

    def test_two_substrates_competitive_mm(self):
        a, b = 2.0, 3.0 / 2.0
        v = _rhs(_two_mm(), SP2).compute_v(_y(SP2, {"A": 2.0, "B": 3.0}))
        np.testing.assert_allclose(v[:2], [VMAX * a / (1 + a + b), VMAX * b / (1 + a + b)], rtol=1e-12)

    def test_plain_mm_fallback_equals_old_patch_form(self):
        """Fallback MM (no product Km): c = (1+x)/(1+x+others), the earlier patch formula."""
        rhs = _rhs(_two_mm(), SP2)
        y = _y(SP2, {"A": 2.0, "B": 3.0})
        a, b = 2.0, 1.5
        np.testing.assert_allclose(rhs.competition_factor(y)[:2], [(1 + a) / (1 + a + b), (1 + b) / (1 + a + b)],
                                   rtol=1e-14)

    def test_edge_reaction_does_not_slow_core(self):
        a, b = 2.0, 1.5
        v = _rhs(_two_mm(), SP2, n_core=1).compute_v(_y(SP2, {"A": 2.0, "B": 3.0}))
        assert v[0] == pytest.approx(VMAX * a / (1 + a), rel=1e-12)
        assert v[1] == pytest.approx(VMAX * b / (1 + a + b), rel=1e-12)

    def test_matrix_matches_vector(self):
        rxns = [_rxn("E1", ["A"], ["P"], {"A": 1.0, "P": 0.5}), _rxn("E1", ["B"], ["Q"], {"B": 2.0, "Q": 0.3})]
        rhs = _rhs(rxns, SP2, options=RateLawOptions(enzyme_competition=True))
        ys = [_y(SP2, {"A": 2.0, "B": 3.0, "P": 0.1}), _y(SP2, {"A": 0.5, "B": 7.0, "Q": 1.0})]
        v_mat = rhs.compute_v(np.column_stack(ys))
        for k, y in enumerate(ys):
            np.testing.assert_allclose(v_mat[:, k], rhs.compute_v(y), rtol=1e-12)

    def test_off_is_identical(self):
        rxns = _two_mm()
        y = _y(SP2, {"A": 2.0, "B": 3.0})
        ref = _rhs(rxns, SP2, options=None).compute_v(y)
        np.testing.assert_array_equal(_rhs(rxns, SP2, options=RateLawOptions()).compute_v(y), ref)

    def test_forward_reverse_duplicates_count_once(self):
        """A<->C written forward and reverse: identical complex sets, so neither is slowed."""
        sp = ["A", "C"]
        rxns = [_rxn("E1", ["A"], ["C"], {"A": 1.0, "C": 2.0}), _rxn("E1", ["C"], ["A"], {"C": 2.0, "A": 1.0})]
        c = _rhs(rxns, sp).competition_factor(_y(sp, {"A": 3.0, "C": 5.0}))
        np.testing.assert_array_equal(c, [1.0, 1.0])

    def test_stoichiometry_two_uses_binomial_multiplicity(self):
        rxns = [_rxn("E1", ["A"], ["P"], {"A": 1.0}, nu={"A": 2}), _rxn("E1", ["B"], ["Q"], {"B": 1.0})]
        a, b = 0.7, 1.3
        c = _rhs(rxns, SP2).competition_factor(_y(SP2, {"A": a, "B": b}))
        assert c[1] == pytest.approx((1 + b) / ((1 + a) ** 2 + b), rel=1e-12)
        assert c[0] == pytest.approx((1 + a) ** 2 / ((1 + a) ** 2 + b), rel=1e-12)

    def test_mixed_group(self):
        """A+B->P and A->Q: the A-only reaction sees the B and A.B forms, c = 1/(1+b)."""
        rxns = [_rxn("E1", ["A", "B"], ["P"], {"A": 1.0, "B": 1.0}), _rxn("E1", ["A"], ["Q"], {"A": 1.0})]
        a, b = 2.0, 0.5
        c = _rhs(rxns, SP2).competition_factor(_y(SP2, {"A": a, "B": b}))
        assert c[0] == 1.0
        assert c[1] == pytest.approx(1.0 / (1.0 + b), rel=1e-12)

    def test_products_switch(self):
        rxns = [_rxn("E1", ["A"], ["P"], {"A": 1.0, "P": 1.0}), _rxn("E1", ["B"], ["Q"], {"B": 1.0, "Q": 1.0})]
        a, b, p, q = 2.0, 0.5, 0.25, 0.75
        y = _y(SP2, {"A": a, "B": b, "P": p, "Q": q})
        d1 = 1 + a + p
        c_on = _rhs(rxns, SP2).competition_factor(y)
        c_off = _rhs(rxns, SP2, options=RateLawOptions(enzyme_competition=True,
                                                       competition_products=False)).competition_factor(y)
        assert c_on[0] == pytest.approx(d1 / (d1 + b + q), rel=1e-12)
        assert c_off[0] == pytest.approx(d1 / (d1 + b), rel=1e-12)

    def test_canonical_km_is_order_independent_and_warns(self, monkeypatch):
        from unittest.mock import MagicMock

        sharing._WARNED.clear()
        log = MagicMock()
        monkeypatch.setattr(sharing, "logger", log)
        sp = ["A", "B", "C", "P", "Q", "R"]
        rxns = [
            _rxn("E1", ["A"], ["P"], {"A": 1.0}),
            _rxn("E1", ["A", "B"], ["Q"], {"A": 4.0, "B": 1.0}),
            _rxn("E1", ["C"], ["R"], {"C": 1.0}),
        ]
        y = _y(sp, {"A": 1.5, "B": 0.4, "C": 2.0})
        c = _rhs(rxns, sp).competition_factor(y)
        assert "spans" in log.warning.call_args[0][0]
        order = [2, 0, 1]
        c_shuffled = _rhs([rxns[i] for i in order], sp).competition_factor(y)
        np.testing.assert_allclose(c_shuffled, c[order], rtol=1e-14)
        # C-reaction sees A (K = geomean(1, 4) = 2), B and A.B with A's canonical K.
        a, bb = 1.5 / 2.0, 0.4
        assert c[2] == pytest.approx(3.0 / (3.0 + a + bb + a * bb), rel=1e-12)

    def test_factor_never_exceeds_one(self):
        rng = np.random.default_rng(0)
        sp = ["A", "B", "C", "P", "Q"]
        rxns = [
            _rxn("E1", ["A", "B"], ["P"], {"A": 0.3, "B": 2.0, "P": 1.0}, keq=50.0),
            _rxn("E1", ["C"], ["Q"], {"C": 0.1, "Q": 0.4}),
            _rxn("E1", ["A"], ["Q"], {"A": 0.5}),
        ]
        rhs = _rhs(rxns, sp)
        for _ in range(50):
            c = rhs.competition_factor(rng.uniform(0, 5, size=len(sp)))
            assert np.all(c <= 1.0 + 1e-15) and np.all(c > 0)

    def test_legacy_form_equals_patch_formula(self):
        sp = ["butanoyl-[ACP]", "hexanoyl-[ACP]", "P1", "P2"]
        rxns = [_rxn("FabI", [sp[0]], ["P1"], {sp[0]: 1.0}), _rxn("FabI", [sp[1]], ["P2"], {sp[1]: 1.0})]
        x1, x2 = 2.0, 3.0
        v = _rhs(rxns, sp, options=RateLawOptions(enzyme_competition=True, competition_form="legacy"),
                 enz={"fabi": E}).compute_v(_y(sp, {sp[0]: x1, sp[1]: x2}))
        np.testing.assert_allclose(v[:2], [VMAX * x1 / (1 + x1 + x2), VMAX * x2 / (1 + x1 + x2)], rtol=1e-14)

    def test_contributions_name_complexes(self):
        rxns = [_rxn("E1", ["A"], ["P"], {"A": 1.0, "P": 0.5}), _rxn("E1", ["B"], ["Q"], {"B": 2.0, "Q": 0.25})]
        rhs = _rhs(rxns, SP2)
        y = _y(SP2, {"A": 1.0, "B": 1.0, "P": 0.5, "Q": 0.5})
        assert rhs.competition_contributions(y) == {
            "e1": pytest.approx({"A": 1.0, "P": 1.0, "B": 0.5, "Q": 2.0})
        }
        br = rhs.competition_breakdown(y)
        assert set(br[0]["others"]) == {"B", "Q"} and br[0]["D"] == pytest.approx(3.0)


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def _model(options=None):
    model = CoreEdgeModel()
    for lab, c in (("A", 2.0), ("B", 3.0), ("P", 0.0), ("Q", 0.0)):
        model.add_core_species(SpeciesData(label=lab, concentration=c, initial_concentration=c))
    model.add_core_species(SpeciesData(label="E1", concentration=E, initial_concentration=E,
                                       is_enzyme=True, constant=True))
    model.core_reactions.extend(_two_mm())
    if options is not None:
        model.rate_law_options = options
    return model


class TestPlumbing:
    def test_simulator_reads_model_options(self):
        r_off = ODESimulator(_model()).simulate(end_time=10.0, time_step=5.0)
        r_on = ODESimulator(_model(ON)).simulate(end_time=10.0, time_step=5.0)
        i = r_on.species_labels.index("P")
        assert r_on.y[i, -1] < r_off.y[i, -1]

    def test_old_pickle_without_options_means_off(self):
        model = _model(ON)
        del model.rate_law_options
        assert rate_law_options_of(model) == RateLawOptions()

    def test_old_env_var_raises(self, monkeypatch):
        monkeypatch.setenv("BEES_ENZYME_COMPETITION", "0")
        with pytest.raises(RuntimeError, match="settings.enzyme_competition"):
            ODESimulator(_model())

    def test_unknown_form_rejected(self):
        with pytest.raises(ValueError):
            RateLawOptions(enzyme_competition=True, competition_form="other")

    def test_settings_default_off(self):
        from bees.schema import Settings
        s = Settings(end_time=10.0)
        assert s.enzyme_competition is True
        assert Settings(end_time=10.0, enzyme_competition=False).enzyme_competition is False

    def test_sbml_and_ode_text_match_simulator(self, tmp_path):
        libsbml = pytest.importorskip("libsbml")
        from unittest.mock import MagicMock
        from bees.exporter import EnlargerExporter

        model = _model(ON)
        logger = MagicMock()
        exporter = EnlargerExporter(
            model=model, profiles=[], output_directory=str(tmp_path), logger=logger,
            reaction_id_by_sig={}, reaction_first_seen_iter={}, reaction_core_enter_iter={},
            reaction_obj_by_sig={}, iteration_summaries=[],
        )
        path = exporter.export_sbml()
        assert path and not (tmp_path / "model.xml.SKIPPED.txt").exists()

        doc = libsbml.readSBML(path)
        m = doc.getModel()
        conc = {"a": 2.0, "b": 3.0, "p": 0.4, "q": 0.2}
        ns = {p.getId(): p.getValue() for p in m.getListOfParameters()}
        ns["compartment1"] = 1.0
        for s in m.getListOfSpecies():
            ns[s.getId()] = conc.get(s.getName().lower(), E)
        rhs = _rhs(model.core_reactions, [s.label for s in model.core_species], options=ON,
                   enz={"e1": E})
        y = _y([s.label for s in model.core_species], {k.upper(): v for k, v in conc.items()} | {"E1": E})
        v = rhs.compute_v(y)
        assert rhs.competition_factor(y)[0] < 1.0
        for i in range(m.getNumReactions()):
            formula = libsbml.formulaToL3String(m.getReaction(i).getKineticLaw().getMath())
            got = eval(formula.replace("^", "**"), {"__builtins__": {}}, ns)
            assert got == pytest.approx(v[i], rel=1e-9), formula
            assert m.getReaction(i).getNumModifiers() >= 1

        ode = ODESimulator(model, logger=logger).export_ode_equations(str(tmp_path / "ode.txt"), 1, 10.0)
        text = open(ode).read()
        assert "Kc(E1, A)" in text and "Kc(E1, B)" in text
