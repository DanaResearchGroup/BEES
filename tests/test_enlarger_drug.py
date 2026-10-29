"""
Tests for the enlarger's inhibition, measured-kinetics and dosing wiring.

To run: pytest -v tests/test_enlarger_drug.py
"""

import tempfile
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from bees.dosing import DoseSchedule
from bees.enlarger import IterativeEnlarger
from bees.inhibition import InhibitorSpec
from bees.reaction_generator import apply_measured_kinetics
from bees.simulator import ODESimulator, SimulationResult


def _reaction(enzyme_label="Enzyme", reactants=("S",), products=("P",), km=1.0, kcat=10.0, km_per=None):
    kin = SimpleNamespace(km=km, kcat=kcat, vmax=None, km_per_substrate=km_per, source="database", kcat_sd=None)
    return SimpleNamespace(
        enzyme_label=enzyme_label, substrate_label=reactants[0], reactant_labels=list(reactants),
        product_labels=list(products),
        stoichiometry={**{r: -1 for r in reactants}, **{p: 1 for p in products}},
        rate_law="Michaelis-Menten", kinetics=kin, ec_number="EC 1.1.1.1",
        template=SimpleNamespace(reversible=False), thermo=None,
    )


def _inh(inhibitor="Drug", type="competitive", ki=0.1, alpha=None, substrate=None):
    return SimpleNamespace(inhibitor=inhibitor, type=type, ki=ki, alpha=alpha, substrate=substrate, source=None)


def _bees(inhibition=None, dosing=None, measured=None, end_time=200.0):
    settings = SimpleNamespace(
        end_time=end_time, time_step=10.0, toleranceMoveToCore=1e-5, toleranceKeepInEdge=0,
        max_iterations=3, max_edge_species=None, termination_conversion=None,
        termination_rate_ratio=None, save_simulation_profiles=False, saveEdgeSpecies=True,
    )
    species = [
        SimpleNamespace(label="S", concentration=10.0, reactive=True, solvent=False, constant=False, smiles=None),
        SimpleNamespace(label="Drug", concentration=0.0, reactive=False, solvent=False, constant=False,
                        smiles="CCN"),
    ]
    enzymes = [SimpleNamespace(label="Enzyme", concentration=0.01, reactive=True, solvent=False,
                               ecnumber="EC 1.1.1.1", constant=True, amino_acid_sequence="ACDEFG",
                               inhibition=inhibition, measured_kinetics=measured)]
    return SimpleNamespace(
        settings=settings, species=species, enzymes=enzymes, dosing=dosing, project="T",
        database=SimpleNamespace(name="db"), environment=SimpleNamespace(temperature=298.15, pH=7.0),
    )


@pytest.fixture
def output_dir():
    with tempfile.TemporaryDirectory() as d:
        yield d


def _enlarger(bees, output_dir, rxn=None, logger=None):
    generator = MagicMock()
    generator._generate_reactions.return_value = [rxn or _reaction()]
    return IterativeEnlarger(bees_object=bees, reaction_generator=generator,
                             logger=logger or MagicMock(), output_directory=output_dir)


class TestApplyInhibition:
    def test_attached_by_enzyme_label_with_lowercase_names(self, output_dir):
        bees = _bees(inhibition=[_inh("Drug", "mixed", 0.2, alpha=3.0, substrate="S")])
        enlarger = _enlarger(bees, output_dir)
        mine, other = _reaction("Enzyme"), _reaction("Other")
        enlarger._apply_inhibition([mine, other])
        assert mine.inhibitors == [InhibitorSpec("drug", "mixed", 0.2, 3.0, "s")]
        assert other.inhibitors is None

    def test_reapplied_each_time_so_it_never_leaks(self, output_dir):
        bees = _bees(inhibition=[_inh()])
        enlarger = _enlarger(bees, output_dir)
        rxn = _reaction()
        enlarger._apply_inhibition([rxn])
        assert rxn.inhibitors
        bees.enzymes[0].inhibition = None
        enlarger._apply_inhibition([rxn])
        assert rxn.inhibitors is None

    def test_omitted_ki_is_estimated_once_per_enzyme_inhibitor_pair(self, output_dir):
        bees = _bees(inhibition=[_inh(ki=None)])
        enlarger = _enlarger(bees, output_dir)
        estimator = MagicMock()
        estimator.estimate_ki.return_value = SimpleNamespace(ki=0.037, ki_sd=None, source="catpred")
        enlarger.reaction_generator.kinetics_estimator = estimator
        rxns = [_reaction(), _reaction()]
        enlarger._apply_inhibition(rxns)
        enlarger._apply_inhibition(rxns)
        estimator.estimate_ki.assert_called_once_with(enzyme_sequence="ACDEFG", inhibitor_smiles="CCN")
        assert rxns[0].inhibitors[0].ki == 0.037

    def test_given_ki_never_calls_the_estimator(self, output_dir):
        enlarger = _enlarger(_bees(inhibition=[_inh(ki=0.5)]), output_dir)
        enlarger.reaction_generator.kinetics_estimator = MagicMock()
        rxn = _reaction()
        enlarger._apply_inhibition([rxn])
        enlarger.reaction_generator.kinetics_estimator.estimate_ki.assert_not_called()
        assert rxn.inhibitors[0].ki == 0.5

    @pytest.mark.parametrize("failure", ["no_estimator", "raises", "no_value"])
    def test_unavailable_estimate_leaves_the_inhibitor_out_with_a_warning(self, output_dir, failure):
        logger = MagicMock()
        enlarger = _enlarger(_bees(inhibition=[_inh(ki=None)]), output_dir, logger=logger)
        estimator = MagicMock()
        estimator.estimate_ki.side_effect = FileNotFoundError("no CatPred") if failure == "raises" else None
        if failure == "no_value":
            estimator.estimate_ki.return_value = SimpleNamespace(ki=None, ki_sd=None, source="catpred(err)")
        enlarger.reaction_generator.kinetics_estimator = None if failure == "no_estimator" else estimator
        rxn = _reaction()
        enlarger._apply_inhibition([rxn])
        assert rxn.inhibitors is None
        assert "left out" in logger.warning.call_args.args[0]


class TestMeasuredKinetics:
    def test_overrides_and_default_km_cover_substrates_without_a_km(self):
        rxn = _reaction(reactants=("A", "ACP-thing"), km=None, kcat=None, km_per={"A": 7.0})
        measured = SimpleNamespace(kcat=30.0, km={"a": 0.5}, default_km=0.02, source="assay")
        apply_measured_kinetics(rxn, measured)
        assert rxn.kinetics.kcat == 30.0
        assert rxn.kinetics.km_per_substrate == {"A": 0.5, "ACP-thing": 0.02}
        assert rxn.kinetics.source == "assay"

    def test_buffered_cofactors_get_no_default_km_and_listed_products_are_kept(self):
        rxn = _reaction(reactants=("A", "H2O"), products=("B",), km=None, kcat=None)
        apply_measured_kinetics(rxn, SimpleNamespace(kcat=None, km={"B": 3.0}, default_km=0.1, source=None))
        assert rxn.kinetics.km_per_substrate == {"B": 3.0, "A": 0.1}

    def test_single_estimated_km_is_kept_for_unlisted_substrates(self):
        rxn = _reaction(reactants=("A", "B"), km=0.4, kcat=1.0)
        apply_measured_kinetics(rxn, SimpleNamespace(kcat=None, km={"A": 9.0}, default_km=None, source=None))
        assert rxn.kinetics.km_per_substrate == {"A": 9.0, "B": 0.4}

    def test_a_reaction_without_a_rate_law_gets_one(self):
        rxn = _reaction(kcat=None, km=None)
        rxn.rate_law = None
        apply_measured_kinetics(rxn, SimpleNamespace(kcat=5.0, km=None, default_km=None, source=None))
        assert rxn.rate_law == "Michaelis-Menten"

    def test_applying_twice_changes_nothing(self):
        rxn = _reaction(reactants=("A", "B"), km=None, kcat=None)
        measured = SimpleNamespace(kcat=5.0, km={"A": 1.0}, default_km=0.2, source=None)
        apply_measured_kinetics(rxn, measured)
        once = (rxn.kinetics.kcat, dict(rxn.kinetics.km_per_substrate))
        apply_measured_kinetics(rxn, measured)
        assert (rxn.kinetics.kcat, rxn.kinetics.km_per_substrate) == once

    def test_reapplied_after_the_rules_layer(self, output_dir):
        """A calibration that rescales kcat must not override the measured value."""
        measured = SimpleNamespace(kcat=42.0, km={"S": 0.3}, default_km=None, source=None)
        enlarger = _enlarger(_bees(measured=measured), output_dir)
        rxn = _reaction(kcat=1.0, km=1.0)

        class Rules:
            def configure_calibrations(self, names): pass
            def by_name(self, name): return SimpleNamespace(params={})
            def apply_all(self, reactions):
                for r in reactions:
                    r.kinetics.kcat *= 1000.0
                    r.kinetics.km_per_substrate = {"S": 999.0}

        enlarger.model.core_reactions.append(rxn)
        with patch("bees.rules.load_rules", return_value=Rules()):
            enlarger._attach_thermo_to_reactions()
        assert rxn.kinetics.kcat == 42.0
        assert rxn.kinetics.km_per_substrate["S"] == 0.3


def test_a_substrate_without_a_km_makes_mass_and_default_km_stops_it():
    """B has no Km, so the law treats it as saturated and the reaction keeps running after B is gone."""
    from bees.core_edge_model import CoreEdgeModel, SpeciesData

    def total_after(default_km):
        model = CoreEdgeModel()
        for lab, c, enz in [("A", 10.0, False), ("B", 0.5, False), ("C", 0.0, False), ("E", 0.01, True)]:
            model.add_core_species(SpeciesData(lab, c, c, enz, enz))
        rxn = _reaction("E", reactants=("A", "B"), products=("C",), km=None, kcat=10.0, km_per={"A": 1.0})
        if default_km:
            apply_measured_kinetics(rxn, SimpleNamespace(kcat=None, km=None, default_km=default_km, source=None))
        model.add_reaction(rxn)
        res = ODESimulator(model).simulate(end_time=200.0, time_step=10.0, rtol=1e-9, atol=1e-12)
        idx = {lab: i for i, lab in enumerate(res.species_labels)}
        return res.y[idx["B"], -1] + res.y[idx["C"], -1]

    assert total_after(None) > 5 * 0.5  # C grows far beyond the B that existed
    assert total_after(0.02) == pytest.approx(0.5, rel=1e-4)


class TestDosingThroughTheEnlarger:
    def test_drug_run_inhibits_and_delivers_the_dose(self, output_dir):
        dose = SimpleNamespace(species="Drug", profile="first_order", dose=0.5, params={"half_life": 20.0},
                               t_start=0.0, elimination_k=None, elimination_half_life=None)
        control = _enlarger(_bees(), output_dir)
        control.run()
        treated = _enlarger(_bees(inhibition=[_inh(ki=0.05)], dosing=[dose]), output_dir)
        treated.run()

        def final(enlarger, label):
            prof = enlarger._profiles[-1]
            return prof.y[prof.species_labels.index(label), -1]

        assert final(treated, "Drug") == pytest.approx(0.5, rel=1e-3)
        assert final(treated, "p") < 0.9 * final(control, "p")

    def test_schedules_follow_the_input_and_reach_every_simulation(self, output_dir):
        dose = SimpleNamespace(species="Drug", profile="bolus", dose=0.5, params=None, t_start=0.0,
                               elimination_k=None, elimination_half_life=30.0)
        enlarger = _enlarger(_bees(dosing=[dose]), output_dir)
        (sch,) = enlarger.dose_schedules
        assert isinstance(sch, DoseSchedule) and sch.elimination_k == pytest.approx(np.log(2) / 30.0)
        with patch.object(ODESimulator, "simulate", autospec=True, side_effect=ODESimulator.simulate) as sim:
            enlarger.run()
        assert sim.call_args_list
        assert all(call.args[0].dosing == [sch] for call in sim.call_args_list)

    def test_resumed_run_passes_the_time_offset(self, output_dir):
        enlarger = _enlarger(_bees(), output_dir)
        flagged = SimulationResult(
            t=np.array([0.0, 40.0]), y=np.zeros((1, 2)), species_labels=["x"], success=True,
            simulation_interrupted=True, interrupt_char_rate=1.0, interrupt_edge_rates={"zz": 0.0},
        )
        finished = SimulationResult(t=np.array([0.0, 160.0]), y=np.zeros((1, 2)), species_labels=["x"], success=True)
        with patch("bees.enlarger.ODESimulator.simulate", side_effect=[flagged, finished, finished, finished]) as sim:
            enlarger.run()
        offsets = [c.kwargs.get("t_offset") for c in sim.call_args_list]
        assert offsets[0] is None and offsets[1] == pytest.approx(40.0)

    def test_dose_of_an_unknown_species_is_reported(self, output_dir):
        dose = SimpleNamespace(species="Ghost", profile="bolus", dose=0.5, params=None, t_start=0.0,
                               elimination_k=None, elimination_half_life=None)
        logger = MagicMock()
        enlarger = _enlarger(_bees(dosing=[dose]), output_dir, logger=logger)
        enlarger.run()
        assert any("Ghost" in c.args[0] for c in logger.warning.call_args_list)

    def test_unresolved_inhibitor_is_reported(self, output_dir):
        logger = MagicMock()
        enlarger = _enlarger(_bees(inhibition=[_inh("Phantom")]), output_dir, logger=logger)
        enlarger.run()
        assert any("phantom" in c.args[0] for c in logger.warning.call_args_list)
