"""
Tests for the IterativeEnlarger module.

To run:  pytest -v tests/test_enlarger.py
"""

import pytest
import tempfile
import numpy as np
from unittest.mock import MagicMock, patch
from types import SimpleNamespace

from bees.enlarger import IterativeEnlarger, EnlargerResult
from bees.simulator import SimulationResult


# ---------------------------------------------------------------------------
# Helper to build a mock reaction
# ---------------------------------------------------------------------------

def _make_reaction(
    enzyme_label="Enzyme",
    substrate_label="S",
    reactant_labels=None,
    product_labels=None,
    stoichiometry=None,
    rate_law="Michaelis-Menten",
    km=1.0,
    kcat=10.0,
):
    kinetics = MagicMock()
    kinetics.km = km
    kinetics.kcat = kcat
    kinetics.vmax = None
    kinetics.km_per_substrate = None

    rxn = MagicMock()
    rxn.enzyme_label = enzyme_label
    rxn.substrate_label = substrate_label
    rxn.reactant_labels = reactant_labels or ["S"]
    rxn.product_labels = product_labels or ["P"]
    rxn.stoichiometry = stoichiometry or {"S": -1, "P": 1}
    rxn.rate_law = rate_law
    rxn.kinetics = kinetics
    # GeneratedReaction-like attributes
    rxn.ec_number = "EC 1.1.1.1"
    rxn.template = MagicMock()
    return rxn


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _interrupt_result(interrupted: bool, edge_label: str = "p") -> SimulationResult:
    """Minimal SimulationResult: interrupted with flux into ``edge_label``, or a full pass."""
    return SimulationResult(
        t=np.array([0.0, 1.0]),
        y=np.zeros((1, 2)),
        species_labels=["x"],
        success=True,
        max_char_rate=1.0,
        final_char_rate=1.0,
        simulation_interrupted=interrupted,
        interrupt_char_rate=1.0 if interrupted else 0.0,
        interrupt_edge_rates={edge_label: 1.0} if interrupted else {},
    )


@pytest.fixture
def mock_bees_object():
    """Minimal bees_object with settings for iterative mode."""
    settings = SimpleNamespace(
        end_time=100.0,
        time_step=10.0,
        toleranceMoveToCore=1e-5,
        toleranceKeepInEdge=0,
        max_iterations=3,
        max_edge_species=None,
        termination_conversion=None,
        termination_rate_ratio=None,
        save_simulation_profiles=False,
        saveEdgeSpecies=True,
    )

    species = [
        SimpleNamespace(label="S", concentration=10.0, reactive=True, solvent=False, constant=False),
    ]
    enzymes = [
        SimpleNamespace(label="Enzyme", concentration=0.01, reactive=True, solvent=False,
                        ecnumber="EC 1.1.1.1", constant=True),
    ]

    bees_obj = SimpleNamespace(
        settings=settings,
        species=species,
        enzymes=enzymes,
        project="TestProject",
        database=SimpleNamespace(name="db"),
        environment=SimpleNamespace(temperature=298.15, pH=7.0),
    )
    return bees_obj


@pytest.fixture
def mock_reaction_generator():
    """ReactionGenerator whose per-enzyme discovery (_generate_reactions) returns S -> P."""
    mg = MagicMock()
    mg._generate_reactions.return_value = [
        _make_reaction(
            enzyme_label="Enzyme",
            substrate_label="S",
            reactant_labels=["S"],
            product_labels=["P"],
            stoichiometry={"S": -1, "P": 1},
            km=1.0,
            kcat=10.0,
        ),
    ]
    return mg


@pytest.fixture
def output_dir():
    with tempfile.TemporaryDirectory() as d:
        yield d


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestEnlargerResult:
    def test_default_values(self):
        r = EnlargerResult()
        assert r.iterations == 0
        assert not r.converged
        assert r.simulation_profiles == []


class TestIterativeEnlarger:
    def test_rate_ratio_termination_rmg_style(self, mock_bees_object, mock_reaction_generator, output_dir):
        """termination_rate_ratio compares final_char_rate / max_char_rate."""
        mock_bees_object.settings.termination_rate_ratio = 0.5
        logger = MagicMock()
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mock_reaction_generator,
            logger=logger,
            output_directory=output_dir,
        )
        sr = SimulationResult(
            t=np.array([0.0, 1.0]),
            y=np.zeros((1, 2)),
            species_labels=["x"],
            success=True,
            max_char_rate=10.0,
            final_char_rate=2.0,
        )
        assert enlarger._check_rate_ratio_termination(sr) is True

        sr2 = SimulationResult(
            t=np.array([0.0, 1.0]),
            y=np.zeros((1, 2)),
            species_labels=["x"],
            success=True,
            max_char_rate=10.0,
            final_char_rate=6.0,
        )
        assert enlarger._check_rate_ratio_termination(sr2) is False

    def test_initialise_model(self, mock_bees_object, mock_reaction_generator, output_dir):
        logger = MagicMock()
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mock_reaction_generator,
            logger=logger,
            output_directory=output_dir,
        )
        enlarger._initialise_model()
        # Should have the substrate and the enzyme in core
        assert enlarger.model.is_core_species("S")
        assert enlarger.model.is_core_species("Enzyme")
        assert len(enlarger.model.core_species) == 2

    def test_ingest_reactions(self, mock_bees_object, mock_reaction_generator, output_dir):
        logger = MagicMock()
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mock_reaction_generator,
            logger=logger,
            output_directory=output_dir,
        )
        enlarger._initialise_model()

        rxn = _make_reaction()
        enlarger._ingest_reactions([rxn])

        # P should have been added to edge (not in core initially)
        assert enlarger.model.is_edge_species("P")
        # S should still be in core
        assert enlarger.model.is_core_species("S")

    def test_run_converges(self, mock_bees_object, mock_reaction_generator, output_dir):
        """The enlarger should run and converge (no new significant species)."""
        logger = MagicMock()
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mock_reaction_generator,
            logger=logger,
            output_directory=output_dir,
        )
        result = enlarger.run()
        assert isinstance(result, EnlargerResult)
        assert result.converged
        assert result.iterations >= 1
        assert result.final_core_species >= 1
        assert "s" in enlarger.model.explored_species_lc

    def test_run_forgets_rule_values_from_a_previous_model(
        self, mock_bees_object, mock_reaction_generator, output_dir
    ):
        """A same-signature reaction from an earlier model must not overwrite this run's kinetics."""
        from bees.reaction_generator import reaction_signature
        from bees.rules import load_rules

        rules = load_rules()
        rxn = mock_reaction_generator._generate_reactions.return_value[0]
        rules._baselines[reaction_signature(rxn)] = {
            "kinetics": {"kcat": 999.0, "km": 1.0, "km_per_substrate": None},
        }
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mock_reaction_generator,
            logger=MagicMock(),
            output_directory=output_dir,
        )
        try:
            enlarger.run()
        finally:
            rules.reset_baselines()
        assert rxn.kinetics.kcat == 10.0

    def test_flux_summary_has_final_iteration_when_run_stops_early(
        self, mock_bees_object, mock_reaction_generator, output_dir
    ):
        failed = SimulationResult(
            t=np.array([]), y=np.zeros((0, 0)), species_labels=[],
            success=False, message="solver failed",
        )
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mock_reaction_generator,
            logger=MagicMock(),
            output_directory=output_dir,
        )
        with patch("bees.enlarger.ODESimulator.simulate", return_value=failed):
            result = enlarger.run()
        assert result.iterations == 1
        assert [row["iteration"] for row in enlarger._iteration_summaries] == [0, 1]

    def test_unmatched_feedback_inhibitor_warns(
        self, mock_bees_object, mock_reaction_generator, output_dir
    ):
        """A misspelled inhibitor silently disables feedback; the final model must warn."""
        logger = MagicMock()
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mock_reaction_generator,
            logger=logger,
            output_directory=output_dir,
        )
        enlarger._initialise_model()
        rxn = _make_reaction(enzyme_label="Enzyme")
        rxn.feedback_inhibitors = {"s": (0.02, 1.0), "dodecanaote": (0.02, 1.0)}
        enlarger.model.core_reactions.append(rxn)

        enlarger._warn_unmatched_feedback_inhibitors()

        warnings = [c.args[0] for c in logger.warning.call_args_list]
        assert len(warnings) == 1
        assert "Enzyme" in warnings[0] and "dodecanaote" in warnings[0]
        assert "'s'" not in warnings[0]

    def test_matched_feedback_inhibitors_do_not_warn(
        self, mock_bees_object, mock_reaction_generator, output_dir
    ):
        logger = MagicMock()
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mock_reaction_generator,
            logger=logger,
            output_directory=output_dir,
        )
        enlarger._initialise_model()
        rxn = _make_reaction(enzyme_label="Enzyme")
        rxn.feedback_inhibitors = {"s": (0.02, 1.0)}
        enlarger.model.core_reactions.append(rxn)

        enlarger._warn_unmatched_feedback_inhibitors()

        logger.warning.assert_not_called()

    def test_iterations_count_simulation_passes(
        self, mock_bees_object, mock_reaction_generator, output_dir
    ):
        """
        Iteration count should reflect simulation passes (interrupt cycles),
        not an outer loop with hidden sub-iterations.
        """
        logger = MagicMock()
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mock_reaction_generator,
            logger=logger,
            output_directory=output_dir,
        )

        # Ensure there is an edge species "p" that can be promoted when the
        # mocked simulator reports interrupt_edge_rates={"p": ...}.
        seed_rxn_p = _make_reaction(
            enzyme_label="Enzyme",
            substrate_label="S",
            reactant_labels=["S"],
            product_labels=["p"],
            stoichiometry={"S": -1, "p": 1},
            km=1.0,
            kcat=10.0,
        )
        seed_rxn_q = _make_reaction(
            enzyme_label="Enzyme",
            substrate_label="S",
            reactant_labels=["S"],
            product_labels=["q"],
            stoichiometry={"S": -1, "q": 1},
            km=1.0,
            kcat=10.0,
        )

        # Initial discovery seeds "p" and "q" into the edge; later calls (promoted species) find nothing.
        enlarger._generate_reactions_for_species = MagicMock(
            side_effect=[[seed_rxn_p, seed_rxn_q]] + [[]] * 10
        )

        # Two interrupted passes then one non-interrupted (converged).
        sims = [_interrupt_result(True, "p"), _interrupt_result(True, "q"), _interrupt_result(False)]

        with patch("bees.enlarger.ODESimulator") as SimCls:
            inst = SimCls.return_value
            inst.simulate.side_effect = sims
            inst.export_ode_equations = MagicMock()
            result = enlarger.run()

        assert result.converged
        assert result.iterations == 3
        assert len(result.simulation_profiles) == 3
        # Iteration summaries include iteration 0 + one row per pass.
        assert len(enlarger._iteration_summaries) == 1 + 3

    def test_max_iterations_termination(self, mock_bees_object, mock_reaction_generator, output_dir):
        """Every pass still promotes something: the loop must stop at max_iterations, unconverged."""
        mock_bees_object.settings.max_iterations = 2
        logger = MagicMock()
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mock_reaction_generator,
            logger=logger,
            output_directory=output_dir,
        )
        seeds = [
            _make_reaction(product_labels=[lab], stoichiometry={"S": -1, lab: 1})
            for lab in ("p", "q", "r")
        ]
        enlarger._generate_reactions_for_species = MagicMock(side_effect=[seeds] + [[]] * 10)

        with patch("bees.enlarger.ODESimulator") as SimCls:
            inst = SimCls.return_value
            inst.simulate.side_effect = [_interrupt_result(True, lab) for lab in ("p", "q", "r")]
            inst.export_ode_equations = MagicMock()
            result = enlarger.run()

        assert result.iterations == 2
        # Two growth passes, then the final model is simulated from t=0.
        assert inst.simulate.call_count == 3
        assert inst.simulate.call_args_list[2].kwargs["interrupt_simulation_tol"] == float("inf")
        assert not result.converged
        assert result.convergence_reason == "Reached max_iterations (2)."

    @patch("bees.enlarger.ODESimulator")
    def test_interrupt_promotes_then_reruns(
        self,
        ode_cls,
        mock_bees_object,
        mock_reaction_generator,
        output_dir,
    ):
        """Interrupt at toleranceInterruptSimulation promotes the edge species; second run completes."""
        from bees.core_edge_model import SpeciesData

        mock_reaction_generator.ensure_estimator_initialized = MagicMock()
        mock_reaction_generator._generate_reactions.return_value = []

        interrupted = SimulationResult(
            t=np.array([0.0, 0.5]),
            y=np.ones((3, 2)),
            species_labels=["S", "Enzyme", "prodx"],
            success=True,
            simulation_interrupted=True,
            interrupt_char_rate=1.0,
            interrupt_edge_rates={"prodx": 0.1},
            max_char_rate=1.0,
            final_char_rate=0.5,
            max_edge_rate_ratio={"prodx": 0.5},
        )
        complete = SimulationResult(
            t=np.array([0.0, 1.0]),
            y=np.ones((3, 2)),
            species_labels=["S", "Enzyme", "prodx"],
            success=True,
            simulation_interrupted=False,
            max_char_rate=1.0,
            final_char_rate=0.5,
            max_edge_rate_ratio={},
        )
        sim_inst = MagicMock()
        sim_inst.simulate.side_effect = [interrupted, complete]
        ode_cls.return_value = sim_inst

        logger = MagicMock()
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mock_reaction_generator,
            logger=logger,
            output_directory=output_dir,
        )
        enlarger._initialise_model()
        enlarger.model.add_edge_species(
            SpeciesData(label="prodx", concentration=0.0, initial_concentration=0.0)
        )
        enlarger.run()
        assert sim_inst.simulate.call_count == 2
        assert enlarger.model.is_core_species("prodx")

    @pytest.mark.parametrize(
        "peak_rr, expect_promoted",
        [(0.5, True), (1e-7, False)],
    )
    @patch("bees.enlarger.ODESimulator")
    def test_full_run_promotes_by_peak_ratio(
        self,
        ode_cls,
        mock_bees_object,
        mock_reaction_generator,
        output_dir,
        peak_rr,
        expect_promoted,
    ):
        """Non-interrupted run: promote edge species whose peak rr >= toleranceMoveToCore (RMG)."""
        from bees.core_edge_model import SpeciesData

        mock_bees_object.settings.toleranceInterruptSimulation = 1.0
        mock_bees_object.settings.end_time = 100.0
        mock_reaction_generator.ensure_estimator_initialized = MagicMock()
        mock_reaction_generator._generate_reactions.return_value = []

        def _complete(ratios):
            return SimulationResult(
                t=np.array([0.0, 100.0]),
                y=np.ones((3, 2)),
                species_labels=["S", "Enzyme", "prodx"],
                success=True,
                simulation_interrupted=False,
                max_char_rate=1.0,
                final_char_rate=0.5,
                max_edge_rate_ratio=ratios,
            )

        sim_inst = MagicMock()
        sim_inst.simulate.side_effect = [_complete({"prodx": peak_rr}), _complete({})]
        ode_cls.return_value = sim_inst

        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mock_reaction_generator,
            logger=MagicMock(),
            output_directory=output_dir,
        )
        enlarger._initialise_model()
        enlarger.model.add_edge_species(
            SpeciesData(label="prodx", concentration=0.0, initial_concentration=0.0)
        )
        result = enlarger.run()

        assert result.converged
        assert enlarger.model.is_core_species("prodx") is expect_promoted
        assert sim_inst.simulate.call_count == (2 if expect_promoted else 1)
        if expect_promoted:
            assert "prodx" not in enlarger._edge_peak_rr
        else:
            assert enlarger._edge_peak_rr["prodx"] == (peak_rr, 1)

    def test_record_edge_peaks_keeps_largest(self, mock_bees_object, mock_reaction_generator, output_dir):
        """The recorded edge peak is the max over growth simulations, with its iteration."""
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mock_reaction_generator,
            logger=MagicMock(),
            output_directory=output_dir,
        )

        def _res(ratios, success=True):
            return SimulationResult(
                t=np.array([0.0]), y=np.zeros((1, 1)), species_labels=["x"],
                success=success, max_edge_rate_ratio=ratios,
            )

        enlarger._record_edge_peaks(_res({"a": 2e-3, "b": 1e-4}), 1)
        enlarger._record_edge_peaks(_res({"a": 5e-3, "b": 5e-5}), 2)
        enlarger._record_edge_peaks(_res({"a": 9.0}, success=False), 3)
        assert enlarger._edge_peak_rr == {"a": (5e-3, 2), "b": (1e-4, 1)}

    @patch("bees.enlarger.ODESimulator")
    def test_terminated_growth_run_is_replaced_by_full_length_profile(
        self,
        ode_cls,
        mock_bees_object,
        mock_reaction_generator,
        output_dir,
    ):
        """Early termination converges; exported profile is re-simulated to end_time (RMG)."""
        from bees.core_edge_model import SpeciesData

        mock_bees_object.settings.toleranceInterruptSimulation = 1.0
        mock_bees_object.settings.end_time = 100.0
        mock_bees_object.settings.termination_rate_ratio = 0.01
        mock_reaction_generator.ensure_estimator_initialized = MagicMock()
        mock_reaction_generator._generate_reactions.return_value = []

        terminated = SimulationResult(
            t=np.array([0.0, 20.0]),
            y=np.ones((3, 2)),
            species_labels=["S", "Enzyme", "prodx"],
            success=True,
            max_char_rate=1.0,
            final_char_rate=0.005,
            max_edge_rate_ratio={"prodx": 1e-7},
            termination_reason="rate_ratio",
        )
        full = SimulationResult(
            t=np.array([0.0, 50.0, 100.0]),
            y=np.ones((3, 3)),
            species_labels=["S", "Enzyme", "prodx"],
            success=True,
        )
        sim_inst = MagicMock()
        sim_inst.simulate.side_effect = [terminated, full]
        ode_cls.return_value = sim_inst

        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mock_reaction_generator,
            logger=MagicMock(),
            output_directory=output_dir,
        )
        enlarger._initialise_model()
        enlarger.model.add_edge_species(
            SpeciesData(label="prodx", concentration=0.0, initial_concentration=0.0)
        )
        result = enlarger.run()

        assert result.converged
        assert result.convergence_reason == "Termination rate ratio reached."
        assert not enlarger.model.is_core_species("prodx")
        assert sim_inst.simulate.call_count == 2
        growth_kwargs = sim_inst.simulate.call_args_list[0].kwargs
        assert growth_kwargs["termination_rate_ratio"] == 0.01
        final_kwargs = sim_inst.simulate.call_args_list[1].kwargs
        assert final_kwargs["end_time"] == 100.0
        assert "termination_rate_ratio" not in final_kwargs
        assert result.simulation_profiles[-1] is full

    @patch("bees.enlarger.ODESimulator")
    def test_no_prune_when_simulation_interrupted(
        self,
        ode_cls,
        mock_bees_object,
        mock_reaction_generator,
        output_dir,
    ):
        """An interrupted pass never prunes the edge, even when nothing is promotable and it resumes."""
        mock_bees_object.settings.toleranceKeepInEdge = 0.01
        mock_bees_object.settings.minEdgeIterationsForPrune = 0
        mock_reaction_generator.ensure_estimator_initialized = MagicMock()
        mock_reaction_generator._generate_reactions.return_value = [
            _make_reaction(
                enzyme_label="Enzyme",
                substrate_label="S",
                reactant_labels=["S"],
                product_labels=["prodx"],
                stoichiometry={"S": -1, "prodx": 1},
                kcat=10.0,
            ),
        ]
        sim_inst = MagicMock()
        sim_inst.simulate.return_value = SimulationResult(
            t=np.array([0.0, 1.0]),
            y=np.ones((2, 2)),
            species_labels=["S", "Enzyme"],
            success=True,
            simulation_interrupted=True,
            max_char_rate=1.0,
            final_char_rate=0.5,
            max_edge_rate_ratio={"prodx": 0.5},
        )
        ode_cls.return_value = sim_inst

        logger = MagicMock()
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mock_reaction_generator,
            logger=logger,
            output_directory=output_dir,
        )
        enlarger.model.prune_edge = MagicMock()
        enlarger.run()
        enlarger.model.prune_edge.assert_not_called()


class TestResumeAfterInterrupt:
    """Interrupt with nothing promotable resumes the pass from the interrupt time."""

    @staticmethod
    def _result(*, t_end, interrupted, success=True, message=""):
        return SimulationResult(
            t=np.array([0.0, t_end]),
            y=np.ones((3, 2)),
            species_labels=["S", "Enzyme", "prodx"],
            success=success,
            message=message,
            simulation_interrupted=interrupted,
            interrupt_char_rate=1.0 if interrupted else 0.0,
            interrupt_edge_rates={"prodx": 1e-12} if interrupted else {},
            max_char_rate=1.0,
            final_char_rate=0.5,
            max_edge_rate_ratio={},
        )

    def _enlarger(self, mock_bees_object, mock_reaction_generator, output_dir):
        from bees.core_edge_model import SpeciesData

        mock_bees_object.settings.toleranceInterruptSimulation = 1.0
        mock_reaction_generator.ensure_estimator_initialized = MagicMock()
        mock_reaction_generator._generate_reactions.return_value = []
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mock_reaction_generator,
            logger=MagicMock(),
            output_directory=output_dir,
        )
        enlarger._initialise_model()
        enlarger.model.add_edge_species(
            SpeciesData(label="prodx", concentration=0.0, initial_concentration=0.0)
        )
        return enlarger

    @patch("bees.enlarger.ODESimulator")
    def test_clean_resume_resimulates_from_zero_and_converges(
        self, ode_cls, mock_bees_object, mock_reaction_generator, output_dir
    ):
        full = self._result(t_end=100.0, interrupted=False)
        sim_inst = MagicMock()
        sim_inst.simulate.side_effect = [
            self._result(t_end=30.0, interrupted=True),
            self._result(t_end=70.0, interrupted=False),
            full,
        ]
        ode_cls.return_value = sim_inst
        enlarger = self._enlarger(mock_bees_object, mock_reaction_generator, output_dir)
        result = enlarger.run()

        assert sim_inst.simulate.call_count == 3
        resume_kwargs = sim_inst.simulate.call_args_list[1].kwargs
        assert resume_kwargs["end_time"] == pytest.approx(70.0)
        for key in ("method", "rtol", "atol", "max_wall_time_s"):
            assert key in resume_kwargs
        final_kwargs = sim_inst.simulate.call_args_list[2].kwargs
        assert final_kwargs["end_time"] == 100.0
        assert final_kwargs["interrupt_simulation_tol"] == float("inf")
        assert result.converged
        assert result.simulation_profiles[-1] is full

    @patch("bees.enlarger.ODESimulator")
    def test_failed_resume_is_reported(
        self, ode_cls, mock_bees_object, mock_reaction_generator, output_dir
    ):
        sim_inst = MagicMock()
        sim_inst.simulate.side_effect = [
            self._result(t_end=30.0, interrupted=True),
            self._result(t_end=40.0, interrupted=False, success=False, message="boom"),
        ]
        ode_cls.return_value = sim_inst
        enlarger = self._enlarger(mock_bees_object, mock_reaction_generator, output_dir)
        result = enlarger.run()

        assert not result.converged
        assert result.convergence_reason == "ODE failure on resume: boom"

    @patch("bees.enlarger.ODESimulator")
    def test_final_profile_resimulated_when_model_grew_after_last_pass(
        self, ode_cls, mock_bees_object, mock_reaction_generator, output_dir
    ):
        """Hitting max_iterations right after a promotion must not export the stale, cut-off pass."""
        mock_bees_object.settings.max_iterations = 1
        interrupted = self._result(t_end=30.0, interrupted=True)
        interrupted.interrupt_edge_rates = {"prodx": 1.0}
        full = self._result(t_end=100.0, interrupted=False)
        sim_inst = MagicMock()
        sim_inst.simulate.side_effect = [interrupted, full]
        ode_cls.return_value = sim_inst
        enlarger = self._enlarger(mock_bees_object, mock_reaction_generator, output_dir)
        result = enlarger.run()

        assert enlarger.model.is_core_species("prodx")
        assert result.convergence_reason == "Reached max_iterations (1)."
        assert sim_inst.simulate.call_count == 2
        assert sim_inst.simulate.call_args_list[1].kwargs["interrupt_simulation_tol"] == float("inf")
        assert result.simulation_profiles[-1] is full

    def test_reaction_level_promotion_keeps_strongest(
        self, mock_bees_object, mock_reaction_generator, output_dir
    ):
        from bees.core_edge_model import SpeciesData

        mock_bees_object.settings.toleranceMoveEdgeReactionToCore = 0.1
        mock_bees_object.settings.max_num_objects_per_iter = 1
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mock_reaction_generator,
            logger=MagicMock(),
            output_directory=output_dir,
        )
        for label in ("weak", "strong"):
            enlarger.model.add_edge_species(SpeciesData(label=label, concentration=0.0))
        enlarger._reaction_obj_by_sig = {
            "sig_weak": _make_reaction(reactant_labels=["S"], product_labels=["weak"]),
            "sig_strong": _make_reaction(reactant_labels=["S"], product_labels=["strong"]),
        }
        sim = SimulationResult(
            t=np.array([0.0]), y=np.zeros((1, 1)), species_labels=["S"], success=True,
        )
        sim.max_edge_reaction_dlnaccum = {"sig_weak": 1.0, "sig_strong": 5.0}
        assert enlarger._reaction_level_promotions(sim) == ["strong"]

    def test_reaction_level_promotion_ignores_label_case(
        self, mock_bees_object, mock_reaction_generator, output_dir
    ):
        from bees.core_edge_model import SpeciesData

        mock_bees_object.settings.toleranceMoveEdgeReactionToCore = 0.1
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mock_reaction_generator,
            logger=MagicMock(),
            output_directory=output_dir,
        )
        enlarger.model.add_edge_species(SpeciesData(label="Malonyl-CoA", concentration=0.0))
        enlarger._reaction_obj_by_sig = {
            "sig": _make_reaction(reactant_labels=["S"], product_labels=["malonyl-coa"]),
        }
        sim = SimulationResult(
            t=np.array([0.0]), y=np.zeros((1, 1)), species_labels=["S"], success=True,
        )
        sim.max_edge_reaction_dlnaccum = {"sig": 5.0}
        assert enlarger._reaction_level_promotions(sim) == ["Malonyl-CoA"]


class TestResetToInitialConsistency:
    """Every outer enlargement iteration must start from initial concentrations."""

    @patch("bees.enlarger.ODESimulator")
    def test_reset_called_before_each_outer_iteration(
        self, ode_cls, mock_bees_object, mock_reaction_generator, output_dir
    ):
        """Each enlargement iteration starts from initial concentrations, even though
        the simulator leaves final-time concentrations in the model."""
        from bees.core_edge_model import SpeciesData

        mock_reaction_generator.ensure_estimator_initialized = MagicMock()
        mock_reaction_generator._generate_reactions.return_value = []
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mock_reaction_generator,
            logger=MagicMock(),
            output_directory=output_dir,
        )
        enlarger._initialise_model()
        enlarger.model.add_edge_species(
            SpeciesData(label="prodx", concentration=0.0, initial_concentration=0.0)
        )

        results = iter([_interrupt_result(True, "prodx"), _interrupt_result(False)])
        at_initial = []

        def simulate(**_kwargs):
            variable = [sp for sp in enlarger.model.core_species if not sp.constant]
            at_initial.append(all(sp.concentration == sp.initial_concentration for sp in variable))
            for sp in variable:
                sp.concentration = sp.initial_concentration + 123.0
            return next(results)

        ode_cls.return_value.simulate.side_effect = simulate
        enlarger.run()

        assert len(at_initial) == 2
        assert all(at_initial), at_initial
