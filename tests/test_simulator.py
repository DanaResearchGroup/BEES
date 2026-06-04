"""
Tests for the ODE simulator module.

To run:  pytest -v tests/test_simulator.py
"""


import pytest
import numpy as np
from itertools import chain, repeat
from unittest.mock import MagicMock, patch

from bees.core_edge_model import CoreEdgeModel, SpeciesData
from bees.simulator import ODESimulator


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
    vmax=None,
    km_per_substrate=None,
):
    kinetics = MagicMock()
    kinetics.km = km
    kinetics.kcat = kcat
    kinetics.vmax = vmax
    kinetics.km_per_substrate = km_per_substrate

    rxn = MagicMock()
    rxn.enzyme_label = enzyme_label
    rxn.substrate_label = substrate_label
    rxn.reactant_labels = reactant_labels or ["S"]
    rxn.product_labels = product_labels or ["P"]
    rxn.stoichiometry = stoichiometry or {"S": -1, "P": 1}
    rxn.rate_law = rate_law
    rxn.kinetics = kinetics
    return rxn


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def simple_model():
    """
    Minimal model: S -> P catalysed by E.
    Uses Michaelis-Menten kinetics with Km=1.0 mM, kcat=10.0 s^-1.
    Initial [S] = 10.0 mM, [E] = 0.01 mM.
    """
    model = CoreEdgeModel()

    model.add_core_species(SpeciesData(
        label="S", concentration=10.0, initial_concentration=10.0,
    ))
    model.add_core_species(SpeciesData(
        label="P", concentration=0.0, initial_concentration=0.0,
    ))
    model.add_core_species(SpeciesData(
        label="Enzyme", concentration=0.01, initial_concentration=0.01,
        is_enzyme=True, constant=True,
    ))

    rxn = _make_reaction(
        enzyme_label="Enzyme",
        substrate_label="S",
        reactant_labels=["S"],
        product_labels=["P"],
        stoichiometry={"S": -1, "P": 1},
        km=1.0,
        kcat=10.0,
    )
    model.core_reactions.append(rxn)
    return model


# ---------------------------------------------------------------------------
# Test ODESimulator
# ---------------------------------------------------------------------------

class TestODESimulator:
    def test_basic_simulation(self, simple_model):
        """S -> P should show S decreasing and P increasing."""
        sim = ODESimulator(simple_model)
        result = sim.simulate(end_time=100.0, time_step=10.0)

        assert result.success
        assert result.t.shape[0] > 1
        assert result.species_labels == ["S", "P", "Enzyme"]
        assert np.all(result.y >= 0.0)

        # S should decrease
        s_idx = 0
        assert result.y[s_idx, -1] < result.y[s_idx, 0]

        # P should increase
        p_idx = 1
        assert result.y[p_idx, -1] > result.y[p_idx, 0]

    def test_mass_conservation(self, simple_model):
        """Total [S] + [P] should be conserved (=10 mM)."""
        sim = ODESimulator(simple_model)
        result = sim.simulate(end_time=100.0, time_step=1.0)

        s_idx = 0
        p_idx = 1
        total = result.y[s_idx, :] + result.y[p_idx, :]
        np.testing.assert_allclose(total, 10.0, atol=0.01)

    def test_enzyme_constant(self, simple_model):
        """Enzyme concentration should remain constant."""
        sim = ODESimulator(simple_model)
        result = sim.simulate(end_time=50.0, time_step=5.0)
        e_idx = 2
        np.testing.assert_allclose(result.y[e_idx, :], 0.01, atol=1e-12)

    def test_empty_model(self):
        """An empty model should simulate without errors."""
        model = CoreEdgeModel()
        sim = ODESimulator(model)
        result = sim.simulate(end_time=10.0)
        assert result.success
        assert result.species_labels == []

    def test_evaluate_core_rates(self, simple_model):
        """Core rates should be non-zero for the simple S->P model."""
        sim = ODESimulator(simple_model)
        result = sim.simulate(end_time=10.0, time_step=1.0)
        core_rates = sim.evaluate_core_rates(result)
        # S should have a negative rate (being consumed)
        assert core_rates.get("s", 0) < 0
        # P should have a positive rate (being produced)
        assert core_rates.get("p", 0) > 0

    def test_flux_trajectory_metrics_attached(self, simple_model):
        """simulate() should attach RMG-style max/final R_char and peak edge dicts."""
        sim = ODESimulator(simple_model)
        result = sim.simulate(end_time=10.0, time_step=1.0)
        assert result.max_char_rate >= 0.0
        assert result.final_char_rate >= 0.0
        assert isinstance(result.max_edge_rate_ratio, dict)
        assert isinstance(result.peak_edge_signed_rate, dict)
        assert result.max_char_rate >= result.final_char_rate or result.max_char_rate > 0

    def test_model_concentrations_updated(self, simple_model):
        """After simulation, model concentrations should reflect the final state."""
        sim = ODESimulator(simple_model)
        sim.simulate(end_time=100.0)
        # S should have decreased from 10.0
        s_sp = simple_model.get_core_species_by_label("S")
        assert s_sp.concentration < 10.0
        assert s_sp.concentration >= 0.0


# ---------------------------------------------------------------------------
# Stepwise interrupt tests
# ---------------------------------------------------------------------------

def _make_core_edge_model_with_immediate_flux():
    """
    Build a model where a core reaction S -> P (core) immediately
    produces edge flux for species E_prod via P -> E_prod (edge).
    Both reactions run from t=0, so edge rr should exceed any
    reasonable tolerance almost immediately.
    """
    model = CoreEdgeModel()

    model.add_core_species(SpeciesData(
        label="S", concentration=10.0, initial_concentration=10.0,
    ))
    model.add_core_species(SpeciesData(
        label="P", concentration=5.0, initial_concentration=5.0,
    ))
    model.add_core_species(SpeciesData(
        label="Enz", concentration=0.01, initial_concentration=0.01,
        is_enzyme=True, constant=True,
    ))

    model.add_edge_species(SpeciesData(
        label="E_prod", concentration=0.0, initial_concentration=0.0,
    ))

    rxn_core = _make_reaction(
        enzyme_label="Enz", substrate_label="S",
        reactant_labels=["S"], product_labels=["P"],
        stoichiometry={"S": -1, "P": 1},
        km=1.0, kcat=10.0,
    )
    model.core_reactions.append(rxn_core)

    rxn_edge = _make_reaction(
        enzyme_label="Enz", substrate_label="P",
        reactant_labels=["P"], product_labels=["E_prod"],
        stoichiometry={"P": -1, "E_prod": 1},
        km=1.0, kcat=10.0,
    )
    model.edge_reactions.append(rxn_edge)
    return model


class TestStepwiseInterrupt:
    """Tests for the RMG-style stepwise interrupt path."""

    def test_interrupt_occurs_at_positive_time(self):
        """With immediate edge flux and tiny tolerance, interrupt must be at t > 0."""
        model = _make_core_edge_model_with_immediate_flux()
        sim = ODESimulator(model)
        result = sim.simulate(
            end_time=100.0,
            interrupt_simulation_tol=1e-8,
        )
        assert result.simulation_interrupted
        assert result.t[-1] > 0.0, "Interrupt must occur after at least one step (t > 0)"

    def test_interrupt_rates_populated(self):
        """When interrupted, interrupt_char_rate and interrupt_edge_rates must be set."""
        model = _make_core_edge_model_with_immediate_flux()
        sim = ODESimulator(model)
        result = sim.simulate(
            end_time=100.0,
            interrupt_simulation_tol=1e-8,
        )
        assert result.simulation_interrupted
        assert result.interrupt_char_rate > 0.0
        assert len(result.interrupt_edge_rates) > 0
        assert "e_prod" in result.interrupt_edge_rates

    def test_no_interrupt_without_tolerance(self):
        """Without interrupt_simulation_tol, simulation should run to end_time."""
        model = _make_core_edge_model_with_immediate_flux()
        sim = ODESimulator(model)
        result = sim.simulate(end_time=1.0, time_step=0.1)
        assert not result.simulation_interrupted
        assert result.t[-1] >= 0.99

    def test_no_interrupt_with_high_tolerance(self):
        """With a very high tolerance, no interrupt should fire."""
        model = _make_core_edge_model_with_immediate_flux()
        sim = ODESimulator(model)
        result = sim.simulate(
            end_time=1.0,
            interrupt_simulation_tol=1e6,
        )
        assert not result.simulation_interrupted

    def test_continuous_path_still_works(self):
        """When no edge species exist, the continuous path should be used."""
        model = CoreEdgeModel()
        model.add_core_species(SpeciesData(
            label="A", concentration=5.0, initial_concentration=5.0,
        ))
        model.add_core_species(SpeciesData(
            label="B", concentration=0.0, initial_concentration=0.0,
        ))
        model.add_core_species(SpeciesData(
            label="Enz", concentration=0.01, initial_concentration=0.01,
            is_enzyme=True, constant=True,
        ))
        rxn = _make_reaction(
            enzyme_label="Enz", substrate_label="A",
            reactant_labels=["A"], product_labels=["B"],
            stoichiometry={"A": -1, "B": 1},
            km=1.0, kcat=5.0,
        )
        model.core_reactions.append(rxn)

        sim = ODESimulator(model)
        result = sim.simulate(
            end_time=10.0, time_step=1.0,
            interrupt_simulation_tol=1e-5,
        )
        assert not result.simulation_interrupted
        assert result.t[-1] >= 9.99

    def test_history_has_multiple_points(self):
        """Stepwise result should contain more than one time point."""
        model = _make_core_edge_model_with_immediate_flux()
        sim = ODESimulator(model)
        result = sim.simulate(
            end_time=100.0,
            interrupt_simulation_tol=1e-8,
        )
        assert len(result.t) >= 2
        assert result.y.shape[1] == len(result.t)

    def test_stepwise_dt_is_capped_at_5s(self):
        """
        The outer-step dt in stepwise mode should never exceed 5 s.

        We check this indirectly via the returned time history, which stores
        the per-step endpoints produced by the outer step loop.
        """
        model = _make_core_edge_model_with_immediate_flux()
        sim = ODESimulator(model)
        result = sim.simulate(
            end_time=50.0,
            interrupt_simulation_tol=1e9,
        )
        assert not result.simulation_interrupted
        assert len(result.t) >= 2
        dts = np.diff(result.t)
        assert np.all(dts <= 5.0 + 1e-12)

    def test_stepwise_wall_clock_limit_stops_pass(self):
        """
        When max_wall_time_s is exceeded after an outer step, the pass should end
        with success=False (enlarger treats this as ODE failure).
        """
        model = _make_core_edge_model_with_immediate_flux()
        sim = ODESimulator(model)
        mono_vals = chain([0.0, 0.0, 100.0], repeat(100.0))
        with patch("bees.simulator.time.monotonic", side_effect=lambda: next(mono_vals)):
            result = sim.simulate(
                end_time=50.0,
                interrupt_simulation_tol=1e9,
                max_wall_time_s=1.0,
                stepwise_heartbeat_interval_s=0.0,
            )
        assert not result.success
        assert "max_wall_time_per_iteration" in result.message
        assert not result.simulation_interrupted

    def test_simulate_accepts_custom_method_and_tolerances(self):
        """Custom method/rtol/atol should not break the stepwise path."""
        model = _make_core_edge_model_with_immediate_flux()
        sim = ODESimulator(model)
        result = sim.simulate(
            end_time=1.0,
            interrupt_simulation_tol=1e9,
            method="BDF",
            rtol=1e-5,
            atol=1e-7,
            stepwise_heartbeat_interval_s=0.0,
        )
        assert result.success
        assert not result.simulation_interrupted
