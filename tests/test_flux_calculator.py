"""
Tests for the flux calculator module.

To run:  pytest -v tests/test_flux_calculator.py
"""


from unittest.mock import MagicMock

from bees.flux_calculator import (
    calculate_characteristic_rate,
    calculate_species_rates,
    compute_mm_rate,
    identify_insignificant_species_from_peak_ratios,
    identify_significant_species_at_interrupt,
)


# ---------------------------------------------------------------------------
# Helper to build a mock GeneratedReaction
# ---------------------------------------------------------------------------

def _make_reaction(
    enzyme_label="Hexokinase",
    substrate_label="Glucose",
    reactant_labels=None,
    product_labels=None,
    stoichiometry=None,
    rate_law="Michaelis-Menten",
    km=0.1,
    kcat=100.0,
    vmax=None,
    km_per_substrate=None,
):
    """Create a minimal mock GeneratedReaction for testing."""
    if reactant_labels is None:
        reactant_labels = ["Glucose", "ATP"]
    if product_labels is None:
        product_labels = ["Glucose-6P", "ADP"]
    if stoichiometry is None:
        stoichiometry = {"Glucose": -1, "ATP": -1, "Glucose-6P": 1, "ADP": 1}

    kinetics = MagicMock()
    kinetics.km = km
    kinetics.kcat = kcat
    kinetics.vmax = vmax
    kinetics.km_per_substrate = km_per_substrate

    rxn = MagicMock()
    rxn.enzyme_label = enzyme_label
    rxn.substrate_label = substrate_label
    rxn.reactant_labels = reactant_labels
    rxn.product_labels = product_labels
    rxn.stoichiometry = stoichiometry
    rxn.rate_law = rate_law
    rxn.kinetics = kinetics
    return rxn


# ---------------------------------------------------------------------------
# compute_mm_rate
# ---------------------------------------------------------------------------

class TestComputeMmRate:
    def test_basic_mm_rate(self):
        rxn = _make_reaction(km=0.1, kcat=100.0)
        conc = {"glucose": 1.0, "atp": 2.0}
        enz_conc = {"hexokinase": 0.001}
        rate = compute_mm_rate(rxn, conc, enz_conc)
        # v = kcat * [E] * [S1]/(Km+[S1]) * [S2]/(Km+[S2])
        expected = 100.0 * 0.001 * (1.0 / (0.1 + 1.0)) * (2.0 / (0.1 + 2.0))
        assert abs(rate - expected) < 1e-10

    def test_vmax_fallback(self):
        rxn = _make_reaction(km=0.5, kcat=None, vmax=10.0)
        conc = {"glucose": 1.0, "atp": 1.0}
        enz_conc = {}  # no enzyme conc needed for vmax mode
        rate = compute_mm_rate(rxn, conc, enz_conc)
        expected = 10.0 * (1.0 / (0.5 + 1.0)) * (1.0 / (0.5 + 1.0))
        assert abs(rate - expected) < 1e-10

    def test_no_kinetics(self):
        rxn = _make_reaction()
        rxn.kinetics = None
        rate = compute_mm_rate(rxn, {}, {})
        assert rate == 0.0

    def test_no_rate_law(self):
        rxn = _make_reaction()
        rxn.rate_law = None
        rate = compute_mm_rate(rxn, {}, {})
        assert rate == 0.0

    def test_zero_substrate(self):
        rxn = _make_reaction(km=0.1, kcat=100.0)
        conc = {"glucose": 0.0, "atp": 1.0}
        enz_conc = {"hexokinase": 0.001}
        rate = compute_mm_rate(rxn, conc, enz_conc)
        assert rate == 0.0

    def test_per_substrate_km(self):
        rxn = _make_reaction(
            km=None,
            kcat=50.0,
            km_per_substrate={"Glucose": 0.2, "ATP": 0.5},
        )
        conc = {"glucose": 1.0, "atp": 2.0}
        enz_conc = {"hexokinase": 0.01}
        rate = compute_mm_rate(rxn, conc, enz_conc)
        expected = 50.0 * 0.01 * (1.0 / (0.2 + 1.0)) * (2.0 / (0.5 + 2.0))
        assert abs(rate - expected) < 1e-10

    def test_negative_concentration_guarded(self):
        rxn = _make_reaction(km=0.1, kcat=100.0)
        conc = {"glucose": -0.5, "atp": 1.0}
        enz_conc = {"hexokinase": 0.001}
        rate = compute_mm_rate(rxn, conc, enz_conc)
        assert rate == 0.0


# ---------------------------------------------------------------------------
# calculate_species_rates
# ---------------------------------------------------------------------------

class TestCalculateSpeciesRates:
    def test_single_reaction(self):
        rxn = _make_reaction(km=0.1, kcat=100.0)
        conc = {"glucose": 1.0, "atp": 2.0}
        enz_conc = {"hexokinase": 0.001}
        rates = calculate_species_rates([rxn], conc, enz_conc)
        v = compute_mm_rate(rxn, conc, enz_conc)
        assert abs(rates.get("glucose", 0) - (-v)) < 1e-10
        assert abs(rates.get("atp", 0) - (-v)) < 1e-10
        assert abs(rates.get("glucose-6p", 0) - v) < 1e-10
        assert abs(rates.get("adp", 0) - v) < 1e-10


# ---------------------------------------------------------------------------
# calculate_characteristic_rate
# ---------------------------------------------------------------------------

class TestCharacteristicRate:
    def test_basic(self):
        rates = {"a": 3.0, "b": 4.0}
        r_char = calculate_characteristic_rate(rates)
        assert abs(r_char - 5.0) < 1e-10  # sqrt(9 + 16)

    def test_empty(self):
        assert calculate_characteristic_rate({}) == 0.0


# ---------------------------------------------------------------------------
# identify_significant_species_at_interrupt
# ---------------------------------------------------------------------------


class TestIdentifySignificantAtInterrupt:
    def test_normal_promotion(self):
        """Species exceeding tol are returned, sorted by rr desc."""
        edge = {"a": 1.0, "b": 0.5, "c": 0.001}
        result = identify_significant_species_at_interrupt(
            edge, char_rate=1.0, tol_move_to_core=0.1
        )
        labels = [sf.label for sf in result]
        assert labels == ["a", "b"]
        assert result[0].normalized_rate > result[1].normalized_rate

    def test_max_objects_truncation(self):
        """Only top max_objects candidates are returned."""
        edge = {"a": 1.0, "b": 0.8, "c": 0.6}
        result = identify_significant_species_at_interrupt(
            edge, char_rate=1.0, tol_move_to_core=0.1, max_objects=2
        )
        assert len(result) == 2
        assert result[0].label == "a"

    def test_flat_core_abs_flux_floor(self):
        """When char_rate == 0, use abs_flux_floor instead of ratio."""
        edge = {"big": 0.5, "tiny": 1e-15, "medium": 0.01}
        result = identify_significant_species_at_interrupt(
            edge, char_rate=0.0, tol_move_to_core=0.1,
            abs_flux_floor=1e-3,
        )
        labels = [sf.label for sf in result]
        assert "big" in labels
        assert "medium" in labels
        assert "tiny" not in labels
        assert result[0].normalized_rate == float("inf")

    def test_flat_core_truncation(self):
        """In flat-core mode, max_objects still caps output."""
        edge = {"a": 1.0, "b": 0.5, "c": 0.3}
        result = identify_significant_species_at_interrupt(
            edge, char_rate=0.0, tol_move_to_core=0.1,
            abs_flux_floor=1e-12, max_objects=1,
        )
        assert len(result) == 1
        assert result[0].label == "a"

    def test_empty_edge(self):
        result = identify_significant_species_at_interrupt(
            {}, char_rate=1.0, tol_move_to_core=0.1
        )
        assert result == []

    def test_none_above_threshold(self):
        edge = {"a": 0.001}
        result = identify_significant_species_at_interrupt(
            edge, char_rate=1.0, tol_move_to_core=0.1
        )
        assert result == []

    def test_flat_core_all_below_floor(self):
        """When char_rate == 0 and all |rates| < floor, nothing promoted."""
        edge = {"a": 1e-20, "b": 1e-18}
        result = identify_significant_species_at_interrupt(
            edge, char_rate=0.0, tol_move_to_core=0.1,
            abs_flux_floor=1e-12,
        )
        assert result == []


# ---------------------------------------------------------------------------
# identify_insignificant_species_from_peak_ratios
# ---------------------------------------------------------------------------


class TestIdentifyInsignificantFromPeakRatios:
    def test_prune_low_peak_ratio(self):
        max_rr = {"a": 0.5, "b": 0.001}
        out = identify_insignificant_species_from_peak_ratios(
            max_rr, max_char_rate=1.0, tol_keep_in_edge=0.01
        )
        assert "b" in out
        assert "a" not in out

    def test_ineligible_skipped(self):
        max_rr = {"a": 0.001}
        out = identify_insignificant_species_from_peak_ratios(
            max_rr,
            max_char_rate=1.0,
            tol_keep_in_edge=0.01,
            ineligible_for_prune={"a"},
        )
        assert out == set()

    def test_zero_max_char(self):
        assert identify_insignificant_species_from_peak_ratios(
            {"a": 0.5}, 0.0, 0.01
        ) == set()
