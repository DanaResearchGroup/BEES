"""
Tests for the kinetics_estimator module.

To run: pytest -v tests/test_kinetics_estimator.py
"""

import pytest

from bees.kinetics_estimator import (
    EstimatedKinetics,
    BaseKineticsEstimator,
    build_estimator,
)


class TestEstimatedKinetics:
    def test_defaults(self):
        k = EstimatedKinetics()
        assert k.km is None
        assert k.kcat is None
        assert k.source == "estimator"

    def test_with_values(self):
        k = EstimatedKinetics(km=0.1, kcat=100.0)
        assert k.km == 0.1
        assert k.kcat == 100.0

    def test_frozen(self):
        k = EstimatedKinetics(km=0.1)
        with pytest.raises(AttributeError):
            k.km = 0.2


class TestBuildEstimator:
    def test_none_returns_none(self):
        assert build_estimator(None) is None

    def test_catpred_returns_estimator(self):
        est = build_estimator("catpred")
        assert est is not None
        assert est.name == "catpred"

    def test_catpred_with_include_sd(self):
        est = build_estimator("catpred", include_sd=True)
        assert est.include_sd is True

    def test_unknown_raises(self):
        with pytest.raises(ValueError, match="Unknown kinetics_estimator"):
            build_estimator("unknown_backend")


class TestBaseKineticsEstimator:
    def test_estimate_not_implemented(self):
        base = BaseKineticsEstimator()
        with pytest.raises(NotImplementedError):
            base.estimate(
                enzyme_sequence="MKTAY",
                reactant_smiles={"S": "C"},
            )
