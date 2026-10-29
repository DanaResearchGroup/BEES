"""
Tests for the kinetics_estimator module.

To run: pytest -v tests/test_kinetics_estimator.py
"""

import pytest

from bees.kinetics_estimator import (
    EstimatedKinetics,
    build_estimator,
)


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


class TestPersistentCache:
    """Persistent prediction cache (default ON; CatPred is the slow step)."""

    def test_on_by_default_global_path(self, monkeypatch, tmp_path):
        from bees.kinetics_estimator import CatPredEstimator
        monkeypatch.delenv("BEES_CATPRED_CACHE", raising=False)
        monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
        est = CatPredEstimator()
        assert est._cache_path == str(tmp_path / "bees" / "catpred_predictions.pkl")

    def test_disabled_when_off(self, monkeypatch):
        from bees.kinetics_estimator import CatPredEstimator
        for val in ("off", "0", "none", "False"):
            monkeypatch.setenv("BEES_CATPRED_CACHE", val)
            est = CatPredEstimator()
            assert est._cache_path is None, f"{val!r} should disable the cache"
            est._save_persistent_cache()  # no-op, must not raise

    def test_roundtrip(self, tmp_path, monkeypatch):
        from bees.kinetics_estimator import CatPredEstimator
        cache = tmp_path / "catpred.pkl"
        monkeypatch.setenv("BEES_CATPRED_CACHE", str(cache))

        est = CatPredEstimator(include_sd=True)
        assert est._cache_path == str(cache)
        key = ("SEQ", (("S", "C"),), None, True)
        est._memo[key] = EstimatedKinetics(km=0.05, kcat=3.14, source="catpred")
        est._save_persistent_cache()
        assert cache.exists()

        est2 = CatPredEstimator(include_sd=True)
        assert key in est2._memo
        assert est2._memo[key].kcat == 3.14
        assert est2._memo[key].km == 0.05

    def test_failed_prediction_not_persisted(self, tmp_path, monkeypatch):
        """A failed CatPred run is reused within the run but never written to disk."""
        import subprocess

        from bees.kinetics_estimator import CatPredEstimator

        cache = tmp_path / "catpred.pkl"
        monkeypatch.setenv("BEES_CATPRED_CACHE", str(cache))
        monkeypatch.chdir(tmp_path)
        monkeypatch.setattr(CatPredEstimator, "CATPRED_DIR", str(tmp_path))
        calls = []

        def _fail(cmd, **kwargs):
            calls.append(cmd)
            raise subprocess.CalledProcessError(1, cmd, stderr="boom")

        monkeypatch.setattr(subprocess, "run", _fail)

        est = CatPredEstimator()
        good_key = ("SEQ0", (("X", "C"),), None, False)
        est._memo[good_key] = EstimatedKinetics(km=1.0, kcat=1.0, source="catpred")
        kwargs = dict(enzyme_sequence="SEQ", reactant_smiles={"S": "CCO"})
        first = est.estimate(**kwargs)
        assert "errors" in first.source
        n_calls = len(calls)
        assert est.estimate(**kwargs).source == first.source
        assert len(calls) == n_calls  # reused in-process, no second subprocess

        est._save_persistent_cache()
        est2 = CatPredEstimator()
        assert good_key in est2._memo
        assert all("errors" not in v.source for v in est2._memo.values())


class TestEstimateKi:
    """CatPred Ki for inhibitors whose YAML omits ki (the CatPred call itself is mocked)."""

    def _estimator(self, monkeypatch, tmp_path):
        from bees.kinetics_estimator import CatPredEstimator
        monkeypatch.setenv("BEES_CATPRED_CACHE", "off")
        monkeypatch.chdir(tmp_path)
        (tmp_path / "catpred").mkdir()
        monkeypatch.setattr(CatPredEstimator, "CATPRED_DIR", str(tmp_path / "catpred"))
        return CatPredEstimator(include_sd=True)

    def _fake_run(self, value, calls):
        import os

        def _run(cmd, cwd=None, **kwargs):
            calls.append(cmd)
            csv_path = cmd[cmd.index("--input_file") + 1]
            base = os.path.splitext(os.path.basename(csv_path))[0]
            out = os.path.join(cwd, "output", base)
            os.makedirs(out)
            with open(os.path.join(out, f"{base}_ki_output.csv"), "w") as fh:
                fh.write(f"SMILES,Prediction_(mM),SD_total\nCCN,{value},0.3\n")
        return _run

    def test_returns_the_predicted_ki_and_memoizes(self, monkeypatch, tmp_path):
        import subprocess
        est = self._estimator(monkeypatch, tmp_path)
        calls = []
        monkeypatch.setattr(subprocess, "run", self._fake_run(0.021, calls))
        first = est.estimate_ki(enzyme_sequence="MKT", inhibitor_smiles="CCN")
        assert first.ki == pytest.approx(0.021) and first.ki_sd > 0 and first.source == "catpred"
        assert est.estimate_ki(enzyme_sequence="MKT", inhibitor_smiles="CCN") is first
        assert len(calls) == 1
        assert "--parameter" in calls[0] and calls[0][calls[0].index("--parameter") + 1] == "ki"

    def test_failed_run_reports_the_error_instead_of_a_value(self, monkeypatch, tmp_path):
        import subprocess
        est = self._estimator(monkeypatch, tmp_path)

        def _fail(cmd, **kwargs):
            raise subprocess.CalledProcessError(1, cmd, stderr="boom")

        monkeypatch.setattr(subprocess, "run", _fail)
        res = est.estimate_ki(enzyme_sequence="MKT", inhibitor_smiles="CCN")
        assert res.ki is None and "boom" in res.source

    def test_needs_sequence_and_smiles(self, monkeypatch, tmp_path):
        est = self._estimator(monkeypatch, tmp_path)
        assert est.estimate_ki(enzyme_sequence="", inhibitor_smiles="CCN").ki is None
        assert est.estimate_ki(enzyme_sequence="MKT", inhibitor_smiles="").ki is None

    def test_missing_catpred_installation_is_an_error(self, monkeypatch, tmp_path):
        from bees.kinetics_estimator import CatPredEstimator
        monkeypatch.setenv("BEES_CATPRED_CACHE", "off")
        monkeypatch.setattr(CatPredEstimator, "CATPRED_DIR", str(tmp_path / "absent"))
        with pytest.raises(FileNotFoundError, match="give 'ki'"):
            CatPredEstimator().estimate_ki(enzyme_sequence="MKT", inhibitor_smiles="CCN")
