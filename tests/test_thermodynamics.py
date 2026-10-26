"""
Tests for bees.thermodynamics.

To run:  pytest -v tests/test_thermodynamics.py

All tests mock equilibrator-api so no live network/database calls are made.
Tests that exercise the real (non-disabled) thermo path are skipped when
BEES_DISABLE_THERMO=1 is set. 
"""

import math
import os
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from bees.common import R

from bees.thermodynamics import (
    CompoundSubstitutor,
    ThermoEngine,
    _DiskCache,
    _is_disabled,
)
from bees.common import smiles_to_inchi as _smiles_to_inchi, smiles_to_inchikey as _smiles_to_inchikey


# Tests that need the real (non-disabled) engine are skipped when the user
# has set BEES_DISABLE_THERMO=1.
_skip_if_disabled = pytest.mark.skipif(
    _is_disabled(),
    reason="BEES_DISABLE_THERMO=1 forces the irreversible fallback path",
)

# ---------------------------------------------------------------------------
# SMILES → InChI / InChIKey helpers
# ---------------------------------------------------------------------------

class TestSmilesConversion:
    def test_smiles_to_inchi_valid(self):
        inchi = _smiles_to_inchi("C")  # methane
        assert inchi is not None
        assert inchi.startswith("InChI=")

    def test_smiles_to_inchi_invalid(self):
        assert _smiles_to_inchi("NOTASMILES!!!") is None
        assert _smiles_to_inchi(None) is None
        assert _smiles_to_inchi("") is None

    def test_smiles_to_inchikey_valid(self):
        ikey = _smiles_to_inchikey("C")
        assert ikey is not None
        assert len(ikey) == 27  # standard InChIKey length XX-XX-X

    def test_smiles_to_inchikey_invalid(self):
        assert _smiles_to_inchikey("NOTASMILES!!!") is None
        assert _smiles_to_inchikey(None) is None


# ---------------------------------------------------------------------------
# Disk cache
# ---------------------------------------------------------------------------

class TestDiskCache:
    def test_roundtrip(self, tmp_path: Path):
        cache = _DiskCache(tmp_path)
        cache.set("ns", "key", {"a": 1})
        assert cache.get("ns", "key") == {"a": 1}

    def test_miss_returns_default(self, tmp_path: Path):
        cache = _DiskCache(tmp_path)
        sentinel = object()
        assert cache.get("ns", "absent", default=sentinel) is sentinel

    def test_persists_none_value(self, tmp_path: Path):
        # Negative caching: None should be stored and distinguished from a miss.
        cache = _DiskCache(tmp_path)
        cache.set("ns", "key", None)
        assert cache.get("ns", "key", default=...) is None

    def test_different_namespaces_independent(self, tmp_path: Path):
        cache = _DiskCache(tmp_path)
        cache.set("ns1", "key", "value1")
        cache.set("ns2", "key", "value2")
        assert cache.get("ns1", "key") == "value1"
        assert cache.get("ns2", "key") == "value2"


# ---------------------------------------------------------------------------
# ThermoEngine — disabled path
# ---------------------------------------------------------------------------

class TestThermoEngineDisabled:
    def test_disabled_via_env(self, tmp_path: Path):
        engine = ThermoEngine(pH=7.4, cache=_DiskCache(tmp_path))
        with patch.dict(os.environ, {"BEES_DISABLE_THERMO": "1"}):
            td = engine.compute_keq(
                stoichiometry={"S": -1, "P": 1},
                smiles_map={"S": "C", "P": "C"},
            )
        assert td.source == "disabled"
        assert td.irreversible is True
        assert td.kcat_rev is None

    @_skip_if_disabled
    def test_fallback_when_cc_unavailable(self, tmp_path: Path):
        engine = ThermoEngine(pH=7.4, cache=_DiskCache(tmp_path))
        with patch.object(engine, "_load_cc", return_value=None):
            td = engine.compute_keq(
                stoichiometry={"S": -1, "P": 1},
                smiles_map={"S": "C", "P": "C"},
            )
        assert td.source == "fallback"
        assert td.irreversible is True


# ---------------------------------------------------------------------------
# ThermoEngine — core logic with mocked equilibrator
# ---------------------------------------------------------------------------

@_skip_if_disabled
class TestThermoEngineWithMockedCC:
    """Mock only _compute_dgr_prime so the Keq/cache logic is tested
    against real code, without any equilibrator-api calls."""

    def _engine(self, tmp_path, dgr_kJmol, sigma_kJmol=1.0):
        engine = ThermoEngine(pH=7.4, T_K=310.15, cache=_DiskCache(tmp_path))
        engine._cc = MagicMock()  # mark CC as loaded
        engine._compute_dgr_prime = MagicMock(return_value=(dgr_kJmol, sigma_kJmol))
        return engine

    def _kwargs(self):
        return dict(
            stoichiometry={"S": -1, "P": 1},
            smiles_map={"S": "C", "P": "C"},
        )

    def test_keq_math(self, tmp_path: Path):
        # Keq = exp(-ΔG°' / RT) at 310.15 K
        engine = self._engine(tmp_path, dgr_kJmol=-10.0)
        td = engine.compute_keq(**self._kwargs())

        assert td.source == "equilibrator"
        assert td.irreversible is False
        assert td.dgr_prime_kJmol == pytest.approx(-10.0)
        assert td.kcat_rev is None  # Haldane is caller's responsibility

        expected_keq = math.exp(10.0 * 1000 / (R * 310.15))
        assert td.keq == pytest.approx(expected_keq, rel=1e-9)

    @pytest.mark.parametrize("dgr", [-50.0, 45.0])
    def test_compute_keq_returns_raw_irreversibility_false(self, tmp_path: Path, dgr):
        # compute_keq returns raw thermo with irreversible=False.
        # The ΔG°' < −cutoff rule has moved to DgrIrreversibility in
        # bees.rules.physics_rules and runs via RULES.apply_all.
        engine = self._engine(tmp_path, dgr_kJmol=dgr)
        td = engine.compute_keq(**self._kwargs())
        assert td.irreversible is False
        assert td.kcat_rev is None
        assert td.dgr_prime_kJmol == pytest.approx(dgr)

    def test_caches_reaction_on_second_call(self, tmp_path: Path):
        engine = self._engine(tmp_path, dgr_kJmol=-10.0)
        engine.compute_keq(**self._kwargs())
        engine.compute_keq(**self._kwargs())
        # _compute_dgr_prime called only once; second call served from cache.
        assert engine._compute_dgr_prime.call_count == 1

    def test_cache_keyed_by_substituted_chemistry(self, tmp_path: Path):
        """Same raw SMILES, different substitutions: the second engine must not reuse
        the first engine's ΔG, which was computed for different effective compounds."""

        class _Sub(CompoundSubstitutor):
            def __init__(self, smiles):
                self.smiles = smiles

            def substitute(self, label, smiles):
                return self.smiles if label == "S" else None

        first = self._engine(tmp_path, dgr_kJmol=-10.0)
        first.substitutor = _Sub("CC")
        second = self._engine(tmp_path, dgr_kJmol=5.0)  # same cache directory
        second.substitutor = _Sub("CCC")

        assert first.compute_keq(**self._kwargs()).dgr_prime_kJmol == pytest.approx(-10.0)
        td = second.compute_keq(**self._kwargs())
        assert second._compute_dgr_prime.call_count == 1
        assert td.dgr_prime_kJmol == pytest.approx(5.0)

        # An identical substitution still hits the cache.
        third = self._engine(tmp_path, dgr_kJmol=99.0)
        third.substitutor = _Sub("CCC")
        assert third.compute_keq(**self._kwargs()).dgr_prime_kJmol == pytest.approx(5.0)
        assert third._compute_dgr_prime.call_count == 0

    def test_fallback_when_dgr_is_none(self, tmp_path: Path):
        engine = ThermoEngine(pH=7.4, T_K=310.15, cache=_DiskCache(tmp_path))
        engine._cc = MagicMock()
        engine._compute_dgr_prime = MagicMock(return_value=(None, None))
        td = engine.compute_keq(**self._kwargs())
        assert td.source == "fallback"
        assert td.irreversible is True


# ---------------------------------------------------------------------------
# Compound resolution (mocked CC)
# ---------------------------------------------------------------------------

@_skip_if_disabled
class TestCompoundResolution:
    def test_inchi_lookup_preferred(self, tmp_path: Path):
        """Primary path: get_compound_by_inchi() is called first."""
        engine = ThermoEngine(pH=7.4, cache=_DiskCache(tmp_path))
        mock_cc = MagicMock()
        mock_cc.get_compound_by_inchi.return_value = MagicMock()
        engine._cc = mock_cc

        cpd = engine._resolve_compound("C")  # methane — valid SMILES
        mock_cc.get_compound_by_inchi.assert_called_once()
        assert cpd is not None

    def test_inchikey_fallback_on_inchi_miss(self, tmp_path: Path):
        """If InChI lookup returns None, fall back to InChIKey prefix search."""
        engine = ThermoEngine(pH=7.4, cache=_DiskCache(tmp_path))
        mock_cc = MagicMock()
        mock_cc.get_compound_by_inchi.return_value = None
        mock_compound = MagicMock()
        mock_compound.inchi_key = _smiles_to_inchikey("C")
        mock_cc.search_compound_by_inchi_key.return_value = [mock_compound]
        engine._cc = mock_cc

        cpd = engine._resolve_compound("C")
        mock_cc.search_compound_by_inchi_key.assert_called_once()
        assert cpd is mock_compound

    def test_none_returned_for_invalid_smiles(self, tmp_path: Path):
        engine = ThermoEngine(pH=7.4, cache=_DiskCache(tmp_path))
        mock_cc = MagicMock()
        mock_cc.get_compound_by_inchi.return_value = None
        mock_cc.search_compound_by_inchi_key.return_value = []
        mock_cc.get_compound.return_value = None
        engine._cc = mock_cc
        cpd = engine._resolve_compound("NOTASMILES!!!")
        assert cpd is None
        mock_cc.get_compound_by_inchi.assert_not_called()
        mock_cc.search_compound_by_inchi_key.assert_not_called()

    def test_cache_miss_is_persisted(self, tmp_path: Path):
        """Negative cache: after a miss, second lookup short-circuits to None.

        Hits are deliberately NOT cached on disk because equilibrator-cache's
        Compound object is bound to a SQLAlchemy session and cannot be safely
        pickled/restored. Only misses are persisted.
        """
        engine = ThermoEngine(pH=7.4, cache=_DiskCache(tmp_path))
        mock_cc = MagicMock()
        mock_cc.get_compound_by_inchi.return_value = None
        mock_cc.search_compound_by_inchi_key.return_value = []
        mock_cc.get_compound.return_value = None
        engine._cc = mock_cc

        first = engine._resolve_compound("C")
        second = engine._resolve_compound("C")
        assert first is None
        assert second is None
        # Second call returns the persisted None without re-hitting CC.
        assert mock_cc.get_compound_by_inchi.call_count == 1

    def test_lookup_error_is_not_persisted(self, tmp_path: Path):
        """A DB error is possibly transient: retry next time, don't cache a miss."""
        engine = ThermoEngine(pH=7.4, cache=_DiskCache(tmp_path))
        mock_cc = MagicMock()
        mock_cc.get_compound_by_inchi.side_effect = RuntimeError("db locked")
        engine._cc = mock_cc

        assert engine._resolve_compound("C") is None
        mock_cc.get_compound_by_inchi.side_effect = None
        mock_cc.get_compound_by_inchi.return_value = MagicMock()
        assert engine._resolve_compound("C") is not None

    def test_ambiguous_inchikey_uses_first_match(self, tmp_path: Path):
        """Multiple InChIKey hits with no exact match → fall back to first.
        """
        engine = ThermoEngine(pH=7.4, cache=_DiskCache(tmp_path))
        mock_cc = MagicMock()
        mock_cc.get_compound_by_inchi.return_value = None
        m1 = MagicMock()
        m1.inchi_key = "AAAAAAAAAAAAAAA-AAAAAAAAAA-A"
        m2 = MagicMock()
        m2.inchi_key = "BBBBBBBBBBBBBBB-BBBBBBBBBB-B"
        mock_cc.search_compound_by_inchi_key.return_value = [m1, m2]
        engine._cc = mock_cc

        cpd = engine._resolve_compound("C")
        assert cpd is m1


@pytest.mark.parametrize("value, disabled", [
    ("", False), ("0", False), ("false", False), ("off", False), ("OFF", False), ("no", False),
    ("1", True), ("true", True), ("yes", True),
])
def test_disable_thermo_env_parsing(monkeypatch, value, disabled):
    monkeypatch.setenv("BEES_DISABLE_THERMO", value)
    assert _is_disabled() is disabled
