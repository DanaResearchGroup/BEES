"""
Tests for the three robustness fixes:
1. ODE failures propagate to execute() and CLI
2. Missing thermodynamics state is explicit
3. Missing database is explicit

To run: pytest -v tests/test_robustness_fixes.py
"""

import os
import pytest
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch
from types import SimpleNamespace

import bees.main as main
from bees.enlarger import IterativeEnlarger, EnlargerResult
from bees.simulator import SimulationResult
from bees.thermodynamics import ThermoEngine, _DiskCache


# ---------------------------------------------------------------------------
# Issue 1: ODE failure propagation
# ---------------------------------------------------------------------------

class TestODEFailurePropagation:
    """ODE failures must propagate to BEES.execute() and return success=False."""
    
    @pytest.fixture
    def mock_bees_object(self):
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
            SimpleNamespace(
                label="Enzyme", concentration=0.01, reactive=True, solvent=False,
                ecnumber="EC 1.1.1.1", constant=True
            ),
        ]
        bees_obj = SimpleNamespace(
            settings=settings,
            species=species,
            enzymes=enzymes,
            project="TestProject",
            database=SimpleNamespace(name="db", missing_policy="fallback_with_warning"),
            environment=SimpleNamespace(temperature=298.15, pH=7.0),
        )
        return bees_obj
    
    def test_enlarger_result_tracks_ode_failure(self, mock_bees_object):
        """EnlargerResult.ode_failed is set when simulation fails."""
        with tempfile.TemporaryDirectory() as output_dir:
            logger = MagicMock()
            mock_reaction_generator = MagicMock()
            mock_reaction_generator._generate_reactions.return_value = []
            mock_reaction_generator.ensure_estimator_initialized = MagicMock()
            
            enlarger = IterativeEnlarger(
                bees_object=mock_bees_object,
                reaction_generator=mock_reaction_generator,
                logger=logger,
                output_directory=output_dir,
            )
            
            failed_sim = SimulationResult(
                t=[],
                y=[],
                species_labels=[],
                success=False,
                message="ODE solver diverged",
            )
            
            with patch("bees.enlarger.ODESimulator") as mock_sim_cls:
                mock_sim_inst = MagicMock()
                mock_sim_inst.simulate.return_value = failed_sim
                mock_sim_cls.return_value = mock_sim_inst
                
                result = enlarger.run()
                
                assert result.ode_failed is True
                assert "ODE failure" in result.convergence_reason
                assert not result.converged
    
    def test_main_execute_returns_failure_on_ode_failure(self, tmp_path, monkeypatch):
        """BEES.execute() returns success=False when ODE fails."""
        # Setup paths
        bees_path = tmp_path / "bees_path"
        (bees_path / "db").mkdir(parents=True)
        (bees_path / "db" / "db.csv").write_text("ecnumber,substrate,vmax,km\n")
        (bees_path / "projects" / "test").mkdir(parents=True)
        
        monkeypatch.setattr(main, "BEES_PATH", str(bees_path))
        monkeypatch.setattr(main.common, "BEES_PATH", str(bees_path))
        monkeypatch.setattr(main.common, "get_git_branch", lambda: "test")
        monkeypatch.setattr(main.common, "get_git_commit", lambda: ("abc", "2025-01-01"))
        monkeypatch.setattr(main.common, "load_and_invert_ontology", lambda p: {})
        
        input_data = {
            "project": "test_ode_failure",
            "species": [{"label": "S", "concentration": 1.0, "reactive": True}],
            "enzymes": [{"label": "E", "concentration": 0.01, "ecnumber": "EC 1.1.1.1", "reactive": True}],
            "environment": {"temperature": 300.0, "pH": 7.0},
            "database": {"name": "db", "missing_policy": "fallback_with_warning"},
            "settings": {"end_time": 10.0, "toleranceMoveToCore": 1e-5},
        }
        
        bees_instance = main.BEES(input_data)
        
        # Mock the enlarger to return an ODE failure
        failed_result = EnlargerResult(
            iterations=1,
            converged=False,
            convergence_reason="ODE failure: solver crashed",
            ode_failed=True,
        )
        
        with patch("bees.main.IterativeEnlarger") as mock_enlarger_cls:
            mock_enlarger = MagicMock()
            mock_enlarger.run.return_value = failed_result
            mock_enlarger._reaction_id_by_sig = {}
            mock_enlarger_cls.return_value = mock_enlarger
            
            with patch("bees.main.ReactionGenerator"):
                result = bees_instance.execute()
        
        assert result["success"] is False
        assert "ODE" in result["convergence_reason"]
        assert "database_path" in result


# ---------------------------------------------------------------------------
# Issue 2: Missing thermodynamics explicit state
# ---------------------------------------------------------------------------

class TestThermodynamicsPolicy:
    """Thermodynamics fallback policy must be configurable and explicit."""
    
    def test_thermo_engine_respects_error_policy(self, tmp_path):
        """ThermoEngine with policy='error' raises on missing data."""
        engine = ThermoEngine(
            pH=7.4,
            cache=_DiskCache(tmp_path),
            missing_data_policy="error",
        )
        # Mock CC to be unavailable
        with patch.object(engine, "_load_cc", return_value=None):
            with pytest.raises(ValueError, match="Thermodynamics computation failed"):
                engine.compute_keq(
                    stoichiometry={"S": -1, "P": 1},
                    smiles_map={"S": "C", "P": "C"},
                )
    
    def test_thermo_engine_warn_policy_logs_warning(self, tmp_path):
        """ThermoEngine with policy='warn_assume_irreversible' logs warning."""
        engine = ThermoEngine(
            pH=7.4,
            cache=_DiskCache(tmp_path),
            missing_data_policy="warn_assume_irreversible",
        )
        
        with patch.object(engine, "_load_cc", return_value=None):
            with patch("bees.thermodynamics.logger") as mock_logger:
                td = engine.compute_keq(
                    stoichiometry={"S": -1, "P": 1},
                    smiles_map={"S": "C", "P": "C"},
                )
                
                assert td.source == "fallback"
                assert td.irreversible is True
                # Should have called warning (not just info/debug)
                assert any(
                    call[0][0].startswith("Thermo fallback")
                    for call in mock_logger.warning.call_args_list
                )
    
    def test_thermo_engine_silent_policy_no_warning(self, tmp_path):
        """ThermoEngine with policy='silent_assume_irreversible' logs only debug."""
        engine = ThermoEngine(
            pH=7.4,
            cache=_DiskCache(tmp_path),
            missing_data_policy="silent_assume_irreversible",
        )
        
        with patch.object(engine, "_load_cc", return_value=None):
            with patch("bees.thermodynamics.logger") as mock_logger:
                td = engine.compute_keq(
                    stoichiometry={"S": -1, "P": 1},
                    smiles_map={"S": "C", "P": "C"},
                )
                
                assert td.source == "fallback"
                # Should NOT have called warning
                mock_logger.warning.assert_not_called()
    
    def test_enlarger_passes_policy_to_thermo_engine(self):
        """IterativeEnlarger passes thermo_missing_data_policy to ThermoEngine."""
        with tempfile.TemporaryDirectory() as output_dir:
            settings = SimpleNamespace(
                end_time=100.0,
                toleranceMoveToCore=1e-5,
                max_iterations=1,
                thermo_missing_data_policy="error",
                thermo_irreversible_cutoff_kJmol=30.0,
            )
            bees_obj = SimpleNamespace(
                settings=settings,
                species=[SimpleNamespace(label="S", concentration=1.0, reactive=True, constant=False)],
                enzymes=[SimpleNamespace(label="E", concentration=0.01, reactive=True, constant=True, ecnumber="EC 1.1.1.1")],
                environment=SimpleNamespace(temperature=300.0, pH=7.0),
            )
            
            mock_reaction_generator = MagicMock()
            enlarger = IterativeEnlarger(
                bees_object=bees_obj,
                reaction_generator=mock_reaction_generator,
                logger=MagicMock(),
                output_directory=output_dir,
            )
            
            engine = enlarger._build_thermo_engine()
            assert engine.missing_data_policy == "error"


# ---------------------------------------------------------------------------
# Issue 3: Missing database explicit state
# ---------------------------------------------------------------------------

class TestDatabaseMissingPolicy:
    """Database resolution must respect missing_policy and record resolved path."""
    
    def test_resolve_db_path_error_policy_raises(self, tmp_path, monkeypatch):
        """resolve_db_path with policy='error' raises when DB is missing."""
        monkeypatch.setattr(main, "BEES_PATH", str(tmp_path))
        (tmp_path / "db").mkdir()
        
        logger = MagicMock()
        
        with pytest.raises(FileNotFoundError, match="not found"):
            main.resolve_db_path("nonexistent", logger, missing_policy="error")
        
        # Should not have logged a warning (it raised instead)
        logger.warning.assert_not_called()
    
    def test_resolve_db_path_fallback_policy_warns(self, tmp_path, monkeypatch):
        """resolve_db_path with policy='fallback_with_warning' logs warning."""
        monkeypatch.setattr(main, "BEES_PATH", str(tmp_path))
        (tmp_path / "db").mkdir()
        (tmp_path / "db" / "db.csv").write_text("")
        
        logger = MagicMock()
        
        resolved = main.resolve_db_path(
            "nonexistent",
            logger,
            missing_policy="fallback_with_warning"
        )
        
        assert resolved == str(tmp_path / "db" / "db.csv")
        # Should have logged a warning
        assert logger.warning.call_count == 1
        assert "nonexistent" in logger.warning.call_args[0][0]
    
    def test_resolve_db_path_succeeds_when_db_exists(self, tmp_path, monkeypatch):
        """resolve_db_path returns correct path when DB exists (both policies)."""
        monkeypatch.setattr(main, "BEES_PATH", str(tmp_path))
        (tmp_path / "db").mkdir()
        (tmp_path / "db" / "ecoli.csv").write_text("")
        
        logger = MagicMock()
        
        # With error policy
        resolved = main.resolve_db_path("ecoli", logger, missing_policy="error")
        assert resolved == str(tmp_path / "db" / "ecoli.csv")
        logger.warning.assert_not_called()
        
        # With fallback policy
        resolved = main.resolve_db_path("ecoli", logger, missing_policy="fallback_with_warning")
        assert resolved == str(tmp_path / "db" / "ecoli.csv")
        logger.warning.assert_not_called()
    
    def test_execute_records_database_path_in_results(self, tmp_path, monkeypatch):
        """BEES.execute() includes database_path in results."""
        bees_path = tmp_path / "bees_path"
        (bees_path / "db").mkdir(parents=True)
        (bees_path / "db" / "testdb.csv").write_text("ecnumber,substrate,vmax,km\n")
        (bees_path / "projects" / "test").mkdir(parents=True)
        
        monkeypatch.setattr(main, "BEES_PATH", str(bees_path))
        monkeypatch.setattr(main.common, "BEES_PATH", str(bees_path))
        monkeypatch.setattr(main.common, "get_git_branch", lambda: "test")
        monkeypatch.setattr(main.common, "get_git_commit", lambda: ("abc", "2025-01-01"))
        monkeypatch.setattr(main.common, "load_and_invert_ontology", lambda p: {})
        
        input_data = {
            "project": "test_db_record",
            "species": [{"label": "S", "concentration": 1.0, "reactive": True}],
            "enzymes": [{"label": "E", "concentration": 0.01, "ecnumber": "EC 1.1.1.1", "reactive": True}],
            "environment": {"temperature": 300.0, "pH": 7.0},
            "database": {"name": "testdb", "missing_policy": "fallback_with_warning"},
            "settings": {"simulation_mode": "batch"},
        }
        
        bees_instance = main.BEES(input_data)
        
        with patch("bees.main.ReactionGenerator") as mock_rg_cls:
            mock_rg = MagicMock()
            mock_rg.generate_reactions.return_value = []
            mock_rg.export_reactions_summary.return_value = None
            mock_rg_cls.return_value = mock_rg
            
            result = bees_instance.execute()
        
        assert result["success"] is True
        assert "database_path" in result
        assert result["database_path"] == str(bees_path / "db" / "testdb.csv")


# ---------------------------------------------------------------------------
# Integration tests
# ---------------------------------------------------------------------------

class TestRobustnessIntegration:
    """End-to-end tests combining multiple fixes."""
    
    def test_missing_db_with_error_policy_fails_before_enlargement(self, tmp_path, monkeypatch):
        """When database.missing_policy='error' and DB is missing, fail early."""
        bees_path = tmp_path / "bees_path"
        (bees_path / "db").mkdir(parents=True)
        (bees_path / "projects" / "test").mkdir(parents=True)
        
        monkeypatch.setattr(main, "BEES_PATH", str(bees_path))
        monkeypatch.setattr(main.common, "BEES_PATH", str(bees_path))
        monkeypatch.setattr(main.common, "get_git_branch", lambda: "test")
        monkeypatch.setattr(main.common, "get_git_commit", lambda: ("abc", "2025-01-01"))
        monkeypatch.setattr(main.common, "load_and_invert_ontology", lambda p: {})
        
        input_data = {
            "project": "test_fail_early",
            "species": [{"label": "S", "concentration": 1.0, "reactive": True}],
            "enzymes": [{"label": "E", "concentration": 0.01, "ecnumber": "EC 1.1.1.1", "reactive": True}],
            "environment": {"temperature": 300.0, "pH": 7.0},
            "database": {"name": "nonexistent", "missing_policy": "error"},
            "settings": {"end_time": 10.0, "toleranceMoveToCore": 1e-5},
        }
        
        bees_instance = main.BEES(input_data)
        
        with pytest.raises(FileNotFoundError, match="not found"):
            bees_instance.execute()
