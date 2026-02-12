"""
Kinetics estimation module.

This module defines a small abstraction layer so BEES can fill missing kinetic
parameters (kcat, Km, Ki) using external tools (e.g., CatPred) without coupling
the core generator to a specific backend.

Current status:
- Provides an interface integrated with CatPred for kcat and Km estimation.
"""

from __future__ import annotations

import logging
import os
import subprocess
import uuid
from dataclasses import dataclass
from typing import Dict, Optional

import pandas as pd

from bees.common import log10_sd_to_linear_sd


@dataclass(frozen=True)
class EstimatedKinetics:
    """
    Estimated kinetic parameters.

    Units:
    - km: mM (single value, for backward compat when only one substrate)
    - km_per_substrate: mM per substrate (CatPred per-substrate Km)
    - ki: mM
    - kcat: 1/s (one per reaction)
    - *_sd: standard deviation of the prediction (same units as the parameter)
    """

    km: Optional[float] = None
    km_per_substrate: Optional[Dict[str, float]] = None
    km_sd: Optional[float] = None
    km_sd_per_substrate: Optional[Dict[str, float]] = None
    kcat: Optional[float] = None
    ki: Optional[float] = None
    kcat_sd: Optional[float] = None
    ki_sd: Optional[float] = None
    source: str = "estimator"


class BaseKineticsEstimator:
    """Adaptor - base class for integrate with external kinetics estimators."""

    name: str = "base"

    def estimate(
        self,
        *,
        enzyme_sequence: str,
        reactant_smiles: Dict[str, str],
        inhibitor_smiles: Optional[str] = None,
    ) -> EstimatedKinetics:
        """
        Estimate kinetics from enzyme and reactant SMILES.

        Args:
            enzyme_sequence: Enzyme amino acid sequence
            reactant_smiles: Map of compound name -> SMILES for all reactants.
                Used for: kcat = concatenated SMILES; Km = one per substrate.
            inhibitor_smiles: Optional inhibitor SMILES for Ki
        """
        raise NotImplementedError


class CatPredEstimator(BaseKineticsEstimator):
    """
    CatPred kinetics estimator.

    Uses a CLI wrapper to call CatPred in its own conda environment.
    This avoids dependency conflicts between BEES and CatPred.
    """

    name = "catpred"

    def __init__(self, include_sd: bool = False):
        """
        Args:
            include_sd: If True, include SD_total (standard deviation) from CatPred output for each parameter.
        """
        self.include_sd = include_sd

    # Paths to CatPred installation and models.
    # Set env vars (CATPRED_DIR, CATPRED_CHECKPOINT_BASE, CATPRED_CONDA_ENV) or
    # edit these defaults to match your installation.
    CATPRED_DIR = os.environ.get("CATPRED_DIR", "/path/to/CatPred")
    CHECKPOINT_BASE = os.environ.get(
        "CATPRED_CHECKPOINT_BASE",
        "/path/to/pretrained/production"
    )
    CONDA_ENV = os.environ.get("CATPRED_CONDA_ENV", "catpred")

    def estimate(
        self,
        *,
        enzyme_sequence: str,
        reactant_smiles: Dict[str, str],
        inhibitor_smiles: Optional[str] = None,
    ) -> EstimatedKinetics:
        """
        Estimate kcat and Km using CatPred.

        - kcat: concatenated SMILES of all reactants (one kcat per reaction)
        - Km: one value per substrate (individual substrate SMILES)
        """
        if not reactant_smiles:
            return EstimatedKinetics(source="catpred(no reactants)")

        # 1. Create a unique ID for this prediction request to avoid file collisions
        run_id = f"bees_{uuid.uuid4().hex[:8]}"
        
        # Ensure a local writable temp directory exists
        local_tmp = os.path.join(os.getcwd(), ".tmp")
        os.makedirs(local_tmp, exist_ok=True)
        
        # Create a writable working directory for CatPred execution
        catpred_work_dir = os.path.join(local_tmp, f"work_{run_id}")
        os.makedirs(catpred_work_dir, exist_ok=True)
        
        # Symlink necessary CatPred files to the writable work dir
        for item in os.listdir(self.CATPRED_DIR):
            if item.startswith('.') or item == "output" or item == "demo":
                continue
            src = os.path.normpath(os.path.join(self.CATPRED_DIR, item))
            dst = os.path.join(catpred_work_dir, item)
            if not os.path.exists(dst):
                try:
                    os.symlink(src, dst)
                except Exception as exc:
                    # Best-effort symlinking: continue if a symlink cannot be created,
                    # but emit a warning so environment or permission issues are visible.
                    print(f"Warning: failed to create symlink from {src} to {dst}: {exc}")
        
        os.makedirs(os.path.join(catpred_work_dir, "output"), exist_ok=True)
        os.makedirs(os.path.join(catpred_work_dir, "demo"), exist_ok=True)

        results = {}

        # --- kcat: one row with concatenated SMILES of all reactants ---
        concatenated_smiles = ".".join(s for s in reactant_smiles.values() if s)
        if concatenated_smiles:
            kcat_csv = os.path.join(catpred_work_dir, f"{run_id}_kcat.csv")
            df_kcat = pd.DataFrame([{
                "SMILES": concatenated_smiles,
                "sequence": enzyme_sequence,
                "pdbpath": "bees_query"
            }])
            try:
                df_kcat.to_csv(kcat_csv, index=False)
            except Exception as e:
                return EstimatedKinetics(source=f"catpred(error writing kcat input: {str(e)})")

            param = "kcat"
            checkpoint_dir = os.path.join(self.CHECKPOINT_BASE, param)
            cmd = [
                "conda", "run", "-n", self.CONDA_ENV,
                "python", "demo_run.py",
                "--parameter", param,
                "--input_file", kcat_csv,
                "--checkpoint_dir", checkpoint_dir
            ]
            try:
                env = os.environ.copy()
                env["TMPDIR"] = local_tmp
                env["TEMP"] = local_tmp
                env["TMP"] = local_tmp
                subprocess.run(cmd, cwd=catpred_work_dir, check=True, capture_output=True, text=True, env=env)
                kcat_basename = os.path.splitext(os.path.basename(kcat_csv))[0]
                output_path = os.path.join(catpred_work_dir, "output", kcat_basename, f"{kcat_basename}_{param}_output.csv")
                if os.path.exists(output_path):
                    df_res = pd.read_csv(output_path)
                    if not df_res.empty:
                        col_name = "Prediction_(s^(-1))"
                        if col_name in df_res.columns:
                            results["kcat"] = float(df_res[col_name].iloc[0])
                        if self.include_sd and "SD_total" in df_res.columns:
                            pred_val = results.get("kcat")
                            sd_log10 = float(df_res["SD_total"].iloc[0])
                            if pred_val is not None and pred_val > 0 and sd_log10 > 0:
                                results["kcat_sd"] = log10_sd_to_linear_sd(pred_val, sd_log10)
            except subprocess.CalledProcessError as e:
                results["kcat_error"] = e.stderr.strip() if e.stderr else str(e)
            except Exception as e:
                results["kcat_error"] = str(e)

        # --- km: one row per reactant (individual substrate SMILES) ---
        reactant_names = list(reactant_smiles.keys())
        reactant_smiles_list = [reactant_smiles.get(n, "") for n in reactant_names]
        if any(reactant_smiles_list):
            km_csv = os.path.join(catpred_work_dir, f"{run_id}_km.csv")
            df_km = pd.DataFrame([
                {"SMILES": smi, "sequence": enzyme_sequence, "pdbpath": "bees_query"}
                for smi in reactant_smiles_list if smi
            ])
            # Track which reactant name corresponds to each row (only rows with SMILES)
            km_reactant_order = [n for n, smi in zip(reactant_names, reactant_smiles_list) if smi]
            try:
                df_km.to_csv(km_csv, index=False)
            except Exception as e:
                pass  # results already populated for kcat
            else:
                param = "km"
                checkpoint_dir = os.path.join(self.CHECKPOINT_BASE, param)
                cmd = [
                    "conda", "run", "-n", self.CONDA_ENV,
                    "python", "demo_run.py",
                    "--parameter", param,
                    "--input_file", km_csv,
                    "--checkpoint_dir", checkpoint_dir
                ]
                try:
                    env = os.environ.copy()
                    env["TMPDIR"] = local_tmp
                    env["TEMP"] = local_tmp
                    env["TMP"] = local_tmp
                    subprocess.run(cmd, cwd=catpred_work_dir, check=True, capture_output=True, text=True, env=env)
                    km_basename = os.path.splitext(os.path.basename(km_csv))[0]
                    output_path = os.path.join(catpred_work_dir, "output", km_basename, f"{km_basename}_{param}_output.csv")
                    if os.path.exists(output_path):
                        df_res = pd.read_csv(output_path)
                        col_name = "Prediction_(mM)"
                        if not df_res.empty and col_name in df_res.columns:
                            km_per_substrate = {}
                            km_sd_per_substrate = {}
                            for i, rname in enumerate(km_reactant_order):
                                if i < len(df_res):
                                    val = float(df_res[col_name].iloc[i])
                                    km_per_substrate[rname] = val
                                    if self.include_sd and "SD_total" in df_res.columns and i < len(df_res):
                                        sd_log10 = float(df_res["SD_total"].iloc[i])
                                        if val > 0 and sd_log10 > 0:
                                            km_sd_per_substrate[rname] = log10_sd_to_linear_sd(val, sd_log10)
                            results["km_per_substrate"] = km_per_substrate
                            if km_sd_per_substrate:
                                results["km_sd_per_substrate"] = km_sd_per_substrate
                            # Backward compat: single km = first substrate
                            if km_per_substrate:
                                first = next(iter(km_per_substrate.values()))
                                results["km"] = first
                                if km_sd_per_substrate:
                                    results["km_sd"] = next(iter(km_sd_per_substrate.values()))
                except subprocess.CalledProcessError as e:
                    results["km_error"] = e.stderr.strip() if e.stderr else str(e)
                except Exception as e:
                    results["km_error"] = str(e)

        # 3. Cleanup
        import shutil
        if os.path.exists(catpred_work_dir):
            try:
                shutil.rmtree(catpred_work_dir)
            except Exception as exc:
                # Best-effort cleanup: ignore failure but log for diagnostics.
                logging.getLogger(__name__).warning(
                    "Failed to remove CatPred work directory %s: %s",
                    catpred_work_dir,
                    exc,
                )

        source = "catpred"
        errors = [results.get(f"{p}_error") for p in ["km", "kcat"] if f"{p}_error" in results]
        if errors:
            source = f"catpred(errors: {'; '.join(e for e in errors if e)})"[:255]

        return EstimatedKinetics(
            km=results.get("km"),
            km_per_substrate=results.get("km_per_substrate"),
            km_sd=results.get("km_sd"),
            km_sd_per_substrate=results.get("km_sd_per_substrate"),
            kcat=results.get("kcat"),
            ki=results.get("ki"),
            kcat_sd=results.get("kcat_sd"),
            ki_sd=results.get("ki_sd"),
            source=source
        )


def build_estimator(
    name: Optional[str],
    include_sd: bool = False,
) -> Optional[BaseKineticsEstimator]:
    """
    Return the appropriate kinetics estimator based on the name.

    Args:
        name: Estimator backend name (e.g. 'catpred')
        include_sd: If True, estimator will include SD (standard deviation) from predictions when available.
    """
    if not name:
        return None
    if name == "catpred":
        return CatPredEstimator(include_sd=include_sd)
    raise ValueError(f"Unknown kinetics_estimator: {name}")

