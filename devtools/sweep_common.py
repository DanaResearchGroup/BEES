#!/usr/bin/env python3

"""Shared helpers of devtools/release_sweep.py and devtools/target_scan.py.

Both build the model through the normal BEES pipeline (iterative enlargement of the project
input, without its dosing) and then re-simulate that fixed model with something changed.
"""

import copy
import math
import os
import sys
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import bees.common as common  # noqa: E402
from bees.dosing import DoseSchedule  # noqa: E402
from bees.simulator import ODESimulator  # noqa: E402

# Okabe-Ito: distinguishable for colour-blind readers and in greyscale.
PALETTE = {
    "control": "#000000",
    "free bolus": "#D55E00",
    "liposome": "#0072B2",
    "slow depot": "#009E73",
}


def load_input(
    path: str,
    out_dir: str,
    *,
    end_time: Optional[float] = None,
    drug: Optional[str] = None,
    ki: Optional[float] = None,
) -> dict:
    """The project's input as a dict, redirected to ``out_dir`` and with its dosing taken out.

    ``ki`` (mM) replaces the Ki of every inhibition entry that names ``drug``.
    """
    data = common.read_yaml_file(path)
    data = copy.deepcopy(data)
    data["project_directory"] = os.path.join(os.path.abspath(out_dir), "model")
    settings = data.setdefault("settings", {})
    settings["output_directory"] = None
    settings.update(save_simulation_plots=False, save_simulation_profiles=False, save_ode_equations=False)
    settings["verbose"] = 30
    if end_time is not None:
        settings["end_time"] = float(end_time)
    if ki is not None:
        if drug is None:
            raise ValueError("--ki needs --drug (the inhibitor whose Ki it replaces)")
        hit = False
        for enz in data.get("enzymes", []):
            for inh in enz.get("inhibition") or []:
                if str(inh.get("inhibitor", "")).lower().strip() == drug.lower().strip():
                    inh["ki"] = float(ki)
                    hit = True
        if not hit:
            raise ValueError(f"--ki: no enzyme lists '{drug}' under inhibition")
    return data


def dosing_entry(data: dict, drug: str) -> Optional[dict]:
    """The input's dosing entry for ``drug`` (first one), if any."""
    for d in data.get("dosing") or []:
        if str(d.get("species", "")).lower().strip() == drug.lower().strip():
            return d
    return None


def build_model(data: dict):
    """Run the enlargement of ``data`` (dosing removed): (bees_object, enlarger, model).

    The inhibitors stay attached to the reactions, so the model can be simulated with any
    dosing afterwards.
    """
    from bees.enlarger import IterativeEnlarger
    from bees.main import BEES, resolve_db_path
    from bees.reaction_generator import ReactionGenerator

    data = copy.deepcopy(data)
    data.pop("dosing", None)
    app = BEES(input_data=data)
    generator = ReactionGenerator(app.bees_object, app.logger, app.output_directory)
    generator.load_kinetic_database(resolve_db_path(app.bees_object.database.name, app.logger), ontology=app.ontology)
    enlarger = IterativeEnlarger(app.bees_object, generator, app.logger, app.output_directory)
    result = enlarger.run()
    enlarger._attach_thermo_to_reactions()
    if not result.model.core_reactions:
        raise RuntimeError("The enlargement produced no reactions; check the input and the log in " + app.output_directory)
    return app.bees_object, enlarger, result.model


def solver_settings(bees_object) -> dict:
    s = bees_object.settings
    return {
        "method": s.ode_method or "BDF",
        "rtol": float(s.ode_rtol) if s.ode_rtol is not None else 1e-8,
        "atol": float(s.ode_atol) if s.ode_atol is not None else 1e-10,
    }


def simulate(model, bees_object, schedules: Sequence[DoseSchedule], end_time: float, n_points: int = 721):
    """Simulate the whole model from its initial concentrations; (t, {label: trajectory})."""
    model.reset_concentrations_to_initial()
    sim = ODESimulator(model, dosing=list(schedules))
    res = sim.simulate(end_time=end_time, time_step=end_time / (n_points - 1), **solver_settings(bees_object))
    if not res.success:
        raise RuntimeError(f"ODE solve failed: {res.message}")
    return res.t, {lab.lower().strip(): res.y[i] for i, lab in enumerate(res.species_labels)}


def parse_labels(spec: str) -> List[str]:
    return [s.strip() for s in spec.split(";") if s.strip()]


def summed_output(traj: Dict[str, np.ndarray], labels: Sequence[str]) -> Tuple[np.ndarray, List[str], List[str]]:
    """Sum of the trajectories of ``labels`` (those in the model); also (used, skipped) labels."""
    used = [lab for lab in labels if lab.lower().strip() in traj]
    skipped = [lab for lab in labels if lab not in used]
    if not used:
        raise ValueError(f"None of the output species {list(labels)} is in the model; it has: {sorted(traj)}")
    total = sum(traj[lab.lower().strip()] for lab in used)
    return np.asarray(total, dtype=float), used, skipped


def percent_lost(control: np.ndarray, other: np.ndarray, index: int) -> float:
    base = float(control[index])
    return 100.0 * (base - float(other[index])) / base if base > 0 else float("nan")


def hours(seconds: float) -> float:
    return seconds / 3600.0


def log_step(message: str) -> None:
    print(message, flush=True)
