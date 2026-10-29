#!/usr/bin/env python3

"""Same total dose, different delivery: how much of a pathway's output does each arm cost?

Builds the model of a project through the normal BEES pipeline (without its dosing), then
simulates four arms of EQUAL total dose of one inhibitor species:

    control      no drug
    free bolus   all of the dose at t = 0
    liposome     biexponential release: 10 % with a 2 h half-life, 90 % with 45 h
    slow depot   Weibull release, tau = 24 h, beta = 0.6

with the same first-order clearance of the free drug in every arm. Per arm it reports the
output lost against the control at the end, the time the output flux spends below half the
control flux, the inhibition AUC and the peak free drug; it writes a CSV and a 3-panel PNG.

    python devtools/release_sweep.py -i projects/methotrexate_dhfr/input.yml \\
        --drug methotrexate --output-species dTMP

The kinetics of the shipped projects are placeholders, see their READMEs.
"""

import argparse
import csv
import math
import os
import sys
from typing import Dict, List, Optional, Sequence

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep_common as sc  # noqa: E402

from bees.dosing import DoseSchedule  # noqa: E402
from bees.simulator import ODESimulator  # noqa: E402

H = 3600.0
ARMS = ("control", "free bolus", "liposome", "slow depot")
# Release shapes of the arms (before --time-scale).
LIPOSOME = {"fraction_fast": 0.1, "half_life_fast": 2 * H, "half_life_slow": 45 * H}
DEPOT = {"tau": 24 * H, "beta": 0.6}


def make_arms(drug: str, dose: float, time_scale: float, elimination_half_life: Optional[float]) -> Dict[str, List[DoseSchedule]]:
    """Dose schedules of every arm; release times are multiplied by ``time_scale``."""
    k_el = math.log(2.0) / elimination_half_life if elimination_half_life else 0.0
    liposome = dict(LIPOSOME, half_life_fast=LIPOSOME["half_life_fast"] * time_scale,
                    half_life_slow=LIPOSOME["half_life_slow"] * time_scale)
    depot = dict(DEPOT, tau=DEPOT["tau"] * time_scale)
    return {
        "control": [],
        "free bolus": [DoseSchedule(drug, "bolus", dose, elimination_k=k_el)],
        "liposome": [DoseSchedule(drug, "biexponential", dose, liposome, elimination_k=k_el)],
        "slow depot": [DoseSchedule(drug, "weibull", dose, depot, elimination_k=k_el)],
    }


def inhibition_auc_hours(model, bees_object, t: np.ndarray, traj_matrix: np.ndarray, labels: Sequence[str]) -> float:
    """Largest, over the inhibited reactions, of the integral of the fractional rate loss (h).

    The loss at each time is 1 - v/v0 for the reaction evaluated at the arm's own state with
    and without its inhibitors, so it measures the drug's action on the enzyme itself, before
    the pathway buffers it. Zero when no reaction is inhibited.
    """
    sim = ODESimulator(model)
    with_inh = sim._build_rhs()
    inhibited = sorted(with_inh._inhibited_rx)
    if not inhibited:
        return 0.0
    saved = {id(r): r.inhibitors for r in model.core_reactions + model.edge_reactions}
    try:
        for r in model.core_reactions + model.edge_reactions:
            r.inhibitors = None
        without = sim._build_rhs()
    finally:
        for r in model.core_reactions + model.edge_reactions:
            r.inhibitors = saved[id(r)]
    v = with_inh.compute_v(traj_matrix)
    v0 = without.compute_v(traj_matrix)
    best = 0.0
    for j in inhibited:
        with np.errstate(divide="ignore", invalid="ignore"):
            loss = np.where(v0[j] > 0, 1.0 - v[j] / v0[j], 0.0)
        area = float(np.sum(0.5 * (loss[1:] + loss[:-1]) * np.diff(t)))
        best = max(best, area)
    return sc.hours(best)


def time_below_half(t: np.ndarray, flux: np.ndarray, flux_control: np.ndarray) -> float:
    """Hours during which the output flux is below half the control flux (where the control has flux)."""
    live = flux_control > 1e-6 * max(float(np.max(flux_control)), 1e-300)
    below = live & (flux < 0.5 * flux_control)
    dt = np.diff(t)
    inside = below[1:] & below[:-1]
    return sc.hours(float(np.sum(dt[inside] + 0.0)))


def run_sweep(model, bees_object, drug, output_labels, dose, time_scale, elimination_half_life, end_time):
    arms = make_arms(drug, dose, time_scale, elimination_half_life)
    runs = {}
    for name, schedules in arms.items():
        t, traj = sc.simulate(model, bees_object, schedules, end_time)
        out, used, skipped = sc.summed_output(traj, output_labels)
        matrix = np.stack([traj[lab.lower().strip()] for lab in _model_labels(model)], axis=0)
        runs[name] = {"t": t, "out": out, "flux": np.gradient(out, t), "drug": traj.get(drug.lower().strip(), np.zeros_like(t)),
                      "matrix": matrix, "used": used, "skipped": skipped}
    ctrl = runs["control"]
    rows = []
    for name in ARMS:
        r = runs[name]
        rows.append({
            "arm": name,
            "total_dose_uM": 1000.0 * dose if name != "control" else 0.0,
            "output_lost_pct": sc.percent_lost(ctrl["out"], r["out"], -1) if name != "control" else 0.0,
            "hours_below_half_control_flux": time_below_half(ctrl["t"], r["flux"], ctrl["flux"]) if name != "control" else 0.0,
            "inhibition_auc_h": inhibition_auc_hours(model, bees_object, r["t"], r["matrix"], _model_labels(model)) if name != "control" else 0.0,
            "peak_free_drug_uM": 1000.0 * float(np.max(r["drug"])) if name != "control" else 0.0,
        })
    return runs, rows


def _model_labels(model) -> List[str]:
    return [lab.lower().strip() for lab in model.get_all_species_labels()]


def write_outputs(runs, rows, out_dir, title, drug, output_note):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, "release_sweep.csv")
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in row.items()})

    t_h = runs["control"]["t"] / H
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
    ctrl = runs["control"]
    for name in ARMS:
        r, color = runs[name], sc.PALETTE[name]
        if name != "control":
            axes[0].plot(t_h, 1000.0 * r["drug"], color=color, lw=2, label=name)
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.where(ctrl["flux"] > 0, r["flux"] / ctrl["flux"], np.nan)
        axes[1].plot(t_h, rel, color=color, lw=2, label=name, ls="--" if name == "control" else "-")
        axes[2].plot(t_h, r["out"], color=color, lw=2, label=name, ls="--" if name == "control" else "-")
    axes[0].set(xlabel="time (h)", ylabel=f"free {drug} (µM)", title="Free drug")
    axes[1].set(xlabel="time (h)", ylabel="output flux / control flux", title="Output flux", ylim=(-0.05, 1.15))
    axes[2].set(xlabel="time (h)", ylabel=f"cumulative {output_note} (mM)", title="Cumulative output")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(alpha=0.25)
    axes[1].legend(frameon=False)
    fig.suptitle(f"{title}: equal total dose, different delivery", fontsize=12)
    png_path = os.path.join(out_dir, "release_sweep.png")
    fig.savefig(png_path, dpi=150)
    plt.close(fig)
    return csv_path, png_path


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, help="BEES project input.yml")
    p.add_argument("--drug", required=True, help="label of the dosed inhibitor species")
    p.add_argument("--output-species", required=True,
                   help="species whose accumulation is the output; several labels separated by ';' are summed "
                        "(labels not in the model are skipped)")
    p.add_argument("--time-scale", type=float, default=1.0, help="multiply the release times of the arms (default 1)")
    p.add_argument("--dose", type=float, help="total dose in mM (default: the input's dosing entry for the drug)")
    p.add_argument("--elimination-half-life", type=float, help="free-drug half-life in s (default: the input's)")
    p.add_argument("--ki", type=float, help="replace the Ki (mM) of every inhibition entry of the drug")
    p.add_argument("--end-time", type=float, help="simulated time in s (default: settings.end_time)")
    p.add_argument("--out-dir", help="output directory (default: <project>/output/release_sweep)")
    return p.parse_args(argv)


def main(argv=None) -> List[dict]:
    args = parse_args(argv)
    raw = sc.common.read_yaml_file(args.input)
    out_dir = args.out_dir or os.path.join(sc.REPO_ROOT, raw.get("project_directory", "."), "output", "release_sweep")
    data = sc.load_input(args.input, out_dir, end_time=args.end_time, drug=args.drug, ki=args.ki)
    entry = sc.dosing_entry(raw, args.drug)

    dose = args.dose if args.dose is not None else (entry or {}).get("dose")
    if dose is None:
        sys.exit(f"No dose: give --dose (mM) or a dosing entry for '{args.drug}' in {args.input}")
    t_half = args.elimination_half_life
    if t_half is None and entry:
        t_half = entry.get("elimination_half_life")
        if t_half is None and entry.get("elimination_k"):
            t_half = math.log(2.0) / entry["elimination_k"]
    end_time = data["settings"].get("end_time")
    if end_time is None:
        sys.exit("The input has no settings.end_time; give --end-time")

    sc.log_step(f"Building the model of {data['project']} ...")
    bees_object, _, model = sc.build_model(data)
    if args.drug.lower().strip() not in _model_labels(model):
        sys.exit(f"'{args.drug}' is not a species of the model")
    labels = sc.parse_labels(args.output_species)
    sc.log_step(f"Simulating {len(ARMS)} arms of {1000 * dose:g} uM {args.drug} over {end_time / H:g} h "
                f"(time scale {args.time_scale:g}, free-drug half-life "
                f"{'none' if not t_half else f'{t_half:g} s'}) ...")
    runs, rows = run_sweep(model, bees_object, args.drug, labels, dose, args.time_scale, t_half, end_time)
    used, skipped = runs["control"]["used"], runs["control"]["skipped"]
    if skipped:
        sc.log_step(f"Output species not in the model, skipped: {', '.join(skipped)}")
    csv_path, png_path = write_outputs(runs, rows, out_dir, data["project"], args.drug,
                                       "+".join(used) if len(used) <= 2 else f"{len(used)} species")

    header = f"{'arm':<12}{'lost %':>9}{'h < 1/2 flux':>14}{'inhib AUC h':>13}{'peak uM':>11}"
    print("\n" + header + "\n" + "-" * len(header))
    for r in rows:
        print(f"{r['arm']:<12}{r['output_lost_pct']:>9.1f}{r['hours_below_half_control_flux']:>14.2f}"
              f"{r['inhibition_auc_h']:>13.2f}{r['peak_free_drug_uM']:>11.4g}")
    print(f"\nwrote {csv_path}\nwrote {png_path}")
    return rows


if __name__ == "__main__":
    main()
