#!/usr/bin/env python3

"""Which enzyme is worth hitting? Knock each enzyme down in turn and see what the output loses.

Builds the model of a project through the normal BEES pipeline, then for every enzyme scales
kcat of all its reactions by --residual (default 0.1, i.e. 90 % inhibition), re-simulates and
reports the percentage of the output lost against the control at 1, 3 and 6 h. Writes a CSV
and a bar chart.

    python devtools/target_scan.py -i projects/fattyAcidSynthesis/fas_inhibitor_delivery/input.yml \\
        --output-species "hexanoate;octanoate;decanoate"

This ranks targets by control of the output flux, not by how druggable they are.
"""

import argparse
import csv
import os
import sys
from typing import Dict, List, Sequence

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep_common as sc  # noqa: E402

H = 3600.0


def scale_enzyme(model, enzyme_label: str, residual: float):
    """Scale kcat (or Vmax) of every reaction of ``enzyme_label``; returns an undo function."""
    saved = []
    for rxn in model.core_reactions + model.edge_reactions:
        if str(rxn.enzyme_label).lower().strip() != enzyme_label.lower().strip():
            continue
        kin = rxn.kinetics
        saved.append((kin, kin.kcat, kin.vmax))
        if kin.kcat is not None:
            kin.kcat = kin.kcat * residual
        if kin.vmax is not None:
            kin.vmax = kin.vmax * residual

    def undo():
        for kin, kcat, vmax in saved:
            kin.kcat, kin.vmax = kcat, vmax

    undo.n_reactions = len(saved)
    return undo


def scan(model, bees_object, enzymes: Sequence[str], output_labels, residual: float, times_h: Sequence[float], end_time: float):
    """Rows {enzyme, n_reactions, lost_pct_<t>h...} sorted by the loss at the last time."""
    t, traj = sc.simulate(model, bees_object, [], end_time)
    control, used, skipped = sc.summed_output(traj, output_labels)
    idx = [int(np.argmin(np.abs(t - h * H))) for h in times_h]
    rows = []
    for name in enzymes:
        undo = scale_enzyme(model, name, residual)
        try:
            if undo.n_reactions == 0:
                row = {"enzyme": name, "n_reactions": 0}
                row.update({f"lost_pct_{h:g}h": float("nan") for h in times_h})
            else:
                _, traj_k = sc.simulate(model, bees_object, [], end_time)
                out, _, _ = sc.summed_output(traj_k, output_labels)
                row = {"enzyme": name, "n_reactions": undo.n_reactions}
                row.update({f"lost_pct_{h:g}h": sc.percent_lost(control, out, i) for h, i in zip(times_h, idx)})
        finally:
            undo()
        rows.append(row)
    key = f"lost_pct_{times_h[-1]:g}h"
    rows.sort(key=lambda r: (np.isnan(r[key]), -(r[key] if not np.isnan(r[key]) else 0.0)))
    return rows, used, skipped


def write_outputs(rows, times_h, out_dir, title, residual):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, "target_scan.csv")
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in row.items()})

    shades = ["#9ecae1", "#4292c6", "#08519c"][-len(times_h):] if len(times_h) <= 3 else None
    fig, ax = plt.subplots(figsize=(max(6, 0.7 * len(rows) + 2), 4.2), constrained_layout=True)
    width = 0.8 / len(times_h)
    x = np.arange(len(rows))
    for k, h in enumerate(times_h):
        vals = [np.nan_to_num(r[f"lost_pct_{h:g}h"]) for r in rows]
        ax.bar(x + (k - (len(times_h) - 1) / 2) * width, vals, width, label=f"{h:g} h",
               color=shades[k] if shades else None)
    ax.set_xticks(x, [r["enzyme"] for r in rows], rotation=45, ha="right")
    ax.set(ylabel="output lost vs control (%)", ylim=(0, 100),
           title=f"{title}: kcat x {residual:g} for one enzyme at a time")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, title="after")
    png_path = os.path.join(out_dir, "target_scan.png")
    fig.savefig(png_path, dpi=150)
    plt.close(fig)
    return csv_path, png_path


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, help="BEES project input.yml")
    p.add_argument("--output-species", required=True,
                   help="species whose accumulation is the output; several labels separated by ';' are summed "
                        "(labels not in the model are skipped)")
    p.add_argument("--residual", type=float, default=0.1, help="fraction of kcat left (default 0.1)")
    p.add_argument("--times-h", default="1,3,6", help="report the loss at these hours (default 1,3,6)")
    p.add_argument("--end-time", type=float, help="simulated time in s (default: settings.end_time)")
    p.add_argument("--out-dir", help="output directory (default: <project>/output/target_scan)")
    return p.parse_args(argv)


def main(argv=None) -> List[dict]:
    args = parse_args(argv)
    if not 0.0 <= args.residual < 1.0:
        sys.exit("--residual must be in [0, 1)")
    raw = sc.common.read_yaml_file(args.input)
    out_dir = args.out_dir or os.path.join(sc.REPO_ROOT, raw.get("project_directory", "."), "output", "target_scan")
    data = sc.load_input(args.input, out_dir, end_time=args.end_time)
    end_time = data["settings"].get("end_time")
    if end_time is None:
        sys.exit("The input has no settings.end_time; give --end-time")
    times_h = [float(x) for x in args.times_h.split(",")]
    late = [h for h in times_h if h * H > end_time * (1 + 1e-9)]
    if late:
        print(f"Skipping times beyond the end of the run ({end_time / H:g} h): {late}")
        times_h = [h for h in times_h if h not in late]
    if not times_h:
        sys.exit("No reporting time is within the simulated span")

    sc.log_step(f"Building the model of {data['project']} ...")
    bees_object, _, model = sc.build_model(data)
    enzymes = [e.label for e in bees_object.enzymes if e.reactive]
    labels = sc.parse_labels(args.output_species)
    sc.log_step(f"Scanning {len(enzymes)} enzymes (kcat x {args.residual:g}) ...")
    rows, used, skipped = scan(model, bees_object, enzymes, labels, args.residual, times_h, end_time)
    if skipped:
        sc.log_step(f"Output species not in the model, skipped: {', '.join(skipped)}")
    csv_path, png_path = write_outputs(rows, times_h, out_dir, data["project"], args.residual)

    cols = [f"lost_pct_{h:g}h" for h in times_h]
    header = f"{'enzyme':<22}{'rxns':>5}" + "".join(f"{c.replace('lost_pct_', '% @ '):>10}" for c in cols)
    print("\n" + header + "\n" + "-" * len(header))
    for r in rows:
        print(f"{r['enzyme']:<22}{r['n_reactions']:>5}" + "".join(f"{r[c]:>10.1f}" for c in cols))
    print(f"\nwrote {csv_path}\nwrote {png_path}")
    return rows


if __name__ == "__main__":
    main()
