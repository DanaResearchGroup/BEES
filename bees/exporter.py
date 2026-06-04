#!/usr/bin/env python3

"""
Enlarger Exporter Module
-----------------------
Holds all export/plot functionality for the iterative enlarger.


"""

from __future__ import annotations

import csv
import hashlib
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple, Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from bees.common import is_general_cofactor_label
from bees.core_edge_model import CoreEdgeModel, SpeciesData
from bees.reaction_generator import GeneratedReaction
from bees.simulator import SimulationResult

try:
    import libsbml as _libsbml  # python-libsbml
    _LIBSBML_AVAILABLE = True
except ImportError:  
    _libsbml = None
    _LIBSBML_AVAILABLE = False


def reaction_signature(
    reaction: GeneratedReaction,
) -> Tuple[str, Tuple[str, ...], Tuple[str, ...]]:
    """
    Canonical reaction signature for stable reaction ID tracking.
    """
    enzyme = str(reaction.enzyme_label).lower().strip()
    reactants = tuple(sorted(str(r).lower().strip() for r in reaction.reactant_labels))
    products = tuple(sorted(str(p).lower().strip() for p in reaction.product_labels))
    return (enzyme, reactants, products)


@dataclass
class EnlargerExporter:
    """
    Exporter for IterativeEnlarger outputs (CSVs, plots, reaction tree).

    Reads state that was produced by the enlarger module and writes artifacts to disk.
    """

    model: CoreEdgeModel
    profiles: List[SimulationResult]
    output_directory: str
    logger: Any

    reaction_id_by_sig: Dict[Tuple[str, Tuple[str, ...], Tuple[str, ...]], int]
    reaction_first_seen_iter: Dict[Tuple[str, Tuple[str, ...], Tuple[str, ...]], int]
    reaction_core_enter_iter: Dict[Tuple[str, Tuple[str, ...], Tuple[str, ...]], Optional[int]]
    reaction_obj_by_sig: Dict[Tuple[str, Tuple[str, ...], Tuple[str, ...]], GeneratedReaction]
    iteration_summaries: List[Dict[str, int]]

    # Plot/export settings
    save_reaction_tree_plots: bool = False
    save_simulation_plots: bool = True
    plot_max_species: Optional[int] = None
    plot_exclude_enzymes: bool = True
    plot_exclude_cofactors: bool = True
    reaction_tree_layout: str = "graphviz"
    reaction_tree_rankdir: str = "TB"
    reaction_tree_fontsize: int = 8

    # Stateful across iterations (used by reaction-tree coloring)
    core_seen_labels: Optional[Set[str]] = None
    # Optional access to original input (used for cofactor exclusion in plots)
    bees_object: Optional[Any] = None

    # ------------------------------------------------------------------
    # Export helpers - flux analysis
    # ------------------------------------------------------------------

    def export_flux_analysis(self, filename: str = "flux_analysis.csv") -> Optional[str]:
        """
        Write per-iteration flux data (core/edge growth summary + reaction history).
        """
        output_path = os.path.join(self.output_directory, filename)
        with open(output_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    "iteration",
                    "core_species",
                    "edge_species",
                    "core_reactions",
                    "edge_reactions",
                ]
            )

            if self.iteration_summaries:
                for row in self.iteration_summaries:
                    writer.writerow(
                        [
                            row["iteration"],
                            row["core_species"],
                            row["edge_species"],
                            row["core_reactions"],
                            row["edge_reactions"],
                        ]
                    )
            else:
                s = self.model.summary()
                writer.writerow(
                    [
                        len(self.profiles),
                        s["core_species"],
                        s["edge_species"],
                        s["core_reactions"],
                        s["edge_reactions"],
                    ]
                )

        self.logger.info(f"Exported flux analysis to {output_path}")

        details_path = os.path.join(self.output_directory, "flux_analysis_reactions.csv")
        core_sigs = {reaction_signature(rxn) for rxn in self.model.core_reactions}
        edge_sigs = {reaction_signature(rxn) for rxn in self.model.edge_reactions}
        all_sigs = sorted(
            set(core_sigs) | set(edge_sigs),
            key=lambda sig: self.reaction_id_by_sig.get(sig, 10**9),
        )
        with open(details_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    "reaction_id",
                    "classification",
                    "first_seen_iteration",
                    "core_enter_iteration",
                ]
            )
            for sig in all_sigs:
                rxn = self.reaction_obj_by_sig.get(sig)
                if rxn is None:
                    continue
                reaction_id = f"R{self.reaction_id_by_sig.get(sig, 0)}"
                classification = "core" if sig in core_sigs else "edge"
                first_seen = self.reaction_first_seen_iter.get(sig, "")
                core_enter = self.reaction_core_enter_iter.get(sig, "")
                writer.writerow(
                    [
                        reaction_id,
                        classification,
                        first_seen,
                        core_enter if core_enter is not None else "",
                    ]
                )
        self.logger.info(f"Exported reaction flux summary to {details_path}")
        return output_path

    # ------------------------------------------------------------------
    # Export helpers - core/edge CSVs
    # ------------------------------------------------------------------

    def export_core_edge_reaction_species_csvs(self) -> Dict[str, str]:
        """
        Export core/edge reaction tables, each followed by a species section.

        Returns:
            Dict with keys "core" and "edge" and absolute output file paths.
        """
        outputs: Dict[str, str] = {}

        def _write_one(
            path: str,
            reactions: List[GeneratedReaction],
            species: List[SpeciesData],
            section_name: str,
        ) -> None:
            with open(path, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(
                    [
                        "index",
                        "reaction_id",
                        "template",
                        "ec_number",
                        "family",
                        "enzyme",
                        "substrate",
                        "reactants",
                        "products",
                    ]
                )
                for idx, rxn in enumerate(reactions, start=1):
                    sig = reaction_signature(rxn)
                    rid_num = self.reaction_id_by_sig.get(sig)
                    rid = f"R{rid_num}" if rid_num is not None else ""
                    template_type = ""
                    family = ""
                    if getattr(rxn, "template", None) is not None:
                        template_type = str(getattr(rxn.template, "template_type", "") or "")
                        ec_class_obj = getattr(rxn.template, "ec_class", None)
                        family = str(getattr(ec_class_obj, "name", "") or "")
                    writer.writerow(
                        [
                            idx,
                            rid,
                            template_type,
                            str(getattr(rxn, "ec_number", "") or ""),
                            family,
                            str(getattr(rxn, "enzyme_label", "") or ""),
                            str(getattr(rxn, "substrate_label", "") or ""),
                            " + ".join(str(r) for r in getattr(rxn, "reactant_labels", []) or []),
                            " + ".join(str(p) for p in getattr(rxn, "product_labels", []) or []),
                        ]
                    )

                writer.writerow([])
                writer.writerow([f"{section_name} species"])
                writer.writerow(["index", "label", "is_enzyme", "constant"])
                for idx, sp in enumerate(species, start=1):
                    writer.writerow(
                        [
                            idx,
                            str(getattr(sp, "label", "") or ""),
                            bool(getattr(sp, "is_enzyme", False)),
                            bool(getattr(sp, "constant", False)),
                        ]
                    )

        core_path = os.path.join(self.output_directory, "core_reactions_species.csv")
        edge_path = os.path.join(self.output_directory, "edge_reactions_species.csv")

        _write_one(
            core_path,
            list(self.model.core_reactions),
            list(self.model.core_species),
            section_name="core",
        )
        _write_one(
            edge_path,
            list(self.model.edge_reactions),
            list(self.model.edge_species),
            section_name="edge",
        )

        self.logger.info(f"Exported core reactions/species CSV to {core_path}")
        self.logger.info(f"Exported edge reactions/species CSV to {edge_path}")
        outputs["core"] = core_path
        outputs["edge"] = edge_path
        return outputs

    # ------------------------------------------------------------------
    # Export helpers - simulation profiles
    # ------------------------------------------------------------------

    def export_simulation_profiles(self, filename: str = "simulation_profiles.csv") -> Optional[str]:
        """
        Write concentration time-series to CSV.

        Concatenates profiles from every iteration into a single file.
        Returns the output path, or None if there are no profiles.
        """
        if not self.profiles:
            return None

        output_path = os.path.join(self.output_directory, filename)
        all_labels: List[str] = []
        seen: Set[str] = set()
        for prof in self.profiles:
            for lab in prof.species_labels:
                if lab not in seen:
                    all_labels.append(lab)
                    seen.add(lab)

        with open(output_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["iteration", "time"] + all_labels)
            for it_idx, prof in enumerate(self.profiles, 1):
                label_to_row = {lab: i for i, lab in enumerate(prof.species_labels)}
                for t_idx in range(prof.t.shape[0]):
                    row = [it_idx, prof.t[t_idx]]
                    for lab in all_labels:
                        ridx = label_to_row.get(lab)
                        if ridx is not None:
                            row.append(prof.y[ridx, t_idx])
                        else:
                            row.append("")
                    writer.writerow(row)

        self.logger.info(f"Exported simulation profiles to {output_path}")
        return output_path

    # ------------------------------------------------------------------
    # Export helpers - reaction tree (visualisation)
    # ------------------------------------------------------------------

    def export_reaction_tree(
        self,
        *,
        iteration: int,
        promoted_labels: Optional[List[str]] = None,
    ) -> None:
        """
        Export a reaction tree plot for this iteration.

        The tree:
        - Includes only core species.
        - Excludes enzymes and general cofactors.
        - Colours newly promoted core species in this iteration differently
          from species that were already in the core.
        """
        if not self.save_reaction_tree_plots or not self.model:
            return

        if self.core_seen_labels is None:
            self.core_seen_labels = set()

        core_nodes: List[SpeciesData] = []
        for sd in self.model.core_species:
            if sd.is_enzyme:
                continue
            if is_general_cofactor_label(sd.label):
                continue
            core_nodes.append(sd)

        if not core_nodes:
            return

        labels = [sd.label for sd in core_nodes]
        label_set = set(labels)

        promoted_set = set(promoted_labels or [])
        new_nodes: Set[str] = set()
        for lab in labels:
            if lab in promoted_set and lab not in self.core_seen_labels:
                new_nodes.add(lab)

        self.core_seen_labels.update(labels)

        edges: List[tuple] = []
        edge_reaction_labels: Dict[tuple, Set[str]] = {}
        for rxn in self.model.core_reactions:
            sig = reaction_signature(rxn)
            reaction_id = self.reaction_id_by_sig.get(sig)
            reaction_tag = f"R{reaction_id}" if reaction_id is not None else ""
            for reactant in rxn.reactant_labels:
                if reactant not in label_set:
                    continue
                for product in rxn.product_labels:
                    if product not in label_set:
                        continue
                    if reactant == product:
                        continue
                    edge = (reactant, product)
                    edges.append(edge)
                    if reaction_tag:
                        edge_reaction_labels.setdefault(edge, set()).add(reaction_tag)

        edges = sorted(set(edges))
        if not edges:
            return

        xs: Dict[str, float] = {}
        ys: Dict[str, float] = {}

        def _simple_layout() -> None:
            old_labels = [lab for lab in labels if lab not in new_nodes]
            new_labels_ordered = [lab for lab in labels if lab in new_nodes]

            def _assign_row(row_labels: List[str], y_val: float) -> None:
                n = len(row_labels)
                if n == 0:
                    return
                if n == 1:
                    xs[row_labels[0]] = 0.5
                    ys[row_labels[0]] = y_val
                    return
                for i, lab in enumerate(row_labels):
                    xs[lab] = i / (n - 1)
                    ys[lab] = y_val

            _assign_row(old_labels, y_val=0.0)
            _assign_row(new_labels_ordered, y_val=-1.0)

        def _graphviz_layout() -> bool:
            dot_exe = shutil.which("dot")
            if not dot_exe:
                return False

            rankdir = str(self.reaction_tree_rankdir or "TB").upper()
            if rankdir not in {"TB", "BT", "LR", "RL"}:
                rankdir = "TB"

            node_id: Dict[str, str] = {}
            for i, lab in enumerate(labels):
                node_id[lab] = f"n{i}"

            dot_lines: List[str] = [
                "digraph ReactionTree {",
                f'  rankdir="{rankdir}";',
                "  splines=true;",
                "  overlap=false;",
                "  nodesep=0.35;",
                "  ranksep=0.6;",
                "  node [shape=circle];",
            ]
            for lab in labels:
                dot_lines.append(f'  {node_id[lab]} [label="{node_id[lab]}"];')
            for src, dst in edges:
                dot_lines.append(f"  {node_id[src]} -> {node_id[dst]};")
            dot_lines.append("}")
            dot = "\n".join(dot_lines)

            try:
                proc = subprocess.run(
                    [dot_exe, "-Tplain"],
                    input=dot.encode("utf-8"),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    check=True,
                )
            except Exception:
                return False

            positions_raw: Dict[str, tuple[float, float]] = {}
            for line in proc.stdout.decode("utf-8", errors="replace").splitlines():
                if not line.startswith("node "):
                    continue
                parts = line.split()
                if len(parts) < 4:
                    continue
                name = parts[1]
                try:
                    x = float(parts[2])
                    y = float(parts[3])
                except ValueError:
                    continue
                positions_raw[name] = (x, y)

            if not positions_raw:
                return False

            xs_vals = [p[0] for p in positions_raw.values()]
            ys_vals = [p[1] for p in positions_raw.values()]
            min_x, max_x = min(xs_vals), max(xs_vals)
            min_y, max_y = min(ys_vals), max(ys_vals)
            span_x = max(1e-9, max_x - min_x)
            span_y = max(1e-9, max_y - min_y)

            inv_node_id = {v: k for k, v in node_id.items()}
            for nid, (x, y) in positions_raw.items():
                lab = inv_node_id.get(nid)
                if not lab:
                    continue
                xs[lab] = (x - min_x) / span_x
                ys[lab] = (y - min_y) / span_y
            return len(xs) > 0 and len(ys) > 0

        used_graphviz = (
            str(self.reaction_tree_layout or "graphviz").lower() == "graphviz"
            and _graphviz_layout()
        )
        if not used_graphviz:
            _simple_layout()

        fig, ax = plt.subplots(figsize=(16, 12), facecolor="white")
        ax.set_facecolor("white")

        for src, dst in edges:
            ax.annotate(
                "",
                xy=(xs[dst], ys[dst]),
                xytext=(xs[src], ys[src]),
                arrowprops=dict(
                    arrowstyle="->",
                    color="#555555",
                    linewidth=1.0,
                    alpha=0.8,
                ),
            )
            rid_set = edge_reaction_labels.get((src, dst), set())
            if rid_set:
                rid_text = ",".join(
                    sorted(
                        rid_set,
                        key=lambda s: int(s[1:])
                        if s.startswith("R") and s[1:].isdigit()
                        else 10**9,
                    )
                )
                mid_x = (xs[src] + xs[dst]) / 2.0
                mid_y = (ys[src] + ys[dst]) / 2.0
                ax.annotate(
                    rid_text,
                    xy=(mid_x, mid_y),
                    xycoords="data",
                    xytext=(0, 7),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=10,
                    color="#222222",
                    zorder=5,
                    bbox=dict(
                        boxstyle="round,pad=0.18",
                        facecolor="white",
                        edgecolor="none",
                        alpha=0.85,
                    ),
                )

        old_color = "#A6CEE3"
        new_color = "#FB9A99"

        def _shorten_label(label: str) -> str:
            s = str(label).strip()
            s = re.sub(r"\s+", " ", s)
            if len(s) <= 10 and " " not in s:
                return s
            tokens = re.split(r"[\s\-_]+", s)
            keep = []
            for t in tokens:
                if not t:
                    continue
                if t.isdigit() or re.fullmatch(r"\d+[A-Za-z]*", t or ""):
                    keep.append(t)
                else:
                    keep.append(t[0].upper())
            base = "".join(keep) or s[:6].upper()
            if len(base) > 12:
                base = base[:12]
            return base

        cumulative_mapping_path = os.path.join(
            self.output_directory,
            "reaction_tree_labels.csv",
        )
        full_to_short: Dict[str, str] = {}
        first_seen_by_full: Dict[str, int] = {}
        used_short_labels: Set[str] = set()

        if os.path.exists(cumulative_mapping_path):
            try:
                with open(cumulative_mapping_path, "r", newline="") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        full = str(row.get("full_label", "")).strip()
                        short = str(row.get("short_label", "")).strip()
                        first_seen_raw = str(row.get("first_seen_iteration", "")).strip()
                        if not full or not short:
                            continue
                        try:
                            first_seen = int(first_seen_raw)
                        except ValueError:
                            first_seen = iteration
                        if full not in full_to_short:
                            full_to_short[full] = short
                            first_seen_by_full[full] = first_seen
                            used_short_labels.add(short)
            except Exception:
                full_to_short = {}
                first_seen_by_full = {}
                used_short_labels = set()

        for lab in labels:
            if lab in full_to_short:
                continue
            base = _shorten_label(lab)
            candidate = base
            if candidate in used_short_labels:
                digest = hashlib.blake2s(
                    str(lab).encode("utf-8"), digest_size=2
                ).hexdigest()
                suffix = digest.upper()
                candidate = f"{base}-{suffix}"
            if candidate in used_short_labels:
                i = 2
                while f"{candidate}{i}" in used_short_labels:
                    i += 1
                candidate = f"{candidate}{i}"
            full_to_short[lab] = candidate
            first_seen_by_full[lab] = iteration
            used_short_labels.add(candidate)

        try:
            with open(cumulative_mapping_path, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["short_label", "full_label", "first_seen_iteration"])
                for full in sorted(
                    full_to_short.keys(),
                    key=lambda x: (first_seen_by_full.get(x, iteration), x.lower()),
                ):
                    w.writerow(
                        [
                            full_to_short[full],
                            full,
                            first_seen_by_full.get(full, iteration),
                        ]
                    )
        except Exception as exc:
            self.logger.warning("Failed to write cumulative label mapping to %s: %s", cumulative_mapping_path, exc)

        for lab in labels:
            display_label = str(full_to_short.get(lab, lab))
            color = new_color if lab in new_nodes else old_color
            lines = display_label.splitlines() if display_label else [""]
            n_lines = max(1, len(lines))
            max_line_len = max((len(line) for line in lines), default=0)
            s = 400 + 45 * (max_line_len**1.15) + 220 * n_lines
            s = max(700, min(s, 8000))
            ax.scatter(
                xs[lab],
                ys[lab],
                s=s,
                c=color,
                edgecolors="#333333",
                linewidths=1.0,
                zorder=3,
            )
            ax.text(
                xs[lab],
                ys[lab],
                display_label,
                ha="center",
                va="center",
                fontsize=max(8, int(self.reaction_tree_fontsize)),
                color="black",
                zorder=4,
                wrap=True,
            )

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim(-0.1, 1.1)
        min_y = min(ys.values())
        max_y = max(ys.values())
        ax.set_ylim(min_y - 0.5, max_y + 0.5)

        ax.set_title(f"Reaction tree – iteration {iteration}", fontsize=14)

        from matplotlib.patches import Patch

        legend_handles = [
            Patch(facecolor=old_color, edgecolor="#333333", label="Existing core species"),
            Patch(facecolor=new_color, edgecolor="#333333", label="New in this iteration"),
        ]
        ax.legend(
            handles=legend_handles,
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
            borderaxespad=0.0,
            frameon=False,
            fontsize=9,
        )

        fig.tight_layout(rect=[0.0, 0.0, 0.78, 1.0])

        out_path = os.path.join(
            self.output_directory,
            f"reaction_tree_iter{iteration}.png",
        )
        fig.savefig(
            out_path,
            dpi=250,
            bbox_inches="tight",
            pad_inches=0.25,
            facecolor="white",
            edgecolor="none",
        )
        plt.close(fig)

        self.logger.info(
            f"Exported reaction tree plot for iteration {iteration} to {out_path}"
        )

    # ------------------------------------------------------------------
    # Export helpers - simulation plots
    # ------------------------------------------------------------------

    def export_simulation_plots(
        self,
        filename_pattern: str = "simulation_plot_iter{}.png",
    ) -> Optional[List[str]]:
        """
        Plot concentration vs time for each iteration and save to PNG.
        """
        if not self.profiles:
            return None

        enzyme_labels: Set[str] = set()
        if self.plot_exclude_enzymes and self.model:
            for sd in self.model.core_species + self.model.edge_species:
                if sd.is_enzyme:
                    enzyme_labels.add(sd.label)

        cofactor_labels: Set[str] = set()
        if self.plot_exclude_cofactors and self.bees_object is not None:
            for sp in getattr(self.bees_object, "species", []) or []:
                if getattr(sp, "reactive", True) is False:
                    cofactor_labels.add(getattr(sp, "label", ""))

        max_species = self.plot_max_species if self.plot_max_species is not None else 12
        exclude = enzyme_labels | cofactor_labels
        paths: List[str] = []

        for it_idx, prof in enumerate(self.profiles, 1):
            candidates = []
            for i, lab in enumerate(prof.species_labels):
                if lab in exclude:
                    continue
                if self.plot_exclude_cofactors and is_general_cofactor_label(lab):
                    continue

                y = prof.y[i, :]
                span = float(y.max() - y.min())
                if span <= 1e-9:
                    continue
                candidates.append((i, lab, span))

            if not candidates:
                continue

            candidates.sort(key=lambda item: item[2], reverse=True)
            if len(candidates) > max_species:
                candidates = candidates[:max_species]

            fig, ax = plt.subplots(figsize=(10, 6), facecolor="white")
            ax.set_facecolor("white")

            colors = [
                "#0173B2",
                "#DE8F05",
                "#029E73",
                "#CC78BC",
                "#CA9161",
                "#FBAFE4",
                "#949494",
                "#ECE133",
                "#56B4E9",
                "#D55E00",
            ]
            ax.set_prop_cycle(color=colors)

            mM_to_uM = 1000.0
            t = prof.t
            for i, lab, _ in candidates:
                ax.plot(t, prof.y[i, :] * mM_to_uM, label=lab)

            ax.set_xlabel("Time (s)", fontsize=14, fontweight="medium")
            ax.set_ylabel(r"Concentration ($\mu$M)", fontsize=14, fontweight="medium")
            ax.set_title(f"Iteration {it_idx}", fontsize=16, fontweight="medium")
            ax.tick_params(axis="both", which="major", labelsize=12)
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:g}"))
            ax.grid(True, alpha=0.35, linestyle="-", linewidth=0.6)
            ax.set_axisbelow(True)
            ax.legend(
                loc="upper left",
                bbox_to_anchor=(1.02, 1.0),
                borderaxespad=0.0,
                frameon=False,
                fontsize=11,
            )

            fig.tight_layout(rect=[0.0, 0.0, 0.78, 1.0])
            out_path = os.path.join(
                self.output_directory,
                filename_pattern.format(it_idx),
            )
            fig.savefig(
                out_path,
                dpi=150,
                bbox_inches="tight",
                pad_inches=0.25,
                facecolor="white",
                edgecolor="none",
            )
            plt.close(fig)
            paths.append(out_path)

        if paths:
            self.logger.info(
                f"Exported {len(paths)} simulation plot(s) to {self.output_directory}"
            )
        return paths if paths else None

    # ------------------------------------------------------------------
    # Export helpers - SBML (for COPASI / any SBML-compatible tool)
    # ------------------------------------------------------------------

    def export_sbml(
        self,
        filename: str = "model.xml",
        core_only: bool = True,
    ) -> Optional[str]:
        """
        Export the reaction network as an SBML Level 3 Version 2 file.

        Structure of the generated SBML:
        - One compartment: ``cytosol`` volume = 1 L, so concentrations in mM
          map directly to amounts in mmol).
        - One ``species`` per model species; ``initialConcentration`` is the
          value from the last iteration (mM).  Enzyme species are set
          ``constant=true, boundaryCondition=true`` so COPASI treats them as
          fixed parameters.
        - One ``parameter`` per unique enzyme holding its concentration (mM).
        - One ``reaction`` per core reaction (or core + edge if
          ``core_only=False``) with:
          - Explicit ``listOfReactants`` / ``listOfProducts`` stoichiometries.
          - A Michaelis-Menten ``kineticLaw`` written as a MathML formula:
            ``kcat * E * prod_i(S_i / (Km_i + S_i))``
          - All kinetic constants stored as local ``parameter`` elements
            inside the kineticLaw.

        Args:
            filename: Output file name inside the project output directory.
            core_only: If True (default) export only core reactions/species.
                       Set False to include edge species/reactions as well.

        Returns:
            Absolute path of the written SBML file, or None (if libsbml is
            not installed or the model is empty)
        """
        if not _LIBSBML_AVAILABLE:
            self.logger.warning(
                "python-libsbml is not installed – SBML export skipped. "
                "Install it with:  pip install python-libsbml"
            )
            return None

        reactions: List[GeneratedReaction] = list(self.model.core_reactions)
        species_list: List[SpeciesData] = list(self.model.core_species)
        if not core_only:
            reactions += list(self.model.edge_reactions)
            species_list += list(self.model.edge_species)

        if not reactions:
            self.logger.warning("SBML export: no reactions in model – skipping.")
            return None

        # ----------------------------------------------------------------
        # 1.  Create SBML document
        # ----------------------------------------------------------------
        doc = _libsbml.SBMLDocument(3, 2)
        model = doc.createModel()
        model.setId("BEES_model")
        model.setName("BEES auto-generated kinetic model")

        # ----------------------------------------------------------------
        # 1.  Unit Definitions
        # ----------------------------------------------------------------
        # Time units = seconds
        tu = model.createUnitDefinition()
        tu.setId("second")
        u = tu.createUnit()
        u.setKind(_libsbml.UNIT_KIND_SECOND)
        u.setExponent(1)
        u.setScale(0)
        u.setMultiplier(1.0)
        model.setTimeUnits("second")

        # Substance units = mmol
        su = model.createUnitDefinition()
        su.setId("mmol")
        u2 = su.createUnit()
        u2.setKind(_libsbml.UNIT_KIND_MOLE)
        u2.setExponent(1)
        u2.setScale(-3)
        u2.setMultiplier(1.0)
        model.setSubstanceUnits("mmol")
        model.setExtentUnits("mmol")

        # Volume units = litre
        vu = model.createUnitDefinition()
        vu.setId("litre")
        u3 = vu.createUnit()
        u3.setKind(_libsbml.UNIT_KIND_LITRE)
        u3.setExponent(1)
        u3.setScale(0)
        u3.setMultiplier(1.0)
        model.setVolumeUnits("litre")

        # ----------------------------------------------------------------
        # 2.  One compartment: cytosol (1 mL = 0.001 L)
        # ----------------------------------------------------------------
        comp = model.createCompartment()
        comp.setId("cytosol")
        comp.setName("Cytosol")
        comp.setSize(0.001)  # 1 mL
        comp.setConstant(True)
        comp.setSpatialDimensions(3)

        # ----------------------------------------------------------------
        # 3.  Species
        # ----------------------------------------------------------------
        # Build a safe SBML id from a label (SBML ids must start with letter/underscore)
        def _sbml_id(label: str) -> str:
            s = re.sub(r"[^A-Za-z0-9_]", "_", str(label).strip())
            if s and s[0].isdigit():
                s = "_" + s
            return s or "_species"

        label_to_id: Dict[str, str] = {}
        used_ids: Set[str] = set()
        for sd in species_list:
            base = _sbml_id(sd.label)
            sid = base
            n = 2
            while sid in used_ids:
                sid = f"{base}_{n}"
                n += 1
            label_to_id[sd.label.lower().strip()] = sid
            used_ids.add(sid)

        # Build initial concentration mapping from original input
        input_concs = {}
        if self.bees_object is not None:
            for sp in getattr(self.bees_object, "species", []) or []:
                input_concs[str(sp.label).lower().strip()] = getattr(sp, "concentration", 0.0)
            for enz in getattr(self.bees_object, "enzymes", []) or []:
                input_concs[str(enz.label).lower().strip()] = getattr(enz, "concentration", 0.0)

        for sd in species_list:
            lc = sd.label.lower().strip()
            sid = label_to_id[lc]
            sp = model.createSpecies()
            sp.setId(sid)
            sp.setName(sd.label)
            sp.setCompartment("cytosol")
            conc = max(input_concs.get(lc, 0.0) or 0.0, 0.0)
            sp.setInitialConcentration(conc)
            sp.setHasOnlySubstanceUnits(False)
            is_enz = bool(getattr(sd, "is_enzyme", False))
            is_const = bool(getattr(sd, "constant", False)) or is_enz
            sp.setConstant(is_const)
            sp.setBoundaryCondition(is_const)

        # ----------------------------------------------------------------
        # 4.  Global parameters: enzyme concentrations
        # ----------------------------------------------------------------
        enzyme_param_ids: Dict[str, str] = {}  # enzyme_label_lc -> param_id
        for sd in species_list:
            if not getattr(sd, "is_enzyme", False):
                continue
            lc = sd.label.lower().strip()
            pid = f"E_{_sbml_id(sd.label)}"
            n = 2
            orig_pid = pid
            while model.getParameter(pid) is not None:
                pid = f"{orig_pid}_{n}"
                n += 1
            p = model.createParameter()
            p.setId(pid)
            p.setName(f"[{sd.label}]")
            p.setValue(max(input_concs.get(lc, 0.0) or 0.0, 0.0))
            p.setConstant(True)
            enzyme_param_ids[lc] = pid

        # ----------------------------------------------------------------
        # 5.  Reactions
        # ----------------------------------------------------------------
        used_rxn_ids: Set[str] = set()

        for rxn_idx, rxn in enumerate(reactions, start=1):
            sig = reaction_signature(rxn)
            rxn_num = self.reaction_id_by_sig.get(sig, rxn_idx)
            rid = f"R{rxn_num}"
            if rid in used_rxn_ids:
                rid = f"R{rxn_num}_{rxn_idx}"
            used_rxn_ids.add(rid)

            sbml_rxn = model.createReaction()
            sbml_rxn.setId(rid)
            sbml_rxn.setName(
                f"{rxn.enzyme_label}: {' + '.join(rxn.reactant_labels)} -> {' + '.join(rxn.product_labels)}"
            )
            sbml_rxn.setReversible(False)

            # Reactants
            for r_label in rxn.reactant_labels:
                lc = r_label.lower().strip()
                sid = label_to_id.get(lc)
                if sid is None:
                    continue
                sr = sbml_rxn.createReactant()
                sr.setSpecies(sid)
                coeff = abs(rxn.stoichiometry.get(r_label, -1))
                sr.setStoichiometry(float(coeff))
                sr.setConstant(True)

            # Products
            for p_label in rxn.product_labels:
                lc = p_label.lower().strip()
                sid = label_to_id.get(lc)
                if sid is None:
                    continue
                sp2 = sbml_rxn.createProduct()
                sp2.setSpecies(sid)
                coeff = abs(rxn.stoichiometry.get(p_label, 1))
                sp2.setStoichiometry(float(coeff))
                sp2.setConstant(True)

            # KineticLaw: kcat * E * prod(S/(Km+S))
            kin = getattr(rxn, "kinetics", None)
            if kin is None or getattr(rxn, "rate_law", None) is None:
                continue

            kl = sbml_rxn.createKineticLaw()

            kcat_val = getattr(kin, "kcat", None)
            km_per = getattr(kin, "km_per_substrate", None) or {}
            km_single = getattr(kin, "km", None)

            # We will make them global parameters so they show up easily in COPASI's parameter list
            if kcat_val is not None:
                pid_kcat = f"kcat_{rid}"
                p_kcat = model.createParameter()
                p_kcat.setId(pid_kcat)
                p_kcat.setName(f"kcat ({rxn.enzyme_label})")
                p_kcat.setValue(float(kcat_val))
                p_kcat.setConstant(True)

            substrate_km_pairs: List[Tuple[str, str]] = []
            for r_label in rxn.reactant_labels:
                r_lc = r_label.lower().strip()
                km_val: Optional[float] = None
                if km_per:
                    km_val = km_per.get(r_label)
                    if km_val is None:
                        km_val = next(
                            (v for k, v in km_per.items() if k.lower().strip() == r_lc),
                            None,
                        )
                    if km_val is None:
                        continue  # saturated – factor = 1, no Km term
                else:
                    km_val = km_single
                if km_val is None or km_val <= 0:
                    continue
                sid = label_to_id.get(r_lc)
                if sid is None:
                    continue
                
                # Make Km a global parameter
                km_pid = f"Km_{rid}_{_sbml_id(r_label)}"
                p_km = model.createParameter()
                
                suffix_n = 2
                orig_km_pid = km_pid
                while model.getParameter(km_pid) is not None:
                    km_pid = f"{orig_km_pid}_{suffix_n}"
                    suffix_n += 1
                    
                p_km.setId(km_pid)
                p_km.setName(f"Km for {r_label} ({rxn.enzyme_label})")
                p_km.setValue(float(km_val))
                p_km.setConstant(True)
                
                substrate_km_pairs.append((sid, km_pid))

            # Build formula string
            e_lc = rxn.enzyme_label.lower().strip()
            e_param = enzyme_param_ids.get(e_lc)
            if e_param is None:
                # Enzyme not in species list – use a fallback numeric value
                e_conc = 0.001  # 1 µM default
                e_token = str(e_conc)
            else:
                e_token = e_param

            if kcat_val is not None and e_param is not None:
                formula_parts = [pid_kcat, e_token]
            elif kcat_val is not None:
                formula_parts = [pid_kcat, e_token]
            else:
                formula_parts = []

            for sid, km_pid in substrate_km_pairs:
                # Use abs(sid) in denominator to prevent division-by-zero if solver overshoots to negative
                formula_parts.append(f"({sid} / ({km_pid} + abs({sid})))")

            if formula_parts:
                formula = " * ".join(formula_parts)
            else:
                formula = "0"

            kl.setFormula(formula)

        # ----------------------------------------------------------------
        # 6.  Validate and write
        # ----------------------------------------------------------------
        doc.setConsistencyChecks(
            _libsbml.LIBSBML_CAT_GENERAL_CONSISTENCY, True
        )
        doc.setConsistencyChecks(
            _libsbml.LIBSBML_CAT_IDENTIFIER_CONSISTENCY, True
        )

        output_path = os.path.join(self.output_directory, filename)
        writer = _libsbml.SBMLWriter()
        writer.setProgramName("BEES")
        writer.setProgramVersion("0.1")
        ok = writer.writeSBMLToFile(doc, output_path)

        if ok:
            n_sp = model.getNumSpecies()
            n_rx = model.getNumReactions()
            self.logger.info(
                f"Exported SBML model to {output_path} "
                f"({n_sp} species, {n_rx} reactions) — load into COPASI or any SBML tool."
            )
            return output_path
        else:
            self.logger.warning(f"SBML export failed (libsbml writer returned error).")
            return None


