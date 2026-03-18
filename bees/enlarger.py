#!/usr/bin/env python3

"""
Iterative Enlarger Module (RMG-style)
--------------------------------------
Orchestrates rate-based iterative model enlargement for biochemical
reaction networks. This model discover and simulation are interleaved so that
reactions are generated only for species whose flux exceeds the promotion threshold.

Algorithm overview:
    1. Initialise core species from user input.
    2. Generate reactions for the initial core species only (single pass).
    3. Run ODE simulation of the current model (core + edge).
    4. Evaluate edge species flux.
    5. Promote edge species whose |R_i| >= epsilon * R_char to the core.
    6. Generate new reactions for the newly promoted species.
    7. Repeat from step 3 until convergence or a termination criterion.

"""

import os
import csv
import re
import shutil
import subprocess
import hashlib
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from bees.common import GENERAL_COFACTORS, get_ontology_equivalents
from bees.core_edge_model import CoreEdgeModel, SpeciesData
from bees.flux_calculator import (
    calculate_characteristic_rate,
    identify_insignificant_species,
    identify_significant_species,
)
from bees.model_generator import GeneratedReaction, ModelGenerator
from bees.simulator import ODESimulator, SimulationResult


@dataclass
class EnlargerResult:
    """
    Summary returned by the iterative enlarger.

    Attributes:
        iterations: Number of enlargement iterations performed.
        converged: Whether convergence was achieved.
        convergence_reason: Human-readable reason for stopping.
        final_core_species: Number of core species at finish.
        final_edge_species: Number of edge species at finish.
        final_core_reactions: Number of core reactions at finish.
        simulation_profiles: List of SimulationResult from each iteration.
        model: The final CoreEdgeModel.
    """
    iterations: int = 0
    converged: bool = False
    convergence_reason: str = ""
    final_core_species: int = 0
    final_edge_species: int = 0
    final_core_reactions: int = 0
    simulation_profiles: list = None  
    model: Optional[CoreEdgeModel] = None

    def __post_init__(self):
        if self.simulation_profiles is None:
            self.simulation_profiles = []


class IterativeEnlarger:
    """
    rate-based iterative model enlargement engine.

    Discovery and simulation are interleaved: reactions are generated
    only for the initial input species and then for each batch of
    newly promoted species.  This avoids the expensive upfront batch
    enumeration of all reachable reactions, focusing kinetics
    estimation (CatPred) on the flux-active part of the network.

    Args:
        bees_object: Validated InputBase from schema.
        model_generator: An already-initialised ModelGenerator with a
                         loaded kinetic database.
        logger: Logger instance.
        output_directory: Path for output files.
    """

    def __init__(
        self,
        bees_object,
        model_generator: ModelGenerator,
        logger,
        output_directory: str,
    ):
        self.bees_object = bees_object
        self.model_generator = model_generator
        self.logger = logger
        self.output_directory = output_directory

        # Settings shortcuts
        settings = bees_object.settings
        self.end_time: float = settings.end_time
        self.time_step: Optional[float] = settings.time_step
        self.tol_move_to_core: float = settings.toleranceMoveToCore
        self.tol_keep_in_edge: float = settings.toleranceKeepInEdge
        self.max_iterations: int = settings.max_iterations
        self.max_edge_species: Optional[int] = settings.max_edge_species
        self.termination_conversion: Optional[Dict[str, float]] = settings.termination_conversion
        self.termination_rate_ratio: Optional[float] = settings.termination_rate_ratio
        self.save_profiles: bool = settings.save_simulation_profiles
        self.save_edge: bool = settings.saveEdgeSpecies
        self.save_ode_equations: bool = getattr(
            settings, "save_ode_equations", False
        )
        self.filter_reactions: bool = settings.filter_reactions
        self.save_reaction_tree_plots: bool = getattr(
            settings, "save_reaction_tree_plots", False
        )
        self.save_simulation_plots: bool = getattr(
            settings, "save_simulation_plots", True
        )
        self.plot_max_species: Optional[int] = getattr(
            settings, "plot_max_species", None
        )
        self.plot_exclude_enzymes: bool = getattr(
            settings, "plot_exclude_enzymes", True
        )
        self.plot_exclude_cofactors: bool = getattr(
            settings, "plot_exclude_cofactors", True
        )
        self.reaction_tree_layout: str = getattr(
            settings, "reaction_tree_layout", "graphviz"
        )
        self.reaction_tree_rankdir: str = getattr(
            settings, "reaction_tree_rankdir", "TB"
        )
        self.reaction_tree_fontsize: int = int(
            getattr(settings, "reaction_tree_fontsize", 8) or 8
        )

        # Internal state
        self.model = CoreEdgeModel()
        self._profiles: List[SimulationResult] = []
        # Track which core species labels have already appeared in previous
        self._core_seen_labels: Set[str] = set()
        
        # 5 lines of code to track the reaction history.
        self._reaction_id_by_sig: Dict[Tuple[str, Tuple[str, ...], Tuple[str, ...]], int] = {}
        self._reaction_first_seen_iter: Dict[Tuple[str, Tuple[str, ...], Tuple[str, ...]], int] = {}
        self._reaction_core_enter_iter: Dict[Tuple[str, Tuple[str, ...], Tuple[str, ...]], Optional[int]] = {}
        self._reaction_obj_by_sig: Dict[Tuple[str, Tuple[str, ...], Tuple[str, ...]], GeneratedReaction] = {}
        self._next_reaction_id: int = 1

        # Per-iteration model size history for exporting summaries.
        # Each entry: {"iteration": int, "core_species": int, "edge_species": int, "core_reactions": int, "edge_reactions": int}
        self._iteration_summaries: List[Dict[str, int]] = []

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def run(self) -> EnlargerResult:
        """
        Execute the iterative rate-based enlargement loop.

        Discovery and simulation are interleaved: reactions are generated
        only for the initial core species and then for each batch of
        newly promoted species, avoiding a full upfront enumeration.

        Returns:
            EnlargerResult summarising the outcome.
        """
        self.logger.info("=" * 60)
        self.logger.info("RATE-BASED ITERATIVE MODEL ENLARGEMENT")
        self.logger.info("=" * 60)

        # Step 0 -- Initialise core/edge from input
        self._initialise_model()
        # Record initial core species labels before any enlargement.
        self._core_seen_labels = {
            s.label for s in self.model.core_species
        }

        # Ensure the kinetics estimator is ready before we start calling _generate_reactions directly.
        self.model_generator.ensure_estimator_initialized()

        # Step 0b -- Generate initial reactions for the input species only
    
        initial_species = [
            sp.label for sp in self.bees_object.species
            if sp.reactive and not getattr(sp, "solvent", False)
        ]
        initial_reactions = self._generate_reactions_for_species(initial_species)
        self._ingest_reactions(initial_reactions)
        self._sync_reaction_tracking(iteration=0)
        self.logger.info(
            f"Initial discovery: {len(initial_reactions)} reaction(s) "
            f"for {len(initial_species)} input species."
        )

        self.logger.info(
            f"Initial model: {self.model.summary()}"
        )
        # Record iteration 0  (after initial discovery)
        s0 = self.model.summary()
        self._iteration_summaries = [{
            "iteration": 0,
            "core_species": s0["core_species"],
            "edge_species": s0["edge_species"],
            "core_reactions": s0["core_reactions"],
            "edge_reactions": s0["edge_reactions"],
        }]

        result = EnlargerResult(model=self.model)

        for iteration in range(1, self.max_iterations + 1):
            self.logger.info("-" * 40)
            self.logger.info(f"Enlargement iteration {iteration}")
            # Track current iteration early so summary is correct even if we break.
            result.iterations = iteration

            # 1. Run ODE simulation
            simulator = ODESimulator(self.model, logger=self.logger)
            sim_result = simulator.simulate(
                end_time=self.end_time,
                time_step=self.time_step,
            )
            self._profiles.append(sim_result)

            if self.save_ode_equations and sim_result.success:
                ode_path = os.path.join(
                    self.output_directory,
                    f"ode_equations_iter{iteration}.txt",
                )
                simulator.export_ode_equations(
                    output_path=ode_path,
                    iteration=iteration,
                    end_time=self.end_time,
                )

            if not sim_result.success:
                self.logger.warning(
                    f"ODE simulation failed: {sim_result.message}"
                )
                result.convergence_reason = f"ODE failure: {sim_result.message}"
                break

            # 2. Evaluate core and edge rates
            core_rates = simulator.evaluate_core_rates(sim_result)
            edge_rates = simulator.evaluate_edge_rates(sim_result)
            r_char = calculate_characteristic_rate(core_rates)

            self.logger.info(f"  R_char = {r_char:.6e} mM/s")

            # 3. Check termination -- conversion
            if self._check_conversion_termination(sim_result):
                result.converged = True
                result.convergence_reason = "Termination conversion target reached."
                self.logger.info(result.convergence_reason)
                break

            # 4. Check termination -- rate ratio
            if self._check_rate_ratio_termination(core_rates, edge_rates, r_char):
                result.converged = True
                result.convergence_reason = "Termination rate ratio reached."
                self.logger.info(result.convergence_reason)
                break

            # 5. Identify significant edge species
            significant = identify_significant_species(
                edge_rates, r_char, self.tol_move_to_core
            )

            if not significant:
                result.converged = True
                result.convergence_reason = (
                    "No edge species exceeded toleranceMoveToCore."
                )
                self.logger.info(result.convergence_reason)
                break

            self.logger.info(
                f"  {len(significant)} edge species exceed threshold:"
            )
            for sf in significant:
                self.logger.info(
                    f"    {sf.label}: |rate|={abs(sf.rate):.4e} mM/s "
                    f"(norm={sf.normalized_rate:.4e})"
                )

            # 6. Promote significant species to core
            promoted_labels = []
            for sf in significant:
                sp = self.model.promote_species_to_core(sf.label)
                if sp is not None:
                    promoted_labels.append(sp.label)

            # 7. Reclassify reactions (edge -> core if all participants now core)
            n_promoted_rxn = self.model.reclassify_reactions()
            self.logger.info(
                f"  Promoted {len(promoted_labels)} species, "
                f"{n_promoted_rxn} reactions moved to core."
            )

            # 8. Generate new reactions for promoted species
            new_reactions = self._generate_reactions_for_species(promoted_labels)
            self._ingest_reactions(new_reactions)
            self._sync_reaction_tracking(iteration=iteration)
            self.logger.info(
                f"  Generated {len(new_reactions)} new reaction(s)."
            )

            # 9. Optional: prune insignificant edge species
            if self.tol_keep_in_edge > 0:
                to_prune = identify_insignificant_species(
                    edge_rates, r_char, self.tol_keep_in_edge
                )
                n_pruned = self.model.prune_edge(to_prune)
                if n_pruned:
                    self.logger.info(f"  Pruned {n_pruned} edge species.")

            # 10. Check max edge species
            if (
                self.max_edge_species is not None
                and len(self.model.edge_species) > self.max_edge_species
            ):
                result.convergence_reason = (
                    f"Exceeded max_edge_species ({self.max_edge_species})."
                )
                self.logger.info(result.convergence_reason)
                self._export_reaction_tree(iteration, promoted_labels)
                break

            self.logger.info(f"  Model size: {self.model.summary()}")
            # Record end-of-iteration snapshot.
            s_it = self.model.summary()
            self._iteration_summaries.append({
                "iteration": iteration,
                "core_species": s_it["core_species"],
                "edge_species": s_it["edge_species"],
                "core_reactions": s_it["core_reactions"],
                "edge_reactions": s_it["edge_reactions"],
            })
            self._export_reaction_tree(iteration, promoted_labels)

        else:
            # Loop exhausted without break
            result.convergence_reason = (
                f"Reached max_iterations ({self.max_iterations})."
            )
            self.logger.info(result.convergence_reason)

        # Finalise result
        summary = self.model.summary()
        result.final_core_species = summary["core_species"]
        result.final_edge_species = summary["edge_species"]
        result.final_core_reactions = summary["core_reactions"]
        result.simulation_profiles = self._profiles
        result.model = self.model
        if result.iterations == 0:
            result.iterations = min(
                self.max_iterations, max(1, len(self._profiles))
            )

        self.logger.info("=" * 60)
        self.logger.info("ENLARGEMENT SUMMARY")
        self.logger.info(f"  Iterations:      {result.iterations}")
        self.logger.info(f"  Converged:       {result.converged}")
        self.logger.info(f"  Reason:          {result.convergence_reason}")
        self.logger.info(f"  Core species:    {result.final_core_species}")
        self.logger.info(f"  Edge species:    {result.final_edge_species}")
        self.logger.info(f"  Core reactions:  {result.final_core_reactions}")
        self.logger.info("=" * 60)

        return result


      # ------------------------------------------------------------------
    # Generate reactions for species
    # ------------------------------------------------------------------

    
    def _generate_reactions_for_species(
        self, species_labels: List[str]
    ) -> List[GeneratedReaction]:
        """
        Use the ModelGenerator to find new reactions involving the given
        species as substrates (crossed with all reactive enzymes).
        """
        from types import SimpleNamespace
        from bees.common import GENERAL_COFACTORS, get_ontology_equivalents

        new_reactions: List[GeneratedReaction] = []
        enzymes = [e for e in self.bees_object.enzymes if e.reactive]
        available_lc = self.model.get_species_labels_lc()

        full_available: Set[str] = set()
        for lc in available_lc:
            full_available.update(get_ontology_equivalents(lc))

        provided_lc = set(available_lc)

        for label in species_labels:
            if label.lower().strip() in GENERAL_COFACTORS:
                continue
            substrate = SimpleNamespace(
                label=label,
                reactive=True,
                solvent=False,
                smiles=None,
            )
            for enzyme in enzymes:
                rxns = self.model_generator._generate_reactions(
                    enzyme,
                    substrate,
                    available_species_labels_lc=full_available,
                    provided_species_labels_lc=provided_lc,
                )
                new_reactions.extend(rxns)

        return new_reactions

    # ------------------------------------------------------------------
    # Termination checks
    # ------------------------------------------------------------------

    def _check_conversion_termination(
        self, sim_result: SimulationResult
    ) -> bool:
        """
        Check if any termination_conversion target has been met.

        Conversion is defined as: X = 1 - C(t) / C(0).
        """
        if not self.termination_conversion:
            return False

        if sim_result.y.shape[1] == 0:
            return False

        label_to_idx = {
            lab.lower().strip(): i
            for i, lab in enumerate(sim_result.species_labels)
        }
        for species_label, target_frac in self.termination_conversion.items():
            lc = species_label.lower().strip()
            idx = label_to_idx.get(lc)
            if idx is None:
                continue
            c0 = sim_result.y[idx, 0]
            if c0 <= 0:
                continue
            c_final = sim_result.y[idx, -1]
            conversion = 1.0 - (max(c_final, 0.0) / c0)
            if conversion >= target_frac:
                self.logger.info(
                    f"  Conversion of '{species_label}': "
                    f"{conversion:.4f} >= {target_frac}"
                )
                return True
        return False

    def _check_rate_ratio_termination(
        self,
        core_rates: Dict[str, float],
        edge_rates: Dict[str, float],
        r_char: float,
    ) -> bool:
        """
        Check if the maximum edge-to-core rate ratio is below the
        termination_rate_ratio threshold.
        """
        if self.termination_rate_ratio is None or r_char <= 0:
            return False

        if not edge_rates:
            return True  # no edge flux at all

        max_edge_rate = max(abs(r) for r in edge_rates.values()) if edge_rates else 0.0
        ratio = max_edge_rate / r_char
        if ratio < self.termination_rate_ratio:
            self.logger.info(
                f"  Max edge/core rate ratio: {ratio:.6e} "
                f"< {self.termination_rate_ratio}"
            )
            return True
        return False


    # ------------------------------------------------------------------
    # Export helpers function - flux analysis
    # ------------------------------------------------------------------
    @staticmethod
    def _is_general_cofactor_label(label: str) -> bool:
        """
        Return True if label appears to be a general cofactor/carrier species.
        """
        normalized = " ".join(
            str(label).lower().strip().replace("_", " ").replace("-", " ").split()
        )
        compact = normalized.replace(" ", "")
        if normalized in GENERAL_COFACTORS or compact in GENERAL_COFACTORS:
            return True

        # Use ontology aliases to catch naming variants (e.g., alpha/beta synonyms).
        equivalents = get_ontology_equivalents(normalized)
        for eq in equivalents:
            eq_norm = " ".join(str(eq).lower().strip().split())
            if eq_norm in GENERAL_COFACTORS or eq_norm.replace(" ", "") in GENERAL_COFACTORS:
                return True
        return False

    def export_flux_analysis(
        self, filename: str = "flux_analysis.csv"
    ) -> Optional[str]:
        """
        Write per-iteration flux data (species rates, R_char, promotion decisions).
        """
       
        output_path = os.path.join(self.output_directory, filename)
        with open(output_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "iteration", "core_species", "edge_species",
                "core_reactions", "edge_reactions",
            ])

            # Prefer recorded history (like RMG: show core/edge growth over time).
            if self._iteration_summaries:
                for row in self._iteration_summaries:
                    writer.writerow([
                        row["iteration"],
                        row["core_species"],
                        row["edge_species"],
                        row["core_reactions"],
                        row["edge_reactions"],
                    ])
            else:
                # Fallback: final snapshot only.
                s = self.model.summary()
                writer.writerow([
                    len(self._profiles),
                    s["core_species"], s["edge_species"],
                    s["core_reactions"], s["edge_reactions"],
                ])

        self.logger.info(f"Exported flux analysis to {output_path}")

        # Detailed per-reaction summary (core/edge membership and iteration history).
        details_path = os.path.join(self.output_directory, "flux_analysis_reactions.csv")
        core_sigs = {
            self._reaction_signature(rxn) for rxn in self.model.core_reactions
        }
        edge_sigs = {
            self._reaction_signature(rxn) for rxn in self.model.edge_reactions
        }
        all_sigs = sorted(
            set(core_sigs) | set(edge_sigs),
            key=lambda sig: self._reaction_id_by_sig.get(sig, 10**9),
        )
        with open(details_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "reaction_id",
                "classification",
                "first_seen_iteration",
                "core_enter_iteration",
            ])
            for sig in all_sigs:
                rxn = self._reaction_obj_by_sig.get(sig)
                if rxn is None:
                    continue
                reaction_id = f"R{self._reaction_id_by_sig.get(sig, 0)}"
                classification = (
                    "core" if sig in core_sigs else "edge"
                )
                first_seen = self._reaction_first_seen_iter.get(sig, "")
                core_enter = self._reaction_core_enter_iter.get(sig, "")
                writer.writerow([
                    reaction_id,
                    classification,
                    first_seen,
                    core_enter if core_enter is not None else "",
                ])
        self.logger.info(f"Exported reaction flux summary to {details_path}")
        return output_path

# ------------------------------------------------------------------
# Export helpers function - simulation profiles
# ------------------------------------------------------------------
    def export_simulation_profiles(
        self, filename: str = "simulation_profiles.csv"
    ) -> Optional[str]:
        """
        Write concentration time-series to CSV.

        Concatenates profiles from every iteration into a single file.
        Returns the output path, or None if there are no profiles.
        """
        if not self._profiles:
            return None

        output_path = os.path.join(self.output_directory, filename)
        # Gather the union of species labels across iterations
        all_labels: List[str] = []
        seen: Set[str] = set()
        for prof in self._profiles:
            for lab in prof.species_labels:
                if lab not in seen:
                    all_labels.append(lab)
                    seen.add(lab)

        with open(output_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["iteration", "time"] + all_labels)
            for it_idx, prof in enumerate(self._profiles, 1):
                label_to_row = {
                    lab: i for i, lab in enumerate(prof.species_labels)
                }
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
    # Export helpers function - visualisation tools
    # ------------------------------------------------------------------

    def _export_reaction_tree(
        self,
        iteration: int,
        promoted_labels: Optional[List[str]] = None,
    ) -> None:
        """
        Export a reaction tree plot for this iteration.

        The tree:
        - Includes **only core species**.
        - Excludes enzymes and general cofactors.
        - Colours newly promoted core species in this iteration differently
          from species that were already in the core.
        """
        if not self.save_reaction_tree_plots or not self.model:
            return

        # Build the set of candidate nodes: core species that are not enzymes
        # and are not general cofactors.
        core_nodes: List[SpeciesData] = []
        for sd in self.model.core_species:
            if sd.is_enzyme:
                continue
            if self._is_general_cofactor_label(sd.label):
                continue
            core_nodes.append(sd)

        if not core_nodes:
            return

        labels = [sd.label for sd in core_nodes]
        label_set = set(labels)

        # Determine which nodes are "new" in this iteration.
        promoted_set = set(promoted_labels or [])
        new_nodes: Set[str] = set()
        for lab in labels:
            if lab in promoted_set and lab not in self._core_seen_labels:
                new_nodes.add(lab)

        # Update the cumulative core history.
        self._core_seen_labels.update(labels)

        # Build directed edges between species based on core reactions:
        # reactants -> products, filtered to core, non-enzyme, non-cofactor species.
        edges: List[tuple] = []
        edge_reaction_labels: Dict[tuple, Set[str]] = {}
        for rxn in self.model.core_reactions:
            sig = self._reaction_signature(rxn)
            reaction_id = self._reaction_id_by_sig.get(sig)
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

        # Draw each species-to-species edge once.
        edges = sorted(set(edges))

        if not edges:
            return

        # Assign coordinates. Prefer Graphviz (structure-aware); fall back to
        # the simple 2-row layout if Graphviz is unavailable.
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

            # Build a DOT graph with safe node IDs (Graphviz IDs cannot contain
            # arbitrary punctuation reliably). Keep mapping for coordinates.
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

            # Parse Graphviz plain output.
            # We care about: `node <name> <x> <y> <w> <h> ...`
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

        # Create the plot.
        fig, ax = plt.subplots(figsize=(16, 12), facecolor="white")
        ax.set_facecolor("white")

        # Draw edges.
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
            # Show reaction IDs on edges (e.g., R2, R43). Aggregate all
            # reactions that map to the same species-to-species connection.
            rid_set = edge_reaction_labels.get((src, dst), set())
            if rid_set:
                rid_text = ",".join(
                    sorted(
                        rid_set,
                        key=lambda s: int(s[1:]) if s.startswith("R") and s[1:].isdigit() else 10**9,
                    )
                )
                mid_x = (xs[src] + xs[dst]) / 2.0
                mid_y = (ys[src] + ys[dst]) / 2.0
                # Place label slightly above the edge so the line doesn't cross the text.
                ax.annotate(
                    rid_text,
                    xy=(mid_x, mid_y),
                    xycoords="data",
                    xytext=(0, 7),  # pixels/points offset upward
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

        # Draw nodes: old vs new in different colours.
        old_color = "#A6CEE3"  # blue
        new_color = "#FB9A99"  # red

        # Short labels for readability (full mapping exported alongside PNG).
        def _shorten_label(label: str) -> str:
            s = str(label).strip()
            s = re.sub(r"\s+", " ", s)
            # Keep common small molecules readable as-is.
            if len(s) <= 10 and " " not in s:
                return s
            # Use a token-based abbreviation.
            tokens = re.split(r"[\s\-_]+", s)
            keep = []
            for t in tokens:
                if not t:
                    continue
                # Keep numbers and short chemical fragments.
                if t.isdigit() or re.fullmatch(r"\d+[A-Za-z]*", t or ""):
                    keep.append(t)
                else:
                    keep.append(t[0].upper())
            base = "".join(keep) or s[:6].upper()
            # If still too long, truncate.
            if len(base) > 12:
                base = base[:12]
            return base

        # Build or update a single cumulative mapping:
        # full_label -> (short_label, first_seen_iteration)
        cumulative_mapping_path = os.path.join(
            self.output_directory,
            "reaction_tree_labels.csv",
        )
        full_to_short: Dict[str, str] = {}
        first_seen_by_full: Dict[str, int] = {}
        used_short_labels: Set[str] = set()

        # Load existing cumulative mapping if present.
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
                        # Keep first occurrence if duplicate rows exist.
                        if full not in full_to_short:
                            full_to_short[full] = short
                            first_seen_by_full[full] = first_seen
                            used_short_labels.add(short)
            except Exception:
                # Continue with a fresh mapping if the file is malformed.
                full_to_short = {}
                first_seen_by_full = {}
                used_short_labels = set()

        # Assign short labels for current iteration nodes.
        for lab in labels:
            if lab in full_to_short:
                continue
            base = _shorten_label(lab)
            candidate = base
            # Ensure uniqueness against all previously assigned short labels.
            if candidate in used_short_labels:
                digest = hashlib.blake2s(str(lab).encode("utf-8"), digest_size=2).hexdigest()
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

        # Persist cumulative mapping.
        try:
            with open(cumulative_mapping_path, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["short_label", "full_label", "first_seen_iteration"])
                for full in sorted(
                    full_to_short.keys(),
                    key=lambda x: (first_seen_by_full.get(x, iteration), x.lower()),
                ):
                    w.writerow([
                        full_to_short[full],
                        full,
                        first_seen_by_full.get(full, iteration),
                    ])
        except Exception:
            # Plot should still export even if mapping write fails.
            pass

        for lab in labels:
            display_label = str(full_to_short.get(lab, lab))
            color = new_color if lab in new_nodes else old_color
            # Adaptive node size: scale marker area by wrapped label size so
            # text fits inside the circle more often (simple heuristic).
            lines = display_label.splitlines() if display_label else [""]
            n_lines = max(1, len(lines))
            max_line_len = max((len(line) for line in lines), default=0)
            # Matplotlib scatter uses marker area in points^2.
            s = 400 + 45 * (max_line_len ** 1.15) + 220 * n_lines
            # Guardrails for readability.
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
        # Add a little vertical padding.
        min_y = min(ys.values())
        max_y = max(ys.values())
        ax.set_ylim(min_y - 0.5, max_y + 0.5)

        ax.set_title(f"Reaction tree – iteration {iteration}", fontsize=14)

        # Legend explaining colours.
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


    def export_simulation_plots(
        self,
        filename_pattern: str = "simulation_plot_iter{}.png",
    ) -> Optional[List[str]]:
        """
        Plot concentration vs time for each iteration and save to PNG.

        Respects plot_exclude_enzymes, plot_exclude_cofactors, and plot_max_species.
        Returns list of written file paths, or None if no profiles.
        """
        if not self._profiles:
            return None

        enzyme_labels: Set[str] = set()
        if self.plot_exclude_enzymes and self.model:
            for sd in self.model.core_species + self.model.edge_species:
                if sd.is_enzyme:
                    enzyme_labels.add(sd.label)

        cofactor_labels: Set[str] = set()
        if self.plot_exclude_cofactors and self.bees_object:
            for sp in getattr(self.bees_object, "species", []) or []:
                if getattr(sp, "reactive", True) is False:
                    cofactor_labels.add(getattr(sp, "label", ""))

        # Keep plots readable by default even when plot_max_species is unset.
        # Users can still override this via plot_max_species in settings.
        max_species = self.plot_max_species if self.plot_max_species is not None else 12

        exclude = enzyme_labels | cofactor_labels
        paths: List[str] = []

        for it_idx, prof in enumerate(self._profiles, 1):
            candidates = []
            for i, lab in enumerate(prof.species_labels):
                if lab in exclude:
                    continue
                if self.plot_exclude_cofactors and self._is_general_cofactor_label(lab):
                    continue

                y = prof.y[i, :]
                # Drop practically flat series to avoid legend clutter.
                span = float(y.max() - y.min())
                if span <= 1e-9:
                    continue
                candidates.append((i, lab, span))

            if not candidates:
                continue

            # Prioritize species with the largest dynamic span.
            candidates.sort(key=lambda item: item[2], reverse=True)
            if len(candidates) > max_species:
                candidates = candidates[:max_species]

            # Publication-ready style: clean white background, large fonts, µM axis
            fig, ax = plt.subplots(figsize=(10, 6), facecolor="white")
            ax.set_facecolor("white")

            # Colorblind-friendly palette (Paul Tol / matplotlib tab10–style)
            colors = [
                "#0173B2", "#DE8F05", "#029E73", "#CC78BC", "#CA9161",
                "#FBAFE4", "#949494", "#ECE133", "#56B4E9", "#D55E00",
            ]
            ax.set_prop_cycle(color=colors)

            # Concentration in µM for cleaner tick labels (data are in mM)
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

            # Clean layout: legend outside, no black border, minimal padding
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
    # Internal helpers
    # ------------------------------------------------------------------

    def _initialise_model(self) -> None:
        """
        Populate the core/edge model from the validated input object.

        Note: All user-provided species and enzymes go into the core.
        """
        # Add substrates / species to core
        for sp in self.bees_object.species:
            conc = sp.concentration
            if isinstance(conc, tuple):
                conc = conc[0]  # use lower bound of range
            sd = SpeciesData(
                label=sp.label,
                concentration=conc,
                initial_concentration=conc,
                is_enzyme=False,
                constant=getattr(sp, "constant", False),
            )
            self.model.add_core_species(sd)

        # Add enzymes to core
        for enz in self.bees_object.enzymes:
            conc = enz.concentration
            if isinstance(conc, tuple):
                conc = conc[0]
            sd = SpeciesData(
                label=enz.label,
                concentration=conc,
                initial_concentration=conc,
                is_enzyme=True,
                constant=True, 
            )
            self.model.add_core_species(sd)

    def _ingest_reactions(self, reactions: List[GeneratedReaction]) -> None:
        """
        Add reactions to the core/edge model, creating SpeciesData entries
        for any previously unknown products (which start in the edge).
        """
        for rxn in reactions:
            # Ensure all participants are tracked
            for label in rxn.reactant_labels + rxn.product_labels:
                label= label.lower().strip()
                if (
                    not self.model.is_core_species(label)
                    and not self.model.is_edge_species(label)
                ):
                    # New species discovered -- add to edge with 0 concentration
                    self.model.add_edge_species(SpeciesData(
                        label=label,
                        concentration=0.0,
                        initial_concentration=0.0,
                    ))
            self.model.add_reaction(rxn)

    @staticmethod
    def _reaction_signature(
        reaction: GeneratedReaction,
    ) -> Tuple[str, Tuple[str, ...], Tuple[str, ...]]:
        """
        Canonical reaction signature for stable reaction ID tracking.
        """
        enzyme = str(reaction.enzyme_label).lower().strip()
        reactants = tuple(
            sorted(str(r).lower().strip() for r in reaction.reactant_labels)
        )
        products = tuple(
            sorted(str(p).lower().strip() for p in reaction.product_labels)
        )
        return (enzyme, reactants, products)

    def _sync_reaction_tracking(self, iteration: int) -> None:
        """
        Synchronize reaction metadata from the current model state.
        Assigns stable IDs and tracks first-seen/core-entry iterations.
        """
        all_reactions = list(self.model.core_reactions) + list(self.model.edge_reactions)
        core_sigs = {
            self._reaction_signature(rxn) for rxn in self.model.core_reactions
        }
        for rxn in all_reactions:
            sig = self._reaction_signature(rxn)
            if sig not in self._reaction_id_by_sig:
                self._reaction_id_by_sig[sig] = self._next_reaction_id
                self._next_reaction_id += 1
            if sig not in self._reaction_first_seen_iter:
                self._reaction_first_seen_iter[sig] = iteration
            if sig not in self._reaction_obj_by_sig:
                self._reaction_obj_by_sig[sig] = rxn
            if sig in core_sigs and sig not in self._reaction_core_enter_iter:
                self._reaction_core_enter_iter[sig] = iteration

    
  