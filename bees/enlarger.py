#!/usr/bin/env python3

"""
Iterative Enlarger Module 
--------------------------------------
Orchestrates rate-based iterative model enlargement for biochemical
reaction networks. 
Algorithm overview:
    1. Initialise core species from user input.
    2. Generate reactions for the initial core species only (single pass).
    3. Run ODE simulation of the current model (core + edge). stops when the termination conditions are met.
    4. Evaluate edge species flux.
    5. Promote edge species whose |R_i| >= epsilon * R_char to the core.
    6. Generate new reactions for the newly promoted species.
    7. Repeat from step 3 until convergence or a termination criterion.
"""

import os
import time
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

from bees.core_edge_model import CoreEdgeModel, SpeciesData
from bees.exporter import EnlargerExporter, reaction_signature
from bees.flux_calculator import (
    identify_insignificant_species_from_peak_ratios,
    identify_significant_species_at_interrupt,
)
from bees.reaction_generator import GeneratedReaction, ModelGenerator
from bees.simulator import ODESimulator, SimulationResult


@dataclass
class EnlargerResult:
    """
    Summary returned by the iterative enlarger.

    Attributes:
        iterations: Number of enlargement iterations performed.
        converged: Whether convergence was achieved.
        convergence_reason: Reason for stopping.
        final_core_species: Number of core species at finish.
        final_edge_species: Number of edge species at finish.
        final_core_reactions: Number of core reactions at finish.
        simulation_profiles: List of SimulationResult from each iteration.
        model: The final CoreEdgeModel.

    Return:
    - Convergence: the model has reached the desired size or the termination conditions are met.
    - Termination: the model has reached the desired size or the termination conditions are met.
    - Error: the model has reached the desired size or the termination conditions are met.
    - Timeout: the model has reached the desired size or the termination conditions are met.
    - Max iterations: the model has reached the desired size or the termination conditions are met.
    - Max wall time: the model has reached the desired size or the termination conditions are met.
    - Max wall time: the model has reached the desired size or the termination conditions are met.
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
    enlargement engine.Discovery and simulation are interleaved: 
    reactions are generated only for the initial input species and then for each 
    batch of newly promoted species.  
   
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
        self.tol_interrupt_simulation: float = getattr(
            settings, "toleranceInterruptSimulation", None
        )
        if self.tol_interrupt_simulation is None:
            self.tol_interrupt_simulation = self.tol_move_to_core
        # RMG-style reaction-level criterion (dlnaccum). None disables it.
        _tol_rxn = getattr(settings, "toleranceMoveEdgeReactionToCore", None)
        self.tol_move_edge_reaction_to_core: Optional[float] = (
            float(_tol_rxn) if _tol_rxn is not None else None
        )
        self.min_edge_iterations_for_prune: int = int(
            getattr(settings, "minEdgeIterationsForPrune", 2) or 0
        )
        self.min_core_species_for_prune: int = int(
            getattr(settings, "minCoreSpeciesForPrune", 0) or 0
        )
        self.max_iterations: int = settings.max_iterations
        self.max_edge_species: Optional[int] = settings.max_edge_species
        self.max_num_objects_per_iter: int = int(
            getattr(settings, "max_num_objects_per_iter", 10)
        )
        self.abs_flux_floor: float = float(
            getattr(settings, "abs_flux_floor", 1e-12)
        )
        self.ode_method: str = getattr(settings, "ode_method", None) or "BDF"
        _ode_rtol = getattr(settings, "ode_rtol", None)
        _ode_atol = getattr(settings, "ode_atol", None)
        self.ode_rtol: float = float(_ode_rtol) if _ode_rtol is not None else 1e-8
        self.ode_atol: float = float(_ode_atol) if _ode_atol is not None else 1e-10
        _mwt = getattr(settings, "max_wall_time_per_iteration", None)
        self.max_wall_time_per_iteration: Optional[float] = (
            float(_mwt) if _mwt is not None else None
        )
        self.stepwise_heartbeat_interval: Optional[float] = getattr(
            settings, "stepwise_heartbeat_interval", None
        )
        self.termination_conversion: Optional[Dict[str, float]] = settings.termination_conversion
        self.termination_rate_ratio: Optional[float] = settings.termination_rate_ratio
        self.save_ode_equations: bool = getattr(
            settings, "save_ode_equations", False
        )
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
        self._edge_species_created_iter: Dict[str, int] = {}
        self._ingest_iteration: int = 0
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

        # Caches / backoff to avoid repeated slow no-hit searches.
        # Keys are normalized (ec_number, substrate_label_lc).
        self._rxn_nohit_cache: Set = set()
        self._enzyme_nohit_counts: Dict = {}
        self._skip_enzyme_ec_numbers: Set = set()

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def run(self) -> EnlargerResult:
        """
        Execute the iterative rate-based enlargement loop.
      
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
        initial_reactions = self._generate_reactions_for_species(
            initial_species
        )
        self._ingest_iteration = 0
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
            self.logger.info(f"--- Enlargement iteration {iteration} ---")
            result.iterations = iteration

            # Every pass starts from t=0 with initial concentrations.
            self.model.reset_concentrations_to_initial()

            simulator = ODESimulator(self.model, logger=self.logger)
            sim_result = simulator.simulate(
                end_time=self.end_time,
                time_step=self.time_step,
                method=self.ode_method,
                rtol=self.ode_rtol,
                atol=self.ode_atol,
                interrupt_simulation_tol=self.tol_interrupt_simulation,
                tol_move_edge_reaction_to_core=self.tol_move_edge_reaction_to_core,
                max_wall_time_s=self.max_wall_time_per_iteration,
                stepwise_heartbeat_interval_s=self.stepwise_heartbeat_interval,
            )
            self._profiles.append(sim_result)

            promoted_labels_combined: List[str] = []

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
                self.logger.warning(f"ODE simulation failed: {sim_result.message}")
                result.convergence_reason = f"ODE failure: {sim_result.message}"
                break

            # Identify edge reactions whose lifetime-max dlnaccum exceeded the
            # reaction-level tolerance (RMG-style chain-branching criterion).
            # We promote *species* that participate in those reactions; a single
            # promotion may activate multiple edge reactions via reclassify.
            rxn_promote_species: List[str] = []
            if (
                self.tol_move_edge_reaction_to_core is not None
                and self.tol_move_edge_reaction_to_core > 0
                and sim_result.max_edge_reaction_dlnaccum
            ):
                # Sort edge reactions by dlnaccum descending; take violators.
                violators = [
                    (sig, dln)
                    for sig, dln in sim_result.max_edge_reaction_dlnaccum.items()
                    if dln > self.tol_move_edge_reaction_to_core
                ]
                violators.sort(key=lambda kv: kv[1], reverse=True)
                edge_label_set = {sp.label for sp in self.model.edge_species}
                for sig, dln in violators[: self.max_num_objects_per_iter]:
                    rxn_obj = self._reaction_obj_by_sig.get(sig)
                    if rxn_obj is None:
                        continue
                    # Promote each edge-side species participating in the reaction.
                    for sp_label in list(rxn_obj.reactant_labels) + list(rxn_obj.product_labels):
                        if sp_label in edge_label_set and sp_label not in rxn_promote_species:
                            rxn_promote_species.append(sp_label)
                if rxn_promote_species and self.logger:
                    self.logger.info(
                        f"  Reaction-level (dlnaccum) promotion: {len(violators)} edge reaction(s) "
                        f"above tol={self.tol_move_edge_reaction_to_core}; "
                        f"promoting {len(rxn_promote_species)} associated species."
                    )

            # Interrupted pass: promote, enlarge, and continue to next iteration.
            if sim_result.simulation_interrupted:
                significant_sub = identify_significant_species_at_interrupt(
                    sim_result.interrupt_edge_rates,
                    sim_result.interrupt_char_rate,
                    self.tol_move_to_core,
                    max_objects=self.max_num_objects_per_iter,
                    abs_flux_floor=self.abs_flux_floor,
                )
                if not significant_sub and not rxn_promote_species:
                    # The interrupt fired but nothing new can be promoted (e.g. a
                    # reaction whose dlnaccum is high but all participants are already
                    # core).  Resume the simulation from the interrupt point rather
                    # than terminating — the interrupt criterion may not fire again
                    # once the system has evolved past this transient.
                    t_resume = float(sim_result.t[-1]) if len(sim_result.t) > 0 else 0.0
                    if t_resume < self.end_time * 0.99:
                        max_edge_abs = max(
                            (abs(r) for r in sim_result.interrupt_edge_rates.values()),
                            default=0.0,
                        )
                        if self.logger:
                            self.logger.info(
                                f"  Interrupt at t={t_resume:.3e} s produced no new promotable "
                                f"species (R_char={sim_result.interrupt_char_rate:.4e}, "
                                f"max |R_edge|={max_edge_abs:.4e}); resuming simulation."
                            )
                        # Re-run from the interrupt point to end_time using current
                        # model concentrations (already updated by the ODE result).
                        remaining_time = self.end_time - t_resume
                        sim_result2 = simulator.simulate(
                            end_time=remaining_time,
                            time_step=self.end_time / 100.0,
                            interrupt_simulation_tol=self.tol_interrupt_simulation,
                            tol_move_edge_reaction_to_core=self.tol_move_edge_reaction_to_core,
                        )
                        if sim_result2.success and not sim_result2.simulation_interrupted:
                            # Completed without another interrupt — treat as a
                            # full convergence-check pass (fall through to the
                            # non-interrupted path below by replacing sim_result).
                            sim_result = sim_result2
                        elif sim_result2.simulation_interrupted:
                            # Another interrupt fired — handle it next iteration.
                            sim_result = sim_result2
                        # If ODE failed on resume, fall through to the break below.
                        if not sim_result.success:
                            result.convergence_reason = f"ODE failure on resume: {sim_result.message}"
                            if self.logger:
                                self.logger.warning(f"  {result.convergence_reason}")
                            break
                        # Re-evaluate with the new sim_result.
                        significant_sub = identify_significant_species_at_interrupt(
                            sim_result.interrupt_edge_rates,
                            sim_result.interrupt_char_rate,
                            self.tol_move_to_core,
                            max_objects=self.max_num_objects_per_iter,
                            abs_flux_floor=self.abs_flux_floor,
                        )
                        rxn_promote_species = []
                        if (
                            self.tol_move_edge_reaction_to_core is not None
                            and self.tol_move_edge_reaction_to_core > 0
                            and sim_result.max_edge_reaction_dlnaccum
                        ):
                            violators2 = [
                                (sig, dln)
                                for sig, dln in sim_result.max_edge_reaction_dlnaccum.items()
                                if dln > self.tol_move_edge_reaction_to_core
                            ]
                            edge_label_set2 = {sp.label for sp in self.model.edge_species}
                            for sig, _dln in violators2[: self.max_num_objects_per_iter]:
                                rxn_obj = self._reaction_obj_by_sig.get(sig)
                                if rxn_obj is None:
                                    continue
                                for sp_label in list(rxn_obj.reactant_labels) + list(rxn_obj.product_labels):
                                    if sp_label in edge_label_set2 and sp_label not in rxn_promote_species:
                                        rxn_promote_species.append(sp_label)
                        if not significant_sub and not rxn_promote_species:
                            max_edge_abs = max(
                                (abs(r) for r in sim_result.interrupt_edge_rates.values()),
                                default=0.0,
                            )
                            result.convergence_reason = (
                                "ODE interrupted but no edge species exceeded "
                                f"threshold (R_char={sim_result.interrupt_char_rate:.4e}, "
                                f"max |R_edge|={max_edge_abs:.4e}, "
                                f"tol={self.tol_move_to_core})."
                            )
                            self.logger.warning(f"  {result.convergence_reason}")
                            break
                    else:
                        max_edge_abs = max(
                            (abs(r) for r in sim_result.interrupt_edge_rates.values()),
                            default=0.0,
                        )
                        result.convergence_reason = (
                            "ODE interrupted but no edge species exceeded "
                            f"threshold (R_char={sim_result.interrupt_char_rate:.4e}, "
                            f"max |R_edge|={max_edge_abs:.4e}, "
                            f"tol={self.tol_move_to_core})."
                        )
                        self.logger.warning(f"  {result.convergence_reason}")
                        break

                t_int = float(sim_result.t[-1]) if len(sim_result.t) > 0 else 0.0
                labels_csv = ", ".join(sf.label for sf in significant_sub[:3])
                if len(significant_sub) > 3:
                    labels_csv += f", ... (+{len(significant_sub) - 3} more)"
                self.logger.info(
                    f"  Interrupt at t={t_int:.6e} s: "
                    f"R_char={sim_result.interrupt_char_rate:.6e}, "
                    f"promoting {len(significant_sub)}: [{labels_csv}]"
                )
                for sf in significant_sub:
                    rr_s = (
                        f"{sf.normalized_rate:.4e}"
                        if sf.normalized_rate != float("inf")
                        else "inf"
                    )
                    self.logger.debug(
                        f"    {sf.label}: rr={rr_s}, |rate|={abs(sf.rate):.4e} mM/s"
                    )

                batch_promoted: List[str] = []
                for sf in significant_sub:
                    sp = self.model.promote_species_to_core(sf.label)
                    if sp is not None:
                        batch_promoted.append(sp.label)
                        promoted_labels_combined.append(sp.label)
                # Also promote species flagged by the reaction-level criterion.
                for sp_label in rxn_promote_species:
                    sp = self.model.promote_species_to_core(sp_label)
                    if sp is not None and sp.label not in batch_promoted:
                        batch_promoted.append(sp.label)
                        promoted_labels_combined.append(sp.label)
                if not batch_promoted:
                    result.convergence_reason = (
                        "Interrupt: significant flux reported but no edge species "
                        "could be promoted."
                    )
                    self.logger.warning(f"  {result.convergence_reason}")
                    break

                n_reclassed = self.model.reclassify_reactions()
                self.logger.info(
                    f"  Promoted {len(batch_promoted)} species to core, "
                    f"{n_reclassed} reactions moved edge->core."
                )

                new_rxns_sub = self._generate_reactions_for_species(
                    batch_promoted
                )
                self._ingest_iteration = iteration
                self._ingest_reactions(new_rxns_sub)
                self._sync_reaction_tracking(iteration=iteration)
                self.logger.info(
                    f"  Generated {len(new_rxns_sub)} new reaction(s); restarting ODE from t=0."
                )

            else:
                # Non-interrupted pass: reached end_time. Apply termination checks
                # and declare convergence.
                r_char_final = sim_result.final_char_rate
                r_char_max = sim_result.max_char_rate
                self.logger.info(
                    f"  R_char at t_end = {r_char_final:.6e} mM/s  |  "
                    f"R_char peak (this run) = {r_char_max:.6e} mM/s"
                )

                if self._check_conversion_termination(sim_result):
                    result.converged = True
                    result.convergence_reason = "Termination conversion target reached."
                    self.logger.info(result.convergence_reason)
                elif self._check_rate_ratio_termination(sim_result):
                    result.converged = True
                    result.convergence_reason = "Termination rate ratio reached."
                    self.logger.info(result.convergence_reason)
                else:
                    result.converged = True
                    result.convergence_reason = (
                        "Simulation reached end_time without exceeding "
                        "toleranceInterruptSimulation."
                    )
                    self.logger.info(result.convergence_reason)

                # Optional pruning (only meaningful on a full trajectory).
                if self.tol_keep_in_edge > 0:
                    n_core = len(self.model.core_species)
                    if n_core < self.min_core_species_for_prune:
                        self.logger.info(
                            "  Skipping edge prune: core species count "
                            f"{n_core} < minCoreSpeciesForPrune "
                            f"({self.min_core_species_for_prune})."
                        )
                    else:
                        ineligible_prune: Set[str] = set()
                        for sp in self.model.edge_species:
                            lc = sp.label.lower().strip()
                            birth = self._edge_species_created_iter.get(lc, 0)
                            if iteration - birth < self.min_edge_iterations_for_prune:
                                ineligible_prune.add(lc)
                        to_prune = identify_insignificant_species_from_peak_ratios(
                            sim_result.max_edge_rate_ratio,
                            sim_result.max_char_rate,
                            self.tol_keep_in_edge,
                            ineligible_for_prune=ineligible_prune,
                        )
                        n_pruned = self.model.prune_edge(to_prune)
                        if n_pruned:
                            self.logger.info(f"  Pruned {n_pruned} edge species.")

            promoted_labels_combined = list(dict.fromkeys(promoted_labels_combined))

            if (
                self.max_edge_species is not None
                and len(self.model.edge_species) > self.max_edge_species
            ):
                result.convergence_reason = (
                    f"Exceeded max_edge_species ({self.max_edge_species})."
                )
                self.logger.info(result.convergence_reason)
                if self.save_reaction_tree_plots:
                    exporter = EnlargerExporter(
                        model=self.model,
                        profiles=self._profiles,
                        output_directory=self.output_directory,
                        logger=self.logger,
                        reaction_id_by_sig=self._reaction_id_by_sig,
                        reaction_first_seen_iter=self._reaction_first_seen_iter,
                        reaction_core_enter_iter=self._reaction_core_enter_iter,
                        reaction_obj_by_sig=self._reaction_obj_by_sig,
                        iteration_summaries=self._iteration_summaries,
                        save_reaction_tree_plots=self.save_reaction_tree_plots,
                        save_simulation_plots=self.save_simulation_plots,
                        plot_max_species=self.plot_max_species,
                        plot_exclude_enzymes=self.plot_exclude_enzymes,
                        plot_exclude_cofactors=self.plot_exclude_cofactors,
                        reaction_tree_layout=self.reaction_tree_layout,
                        reaction_tree_rankdir=self.reaction_tree_rankdir,
                        reaction_tree_fontsize=self.reaction_tree_fontsize,
                        core_seen_labels=self._core_seen_labels,
                        bees_object=self.bees_object,
                    )
                    exporter.export_reaction_tree(
                        iteration=iteration,
                        promoted_labels=promoted_labels_combined,
                    )
                break

            self.logger.info(f"  Model status: {self.model.summary()}")
            s_it = self.model.summary()
            self._iteration_summaries.append({
                "iteration": iteration,
                "core_species": s_it["core_species"],
                "edge_species": s_it["edge_species"],
                "core_reactions": s_it["core_reactions"],
                "edge_reactions": s_it["edge_reactions"],
            })
            if self.save_reaction_tree_plots:
                exporter = EnlargerExporter(
                    model=self.model,
                    profiles=self._profiles,
                    output_directory=self.output_directory,
                    logger=self.logger,
                    reaction_id_by_sig=self._reaction_id_by_sig,
                    reaction_first_seen_iter=self._reaction_first_seen_iter,
                    reaction_core_enter_iter=self._reaction_core_enter_iter,
                    reaction_obj_by_sig=self._reaction_obj_by_sig,
                    iteration_summaries=self._iteration_summaries,
                    save_reaction_tree_plots=self.save_reaction_tree_plots,
                    save_simulation_plots=self.save_simulation_plots,
                    plot_max_species=self.plot_max_species,
                    plot_exclude_enzymes=self.plot_exclude_enzymes,
                    plot_exclude_cofactors=self.plot_exclude_cofactors,
                    reaction_tree_layout=self.reaction_tree_layout,
                    reaction_tree_rankdir=self.reaction_tree_rankdir,
                    reaction_tree_fontsize=self.reaction_tree_fontsize,
                    core_seen_labels=self._core_seen_labels,
                    bees_object=self.bees_object,
                )
                exporter.export_reaction_tree(
                    iteration=iteration,
                    promoted_labels=promoted_labels_combined,
                )

            if result.converged:
                break

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
        self,
        species_labels: List[str],
    ) -> List[GeneratedReaction]:
        """
        Use the ModelGenerator mudule to find new reactions involving the given
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
                _ecval = getattr(enzyme, "ecnumber", None) or ""
                if isinstance(_ecval, list):
                    ec = tuple(e.strip() for e in _ecval)
                else:
                    ec = _ecval.strip()
                if ec and ec in self._skip_enzyme_ec_numbers:
                    continue
                key = (ec, label.lower().strip())
                if ec and key in self._rxn_nohit_cache:
                    continue

                t0 = time.time()
                rxns = self.model_generator._generate_reactions(
                    enzyme,
                    substrate,
                    available_species_labels_lc=full_available,
                    provided_species_labels_lc=provided_lc,
                )
                dt = time.time() - t0
                new_reactions.extend(rxns)

                # If we got no reactions, remember it to avoid repeating expensive queries.
                if ec and not rxns:
                    self._rxn_nohit_cache.add(key)
                    self._enzyme_nohit_counts[ec] = self._enzyme_nohit_counts.get(ec, 0) + 1

                    # Backoff policy for extremely slow, consistently-unproductive enzymes.
                    # FabA (EC 4.2.1.60) is observed to be very expensive in fatty-acid projects.
                    if ec == "EC 4.2.1.60" and self._enzyme_nohit_counts[ec] >= 1 and dt > 30.0:
                        self._skip_enzyme_ec_numbers.add(ec)
                        self.logger.info(
                            f"Skipping further reaction-generation calls for {ec} (repeated no-hits; last call {dt:.1f}s)."
                        )

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
        sim_result: SimulationResult,
    ) -> bool:
        """
        RMG-style TerminationRateRatio: stop when R_char at t_end has fallen
        below a fraction of R_char peak over this simulation (system has
        slowed relative to peak activity).
        """
        if self.termination_rate_ratio is None:
            return False

        max_cr = sim_result.max_char_rate
        final_cr = sim_result.final_char_rate
        if max_cr <= 0.0:
            return False

        ratio = final_cr / max_cr
        if ratio < self.termination_rate_ratio:
            self.logger.info(
                f"  Char rate ratio (R_char at t_end / peak): {ratio:.6e} "
                f"< {self.termination_rate_ratio}"
            )
            return True
        return False


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
                    self._edge_species_created_iter[label] = self._ingest_iteration
            self.model.add_reaction(rxn)

    def _sync_reaction_tracking(self, iteration: int) -> None:
        """
        Synchronize reaction metadata from the current model state.
        Assigns stable IDs and tracks first-seen/core-entry iterations.
        """
        all_reactions = list(self.model.core_reactions) + list(self.model.edge_reactions)
        core_sigs = {
            reaction_signature(rxn) for rxn in self.model.core_reactions
        }
        for rxn in all_reactions:
            sig = reaction_signature(rxn)
            if sig not in self._reaction_id_by_sig:
                self._reaction_id_by_sig[sig] = self._next_reaction_id
                self._next_reaction_id += 1
            if sig not in self._reaction_first_seen_iter:
                self._reaction_first_seen_iter[sig] = iteration
            if sig not in self._reaction_obj_by_sig:
                self._reaction_obj_by_sig[sig] = rxn
            if sig in core_sigs and sig not in self._reaction_core_enter_iter:
                self._reaction_core_enter_iter[sig] = iteration

    
  