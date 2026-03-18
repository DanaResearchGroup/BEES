#!/usr/bin/env python3

"""
ODE Simulator Module
--------------------
Integrates the reaction network over time using scipy's ODE solvers.

Supports Michaelis-Menten and Reversible-Michaelis-Menten rate laws.
The simulator integrates both core and edge species and reactions so
that flux can be computed for edge species (enabling promotion decisions).
Edge species start at zero concentration.

All concentrations are in mM, time in seconds, rates in mM/s.
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import os

import numpy as np
from scipy.integrate import solve_ivp

from bees.common import get_ontology_equivalents
from bees.core_edge_model import CoreEdgeModel
from bees.flux_calculator import compute_mm_rate


@dataclass
class SimulationResult:
    """
    Container for ODE simulation output.

    Attributes:
        t: Array of time points (seconds).
        y: 2-D array of shape (n_species, n_timepoints).
        species_labels: Ordered list of species labels matching rows of y.
        success: Whether the integration succeeded.
        message: Solver message.
    """
    t: np.ndarray
    y: np.ndarray
    species_labels: List[str]
    success: bool = True
    message: str = ""


class ODESimulator:
    """
    Integrate the biochemical reaction network for a CoreEdgeModel.

    Simulates both core and edge species and reactions so that flux can
    be computed for edge species. Edge species start at zero concentration.

    Parameters
    ----------
    model : CoreEdgeModel
        The core/edge model to simulate.
    logger : optional
        Logger instance for diagnostic messages.
    """

    def __init__(self, model: CoreEdgeModel, logger=None):
        self.model = model
        self.logger = logger

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def simulate(
        self,
        end_time: float,
        time_step: Optional[float] = None,
        method: str = "BDF",
        rtol: float = 1e-8,
        atol: float = 1e-10,
    ) -> SimulationResult:
        """
        Run an ODE simulation of the full model (core + edge species and reactions).

        Edge species start at zero concentration. This enables flux to be
        computed for edge species so they can be promoted to the core.

        Args:
            end_time: Simulation end time in seconds.
            time_step: If given, store solution at these intervals.
            method: Integration method for solve_ivp (default BDF,
                    good for stiff biochemical systems).
            rtol: Relative tolerance for the solver.
            atol: Absolute tolerance for the solver.

        Returns:
            SimulationResult with time series data.
        """
        species_labels = self.model.get_all_species_labels()
        n_species = len(species_labels)

        if n_species == 0:
            return SimulationResult(
                t=np.array([0.0]),
                y=np.empty((0, 1)),
                species_labels=[],
                success=True,
                message="No species to simulate.",
            )

        # Build index maps used by the RHS function
        label_to_idx = {lab.lower().strip(): i for i, lab in enumerate(species_labels)}
        alias_to_model_label = self._build_alias_to_model_label(species_labels)
        enzyme_conc_map = self._build_enzyme_concentration_map()

        # Constant mask: core species can be constant, edge species never
        constant_mask = np.array(
            [s.constant for s in self.model.core_species]
            + [False] * len(self.model.edge_species),
            dtype=bool,
        )

        y0 = np.array(self.model.get_all_concentration_vector(), dtype=float)

        # Time span / evaluation points
        t_span = (0.0, end_time)
        t_eval = None
        if time_step is not None and time_step > 0:
            t_eval = np.arange(0.0, end_time + time_step * 0.5, time_step)
            t_eval = t_eval[t_eval <= end_time]

        # Include both core and edge reactions so edge species get flux
        reactions = self.model.core_reactions + self.model.edge_reactions

        def rhs(t, y):
            dydt = np.zeros(n_species)
            conc_dict = self._build_conc_dict_with_ontology(species_labels, y)

            for rxn in reactions:
                v = compute_mm_rate(rxn, conc_dict, enzyme_conc_map)
                if v == 0.0:
                    continue
                for species_label, coeff in rxn.stoichiometry.items():
                    lc = species_label.lower().strip()
                    model_lc = alias_to_model_label.get(lc, lc)
                    idx = label_to_idx.get(model_lc)
                    if idx is not None:
                        dydt[idx] += coeff * v

            # Zero out derivatives for constant species
            dydt[constant_mask] = 0.0
            return dydt

        if self.logger:
            self.logger.info(
                f"Running ODE simulation: {n_species} species, "
                f"{len(reactions)} reactions, t=[0, {end_time}] s"
            )

        sol = solve_ivp(
            rhs,
            t_span,
            y0,
            method=method,
            t_eval=t_eval,
            rtol=rtol,
            atol=atol,
            dense_output=True,
            max_step=end_time / 10 if end_time > 0 else np.inf,
        )

        result = SimulationResult(
            t=sol.t,
            y=sol.y,
            species_labels=species_labels,
            success=sol.success,
            message=sol.message if hasattr(sol, "message") else "",
        )

        # Update model concentrations to final state (both core and edge)
        if sol.success and sol.y.shape[1] > 0:
            final_conc = sol.y[:, -1].tolist()
            self.model.set_all_concentrations(final_conc)

        if self.logger:
            status = "succeeded" if sol.success else "FAILED"
            self.logger.info(f"ODE simulation {status} ({sol.t.shape[0]} time points)")

        return result

    # ------------------------------------------------------------------
    # ODE export
    # ------------------------------------------------------------------

    def export_ode_equations(
        self,
        output_path: str,
        iteration: int,
        end_time: float,
    ) -> str:
        """
        Export the ODE equations solved at this iteration to a text file.

        Writes reaction rate laws (with Km, kcat) and dC/dt for each species.

        Args:
            output_path: Path for the output file.
            iteration: Enlargement iteration number.
            end_time: Simulation end time (seconds).

        Returns:
            The output path written.
        """

        species_labels = self.model.get_all_species_labels()
        reactions = self.model.core_reactions + self.model.edge_reactions
        label_to_idx = {lab.lower().strip(): i for i, lab in enumerate(species_labels)}
        edge_labels_lc = {sp.label.lower().strip() for sp in self.model.edge_species}
        enzyme_map = self._build_enzyme_concentration_map()

        # Build species (lc) -> list of (coeff, rxn_idx) for ODE terms
        species_terms: Dict[str, List[Tuple[int, int]]] = {
            lab.lower().strip(): [] for lab in species_labels
        }
        constant_mask = {}
        for i, sp in enumerate(self.model.core_species):
            constant_mask[sp.label.lower().strip()] = sp.constant
        for sp in self.model.edge_species:
            constant_mask[sp.label.lower().strip()] = False

        lines: List[str] = []
        lines.append("=" * 80)
        lines.append(f"ODE EQUATIONS - Enlargement Iteration {iteration}")
        lines.append("=" * 80)
        lines.append(
            f"Species: {len(species_labels)}  |  Reactions: {len(reactions)}  |  t in [0, {end_time}] s"
        )
        lines.append("")

        # --- Reaction rate laws ---
        lines.append("-" * 80)
        lines.append("Reactions (rate laws with parameters)")
        lines.append("-" * 80)

        for r_idx, rxn in enumerate(reactions, 1):
            reactants = " + ".join(
                s for s, c in rxn.stoichiometry.items() if c < 0
            )
            products = " + ".join(
                s for s, c in rxn.stoichiometry.items() if c > 0
            )
            lines.append(f"R{r_idx}: {reactants} -> {products}")
            lines.append(f"    Enzyme: {rxn.enzyme_label}  ({rxn.ec_number or 'N/A'})")

            kin = rxn.kinetics
            if kin is None or rxn.rate_law is None:
                lines.append("    v = (no kinetics - rate = 0)")
            else:
                e_conc = enzyme_map.get(rxn.enzyme_label.lower().strip(), 0.0)
                kcat = kin.kcat
                vmax = kin.vmax
                km_per = getattr(kin, "km_per_substrate", None) or {}
                km_single = kin.km

                if kcat is not None and e_conc > 0:
                    rate_pre = f"v{r_idx} = kcat * [{rxn.enzyme_label}] * "
                elif vmax is not None:
                    rate_pre = f"v{r_idx} = Vmax * "
                else:
                    rate_pre = f"v{r_idx} = (missing kcat/Vmax)"

                sat_parts = []
                for r in rxn.reactant_labels:
                    km_val = km_per.get(r) if km_per else None
                    if km_val is None:
                        km_val = km_single
                    if km_val is not None and km_val > 0:
                        sat_parts.append(f"[{r}]/(Km_{r}+[{r}])")
                    else:
                        sat_parts.append(f"[{r}]")

                if sat_parts:
                    rate_str = rate_pre + " * ".join(sat_parts)
                else:
                    rate_str = rate_pre.rstrip(" * ")

                lines.append(f"    {rate_str}")

                params = []
                if kcat is not None:
                    params.append(f"kcat={kcat:.6g} 1/s")
                    if getattr(kin, "kcat_sd", None) is not None:
                        params.append(f"kcat_sd={kin.kcat_sd:.6g} 1/s")
                if vmax is not None:
                    params.append(f"Vmax={vmax:.6g} mM/s")
                if km_per:
                    km_sd_per = getattr(kin, "km_sd_per_substrate", None) or {}
                    for r, k in km_per.items():
                        if k is not None:
                            params.append(f"Km({r})={k:.6g} mM")
                            if km_sd_per and r in km_sd_per:
                                params.append(f"Km_sd({r})={km_sd_per[r]:.6g} mM")
                elif km_single is not None:
                    params.append(f"Km={km_single:.6g} mM")
                    if getattr(kin, "km_sd", None) is not None:
                        params.append(f"Km_sd={kin.km_sd:.6g} mM")
                if params:
                    lines.append(f"    Parameters: {', '.join(params)}")
            lines.append("")

            for species_label, coeff in rxn.stoichiometry.items():
                lc = species_label.lower().strip()
                if lc in label_to_idx and lc in species_terms:
                    species_terms[lc].append((coeff, r_idx))

        # --- ODEs for each species ---
        lines.append("-" * 80)
        lines.append("ODEs (dC/dt for each species)")
        lines.append("-" * 80)

        for lab in species_labels:
            lc = lab.lower().strip()
            terms = species_terms.get(lc, [])
            is_edge = lc in edge_labels_lc
            is_const = constant_mask.get(lc, False)

            if is_const:
                lines.append(f"d[{lab}]/dt = 0  (constant)")
            elif not terms:
                lines.append(f"d[{lab}]/dt = 0")
            else:
                term_strs = []
                for coeff, r_idx in terms:
                    if coeff > 0:
                        term_strs.append(f"+{coeff}*v{r_idx}")
                    else:
                        term_strs.append(f"{coeff}*v{r_idx}")
                rhs = " ".join(term_strs).lstrip("+") or "0"
                suffix = "  (edge)" if is_edge else ""
                lines.append(f"d[{lab}]/dt = {rhs}{suffix}")
        lines.append("")
        lines.append("=" * 80)

        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        with open(output_path, "w") as f:
            f.write("\n".join(lines))

        if self.logger:
            self.logger.info(f"Exported ODE equations to {output_path}")

        return output_path

    # ------------------------------------------------------------------
    # Edge species rate evaluation
    # ------------------------------------------------------------------

    def evaluate_edge_rates(
        self,
        sim_result: SimulationResult,
    ) -> Dict[str, float]:
        """
        Evaluate net production rates for *edge* species using peak flux over
        the whole simulation (RMG-style).

        Uses the maximum |flux| over all time points so that species with
        significant transient flux are promoted even when near steady state
        at the end. This matches RMG: "flux at some point" for promotion.

        Returns:
            Dictionary mapping edge species label (lc) -> signed rate (mM/s)
            at the time of peak |flux|.
        """
        if not sim_result.success or sim_result.y.shape[1] == 0:
            return {}

        species_labels = sim_result.species_labels
        alias_to_model_label = self._build_alias_to_model_label(species_labels)
        enzyme_conc_map = self._build_enzyme_concentration_map()
        edge_labels_lc = {sp.label.lower().strip() for sp in self.model.edge_species}
        all_reactions = self.model.core_reactions + self.model.edge_reactions

        # Track peak |flux| and signed flux at that time for each edge species
        peak_abs: Dict[str, float] = {}
        flux_at_peak: Dict[str, float] = {}

        for t in range(sim_result.y.shape[1]):
            y_t = sim_result.y[:, t]
            conc_dict = self._build_conc_dict_with_ontology(species_labels, y_t)

            instant_rates: Dict[str, float] = {}
            for rxn in all_reactions:
                v = compute_mm_rate(rxn, conc_dict, enzyme_conc_map)
                if v == 0.0:
                    continue
                for species_label, coeff in rxn.stoichiometry.items():
                    lc = species_label.lower().strip()
                    model_lc = alias_to_model_label.get(lc, lc)
                    if model_lc in edge_labels_lc:
                        instant_rates[model_lc] = instant_rates.get(model_lc, 0.0) + coeff * v

            for lc, rate in instant_rates.items():
                abs_r = abs(rate)
                if abs_r > peak_abs.get(lc, 0.0):
                    peak_abs[lc] = abs_r
                    flux_at_peak[lc] = rate

        return flux_at_peak

    def evaluate_core_rates(
        self,
        sim_result: SimulationResult,
    ) -> Dict[str, float]:
        """
        Evaluate net production rates for *core* species at the final
        simulation time point. Uses all reactions (core + edge) since
        core species can participate in edge reactions.

        Returns:
            Dictionary mapping core species label (lc) -> net rate (mM/s).
        """
        if not sim_result.success or sim_result.y.shape[1] == 0:
            return {}

        species_labels = sim_result.species_labels
        final_y = sim_result.y[:, -1]
        conc_dict = self._build_conc_dict_with_ontology(species_labels, final_y)
        alias_to_model_label = self._build_alias_to_model_label(species_labels)
        enzyme_conc_map = self._build_enzyme_concentration_map()

        core_labels_lc = {sp.label.lower().strip() for sp in self.model.core_species}
        core_rates: Dict[str, float] = {}

        # Use all reactions: core species can be consumed/produced by edge reactions
        all_reactions = self.model.core_reactions + self.model.edge_reactions
        for rxn in all_reactions:
            v = compute_mm_rate(rxn, conc_dict, enzyme_conc_map)
            if v == 0.0:
                continue
            for species_label, coeff in rxn.stoichiometry.items():
                lc = species_label.lower().strip()
                model_lc = alias_to_model_label.get(lc, lc)
                if model_lc in core_labels_lc:
                    core_rates[model_lc] = core_rates.get(model_lc, 0.0) + coeff * v

        return core_rates

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _build_conc_dict_with_ontology(
        self,
        species_labels: List[str],
        y_vector: np.ndarray,
    ) -> Dict[str, float]:
        """
        Build concentration dict (label_lc -> mM) including ontology equivalents
        so that reactions using DB names (e.g. D-fructose 1,6-bisphosphate)
        resolve to the model species concentration (e.g. beta-d-fructofuranose
        1,6-bisphosphate).
        """
        conc_dict: Dict[str, float] = {}
        for i, lab in enumerate(species_labels):
            val = max(y_vector[i], 0.0) if i < len(y_vector) else 0.0
            lc = lab.lower().strip()
            conc_dict[lc] = val
            for equiv in get_ontology_equivalents(lab):
                conc_dict[equiv.lower().strip()] = val
        return conc_dict

    def _build_alias_to_model_label(
        self,
        species_labels: List[str],
    ) -> Dict[str, str]:
        """
        Map any ontology alias (lowercase) to the model's species label (lc).
        Used to apply flux to the correct state variable when reaction
        stoichiometry uses a different name for the same compound.
        """
        alias_to_model: Dict[str, str] = {}
        for lab in species_labels:
            model_lc = lab.lower().strip()
            alias_to_model[model_lc] = model_lc
            for equiv in get_ontology_equivalents(lab):
                alias_to_model[equiv.lower().strip()] = model_lc
        return alias_to_model

    def _build_enzyme_concentration_map(self) -> Dict[str, float]:
        """
        Build a mapping of enzyme label (lowercase) -> concentration (mM)
        from the core species that are flagged as enzymes.
        """
        enzyme_map: Dict[str, float] = {}
        for sp in self.model.core_species:
            if sp.is_enzyme:
                enzyme_map[sp.label.lower().strip()] = sp.concentration
        return enzyme_map
