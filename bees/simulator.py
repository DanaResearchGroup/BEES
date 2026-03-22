#!/usr/bin/env python3

"""
ODE Simulator Module
--------------------
Integrates the reaction network over time using scipy's ODE solvers.

All concentrations are in mM, time in seconds, rates in mM/s.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import math
import os
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import csc_matrix

from bees.common import get_ontology_equivalents
from bees.core_edge_model import CoreEdgeModel
from bees.flux_calculator import compute_mm_rate


@dataclass
class SimulationResult:
    """
    ODE simulation output.

   
    - **y shape**: (n_species, n_timepoints). Concentrations are clipped to >= 0 on success.

    Key fields:
    - `max_char_rate` / `final_char_rate`: peak and final core characteristic rate (R_char).
    - `max_edge_rate_ratio`: per edge species, max over time of |edge_rate| / R_char (when R_char > 0).
    - `peak_edge_signed_rate`: signed edge rate at the time of peak |flux|.
    - `simulation_interrupted` + `interrupt_*`: snapshot of rates at the interrupt time (if used).
    """
    t: np.ndarray
    y: np.ndarray
    species_labels: List[str]
    success: bool = True
    message: str = ""
    max_char_rate: float = 0.0
    final_char_rate: float = 0.0
    max_edge_rate_ratio: Dict[str, float] = field(default_factory=dict)
    peak_edge_signed_rate: Dict[str, float] = field(default_factory=dict)
    simulation_interrupted: bool = False
    interrupt_char_rate: float = 0.0
    interrupt_edge_rates: Dict[str, float] = field(default_factory=dict)
    # RMG-style reaction-level dynamics metric (dlnaccum) per edge reaction.
    # Keyed by reaction signature ((reactants_sorted_tuple, products_sorted_tuple)).
    # max_edge_reaction_dlnaccum is the maximum dlnaccum_j seen during the run.
    max_edge_reaction_dlnaccum: Dict[Tuple, float] = field(default_factory=dict)
    # Snapshot of dlnaccum at the moment of interrupt (only set if interrupt fires
    # via the dlnaccum criterion).
    interrupt_edge_reaction_dlnaccum: Dict[Tuple, float] = field(default_factory=dict)


class _VectorizedRHS:
    """
    Pre-compiled, NumPy-vectorized ODE right-hand side for a fixed reaction
    network.

    Builds all matrices and index arrays once at construction time so that
    each call to ``__call__(t, y)`` runs without any Python-level loops or
    dict lookups.

    Parameters
    ----------
    species_labels : list of str
        Ordered list of all species labels (core + edge).
    reactions : list of GeneratedReaction
        All reactions (core + edge) to include.
    alias_to_model_label : dict
        Mapping from any alias (lc) to the canonical model label (lc).
    enzyme_conc_map : dict
        Enzyme label (lc) -> concentration (mM).
    constant_mask : np.ndarray of bool
        True for species whose concentration must be held constant (dydt=0).
    """

    def __init__(
        self,
        species_labels: List[str],
        reactions: list,
        alias_to_model_label: Dict[str, str],
        enzyme_conc_map: Dict[str, float],
        constant_mask: np.ndarray,
    ):
        n_sp = len(species_labels)
        n_rx = len(reactions)
        label_to_idx = {lab.lower().strip(): i for i, lab in enumerate(species_labels)}

        # ---- Per-reaction Vmax ----------------------------------------
        vmax = np.zeros(n_rx, dtype=np.float64)
        for j, rxn in enumerate(reactions):
            kin = rxn.kinetics
            if kin is None or rxn.rate_law is None:
                continue
            e_key = rxn.enzyme_label.lower().strip()
            e_conc = enzyme_conc_map.get(e_key, 0.0)
            if kin.kcat is not None and e_conc > 0:
                vmax[j] = kin.kcat * e_conc
            elif kin.vmax is not None:
                vmax[j] = kin.vmax

        # ---- Saturation substrate lists (ragged -> matrix) -------------
        # For each reaction we store a list of (species_idx, km_value) pairs.
        # We pack these into matrices (n_rx, K) where K is max substrates.
        sub_idx_list: List[List[int]] = []
        km_list:      List[List[float]] = []

        for rxn in reactions:
            kin = rxn.kinetics
            pairs_idx: List[int] = []
            pairs_km:  List[float] = []
            if kin is not None and rxn.rate_law is not None:
                km_per = getattr(kin, "km_per_substrate", None) or {}
                km_single = kin.km
                for reactant in rxn.reactant_labels:
                    r_lc = reactant.lower().strip()
                    km_val = None
                    if km_per:
                        km_val = km_per.get(reactant)
                        if km_val is None:
                            km_val = next(
                                (v for k, v in km_per.items()
                                 if k.lower().strip() == r_lc),
                                None,
                            )
                        if km_val is None:
                            continue
                    else:
                        km_val = km_single
                    if km_val is None or km_val <= 0:
                        continue
                    model_lc = alias_to_model_label.get(r_lc, r_lc)
                    idx = label_to_idx.get(model_lc)
                    if idx is None:
                        idx = label_to_idx.get(r_lc)
                    if idx is not None:
                        pairs_idx.append(idx)
                        pairs_km.append(float(km_val))
            sub_idx_list.append(pairs_idx)
            km_list.append(pairs_km)

        # Build (n_rx, K) matrices
        K = max((len(idxs) for idxs in sub_idx_list), default=0)
        # Pad with dummy index pointing to a constant 1.0 concentration
        self._sub_idx_mat = np.full((n_rx, K), fill_value=n_sp, dtype=np.int32)
        self._km_mat = np.zeros((n_rx, K), dtype=np.float64)
        for j, (idxs, kms) in enumerate(zip(sub_idx_list, km_list)):
            nj = len(idxs)
            if nj > 0:
                self._sub_idx_mat[j, :nj] = idxs
                self._km_mat[j, :nj] = kms

        # ---- Stoichiometry matrix S  (n_species x n_reactions) ---------
        rows, cols, data = [], [], []
        for j, rxn in enumerate(reactions):
            for sp_label, coeff in rxn.stoichiometry.items():
                lc = sp_label.lower().strip()
                model_lc = alias_to_model_label.get(lc, lc)
                idx = label_to_idx.get(model_lc)
                if idx is not None:
                    rows.append(idx)
                    cols.append(j)
                    data.append(float(coeff))
        S = csc_matrix((data, (rows, cols)), shape=(n_sp, n_rx), dtype=np.float64)

        # Store as dense if small, sparse otherwise
        if n_rx <= 200:
            self._S = S.toarray()
            self._sparse = False
        else:
            self._S = S
            self._sparse = True

        # Pre-compute production and consumption matrices for metrics
        S_prod = S.copy()
        S_prod.data[S_prod.data < 0] = 0.0
        S_prod.eliminate_zeros()
        S_cons = S.copy()
        S_cons.data[S_cons.data > 0] = 0.0
        S_cons.data = np.abs(S_cons.data)
        S_cons.eliminate_zeros()

        if n_rx <= 200:
            self._S_prod = S_prod.toarray()
            self._S_cons = S_cons.toarray()
        else:
            self._S_prod = S_prod
            self._S_cons = S_cons

        self._vmax = vmax
        self._n_rx = n_rx
        self._constant_mask = constant_mask
        self._n_sp = n_sp
        self.label_to_idx = label_to_idx
        self.alias_to_model_label = alias_to_model_label

    def __call__(self, t: float, y: np.ndarray) -> np.ndarray:
        """Evaluate dydt = S @ v(y)."""
        return self.compute_dydt(y)

    def compute_v(self, y: np.ndarray) -> np.ndarray:
        """
        Compute reaction rate vector v(y).
        Supports y as (n_species,) or (n_species, n_timepoints).
        """
        is_mat = (y.ndim > 1)
        y_nn = np.maximum(y, 0.0)

        # Append dummy 1.0 for padding
        if is_mat:
            ones = np.ones((1, y.shape[1]), dtype=y.dtype)
            y_ext = np.vstack([y_nn, ones])
        else:
            y_ext = np.append(y_nn, 1.0)

        # Saturation: s / (Km + s)
        # s shape: (n_rx, K) if vector, (n_rx, K, n_t) if matrix
        s = y_ext[self._sub_idx_mat]
        
        if is_mat:
            km = self._km_mat[:, :, np.newaxis]
            sat = s / (km + s)
            # product over K dimension (axis 1)
            v = self._vmax[:, np.newaxis] * np.prod(sat, axis=1)
        else:
            sat = s / (self._km_mat + s)
            v = self._vmax * np.prod(sat, axis=1)
        
        return v

    def compute_dydt(self, y: np.ndarray) -> np.ndarray:
        """Compute dydt = S @ v(y). Supports matrix y."""
        v = self.compute_v(y)
        if self._sparse:
            dydt = self._S.dot(v)
            if hasattr(dydt, "toarray"):
                dydt = dydt.toarray()
            if dydt.ndim > 1 and v.ndim == 1:
                dydt = dydt.ravel()
        else:
            dydt = self._S @ v
        
        if dydt.ndim > 1:
            dydt[self._constant_mask, :] = 0.0
        else:
            dydt[self._constant_mask] = 0.0
        return dydt

    def compute_prod_cons(self, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Compute production and consumption rates for all species."""
        v = self.compute_v(y)
        if self._sparse:
            p = self._S_prod.dot(v)
            c = self._S_cons.dot(v)
            if hasattr(p, "toarray"):
                p = p.toarray()
            if hasattr(c, "toarray"):
                c = c.toarray()
            if p.ndim > 1 and v.ndim == 1:
                p = p.ravel()
            if c.ndim > 1 and v.ndim == 1:
                c = c.ravel()
        else:
            p = self._S_prod @ v
            c = self._S_cons @ v
        return p, c


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

    @staticmethod
    def _concentrations_nonnegative(y: np.ndarray) -> np.ndarray:
        """Project concentrations (mM) to non-negative values, in-place safe."""
        out = np.asarray(y, dtype=float)
        return np.maximum(out, 0.0)

    # ------------------------------------------------------------------
     


    def simulate(
        self,
        end_time: float,
        time_step: Optional[float] = None,
        method: str = "BDF",
        rtol: float = 1e-8,
        atol: float = 1e-10,
        interrupt_simulation_tol: Optional[float] = None,
        tol_move_edge_reaction_to_core: Optional[float] = None,
        max_wall_time_s: Optional[float] = None,
        stepwise_heartbeat_interval_s: Optional[float] = None,
    ) -> SimulationResult:
        """
        Run an ODE simulation of the full model (core + edge species and reactions).

        Args:
            end_time: Simulation end time in seconds.
            time_step: If given, store solution at these intervals.
            method: Integration method for solve_ivp (default BDF,
                    good for stiff biochemical systems).
            rtol: Relative tolerance for the solver.
            atol: Absolute tolerance for the solver.
            interrupt_simulation_tol: If set, stop integration when any edge species
                flux ratio |R_edge|/R_char exceeds this value.
                Iterative enlarger passes ``Settings.toleranceInterruptSimulation``
                (defaulting to ``toleranceMoveToCore`` when unset).
            max_wall_time_s: Stepwise mode only. Stop the pass after this many seconds
                of wall-clock time (``success=False``).
            stepwise_heartbeat_interval_s: Stepwise mode only. Log INFO progress
                periodically. None defaults to 30 s. 0 disables.

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

        edge_labels_lc = {sp.label.lower().strip() for sp in self.model.edge_species}

        use_stepwise = (
            interrupt_simulation_tol is not None
            and interrupt_simulation_tol > 0
            and len(edge_labels_lc) > 0
        )
        if use_stepwise:
            return self._simulate_stepwise(
                end_time=end_time,
                time_step=time_step,
                method=method,
                rtol=rtol,
                atol=atol,
                interrupt_simulation_tol=float(interrupt_simulation_tol),
                tol_move_edge_reaction_to_core=tol_move_edge_reaction_to_core,
                max_wall_time_s=max_wall_time_s,
                stepwise_heartbeat_interval_s=stepwise_heartbeat_interval_s,
            )
        return self._simulate_continuous(
            end_time=end_time,
            time_step=time_step,
            method=method,
            rtol=rtol,
            atol=atol,
        )

    # ------------------------------------------------------------------
    # Continuous (non-interrupt) simulation path
    # ------------------------------------------------------------------

    def _simulate_continuous(
        self,
        end_time: float,
        time_step: Optional[float],
        method: str,
        rtol: float,
        atol: float,
    ) -> SimulationResult:
        """Full-span solve_ivp without flux-interrupt events."""
        species_labels = self.model.get_all_species_labels()
        n_species = len(species_labels)
        label_to_idx = {lab.lower().strip(): i for i, lab in enumerate(species_labels)}
        alias_to_model_label = self._build_alias_to_model_label(species_labels)
        enzyme_conc_map = self._build_enzyme_concentration_map()
        constant_mask = np.array(
            [s.constant for s in self.model.core_species]
            + [False] * len(self.model.edge_species),
            dtype=bool,
        )
        y0 = self._concentrations_nonnegative(
            np.array(self.model.get_all_concentration_vector(), dtype=float)
        )
        reactions = self.model.core_reactions + self.model.edge_reactions

        t_eval = None
        if time_step is not None and time_step > 0:
            t_eval = np.arange(0.0, end_time + time_step * 0.5, time_step)
            t_eval = t_eval[t_eval <= end_time]

        rhs = _VectorizedRHS(
            species_labels=species_labels,
            reactions=reactions,
            alias_to_model_label=alias_to_model_label,
            enzyme_conc_map=enzyme_conc_map,
            constant_mask=constant_mask,
        )

        if self.logger:
            self.logger.debug(
                f"Running ODE simulation (continuous): {n_species} species, "
                f"{len(reactions)} reactions, t=[0, {end_time}] s"
            )

        sol = solve_ivp(
            rhs, (0.0, end_time), y0,
            method=method, t_eval=t_eval,
            rtol=rtol, atol=atol,
            dense_output=True,
            max_step=end_time / 10 if end_time > 0 else np.inf,
        )

        t_out = sol.t
        y_raw = sol.y
        extra_msg = ""
        ok = sol.success
        if y_raw.size > 0 and ok:
            if not np.isfinite(y_raw).all():
                ok = False
                extra_msg = "Non-finite values (NaN or inf) in concentration state."
                if self.logger:
                    self.logger.error(f"ODE state invalid: {extra_msg}")
        y_out = (
            self._concentrations_nonnegative(y_raw) if y_raw.size else y_raw
        )

        result = SimulationResult(
            t=t_out, y=y_out,
            species_labels=species_labels,
            success=ok,
            message=(
                (sol.message if hasattr(sol, "message") else "")
                + (f" {extra_msg}" if extra_msg else "")
            ).strip(),
        )

        if ok and y_out.shape[1] > 0:
            self.model.set_all_concentrations(y_out[:, -1].tolist())

        self._attach_flux_metrics(result, rhs=rhs)

        if self.logger:
            if not sol.success:
                self.logger.warning(
                    f"ODE solver FAILED: {getattr(sol, 'message', 'unknown')}"
                )
            elif ok:
                self.logger.debug(
                    f"ODE simulation completed "
                    f"({t_out.shape[0]} time points)"
                )
        return result

    # ------------------------------------------------------------------
    # Simulation with flux-ratio interrupt
    # ------------------------------------------------------------------

    def _simulate_stepwise(
        self,
        end_time: float,
        time_step: Optional[float],
        method: str,
        rtol: float,
        atol: float,
        interrupt_simulation_tol: float,
        tol_move_edge_reaction_to_core: Optional[float] = None,
        max_wall_time_s: Optional[float] = None,
        stepwise_heartbeat_interval_s: Optional[float] = None,
    ) -> SimulationResult:
        """
        advance by dt, compute R_char and edge rr_i
        after each step, interrupt via boolean when any rr_i exceeds the
        tolerance. At least one step is always taken before an interrupt
        can fire (earliest interrupt at t > 0).
        """
        species_labels = self.model.get_all_species_labels()
        n_species = len(species_labels)
        label_to_idx = {lab.lower().strip(): i for i, lab in enumerate(species_labels)}
        alias_to_model_label = self._build_alias_to_model_label(species_labels)
        enzyme_conc_map = self._build_enzyme_concentration_map()
        constant_mask = np.array(
            [s.constant for s in self.model.core_species]
            + [False] * len(self.model.edge_species),
            dtype=bool,
        )
        y0 = self._concentrations_nonnegative(
            np.array(self.model.get_all_concentration_vector(), dtype=float)
        )
        reactions = self.model.core_reactions + self.model.edge_reactions
        edge_labels_lc = {sp.label.lower().strip() for sp in self.model.edge_species}
        edge_labels_lc_list = [sp.label.lower().strip() for sp in self.model.edge_species]
        n_core = len(self.model.core_species)
        n_core_rxns = len(self.model.core_reactions)

        rhs = _VectorizedRHS(
            species_labels=species_labels,
            reactions=reactions,
            alias_to_model_label=alias_to_model_label,
            enzyme_conc_map=enzyme_conc_map,
            constant_mask=constant_mask,
        )

        # ---- Edge-reaction setup for dlnaccum (RMG reaction-level criterion).
        # For each edge reaction j, we need the species indices (reactants and
        # products separately) to look up consumption/production at each step.
        # Signature must match `bees.exporter.reaction_signature`:
        # (enzyme_lc, reactants_sorted_tuple, products_sorted_tuple).
        edge_rxn_sigs: List[Tuple[str, Tuple[str, ...], Tuple[str, ...]]] = []
        edge_rxn_reactant_idx: List[List[int]] = []
        edge_rxn_product_idx: List[List[int]] = []
        edge_rxn_global_idx: List[int] = []
        for j_local, rxn in enumerate(self.model.edge_reactions):
            j_global = n_core_rxns + j_local
            r_idx: List[int] = []
            p_idx: List[int] = []
            for sp_label, coeff in rxn.stoichiometry.items():
                lc = sp_label.lower().strip()
                model_lc = alias_to_model_label.get(lc, lc)
                idx = label_to_idx.get(model_lc)
                if idx is None:
                    idx = label_to_idx.get(lc)
                if idx is None:
                    continue
                # Stoichiometric multiplicity: emulate RMG, which loops over
                # reactant/product *index lists* (a coefficient of 2 contributes twice).
                mult = int(round(abs(float(coeff))))
                if mult <= 0:
                    mult = 1
                if float(coeff) < 0:
                    r_idx.extend([idx] * mult)
                elif float(coeff) > 0:
                    p_idx.extend([idx] * mult)
            sig = (
                str(getattr(rxn, "enzyme_label", "")).lower().strip(),
                tuple(sorted(str(r).lower().strip() for r in rxn.reactant_labels)),
                tuple(sorted(str(p).lower().strip() for p in rxn.product_labels)),
            )
            edge_rxn_sigs.append(sig)
            edge_rxn_reactant_idx.append(r_idx)
            edge_rxn_product_idx.append(p_idx)
            edge_rxn_global_idx.append(j_global)

        max_edge_rxn_dlnaccum: Dict[Tuple, float] = {sig: 0.0 for sig in edge_rxn_sigs}
        interrupt_edge_rxn_dlnaccum: Dict[Tuple, float] = {}
        # Small floor to keep ln(1 + v/R) finite when R ~ 0.
        _RATE_FLOOR = 1e-30

        if stepwise_heartbeat_interval_s is None:
            heartbeat_interval = 30.0
        else:
            heartbeat_interval = float(stepwise_heartbeat_interval_s)

        if self.logger:
            self.logger.debug(
                f"Running ODE simulation (stepwise): {n_species} species, "
                f"{len(reactions)} reactions, t=[0, {end_time}] s, "
                f"tol_interrupt={interrupt_simulation_tol}"
            )
            if max_wall_time_s is not None:
                self.logger.debug(
                    f"Stepwise wall-clock limit: {max_wall_time_s} s per pass"
                )
            if heartbeat_interval > 0:
                self.logger.debug(
                    f"Stepwise heartbeat interval: {heartbeat_interval} s"
                )

        t_history: List[float] = [0.0]
        y_history: List[np.ndarray] = [y0.copy()]

        wall_start = time.monotonic()
        last_heartbeat = wall_start
        wall_time_exceeded = False

        current_t = 0.0
        current_y = y0.copy()
        # Outer-step size bounds. A floor avoids spending wall time on 1e-12 s
        # segments when max_rr sits below interrupt_tol (no interrupt possible
        # yet) but old backoff logic kept shrinking dt.
        dt_max = 5.0
        if end_time > 0.0:
            dt_floor_outer = max(1e-12, min(1e-6, end_time * 1e-8))
        else:
            dt_floor_outer = 1e-12
        dt = dt_floor_outer
        dt_min = dt_floor_outer
        dt_growth = 2.0
        near_frac = 0.5
        max_char_rate = 0.0
        interrupted = False
        interrupt_char = 0.0
        interrupt_edge: Dict[str, float] = {}
        first_step = True
        outer_step_idx = 0
        state_invalid = False
        invalid_message = ""

        while current_t < end_time:
            now_wall = time.monotonic()
            if max_wall_time_s is not None and (now_wall - wall_start) > max_wall_time_s:
                wall_time_exceeded = True
                if self.logger:
                    self.logger.warning(
                        "ODE stepwise pass stopped: exceeded "
                        f"max_wall_time_per_iteration ({max_wall_time_s} s wall-clock) "
                        f"at simulation t={current_t:.6e} s."
                    )
                break

            dt = min(max(dt, dt_min), dt_max)
            target_t = min(current_t + dt, end_time)
            if target_t <= current_t:
                break

            try:
                sol = solve_ivp(
                    rhs, (current_t, target_t), current_y,
                    method=method, rtol=rtol, atol=atol,
                    dense_output=False,
                    max_step=target_t - current_t,
                )
            except Exception:
                if self.logger:
                    self.logger.warning(
                        f"ODE solver exception at t={current_t:.6e} s"
                    )
                break

            if not sol.success or sol.y.shape[1] == 0:
                if self.logger:
                    self.logger.warning(
                        f"ODE solver failed at t={current_t:.6e} s: "
                        f"{getattr(sol, 'message', 'unknown')}"
                    )
                break

            current_t = float(sol.t[-1])
            current_y = sol.y[:, -1].copy()
            outer_step_idx += 1
            if not np.isfinite(current_y).all():
                state_invalid = True
                invalid_message = (
                    "Non-finite values (NaN or inf) in concentration state."
                )
                if self.logger:
                    self.logger.error(
                        f"ODE state invalid at t={current_t:.6e} s after outer step. "
                        f"{invalid_message}"
                    )
                break
            current_y = self._concentrations_nonnegative(current_y)

            t_history.append(current_t)
            y_history.append(current_y.copy())

            dydt_all = rhs.compute_dydt(current_y)
            char_rate = float(np.linalg.norm(dydt_all[:n_core]))
            edge_rates = {
                lab: float(dydt_all[n_core + i])
                for i, lab in enumerate(edge_labels_lc_list)
            }
            max_char_rate = max(max_char_rate, char_rate)

            max_rr = 0.0
            if char_rate > 0.0 and edge_rates:
                max_rr = max((abs(r) / char_rate for r in edge_rates.values()), default=0.0)

            # ---- RMG-style reaction-level dlnaccum (instantaneous).
            # dlnaccum_j = Σ_i ln(1 + v_j / R_i) over species i touched by edge
            # reaction j; R_i is consumption (L_i) for reactants, production
            # (P_i) for products. Tracks per-step max over the run.
            max_rxn_dlnaccum_step = 0.0
            max_rxn_dlnaccum_sig: Optional[Tuple] = None
            if edge_rxn_sigs:
                v_all = rhs.compute_v(current_y)
                prod_all, cons_all = rhs.compute_prod_cons(current_y)
                # prod_all/cons_all are species production/consumption magnitudes.
                for k, sig in enumerate(edge_rxn_sigs):
                    j_global = edge_rxn_global_idx[k]
                    v_j = float(v_all[j_global])
                    if v_j <= 0.0:
                        continue
                    dln = 0.0
                    for i_idx in edge_rxn_reactant_idx[k]:
                        L_i = float(cons_all[i_idx])
                        dln += math.log1p(v_j / max(L_i, _RATE_FLOOR))
                    for i_idx in edge_rxn_product_idx[k]:
                        P_i = float(prod_all[i_idx])
                        dln += math.log1p(v_j / max(P_i, _RATE_FLOOR))
                    if dln > max_edge_rxn_dlnaccum[sig]:
                        max_edge_rxn_dlnaccum[sig] = dln
                    if dln > max_rxn_dlnaccum_step:
                        max_rxn_dlnaccum_step = dln
                        max_rxn_dlnaccum_sig = sig

            now_wall = time.monotonic()
            if heartbeat_interval > 0 and self.logger and (
                now_wall - last_heartbeat >= heartbeat_interval
            ):
                elapsed = now_wall - wall_start
                self.logger.info(
                    "ODE stepwise progress: "
                    f"sim_t={current_t:.6e} s, wall_elapsed={elapsed:.1f} s, "
                    f"dt={dt:.4e} s, "
                    f"max_rr={max_rr:.4e}, R_char={char_rate:.4e} mM/s"
                )
                last_heartbeat = now_wall

            if max_wall_time_s is not None and (now_wall - wall_start) > max_wall_time_s:
                wall_time_exceeded = True
                if self.logger:
                    self.logger.warning(
                        "ODE stepwise pass stopped: exceeded "
                        f"max_wall_time_per_iteration ({max_wall_time_s} s wall-clock) "
                        f"after completing outer step at simulation t={current_t:.6e} s "
                                            )
                break

            if not first_step and max_rr > interrupt_simulation_tol:
                interrupted = True
                interrupt_char = char_rate
                interrupt_edge = edge_rates
                # Snapshot current dlnaccum at the interrupt — these are
                # available as candidates for reaction-level promotion too.
                interrupt_edge_rxn_dlnaccum = dict(max_edge_rxn_dlnaccum)
                if self.logger:
                    self.logger.info(
                        f"ODE interrupted at t={current_t:.6e} s "
                        f"(max_rr={max_rr:.4e}, "
                        f"tol={interrupt_simulation_tol}, "
                        f"{len(t_history)} time points)"
                    )
                break

            # Reaction-level interrupt: if any edge reaction's dlnaccum exceeds
            # the threshold, halt and record candidates for promotion.
            if (
                not first_step
                and tol_move_edge_reaction_to_core is not None
                and tol_move_edge_reaction_to_core > 0
                and max_rxn_dlnaccum_step > tol_move_edge_reaction_to_core
            ):
                interrupted = True
                interrupt_char = char_rate
                interrupt_edge = edge_rates
                interrupt_edge_rxn_dlnaccum = dict(max_edge_rxn_dlnaccum)
                if self.logger:
                    self.logger.info(
                        f"ODE interrupted (reaction dlnaccum) at t={current_t:.6e} s "
                        f"(max_dlnaccum={max_rxn_dlnaccum_step:.4e}, "
                        f"tol={tol_move_edge_reaction_to_core}, "
                        f"sig={max_rxn_dlnaccum_sig})"
                    )
                break

            first_step = False

            # Outer-step adaptation (no shrink while max_rr <= tol):
            # Interrupt fires only when max_rr > interrupt_simulation_tol. If
            # max_rr is below tol, shrinking dt cannot create an interrupt; it
            # only traps the integrator at dt_min. Hold dt in the "near"
            # band; grow when comfortably below.
            dt = max(dt, dt_floor_outer)
            if max_rr >= near_frac * interrupt_simulation_tol:
                pass
            else:
                dt = min(dt * dt_growth, dt_max)

        t_arr = np.array(t_history)
        y_arr = np.column_stack(y_history)

        fail_msg = ""
        if wall_time_exceeded:
            fail_msg = (
                f"Exceeded max_wall_time_per_iteration ({max_wall_time_s} s wall-clock)."
            )
        if state_invalid:
            fail_msg = (
                f"{fail_msg + ' ' if fail_msg else ''}"
                f"Invalid concentration state: {invalid_message}"
            ).strip()

        result = SimulationResult(
            t=t_arr,
            y=y_arr,
            species_labels=species_labels,
            success=(not wall_time_exceeded and not state_invalid),
            message=fail_msg,
            simulation_interrupted=interrupted,
            max_char_rate=max_char_rate,
        )

        if result.success and y_arr.shape[1] > 0:
            self.model.set_all_concentrations(y_arr[:, -1].tolist())

        self._attach_flux_metrics(result, rhs=rhs)
        # Stepwise path already computes per-step maxima; continuous path computes here.

        if interrupted:
            result.interrupt_char_rate = interrupt_char
            result.interrupt_edge_rates = interrupt_edge
            result.interrupt_edge_reaction_dlnaccum = interrupt_edge_rxn_dlnaccum
        else:
            if self.logger and not wall_time_exceeded and not state_invalid:
                self.logger.debug(
                    f"ODE simulation completed (stepwise) "
                    f"({len(t_history)} time points)"
                )

        # Always expose the lifetime-max dlnaccum so the enlarger can promote
        # reactions discovered above tolerance even on a non-interrupted pass.
        result.max_edge_reaction_dlnaccum = max_edge_rxn_dlnaccum

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
                    # IMPORTANT: when per-substrate Km exists, only apply Km
                    # saturation terms to those explicitly provided. Do not
                    # apply a global/single Km fallback to other reactants.
                    km_val = None
                    if km_per:
                        km_val = km_per.get(r)
                        if km_val is None:
                            # No per-substrate Km for this reactant -> print as a
                            # plain concentration factor.
                            km_val = None
                    else:
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
    # RMG-style flux trajectory metrics
    # ------------------------------------------------------------------

    def _attach_flux_metrics(self, sim_result: SimulationResult, rhs: _VectorizedRHS) -> None:
        """
        Populate max_char_rate, final_char_rate, max_edge_rate_ratio, and
        peak_edge_signed_rate by scanning the stored trajectory using vectorized RHS.
        """
        if not sim_result.success or sim_result.y.size == 0:
            sim_result.max_char_rate = 0.0
            sim_result.final_char_rate = 0.0
            sim_result.max_edge_rate_ratio = {}
            sim_result.peak_edge_signed_rate = {}
            return

        n_core = len(self.model.core_species)
        n_edge = len(self.model.edge_species)
        edge_labels = [sp.label.lower().strip() for sp in self.model.edge_species]
        
        # Vectorized calculation over all time points
        dydt_all_mat = rhs.compute_dydt(sim_result.y)
        # dydt_all_mat is (n_sp, n_t)
        
        dydt_core_mat = dydt_all_mat[:n_core, :]
        char_rates = np.linalg.norm(dydt_core_mat, axis=0)
        
        sim_result.max_char_rate = float(np.max(char_rates))
        sim_result.final_char_rate = float(char_rates[-1])
        
        if n_edge > 0:
            dydt_edge_mat = dydt_all_mat[n_core:, :]
            abs_edge_mat = np.abs(dydt_edge_mat)
            
            # Peak flux metrics
            peak_indices = np.argmax(abs_edge_mat, axis=1)
            flux_at_peak = dydt_edge_mat[np.arange(n_edge), peak_indices]
            sim_result.peak_edge_signed_rate = {
                edge_labels[i]: float(flux_at_peak[i]) for i in range(n_edge)
            }
            
            # Max ratio metrics
            # rr_mat = abs_edge_mat / char_rates
            # avoid divide by zero
            safe_char = char_rates.copy()
            safe_char[safe_char == 0] = 1e-100 # essentially zero but avoids NaN
            rr_mat = abs_edge_mat / safe_char
            max_ratios = np.max(rr_mat, axis=1)
            sim_result.max_edge_rate_ratio = {
                edge_labels[i]: float(max_ratios[i]) for i in range(n_edge)
            }
        else:
            sim_result.peak_edge_signed_rate = {}
            sim_result.max_edge_rate_ratio = {}

    def _compute_instantaneous_rates(
        self,
        y_vector: np.ndarray,
        species_labels: List[str],
        reactions: list,
        alias_to_model_label: Dict[str, str],
        enzyme_conc_map: Dict[str, float],
        core_labels_lc: set,
        edge_labels_lc: set,
    ) -> Tuple[float, Dict[str, float]]:
        """
        Compute instantaneous R_char and per-edge-species net rates at an
        arbitrary state vector.

        Returns:
            (char_rate, edge_rates) where char_rate is the L2 norm of core
            species net rates (mM/s) and edge_rates maps edge label (lc) to
            its instantaneous net rate (mM/s).
        """
        conc_dict = self._build_conc_dict_with_ontology(species_labels, y_vector)

        instant_core: Dict[str, float] = {}
        instant_edge: Dict[str, float] = {}

        for rxn in reactions:
            v = compute_mm_rate(rxn, conc_dict, enzyme_conc_map)
            if v == 0.0:
                continue
            for species_label, coeff in rxn.stoichiometry.items():
                lc = species_label.lower().strip()
                model_lc = alias_to_model_label.get(lc, lc)
                if model_lc in core_labels_lc:
                    instant_core[model_lc] = (
                        instant_core.get(model_lc, 0.0) + coeff * v
                    )
                elif model_lc in edge_labels_lc:
                    instant_edge[model_lc] = (
                        instant_edge.get(model_lc, 0.0) + coeff * v
                    )

        char_rate = (
            math.sqrt(sum(r * r for r in instant_core.values()))
            if instant_core
            else 0.0
        )
        return char_rate, instant_edge

    def _attach_interrupt_rates(
        self,
        sim_result: SimulationResult,
        species_labels: List[str],
        reactions: list,
        alias_to_model_label: Dict[str, str],
        enzyme_conc_map: Dict[str, float],
        core_labels_lc: set,
        edge_labels_lc: set,
    ) -> None:
        """
        Compute instantaneous core and edge rates at the interrupt time
        (last stored time point) and populate ``interrupt_char_rate`` and
        ``interrupt_edge_rates`` on the result.

        RMG promotes species based on the instantaneous flux ratio at the
        exact interrupt time, not trajectory peaks.
        """
        char_rate, instant_edge = self._compute_instantaneous_rates(
            sim_result.y[:, -1],
            species_labels, reactions, alias_to_model_label,
            enzyme_conc_map, core_labels_lc, edge_labels_lc,
        )
        sim_result.interrupt_char_rate = char_rate
        sim_result.interrupt_edge_rates = instant_edge

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

        if sim_result.peak_edge_signed_rate:
            return dict(sim_result.peak_edge_signed_rate)

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
            conc_dict[lc] = max(conc_dict.get(lc, 0.0), val)
            for equiv in get_ontology_equivalents(lab):
                eqlc = equiv.lower().strip()
                conc_dict[eqlc] = max(conc_dict.get(eqlc, 0.0), val)
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
