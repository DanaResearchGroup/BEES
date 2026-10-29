#!/usr/bin/env python3

"""ODE simulator: scipy solve_ivp with a vectorized RHS.

Concentrations mM, time s, rates mM/s.
"""

from dataclasses import dataclass, field
import copy
from types import SimpleNamespace
from typing import Dict, List, Optional, Tuple
import math
import os
import time
import numpy as np
from scipy.integrate import BDF, DOP853, RK23, RK45, Radau, solve_ivp
from scipy.sparse import csc_matrix

from bees.common import get_ontology_equivalents
from bees.core_edge_model import CoreEdgeModel, RateLawOptions
from bees.enzyme_sharing import PartitionCompetition, format_sharing_terms
from bees.flux_calculator import (
    compute_mm_rate_and_den,
    compute_reversible_mm_rate_and_den,
    format_rate_law,
    has_complete_explicit_product_kms,
    rate_law_spec,
)
from bees.conservator import Conservator
from bees.dosing import DosePlan, DoseSchedule
from bees.inhibition import (
    SUBSTRATE_SPECIFIC_TYPES,
    InhibitionTable,
    format_inhibition_terms,
    morrison_expression,
    resolve_slot,
)
from bees.reaction_generator import reaction_signature

_OLD_MODEL_NOTED = False

# Solvers that clip their step to t_bound, so one instance can be moved from check
# point to check point (keeps its Jacobian and step size). LSODA cannot.
_CONTINUABLE_SOLVERS = {"BDF": BDF, "Radau": Radau, "RK45": RK45, "RK23": RK23, "DOP853": DOP853}


def rate_law_options_of(model, logger=None) -> RateLawOptions:
    """The model's RateLawOptions; models pickled before the options existed mean all off."""
    global _OLD_MODEL_NOTED
    opts = getattr(model, "rate_law_options", None)
    if not isinstance(opts, RateLawOptions):
        if not _OLD_MODEL_NOTED and logger:
            logger.info("Model has no rate_law_options (older pickle): enzyme competition off.")
        _OLD_MODEL_NOTED = True
        return RateLawOptions()
    return opts


@dataclass
class SimulationResult:
    """ODE simulation output."""
    t: np.ndarray
    y: np.ndarray
    species_labels: List[str]
    success: bool = True
    message: str = ""
    max_char_rate: float = 0.0
    final_char_rate: float = 0.0
    max_edge_rate_ratio: Dict[str, float] = field(default_factory=dict)
    simulation_interrupted: bool = False
    interrupt_char_rate: float = 0.0
    interrupt_edge_rates: Dict[str, float] = field(default_factory=dict)
    # Signature is reaction_signature (includes enzyme).
    max_edge_reaction_dlnaccum: Dict[Tuple, float] = field(default_factory=dict)
    # "rate_ratio" or "conversion" when the stepwise pass stopped at a termination criterion.
    termination_reason: str = ""

    @property
    def terminated_early(self) -> bool:
        return bool(self.termination_reason)

def build_alias_to_model_label(species_labels: List[str]) -> Dict[str, str]:
    """Map ontology aliases to model species labels so stoichiometry lookups resolve correctly."""
    alias_to_model: Dict[str, str] = {}
    for lab in species_labels:
        model_lc = lab.lower().strip()
        alias_to_model[model_lc] = model_lc
        for equiv in get_ontology_equivalents(lab):
            alias_to_model[equiv.lower().strip()] = model_lc
    return alias_to_model


class _VectorizedRHS:
    """NumPy-vectorized ODE RHS compiled once (no per-call Python loops)."""

    def __init__(
        self,
        species_labels: List[str],
        reactions: list,
        alias_to_model_label: Dict[str, str],
        enzyme_conc_map: Dict[str, float],
        constant_mask: np.ndarray,
        n_core_species: Optional[int] = None,
        n_core_reactions: Optional[int] = None,
        options: Optional[RateLawOptions] = None,
        dose_plan: Optional[DosePlan] = None,
    ):
        n_sp = len(species_labels)
        n_rx = len(reactions)
        label_to_idx = {lab.lower().strip(): i for i, lab in enumerate(species_labels)}

        vmax = np.zeros(n_rx, dtype=np.float64)
        reversible_indices: List[int] = []
        product_inhibited_indices: List[int] = []
        for j, rxn in enumerate(reactions):
            is_reversible = (
                getattr(rxn.template, "reversible", False)
                and getattr(rxn, "thermo", None) is not None
                and not rxn.thermo.irreversible
            )
            if is_reversible:
                reversible_indices.append(j)
                continue
            if has_complete_explicit_product_kms(rxn):
                product_inhibited_indices.append(j)
                continue  # leave vmax[j] = 0 so the overlay owns this column
            kin = rxn.kinetics
            if kin is None or rxn.rate_law is None:
                continue
            e_key = rxn.enzyme_label.lower().strip()
            e_conc = enzyme_conc_map.get(e_key, 0.0)
            if kin.kcat is not None and e_conc > 0:
                vmax[j] = kin.kcat * e_conc
            elif kin.vmax is not None:
                vmax[j] = kin.vmax

        sub_idx_list: List[List[int]] = []
        km_list:      List[List[float]] = []
        nu_list:      List[List[float]] = []
        slot_labels:  List[List[str]] = []

        for rxn in reactions:
            kin = rxn.kinetics
            pairs_idx: List[int] = []
            pairs_km:  List[float] = []
            pairs_nu:  List[float] = []
            pairs_lab: List[str] = []
            if kin is not None and rxn.rate_law is not None:
                km_per = getattr(kin, "km_per_substrate", None) or {}
                km_single = kin.km
                stoich = rxn.stoichiometry
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
                        pairs_nu.append(float(abs(stoich.get(reactant, 1))))
                        pairs_lab.append(reactant)
            slot_labels.append(pairs_lab)
            sub_idx_list.append(pairs_idx)
            km_list.append(pairs_km)
            nu_list.append(pairs_nu)

        K = max((len(idxs) for idxs in sub_idx_list), default=0)
        # Pad dummy conc 1.0 so (s/(Km+s))**0 = 1.
        self._sub_idx_mat = np.full((n_rx, K), fill_value=n_sp, dtype=np.int32)
        self._km_mat = np.zeros((n_rx, K), dtype=np.float64)
        self._nu_mat = np.zeros((n_rx, K), dtype=np.float64)
        for j, (idxs, kms, nus) in enumerate(zip(sub_idx_list, km_list, nu_list)):
            nj = len(idxs)
            if nj > 0:
                self._sub_idx_mat[j, :nj] = idxs
                self._km_mat[j, :nj] = kms
                self._nu_mat[j, :nj] = nus

        fb_idx_list: List[List[int]] = []
        fb_ki_list:  List[List[float]] = []
        fb_h_list:   List[List[float]] = []
        for rxn in reactions:
            inh = getattr(rxn, "feedback_inhibitors", None)
            if not isinstance(inh, dict):  # None, or MagicMock in tests
                inh = {}
            pi: List[int] = []
            pk: List[float] = []
            ph: List[float] = []
            for label_lc, params in inh.items():
                try:
                    ki, hill = float(params[0]), float(params[1])
                except (TypeError, IndexError, ValueError):
                    continue
                if ki <= 0:
                    continue
                model_lc = alias_to_model_label.get(label_lc, label_lc)
                idx = label_to_idx.get(model_lc, label_to_idx.get(label_lc))
                if idx is not None:
                    pi.append(idx)
                    pk.append(ki)
                    ph.append(hill)
            fb_idx_list.append(pi)
            fb_ki_list.append(pk)
            fb_h_list.append(ph)

        M = max((len(p) for p in fb_idx_list), default=0)
        self._has_feedback = M > 0 and any(fb_idx_list)
        # Feedback pad Ki=inf so (1.0/inf)^h = 0 -> factor 1.
        self._fb_idx_mat = np.full((n_rx, M), fill_value=n_sp, dtype=np.int32)
        self._fb_ki_mat = np.full((n_rx, M), fill_value=np.inf, dtype=np.float64)
        self._fb_hill_mat = np.ones((n_rx, M), dtype=np.float64)
        for j, (idxs, kis, hs) in enumerate(zip(fb_idx_list, fb_ki_list, fb_h_list)):
            nj = len(idxs)
            if nj > 0:
                self._fb_idx_mat[j, :nj] = idxs
                self._fb_ki_mat[j, :nj] = kis
                self._fb_hill_mat[j, :nj] = hs

        # Reversible inhibitors (rxn.inhibitors): each acts on one substrate slot of the law.
        self._slot_labels = slot_labels
        entries = []
        self.unresolved_inhibitors: List[Tuple[str, str, str]] = []
        for j, rxn in enumerate(reactions):
            specs = getattr(rxn, "inhibitors", None)
            if not isinstance(specs, (list, tuple)):  # None, or MagicMock in tests
                continue
            e_conc = enzyme_conc_map.get(str(rxn.enzyme_label).lower().strip(), 0.0)
            for sp in specs:
                model_lc = alias_to_model_label.get(sp.inhibitor, sp.inhibitor)
                idx = label_to_idx.get(model_lc, label_to_idx.get(sp.inhibitor))
                if idx is None:
                    self.unresolved_inhibitors.append((rxn.enzyme_label, sp.inhibitor, "not a model species"))
                    continue
                if sp.type in SUBSTRATE_SPECIFIC_TYPES:
                    slot = resolve_slot(sp.substrate, slot_labels[j], alias_to_model_label)
                    if slot is None:
                        if sp.substrate is not None and slot_labels[j]:
                            reason = f"{sp.substrate!r} is not a substrate with a Km in {rxn!r}"
                        else:
                            reason = f"no substrate with a Km in {rxn!r}"
                        self.unresolved_inhibitors.append((rxn.enzyme_label, sp.inhibitor, reason))
                        continue
                else:
                    slot = K
                entries.append((j, idx, slot, sp, e_conc))
        self._inhibition = InhibitionTable(n_rx, K, entries)
        self._inhibited_rx = set(self._inhibition.inhibited_reactions)

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

        if n_rx <= 200:
            self._S = S.toarray()
            self._sparse = False
        else:
            self._S = S
            self._sparse = True

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

        # Edge isolation: species rows 0, reaction fluxes zeroed; omit n_core = all core.
        self._n_core_species = n_sp if n_core_species is None else int(n_core_species)
        self._n_core_reactions = n_rx if n_core_reactions is None else int(n_core_reactions)
        self._edge_species_mask = np.zeros(n_sp, dtype=bool)
        if self._n_core_species < n_sp:
            self._edge_species_mask[self._n_core_species:] = True
        self._integration_mask = constant_mask | self._edge_species_mask

        self._reversible_indices = reversible_indices
        self._reversible_reactions = [reactions[j] for j in reversible_indices]
        self._product_inhibited_indices = product_inhibited_indices
        self._product_inhibited_reactions = [reactions[j] for j in product_inhibited_indices]
        self._species_labels = list(species_labels)
        self._enzyme_conc_map = dict(enzyme_conc_map)

        self._reactions = reactions
        self._dose_plan = dose_plan
        self.options = options or RateLawOptions()
        self._competition = None
        if self.options.enzyme_competition:
            self._competition = PartitionCompetition(
                reactions, self._n_core_reactions, species_labels, label_to_idx,
                alias_to_model_label, enzyme_conc_map,
                include_products=self.options.competition_products,
            )

    def __call__(self, t: float, y: np.ndarray) -> np.ndarray:
        """Evaluate dydt = S @ v(y), plus the dosing sources and eliminations at time t."""
        dydt = self.compute_dydt(y)
        if self._dose_plan is not None and self._dose_plan.active:
            dose = self._dose_plan.dydt(t, y)
            dose[self._integration_mask] = 0.0
            dydt = dydt + dose
        return dydt

    def compute_v(self, y: np.ndarray) -> np.ndarray:
        """Compute rate vector v(y). Accepts (n_species,) or (n_species, n_timepoints)."""
        return self._rates(y)[0]

    def compute_v_and_den(self, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Rates v(y) and each reaction's rate-law denominator D (same shape as v)."""
        return self._rates(y, want_den=True)

    def competition_factor(self, y: np.ndarray) -> np.ndarray:
        """c_r per reaction (ones when competition is off)."""
        y_ext, den, shape = self._ext_den(y)
        return self._competition_c(y_ext, den, shape)

    def competition_contributions(self, y: np.ndarray) -> Dict[str, Dict[str, float]]:
        """{enzyme: {complex: weight in Z}} for a 1-D state (partition competition only)."""
        if self._competition is None:
            return {}
        return self._competition.contributions(np.append(y, 1.0))

    def competition_breakdown(self, y: np.ndarray) -> Dict[int, dict]:
        """Per competing reaction: enzyme, own D_r and the other complexes' weights in Z_r."""
        if self._competition is None:
            return {}
        y_ext, den, _ = self._ext_den(y)
        return self._competition.breakdown(y_ext, den)

    def _ext_den(self, y: np.ndarray):
        y_ext = np.vstack([y, np.ones((1, y.shape[1]))]) if y.ndim > 1 else np.append(y, 1.0)
        v, den = self._rates(y, want_den=True, apply_sharing=False)
        return y_ext, den, v.shape

    def _competition_c(self, y_ext, den, shape) -> np.ndarray:
        if self._competition is not None:
            return self._competition.factor(y_ext, den)
        return np.ones(shape)

    def _rates(
        self, y: np.ndarray, want_den: bool = False, apply_sharing: bool = True
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        is_mat = (y.ndim > 1)
        sharing = apply_sharing and self._competition is not None
        want_den = want_den or sharing

        # NEVER clip y globally (breaks S@v mass balance). Dummy 1.0 is the pad sentinel.
        if is_mat:
            ones = np.ones((1, y.shape[1]), dtype=y.dtype)
            y_ext = np.vstack([y, ones])
        else:
            y_ext = np.append(y, 1.0)

        # Inhibitors rescale each reaction's Vmax and the Km of the substrate they act on.
        v_scale = km_scale = None
        if self._inhibition.active:
            v_scale, km_scale = self._inhibition.scales(y_ext)

        # Clip only inside saturation so rates stay non-negative.
        s = y_ext[self._sub_idx_mat]
        s_pos = np.maximum(s, 0.0)

        if is_mat:
            km = self._km_mat[:, :, np.newaxis]
            nu = self._nu_mat[:, :, np.newaxis]
            vmax = self._vmax[:, np.newaxis]
        else:
            km = self._km_mat
            nu = self._nu_mat
            vmax = self._vmax
        if v_scale is not None:
            km = km * km_scale
            vmax = vmax * v_scale

        sat = s_pos / (km + s_pos)
        v = vmax * np.prod(sat ** nu, axis=1)

        den = None
        if want_den:
            # Pad columns have Km = 0, nu = 0: (1 + 1/0)**0 = 1.
            with np.errstate(divide="ignore"):
                den = np.prod((1.0 + s_pos / km) ** nu, axis=1)

        # Overlays apply max(c, 0) themselves.
        _need_overlay = self._reversible_indices or self._product_inhibited_indices
        if _need_overlay and not is_mat:
            conc = {
                lab.lower().strip(): float(y[i])
                for i, lab in enumerate(self._species_labels)
            }

        if self._reversible_indices:
            if is_mat:
                n_t = y.shape[1]
                for col, rxn in zip(self._reversible_indices, self._reversible_reactions):
                    for k in range(n_t):
                        conc_k = {
                            lab.lower().strip(): float(y[i, k])
                            for i, lab in enumerate(self._species_labels)
                        }
                        v[col, k], d = compute_reversible_mm_rate_and_den(
                            self._inhibited_view(col, rxn, v_scale, km_scale, k), conc_k, self._enzyme_conc_map
                        )
                        if den is not None:
                            den[col, k] = d
            else:
                for col, rxn in zip(self._reversible_indices, self._reversible_reactions):
                    v[col], d = compute_reversible_mm_rate_and_den(
                        self._inhibited_view(col, rxn, v_scale, km_scale), conc, self._enzyme_conc_map
                    )
                    if den is not None:
                        den[col] = d

        if self._product_inhibited_indices:
            if is_mat:
                n_t = y.shape[1]
                for col, rxn in zip(self._product_inhibited_indices, self._product_inhibited_reactions):
                    for k in range(n_t):
                        conc_k = {
                            lab.lower().strip(): float(y[i, k])
                            for i, lab in enumerate(self._species_labels)
                        }
                        v[col, k], d = compute_mm_rate_and_den(
                            self._inhibited_view(col, rxn, v_scale, km_scale, k), conc_k, self._enzyme_conc_map
                        )
                        if den is not None:
                            den[col, k] = d
            else:
                for col, rxn in zip(self._product_inhibited_indices, self._product_inhibited_reactions):
                    v[col], d = compute_mm_rate_and_den(
                        self._inhibited_view(col, rxn, v_scale, km_scale), conc, self._enzyme_conc_map
                    )
                    if den is not None:
                        den[col] = d

        # Feedback inhibitors compete for one site: 1 / (1 + sum_j ([I_j]/Ki_j)^h_j).
        if self._has_feedback:
            c = np.maximum(y_ext[self._fb_idx_mat], 0.0)
            if is_mat:
                ki = self._fb_ki_mat[:, :, np.newaxis]
                hill = self._fb_hill_mat[:, :, np.newaxis]
            else:
                ki = self._fb_ki_mat
                hill = self._fb_hill_mat
            v = v / (1.0 + np.sum((c / ki) ** hill, axis=1))

        if sharing:
            v = v * self._competition_c(y_ext, den, v.shape)

        return v, den

    def _inhibited_view(self, col: int, rxn, v_scale, km_scale, k: Optional[int] = None):
        """``rxn`` with its Vmax and substrate Kms replaced by the inhibited ones (the same
        reaction when nothing inhibits it), for rate laws that read them off the reaction."""
        if v_scale is None or col not in self._inhibited_rx:
            return rxn
        v = float(v_scale[col] if k is None else v_scale[col, k])
        km_row = km_scale[col] if k is None else km_scale[col, :, k]
        kin = copy.copy(rxn.kinetics)
        if kin.kcat is not None:
            kin.kcat = kin.kcat * v
        if kin.vmax is not None:
            kin.vmax = kin.vmax * v
        km_per = dict(getattr(kin, "km_per_substrate", None) or {})
        for slot, lab in enumerate(self._slot_labels[col]):
            km_per[lab] = float(self._km_mat[col, slot] * km_row[slot])
        kin.km_per_substrate = km_per
        view = copy.copy(rxn)
        view.kinetics = kin
        return view

    def compute_dydt(self, y: np.ndarray) -> np.ndarray:
        """dydt = S @ v(y) with edge isolation for the integrator."""
        v = self.compute_v(y)
        if self._n_core_reactions < self._n_rx:
            # Zero edge reactions on a COPY so compute_v's result stays intact.
            v = v.copy()
            if v.ndim == 1:
                v[self._n_core_reactions:] = 0.0
            else:
                v[self._n_core_reactions:, :] = 0.0
        if self._sparse:
            dydt = self._S.dot(v)
            if hasattr(dydt, "toarray"):
                dydt = dydt.toarray()
            if dydt.ndim > 1 and v.ndim == 1:
                dydt = dydt.ravel()
        else:
            dydt = self._S @ v

        if dydt.ndim > 1:
            dydt[self._integration_mask, :] = 0.0
        else:
            dydt[self._integration_mask] = 0.0
        return dydt

    def compute_dydt_unmasked(self, y: np.ndarray) -> np.ndarray:
        """dydt = S @ v(y) without edge isolation. Constant-mask still applied."""
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
    """Integrate the biochemical reaction network for a CoreEdgeModel."""

    def __init__(self, model: CoreEdgeModel, logger=None, dosing: Optional[List[DoseSchedule]] = None):
        self.model = model
        self.logger = logger
        self.options = rate_law_options_of(model, logger)
        self.dosing: List[DoseSchedule] = list(dosing or [])

    def _dose_plan(self, species_labels: List[str], t_offset: float) -> Optional[DosePlan]:
        """The dosing schedules laid over this run (run time 0 = absolute time ``t_offset``)."""
        if not self.dosing:
            return None
        label_to_idx = {lab.lower().strip(): i for i, lab in enumerate(species_labels)}
        plan = DosePlan(self.dosing, label_to_idx, len(species_labels), t_offset)
        if plan.missing and self.logger:
            self.logger.warning(
                f"Dosing: species not in the model, so their doses are ignored: {', '.join(plan.missing)}"
            )
        return plan if plan.active else None

    @staticmethod
    def _concentrations_nonnegative(y: np.ndarray) -> np.ndarray:
        """Project concentrations (mM) to non-negative values, in-place safe."""
        out = np.asarray(y, dtype=float)
        return np.maximum(out, 0.0)

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
        termination_rate_ratio: Optional[float] = None,
        termination_conversion: Optional[Dict[str, float]] = None,
        t_offset: float = 0.0,
    ) -> SimulationResult:
        """Run an ODE simulation of the full model (core + edge species and reactions).

        ``t_offset`` is the absolute time at which this run starts (a run resumed from the
        model's current concentrations); dosing schedules are laid over absolute time.

        termination_rate_ratio / termination_conversion stop the stepwise pass as soon
        as R_char/peak drops below the ratio or a conversion target is met;
        they have no effect on the continuous path.
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
                termination_rate_ratio=termination_rate_ratio,
                termination_conversion=termination_conversion,
                t_offset=t_offset,
            )
        return self._simulate_continuous(
            end_time=end_time,
            time_step=time_step,
            method=method,
            rtol=rtol,
            atol=atol,
            t_offset=t_offset,
        )

    def _simulate_continuous(
        self,
        end_time: float,
        time_step: Optional[float],
        method: str,
        rtol: float,
        atol: float,
        t_offset: float = 0.0,
    ) -> SimulationResult:
        """Full-span solve_ivp without flux-interrupt events."""
        species_labels = self.model.get_all_species_labels()
        n_species = len(species_labels)
        alias_to_model_label = build_alias_to_model_label(species_labels)
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
        plan = self._dose_plan(species_labels, t_offset)
        if plan is not None:
            y0 = y0 + plan.y0_add

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
            n_core_species=len(self.model.core_species),
            n_core_reactions=len(self.model.core_reactions),
            options=self.options,
            dose_plan=plan,
        )
        conservator = Conservator(species_labels)

        if self.logger:
            self.logger.debug(
                f"Running ODE simulation (continuous): {n_species} species, "
                f"{len(reactions)} reactions, t=[0, {end_time}] s"
            )

        max_step = end_time / 10 if end_time > 0 else np.inf
        if plan is None or not (plan.breakpoints or plan.events):
            sol = solve_ivp(
                rhs, (0.0, end_time), y0,
                method=method, t_eval=t_eval,
                rtol=rtol, atol=atol,
                max_step=max_step,
            )
        else:
            sol = self._solve_segments(rhs, plan, y0, end_time, t_eval, method, rtol, atol, max_step)

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
            conservator.clip_2d(y_raw) if y_raw.size else y_raw
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

    @staticmethod
    def _solve_segments(rhs, plan: DosePlan, y0, end_time, t_eval, method, rtol, atol, max_step):
        """solve_ivp in pieces that end at every dosing breakpoint, so the integrator never
        steps over a change of the source; impulses (boluses) are added between pieces."""
        cuts = [0.0] + [p for p in plan.breakpoints if p < end_time] + [end_time]
        t_all: List[np.ndarray] = []
        y_all: List[np.ndarray] = []
        y = y0.copy()
        success, message = True, ""
        for a, b in zip(cuts[:-1], cuts[1:]):
            for idx, amount in plan.impulses_at(a):
                y[idx] += amount
            if t_all and t_all[-1].size and t_all[-1][-1] == a and plan.impulses_at(a):
                y_all[-1][:, -1] = y  # the state at a is the post-impulse one
            sol = solve_ivp(
                rhs, (a, b), y, method=method, rtol=rtol, atol=atol,
                max_step=max_step, dense_output=t_eval is not None,
            )
            if not sol.success:
                success, message = False, sol.message
                t_all.append(sol.t)
                y_all.append(sol.y)
                break
            if t_eval is None:
                t_seg, y_seg = sol.t, sol.y
            else:
                grid = t_eval[(t_eval >= a) & (t_eval <= b)] if a == 0.0 else t_eval[(t_eval > a) & (t_eval <= b)]
                t_seg, y_seg = grid, (sol.sol(grid) if len(grid) else np.empty((len(y), 0)))
            if t_all and t_eval is None:
                t_seg, y_seg = t_seg[1:], y_seg[:, 1:]
            t_all.append(t_seg)
            y_all.append(y_seg)
            y = sol.y[:, -1].copy()
        return SimpleNamespace(
            t=np.concatenate(t_all), y=np.concatenate(y_all, axis=1), success=success, message=message,
        )

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
        termination_rate_ratio: Optional[float] = None,
        termination_conversion: Optional[Dict[str, float]] = None,
        t_offset: float = 0.0,
    ) -> SimulationResult:
        """Advance by dt; interrupt when any edge flux ratio exceeds tolerance (earliest at t > 0).

        Termination criteria are checked after the interrupt checks at every outer step.
        """
        species_labels = self.model.get_all_species_labels()
        n_species = len(species_labels)
        label_to_idx = {lab.lower().strip(): i for i, lab in enumerate(species_labels)}
        alias_to_model_label = build_alias_to_model_label(species_labels)
        enzyme_conc_map = self._build_enzyme_concentration_map()
        constant_mask = np.array(
            [s.constant for s in self.model.core_species]
            + [False] * len(self.model.edge_species),
            dtype=bool,
        )
        y0 = self._concentrations_nonnegative(
            np.array(self.model.get_all_concentration_vector(), dtype=float)
        )
        plan = self._dose_plan(species_labels, t_offset)
        if plan is not None:
            y0 = y0 + plan.y0_add
        reactions = self.model.core_reactions + self.model.edge_reactions
        edge_labels_lc_list = [sp.label.lower().strip() for sp in self.model.edge_species]
        n_core = len(self.model.core_species)
        n_core_rxns = len(self.model.core_reactions)

        rhs = _VectorizedRHS(
            species_labels=species_labels,
            reactions=reactions,
            alias_to_model_label=alias_to_model_label,
            enzyme_conc_map=enzyme_conc_map,
            constant_mask=constant_mask,
            n_core_species=n_core,
            n_core_reactions=n_core_rxns,
            options=self.options,
            dose_plan=plan,
        )
        conservator = Conservator(species_labels)

        # dlnaccum signature is reaction_signature (includes enzyme).
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
                mult = int(round(abs(float(coeff))))
                if mult <= 0:
                    mult = 1
                if float(coeff) < 0:
                    r_idx.extend([idx] * mult)
                elif float(coeff) > 0:
                    p_idx.extend([idx] * mult)
            sig = reaction_signature(rxn)
            edge_rxn_sigs.append(sig)
            edge_rxn_reactant_idx.append(r_idx)
            edge_rxn_product_idx.append(p_idx)
            edge_rxn_global_idx.append(j_global)

        max_edge_rxn_dlnaccum: Dict[Tuple, float] = {sig: 0.0 for sig in edge_rxn_sigs}

        conversion_targets: List[Tuple[str, int, float, float]] = []
        for sp_label, target_frac in (termination_conversion or {}).items():
            idx = label_to_idx.get(sp_label.lower().strip())
            if idx is not None and y0[idx] > 0.0:
                conversion_targets.append((sp_label, idx, float(y0[idx]), float(target_frac)))
        termination_reason = ""

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
        dt_max = 5.0
        if end_time > 0.0:
            dt_floor_outer = max(1e-12, min(1e-6, end_time * 1e-8))
        else:
            dt_floor_outer = 1e-12
        dt = dt_floor_outer
        dt_min = dt_floor_outer
        dt_growth = 2.0
        near_frac = 0.5
        # Peak includes t=0 (often the peak), like _attach_flux_metrics over the trajectory.
        max_char_rate = float(np.linalg.norm(rhs.compute_dydt_unmasked(y0)[:n_core]))
        interrupted = False
        interrupt_char = 0.0
        interrupt_edge: Dict[str, float] = {}
        first_step = True
        outer_step_idx = 0
        state_invalid = False
        invalid_message = ""
        solver_failed = False
        solver_message = ""
        solver_cls = _CONTINUABLE_SOLVERS.get(method)
        solver = None

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
            # Never step over a change of the dosing source; restart the solver after it.
            at_dose_break = False
            if plan is not None:
                t_break = plan.next_break(current_t, end_time)
                if t_break <= target_t:
                    target_t, at_dose_break = t_break, t_break < end_time
            if target_t <= current_t:
                break

            try:
                if solver_cls is not None:
                    if solver is None:
                        solver = solver_cls(
                            rhs, current_t, current_y, t_bound=target_t, rtol=rtol, atol=atol,
                        )
                    else:
                        solver.t_bound = target_t
                        solver.status = "running"
                    message = None
                    while solver.status == "running":
                        message = solver.step()
                    ok = solver.status == "finished"
                    new_t, new_y = solver.t, solver.y
                else:
                    sol = solve_ivp(
                        rhs, (current_t, target_t), current_y,
                        method=method, rtol=rtol, atol=atol,
                        dense_output=False,
                        max_step=target_t - current_t,
                    )
                    ok = bool(sol.success) and sol.y.shape[1] > 0
                    message = sol.message
                    if ok:
                        new_t, new_y = sol.t[-1], sol.y[:, -1]
            except Exception as exc:
                solver_failed = True
                solver_message = f"ODE solver exception at t={current_t:.6e} s: {exc!r}"
                if self.logger:
                    self.logger.warning(solver_message)
                break

            if not ok:
                solver_failed = True
                solver_message = (
                    f"ODE solver failed at t={current_t:.6e} s: {message or 'unknown'}"
                )
                if self.logger:
                    self.logger.warning(solver_message)
                break

            current_t = float(new_t)
            current_y = np.array(new_y, dtype=float)
            if at_dose_break:
                for idx, amount in plan.impulses_at(current_t):
                    current_y[idx] += amount
                solver = None
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
            # Do NOT clip current_y mid-integration (BDF conservation drift).

            t_history.append(current_t)
            y_history.append(conservator.clip(current_y).copy())

            # Enlarger uses UNMASKED residual (bootstrap: isolated R_char=0).
            dydt_unmasked_step = rhs.compute_dydt_unmasked(current_y)
            char_rate = float(np.linalg.norm(dydt_unmasked_step[:n_core]))
            if edge_labels_lc_list:
                edge_rates = {
                    lab: float(dydt_unmasked_step[n_core + i])
                    for i, lab in enumerate(edge_labels_lc_list)
                }
            else:
                edge_rates = {}
            max_char_rate = max(max_char_rate, char_rate)

            max_rr = 0.0
            if char_rate > 0.0 and edge_rates:
                max_rr = max((abs(r) / char_rate for r in edge_rates.values()), default=0.0)

            max_rxn_dlnaccum_step = 0.0
            max_rxn_dlnaccum_sig: Optional[Tuple] = None
            if edge_rxn_sigs:
                v_all = rhs.compute_v(current_y)
                prod_all, cons_all = rhs.compute_prod_cons(current_y)
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
                        f"after completing outer step at simulation t={current_t:.6e} s."
                    )
                break

            if not first_step and max_rr > interrupt_simulation_tol:
                interrupted = True
                interrupt_char = char_rate
                interrupt_edge = edge_rates
                if self.logger:
                    self.logger.info(
                        f"ODE interrupted at t={current_t:.6e} s "
                        f"(max_rr={max_rr:.4e}, "
                        f"tol={interrupt_simulation_tol}, "
                        f"{len(t_history)} time points)"
                    )
                break

            if (
                not first_step
                and tol_move_edge_reaction_to_core is not None
                and tol_move_edge_reaction_to_core > 0
                and max_rxn_dlnaccum_step > tol_move_edge_reaction_to_core
            ):
                interrupted = True
                interrupt_char = char_rate
                interrupt_edge = edge_rates
                if self.logger:
                    self.logger.info(
                        f"ODE interrupted (reaction dlnaccum) at t={current_t:.6e} s "
                        f"(max_dlnaccum={max_rxn_dlnaccum_step:.4e}, "
                        f"tol={tol_move_edge_reaction_to_core}, "
                        f"sig={max_rxn_dlnaccum_sig})"
                    )
                break

            if not first_step:
                if (
                    termination_rate_ratio is not None
                    and max_char_rate > 0.0
                    and char_rate / max_char_rate < termination_rate_ratio
                ):
                    termination_reason = "rate_ratio"
                    detail = (
                        f"R_char/peak={char_rate / max_char_rate:.4e} "
                        f"< {termination_rate_ratio}"
                    )
                else:
                    for sp_label, idx, c0, target_frac in conversion_targets:
                        conversion = 1.0 - max(float(current_y[idx]), 0.0) / c0
                        if conversion >= target_frac:
                            termination_reason = "conversion"
                            detail = (
                                f"conversion of '{sp_label}'={conversion:.4f} "
                                f">= {target_frac}"
                            )
                            break
                if termination_reason:
                    if self.logger:
                        self.logger.info(
                            f"ODE terminated at t={current_t:.6e} s ({detail})"
                        )
                    break

            first_step = False

            # Do not shrink dt while max_rr <= tol.
            dt = max(dt, dt_floor_outer)
            if not max_rr >= near_frac * interrupt_simulation_tol:
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
        if solver_failed:
            fail_msg = f"{fail_msg + ' ' if fail_msg else ''}{solver_message}".strip()

        result = SimulationResult(
            t=t_arr,
            y=y_arr,
            species_labels=species_labels,
            success=(not wall_time_exceeded and not state_invalid and not solver_failed),
            message=fail_msg,
            simulation_interrupted=interrupted,
            max_char_rate=max_char_rate,
            termination_reason=termination_reason,
        )

        if result.success and y_arr.shape[1] > 0:
            self.model.set_all_concentrations(
                conservator.clip(current_y).tolist()
            )
        self._attach_flux_metrics(result, rhs=rhs)

        if interrupted:
            result.interrupt_char_rate = interrupt_char
            result.interrupt_edge_rates = interrupt_edge
        else:
            if self.logger and result.success:
                self.logger.debug(
                    f"ODE simulation completed (stepwise) "
                    f"({len(t_history)} time points)"
                )

        result.max_edge_reaction_dlnaccum = max_edge_rxn_dlnaccum

        return result

    def export_ode_equations(
        self,
        output_path: str,
        iteration: int,
        end_time: float,
    ) -> str:
        """Export ODE equations solved at this iteration to a text file."""

        species_labels = self.model.get_all_species_labels()
        reactions = self.model.core_reactions + self.model.edge_reactions
        label_to_idx = {lab.lower().strip(): i for i, lab in enumerate(species_labels)}
        edge_labels_lc = {sp.label.lower().strip() for sp in self.model.edge_species}
        enzyme_map = self._build_enzyme_concentration_map()
        shown = {lab.lower().strip(): lab for lab in species_labels}

        species_terms: Dict[str, List[Tuple[int, int]]] = {
            lab.lower().strip(): [] for lab in species_labels
        }
        constant_mask = {}
        for sp in self.model.core_species:
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

        lines.append("-" * 80)
        lines.append("Reactions (rate laws with parameters)")
        lines.append("-" * 80)

        _ode_sharing = None
        if self.options.enzyme_competition:
            _ode_sharing = PartitionCompetition(
                reactions,
                len(self.model.core_reactions),
                species_labels,
                label_to_idx,
                {lab.lower().strip(): lab.lower().strip() for lab in species_labels},
                enzyme_map,
                include_products=self.options.competition_products,
            )

        for r_idx, rxn in enumerate(reactions, 1):
            reactants = " + ".join(
                s for s, c in rxn.stoichiometry.items() if c < 0
            )
            products = " + ".join(
                s for s, c in rxn.stoichiometry.items() if c > 0
            )
            td = getattr(rxn, "thermo", None)
            spec = rate_law_spec(rxn, enzyme_map)
            _is_rev = spec.form == "reversible"
            arrow = "<=>" if _is_rev else "->"
            lines.append(f"R{r_idx}: {reactants} {arrow} {products}")
            lines.append(f"    Enzyme: {rxn.enzyme_label}  ({rxn.ec_number or 'N/A'})")

            kin = rxn.kinetics
            if kin is None or rxn.rate_law is None:
                lines.append("    v = (no kinetics - rate = 0)")
            else:
                kcat = kin.kcat
                vmax = kin.vmax
                km_per = getattr(kin, "km_per_substrate", None) or {}
                km_single = kin.km

                if spec.kcat is not None:
                    vf = f"kcat * [{rxn.enzyme_label}]"
                elif spec.vmax is not None:
                    vf = "Vmax"
                else:
                    vf = "0"
                others, k_of = None, {}
                if _ode_sharing is not None:
                    terms = _ode_sharing.other_terms(r_idx - 1)
                    k_of = {lab: k for _, ligands in terms for lab, _, k in ligands}
                    enzyme = rxn.enzyme_label
                    others = format_sharing_terms(
                        terms,
                        lambda lab: f"[{lab}]" if lab.lower().strip() in label_to_idx else None,
                        lambda lab, _e=enzyme: f"Kc({_e}, {lab})",
                    )
                inh_specs = getattr(rxn, "inhibitors", None)
                inh_specs = inh_specs if isinstance(inh_specs, (list, tuple)) else []
                inh_params: List[str] = []
                inh_v, inh_km = "1", {}
                if inh_specs and spec.form != "zero":
                    def _inh_param(kind, sp, value, _p=inh_params):
                        inh = shown.get(sp.inhibitor, sp.inhibitor)
                        name = f"{kind}({inh})" if sp.substrate is None else f"{kind}({inh}, {sp.substrate})"
                        _p.append(f"{name}={value:.6g} mM")
                        return name
                    first_sub = spec.substrates[0][0] if spec.substrates else None
                    inh_v, inh_km = format_inhibition_terms(
                        inh_specs,
                        lambda lab: f"[{shown[lab]}]" if lab in label_to_idx else None,
                        _inh_param,
                        first_sub,
                    )
                    for sp in inh_specs:
                        if sp.type == "tight_binding" and sp.inhibitor in label_to_idx:
                            ki_name = f"Ki({shown[sp.inhibitor]})"
                            inh_params.append(f"{ki_name}={sp.ki:.6g} mM")
                            inh_v += " * " + morrison_expression(f"[{shown[sp.inhibitor]}]", f"[{rxn.enzyme_label}]", ki_name)
                    if inh_v.startswith("1 * "):
                        inh_v = inh_v[4:]
                    if inh_v != "1":
                        vf = f"{vf} * {inh_v}"
                expr = format_rate_law(
                    spec,
                    conc=lambda lab: f"[{lab}]" if lab.lower().strip() in label_to_idx else None,
                    km=lambda lab, role: (
                        f"(Km_{lab} * {inh_km[lab.lower().strip()]})"
                        if role == "S" and lab.lower().strip() in inh_km else f"Km_{lab}"
                    ),
                    vf=vf,
                    others=others,
                )
                rate_str = f"v{r_idx} = {expr}"
                fb = getattr(rxn, "feedback_inhibitors", None)
                if expr != "0" and isinstance(fb, dict) and fb:
                    fb_terms = " + ".join(f"([{lab}]/Ki_{lab})^h_{lab}" for lab in fb)
                    rate_str = f"v{r_idx} = ({expr}) / (1 + {fb_terms})"

                lines.append(f"    {rate_str}")

                params = []
                if kcat is not None:
                    params.append(f"kcat={kcat:.6g} 1/s")
                    if getattr(kin, "kcat_sd", None) is not None:
                        params.append(f"kcat_sd={kin.kcat_sd:.6g} 1/s")
                if _is_rev:
                    params.append(f"Keq={spec.keq:.6g} mM^{spec.delta_nu:g}")
                if isinstance(fb, dict):
                    for lab, (ki_val, hill_val) in fb.items():
                        params.append(f"Ki({lab})={float(ki_val):.6g} mM, h({lab})={float(hill_val):g}")
                if vmax is not None:
                    params.append(f"Vmax={vmax:.6g} mM/s")
                for lab, k_can in k_of.items():
                    params.append(f"Kc({rxn.enzyme_label}, {lab})={k_can:.6g} mM")
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
                params.extend(inh_params)
                if params:
                    lines.append(f"    Parameters: {', '.join(params)}")
            if td is not None:
                thermo_parts = []
                if td.dgr_prime_kJmol is not None:
                    thermo_parts.append(f"ΔG°'={td.dgr_prime_kJmol:.2f} kJ/mol")
                    if td.sigma_kJmol is not None and math.isfinite(td.sigma_kJmol):
                        thermo_parts.append(f"σ={td.sigma_kJmol:.2f} kJ/mol")
                thermo_parts.append(f"source={td.source}")
                lines.append(f"    Thermodynamics: {', '.join(thermo_parts)}")
            lines.append("")

            for species_label, coeff in rxn.stoichiometry.items():
                lc = species_label.lower().strip()
                if lc in label_to_idx and lc in species_terms:
                    species_terms[lc].append((coeff, r_idx))

        lines.append("-" * 80)
        lines.append("ODEs (dC/dt for each species)")
        lines.append("-" * 80)

        dose_terms: Dict[str, List[str]] = {}
        eliminations: Dict[str, float] = {}
        for n, sch in enumerate(self.dosing, 1):
            lc = sch.species.lower().strip()
            if lc not in label_to_idx:
                continue
            if sch.profile != "bolus":
                dose_terms.setdefault(lc, []).append(f"+ u{n}(t)")
            if sch.elimination_k:
                eliminations[lc] = sch.elimination_k
        for lc, k_el in eliminations.items():
            dose_terms.setdefault(lc, []).append(f"- kel_{shown[lc]}*[{shown[lc]}]")

        for lab in species_labels:
            lc = lab.lower().strip()
            terms = species_terms.get(lc, [])
            is_edge = lc in edge_labels_lc
            is_const = constant_mask.get(lc, False)
            extra = " ".join(dose_terms.get(lc, []))

            if is_const:
                lines.append(f"d[{lab}]/dt = 0  (constant)")
            elif not terms and not extra:
                lines.append(f"d[{lab}]/dt = 0")
            else:
                term_strs = []
                for coeff, r_idx in terms:
                    if coeff > 0:
                        term_strs.append(f"+{coeff}*v{r_idx}")
                    else:
                        term_strs.append(f"{coeff}*v{r_idx}")
                rhs = " ".join(term_strs).lstrip("+") or "0"
                if extra:
                    rhs = f"{rhs} {extra}" if terms else extra.lstrip("+ ")
                suffix = "  (edge)" if is_edge else ""
                lines.append(f"d[{lab}]/dt = {rhs}{suffix}")

        if self.dosing:
            lines.append("")
            lines.append("-" * 80)
            lines.append("Dosing (additive sources u(t) in mM/s and first-order elimination)")
            lines.append("-" * 80)
            for n, sch in enumerate(self.dosing, 1):
                if sch.species.lower().strip() not in label_to_idx:
                    lines.append(f"u{n}: {sch.species} is not in the model - ignored")
                    continue
                params = ", ".join(
                    f"{k}={v:.6g}" if isinstance(v, float) else f"{k}=<{len(v)} rows>"
                    for k, v in sch.params.items()
                )
                lines.append(
                    f"u{n}(t): {sch.profile} release of {sch.dose:.6g} mM {sch.species} from "
                    f"t_start={sch.t_start:g} s" + (f" ({params})" if params else "")
                    + (" - an impulse at t_start, no rate term" if sch.profile == "bolus" else "")
                )
                if sch.elimination_k:
                    lines.append(
                        f"kel_{shown[sch.species.lower().strip()]} = {sch.elimination_k:.6g} 1/s "
                        f"(half-life {math.log(2.0) / sch.elimination_k:.6g} s)"
                    )
        lines.append("")
        lines.append("=" * 80)

        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        with open(output_path, "w") as f:
            f.write("\n".join(lines))

        if self.logger:
            self.logger.info(f"Exported ODE equations to {output_path}")

        return output_path

    def _attach_flux_metrics(self, sim_result: SimulationResult, rhs: _VectorizedRHS) -> None:
        """Populate max_char_rate, final_char_rate, and max_edge_rate_ratio from the trajectory."""
        if not sim_result.success or sim_result.y.size == 0:
            sim_result.max_char_rate = 0.0
            sim_result.final_char_rate = 0.0
            sim_result.max_edge_rate_ratio = {}
            return

        n_core = len(self.model.core_species)
        n_edge = len(self.model.edge_species)
        edge_labels = [sp.label.lower().strip() for sp in self.model.edge_species]

        # Enlarger uses UNMASKED residual (bootstrap: isolated R_char=0).
        dydt_unmasked = rhs.compute_dydt_unmasked(sim_result.y)
        dydt_core_mat = dydt_unmasked[:n_core, :]
        char_rates = np.linalg.norm(dydt_core_mat, axis=0)

        sim_result.max_char_rate = float(np.max(char_rates))
        sim_result.final_char_rate = float(char_rates[-1])

        if n_edge > 0:
            dydt_edge_mat = dydt_unmasked[n_core:, :]
            abs_edge_mat = np.abs(dydt_edge_mat)
            safe_char = char_rates.copy()
            safe_char[safe_char == 0] = 1e-100
            rr_mat = abs_edge_mat / safe_char
            max_ratios = np.max(rr_mat, axis=1)
            sim_result.max_edge_rate_ratio = {
                edge_labels[i]: float(max_ratios[i]) for i in range(n_edge)
            }
        else:
            sim_result.max_edge_rate_ratio = {}

    def _build_rhs(self, core_only: bool = False) -> "_VectorizedRHS":
        """The RHS of the model (core + edge, or the core reactions only)."""
        species_labels = self.model.get_all_species_labels()
        return _VectorizedRHS(
            species_labels=species_labels,
            reactions=self.model.core_reactions + ([] if core_only else self.model.edge_reactions),
            alias_to_model_label=build_alias_to_model_label(species_labels),
            enzyme_conc_map=self._build_enzyme_concentration_map(),
            constant_mask=np.array(
                [s.constant for s in self.model.core_species] + [False] * len(self.model.edge_species),
                dtype=bool,
            ),
            n_core_species=len(self.model.core_species),
            n_core_reactions=len(self.model.core_reactions),
            options=self.options,
        )

    def unmatched_feedback_inhibitors(self) -> Dict[str, List[str]]:
        """{enzyme: [inhibitor labels]} of core reactions whose inhibitor is not a model species.

        The RHS drops such inhibitors from the feedback term, so the run would otherwise
        proceed with that feedback silently off.
        """
        species_labels = self.model.get_all_species_labels()
        label_to_idx = {lab.lower().strip(): i for i, lab in enumerate(species_labels)}
        alias = build_alias_to_model_label(species_labels)
        missing: Dict[str, List[str]] = {}
        for rxn in self.model.core_reactions:
            inh = getattr(rxn, "feedback_inhibitors", None)
            if not isinstance(inh, dict):
                continue
            for label_lc in inh:
                if alias.get(label_lc, label_lc) in label_to_idx or label_lc in label_to_idx:
                    continue
                labels = missing.setdefault(str(rxn.enzyme_label), [])
                if label_lc not in labels:
                    labels.append(label_lc)
        return missing

    def _build_enzyme_concentration_map(self) -> Dict[str, float]:
        """Return enzyme label (lc) -> concentration (mM) for all enzyme core species."""
        enzyme_map: Dict[str, float] = {}
        for sp in self.model.core_species:
            if sp.is_enzyme:
                enzyme_map[sp.label.lower().strip()] = sp.concentration
        return enzyme_map
