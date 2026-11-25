"""
Tests for FAS-II parameter sampling (inventory, tied Km, and draw indexing).

To run:  pytest -v tests/test_uncertainty.py

Covers: exact log10 SD recovery, reverse-query Km SD retention, joint ΔG°′ covariance
assembly, rule-baseline isolation, and index-based deterministic sampling.
"""

from __future__ import annotations

import math
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace
from typing import Optional
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from bees.common import linear_sd_to_log10_sd, log10_sd_to_linear_sd
from bees.reaction_generator import ReactionGenerator
from bees.rules import Reference, RuleRegistry
from bees.rules.calibrations.fas import MeasuredKcatAnchor, TesALongChainPreference
from bees.rules.physics_rules import DgrIrreversibility
from bees.thermodynamics import ThermoEngine, _DiskCache


# ---------------------------------------------------------------------------
# log10 SD recovery
# ---------------------------------------------------------------------------

class TestLog10SdRecovery:
    @pytest.mark.parametrize("mean", [1e-6, 0.0356, 1.0, 15.67, 4.2e3])
    @pytest.mark.parametrize("sd_log10", [0.05, 0.35, 0.7, 1.06, 2.0])
    def test_round_trip_exact(self, mean, sd_log10):
        sd_lin = log10_sd_to_linear_sd(mean, sd_log10)
        assert linear_sd_to_log10_sd(mean, sd_lin) == pytest.approx(sd_log10, rel=1e-10)

    def test_degenerate_inputs_return_zero(self):
        assert linear_sd_to_log10_sd(0.0, 1.0) == 0.0
        assert linear_sd_to_log10_sd(1.0, 0.0) == 0.0
        assert linear_sd_to_log10_sd(-1.0, 1.0) == 0.0

    def test_post_rule_mean_gives_wrong_sd(self):
        # The inverse needs the mean the SD was computed with (pre-rule CatPred value).
        sd_lin = log10_sd_to_linear_sd(0.04, 0.6)
        assert linear_sd_to_log10_sd(0.04 * 1e-3, sd_lin) > 1.0


# ---------------------------------------------------------------------------
# Reverse-query Km SD retention
# ---------------------------------------------------------------------------

def _generator(est_km, est_sd):
    bees = MagicMock()
    bees.enzymes = []
    bees.species = []
    bees.settings = MagicMock()
    gen = ReactionGenerator(bees, MagicMock(), output_directory="/tmp")
    est = SimpleNamespace(
        km_per_substrate=dict(est_km),
        km_sd_per_substrate=dict(est_sd),
        compound_smiles={},
    )
    gen.kinetics_estimator = MagicMock()
    gen.kinetics_estimator.estimate = MagicMock(return_value=est)
    return gen, est


def _fabg_like():
    return SimpleNamespace(
        enzyme_label="FabG",
        substrate_label="3-oxobutanoyl-[ACP]",
        product_labels=["(3R)-hydroxybutanoyl-[ACP]", "NADP"],
        reactant_labels=["3-oxobutanoyl-[ACP]", "NADPH", "H+"],
        stoichiometry={
            "3-oxobutanoyl-[ACP]": -1, "NADPH": -1, "H+": -1,
            "(3R)-hydroxybutanoyl-[ACP]": 1, "NADP": 1,
        },
        ec_number="EC 1.1.1.100",
    )


_SMILES = {
    "(3R)-hydroxybutanoyl-[ACP]": "CC(O)CC(=O)S",
    "NADP": "NC(=O)C1=C[N+](=CC=C1)C1OC(COP(=O)(O)OP(=O)(O)OCC2OC("
            "N3C=NC4=C(N)N=CN=C43)C(OP(=O)(O)O)C2O)C(O)C1O",
    "H+": "[H+]",
}


class TestReverseKmSd:
    def test_lookup_fills_product_sds(self):
        gen, _ = _generator(
            {"(3R)-hydroxybutanoyl-[ACP]": 0.2, "NADP": 0.26},
            {"(3R)-hydroxybutanoyl-[ACP]": 1.0, "NADP": 0.84},
        )
        rxn = _fabg_like()
        sds: dict = {}
        kms, _ = gen._lookup_product_kms_via_reverse_query(
            reaction=rxn, stoich=rxn.stoichiometry, ec_numbers_to_try=["EC 1.1.1.100"],
            temp_range=None, ph_range=None, provided_species_labels_lc=set(),
            enzyme_sequence="MNF", smiles_map=_SMILES, substitutor=None,
            product_km_sds=sds,
        )
        assert kms == {"(3R)-hydroxybutanoyl-[ACP]": 0.2, "NADP": 0.26}
        assert sds == {"(3R)-hydroxybutanoyl-[ACP]": 1.0, "NADP": 0.84}

    def test_lookup_without_out_param_keeps_two_tuple(self):
        gen, _ = _generator({"NADP": 0.26}, {"NADP": 0.84})
        rxn = _fabg_like()
        out = gen._lookup_product_kms_via_reverse_query(
            reaction=rxn, stoich=rxn.stoichiometry, ec_numbers_to_try=["EC 1.1.1.100"],
            temp_range=None, ph_range=None, provided_species_labels_lc=set(),
            enzyme_sequence="MNF", smiles_map=_SMILES, substitutor=None,
        )
        assert len(out) == 2

    def test_attach_merges_sds_into_new_dict(self):
        gen, _ = _generator({}, {})
        rxn = _fabg_like()
        fwd_km = {"3-oxobutanoyl-[ACP]": 0.087, "NADPH": 0.12}
        fwd_sd = {"3-oxobutanoyl-[ACP]": 0.94, "NADPH": 0.49}
        rxn.kinetics = SimpleNamespace(
            km_per_substrate=fwd_km, km_sd_per_substrate=fwd_sd, compound_smiles={},
        )

        def _fake_lookup(**kwargs):
            kwargs["product_km_sds"].update({"NADP": 0.84, "(3R)-hydroxybutanoyl-[ACP]": 1.0})
            return {"NADP": 0.26, "(3R)-hydroxybutanoyl-[ACP]": 0.2}, {}

        gen._lookup_product_kms_via_reverse_query = _fake_lookup
        engine = MagicMock()
        engine.compute_keq = MagicMock(return_value=SimpleNamespace(irreversible=False))
        gen._attach_raw_thermo_eagerly(
            rxn, engine, global_smiles_map={}, ec_numbers_to_try=["EC 1.1.1.100"],
            temp_range=None, ph_range=None, provided_species_labels_lc=set(),
        )
        kin = rxn.kinetics
        assert kin.km_per_substrate["NADP"] == 0.26
        assert kin.km_sd_per_substrate["NADP"] == 0.84
        assert kin.km_sd_per_substrate["(3R)-hydroxybutanoyl-[ACP]"] == 1.0
        assert kin.km_sd_per_substrate["NADPH"] == 0.49
        # The forward SD dict may be the CatPred memo's object: must not be mutated.
        assert fwd_sd == {"3-oxobutanoyl-[ACP]": 0.94, "NADPH": 0.49}
        assert kin.km_sd_per_substrate is not fwd_sd

    def test_attach_does_not_borrow_sd_for_existing_km(self):
        gen, _ = _generator({}, {})
        rxn = _fabg_like()
        rxn.kinetics = SimpleNamespace(
            km_per_substrate={"NADP": 0.5}, km_sd_per_substrate={}, compound_smiles={},
        )

        def _fake_lookup(**kwargs):
            kwargs["product_km_sds"]["NADP"] = 0.84
            return {"NADP": 0.26}, {}

        gen._lookup_product_kms_via_reverse_query = _fake_lookup
        engine = MagicMock()
        gen._attach_raw_thermo_eagerly(
            rxn, engine, global_smiles_map={}, ec_numbers_to_try=[],
            temp_range=None, ph_range=None, provided_species_labels_lc=set(),
        )
        assert rxn.kinetics.km_per_substrate["NADP"] == 0.5
        assert "NADP" not in (rxn.kinetics.km_sd_per_substrate or {})


# ---------------------------------------------------------------------------
# Joint ΔG°′ covariance
# ---------------------------------------------------------------------------

class _Q:
    def __init__(self, arr):
        self._arr = np.asarray(arr, dtype=float)

    def m_as(self, _unit):
        return self._arr


class TestJointDgr:
    def _engine(self, tmp_path, statuses, values, cov):
        engine = ThermoEngine(pH=7.4, cache=_DiskCache(tmp_path))
        cc = MagicMock()
        cc.standard_dg_prime_multi = MagicMock(return_value=(_Q(values), _Q(cov)))
        engine._cc = cc
        it = iter(statuses)
        engine._build_cc_reaction = MagicMock(side_effect=lambda s, m: next(it))
        return engine, cc

    def test_shapes_and_placement(self, tmp_path):
        statuses = [("ok", "r0"), ("null", None), (None, None), ("ok", "r3")]
        cov = [[4.0, 1.5], [1.5, 9.0]]
        engine, cc = self._engine(tmp_path, statuses, [-10.0, 5.0], cov)
        with patch.dict(os.environ, {"BEES_DISABLE_THERMO": "0"}):
            dg, c, ok = engine.joint_dgr_prime([{}] * 4, [{}] * 4)
        assert dg.shape == (4,) and c.shape == (4, 4) and ok.shape == (4,)
        assert dg[0] == -10.0 and dg[3] == 5.0
        assert dg[1] == 0.0 and math.isnan(dg[2])
        assert list(ok) == [True, True, False, True]
        assert c[0, 0] == 4.0 and c[3, 3] == 9.0 and c[0, 3] == 1.5 and c[3, 0] == 1.5
        assert np.all(c[1] == 0) and np.all(c[2] == 0)
        args, kwargs = cc.standard_dg_prime_multi.call_args
        assert args[0] == ["r0", "r3"]
        assert kwargs.get("uncertainty_representation") == "cov"

    def test_disabled_returns_none(self, tmp_path):
        engine = ThermoEngine(pH=7.4, cache=_DiskCache(tmp_path))
        with patch.dict(os.environ, {"BEES_DISABLE_THERMO": "1"}):
            assert engine.joint_dgr_prime([{}], [{}]) is None

    def test_scalar_path_unchanged_for_null_reaction(self, tmp_path):
        engine = ThermoEngine(pH=7.4, cache=_DiskCache(tmp_path))
        engine._cc = MagicMock()
        engine._build_cc_reaction = MagicMock(return_value=("null", None))
        assert engine._compute_dgr_prime({"A": -1, "B": 1}, {}) == (0.0, 0.0)
        engine._build_cc_reaction = MagicMock(return_value=(None, None))
        assert engine._compute_dgr_prime({"A": -1, "B": 1}, {}) == (None, None)


# ---------------------------------------------------------------------------
# Rule-baseline isolation
# ---------------------------------------------------------------------------

@dataclass
class _Kin:
    kcat: Optional[float] = None
    km: Optional[float] = None
    km_per_substrate: Optional[dict] = None
    compound_smiles: Optional[dict] = None


@dataclass
class _Thermo:
    dgr_prime_kJmol: float = 0.0
    sigma_kJmol: float = 1.0
    keq: float = 1.0
    kcat_rev: Optional[float] = None
    irreversible: bool = False
    source: str = "equilibrator"


@dataclass
class _Rxn:
    enzyme_label: str = "E1"
    ec_number: Optional[str] = None
    reactant_labels: list = field(default_factory=list)
    product_labels: list = field(default_factory=list)
    stoichiometry: dict = field(default_factory=dict)
    kinetics: Optional[_Kin] = None
    thermo: Optional[_Thermo] = None


def _ref():
    return Reference(authors=("Test",), title="t", year="2024")


def _registry():
    reg = RuleRegistry()
    reg.register(DgrIrreversibility(
        name="dgr", description="t", reference=_ref(), reference_type="theoretical",
        params={"dgr_kjmol_cutoff": 30.0}))
    reg.register(MeasuredKcatAnchor(
        name="anchor", description="t", reference=_ref(), reference_type="experimental",
        params={"ec_kcat": {"1.3.1.9": 15.0}}))
    reg.register(TesALongChainPreference(
        name="tesa", description="t", reference=_ref(), reference_type="experimental",
        params={"ec_numbers": frozenset({"3.1.2.14"}), "n_half": 15.3, "sharpness_k": 0.64}))
    return reg


def _tesa_c16():
    return _Rxn(
        enzyme_label="TesA", ec_number="EC 3.1.2.14",
        reactant_labels=["hexadecanoyl-[ACP]", "H2O"],
        product_labels=["hexadecanoate", "holo-[ACP]"],
        stoichiometry={"hexadecanoyl-[ACP]": -1, "H2O": -1, "hexadecanoate": 1,
                       "holo-[ACP]": 1},
        kinetics=_Kin(kcat=10.0, km_per_substrate={"hexadecanoyl-[ACP]": 0.1},
                      compound_smiles={"hexadecanoyl-[ACP]": "CCCCCCCCCCCCCCCC(=O)S"}),
        thermo=_Thermo(dgr_prime_kJmol=-29.0, keq=1e5),
    )


def _fabi():
    return _Rxn(
        enzyme_label="FabI", ec_number="EC 1.3.1.9",
        reactant_labels=["S"], product_labels=["P"], stoichiometry={"S": -1, "P": 1},
        kinetics=_Kin(kcat=3.0, km_per_substrate={"S": 0.05}),
        thermo=_Thermo(dgr_prime_kJmol=-51.0, keq=1e9),
    )


class TestBaselineIsolation:
    def test_seeded_sample_goes_through_rules_and_does_not_leak(self):
        reg = _registry()
        tesa, fabi = _tesa_c16(), _fabi()
        reg.apply_all([tesa, fabi])
        tesa_kcat = tesa.kinetics.kcat
        assert tesa_kcat < 10.0
        assert fabi.kinetics.kcat == 15.0
        assert tesa.thermo.irreversible is False

        base_tesa = reg.baseline_snapshot(tesa)
        assert base_tesa["kinetics"]["kcat"] == 10.0
        sample = reg.baseline_snapshot(tesa)
        sample["kinetics"]["kcat"] = 20.0
        sample["thermo"]["dgr_prime_kJmol"] = -31.0
        sample_fabi = reg.baseline_snapshot(fabi)
        sample_fabi["kinetics"]["kcat"] = 999.0

        anchor = reg.by_name("anchor")
        saved = anchor.params["ec_kcat"]
        with reg.isolated_baselines():
            anchor.params["ec_kcat"] = {"1.3.1.9": 15.3}
            reg.seed_baseline(tesa, sample)
            reg.seed_baseline(fabi, sample_fabi)
            reg.apply_all([tesa, fabi])
            # TesA multiplier acts on the sampled value; the anchor overrides FabI's draw.
            assert tesa.kinetics.kcat == pytest.approx(2.0 * tesa_kcat)
            assert fabi.kinetics.kcat == 15.3
            # Sampled ΔG°′ crossed the cutoff: irreversibility re-decided.
            assert tesa.thermo.irreversible is True
            anchor.params["ec_kcat"] = saved

        reg.apply_all([tesa, fabi])
        assert tesa.kinetics.kcat == pytest.approx(tesa_kcat)
        assert fabi.kinetics.kcat == 15.0
        assert tesa.thermo.irreversible is False
        assert tesa.thermo.dgr_prime_kJmol == -29.0

    def test_repeated_samples_do_not_compound(self):
        reg = _registry()
        tesa = _tesa_c16()
        reg.apply_all([tesa])
        once = tesa.kinetics.kcat
        snap = reg.baseline_snapshot(tesa)
        for _ in range(3):
            with reg.isolated_baselines():
                reg.seed_baseline(tesa, snap)
                reg.apply_all([tesa])
                reg.apply_all([tesa])
            assert tesa.kinetics.kcat == pytest.approx(once)

    def test_isolated_store_restored_on_error(self):
        reg = _registry()
        tesa = _tesa_c16()
        reg.apply_all([tesa])
        with pytest.raises(RuntimeError):
            with reg.isolated_baselines():
                raise RuntimeError("boom")
        assert reg.baseline_snapshot(tesa)["kinetics"]["kcat"] == 10.0

    def test_unseen_reaction_has_no_snapshot(self):
        assert RuleRegistry().baseline_snapshot(_tesa_c16()) is None


# ---------------------------------------------------------------------------
# Parameter sampling helpers (project script)
# ---------------------------------------------------------------------------

_PROJECT = (Path(__file__).resolve().parents[1]
            / "projects/fattyAcidSynthesis/fattyAcidSynthesis_ecoli")


@pytest.fixture(scope="module")
def sampling():
    sys.path.insert(0, str(_PROJECT))
    try:
        import uncertainty as mod
    except FileNotFoundError as exc:  # S2A experiment CSV lives in local-only knowledge/
        pytest.skip(f"FAS sampling script data unavailable: {exc}")
    return mod


def _inventory(sampling):
    D = sampling.Dim
    kcat = [D("kcat", 0, None, 5.0, 0.7, "FabB", "r0", rxns=(0,)),
            D("kcat", 1, None, 2.0, 0.3, "FabF", "r1", rxns=(1, 2))]
    km = [D("km", 0, "S", 0.04, 0.6, "FabB", "r0", "reactant", rxns=(0,)),
          D("km", 1, "S", 0.01, 0.5, "FabF", "r1", "reactant", rxns=(1, 2))]
    dg = [D("dg", 0, None, -10.0, 2.0, "FabB", "r0", rxns=(0,)),
          D("dg", 1, None, -29.0, 1.5, "FabF", "r1", rxns=(1,))]
    chol = sampling._psd_factor(np.array([[4.0, 1.0], [1.0, 2.25]]))
    return sampling.Inventory(kcat, km, dg, chol, 0.021, [])


class TestSampling:
    def test_tie_merges_batch_noise_repeats_only(self, sampling):
        D = sampling.Dim
        dims = [D("km", 0, "NADPH", 0.1198231657, 0.536, "FabG", "a", "reactant"),
                D("km", 1, "nadph", 0.1198231764, 0.536, "FabG", "b", "reactant"),
                D("km", 2, "NADPH", 0.1198231700, 0.536, "FabI", "c", "reactant"),
                D("km", 3, "holo-[ACP]", 0.039, 0.3, "FabD", "d", "product"),
                D("km", 4, "holo-[ACP]", 0.076, 0.3, "FabD", "e", "product"),
                D("km", 5, None, 0.5, 0.2, "FabB", "f", "scalar")]
        tied = sampling.tie_identical_km(dims)
        assert len(tied) == 5
        assert tied[0].km_targets() == ((0, "NADPH"), (1, "nadph"))
        assert [d.km_targets() for d in tied[1:]] == [
            ((2, "NADPH"),), ((3, "holo-[ACP]"),), ((4, "holo-[ACP]"),), ((5, None),)]

    def test_perturbed_snapshots_write_tied_km_under_each_label(self, sampling):
        D = sampling.Dim
        km = sampling.tie_identical_km([D("km", 0, "NADPH", 0.12, 0.5, "FabG", "a"),
                                        D("km", 1, "nadph", 0.12, 0.5, "FabG", "b")])
        inv = sampling.Inventory([], km, [], np.zeros((0, 0)), 0.021, [])
        base = [{"kinetics": {"kcat": 1.0, "km": 0.12, "km_per_substrate": {"NADPH": 0.12}}},
                {"kinetics": {"kcat": 1.0, "km": 0.12, "km_per_substrate": {"nadph": 0.12}}}]
        vals = {"kcat": np.array([]), "km": np.array([0.5]), "dg": np.array([])}
        snaps = sampling.perturbed_snapshots(base, inv, vals, 298.15)
        assert snaps[0]["kinetics"]["km_per_substrate"]["NADPH"] == vals["km"][0]
        assert snaps[1]["kinetics"]["km_per_substrate"]["nadph"] == vals["km"][0]

    def test_perturbed_snapshots_write_shared_groups_and_keep_baseline(self, sampling):
        inv = _inventory(sampling)
        base = [
            {"kinetics": {"kcat": 5.0, "km": 0.04, "km_per_substrate": {"S": 0.04}},
             "thermo": {"dgr_prime_kJmol": -10.0, "keq": 50.0}},
            {"kinetics": {"kcat": 2.0, "km": 0.01, "km_per_substrate": {"S": 0.01}},
             "thermo": {"dgr_prime_kJmol": -29.0, "keq": 1e5}},
            {"kinetics": {"kcat": 2.0, "km": 0.01, "km_per_substrate": {"S": 0.01}},
             "thermo": {"dgr_prime_kJmol": -29.0, "keq": 1e5}},
        ]
        zero = sampling.perturbed_snapshots(base, inv, {
            "kcat": np.array([5.0, 2.0]),
            "km": np.array([0.04, 0.01]),
            "dg": np.array([-10.0, -29.0]),
        }, 298.15)
        assert zero == base
        vals = {
            "kcat": np.array([9.0, 4.0]),
            "km": np.array([0.2, 0.08]),
            "dg": np.array([-12.0, -20.0]),
        }
        snaps = sampling.perturbed_snapshots(base, inv, vals, 298.15)
        assert base[1]["kinetics"]["kcat"] == 2.0
        assert snaps[1]["kinetics"]["kcat"] == snaps[2]["kinetics"]["kcat"] == vals["kcat"][1]
        assert snaps[1]["kinetics"]["km_per_substrate"]["S"] == vals["km"][1]
        assert snaps[2]["kinetics"]["km_per_substrate"]["S"] == vals["km"][1]
        assert snaps[1]["kinetics"]["km"] == vals["km"][1]
        assert snaps[0]["thermo"]["dgr_prime_kJmol"] == vals["dg"][0]
        assert snaps[0]["thermo"]["keq"] == pytest.approx(
            math.exp(-vals["dg"][0] * 1000 / (8.314462618 * 298.15)), rel=1e-6)


# ---------------------------------------------------------------------------
# Published ranking and band statistics
# ---------------------------------------------------------------------------

_PROJECT = os.path.join(
    os.path.dirname(__file__),
    "..",
    "projects",
    "fattyAcidSynthesis",
    "fattyAcidSynthesis_ecoli",
)
sys.path.insert(0, os.path.abspath(_PROJECT))

import uncertainty as unc  # noqa: E402


def test_baseline_values_are_the_unperturbed_means():
    D = unc.Dim
    inv = unc.Inventory(
        [D("kcat", 0, None, 5.0, 0.7, "FabB", "r0", rxns=(0,))],
        [D("km", 0, "S", 0.04, 0.6, "FabB", "r0", "reactant", rxns=(0,))],
        [D("dg", 0, None, -10.0, 2.0, "FabB", "r0", rxns=(0,))],
        np.zeros((0, 0)),
        0.021,
        [],
    )
    base = unc._baseline_values(inv)
    np.testing.assert_array_equal(base["kcat"], [5.0])
    np.testing.assert_array_equal(base["km"], [0.04])
    np.testing.assert_array_equal(base["dg"], [-10.0])
    assert float(base["fabi_kcat"]) == unc.FABI_KCAT_MEAN
    assert float(base["fabh_ki_factor"]) == 1.0
    assert float(base["fabi_ki_factor"]) == 1.0


def test_ln_sd_converts_log10():
    assert unc.ln_sd(1.0) == pytest.approx(math.log(10.0))
    assert unc.ln_sd(0.5) == pytest.approx(0.5 * math.log(10.0))


def test_coeff_of_sums_slope_times_sigma():
    assert unc.coeff_of([(2.0, 0.5), (-1.0, 1.0)]) == pytest.approx(0.0)
    assert unc.coeff_of([(0.2, unc.ln_sd(1.0))]) == pytest.approx(0.2 * math.log(10.0))


def test_zero_coeff_gives_zero_variance_share():
    assert unc.coeff_of([]) == 0.0
    groups = [
        unc.Group("a", "kcat", 1.0, 1.0),
        unc.Group("b", "km", 0.0, 0.0),
    ]
    chosen, ridge = unc.select_groups(groups)
    assert [g.name for g in chosen] == ["a"]
    assert ridge is False


def test_dg_variance_is_quadratic_form():
    s = np.array([1.0, 0.0])
    cov = np.array([[4.0, 0.0], [0.0, 9.0]])
    assert unc.dg_variance(s, cov) == pytest.approx(4.0)


def test_select_top_ode_skips_dg_and_fabi_kcat():
    groups = [
        unc.Group("FabI Km", "km", 1.0, 1.0, sd_log10=0.5),
        unc.Group("dG joint", "dg", 0.5, 0.25, sd_log10=0.0),
        unc.Group("FabG kcat", "kcat", 0.4, 0.16, sd_log10=0.6),
        unc.Group("FabI kcat", "kcat", 0.01, 1e-4, sd_log10=0.0, is_fabi_kcat=True),
        unc.Group("FabF kcat", "kcat", 0.3, 0.09, sd_log10=0.4),
    ]
    chosen = unc.select_top_ode(groups, n=3)
    assert [g.name for g in chosen] == ["FabI Km", "FabG kcat", "FabF kcat"]


def test_design_loguniform_in_unit_interval():
    u = unc.design_loguniform(32, 3, seed=1)
    assert u.shape == (32, 3)
    assert np.all(u >= -1.0) and np.all(u <= 1.0)


def test_bootstrap_percentile_ci_covers_point_estimate():
    rng = np.random.default_rng(0)
    x = rng.normal(10.0, 2.0, size=200)
    q = 5.0
    point = float(np.percentile(x, q))
    lo, hi = unc.bootstrap_percentile_ci(x, q, n_boot=500, seed=1)
    assert lo <= point <= hi


def test_stability_ok_relative_and_absolute():
    assert unc.stability_ok(10.0, 10.5) is True
    assert unc.stability_ok(10.0, 13.0) is False
    assert unc.stability_ok(0.2, 0.5) is True  # both < 1 µM, |Δ| < 1
    assert unc.stability_ok(0.1, 2.0) is False


def test_stability_at_12min_on_identical_halves():
    x = np.linspace(5.0, 50.0, 128)
    stab = unc.stability_at_12min(x, split=64)
    assert stab["n_first"] == 64
    assert stab["n_all"] == 128
    # The first half is systematically lower, so mean and percentiles should fail.
    assert stab["pass"] is False


def test_band_label_is_the_sd_box():
    assert unc.BAND_LABEL == "±1 CatPred SD box, top 4 groups"
    assert "5–95% of CatPred" not in unc.BAND_LABEL
