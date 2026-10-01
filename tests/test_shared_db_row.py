"""Isozymes that share an EC must not share one database row object.

The fake estimator returns the same EstimatedKinetics instance and the same
nested Km dicts on every matching call, the way the CatPred memo does.
"""

from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from bees.kinetics_estimator import EstimatedKinetics
from bees.reaction_generator import ReactionGenerator
from bees.rules import load_rules
from bees.rules.base import RULES
from bees.rules.physics_rules import (
    _DG_PER_CH2_KJMOL,
    _HYDROPHOBIC_REFERENCE_CHAIN,
    _HYDROPHOBIC_TEMPERATURE_K,
    _R_KJ,
)
from db.reaction_database import KineticData

_SUB = "(3R)-hydroxytetradecanoyl-[ACP]"
_PROD = "(2E)-tetradecenoyl-[ACP]"
_SMI = {
    _SUB: "CCCCCCCCCCC[C@@H](O)CC(=O)S",
    _PROD: "CCCCCCCCCCC=CC(=O)S",
}
_EC = "EC 4.2.1.59"


class _SharedMemoEstimator:
    """One cached EstimatedKinetics object, and its dicts, per matching call."""

    def __init__(self):
        self.memo = {}

    def estimate(
        self,
        *,
        enzyme_sequence,
        reactant_smiles,
        inhibitor_smiles=None,
        ec_number=None,
    ):
        key = (enzyme_sequence, tuple(sorted(reactant_smiles)))
        hit = self.memo.get(key)
        if hit is None:
            km_value = 0.040 if str(enzyme_sequence).endswith("A") else 0.050
            km = {lab: km_value for lab in reactant_smiles}
            sd = {lab: 0.01 for lab in reactant_smiles}
            hit = EstimatedKinetics(
                km=km_value,
                km_per_substrate=km,
                km_sd_per_substrate=sd,
                kcat=11.5 if str(enzyme_sequence).endswith("A") else 7.36,
                kcat_sd=1.0,
                source="fake-cache",
            )
            self.memo[key] = hit
        return hit


class _OneRowDB:
    """Every matching query returns the same KineticData object."""

    def __init__(self, row):
        self.row = row

    def query_by_enzyme_substrate(self, **kwargs):
        if kwargs.get("substrate_label") != _SUB:
            return []
        return [self.row]


def _row():
    return KineticData(
        ec_number=_EC,
        enzyme_name="shared",
        reaction_string=f"{_SUB} = {_PROD}",
        stoichiometry={_SUB: -1, _PROD: 1},
        compound_smiles=dict(_SMI),
        source="database",
        kcat=1.0,
        km_per_substrate={_SUB: 9.0},
    )


def _enzyme(label, sequence):
    return SimpleNamespace(
        label=label,
        ecnumber=_EC,
        amino_acid_sequence=sequence,
        concentration=0.001,
        reactive=True,
    )


def _substrate():
    return SimpleNamespace(label=_SUB, smiles=_SMI[_SUB], reactive=True)


def _generator(enzymes, estimator):
    species = [
        SimpleNamespace(label=_SUB, smiles=_SMI[_SUB], reactive=True),
        SimpleNamespace(label=_PROD, smiles=_SMI[_PROD], reactive=True),
    ]
    bees = SimpleNamespace(
        species=species,
        enzymes=enzymes,
        environment=SimpleNamespace(temperature=298.15, pH=7.4),
        settings=SimpleNamespace(smiles_mode="auto", estimate_kinetics=True),
    )
    gen = ReactionGenerator(bees, MagicMock(), output_directory="/tmp")
    gen.kinetic_db = _OneRowDB(_row())
    gen.kinetics_estimator = estimator
    return gen


def _generate_pair(with_thermo=False):
    enzymes = [_enzyme("FabA", "SEQA"), _enzyme("FabZ", "SEQZ")]
    estimator = _SharedMemoEstimator()
    gen = _generator(enzymes, estimator)
    available = {_SUB.lower(), _PROD.lower()}
    thermo = None
    if with_thermo:
        thermo = MagicMock()
        thermo.compute_keq.return_value = SimpleNamespace(
            dgr_prime_kJmol=-10.0,
            sigma_kJmol=1.0,
            keq=50.0,
            irreversible=False,
            source="fake",
        )
    made = []
    for enzyme in enzymes:
        made.extend(
            gen._generate_reactions(
                enzyme,
                _substrate(),
                available_species_labels_lc=available,
                provided_species_labels_lc=available,
                thermo_engine=thermo,
                global_smiles_map=dict(_SMI),
            )
        )
    by_enzyme = {rxn.enzyme_label: rxn for rxn in made}
    return estimator, by_enzyme


def _factor(n=14):
    import math
    return math.exp(
        _DG_PER_CH2_KJMOL * (n - _HYDROPHOBIC_REFERENCE_CHAIN)
        / (_R_KJ * _HYDROPHOBIC_TEMPERATURE_K)
    )


def test_isozymes_keep_separate_kinetics():
    _estimator, by_enzyme = _generate_pair()
    fab_a = by_enzyme["FabA"]
    fab_z = by_enzyme["FabZ"]
    assert fab_a.kinetics is not fab_z.kinetics
    assert fab_a.kinetics.kcat == pytest.approx(11.5)
    assert fab_z.kinetics.kcat == pytest.approx(7.36)
    assert fab_a.kinetics.km_per_substrate[_SUB] == pytest.approx(0.040)
    assert fab_z.kinetics.km_per_substrate[_SUB] == pytest.approx(0.050)


def test_chain_length_km_scales_once_per_reaction():
    _estimator, by_enzyme = _generate_pair()
    load_rules()
    RULES._baselines.clear()
    RULES.apply_all([by_enzyme["FabA"], by_enzyme["FabZ"]])
    once = _factor()
    assert by_enzyme["FabA"].kinetics.km_per_substrate[_SUB] == pytest.approx(0.040 * once)
    assert by_enzyme["FabZ"].kinetics.km_per_substrate[_SUB] == pytest.approx(0.050 * once)


def test_generation_does_not_mutate_estimator_memo():
    estimator, by_enzyme = _generate_pair(with_thermo=True)
    assert "FabA" in by_enzyme and "FabZ" in by_enzyme
    forward = estimator.memo[("SEQA", (_SUB,))]
    assert set(forward.km_per_substrate) == {_SUB}
    assert forward.km_per_substrate[_SUB] == pytest.approx(0.040)
    reverse = estimator.memo[("SEQA", (_PROD,))]
    assert reverse is not forward
    assert set(reverse.km_per_substrate) == {_PROD}
