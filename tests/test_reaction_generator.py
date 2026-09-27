"""Tests for ReactionGenerator and GeneratedReaction."""

from types import SimpleNamespace
from unittest.mock import MagicMock

from bees.reaction_generator import GeneratedReaction, ReactionGenerator
from bees.reaction_template import ECClass, ReactionTemplate


def _make_template(ec_class=ECClass.TRANSFERASE):
    return ReactionTemplate(
        template_type="phosphorylation",
        ec_class=ec_class,
    )


class TestGeneratedReaction:
    def test_repr(self):
        rxn = GeneratedReaction(
            enzyme_label="Hexokinase",
            substrate_label="Glucose",
            ec_number="EC 2.7.1.1",
            template=_make_template(),
            kinetics=None,
            reactant_labels=["Glucose", "ATP"],
            product_labels=["Glucose-6P", "ADP"],
            stoichiometry={"Glucose": -1, "ATP": -1, "Glucose-6P": 1, "ADP": 1},
        )
        s = repr(rxn)
        assert "Glucose" in s
        assert "ATP" in s
        assert "Glucose-6P" in s
        assert "ADP" in s

    def test_with_rate_law(self):
        kin = MagicMock()
        rxn = GeneratedReaction(
            enzyme_label="E",
            substrate_label="S",
            ec_number="EC 1.1.1.1",
            template=_make_template(),
            kinetics=kin,
            reactant_labels=["S"],
            product_labels=["P"],
            stoichiometry={"S": -1, "P": 1},
            rate_law="Michaelis-Menten",
        )
        assert rxn.rate_law == "Michaelis-Menten"


class TestReactionGenerator:
    def test_init(self):
        bees_obj = MagicMock()
        logger = MagicMock()
        gen = ReactionGenerator(bees_obj, logger, "/tmp/out")
        assert gen.bees_object is bees_obj
        assert gen.logger is logger
        assert gen.output_directory == "/tmp/out"
        assert gen.reactions == []
        assert gen.kinetic_db is None
        assert gen.kinetics_estimator is None

    def test_check_heavy_atom_balance_balanced(self):
        """Stoichiometry is balanced if total heavy atoms on both sides match."""
        gen = ReactionGenerator(MagicMock(), MagicMock(), "/tmp/out")

        smiles_map = {
            "A": "CCO",  # 3 heavy atoms
            "B": "CCO",  # 3 heavy atoms
        }

        gen._resolve_smiles_for_compound = MagicMock(  # type: ignore[method-assign]
            side_effect=lambda compound_label, kinetic_data, substrate_label: smiles_map.get(compound_label)
        )

        stoich = {"A": -1, "B": 1}
        assert gen._check_heavy_atom_balance(stoich, kinetic_data=None, substrate_label="A") is True

    def test_check_heavy_atom_balance_unbalanced(self):
        gen = ReactionGenerator(MagicMock(), MagicMock(), "/tmp/out")

        smiles_map = {
            "A": "CCO",  # 3
            "B": "CC",  # 2
        }
        gen._resolve_smiles_for_compound = MagicMock(  # type: ignore[method-assign]
            side_effect=lambda compound_label, kinetic_data, substrate_label: smiles_map.get(compound_label)
        )

        stoich = {"A": -1, "B": 1}
        assert gen._check_heavy_atom_balance(stoich, kinetic_data=None, substrate_label="A") is False

    def test_check_heavy_atom_balance_missing_smiles_returns_none(self):
        gen = ReactionGenerator(MagicMock(), MagicMock(), "/tmp/out")
        gen._resolve_smiles_for_compound = MagicMock(return_value=None)  # type: ignore[method-assign]

        stoich = {"A": -1, "B": 1}
        assert gen._check_heavy_atom_balance(stoich, kinetic_data=None, substrate_label="A") is None

    def test_check_heavy_atom_balance_invalid_coeff_returns_none(self):
        gen = ReactionGenerator(MagicMock(), MagicMock(), "/tmp/out")
        gen._resolve_smiles_for_compound = MagicMock(return_value="CC")  # type: ignore[method-assign]

        stoich = {"A": "not_a_number", "B": 1}
        assert gen._check_heavy_atom_balance(stoich, kinetic_data=None, substrate_label="A") is None


def _make_generator_with_estimator(captured: dict) -> ReactionGenerator:
    bees = MagicMock()
    bees.enzymes = []
    bees.species = []
    bees.settings = MagicMock()
    logger = MagicMock()
    gen = ReactionGenerator(bees, logger, output_directory="/tmp")

    def _estimate(**kwargs):
        captured["reactant_smiles"] = dict(kwargs.get("reactant_smiles") or {})
        est = MagicMock()
        est.km_per_substrate = {lab: 0.05 for lab in captured["reactant_smiles"]}
        est.compound_smiles = dict(captured["reactant_smiles"])
        return est

    gen.kinetics_estimator = MagicMock()
    gen.kinetics_estimator.estimate = _estimate
    return gen


def test_regulatory_cofactor_product_queried_buffered_skipped():
    """CoA gets a reverse Km query; H2O and CO2 do not."""
    captured: dict = {}
    gen = _make_generator_with_estimator(captured)

    rxn = SimpleNamespace(
        enzyme_label="FabD",
        substrate_label="Malonyl-CoA",
        product_labels=["malonyl-[ACP]", "Coenzyme A", "H2O", "Carbon dioxide"],
        stoichiometry={
            "holo-[ACP]": -1,
            "Malonyl-CoA": -1,
            "malonyl-[ACP]": 1,
            "Coenzyme A": 1,
            "H2O": 1,
            "Carbon dioxide": 1,
        },
        ec_number="EC 2.3.1.39",
    )
    smiles_map = {
        "malonyl-[ACP]": "CC(=O)CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)O",
        "Coenzyme A": (
            "CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)"
            "N1C=NC2=C1N=CN=C2N)C(O)C(=O)NCCC(=O)NCCS"
        ),
        "H2O": "O",
        "Carbon dioxide": "O=C=O",
    }

    product_kms, _ = gen._lookup_product_kms_via_reverse_query(
        reaction=rxn,
        stoich=rxn.stoichiometry,
        ec_numbers_to_try=["EC 2.3.1.39"],
        temp_range=None,
        ph_range=None,
        provided_species_labels_lc=set(),
        enzyme_sequence="MKT",
        smiles_map=smiles_map,
        substitutor=None,
    )

    queried = set(captured.get("reactant_smiles", {}))
    assert "Coenzyme A" in queried
    assert "malonyl-[ACP]" in queried
    assert "H2O" not in queried
    assert "Carbon dioxide" not in queried
    assert "Coenzyme A" in product_kms
    assert "H2O" not in product_kms


def test_nadp_product_queried_for_fabg_like():
    """NADP (regulatory) is queried; H+ (buffered) is not."""
    captured: dict = {}
    gen = _make_generator_with_estimator(captured)
    rxn = SimpleNamespace(
        enzyme_label="FabG",
        substrate_label="3-oxobutanoyl-[ACP]",
        product_labels=["(3R)-hydroxybutanoyl-[ACP]", "NADP"],
        stoichiometry={
            "3-oxobutanoyl-[ACP]": -1,
            "NADPH": -1,
            "H+": -1,
            "(3R)-hydroxybutanoyl-[ACP]": 1,
            "NADP": 1,
        },
        ec_number="EC 1.1.1.100",
    )
    smiles_map = {
        "(3R)-hydroxybutanoyl-[ACP]": "CC(O)CC(=O)S",
        "NADP": (
            "NC(=O)C1=C[N+](=CC=C1)C1OC(COP(=O)(O)OP(=O)(O)OCC2OC("
            "N3C=NC4=C(N)N=CN=C43)C(OP(=O)(O)O)C2O)C(O)C1O"
        ),
        "H+": "[H+]",
    }
    product_kms, _ = gen._lookup_product_kms_via_reverse_query(
        reaction=rxn,
        stoich=rxn.stoichiometry,
        ec_numbers_to_try=["EC 1.1.1.100"],
        temp_range=None,
        ph_range=None,
        provided_species_labels_lc=set(),
        enzyme_sequence="MNF",
        smiles_map=smiles_map,
        substitutor=None,
    )
    assert "NADP" in captured["reactant_smiles"]
    assert "H+" not in captured["reactant_smiles"]
    assert product_kms.get("NADP", 0) > 0
