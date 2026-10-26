"""Acyl-chain-length detector validation corpus.

Gates the HydrophobicChainLengthKm rule: the detector must return the correct
carbon count for every acyl species the rule will encounter, and None for
non-acyl species, BEFORE the rule is allowed to act on Km. Kept in its own file
so the detector can be validated independently of the rule.

Chain length comes from structure only. The critical regression is
`test_acyl_coa_returns_acyl_length_not_total`: the detector must return the ACYL
chain length, never the total molecule carbon count (a full acyl-CoA carries
~20 extra carbons in the pantetheine tail).
"""

from __future__ import annotations

import pytest

from bees.cofactors import COA_TAIL, PPANT_HANDLE
from bees.rules.physics_rules import detect_acyl_chain_length

_SNAC = "SCCNC(C)=O"
_CARRIERS = {"free acid": "O", "acyl-carrier handle": PPANT_HANDLE, "CoA": COA_TAIL, "SNAC": _SNAC}


def _saturated(n: int) -> str:
    return "C" * (n - 1) + "C(=O)"


# ---------------------------------------------------------------------------
# Saturated acyl, every carrier
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("carrier", sorted(_CARRIERS))
@pytest.mark.parametrize("n", [2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
def test_saturated_chain_lengths(n, carrier):
    assert detect_acyl_chain_length("x", _saturated(n) + _CARRIERS[carrier]) == n


# ---------------------------------------------------------------------------
# Modifiers do not change the carbon count
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("acyl,expected", [
    ("CC(=O)CC(=O)", 4),                       # 3-oxobutanoyl
    ("C[C@@H](O)CC(=O)", 4),                   # (3R)-hydroxybutanoyl
    ("C/C=C/C(=O)", 4),                        # (2E)-butenoyl
    ("CCCC(=O)CC(=O)", 6),                     # 3-oxohexanoyl
    ("CCCCCCCCCCC[C@@H](O)CC(=O)", 14),        # (3R)-hydroxytetradecanoyl
    ("CCCCCC/C=C\\CCCCCCCC(=O)", 16),          # (9Z)-hexadecenoyl
    ("OC(=O)CC(=O)", 3),                       # malonyl
])
def test_modifier_chain_lengths(acyl, expected):
    assert detect_acyl_chain_length("x", acyl + PPANT_HANDLE) == expected


# ---------------------------------------------------------------------------
# Non-acyl / cofactors / bare carriers / no structure → None (rule must skip)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("label,smiles", [
    ("NADPH", None),
    ("H2O", "O"),
    ("CO2", "O=C=O"),
    ("carbon dioxide", "O=C=O"),
    ("holo-carrier", PPANT_HANDLE),
    ("CoA", COA_TAIL),
    ("acyl without structure", None),
    ("", None),
])
def test_non_acyl_returns_none(label, smiles):
    assert detect_acyl_chain_length(label, smiles) is None


# ---------------------------------------------------------------------------
# Acyl-moiety extraction, NOT total carbon count
# ---------------------------------------------------------------------------

class TestAcylMoietyExtraction:
    def test_acyl_coa_returns_acyl_length_not_total(self):
        coa_smiles = _saturated(4) + COA_TAIL
        assert coa_smiles.count("C") > 15  # the trap: total carbons are many
        n = detect_acyl_chain_length("x", coa_smiles)
        assert n == 4, f"expected acyl length 4, got {n} (must not count CoA tail)"

    def test_unparseable_smiles_returns_none(self):
        assert detect_acyl_chain_length("x", "not_a_smiles") is None

    def test_no_anchor_returns_none(self):
        # A plain alkane has no thioester/carboxyl carbonyl → no acyl anchor.
        assert detect_acyl_chain_length("x", "CCCCCC") is None

    def test_label_is_not_parsed_for_chain_length(self):
        # A misleading name does not override the structure.
        assert detect_acyl_chain_length("butanoyl-[ACP]", _saturated(12) + "O") == 12


# ---------------------------------------------------------------------------
# detect_acyl_carrier — what the acyl carbonyl is bound to, from structure only
# ---------------------------------------------------------------------------

class TestAcylCarrier:
    def test_free_acid_same_key_charged_or_not(self):
        from bees.rules.helpers import detect_acyl_carrier
        assert detect_acyl_carrier("CCCC(=O)O") == detect_acyl_carrier("CCCC(=O)[O-]") == "O"

    def test_same_carrier_across_chain_lengths(self):
        from bees.rules.helpers import detect_acyl_carrier
        assert detect_acyl_carrier("CCCC(=O)" + COA_TAIL) == detect_acyl_carrier(
            "CCCCCCCCCCCCCCCC(=O)" + COA_TAIL)

    def test_different_carriers_differ(self):
        from bees.rules.helpers import detect_acyl_carrier
        keys = {detect_acyl_carrier(s) for s in
                ("CCCC(=O)O", "CCCC(=O)" + _SNAC, "CCCC(=O)" + COA_TAIL)}
        assert len(keys) == 3

    def test_malonyl_thioester_anchors_on_carrier_not_free_carboxyl(self):
        from bees.rules.helpers import detect_acyl_carrier
        assert detect_acyl_carrier("OC(=O)CC(=O)" + _SNAC) == detect_acyl_carrier("CCCC(=O)" + _SNAC)
        assert detect_acyl_chain_length("x", "OC(=O)CC(=O)" + _SNAC) == 3

    def test_no_structure_or_no_acyl_group(self):
        from bees.rules.helpers import detect_acyl_carrier
        assert detect_acyl_carrier(None) is None
        assert detect_acyl_carrier("CCCCCC") is None
        assert detect_acyl_carrier("not_a_smiles") is None
