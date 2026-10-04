"""When two database rows share an EC and a reaction, keep the named enzyme's row."""

from bees.logger import Logger
from db.reaction_database import KineticData, ReactionDatabase

_EC = "EC 4.2.1.59"
_C4 = {
    "(3R)-hydroxybutanoyl-[ACP]": -1,
    "(2E)-butenoyl-[ACP]": 1,
    "Water": 1,
}
_C18 = {
    "(3R)-hydroxyoctadecanoyl-[ACP]": -1,
    "(2E)-octadecenoyl-[ACP]": 1,
    "Water": 1,
}


def _row(name, stoich, kcat=None):
    return KineticData(
        ec_number=_EC,
        enzyme_name=name,
        reaction_string=name,
        stoichiometry=dict(stoich),
        kcat=kcat,
        meta={"uniprot_entry": name},
    )


def _db(tmp_path, rows):
    db = ReactionDatabase(Logger(str(tmp_path), None, 0.0))
    db._ec_index[_EC] = rows
    return db


def _only(db, substrate, enzyme_label=None):
    found = db.query_by_enzyme_substrate(
        ec_number=_EC,
        substrate_label=substrate,
        return_all=True,
        enzyme_label=enzyme_label,
    )
    assert len(found) == 1
    return found[0]


def test_named_enzyme_wins_over_higher_priority_row(tmp_path):
    # FabA has a kcat, so it scores higher and wins when no name is given.
    fab_a = _row("FabA", _C4, kcat=1.0)
    fab_z = _row("FabZ", _C4)
    db = _db(tmp_path, [fab_a, fab_z])
    assert _only(db, "(3R)-hydroxybutanoyl-[ACP]", "FabZ") is fab_z
    assert _only(db, "(3R)-hydroxybutanoyl-[ACP]", "FabA") is fab_a


def test_name_match_ignores_case_and_trailing_space(tmp_path):
    fab_a = _row("FabA", _C4, kcat=1.0)
    fab_z = _row("FabZ", _C4)
    db = _db(tmp_path, [fab_a, fab_z])
    assert _only(db, "(3R)-hydroxybutanoyl-[ACP]", "fabz ") is fab_z


def test_unknown_label_matches_no_label(tmp_path):
    fab_a = _row("FabA", _C4, kcat=1.0)
    fab_z = _row("FabZ", _C4)
    db = _db(tmp_path, [fab_a, fab_z])
    unnamed = _only(db, "(3R)-hydroxybutanoyl-[ACP]")
    custom = _only(db, "(3R)-hydroxybutanoyl-[ACP]", "MyDehydratase")
    assert custom is unnamed is fab_a


def test_falls_back_when_named_enzyme_has_no_row(tmp_path):
    fab_a = _row("FabA", _C18, kcat=1.0)
    db = _db(tmp_path, [fab_a])
    assert _only(db, "(3R)-hydroxyoctadecanoyl-[ACP]", "FabZ") is fab_a
