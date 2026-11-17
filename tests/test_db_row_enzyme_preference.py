"""When two database rows share an EC and a reaction, keep the named enzyme's row."""

from bees.logger import Logger
from db.reaction_database import KineticData, ReactionDatabase, _ROW_ID_KEYS

_COMPLETE = "EC 1.2.3.4"
_SHORT = {
    "short reactant": -1,
    "short product": 1,
    "Water": 1,
}
_LONG = {
    "long reactant": -1,
    "long product": 1,
    "Water": 1,
}


def _row(name, stoich, kcat=None):
    return KineticData(
        ec_number=_COMPLETE,
        enzyme_name=name,
        reaction_string=name,
        stoichiometry=dict(stoich),
        kcat=kcat,
        meta={"entry_name": name},
    )


def _db(tmp_path, rows):
    db = ReactionDatabase(Logger(str(tmp_path), None, 0.0))
    db._ec_index[_COMPLETE] = rows
    return db


def _only(db, substrate, enzyme_label=None):
    found = db.query_by_enzyme_substrate(
        ec_number=_COMPLETE,
        substrate_label=substrate,
        return_all=True,
        enzyme_label=enzyme_label,
    )
    assert len(found) == 1
    return found[0]


def test_named_enzyme_wins_over_higher_priority_row(tmp_path):
    # The first row has a kcat, so it scores higher when no name is given.
    first = _row("EnzymeA", _SHORT, kcat=1.0)
    second = _row("EnzymeB", _SHORT)
    db = _db(tmp_path, [first, second])
    assert _only(db, "short reactant", "EnzymeB") is second
    assert _only(db, "short reactant", "EnzymeA") is first


def test_name_match_ignores_case_and_trailing_space(tmp_path):
    first = _row("EnzymeA", _SHORT, kcat=1.0)
    second = _row("EnzymeB", _SHORT)
    db = _db(tmp_path, [first, second])
    assert _only(db, "short reactant", "enzymeb ") is second


def test_unknown_label_matches_no_label(tmp_path):
    first = _row("EnzymeA", _SHORT, kcat=1.0)
    second = _row("EnzymeB", _SHORT)
    db = _db(tmp_path, [first, second])
    unnamed = _only(db, "short reactant")
    custom = _only(db, "short reactant", "MyEnzyme")
    assert custom is unnamed is first


def test_falls_back_when_named_enzyme_has_no_row(tmp_path):
    first = _row("EnzymeA", _LONG, kcat=1.0)
    db = _db(tmp_path, [first])
    assert _only(db, "long reactant", "EnzymeB") is first


_INCOMPLETE = "EC 1.2.3.-"
_SHARED = {"shared reactant": -1, "Water": -1, "shared product": 1}
_OWN = {"own reactant": -1, "Water": -1, "own product": 1}
_ACCESSION_KEY = _ROW_ID_KEYS[1]


def test_incomplete_ec_does_not_borrow_another_enzyme(tmp_path):
    """An incomplete EC must not hand one enzyme a row that belongs to another."""
    other = KineticData(
        ec_number=_INCOMPLETE,
        enzyme_name="Other enzyme",
        reaction_string="shared reactant + Water = shared product",
        stoichiometry=dict(_SHARED),
        meta={"entry_name": "OTHER_ENTRY", _ACCESSION_KEY: "ACC1"},
    )
    owner = KineticData(
        ec_number=_INCOMPLETE,
        enzyme_name="Owner enzyme",
        reaction_string="own reactant + Water",
        stoichiometry=dict(_OWN),
        meta={"entry_name": "OWNER_ENTRY"},
    )
    db = ReactionDatabase(Logger(str(tmp_path), None, 0.0))
    db._ec_index[_INCOMPLETE] = [other, owner]

    borrowed = db.query_by_enzyme_substrate(
        ec_number=_INCOMPLETE,
        substrate_label="shared reactant",
        return_all=True,
        enzyme_label="OWNER_ENTRY",
    )
    assert borrowed == []

    own = db.query_by_enzyme_substrate(
        ec_number=_INCOMPLETE,
        substrate_label="own reactant",
        return_all=True,
        enzyme_label="OWNER_ENTRY",
    )
    assert own == [owner]


def test_incomplete_ec_keeps_own_copy_when_the_equation_is_shared(tmp_path):
    """An enzyme keeps its own row when another enzyme has the same equation."""
    def shared_row(entry):
        return KineticData(
            ec_number=_INCOMPLETE,
            enzyme_name=entry,
            reaction_string="shared reactant + Water = shared product",
            stoichiometry=dict(_SHARED),
            meta={"entry_name": entry},
        )
    other = shared_row("OTHER_ENTRY")
    owner = shared_row("OWNER_ENTRY")
    db = ReactionDatabase(Logger(str(tmp_path), None, 0.0))
    db._ec_index[_INCOMPLETE] = [other, owner]
    found = db.query_by_enzyme_substrate(
        ec_number=_INCOMPLETE,
        substrate_label="shared reactant",
        return_all=True,
        enzyme_label="OWNER_ENTRY",
    )
    assert found == [owner]
