"""
Every tracked project input must validate against the schema and name only known calibrations.

To run: pytest -v tests/test_project_inputs.py
"""

import glob
import os

import pytest
import yaml

from bees.common import BEES_PATH
from bees.rules import load_rules
from bees.schema import InputBase

_INPUTS = sorted(
    p for p in glob.glob(os.path.join(BEES_PATH, "projects", "**", "input*.yml"), recursive=True)
    if f"{os.sep}output{os.sep}" not in p
)


def test_inputs_found():
    assert _INPUTS


@pytest.mark.parametrize("path", _INPUTS, ids=lambda p: os.path.relpath(p, BEES_PATH))
def test_project_input_validates(path):
    with open(path) as fh:
        data = yaml.safe_load(fh)
    parsed = InputBase(**data)
    rules = load_rules()
    for name in parsed.settings.calibrations or []:
        rules.by_name(name)


_COMMENTED = os.path.join(BEES_PATH, "projects", "commented", "input.yml")


def test_commented_template_validates_and_uses_the_drug_keys():
    with open(_COMMENTED) as fh:
        data = yaml.safe_load(fh)
    parsed = InputBase(**data)
    enzyme = parsed.enzymes[0]
    assert enzyme.inhibition and enzyme.measured_kinetics and parsed.dosing
    assert any(sp.concentration == 0 for sp in parsed.species)


def test_commented_template_documents_every_schema_key():
    """A key the schema accepts but the template never mentions is undocumented."""
    from bees import schema

    with open(_COMMENTED) as fh:
        text = fh.read()
    models = [
        schema.Species, schema.Enzyme, schema.FeedbackInhibition, schema.Inhibition,
        schema.MeasuredKinetics, schema.Environment, schema.Settings, schema.Database, schema.Dosing,
        schema.InputBase,
    ]
    missing = sorted(
        f"{m.__name__}.{name}" for m in models for name in m.model_fields if f"{name}" not in text
    )
    assert not missing, f"keys missing from projects/commented/input.yml: {missing}"
