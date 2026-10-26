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
