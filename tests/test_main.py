"""
Tests for BEES main module (parser + BEES init)

Run: pytest -v tests/test_main.py
"""

import os
import yaml
import pytest
from unittest import mock

import main
from bees.common import InputError  # raised by read_yaml_file on missing file


# --------------------------
# Shared fixtures/utilities
# --------------------------

@pytest.fixture(autouse=True)
def patch_paths_and_deps(tmp_path, monkeypatch):
    """
    - Redirect BEES_PATH and PROJECTS_BASE_PATH to temp dirs.
    - Create examples/minimal/ for the parser's default path.
    - Stub git helpers for deterministic logs.
    """
    bees_path = tmp_path / "bees_path"
    projects_base = tmp_path / "projects_base"

    (bees_path / 'examples' / 'minimal').mkdir(parents=True)
    projects_base.mkdir(parents=True)

    # Patch globals used by main.py and its imported common
    monkeypatch.setattr(main, 'BEES_PATH', str(bees_path))
    monkeypatch.setattr(main, 'PROJECTS_BASE_PATH', str(projects_base))
    monkeypatch.setattr(main.common, 'BEES_PATH', str(bees_path))
    monkeypatch.setattr(main.common, 'PROJECTS_BASE_PATH', str(projects_base))

    # Stable git info
    monkeypatch.setattr(main.common, 'get_git_branch', lambda: 'test-branch')
    monkeypatch.setattr(main.common, 'get_git_commit', lambda: ('abcdef', '2025-07-27'))

    yield


def _write_yaml(path, data):
    with open(path, "w") as f:
        yaml.dump(data, f)


# --------------------------
# Tests: parse_and_load_input
# --------------------------

def test_parser_uses_default_input_when_none_provided(tmp_path, monkeypatch):
    """
    No --input_file -> load BEES_PATH/examples/minimal/input.yml
    """
    default_dir = os.path.join(main.BEES_PATH, 'examples', 'minimal')
    default_yaml = os.path.join(default_dir, 'input.yml')
    _write_yaml(default_yaml, {
        'project': 'DefaultProj',
        'species': [],
        'enzymes': [],
        'environment': {'temperature': 300},
        'database': {'name': 'db', 'solver': 'odeint'},
        'settings': {'end_time': 10, 'verbose': 20},
    })

    import sys
    old_argv = sys.argv
    sys.argv = ['main.py']  # simulate: python main.py
    try:
        out = main.parse_and_load_input()
    finally:
        sys.argv = old_argv

    assert out['project'] == 'DefaultProj'
    assert out['settings']['verbose'] == 20
    assert out['settings']['end_time'] == 10


def test_parser_cli_overrides_yaml(tmp_path):
    """
    --project, --verbose, --output_directory override YAML.
    """
    yml = tmp_path / "custom.yml"
    _write_yaml(yml, {
        'project': 'YAMLProj',
        'species': [],
        'enzymes': [],
        'environment': {'temperature': 310},
        'database': {'name': 'db', 'solver': 'odeint'},
        'settings': {'end_time': 100, 'time_step': 0.5, 'verbose': 30}
    })

    import sys
    old_argv = sys.argv
    sys.argv = [
        'main.py',
        '--input_file', str(yml),
        '--project', 'CLIProj',
        '--verbose', '10',
        '--output_directory', 'results/run_01',
    ]
    try:
        out = main.parse_and_load_input()
    finally:
        sys.argv = old_argv

    assert out['project'] == 'CLIProj'
    assert out['settings']['verbose'] == 10
    assert out['settings']['output_directory'] == 'results/run_01'
    # unchanged YAML bits
    assert out['settings']['end_time'] == 100
    assert out['settings']['time_step'] == 0.5
    assert out['environment']['temperature'] == 310


def test_parser_missing_yaml_raises_inputerror(tmp_path):
    """
    Non-existent --input_file -> bees.common.read_yaml_file raises InputError.
    """
    import sys
    old_argv = sys.argv
    sys.argv = ['main.py', '--input_file', str(tmp_path / 'does_not_exist.yml')]
    try:
        with pytest.raises(InputError):
            main.parse_and_load_input()
    finally:
        sys.argv = old_argv


def test_parser_requires_project_when_yaml_omits_it(tmp_path):
    """
    YAML without 'project' and no --project -> ValueError (guard in parser).
    """
    yml = tmp_path / "no_project.yml"
    _write_yaml(yml, {
        'species': [],
        'enzymes': [],
        'environment': {'temperature': 300},
        'database': {'name': 'db', 'solver': 'odeint'},
        'settings': {'end_time': 10, 'verbose': 20},
    })

    import sys
    old_argv = sys.argv
    sys.argv = ['main.py', '--input_file', str(yml)]
    try:
        with pytest.raises(ValueError, match="Project name is required"):
            main.parse_and_load_input()
    finally:
        sys.argv = old_argv


def test_parser_injects_settings_when_missing(tmp_path):
    """
    If YAML lacks 'settings', parser should create it so overrides can apply.
    """
    yml = tmp_path / "no_settings.yml"
    _write_yaml(yml, {
        'project': 'P',
        'species': [],
        'enzymes': [],
        'environment': {'temperature': 300},
        'database': {'name': 'db', 'solver': 'odeint'},
        # no 'settings'
    })

    import sys
    old_argv = sys.argv
    sys.argv = ['main.py', '--input_file', str(yml), '--verbose', '40', '--output_directory', 'out']
    try:
        out = main.parse_and_load_input()
    finally:
        sys.argv = old_argv

    assert out['project'] == 'P'
    # parser created settings and applied overrides
    assert out['settings']['verbose'] == 40
    assert out['settings']['output_directory'] == 'out'


# --------------------------
# Tests: BEES.__init__
# --------------------------

class DummySchema:
    """
    Minimal stand-in for bees.schema.InputBase so we don't rely on full validation here.
    """
    def __init__(self, **kwargs):
        self.project = kwargs.get('project')
        self.species = kwargs.get('species', [])
        self.enzymes = kwargs.get('enzymes', [])
        class E: pass
        self.environment = E()
        env = kwargs.get('environment', {})
        self.environment.temperature = env.get('temperature')
        self.environment.pH = env.get('pH', 7.0)
        class D: pass
        self.database = D()
        db = kwargs.get('database', {})
        self.database.name = db.get('name')
        self.database.solver = db.get('solver')
        self.database.rate_law = db.get('rate_law', None)
        self.database.parameter_estimator = db.get('parameter_estimator', None)
        class S: pass
        self.settings = S()
        st = kwargs.get('settings', {})
        self.settings.end_time = st.get('end_time')
        self.settings.verbose = st.get('verbose')
        self.settings.time_step = st.get('time_step', 1.0)
    def model_dump(self, exclude_unset):
        return {}

@mock.patch('main.Logger')
@mock.patch('main.InputBase', new=DummySchema)
def test_BEES_init_success(MockLogger):
    """
    Creates directories, configures logger, and stores schema.
    """
    data = {
        'project': 'testproj',
        'species': [{'label': 'A'}],
        'enzymes': [{'name': 'E'}],
        'environment': {'temperature': 310, 'pH': 7.4},
        'database': {'name': 'db', 'solver': 'odeint'},
        'settings': {'end_time': 15, 'verbose': 10, 'time_step': 0.5},
    }
    bees = main.BEES(data)

    assert os.path.isdir(bees.project_directory)
    assert bees.input_schema.project == 'testproj'
    MockLogger.assert_called_once()


@mock.patch('main.Logger')
def test_BEES_init_dir_fail(MockLogger, monkeypatch):
    """
    If directory creation fails, BEES raises FileNotFoundError.
    """
    data = {'project': 'proj'}
    monkeypatch.setattr(os, 'makedirs', lambda *a, **k: (_ for _ in ()).throw(OSError("deny")))
    with pytest.raises(FileNotFoundError):
        main.BEES(data)


@mock.patch('main.Logger')
@mock.patch('main.InputBase', side_effect=Exception("schema error"))
def test_BEES_init_schema_fail(MockInputBase, MockLogger):
    """
    If schema validation fails, BEES raises ValueError.
    """
    data = {
        'project': 'x',
        'species': [{'label': 'A'}],
        'enzymes': [{'name': 'E'}],
        'environment': {'temperature': 300, 'pH': 7.0},
        'database': {'name': 'db', 'solver': 'odeint'},
        'settings': {'end_time': 5, 'verbose': 20, 'time_step': 1.0},
    }
    with pytest.raises(ValueError):
        main.BEES(data)
