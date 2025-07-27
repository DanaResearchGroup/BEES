import os
import yaml
import pytest
import main
from unittest import mock

# Autouse fixture to patch paths and dependencies
@pytest.fixture(autouse=True)
def patch_paths_and_deps(tmp_path, monkeypatch):
    # Patch BEES_PATH and PROJECTS_BASE_PATH globals in main and common modules
    bees_path = tmp_path / "bees_path"
    projects_base = tmp_path / "projects_base"
    # Create necessary directories for fixed input
    (bees_path / 'examples' / 'minimal').mkdir(parents=True)
    projects_base.mkdir(parents=True)

    # Monkey-patch global constants
    monkeypatch.setattr(main, 'BEES_PATH', str(bees_path))
    monkeypatch.setattr(main, 'PROJECTS_BASE_PATH', str(projects_base))
    monkeypatch.setattr(main.common, 'BEES_PATH', str(bees_path))
    monkeypatch.setattr(main.common, 'PROJECTS_BASE_PATH', str(projects_base))

    # Monkey-patch Git utilities
    monkeypatch.setattr(main.common, 'get_git_branch', lambda: 'test-branch')
    monkeypatch.setattr(main.common, 'get_git_commit', lambda: ('abcdef', '2025-07-27'))

    yield

# -------- Tests for load_input_from_fixed_path() -------- #

def test_load_input_from_fixed_path_success():
    """Returns dict when input.yml is valid YAML"""
    fixed_dir = os.path.join(main.BEES_PATH, 'examples', 'minimal')
    input_file = os.path.join(fixed_dir, 'input.yml')
    valid_data = {
        'project': 'p',
        'species': [],
        'enzymes': [],
        'environment': {'temperature': 300},
        'database': {'name': 'db', 'solver': 'odeint'},
        'settings': {'end_time': 10, 'verbose': 20}
    }
    with open(input_file, 'w') as f:
        yaml.dump(valid_data, f)

    result = main.load_input_from_fixed_path()
    assert isinstance(result, dict)
    assert result['project'] == 'p'


def test_load_input_from_fixed_path_missing():
    """Raises FileNotFoundError when input.yml does not exist"""
    fixed_dir = os.path.join(main.BEES_PATH, 'examples', 'minimal')
    input_file = os.path.join(fixed_dir, 'input.yml')
    if os.path.exists(input_file):
        os.remove(input_file)
    with pytest.raises(FileNotFoundError):
        main.load_input_from_fixed_path()


def test_load_input_from_fixed_path_invalid_yaml():
    """Raises ValueError when input.yml is invalid YAML"""
    fixed_dir = os.path.join(main.BEES_PATH, 'examples', 'minimal')
    input_file = os.path.join(fixed_dir, 'input.yml')
    with open(input_file, 'w') as f:
        f.write("invalid: [unclosed list")
    with pytest.raises(ValueError):
        main.load_input_from_fixed_path()

# -------- Tests for BEES.__init__() -------- #

# Dummy schema to bypass real validation
class DummySchema:
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
    """Initialization succeeds: directories, logger, schema are set up"""
    input_data = {
        'project': 'testproj',
        'species': [{'label':'A'}],
        'enzymes': [{'name':'E'}],
        'environment': {'temperature': 310, 'pH': 7.4},
        'database': {'name': 'db', 'solver': 'odeint'},
        'settings': {'end_time': 15, 'verbose': 10, 'time_step': 0.5}
    }
    bees = main.BEES(input_data)

    # Directories exist
    assert os.path.isdir(bees.project_directory)
    # Schema picked up
    assert bees.input_schema.project == 'testproj'
    # Logger was instantiated
    MockLogger.assert_called_once()

@mock.patch('main.Logger')
def test_BEES_init_dir_fail(MockLogger, monkeypatch):
    """Raises FileNotFoundError if os.makedirs fails"""
    input_data = {'project':'proj'}
    # Force os.makedirs to error
    monkeypatch.setattr(os, 'makedirs', lambda *args, **kwargs: (_ for _ in ()).throw(OSError("deny")))
    with pytest.raises(FileNotFoundError):
        main.BEES(input_data)

@mock.patch('main.Logger')
@mock.patch('main.InputBase', side_effect=Exception("schema error"))
def test_BEES_init_schema_fail(MockInputBase, MockLogger):
    """Raises ValueError when schema validation fails"""
    data = {
        'project': 'x',
        'species': [{'label':'A'}],
        'enzymes': [{'name':'E'}],
        'environment': {'temperature':300, 'pH':7.0},
        'database': {'name':'db', 'solver':'odeint'},
        'settings': {'end_time':5, 'verbose':20, 'time_step':1.0}
    }
    with pytest.raises(ValueError):
        main.BEES(data)
