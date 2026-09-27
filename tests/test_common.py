""""
Test BEES common module


This module Test the common module which contains functions which are shared across multiple  modules.
To run the tests, use pytest and the command line: pytest -v tests/test_common.py


"""


import pytest
import os

import time
import yaml
from unittest.mock import patch
import re


# Import functions and constants from the common module
from bees.common import (
    get_git_branch,
    get_git_commit,
    InputError,
    read_yaml_file,
    save_yaml_file,
    to_yaml,
    time_lapse,
    dict_to_str,
    heavy_atom_count,
    get_ontology_equivalents,
)

def _reset_ontology_caches():
    import bees.common as common
    import bees.reaction_utils as reaction_utils
    common._CHEMICAL_ONTOLOGY = None
    common._ONTOLOGY_CATEGORIES_RAW = None
    common.get_ontology_equivalents.cache_clear()
    reaction_utils._EC_ALIASES_CACHE = None
    reaction_utils._ENZYME_DOMAIN_COFACTORS_CACHE = None


# Mock BEES_PATH for isolated testing
@pytest.fixture(autouse=True)
def mock_bees_paths(tmp_path):
    """
    Fixture to mock BEES_PATH to temporary directory
    for isolated testing of file operations.
    Also creates a dummy .git directory for git-related tests.
    """
    import shutil
    import bees.common as common

    real_root = common.BEES_PATH
    mock_bees_root = tmp_path / 'BEES_ROOT'
    mock_git_dir = mock_bees_root / '.git'

    os.makedirs(mock_bees_root, exist_ok=True)
    mock_git_dir.mkdir(exist_ok=True) # Create dummy .git directory for git tests

    # Ontology loaders cache by first hit. Without a copy, this fixture
    # caches an empty table and later test files see no acyl-CoA categories.
    db_src = os.path.join(real_root, 'db', 'ontology.yaml')
    db_dst = mock_bees_root / 'db'
    db_dst.mkdir(exist_ok=True)
    if os.path.exists(db_src):
        shutil.copy(db_src, db_dst / 'ontology.yaml')

    _reset_ontology_caches()
    with patch('bees.common.BEES_PATH', str(mock_bees_root)):
        yield
    _reset_ontology_caches()

def test_get_git_branch():
    """Test get_git_branch function."""
    # Mock subprocess.check_output for a controlled test environment
    with patch('subprocess.check_output') as mock_subproc:
        # Mocking the command to simulate git branch output
        mock_subproc.return_value = b'  main\n* develop\n'
        branch = get_git_branch()
        assert branch == 'develop'

        mock_subproc.return_value = b'* master\n'
        branch = get_git_branch()
        assert branch == 'master'

    # Test when .git directory does not exist (mock os.path.exists)
    with patch('os.path.exists', return_value=False):
        branch = get_git_branch()
        assert branch == ''

def test_get_git_commit():
    """Test get_git_commit function."""
    with patch('subprocess.check_output') as mock_subproc:
        # Mocking the command to simulate git log output
        mock_subproc.return_value = b'abcdef1234567890\nThu Jan 1 00:00:00 1970 +0000\n'
        commit, date = get_git_commit()
        assert commit == 'abcdef1234567890'
        assert date == 'Thu Jan 1 00:00:00 1970 +0000'

    # Test when .git directory does not exist (mock os.path.exists)
    with patch('os.path.exists', return_value=False):
        commit, date = get_git_commit()
        assert commit == ''
        assert date == ''

def test_InputError():
    """Test InputError exception."""
    with pytest.raises(InputError, match="This is a test error"):
        raise InputError("This is a test error")

def test_read_yaml_file(tmp_path):
    """Test read_yaml_file function."""
    # Test reading a valid YAML file
    test_yaml_content = {"key": "value", "number": 123}
    test_yaml_path = tmp_path / "test.yml"
    with open(test_yaml_path, "w") as f:
        yaml.dump(test_yaml_content, f)

    content = read_yaml_file(str(test_yaml_path))
    assert content == test_yaml_content

    # Test with project_directory for relative path resolution (file path itself)
    nested_dir = tmp_path / "project_dir"
    nested_dir.mkdir()
    nested_yaml_path = nested_dir / "nested.yml"
    with open(nested_yaml_path, "w") as f:
        yaml.dump({"nested_key": "nested_value"}, f)
    
    content = read_yaml_file("nested.yml", project_directory=str(nested_dir))
    assert content == {"nested_key": "nested_value"}

    # Test non-existent file
    with pytest.raises(InputError, match="Could not find the YAML file"):
        read_yaml_file(str(tmp_path / "non_existent.yml"))

    # Test invalid YAML content
    invalid_yaml_path = tmp_path / "invalid.yml"
    with open(invalid_yaml_path, "w") as f:
        f.write("key: - value") # Invalid YAML syntax (sequence entry where mapping key is expected)

    # The error message can vary slightly between PyYAML versions,
    # so a more general regex for YAMLError is appropriate.
    with pytest.raises(yaml.YAMLError):
        read_yaml_file(str(invalid_yaml_path))

    # Test invalid path type
    with pytest.raises(InputError, match="path must be a string"):
        read_yaml_file(123)


def test_get_ontology_equivalents_acp_thioester_permutations():
    eq = get_ontology_equivalents("3-oxo-(5Z)-dodecenoyl-[ACP]")
    assert "(5z)-3-oxododecenoyl-[acp]" in eq

    eq2 = get_ontology_equivalents("(5Z)-3-oxododecenoyl-[ACP]")
    assert "3-oxo-(5z)-dodecenoyl-[acp]" in eq2

def test_save_yaml_file(tmp_path):
    """Test save_yaml_file function."""
    output_path = tmp_path / "output.yml"
    data_to_save = {"data_key": "data_value", "list_data": [1, 2, 3]}
    save_yaml_file(str(output_path), data_to_save)

    assert output_path.exists()
    with open(output_path, "r") as f:
        loaded_data = yaml.safe_load(f)
    assert loaded_data == data_to_save

    # Test saving to a non-existent directory
    nested_output_dir = tmp_path / "new_dir" / "sub_dir"
    nested_output_path = nested_output_dir / "nested_output.yml"
    save_yaml_file(str(nested_output_path), {"nested": True})
    assert nested_output_path.exists()

    # Test invalid path type
    with pytest.raises(InputError, match="path must be a string"):
        save_yaml_file(123, {})

def test_to_yaml():
    """Test to_yaml function."""
    data = {"name": "Test", "value": 42}
    yaml_string = to_yaml(data)
    # Use safe_load to verify the output YAML string
    loaded_data = yaml.safe_load(yaml_string)
    assert loaded_data == data

    # Test with multiline string representation
    data_multiline = {"description": "This is a\nmultiline\nstring."}
    yaml_string_multiline = to_yaml(data_multiline) # Use to_yaml to apply the custom representer
    # The string_representer uses style='|'. PyYAML often adds '|- ' for block literals.
    # The regex should be flexible for leading spaces and the exact chomping indicator.
    # It should match "description: " followed by optional " |-" and then the multiline content.
    assert re.search(r"description:\s*\|-?\s*\n\s*This is a\n\s*multiline\n\s*string\.?", yaml_string_multiline) is not None


def test_string_representer():
    """Test string_representer for YAML multiline strings."""
    # This function is typically used by yaml.dump, so we'll test it indirectly
    # by dumping a dictionary with a multiline string.
    data = {"key": "single line string"}
    yaml_str = yaml.dump(data, Dumper=yaml.Dumper, default_flow_style=False)
    assert "single line string" in yaml_str
    assert "|" not in yaml_str # Should not use block literal for single line

    data_multiline = {"description": "This is a\nmultiline\nstring."}
    yaml_str_multiline = to_yaml(data_multiline) # Use to_yaml to apply the custom representer
    # The regex checks for the key, optional chomping indicator and the multiline content
    assert re.search(r"description:\s*\|-?\s*\n\s*This is a\n\s*multiline\n\s*string\.?", yaml_str_multiline) is not None


"""
Bond-length and distance-matrix utilities were removed from `bees.common` since they are not
used by BEES' core pipeline. The associated tests were intentionally removed.
"""

def test_time_lapse():
    """Test time_lapse function."""
    t0 = time.time()
    time.sleep(1.1) # Sleep for a bit to ensure non-zero time
    elapsed_time_str = time_lapse(t0)
    # Check format, can't check exact time due to execution variations
    assert re.match(r'(\d+ days, )?\d{2}:\d{2}:\d{2}', elapsed_time_str) is not None

    # Test with a known time difference for specific formatting
    mock_start_time = time.time() - (2 * 24 * 3600 + 3 * 3600 + 4 * 60 + 5) # 2 days, 03:04:05
    with patch('time.time', return_value=mock_start_time + (2 * 24 * 3600 + 3 * 3600 + 4 * 60 + 5)):
        assert time_lapse(mock_start_time) == "2 days, 03:04:05"

    mock_start_time_no_days = time.time() - (3 * 3600 + 4 * 60 + 5) # 03:04:05
    with patch('time.time', return_value=mock_start_time_no_days + (3 * 3600 + 4 * 60 + 5)):
        assert time_lapse(mock_start_time_no_days) == "03:04:05"


def test_dict_to_str():
    """Test dict_to_str function."""
    test_dict = {
        "key1": "value1",
        "key2": {
            "nested_key1": 123,
            "nested_key2": "abc"
        },
        "key3": True
    }
    expected_str = (
        "key1: value1\n"
        "key2:\n"
        "  nested_key1: 123\n"
        "  nested_key2: abc\n"
        "key3: True\n"
    )
    assert dict_to_str(test_dict) == expected_str

    # Test with empty dict
    assert dict_to_str({}) == ""

    # Test with different level
    test_dict_level = {"outer": {"inner": "val"}}
    expected_str_level = (
        "outer:\n"
        "  inner: val\n"
    )
    assert dict_to_str(test_dict_level, level=0) == expected_str_level

def test_heavy_atom_count_valid_smiles():
    # ethanol: C-C-O = 3 heavy atoms
    assert heavy_atom_count("CCO") == 3


def test_heavy_atom_count_empty_or_invalid_returns_none():
    assert heavy_atom_count("") is None
    assert heavy_atom_count(None) is None
    assert heavy_atom_count("not_a_smiles") is None

