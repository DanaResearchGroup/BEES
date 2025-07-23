"""
BEES common module
"""
"""
This module contains functions which are shared across multiple  modules.
As such, it should not import any other BEES module (specifically ones that use the logger defined here)
to avoid circular imports.

This VERSION based on is the full ARC version, using `semantic versioning <https://semver.org/>`_.
"""

import os  # Essential for path manipulations.
import subprocess
from typing import TYPE_CHECKING, Any, List, Optional, Tuple, Union
import datetime
import shutil

import time

import yaml
import numpy as np
import math # Added this import for math.exp
import re
from collections import deque # Keep deque for potential future use or if other commented functions use it


"""
#TODO:
1. learn each part of the code how its work and if its neccery in our code
2. delete/adopt which module that is neccery to BEES
3. writing test.

"""
 
# Absolute path to the BEES folder.

BEES_PATH = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))

# Define base paths for projects and data relative to BEES_PATH
# These paths are used to store project files and data files.

PROJECTS_BASE_PATH = os.path.join(BEES_PATH, 'projects')
DATA_BASE_PATH = os.path.join(BEES_PATH, 'data')

VERSION = '1.1.0'



# Constants
R = 8.31446261815324  # J/(mol*K)
EA_UNIT_CONVERSION = {'J/mol': 1, 'kJ/mol': 1e+3, 'cal/mol': 4.184, 'kcal/mol': 4.184e+3}



#All the functions in the common module



def get_git_branch(path: Optional[str] = None) -> str:
    """
    Get the git branch to be logged.

    Args:
        path (str, optional): The path to check.

    Returns: str
        The git branch name.
    """
    path = path or BEES_PATH
    if os.path.exists(os.path.join(path, '.git')):
        try:
            branch_list = subprocess.check_output(['git', 'branch'], cwd=path).splitlines()
        except (subprocess.CalledProcessError, OSError):
            return ''
        for branch_name in branch_list:
            if '*' in branch_name.decode():
                return branch_name.decode()[2:]
    else:
        return ''


def get_git_commit(path: Optional[str] = None) -> Tuple[str, str]:
    """
    Get the recent git commit to be logged.

    Note:
        Returns empty strings if hash and date cannot be determined.

    Args:
        path (str, optional): The path to check.

    Returns: tuple
        The git HEAD commit hash and the git HEAD commit date, each as a string.
    """
    path = path or BEES_PATH
    head, date = '', ''
    if os.path.exists(os.path.join(path, '.git')):
        try:
            head, date = subprocess.check_output(['git', 'log', '--format=%H%n%cd', '-1'], cwd=path).splitlines()
            head, date = head.decode(), date.decode()
        except (subprocess.CalledProcessError, OSError):
            return head, date
    return head, date

class InputError(Exception):
    """An exception class for reporting errors in BEES input files or parameters. 
    TODO : we might move this one out of the common module to a separate file in the bees module."""
    pass


def globalize_paths(file_path: str,
                    project_directory: str,
                    ) -> str:
    """
    Rebase all file paths in the contents of the given file on the current project path.
    Useful when restarting an BEES project in a different folder or on a different machine.

    Args:
        file_path (str): A path to the file to check.
                         The contents of this file will be changed and saved as a different file.
        project_directory (str): The current project directory to rebase upon.

    Returns: str
        A path to the respective file with rebased absolute file paths.
    """
    modified = False
    new_lines = []
    # Ensure project_directory has a trailing slash for consistent path construction
    # and normalize path to handle different OS separators
    normalized_project_directory = os.path.normpath(project_directory).rstrip(os.sep) + os.sep

    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        # Pass the original line to globalize_path, it will handle internal stripping/normalization
        rebased_line = globalize_path(line, normalized_project_directory)

        if line != rebased_line: # Check if the line was actually changed by globalize_path
            modified = True
        new_lines.append(rebased_line)
            
    if modified:
        base_name, file_name = os.path.split(file_path)
        file_name_splits = file_name.split('.')
        new_file_name = '.'.join(file_name_splits[:-1]) + '_globalized.' + str(file_name_splits[-1])
        new_path = os.path.join(base_name, new_file_name)
        with open(new_path, 'w') as f:
            f.writelines(new_lines)
        return new_path
    else:
        return file_path


def globalize_path(string: str,
                   project_directory: str, # This should already be normalized and end with a slash from globalize_paths
                   ) -> str:
    """
    Rebase an absolute file path on the current project path.
    Useful when restarting an BEES project in a different folder or on a different machine.

    Args:
        string (str): A string containing a path to rebase.
        project_directory (str): The current project directory to rebase upon.
                                 Expected to be normalized and end with a slash.

    Returns: str
        A string with the rebased path, or the original string if no change is needed.
    """
    # project_directory is expected to be normalized and end with a slash already
    normalized_project_directory = project_directory 

    # Regex to capture:
    # Group 'leading_ws': Optional leading whitespace
    # Group 'key_prefix': Optional YAML key and colon (e.g., "key: ")
    # Group 'path_value': The path value itself (non-greedy, matches until end of line or before trailing whitespace)
    # Group 'trailing_ws_before_newline': Any whitespace characters immediately before the final newline (if present)
    # Group 'newline': The actual newline character (if present)
    path_pattern = r'^(?P<leading_ws>\s*)(?P<key_prefix>\w+\s*:\s*)?(?P<path_value>.*?)(?P<trailing_ws_before_newline>\s*)(?P<newline>\n?)$'
    match = re.match(path_pattern, string)

    if not match:
        return string # Return original string if it doesn't match the expected pattern

    leading_ws = match.group('leading_ws')
    key_prefix = match.group('key_prefix') if match.group('key_prefix') else ""
    raw_path_value = match.group('path_value')
    trailing_ws_before_newline = match.group('trailing_ws_before_newline')
    newline = match.group('newline') # This will be '\n' or ''

    # Normalize the raw path value for comparison
    normalized_raw_path_value = os.path.normpath(raw_path_value).rstrip(os.sep)

    sub_path_marker = None
    if '/calcs/Species/' in normalized_raw_path_value:
        sub_path_marker = '/calcs/Species/'
    elif '/calcs/TSs/' in normalized_raw_path_value:
        sub_path_marker = '/calcs/TSs/'

    if sub_path_marker:
        # Split the normalized raw path value to get the old root and the specific path
        parts_of_path = normalized_raw_path_value.split(sub_path_marker, 1)
        
        # Construct the path that *would* be rebased from the old root
        old_root_candidate = parts_of_path[0]
        
        # Normalize the old root candidate for comparison
        normalized_old_root_candidate = os.path.normpath(old_root_candidate).rstrip(os.sep) + os.sep if old_root_candidate else ''

        # If the old root is already the target project directory, no change is needed.
        if normalized_old_root_candidate == normalized_project_directory:
            return string # Return original string as no effective change is required.

        # If it's not already matching the target project_directory, then rebase.
        new_path_value = normalized_project_directory + sub_path_marker[1:] + parts_of_path[1]
        result = leading_ws + key_prefix + new_path_value + trailing_ws_before_newline + newline
        return result

    else:
        # Only rebase if it's explicitly the 'project_directory: ' field, or a standalone absolute path that *is* the old project root
        if key_prefix == 'project_directory: ':
            current_project_dir_in_string = raw_path_value
            normalized_current_in_string = os.path.normpath(current_project_dir_in_string).rstrip(os.sep) + os.sep

            if normalized_current_in_string == normalized_project_directory:
                return string # Already matches, no change needed
            else:
                result = leading_ws + key_prefix + normalized_project_directory + trailing_ws_before_newline + newline
                return result
        elif os.path.isabs(normalized_raw_path_value) and \
             normalized_raw_path_value + os.sep == os.path.normpath('/old/project/root/').rstrip(os.sep) + os.sep: # Check if it's the specific old root
            # This handles the case where the line is just the old project root path, e.g., "/old/project/root/\n"
            new_path_value = normalized_project_directory.rstrip(os.sep) + os.sep
            result = leading_ws + key_prefix + new_path_value + trailing_ws_before_newline + newline
            return result
        
        # If it's not a path with a sub_path_marker, not a 'project_directory:' line,
        # and not the specific old project root path, return original string.
        return string


def read_yaml_file(path: str,
                   project_directory: Optional[str] = None,
                   ) -> Union[dict, list]:
    """
    Read a YAML file (usually an input / restart file, but also conformers file)
    and return the parameters as python variables.

    Args:
        path (str): The YAML file path to read.
        project_directory (str, optional): The current project directory to rebase upon.
                                           Used to resolve relative 'path' if provided.

    Returns: Union[dict, list]
        The content read from the file.
    """
    if not isinstance(path, str):
        raise InputError(f'path must be a string, got {path} which is a {type(path)}')
    
    original_path_for_error = path # Keep original path for error message

    # --- FIX: Resolve the file path if it's relative and project_directory is provided ---
    if project_directory is not None and not os.path.isabs(path):
        path = os.path.join(project_directory, path)
    # --- END FIX ---

    if not os.path.isfile(path):
        raise InputError(f'Could not find the YAML file {original_path_for_error} (resolved to {path})')
    
    with open(path, 'r') as f:
        content = yaml.load(stream=f, Loader=yaml.FullLoader)
    
    # No integrated globalization logic here.
    # Globalization of content should be handled by calling globalize_paths separately if needed.

    return content # <--- FIXED: Ensure content is returned


def save_yaml_file(path: str,
                   content: Union[list, dict],
                   ) -> None:
    """
    Save a YAML file (usually an input / restart file, but also conformers file).

    Args:
        path (str): The YAML file path to save.
        content (list, dict): The content to save.
    """
    if not isinstance(path, str):
        raise InputError(f'path must be a string, got {path} which is a {type(path)}')
    yaml_str = to_yaml(py_content=content)
    if '/' in path and os.path.dirname(path) and not os.path.exists(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))
    with open(path, 'w') as f:
        f.write(yaml_str)


def from_yaml(yaml_str: str) -> Union[dict, list]:
    """
    Read a YAML string and decode to the respective Python object.
    Args:
        yaml_str (str): The YAML string content.
    Returns: Union[dict, list]
        The respective Python object.
    """
    return yaml.load(stream=yaml_str, Loader=yaml.FullLoader)


def to_yaml(py_content: Union[list, dict]) -> str:
    """
    Convert a Python list or dictionary to a YAML string format.

    Args:
        py_content (list, dict): The Python content to save.

    Returns: str
        The corresponding YAML representation.
    """
    yaml.add_representer(str, string_representer)
    yaml_str = yaml.dump(data=py_content)
    return yaml_str


def string_representer(dumper, data):
    """
    Add a custom string representer to use block literals for multiline strings.
    """
    if len(data.splitlines()) > 1:
        return dumper.represent_scalar(tag='tag:yaml.org,2002:str', value=data, style='|')
    return dumper.represent_scalar(tag='tag:yaml.org,2002:str', value=data)


def get_ordinal_indicator(number: int) -> str:
    """
    Returns the ordinal indicator for an integer.

    Args:
        number (int): An integer for which the ordinal indicator will be determined.

    Returns: str
        The integer's ordinal indicator.
    """
    ordinal_dict = {1: 'st', 2: 'nd', 3: 'rd'}
    if number > 13:
        number %= 10
    if number in list(ordinal_dict.keys()):
        return ordinal_dict[number]
    return 'th'


def get_number_with_ordinal_indicator(number: int) -> str:
    """
    Returns the number as a string with the ordinal indicator.

    Args:
        number (int): An integer for which the ordinal indicator will be determined.

    Returns: str
        The number with the respective ordinal indicator.
    """
    return f'{number}{get_ordinal_indicator(number)}'





# A bond length dictionary of single bonds in Angstrom.
# https://sites.google.com/site/chempendix/bond-lengths
# https://courses.lumenlearning.com/suny-potsdam-organicchemistry/chapter/1-3-basics-of-bonding/
# 'N-O' is taken from NH2OH, 'N+_N+' and 'N+_O-' are taken from N2O4.
# 'H_H' was artificially modified from 0.74 to 1.0 since it collides quickly at 0.55.
SINGLE_BOND_LENGTH = {'Br_Br': 2.29, 'Br_Cr': 1.94, 'Br_H': 1.41,
                      'C_C': 1.54, 'C_Cl': 1.77, 'C_F': 1.35, 'C_H': 1.09, 'C_I': 2.13,
                      'C_N': 1.47, 'C_O': 1.43, 'C_P': 1.87, 'C_S': 1.81, 'C_Si': 1.86,
                      'Cl_Cl': 1.99, 'Cl_H': 1.27, 'Cl_N': 1.75, 'Cl_Si': 2.03, 'Cl_P': 2.03, 'Cl_S': 2.07,
                      'F_F': 1.42, 'F_H': 0.92, 'F_P': 1.57, 'F_S': 1.56, 'F_Si': 1.56, 'F_Xe': 1.90,
                      'H_H': 1.0, 'H_I': 1.61, 'H_N': 1.04, 'H_O': 0.96, 'H_P': 1.42, 'H_S': 1.34, 'H_Si': 1.48,
                      'I_I': 2.66,
                      'N_N': 1.45, 'N+1_N+1': 1.81, 'N_O': 1.44, 'N+1_O-1': 1.2,
                      'O_O': 1.48, 'O_P': 1.63, 'O_S': 1.58, 'O_Si': 1.66,
                      'P_P': 2.21,
                      'S_S': 2.05,
                      'Si_Si': 2.35,
                      }


def get_single_bond_length(symbol_1: str,
                           symbol_2: str,
                           charge_1: int = 0,
                           charge_2: int = 0,
                           ) -> float:
    """
    Get an approximate for a single bond length between two elements.

    Args:
        symbol_1 (str): Symbol 1.
        symbol_2 (str): Symbol 2.
        charge_1 (int, optional): The partial charge of the atom represented by ``symbol_1``.
        charge_2 (int, optional): The partial charge of the atom represented by ``symbol_2``.

    Returns: float
        The estimated single bond length in Angstrom.
    """
    if charge_1 and charge_2:
        symbol_1 = f"{symbol_1}{'+' if charge_1 > 0 else ''}{charge_1}"
        symbol_2 = f"{symbol_2}{'+' if charge_2 > 0 else ''}{charge_2}"
    bond1, bond2 = '_'.join([symbol_1, symbol_2]), '_'.join([symbol_2, symbol_1])
    if bond1 in SINGLE_BOND_LENGTH.keys():
        return SINGLE_BOND_LENGTH[bond1]
    if bond2 in SINGLE_BOND_LENGTH.keys():
        return SINGLE_BOND_LENGTH[bond2]
    return 2.5


def get_bonds_from_dmat(dmat: np.ndarray,
                        elements: Union[Tuple[str, ...], List[str]],
                        charges: Optional[List[int]] = None,
                        tolerance: float = 1.2,
                        bond_lone_hydrogens: bool = True,
                        ) -> List[Tuple[int, int]]:
    """
    Guess the connectivity of a molecule from its distance matrix representation.

    Args:
        dmat (np.ndarray): An NxN matrix of atom distances in Angstrom.
        elements (List[str]): The corresponding element list in the same atomic order.
        charges (List[int], optional): A corresponding list of formal atomic charges.
        tolerance (float, optional): A factor by which the single bond length threshold is multiplied for the check.
        bond_lone_hydrogens (bool, optional): Whether to assign a bond to hydrogen atoms which were not identified
                                              as bonded. If so, the closest atom will be considered.

    Returns:
        List[Tuple[int, int]]:
            A list of tuple entries, each represents a bond and contains sorted atom indices.
    """
    if len(elements) != dmat.shape[0] or len(elements) != dmat.shape[1] or len(dmat.shape) != 2:
        raise ValueError(f'The dimensions of the DMat {dmat.shape} must be equal to the number of elements {len(elements)}')
    bonds, bonded_hydrogens = list(), list()
    charges = charges or [0] * len(elements)
    # Heavy atoms
    for i, e_1 in enumerate(elements):
        for j, e_2 in enumerate(elements):
            if i > j and e_1 != 'H' and e_2 != 'H' and dmat[i, j] < tolerance * \
                    get_single_bond_length(symbol_1=e_1,
                                           symbol_2=e_2,
                                           charge_1=charges[i],
                                           charge_2=charges[j]):
                bonds.append(tuple(sorted([i, j])))
    # Hydrogen atoms except for H--H bonds, make sure each H has only one bond (the closest heavy atom).
    for i, e_1 in enumerate(elements):
        if e_1 == 'H':
            j = get_extremum_index(lst=dmat[i], return_min=True, skip_values=[0])
            if i != j and (elements[j] != 'H'):
                bonds.append(tuple(sorted([i, j])))
                bonded_hydrogens.append(i)
    # Lone hydrogens, also important for the H2 molecule.
    if bond_lone_hydrogens:
        for i, e_1 in enumerate(elements):
            j = get_extremum_index(lst=dmat[i], return_min=True, skip_values=[0])
            bond = tuple(sorted([i, j]))
            if i != j and e_1 == 'H' and i not in bonded_hydrogens and j not in bonded_hydrogens and bond not in bonds:
                bonds.append(bond)
    return bonds




def extremum_list(lst: list,
                  return_min: bool = True,
                  ) -> Optional[Union[int, None]]:
    """
    A helper function for finding the extremum (either minimum or maximum) of a list of numbers (int/float)
    where some entries could be ``None``.

    Args:
        lst (list): The list.
        return_min (bool, optional): Whether to return the minimum or the maximum.
                                    ``True`` for minimum, ``False`` for maximum, ``True`` by default.

    Returns: Optional[Union[int, None]]
        The entry with the minimal/maximal value.
    """
    if lst is None or len(lst) == 0 or all([entry is None for entry in lst]):
        return None
    elif len(lst) == 1:
        return lst[0]
    if return_min:
        return min([entry for entry in lst if entry is not None])
    else:
        return max([entry for entry in lst if entry is not None])


def get_extremum_index(lst: list,
                       return_min: bool = True,
                       skip_values: Optional[list] = None
                       ) -> Optional[Union[int, None]]:
    """
    A helper function for finding the extremum (either minimum or maximum) of a list of numbers (int/float)
    where some entries could be ``None``.

    Args:
        lst (list): The list.
        return_min (bool, optional): Whether to return the minimum or the maximum.
                                    ``True`` for minimum, ``False`` for maximum, ``True`` by default.
        skip_values (list, optional): Values to skip when checking for extermum, e.g., 0.

    Returns: Optional[Union[int, None]]
        The index of an entry with the minimal/maximal value.
    """
    if len(lst) == 0:
        return None
    elif all([entry is None for entry in lst]):
        return None
    elif len(lst) == 1:
        return 0
    
    skip_values_set = set(skip_values) if skip_values is not None else set()
    skip_values_set.add(None) # Always skip None values

    # Initialize extremum_index with the first valid (non-skipped) index
    extremum_index = -1
    for i, entry in enumerate(lst):
        if entry not in skip_values_set:
            extremum_index = i
            break
    
    if extremum_index == -1: # No valid entries found
        return None

    # Iterate from the first valid entry
    for i in range(extremum_index + 1, len(lst)):
        entry = lst[i]
        if entry in skip_values_set:
            continue
        # Ensure that lst[extremum_index] is not None before comparison
        if return_min:
            if lst[extremum_index] is None or (entry is not None and entry < lst[extremum_index]):
                extremum_index = i
        else:
            if lst[extremum_index] is None or (entry is not None and entry > lst[extremum_index]):
                extremum_index = i
    return extremum_index


def sum_list_entries(lst: List[Union[int, float]],
                     multipliers: Optional[List[Union[int, float]]] = None,
                     ) -> Optional[float]:
    """
    Sum all entries in a list. If any entry is ``None``, return ``None``.
    If ``multipliers`` is given, multiply each entry in ``lst`` by the respective multiplier entry.

    Args:
        lst (list): The list to process.
        multipliers (list, optional): A list of multipliers.

    Returns:
        Optional[float]: The result.
    """
    if any(entry is None or not isinstance(entry, (int, float)) for entry in lst):
        return None
    if multipliers is None:
        return sum(lst)
    return float(np.dot(lst, multipliers + [1] * (len(lst) - len(multipliers))))


def sort_two_lists_by_the_first(list1: List[Union[float, int, None]],
                                list2: List[Union[float, int, None]],
                                ) -> Tuple[List[Union[float, int]], List[Union[float, int]]]:
    """
    Sort two lists in increasing order by the values of the first list.
    Ignoring None entries from list1 and their respective entries in list2.
    The function was written in this format rather the more pytonic ``zip(*sorted(zip(list1, list2)))`` style
    to accommodate for dictionaries as entries of list2, otherwise a
    ``TypeError: '<' not supported between instances of 'dict' and 'dict'`` error is raised.

    Args:
        list1 (list, tuple): Entries are floats or ints (could also be None).
        list2 (list, tuple): Entries could be anything.

    Raises:
        InputError: If types are wrong, or lists are not the same length.

    Returns: Tuple[list, list]
        - Sorted values from list1, ignoring None entries.
        - Respective entries from list2.
    """
    if not isinstance(list1, (list, tuple)) or not isinstance(list2, (list, tuple)):
        raise InputError(f'Arguments must be lists, got: {type(list1)} and {type(list2)}')
    for entry in list1:
        if not isinstance(entry, (float, int)) and entry is not None:
            raise InputError(f'Entries of list1 must be either floats or integers, got: {type(entry)}.')
    if len(list1) != len(list2):
        raise InputError(f'Both lists must be the same length, got {len(list1)} and {len(list2)}')

    # remove None entries from list1 and their respective entries from list2:
    new_list1, new_list2 = list(), list()
    for entry1, entry2 in zip(list1, list2):
        if entry1 is not None:
            new_list1.append(entry1)
            new_list2.append(entry2)
    indices = list(range(len(new_list1)))

    zipped_lists = zip(new_list1, indices)
    sorted_lists = sorted(zipped_lists)
    sorted_list1 = [x for x, _ in sorted_lists]
    sorted_indices = [x for _, x in sorted_lists]
    sorted_list2 = [0] * len(new_list2)
    for counter, index in enumerate(sorted_indices):
        sorted_list2[counter] = new_list2[index]
    return sorted_list1, sorted_list2


def check_that_all_entries_are_in_list(list_1: Union[list, tuple], list_2: Union[list, tuple]) -> bool:
    """
    Check that all entries from ``list_2`` are in ``list_1``, and that the lists are the same length.
    Useful for testing that two lists are equal regardless of entry order.

    Args:
        list_1 (list, tuple): Entries are floats or ints (could also be None).
        list_2 (list, tuple): Entries could be anything.

    Returns: bool
        Whether all entries from ``list_2`` are in ``list_1`` and the lists are the same length.
    """
    if len(list_1) != len(list_2):
        return False
    for entry in list_2:
        if entry not in list_1:
            return False
    return True


def key_by_val(dictionary: dict,
               value: Any,
               ) -> Any:
    """
    A helper function for getting a key from a dictionary corresponding to a certain value.
    Does not check for value unicity.

    Args:
        dictionary (dict): The dictionary.
        value: The value.

    Raises:
        ValueError: If the value could not be found in the dictionary.

    Returns: Any
        The key.
    """
    for key, val in dictionary.items():
        # The 'X' prefix handling is only for integer values, where val is a string like 'X2'
        if val == value or (isinstance(value, int) and isinstance(val, str) and val == f'X{value}'):
            return key
    raise ValueError(f'Could not find value {value} in the dictionary\n{dictionary}')


def is_str_float(value: Optional[str]) -> bool:
    """
    Check whether a string can be converted to a floating number.

    Args:
        value (str): The string to check.

    Returns: bool
        ``True`` if it can, ``False`` otherwise.
    """
    try:
        float(value)
        return True
    except (ValueError, TypeError):
        return False


def is_str_int(value: Optional[str]) -> bool:
    """
    Check whether a string can be converted to an integer.

    Args:
        value (str): The string to check.

    Returns: bool
        ``True`` if it can, ``False`` otherwise.
    """
    try:
        int(value)
        return True
    except (ValueError, TypeError):
        return False


def clean_text(text: str) -> str:
    """
    Clean a text string from leading and trailing whitespaces, newline characters, and double quotes.
    Also removes trailing commas.

    Args:
        text (str): The text to clean.

    Returns:
        str: The cleaned text.
    """
    # Remove leading/trailing quotes and newlines, then strip spaces
    text = re.sub(r'^\s*["\n]*|["\n]*\s*$', '', text)
    # Remove trailing commas
    text = text.rstrip(',')
    # Final strip to remove any remaining whitespace (e.g., from original inner newlines)
    text = text.strip()
    return text


def time_lapse(t0) -> str:
    """
    A helper function returning the elapsed time since t0.

    Args:
        t0 (float): The initial time the count starts from.

    Returns: str
        A "D HH:MM:SS" formatted time difference between now and t0.
    """
    t = time.time() - t0
    m, s = divmod(t, 60)
    h, m = divmod(m, 60)
    d, h = divmod(h, 24)
    if d > 0:
        d_str = str(int(d)) + ' days, '
    else:
        d_str = ''
    return f'{d_str}{int(h):02.0f}:{int(m):02.0f}:{int(s):02.0f}'

def dict_to_str(dictionary: dict,
                level: int = 0,
                ) -> str:
    """
    A helper function to log dictionaries in a pretty way.

    Args:
        dictionary (dict): A general python dictionary.
        level (int): A recursion level counter, sets the visual indentation.

    Returns:
        str: A text representation for the dictionary.
    """
    message = ''
    for key, value in dictionary.items():
        if isinstance(value, dict):
            message += ' ' * level * 2 + str(key) + ':\n' + dict_to_str(value, level + 1)
        else:
            message += ' ' * level * 2 + str(key) + ': ' + str(value) + '\n'
    return message


def timedelta_from_str(time_str: str):
    """
    Get a datetime.timedelta object from its str() representation.
    Supports formats like "1hr", "30m", "15s", "1hr30m", "1hr30m15s".

    Args:
        time_str (str): The string representation of a time duration.

    Returns:
        datetime.timedelta: The corresponding timedelta object, or timedelta(0) if empty, None if invalid.
    """
    if not time_str:
        return datetime.timedelta(0)

    regex = re.compile(r'((?P<hours>\d+?)hr)?((?P<minutes>\d+?)m)?((?P<seconds>\d+?)s)?')
    parts = regex.match(time_str)
    
    if not parts or not any(parts.groupdict().values()):
        return None
    
    time_params = {}
    for (name, param) in parts.groupdict().items():
        if param:
            time_params[name] = int(param)
    
    return datetime.timedelta(**time_params)

    Args:
        time_str (str): The string representation of a datetime.timedelta object.

def convert_list_index_0_to_1(_list: Union[list, tuple], direction: int = 1) -> Union[list, tuple]:
    """
    Convert a list from 0-indexed to 1-indexed, or vice versa.
    Ensures positive values in the resulting list.

    Args:
        _list (list): The list to be converted.
        direction (int, optional): Either 1 or -1 to convert 0-indexed to 1-indexed or vice versa, respectively.

    Raises:
        ValueError: If the new list contains negative values.

    Returns:
        Union[list, tuple]: The converted indices.
    """
    new_list = [item + direction for item in _list]
    if any(val < 0 for val in new_list):
        raise ValueError(f'The resulting list from converting {_list} has negative values:\n{new_list}')
    if isinstance(_list, tuple):
        new_list = tuple(new_list)
    return new_list


# def bees_mol_to_dict_repr(mol: Molecule,
#                          reset_atom_ids: bool = False,
#                          testing: bool = False,
#                          ) -> dict:
#     """
#     based on the Generate a dict representation of an RMG ``Molecule`` object instance.

#     Args:
#         mol (Molecule): The RMG ``Molecule`` object instance.
#         reset_atom_ids (bool, optional): Whether to reset the atom IDs in the .mol Molecule attribute.
#                                          Useful when copying the object to avoid duplicate atom IDs between
#                                          different object instances.
#         testing (bool, optional): Whether this is called during a test, in which case atom IDs should be deterministic.

#     Returns:
#         dict: The corresponding dict representation.
#     """
#     mol = mol.copy(deep=True)
#     if testing:
#         counter = 0
#         for atom in mol.atoms:
#             atom.id = counter
#             counter += 1
#     elif len(mol.atoms) > 1 and mol.atoms[0].id == mol.atoms[1].id or reset_atom_ids:
#         mol.assign_atom_ids()
#     return {'atoms': [{'element': {'number': atom.element.number,
#                                    'isotope': atom.element.isotope,
#                                    },
#                        'radical_electrons': atom.radical_electrons,
#                        'charge': atom.charge,
#                        'label': atom.label,
#                        'lone_pairs': atom.lone_pairs,
#                        'id': atom.id,
#                        'props': atom.props,
#                        'atomtype': atom.atomtype.label,
#                        'edges': {atom_2.id: bond.order
#                                  for atom_2, bond in atom.edges.items()},
#                        } for atom in mol.atoms],
#             'multiplicity': mol.multiplicity,
#             'props': mol.props,
#             'atom_order': [atom.id for atom in mol.atoms]
#             }


# def rmg_mol_from_dict_repr(representation: dict,
#                            is_ts: bool = False,
#                            ) -> Optional[Molecule]:
#     """
#     Generate a dict representation of an RMG ``Molecule`` object instance.

#     Args:
#         representation (dict): A dict representation of an RMG ``Molecule`` object instance.
#         is_ts (bool, optional): Whether the ``Molecule`` represents a TS.

#     Returns:
#         ``Molecule``: The corresponding RMG ``Molecule`` object instance.

#     """
#     mol = Molecule(multiplicity=representation['multiplicity'],
#                    props=representation['props'])
#     atoms = {atom_dict['id']: Atom(element=get_element(value=atom_dict['element']['number'],
#                                                        isotope=atom_dict['element']['isotope']),
#                                    radical_electrons=atom_dict['radical_electrons'],
#                                    charge=atom_dict['charge'],
#                                    lone_pairs=atom_dict['lone_pairs'],
#                                    id=atom_dict['id'],
#                                    props=atom_dict['props'],
#                                    ) for atom_dict in representation['atoms']}
#     for atom_dict in representation['atoms']:
#         atoms[atom_dict['id']].atomtype = ATOMTYPES[atom_dict['atomtype']]
#     mol.atoms = list(atoms[atom_id] for atom_id in representation['atom_order'])
#     for i, atom_1 in enumerate(atoms.values()):
#         for atom_2_id, bond_order in representation['atoms'][i]['edges'].items():
#             bond = Bond(atom_1, atoms[atom_2_id], bond_order)
#             mol.add_bond(bond)
#     mol.update_atomtypes(raise_exception=False)
#     mol.update_multiplicity()
#     if not is_ts:
#         mol.identify_ring_membership()
#         mol.update_connectivity_values()
#     return mol


def generate_resonance_structures(object_: Union['Species', Any],
                                  keep_isomorphic: bool = False,
                                  filter_structures: bool = True,
                                  save_order: bool = True,
                                  ) -> Optional[List[Any]]:
    """
    Safely generate resonance structures for either an  BEES speicies (based on RMG Molecule or an RMG Species object instances).

    Args:
        object_ (Species, Molecule): The object to generate resonance structures for.
        keep_isomorphic (bool, optional): Whether to keep isomorphic isomers.
        filter_structures (bool, optional): Whether to filter resonance structures.
        save_order (bool, optional): Whether to make sure atom order is preserved.

    Returns:
        Optional[List[Molecule]]: If a ``Molecule`` object instance was given, the function returns a list of resonance
                                  structures (each is a ``Molecule`` object instance). If a ``Species`` object instance
                                  is given, the resonance structures are stored within the given object
                                  (in a .molecule attribute), and the function returns ``None``.
    """
    result = None
    try:
        # This function assumes 'object_' has a 'generate_resonance_structures' method.
        # This part of the code is commented out in common.py, so it's a placeholder.
        # result = object_.generate_resonance_structures(keep_isomorphic=keep_isomorphic,
        #                                                filter_structures=filter_structures,
        #                                                save_order=save_order,
        #                                                )
        pass # Placeholder for now
    except (TypeError, ValueError):
        pass
    return result


def calc_rmsd(x: Union[list, np.array],
              y: Union[list, np.array],
              ) -> float:
    """
    Compute the root-mean-square deviation between two matrices.

    Args:
        x (np.array): Matrix 1.
        y (np.array): Matrix 2.

    Returns:
        float: The RMSD score of two matrices.
    """
    x = np.array(x) if isinstance(x, list) else x
    y = np.array(y) if isinstance(y, list) else y
    d = x - y
    n = x.shape[0]
    sqr_sum = (d ** 2).sum()
    rmsd = np.sqrt(sqr_sum / n)
    return float(rmsd)


def safe_copy_file(source: str,
                   destination: str,
                   wait: int = 10,
                   max_cycles: int = 15,
                   ):
    """
    Copy a file safely.

    Args:
        source (str): The full path to the file to be copied.
        destination (str): The full path to the file destination.
        wait (int, optional): The number of seconds to wait between cycles.
        max_cycles (int, optional): The maximum number of cycles to try.
    """
    for i in range(max_cycles):
        try:
            shutil.copyfile(src=source, dst=destination)
        except shutil.SameFileError:
            break
        except NotADirectoryError:
            print(f'Expected "source" and "destination" to be directories. Trying to extract base paths.')
            if os.path.isfile(source):
                source = os.path.dirname(source)
            if os.path.isfile(destination):
                destination = os.path.dirname(destination)
        except OSError:
            time.sleep(wait)
        else:
            break
        if i >= max_cycles:
            break


# def dfs(mol: Molecule,
#         start: int,
#         sort_result: bool = True,
#         ) -> List[int]:
#     """
#     A depth-first search algorithm for graph traversal of a Molecule object instance.

#     Args:
#         mol (Molecule): The Molecule to search.
#         start (int): The index of the first atom in the search.
#         sort_result (bool, optional): Whether to sort the returned visited indices.

#     Returns:
#         List[int]: Indices of all atoms connected to the starting atom.
#     """
#     if start >= len(mol.atoms):
#         raise ValueError(f'Got a wrong start number {start} for a molecule with only {len(mol.atoms)} atoms.')
#     visited = list()
#     stack = deque()
#     stack.append(start)
#     while stack:
#         key = stack.pop()
#         if key in visited:
#             continue
#         visited.append(key)
#         for atom in mol.atoms[key].edges.keys():
#             stack.append(mol.atoms.index(atom))
#     visited = sorted(visited) if sort_result else visited
#     return visited



def convert_to_hours(time_str:str) -> float:
    """Convert walltime string in format HH:MM:SS to hours.

    Args:
        time_str (str): A time string in format HH:MM:SS
    
    Returns:
        float: The time in hours
    """
    h, m, s = map(int, time_str.split(':'))
    return h + m / 60 + s / 3600


def calculate_arrhenius_rate_coefficient(A: float, n: float, Ea: float, T: float, Ea_units: str = 'kJ/mol') -> float:
    """
    Calculate the Arrhenius rate coefficient.

    Args:
        A (float): Pre-exponential factor in cm^3, mol, s units.
        n (float): Temperature exponent.
        Ea (float): Activation energy in J/mol.
        T (float): Temperature in Kelvin.
        Ea_units (str): Units of the rate coefficient.

    Returns:
        float: The rate coefficient at the specified temperature.
    """
    if Ea_units not in EA_UNIT_CONVERSION:
        raise ValueError(f"Unsupported Ea units: {Ea_units}")
    return A * (T ** n) * math.exp(-1 * (Ea * EA_UNIT_CONVERSION[Ea_units]) / (R * T))

