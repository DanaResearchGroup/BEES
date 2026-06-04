"""
BEES common module

This module contains functions which are shared across multiple  modules.


"""

import logging
import os
import subprocess
from typing import List, Optional, Tuple, Union, Dict
import time
import yaml
import math
import re
import functools
from rdkit import Chem

 
# Absolute path to the BEES folder.

BEES_PATH = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))

# Define base paths for projects and data relative to BEES_PATH


PROJECTS_BASE_PATH = os.path.join(BEES_PATH, 'projects')
DATA_BASE_PATH = os.path.join(BEES_PATH, 'data')

VERSION = '0.1.0'  



# General cofactors (energy carriers, redox agents, vitamin-derived, metals).
# IMPORTANT: This set is used to identify cofactors that should NOT trigger new reaction discovery or ontology expansion
GENERAL_COFACTORS = {
    'atp', 'adp', 'amp', 'gtp', 'gdp', 'gmp', 'utp', 'udp', 'ump', 'ctp', 'cdp', 'cmp',
    'nadh', 'nad+', 'nad', 'nadh2', 'nadph', 'nadp+', 'nadp', 'nadph2',
    'fad', 'fadh2', 'coa', 'coenzyme a', 'coash',
    'h2o', 'water', 'h+', 'proton',
    'phosphate', 'pi', 'orthophosphate', 'pyrophosphate', 'ppi', 'diphosphate',
    'co2', 'carbon dioxide', 'hco3-', 'bicarbonate', 'o2', 'oxygen', 'h2', 'hydrogen',
    'tpp', 'thiamine pyrophosphate', 'thiamin pyrophosphate',
    'plp', 'pyridoxal phosphate', 'pyridoxal 5-phosphate',
    'thf', 'tetrahydrofolate', 'tetrahydrofolic acid', 'h4folate',
    'sam', 's-adenosylmethionine', 's-adenosyl-l-methionine',
    'lipoic acid', 'lipoamide', 'lipoyl',
    'cobalamin', 'vitamin b12', 'adenosylcobalamin', 'methylcobalamin',
    'mg2+', 'magnesium', 'mn2+', 'manganese', 'zn2+', 'zinc',
    'fe2+', 'fe3+', 'iron', 'cu2+', 'copper', 'ca2+', 'calcium',
}

# Subset of cofactors that are treated as "always available" for the purposes of reactant admission (i.e. they do not need to appear in the input species
# list). This intentionally excludes high‑energy carriers like ATP / NAD(H)/NADP(H);
# those must be provided explicitly or produced in the network.
COFACTORS_ALWAYS_AVAILABLE = {
    'h2o', 'water',
    'h+', 'proton',
    'phosphate', 'pi', 'orthophosphate', 'pyrophosphate', 'ppi', 'diphosphate',
    'co2', 'carbon dioxide', 'hco3-', 'bicarbonate',
    'o2', 'oxygen',
    'h2', 'hydrogen',
    'mg2+', 'magnesium',
    'mn2+', 'manganese',
    'zn2+', 'zinc',
    'fe2+', 'fe3+', 'iron',
    'cu2+', 'copper',
    'ca2+', 'calcium',
}

# EC number aliases
EC_ALIASES = {
    'EC 2.3.1.85': ['EC 2.3.1.86'],
    'EC 2.3.1.86': ['EC 2.3.1.85'],
}

# Enzyme domain cofactors
ENZYME_DOMAIN_COFACTORS = {
    'synthase': ['acp', 'acyl carrier protein'],
    'carboxylase': ['biotin'],
}

# Chemical ontology lookup table (populated by load_chemical_ontology at runtime)
_CHEMICAL_ONTOLOGY: Optional[Dict[str, List[str]]] = None
_ONTOLOGY_CATEGORIES_RAW: Optional[Dict[str, List[str]]] = None

# Acyl prefixes for ACP metabolites
_ACP_ACYL_PREFIXES = frozenset({
    "acetyl", "malonyl", "propionyl", "butyryl", "butanoyl", "acetoacetyl",
    "hexanoyl", "octanoyl", "decanoyl", "dodecanoyl", "tetradecanoyl",
    "hexadecanoyl", "octadecanoyl", "icosanoyl", "docosanoyl", "tetracosanoyl",
    "hexacosanoyl", "lauroyl", "myristoyl", "palmitoyl", "stearoyl",
    "cerotoyl", "lignoceroyl", "hydroxybutyryl", "crotonyl",
})

# Energy carriers for strict matching
ENERGY_CARRIERS = {
    "atp", "adp", "amp", "gtp", "gdp", "gmp", "itp", "idp", "imp",
    "utp", "udp", "ump", "ctp", "cdp", "cmp",
    "dttp", "dtdp", "dtmp", "datp", "dadp", "damp",
    "dgtp", "dgdp", "dgmp", "dctp", "dcdp", "dcmp"
}
# Constants
R = 8.31446261815324  # J/(mol*K)
EA_UNIT_CONVERSION = {'J/mol': 1, 'kJ/mol': 1e+3, 'cal/mol': 4.184, 'kcal/mol': 4.184e+3}


########################################################
#All the functions in the common module
########################################################

def canonical_smiles(smiles: Optional[str]) -> Optional[str]:
    """
    Return canonical SMILES for the given string, preserving stereochemistry.
    Returns None if input is empty; on parse/canonicalization failure returns
    the stripped original input so callers can still use it for exact matching.
    """
    s = str(smiles).strip() if smiles else ""
    if not s:
        return None
    try:
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            return s  # fallback: return original so callers can use for exact match
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception:
        return s  # fallback: return original


def heavy_atom_count(smiles: Optional[str]) -> Optional[int]:
    """
    Return the number of heavy (non-hydrogen) atoms in a SMILES string.
    Returns None if the input is empty, unparseable, or otherwise invalid.
    """
    s = str(smiles).strip() if smiles else ""
    if not s:
        return None
    try:
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            return None
        return mol.GetNumHeavyAtoms()
    except Exception:
        return None


def _normalize_compound_label(label: str) -> str:
    """Normalize compound label for comparison."""
    return str(label).lower().strip().replace("-", " ").replace("_", " ")

def _has_acyl_attachment_to_acp(label: str) -> bool:
    """True if label describes acyl-ACP (e.g. acetyl-ACP), not carrier-only."""
    normalized = _normalize_compound_label(label)
    for prefix in _ACP_ACYL_PREFIXES:
        if prefix in normalized:
            return True
    if "3 oxoacyl" in normalized or "3 oxo" in normalized and "acp" in normalized:
        return True
    if "3 hydroxyacyl" in normalized or "3 hydroxy" in normalized and "acp" in normalized:
        return True
    if "enoyl" in normalized or "dehydroacyl" in normalized:
        return True
    if normalized.startswith("acyl ") or " acyl " in normalized:
        return True
    return False

def get_coenzyme_like_flags(label: str) -> Dict[str, bool]:
    """Get flags indicating if compound is a coenzyme-like molecule."""
    normalized = _normalize_compound_label(label)
    is_free_coa = normalized in {"coa", "coenzyme a", "co enzyme a", "coenzyme a (coa)"}
    is_acp_like = ("acyl carrier protein" in normalized or "acp" in normalized.split())
    is_acyl_acp = is_acp_like and _has_acyl_attachment_to_acp(label)
    is_acp_carrier_only = is_acp_like and not is_acyl_acp
    is_biotin_like = ("biotin" in normalized or "carboxyl carrier protein" in normalized)
    is_coenzyme_like = is_free_coa or is_acp_carrier_only
    return {
        "normalized": normalized,
        "is_free_coa": is_free_coa,
        "is_acp_like": is_acp_like,
        "is_acyl_acp": is_acyl_acp,
        "is_acp_carrier_only": is_acp_carrier_only,
        "is_biotin_like": is_biotin_like,
        "is_coenzyme_like": is_coenzyme_like,
    }

def load_chemical_ontology(ontology_path: Optional[str] = None) -> Dict[str, List[str]]:
    """Load and invert chemical ontology from YAML (material -> [categories])."""
    global _CHEMICAL_ONTOLOGY
    if _CHEMICAL_ONTOLOGY is not None:
        return _CHEMICAL_ONTOLOGY
    if ontology_path is None:
        ontology_path = os.path.join(BEES_PATH, 'db', 'ontology.yaml')
    if not os.path.exists(ontology_path):
        _CHEMICAL_ONTOLOGY = {}
        return _CHEMICAL_ONTOLOGY
    try:
        categories = read_yaml_file(ontology_path)
    except Exception:
        _CHEMICAL_ONTOLOGY = {}
        return _CHEMICAL_ONTOLOGY
    inverted_ontology: Dict[str, List[str]] = {}
    if isinstance(categories, dict):
        for category, materials in categories.items():
            if isinstance(materials, list):
                for material in materials:
                    material_lower = material.lower().strip()
                    if material_lower not in inverted_ontology:
                        inverted_ontology[material_lower] = []
                    inverted_ontology[material_lower].append(category)
    _CHEMICAL_ONTOLOGY = inverted_ontology
    return _CHEMICAL_ONTOLOGY

def load_ontology_categories(ontology_path: Optional[str] = None) -> Dict[str, List[str]]:
    """Load raw ontology mapping (category -> [materials])."""
    global _ONTOLOGY_CATEGORIES_RAW
    if _ONTOLOGY_CATEGORIES_RAW is not None:
        return _ONTOLOGY_CATEGORIES_RAW
    if ontology_path is None:
        ontology_path = os.path.join(BEES_PATH, 'db', 'ontology.yaml')
    if not os.path.exists(ontology_path):
        _ONTOLOGY_CATEGORIES_RAW = {}
        return _ONTOLOGY_CATEGORIES_RAW
    try:
        categories = read_yaml_file(ontology_path)
        if not isinstance(categories, dict):
            _ONTOLOGY_CATEGORIES_RAW = {}
            return _ONTOLOGY_CATEGORIES_RAW
        normalized: Dict[str, List[str]] = {}
        for category, materials in categories.items():
            if not isinstance(materials, list):
                continue
            cat_lc = str(category).lower().strip()
            normalized[cat_lc] = [str(m).lower().strip() for m in materials]
        _ONTOLOGY_CATEGORIES_RAW = normalized
        return _ONTOLOGY_CATEGORIES_RAW
    except Exception:
        _ONTOLOGY_CATEGORIES_RAW = {}
        return _ONTOLOGY_CATEGORIES_RAW

@functools.lru_cache(maxsize=4096)
def get_ontology_equivalents(label: str) -> List[str]:
    """Return equivalent labels based on ontology categories.

    Results are cached (lru_cache) so repeated calls with the same label
    (e.g. during ODE integration) are free after the first lookup.
    """
    l_lc = str(label).lower().strip()
    out: List[str] = [l_lc]
    if l_lc in ENERGY_CARRIERS:
        return out
    out.extend(_acp_thioester_label_permutations(l_lc))
    aliases = get_chemical_aliases(l_lc)
    out.extend(aliases)
    categories_raw = load_ontology_categories()
    # IMPORTANT: only expand category members when the input label *is itself* a category.
    # This prevents concrete molecules (e.g. "hexanoate") from becoming equivalent to
    # other members in the same category (e.g. "octanoate") while still allowing
    # concrete molecules to match generic reaction templates via their category alias.
    if l_lc in categories_raw:
        out.extend(categories_raw[l_lc])
    return list(dict.fromkeys(out))


_ACP_SUFFIX_RE = re.compile(r"-\[(?:acp)\]$", re.IGNORECASE)


def _acp_thioester_label_permutations(label_lc: str) -> List[str]:
    """
    helper funaction for get_ontology_equivalents, for ACP-thioester labels.
    Generate narrowly-scoped equivalence permutations for ACP-thioester labels.

    Note: This is intentionally limited to labels ending in `-[ACP]` (case-insensitive) to
    avoid broad/unsafe ontology expansion. It targets common UniProt/Rhea naming
    re-orderings such as:

      - `3-oxo-(5z)-dodecenoyl-[acp]` <-> `(5z)-3-oxododecenoyl-[acp]`
      - `3-hydroxy-(5z)-dodecenoyl-[acp]` <-> `(5z)-3-hydroxydodecenoyl-[acp]`
    """
    if not label_lc:
        return []
    s = str(label_lc).lower().strip()
    if not _ACP_SUFFIX_RE.search(s):
        return []

    out: List[str] = []

    # Pattern A (forward): 3-oxo-(5z)-Xoyl-[acp]  ->  (5z)-3-oxoXoyl-[acp]
    #                       3-hydroxy-(5z)-Xoyl-[acp] -> (5z)-3-hydroxyXoyl-[acp]
    m = re.match(r"^(3-(?:oxo|hydroxy))-\(([^)]+)\)-(.+)$", s)
    if m:
        group = m.group(1)  # "3-oxo" or "3-hydroxy"
        stereo = m.group(2)  # e.g. "5z" or "3r,5z"
        tail = m.group(3)    # rest, typically "...oyl-[acp]"
        out.append(f"({stereo})-{group}{tail}")

    # Pattern A (reverse): (5z)-3-oxoXoyl-[acp]  ->  3-oxo-(5z)-Xoyl-[acp]
    m = re.match(r"^\(([^)]+)\)-(3-(?:oxo|hydroxy))(.+)$", s)
    if m:
        stereo = m.group(1)
        group = m.group(2)
        tail = m.group(3)
        out.append(f"{group}-({stereo})-{tail.lstrip('-')}")

    return list(dict.fromkeys(out))


def is_general_cofactor_label(label: str) -> bool:
    """
    Return True if label appears to be a general cofactor/carrier species.

    Uses `GENERAL_COFACTORS` plus ontology equivalents to catch naming variants.
    """
    normalized = " ".join(
        str(label).lower().strip().replace("_", " ").replace("-", " ").split()
    )
    compact = normalized.replace(" ", "")
    if normalized in GENERAL_COFACTORS or compact in GENERAL_COFACTORS:
        return True

    equivalents = get_ontology_equivalents(normalized)
    for eq in equivalents:
        eq_norm = " ".join(str(eq).lower().strip().split())
        if eq_norm in GENERAL_COFACTORS or eq_norm.replace(" ", "") in GENERAL_COFACTORS:
            return True
    return False

def get_chemical_aliases(molecule_label: str) -> List[str]:
    """Get all chemical class aliases for a molecule label."""
    molecule_lower = molecule_label.lower().strip()
    aliases = [molecule_lower]
    ontology = load_chemical_ontology()
    if molecule_lower in ontology:
        aliases.extend(ontology[molecule_lower])
    ontology_path = os.path.join(BEES_PATH, 'db', 'ontology.yaml')
    if os.path.exists(ontology_path):
        try:
            categories = read_yaml_file(ontology_path)
            if isinstance(categories, dict):
                for category in categories.keys():
                    if category.lower() == molecule_lower:
                        aliases.append(category)
                        break
        except Exception:
            # Non-fatal: aliases from the main ontology are still returned.
            logging.getLogger(__name__).debug("Failed to read ontology.yaml for category aliases", exc_info=True)
    return aliases

def log10_sd_to_linear_sd(linear_mean: float, sd_log10: float) -> float:
    """Convert standard deviation from log10 space to linear space."""
    ln10 = math.log(10)
    sigma_ln = ln10 * sd_log10
    exp_s2 = math.exp(sigma_ln * sigma_ln)
    var_linear = (exp_s2 - 1) * (linear_mean * linear_mean) * exp_s2
    return math.sqrt(max(0, var_linear))

def load_and_invert_ontology(path: str) -> Dict[str, List[str]]:
    """Read category-based YAML ontology and invert for lookups."""
    if not os.path.exists(path):
        return {}
    try:
        categories = read_yaml_file(path)
        if not categories or not isinstance(categories, dict):
            return {}
        lookup_table = {}
        for category, materials in categories.items():
            if not isinstance(materials, list):
                continue
            for material in materials:
                mat_lower = material.lower().strip()
                if mat_lower not in lookup_table:
                    lookup_table[mat_lower] = []
                cat_name = category.strip()
                if cat_name not in lookup_table[mat_lower]:
                    lookup_table[mat_lower].append(cat_name)
        return lookup_table
    except Exception:
        return {}




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
    """An exception class for reporting errors in BEES input files or parameters."""
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
    # Ensure project_directory has a trailing slash for consistent path construction and normalize path to handle different OS separators
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
    if 'calcs' in normalized_raw_path_value.lower() and 'species' in normalized_raw_path_value.lower():
        sub_path_marker = os.path.join('calcs', 'Species')

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
                return string 
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

   
    if project_directory is not None and not os.path.isabs(path):
        path = os.path.join(project_directory, path)


    if not os.path.isfile(path):
        raise InputError(f'Could not find the YAML file {original_path_for_error} (resolved to {path})')
    
    with open(path, 'r') as f:
        content = yaml.load(stream=f, Loader=yaml.FullLoader)
    
 
    return content 


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
    dir_path = os.path.dirname(path)
    if dir_path and not os.path.exists(dir_path):
        os.makedirs(dir_path, exist_ok=True)
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

