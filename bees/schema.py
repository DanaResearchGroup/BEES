"""
BEES schema module
used for input validation

This module defines the schema for BEES input validation using Pydantic.

"""

from enum import Enum
from typing import Dict, List, Optional, Tuple, Union, Literal, Annotated
from pydantic import BaseModel, conint, confloat, constr, field_validator, ValidationInfo, Field
from rdkit import Chem


class TerminationTimeEnum(str, Enum):
    micro_s = 'micro-s'
    ms = 'ms'
    s = 's'
    hrs = 'hrs'
    minutes = 'min'
    hours = 'hours'
    days = 'days'


class Species(BaseModel):
    """
    A class for validate input.BEES.Species arguments.
    """
    label: str
    concentration: Union[Annotated[float, Field(gt=0)], Tuple[Annotated[float, Field(gt=0)], Annotated[float, Field(gt=0)]]] = None
    smiles: Optional[str] = None
    inchi: Optional[str] = None
    adjlist: Optional[str] = None
    charge: Optional[float] = 0
    constant: bool = False
    reactive: bool = True
    observable: bool = True
    solvent: bool = False
    xyz: Optional[Union[dict, str]] = None  # PLACEHOLDER: XYZ coordinates for future structural analysis (inherited from RMG, not currently used)

    class Config:
        extra = "forbid"

    @classmethod
    @field_validator('concentration')
    def check_concentration_range_order(cls, value, info: ValidationInfo):
        label = info.data.get('label')
        if value is None:
            raise ValueError(f"Concentration must be specified for '{label}'")
        if isinstance(value, float):
            if value < 0:
                raise ValueError(f"Concentration cannot be negative. Got {value} for '{label}'")
            return value
        if isinstance(value, tuple):
            if value[0] == value[1]:
                raise ValueError("Concentration range cannot have identical values")
            if value[0] < 0 or value[1] < 0:
                raise ValueError(f"Concentration cannot be negative. Got {value} for '{label}'")
            if value[0] > value[1]:
                raise ValueError(f"Concentration range min value ({value[0]}) cannot be greater than max value ({value[1]}) for '{label}'")
            return value
        return value

    @classmethod
    @field_validator('constant')
    def check_constant_species(cls, value, info: ValidationInfo):
        if value:
            if info.data.get('concentration') and isinstance(info.data['concentration'], tuple):
                raise ValueError("Constant species cannot have a concentration range")
            if info.data.get('reactive'):
                raise ValueError("Reactive species cannot be constant")
            if info.data.get('observable'):
                raise ValueError("Observable species cannot be constant")
        return value

    @classmethod
    @field_validator('smiles')
    def validate_smiles(cls, value):
        if value:
            try:
                mol = Chem.MolFromSmiles(value)
                if not mol:
                    raise ValueError("Invalid SMILES string")
            except Exception:
                raise ValueError("Invalid SMILES string")
        return value

    @classmethod
    @field_validator('inchi')
    def validate_inchi(cls, value):
        if value:
            try:
                mol = Chem.MolFromInchi(value)
                if not mol:
                    raise ValueError("Invalid InChI string")
            except Exception:
                raise ValueError("Invalid InChI string")
        return value

    @classmethod
    @field_validator('adjlist')
    def validate_adjlist(cls, value):
        if value:
            try:
                # Assuming adjlist is in a format parsable by MolFromMolBlock (e.g., MDL Molfile)
                mol = Chem.MolFromMolBlock(value)
                if not mol:
                    raise ValueError("Invalid adjacency list (MolBlock format expected)")
            except Exception:
                raise ValueError("Invalid adjacency list (MolBlock format expected)")
        return value


class Enzyme(Species):
    """
    A class for validate input.BEES.Enzyme arguments if there are any.
    Inherits from Species, adding specific fields for enzymes.
    """
    ecnumber: Optional[constr(pattern=r"^EC \d+\.\d+\.\d+\.\d+$")] = None
    amino_acid_sequence: Optional[str] = None

    @classmethod
    @field_validator('label')
    def check_label_not_empty(cls, value):
        if not value.strip():
            raise ValueError("Label cannot be empty")
        return value

    @classmethod
    @field_validator('ecnumber')
    def validate_ecnumber(cls, value):
        if value and not value.startswith("EC "):
            raise ValueError("Invalid EC number: must start with 'EC '")
        return value

    @classmethod
    @field_validator('amino_acid_sequence')
    def validate_amino_acid_sequence(cls, value):
        if value is None:
            return value
        
        # Check for spaces
        if ' ' in value:
            raise ValueError("Amino acid sequence cannot contain spaces")
        
        # Check if all characters are uppercase
        if not value.isupper():
            raise ValueError("Amino acid sequence must be in capital letters")
        
        # Valid amino acid single-letter codes
        valid_aa_codes = set('ACDEFGHIKLMNPQRSTVWY')
        
        # Check if all characters are valid amino acid codes
        invalid_chars = set(value) - valid_aa_codes
        if invalid_chars:
            raise ValueError(f"Amino acid sequence contains invalid characters: {', '.join(sorted(invalid_chars))}. Valid codes are: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y")
        
        return value


class Environment(BaseModel):
    temperature: Union[confloat(gt=0), Tuple[confloat(gt=0), confloat(gt=0)]]  # Currently used: Single value or range
    pH: Union[confloat(ge=0, le=14), Tuple[confloat(ge=0, le=14), confloat(ge=0, le=14)]]  # Currently used: Single value or range
    # PLACEHOLDERS for future functionality - not currently used
    ionic_strength: Optional[confloat(ge=0)] = None  # PLACEHOLDER: Ionic strength for future calculations
    oxygen_level: Optional[confloat(ge=0, le=1)] = None  # PLACEHOLDER: Oxygen level for future aerobic/anaerobic conditions
    seed_mechanisms: Optional[List[str]] = None  # PLACEHOLDER: Pre-defined reaction mechanisms for future use

    class Config:
        extra = "forbid"

    @classmethod
    @field_validator('temperature')
    def validate_temperature_range(cls, value):
        if isinstance(value, tuple) and len(value) != 2:
            raise ValueError("Temperature as list must have exactly 2 values (min, max)")
        if isinstance(value, tuple) and value[0] > value[1]:
            raise ValueError("Temperature range min value cannot be greater than max value")
        return value

    @classmethod
    @field_validator('pH')
    def validate_pH_range(cls, value):
        if isinstance(value, tuple) and len(value) != 2:
            raise ValueError("pH as list must have exactly 2 values (min, max)")
        if isinstance(value, tuple) and value[0] > value[1]:
            raise ValueError("pH range min value cannot be greater than max value")
        return value


class Settings(BaseModel):
    """
    A class for validate input.BEES.Settings arguments.
    """

    end_time: confloat(gt=0)  # PLACEHOLDER: Simulation end time (ODE simulation not yet implemented)
    time_step: confloat(gt=0)  # PLACEHOLDER: Simulation time step (ODE simulation not yet implemented)
    time_units : TerminationTimeEnum = TerminationTimeEnum.s  # PLACEHOLDER: Time units for simulation (ODE simulation not yet implemented)

    # PLACEHOLDERS for future iterative network refinement functionality (not yet implemented)
    toleranceKeepInEdge: confloat(gt=0) = 0  # PLACEHOLDER: Threshold for keeping species in edge during iterative refinement
    toleranceMoveToCore: confloat(gt=0) = 1e-5  # PLACEHOLDER: Threshold for moving species to core during iterative refinement
    termination_conversion: Optional[Dict[str, confloat(gt=0, lt=1)]] = None  # PLACEHOLDER: Species conversion termination criteria
    termination_rate_ratio: Optional[confloat(gt=0, lt=1)] = None  # PLACEHOLDER: Rate ratio termination criteria
    max_edge_species: Optional[conint(gt=0)] = None  # PLACEHOLDER: Maximum edge species during iterative refinement
    filter_reactions: bool = True  # PLACEHOLDER: Reaction filtering during iterative refinement
    modify_concentration_ranges_together: bool = True  # PLACEHOLDER: Concentration range modification behavior

    max_iterations: conint(gt=0) = 50  # PLACEHOLDER: Maximum iterations for iterative refinement (not yet implemented)
    verbose: Optional[conint(ge=10, le=50)] = 20  # Currently used: Logging verbosity level
    saveEdgeSpecies: bool = True  # PLACEHOLDER: Whether to save edge species during iterative refinement
    output_directory: Optional[str] = None  # Currently used: Custom output directory path
    generate_plots: bool = False  # PLACEHOLDER: Plot generation functionality (not yet implemented)
    save_simulation_profiles: bool = False  # PLACEHOLDER: Simulation profile saving (ODE simulation not yet implemented)

    class Config:
        extra = "forbid"

    @classmethod
    @field_validator('time_step')
    def validate_time_step(cls, value, info: ValidationInfo):
        end_time = info.data.get('end_time')
        if end_time is not None and value >= end_time:
            raise ValueError(f"'time_step' must be smaller than 'end_time' ({end_time}). Got: {value}")
        return value

    @classmethod
    @field_validator('termination_conversion')
    def validate_termination_conversion(cls, value):
        if value:
            for species, frac in value.items():
                if not (0 < frac < 1):
                    raise ValueError(f"termination_conversion values must be between 0 and 1. Got: {species}: {frac}")
        return value

    @classmethod
    @field_validator('termination_rate_ratio')
    def validate_rate_ratio(cls, value):
        if value and not (0 < value < 1):
            raise ValueError("termination_rate_ratio must be between 0 and 1 (exclusive).")
        return value

    @classmethod
    @field_validator('verbose')
    def validate_verbose_level(cls, value):
        if value is not None and value not in [10, 20, 30, 40, 50]:
            raise ValueError("Verbose level must be 10, 20, 30, 40, or 50")
        return value


class SpeciesConstraints(BaseModel):
    """
    PLACEHOLDER: Species constraints for future iterative network refinement functionality.
    This class is defined but not currently used in the InputBase schema.
    All fields are placeholders for future functionality.
    """
    allowed: List[Literal['input species', 'seed mechanisms', 'reaction libraries']] = ['input species', 'seed mechanisms', 'reaction libraries']  # PLACEHOLDER
    tolerance_thermo_keep_species_in_edge: Optional[confloat(gt=0)] = None  # PLACEHOLDER: Tolerance for thermodynamic properties to keep species in the edge
    max_C_atoms: Optional[conint(gt=0)] = None  # PLACEHOLDER: Maximum number of Carbon atoms allowed in generated species
    max_O_atoms: Optional[conint(gt=0)] = None  # PLACEHOLDER: Maximum number of Oxygen atoms allowed in generated species
    max_N_atoms: Optional[conint(gt=0)] = None  # PLACEHOLDER: Maximum number of Nitrogen atoms allowed in generated species
    max_Si_atoms: Optional[conint(gt=0)] = None  # PLACEHOLDER: Maximum number of Silicon atoms allowed in generated species
    max_S_atoms: Optional[conint(gt=0)] = None  # PLACEHOLDER: Maximum number of Sulfur atoms allowed in generated species
    max_heavy_atoms: Optional[conint(gt=0)] = None  # PLACEHOLDER: Maximum number of heavy (non-hydrogen) atoms allowed in generated species
    max_radical_electrons: Optional[conint(ge=0)] = None  # PLACEHOLDER: Maximum number of unpaired electrons (radicals) allowed in generated species
    max_singlet_carbenes: Optional[conint(ge=0)] = 1  # PLACEHOLDER: Maximum number of singlet carbenes allowed in generated species
    max_carbene_radicals: Optional[conint(ge=0)] = 0  # PLACEHOLDER: Maximum number of carbene radicals allowed in generated species
    allow_singlet_O2: bool = True  # PLACEHOLDER: Whether to allow singlet oxygen (O2(a1Δg)) as an input species

    class Config:
        extra = "forbid"

    @classmethod
    @field_validator('allowed')
    def check_allowed_not_empty(cls, value):
        if not value:
            raise ValueError("'allowed' list cannot be empty")
        return value


class Database(BaseModel):
    name: constr(min_length=1)
    parameters: Optional[Dict[str, float]] = None
    # Placeholder for future functionality - not currently used
    thermo_libraries: Optional[List[str]] = None  # PLACEHOLDER: Thermodynamic libraries for future use
    kinetics_libraries: Optional[List[str]] = None  # PLACEHOLDER: Kinetics libraries for future use
    chemistry_sets: Optional[List[str]] = None  # PLACEHOLDER: Pre-defined collections of species and reactions for future use
    use_low_credence_libraries: bool = False  # PLACEHOLDER: For future library filtering functionality
    seed_mechanism: Optional[List[str]] = None  # PLACEHOLDER: Pre-defined reaction mechanisms for future use
    kinetics_depositories: Union[List[str], Literal['default']] = 'default'  # PLACEHOLDER: Kinetics data sources for future use
    kinetics_families: Union[List[str], Literal['default']] = 'default'  # PLACEHOLDER: Reaction families for future use
    # Currently used fields
    solver: Literal['odeint', 'CVODE', 'BDF'] = 'odeint'  # PLACEHOLDER: ODE solver name (ODE simulation not yet implemented)
    kinetics_estimator: str = 'rate rules'  # PLACEHOLDER: Kinetics estimation method name (not yet fully implemented)

    class Config:
        extra = "forbid"

    @classmethod
    @field_validator('name')
    def check_name_not_empty(cls, value):
        if not value.strip():
            raise ValueError("Name cannot be empty")
        return value

    @classmethod
    @field_validator('thermo_libraries')
    def validate_thermo_libraries(cls, value):
        if value:
            if not isinstance(value, list):
                raise ValueError("thermo_libraries must be a list")
            for item in value:
                if not isinstance(item, str):
                    raise ValueError("Each thermo library must be a string")
        return value

    @classmethod
    @field_validator('kinetics_libraries')
    def validate_kinetics_libraries(cls, value):
        if value:
            if not isinstance(value, list):
                raise ValueError("kinetics_libraries must be a list")
            for item in value:
                if not isinstance(item, str):
                    raise ValueError("Each kinetics library must be a string")
        return value

    @classmethod
    @field_validator('kinetics_depositories')
    def validate_kinetics_depositories(cls, value):
        if isinstance(value, list):
            for item in value:
                if not isinstance(item, str):
                    raise ValueError("Each kinetics depository must be a string")
        return value

    @classmethod
    @field_validator('kinetics_families')
    def validate_kinetics_families(cls, value):
        if isinstance(value, list):
            for item in value:
                if not isinstance(item, str):
                    raise ValueError("Each kinetics family must be a string")
        return value


class InputBase(BaseModel):
    project: constr(max_length=255)
    project_directory: Optional[constr(max_length=255)] = None
    verbose: Optional[conint(ge=10, le=50)] = 20 # Moved to Settings
    species: List[Species]
    enzymes: List[Enzyme]
    environment: Environment
    settings: Settings
    database: Database

    class Config:
        extra = "forbid"

    @classmethod
    @field_validator('project')
    def check_project_not_empty(cls, value):
        if not value.strip():
            raise ValueError("Project name cannot be empty")
        return value

    @classmethod
    @field_validator('species')
    def check_species_list_not_empty(cls, value):
        if not value:
            raise ValueError("Species list cannot be empty")
        return value

    @classmethod
    @field_validator('enzymes')
    def check_enzymes_list_not_empty(cls, value):
        if not value:
            raise ValueError("Enzymes list cannot be empty")
        return value
