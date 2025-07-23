"""
BEES schema module
used for input validation
"""

from enum import Enum
from typing import Dict, List, Optional, Tuple, Union, Literal
from pydantic import BaseModel, conint, confloat, constr, field_validator, ValidationInfo
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
    concentration: Union[confloat(gt=0), Tuple[confloat(gt=0), confloat(gt=0)]] = None
    smiles: Optional[str] = None
    inchi: Optional[str] = None
    adjlist: Optional[str] = None
    charge: Optional[float] = 0
    constant: bool = False
    reactive: bool = True
    observable: bool = True
    solvent: bool = False
    xyz: Optional[Union[dict, str]] = None

    class Config:
        extra = "forbid"

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

    @field_validator('label')
    def check_label_not_empty(cls, value):
        if not value.strip():
            raise ValueError("Label cannot be empty")
        return value

    @field_validator('ecnumber')
    def validate_ecnumber(cls, value):
        if value and not value.startswith("EC "):
            raise ValueError("Invalid EC number: must start with 'EC '")
        return value


class Environment(BaseModel):
    temperature: Union[confloat(gt=0), Tuple[confloat(gt=0), confloat(gt=0)]]
    pH: Union[confloat(ge=0, le=14), Tuple[confloat(ge=0, le=14), confloat(ge=0, le=14)]]
    ionic_strength: Optional[confloat(ge=0)] = None
    oxygen_level: Optional[confloat(ge=0, le=1)] = None
    seed_mechanisms: Optional[List[str]] = None # pre-defined known or hypothesized reaction that we want to include from the beginning.

    class Config:
        extra = "forbid"

    @field_validator('temperature')
    def validate_temperature_range(cls, value):
        if isinstance(value, tuple) and len(value) != 2:
            raise ValueError("Temperature as list must have exactly 2 values (min, max)")
        if isinstance(value, tuple) and value[0] > value[1]:
            raise ValueError("Temperature range min value cannot be greater than max value")
        return value

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

    end_time: confloat(gt=0)
    time_step: confloat(gt=0)
    time_units : TerminationTimeEnum = TerminationTimeEnum.s # time units for the simulation, default is seconds

    #rate_law: Literal['Michaelis-Menten', 'Hill', 'MassAction'] # This line is commented out, so it's not part of the schema
    toleranceKeepInEdge: confloat(gt=0) = 0 #threshold for
    toleranceMoveToCore: confloat(gt=0) = 1e-5 # threshold for moving species from edge to core
    termination_conversion: Optional[Dict[str, confloat(gt=0, lt=1)]] = None # species: fraction of species that should be converted to termination species
    termination_rate_ratio: Optional[confloat(gt=0, lt=1)] = None # ADDED THIS LINE: termination_rate_ratio field
    max_edge_species: Optional[conint(gt=0)] = None # maximum number of species in the edge
    filter_reactions: bool = True # whether to filter reactions based on species constraints
    modify_concentration_ranges_together: bool = True

    max_iterations: conint(gt=0) = 50
    verbose: Optional[conint(ge=10, le=50)] = 20
    saveEdgeSpecies: bool = True # Corrected syntax: removed '=' before bool and trailing comma
    output_directory: Optional[str] = None
    generate_plots: bool = False
    save_simulation_profiles: bool = False

    class Config:
        extra = "forbid"

    @field_validator('time_step')
    def validate_time_step(cls, value, info: ValidationInfo):
        end_time = info.data.get('end_time')
        if end_time is not None and value >= end_time:
            raise ValueError(f"'time_step' must be smaller than 'end_time' ({end_time}). Got: {value}")
        return value

    @field_validator('termination_conversion')
    def validate_termination_conversion(cls, value):
        if value:
            for species, frac in value.items():
                if not (0 < frac < 1):
                    raise ValueError(f"termination_conversion values must be between 0 and 1. Got: {species}: {frac}")
        return value

    @field_validator('termination_rate_ratio')
    def validate_rate_ratio(cls, value):
        if value and not (0 < value < 1):
            raise ValueError("termination_rate_ratio must be between 0 and 1 (exclusive).")
        return value

    @field_validator('verbose')
    def validate_verbose_level(cls, value):
        if value is not None and value not in [10, 20, 30, 40, 50]:
            raise ValueError("Verbose level must be 10, 20, 30, 40, or 50")
        return value


class SpeciesConstraints(BaseModel):
    allowed: List[Literal['input species', 'seed mechanisms', 'reaction libraries']] = ['input species', 'seed mechanisms', 'reaction libraries']
    # Removed 'termination_conversion' from here to avoid duplication, it's in Settings now.
    tolerance_thermo_keep_species_in_edge: Optional[confloat(gt=0)] = None # Tolerance for thermodynamic properties to keep species in the edge (e.g., max deviation from a threshold)
    max_C_atoms: Optional[conint(gt=0)] = None # Maximum number of Carbon atoms allowed in generated species
    max_O_atoms: Optional[conint(gt=0)] = None # Maximum number of Oxygen atoms allowed in generated species
    max_N_atoms: Optional[conint(gt=0)] = None # Maximum number of Nitrogen atoms allowed in generated species
    max_Si_atoms: Optional[conint(gt=0)] = None # Maximum number of Silicon atoms allowed in generated species
    max_S_atoms: Optional[conint(gt=0)] = None # Maximum number of Sulfur atoms allowed in generated species
    max_heavy_atoms: Optional[conint(gt=0)] = None # Maximum number of heavy (non-hydrogen) atoms allowed in generated species
    max_radical_electrons: Optional[conint(ge=0)] = None # Maximum number of unpaired electrons (radicals) allowed in generated species
    max_singlet_carbenes: Optional[conint(ge=0)] = 1 # Maximum number of singlet carbenes allowed in generated species
    max_carbene_radicals: Optional[conint(ge=0)] = 0 # Maximum number of carbene radicals allowed in generated species
    allow_singlet_O2: bool = True # Whether to allow singlet oxygen (O2(a1Δg)) as an input species

    class Config:
        extra = "forbid"

    @field_validator('allowed')
    def check_allowed_not_empty(cls, value):
        if not value:
            raise ValueError("'allowed' list cannot be empty")
        return value


class Database(BaseModel):
    name: constr(min_length=1)
    parameters: Optional[Dict[str, float]] = None
    thermo_libraries: Optional[List[str]] = None
    kinetics_libraries: Optional[List[str]] = None
    chemistry_sets: Optional[List[str]] = None # Pre-defined collections of species and reactions to include
    use_low_credence_libraries: bool = False
    seed_mechanism: Optional[List[str]] = None # Corrected syntax: removed '=' before Optional
    kinetics_depositories: Union[List[str], Literal['default']] = 'default' # Pre-defined kinetics depositories to use. 'default' will be a list of common depositories.
    kinetics_families: Union[List[str], Literal['default']] = 'default' # Pre-defined kinetics families to use. 'default' will be a list of common families.
    solver: Literal['odeint', 'CVODE', 'BDF'] = 'odeint' # Name of the ODE solver to use for kinetic simulations.
    kinetics_estimator: str = 'rate rules' # Name of the kinetics estimator to use, default is 'rate rules', but more might be added in the future.

    class Config:
        extra = "forbid"

    @field_validator('name')
    def check_name_not_empty(cls, value):
        if not value.strip():
            raise ValueError("Name cannot be empty")
        return value

    @field_validator('thermo_libraries')
    def validate_thermo_libraries(cls, value):
        if value:
            if not isinstance(value, list):
                raise ValueError("thermo_libraries must be a list")
            for item in value:
                if not isinstance(item, str):
                    raise ValueError("Each thermo library must be a string")
        return value

    @field_validator('kinetics_libraries')
    def validate_kinetics_libraries(cls, value):
        if value:
            if not isinstance(value, list):
                raise ValueError("kinetics_libraries must be a list")
            for item in value:
                if not isinstance(item, str):
                    raise ValueError("Each kinetics library must be a string")
        return value

    @field_validator('kinetics_depositories')
    def validate_kinetics_depositories(cls, value):
        if isinstance(value, list):
            for item in value:
                if not isinstance(item, str):
                    raise ValueError("Each kinetics depository must be a string")
        return value

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

    @field_validator('project')
    def check_project_not_empty(cls, value):
        if not value.strip():
            raise ValueError("Project name cannot be empty")
        return value

    @field_validator('species')
    def check_species_list_not_empty(cls, value):
        if not value:
            raise ValueError("Species list cannot be empty")
        return value

    @field_validator('enzymes')
    def check_enzymes_list_not_empty(cls, value):
        if not value:
            raise ValueError("Enzymes list cannot be empty")
        return value
