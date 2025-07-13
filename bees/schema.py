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
    label: str
    concentration: Union[confloat(gt=0), Tuple[confloat(gt=0), confloat(gt=0)]] = None
    smiles: Optional[str] = None
    inchi: Optional[str] = None
    charge: Optional[float] = 0
    constant: bool = False
    reactive: bool = True
    observable: bool = False
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


class Enzyme(Species):
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
    seed_mechanisms: Optional[List[str]] = None

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
    end_time: confloat(gt=0)
    time_step: confloat(gt=0)
    solver: Literal['odeint', 'CVODE', 'BDF'] = 'odeint'
    use_core_edge: bool = True
    threshold: confloat(gt=0) = 1e-3
    max_iterations: conint(gt=0) = 50
    ml_models_enabled: bool = True
    stop_at_steady_state: bool = True
    termination_conversion: Optional[Dict[str, confloat(gt=0, lt=1)]] = None
    termination_rate_ratio: Optional[confloat(gt=0, lt=1)] = None
    verbose: Optional[conint(ge=10, le=50)] = 20 # Added this field
    output_directory: Optional[str] = None # Added this field
    flux_adapter: Literal['RMG', 'Cantera', 'Custom'] = 'RMG' # Added this field
    profiles_adapter: Literal['RMG', 'Cantera', 'Custom'] = 'RMG' # Added this field
    generate_plots: bool = False # Added this field
    save_simulation_profiles: bool = False # Added this field

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

    @field_validator('verbose') # Re-added this custom validator
    def validate_verbose_level(cls, value):
        if value is not None and value not in [10, 20, 30, 40, 50]:
            raise ValueError("Verbose level must be 10, 20, 30, 40, or 50")
        return value


class SpeciesConstraints(BaseModel):
    allowed: List[Literal['input species', 'seed mechanisms', 'reaction libraries']] = ['input species', 'seed mechanisms', 'reaction libraries']
    max_C_atoms: Optional[conint(gt=0)] = None
    max_O_atoms: Optional[conint(gt=0)] = None
    max_N_atoms: Optional[conint(gt=0)] = None
    max_Si_atoms: Optional[conint(gt=0)] = None
    max_S_atoms: Optional[conint(gt=0)] = None
    max_heavy_atoms: Optional[conint(gt=0)] = None
    max_radical_electrons: Optional[conint(ge=0)] = None
    max_singlet_carbenes: Optional[conint(ge=0)] = 1
    max_carbene_radicals: Optional[conint(ge=0)] = 0
    allow_singlet_O2: bool = True

    class Config:
        extra = "forbid"

    @field_validator('allowed')
    def check_allowed_not_empty(cls, value):
        if not value:
            raise ValueError("'allowed' list cannot be empty")
        return value


class Database(BaseModel):
    name: constr(min_length=1)
    rate_law: Literal['Michaelis-Menten', 'Hill', 'MassAction']
    parameter_estimator: Literal['ML', 'group_contribution', 'fixed_defaults'] = 'ML'
    parameters: Optional[Dict[str, float]] = None
    thermo_libraries: Optional[List[str]] = None
    kinetics_libraries: Optional[List[str]] = None
    use_low_credence_libraries: bool = False
    kinetics_depositories: Union[List[str], Literal['default']] = 'default'
    kinetics_families: Union[List[str], Literal['default']] = 'default'
    species_constraints: Optional[SpeciesConstraints] = None

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

    @field_validator('species_constraints')
    def validate_species_constraints_type(cls, value):
        if value is not None and not isinstance(value, SpeciesConstraints):
            raise ValueError("species_constraints must be an instance of SpeciesConstraints")
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