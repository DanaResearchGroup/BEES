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
                raise ValueError(f"Concentration range cannot have identical values. Got {value} for '{label}'")
            if value[0] > value[1]:
                value = (value[1], value[0])
            if value[0] < 0 or value[1] < 0:
                raise ValueError("Concentration cannot be negative")
        return value

    @field_validator('constant')
    def check_constant_species(cls, value, info: ValidationInfo):
        if value and isinstance(info.data.get('concentration'), tuple):
            raise ValueError("Constant species cannot have a concentration range.")
        return value

    @field_validator('reactive')
    def check_reactive(cls, value, info: ValidationInfo):
        if value and info.data.get('constant'):
            raise ValueError("Reactive species cannot be constant.")
        return value

    @field_validator('observable')
    def check_observable(cls, value, info: ValidationInfo):
        if value and info.data.get('constant'):
            raise ValueError("Observable species cannot be constant.")
        return value

    @field_validator('smiles')
    def check_smiles(cls, value):
        if value and not isinstance(value, str):
            raise ValueError("SMILES must be a string.")
        if value and Chem.MolFromSmiles(value) is None:
            raise ValueError(f"Invalid SMILES string: {value}")
        return value

    @field_validator('inchi')
    def validate_inchi(cls, value):
        if value is None:
            return value
        try:
            from rdkit.Chem import inchi
            mol = inchi.MolFromInchi(value)
            if mol is None:
                raise ValueError
        except Exception:
            raise ValueError(f"Invalid InChI string: {value}")
        return value


class Enzyme(BaseModel):
    label: str
    concentration: Union[confloat(gt=0), Tuple[confloat(gt=0), confloat(gt=0)]] = None
    ecnumber: Optional[str] = None
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
                raise ValueError(f"Concentration range cannot have identical values. Got {value} for '{label}'")
            if value[0] > value[1]:
                value = (value[1], value[0])
            if value[0] < 0 or value[1] < 0:
                raise ValueError("Concentration cannot be negative")
        return value
    
    @field_validator('label')
    def check_label(cls, value):
        if not isinstance(value, str):
            raise ValueError("Label must be a string.")
        if not value:
            raise ValueError("Label cannot be empty.")
        return value    
    



    @field_validator('constant')
    def check_constant_species(cls, value, info: ValidationInfo):
        if value and isinstance(info.data.get('concentration'), tuple):
            raise ValueError("Constant species cannot have a concentration range.")
        return value

    @field_validator('reactive')
    def check_reactive(cls, value, info: ValidationInfo):
        if value and info.data.get('constant'):
            raise ValueError("Reactive species cannot be constant.")
        return value

    @field_validator('observable')
    def check_observable(cls, value, info: ValidationInfo):
        if value and info.data.get('constant'):
            raise ValueError("Observable species cannot be constant.")
        return value

    @field_validator('ecnumber')
    def check_ecnumber(cls, value):
        if value is None:
            return value
        if not isinstance(value, str):
            raise ValueError("EC number must be a string.")
        if not value.startswith("EC"):
            raise ValueError(f"Invalid EC number: {value}. Must start with 'EC'.")
        return value



class SpeciesConstraints(BaseModel):
    allowed: List[str] = ['input species', 'seed mechanisms', 'reaction libraries']

    class Config:
        extra = "forbid"

    @field_validator('allowed')
    def check_allowed(cls, value):
        valid = ['input species', 'seed mechanisms', 'reaction libraries']
        if not value:
            raise ValueError("'allowed' list cannot be empty.")
        for val in value:
            if val not in valid:
                raise ValueError(f"Each 'allowed' value must be one of {valid}. Got: '{val}'")
        return value


class Database(BaseModel):
    name: str
    thermo_libraries: Optional[List[str]] = None
    kinetics_libraries: Optional[List[str]] = None
    use_low_credence_libraries: bool = False
    kinetics_depositories: Union[List[str], str] = 'default'
    kinetics_families: Union[str, List[str]] = 'default'
    parameter_estimator: Literal["ML", "group_contribution", "fixed_defaults"] = "ML"
    rate_law: Literal["Michaelis-Menten", "Hill", "MassAction"]
    parameters: Optional[Dict[str, float]] = None
    species_constraints: Optional[SpeciesConstraints] = None

    class Config:
        extra = "forbid"

    @field_validator('name')
    def check_name(cls, value):
        if not isinstance(value, str):
            raise ValueError("Name must be a string.")
        if not value:
            raise ValueError("Name cannot be empty.")
        return value    
    
    @field_validator('rate_law')
    def check_rate_law(cls, value):
        valid_rate_laws = ["Michaelis-Menten", "Hill", "MassAction"]
        if value not in valid_rate_laws:
            raise ValueError(f"Invalid rate law: {value}. Must be one of {valid_rate_laws}.")
        return value
    
    @field_validator('parameter_estimator')
    def check_parameter_estimator(cls, value):
        valid_estimators = ["ML", "group_contribution", "fixed_defaults"]
        if value not in valid_estimators:
            raise ValueError(f"Invalid parameter estimator: {value}. Must be one of {valid_estimators}.")
        return value
    
    @field_validator('kinetics_families')
    def check_kinetics_families(cls, value):
        if isinstance(value, str) and value != 'default':
            raise ValueError("kinetics_families must be a list or 'default'.")
        if isinstance(value, list):
            for family in value:
                if not isinstance(family, str):
                    raise ValueError(f"Each kinetics family must be a string. Got: {family}")
        return value
    
    @field_validator('kinetics_depositories')
    def check_kinetics_depositories(cls, value):
        if isinstance(value, str) and value != 'default':
            raise ValueError("kinetics_depositories must be a list or 'default'.")
        if isinstance(value, list):
            for depo in value:
                if not isinstance(depo, str):
                    raise ValueError(f"Each kinetics depository must be a string. Got: {depo}")
        return value
    
    @field_validator('thermo_libraries')
    def check_thermo_libraries(cls, value):
        if value is not None:
            if not isinstance(value, list):
                raise ValueError("thermo_libraries must be a list.")
            for lib in value:
                if not isinstance(lib, str):
                    raise ValueError(f"Each thermo library must be a string. Got: {lib}")
        return value
    
    @field_validator('kinetics_libraries')
    def check_kinetics_libraries(cls, value):
        if value is not None:
            if not isinstance(value, list):
                raise ValueError("kinetics_libraries must be a list.")
            for lib in value:
                if not isinstance(lib, str):
                    raise ValueError(f"Each kinetics library must be a string. Got: {lib}")
        return value

    @field_validator('species_constraints')
    def check_species_constraints(cls, value):
        if value is not None:
            if not isinstance(value, SpeciesConstraints):
                raise ValueError("species_constraints must be an instance of SpeciesConstraints.")
            if not value.allowed:
                raise ValueError("species_constraints must have at least one allowed entry.")
        return value
    


class RadicalTypeEnum(str, Enum):
    radical = 'radical'
    alkoxyl = 'alkoxyl'
    peroxyl = 'peroxyl'


class Environment(BaseModel):
    temperature: Union[confloat(gt=0), List[confloat(gt=0)]]
    pH: Union[confloat(gt=0), List[confloat(gt=0)]]
    ionic_strength: Optional[float] = None
    oxygen_level: Optional[float] = None
    seed_mechanisms: List[str] = list()

    class Config:          
        extra = "forbid"


    @field_validator('temperature')
    def check_t(cls, value):
        if isinstance(value, list) and len(value) != 2:
            raise ValueError(f'Temperature as list must have exactly 2 values (min, max). Got: {value}')
        return value

    @field_validator('pH')
    def validate_pH(cls, value):
        if isinstance(value, (int, float)) and not (0 <= value <= 14):
            raise ValueError(f'pH must be between 0 and 14. Got: {value}')
        return value

    @field_validator('ionic_strength')
    def validate_ionic_strength(cls, value):
        if value is not None and value < 0:
            raise ValueError(f'Ionic strength cannot be negative. Got: {value}')
        return value

    @field_validator('oxygen_level')
    def validate_oxygen_level(cls, value):
        if value is not None and not (0 <= value <= 1):
            raise ValueError(f'Oxygen level must be between 0 and 1. Got: {value}')
        return value


class Settings(BaseModel):
    end_time: confloat(gt=0)  
    time_step: confloat(gt=0)
    solver: Literal["odeint", "CVODE", "custom"]
    use_core_edge: bool = True
    threshold: confloat(gt=0) = 1e-3
    max_iterations: confloat(gt=0) = 50
    ml_models_enabled: bool = True
    stop_at_steady_state: bool = True
    termination_conversion: Optional[Dict[str, confloat(gt=0, lt=1)]] = None
    termination_rate_ratio: Optional[confloat(gt=0, lt=1)] = None

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


class InputBase(BaseModel):
    project: constr(max_length=255)
    project_directory: Optional[constr(max_length=255)] = None
    species: List[Species]
    enzymes: List[Enzyme]
    environment: Environment
    database: Database
    settings: Settings

    class Config:
        extra = "forbid"

