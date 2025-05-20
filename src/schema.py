"""
BEES schema module
used for input validation
"""

from enum import Enum
from typing import Dict, List, Optional, Tuple, Union, Literal 
import pydantic 
from pydantic import BaseModel, conint, confloat, constr, model_validator, field_validator
from rdkit import Chem


class TerminationTimeEnum(str, Enum):
    """
    The supported termination type units in an "reactor".
    """
    micro_s = 'micro-s'
    ms = 'ms'
    s = 's'
    hrs = 'hrs'
    minutes = 'min'
    hours = 'hours'
    days = 'days'


class Species(BaseModel):
    """
    A class for validating input.species arguments
    """
    label: str
    concentration: Union[confloat(gt=0), Tuple[confloat(gt=0), confloat(gt=0)]] = None # concentration in mol/L
    smiles: Optional[str] = None
    inchi: Optional[str] = None
    charge: Optional[float] = 0
    constant: bool = False
    reactive: bool = True
    observable: bool = False  # Whether this is an observable for SA
    solvent: bool = False # True if the compound is a solvent
    xyz: Optional[Union[dict, str]] = None

    class Config:
        extra = "forbid"

    @field_validator('concentration')
    def check_concentration_range_order(cls, value, values):
        """Ensure range is valid and ordered"""
        label = values.get('label')
        if value is None:   
            raise ValueError(f"Concentration must be specified for '{label}'")
        if isinstance(value, float):       
            if value < 0:
                raise ValueError(f"Concentration cannot be negative. Got {value} for '{label}'")
            return value
        # If it's a tuple, check the range                       
        if isinstance(value, tuple):
            if value[0] == value[1]:
                raise ValueError(f"Concentration range cannot have identical values. Got {value} for '{label}'")
            if value[0] > value[1]:
                value = (value[1], value[0])
            if value[0] < 0 or value[1] < 0:
                raise ValueError(f"Concentration cannot be negative")
        return value

    @field_validator('constant')
    def check_constant_species(cls, value, values):
        if value and isinstance(values.get('concentration'), tuple):
            raise ValueError("Constant species cannot have a concentration range.")
        return value
    
    @field_validator('reactive')
    def check_reactive(cls, value, values):
        """Check if the compound is reactive"""
        if value and values.get('constant'):
            raise ValueError("Reactive species cannot be constant.")
        return value

    @field_validator('observable')
    def check_observable(cls, value, values):
        """Check if the compound is observable"""
        if value and values.get('constant'):
            raise ValueError("Observable species cannot be constant.")
        return value

    @field_validator('smiles')
    def check_smiles(cls, value):
        """Check if the smiles is valid"""
        if value and not isinstance(value, str):
            raise ValueError("SMILES must be a string.")
        if Chem.MolFromSmiles(value) is None:
            raise ValueError(f"Invalid SMILES string: {value}")
        return value

    @field_validator('inchi')
    def validate_inchi(cls, value):
        """Check that InChI string is chemically valid"""
        if value is +None:
            return value
        try:
            from rdkit.Chem import inchi
            mol = inchi.MolFromInchi(value)
            if mol is None:
                raise ValueError
        except Exception:
            raise ValueError(f"Invalid InChI string: {value}")
        return value


class Enzyme(BaseModel): # TODO: add enzyme specific fields
    """
    A class for validating input/species arguments
    """
    label: str
    concentration: Union[confloat(gt=0), Tuple[confloat(gt=0), confloat(gt=0)]] = None # concentration in mol/L
    ecnumber: Optional[str] = None
    charge: Optional[float] = 0
    constant: bool = False
    reactive: bool = True
    observable: bool = False  # Whether this is an observable for SA
    solvent: bool = False # True if the compound is a solvent
    xyz: Optional[Union[dict, str]] = None

    class Config:
        extra = "forbid"

    @field_validator('concentration')
    def check_concentration_range_order(cls, value, values):
        """Ensure range is valid and ordered"""
        label = values.get('label')
        if value is None:
            raise ValueError(f"Concentration must be specified for '{label}'")
        if isinstance(value, float):
            if value < 0:
                raise ValueError(f"Concentration cannot be negative. Got {value} for '{label}'")
            return value
        # If it's a tuple, check the range
        if isinstance(value, tuple):
            if value[0] == value[1]:
                raise ValueError(f"Concentration range cannot have identical values. Got {value} for '{label}'")
            if value[0] > value[1]:
                value = (value[1], value[0])
            if value[0] < 0 or value[1] < 0:
                raise ValueError(f"Concentration cannot be negative")
        return value

    @field_validator('constant')
    def check_constant_species(cls, value, values):
        if value and isinstance(values.get('concentration'), tuple):
            raise ValueError("Constant species cannot have a concentration range.")
        return value

    @field_validator('reactive')
    def check_reactive(cls, value, values):
        """Check if the compound is reactive"""
        if value and values.get('constant'):
            raise ValueError("Reactive species cannot be constant.")
        return value

    @field_validator('observable')
    def check_observable(cls, value, values):
        """Check if the compound is observable"""
        if value and values.get('constant'):
            raise ValueError("Observable species cannot be constant.")
        return value
    
    @field_validator('ecnumber')
    def check_ecnumber(cls, value):
        """Check if the EC number is valid"""
        if value is None:
            return value
        if not isinstance(value, str):
            raise ValueError("EC number must be a string.")
        if not value.startswith("EC"):
            raise ValueError(f"Invalid EC number: {value}. Must start with 'EC'.")

    

    


class SpeciesConstraints(BaseModel):
    """
    A class for validating input.species_constraints arguments
    """
    allowed: List[str] = ['input species', 'seed mechanisms', 'reaction libraries']
    #max_C_atoms: conint(ge=0)
    #max_O_atoms: conint(ge=0)
    #max_N_atoms: conint(ge=0)
    #max_Si_atoms: conint(ge=0)
    #max_S_atoms: conint(ge=0)
    #max_heavy_atoms: conint(ge=0)
    #max_radical_electrons: conint(ge=0)
    #max_singlet_carbenes: conint(ge=0) = 1
    #max_carbene_radicals: conint(ge=0) = 0
    #allow_singlet_O2: bool = True

    class Config:
        extra = "forbid"

    @field_validator('allowed')
    def check_allowed(cls, value):
        """SpeciesConstraints.allowed field_validator"""
        for val in value:
            if val not in ['input species', 'seed mechanisms', 'reaction libraries']:
                raise ValueError(f"The allowed species in the  species constraints list must be in\n"
                                 f"['input species', 'seed mechanisms', 'reaction libraries'].\n"
                                 f"Got: {val} in {value}")
        return value

class Database(BaseModel):
    name: str
    thermo_libraries: Optional[List[str]] = None
    kinetics_libraries: Optional[List[str]] = None
    use_low_credence_libraries: bool = False # True if low credence libraries should be used
    #transport_libraries: List[str] = ['OneDMinN2', 'PrimaryTransportLibrary', 'NOx2018', 'GRI-Mech']
    kinetics_depositories: Union[List[str], str] = 'default' # depositories are used to estimate the rate law 
    kinetics_families: Union[str, List[str]] = 'default'
    parameter_estimator: Literal["ML", "group_contribution", "fixed_defaults"] = "ML" # how to estimate 
    rate_law: Literal["Michaelis-Menten", "Hill", "MassAction"]
    parameters: Optional[Dict[str, float]]  # e.g., {"Km": 0.1, "Vmax": 2.0}
    species_constraints: Optional[SpeciesConstraints] = None
    

class RadicalTypeEnum(str, Enum):
    """
    The supported radical ``types`` entries for ``generate_radicals()``.
    """
    radical = 'radical'
    alkoxyl = 'alkoxyl'
    peroxyl = 'peroxyl'


class Environment(BaseModel):
    """
    A class for validating input.Enviornment arguments
    """
    temperature: Union[confloat(gt=0), List[confloat(gt=0)]]  # in Celsius
    pH: Union[confloat(gt=0), List[confloat(gt=0)]]
    ionic_strength: Optional[float] = None
    oxygen_level: Optional[float] = None
    seed_mechanisms: List[str] = list()   # todo: database

    @field_validator('temperature')
    def check_t(cls, value):
        """Environment.temperature field_validator"""
        if isinstance(value, list) and len(value) != 2:
            raise ValueError(f'When specifying the temperature as a list, only two values are allowed (T min, T max),\n'
                             f'got {len(value)} values: {value}.')
        return value
    @field_validator('pH')
    def validate_pH(cls, value):
        if not (0 <= value <= 14):
            raise ValueError(f'pH must be between 0 and 14. Got: {value}')
        return value
    

    @field_validator('ionic_strength')
    def validate_ionic_strength(cls, value):
        """Environment.ionic_strength field_validator"""
        if value is not None and value < 0:
            raise ValueError(f'Ionic strength cannot be negative. Got: {value}')
        return value

    @field_validator('oxygen_level')
    def validate_oxygen_level(cls, value):
        """Environment.oxgen_level field_validator"""
        if value is not None and not (0 <= value <= 1):
            raise ValueError(f'Oxygen level must be between 0 and 1. Got: {value}')
        return value
    
    
class Settings(BaseModel):
    """
    A class for validating input.settings arguments     
    """

    end_time: confloat(gt=0)  
    time_step: confloat(gt=0)
    solver: Literal["odeint", "CVODE", "custom"]
    use_core_edge: bool = True
    threshold: confloat(gt=0) = 1e-3 #epsilon
    max_iterations: confloat(gt=0) = 50
    ml_models_enabled: bool = True
    stop_at_steady_state: bool = True
    termination_conversion: Optional[Dict[str, confloat(gt=0, lt=1)]] = None
    termination_rate_ratio: Optional[confloat(gt=0, lt=1)] = None  # Allows convergence-based logic

    class Config:
        extra = "forbid"
    @field_validator('time_step')
    def validate_time_step(cls, value, values):
        if 'end_time' in values and value >= values['end_time']:
            raise ValueError(f"'time_step' must be smaller than total simulation time {values['end_time']}. Got: {value}")
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
    """
    An InputBase class for validating input arguments
    """
    project: constr(max_length=255)
    project_directory: Optional[constr(max_length=255)] = None
    species: List[Species]
    enzymes: List[Enzyme]
    environment: Environment
    database: Database
    settings: Settings

    class Config:
        extra = "forbid"
