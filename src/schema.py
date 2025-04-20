"""
BEES schema module
used for input validation
"""

from enum import Enum
from typing import Dict, List, Optional, Tuple, Union, Literal 

from pydantic import BaseModel, conint, confloat, constr, root_validator, validator
from rdkit import Chem


class TerminationTimeEnum(str, Enum):
    """
    The supported termination type units in an BEES "reactor".
    """
    micro_s = 'micro-s'
    ms = 'ms'
    s = 's'
    hrs = 'hrs'
    minutes = 'min'
    hours = 'hours'
    days = 'days'

class BEESCompounds(BaseModel):
    """
    A class for validating input.BEES.compound arguments
    """
    label: str
    type: Literal["substrate", "enzyme", "cofactor"]
    concentration: Union[confloat(gt=0), Tuple[confloat(gt=0), confloat(gt=0)]] = None # concentration in mol/L
    structure_smiles: Optional[str] = None
    structure_inchi: Optional[str] = None 
    charge: Optional[float] = None
    compartment: Optional[str] = None  # e.g., 'cytoplasm', 
    kinetic_params: Optional[Dict[str, float]] = None # e.g., {"Km": 0.1, "Vmax": 2.0}
    reactive: bool = True
    constant: bool = False
    observable: bool = False # True if the compound is observable, meaning it can be measured and tracked in the simulation
    SA_observable: bool = False # True if the compound is SA observable, meaning it can be measured and tracked in the sanstivity analysis simulation
    UA_observable: bool = False #  True if the compound is UA observable, meaning it can be measured and tracked in the uncertainty analysis simulation
    solvent: bool = False # True if the compound is a solvent
    xyz: Optional[Union[List[Union[dict, str]], dict, str]] = None # Holds 3D coordinates (geometry)
    #seed_all_rads: Optional[List[RadicalTypeEnum]] = None
    
    
    class Config:
        extra = "forbid"


    @validator('concentrations')
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

    @validator('constant')
    def check_constant_species(cls, value, values):
        if value and isinstance(values.get('concentration'), tuple):
            raise ValueError("Constant species cannot have a concentration range.")
        return value
    
    #@validator('reactive') 
    #def check_reactive(cls, value, values):
        #"""Check if the compound is reactive"""
        if value and values.get('constant'):
            #raise ValueError("Reactive species cannot be constant.")
        #return value

    #@validator('observable')       
    #def check_observable(cls, value, values):  
        #"""Check if the compound is observable"""
        #if value and values.get('constant'):
            #raise ValueError("Observable species cannot be constant.")
        #return value

    #@validator('structure_smiles')
    #def check_structure_smiles(cls, value):
            """Check if the structure is valid"""
        if value and not isinstance(value, str):
            raise ValueError("Structure must be a string.")
        if Chem.MolFromSmiles(value) is None:
            raise ValueError(f"Invalid SMILES string: {value}")
        #return value
    #@validator('structure_inchi')
    #def validate_inchi(cls, value):
    #"""Check that InChI string is chemically valid"""
    #if value is None:
        return value
    #try:
        mol = inchi.MolFromInchi(value)
        if mol is None:
            raise ValueError
    #except Exception:
        raise ValueError(f"Invalid InChI string: {value}")
    #return value

class BEESSpeciesConstraints(BaseModel):
    """
    A class for validating input.BEES.species_constraints arguments
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

    @validator('allowed')
    def check_allowed(cls, value):
        """BEESSpeciesConstraints.allowed validator"""
        for val in value:
            if val not in ['input species', 'seed mechanisms', 'reaction libraries']:
                raise ValueError(f"The allowed species in the BEES species constraints list must be in\n"
                                 f"['input species', 'seed mechanisms', 'reaction libraries'].\n"
                                 f"Got: {val} in {value}")
        return value

class BEESReactionRule(BaseModel):
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
    species_constraints: Optional[BEESSpeciesConstraints] = None
    

class RadicalTypeEnum(str, Enum):
    """
    The supported radical ``types`` entries for ``generate_radicals()``.
    """
    radical = 'radical'
    alkoxyl = 'alkoxyl'
    peroxyl = 'peroxyl'


class BEESEnvironment(BaseModel):
    """
    A class for validating input.BEES.Enviornment arguments
    """
    
    temperature: Union[confloat(gt=0), List[confloat(gt=0)]]  # in Celsius or Kelvin ?
    pH: Union[confloat(gt=0), List[confloat(gt=0)]]
   #Volume: Optional[Union[confloat(gt=0), List[confloat(gt=0)]]] = None
    ionic_strength: Optional[float] = None
    Oxygen_level: Optional[float] = None
    seed_mechanisms: List[str] = list()


    @validator('temperature')
    def check_t(cls, value):
        """BEESEnvironment.temperature validator"""
        if isinstance(value, list) and len(value) != 2:
            raise ValueError(f'When specifying the temperature as a list, only two values are allowed (T min, T max),\n'
                             f'got {len(value)} values: {value}.')
        return value
    @validator('pH')
    def validate_pH(cls, value):
        if not (0 <= value <= 14):
            raise ValueError(f'pH must be between 0 and 14. Got: {value}')
        return value
    

    @validator('ionic_strength')
    def validate_ionic_strength(cls, value):
        """BEESEnvironment.ionic_strength validator"""
        if value is not None and value < 0:
            raise ValueError(f'Ionic strength cannot be negative. Got: {value}')
        return value

    @validator('Oxygen_level')
    def validate_oxygen_level(cls, value):
        """BEESEnvironment.Oxgen_level validator"""
        if value is not None and not (0 <= value <= 1):
            raise ValueError(f'Oxygen level must be between 0 and 1. Got: {value}')
        return value
    
    #@validator('V', always=True)
    #def check_v(cls, value, values):
        """BEESEnvironment.volume validator"""
        if isinstance(value, list) and len(value) != 2:
            raise ValueError(f'When specifying the volume as a list, only two values are allowed (V min, V max),\n'
                             f'got {len(value)} values: {value}.')
        if 'type' in values and 'liquid' in values['type'] and value is None:
            raise ValueError('The reactor volume must be specified for a liquid-phase reactor.')
    
    
class BEESModelSettings(BaseModel):
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
    @validator('time_step')
    def validate_time_step(cls, value, values):
        if 'end_time' in values and value >= values['end_time']:
            raise ValueError(f"'time_step' must be smaller than total simulation time {values['end_time']}. Got: {value}")
        return value

    @validator('termination_conversion')
    def validate_termination_conversion(cls, value):
        if value:
            for species, frac in value.items():
                if not (0 < frac < 1):
                    raise ValueError(f"termination_conversion values must be between 0 and 1. Got: {species}: {frac}")
        return value

    @validator('termination_rate_ratio')
    def validate_rate_ratio(cls, value):
        if value and not (0 < value < 1):
            raise ValueError("termination_rate_ratio must be between 0 and 1 (exclusive).")
        return value

class BEESInputBase(BaseModel):
    """
    An InputBase class for validating input arguments
    """
    BEES_execution_type: Optional[str] = None
    memory: Optional[conint(ge=0)] = None
    cpus: Optional[conint(gt=0)] = None
    project: constr(max_length=255)
    project_directory: Optional[constr(max_length=255)] = None
    compounds: List[BEESCompounds]
    environment: BEESEnvironment
    rules: List[BEESReactionRule]
    model_settings: BEESModelSettings

    #verbose: conint(ge=10, le=30, multiple_of=10) = 20     #didn't get what this line doing 
    
    #qm: Optional[QM] = None

    class Config:
        extra = "forbid"

    #@validator('t3', always=True)
    #def check_t3(cls, value):
        #"""InputBase.t3 validator"""
        #return value or T3()

    #@validator('qm', always=True)
    #def check_qm(cls, value):
        """InputBase.qm validator"""
        # Removed invalid return statement
        pass

    #@root_validator(pre=True)
    #def validate_rmg_t3(cls, values):
       # """InputBase.validate_rmg_t3"""
        #if 'rmg' in values and 't3' in values and values['t3']:
            # check termination time for global UA
         #   if 'uncertainty' in values['t3'] and values['t3']['uncertainty'] is not None:
          #      ua_termination_time = values['t3']['uncertainty']['termination_time']
           #     rmg_reactor_termination_times = [reactor['termination_time'] for reactor in values['rmg']['reactors']]
            #    if all([termination_time is None for termination_time in rmg_reactor_termination_times]) \
             #           and ua_termination_time is None and values['t3']['uncertainty']['global_analysis']:
              #      raise ValueError('If a global uncertainty analysis is requested, a termination time must be '
               #                      'specified either under "t3.uncertainty.termination_time" '
                #                     'or in at least one RMG reactor.')
            # check solvent for liquid phase
            #reactor_types = set([reactor['type'] for reactor in values['rmg']['reactors']])
            #solvents = list()
            #for species in values['rmg']['species']:
                #if 'solvent' in species and species['solvent']:
                    #solvents.append(species['label'])
            #if any(['liquid' in reactor_type for reactor_type in reactor_types]):
                #if not len(solvents):
                    #raise ValueError('One species must be defined as the solvent when using liquid phase reactors.')
                #if len(solvents) > 1:
                    #raise ValueError(f'Only one solvent can be specified, got: {solvents}')
            #else:
                #if len(solvents):
                 #   raise ValueError(f'No solvent species are allowed for gas phase reactors, got: {solvents}')
            # check core_tolerance and max_T3_iterations
            #if 'model' in values['rmg'] and 'core_tolerance' in values['rmg']['model'] \
             #       and 'options' in values['t3'] and 'max_T3_iterations' in values['t3']['options'] \
              #      and not isinstance(values['rmg']['model']['core_tolerance'], float) \
               #     and len(values['rmg']['model']['core_tolerance']) > values['t3']['options']['max_T3_iterations']:
                #raise ValueError(f'The number of RMG core tolerances ({len(values["rmg"]["model"]["core_tolerance"])}) '
                 #                f'cannot be greater than the max number of T3 iterations '
                   #              f'({values["t3"]["options"]["max_T3_iterations"]}).')
        #return values
