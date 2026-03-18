"""
BEES schema module
used for input validation
"""

from typing import Dict, List, Optional, Tuple, Union, Literal, Annotated
from pydantic import BaseModel, conint, confloat, constr, field_validator, ValidationInfo, Field, ConfigDict
from rdkit import Chem


class Species(BaseModel):
    """
    A class for validate input.BEES.Species arguments.
    """
    label: str
    concentration: Union[Annotated[float, Field(gt=0)], Tuple[Annotated[float, Field(gt=0)], Annotated[float, Field(gt=0)]]] = None
    smiles: Optional[str] = None
    constant: bool = False
    reactive: bool = True
    solvent: bool = False
  
    model_config = ConfigDict(extra="forbid")


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

class Enzyme(Species):
    """
    A class for validate input.BEES.Enzyme arguments if there are any.
    Inherits from Species, adding specific fields for enzymes.
    """
    ecnumber: Optional[Union[constr(pattern=r"^EC \d+\.\d+\.\d+\.\d+$"), List[constr(pattern=r"^EC \d+\.\d+\.\d+\.\d+$")]]] = None
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
        if value is None:
            return value
        entries = [value] if isinstance(value, str) else value
        for ec in entries:
            if not ec.startswith("EC "):
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
    temperature: Union[confloat(gt=0), Tuple[confloat(gt=0), confloat(gt=0)]]  
    pH: Union[confloat(ge=0, le=14), Tuple[confloat(ge=0, le=14), confloat(ge=0, le=14)]]  
    

    model_config = ConfigDict(extra="forbid")

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

    end_time: Optional[confloat(gt=0)] = None
    time_step: Optional[confloat(gt=0)] = None

    estimate_kinetics: bool = False
    kinetics_estimator: Optional[Literal['catpred']] = None
    kinetics_include_sd: bool = True
    smiles_mode: Literal['auto', 'interactive'] = 'auto'

    toleranceKeepInEdge: confloat(ge=0) = 0
    toleranceMoveToCore: confloat(gt=0) = 1e-5
    # ODE early-stop (flux ratio |R_edge|/R_char); defaults to toleranceMoveToCore if unset.
    toleranceInterruptSimulation: Optional[confloat(gt=0)] = None
    # reaction-level promotion criterion. dlnaccum_j = Σ_i ln(1 + v_j / R_i)
    # over species i touched by edge reaction j; R_i is consumption for reactants
    # and production for products. None disables this criterion.
    toleranceMoveEdgeReactionToCore: Optional[confloat(gt=0)] = None
    minEdgeIterationsForPrune: conint(ge=0) = 2
    minCoreSpeciesForPrune: conint(ge=0) = 0
    termination_conversion: Optional[Dict[str, confloat(gt=0, lt=1)]] = None
    termination_rate_ratio: Optional[confloat(gt=0, lt=1)] = None
    max_edge_species: Optional[conint(gt=0)] = None
    max_num_objects_per_iter: conint(gt=0) = 10
    # When R_char == 0 (flat core), promote edge species only if |R_i| > this
    # absolute floor (mM/s). Prevents spurious promotions from numerical noise
    # while still allowing species with real flux to enter the core.
    abs_flux_floor: confloat(gt=0) = 1e-12
    max_iterations: conint(gt=0) = 50

    # ODE solver (scipy.integrate.solve_ivp). Unset fields use simulator defaults.
    ode_method: Optional[str] = None
    ode_rtol: Optional[confloat(gt=0)] = None
    ode_atol: Optional[confloat(gt=0)] = None
    # Wall-clock limit for one enlargement pass (stepwise ODE only). None = no limit.
    max_wall_time_per_iteration: Optional[confloat(gt=0)] = None
    # INFO heartbeat interval (seconds) during stepwise integration. 0 = disabled.
    # None uses the simulator default (30 s).
    stepwise_heartbeat_interval: Optional[confloat(ge=0)] = None

    verbose: Optional[conint(ge=10, le=50)] = 20
    saveEdgeSpecies: bool = True
    output_directory: Optional[str] = None
    save_simulation_profiles: bool = True
    save_simulation_plots: bool = True 
    save_reaction_tree_plots: bool = False
    reaction_tree_layout: Literal["graphviz", "simple"] = "graphviz"
    reaction_tree_rankdir: Literal["TB", "BT", "LR", "RL"] = "TB"
    reaction_tree_fontsize: conint(ge=8, le=24) = 8
    plot_max_species: Optional[conint(gt=0)] = None
    plot_exclude_enzymes: bool = True
    plot_exclude_cofactors: bool = True
    save_ode_equations: bool = True

    model_config = ConfigDict(extra="forbid")

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


class Database(BaseModel):
    name: constr(min_length=1)

    model_config = ConfigDict(extra="forbid")

    @classmethod
    @field_validator('name')
    def check_name_not_empty(cls, value):
        if not value.strip():
            raise ValueError("Name cannot be empty")
        return value


class InputBase(BaseModel):
    project: constr(max_length=255)
    project_directory: Optional[constr(max_length=255)] = None
    species: List[Species]
    enzymes: List[Enzyme]
    environment: Environment
    settings: Settings
    database: Database

    model_config = ConfigDict(extra="forbid")

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

