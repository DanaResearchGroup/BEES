"""Pydantic input schema for BEES YAML."""

import math
from typing import Any, Dict, List, Optional, Tuple, Union, Literal, Annotated
from pydantic import BaseModel, ConfigDict, conint, confloat, constr, field_validator, model_validator, ValidationInfo, Field
from rdkit import Chem

from bees.dosing import normalize_params

class Species(BaseModel):
    label: str
    # 0 is allowed: a drug species that only enters through `dosing`.
    concentration: Union[Annotated[float, Field(ge=0)], Tuple[Annotated[float, Field(ge=0)], Annotated[float, Field(ge=0)]]] = Field(
        default=None, validate_default=True,
    )
    smiles: Optional[str] = None
    constant: bool = False
    reactive: bool = True
    solvent: bool = False
  
    model_config = ConfigDict(extra="forbid")

    @field_validator('concentration', mode='before')
    @classmethod
    def check_concentration_given(cls, value, info: ValidationInfo):
        if value is None:
            raise ValueError(f"Concentration must be specified for '{info.data.get('label')}'")
        return value

    @field_validator('concentration')
    @classmethod
    def check_concentration_range_order(cls, value, info: ValidationInfo):
        label = info.data.get('label')
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
    @classmethod
    def check_constant_species(cls, value, info: ValidationInfo):
        if value:
            if info.data.get('concentration') and isinstance(info.data['concentration'], tuple):
                raise ValueError("Constant species cannot have a concentration range")
        return value

    @field_validator('smiles')
    @classmethod
    def validate_smiles(cls, value):
        if value:
            try:
                mol = Chem.MolFromSmiles(value)
                if not mol:
                    raise ValueError("Invalid SMILES string")
            except Exception:
                raise ValueError("Invalid SMILES string")
        return value

class FeedbackInhibition(BaseModel):
    """End-product/feedback inhibition on an enzyme (by label, not EC)."""
    inhibitors: Annotated[List[constr(min_length=1)], Field(min_length=1)]
    ki: Optional[Annotated[float, Field(gt=0)]] = None  # mM
    hill: Annotated[float, Field(gt=0)] = 1.0

    model_config = ConfigDict(extra="forbid")

class Inhibition(BaseModel):
    """Reversible inhibition of an enzyme by a model species (see bees/inhibition.py)."""
    inhibitor: constr(min_length=1)  # label of a species
    type: Literal["competitive", "uncompetitive", "noncompetitive", "mixed", "tight_binding"]
    ki: Optional[Annotated[float, Field(gt=0)]] = None  # mM; None -> estimated (CatPred)
    alpha: Optional[Annotated[float, Field(gt=0)]] = None  # mixed only: Ki' = alpha * Ki
    substrate: Optional[constr(min_length=1)] = None  # which substrate it acts against
    source: Optional[str] = None  # where Ki comes from (free text)

    model_config = ConfigDict(extra="forbid")

    @model_validator(mode="after")
    def _check_type_fields(self):
        if self.type == "mixed" and self.alpha is None:
            raise ValueError("Mixed inhibition needs 'alpha' (Ki' = alpha * Ki)")
        if self.type != "mixed" and self.alpha is not None:
            raise ValueError(f"'alpha' only applies to mixed inhibition, not '{self.type}'")
        if self.substrate is not None and self.type in ("noncompetitive", "tight_binding"):
            raise ValueError(f"'substrate' does not apply to {self.type} inhibition (it acts on the whole enzyme)")
        return self

class MeasuredKinetics(BaseModel):
    """Experimental kinetics of an enzyme; override database/CatPred values."""
    kcat: Optional[Annotated[float, Field(gt=0)]] = None  # 1/s
    km: Optional[Dict[str, Annotated[float, Field(gt=0)]]] = None  # substrate label -> mM
    # Km for every substrate not listed in `km`. A substrate without a Km is left out of the
    # rate law, i.e. treated as saturated, which lets the reaction run after it is used up.
    default_km: Optional[Annotated[float, Field(gt=0)]] = None  # mM
    source: Optional[str] = None

    model_config = ConfigDict(extra="forbid")

    @model_validator(mode="after")
    def _check_something_given(self):
        if self.kcat is None and not self.km and self.default_km is None:
            raise ValueError("measured_kinetics needs at least one of kcat, km, default_km")
        return self

class Enzyme(Species):
    ecnumber: Optional[Union[constr(pattern=r"^EC \d+\.\d+\.\d+\.\d+$"), List[constr(pattern=r"^EC \d+\.\d+\.\d+\.\d+$")]]] = None
    amino_acid_sequence: Optional[str] = None
    feedback_inhibition: Optional[FeedbackInhibition] = None
    inhibition: Optional[List[Inhibition]] = None
    measured_kinetics: Optional[MeasuredKinetics] = None

    @field_validator('label')
    @classmethod
    def check_label_not_empty(cls, value):
        if not value.strip():
            raise ValueError("Label cannot be empty")
        return value

    @field_validator('ecnumber')
    @classmethod
    def validate_ecnumber(cls, value):
        if value is None:
            return value
        entries = [value] if isinstance(value, str) else value
        for ec in entries:
            if not ec.startswith("EC "):
                raise ValueError("Invalid EC number: must start with 'EC '")
        return value

    @field_validator('amino_acid_sequence')
    @classmethod
    def validate_amino_acid_sequence(cls, value):
        if value is None:
            return value
        
        if ' ' in value:
            raise ValueError("Amino acid sequence cannot contain spaces")
        if not value.isupper():
            raise ValueError("Amino acid sequence must be in capital letters")
        valid_aa_codes = set('ACDEFGHIKLMNPQRSTVWY')
        invalid_chars = set(value) - valid_aa_codes
        if invalid_chars:
            raise ValueError(f"Amino acid sequence contains invalid characters: {', '.join(sorted(invalid_chars))}. Valid codes are: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y")
        
        return value

class Environment(BaseModel):
    temperature: Union[confloat(gt=0), Tuple[confloat(gt=0), confloat(gt=0)]]
    pH: Union[confloat(ge=0, le=14), Tuple[confloat(ge=0, le=14), confloat(ge=0, le=14)]]
    ionic_strength_M: confloat(ge=0, le=2.0) = 0.25
    pMg: confloat(ge=0, le=10) = 3.0

    model_config = ConfigDict(extra="forbid")

    @field_validator('temperature')
    @classmethod
    def validate_temperature_range(cls, value):
        if isinstance(value, tuple) and len(value) != 2:
            raise ValueError("Temperature as list must have exactly 2 values (min, max)")
        if isinstance(value, tuple) and value[0] > value[1]:
            raise ValueError("Temperature range min value cannot be greater than max value")
        return value

    @field_validator('pH')
    @classmethod
    def validate_pH_range(cls, value):
        if isinstance(value, tuple) and len(value) != 2:
            raise ValueError("pH as list must have exactly 2 values (min, max)")
        if isinstance(value, tuple) and value[0] > value[1]:
            raise ValueError("pH range min value cannot be greater than max value")
        return value

class Settings(BaseModel):
    end_time: Optional[confloat(gt=0)] = None
    time_step: Optional[confloat(gt=0)] = None
    simulation_mode: Optional[Literal["iterative", "batch"]] = None

    estimate_kinetics: bool = False
    kinetics_estimator: Optional[Literal['catpred']] = None
    kinetics_include_sd: bool = True
    # DEPRECATED: YAML kcat scale; use settings.calibrations.
    kinetics_ec_kcat_scale: Optional[Dict[str, Annotated[float, Field(gt=0)]]] = None
    calibrations: Optional[List[str]] = None  
    smiles_mode: Literal['auto', 'interactive'] = 'auto'

    toleranceKeepInEdge: confloat(ge=0) = 0
    toleranceMoveToCore: confloat(gt=0) = 1e-5
    toleranceInterruptSimulation: Optional[confloat(gt=0)] = None  # |R_edge|/R_char; default = toleranceMoveToCore
    toleranceMoveEdgeReactionToCore: Optional[confloat(gt=0)] = None  # dlnaccum; None disables
    minEdgeIterationsForPrune: conint(ge=0) = 2
    minCoreSpeciesForPrune: conint(ge=0) = 0
    termination_conversion: Optional[Dict[str, confloat(gt=0, lt=1)]] = None
    termination_rate_ratio: Optional[confloat(gt=0, lt=1)] = None
    max_edge_species: Optional[conint(gt=0)] = None
    max_num_objects_per_iter: conint(gt=0) = 10
    abs_flux_floor: confloat(gt=0) = 1e-12  # mM/s noise floor when R_char == 0
    max_iterations: conint(gt=0) = 50

    ode_method: Optional[str] = None
    ode_rtol: Optional[confloat(gt=0)] = None
    ode_atol: Optional[confloat(gt=0)] = None
    max_wall_time_per_iteration: Optional[confloat(gt=0)] = None  # stepwise ODE only

    verbose: Optional[conint(ge=10, le=50)] = 20
    saveEdgeSpecies: bool = True
    output_directory: Optional[str] = None
    save_simulation_profiles: bool = True
    save_simulation_plots: bool = True 
    plot_max_species: Optional[conint(gt=0)] = None
    plot_exclude_enzymes: bool = True
    plot_exclude_cofactors: bool = True
    save_ode_equations: bool = True

    thermo_smiles_substitutions: Optional[Dict[str, str]] = None
    thermo_irreversible_cutoff_kJmol: Annotated[float, Field(gt=0)] = 30.0  # ΔG°' < −this kJ/mol

    # Reactions of one enzyme share its free enzyme (partition function over CM bound forms).
    enzyme_competition: bool = True

    model_config = ConfigDict(extra="forbid")

    @field_validator('time_step')
    @classmethod
    def validate_time_step(cls, value, info: ValidationInfo):
        end_time = info.data.get('end_time')
        if end_time is not None and value >= end_time:
            raise ValueError(f"'time_step' must be smaller than 'end_time' ({end_time}). Got: {value}")
        return value

    @field_validator('termination_conversion')
    @classmethod
    def validate_termination_conversion(cls, value):
        if value:
            for species, frac in value.items():
                if not (0 < frac < 1):
                    raise ValueError(f"termination_conversion values must be between 0 and 1. Got: {species}: {frac}")
        return value

    @field_validator('termination_rate_ratio')
    @classmethod
    def validate_rate_ratio(cls, value):
        if value and not (0 < value < 1):
            raise ValueError("termination_rate_ratio must be between 0 and 1 (exclusive).")
        return value

    @model_validator(mode="after")
    def _require_at_least_one_termination(self):
        if self.simulation_mode == "batch":
            return self
        if (
            self.end_time is None
            and not self.termination_conversion
            and self.termination_rate_ratio is None
        ):
            raise ValueError(
                "At least one termination criterion must be set for iterative mode: "
                "end_time, termination_conversion, or termination_rate_ratio."
            )
        return self

    @field_validator('verbose')
    @classmethod
    def validate_verbose_level(cls, value):
        if value is not None and value not in [10, 20, 30, 40, 50]:
            raise ValueError("Verbose level must be 10, 20, 30, 40, or 50")
        return value

class Database(BaseModel):
    name: constr(min_length=1)

    model_config = ConfigDict(extra="forbid")

    @field_validator('name')
    @classmethod
    def check_name_not_empty(cls, value):
        if not value.strip():
            raise ValueError("Name cannot be empty")
        return value

class Dosing(BaseModel):
    """Time-dependent input of a species (a drug released by its carrier); see bees/dosing.py."""
    species: constr(min_length=1)
    profile: Literal[
        "bolus", "constant", "first_order", "biexponential", "weibull", "korsmeyer_peppas", "table",
    ]
    dose: Annotated[float, Field(gt=0)]  # mM released in total
    t_start: Annotated[float, Field(ge=0)] = 0.0  # s
    params: Optional[Dict[str, Any]] = None  # profile parameters, in seconds
    elimination_k: Optional[Annotated[float, Field(gt=0)]] = None  # 1/s
    elimination_half_life: Optional[Annotated[float, Field(gt=0)]] = None  # s

    model_config = ConfigDict(extra="forbid")

    @model_validator(mode="after")
    def _check_profile(self):
        if self.elimination_k is not None and self.elimination_half_life is not None:
            raise ValueError("Give elimination_k or elimination_half_life, not both")
        normalize_params(self.profile, self.params)  # raises with a message naming the profile
        return self

class InputBase(BaseModel):
    project: constr(max_length=255)
    project_directory: Optional[constr(max_length=255)] = None
    species: List[Species]
    enzymes: List[Enzyme]
    environment: Environment
    settings: Settings
    database: Database
    dosing: Optional[List[Dosing]] = None

    model_config = ConfigDict(extra="forbid")

    @field_validator('project')
    @classmethod
    def check_project_not_empty(cls, value):
        if not value.strip():
            raise ValueError("Project name cannot be empty")
        return value

    @field_validator('species')
    @classmethod
    def check_species_list_not_empty(cls, value):
        if not value:
            raise ValueError("Species list cannot be empty")
        return value

    @field_validator('enzymes')
    @classmethod
    def check_enzymes_list_not_empty(cls, value):
        if not value:
            raise ValueError("Enzymes list cannot be empty")
        return value

    @model_validator(mode="after")
    def _check_cross_references(self):
        species = {sp.label.lower().strip(): sp for sp in self.species}
        enzyme_labels = {e.label.lower().strip() for e in self.enzymes}
        for enz in self.enzymes:
            for inh in enz.inhibition or []:
                key = inh.inhibitor.lower().strip()
                sp = species.get(key)
                if sp is None:
                    raise ValueError(
                        f"Enzyme '{enz.label}': inhibitor '{inh.inhibitor}' is not in 'species'"
                    )
                if inh.ki is None:
                    if not self.settings.estimate_kinetics:
                        raise ValueError(
                            f"Enzyme '{enz.label}': Ki for '{inh.inhibitor}' is omitted, which asks for an "
                            "estimate; set settings.estimate_kinetics or give 'ki'"
                        )
                    if not sp.smiles or not enz.amino_acid_sequence:
                        raise ValueError(
                            f"Enzyme '{enz.label}': estimating Ki for '{inh.inhibitor}' needs the inhibitor's "
                            "'smiles' and the enzyme's 'amino_acid_sequence'; otherwise give 'ki'"
                        )
        elimination: Dict[str, float] = {}
        for d in self.dosing or []:
            key = d.species.lower().strip()
            sp = species.get(key)
            if sp is None:
                raise ValueError(f"Dosing: species '{d.species}' is not in 'species'")
            if sp.constant:
                raise ValueError(f"Dosing: species '{d.species}' is constant, so it cannot be dosed")
            if isinstance(sp.concentration, tuple):
                raise ValueError(f"Dosing: species '{d.species}' must have a single concentration")
            k = d.elimination_k
            if k is None and d.elimination_half_life is not None:
                k = math.log(2.0) / d.elimination_half_life
            if k is not None:
                if key in elimination and not math.isclose(elimination[key], k, rel_tol=1e-9):
                    raise ValueError(f"Dosing: conflicting elimination rates for species '{d.species}'")
                elimination[key] = k
        return self
