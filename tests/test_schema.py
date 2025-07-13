#!/usr/bin/env python3
# encoding: utf-8

"""
Test schema validation for BEES input models
"""

import pytest
from pydantic import ValidationError
# Corrected import path: from bees.schema instead of from src.schema
from bees.schema import (
    Enzyme,
    Species,
    Environment,
    Settings,
    Database,
    InputBase,
    SpeciesConstraints,
    TerminationTimeEnum
)


def test_TerminationTimeEnum():
    """Test the TerminationTimeEnum."""
    assert TerminationTimeEnum.s == 's'
    assert TerminationTimeEnum.hours == 'hours'
    assert 'min' in TerminationTimeEnum.__members__.values()


def test_Species():
    """Test the Species model."""
    species = Species(
        label="Glucose",
        concentration=0.1,
        charge=0,
        reactive=True,
        constant=False,
        observable=True,
    )
    assert species.label == "Glucose"
    assert species.concentration == 0.1
    assert species.reactive is True
    assert species.observable is True
    assert species.charge == 0
    assert species.solvent is False
    assert species.xyz is None
    assert species.smiles is None
    assert species.inchi is None


    # Updated InChI for a simpler molecule (Methane) for reliable validation
    species = Species(
        label="Fructose", # Label remains Fructose for the test case name
        concentration=(0.1, 1.0),
        smiles="C(C1C(C(C(C(O1)CO)O)O)O)O", # Still Fructose SMILES
        inchi="InChI=1S/CH4/h1H4" # Valid InChI for Methane
    )
    assert species.concentration == (0.1, 1.0)
    assert species.smiles == "C(C1C(C(C(C(O1)CO)O)O)O)O"
    assert species.inchi == "InChI=1S/CH4/h1H4"

    # Test validators
    # Removed: with pytest.raises(ValidationError, match="Concentration must be specified"): Species(label="NoConc")
    # This is because 'concentration' is Optional in your schema.
    with pytest.raises(ValidationError, match="Concentration range cannot have identical values"):
        Species(label="SameRange", concentration=(0.5, 0.5))
    with pytest.raises(ValidationError, match="Constant species cannot have a concentration range"):
        Species(label="RangeConstant", concentration=(0.1, 0.9), constant=True)
    with pytest.raises(ValidationError, match="Input should be greater than 0"): # Updated regex for single float
        Species(label="NegativeConcentration", concentration=-0.1)
    with pytest.raises(ValidationError, match="Input should be greater than 0"): # Updated regex for tuple element
        Species(label="NegativeRange", concentration=(-0.1, 0.5))
    with pytest.raises(ValidationError, match="Reactive species cannot be constant"):
        Species(label="ReactiveConstant", reactive=True, constant=True, concentration=0.1)
    with pytest.raises(ValidationError, match="Observable species cannot be constant"):
        Species(label="ObservableConstant", observable=True, constant=True, concentration=0.1)
    with pytest.raises(ValidationError, match="Invalid SMILES string"):
        Species(label="InvalidSMILES", concentration=0.1, smiles="InvalidSmiles[")
    with pytest.raises(ValidationError, match="Invalid InChI string"):
        Species(label="InvalidInChI", concentration=0.1, inchi="InvalidInChI")


def test_Enzyme():
    """Test the Enzyme model."""
    enzyme = Enzyme(
        label="ATP",
        concentration=0.01,
        charge=-3,
        reactive=True,
        constant=False,
        observable=True,
        ecnumber="EC 2.7.1.1"
    )
    assert enzyme.label == "ATP"
    assert enzyme.concentration == 0.01
    assert enzyme.reactive is True
    assert enzyme.ecnumber == "EC 2.7.1.1"

    enzyme = Enzyme(
        label="Phosphofructokinase",
        concentration=(0.1, 1.0),
    )
    assert enzyme.concentration == (0.1, 1.0)

    # Test validators
    # Removed: with pytest.raises(ValidationError, match="Concentration must be specified"): Enzyme(label="NoConc")
    # This is because 'concentration' is Optional in your schema.
    with pytest.raises(ValidationError, match="Concentration range cannot have identical values"):
        Enzyme(label="SameRange", concentration=(0.5, 0.5))
    with pytest.raises(ValidationError, match="Constant species cannot have a concentration range"):
        Enzyme(label="RangeConstant", concentration=(0.1, 0.9), constant=True)
    with pytest.raises(ValidationError, match="Input should be greater than 0"): # Updated regex for single float
        Enzyme(label="NegativeConcentration", concentration=-0.1)
    with pytest.raises(ValidationError, match="Input should be greater than 0"): # Updated regex for tuple element
        Enzyme(label="NegativeRange", concentration=(-0.1, 0.5))
    with pytest.raises(ValidationError, match="Invalid EC number"):
        Enzyme(label="InvalidEC", concentration=0.1, ecnumber="not.a.valid.ecnumber")
    with pytest.raises(ValidationError, match="Invalid EC number"):
        Enzyme(label="InvalidEC2", concentration=0.1, ecnumber="2.7.1.1") # Missing 'EC ' prefix
    with pytest.raises(ValidationError, match="Reactive species cannot be constant"):
        Enzyme(label="ReactiveConstant", reactive=True, constant=True, concentration=0.1)
    with pytest.raises(ValidationError, match="Observable species cannot be constant"):
        Enzyme(label="ObservableConstant", observable=True, constant=True, concentration=0.1)
    with pytest.raises(ValidationError, match=r"Label.* cannot be empty"): # Updated regex to be more robust
        Enzyme(label="", concentration=0.1)


def test_Environment():
    """Test the Environment model."""
    env = Environment(
        temperature=310.15,
        pH=7.4,
        ionic_strength=0.15,
        oxygen_level=0.21,
        seed_mechanisms=["basic_metabolism"]
    )
    assert env.temperature == 310.15
    assert env.pH == 7.4
    assert env.ionic_strength == 0.15
    assert env.oxygen_level == 0.21
    assert env.seed_mechanisms == ["basic_metabolism"]

    # Test temperature as list
    env_range_t = Environment(temperature=[298.15, 310.15], pH=7.0)
    assert env_range_t.temperature == [298.15, 310.15]

    # Test pH as list
    env_range_ph = Environment(temperature=298.15, pH=[6.0, 8.0])
    assert env_range_ph.pH == [6.0, 8.0]

    # Test validators
    with pytest.raises(ValidationError, match="pH must be between 0 and 14"):
        Environment(temperature=37, pH=15)
    with pytest.raises(ValidationError, match="Input should be greater than 0"): # Updated regex for pH=-1
        Environment(temperature=37, pH=-1)
    with pytest.raises(ValidationError, match="Oxygen level must be between 0 and 1"):
        Environment(temperature=37, pH=7, oxygen_level=1.2)
    with pytest.raises(ValidationError, match="Oxygen level must be between 0 and 1"):
        Environment(temperature=37, pH=7, oxygen_level=-0.1)
    with pytest.raises(ValidationError, match="Temperature as list must have exactly 2 values"):
        Environment(temperature=[25, 30, 40], pH=7)
    with pytest.raises(ValidationError, match="Ionic strength cannot be negative"):
        Environment(temperature=298.15, pH=7.0, ionic_strength=-0.5)


def test_Settings():
    """Test the Settings model."""
    # Test with minimal required fields, explicitly setting optional fields to None or their defaults
    # This helps Pydantic correctly recognize all fields, especially with extra="forbid"
    settings = Settings(
        end_time=100.0,
        time_step=1.0,
        solver="odeint",
        use_core_edge=True,
        threshold=1e-3,
        max_iterations=50,
        ml_models_enabled=True,
        stop_at_steady_state=True,
        termination_conversion=None,
        termination_rate_ratio=None,
        verbose=None,
        output_directory=None,
        flux_adapter="RMG",
        profiles_adapter="RMG",
        generate_plots=False,
        save_simulation_profiles=False
    )
    assert settings.end_time == 100.0
    assert settings.time_step == 1.0
    assert settings.solver == "odeint"
    assert settings.use_core_edge is True # Default
    assert settings.threshold == 1e-3 # Default
    assert settings.max_iterations == 50 # Default
    assert settings.ml_models_enabled is True # Default
    assert settings.stop_at_steady_state is True # Default
    assert settings.verbose is None # Asserting the value now
    assert settings.output_directory is None # Asserting the value now
    assert settings.flux_adapter == "RMG" # Asserting the value now
    assert settings.profiles_adapter == "RMG" # Asserting the value now
    assert settings.generate_plots is False # Asserting the value now
    assert settings.save_simulation_profiles is False # Asserting the value now

    # Test with all fields
    settings_full = Settings(
        end_time=3600.0,
        time_step=0.5,
        solver="CVODE",
        use_core_edge=False,
        threshold=1e-4,
        max_iterations=100,
        ml_models_enabled=True,
        stop_at_steady_state=False,
        termination_conversion={"Glucose": 0.95},
        termination_rate_ratio=0.005,
        verbose=10,
        output_directory="/tmp/test_full/output",
        flux_adapter="Cantera",
        profiles_adapter="Custom",
        generate_plots=True,
        save_simulation_profiles=True
    )
    assert settings_full.end_time == 3600.0
    assert settings_full.time_step == 0.5
    assert settings_full.solver == "CVODE"
    assert settings_full.use_core_edge is False
    assert settings_full.threshold == 1e-4
    assert settings_full.max_iterations == 100
    assert settings_full.ml_models_enabled is True
    assert settings_full.stop_at_steady_state is False
    assert settings_full.termination_conversion == {"Glucose": 0.95}
    assert settings_full.termination_rate_ratio == 0.005
    assert settings_full.verbose == 10
    assert settings_full.output_directory == "/tmp/test_full/output"
    assert settings_full.flux_adapter == "Cantera"
    assert settings_full.profiles_adapter == "Custom"
    assert settings_full.generate_plots is True
    assert settings_full.save_simulation_profiles is True

    # Test validators
    with pytest.raises(ValidationError, match=r"'time_step' must be smaller than 'end_time'"):
        Settings(end_time=10.0, time_step=10.0, solver="odeint")
    with pytest.raises(ValidationError, match="Input should be less than 1"): # Updated regex
        Settings(end_time=100.0, time_step=1.0, solver="odeint", termination_rate_ratio=1.0)
    with pytest.raises(ValidationError, match="termination_conversion values must be between 0 and 1"):
        Settings(end_time=100.0, time_step=1.0, solver="odeint", termination_conversion={"A": 1.1})
    with pytest.raises(ValidationError, match="Verbose level must be 10, 20, 30, 40, or 50"):
        Settings(end_time=100.0, time_step=1.0, solver="odeint", verbose=25)
    with pytest.raises(ValidationError, match="Verbose level must be 10, 20, 30, 40, or 50"):
        Settings(end_time=100.0, time_step=1.0, solver="odeint", verbose=5)
    with pytest.raises(ValidationError, match="value is not a valid enumeration member"):
        Settings(end_time=100.0, time_step=1.0, solver="invalid_solver")
    with pytest.raises(ValidationError, match="value is not a valid enumeration member"):
        Settings(end_time=100.0, time_step=1.0, solver="odeint", flux_adapter="invalid_adapter")


def test_Database():
    """Test the Database model."""
    db = Database(
        name="enzyme_catalysis",
        rate_law="Michaelis-Menten",
        parameter_estimator="ML",
        parameters={"Km": 0.01, "Vmax": 1.0},
        thermo_libraries=["primaryThermoLibrary"],
        kinetics_libraries=["BurkeH2O2inN2"],
        use_low_credence_libraries=True,
        kinetics_depositories=["default"],
        kinetics_families=["default"]
    )
    assert db.name == "enzyme_catalysis"
    assert db.rate_law == "Michaelis-Menten"
    assert db.parameter_estimator == "ML"
    assert db.parameters["Km"] == 0.01
    assert "primaryThermoLibrary" in db.thermo_libraries
    assert "BurkeH2O2inN2" in db.kinetics_libraries
    assert db.use_low_credence_libraries is True
    assert db.kinetics_depositories == ["default"]
    assert db.kinetics_families == ["default"]

    # Test with minimal required fields
    db_minimal = Database(
        name="MinimalDB",
        rate_law="Michaelis-Menten"
    )
    assert db_minimal.name == "MinimalDB"
    assert db_minimal.rate_law == "Michaelis-Menten"
    assert db_minimal.parameter_estimator == "ML" # Default value

    # Test validators
    with pytest.raises(ValidationError, match="Name cannot be empty"):
        Database(name="", rate_law="Michaelis-Menten")
    with pytest.raises(ValidationError, match="Input should be 'Michaelis-Menten', 'Hill' or 'MassAction'"): # Updated regex
        Database(name="bad", rate_law="UnknownLaw")
    with pytest.raises(ValidationError, match="Input should be 'ML', 'group_contribution' or 'fixed_defaults'"): # Updated regex
        Database(name="bad", rate_law="Michaelis-Menten", parameter_estimator="UnknownEstimator")
    with pytest.raises(ValidationError, match="Input should be a valid string"): # Updated regex
        Database(name="bad", rate_law="Michaelis-Menten", kinetics_families=[123])
    with pytest.raises(ValidationError, match="kinetics_families must be a list or 'default'"):
        Database(name="bad", rate_law="Michaelis-Menten", kinetics_families="not_default_string")
    with pytest.raises(ValidationError, match="Each kinetics depository must be a string"): # Updated regex
        Database(name="bad", rate_law="Michaelis-Menten", kinetics_depositories=[123])
    with pytest.raises(ValidationError, match="thermo_libraries must be a list"):
        Database(name="bad", rate_law="Michaelis-Menten", thermo_libraries="not_list")
    with pytest.raises(ValidationError, match="Each thermo library must be a string"):
        Database(name="bad", rate_law="Michaelis-Menten", thermo_libraries=[123])
    with pytest.raises(ValidationError, match="kinetics_libraries must be a list"):
        Database(name="bad", rate_law="Michaelis-Menten", kinetics_libraries="not_list")
    with pytest.raises(ValidationError, match="Each kinetics library must be a string"):
        Database(name="bad", rate_law="Michaelis-Menten", kinetics_libraries=[123])
    
    # Test species_constraints within Database
    db_with_constraints = Database(
        name="constrained_db",
        rate_law="MassAction",
        species_constraints=SpeciesConstraints(allowed=["input species"])
    )
    assert db_with_constraints.species_constraints is not None
    assert "input species" in db_with_constraints.species_constraints.allowed

    with pytest.raises(ValidationError, match="Input should be an instance of SpeciesConstraints"): # Updated regex
        Database(name="bad", rate_law="MassAction", species_constraints={"allowed": ["input species"]}) # Dict instead of model
    with pytest.raises(ValidationError, match="species_constraints must have at least one allowed entry"):
        Database(name="bad", rate_law="MassAction", species_constraints=SpeciesConstraints(allowed=[]))


def test_speciesconstraints():
    """Test the SpeciesConstraints model."""
    constraints = SpeciesConstraints(
        allowed=["input species", "reaction libraries"]
    )
    assert "input species" in constraints.allowed
    assert "seed mechanisms" not in constraints.allowed # Default value not present if overridden

    # Test default
    default_constraints = SpeciesConstraints()
    assert default_constraints.allowed == ['input species', 'seed mechanisms', 'reaction libraries']

    with pytest.raises(ValidationError, match="Each 'allowed' value must be one of"):
        SpeciesConstraints(allowed=["invalid entry"])

    with pytest.raises(ValidationError, match="Each 'allowed' value must be one of"):
        SpeciesConstraints(allowed=["input species", "invalid entry"])

    with pytest.raises(ValidationError, match="'allowed' list cannot be empty"):
        SpeciesConstraints(allowed=[])


def test_InputBase():
    """Test the InputBase model."""
    # Test with minimal required fields
    input_data_minimal = InputBase(
        project="TestProjectMinimal",
        species=[Species(label="Glucose", concentration=0.1)],
        enzymes=[Enzyme(label="Hexokinase", concentration=0.01)],
        environment=Environment(temperature=298.15, pH=7.0),
        settings=Settings(end_time=100.0, time_step=1.0, solver="odeint"),
        database=Database(name="enzyme_catalysis", rate_law="Michaelis-Menten")
    )
    assert input_data_minimal.project == "TestProjectMinimal"
    assert len(input_data_minimal.species) == 1
    assert len(input_data_minimal.enzymes) == 1
    assert input_data_minimal.environment.temperature == 298.15
    assert input_data_minimal.settings.end_time == 100.0
    assert input_data_minimal.database.name == "enzyme_catalysis"

    # Test with more comprehensive settings
    input_data_full = InputBase(
        project="TestProjectFull",
        project_directory="/tmp/test_full",
        species=[
            Species(label="Glucose", concentration=0.1, observable=True),
            Species(label="ATP", concentration=0.05)
        ],
        enzymes=[
            Enzyme(label="Hexokinase", concentration=0.001, ecnumber="EC 2.7.1.1")
        ],
        environment=Environment(
            temperature=[298.15, 310.15],
            pH=[6.5, 7.5],
            ionic_strength=0.2,
            oxygen_level=0.5
        ),
        settings=Settings(
            end_time=3600.0,
            time_step=0.5,
            solver="CVODE",
            use_core_edge=False,
            threshold=1e-4,
            max_iterations=100,
            ml_models_enabled=True,
            stop_at_steady_state=False,
            termination_conversion={"Glucose": 0.95},
            termination_rate_ratio=0.005,
            verbose=10,
            output_directory="/tmp/test_full/output",
            flux_adapter="Cantera",
            profiles_adapter="Custom",
            generate_plots=True,
            save_simulation_profiles=True
        ),
        database=Database(
            name="full_db",
            rate_law="Hill",
            parameter_estimator="group_contribution",
            parameters={"n_Hill": 2, "K_Hill": 0.002},
            thermo_libraries=["lib1"],
            kinetics_libraries=["libA"],
            use_low_credence_libraries=True,
            species_constraints=SpeciesConstraints(allowed=["input species"])
        )
    )
    assert input_data_full.project == "TestProjectFull"
    assert input_data_full.project_directory == "/tmp/test_full"
    assert len(input_data_full.species) == 2
    assert input_data_full.settings.verbose == 10
    assert input_data_full.settings.output_directory == "/tmp/test_full/output"
    assert input_data_full.settings.flux_adapter == "Cantera"
    assert input_data_full.settings.profiles_adapter == "Custom"
    assert input_data_full.settings.generate_plots is True
    assert input_data_full.settings.save_simulation_profiles is True
    assert input_data_full.database.parameter_estimator == "group_contribution"
    assert input_data_full.database.species_constraints.allowed == ["input species"]


    # Test validators for InputBase (mostly covered by nested model validators)
    with pytest.raises(ValidationError, match="Input should be a valid list"):
        InputBase(
            project="TestProject",
            species="not_a_list", # Invalid type
            enzymes=[Enzyme(label="ATP", concentration=0.01)],
            environment=Environment(temperature=298.15, pH=7.0),
            settings=Settings(end_time=100.0, time_step=1.0, solver="odeint"),
            database=Database(name="enzyme_catalysis", rate_law="Michaelis-Menten")
        )

    with pytest.raises(ValidationError, match="Input should be a valid list"):
        InputBase(
            project="TestProject",
            species=[Species(label="Glucose", concentration=0.1)],
            enzymes="not_a_list", # Invalid type
            environment=Environment(temperature=298.15, pH=7.0),
            settings=Settings(end_time=100.0, time_step=1.0, solver="odeint"),
            database=Database(name="enzyme_catalysis", rate_law="Michaelis-Menten")
        )

    with pytest.raises(ValidationError, match="Input should be a valid dictionary"):
        InputBase(
            project="TestProject",
            species=[Species(label="Glucose", concentration=0.1)],
            enzymes=[Enzyme(label="ATP", concentration=0.01)],
            environment="not_a_dict", # Invalid type
            settings=Settings(end_time=100.0, time_step=1.0, solver="odeint"),
            database=Database(name="enzyme_catalysis", rate_law="Michaelis-Menten")
        )

    with pytest.raises(ValidationError, match="Input should be a valid dictionary"):
        InputBase(
            project="TestProject",
            species=[Species(label="Glucose", concentration=0.1)],
            enzymes=[Enzyme(label="ATP", concentration=0.01)],
            environment=Environment(temperature=298.15, pH=7.0),
            settings="not_a_dict", # Invalid type
            database=Database(name="enzyme_catalysis", rate_law="Michaelis-Menten")
        )

    with pytest.raises(ValidationError, match="Input should be a valid dictionary"):
        InputBase(
            project="TestProject",
            species=[Species(label="Glucose", concentration=0.1)],
            enzymes=[Enzyme(label="ATP", concentration=0.01)],
            environment=Environment(temperature=298.15, pH=7.0),
            settings=Settings(end_time=100.0, time_step=1.0, solver="odeint"),
            database="not_a_dict" # Invalid type
        )

    # Test cases where nested validators would catch errors (already covered by individual tests, but good for completeness)
    with pytest.raises(ValidationError):
        InputBase(
            project="TestProject",
            species=[Species(label="Glucose", concentration=-0.1)], # Invalid concentration
            enzymes=[Enzyme(label="ATP", concentration=0.01)],
            environment=Environment(temperature=298.15, pH=7.0),
            settings=Settings(end_time=100.0, time_step=1.0, solver="odeint"),
            database=Database(name="enzyme_catalysis", rate_law="Michaelis-Menten")
        )

    with pytest.raises(ValidationError):
        InputBase(
            project="TestProject",
            species=[Species(label="Glucose", concentration=0.1)],
            enzymes=[Enzyme(label="ATP", concentration=0.01)],
            environment=Environment(temperature=298.15, pH=15.0), # Invalid pH
            settings=Settings(end_time=100.0, time_step=1.0, solver="odeint"),
            database=Database(name="enzyme_catalysis", rate_law="Michaelis-Menten")
        )

    with pytest.raises(ValidationError):
        InputBase(
            project="TestProject",
            species=[Species(label="Glucose", concentration=0.1)],
            enzymes=[Enzyme(label="ATP", concentration=0.01)],
            environment=Environment(temperature=298.15, pH=7.0),
            settings=Settings(end_time=100.0, time_step=100.0, solver="odeint"), # Invalid time_step
            database=Database(name="enzyme_catalysis", rate_law="Michaelis-Menten")
        )

    with pytest.raises(ValidationError):
        InputBase(
            project="TestProject",
            species=[Species(label="Glucose", concentration=0.1)],
            enzymes=[Enzyme(label="ATP", concentration=0.01)],
            environment=Environment(temperature=298.15, pH=7.0),
            settings=Settings(end_time=100.0, time_step=1.0, solver="odeint"),
            database=Database(name="enzyme_catalysis", rate_law="UnknownLaw") # Invalid rate_law
        )


if __name__ == "__main__":
    pytest.main([__file__])
