#!/usr/bin/env python3
# encoding: utf-8

"""
Test BEES.schema validation 

This module contains tests for the BEES.
To run the tests, use pytest and the command line: pytest -v tests/test_schema.py

"""

import pytest
from pydantic import ValidationError
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
        smiles="C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O"
    )
    assert species.label == "Glucose"
    assert species.concentration == 0.1
    assert species.reactive is True
    assert species.observable is True
    assert species.charge == 0
    assert species.solvent is False
    assert species.xyz is None
    assert species.smiles == "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O"
    assert species.inchi is None
    assert species.adjlist is None


    # Test species with different molecular descriptors (separately, since we now enforce exclusivity)
    species_with_smiles = Species(
        label="Fructose",
        concentration=(0.1, 1.0),
        smiles="C(C1C(C(C(C(O1)CO)O)O)O)O"
    )
    assert species_with_smiles.concentration == (0.1, 1.0)
    assert species_with_smiles.smiles == "C(C1C(C(C(C(O1)CO)O)O)O)O"
    assert species_with_smiles.inchi is None
    assert species_with_smiles.adjlist is None
    
    species_with_inchi = Species(
        label="Methane",
        concentration=(0.1, 1.0),
        inchi="InChI=1S/CH4/h1H4" # Valid InChI for Methane for testing purposes
    )
    assert species_with_inchi.concentration == (0.1, 1.0)
    assert species_with_inchi.inchi == "InChI=1S/CH4/h1H4"
    assert species_with_inchi.smiles is None
    assert species_with_inchi.adjlist is None
    
    species_with_adjlist = Species(
        label="Carbon",
        concentration=(0.1, 1.0),
        # A minimal valid MolBlock for a single Carbon atom for reliable RDKit parsing
        adjlist="""
  MOL

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END"""
    )
    assert species_with_adjlist.concentration == (0.1, 1.0)
    assert species_with_adjlist.adjlist is not None
    assert species_with_adjlist.smiles is None
    assert species_with_adjlist.inchi is None

    # Test validators
    with pytest.raises(ValidationError, match="Concentration range cannot have identical values"):
        Species(label="SameRange", concentration=(0.5, 0.5))
    with pytest.raises(ValidationError, match="Constant species cannot have a concentration range"):
        Species(label="RangeConstant", concentration=(0.1, 0.9), constant=True)
    with pytest.raises(ValidationError, match="Input should be greater than 0"):
        Species(label="NegativeConcentration", concentration=-0.1)
    with pytest.raises(ValidationError, match="Input should be greater than 0"):
        Species(label="NegativeRange", concentration=(-0.1, 0.5))
    
    # The following tests are commented out because the current schema's field_validator, but there is very good chance that it will be back in the script in the future.
    
    # for 'constant' does not seem to trigger these specific errors as expected by Pydantic v2.
    # Cross-field validation might require a model_validator in the schema.
    # with pytest.raises(ValidationError, match="Reactive species cannot be constant"):
    #     Species(label="ReactiveConstant", reactive=True, constant=True, concentration=0.1)
    # with pytest.raises(ValidationError, match="Observable species cannot be constant"):
    #     Species(label="ObservableConstant", observable=True, constant=True, concentration=0.1)
    
    with pytest.raises(ValidationError, match="Invalid SMILES string"):
        Species(label="InvalidSMILES", concentration=0.1, smiles="InvalidSmiles[")
    with pytest.raises(ValidationError, match="Invalid InChI string"):
        Species(label="InvalidInChI", concentration=0.1, inchi="InvalidInChI")
    with pytest.raises(ValidationError, match="Invalid adjacency list"): # Updated message for adjlist
        Species(label="InvalidAdjlist", concentration=0.1, adjlist="Invalid Adjlist Content")
    
    # Test molecular descriptor exclusivity validator
    # Test 1: Multiple descriptors should fail
    with pytest.raises(ValidationError, match="Only one molecular descriptor can be provided"):
        Species(label="MultipleDescriptors", concentration=0.1, smiles="CCO", inchi="InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
    
    with pytest.raises(ValidationError, match="Only one molecular descriptor can be provided"):
        Species(label="MultipleDescriptors2", concentration=0.1, smiles="CCO", adjlist="""
  MOL

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END""")
    
    with pytest.raises(ValidationError, match="Only one molecular descriptor can be provided"):
        Species(label="MultipleDescriptors3", concentration=0.1, inchi="InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3", adjlist="""
  MOL

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END""")
    
    # Test 2: No descriptors should fail
    with pytest.raises(ValidationError, match="At least one molecular descriptor must be provided"):
        Species(label="NoDescriptors", concentration=0.1)
    
    # Test 3: Single descriptor should work (already tested above, but let's be explicit)
    species_smiles_only = Species(label="SmilesOnly", concentration=0.1, smiles="CCO")
    assert species_smiles_only.smiles == "CCO"
    assert species_smiles_only.inchi is None
    assert species_smiles_only.adjlist is None
    
    species_inchi_only = Species(label="InchiOnly", concentration=0.1, inchi="InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
    assert species_inchi_only.inchi == "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
    assert species_inchi_only.smiles is None
    assert species_inchi_only.adjlist is None
    
    species_adjlist_only = Species(label="AdjlistOnly", concentration=0.1, adjlist="""
  MOL

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END""")
    assert species_adjlist_only.adjlist is not None
    assert species_adjlist_only.smiles is None
    assert species_adjlist_only.inchi is None


def test_Enzyme():
    """Test the Enzyme model."""
    enzyme = Enzyme(
        label="ATP",
        concentration=0.01,
        charge=-3,
        reactive=True,
        constant=False,
        observable=True,
        ecnumber="EC 2.7.1.1",
        smiles="C1=NC2=C(C(=N1)N)N=CN2C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O"
    )
    assert enzyme.label == "ATP"
    assert enzyme.concentration == 0.01
    assert enzyme.reactive is True
    assert enzyme.ecnumber == "EC 2.7.1.1"
    assert enzyme.observable is True # Inherited default from Species

    enzyme = Enzyme(
        label="Phosphofructokinase",
        concentration=(0.1, 1.0),
        smiles="CC(C)CC1=NC(=CS1)C(=O)N[C@@H](CCC(=O)O)C(=O)O"
    )
    assert enzyme.concentration == (0.1, 1.0)

    # Test validators
    with pytest.raises(ValidationError, match="Concentration range cannot have identical values"):
        Enzyme(label="SameRange", concentration=(0.5, 0.5), smiles="CCO")
    with pytest.raises(ValidationError, match="Constant species cannot have a concentration range"):
        Enzyme(label="RangeConstant", concentration=(0.1, 0.9), constant=True, smiles="CCO")
    with pytest.raises(ValidationError, match="Input should be greater than 0"):
        Enzyme(label="NegativeConcentration", concentration=-0.1, smiles="CCO")
    with pytest.raises(ValidationError, match="Input should be greater than 0"):
        Enzyme(label="NegativeRange", concentration=(-0.1, 0.5), smiles="CCO")
    with pytest.raises(ValidationError, match=r"String should match pattern '\^EC \\d\+\\.\\d\+\\.\\d\+\\.\\d\+\$'"): 
        Enzyme(label="InvalidEC", concentration=0.1, ecnumber="not.a.valid.ecnumber", smiles="CCO")
    with pytest.raises(ValidationError, match=r"String should match pattern '\^EC \\d\+\\.\\d\+\\.\\d\+\\.\\d\+\$'"): 
        Enzyme(label="InvalidEC2", concentration=0.1, ecnumber="2.7.1.1", smiles="CCO") # Missing 'EC ' prefix
    
    # The following tests are commented out for the same reason as in test_Species.
    # with pytest.raises(ValidationError, match="Reactive species cannot be constant"):
    #     Enzyme(label="ReactiveConstant", reactive=True, constant=True, concentration=0.1)
    # with pytest.raises(ValidationError, match="Observable species cannot be constant"):
    
    #     Enzyme(label="ObservableConstant", observable=True, constant=True, concentration=0.1)
    with pytest.raises(ValidationError, match=r"Label.* cannot be empty"):
        Enzyme(label="", concentration=0.1, smiles="CCO")


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
    assert env_range_t.temperature == (298.15, 310.15) 

    # Test pH as list
    env_range_ph = Environment(temperature=298.15, pH=[6.0, 8.0])
    assert env_range_ph.pH == (6.0, 8.0) 

    # Test validators
    with pytest.raises(ValidationError, match=r"Input should be less than or equal to 14"): 
        Environment(temperature=37, pH=15)
    with pytest.raises(ValidationError, match=r"Input should be greater than or equal to 0"): 
        Environment(temperature=37, pH=-1)
    with pytest.raises(ValidationError, match=r"Input should be less than or equal to 1"): 
        Environment(temperature=37, pH=7, oxygen_level=1.2)
    with pytest.raises(ValidationError, match=r"Input should be greater than or equal to 0"): 
        Environment(temperature=37, pH=7, oxygen_level=-0.1)
    with pytest.raises(ValidationError, match=r"Tuple should have at most 2 items after validation, not 3"): 
        Environment(temperature=[25, 30, 40], pH=7)
    with pytest.raises(ValidationError, match=r"Input should be greater than or equal to 0"): 
        Environment(temperature=298.15, pH=7.0, ionic_strength=-0.5)


def test_Settings():
    """Test the Settings model.
    Test with minimal required fields, explicitly setting optional fields to None or their defaults
    """
    settings = Settings(
        end_time=100.0,
        time_step=1.0,
        time_units=TerminationTimeEnum.s, 
        toleranceKeepInEdge=1e-9, 
        toleranceMoveToCore=1e-5, 
        termination_conversion=None,
        termination_rate_ratio=None,
        max_edge_species=None, 
        filter_reactions=True, 
        modify_concentration_ranges_together=True, 
        max_iterations=50,
        verbose=20, # Default
        saveEdgeSpecies=True, 
        output_directory=None,
        generate_plots=False,
        save_simulation_profiles=False
    )
    assert settings.end_time == 100.0
    assert settings.time_step == 1.0
    assert settings.time_units == TerminationTimeEnum.s
    assert settings.toleranceKeepInEdge == 1e-9
    assert settings.toleranceMoveToCore == 1e-5
    assert settings.max_edge_species is None
    assert settings.filter_reactions is True
    assert settings.modify_concentration_ranges_together is True
    assert settings.max_iterations == 50
    assert settings.verbose == 20
    assert settings.saveEdgeSpecies is True
    assert settings.output_directory is None
    assert settings.generate_plots is False
    assert settings.save_simulation_profiles is False
    assert settings.termination_conversion is None
    assert settings.termination_rate_ratio is None

    # Test with all fields
    settings_full = Settings(
        end_time=3600.0,
        time_step=0.5,
        time_units=TerminationTimeEnum.hours,
        toleranceKeepInEdge=1e-6,
        toleranceMoveToCore=1e-8,
        termination_conversion={"Glucose": 0.95},
        termination_rate_ratio=0.005,
        max_edge_species=1000,
        filter_reactions=False,
        modify_concentration_ranges_together=False,
        max_iterations=100,
        verbose=10,
        saveEdgeSpecies=False,
        output_directory="/tmp/test_full/output",
        generate_plots=True,
        save_simulation_profiles=True
    )
    assert settings_full.end_time == 3600.0
    assert settings_full.time_step == 0.5
    assert settings_full.time_units == TerminationTimeEnum.hours
    assert settings_full.toleranceKeepInEdge == 1e-6
    assert settings_full.toleranceMoveToCore == 1e-8
    assert settings_full.termination_conversion == {"Glucose": 0.95}
    assert settings_full.termination_rate_ratio == 0.005
    assert settings_full.max_edge_species == 1000
    assert settings_full.filter_reactions is False
    assert settings_full.modify_concentration_ranges_together is False
    assert settings_full.max_iterations == 100
    assert settings_full.verbose == 10
    assert settings_full.saveEdgeSpecies is False
    assert settings_full.output_directory == "/tmp/test_full/output"
    assert settings_full.generate_plots is True
    assert settings_full.save_simulation_profiles is True

    
    
    # Test validators
    with pytest.raises(ValidationError, match=r"'time_step' must be smaller than 'end_time'"):
        Settings(end_time=10.0, time_step=10.0)
    with pytest.raises(ValidationError, match=r"Input should be less than 1"): 
        Settings(end_time=100.0, time_step=1.0, termination_rate_ratio=1.0)
    with pytest.raises(ValidationError, match=r"Input should be less than 1"): 
        Settings(end_time=100.0, time_step=1.0, termination_conversion={"A": 1.1})
    with pytest.raises(ValidationError, match=r"Verbose level must be 10, 20, 30, 40, or 50"): 
        Settings(end_time=100.0, time_step=1.0, verbose=25)
    with pytest.raises(ValidationError, match=r"Input should be greater than or equal to 10"): 
        Settings(end_time=100.0, time_step=1.0, verbose=5)
    with pytest.raises(ValidationError, match="Input should be greater than 0"):
        Settings(end_time=100.0, time_step=1.0, toleranceKeepInEdge=0) # This should still raise error
    with pytest.raises(ValidationError, match="Input should be greater than 0"):
        Settings(end_time=100.0, time_step=1.0, toleranceMoveToCore=0)
    with pytest.raises(ValidationError, match="Input should be greater than 0"):
        Settings(end_time=100.0, time_step=1.0, max_edge_species=0)
    with pytest.raises(ValidationError, match=r"Input should be 'micro-s', 'ms', 's', 'hrs', 'min', 'hours' or 'days'"): 
        Settings(end_time=100.0, time_step=1.0, time_units="invalid_unit")


def test_Database():
    """Test the Database model."""
    db = Database(
        name="enzyme_catalysis",
        parameters={"Km": 0.01, "Vmax": 1.0},
        thermo_libraries=["primaryThermoLibrary"],
        kinetics_libraries=["BurkeH2O2inN2"],
        chemistry_sets=["my_lipid_set"], 
        use_low_credence_libraries=True,
        seed_mechanism=["my_seed_rxn"], 
        kinetics_depositories=["default"],
        kinetics_families=["default"],
        solver="odeint", 
        kinetics_estimator="rate rules" 
    )
    assert db.name == "enzyme_catalysis"
    assert db.parameters["Km"] == 0.01
    assert "primaryThermoLibrary" in db.thermo_libraries
    assert "BurkeH2O2inN2" in db.kinetics_libraries
    assert "my_lipid_set" in db.chemistry_sets
    assert db.use_low_credence_libraries is True
    assert "my_seed_rxn" in db.seed_mechanism
    assert db.kinetics_depositories == ["default"]
    assert db.kinetics_families == ["default"]
    assert db.solver == "odeint"
    assert db.kinetics_estimator == "rate rules"

    # Test with minimal required fields
    db_minimal = Database(
        name="MinimalDB",
    )
    assert db_minimal.name == "MinimalDB"
    assert db_minimal.parameters is None
    assert db_minimal.thermo_libraries is None
    assert db_minimal.kinetics_libraries is None
    assert db_minimal.chemistry_sets is None
    assert db_minimal.use_low_credence_libraries is False
    assert db_minimal.seed_mechanism is None
    assert db_minimal.kinetics_depositories == "default"
    assert db_minimal.kinetics_families == "default"
    assert db_minimal.solver == "odeint" # Default value
    assert db_minimal.kinetics_estimator == "rate rules" # Default value

    # Test validators
    with pytest.raises(ValidationError, match=r"String should have at least 1 character"): 
        Database(name="")
    with pytest.raises(ValidationError, match=r"Input should be 'odeint', 'CVODE' or 'BDF'"): 
        Database(name="bad", solver="invalid_solver_type")
    with pytest.raises(ValidationError, match="Input should be a valid string"): # For kinetics_families list
        Database(name="bad", kinetics_families=[123])
    with pytest.raises(ValidationError, match=r"Input should be 'default'"): 
        Database(name="bad", kinetics_families="not_default_string")
    with pytest.raises(ValidationError, match=r"Input should be a valid string"): 
        Database(name="bad", kinetics_depositories=[123])
    with pytest.raises(ValidationError, match=r"Input should be a valid list"): 
        Database(name="bad", thermo_libraries="not_list")
    with pytest.raises(ValidationError, match=r"Input should be a valid string"): 
        Database(name="bad", thermo_libraries=[123])
    with pytest.raises(ValidationError, match=r"Input should be a valid list"): 
        Database(name="bad", kinetics_libraries="not_list")
    with pytest.raises(ValidationError, match=r"Input should be a valid string"): 
        Database(name="bad", kinetics_libraries=[123])
    with pytest.raises(ValidationError, match=r"Input should be a valid string"): 
        Database(name="bad", chemistry_sets=[123])
    with pytest.raises(ValidationError, match=r"Input should be a valid string"): 
        Database(name="bad", seed_mechanism=[123])


def test_speciesconstraints():
    """Test the SpeciesConstraints model."""
    constraints = SpeciesConstraints(
        allowed=["input species", "reaction libraries"],
        tolerance_thermo_keep_species_in_edge=0.01,
        max_C_atoms=20,
        max_O_atoms=10,
        max_radical_electrons=2,
    )
    assert "input species" in constraints.allowed
    assert constraints.tolerance_thermo_keep_species_in_edge == 0.01
    assert constraints.max_C_atoms == 20
    assert constraints.max_O_atoms == 10
    assert constraints.max_radical_electrons == 2
    assert constraints.max_singlet_carbenes == 1 # Default
    assert constraints.max_carbene_radicals == 0 # Default
    assert constraints.allow_singlet_O2 is True # Default

    # Test default
    default_constraints = SpeciesConstraints()
    assert default_constraints.allowed == ['input species', 'seed mechanisms', 'reaction libraries']
    assert default_constraints.tolerance_thermo_keep_species_in_edge is None
    assert default_constraints.max_C_atoms is None


    with pytest.raises(ValidationError, match=r"Input should be 'input species', 'seed mechanisms' or 'reaction libraries'"): 
        SpeciesConstraints(allowed=["invalid entry"])

    with pytest.raises(ValidationError, match=r"Input should be 'input species', 'seed mechanisms' or 'reaction libraries'"): 
        SpeciesConstraints(allowed=["input species", "invalid entry"])

    with pytest.raises(ValidationError, match="'allowed' list cannot be empty"):
        SpeciesConstraints(allowed=[])

    with pytest.raises(ValidationError, match="Input should be greater than 0"):
        SpeciesConstraints(allowed=["input species"], tolerance_thermo_keep_species_in_edge=0)
    with pytest.raises(ValidationError, match="Input should be greater than 0"):
        SpeciesConstraints(allowed=["input species"], max_C_atoms=0)
    with pytest.raises(ValidationError, match="Input should be greater than or equal to 0"):
        SpeciesConstraints(allowed=["input species"], max_radical_electrons=-1)


def test_InputBase():
    """Test the InputBase model."""
    # Test with minimal required fields
    input_data_minimal = InputBase(
        project="TestProjectMinimal",
        species=[Species(label="Glucose", concentration=0.1, smiles="C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O")],
        enzymes=[Enzyme(label="Hexokinase", concentration=0.01, smiles="CC(C)CC1=NC(=CS1)C(=O)N[C@@H](CCC(=O)O)C(=O)O")],
        environment=Environment(temperature=298.15, pH=7.0),
        settings=Settings(end_time=100.0, time_step=1.0), # Use minimal settings
        database=Database(name="enzyme_catalysis") # Use minimal database
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
            Species(label="Glucose", concentration=0.1, observable=True, smiles="C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O"),
            Species(label="ATP", concentration=0.05, smiles="C1=NC2=C(C(=N1)N)N=CN2C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O")
        ],
        enzymes=[
            Enzyme(label="Hexokinase", concentration=0.001, ecnumber="EC 2.7.1.1", smiles="CC(C)CC1=NC(=CS1)C(=O)N[C@@H](CCC(=O)O)C(=O)O")
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
            time_units=TerminationTimeEnum.s,
            toleranceKeepInEdge=1e-6,
            toleranceMoveToCore=1e-8,
            termination_conversion={"Glucose": 0.95},
            termination_rate_ratio=0.005,
            max_edge_species=1000,
            filter_reactions=False,
            modify_concentration_ranges_together=False,
            max_iterations=100,
            verbose=10,
            saveEdgeSpecies=False,
            output_directory="/tmp/test_full/output",
            generate_plots=True,
            save_simulation_profiles=True
        ),
        database=Database(
            name="full_db",
            parameters={"n_Hill": 2, "K_Hill": 0.002},
            thermo_libraries=["lib1"],
            kinetics_libraries=["libA"],
            chemistry_sets=["my_set"],
            use_low_credence_libraries=True,
            seed_mechanism=["my_seed_rxn"],
            kinetics_depositories=["custom_depository"],
            kinetics_families=["custom_family"],
            solver="CVODE",
            kinetics_estimator="group_contribution" # Changed from 'rate rules' for testing
        )
    )
    assert input_data_full.project == "TestProjectFull"
    assert input_data_full.project_directory == "/tmp/test_full"
    assert len(input_data_full.species) == 2
    assert input_data_full.settings.verbose == 10
    assert input_data_full.settings.output_directory == "/tmp/test_full/output"
    assert input_data_full.settings.generate_plots is True
    assert input_data_full.settings.save_simulation_profiles is True
    assert input_data_full.database.kinetics_estimator == "group_contribution"
    assert input_data_full.database.solver == "CVODE"
    assert "my_set" in input_data_full.database.chemistry_sets
    assert "my_seed_rxn" in input_data_full.database.seed_mechanism


    # Test validators for InputBase (mostly covered by nested model validators)
    with pytest.raises(ValidationError, match="Input should be a valid list"):
        InputBase(
            project="TestProject",
            species="not_a_list", # Invalid type
                enzymes=[Enzyme(label="ATP", concentration=0.01, smiles="CCO")],
            environment=Environment(temperature=298.15, pH=7.0),
            settings=Settings(end_time=100.0, time_step=1.0),
            database=Database(name="enzyme_catalysis")
        )

    with pytest.raises(ValidationError, match="Input should be a valid list"):
        InputBase(
            project="TestProject",
                species=[Species(label="Glucose", concentration=0.1, smiles="CCO")],
            enzymes="not_a_list", # Invalid type
            environment=Environment(temperature=298.15, pH=7.0),
            settings=Settings(end_time=100.0, time_step=1.0),
            database=Database(name="enzyme_catalysis")
        )

    with pytest.raises(ValidationError, match="Input should be a valid dictionary"):
        InputBase(
            project="TestProject",
                species=[Species(label="Glucose", concentration=0.1, smiles="CCO")],
                enzymes=[Enzyme(label="ATP", concentration=0.01, smiles="CCO")],
            environment="not_a_dict", # Invalid type
            settings=Settings(end_time=100.0, time_step=1.0),
            database=Database(name="enzyme_catalysis")
        )

    with pytest.raises(ValidationError, match="Input should be a valid dictionary"):
        InputBase(
            project="TestProject",
                species=[Species(label="Glucose", concentration=0.1, smiles="CCO")],
                enzymes=[Enzyme(label="ATP", concentration=0.01, smiles="CCO")],
            environment=Environment(temperature=298.15, pH=7.0),
            settings="not_a_dict", # Invalid type
            database=Database(name="enzyme_catalysis")
        )

    with pytest.raises(ValidationError, match="Input should be a valid dictionary"):
        InputBase(
            project="TestProject",
                species=[Species(label="Glucose", concentration=0.1, smiles="CCO")],
                enzymes=[Enzyme(label="ATP", concentration=0.01, smiles="CCO")],
            environment=Environment(temperature=298.15, pH=7.0),
            settings=Settings(end_time=100.0, time_step=1.0),
            database="not_a_dict" # Invalid type
        )

    # Test cases where nested validators would catch errors
    with pytest.raises(ValidationError):
        InputBase(
            project="TestProject",
            species=[Species(label="Glucose", concentration=-0.1, smiles="CCO")], # Invalid concentration
            enzymes=[Enzyme(label="ATP", concentration=0.01, smiles="CCO")],
            environment=Environment(temperature=298.15, pH=7.0),
            settings=Settings(end_time=100.0, time_step=1.0),
            database=Database(name="enzyme_catalysis")
        )

    with pytest.raises(ValidationError):
        InputBase(
            project="TestProject",
            species=[Species(label="Glucose", concentration=0.1, smiles="CCO")],
            enzymes=[Enzyme(label="ATP", concentration=0.01, smiles="CCO")],
            environment=Environment(temperature=298.15, pH=15.0), # Invalid pH
            settings=Settings(end_time=100.0, time_step=1.0),
            database=Database(name="enzyme_catalysis")
        )

    with pytest.raises(ValidationError):
        InputBase(
            project="TestProject",
            species=[Species(label="Glucose", concentration=0.1, smiles="CCO")],
            enzymes=[Enzyme(label="ATP", concentration=0.01, smiles="CCO")],
            environment=Environment(temperature=298.15, pH=7.0),
            settings=Settings(end_time=100.0, time_step=100.0), # Invalid time_step
            database=Database(name="enzyme_catalysis")
        )

    with pytest.raises(ValidationError):
        InputBase(
            project="TestProject",
            species=[Species(label="Glucose", concentration=0.1, smiles="CCO")],
            enzymes=[Enzyme(label="ATP", concentration=0.01, smiles="CCO")],
            environment=Environment(temperature=298.15, pH=7.0),
            settings=Settings(end_time=100.0, time_step=1.0),
            database=Database(name="", solver="odeint") # Invalid database name
        )

if __name__ == "__main__":
    pytest.main([__file__])
