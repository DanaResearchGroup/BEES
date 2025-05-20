#!/usr/bin/env python3
# encoding: utf-8

"""
Test schema validation for BEES input models
"""

import pytest
from pydantic import BaseModel, Field
from typing import List, Dict, Union, Optional
from pydantic import ValidationError
from src.schema import (
    Enzyme,
    Species,
    Environment,
    Settings,
    Database,
    InputBase,
    SpeciesConstraints
)


def test_Species(): # TODO: add test for species
    """Test Species schema"""       
    species = Species(
        label="Glucose",
        type="substrate",
        concentration=0.1,
        charge=0,
        reactive=True,
        constant=False,
        observable=True,
    )
    assert species.label == "Glucose"
    assert species.type == "substrate"
    assert species.concentration == 0.1
    assert species.reactive is True

    # Tuple concentration range is allowed if valid
    species = Species(
        label="Fructose",
        type="substrate",
        concentration=(0.1, 1.0),
    )
    assert species.concentration == (0.1, 1.0)

def test_Enzyme():  
    """Test Enzyme schema, Creates a compound with a single float concentration: 0.01"""
    compound = Enzyme(
        label="ATP",
        type="cofactor",
        concentration=0.01,
        charge=-3,
        reactive=True,
        constant=False,
        observable=True,
    )
    assert compound.label == "ATP"
    assert compound.type == "cofactor"
    assert compound.concentration == 0.01
    assert compound.reactive is True

    # Tuple concentration range is allowed if valid
    compound = Enzyme(
        label="Glucose",
        type="substrate",
        concentration=(0.1, 1.0),
    )
    assert compound.concentration == (0.1, 1.0)

    # Invalid range: same values
    with pytest.raises(ValidationError):
        Enzyme(
            label="Test",
            type="substrate",
            concentration=(0.5, 0.5),
        )

    # Constant with range should fail
    with pytest.raises(ValidationError):
        Enzyme(
            label="Test",
            type="substrate",
            concentration=(0.1, 0.9),
            constant=True,
            observable=True,
        )


        # Valid SMILES and InChI
    compound = Enzyme(
        label="Ethanol",
        type="substrate",
        concentration=0.1,
        structure_smiles="CCO",
        structure_inchi="InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
    )
    assert compound.structure_smiles == "CCO"
    assert compound.structure_inchi.startswith("InChI=")

    # Invalid SMILES
    with pytest.raises(ValidationError):
        Enzyme(
            label="BadSMILES",
            type="substrate",
            concentration=0.1,
            structure_smiles="not-a-smiles"
        )

    # updtaetInvalid InChI
    with pytest.raises(ValidationError):
        Enzyme(
            label="BadInChI",
            type="substrate",
            concentration=0.1,
            structure_inchi="not-an-inchi"
        )



def test_Environment():
    """Test Environment schema"""
    env = Environment(
        temperature=[25, 37],
        pH=7.0,
        ionic_strength=0.1,
        Oxygen_level=0.8,
        seed_mechanisms=["basic_metabolism"]
    )
    assert env.temperature == [25, 37]
    assert env.pH == 7.0

    with pytest.raises(ValidationError):
        # Invalid pH
        Environment(temperature=37, pH=15, seed_mechanisms=[])

    with pytest.raises(ValidationError):
        # Invalid oxygen level
        Environment(temperature=37, pH=7, Oxygen_level=1.2, seed_mechanisms=[])

    with pytest.raises(ValidationError):
        # Invalid temperature list size
        Environment(temperature=[25, 30, 40], pH=7, seed_mechanisms=[])

    


def test_settings():
    """Test Settings schema"""
    model = Settings(
        end_time=100.0,
        time_step=1.0,
        solver="CVODE",
        threshold=0.001,
        max_iterations=50,
        stop_at_steady_state=True
    )
    assert model.solver == "CVODE"

    # time_step must be < end_time
    with pytest.raises(ValidationError):
        Settings(end_time=10, time_step=10, solver="CVODE")

    # termination rate must be < 1
    with pytest.raises(ValidationError):
        Settings(
            end_time=100,
            time_step=1,
            solver="odeint",
            termination_rate_ratio=1.0
        )


def test_Database():
    """Test Database schema"""
    rule = Database(
        name="enzyme_catalysis",
        rate_law="Michaelis-Menten",
        parameter_estimator="ML",
        parameters={"Km": 0.01, "Vmax": 1.0}
    )
    assert rule.rate_law == "Michaelis-Menten"
    assert rule.parameter_estimator == "ML"
    assert rule.parameters["Km"] == 0.01

    with pytest.raises(ValidationError):
        # Invalid rate law
        Database(
            name="bad",
            rate_law="UnknownLaw",
            parameter_estimator="ML"
        )


def test_speciesconstraints():
    """Test SpeciesConstraints schema"""
    constraints = SpeciesConstraints(
        allowed=["input species", "reaction libraries"]
    )
    assert "input species" in constraints.allowed

    with pytest.raises(ValidationError):
        SpeciesConstraints(
            allowed=["invalid entry"]
        )
    # Invalid entry in allowed
    with pytest.raises(ValidationError):
        SpeciesConstraints(
            allowed=["input species", "invalid entry"]
        )
    # Empty allowed list
    with pytest.raises(ValidationError):
        SpeciesConstraints(
            allowed=[]
        )
    # Invalid entry in allowed              
    with pytest.raises(ValidationError):
        SpeciesConstraints(
            allowed=["input species", "invalid entry"]
        )
    # Empty allowed list
    with pytest.raises(ValidationError):
        SpeciesConstraints(
            allowed=[]
        )
    # Invalid entry in allowed  





def test_InputBase():
    """Test InputBase with minimal valid input"""
    env = Environment(temperature=37, pH=7, seed_mechanisms=["base"])
    compound = Enzyme(label="H2O", type="solvent", concentration=1.0)
    rule = Database(name="test", rate_law="MassAction", parameter_estimator="ML")
    model = Settings(end_time=100.0, time_step=1.0, solver="odeint")
    
    input = InputBase(
        project="TestProject",
        compounds=[compound],
        environment=env,
        rules=[rule],
        model_settings=model,

    )
    assert input.project == "TestProject"
    assert input.environment.pH == 7
    assert input.compounds[0].label == "H2O"
    assert input.rules[0].name == "test"
    assert input.model_settings.end_time == 100.0

