#!/usr/bin/env python3
# encoding: utf-8

"""
Test schema validation for BEES input models
"""

import pytest
import schema
from pydantic import BaseModel, Field
from typing import List, Dict, Union, Optional
from pydantic import ValidationError
from schema import (
    BEESCompounds,
    BEESEnvironment,
    BEESModelSettings,
    BEESReactionRule,
    BEESInputBase,
    BEESSpeciesConstraints
)


def test_BEESCompounds():  
    """Test BEESCompounds schema, Creates a compound with a single float concentration: 0.01"""
    compound = BEESCompounds(
        label="ATP",
        type="cofactor",
        initial_concentration=0.01,
        charge=-3,
        reactive=True,
        constant=False,
        observable=True,
    )
    assert compound.label == "ATP"
    assert compound.type == "cofactor"
    assert compound.initial_concentration == 0.01
    assert compound.reactive is True

    # Tuple concentration range is allowed if valid
    compound = BEESCompounds(
        label="Glucose",
        type="substrate",
        initial_concentration=(0.1, 1.0),
    )
    assert compound.initial_concentration == (0.1, 1.0)

    # Invalid range: same values
    with pytest.raises(ValidationError):
        BEESCompounds(
            label="Test",
            type="substrate",
            initial_concentration=(0.5, 0.5),
        )

    # Constant with range should fail
    with pytest.raises(ValidationError):
        BEESCompounds(
            label="Test",
            type="substrate",
            initial_concentration=(0.1, 0.9),
            constant=True,
            observable=True,
        )


        # Valid SMILES and InChI
    compound = BEESCompounds(
        label="Ethanol",
        type="substrate",
        initial_concentration=0.1,
        structure_smiles="CCO",
        structure_inchi="InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
    )
    assert compound.structure_smiles == "CCO"
    assert compound.structure_inchi.startswith("InChI=")

    # Invalid SMILES
    with pytest.raises(ValidationError):
        BEESCompounds(
            label="BadSMILES",
            type="substrate",
            initial_concentration=0.1,
            structure_smiles="not-a-smiles"
        )

    # updtaetInvalid InChI
    with pytest.raises(ValidationError):
        BEESCompounds(
            label="BadInChI",
            type="substrate",
            initial_concentration=0.1,
            structure_inchi="not-an-inchi"
        )



def test_BEESEnvironment():
    """Test BEESEnvironment schema"""
    env = BEESEnvironment(
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
        BEESEnvironment(temperature=37, pH=15, seed_mechanisms=[])

    with pytest.raises(ValidationError):
        # Invalid oxygen level
        BEESEnvironment(temperature=37, pH=7, Oxygen_level=1.2, seed_mechanisms=[])

    with pytest.raises(ValidationError):
        # Invalid temperature list size
        BEESEnvironment(temperature=[25, 30, 40], pH=7, seed_mechanisms=[])

    


def test_bees_Modelsettings():
    """Test BEESModelSettings schema"""
    model = BEESModelSettings(
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
        BEESModelSettings(end_time=10, time_step=10, solver="CVODE")

    # termination rate must be < 1
    with pytest.raises(ValidationError):
        BEESModelSettings(
            end_time=100,
            time_step=1,
            solver="odeint",
            termination_rate_ratio=1.0
        )


def test_BeesReactionRule():
    """Test BEESReactionRule schema"""
    rule = BEESReactionRule(
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
        BEESReactionRule(
            name="bad",
            rate_law="UnknownLaw",
            parameter_estimator="ML"
        )


def test_speciesconstraints():
    """Test BEESSpeciesConstraints schema"""
    constraints = BEESSpeciesConstraints(
        allowed=["input species", "reaction libraries"]
    )
    assert "input species" in constraints.allowed

    with pytest.raises(ValidationError):
        BEESSpeciesConstraints(
            allowed=["invalid entry"]
        )
    # Invalid entry in allowed
    with pytest.raises(ValidationError):
        BEESSpeciesConstraints(
            allowed=["input species", "invalid entry"]
        )
    # Empty allowed list
    with pytest.raises(ValidationError):
        BEESSpeciesConstraints(
            allowed=[]
        )
    # Invalid entry in allowed              
    with pytest.raises(ValidationError):
        BEESSpeciesConstraints(
            allowed=["input species", "invalid entry"]
        )
    # Empty allowed list
    with pytest.raises(ValidationError):
        BEESSpeciesConstraints(
            allowed=[]
        )
    # Invalid entry in allowed  



def test_BEESInputBase():
    """Test BEESInputBase with minimal valid input"""
    env = BEESEnvironment(temperature=37, pH=7, seed_mechanisms=["base"])
    compound = BEESCompounds(label="H2O", type="solvent", initial_concentration=1.0)
    rule = BEESReactionRule(name="test", rate_law="MassAction", parameter_estimator="ML")
    model = BEESModelSettings(end_time=100.0, time_step=1.0, solver="odeint")
    
    bees_input = BEESInputBase(
        project="TestProject",
        compounds=[compound],
        environment=env,
        rules=[rule],
        model_settings=model,

    )
    assert bees_input.project == "TestProject"
    assert bees_input.environment.pH == 7
    assert bees_input.compounds[0].label == "H2O"
    assert bees_input.rules[0].name == "test"
    assert bees_input.model_settings.end_time == 100.0

