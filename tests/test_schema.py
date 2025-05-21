#!/usr/bin/env python3
# encoding: utf-8

"""
Test schema validation for BEES input models
"""

import pytest
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


def test_Species():
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

    species = Species(
        label="Fructose",
        concentration=(0.1, 1.0),
    )
    assert species.concentration == (0.1, 1.0)


def test_Enzyme():
    enzyme = Enzyme(
        label="ATP",
        concentration=0.01,
        charge=-3,
        reactive=True,
        constant=False,
        observable=True,
    )
    assert enzyme.label == "ATP"
    assert enzyme.concentration == 0.01
    assert enzyme.reactive is True

    enzyme = Enzyme(
        label="Glucose",
        concentration=(0.1, 1.0),
    )
    assert enzyme.concentration == (0.1, 1.0)

    with pytest.raises(ValidationError):
        Enzyme(label="SameRange", concentration=(0.5, 0.5))

    with pytest.raises(ValidationError):
        Enzyme(label="RangeConstant", concentration=(0.1, 0.9), constant=True)

    with pytest.raises(ValidationError):
        Enzyme(label="NegativeConcentration", concentration=-0.1)

    with pytest.raises(ValidationError):
        Enzyme(label="NegativeRange", concentration=(-0.1, 0.5))

    with pytest.raises(ValidationError):
        Enzyme(label="InvalidEC", ecnumber="not.a.valid.ecnumber")


def test_Environment():
    env = Environment(
        temperature=[25, 37],
        pH=7.0,
        ionic_strength=0.1,
        oxygen_level=0.8,
        seed_mechanisms=["basic_metabolism"]
    )
    assert env.temperature == [25, 37]
    assert env.pH == 7.0

    with pytest.raises(ValidationError):
        Environment(temperature=37, pH=15, seed_mechanisms=[])

    with pytest.raises(ValidationError):
        Environment(temperature=37, pH=7, oxygen_level=1.2, seed_mechanisms=[])

    with pytest.raises(ValidationError):
        Environment(temperature=[25, 30, 40], pH=7, seed_mechanisms=[])


def test_settings():
    model = Settings(
        end_time=100.0,
        time_step=1.0,
        solver="CVODE",
        threshold=0.001,
        max_iterations=50,
        stop_at_steady_state=True
    )
    assert model.solver == "CVODE"

    with pytest.raises(ValidationError):
        Settings(end_time=10, time_step=10, solver="CVODE")

    with pytest.raises(ValidationError):
        Settings(
            end_time=100,
            time_step=1,
            solver="odeint",
            termination_rate_ratio=1.0
        )


def test_Database():
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
        Database(
            name="bad",
            rate_law="UnknownLaw",
            parameter_estimator="ML"
        )


def test_speciesconstraints():
    constraints = SpeciesConstraints(
        allowed=["input species", "reaction libraries"]
    )
    assert "input species" in constraints.allowed

    with pytest.raises(ValidationError):
        SpeciesConstraints(allowed=["invalid entry"])

    with pytest.raises(ValidationError):
        SpeciesConstraints(allowed=["input species", "invalid entry"])

    with pytest.raises(ValidationError):
        SpeciesConstraints(allowed=[])


def test_InputBase():
    input_data = InputBase(
        project="TestProject",
        species=[Species(label="Glucose", concentration=0.1)],
        enzymes=[Enzyme(label="ATP", concentration=0.01)],
        environment=Environment(temperature=[25, 37], pH=7.0),
        settings=Settings(end_time=100.0, time_step=1.0, solver="odeint"),
        database=Database(name="enzyme_catalysis", rate_law="Michaelis-Menten")
    )
    assert input_data.project == "TestProject"

    with pytest.raises(ValidationError):
        InputBase(
            project="TestProject",
            species=[Species(label="Glucose", concentration=-0.1)],
            enzymes=[Enzyme(label="ATP", concentration=0.01)],
            environment=Environment(temperature=[25, 37], pH=7.0),
            settings=Settings(end_time=100.0, time_step=1.0, solver="odeint"),
            database=Database(name="enzyme_catalysis", rate_law="Michaelis-Menten")
        )

    with pytest.raises(ValidationError):
        InputBase(
            project="TestProject",
            species=[Species(label="Glucose", concentration=0.1)],
            enzymes=[Enzyme(label="ATP", concentration=-0.01)],
            environment=Environment(temperature=[25, 37], pH=7.0),
            settings=Settings(end_time=100.0, time_step=1.0, solver="odeint"),
            database=Database(name="enzyme_catalysis", rate_law="Michaelis-Menten")
        )

    with pytest.raises(ValidationError):
        InputBase(
            project="TestProject",
            species=[Species(label="Glucose", concentration=0.1)],
            enzymes=[Enzyme(label="ATP", concentration=0.01)],
            environment=Environment(temperature=[25, 37], pH=15.0),
            settings=Settings(end_time=100.0, time_step=1.0, solver="odeint"),
            database=Database(name="enzyme_catalysis", rate_law="Michaelis-Menten")
        )

    with pytest.raises(ValidationError):
        InputBase(
            project="TestProject",
            species=[Species(label="Glucose", concentration=0.1)],
            enzymes=[Enzyme(label="ATP", concentration=0.01)],
            environment=Environment(temperature=[25, 37], pH=7.0),
            settings=Settings(end_time=100.0, time_step=100.0, solver="odeint"),
            database=Database(name="enzyme_catalysis", rate_law="Michaelis-Menten")
        )

    with pytest.raises(ValidationError):
        InputBase(
            project="TestProject",
            species=[Species(label="Glucose", concentration=0.1)],
            enzymes=[Enzyme(label="ATP", concentration=0.01)],
            environment=Environment(temperature=[25, 37], pH=7.0),
            settings=Settings(end_time=100.0, time_step=1.0, solver="odeint"),
            database=Database(name="enzyme_catalysis", rate_law="UnknownLaw")
        )


if __name__ == "__main__":
    pytest.main([__file__])
