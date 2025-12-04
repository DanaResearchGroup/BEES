"""
Main BEES application module for kinetic model generation and refinement.

This module provides the BEES class, which orchestrates the kinetic model generation pipeline.
It handles input validation, project initialization, logging setup, and execution coordination.

The module processes YAML input files containing species, enzymes, environmental conditions,
and simulation settings. It validates inputs against Pydantic schemas, initializes project
directories and logging, and executes the model generation workflow.

Required inputs:
    - Input dictionary with 'project', 'species', 'enzymes', 'environment', 'database', and 'settings'
    - Valid project directory path (created automatically if missing)

This is probably the most important module in the code.
 run the code by executing this script directly python bees.py --input_file + directory
(decribed in BEES.py as well)
"""


import os
import time
from typing import Any, Dict

import bees.common as common
from bees.logger import Logger
from bees.schema import InputBase

# Base paths
BEES_PATH = common.BEES_PATH
PROJECTS_BASE_PATH = common.PROJECTS_BASE_PATH


class BEES():
    """
    The main BEES application class.
    Orchestrates: setup, schema validation, logging, execution.
    Also, here is where the model generation will be started. right now it is just a placeholder.
    """

    def __init__(self, input_data: Dict[str, Any]):
        self.t0 = time.time()
        self.input_data = input_data
        self.project: str = input_data.get("project", "default_project")

        # Determine project/output directories
        if "project_directory" in input_data:
            self.project_directory = input_data["project_directory"]
            if not os.path.isabs(self.project_directory):
                self.project_directory = os.path.join(PROJECTS_BASE_PATH, self.project_directory)
            elif not self.project_directory.startswith(PROJECTS_BASE_PATH):
                self.project_directory = os.path.join(PROJECTS_BASE_PATH, self.project)
        else:
            self.project_directory = os.path.join(PROJECTS_BASE_PATH, self.project)

        if "output_directory" in input_data.get("settings", {}):
            specified_output_dir = input_data["settings"]["output_directory"]
            if os.path.isabs(specified_output_dir):
                if not specified_output_dir.startswith(self.project_directory):
                    raise ValueError(
                        f"Absolute output path '{specified_output_dir}' is not within "
                        f"the project directory '{self.project_directory}'."
                    )
                self.output_directory = specified_output_dir
            else:
                self.output_directory = os.path.join(self.project_directory, specified_output_dir)
        else:
            self.output_directory = self.project_directory

        # Create directories
        try:
            os.makedirs(self.project_directory, exist_ok=True)
            if self.output_directory != self.project_directory:
                os.makedirs(self.output_directory, exist_ok=True)
        except OSError as e:
            raise FileNotFoundError(f"Failed to create project/output directory: {e}")

        # Initialize logger
        self.logger = Logger(
            project_directory=self.project_directory,
            verbose=input_data.get("settings", {}).get("verbose"),
            t0=self.t0,
        )
        self.logger.info(f"BEES version: {common.VERSION}")
        git_branch = common.get_git_branch()
        commit_hash, commit_date = common.get_git_commit()
        self.logger.info(f"Git branch: {git_branch}")
        self.logger.info(f"Git commit: {commit_hash} ({commit_date})")
        self.logger.log_args(input_data)

        # Validate schema
        try:
            self.bees_object = InputBase(**input_data)
            self.logger.info("Input validated successfully against the schema.")
            verified_input_path = os.path.join(self.project_directory, "input.yml")
            common.save_yaml_file(
                verified_input_path, self.bees_object.model_dump(exclude_unset=True)
            )
            self.logger.info(f"Saving validated input to {verified_input_path}")
        except Exception as e:
            self.logger.error(f"Input validation error: {e}")
            self.logger.log_footer(success=False)
            raise ValueError(f"Invalid input parameters provided: {e}")

        self.logger.info(
            f"BEES project {self.project} initialized successfully in {common.time_lapse(self.t0)}."
        )

    def execute(self):
        """
        Execute the BEES kinetic model generation pipeline.
        
        Current functionality:
        - Logs project initialization and input summary (species count, enzymes count, temperature, database info)
        - Validates solver configuration and simulation end time settings
        - Placeholder for future model generation (currently just sleeps for 1 second)
        
        Future planned functionality:
        - Generate reaction networks from input species and enzymes
        - Apply kinetic models and rate laws from database
        - Perform parameter estimation for unknown kinetic parameters
        - Run kinetic simulations using specified solver
        - Generate output files and plots
        
        Returns:
            dict: Execution results (currently empty, will contain model data in future)
        """

        self.logger.info(f"Starting BEES execution for project '{self.project}'...")

        self.logger.info(
            f"Input has {len(self.bees_object.species)} species and {len(self.bees_object.enzymes)} enzymes."
        )
        self.logger.info(f"Environment Temperature: {self.bees_object.environment.temperature} K")

        if hasattr(self.bees_object.database, "rate_law"):
            self.logger.info(
                f"Using '{self.bees_object.database.rate_law}' rate law "
                f"from database '{self.bees_object.database.name}'."
            )
        else:
            self.logger.info(
                f"Database name: '{self.bees_object.database.name}'. No specific rate law defined."
            )

        if hasattr(self.bees_object.database, "parameter_estimator"):
            self.logger.info(
                f"Using parameter estimator: '{self.bees_object.database.parameter_estimator}'."
            )
        else:
            self.logger.info(
                f"No specific parameter estimator defined for database '{self.bees_object.database.name}'."
            )

        if hasattr(self.bees_object, "settings") and hasattr(self.bees_object.settings, "end_time"):
            self.logger.info(f"The simulation will run until {self.bees_object.settings.end_time} time units.")
        else:
            self.logger.warning("End time setting is not defined.")
        self.logger.info(f"Using solver: {self.bees_object.database.solver}")
        
        time.sleep(1)  # Simulated work

        self.logger.info(
            f"BEES execution for project '{self.project}' completed in {common.time_lapse(self.t0)}."
        )
        self.logger.log_footer(success=True)