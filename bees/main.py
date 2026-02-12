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
from bees.model_generator import ModelGenerator

# Base paths
BEES_PATH = common.BEES_PATH


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

        # 1. Setup Directories
        self._setup_directories()

        # 2. Initialize Logger
        self._setup_logging()

        # 3. Log System Info
        self._log_system_info()

        # 4. Load Chemical Ontology
        self._load_ontology()

        # 5. Fetch missing data from DB
        self._fetch_missing_enzyme_data()

        # 6. Validate Schema
        self._validate_input()

    def _setup_directories(self):
        """Determine and create project/output directories."""
        # 1. Determine the base directory for the project
        if "project_directory" in self.input_data:
            self.base_directory = self.input_data["project_directory"]
            if not os.path.isabs(self.base_directory):
                # If relative, it's relative to the installation root
                self.base_directory = os.path.normpath(os.path.join(common.BEES_PATH, self.base_directory))
        else:
            # Default to projects/<project_name>
            self.base_directory = os.path.join(common.BEES_PATH, "projects", self.project)

        # 2. Set project_directory to 'output' inside base_directory
        # This is where logs and output files will go by default
        self.project_directory = os.path.join(self.base_directory, "output")

        specified_output_dir = self.input_data.get("settings", {}).get("output_directory")
        if specified_output_dir:
            if os.path.isabs(specified_output_dir):
                if not specified_output_dir.startswith(self.base_directory):
                    # We still want to allow absolute paths outside, but maybe warn?
                    # For now, let's just use it if absolute.
                    self.output_directory = specified_output_dir
                else:
                    self.output_directory = specified_output_dir
            else:
                self.output_directory = os.path.join(self.base_directory, specified_output_dir)
        else:
            self.output_directory = self.project_directory

        try:
            os.makedirs(self.project_directory, exist_ok=True)
            if self.output_directory != self.project_directory:
                os.makedirs(self.output_directory, exist_ok=True)
        except OSError as e:
            raise FileNotFoundError(f"Failed to create output directory: {e}")

    def _setup_logging(self):
        """Initialize the BEES logger."""
        self.logger = Logger(
            project_directory=self.project_directory,
            verbose=self.input_data.get("settings", {}).get("verbose"),
            t0=self.t0,
        )

    def _log_system_info(self):
        """Log BEES version and Git information."""
        self.logger.info(f"BEES version: {common.VERSION}")
        git_branch = common.get_git_branch()
        commit_hash, commit_date = common.get_git_commit()
        self.logger.info(f"Git branch: {git_branch}")
        self.logger.info(f"Git commit: {commit_hash} ({commit_date})")
        self.logger.log_args(self.input_data)

    def _load_ontology(self):
        """Load and invert the chemical ontology for alias matching."""
        ontology_path = os.path.join(BEES_PATH, "db", "ontology.yaml")
        self.ontology = common.load_and_invert_ontology(ontology_path)
        if self.ontology:
            self.logger.info(f"Loaded chemical ontology from {ontology_path} ({len(self.ontology)} aliases)")

    def _fetch_missing_enzyme_data(self):
        """Fill in missing enzyme attributes from the database before validation."""
        # Function-first mode: reaction generation is based on EC numbers + substrates.
        # Amino acid sequences (if provided) are used only for kinetic estimation, not for discovery.
        return

    def _validate_input(self):
        """Validate input data against Pydantic schema."""
        try:
            self.bees_object = InputBase(**self.input_data)
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

        # Log input summary
        reactive_species = [s for s in self.bees_object.species if s.reactive]
        cofactor_species = [s for s in self.bees_object.species if not s.reactive]
        reactive_enzymes = [e for e in self.bees_object.enzymes if e.reactive]
        self.logger.info(
            f"Input: {len(reactive_species)} reactive species, "
            f"{len(cofactor_species)} cofactors, "
            f"{len(reactive_enzymes)} enzymes"
        )
        self.logger.info(f"Environment: T={self.bees_object.environment.temperature} K, "
                        f"pH={self.bees_object.environment.pH}")
        self.logger.info(f"Database: '{self.bees_object.database.name}'")

        # Log estimation settings
        est_enabled = getattr(self.bees_object.settings, "estimate_kinetics", False)
        est_backend = getattr(self.bees_object.settings, "kinetics_estimator", None)
        smiles_mode = getattr(self.bees_object.settings, "smiles_mode", "auto")
        if est_enabled and est_backend:
            self.logger.info(f"Kinetics estimation: {est_backend} (smiles_mode={smiles_mode})")
        elif est_enabled:
            self.logger.info("Kinetics estimation: enabled (no backend specified)")
        else:
            self.logger.debug("Kinetics estimation: disabled")
        # Removed placeholder solver log to avoid noise
        
        # Generate reactions using ModelGenerator
        # Resolve DB path from database.name (e.g. db -> db/db.csv)
        db_name = self.bees_object.database.name.strip()
        db_path_from_name = os.path.join(BEES_PATH, "db", f"{db_name}.csv")
        default_db = os.path.join(BEES_PATH, "db", "db.csv")
        if os.path.exists(db_path_from_name):
            db_path = db_path_from_name
        else:
            self.logger.info(
                f"No specific parameter estimator defined for database '{self.bees_object.database.name}'."
            )

        if hasattr(self.bees_object, "settings") and hasattr(self.bees_object.settings, "end_time") and self.bees_object.settings.end_time is not None:
            self.logger.info(f"The simulation will run until {self.bees_object.settings.end_time} time units.")
        else:
            self.logger.info("No simulation end time specified (reaction network generation only).")
        # Removed placeholder solver log to avoid noise
        
        # Generate reactions using ModelGenerator
        # Resolve DB path from database.name (e.g. db -> db/db.csv)
        db_name = self.bees_object.database.name.strip()
        db_path_from_name = os.path.join(BEES_PATH, "db", f"{db_name}.csv")
        default_db = os.path.join(BEES_PATH, "db", "db.csv")
        if os.path.exists(db_path_from_name):
            db_path = db_path_from_name
        else:
            db_path = default_db
        if getattr(self.bees_object, "seed_model", None):
            self.logger.info(f"Using seed model: {self.bees_object.seed_model}")
        model_generator = ModelGenerator(
            bees_object=self.bees_object,
            logger=self.logger,
            output_directory=self.output_directory
        )
        
        # Load kinetic database with ontology
        model_generator.load_kinetic_database(db_path, ontology=self.ontology)
        
        # Generate reactions
        reactions = model_generator.generate_reactions()
        
        # Export reactions summary
        summary_path = None
        if reactions:
            summary_path = model_generator.export_reactions_summary()
            self.logger.info(f"✓ Exported reaction summary to: {summary_path}")

        execution_time = common.time_lapse(self.t0)
        self.logger.info(
            f"BEES execution for project '{self.project}' completed in {execution_time}."
        )
        self.logger.log_footer(success=True)
        
        # Return results dictionary
        return {
            'success': True,
            'project': self.project,
            'n_reactions': len(reactions),
            'execution_time': execution_time,
            'summary_path': summary_path,
            'sbml_path': None,  # SBML export not yet implemented
            'message': f'Execution completed successfully. Generated {len(reactions)} reaction(s).'
        }