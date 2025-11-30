"""
The Biochemical Engine for Enzymatic kinetic modelS (BEES) for iterative kinetic model generation and refinement

This is probably the most important module in the code.
 run the code by executing this script directly python bees.py --input_file + directory
(decribed in BEES.py as well)
"""


import os
import time
from typing import Any, Dict, List

import bees.common as common
from bees.logger import Logger
from bees.schema import InputBase
from bees.model_generator import ModelGenerator

# Base paths
BEES_PATH = common.BEES_PATH
PROJECTS_BASE_PATH = common.PROJECTS_BASE_PATH


class BEES(object):
    """
    The main BEES application class.
    Orchestrates: setup, schema validation, logging, execution.
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

    def execute(self) -> Dict[str, Any]:
        """
        Execute BEES workflow: generate biochemical reaction model.
        
        This method orchestrates:
        1. Logging input configuration
        2. Loading kinetic database
        3. Generating reactions from enzyme-substrate pairs
        4. Exporting results
        
        Returns:
            dict: Execution results containing:
                - success (bool): Whether execution completed successfully
                - project (str): Project name
                - n_reactions (int): Number of reactions generated
                - reactions (List[GeneratedReaction]): Generated reactions
                - summary_path (str): Path to reactions summary file
                - execution_time (str): Total execution time
        """

        self.logger.info(f"Starting BEES execution for project '{self.project}'...")

        self.logger.info(
            f"Input has {len(self.bees_object.species)} species and {len(self.bees_object.enzymes)} enzymes."
        )
        self.logger.info(f"Environment Temperature: {self.bees_object.environment.temperature} K")

        if hasattr(self.bees_object.database, "solver"):
            self.logger.info(
                f"Using solver: '{self.bees_object.database.solver}' "
                f"with database '{self.bees_object.database.name}'."
            )
        else:
            self.logger.info(
                f"Database name: '{self.bees_object.database.name}'. Using default solver."
            )

        if hasattr(self.bees_object.database, "kinetics_estimator"):
            self.logger.info(
                f"Using kinetics estimator: '{self.bees_object.database.kinetics_estimator}'."
            )
        else:
            self.logger.info(
                f"No specific kinetics estimator defined for database '{self.bees_object.database.name}'."
            )

        if hasattr(self.bees_object, "settings") and hasattr(self.bees_object.settings, "end_time"):
            self.logger.info(f"The simulation will run until {self.bees_object.settings.end_time} time units.")
        else:
            self.logger.warning("End time setting is not defined.")
        self.logger.info(f"Using solver: {self.bees_object.database.solver}")
        
        """
        Model Generation Phase
        """
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("STARTING MODEL GENERATION")
        self.logger.info("=" * 60)
        
        # Initialize model generator
        model_gen = ModelGenerator(
            bees_object=self.bees_object,
            logger=self.logger,
            output_directory=self.output_directory
        )
        
        # Load kinetic database
        db_path = os.path.join(common.BEES_PATH, 'db', 'db.csv')
        try:
            model_gen.load_kinetic_database(db_path)
        except FileNotFoundError:
            self.logger.warning(f"Kinetic database not found at {db_path}. "
                              f"Will proceed without database parameters.")
        except Exception as e:
            self.logger.error(f"Error loading kinetic database: {e}")
            self.logger.log_footer(success=False)
            raise
        
        # Generate reactions
        try:
            reactions = model_gen.generate_reactions()
            
            if not reactions:
                self.logger.warning("No reactions were generated! Check your input configuration.")
                # Store empty results as instance variables
                self.reactions = []
                self.model_generator = model_gen
                
                results = {
                    'success': False,
                    'project': self.project,
                    'n_reactions': 0,
                    'reactions': [],
                    'summary_path': None,
                    'execution_time': common.time_lapse(self.t0),
                    'message': 'No reactions generated. Check enzyme-substrate pairs and database.'
                }
            else:
                self.logger.info("")
                self.logger.info(f"✓ Successfully generated {len(reactions)} reaction(s)")
                
                # Export summary
                summary_path = model_gen.export_reactions_summary()
                self.logger.info(f"✓ Exported reaction summary to: {summary_path}")
                
                # Store results as instance variables (Option 3)
                self.reactions = reactions
                self.model_generator = model_gen
                
                # Prepare return dictionary (Option 1)
                results = {
                    'success': True,
                    'project': self.project,
                    'n_reactions': len(reactions),
                    'reactions': reactions,
                    'summary_path': summary_path,
                    'execution_time': common.time_lapse(self.t0),
                    'database_path': db_path
                }
                
        except Exception as e:
            self.logger.error(f"Error during reaction generation: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            self.logger.log_footer(success=False)
            raise

        self.logger.info("")
        self.logger.info(
            f"BEES execution for project '{self.project}' completed in {results['execution_time']}."
        )
        self.logger.log_footer(success=results['success'])
        
        return results