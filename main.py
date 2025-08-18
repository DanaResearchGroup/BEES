"""
The Biochemical Engine for Enzymatic kinetic modelS (BEES) for iterative kinetic model generation and refinement

This is probably the most important module in the code.
 run the code by executing this script directly python main.py --input_file + directory

# TODO:
# 1. Learn from T3 and other papers: Continue researching T3 or similar projects for architectural patterns
#    and specific functionalities that might be relevant for BEES's kinetic model generation and refinement.
# 2. Refine imports: Continuously evaluate and add/remove imports as BEES's functionality evolves.
# 3. Implement core logic for kinetic modeling: This will involve orchestrating different steps like
#    data processing, model building, simulation, and refinement.
# 4. Write comprehensive tests for the main module.
"""


import sys
import os
import time
from typing import Any, Dict, List, Optional
import argparse

# Import necessary modules from BEES
import bees.common as common
from bees.logger import Logger 
from bees.schema import InputBase 

# Define global paths from common
BEES_PATH = common.BEES_PATH
PROJECTS_BASE_PATH = common.PROJECTS_BASE_PATH


#TODO: work on Parsing function 

def parse_and_load_input() -> dict:
    parser = argparse.ArgumentParser(description='BEES')
    parser.add_argument('-i','--input_file', type=str)
    parser.add_argument('-p','--project', type=str)
    parser.add_argument('-v','--verbose', type=int, choices=[10,20,30,40,50])
    parser.add_argument('-o','--output_directory', type=str)
    args = parser.parse_args()

    default_input_path = os.path.join(BEES_PATH, 'examples', 'minimal', 'input.yml')
    input_path = args.input_file or default_input_path

    # Single source of truth for YAML loading:
    input_data = common.read_yaml_file(input_path) 

    # Ensure nested dict exists before overrides
    input_data.setdefault('settings', {})

    if args.project: input_data['project'] = args.project
    if args.verbose is not None: input_data['settings']['verbose'] = args.verbose
    if args.output_directory: input_data['settings']['output_directory'] = args.output_directory

    if not input_data.get('project'):
        raise ValueError("Project name is required (YAML or --project).")
    return input_data


class BEES(object):
    """
    The main BEES application class.

    This class serves as the central orchestrator (manger) for a BEES project to run.
    It is responsible for:

    1.  Initializing the project environment (directories, logger).
    2.  Validating the raw user input against the defined schema
    3.  Converting validated input data into internal BEES objects (species, enzymes and more).
    4.  Executing the core kinetic modeling workflow.

    """

    def __init__(self, input_data: Dict[str, Any]):
        """
        Initializes the BEES application.

        Args:
            input_data (Dict[str, Any]): The raw input data dictionary loaded from the input file.
        
        Raises:
            ValueError: If input parameters are invalid according to the schema.
            FileNotFoundError: If project directory cannot be created.
        """
        self.t0 = time.time() # Start time
        self.input_data = input_data
        self.project: str = input_data.get("project", "default_project") # Default if missing
        
        # Determine project directory and output directory
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
                    raise ValueError(f"Absolute output path '{specified_output_dir}' is not within the project directory '{self.project_directory}'.")
                self.output_directory = specified_output_dir
            else:
                self.output_directory = os.path.join(self.project_directory, specified_output_dir)
        else:
            self.output_directory = self.project_directory

        # Create project directory if it does not exist
        try:
            os.makedirs(self.project_directory, exist_ok=True)
            if self.output_directory != self.project_directory:
                os.makedirs(self.output_directory, exist_ok=True)
        except OSError as e:
            raise FileNotFoundError(f"Failed to create project/output directory: {e}")

        # Initialize the logger - only after determining project_directory
        self.logger = Logger(
            project_directory=self.project_directory,
            verbose=input_data.get("settings", {}).get("verbose"),
            t0=self.t0
        )
        self.logger.info(f'BEES version: {common.VERSION}')
        git_branch = common.get_git_branch()
        commit_hash, commit_date = common.get_git_commit()
        self.logger.info(f'Git branch: {git_branch}')
        self.logger.info(f'Git commit: {commit_hash} ({commit_date})')
        self.logger.log_args(input_data)

        # Validate input using Pydantic schema
        try:
            self.input_schema = InputBase(**input_data)
            self.logger.info("Input validated successfully against the schema.")
            verified_input_path = os.path.join(self.project_directory, 'input.yml')
            common.save_yaml_file(verified_input_path, self.input_schema.model_dump(exclude_unset=True))
            self.logger.info(f"Saving validated input to {verified_input_path}")

        except Exception as e:
            self.logger.error(f"Input validation error: {e}")
            self.logger.log_footer(success=False)
            raise ValueError(f"Invalid input parameters provided: {e}")

        self.species_objects = self._create_species_objects(self.input_schema.species)
        self.enzyme_objects = self._create_enzyme_objects(self.input_schema.enzymes)

        self.logger.info(f"BEES project {self.project} initialized successfully in {common.time_lapse(self.t0)}.")

    def _create_species_objects(self, species_data: List[Any]) -> List[Any]:
        """
        Placeholder method to create species objects from input data.
        """
        self.logger.debug(f"Creating {len(species_data)} species objects...")
        return species_data

    def _create_enzyme_objects(self, enzyme_data: List[Any]) -> List[Any]:
        """
        Placeholder method to create enzyme objects from input data.
        """
        self.logger.debug(f"Creating {len(enzyme_data)} enzyme objects...")
        return enzyme_data
    
    def execute(self):
        """
        Executes the BEES simulation based on the loaded input schema.
        """
        self.logger.info(f"Starting BEES execution for project '{self.project}'...")

        self.logger.info(f"Input has {len(self.input_schema.species)} species and {len(self.input_schema.enzymes)} enzymes.")
        self.logger.info(f"Environment Temperature: {self.input_schema.environment.temperature} K")

        if hasattr(self.input_schema.database, 'rate_law'):
            self.logger.info(f"Using '{self.input_schema.database.rate_law}' rate law from database '{self.input_schema.database.name}'.")
        else:
            self.logger.info(f"Database name: '{self.input_schema.database.name}'. No specific rate law defined.")

        if hasattr(self.input_schema.database, 'parameter_estimator'):
            self.logger.info(f"Using parameter estimator: '{self.input_schema.database.parameter_estimator}'.")
        else:
            self.logger.info(f"No specific parameter estimator defined for database '{self.input_schema.database.name}'.")
            
        self.logger.info(f"The simulation will run until {self.input_schema.settings.end_time} time units.")
        self.logger.info(f"Using solver: {self.input_schema.database.solver}") 
        self.logger.debug("Simulating main BEES process (this will be replaced by actual logic)...")
        time.sleep(1) # Simulate some work

        self.logger.info(f"BEES execution for project '{self.project}' completed in {common.time_lapse(self.t0)}.")
        self.logger.log_footer(success=True)

def main():
    """
       Main entry point for the BEES application.

    This function is the first to be executed when the BEES script is run.
    
    It handles the higher -level orchestration of the application, including:
    1.  Loading the user's input configuration.
    2.  Instantiating the core BEES application class.
    3.  Executing the main BEES workflow.
    4.  Providing error handling and logging for critical failures.
    """
    
    bees_instance = None # Initialize to None for error handling in finally block
    try:
         # NEW: parse CLI + load/merge YAML
        input_data_for_bees = parse_and_load_input()

        # Initialize BEES with the loaded+merged input
        bees_instance = BEES(input_data=input_data_for_bees)

        # Execute workflow
        bees_instance.execute()

    except FileNotFoundError as e:
        if bees_instance is None or not Logger._initialized:
            print(f"CRITICAL: An unhandled error occurred before BEES logger could be fully initialized: {e}", file=sys.stderr)
        else:
            bees_instance.logger.critical(f"An unhandled error occurred: {e}")
            bees_instance.logger.log_footer(success=False)
        sys.exit(1)
    except Exception as e:
        if bees_instance is None or not Logger._initialized:
            print(f"CRITICAL: An unhandled error occurred before BEES logger could be fully initialized: {e}", file=sys.stderr)
        else:
            bees_instance.logger.critical(f"An unhandled error occurred during BEES execution: {e}")
            bees_instance.logger.log_footer(success=False)
        sys.exit(1)


if __name__ == '__main__':
    main()




