"""
The Biochemical Engine for Enzymatic kinetic modelS (BEES) for iterative kinetic model generation and refinement
# This is probably the most important module in the code.

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
import yaml # Direct import for YAML file handling
import logging

# Import internal modules
import bees.common as common
from bees.logger import Logger # Direct import of Logger class
from bees.schema import InputBase # Direct import of InputBase class

# Define global paths from common
BEES_PATH = common.BEES_PATH
PROJECTS_BASE_PATH = common.PROJECTS_BASE_PATH



# Load YAML utility function
# 

def load_yaml(file_path: str) -> Dict[str, Any]:
    """
    Loads content from a YAML file.
    Raises FileNotFoundError or yaml.YAMLError on failure.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    with open(file_path, 'r') as f:
        try:
            return yaml.safe_load(f)
        except yaml.YAMLError as e:
            raise yaml.YAMLError(f"YAML parsing error in file {file_path}: {e}")

def load_input_from_fixed_path() -> Dict[str, Any]:
    """
    Loads input data from a fixed. We will proably after it will be a command line argument with parsing.
   
    """
    fixed_input_path = os.path.join(BEES_PATH, 'examples', 'minimal', 'input.yml')
    
    if not os.path.exists(fixed_input_path):
        raise FileNotFoundError(f"Fixed input file not found at: {fixed_input_path}")
    try:
        input_data = load_yaml(fixed_input_path)
        return input_data
    except yaml.YAMLError as e:
        raise ValueError(f"Error parsing fixed YAML input file at {fixed_input_path}: {e}")
    

# Main BEES class
# This class will handle the initialization and execution of the BEES workflow.
class BEES(object):
    """
    The main BEES application class.
    Handles initialization, input validation and execution of the workflow.
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
        self.logger.log_header()
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

        self.species_objects = self.create_species_objects(self.input_schema.species)
        self.enzyme_objects = self.create_enzyme_objects(self.input_schema.enzymes)

        self.logger.info(f"BEES project {self.project} initialized successfully in {common.time_lapse(self.t0)}.")

    def create_species_objects(self, species_data: List[Any]) -> List[Any]:
        """
        Placeholder method to create species objects from input data.
        """
        self.logger.debug(f"Creating {len(species_data)} species objects...")
        return species_data

    def create_enzyme_objects(self, enzyme_data: List[Any]) -> List[Any]:
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
    Handles input loading, BEES initialization, and execution.
    """
    bees_instance = None # Initialize to None for error handling in finally block
    try:
        # Load input from a fixed path (temporarily bypassing command-line argument parsing)
        input_data_for_bees = load_input_from_fixed_path()
        
        # Initialize BEES with the loaded input data
        bees_instance = BEES(input_data=input_data_for_bees)
        
        # Execute the BEES workflow
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

# You can run the code by executing this script directly or with /opt/miniforge/envs/bees_env/bin/python /home/omerkfir/BEES/main.py



#TODO: work on Parsing function 

# def parse_and_load_input() -> Dict[str, Any]:
#     """
#     Parses arguments and loads/merges input data from a YAML file.
#     Arguments will override values from the YAML file.

#     Returns:
#         Dict[str, Any]: A dictionary containing the merged input parameters.
#     Raises:
#         FileNotFoundError: If the specified input file does not exist.
#         ValueError: If there's an error parsing the YAML file.
#     """
    
#     parser = argparse.ArgumentParser(description='BEES: Biochemical Engine for Enzymatic modelS')
#     parser.add_argument('--project', type=str,
#                         help='Name of the BEES project. If provided, overrides project name in input_file.')
#     parser.add_argument('--input_file', type=str,
#                         help='Path to a BEES input YAML file. If not provided, project and verbose must be CLI arguments.')
#     parser.add_argument('--verbose', type=int,
#                         choices=[logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL],
#                         help='Verbosity level (10=DEBUG, 20=INFO, 30=WARNING, 40=ERROR, 50=CRITICAL). If provided, overrides verbose in input_file.')
#     parser.add_argument('--output_directory', type=str,
#                         help='Path to the output directory. If provided, overrides output_directory in input_file.')

#     args = parser.parse_args()

#     input_data = {}
#     if args.input_file:
#         try:
#             input_data = read_yaml_file(args.input_file)
#         except FileNotFoundError:
#             raise FileNotFoundError(f"Input file not found: {args.input_file}")
#         except yaml.YAMLError as e:
#             raise ValueError(f"Error parsing YAML input file {args.input_file}: {e}")

#     # Override file-loaded parameters with any direct CLI arguments
#     # Special handling for 'project', 'verbose' (nested under settings), and 'output_directory' (nested under settings)
#     if args.project:
#         input_data['project'] = args.project

#     if args.verbose is not None:
#         if 'settings' not in input_data:
#             input_data['settings'] = {}
#         input_data['settings']['verbose'] = args.verbose

#     if args.output_directory:
#         if 'settings' not in input_data:
#             input_data['settings'] = {}
#         input_data['settings']['output_directory'] = args.output_directory
    
#     # Basic check if no input file and no project name provided via CLI
#     if not args.input_file and not args.project:
#         parser.error("Either '--input_file' must be provided, or '--project' must be provided")

#     return input_data

