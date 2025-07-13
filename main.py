"""
The Biochemical Engine for Enzymatic kinetic modelS (bees) for iterative kinetic model generation and refinement
# This is probably the most important module in the code.

# TODO:
# 1. Learn from T3 and other papers: Continue researching T3 or similar projects for architectural patterns
#    and specific functionalities that might be relevant for BEES's kinetic model generation and refinement.
# 2. Refine imports: Continuously evaluate and add/remove imports as BEES's functionality evolves.
# 3. Implement core logic for kinetic modeling: This will involve orchestrating different steps like
#    data processing, model building, simulation, and refinement.
# 4. Write comprehensive tests for the main module.
"""


import os  # Essential for path manipulations.
import yaml
import argparse 
from typing import Optional, Union, List, Any, Dict  # For type hinting, improving code clarity and maintainability.
import argparse
import datetime
import logging
import shutil
import sys
import time


from bees.logger import Logger 
from bees.schema import InputBase 
from bees.common import (
    PROJECTS_BASE_PATH,
    DATA_BASE_PATH, # Included for completeness if used later
    save_yaml_file,
    dict_to_str,
    time_lapse,
    get_git_branch,
    get_git_commit,
    VERSION, # Assuming VERSION is defined in common.py
    read_yaml_file, # For loading input YAML
    InputError # For custom exceptions
)


# Commented out imports - keep them commented for now, if their functionality is not yet implemented.
# from bees.runners.rmg_runner import rmg_runner
# from bees.simulate.factory import simulate_factory
# from bees.utils.libraries import add_to_rmg_libraries
# from bees.utils.writer import write_pdep_network_file, write_rmg_input_file


# RMG_THERMO_LIB_BASE_PATH and RMG_KINETICS_LIB_BASE_PATH:
# These depend on 'rmg_settings'. If you plan to integrate RMG,
# you'll need to define how 'rmg_settings' are loaded/accessed.
# RMG_THERMO_LIB_BASE_PATH = os.path.join(rmg_settings['database.directory'], 'thermo', 'libraries')
# RMG_KINETICS_LIB_BASE_PATH = os.path.join(rmg_settings['database.directory'], 'kinetics', 'libraries')


# The load_yaml function is fine as is.
def load_yaml(file_path):
    with open(file_path, 'r') as stream:
        return yaml.safe_load(stream)

class BEES(object):
    """
    The main class for the BEES platform.
    This class orchestrates the entire workflow from input parsing to execution
    of kinetic model generation and refinement.

    Attributes:
        t0 (float): The timestamp when the BEES instance was initialized (for timing).
        project (str): The name of the current BEES project.
        project_directory (str): The full path to the project's working directory.
        output_directory (str): The designated output directory for results.
        input_schema (InputBase): An instance of the InputBase schema, holding all validated input parameters.
        logger (Logger): The BEES logger instance for logging messages.
    """

    def __init__(self, input_data: Dict[str, Any]):
        """
        Initialize the BEES class with the given input data.

        This constructor validates the input against the schema, sets up the project
        directory, and initializes the logging system.

        Args:
            input_data (Dict[str, Any]): A dictionary containing all BEES input parameters.
                                         This dictionary is expected to conform to the
                                         `InputBase` schema.
        Raises:
            ValueError: If the input data is invalid according to the `InputBase` schema.
        """
        self.t0 = time.time() # Record start time for overall execution timing

        # 1. Validate input data against the schema using Pydantic
        try:
            self.input_schema = InputBase(**input_data)
        except Exception as e:
            # Catch Pydantic validation errors and re-raise as a more general ValueError
            raise ValueError(f"Invalid input parameters provided: {e}")

        self.project = self.input_schema.project
        
        # Determine project directory based on common.PROJECTS_BASE_PATH
        self.project_directory = os.path.join(PROJECTS_BASE_PATH, self.project)

        # Determine output directory: user-specified or default to project directory
        self.output_directory = self.input_schema.settings.output_directory \
                                if hasattr(self.input_schema.settings, 'output_directory') and self.input_schema.settings.output_directory \
                                else self.project_directory



        # 2. Set up project directory (and output directory if different)
        # Create project directory if it doesn't exist after validating input
        if not os.path.exists(self.project_directory):
            os.makedirs(self.project_directory)
            # Use print here as the logger might not be fully initialized yet
            print(f"INFO: Created project directory: {self.project_directory}")
        else:
            print(f"INFO: Project directory already exists: {self.project_directory}")

        # Create output directory if it's different from project directory
        if self.output_directory != self.project_directory and not os.path.exists(self.output_directory):
            os.makedirs(self.output_directory)
            print(f"INFO: Created output directory: {self.output_directory}")

        
        # 3. Initialize the Logger
        # Get verbose level from schema.settings.verbose (default to INFO if not specified)
        verbose_level = self.input_schema.settings.verbose \
                        if hasattr(self.input_schema.settings, 'verbose') and self.input_schema.settings.verbose is not None \
                        else logging.INFO
        
        self.logger = Logger(
            project=self.project,
            project_directory=self.project_directory,
            verbose=verbose_level,
            t0=self.t0 # Pass the initial timestamp for logging header/footer
        )

        # Log header and system information using the initialized logger
        self.logger.log_header()
        self.logger.info(f'BEES version: {VERSION}')
        git_branch = get_git_branch()
        git_commit, git_date = get_git_commit()
        if git_branch:
            self.logger.info(f'Git branch: {git_branch}')
        if git_commit:
            self.logger.info(f'Git commit: {git_commit} ({git_date})')

        self.logger.info(f'BEES project {self.project} initialized successfully in {time_lapse(self.t0)}.')
        # Log the validated input parameters, excluding default values for conciseness
        self.logger.log_args(self.input_schema.model_dump(exclude_defaults=True))


        # 4. Save the validated input to the project directory for reproducibility
        # This creates a canonical input.yml in the project folder
        self.logger.info(f"Saving validated input to {os.path.join(self.project_directory, 'input.yml')}")
        save_yaml_file(os.path.join(self.project_directory, 'input.yml'), self.input_schema.model_dump())

        
        
        # TODO: Initialize other core components/attributes based on self.input_schema that we might need later
        # For example, if you have modules for generating reaction networks, simulations, etc., initialize
        # self.species_objects = self._create_species_objects()
        # self.reaction_network = None

    
    #Temporary placeholder for the execute method
    # This method will be expanded in the next steps

    def execute(self):
        """
        The main execution method for the BEES project.
        This method orchestrates the various steps of kinetic model generation and refinement.
        
        #TODO: Implement the other functios and the necssary logic for the BEES workflow. 
        
        For now there no much to execute, but this will be the entry point for the core logic of BEES.
        
        """
        self.logger.info(f"Starting BEES execution for project '{self.project}'...")

        # --- Placeholder for future core logic ---
        # This is where you would call the main computational modules of BEES.
        # Example of accessing validated input parameters:
        self.logger.info(f"Input has {len(self.input_schema.species)} species and {len(self.input_schema.enzymes)} enzymes.")
        self.logger.info(f"Environment Temperature: {self.input_schema.environment.temperature} K")
        self.logger.info(f"Using '{self.input_schema.database.rate_law}' rate law from database '{self.input_schema.database.name}'.")
        self.logger.info(f"Simulation will run until {self.input_schema.settings.end_time} time units.")

        # Example: Call a hypothetical generator module
        # self.logger.info("Generating initial reaction network...")
        # self.reaction_network = self.generator.generate_network(self.input_schema)
        # self.logger.info(f"Generated network with {len(self.reaction_network.reactions)} reactions.")

        # Example: Call a hypothetical simulation module
        # self.logger.info("Running kinetic simulation...")
        # simulation_results = self.simulator.run_simulation(self.reaction_network, self.input_schema)
        # self.logger.info("Simulation completed.")

        # Example: Call a hypothetical refinement module
        # self.logger.info("Refining model based on simulation results...")
        # self.refiner.refine_model(self.reaction_network, simulation_results, self.input_schema)
        # self.logger.info("Model refinement completed.")

        self.logger.info(f"BEES execution for project '{self.project}' completed in {time_lapse(self.t0)}.")
        self.logger.log_footer() # Log footer at the end of successful execution


def parse_and_load_input() -> Dict[str, Any]:
    """
    Parses command-line (CLI) arguments and loads/merges input data from a YAML file.
    CLI arguments will override values from the YAML file.

    Returns:
        Dict[str, Any]: A dictionary containing the merged input parameters.
    Raises:
        FileNotFoundError: If the specified input file does not exist.
        ValueError: If there's an error parsing the YAML file.
    """
    
    parser = argparse.ArgumentParser(description='BEES: Biochemical Engine for Enzymatic modelS')
    parser.add_argument('--project', type=str,
                        help='Name of the BEES project. If provided, overrides project name in input_file.')
    parser.add_argument('--input_file', type=str,
                        help='Path to a BEES input YAML file. If not provided, project and verbose must be CLI arguments.')
    parser.add_argument('--verbose', type=int,
                        choices=[logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL],
                        help='Verbosity level (10=DEBUG, 20=INFO, 30=WARNING, 40=ERROR, 50=CRITICAL). If provided, overrides verbose in input_file.')
    parser.add_argument('--output_directory', type=str,
                        help='Path to the output directory. If provided, overrides output_directory in input_file.')

    args = parser.parse_args()

    input_data = {}
    if args.input_file:
        try:
            input_data = read_yaml_file(args.input_file)
        except FileNotFoundError:
            raise FileNotFoundError(f"Input file not found: {args.input_file}")
        except yaml.YAMLError as e:
            raise ValueError(f"Error parsing YAML input file {args.input_file}: {e}")

    # Override file-loaded parameters with any direct CLI arguments
    # Special handling for 'project', 'verbose' (nested under settings), and 'output_directory' (nested under settings)
    if args.project:
        input_data['project'] = args.project

    if args.verbose is not None:
        if 'settings' not in input_data:
            input_data['settings'] = {}
        input_data['settings']['verbose'] = args.verbose

    if args.output_directory:
        if 'settings' not in input_data:
            input_data['settings'] = {}
        input_data['settings']['output_directory'] = args.output_directory
    
    # Basic check if no input file and no project name provided via CLI
    if not args.input_file and not args.project:
        parser.error("Either '--input_file' must be provided, or '--project' must be provided via command line.")

    return input_data


if __name__ == '__main__':

    # Initialize a basic logger for early errors that occur before BEES logger is fully set up.
    # This ensures some logging even if the main BEES initialization fails.
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    bees_instance: Optional[BEES] = None # Initialize to None for error handling in case of early failure
    try:
        # Parse arguments and load/merge input data
        input_data_for_bees = parse_and_load_input()

        # Create an instance of the BEES application
        bees_instance = BEES(input_data=input_data_for_bees)

        # Execute the main BEES workflow
        bees_instance.execute()

    except Exception as e:
        # Centralized error handling for the main execution block
        if bees_instance and hasattr(bees_instance, 'logger'):
            # If the BEES logger was successfully initialized, use it for critical errors
            bees_instance.logger.critical(f"An unhandled error occurred during BEES execution: {e}", exc_info=True)
            bees_instance.logger.log_footer(success=False) # Log footer indicating failure
        else:
            # Fallback to basic logging if BEES logger failed to initialize
            logging.critical(f"An unhandled error occurred before BEES logger could be fully initialized: {e}", exc_info=True)
        sys.exit(1) # Exit with a non-zero status code to indicate an error

