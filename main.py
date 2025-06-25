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

import datetime

import os  # Essential for path manipulations.
import shutil # Useful for file operations like copying/deleting directories.
import yaml
import argparse # Crucial for command-line argument parsing if this is the main entry point.
from typing import List, Optional, Tuple, Union
import inspect  # Potentially useful for introspection, keep for now.
import re  # Potentially useful for string parsing, keep for now.


#from bees.common import (DATA_BASE_PATH, PROJECTS_BASE_PATH,VALID_CHARS, delete_root_rmg_log, # Keep if RMG integration is planned.
                         get_species_by_label, # Keep if species management is relevant.
                         time_lapse,
                         save_yaml_file # Add this import as it's used in write_bees_input_file
                         )
from bees.logger import Logger # Importing the Logger class
from bees.schema import InputBase # Importing the Pydantic schema for input validation


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


class BEES(object): # Changed class name to BEES (capitalized) for consistency and convention.
    """
    The main class for the BEES platform.

    Attributes:
        project (str): The name of the current BEES project.
        input_schema (InputBase): An instance of the InputBase schema, holding all input parameters.
        logger (Logger): The BEES logger instance for logging messages.
        project_directory (str): The full path to the project's working directory.
        t0 (float): The timestamp when the BEES instance was initialized.
        # Add other relevant attributes here as the project develops.
    """

    def __init__(self, 
                 input_file_path: Optional[str] = None, # Path to a YAML input file
                 **kwargs # Allows for direct passing of input parameters
                 ):
        """
        Initialize the BEES class.

        This constructor can accept parameters either from a YAML input file
        or directly as keyword arguments. Keyword arguments will override
        parameters loaded from the YAML file.

        Args:
            input_file_path (Optional[str]): Path to a YAML file containing BEES input parameters.
                                              Defaults to None.
            **kwargs: Direct input parameters that conform to the InputBase schema.
                      These will override values from the input file.
        """
        self.t0 = time.time()  # Initialize time at the very beginning of the constructor.

        input_data = {}
        if input_file_path:
            # Load parameters from the specified YAML file
            try:
                input_data = load_yaml(input_file_path)
            except FileNotFoundError:
                raise FileNotFoundError(f"Input file not found: {input_file_path}")
            except yaml.YAMLError as e:
                raise ValueError(f"Error parsing YAML input file {input_file_path}: {e}")

        # Override file-loaded parameters with any direct kwargs
        input_data.update(kwargs)

        # Validate input data against the schema
        try:
            self.input_schema = InputBase(**input_data)
        except Exception as e:
            # Pydantic validation errors are typically more specific, but a general catch is fine here.
            raise ValueError(f"Invalid input parameters provided: {e}")

        self.project = self.input_schema.project # Project name from the validated input.

        # Initialize the logger using the BEES Logger class.
        # The Logger class should handle setting up handlers (console/file) and verbosity.
        self.logger = Logger(
            project_name=self.project,
            project_directory=os.path.join(PROJECTS_BASE_PATH, self.project), # Pass relevant paths
            verbose=self.input_schema.verbose # Get verbose level from schema
        )
        self.logger.log_args(self.input_schema.dict(exclude_defaults=True)) # Log provided arguments

        # Set up project directory
        self.project_directory = os.path.join(PROJECTS_BASE_PATH, self.project)
        if not os.path.exists(self.project_directory):
            os.makedirs(self.project_directory)
            self.logger.info(f"Created project directory: {self.project_directory}")
        else:
            self.logger.info(f"Project directory already exists: {self.project_directory}")

        # Now that input is validated and project directory is set up,
        # you can initialize other attributes based on self.input_schema
        # Example:
        # self.species = [] # If species are managed by BEES
        # self.reactions = [] # If reactions are managed by BEES


    def as_dict(self) -> dict:
        """
        Convert the BEES object's current state to a dictionary representation,
        primarily by serializing its input schema.

        Returns:
            dict: The dictionary representation of the BEES object.
        """
        # The input_schema itself can be converted to a dictionary.
        # You might choose to include other runtime attributes if they are
        # important for a "restart" dictionary.
        return {
            'project': self.project,
            'input_schema': self.input_schema.dict(), # Convert the Pydantic schema to a dict
            # Add other runtime attributes here if needed for state saving
        }
    

    def write_bees_input_file(self, path: Optional[str] = None, all_args: bool = False) -> None:
        """
        Write the BEES input file based on the provided input object.

        Args:
            path (Optional[str]): The path to save the input file. Defaults to
                                  'bees_auto_saved_input.yml' within the project directory.
            all_args (bool): Whether to include all arguments (including defaults)
                             or only those explicitly set by the user. Defaults to False.
        """
        if path is None:
            path = os.path.join(self.project_directory, 'bees_auto_saved_input.yml')
        
        # Ensure the path is a file path, not a directory.
        if os.path.isdir(path):
            path = os.path.join(path, 'bees_auto_saved_input.yml')
        
        base_path = os.path.dirname(path)
        if not os.path.isdir(base_path):
            os.makedirs(base_path)
            self.logger.debug(f"Created directory for input file: {base_path}") # Use debug for directory creation

        self.logger.info(f'Writing input file to {path}')
        
        # Decide which dictionary representation of the schema to save
        content_to_save = self.input_schema.dict() if all_args else self.input_schema.dict(exclude_unset=True)
        
        # Ensure save_yaml_file is imported from bees.common
        save_yaml_file(path=path, content=content_to_save)


# Example of how BEES might be instantiated (for testing/demonstration)
if __name__ == '__main__':
    # This block will run when main.py is executed directly
    # You can add argument parsing here using argparse if you want to run from command line
    parser = argparse.ArgumentParser(description="Run BEES project.")
    parser.add_argument("--project", type=str, required=True, help="Name of the BEES project.")
    parser.add_argument("--verbose", type=int, default=logging.INFO,
                        help="Logging level (10=debug, 20=info, 30=warning).")
    parser.add_argument("--input_file", type=str, help="Path to a BEES input YAML file.")

    args = parser.parse_args()

    try:
        bees_instance = BEES(
            input_file_path=args.input_file,
            project=args.project,
            verbose=args.verbose
        )
        bees_instance.logger.info(f"BEES project '{bees_instance.project}' initialized successfully.")
        
        # Example: Save the initialized input for demonstration
        bees_instance.write_bees_input_file(all_args=True)

        # Here you would call the main execution method, e.g., bees_instance.execute()
        # but that method doesn't exist yet.
        bees_instance.logger.info("BEES setup complete. Ready for execution.")

    except Exception as e:
        # A basic logger might not be fully set up yet for all errors at this stage
        print(f"Error during BEES initialization: {e}", file=sys.stderr)
        # If logger is already setup:
        # logging.getLogger('BEES').critical(f"Error during BEES initialization: {e}")
        exit(1)