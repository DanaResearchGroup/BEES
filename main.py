"""
The Biochemical Engine for Enzymatic kinetic modelS (bees) for iterative kinetic model generation and refinement
# this probably the most important module in the code

#TODO: 1. learn from T3 what they do and what modifications need to be done
        2. check from other pappers what they do and what modifications need to be done
   
"""

""""
spesific tasks for importing:

#TODO: 1. check if this import is needed or not 
#      2. see which import is needed to be add 
#      3. writing the main tests
"""


import datetime
import inspect
import os
import re
import shutil
import yaml
import argparse
from typing import List, Optional, Tuple, Union



from bees.common import (DATA_BASE_PATH,
                       PROJECTS_BASE_PATH,
                       VALID_CHARS,
                       delete_root_rmg_log,
                       get_species_by_label,
                       time_lapse)
from bees.logger import Logger
from bees.schema import InputBase




# from bees.runners.rmg_runner import rmg_runner
# from bees.simulate.factory import simulate_factory
# from bees.utils.libraries import add_to_rmg_libraries
# from bees.utils.writer import write_pdep_network_file, write_rmg_input_file

"""""
#TODO: 1. check if OS is needed or not
#      2. if so, learn how to use it 
"""


#RMG_THERMO_LIB_BASE_PATH = os.path.join(rmg_settings['database.directory'], 'thermo', 'libraries')
#RMG_KINETICS_LIB_BASE_PATH = os.path.join(rmg_settings['database.directory'], 'kinetics', 'libraries')


def load_yaml(file_path):
    with open(file_path, 'r') as stream:
        return yaml.safe_load(stream)


class bees(object):
    """
    The main class for the bees platform""
    
    #TODO:  1.learn fron T3 what they do and what modifications need to be done
            2. fill init objects with the correct attributes.
    """


    def __init__(self, 
                 project: str,
                 #TODO: fill this with the correct type):
        ):
        
        self.schema = InputBase.schema()
        self.logger = Logger()
        
        
        
        """
        Initialize the bees class with the given input.
        
        Args:
            input (InputBase): The input object containing the necessary parameters.
        """


    def as_dict(self) -> dict:
        """
        Convert the bees object to a dictionary representation.
        
        Returns:
            dict: The dictionary representation of the bees object.
        """
        return {
            'project': self.project,
          #TODO: fill this with the correct type
        }
    


    def write_bees_input_file(self, input: InputBase, path: Optional[str] = None, all_args: bool = True) -> None:
        """
        Write the bees input file based on the provided input object.

        Args:
            input (InputBase): The input object containing the necessary parameters.
            path (Optional[str]): The path to save the input file. Defaults to None.
            all_args (bool): Whether to include all arguments or only those set. Defaults to True.
        """
        if path is None:
            path = os.path.join(self.project_directory, 'bees_auto_saved_input.yml')
        if os.path.isdir(path):
            path += '/' if path[-1] != '/' else ''
            path += 'bees_auto_saved_input.yml'
        base_path = os.path.dirname(path)
        if not os.path.isdir(base_path):
            os.makedirs(base_path)
        self.logger.info(f'\n\nWriting input file to {path}')
        save_yaml_file(path=path, content=self.schema if all_args else self.schema_exclude_unset)
    

    def execute(self):
        """
        Execute the bees platform with the provided input object.
        """



#This File isT