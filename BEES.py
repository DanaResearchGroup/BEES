#!/usr/bin/env python3
r"""
Executable wrapper for BEES (Biochemical Engine for Enzymatic kinetic modelS)


This simply delegates to main.py so users don't need to call `python main.py`. 

For linux users:
Make sure this file is executable: chmod +x BEES.py

To run BEES from anywhere, add to your PATH (in ~/.bashrc or ~/.zshrc):
    export PATH="$PATH:$HOME/BEES"

Then you can run BEES like this:
cd ~/BEES
for linux users:
./BEES.py -i ~/BEES/examples/minimal/input.yml -p MyProject
for windows users:
.\BEES.bat -i examples\\minimal\\input.yml -p MyProject  

#note that some tests are still run good only on linux (for now)

"""

import sys
import os
import argparse

# Check if BEES modules can be imported (more flexible than checking environment name)
try:
    import bees.common as common
    import bees.schema
    import bees.main
except ImportError as e:
    print(f"BEES modules not found. Please ensure BEES is properly installed in your Python environment.")
    print(f"Error: {e}")
    sys.exit(1)

# Import necessary modules from BEES
from bees.main import BEES

BEES_PATH = common.BEES_PATH

if "bees_env" not in sys.executable:
    print("Please activate the 'bees_env' environment before running BEES.")
    sys.exit(1)

def parse_and_load_input() -> dict:
    parser = argparse.ArgumentParser(description="BEES")
    parser.add_argument("-i", "--input_file", type=str, help="Path to input YAML")
    parser.add_argument("-p", "--project", type=str, help="Project name")
    parser.add_argument("-v", "--verbose", type=int, choices=[10, 20, 30, 40, 50], help="Log level")
    parser.add_argument("-o", "--output_directory", type=str, help="Output directory")
    args = parser.parse_args()

    default_input_path = os.path.join(BEES_PATH, "examples", "minimal", "input.yml")
    input_path = args.input_file or default_input_path

    input_data = common.read_yaml_file(input_path)
    input_data.setdefault("settings", {})

    if args.project:
        input_data["project"] = args.project
    if args.verbose is not None:
        input_data["settings"]["verbose"] = args.verbose
    if args.output_directory:
        input_data["settings"]["output_directory"] = args.output_directory

    if not input_data.get("project"):
        raise ValueError("Project name is required!")

    return input_data

def main():
    """
    the main BEES excutable function
    Entrypoint without try/except — exceptions bubble up.
    """
    input_data = parse_and_load_input()
    bees_instance = BEES(input_data=input_data)
    bees_instance.execute()


if __name__ == "__main__":
    main()
