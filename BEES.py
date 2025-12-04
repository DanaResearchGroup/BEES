#!/usr/bin/env python3
r"""
Executable wrapper for BEES (Biochemical Engine for Enzymatic kinetic modelS)


This simply delegates to main.py so users don’t need to call `python main.py`. 

For linux users:
Make sure this file is executable: chmod +x BEES.py

To run BEES from anywhere, add to your PATH (in ~/.bashrc or ~/.zshrc):
    export PATH="$PATH:$HOME/BEES"

Then you can run BEES like this:
cd ~/BEES
for linux users:
./BEES.py -i ~/BEES/examples/minimal/input.yml 
for windows users:
.\BEES.bat -i examples\\minimal\\input.yml 

# Note that some tests are still run good only on linux (for now)
# More examples can bo fonded on projects folder.

"""

import sys
import os
import argparse

# Import necessary modules from BEES
import bees.common as common

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
    results = bees_instance.execute()
    
    # Print summary of results
    if results and results.get('success'):
        print(f"\n{'='*60}")
        print(f"✓ BEES Execution Summary:")
        print(f"  Project: {results['project']}")
        print(f"  Reactions Generated: {results['n_reactions']}")
        print(f"  Execution Time: {results['execution_time']}")
        if results.get('summary_path'):
            print(f"  Summary: {results['summary_path']}")
        print(f"{'='*60}\n")
    else:
        print(f"\n⚠ BEES execution completed with warnings or errors.")
        if results and results.get('message'):
            print(f"  {results['message']}\n")
    
    return results


if __name__ == "__main__":
    main()
