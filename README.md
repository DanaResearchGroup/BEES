# BEES
Biochemical Engine for Enzymatic kinetic modelS

## Overview

BEES is a tool for generating biochemical reaction networks from enzyme and substrate inputs. It automatically creates reaction networks by combining enzyme-substrate pairs, retrieves kinetic parameters from a database when available and take envriomantal factors into consideration.  

## Purpose

BEES will helps researchers and engineering to model biochemical systems by:
- Generating reaction networks from enzyme and substrate specifications
- Retrieving kinetic parameters from a database
- Creating structured reaction summaries with stoichiometry and rate laws

## Current Capabilities

This version of BEES provides two main capabilities:

1. **Network Generation**: Automatically generates biochemical reaction networks from enzyme-substrate pairs. Given a set of enzymes and substrates, BEES creates all possible reactions, determines products, and establishes stoichiometry based on enzyme EC numbers.

2. **Kinetic Parameter Retrieval**: Queries a kinetic database to find kinetic parameters (Km, kcat, Vmax, ΔG) for enzyme-substrate reactions. When parameters are found, they are included in the reaction model. When not found, the system falls back to template-based generation without kinetic parameters.

## Installation

### Prerequisites

- Python 3.12 or higher
- Conda (for environment management)

### Setup

1. Clone the repository:
```bash
git clone <repository-url>
cd BEES
```

2. Create the conda environment:
```bash
conda env create -f environment.yaml
```

3. Activate the environment:
```bash
conda activate bees_env
```

4. Verify installation:
```bash
python BEES.py --help
```

## Quick Start

1. Create an input YAML file. You can start with the minimal example:
```bash
cp examples/minimal/input.yml my_input.yml
```

2. Edit `my_input.yml` to specify your species and enzymes. At minimum, you need:
   - A project name
   - At least one species (substrate)
   - At least one enzyme with an EC number
   - Environment settings (temperature, pH)
   - Database settings
   - Settings (end_time, time_step)

3. Run BEES:
```bash
python BEES.py --input_file my_input.yml
```

4. Check the output. BEES will create a project directory in `projects/` containing:
   - `reactions_summary.txt` - Summary of all generated reactions
   - `<ProjectName>.log` - Execution log
   - `<ProjectName>_errors.log` - Error log (if any)

## Running Examples

BEES includes several example configurations that demonstrate different capabilities. you will find all the explantion in the example file in " Test_project_summeary.md"

## Database

BEES uses a CSV database located at `db/db.csv`. The database contains kinetic parameters for various enzyme-substrate reactions. See `docs/kinetic_database_schema.md` for database structure details.

## Documentation

- `examples/TEST_PROJECTS_SUMMARY.md` - Complete guide to all example projects
- `docs/kinetic_database_schema.md` - Database schema documentation
- `examples/example_output/` - Example output files showing expected format

## Development

To run tests:
```bash
pytest tests/
```

To run a specific test:
```bash
pytest tests/test_comprehensive_demo.py -v
```

## License

See LICENSE file for details.
