# BEES
Biochemical Engine for Enzymatic kinetic modelS

## What BEES does

- Builds biochemical reaction networks from enzyme-substrate pairs (EC-based stoichiometry).
- Pulls kinetic parameters (Km, kcat, Vmax, delta-G) from a local database (37K+ reactions from BKMS).
- Can estimate missing kinetics via [CatPred](https://github.com/DanaResearchGroup/CatPred) when the database has no match.

## Installation

**Prerequisites:** Python 3.12+, Conda. For CatPred: `git` and `wget` (or `curl`) on PATH.

From the repo root (cd ~/BEES):

```bash
./install.sh
```

Or `make install`. This creates the `bees_env` conda environment, installs CatPred in the **parent directory** of BEES (sibling folders `CatPred` and `catpred_pipeline`, ARC-style), and writes `.env.bees` in the BEES root. For BEES only (no CatPred): `make install-bees` or create the env manually with `conda env create -f environment.yaml`.

- **BEES only:** `conda activate bees_env` then `python BEES.py -i projects/minimal/input.yml`
- **With CatPred:** `conda activate bees_env`, `source .env.bees`, then run BEES with an input that has `estimate_kinetics: true`

`.env.bees` is gitignored.

## Quick Start

1. Copy and edit an input file: `cp projects/minimal/input.yml my_input.yml` — set project name, at least one species and one enzyme (with EC number), environment (temperature, pH), `database.name: db`, and settings (end_time, time_step).
2. Run: `python BEES.py --input_file my_input.yml`
3. Output appears in `output/` in your project folder: `reactions_summary.txt`, logs.

## Examples

Example configs in `projects/`: `minimal/`, `Glycolysis/`, `fattyAcidSynthesis/FattyAcidSynthesisDemo/`, `ComprehensiveDemo/`, `commented/` (commented template). See `projects/Project_folder_README.md` for descriptions.

## Database

Kinetics are read from `db/db.csv`. In your input YAML set `database.name: db`.

## Kinetics estimation (CatPred)

CatPred runs in a separate conda env. If you used `./install.sh`, run `source .env.bees` after `conda activate bees_env` when using estimation.

Enable in your input:

```yaml
settings:
  estimate_kinetics: true
  kinetics_estimator: catpred
  kinetics_include_sd: false
  smiles_mode: auto   # or 'interactive' to prompt for missing SMILES
```

Enzymes that need estimation must have `amino_acid_sequence` (or BEES will try to resolve it by EC). When the database has no kinetics, BEES calls CatPred as a subprocess; predictions go into the reaction summary.

**Manual CatPred setup (without install.sh):** clone CatPred, create a `catpred` conda env, download production checkpoints for kcat/km (and optionally ki). Set `CATPRED_DIR`, `CATPRED_CHECKPOINT_BASE`, and `CATPRED_CONDA_ENV` (or use the `.env.bees` written by install.sh). Defaults are in `bees/kinetics_estimator.py`.

## Project structure

```
BEES/
  BEES.py       # CLI
  install.sh    # One-command setup
  bees/         # Core package (main, schema, model_generator, kinetics_estimator, ...)
  db/           # db.csv, reaction_database.py, ontology.yaml
  projects/     # Example configs
  tests/
```

## Development

Run tests: `pytest tests/ -v` (or `make test`).

## License

See LICENSE.
