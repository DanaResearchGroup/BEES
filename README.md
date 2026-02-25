# BEES
Biochemical Engine for Enzymatic kinetic modelS

## What BEES does

- Builds biochemical reaction networks from enzyme-substrate pairs (EC-based stoichiometry).
- Pulls kinetic parameters (Km, kcat, Vmax, delta-G) from a local database (37K+ reactions from BKMS).
- Can estimate missing kinetics via [CatPred](https://github.com/DanaResearchGroup/CatPred) when the database has no match.

## Installation

**Prerequisites:** Python 3.12+, Conda. For CatPred: `git` and `wget` (or `curl`) on PATH. Use a normal directory (not e.g. inside Trash).

### Clone and install this branch (from zero)

1. **Clone and enter the repo** (use your branch if different from `dev-for-pr`):
   ```bash
   git clone <repo-url> BEES
   cd BEES
   git checkout dev-for-pr
   ```

2. **One-shot install (BEES + CatPred):**
   ```bash
   ./install.sh
   ```
   This creates the `bees_env` conda env, clones CatPred into the **parent** of BEES (`../CatPred`, `../catpred_pipeline`), downloads and extracts the pretrained data (~1 GB), and writes `.env.bees` in the BEES root. If the download or extraction fails, the script exits with a clear error.

3. **Run BEES:**
   ```bash
   conda activate bees_env
   python BEES.py -i projects/Glycolysis/input.yml
   ```
   No need to `source .env.bees`; BEES loads it automatically when present.

4. **If the install script reports that it could not find `kcat/` and `km/`:** set `CATPRED_CHECKPOINT_BASE` manually to the directory that contains those folders (often `.../catpred_pipeline/data/pretrained/production`). Put it in `.env.bees` or export it before running; see [Manual CatPred setup](#kinetics-estimation-catpred) below.

**BEES only (no CatPred):** run `make install-bees` (or `conda env create -f environment.yaml -n bees_env`) instead of `./install.sh`. Then run as above; kinetics come from the database only.

`.env.bees` is gitignored and is created by `install.sh` when CatPred is installed.

## Quick Start

1. Copy and edit an input file: `cp projects/minimal/input.yml my_input.yml` — set project name, at least one species and one enzyme (with EC number), environment (temperature, pH), `database.name: db`, and settings (end_time, time_step).
2. Run: `python BEES.py --input_file my_input.yml`
3. Output appears in `output/` in your project folder: `reactions_summary.txt`, logs.

## Examples

Example configs in `projects/`: `minimal/`, `Glycolysis/`, `fattyAcidSynthesis/FattyAcidSynthesisDemo/`, `ComprehensiveDemo/`, `commented/` (commented template). See `projects/Project_folder_README.md` for descriptions.

## Database

Kinetics are read from `db/db.csv`. In your input YAML set `database.name: db`.

## Kinetics estimation (CatPred)

CatPred runs in a separate conda env. If you used `./install.sh`, BEES auto-loads `.env.bees` when present, so just `conda activate bees_env` and run.

Enable in your input:

```yaml
settings:
  estimate_kinetics: true
  kinetics_estimator: catpred
  kinetics_include_sd: false
  smiles_mode: auto   # or 'interactive' to prompt for missing SMILES
```

Enzymes that need estimation must have `amino_acid_sequence` (or BEES will try to resolve it by EC). When the database has no kinetics, BEES calls CatPred as a subprocess; predictions go into the reaction summary.

**Manual CatPred setup (without install.sh):** clone CatPred into a sibling of BEES (e.g. `CatPred` and `catpred_pipeline`), download the pretrained archive into `catpred_pipeline`, extract it, then create the `catpred` conda env. Set `CATPRED_DIR` to the CatPred clone, `CATPRED_CHECKPOINT_BASE` to the directory that contains `kcat/` and `km/` (often `.../catpred_pipeline/data/pretrained/production`), and `CATPRED_CONDA_ENV=catpred`. BEES auto-loads `.env.bees` when present; otherwise export these in your shell before running.

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
