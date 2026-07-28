# BEES
Biochemical Engine for Enzymatic kinetic modelS

## What BEES does

- Builds biochemical reaction networks from species + enzymes (EC numbers).
- Pulls reaction templates and stoichiometry from a local database (37K+ reactions).
- Estimates missing kinetics via [CatPred](https://github.com/DanaResearchGroup/CatPred) when the database has no match.
- Computes thermodynamics (ΔG°′, Keq, reverse kcat) via [equilibrator-api](https://gitlab.com/equilibrator/equilibrator-api).
- Applies a small rules layer after kinetics and thermo: physics corrections always run, system-specific calibrations are optional.
- Grows the network with a core–edge iterative enlargement algorithm driven by ODE simulation.
- Exports a reaction summary, tables, plots, flux analysis, and an SBML Level 3 `model.xml`.

## Installation

**Prerequisites:** Python 3.12+, Conda. For CatPred: `git` and `wget` (or `curl`) on PATH.

1. **Clone and enter the repo**
   ```bash
   git clone https://github.com/DanaResearchGroup/BEES.git
   cd BEES
   ```

2. **Install BEES + CatPred**
   ```bash
   ./install.sh
   ```
   This creates the `bees_env` conda env, clones CatPred into the **parent** of BEES (`../CatPred`, `../catpred_pipeline`), downloads pretrained data (~1 GB), and writes `.env.bees` in the BEES root. BEES loads `.env.bees` automatically (it is gitignored).

   BEES only, no CatPred: `make install-bees`. No option to estimate kinetic in that method.
3. **Run**
   ```bash
   conda activate bees_env
   python BEES.py -i projects/minimal/input.yml
   ```

4. **If install cannot find `kcat/` and `km/`:** set `CATPRED_CHECKPOINT_BASE` to the directory that contains those folders (often `.../catpred_pipeline/data/pretrained/production`). Put it in `.env.bees` or export it before running. See [Kinetics estimation (CatPred)](#kinetics-estimation-catpred).

## Quick start

1. Copy an input YAML and edit it: `cp projects/minimal/input.yml my_input.yml`
2. Set a project name, at least one species and one enzyme (with EC number), `environment` (temperature, pH), `database.name`, and `settings` (at least `end_time` or another termination criterion).
3. Run: `python BEES.py -i my_input.yml`  
   Other flags: `-v 10` (log level 10/20/30/40/50), `-o <dir>` (output directory), `-p <name>` (project name).
4. Output lands in `<project_directory>/output/` (or `settings.output_directory`):
   - `reactions_summary.txt`, `output.log`
   - `simulation_profiles.csv`, plots
   - `flux_analysis.csv`, `core_reactions_species.csv`, `edge_reactions_species.csv`
   - `model.xml` (SBML Level 3)

Example configs in `projects/`:

- `minimal/` — first three glycolysis steps
- `fattyAcidSynthesis/fattyAcidSynthesis_ecoli/` — *E. coli* FAS II (uses `database.name: ecoli`)
- `commented/` — annotated YAML template (every allowed key)

## Input YAML

Unknown keys are rejected. Full key list: `projects/commented/input.yml` 

Required top-level blocks: `project`, `species`, `enzymes`, `environment`, `settings`, `database`.

Enzymes that need CatPred must include `amino_acid_sequence` (uppercase).

Two run modes (`settings.simulation_mode`):

- `iterative` — core–edge enlargement + ODE (default when a termination criterion is set and `toleranceMoveToCore > 0`)
- `batch` — reaction discovery only, no kinetics-driven pruning

At least one of `end_time`, `termination_conversion`, or `termination_rate_ratio` is required in iterative mode.

## Database

Set `database.name` to a file in `db/` (without `.csv`):

- `db` → `db/db.csv` (general BKMS-derived set)
- `ecoli` → `db/ecoli.csv`

```yaml
database:
  name: db
```

## Kinetics estimation (CatPred)

CatPred runs in a separate conda env (`catpred`). After `./install.sh`, just activate `bees_env` and run.

Enable in the input:

```yaml
settings:
  estimate_kinetics: true
  kinetics_estimator: catpred
  kinetics_include_sd: false
  smiles_mode: auto   # or 'interactive' to prompt for missing SMILES
```

Predictions are cached on disk by default (`~/.cache/bees/catpred_predictions.pkl`), so the first run of a network is slow and later runs skip the ML step. Disable with `BEES_CATPRED_CACHE=off`.

**Manual CatPred setup (without install.sh):** clone CatPred into a sibling of BEES (`CatPred` and `catpred_pipeline`), download the pretrained archive into `catpred_pipeline`, extract it, then create the `catpred` conda env. Set `CATPRED_DIR` to the CatPred clone, `CATPRED_CHECKPOINT_BASE` to the directory that contains `kcat/` and `km/`, and `CATPRED_CONDA_ENV=catpred`.

## Thermodynamics (equilibrator-api)

BEES computes ΔG°′ / Keq at the model's pH, ionic strength, pMg, and temperature, then gets reverse kcat from the Haldane relation. The first run downloads Component-Contribution data (~hundreds of MB) into `~/.cache/equilibrator`.

Set `BEES_DISABLE_THERMO=1` to skip this layer (every reaction is then treated as irreversible).

## Rules and calibrations

Physics corrections (irreversibility cutoffs, Haldane reverse kcat, and similar) always run.

System-specific calibrations are **off by default**. 

This version inculde only small proof of conecpt calibration and rules.

```yaml
settings:
  calibrations:
    - fabI_enoyl_reductase_measured_kcat
    - tesa_long_chain_preference
```


## Development

```bash
conda activate bees_env
pytest tests/ -v    # or: make test
```

## License

See LICENSE.
