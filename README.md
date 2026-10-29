# BEES
Biochemical Engine for Enzymatic kinetic modelS

## What BEES does

- Builds biochemical reaction networks from species + enzymes using local db.
- Estimates missing kinetics(Km and Kcat) via integrated Deep learning estimator CatPred.
- Computes thermodynamics (ΔG°′, Keq,) via equilibrator-api package.
- Applies a small proof of conecpt rules layer after kinetics and thermo: physics corrections always run, system-specific calibrations are optional.
- Implements the flux based iterative algorithm
- Exports results a reaction summary, tables, plots, flux analysis and an SBML Level 3 (`model.xml`).

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

   BEES only, no CatPred: `make install-bees`.
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

- `minimal/` — first three glycolysis steps using the general db.
- `fattyAcidSynthesis/fattyAcidSynthesis_ecoli/` — *E. coli* FAS II (uses `database.name: ecoli`)
- `methotrexate_dhfr/`, `glycolysis_2dg/`, `fattyAcidSynthesis/fas_inhibitor_delivery/` — drug-delivery showcases (see below); `database.name` `folate` / `glycolysis` / `ecoli`
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

Enable in the input

```yaml
database:
  name: db
```

## Kinetics estimation (CatPred)

Enable in the input:

```yaml
settings:
  estimate_kinetics: true
  kinetics_estimator: catpred
  kinetics_include_sd: false
  smiles_mode: auto   # or 'interactive' to prompt for missing SMILES
```

Predictions able to be cached on disk by default (`~/.cache/bees/catpred_predictions.pkl`), so the first run of a network is slow and later runs skip the ML step. Disable with `BEES_CATPRED_CACHE=off`.

**Manual CatPred setup (without install.sh):** clone CatPred into a sibling of BEES (`CatPred` and `catpred_pipeline`), download the pretrained archive into `catpred_pipeline`, extract it, then create the `catpred` conda env. Set `CATPRED_DIR` to the CatPred clone, `CATPRED_CHECKPOINT_BASE` to the directory that contains `kcat/` and `km/`, and `CATPRED_CONDA_ENV=catpred`.

## Thermodynamics (equilibrator-api)

BEES computes ΔG°′ / Keq at the model's pH, ionic strength, pMg, and temperature. The first run downloads Component-Contribution data (~hundreds of MB) into `~/.cache/equilibrator`.

Set `BEES_DISABLE_THERMO=1` to skip this layer (every reaction is then treated as irreversible).

## Rules and calibrations

Enable in the input:

```yaml
settings:
  calibrations:
    - fabI_enoyl_reductase_measured_kcat
    - tesa_long_chain_preference
```


## Inhibitors and time-dependent dosing

BEES computes what a drug does to an enzyme pathway; the release profile of its carrier decides when it does it.

```yaml
species:
  - {label: methotrexate, concentration: 0, smiles: "...", reactive: false}   # enters only by dosing
enzymes:
  - label: DHFR
    measured_kinetics:            # optional experimental values, override the database / CatPred
      kcat: 12.0
      km: {FH2: 2.0e-4}           # mM per substrate
      default_km: 0.005           # Km of every substrate not listed (a substrate without Km counts as saturated)
    inhibition:
      - inhibitor: methotrexate
        type: tight_binding       # competitive | uncompetitive | noncompetitive | mixed (alpha:) | tight_binding
        ki: 1.0e-7                # mM; omit to estimate with CatPred (needs smiles + amino_acid_sequence)
        substrate: FH2            # competitive / uncompetitive / mixed only; default: first substrate with a Km
dosing:
  - species: methotrexate
    profile: biexponential        # bolus | constant | first_order | biexponential | weibull | korsmeyer_peppas | table
    dose: 0.002                   # mM released in total
    params: {fraction_fast: 0.1, half_life_fast: 7200, half_life_slow: 162000}   # seconds
    t_start: 0
    elimination_half_life: 10000  # or elimination_k (1/s): first-order loss of the free drug
```

Inhibition changes the rate law: competitive scales the substrate Km by `1 + I/Ki`, noncompetitive scales Vmax by `1/(1 + I/Ki)`, mixed and uncompetitive both (`bees/inhibition.py`); tight binding uses the Morrison quadratic with the enzyme concentration. The older `feedback_inhibition` is unchanged. Dosing adds a source `u(t)` (and the elimination) to the drug's equation, integrated in segments at the profile's breakpoints; both go into the ODE text and the SBML export (`DOSE_`/`ELIM_` pseudo-reactions, Ki/alpha parameters, inhibitors as modifiers).

Tools in `devtools/`:

- `release_sweep.py` — same total dose as free bolus / liposome / slow depot: output lost, time below half the control flux, inhibition AUC, peak free drug (CSV + PNG)
- `target_scan.py` — scale kcat of one enzyme at a time, output lost at 1/3/6 h
- `fetch_inhibitors.py` — median Ki per target/compound from ChEMBL (needs network access to ebi.ac.uk)

```bash
python devtools/release_sweep.py -i projects/methotrexate_dhfr/input.yml --drug methotrexate --output-species dTMP
python devtools/target_scan.py   -i projects/methotrexate_dhfr/input.yml --output-species dTMP
```

The kinetics in the three showcase projects are order-of-magnitude placeholders; see their READMEs.

## Development

```bash
conda activate bees_env
pytest tests/ -v    # or: make test
```

## License

See LICENSE.
