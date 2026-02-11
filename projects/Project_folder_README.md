# BEES Project Examples Summary

This document lists the current example projects under `projects/`.

## How to Run Any Example

From the repo root:
```bash
python BEES.py -i projects/<PROJECT_NAME>/input.yml
```

## Current Example Projects

### 1. Minimal
**Purpose:** Smallest viable example for quick validation.

**Path:** `projects/minimal/input.yml`

**Highlights:**
- Single species and single enzyme
- Uses `db` database
- CatPred kinetics estimation enabled

---

### 2. Glycolysis
**Purpose:** Multi-enzyme pathway example with realistic inputs.

**Path:** `projects/Glycolysis/input.yml`

**Highlights:**
- Many enzymes with amino acid sequences
- CatPred estimation enabled
- Demonstrates larger reaction network generation

---

### 3. FattyAcidSynthesisDemo
**Purpose:** Larger pathway example with long enzyme sequences.

**Path:** `projects/fattyAcidSynthesis/FattyAcidSynthesisDemo/input.yml`

**Highlights:**
- Multiple enzymes with long amino acid sequences
- CatPred estimation enabled
- Shows complex substrate/enzyme combinations

---

## Output Files

Each run generates output files inside the project directory:

1. `<ProjectName>.log` - main execution log
2. `<ProjectName>_errors.log` - error log (if any)
3. `reactions_summary.txt` - summary of generated reactions
4. `input.yml` - validated input snapshot

## Database Information

BEES uses the `db/db.csv` database. In input files:

- Use `database.name: db`
- Kinetics estimation is configured under `settings`, not `database`
- `settings.kinetics_estimator` supports `catpred`

## Notes

- `amino_acid_sequence` is optional for enzymes, but CatPred requires it.
- If missing, BEES attempts to fetch the sequence from the database by EC number.

**Last Updated:** Feb 4, 2026
# BEES Test Projects Summary

*This document provides an overview of all test projects in the `/projects/` folder and how to run them*

***How to Run Any Test***
1. Open the terminal
2. Run the following command (from the project root):
```bash
python BEES.py -i projects/<TEST_NAME>/input.yml
```



***Output Files***

Each test generates the following files in an `output/` folder inside its project directory:

1. **`reactions_summary.txt`** - Human-readable summary of the network, including kinetic parameters and connectivity.
2. **`<ProjectName>.log`** - Detailed execution log.
3. **`<ProjectName>_errors.log`** - Log for tracking any warnings or errors.
4. **`input.yml`** - The validated version of the input file used for the run.

---


## Database Information

BEES uses a comprehensive kinetic database located at `db/db.csv` (tracked via Git LFS), currently containing **37,119 reactions** originally imported from the **BKMS** dataset (BRENDA, KEGG, MetaCyc, SABIO-RK) and passed the filtration process.

**Features:**
- **Full Stoichiometry:** Parsed from reaction templates (e.g., `ATP + (R)-pantoate + beta-alanine`).
- **SMILES Mapping:** Includes primary and candidate SMILES for reactants and products.
- **Trigger Logic:** Smart search ignores general cofactors (ATP, H2O, etc.) unless they are specifically matched in the database.

---

**Important Notes**

1. **Cofactor Handling:** BEES automatically adds cofactors required by the database stoichiometry. It also adds standard products for cofactors (e.g., ATP consumed → ADP + Pi produced, except in phosphorylation where Pi is part of the product).
2. **Placeholders:** Fields like `end_time`, `time_step`, and `kinetics_estimator` are currently placeholders for future ODE simulation and estimation modules.
3. **SMILES:** SMILES are required for all non-enzyme species to ensure chemical identity.

---

**Last Updated:** January 15, 2026
