# BEES Project Examples Summary

This document lists the example projects under `projects/` and how to run them.

## How to Run Any Example

From the repo root:

```bash
python BEES.py -i projects/<PROJECT_PATH>/input.yml
```

Replace `<PROJECT_PATH>` with the folder path to the project (e.g. `minimal`, `Glycolysis`, `fattyAcidSynthesis/fattyAcid_ecoli`).

---


## Output Files

Each run generates output files in an `output/` folder inside the project directory (or in `settings.output_directory` if specified):

1. `output.log` - main execution log
2. `output_errors.log` - error log (if any)
3. `reactions_summary.txt` - summary of generated reactions
4. `input.yml` - validated input

For iterative mode (when `end_time` and `toleranceMoveToCore > 0`):

5. `ode_equations_iterN.txt` - ODE equations and parameters (including SD when available) per iteration
6. `simulation_profiles.csv` - concentration time-series (if `save_simulation_profiles: true`)
7. `flux_analysis.csv` 

---

## Database Information

BEES uses the database currently containing **37,119 reactions** imported from the **BKMS** dataset (BRENDA, KEGG, MetaCyc, SABIO-RK) and passed the filtration process.

**Features:**
- **Full Stoichiometry:** Parsed from reaction templates (e.g., `ATP + (R)-pantoate + beta-alanine`).
- **SMILES Mapping:** Includes primary and candidate SMILES for reactants and products.
- **Trigger Logic:** Smart search ignores general cofactors (ATP, H2O, etc.) unless they are specifically matched in the database.

**In input files:**
- Use `database.name: db`
- `settings.kinetics_estimator` supports `catpred` only

---

## Important Notes

1. **Cofactor Handling:** BEES  can automatically adds some genral cofactors required by the database stoichiometry. It also adds standard products for cofactors (e.g., ATP consumed -> ADP + product-Pi produced)

2. **Kinetics Estimation:** `kinetics_estimator: catpred` is implemented. `amino_acid_sequence` is optional for enzymes, but CatPred requires it for estimation. If missing, kinetics estimation will fail for that enzyme.

3. **Simulation Parameters:** `end_time` and `time_step` are used by the enlarger and simulator when iterative or ODE simulation is run. Have to be in the inuput.

4. **SMILES:** SMILES are required for all non-enzyme species to ensure chemical identity. the intial SMILES should privide from the user.

5. **Iterative enlargement stall fallback:** In some mechanisms you can see `R_char > 0` while `max_rr` stays near 0 for long periods, so no edge species get promoted and iterations run for a long time. BEES supports an optional reaction-level fallback (inspired by RMG’s “dynamics number” (dlnaccum)) to trigger promotion even when edge **net species rates** are uninformative.
   - Enable by setting `settings.toleranceMoveEdgeReactionToCore` to a finite value (it is **disabled by default** when unset).
   - Optional knobs:
     - `settings.dynamics_time_scale`: don’t evaluate dynamics before this simulation time (seconds). Default `0.0`.
     - `settings.dynamics_outer_step_interval`: evaluate dynamics every N stepwise outer-steps (performance knob). Default `1`.

Example:

```yaml
settings:
  toleranceMoveToCore: 1e-5
  toleranceInterruptSimulation: 1e-5
  toleranceMoveEdgeReactionToCore: 0.1
  dynamics_time_scale: 0.0
  dynamics_outer_step_interval: 1
```

---


