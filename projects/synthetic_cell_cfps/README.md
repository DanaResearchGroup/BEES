# Synthetic cell (cell-free protein synthesis): energy and protein output

**The recipe is real, the enzyme side is a guess.** The inner-solution concentrations come from Table S.2 of the
final-year project report (final requested inner solution volume for the synthetic cell, lab table of 23/03/2024):
HEPES-KOH pH 8 55 mM, Mg-acetate 14, K-acetate 50, NH4-acetate 155, PEG 6000 3 %, 3-PGA 40, amino acids (mixtures I
and II) 2.5 + 2.5, ATP 1.2, GTP 1, UTP 0.8, IPTG 1, sucrose 200 mM, S30-T7 lysate 34 %, plasmid 10 ug/mL (the NH4-acetate
and amino-acid stock concentrations were hard to read in the photographed table). The report gives no enzyme content of
the lysate, so **every enzyme concentration, kcat and Km here is a placeholder**. Nothing below is a result about
the real system.

## What is modelled (`db/cfps.csv`)

3-PGA -> 2-PGA (PGM) -> PEP (ENO) -> pyruvate + ATP from ADP (PYK), and one lumped translation step:
`amino acids + 4 ATP -> protein residue + 4 ADP + 4 Pi` (4 ATP-equivalents per residue, GTP not separate). The output
"protein residue" is in mM of residues. Buffer, salts, PEG, sucrose, IPTG and DNA are not modelled, and there is no ATP
regeneration from AMP (no adenylate kinase).

The recipe has 40 mM 3-PGA at t = 0. In `input.yml` 10 mM is inside at t = 0 and 30 mM is supplied through `dosing`
(first-order, half-life 1 h), so that the sweep can compare supplying all of it at once (the recipe) with slower supplies.

## Run

```bash
python BEES.py -i projects/synthetic_cell_cfps/input.yml
python devtools/release_sweep.py -i projects/synthetic_cell_cfps/input.yml \
  --drug 3-phospho-D-glycerate --output-species "protein residue" --dose 30 --time-scale 0.5
```

(In the sweep the "drug" is the supplied 3-PGA; "output lost" is negative because the control has only the 10 mM inside,
and "peak free drug" is the 3-PGA peak. Read the numbers as gains.)

## Result (placeholders), 8 h, 30 mM extra 3-PGA

| arm | ATP at 1 / 2 / 4 / 8 h (mM) | protein residue at 1 / 2 / 4 / 8 h (mM) |
|---|---|---|
| 10 mM only (no supply) | 0.27 / 0.11 / 0.07 / 0.05 | 2.7 / 2.8 / 2.8 / 2.8 |
| all 30 mM at t = 0 (the recipe) | 1.3 / 1.3 / 1.3 / 1.3 | 3.3 / 5.0 / 5.0 / 5.0 |
| liposome-like slow release | 0.65 / 0.37 / 0.33 / 0.47 | 3.2 / 3.7 / 4.2 / 4.9 |
| Weibull depot | 1.3 / 0.69 / 1.3 / 1.3 | 3.3 / 4.8 / 5.0 / 5.0 |

With these placeholders the amino acids (5 mM) are the limit: all supply profiles that deliver the 40 mM end at the
same protein yield (about 5 mM residues) and differ only in how fast they get there, and the 10 mM start alone is
energy-limited (ATP exhausted, 2.8 mM). So the supply profile changes the *timing* here, not the final yield.
Whether that holds depends on the real ratio of energy supply to amino acids and on the lysate's ATP-consuming side
reactions, which this model does not have.

## What would make it real

The lysate's enzyme concentrations (total protein mg/mL and the PGM / ENO / PYK / ribosome content), the measured ATP
and protein (e.g. GFP) time courses of your experiment, and the other ATP sinks of the lysate.
