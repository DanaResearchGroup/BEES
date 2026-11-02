# Synthetic cell (cell-free protein synthesis): energy and protein output

**The recipe is real, the enzyme side is a guess.** The inner-solution concentrations come from Table S.2 of the
final-year project report (final requested inner solution volume for the synthetic cell, lab table of 23/03/2024):
HEPES-KOH pH 8 55 mM, Mg-acetate 14, K-acetate 50, NH4-acetate 155, PEG 6000 3 %, 3-PGA 40, 20 amino acids at 2.5 mM each (mixtures I and II, 50 mM in all), ATP 1.2, GTP 1, UTP 0.8, IPTG 1, sucrose 200 mM, S30-T7 lysate 34 %, plasmid 10 ug/mL (the NH4-acetate
and amino-acid stock concentrations were hard to read in the photographed table). The report gives no enzyme content of
the lysate, so **every enzyme concentration, kcat and Km here is a placeholder**. Nothing below is a result about
the real system.

## What is modelled (`db/cfps.csv`)

3-PGA -> 2-PGA (PGM) -> PEP (ENO) -> pyruvate + ATP from ADP (PYK), and one lumped translation step:
`amino acids + 5 ATP -> protein residue + 5 ADP + 5 Pi` (about 5 ATP-equivalents per residue, Ackerman et al. 2026; GTP not separate). The output
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

(An earlier version of this README took the amino acids as 5 mM in all and concluded that they limit the yield. That was
wrong: "2.5 mM" is per amino acid, 50 mM in all, so the energy is what limits.)

| arm | ATP at 0.5 / 1 / 2 / 4 / 8 h (mM) | protein residue at 0.5 / 1 / 2 / 4 / 8 h (mM) |
|---|---|---|
| 10 mM only (no supply) | 1.29 / 0.32 / 0.15 / 0.11 / 0.08 | 1.4 / 2.2 / 2.2 / 2.2 / 2.2 |
| all 30 mM at t = 0 (the recipe) | 1.29 / 1.29 / 1.29 / 0.16 / 0.09 | 1.4 / 2.8 / 5.6 / 8.2 / 8.2 |
| liposome-like slow release | 1.29 / 0.60 / 0.41 / 0.36 / 0.32 | 1.4 / 2.6 / 2.9 / 3.4 / 4.0 |
| Weibull depot | 1.29 / 1.29 / 0.50 / 0.42 / 0.35 | 1.4 / 2.8 / 3.9 / 4.6 / 5.4 |

With these placeholders the energy supply limits the yield: delivering the whole 40 mM at once gives the most protein and the
fastest, and slower supply gives less within 8 h because the ATP stays low and translation is slow. (In the sweep "output lost"
is negative = a gain.)

## Compared with the published experiments (Schroeder lab)

- Same recipe, lysate and yield: Krinsky et al., "A simple and rapid method for preparing a cell-free bacterial lysate for
  protein synthesis" (S30-T7, lysate total protein set to 22 mg/mL and used at 30 % v/v, so about 6.6 mg/mL in the reaction): 40-150 ug/mL
  protein in 2 h at 37 C. At about 110 Da per residue that is **about 0.4-1.4 mM of residues** (my conversion).
  The model gives 5.6 mM at 2 h for the recipe, about 4x too much: the real lysate wastes ATP elsewhere, which this model
  does not have.
- Ackerman et al. 2026, *ACS Synth Biol* ([DOI 10.1021/acssynbio.5c00649](https://doi.org/10.1021/acssynbio.5c00649)): translation costs the
  equivalent of about 5 ATP per amino acid (used here); ATP falls as the reaction proceeds and protein production plateaus because
  ATP runs out (adding ATP gives +14 %, adding 3-PGA +82 %); activity lasts about 190 min in synthetic cells at 37 C and the
  optimum is 2.4 mM ATP with 70 mM 3-PGA in bulk; higher phosphoserine extends the activity. The ATP time courses are in the
  paper's figures (luminescence), not in its text, so I could not take numbers from them.
- Not found: your own ATP and protein measurements (see below).

## What would make it real

The ATP and protein (GFP) time courses of your experiment, the lysate's PGM / ENO / PYK / ribosome content, and the lysate's other ATP sinks.
With the published 0.4-1.4 mM residues in 2 h and an ATP plateau after about 3 h, a lumped ATP-consuming step could be calibrated.
