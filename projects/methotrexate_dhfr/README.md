# Methotrexate on DHFR: same dose, different release

**Every kinetic parameter here is an order-of-magnitude placeholder** (typed into `measured_kinetics`;
CatPred is not used). The numbers below show how the model behaves, not what a cell does.

## System

`db/folate.csv` (curated): SHMT `serine + FH4 -> glycine + 5,10-methylene-FH4`, TS
`dUMP + 5,10-methylene-FH4 -> dTMP + FH2`, DHFR `FH2 + NADPH -> FH4 + NADP`. The folate pool
(0.01 mM, all FH4 at t = 0) is conserved and cycles; dTMP accumulation is the output flux. dUMP, serine,
glycine, NADPH, NADP, H+ and H2O are `constant: true`. (The folate is called FH4, not THF: BEES lists
THF as a general cofactor, which cannot seed reaction discovery.)

Methotrexate: `tight_binding` on DHFR (Morrison, uses [DHFR] = 1e-4 mM), Ki = 1e-7 mM (0.1 nM).
ChEMBL, target CHEMBL202, molecule CHEMBL34259, Ki, binding assays, relation "=", no validity flag:
activities 3055901 (0.077 nM), 939069 (0.179 nM), 18251252 (0.01 nM), 228406 (0.11 nM, a Leu22-Phe
mutant), 939068 (0.038 nM, listed for *T. cruzi* DHFR) and 2219757 (4420 nM). The median of the first
three is 0.077 nM and of all six 0.094 nM; 0.1 nM is used. Checked in ChEMBL v34 through the ChEMBL
tools; `devtools/fetch_inhibitors.py` reproduces the table when ebi.ac.uk is reachable.

Dose 2 uM (0.002 mM) over 72 h, free-drug half-life 10000 s (the same in every arm).
The Biomodels model BIOMD0000000018 (Morrison & Allegra 1989) was not compared (not reachable here).

## Run

```bash
python BEES.py -i projects/methotrexate_dhfr/input.yml        # liposome arm of the input
python devtools/release_sweep.py -i projects/methotrexate_dhfr/input.yml --drug methotrexate --output-species dTMP
python devtools/target_scan.py   -i projects/methotrexate_dhfr/input.yml --output-species dTMP --times-h 6,24,72
```

## Result (these placeholders)

`release_sweep.py`, 72 h, equal 2 uM total dose:

| arm | dTMP lost vs control | h below half control flux | inhibition AUC (h) | peak free drug |
|---|---|---|---|---|
| free bolus | 16.7 % | 12.1 | 15.9 | 2.00 uM |
| liposome (biexponential, 10 % at 2 h / 90 % at 45 h) | 22.1 % | 16.0 | 50.2 | 0.15 uM |
| slow depot (Weibull, tau 24 h, beta 0.6) | 31.5 % | 22.9 | 44.5 | 0.32 uM |

The liposome costs more output than the bolus at 1/13 of the peak free drug. An earlier build of this
example (different placeholders) gave about 58 % for the liposome against about 19 % for the bolus at
the same 1/13 peak; the ordering and the peak ratio are reproduced, the loss of the liposome is not.
That number depends strongly on one parameter, the enzyme concentration, because the liposome holds the free
drug at about the same level as DHFR (0.1 uM) and tight binding is switch-like there:

| [DHFR] (mM) | bolus | liposome | slow depot |
|---|---|---|---|
| 3e-5 | 24.4 % | 98.6 % | 87.2 % |
| 5e-5 | 21.0 % | 83.7 % | 58.0 % |
| 1e-4 (input) | 16.7 % | 22.1 % | 31.5 % |

With TS and SHMT at 1e-5 mM instead of 5e-4 (DHFR then has far more capacity than the pathway needs) a steady
low exposure does nothing: bolus 11.1 %, liposome 1.3 %, depot 13.0 %. Whether the release profile matters thus
depends on how close the target enzyme is to controlling the flux, which is what `target_scan.py` shows:
kcat x 0.1 for TS or SHMT loses about 84 % of the dTMP, for DHFR about 14 %.

## Exploring enzyme levels (healthy-like vs cancer-like)

`devtools/enzyme_level_scan.py` scales DHFR (x0.3 ... x30 of 1e-4 mM) and TS/SHMT (x0.2, x1, x3). dTMP lost over 72 h:

| DHFR x | TS/SHMT x | bolus | liposome | slow depot |
|---|---|---|---|---|
| 0.3 | 1 | 24.4 | 98.6 | 87.2 |
| 1 | 1 | 16.7 | 22.1 | 31.5 |
| 3 | 1 | 10.4 | 0.0 | 3.9 |
| 10 | 1 | 3.6 | 0.0 | 0.0 |
| 30 | 1 | 0.0 | 0.0 | 0.0 |
| 1 | 0.2 / 3 | 15.7 / 17.7 | 15.0 / 33.6 | 26.9 / 36.7 |

The DHFR level sets the outcome: a tight binder is titrated away by the enzyme, so the low free drug of a slow release
stops inhibiting first when DHFR is raised (as in DHFR-amplified, methotrexate-resistant cells), whereas the
bolus keeps a short window of full inhibition. With little DHFR the slow release wins. TS/SHMT matter much less.
Literature found through PubMed (not absolute numbers): DHFR is up-regulated in breast tumours against adjacent
normal tissue (A-to-I editing, [10.1074/jbc.M117.775684](https://doi.org/10.1074/jbc.M117.775684)) and amplified in methotrexate-resistant colon
cancer cells ([10.1038/s41416-024-02664-0](https://doi.org/10.1038/s41416-024-02664-0)). I did not find absolute enzyme or metabolite
concentrations for healthy vs tumour tissue with the search tools available; the fold changes above are a scan,
not measured values. Next step: take absolute levels from a proteomics resource (for example copy-number datasets) and
metabolomics tables for the tissues of interest.

## Still a placeholder

All kcat / Km values and enzyme / metabolite concentrations; the Ki is a literature median across assays
of different quality (see above). Absolute dTMP amounts are meaningless.
