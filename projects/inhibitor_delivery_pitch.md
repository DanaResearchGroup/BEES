# BEES for drug-delivery design: what the drug does, once you know when it arrives

For Prof. Avi Schroeder (Technion, Chemical Engineering). One page.

## The gap

Your lab controls **when** a drug reaches a cell: carrier, release profile, payload. What a given exposure
**does** to an enzyme pathway is usually left to a cell assay or to a hand-written ODE model for one pathway.
BEES (rate-based model generation for enzymatic networks, `github.com/DanaResearchGroup/BEES`) builds the pathway
model from a reaction database and kinetic estimates, and now takes two new inputs: reversible **enzyme inhibition**
(competitive, uncompetitive, noncompetitive, mixed, tight-binding) and a time-dependent **dosing** term (bolus,
constant, first-order, biexponential, Weibull, Korsmeyer-Peppas, tabulated rate; with clearance of the free drug).
The question it can answer: *for the same total dose, which release profile costs the pathway most - or least?*

## Three examples, one tool (`devtools/release_sweep.py`)

All kinetics are order-of-magnitude placeholders, not results about biology (each project README says what is a placeholder).
Same total dose in every arm; liposome = biexponential (10 % with a 2 h half-life, 90 % with 45 h; x0.1 for the 6 h example);
slow depot = Weibull (tau 24 h, beta 0.6).

| system | drug, dose | free bolus | liposome | peak free drug (bolus : liposome) |
|---|---|---|---|---|
| thymidylate cycle (SHMT, TS, DHFR), 72 h | methotrexate, 2 uM, tight binding Ki 0.1 nM (ChEMBL median) | 16.7 % dTMP lost | 22.1 % | 2.0 : 0.15 uM (13 : 1) |
| upper glycolysis (HK, PGI, PFK), 2 h | 2-deoxyglucose, 5 mM, competitive Ki 0.3 mM (placeholder) | 0.1 % FBP lost | 0.0 % | 5000 : 121 uM |
| E. coli FAS II, 6 h | thiolactomycin, 0.2 mM, competitive vs malonyl-ACP (placeholder Kis), free-drug half-life 30 min | 14.0 % fatty acids lost | 11.2 % | 200 : 21 uM (about 10 : 1) |

What the three show:

- **Methotrexate:** the liposome loses *more* output than the bolus at 1/13 of the peak drug, because tight
  binding makes the effect depend on whether free drug stays above the enzyme concentration (0.1 uM), not on the peak.
  The size of the effect depends strongly on [DHFR] (liposome 22 % at 0.1 uM, 84 % at 0.05 uM), so the enzyme level is the number to measure.
- **2-DG:** a competitive drug against a high-affinity substrate is nearly harmless over the run; a 5 mM bolus delays
  the first minutes by about 80 % (about 20 % at 5 mM glucose), the liposome never gets there. (2-DG is also a substrate of
  hexokinase in reality; the model only has the competition.)
- **FAS:** without clearance the profile matters a lot (thiolactomycin: bolus 72 %, depot 53 %, liposome 36 %); with a 30 min clearance
  the bolus and the liposome cost similar output, the liposome with a tenth of the peak. Triclosan below the FabI concentration
  does almost nothing (0.1-0.4 %). `target_scan.py` ranks enzymes by control of the output (here FabB and FabH above FabI;
  FabA/FabZ about nothing).

## Hook: your release model

Enden & Schroeder, *Ann Biomed Eng* 2009 (PMID 19731036) models release from liposomes; its release curves are the
kind of input BEES takes (`dosing: profile: table` takes any measured release curve, scaled to the dose).
A measured release profile for your formulation, plus a Ki for the payload from ChEMBL, gives the predicted
pathway cost in minutes.

## Follow-ups

- Collagozome model (PMID 31503443): release inside a matrix as a dosing curve.
- Synthetic-cell energy metabolism (PMID 42312441): a defined small network to test predicted vs measured output.
- Barcoded-nanoparticle screens (PMID 27830705): many formulations at once - `release_sweep.py` is one formulation per row of a table.
- Replace the placeholders: CatPred kcat/Km/Ki or measured values, then compare with an assay in cells or a cell-free extract.

## Data and licences

| source | used for | licence |
|---|---|---|
| ChEMBL | inhibition constants (`devtools/fetch_inhibitors.py`) | CC BY-SA 3.0 (attribution + share-alike) |
| BRENDA | kinetic parameters, reaction database (BKMS-react) | CC BY 4.0 |
| SABIO-RK | kinetic parameters | academic use only |

Methotrexate Ki (0.1 nM) is the median of unflagged ChEMBL binding records for human DHFR; the triclosan Ki (7 nM)
is from the literature and is marked to be verified; the 2-DG and thiolactomycin constants are placeholders.
