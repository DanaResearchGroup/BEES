# 2-deoxyglucose against hexokinase

**Placeholder kinetics, and a placeholder Ki.** ChEMBL has no hexokinase Ki for 2-deoxyglucose, so
Ki = 0.3 mM is a stand-in; kcat / Km values are order-of-magnitude numbers typed into `measured_kinetics`
(CatPred is not used). This is an example of the machinery, not a result about biology.

## System

`db/glycolysis.csv` (curated): HK `D-glucose + ATP -> G6P + ADP`, PGI `G6P <=> F6P`,
PFK `F6P + ATP -> FBP + ADP`. Glucose starts at 0.2 mM; ATP and ADP are constant; FBP is the output.
2-DG is a competitive inhibitor of HK against glucose (Km_app = Km (1 + I/Ki)), dose 5 mM,
free-drug half-life 30 min.

2-DG is really *phosphorylated* by hexokinase and 2-DG-6-phosphate inhibits further steps; this model
captures only the competition for the glucose site.

## Run

```bash
python BEES.py -i projects/glycolysis_2dg/input.yml
python devtools/release_sweep.py -i projects/glycolysis_2dg/input.yml --drug 2-deoxyglucose --output-species "D-fructose 1,6-bisphosphate"
python devtools/target_scan.py   -i projects/glycolysis_2dg/input.yml --output-species "D-fructose 1,6-bisphosphate" --times-h 0.25,0.5,1
```

## Result (these placeholders)

Output lost (FBP, vs control), equal 5 mM dose:

| glucose, horizon | free bolus | liposome |
|---|---|---|
| 0.2 mM, 2 h (input) | 0.1 % | 0.0 % |
| 0.2 mM, 5 min | 81.9 % | 0.6 % |
| 0.2 mM, 15 min | 77.8 % | 2.1 % |
| 5 mM, 5 min | 20.4 % | 0.0 % |
| 5 mM, 15 min | 19.9 % | 0.1 % |

Over the 2 h run the effect is small because the glucose (0.2 mM) is converted in the end anyway; 2-DG only delays it,
and only a bolus (5 mM = 17 Ki) delays it noticeably. At 5 mM glucose the early bolus effect falls from about 80 % to about 20 % (Km_app 1.8 mM against 5 mM substrate): it weakens
but is not gone for a 5 mM bolus; the liposome never reaches an effective concentration. The target scan (kcat x 0.1, FBP lost at
15 min / 30 min / 1 h): HK 86 / 76 / 53 %, PFK 77 / 57 / 20 %, PGI 55 / 21 / 1 %.

## Still a placeholder

Everything quantitative: kcat, Km, enzyme levels, the 2-DG Ki.
