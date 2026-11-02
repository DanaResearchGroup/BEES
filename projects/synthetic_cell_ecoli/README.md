# Synthetic-cell energy metabolism on the real E. coli database

Nothing curated: `database.name: ecoli` (db/ecoli.csv, 2332 rows). The input gives the CFPS recipe (3-PGA 40 mM, ATP 1.2 mM,
amino acids 2.5 mM each) and seven E. coli enzymes (GpmM, Eno, Pyk, Pgk, SerA, SerC, SerB) with **placeholder** enzyme levels and
kinetics (the lysate's enzyme content is unknown; NAD, NADH and 2-oxoglutarate are guesses). The core/edge enlargement finds
the reactions itself.

## What the enlargement did (python BEES.py -i projects/synthetic_cell_ecoli/input.yml)

Seeded with 3-PGA alone it found only 2 reactions (3-PGA + ATP -> 1,3-BPG + ADP by Pgk, and 3-PGA + NAD -> 3-phosphooxypyruvate by SerA),
because db/ecoli.csv lists phosphoglycerate mutase as 2-PGA -> 3-PGA and pyruvate kinase as pyruvate + ATP -> PEP,
and BEES finds a row from its reactant side only (reversed matches are skipped by design).
With 2-PGA and pyruvate also seeded (concentration 0, allowed now), it found 5 reactions in 3 iterations, promoted the products to
the core, resolved all thermodynamics (0 fallbacks, equilibrator) and ran the reversible ones both ways:

| reaction (database direction) | rate law | Keq |
|---|---|---|
| 2-PGA <=> 3-PGA (GpmM) | reversible | 5.7 |
| 2-PGA <=> PEP + H2O (Eno) | reversible | 4.3 |
| 3-PGA + ATP <=> 1,3-BPG + ADP (Pgk) | reversible | 4.4e-4 |
| 3-PGA + NAD -> 3-phosphooxypyruvate + NADH (SerA) | forward | |
| pyruvate + ATP -> PEP + ADP (Pyk) | forward | |

## Result and limits (placeholders)

After 2 h: 3-PGA falls from 40 to 20 mM and PEP rises to 15 mM, 2-PGA 3.5 mM, but **ATP stays at about 1.1 mM, no pyruvate forms**.
The step that regenerates ATP, PEP + ADP -> pyruvate + ATP, is the reverse of the database row, and the reaction BEES builds from it is
forward-only (uphill) instead of running backwards, so the energy route stops at PEP. The same direction convention hides
phosphoserine aminotransferase (O-phosphoserine is only a reactant in the database), so the serine route (the phosphoserine
pathway of Ackerman et al. 2026) is not reached from 3-PGA. There is no translation row in the database either, so no protein is made.

The enlargement mechanism works; what blocks the synthetic-cell energy model is how the database rows are oriented.
Options: let BEES also generate the reverse of a row when its thermodynamics make the database direction strongly uphill, or add the missing
rows (pyruvate kinase PEP -> pyruvate, translation) to a small supplementary table, as done in `../synthetic_cell_cfps`.
