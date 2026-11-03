# Synthetic-cell energy metabolism on the real E. coli database

Nothing curated: `database.name: ecoli` (db/ecoli.csv, 2332 rows). The input gives the CFPS recipe (3-PGA 40 mM, ATP 1.2 mM,
amino acids 2.5 mM each) and seven E. coli enzymes (GpmM, Eno, Pyk, Pgk, SerA, SerC, SerB) with **placeholder** enzyme levels and
kinetics (the lysate's enzyme content is unknown; NAD, NADH and 2-oxoglutarate are guesses). The core/edge enlargement finds
the reactions itself.

## What the enlargement does (python BEES.py -i projects/synthetic_cell_ecoli/input.yml)

Seeded with 3-PGA alone it finds only 2 reactions (3-PGA + ATP -> 1,3-BPG + ADP by Pgk, 3-PGA + NAD -> 3-phosphooxypyruvate + NADH by SerA):
db/ecoli.csv lists phosphoglycerate mutase as 2-PGA -> 3-PGA and pyruvate kinase as pyruvate + ATP -> PEP, and BEES finds a row
from its reactant side only (reversed matches are skipped by design, the reaction then runs both ways through the thermodynamics).
So 2-PGA and pyruvate are also seeded, at concentration 0 (allowed now), in the direction of the database rows. Then it finds 5 reactions in 3
iterations, promotes their products to the core and resolves all thermodynamics (0 fallbacks, equilibrator):

| reaction (database direction) | rate law | Keq |
|---|---|---|
| 2-PGA <=> 3-PGA (GpmM) | reversible | 5.7 |
| 2-PGA <=> PEP + H2O (Eno) | reversible | 4.3 |
| 3-PGA + ATP <=> 1,3-BPG + ADP (Pgk) | reversible | 4.4e-4 |
| 3-PGA + NAD <=> 3-phosphooxypyruvate + NADH + H+ (SerA) | reversible | |
| pyruvate + ATP <=> PEP + ADP + H+ (Pyk) | reversible | 4.4e-4 |

**Species names matter.** The database spells the proton `H(+)` and water `Water`. BEES treats `H+` and `H2O` as buffered species that do not need a
product Km, but not `H(+)`; with `H(+)` as the species label every reaction that releases a proton (here pyruvate kinase, SerA)
was silently forced to forward-only, and the ATP-regenerating step PEP + ADP -> pyruvate + ATP never ran. Name them `H+` and `H2O` in the input.

## Result (placeholders), 2 h

3-PGA falls from 40 to 10 mM, PEP rises to 7.9 mM, pyruvate to 10 mM, 2-PGA 1.8 mM (the reversible steps settle at their
equilibria), and **ATP stays at 1.3 mM while ADP drops to about 0.001 mM**: the 3-PGA -> pyruvate route regenerates ATP as in the lab recipe.
Nothing consumes ATP in this network, so ATP does not fall; there is no translation row in the database, hence no protein.

## Still open

- The serine route stops at 3-phosphooxypyruvate: phosphoserine aminotransferase lists O-phosphoserine as its reactant, so it is not reached from 3-PGA.
- Letting BEES also generate the reverse of a database row when it is triggered from the product side would remove the need for the zero-concentration seeds;
  not implemented (it would change which reactions are found in every project).
- A translation (ATP-consuming) step and real enzyme levels: see `../synthetic_cell_cfps`.
