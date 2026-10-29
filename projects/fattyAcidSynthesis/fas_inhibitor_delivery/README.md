# E. coli FAS II: two inhibitors, three delivery profiles

**Placeholder kinetics.** This is a copy of `../fattyAcidSynthesis_ecoli` (same enzymes, `db/ecoli.csv`) with
CatPred replaced by `measured_kinetics`: one rounded kcat per enzyme (FabD 100, FabH 2, FabB 3, FabF 3,
FabG 50, FabA 30, FabZ 30, FabI 20, TesA 1 s^-1), a few cofactor / primer Kms, `default_km` 0.02 mM for acyl-ACPs,
and a high Km for TesA on short acyl-ACPs (butanoyl 20, hexanoyl 20, octanoyl 10, decanoyl 5, (3Z)-decenoyl 5 mM)
to imitate its long-chain preference. The calibrations of the original input are not used.
Nothing here is a result about biology. If CatPred is installed, delete the `measured_kinetics` blocks and set
`estimate_kinetics: true` to use it.

NADPH, NADH, malonyl-CoA and acetyl-CoA are `constant: true`. Without that, malonyl-CoA (0.5 mM) is used up within
half an hour with these kcats, the 6 h output is a batch endpoint, and an inhibitor that shifts the chain-length
distribution even *raises* the number of fatty-acid molecules.

Drugs (species with concentration 0, entering only through `dosing`):

| drug | target | type | Ki | status |
|---|---|---|---|---|
| triclosan | FabI | tight binding | 7e-9 mM | Ward et al. 1999 - **VERIFY** |
| thiolactomycin | FabH / FabB / FabF | competitive vs malonyl-[ACP] | 0.05 / 0.02 / 0.02 mM | **placeholders** |

Cerulenin (covalent) is out of scope. `end_time` 21600 s (6 h); sweeps run with `--time-scale 0.1`
(release times x 0.1: liposome half-lives 12 min / 4.5 h, depot tau 2.4 h); output = the sum of all free
fatty acids (butanoate ... icosanoate).

## Run (the build takes about 10 minutes)

```bash
python devtools/release_sweep.py -i projects/fattyAcidSynthesis/fas_inhibitor_delivery/input.yml \
  --drug thiolactomycin --dose 0.2 --elimination-half-life 1800 --time-scale 0.1 \
  --output-species "butanoate;hexanoate;octanoate;decanoate;(5z)-dodecenoate;dodecanoate;(7z)-tetradecenoate;tetradecanoate;(9z)-hexadecenoate;hexadecanoate;(11z)-octadecenoate;octadecanoate;(13z)-icosenoate;icosanoate"
python devtools/target_scan.py -i ... --output-species "<same>"
```

## Results (these placeholders; output lost over 6 h vs control)

Thiolactomycin 0.2 mM:

| free-drug clearance | free bolus | liposome | slow depot |
|---|---|---|---|
| none | 71.7 % (peak 200 uM) | 36.4 % (129 uM) | 53.2 % (165 uM) |
| half-life 1800 s | 14.0 % (peak 200 uM) | 11.2 % (peak 21 uM) | 15.8 % (42 uM) |

Without clearance the slow release is much better than the bolus, since the drug stays and the earlier it is there the
more it costs. With a 30 min half-life the bolus and the liposome cost about the same output (14.0 vs 11.2 %), the
liposome at one tenth of the peak free drug: delivery then buys a lower peak, not less inhibition of the output.

Triclosan (FabI is 1 uM, Ki 7 nM, so the dose is compared with the enzyme):

| dose | free bolus | liposome | slow depot |
|---|---|---|---|
| 0.0005 mM (sub-stoichiometric) | 0.4 % | 0.1 % | 0.2 % |
| 0.01 mM | 20.0 % | 20.3 % | 27.9 % |

Sub-stoichiometric triclosan barely matters. At 0.01 mM it costs 20-28 % here; in the earlier build every arm lost most
of the output, which these placeholders do not reproduce (FabI controls only part of the flux, below).

Target scan (kcat x 0.1; output lost at 1 / 3 / 6 h):

| enzyme | 1 h | 3 h | 6 h |
|---|---|---|---|
| FabB | 92.7 | 90.0 | 87.7 |
| FabH | 67.1 | 68.7 | 68.7 |
| FabI | 46.6 | 40.2 | 39.0 |
| FabF | 90.2 | 52.4 | 17.7 |
| TesA | 96.5 | 40.0 | 17.0 |
| FabD, FabG, FabZ, FabA | below 7 | | |

FabH and FabB dominate here, FabI is third, FabA / FabZ / FabD / FabG about 0; an earlier build found FabH dominant
(about 88 % at 90 % inhibition), FabI and TesA 32-37 % and FabB / FabF about 0. The ranking depends on the
placeholder kcats (FabB and FabF share 3 s^-1 here), which is the reason to replace them with CatPred or measured values
before reading anything into it. Figures are in `figures/`.

## Still a placeholder

All kcat / Km / enzyme concentrations, both thiolactomycin and the triclosan Ki (to be checked against the primary
source), the doses.
