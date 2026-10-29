# Methotrexate cellular and in-vivo pharmacology quantities for modelling delivery profiles

Research limits: PubMed abstracts were reachable; PMC, ncbi.nlm.nih.gov and nature.com full text were NOT (host not in the session's allowed domains). Numbers below come mainly from abstracts. Items marked "secondary" come from web-search summaries and were not checked against the primary text.

## What intracellular MTX concentration and exposure duration inhibit DHFR and deplete thymidylate (cancer vs normal)?

### Takeaway
Inhibition needs intracellular MTX at or above the DHFR binding capacity, plus an excess of FREE drug, because the MTX-DHFR complex dissociates slowly. In L1210 cells, DNA synthesis (deoxyuridine incorporation) was suppressed to <2% of control when extracellular MTX was >1 uM and DHFR activity was undetectable above 0.5 uM extracellular MTX. In breast cancer cells, 2 uM MTX for 24 h produced polyglutamates in excess of DHFR capacity, and thymidylate synthesis stayed inhibited 24 h after drug removal. No quantitative normal-cell vs cancer-cell comparison was found.

### Cited Findings
- L1210 cells, MTX-DHFR equivalence study: extracellular MTX >1 uM depressed deoxyuridine incorporation into DNA to <2% of control. Incorporation recovered to 62% of control at the point of MTX-DHFR equivalence (all high-affinity sites just saturated). DHFR activity was undetectable at extracellular MTX >0.50 uM. DHFR activity never exceeded 6.09% of control at the equivalence point. 5% of DHFR activity sustained 50% of control UdR incorporation, and 20% of DHFR activity restored control incorporation. The cellular MTX-vs-DNA-synthesis relationship was sigmoid; the authors conclude "free" intracellular MTX is needed for maximal effect. — Bender & Makula, Br J Cancer 1978, PMID 638018, DOI 10.1038/bjc.1978.60
- Because the drug-enzyme complex is slowly dissociable, drug (free) in excess of that needed to bind all DHFR sites is required to compete with substrate. — Jolivet & Chabner, J Clin Invest 1983, PMID 6193143, DOI 10.1172/JCI111048
- ZR-75-B cells, 2 uM MTX for 24 h: polyglutamates (Glu2-5) replaced the parent drug on DHFR and made up 85% of both free and bound intracellular drug at 24 h. — PMID 6193143
- MCF-7 and ZR-75-B, 2 uM MTX for 24 h: polyglutamates were 55.4% (51.9 nmol/g) and 87.6% (62.4 nmol/g) of intracellular drug. MDA-231 (poor accumulator): 32% (4.06 nmol/g) polyglutamates. Retention of polyglutamates in excess of DHFR binding capacity gave prolonged inhibition; after 24 h in 2 uM MTX plus 24 h drug-free, deoxyuridine incorporation was still 30% of control (MCF-7) and 34.7% (ZR-75-B), with a 30% fall in cell number. MDA-231 recovered quickly after 2 uM; it needed 24 h in 10 uM MTX (11.3 nmol/g polyglutamates) for prolonged inhibition. — Jolivet et al., J Clin Invest 1982, PMID 6178756, DOI 10.1172/jci110624
- Folate-cycle kinetic model of MCF-7 (1 uM MTX over 21 h); results depend strongly on the assumption that MTX di/triglutamates inhibit thymidylate synthase in intact cells as much as with purified enzyme. Useful as a ready model structure with fitted pool sizes. — Morrison & Allegra, J Biol Chem 1989, PMID 2732237 (no DOI listed in PubMed record)

### Inferences
- A delivery model should treat the effective threshold as [MTX-free]_intracellular > 0 after total intracellular MTX exceeds total DHFR (titration). Below equivalence, DHFR is not fully inhibited; 95% DHFR inhibition is needed to halve DNA synthesis (from PMID 638018), i.e. DHFR is in large excess over what is rate-limiting.
- Polyglutamation extends the effective exposure time: dissociation from DHFR slows with chain length (see next section), so exposure duration of ~24 h at ~1-2 uM extracellular produces sustained effect in cells that polyglutamate efficiently.

### Gaps
- No primary-source numbers found for a normal cell vs tumour cell comparison of intracellular MTX threshold, polyglutamation or recovery (the MDA-231 vs MCF-7/ZR-75-B comparison is cancer vs cancer only).
- Allegra 1985-87 papers on intracellular free MTX/DHFR titration were not retrieved (search by author returned only PMID 3860829, 2413074 and 2732237). Allegra et al. J Biol Chem 1985;260:9720 (TS inhibition by MTX polyglutamates) was cited in a secondary summary but its PMID/abstract was not verified.

## Intracellular free MTX vs polyglutamates, and polyglutamate dissociation kinetics

### Takeaway
Polyglutamates dominate intracellular drug after 24 h at 1-2 uM in sensitive breast cancer lines (55-88% of total, ~50-60 nmol/g wet weight), and longer chains dissociate more slowly from DHFR (t1/2 12 min for Glu1 to ~120 min for Glu5) and are retained far longer.

### Cited Findings
- Dissociation t1/2 of bound MTX-Glu1 12 min, Glu2 30 min, Glu3 102 min, Glu4 108 min, Glu5 120 min (displacement by >=5-fold excess unlabeled MTX). After transfer to drug-free medium: >90% of Glu1 left within 1 h, >90% of Glu2 within 6 h, >90% of Glu3 within 24 h; free Glu4 fell only 63% and bound Glu4 only 23% in 24 h. — PMID 6193143, DOI 10.1172/JCI111048
- MCF-7 and ZR-75-B retained 22% and 51% of total intracellular drug after 24 h drug-free. — PMID 6178756
- Polyglutamate potency on other enzymes: MTX and DHF monoglutamates are weak AICAR transformylase inhibitors (Ki 143 uM and 63 uM); MTX Glu4/Glu5 Ki 5.6e-8 M (about 2500-fold more potent than MTX); DHF-Glu5 Ki 4.3e-8 M. — Allegra, Drake, Jolivet, Chabner, PNAS 1985, PMID 3860829, DOI 10.1073/pnas.82.15.4881
- Review framing "Is MTX a prodrug?" (no abstract available in PubMed). — Chabner et al., J Clin Invest 1985, PMID 2413074, DOI 10.1172/JCI112088
- Secondary (grant progress report, not verified): several-hundred to thousand-fold greater inhibition of thymidylate synthase and AICAR transformylase by higher MTX polyglutamates. — [Grantome Z01 CM006513](https://grantome.com/grant/NIH/Z01-CM006513-09)

### Inferences
- Units are nmol per g (wet weight) in these papers; with ~1 g/mL cell density and ~70% water, 50-60 nmol/g is roughly 70-90 uM total intracellular MTX species (inference, not stated in source).

### Gaps
- Free MTX (unbound) concentration in uM was not reported in the retrieved abstracts.

## Cellular uptake (RFC / SLC19A1) kinetics

### Takeaway
Only relative kinetics were retrieved; absolute Km (commonly cited ~5-7 uM for L1210) and Vmax values were not found in the accessible sources and are NOT confirmed here.

### Cited Findings
- RFC1 transfected into carrier-deficient L1210 (MTXrA-R16): influx increased ~9-fold without change in Km; efflux rate constant up 5.1-fold; steady-state free intracellular MTX only 2.1-fold higher at 1 uM extracellular MTX. Concentrative advantage ([MTX]i in R16 / L1210) 1.8 at 0.1 uM and 3.8 at 30 uM external MTX. Polyglutamates 2-fold higher; MTX IC50 about 50% lower. — Zhao...Goldman, J Biol Chem 1997, PMID 9261128, DOI 10.1074/jbc.272.34.21207
- Folate receptor-beta mediated MTX influx in F2-MTXrA cells: influx Kt one-third and Vmax one-fourth that of RFC; 0.2 uM MTX influx one-third of parental L1210. — Sierra...Goldman, Biochem Pharmacol 1995, PMID 7488246, DOI 10.1016/0006-2952(95)94097-y
- Secondary (not verified): S180 cell influx Km and transporter KD ~4-fold higher than L1210 due to a codon-297 RFC1 change; L1210 variants with 3-14-fold higher influx had higher Vmax with unchanged Km (Yang et al., J Biol Chem 1984;259:13139, per search summary).

### Inferences
- Because RFC is bidirectional with a large efflux, steady-state free MTX is only a few-fold above extracellular at uM levels, so free MTX should be modelled as influx minus efflux minus DHFR/polyglutamate sinks, not as accumulation against a large gradient.

### Gaps
- Absolute Km (uM), Vmax (pmol/min/mg or per 10^? cells) for RFC-MTX influx not obtained; PMC/full text unreachable.

## DHFR concentration in cells (uM) and drug:enzyme ratio

### Takeaway
No direct uM value for cellular DHFR was found. Only activity-based fold changes and the qualitative result that MTX-DHFR equivalence is reached at roughly sub-uM to ~1 uM extracellular MTX in L1210.

### Cited Findings
- L1210: at MTX-DHFR equivalence, DHFR activity <=6.09% of control and extracellular MTX about 0.5-1 uM gave full DHFR inhibition. — PMID 638018
- Secondary (unverified, from search summary): DHFR activity reported in pmol/mg protein for human melanoma lines; conversion to uM needs cell volume (not given).

### Inferences
- Total DHFR can be estimated from the MCF-7/ZR-75-B polyglutamate totals above (tens of nmol/g) only as an upper bound, since polyglutamates exceed DHFR capacity by design. A modeller should treat DHFR as an adjustable parameter in the roughly 0.1-few uM range (inference only, not sourced); this must be flagged as assumed.

### Gaps
- DHFR content per cell or in uM (e.g. from stoichiometric titration with [3H]MTX in intact cells) was not found. Harrap/Jackson-type binding-capacity papers (nmol per 10^9 cells) were not retrieved.

## DHFR-MTX Ki, tight-binding analysis, and IC50/GI50

### Takeaway
Human DHFR Ki for MTX is about 3.4 pM (slow, tight-binding competitive). IC50/GI50 values are only available from vendor-level summaries (HCT116 GI50 0.015 uM at 72 h).

### Cited Findings
- Human DHFR: MTX is a slow, tight-binding competitive inhibitor, KiMTX = 3.4 pM. — Volpato...Pelletier, J Biol Chem 2009, PMID 19478082, DOI 10.1074/jbc.M109.018010 (value read from the paper's introduction via search summary; not independently verified)
- Other species (secondary, summary only): S. faecium DHFR overall Ki about 58 pM; N. gonorrhoeae Ki 13 pM. Not for human modelling.
- HCT116 GI50 0.015 uM (72 h, MTT) — vendor listing, original paper not identified: [Selleck methotrexate](https://selleck.cn/products/Abitrexate.html)
- A549 IC50 0.25 uM (72 h MTT) and 0.99 uM with leucovorin — compiled listing, low quality: same search summary.
- Resistant CCRF-CEM R30dm EC50 0.0079-0.0145 uM at 120 h (vendor listing, ambiguous).

### Inferences
- With Ki ~ pM and enzyme in sub-uM to uM range, MTX behaves as a stoichiometric (titration) inhibitor: Morrison quadratic regime, so [I]total must exceed [E]total for inhibition. Effective potency in cells is governed by [E]total and by uptake, not by Ki.

### Gaps
- Morrison & Allegra 1989 retrieved (PMID 2732237) is a folate-cycle model, not an enzyme Ki paper; the specific Morrison-equation analysis papers and enzyme concentrations used in Ki assays were not found.
- No normal-cell (fibroblast, marrow, intestinal epithelium) IC50 values found. No MCF-7 or CCRF-CEM 72 h IC50 found.

## Plasma pharmacokinetics (low vs high dose)

### Takeaway
High-dose infusion gives biphasic decline with distribution t1/2 ~2.4 h and terminal t1/2 ~20 h; low-dose terminal t1/2 is shorter (3-10 h per a review summary).

### Cited Findings
- 34 patients (1-25 y), 12 g/m2 over 4 h: peak serum MTX 674-1778 uM (mean 1201 +/- 293 uM); t1/2 alpha 2.42 +/- 0.45 h; t1/2 beta 19.9 +/- 7.6 h; Vss 0.56 +/- 0.18 L/kg; clearance 71 +/- 20 mL/min/m2. Low-risk thresholds (97.5th percentile): <=6.3 uM at 24 h, <=0.77 uM at 48 h, <=0.33 uM at 72 h after end of infusion. — Raude et al., Int J Clin Pharmacol Ther Toxicol 1988, PMID 3209285 (no DOI in PubMed record)
- Secondary (not primary-checked): high-dose t1/2 8-15 h vs 3-10 h at low dose; 5 g/m2 over 24 h gave steady-state ~66 uM with three half-lives 1.8 h, 6.4 h and ~15 h; Japanese infusion study mean 330 uM at 6 h, 72 uM at 12 h, 16 uM at 24 h (low clearance) vs 210, 18 and 1.0 uM (high clearance); plasma >10 uM at 48 h or >0.1 uM at 48 h trigger intensified leucovorin rescue. Sources: [Hiroshima repository](https://ir.lib.hiroshima-u.ac.jp/37841/files/33706); PMC1885068 (adult cohort, 300 mg/m2-12 g/m2; unread because PMC unreachable).

### Gaps
- Low-dose oral/IM plasma concentrations (e.g. 7.5-25 mg weekly, peak ~0.5-1 uM) and FDA label values were not retrieved.
- Tissue/tumour interstitial MTX levels not found.

## Liposomal MTX release profiles and half-lives

### Takeaway
Only fragmentary, secondary data: DPPC liposomes show biphasic first-order release with first-phase time ~40 h (non-lyophilised) or ~15 h (lyophilised); half-lives lengthen ~2.5-fold when complexed with DEAE-dextran. Temperature and medium for that figure were not stated. Schroeder-lab data not found.

### Cited Findings
- DPPC liposomes: biphasic first-order MTX release, first-phase ~40 h non-lyophilised vs ~15 h lyophilised — [Purdue dissertation AAI8210268](https://docs.lib.purdue.edu/dissertations/AAI8210268) (abstract only, search summary)
- Sodium MTX liposomes: 33% MTX lost after 1 month at 25 C; ~70% retention after freeze-drying with trehalose/PEG-20000 — [PDA J Pharm Sci Technol 52(1):23](https://journal.pda.org/content/52/1/23) (search summary)
- Intravitreal MTX liposomes in pig: therapeutic levels sustained >6 weeks (in vivo, Adv Healthc Mater 2025) — URL not captured.

### Gaps
- No serum 37 C release t1/2 for modern (e.g. PEGylated, remote-loaded) liposomal MTX; no PMID/DOI for the sources above. Plasma t1/2 of liposomal MTX in vivo not found.

## DHFR amplification fold in resistant lines

### Takeaway
Typical resistant sublines show about 10-20-fold higher DHFR activity with ~18 gene copies, with 75-fold MTX resistance in one CCRF-CEM subline; all from secondary summaries.

### Cited Findings
- CCRF-CEM/Ri: ~20-fold DHFR activity, 75-fold MTX resistance, ~18 DHFR gene copies; CEM/R3: DHFR activity 10-20x parent. — search summary only (underlying papers not identified; see query results at [Rutgers record](https://www.researchwithrutgers.org/en/publications/molecular-and-karyological-analysis-of-methotrexate-resistant-and/))

### Gaps
- No PMID/DOI verified; classic Alt/Schimke L1210 and 3T6 amplification numbers (up to several-hundred-fold) were not retrieved.
