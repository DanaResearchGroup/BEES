# Absolute abundance and activity of DHFR, TYMS, SHMT1/SHMT2 in healthy vs cancerous human tissue/cells

Scope note: the session could not fetch NCBI/PMC, EuropePMC or Spandidos full text (hosts not in the allowed network list). Numbers below come from PubMed abstracts retrieved through the PubMed tool and from search-result summaries. The big proteomics datasets (Wisniewski/Mann proteomic ruler, Nagaraj/Geiger, CPTAC, PaxDb) could NOT be queried, so per-cell copy numbers for DHFR/TYMS/SHMT are a major gap. Nothing below is guessed; unsourced arithmetic is labelled as an assumption-based estimate.

## 1. What are typical absolute concentrations (uM) or copy numbers per cell of DHFR, TYMS, SHMT in normal vs tumor cells?

### Takeaway
I found no citable per-cell copy number or uM value for DHFR, TYMS, SHMT1 or SHMT2 in any human cell. Only tissue-level TS "binding site" values (fmol/mg protein, pmol/g tissue) are available, and these convert to roughly 1-35 nM TS. The proteomic-ruler datasets that would give copy numbers were not reachable.

### Cited Findings
- TS measured by FdUMP-binding assay in resected gastric cancer ranged 7.1-176.6 fmol/mg protein; in human tumor xenografts (gastric, colon, pancreatic, breast) 84.7-775.5 fmol/mg protein — [Ishikawa 1999 Anticancer Res, PMID 10697632](https://pubmed.ncbi.nlm.nih.gov/10697632/); [Ishikawa 1999 Clin Cancer Res, PMID 10213225](https://pubmed.ncbi.nlm.nih.gov/10213225/)
- TS by chemiluminescent ELISA (TS106 antibody), 7 breast cancer/sarcoma samples: mean 131 fmol/mg cytosolic protein (range 60-240); ELISA read uniformly higher than FdUMP binding — search-result summary of the 1993 ELISA method paper (original paper not identified/verified); see [search summary context via Merck/Millipore TS106 page](https://merckmillipore.com/SE/en/product/Anti-Thymidylate-Synthase-Antibody-clone-TS106,MM_NF-MAB4130) (weak source, low confidence on attribution)
- TS in pmol/g tissue (FdUMP binding), 13 colorectal cancer patients: primary tumor 6.6 +/- 1.40, liver metastasis 3.2 +/- 0.83, normal colonic mucosa 2.5 +/- 0.83 (mean +/- SE) — [Umekita 2000 Gan To Kagaku Ryoho, PMID 11086436](https://pubmed.ncbi.nlm.nih.gov/11086436/)
- Larsson 1996 (cited only via a search summary, not verified): primary colorectal cancer 4.2 pmol/g vs normal liver 0.58 pmol/g — no PMID obtained; treat as unverified.
- Total protein density context: cells contain roughly 1-3 million protein molecules per um^3, and published per-cell total-protein copy numbers in HeLa differ up to ~10-fold between studies — [Dolgalev & Poverennaya 2023 review, PMC9916689 (as listed in search results)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9916689). Per-enzyme values not given there.
- DHFR2 protein is undetectable in most adult cell types; mitochondrial DHFR activity is attributed to DHFR itself — [PMC10884974](https://pmc.ncbi.nlm.nih.gov/articles/PMC10884974)
- HeLa proteomic ruler method (histone-normalised copy numbers): [Wisniewski et al. 2014, PMID 25225357](https://pubmed.ncbi.nlm.nih.gov/25225357/) — method only; enzyme values not retrieved.

### Inferences
- Conversion rule (pure arithmetic): 1 uM = ~602 molecules/um^3, so a cell of 2,000 um^3 (ASSUMED typical adherent line volume, uncited) has ~1.2 million copies per uM.
- TS FdUMP-binding values converted to concentration, ASSUMING intracellular protein of ~200 mg/mL (uncited typical value): 7-177 fmol/mg -> ~1.4-35 nM; 131 fmol/mg -> ~26 nM. If tissue density ~1 g/mL, 2.5-6.6 pmol/g -> ~2.5-6.6 nM (assumes all TS sites in tissue water volume). So TS is roughly low-nM in tissue, i.e. ~1,000-20,000 copies per um^3-equivalent scale is not reached; per-cell copies at 2,000 um^3 would be ~2-40 x 10^3 per nM... (not computed further; treat as order-of-magnitude 10^4-10^5 copies/cell, unverified). Cytosolic-extract protein is more dilute than cellular protein, so these are bounds, not measurements.

### Gaps
- No per-cell copy number or uM value for DHFR, TYMS, SHMT1, SHMT2 in any human cell line or tissue (Wisniewski/Mann, Nagaraj, Geiger, PaxDb, CPTAC supplementary tables not accessible in this session). Recommended next step: pull the proteomic-ruler supplementary table for HeLa/other lines and look up DHFR, TYMS, SHMT1, SHMT2 directly.
- No SHMT1/SHMT2 concentration of any kind found.
- No normal-vs-tumor copy-number comparison found.

## 2. Reported tumor/normal fold changes for each enzyme, by tissue

### Takeaway
TS has the best quantitative data: tumor/normal activity is ~1.6x (colorectal) to ~10x (head and neck), and mRNA ~3-3.5x. For DHFR and SHMT1/2 I found only qualitative statements (higher in tumor) and no numeric tumor/normal ratios.

### Cited Findings
- Colorectal (n total 43 GI patients): TS activity median 49.8 tumor vs 30.8 normal mucosa "pmol/mg protein" (ratio 1.62, p=0.023); stomach 81.0 vs 38.0 (ratio 2.13, p=0.012); TS mRNA (TS/GAPDH) colorectum 3.20 vs 0.91 (3.5x, p=0.001), stomach 6.0 vs 2.0 (3.0x, p=0.009). Authors found no correlation between protein, mRNA and activity — [Miyamoto 2001 Int J Oncol, PMID 11251164, DOI 10.3892/ijo.18.4.705](https://doi.org/10.3892/ijo.18.4.705)
- Colorectal (13 pts): TS pmol/g tissue primary tumor 6.6 vs normal mucosa 2.5 (ratio 2.6); liver metastasis 3.2 (1.3x normal mucosa) — [Umekita 2000, PMID 11086436](https://pubmed.ncbi.nlm.nih.gov/11086436/)
- Head and neck SCC (n=37): TS activity 9.2 +/- 21.5 pmol/min/mg protein in tumor vs 0.9 +/- 1.2 in adjacent mucosa (mean ratio ~10x, p=0.0001; huge SD, so ratio of means is skewed); folate pools similar 1.2 vs 0.8 pmol/mg protein (NS) — [Parise 1994 Head Neck, PMID 8021136, DOI 10.1002/hed.2880160210](https://doi.org/10.1002/hed.2880160210)
- Pediatric/biopsy series (7 pts): normal mucosa TS 44-107 pmol/hr/mg protein (=0.73-1.8 pmol/min/mg) at 10 uM dUMP; paired tumors always higher, 52-3000 pmol/hr/mg (~0.9-50 pmol/min/mg; up to ~30x) — Peters 1988 Pediatr Res, as reported in a search summary (nature.com/articles/pr1988508); not independently verified.
- Gastric (n=41): TS activity and mRNA significantly higher in tumor than adjacent mucosa (no numbers in abstract) — [PMID 10697632](https://pubmed.ncbi.nlm.nih.gov/10697632/)
- Lung (168 pts, paired ELISA): TS higher in tumor than non-tumor tissue; lower in adenocarcinoma than other subtypes; numeric values not retrieved — [Shiina 2017 Tohoku J Exp Med](https://www.jstage.jst.go.jp/article/tjem/242/4/242_303/_html/-char/en)
- CRC mRNA by laser capture (76 pts): 30/76 (39%) had higher TS mRNA in tumor vs normal mucosa; no fold values — [Liu 2011 J Hematol Oncol, PMID 21824439](https://pubmed.ncbi.nlm.nih.gov/21824439/)
- DHFR (qualitative): partially purified DHFR from human GI, kidney, lung and breast cancers was higher than in normal tissue in almost all cases; reported MTX-titration values ~5.09 x 10^-14 mol/mg protein (kidney) and 1.07 x 10^-10 mol/mg (stomach) — [Chulalongkorn repository record](https://cuir.car.chula.ac.th/handle/123456789/2154?mode=full). The two values differ by ~2000-fold and units look inconsistent; treat as unreliable until full text is checked.
- SHMT2 (qualitative): mRNA significantly elevated in >85% of cancer types vs normal (TCGA); protein higher in lung adenocarcinoma tumor vs normal (CPTAC); breast IHC: no SHMT2 in distant non-cancerous cells, stained in all 76 cancers — [PMC8408114](https://pmc.ncbi.nlm.nih.gov/articles/PMC8408114/); [PMC8049788](https://pmc.ncbi.nlm.nih.gov/articles/PMC8049788)
- SHMT1 and SHMT2 up-regulated in patient lung cancer tissue (review statement) — see search results above, no values.
- Rat (not human) serine metabolism enzymes: SHMT activity coordinately increased in tissues with high cell turnover and in proliferating human lymphocytes after PHA stimulation (no numbers in abstract) — [Snell 1984 Adv Enzyme Regul, PMID 6089514](https://pubmed.ncbi.nlm.nih.gov/6089514/)

### Inferences
- Activity ratios are modest (1.6-2.6x) in GI tissue but can be ~10x in head and neck; ratios of means are inflated by skewed distributions.
- Methods disagree: Miyamoto reports 49.8 pmol/mg for colorectal TS whereas FdUMP-binding studies report fmol/mg (e.g. 7-177 fmol/mg in gastric). The ~1000-fold gap suggests different assays/units (possibly a unit typo in the abstract); do not combine across studies.

### Gaps
- No numeric tumor/normal ratios for DHFR (activity or protein) in human tissue, by tissue type.
- No numeric SHMT1/SHMT2 protein or activity ratios (only qualitative or mRNA).
- No CPTAC/Human Protein Atlas log2 fold changes retrieved for DHFR, TYMS, SHMT1, SHMT2.
- No numeric TS data for breast, liver, leukemia tissue (only the 7-sample breast/sarcoma ELISA value above).

## 3. Enzyme activity values per mg protein, and conversion to approximate uM

### Takeaway
TS activity in tissue is ~1-50 pmol/min/mg protein (normal mucosa ~0.7-1.8; tumors up to ~9 mean). I found no DHFR or SHMT specific-activity values in human tissue. Using stated assumptions, TS at ~1-10 pmol/min/mg corresponds to ~nM-level enzyme.

### Cited Findings
- TS activity, head and neck: tumor 9.2 +/- 21.5, mucosa 0.9 +/- 1.2 pmol/min/mg protein — [PMID 8021136](https://pubmed.ncbi.nlm.nih.gov/8021136/)
- TS activity at 10 uM dUMP, colon mucosa 44-107 pmol/hr/mg protein; tumors 52-3000 (Peters 1988, via summary; not verified).
- DPD (not in scope but same assays) 31.7-82.1 pmol/min/mg protein in colon/liver metastasis — [PMID 11086436](https://pubmed.ncbi.nlm.nih.gov/11086436/) (shows assay unit conventions only)
- Recombinant human DHFR commercial standard specification: >=1 unit/mg protein (1 U = 1 umol NADPH oxidised/min) — [Sigma D6566 product page](https://sigmaaldrich.com/DM/en/product/sigma/d6566) (not a tissue value)
- Gene-amplified DHFR in mouse L5178Y MTX-resistant cells made up ~10-11% of total soluble cell protein (300-fold over parent) — search summary of old primary literature (https://pubannotation.org/docs/sourcedb/PubMed/sourceid/3856428 and related), mouse not human.

### Inferences
Conversion method (all constants are ASSUMPTIONS, uncited):
1. Cell protein density ~200 mg/mL (range 150-300).
2. [E] (uM) = specific activity (umol/min/mg) x protein density (mg/mL) / kcat (per min).
3. Example, TS in tumor, 9.2 pmol/min/mg = 9.2e-6 umol/min/mg; x 200 mg/mL = 1.8e-3 umol/min/mL (= 1.8e-3 mM-equivalent flux /min... i.e. 1.8 uM/min). If kcat were ~1 s^-1 = 60 min^-1 (ASSUMED, uncited; measured in vitro rates differ and assays at 10 uM dUMP are not Vmax), [TS] ~ 0.03 uM = ~30 nM. For normal mucosa 0.9 pmol/min/mg -> ~3 nM. These are consistent in order of magnitude with the FdUMP-binding conversions in section 1 (1-35 nM). Because assay conditions and kcat are uncertain, treat as +/- 10x.
4. Mouse L5178Y DHFR: if soluble protein ~100 mg/mL (ASSUMED) and DHFR 10.5% with MW ~21.5 kDa (known approximate), DHFR ~10 mg/mL -> ~0.5 mM. This is an extreme amplified, non-human upper bound.
5. Per-cell copies = uM x 602 molecules/um^3 x cell volume (um^3); e.g. 0.03 uM x 602 x 2000 = ~36,000 copies.

### Gaps
- No human DHFR or SHMT specific activity (nmol/min/mg or U/mg) in normal or tumor tissue/cells found; the DHFR kit literature is rat tissue only.
- No cell-line TS/DHFR activity table retrieved.
- kcat, cell volume and protein density values used above are my assumptions, not from retrieved sources.

## 4. DHFR amplification fold in methotrexate-resistant cells

### Takeaway
Human MTX-resistant lines show ~10-25x DHFR activity increases with ~18 gene copies in CCRF-CEM; mouse lines can reach 160-300x amplification; some resistant lines (HL-60) raise activity without more protein.

### Cited Findings
- Human CCRF-CEM/R1 (T-leukemia), stepwise to 0.2 uM MTX: ~75-fold resistant, ~20-fold DHFR activity, ~18 DHFR gene copies per cell (dot-blot), elongated marker chromosome with homogeneously staining region — [Rutgers record "molecular and karyological analysis of methotrexate resistant and..."](https://www.researchwithrutgers.org/en/publications/molecular-and-karyological-analysis-of-methotrexate-resistant-and/); [BVS record mdl-3855283](https://pesquisa.bvsalud.org/portal/resource/en/mdl-3855283)
- CEM/R3 (2.0 uM MTX) and R1: DHFR activity 10-20x parent — [Cellosaurus publication CLPUB00028](https://cellosaurus.org/cellopub/CLPUB00028)
- CEM subline whose resistance came only from ~20x DHFR increase: MTX EC50 915 nM vs 12.7 nM parental (~72x) — [PMC3850758 Table 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC3850758/table/T2)
- Human colon HCT-8R4: 25-fold increase in DHFR (activity/protein) with more gene copies and transcripts plus Phe-to-Ser mutation raising MTX Kd ~8x — [PubAnnotation PMID 2914962](https://pubannotation.org/docs/sourcedb/PubMed/sourceid/2914962)
- HL-60 MTX-resistant: 20x DHFR activity without matching rise in DHFR protein (more active form) — [PubAnnotation PMID 3856428](https://pubannotation.org/docs/sourcedb/PubMed/sourceid/3856428) (attribution of this finding to this ID is from a search summary; verify)
- CEM/MTX (15.8x resistant, folinic acid): DHFR amplification plus reduced MTX uptake — [UTMB record](https://researchexperts.utmb.edu/en/publications/collateral-sensitivity-to-azidothymidine-in-methotrexate-resistan/)
- HeLa BU-25 parental: low DHFR activity vs very high in MTX-resistant amplified variant (no fold given) — [PMID 1529671](https://www.ncbi.nlm.nih.gov/pubmed/1529671)
- Mouse (non-human): L5178Y 300x DHFR, 10-11% of soluble protein; MELC clone 160x gene copy; carrot cells ~10x — summaries via pubannotation/other search results.

### Inferences
- Human lines typically show 10-25x activity from modest (~18 copy) amplification; extreme 100-300x amplification is mouse data.
- Combined with section 3 conversion, a 20x increase on a baseline of unknown DHFR uM cannot be put in absolute uM without a parental baseline.

### Gaps
- No MCF-7, KB or other epithelial MTX-resistant line fold data retrieved.
- No baseline absolute DHFR copy number in parental lines, so resistant-cell absolute abundance cannot be given.
- Several items here come from abstract-level search summaries; PMIDs for CEM entries were not retrieved.
