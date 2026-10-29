# Folate-cycle and nucleotide metabolite concentrations: healthy vs tumor

Access note: in this session WebFetch could not reach nature.com, ncbi.nlm.nih.gov, pubmed.ncbi.nlm.nih.gov or bionumbers.hms.harvard.edu (DNS/network allow-list). Only the PubMed MCP tool (full text of PMC articles, metadata) and WebSearch summaries worked. Supplementary tables (where most absolute concentrations live) could not be opened. Many requested numbers are therefore reported as NOT FOUND. Items marked "[search-summary only]" come from WebSearch snippets that I could not open and verify, so treat them as provisional.

## What is the typical total intracellular folate concentration in normal vs tumor cells/tissue?

### Takeaway
I did not find a verified total-folate pool in mM (or a clean paired normal vs tumor total-folate comparison). The only quantitative tumor-vs-normal tissue data I could verify are colorectal LC-MS/MS values in pmol/g wet weight (low nmol/kg range, roughly 1-4 uM tissue-equivalent after leucovorin; baseline lower). These are from patients, many given leucovorin, so they are not true baseline.

### Cited Findings
- Colorectal cancer, 76 analysed patients (80 enrolled), LC-MS/MS of THF, 5,10-methylene-THF and 5-methyl-THF in tumour and macroscopically normal mucosa 10 cm from tumour, reported in pmol/g wet weight. Mean 5,10-methyleneTHF in tumour for the four groups (0, 60, 200, 500 mg/m2 leucovorin IV at surgery) was 880, 1,769, 3,024 and 3,723 pmol/g. Highest control (no leucovorin) tumour methyleneTHF value was 1,714 pmol/g. — [Taflin et al. 2014, Cancer Chemother Pharmacol, PMID 25238909, DOI 10.1007/s00280-014-2591-9](https://doi.org/10.1007/s00280-014-2591-9) (full text read via PMC4236605)
- Same study: in untreated controls, methyleneTHF and THF were significantly higher in tumour than in normal-appearing mucosa; mucosa folates were generally lower in the control group. 5-methylTHF did not differ significantly between tumour and mucosa in any leucovorin group, but in the control group mucosa 5-methylTHF was significantly higher than tumour. Rectal tumours had significantly lower methyleneTHF than colon tumours. Large inter-patient variability. Exact table values for THF and 5-methylTHF were not in the extracted text. — [Taflin 2014](https://doi.org/10.1007/s00280-014-2591-9)
- Method caveat from that paper: the method deconjugates polyglutamates, so values are monoglutamate-equivalent totals per form, not a polyglutamate distribution. Tissue homogenised in 10x volume buffer with ascorbate; recovery 98%, 87%, 93% for THF, methyleneTHF, methylTHF. — [Taflin 2014](https://doi.org/10.1007/s00280-014-2591-9)
- Cultured human cells: cytoplasmic 5-methyl-THF accumulation was markedly inhibited once intracellular 5-methyl-THF reached about 5-7 pmol/10^6 cells. In normal human fibroblasts most 5-methyl-THF was polyglutamated. [search-summary only; JCI/PMC links not opened] — [JCI 114310](https://jci.org/articles/view/114310); [PMC372044](https://pmc.ncbi.nlm.nih.gov/articles/PMC372044/?page=-1)
- In 14 human cancer cell lines, basal intracellular folate correlated positively with FPGS activity (numbers not retrieved). [search-summary only] — [BJC 2013 pdf](https://www.nature.com/articles/bjc2013579.pdf)

### Inferences
- Conversion (my assumption, not from a source): tissue density ~1 g/mL and ~70-80% water, so 1,000-4,000 pmol/g ≈ 1-4 nmol/mL tissue ≈ roughly 1.3-5 uM per tissue water, and higher per cell water if only the cellular compartment is considered. 5,10-methyleneTHF alone is therefore low-uM in colorectal tumour tissue after leucovorin, and below ~1.7 uM-equivalent without it.
- 5 pmol/10^6 cells with a typical ~2 pL cell volume (assumption: 2e-12 L, HeLa-like) ≈ 2.5 uM; with 4 pL ≈ 1.2 uM. So total cellular folate of order 1-10 uM is the plausible magnitude, but this rests on my volume assumption and an unverified cell-culture number.

### Gaps
- No verified total folate (sum of all forms plus polyglutamates) in normal vs tumor tissue in a single paper; no polyglutamate chain-length distribution values retrieved (likely in literature on FPGS and in the Liu 2012 Cancer Prev Res colorectal mucosa study, not opened).
- Verified numbers for DHF, and a normal-vs-tumour 5-methyl-THF value in pmol/g, not obtained (Taflin tables not extracted).
- Cheradame 1997 (Eur J Cancer, tumour folate) and Nature Metabolism 2021 (A549 LC-MS folate species per cell) were suggested by search results but I could not open them; no PMIDs verified.

## What are the dUMP, dTMP, serine, glycine, NADPH, NADP+ concentrations in cells, normal vs cancer?

### Takeaway
I could not retrieve any verified absolute concentrations for serine, glycine, NADP+, dTMP or the NADPH/NADP+ ratio; they sit in supplementary tables I could not access. Provisional values: cytosolic free NADPH about 3 uM and mitochondrial about 37 uM (sensor-based, HeLa); baseline dUMP roughly 2-5 pmol/10^6 cells in tumour cell lines; dTTP about 20 pmol/10^6 cells in untreated HCT116. No direct normal-vs-cancer comparison was found for any of these.

### Cited Findings
- Free (not total) NADPH measured with genetically encoded sensors: cytosol 3.1 +/- 0.3 uM, mitochondrial matrix about 37 +/- 2 uM; work in HeLa cells (Tao et al. 2017, Nature Methods, nmeth.4306; PMID and exact DOI not verified by me). [search-summary only] — [Nature Methods pdf](https://www.nature.com/articles/nmeth.4306.pdf)
- Park/Rabinowitz iBMK data (mouse, immortalized, not human): NADPH, ATP, NAD+ almost always saturating at enzyme sites; NADP+ not reliably saturating, which the authors interpret as prioritising a high NADPH/NADP+ ratio. Sum of measured metabolites about 200 mM in all organisms; glutamate the most abundant metabolite; amino acids the largest class. Numeric NADPH, NADP+, serine, glycine values are in the supplementary compilation, which I could not open. — [Park et al. 2016 Nat Chem Biol, PMID 27159581, DOI 10.1038/nchembio.2077](https://doi.org/10.1038/nchembio.2077) (full text PMC4912430)
- Conversion assumptions in that paper: for iBMK cells, 70% of packed cell volume taken as aqueous volume; NADH and NADPH and nucleotides were extracted with acidic solvent (0.1 M formic acid) then neutralised. Inorganic phosphate input 5 mM for iBMK. — [Park 2016](https://doi.org/10.1038/nchembio.2077)
- Baseline dUMP: older study reported dUMP of 2-5 nmol per 10^9 cells (= 2-5 pmol/10^6 cells) in logarithmically growing tumour cell lines, rising sharply after 5-FUra. [search-summary only; original (PMC348628 or PMC383271) not opened] — [PMC348628](https://pmc.ncbi.nlm.nih.gov/articles/PMC348628); [PMC383271](https://pmc.ncbi.nlm.nih.gov/articles/PMC383271)
- dTTP: HCT116 colon cancer cells, vehicle-treated about 20.4 pmol dTTP per 10^6 cells; 1 uM FdUrd reduced it to about 4.5. This comes from a patent text, not peer-reviewed. [search-summary only] — [US patent 8815509](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/8815509)
- FdUrd-treated lymphoblasts reached dUMP about 1000 pmol/10^6 cells (drug-perturbed, not baseline). [search-summary only] — [Mol Pharmacol 21:211](https://molpharm.aspetjournals.org/content/molpharm/21/1/211.full.pdf)
- LC-MS/MS methods for dUMP/TMP in cultured cells exist (J Chromatogr B 2009; Biochem Pharmacol 2011 reports dUMP about 2-fold higher in hypoxic PANC1). [search-summary only] — [catalog record](https://katalog.lib.cas.cz/KNAV/EdsRecord/edsair,edsair.doi.dedup.....0a22ad7f5e8e6abe474946cb126cf42e)

### Inferences
- Conversion (assumption): 2-5 pmol/10^6 cells at ~2 pL/cell ≈ 1-2.5 uM dUMP; 20 pmol/10^6 cells dTTP at ~2 pL ≈ 10 uM (about 5 uM at 4 pL). These are rough and depend on cell volume.
- The sensor-based free NADPH (uM) is far below total NADPH from extraction methods (typically hundreds of uM to low mM range in literature), because most is protein-bound; the two are not interchangeable. This is general knowledge, not sourced here.
- The Rabinowitz statement that NADP+ is often not saturating, while NADPH is, implies a strongly reduced NADPH/NADP+ pool, but no ratio value was retrieved.

### Gaps
- Serine, glycine, NADP+, dTMP, and NADPH/NADP+ ratio numbers: not found (supplementary tables inaccessible). No verified normal-vs-tumour paired values for any of these. Candidate sources to open once network access allows: Park 2016 Supplementary Table (iBMK compilation; BioNumbers entries 112692, 112694, 112689 summarise it), Bennett 2009 (E. coli), Lu/Rabinowitz-type papers on NADPH/NADP+ in mammalian cells, tumour-vs-normal tissue metabolomics (e.g. Hu 2013 colon / Kami 2013 type studies; not verified here).
- Human-specific absolute data: none verified; Park is mouse iBMK.

## Which absolute-quantification metabolomics datasets exist?

### Takeaway
Key absolute-concentration datasets verified: Bennett 2009 (E. coli, more than 100 metabolites, about 300 mM total) and Park 2016 (E. coli, yeast, mouse iBMK cells, about 200 mM total). Both use 13C-labelled internal standards and LC-MS. No verified human tumour-tissue absolute dataset covering the requested folate/nucleotide set was found.

### Cited Findings
- Bennett et al. 2009: LC-MS/MS quantified more than 100 metabolites in aerobic, exponentially growing E. coli on glucose, glycerol or acetate; total pool about 300 mM; glutamate most abundant; most substrate concentrations exceed Km. — [Nat Chem Biol 5:593, PMID 19561621, DOI 10.1038/nchembio.186](https://doi.org/10.1038/nchembio.186) (PMC2754216)
- A BioNumbers entry (via search summary) gives glucose-fed E. coli glutamate about 96 mM from Bennett Table 1; I could not confirm it in the paper itself. [search-summary only] — [BioNumbers 104691](https://bionumbers.hms.harvard.edu/bionumber.aspx?id=104691)
- Park 2016: absolute concentrations and fluxes in E. coli, yeast and mouse iBMK cells (integrated with thermodynamics); BioNumbers summarises iBMK measured-metabolite total as about 176 mM vs 238 mM yeast and 243 mM E. coli. [BioNumbers figures search-summary only] — [Park 2016](https://doi.org/10.1038/nchembio.2077); [BioNumbers 112694](https://bionumbers.hms.harvard.edu/bionumber.aspx?id=112694)
- Method for absolute quantitation in Park 2016: grow cells in U-13C glucose/glutamine, extract with known concentrations of authentic standards; amino acids via unlabeled DMEM plus 13C,15N amino acid internal standards; concentration normalised to packed cell volume (iBMK) with 70% aqueous fraction. — [Park 2016](https://doi.org/10.1038/nchembio.2077)

### Inferences
- For converting per-cell data, Park's approach (packed cell volume x 0.7) is a defensible template; for tissue, divide pmol/g by ~0.7-0.8 g water/g to approximate per-water concentration (my assumption).

### Gaps
- HMDB/ECMDB normal-concentration entries were not queried (no accessible fetch); no data on tumour tissue absolute quantification (e.g. LC-MS with isotope standards for NADPH/NADP+, serine, glycine) retrieved.
- Lu/Rabinowitz folate-species work on mammalian cells (e.g. compartmentalised folate measurements) not located in this session.
