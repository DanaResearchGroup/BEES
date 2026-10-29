# Glucose and upper-glycolysis metabolite concentrations: healthy vs tumor tissue and tumor microenvironment

NOTE ON COVERAGE: In this session, WebFetch could not reach elifesciences.org, aacrjournals.org or pmc.ncbi.nlm.nih.gov (DNS/network allowlist), and the PubMed MCP returns abstracts and PMC full text (text only, no figure data or supplementary tables). So most absolute numbers (mM, nmol/g) could NOT be retrieved. Only numbers visible in abstracts or text are reported below. Everything else is listed as a gap. To close the gaps, those hosts need to be allowed, or someone needs to read the figures and tables directly.

## 1. Interstitial/extracellular glucose in tumors vs normal tissue and plasma (mM)

### Takeaway
Qualitatively, tumor glucose is lower than plasma in several tumor types (Hirayama 2009 human colon and stomach tissue; Sullivan 2019 mouse PDAC tumor interstitial fluid), while normal skin interstitial glucose is about the same as plasma. I could not retrieve the actual mM values from any of these papers.

### Cited Findings
- Hirayama 2009 measured 94 metabolites in colon and 95 in stomach, in tumor and grossly normal tissue from 16 colon and 12 stomach cancer patients. The abstract reports "extremely low glucose and high lactate and glycolytic intermediate concentrations" in both tumor types, and calls the tumor microenvironment "unexpectedly poor" nutritionally. The abstract gives no numbers. — [PubMed 19458066](https://pubmed.ncbi.nlm.nih.gov/19458066/), DOI 10.1158/0008-5472.CAN-08-4806
- Sullivan et al. 2019 (eLife 8:e44235) measured about 118-131 metabolites by LC/MS in tumor interstitial fluid (TIF) and plasma from mouse PDAC and LUAD models. The full text states that "glucose levels in the TIF of some tumors are lower than in the circulation, while lactate levels are higher" (citing earlier work), and that "the IF glucose concentration in healthy skin is very similar to that of circulating plasma". Glucose and lactate are quantified in PDAC TIF vs plasma (paired mice) and in Keap1-wild-type vs Keap1-null LUAD allograft TIF, but the values are in figures only and were not in the text I retrieved. — [PMC6510537](https://pmc.ncbi.nlm.nih.gov/articles/PMC6510537), DOI 10.7554/eLife.44235
- Chang et al. 2015 (Cell 162:1229) report that, in a mouse sarcoma model, tumor glucose consumption restricts T cells, and that anti-CTLA-4, anti-PD-1 and anti-PD-L1 antibodies restore glucose in the tumor microenvironment. The abstract gives no mM values. — [PubMed 26321679](https://pubmed.ncbi.nlm.nih.gov/26321679/), DOI 10.1016/j.cell.2015.08.016
- Ho et al. 2015 (Cell 162:1217, phosphoenolpyruvate as a T-cell checkpoint): I found only a secondary commentary, which does not give tumor glucose mM. — [PMC6327313](https://pmc.ncbi.nlm.nih.gov/articles/PMC6327313) (commentary only)

### Inferences
- Normal plasma glucose (~5 mM fasting) is the reference. I did not retrieve a citable source for that number in this session.
- Direction of effect: tumor glucose < plasma in colon, stomach and PDAC. The magnitude is unknown from my sources.

### Gaps
- Exact mM glucose in tumor interstitial fluid (Sullivan 2019 Fig. 1/2 data), in Hirayama tissue (nmol/g), and in Chang 2015 or Ho 2015: NOT FOUND. These are in figures, tables or supplements that I could not access.
- No measurement of normal-organ interstitial glucose (non-skin) was found.
- I did not find early direct tumor-glucose measurements (for example Gullino-era work). A PubMed search for in vivo extracellular tumor glucose returned nothing usable.

## 2. Intracellular/tissue G6P, F6P, FBP, ATP, ADP, AMP, lactate, glucose: normal vs tumor

### Takeaway
The only human tissue dataset identified is Hirayama 2009 (CE-TOFMS, nmol/g). Its abstract confirms elevated lactate and glycolytic intermediates and very low glucose in tumors, but the per-metabolite concentrations were not retrievable. I found no mammalian cell-level mM values.

### Cited Findings
- Hirayama 2009: tumor tissue had extremely low glucose and high lactate and glycolytic intermediate concentrations relative to paired normal tissue, in both colon and stomach. The abstract also reports accumulation of all amino acids except glutamine in tumors. — [PubMed 19458066](https://pubmed.ncbi.nlm.nih.gov/19458066/)
- Park, Rabinowitz et al. 2016 (Nat Chem Biol), "Metabolite concentrations, fluxes and free energies imply efficient enzyme usage", is PMID 27159581. I located it in PubMed but did not read it. It is, to my recollection, an E. coli study and so is not a mammalian source. — [PubMed 27159581](https://pubmed.ncbi.nlm.nih.gov/27159581/)

### Inferences
- Hirayama reports tissue content (nmol/g), not concentrations. Converting to mM would need an assumed tissue water or cell volume (~1 g tissue is roughly 1 mL, with an intracellular fraction of ~0.5-0.8). I did not apply this, because no nmol/g values were retrieved.
- The "low glucose and high glycolytic intermediates" pattern fits the idea that tissue glucose is limiting in tumors even while upstream intermediates accumulate. This is interpretation only.

### Gaps
- G6P, F6P, FBP, ATP, ADP, AMP, lactate and intracellular glucose values for normal vs tumor tissue: NOT FOUND.
- Cell-line intracellular concentrations (for example Bennett 2009, which is E. coli, or mammalian work such as Park 2016): NOT FOUND for mammalian cells.
- HMDB normal reference values: not reached in this session.

## 3. 2-deoxyglucose (2-DG): in vivo/clinical concentrations, hexokinase Ki, Km

### Takeaway
Clinical dosing is documented (Raez 2013: 63 mg/kg/day recommended). Plasma 2-DG concentrations were not available in the abstract. Hexokinase data found: 2-DG Ki 0.074 and 0.15 mM against human brain mitochondrial and soluble hexokinase (competitive with glucose). I found no HK1 vs HK2 isoform-specific Km for glucose.

### Cited Findings
- Phase I trial (n=34), oral 2-DG daily for 7 days every other week, starting at 2 mg/kg, with docetaxel 30 mg/m2. 63 mg/kg was selected as the clinically tolerable dose. At 63-88 mg/kg the effects were reversible hyperglycemia (100%), GI bleeding (6%) and reversible grade 3 QTc prolongation (22%). There was no PK interaction between 2-DG and docetaxel. The abstract gives no plasma concentrations. — [PubMed 23228990](https://pubmed.ncbi.nlm.nih.gov/23228990/), DOI 10.1007/s00280-012-2045-1
- Human brain hexokinase: 2-DG is a competitive inhibitor with respect to glucose, with Ki = 0.074 mM (mitochondrial form) and 0.15 mM (soluble form). Glucosamine Ki was 0.58 and 0.29 mM. Isoform not specified; human brain is mainly HK1. — [PubMed 6693884](https://pubmed.ncbi.nlm.nih.gov/6693884/), DOI 10.1111/j.1471-4159.1984.tb02707.x
- Rat brain hexokinase, reverse reaction: apparent Km for 2-DG-6-phosphate is 2 mM, versus Km and Ki of about 0.2 mM for G6P (allosteric inhibition by G6P). Km for MgADP is 16 mM. — [PubMed 4004267](https://pubmed.ncbi.nlm.nih.gov/4004267/), DOI 10.1016/0003-9861(85)90693-9
- Vendor page (not primary; unverified): 2-DG-6-phosphate Ki = 1.45 mM against bovine heart hexokinase (ATP-competitive). — [MedChemExpress](https://www.medchemexpress.com/2-deoxy-d-glucose-6-phosphate.html)
- In rat macrophages, free 2-dGlc accumulating in the cytosol saturates hexokinase. — [PubMed 8388727](https://pubmed.ncbi.nlm.nih.gov/8388727/)

### Inferences
- A Ki of about 0.07-0.15 mM for 2-DG means that, against a glucose Km of order 0.05-0.1 mM for HK1/HK2 (typical literature value, which I did not verify here), 2-DG would be a potent competitor at micromolar to tens of micromolar concentrations. Intracellular glucose competes, and tumor intracellular glucose is reported to be low, which favors inhibition. This is inference only.

### Gaps
- Plasma 2-DG Cmax/AUC in humans: NOT FOUND in the abstract. The full text of Raez 2013 would have it.
- Km of glucose for HK1 vs HK2 (human, recombinant): NOT FOUND. My PubMed query returned nothing.
- 2-DG Ki specific to HK2: NOT FOUND.
- Intracellular 2-DG-6-phosphate concentrations in tumors: NOT FOUND.
- Other clinical 2-DG trials (for example Singh 2005 radiotherapy trials in glioma): not searched.
