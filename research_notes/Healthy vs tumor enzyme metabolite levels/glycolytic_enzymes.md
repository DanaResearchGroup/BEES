# Upper glycolysis enzymes (HK1/HK2, GPI, PFKM/PFKL/PFKP): abundance and activity, healthy vs tumor

Scope note: this session could reach only WebSearch, the PubMed MCP (abstracts) and search snippets. Direct fetches of pubmed.ncbi.nlm.nih.gov, pmc.ncbi.nlm.nih.gov, proteinatlas.org and spandidos-publications.com failed with DNS errors (hosts not in the session's allowed network domains). Supplementary tables (proteomic ruler, CPTAC, HPA, PaxDb) and full-text tables were therefore NOT readable. Only numbers visible in abstracts or snippets are reported below. Anything else is listed as a gap. Do not treat the gaps as "no data exist".

## Absolute concentrations (uM) or copy numbers per cell of HK, GPI, PFK in normal vs tumor

### Takeaway
No per-enzyme copy number or uM value for HK1/HK2, GPI or PFKM/PFKL/PFKP was retrievable in this session. The datasets that contain them (Nagaraj 2011, Wisniewski 2014, Bekker-Jensen HeLa, Albe 1990 tables) were identified, but their numeric tables were behind blocked hosts.

### Cited Findings
- Proteomic ruler method (Wisniewski, Hein, Cox, Mann 2014, MCP; PMID 25225357): uses histone MS signal as a proxy for DNA, hence cell number, to give copies per cell and concentration without spike-ins. It agreed well with SILAC-PrEST precision values and was applied to mouse organ proteomes. The abstract gives no values for glycolytic enzymes. — [DOI 10.1074/mcp.M113.037309](https://doi.org/10.1074/mcp.M113.037309), [PubMed 25225357](https://pubmed.ncbi.nlm.nih.gov/25225357/)
- Nagaraj et al. 2011 (Mol Syst Biol, HeLa; DOI 10.1038/msb.2011.81): 10,255 proteins from 9,207 genes. More than 90% of abundances lie within a factor of 60 of the median. A secondary summary (BioNumbers BNID 108425) lists a HeLa median of about 18,000 copies/cell. This is a proteome-wide median, not a glycolytic enzyme value. — [Search result for the paper](https://research.uni-luebeck.de/en/publications/deep-proteome-and-transcriptome-mapping-of-a-human-cancer-cell-li/), [BioNumbers](https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=0&id=108425)
- HeLa total protein per cell varies widely between studies: roughly 200-300 pg/cell is reported, up to a 10-fold spread across studies. HeLa stocks from different labs also differ in copy numbers. Treat any HeLa copy number as uncertain by several-fold. — [MDPI 2023 review](https://www.mdpi.com/1422-0067/24/3/2081), [PMC8293326](https://pmc.ncbi.nlm.nih.gov/articles/PMC8293326)
- Albe, Butler, Wright 1990 (J Theor Biol 143:163-195; PMID 2200929) calculated enzyme site concentrations from specific activities and molecular weights across organisms (including human data). Of 140 metabolite:enzyme ratios, 88% were >=1. Of the 17 cases where enzyme exceeds metabolite, 16 were in glycolysis, so glycolytic enzymes are present at concentrations comparable to or above their substrates. Individual HK/GPI/PFK uM values are in its tables, which I could not read. — [DOI 10.1016/s0022-5193(05)80266-8](https://doi.org/10.1016/s0022-5193(05)80266-8)
- A review snippet gives an older PFK estimate of about 5 x 10^4 molecules per cell, with no organism or tissue attribution in the snippet. Low confidence. — [Springer review, "Glycolytic Enzymes as a Multi-Enzyme System"](https://springerlink.bibliotecabuap.elogim.com/content/pdf/10.1007/BF01648926.pdf)
- BioNumbers (citing Albe 1990 tables, rabbit muscle): hexokinase kcat about 12,000 /min (200 /s); phosphoglucose isomerase kcat about 125,400 /min (about 2,090 /s). These are kcat values, not concentrations. Rabbit, not human. — [HK](https://bionumbers.hms.harvard.edu/bionumber.aspx?id=100963), [PGI](https://bionumbers.hms.harvard.edu/bionumber.aspx?id=100993)

### Inferences
- Conversion recipe (my arithmetic, assumptions stated, not from a source). Copies/cell to uM: conc (uM) = copies / (6.022e23 x V_cell) x 1e6. With an assumed HeLa-like cell volume of 2 pL, 1 x 10^6 copies = 1.66e-18 mol / 2e-12 L = about 0.83 uM. So 10^5-10^7 copies/cell corresponds to about 0.08-8 uM. The factor of 2-3 uncertainty in cell volume and the roughly 10-fold uncertainty in total protein per cell propagate directly.
- Activity to uM: [E] (uM) = V (umol/min/mg protein) x protein concentration (mg/mL, assume about 150-200 mg/mL of cytosol-equivalent cell protein, an assumption) / kcat (/min) x 1000 / (units conversion). Example for HK using the rabbit muscle kcat of 12,000 /min: 1 nmol/min/mg protein at 150 mg/mL gives 150 nmol/min/mL = 150 uM-substrate/min, divided by 12,000 /min = about 0.0125 uM (12 nM). Human tumor HK kcat was not found, so this is illustrative only. Real HK activities in cell lines (see below) of 1-20 nmol/min/mg would then correspond to roughly 0.01-0.25 uM. Check the arithmetic before using, since the kcat basis is a different species and tissue.
- Glycolytic enzymes being the main cases where enzyme exceeds substrate (Albe 1990) suggests uM-range enzyme concentrations are physiologically typical, not nM.

### Gaps
- Per-cell copy numbers for HK1, HK2, GPI, PFKL, PFKM, PFKP in HeLa, other cell lines (Geiger 2012, Nagaraj 2011, Wisniewski 2014, Bekker-Jensen 2017 / Cell Systems) were not obtained: supplementary tables not accessible.
- CPTAC, PaxDb and HPA absolute values (ppm or normalized units) not obtained. HPA and PaxDb hosts were blocked or not surfaced.
- Absolute normal-tissue concentrations in uM for these enzymes (Albe 1990 tables) not obtained.
- A human-specific kcat for HK1/HK2, GPI and PFK was not found.

## Tumor/normal fold changes and activity values per tissue

### Takeaway
Classical literature shows higher HK, PFK and PK activity in many malignancies, but the direction is tissue-dependent. Brain tumors (gliomas) show LOWER total HK than normal brain (2-7 fold lower), because normal brain is already an extremely HK-rich tissue. Quantitative per-tissue U/g values from the Weber/Pedersen/Warburg-era tables were not obtainable in this session.

### Cited Findings
- Human normal vs malignant tissues (gastric mucosa vs stomach cancer, and other tissues): hexokinase activity was "significantly higher" in tumors (no numeric ratio in the abstract). The tumor isoenzyme pattern is simplified by loss of the slowest (type IV-like, glucokinase-like) isoenzyme. Km for glucose: tumor isoenzymes about 10^-6 M, normal counterparts about 10^-5 M (about 10-fold higher affinity in tumor). — Monakhov et al. 1978, J Natl Cancer Inst 61:27-34, PMID 276635, [DOI 10.1093/jnci/61.1.27](https://doi.org/10.1093/jnci/61.1.27)
- Human PFK (isozymes M, L, P) in human leukemias, lymphomas and transformed cell lines vs normal counterparts: PFK shows quantitative increases in total activity and isozyme changes. The regulated liver-type (L) isozyme persists and increases, followed by an increase in platelet-type (P). Total PFK activity tracks replication rate. No numbers in the abstract. — Vora, Halper, Knowles 1985, Cancer Res 45:2993-3001, PMID 3159473 (no DOI listed in PubMed record); [PubMed record via search tool]
- Maximal activities of key enzymes in homogenates of normal, neoplastic and suppressed cells (human lymphocytes, macrophages and tumor lines among them): HK activity is consistently lower than PFK activity in tumor cells, so HK may approach saturation by intracellular glucose and HK plus PFK can co-control flux. Pyruvate kinase: non-tumorigenic 28.4-574, tumorigenic 899-1280, metastatic 1590-1627 nmol/min/mg protein. PK : 2xPFK ratio is 22.4 (neoplastic mean) vs 3.8 (muscle from 60 animals). HK and PFK individual values are in the paper's tables (not readable here). — Board, Humm, Newsholme 1990, Biochem J 265:503-509, PMID 2302181, [DOI 10.1042/bj2650503](https://doi.org/10.1042/bj2650503); review version Newsholme & Board 1991, PMID 1877389, [DOI 10.1016/0065-2571(91)90015-e](https://doi.org/10.1016/0065-2571(91)90015-e)
- Glucokinase-like (high Km 8-11.8 mM) glucose-phosphorylating activity in 8 human tumor cell lines, 1 human tumor tissue and 4 rat tumor lines: 1.1-20 nmol/min/mg cell protein, consistently absent from untransformed cell lines/tissues except rat liver (18 nmol/min/mg). — Board, Colquhoun, Newsholme 1995, Cancer Res 55:3278-3285, PMID 7614462
- HT29 human colon carcinoma: HK, PFK-1 and PK activities were about 2-fold higher in poorly differentiated, highly glycolytic Glc+ cells than in differentiated Glc- cells. — Gauthier et al. 1989, Biochim Biophys Acta 975:231-238, PMID 2525930, [DOI 10.1016/s0005-2728(89)80253-1](https://doi.org/10.1016/s0005-2728(89)80253-1)
- Gliomas: total HK activity 2-4 fold lower than normal rat brain in 4 human xenografted gliomas; mitochondria-bound HK 9-69% of the total; lactate/pyruvate 3-4 fold higher. — Oudard et al. 1997, Anticancer Res 17:1903-1911, PMID 9216643
- Gliomas (human): mean HK activity about 3-fold lower than normal tissue and mostly mitochondria-bound (Oudard 1996, Oncol Rep 3:165, via search snippet; the same group reports a mean 7-fold lower HK in individual glioblastomas in Br J Cancer, via snippet). The Oudard Oncol Rep paper could not be opened (host blocked); exact U values are not available. — [Oncol Rep snippet](https://spandidos-publications.com/10.3892/or.3.1.165)
- Cultured glioma lines: HK and PFK not different between low-grade glioma and normal astrocytes, but consistently and significantly increased in high-grade glioma lines (search snippet; magnitude not given). — [Oncol Rep snippet](https://spandidos-publications.com/10.3892/or.3.1.165)
- Modern expression-level data: a pan-cancer TCGA/TIMER2 analysis reports PFKP higher in most tumor types, with a CPTAC lung adenocarcinoma protein z-score showing higher PFKP in tumor; fold changes not read (figure only). — [Cancer Cell Int 2024 / pan-cancer PFKP](https://link.springer.com/article/10.1186/s12935-024-03497-w), [Sci Rep 2023](https://www.nature.com/articles/s41598-023-43982-2.pdf)
- HK2 pan-cancer TIMER2 (33 tumor types): tumor > normal in CHOL, HNSC, KIRC, KIRP, LUSC, PRAD, STAD, THCA, UCEC (mRNA, no fold values in snippet). — [Sci Rep 2022](https://www.nature.com/articles/s41598-022-23598-8.pdf)
- A ccRCC study notes that in tumors and cancer cell lines, PFK-P or PFK-L, or both, are the most abundant isoforms (qualitative). — [Cell Oncol 2020](https://springerlink.bibliotecabuap.elogim.com/10.1007/s13402-020-00508-6)

### Inferences
- A consistent picture: epithelial and hematologic tumors tend to have about 1.5-3x higher HK/PFK activity than their normal counterpart or less-glycolytic sublines (HT29: about 2x), but this is not universal. Glioma shows the opposite for total HK, because the comparator (brain) has high baseline HK.
- Tumor HK affinity for glucose shifts (Km 10^-5 to 10^-6 M in human stomach), consistent with the loss of low-affinity isoforms and gain of HK2/HK1-type.
- The Weber, Pedersen and Warburg-era U/g tissue tables were not obtained; the correct way to fill this is from those primary papers, not from this note.

### Gaps
- Classic Weber (Weber G, et al.; "Biochemical strategy of cancer cells and the design of chemotherapy", Cancer Res 1983; NEJM 296:486, 1977) hexokinase/PFK activities in human tumors in umol/min/g tissue: not found in an accessible form. My query "hexokinase activity human tumors Weber" returned only an unrelated (Wolfgang Weber) PET paper.
- Numeric tumor/normal fold change per tissue for HK1, HK2, GPI, PFKL, PFKM, PFKP from CPTAC proteomics: not obtained (UALCAN/LinkedOmics not reachable).
- GPI activity in tumor vs normal tissue: a PubMed search returned 61 hits (e.g. PMIDs 38950696, 36915794, 36401990, 36181462, 35820598, 35718053) but I did not open them; no activity value is claimed here.
- Specific activity values (U/g, umol/min/g) for human tumor vs normal tissue pairs: not retrieved.

## HK2 specifics (mitochondrial binding, overexpression in tumors)

### Takeaway
HK2 is highly expressed in many tumors but in few normal adult tissues, and it binds the outer mitochondrial membrane through VDAC via its N-terminus, where it supports glycolysis and suppresses apoptosis. Quantitative bound/free fractions in human tumors were found only for gliomas (9-69% bound).

### Cited Findings
- HK2 is overexpressed in cancer cells and bound to VDAC on the outer mitochondrial membrane, supporting the glycolytic phenotype and suppressing mitochondrial apoptosis. — Mathupala, Ko, Pedersen 2006, Oncogene, [PMC3385868](https://pmc.ncbi.nlm.nih.gov/articles/PMC3385868)
- The mitochondrial binding determinant lies in the first ten N-terminal residues of HK2. — [PMC4613670](https://pmc.ncbi.nlm.nih.gov/articles/PMC4613670)
- HK2 expression in normal adult tissue is mainly skeletal muscle, heart and adipose tissue (insulin-sensitive). HK2 is highly expressed in many cancer cell types. — [Sci Rep 2022 pan-cancer](https://www.nature.com/articles/s41598-022-23598-8.pdf); a Japanese review snippet also states that type II is only slightly expressed in normal heart, muscle and adipose and is remarkably expressed in malignant tumor cells ([Yakugaku Zasshi](https://www.jstage.jst.go.jp/browse/yakushi/120/8/_contents/-char/ja)).
- Mitochondria-bound HK in HT29 colon cancer cells: a large fraction of HK is mitochondrially bound, and binding efficiency of ATP channeling differs by differentiation state (exogenous ATP 6.5x more efficient than endogenous in Glc+, only 1.8x in Glc-). — PMID 2525930 (above)
- Mitochondria-bound HK fraction in xenografted human gliomas: 9-69% of total HK. — PMID 9216643 (above)

### Inferences
- HK2's tumor relevance comes more from relocalization and isoform switching than from total activity. In glioma, total HK drops yet mitochondrial binding rises.

### Gaps
- HK2 tumor/normal protein fold changes by tissue and copy number per cell: not obtained.
- HK1 vs HK2 split in specific human tumor types in absolute units: not obtained.
- Percentage of HK2 mitochondria-bound in non-glioma human tumors: not obtained.
