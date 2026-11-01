# Literature Verification Report — `feat/drug-delivery-inhibitors`

**Scope:** This report follows up on the Claude Code session that built the `feat/drug-delivery-inhibitors` branch (base commit `ad48a1d`) but was blocked from PubMed, PMC, EuropePMC, EBI, Nature, and Protein Atlas by its network sandbox. All citations and placeholder kinetic values that session flagged `VERIFY`, `unverified`, or left as literature gaps were re-checked here against PubMed, PMC, ChEMBL, and the Human Protein Atlas, which were reachable from this session.

**Repo state note:** the local clone at `/home/omerkfir/BEES` is mounted read-only in this session, so none of the fixes below have been committed. Section "Action items for the next Claude Code session" lists exactly what to change and where.

---

## 1. Pitch-doc follow-up citations (`projects/inhibitor_delivery_pitch.md`)

All four PMIDs cited in the pitch document's "follow-up reading" section were re-fetched directly from PubMed and all check out.

| PMID | Claimed paper | Verified title / journal / authors | Status |
|---|---|---|---|
| 19731036 | Enden & Schroeder 2009, liposome ultrasound release model | "A mathematical model of drug release from liposomes by low frequency ultrasound." Ann Biomed Eng, 2009, Enden & Schroeder | **Confirmed** |
| 31503443 | Zinger et al. 2019, collagozome model | "Collagenase Nanoparticles Enhance the Penetration of Drugs into Pancreatic Tumors." ACS Nano, 2019-09-20, Zinger A. et al. (Technion/Schroeder lab) | **Confirmed** — abstract matches pitch doc's collagen %vol and tumor-shrinkage figures |
| 42312441 | Ackerman et al. 2026, synthetic-cell energy metabolism | "Phosphoserine as an Alternative Energy Source for Cell-Free Protein Synthesis with Increased Yield and Prolonged Activity." ACS Synth Biol, 2026-06-18, Ackerman S. et al. (Schroeder lab) | **Confirmed** — this PMID is legitimately a 2026 publication, not a fabricated/out-of-range number as initially suspected |
| 27830705 | Yaari et al. 2016, barcoded-nanoparticle screens | "Theranostic barcoded nanoparticles for personalized cancer medicine." Nat Commun, 2016-11-10, Yaari Z. et al. (Schroeder lab) | **Confirmed** |

**No changes needed** to the pitch doc's citation list.

---

## 2. Triclosan → FabI Ki (Ward 1999, flagged `VERIFY`)

The showcase `projects/fattyAcidSynthesis/fas_inhibitor_delivery/` uses a tight-binding Ki of **7×10⁻⁹ mM (7 pM)** for triclosan inhibiting FabI, citing "Ward 1999, VERIFY."

**Found:** Ward WH et al., *Biochemistry*, 1999-09-21, PMID 10493822, "Kinetic and structural characteristics of the inhibition of enoyl (acyl carrier protein) reductase by triclosan."

The paper reports multiple related constants for *E. coli* EACPR/FabI, not a single Ki:
- Triclosan binds to E-NAD(+) complex, with a dissociation constant around 20-40 pM.
- Triclosan follows competitive kinetics with respect to NADH, giving an inhibition constant of 38 pM at zero NADH and saturating NAD(+). Uncompetitive kinetics are observed when NAD(+) is varied, giving an inhibition constant of 22 pM at saturating NAD(+).
- The association rate constant (k(on)) is estimated as 2.6 x 10(7) s(-1) M(-1)... the ratio k(off)/k(on) = 7.1 pM is similar to the inhibition constants from the steady-state studies.

**Verdict:** The citation ("Ward 1999") is correct, and the BEES value (7 pM) is well-supported — it matches the koff/kon-derived constant (7.1 pM) almost exactly — but it sits slightly below the two directly-measured steady-state Ki values (22 pM uncompetitive, 38 pM competitive). The "VERIFY" flag can be removed; the number is real and traceable, but the showcase README should state explicitly which of the three reported constants (koff/kon vs. competitive-Ki vs. uncompetitive-Ki) is being used, since they differ by up to ~5x.

---

## 3. Methotrexate → DHFR Ki (ChEMBL median vs. Volpato 2009)

The showcase `projects/methotrexate_dhfr/` uses **Ki = 1×10⁻⁷ mM (0.1 nM)**, described as the median of three "unflagged" ChEMBL activity records for target CHEMBL202 / molecule CHEMBL34259.

**Re-fetched directly from the live ChEMBL REST API:**

| Activity ID | Ki (nM) | Source document | Target species |
|---|---|---|---|
| 3055901 | 0.077 | CHEMBL1152608 — Bioorg Med Chem 2009, PMID 19748785 | Human |
| 939069 | 0.179 | CHEMBL1131724 — Bioorg Med Chem Lett 1999, PMID 10360757 | ***Trypanosoma cruzi* DHFR, not human** |
| 18251252 | 0.01 | CHEMBL4052640 — Eur J Med Chem 2017, PMID 28477572 | Human |

**Problem found:** activity 939069 is methotrexate's Ki against *T. cruzi* DHFR from a paper on selective anti-trypanosomal DHFR inhibitors — not human DHFR. Mixing it into a "human DHFR Ki" median is a species error. Excluding it, the two human-DHFR values are 0.077 and 0.01 nM, median 0.0435 nM — still not matching the 0.1 nM figure used.

**Primary-literature cross-check — Volpato et al. 2009** (PMID 19478082, PMC2740434, DOI 10.1074/jbc.M109.018010, *J Biol Chem*): "Methotrexate is a slow, tight-binding, competitive inhibitor of human dihydrofolate reductase (hDHFR)." The paper's introduction (confirmed via the jbc.org full-text mirror, not retrievable from PMC's XML which this publisher blocks) states the wild-type Ki for methotrexate against human DHFR is **3.4 pM (0.0034 nM)** — roughly 23× tighter than the ChEMBL human-only median and 29× tighter than the value BEES currently uses.

| Source | Ki (mM) | Ki (nM) |
|---|---|---|
| BEES current value | 1e-7 | 0.1 |
| ChEMBL median (human-only, 2 records) | 4.35e-8 | 0.0435 |
| ChEMBL median (all 3 records, species-mixed) | 7.7e-8 | 0.077 |
| Volpato et al. 2009 (primary literature) | 3.4e-9 | 0.0034 |

**Verdict:** two issues, both actionable:
1. Activity 939069 should be dropped from any "human DHFR" ChEMBL query — it's *T. cruzi* DHFR.
2. The literature value (Volpato 2009, 3.4 pM) is an order of magnitude tighter than any ChEMBL assay value used. ChEMBL binding assays can differ from kinetically-determined Ki by this much depending on assay conditions (pre-incubation time matters enormously for tight-binding inhibitors — this is the same phenomenon Ward et al. 1999 flags for triclosan above). Recommend citing Volpato 2009 directly for this showcase's Ki, since the showcase already frames methotrexate as a literal tight-binding inhibitor case and Volpato's value is purpose-measured for exactly that regime.

---

## 4. Smaller unverified claims (`research_notes/`)

| Claim | PMID found | Verdict |
|---|---|---|
| Larsson 1996 TS tumor-vs-normal-liver (previously: no PMID found) | **8695163** — Larsson et al., *Acta Oncologica* 1996, "Thymidylate synthase in advanced gastrointestinal and breast cancers." | **Resolved.** TS levels were 2.4 +/- 0.31 pmol/g in liver metastases of colorectal cancer (n = 87), 4.2 +/- 1.0 pmol/g in primary colorectal cancer (n = 13)... and 0.58 +/- 0.075 pmol/g in normal liver parenchyma (n = 24). TS levels were significantly higher in malignant tumor tissues compared to normal liver parenchyma. Direction is consistent with the Umekita 2000 comparison already in the notes. |
| PMID 3856428 attributed to "20-fold DHFR activity increase, MTX-resistant HL-60" (flagged: verify this is not a hallucinated search-summary attribution) | 3856428 | **Confirmed accurate.** Dedhar, Hartley & Goldie, *Biochem J* 1985: "Methotrexate(MTX)-resistant human promyelocytic-leukaemia cells (HL-60)... have a 20-fold increase in dihydrofolate reductase (DHFR) activity as compared with the sensitive cells. This increase is not associated with a concomitant increase in DHFR protein." The attribution was correct; flag can be removed. |
| Monakhov et al. 1978 hexokinase tumor-vs-normal, PMID 276635 | 276635 | **Confirmed.** "Physicochemical properties and isoenzyme composition of hexokinase from normal and malignant human tissues." *J Natl Cancer Inst* 1978 — matches the described elevated tumor HK activity and isoenzyme-pattern simplification. |
| Volpato et al. 2009 metadata spot-check | 19478082 | **Confirmed.** PMCID PMC2740434, DOI 10.1074/jbc.M109.018010 all resolve correctly and match. |

---

## 5. Proteomic copy-number gap (folate and glycolytic enzymes) — still open

The blocked Claude Code session could not get absolute per-cell copy numbers for DHFR, TYMS, SHMT1, SHMT2, HK1, HK2, GPI, PFKL, PFKM, or PFKP from any proteomic-ruler source (Wisniewski 2014, Nagaraj 2011, CPTAC, PaxDb). This session had broader network access but **this specific gap remains unresolved**:

- **PaxDb** (pax-db.org, v5.pax-db.org): reachable (after allowlisting), but both are JavaScript single-page apps with no discoverable REST/JSON endpoint for per-protein values from this sandbox — data loads client-side from binary/bundled assets this session could not locate a URL for.
- **ProteomicsDB** (proteomicsdb.org): reachable (after allowlisting), but its OData-style API rejected every parameter schema tried (`proteinpertissue.xsodata`, `proteinexpression.xsodata`) with "no key property" errors — the exact query signature wasn't discoverable without its web UI.
- **Human Protein Atlas** (proteinatlas.org API): **reachable without any allowlist request** and returned results for all 10 genes, but only qualitative IHC-reliability and subcellular-localization fields — not nTPM or copy-number columns (those returned `null` for every gene tested, likely because the cell-line-specific endpoint needs a different column key than `rnaclsm`).
- **Nagaraj et al. 2011** (PMID 21988835, *Mol Syst Biol*, HeLa deep proteome, 10255 proteins) and **Wisniewski et al. 2014** (PMID 25225357, *Mol Cell Proteomics*, proteomic-ruler method) full texts were fetched and confirm the methodology exists and per-protein copy numbers were computed in both studies, but the actual per-gene values are in supplementary Excel tables, not in the article body, and those supplementary files were not retrievable through any API tried.

**This gap is not resolved.** It should remain flagged as an open item, not be treated as "checked and absent."

---

## 6. Summary table

See [literature_verification_status.csv](literature_verification_status.csv) for the full machine-readable table. Status counts:

| Status | Count |
|---|---|
| Resolved / confirmed as-is | 8 |
| Resolved with a correction needed | 2 (triclosan Ki — annotate which constant; methotrexate Ki — fix species mix + consider literature value) |
| Still open (proteomic copy numbers) | 2 |

---

## Action items for the next Claude Code session

1. **Triclosan Ki** (`projects/fattyAcidSynthesis/fas_inhibitor_delivery/`): remove the "VERIFY" flag on Ward 1999 — citation is correct (PMID 10493822). Optionally annotate in the README that 7e-9 mM corresponds to the koff/kon-derived constant (7.1 pM) rather than either of the two directly-measured steady-state Ki values (22 pM, 38 pM) reported in the same paper.
2. **Methotrexate Ki** (`projects/methotrexate_dhfr/`):
   - Drop ChEMBL activity 939069 from the "unflagged human DHFR" set — it is methotrexate vs. ***Trypanosoma cruzi*** DHFR (CHEMBL1131724), not human.
   - Consider replacing the current Ki (1e-7 mM / 0.1 nM) with the Volpato et al. 2009 value, **Ki = 3.4e-9 mM (3.4 pM)**, PMID 19478082 — this is a primary-literature kinetically-determined value for the exact tight-binding regime the showcase models, and is ~29x tighter than the current figure. If keeping the ChEMBL-derived value instead, cite it as the median of the two human-only records (0.077 and 0.01 nM, median 0.0435 nM) rather than three.
3. **Research notes**:
   - Add PMID 8695163 (Larsson et al. 1996, *Acta Oncologica*) as the citation for the previously-unsourced TS tumor-vs-normal-liver comparison in `folate_enzymes.md`.
   - Remove the "verify this attribution" flag on PMID 3856428 in `folate_enzymes.md` — the 20-fold DHFR activity increase in MTX-resistant HL-60 cells is accurately attributed.
4. **Proteomic copy-number gap** (DHFR/TYMS/SHMT1/SHMT2/HK1/HK2/GPI/PFKL/PFKM/PFKP): still open. Next attempt should target the supplementary tables directly — e.g. the Nagaraj 2011 Supplementary Table S1 (hosted via the journal's own supplementary-material links, often a separate domain from the main article) or PaxDb's static bulk-download `.txt` files (not its web API) if those specific file paths can be allowlisted.
5. **2-deoxyglucose and thiolactomycin Ki values** remain explicit placeholders per the original brief and were out of scope for this verification pass (not flagged with a specific citation to check); a literature search for real values was not attempted here since no citation was given to verify against — this would need a fresh literature search, not a verification, and should be scoped as new work if desired.
