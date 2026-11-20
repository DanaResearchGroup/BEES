# Rebuild `db/ecoli.csv` (BEES E. coli K-12 reaction database)

`build_ecoli_csv.py` rebuilds `db/ecoli.csv` of BEES v1.0.0 from the UniProt download plus
the curation tables in `data/`. Standard-library Python 3.10+, no internet needed for the build.

Result with the 2026-04-20 UniProt file (checked 2026-11-20, from a clean copy of this folder):
- Byte-identical to the released `db/ecoli.csv`, except for 2 lines (the last two rows) that
  have Windows line endings (CRLF) in the released file and LF in the rebuild. Every row and
  every column is equal; JSON cells are equal as data.

## 1. Download the raw input (not shipped)

Only UniProt is needed for the build. The file must be the **2026-04-20** snapshot, see the warning below.

- Source: UniProt proteome UP000000625 (E. coli K-12), https://www.uniprot.org/proteomes/UP000000625
- Licence: CC BY 4.0; cite The UniProt Consortium.
- Download (the same query, as a gzipped TSV with these 10 columns: Entry, Reviewed, Entry Name,
  Protein names, Gene Names, Organism, Length, Catalytic activity, Sequence, EC number):

```bash
curl -s "https://rest.uniprot.org/uniprotkb/stream?format=tsv&query=proteome:UP000000625&fields=accession,reviewed,id,protein_name,gene_names,organism_name,length,cc_catalytic_activity,sequence,ec" \
  | gzip > uniprot_UP000000625.tsv.gz
```

- SHA-256 of the 2026-04-20 file the release was built from:
  `b6809a45c81de03354856429441a29fbbb6919a24ca52b2e657ed5d63e213d99`
  (check with `sha256sum uniprot_UP000000625.tsv.gz`). A download made today does not give this
  file; see the warning below.

> **Warning.** UniProt changes between releases. A download made today (2026-11-20) has
> 2,417 reaction rows instead of 2,327, and the build then no longer matches `db/ecoli.csv`.
> The hand-curation tables (`dropped_rows.csv`, `row_exceptions.csv`, `enzyme_names.csv`,
> `added_rows.csv`) refer to rows by their position in the 2026-04-20 file, so they do not
> line up with a newer download. Use the 2026-04-20 file to reproduce the release.

## 2. Build and compare

Run from this folder (`db/db_build/ecoli_csv/`):

```bash
python build_ecoli_csv.py --uniprot uniprot_UP000000625.tsv.gz --out ecoli.csv
python build_ecoli_csv.py --uniprot uniprot_UP000000625.tsv.gz --out ecoli.csv --check ../../ecoli.csv
```

The second command prints `Same content as the released file` and `Not byte-identical`
(the CRLF lines above). Use `cmp ecoli.csv ../../ecoli.csv` to see the byte difference.

## What the build does
1. **UniProt** (the TSV from step 1):
   - Gives one row per catalytic-activity reaction. The EC number is the inline `EC=`, else
     the protein's single EC, else ECs written in the protein name; multi-EC reactions get
     one row per EC.
   - Removes EC numbers from the protein name and moves synonyms in parentheses to
     `Protein aliases`.
   - Skips reactions without an EC number or a parseable equation.
2. **ECMDB names**:
   - Compounds that match ECMDB (2015) are renamed to the ECMDB name in the stoichiometry,
     and also inside the reaction text (word-boundary replacement, longest names first).
   - This is why some names are capitalised differently in different reactions.
3. **SMILES**: taken per compound name from `data/compound_smiles.csv`, with the source
   label from `data/smiles_source_labels.csv`.
4. **Hand curation** (recorded in git between 2026-05 and 2026-11; see "Provenance of the tables"):
   - 70 short enzyme names (FabA, FabZ, TesA, ...).
   - 14 dropped rows: 8 FabA dehydration rows, the generic TesA acyl-CoA row, and 5 SurE
     rows under EC 3.6.1.11 (removed in commit `02624a2`, 2026-11-17).
   - 54 added FAS rows: 38 chain-length rows inferred from generic UniProt reactions,
     12 from BKMS, 2 unsaturated-termination rows and 2 extra UniProt-style rows.
   - 18 single cells (`row_exceptions.csv`).
5. **Ordering**: the released file orders some JSON keys in a fixed way that the script
   reproduces: compound names sorted case-insensitively; ACP-proxy entries in order of
   appearance in the reaction, with `malonyl-[ACP]` last.

## Files
| File | Content |
|---|---|
| `build_ecoli_csv.py` | the build (`--uniprot` = the TSV from step 1) |
| `data/ecmdb_name_map.csv` | 414 UniProt names -> ECMDB name (exact name, short alias list, or name without bracketed parts; case-insensitive) |
| `data/compound_smiles.csv`, `data/smiles_source_labels.csv` | per-name SMILES and source label (ECMDB, DB_CSV, PUBCHEM, ACP_PPANT_PROXY, HARDCODED, CHEBI_SCOPED, CORRECTED_DOUBLE_BOND, manual fix 2026-11, NOT_FOUND, N_A); ACP proxy provenance in the last column |
| `data/enzyme_names.csv`, `data/dropped_rows.csv`, `data/added_rows.csv`, `data/row_exceptions.csv` | hand curation (rows keyed by position in the 2026-04-20 UniProt file) |
| `lookup_smiles.py` | the SMILES lookup method (generic -> hard-coded -> ECMDB -> db.csv -> PubChem); not needed for the build. Needs `ecmdb.json` (not shipped, see below), the UniProt TSV (`--uniprot`) and `../db_csv/` (for its SMILES table and PubChem helper). Offline it agrees with the frozen table for 406/417 ECMDB, 370/381 DB_CSV and all HARDCODED/N_A names; `--pubchem` adds the PubChem tier |
| `extract_tables.py` | one-time script that cut the data tables out of the released file (needs `ecmdb.json`, the UniProt TSV and the released `ecoli.csv`) |
| `audit_curation.py` | checks every curation entry against the git history of `db/ecoli.csv` (needs a BEES clone) |

`ecmdb.json` (for `lookup_smiles.py` and `extract_tables.py` only, not for the build) is the
ECMDB 2015 JSON download: https://ecmdb.ca/downloads

## Provenance of the tables
The SMILES and curation tables were cut from the released file, so on their own they only
restate it. `audit_curation.py --repo <BEES clone>` checks each curation entry against git
independently. Result on 2026-11-09:
- Every entry first appears in a real commit of `db/ecoli.csv`.
- 70 enzyme names, 14 exceptions and 1 dropped row date from the first commit (`1f5b83fd`, 2026-05-05).
- 52 added rows: `2c0e0ac6` (2026-07-28). 2 added rows: `bc4e1ee8` (2026-10-07).
- 8 dropped rows and 4 exceptions: `a329b35a` (2026-08-17).
- 5 dropped SurE rows: `02624a2` (2026-11-17).
- `added_rows.csv` `row_number` values are positions in the released file. They were
  re-based on 2026-11-20 after the SurE removal (each added row moved up 5 places).

Earlier simulation check: on 2026-11-09 (before the SurE removal) the production FAS model run
with the rebuilt files gave byte-identical `reactions_summary.txt`, core/edge CSVs and SBML
(date line aside). Not re-run since.

## Data sources and licences
| Source | Used for | Terms |
|---|---|---|
| UniProtKB (UP000000625, 2026-04-20) | proteins, reactions, EC numbers (not shipped; step 1) | CC BY 4.0; cite The UniProt Consortium |
| ECMDB 2015 (Sajed et al. 2016, NAR; https://ecmdb.ca/downloads) | compound names (`ecmdb_name_map.csv`) and SMILES labelled ECMDB | free for academic use with citation; check ECMDB terms before commercial use |
| PubChem (https://pubchem.ncbi.nlm.nih.gov) | SMILES labelled PUBCHEM and DB_CSV (via db.csv) | NCBI public data; cite Kim et al. (PubChem) |
| ChEBI | 4 SMILES labelled CHEBI_SCOPED, 3 of the 5 manual fixes | CC BY 4.0 |
| BKMS-react | 12 added rows inferred from BKMS reactions (`bkms_chainlength_inferred`) | BKMS licence (https://bkms.brenda-enzymes.org/license.php); same permission question as db.csv |

## Known data issues (kept, because the paper used this file)
- 39 SMILES entries (35 rows, 9 compounds) have doubled backslashes (`\\`) from an old TSV escape step.
  One of the compounds is (9Z)-octadecenoate (oleate) in the TesA rows, which is part of FAS.
  RDKit reads the same molecule (same InChIKey) with the single or doubled backslash for all
  9 compounds, so the chemistry is unaffected.
- ECMDB ion entries are coarse: Fe(3+) carries the Fe(2+) SMILES, and Zn(2+) is neutral `[Zn]`.
- Wrong ECMDB structures for CMP, CO and two sugar phosphates were fixed on 2026-11-03 (label "manual fix 2026-11").
- Two released rows have CRLF line endings (see the top of this file).
