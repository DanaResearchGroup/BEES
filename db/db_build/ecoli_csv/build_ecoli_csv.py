#!/usr/bin/env python3
"""Build the BEES E. coli reaction database (db/ecoli.csv) from the UniProt download.

Source data:
  UniProt proteome UP000000625 (E. coli K-12, https://www.uniprot.org/proteomes/UP000000625),
  as a TSV with the columns Entry, Reviewed, Entry Name, Protein names, Gene Names,
  Organism, Length, Catalytic activity, Sequence, EC number. Not shipped (CC BY 4.0, but
  UniProt changes between releases): download it as described in REBUILD_INSTRUCTIONS.md and pass it
  with --uniprot. The released db/ecoli.csv comes from the 2026-04-20 download.

Data files shipped with this script (see REBUILD_INSTRUCTIONS.md for how each was made):
  data/ecmdb_name_map.csv        UniProt compound name -> ECMDB compound name
  data/compound_smiles.csv       compound name -> SMILES
  data/smiles_source_labels.csv  compound name -> where its SMILES came from
  data/enzyme_names.csv          short enzyme names (FabA, FabZ, ...) for 70 rows
  data/dropped_rows.csv          14 UniProt rows removed by hand
  data/added_rows.csv            54 hand-curated / inferred FAS rows, with their row number
  data/row_exceptions.csv        18 cells that differ from the per-name tables

Steps:
  1. UniProt: one row per catalytic-activity reaction and EC number; protein names
     cleaned (EC numbers removed, parenthesised synonyms moved to "Protein aliases").
  2. Rename compounds that exist in ECMDB to their ECMDB names (stoichiometry and text).
  3. Attach SMILES and their source labels.
  4. Hand curation: enzyme names, dropped and added rows, row exceptions.
  5. Write ecoli.csv; with --check, compare it with a released file.

Usage:
  python build_ecoli_csv.py --uniprot uniprotkb_proteome_UP000000625_2026_04_20_with_catalytic_activity.tsv.gz --out ecoli.csv
  python build_ecoli_csv.py --uniprot <same file> --out ecoli.csv --check ../../ecoli.csv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
from pathlib import Path

csv.field_size_limit(10**9)

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"

COLUMNS = ["Entry", "Entry Name", "Protein aliases", "Sequence", "index", "enzyme",
           "ec_number", "reaction", "stoichiometry", "meta", "compound_smiles"]
JSON_COLUMNS = {"stoichiometry", "meta", "compound_smiles"}

CATALYTIC_ACTIVITY = re.compile(r"\bCATALYTIC ACTIVITY:\s*")
INLINE_EC = re.compile(r"(?:^|[;\s])EC=([^;]+?)\s*(?:;|$)")
SUBUNIT_PREFIX = re.compile(r"^\s*(\[[^\]]+\]):\s*")
EC_IN_NAME = re.compile(r"\(EC\s+([^)]+?)\)")
PARENTHESES = re.compile(r"\(([^()]*)\)")
BRACKETS = re.compile(r"\[([^\[\]]*?)\]")
SPACES = re.compile(r"\s+")
PLUS = re.compile(r"\s+\+\s+")
COEFFICIENT = re.compile(r"^\s*(\d+)\s+(.+?)\s*$")
LABEL_COEFFICIENT = re.compile(r"^(?:\d+|n)\s+")


# ---------------------------------------------------------------- step 1: UniProt

def read_uniprot(path: Path) -> list[dict]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE, escapechar="\\"))


def clean_protein_name(name: str) -> tuple[str, list[str], list[str]]:
    """'Foo synthase (EC 1.2.3.4) (Bar)' -> ('Foo synthase', ['1.2.3.4'], ['Bar'])."""
    ecs = [SPACES.sub(" ", m.group(1).strip()).replace("EC ", "").strip() for m in EC_IN_NAME.finditer(name)]
    rest = EC_IN_NAME.sub("", name)
    aliases = [SPACES.sub(" ", m.group(1).strip()) for m in PARENTHESES.finditer(rest)]
    aliases += [SPACES.sub(" ", m.group(1).strip()) for m in BRACKETS.finditer(rest)]
    clean = SPACES.sub(" ", BRACKETS.sub("", PARENTHESES.sub("", rest))).strip()
    return clean, [e for e in ecs if e], [a for a in aliases if a]


def split_catalytic_activity(cell: str) -> list[tuple[str | None, str, str | None]]:
    """Catalytic activity cell -> [(reaction, rest of the annotation, inline EC), ...]."""
    out = []
    for block in (b.strip() for b in CATALYTIC_ACTIVITY.split(cell or "")[1:]):
        if not block:
            continue
        prefix = ""
        m = SUBUNIT_PREFIX.match(block)
        if m:
            prefix, block = m.group(1).strip(), block[m.end():].lstrip()
        m_ec = INLINE_EC.search(block)
        inline_ec = m_ec.group(1).strip() if m_ec else None
        if "Reaction=" not in block:
            reaction, source = None, block.strip().strip(";").strip()
        else:
            after = block.split("Reaction=", 1)[1].strip()
            reaction, _, rest = after.partition(";")
            reaction, source = reaction.strip() or None, rest.strip().strip(";").strip()
        out.append((reaction, f"{prefix} {source}".strip() if prefix else source, inline_ec))
    return out


def parse_stoichiometry(reaction: str) -> dict[str, int] | None:
    """'2 A + B = C' -> {'A': -2, 'B': -1, 'C': 1}; a leading integer is a coefficient."""
    left, sep, right = reaction.partition("=")
    if not sep or not left.strip() or not right.strip():
        return None
    stoichiometry: dict[str, int] = {}
    for side, sign in ((left, -1), (right, 1)):
        for term in PLUS.split(side.strip()):
            term = SPACES.sub(" ", term.strip())
            if not term:
                continue
            m = COEFFICIENT.match(term)
            n, name = (int(m.group(1)), m.group(2).strip()) if m else (1, term)
            stoichiometry[name] = stoichiometry.get(name, 0) + sign * n
    return {k: v for k, v in stoichiometry.items() if v != 0} or None


def reaction_names(reaction: str) -> list[str]:
    """Compound names of a reaction in order of first appearance; '2 ' / 'n ' coefficients dropped."""
    names = []
    for side in reaction.split(" = ", 1):
        names += [LABEL_COEFFICIENT.sub("", t.strip()) for t in PLUS.split(side.strip()) if t.strip()]
    return list(dict.fromkeys(names))


def label_names(reaction: str) -> list[str]:
    """Compound names used for SMILES look-ups, sorted case-insensitively (the released order)."""
    return sorted(set(reaction_names(reaction)), key=lambda s: (s.lower(), s))


def uniprot_reactions(proteins: list[dict]) -> list[dict]:
    """One row per (protein, catalytic-activity reaction, EC number)."""
    rows = []
    for p in proteins:
        ec_column = [e.strip() for e in (p["EC number"] or "").split(";") if e.strip()]
        name, name_ecs, aliases = clean_protein_name(p["Protein names"] or "")
        sequence = p["Sequence"].strip()
        for reaction, source, inline_ec in split_catalytic_activity((p["Catalytic activity"] or "").strip()):
            ec = inline_ec or (ec_column[0] if len(ec_column) == 1 else "")
            if not ec and name_ecs:
                ec = "; ".join(dict.fromkeys(name_ecs))
            ecs = list(dict.fromkeys(e.replace("EC", "").strip() for e in ec.split(";") if e.replace("EC", "").strip()))
            stoichiometry = parse_stoichiometry(reaction) if reaction else None
            if not reaction or not ecs or not stoichiometry:
                continue
            for one_ec in ecs:
                rows.append({
                    "Entry": p["Entry"].strip(),
                    "Entry Name": p["Entry Name"].strip(),
                    "Protein aliases": "; ".join(dict.fromkeys(aliases)),
                    "Sequence": sequence,
                    "index": str(len(rows) + 1),
                    "enzyme": name or p["Entry Name"].strip() or p["Entry"].strip(),
                    "ec_number": f"EC {one_ec}",
                    "reaction": reaction,
                    "stoichiometry": stoichiometry,
                    "gene_names": p["Gene Names"].strip(),
                    "reaction_source": source,
                    "lookup_names": label_names(reaction),
                    "reaction_names": reaction_names(reaction),
                })
    return rows


# ---------------------------------------------------------------- step 2: ECMDB names

def rename_to_ecmdb(row: dict, name_map: dict[str, str]) -> list[str]:
    """Rename compounds in place; return the ECMDB names used (in stoichiometry order)."""
    stoichiometry, matched, renamed = {}, [], {}
    for name, n in row["stoichiometry"].items():
        new = name_map.get(name)
        if new:
            matched.append(new)
            if new != name:
                renamed[name] = new
        stoichiometry[new or name] = n
    reaction = row["reaction"]
    for old in sorted(renamed, key=len, reverse=True):
        reaction = re.sub(rf"\b{re.escape(old)}\b", renamed[old], reaction)
    row["stoichiometry"], row["reaction"] = stoichiometry, reaction
    return list(dict.fromkeys(matched))


# ---------------------------------------------------------------- step 3: SMILES

def attach_smiles(row: dict, ecmdb_names: list[str], smiles: dict, labels: dict, acp: dict) -> None:
    sources = {name: labels[name] for name in row["lookup_names"] if name in labels}
    for name in ecmdb_names:
        sources[name] = labels.get(name, "ECMDB")
    # Order of the proxy entries as in the released file: order of appearance in the
    # reaction, with malonyl-[ACP] last (it was added last in the original script).
    proxy_order = [n for n in row["reaction_names"] if n != "malonyl-[ACP]"] + \
                  [n for n in row["reaction_names"] if n == "malonyl-[ACP]"]
    proxies = {name: acp[name] for name in proxy_order
               if sources.get(name) == "ACP_PPANT_PROXY" and name in acp}

    meta = {
        "source": "uniprot_tsv",
        "uniprot_entry": row["Entry"],
        "entry_name": row["Entry Name"],
        "gene_names": row["gene_names"],
        "protein_aliases": row["Protein aliases"],
        "reaction_source": row["reaction_source"],
        "compound_smiles_source": sources,
        "sequence_length": len(row["Sequence"]),
        "sequence_sha1": hashlib.sha1(row["Sequence"].encode()).hexdigest() if row["Sequence"] else "",
    }
    if proxies:
        meta["acp_ppant_proxy"] = proxies
    row["meta"] = meta
    row["compound_smiles"] = {name: smiles[name] for name in row["stoichiometry"] if name in smiles}


# ---------------------------------------------------------------- step 4: curation

def apply_exceptions(row: dict, exceptions: list[dict]) -> None:
    for ex in exceptions:
        value = json.loads(ex["value"]) if ex["value"] else None
        target = {"compound_smiles": row["compound_smiles"],
                  "compound_smiles_source": row["meta"]["compound_smiles_source"],
                  "meta": row["meta"]}[ex["field"]]
        if value is None:
            target.pop(ex["key"], None)
        else:
            target[ex["key"]] = value


def to_csv_row(row: dict) -> dict:
    out = {c: row[c] for c in COLUMNS if c not in JSON_COLUMNS}
    out["stoichiometry"] = json.dumps(row["stoichiometry"], ensure_ascii=False)
    out["meta"] = json.dumps(row["meta"], ensure_ascii=False)
    out["compound_smiles"] = json.dumps(row["compound_smiles"], ensure_ascii=False)
    return out


# ---------------------------------------------------------------- tables

def read_table(name: str) -> list[dict]:
    with (DATA / name).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def load_tables() -> dict:
    labels = read_table("smiles_source_labels.csv")
    exceptions: dict[str, list[dict]] = {}
    for ex in read_table("row_exceptions.csv"):
        exceptions.setdefault(ex["index"], []).append(ex)
    return {
        "name_map": {r["name"]: r["ecmdb_name"] for r in read_table("ecmdb_name_map.csv")},
        "smiles": {r["name"]: r["smiles"] for r in read_table("compound_smiles.csv")},
        "labels": {r["name"]: r["label"] for r in labels},
        "acp": {r["name"]: json.loads(r["acp_ppant_proxy"]) for r in labels if r["acp_ppant_proxy"]},
        "enzyme": {r["index"]: r["enzyme"] for r in read_table("enzyme_names.csv")},
        "dropped": {r["index"] for r in read_table("dropped_rows.csv")},
        "added": read_table("added_rows.csv"),
        "exceptions": exceptions,
    }


def build(tables: dict, uniprot: Path) -> list[dict]:
    rows = uniprot_reactions(read_uniprot(uniprot))
    out = []
    for row in rows:
        ecmdb_names = rename_to_ecmdb(row, tables["name_map"])
        attach_smiles(row, ecmdb_names, tables["smiles"], tables["labels"], tables["acp"])
        apply_exceptions(row, tables["exceptions"].get(row["index"], []))
        row["enzyme"] = tables["enzyme"].get(row["index"], row["enzyme"])
        if row["index"] not in tables["dropped"]:
            out.append(to_csv_row(row))

    for added in sorted(tables["added"], key=lambda r: int(r["row_number"])):
        out.insert(int(added["row_number"]) - 1, {c: added[c] for c in COLUMNS})
    return out


# ---------------------------------------------------------------- check

def compare(rows: list[dict], released: Path) -> None:
    with released.open(newline="", encoding="utf-8") as f:
        ref = list(csv.DictReader(f))
    print(f"Rows: built {len(rows)}, released {len(ref)}")
    diffs: dict[str, int] = {}
    for a, b in zip(rows, ref):
        for c in COLUMNS:
            same = json.loads(a[c]) == json.loads(b[c]) if c in JSON_COLUMNS else a[c] == b[c]
            if not same:
                diffs[c] = diffs.get(c, 0) + 1
    if len(rows) == len(ref) and not diffs:
        print("Same content as the released file (every row and column; JSON compared as data)")
    else:
        print(f"Differences per column: {diffs}")


def main() -> None:
    ap = argparse.ArgumentParser(description="Build BEES ecoli.csv from UniProt.")
    ap.add_argument("--uniprot", type=Path, required=True,
                    help="UniProt TSV.gz with catalytic activity (see REBUILD_INSTRUCTIONS.md for the download)")
    ap.add_argument("--out", type=Path, default=Path("ecoli.csv"))
    ap.add_argument("--check", type=Path, help="released ecoli.csv to compare with")
    args = ap.parse_args()

    rows = build(load_tables(), args.uniprot)
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=COLUMNS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    data = buffer.getvalue().encode("utf-8")
    args.out.write_bytes(data)
    print(f"Wrote {args.out} ({len(rows)} rows), SHA-256 {hashlib.sha256(data).hexdigest()}")

    if args.check:
        compare(rows, args.check)
        print("Byte-identical" if data == args.check.read_bytes() else "Not byte-identical")


if __name__ == "__main__":
    main()
