#!/usr/bin/env python3
"""Look up SMILES for the compounds in ecoli.csv (the method behind data/compound_smiles.csv).

NOT needed to build ecoli.csv: build_ecoli_csv.py uses the frozen tables. This script
shows how the SMILES were found and reports how often a fresh lookup agrees.

For each compound name in the UniProt reactions (before and after ECMDB renaming),
the first tier that answers wins:
  1. generic names ("an acceptor", "...-[protein]", ...)        -> N_A, no SMILES
  2. a short hard-coded list (water, H+, O2, CO2, phosphate, ...) -> HARDCODED
  3. ECMDB 2015 (name / IUPAC / traditional name), trying the name
     without "(in)"/"(out)" and without a charge such as "(+)"   -> ECMDB
  4. the db.csv SMILES table (../db_csv/data/compound_smiles.csv) -> DB_CSV
  5. with --pubchem: PubChem name search                          -> PUBCHEM
  otherwise NOT_FOUND.

Not reproduced here (kept as-is in the frozen tables): 4 ChEBI matches scoped to the
reaction's ChEBI cross-references, 50 acyl-ACP structures built by a
phosphopantetheine proxy, 3 corrected double bonds and 5 manual fixes (2026-11).

Input not shipped: ecmdb.json, the ECMDB 2015 JSON download (https://ecmdb.ca/downloads).

Usage:
  python lookup_smiles.py --ecmdb ecmdb.json --uniprot <UniProt TSV.gz> [--pubchem] [--out compound_smiles_lookup.csv]
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import re
import sys
import time
from pathlib import Path

import build_ecoli_csv as b

sys.path.insert(0, str(b.HERE.parent / "db_csv"))
import fetch_smiles as pubchem_lookup  # noqa: E402

DB_CSV_SMILES = b.HERE.parent / "db_csv" / "data" / "compound_smiles.csv"

GENERIC = re.compile(r"\bacceptor\b|\bdonor\b|\bdonated\s+group\b|\ban\s+aldehyde\b|\ban\s+alcohol\b"
                     r"|\bside\s+chain\b|\br-?group\b|-r$|\[protein\]|\[enzyme\]", re.IGNORECASE)
HARDCODED = {
    "h2o": "O", "water": "O", "h(+)": "[H+]", "h+": "[H+]", "proton": "[H+]",
    "o2": "O=O", "dioxygen": "O=O", "co2": "O=C=O", "carbon dioxide": "O=C=O",
    "pi": "OP(=O)(O)O", "phosphate": "OP(=O)(O)O", "orthophosphate": "OP(=O)(O)O",
    "hydrogenphosphate": "OP(=O)(O)O", "nh3": "N", "ammonia": "N",
    "hco3-": "O=C([O-])O", "bicarbonate": "O=C([O-])O",
    "diphosphate": "[O-]P(=O)([O-])OP(=O)([O-])[O-]", "pyrophosphate": "[O-]P(=O)([O-])OP(=O)([O-])[O-]",
    "ppi": "[O-]P(=O)([O-])OP(=O)([O-])[O-]",
}
COMPARTMENT = re.compile(r"\((?:in|out)\)\s*$", re.IGNORECASE)
CHARGE = re.compile(r"\((?:\+|-|\d+[+-])\)\s*$")


def key(name: str) -> str:
    return " ".join(name.strip().lower().split())


def candidates(name: str) -> list[str]:
    no_compartment = COMPARTMENT.sub("", name).strip()
    forms = [name, no_compartment, CHARGE.sub("", name).strip(), CHARGE.sub("", no_compartment).strip()]
    return list(dict.fromkeys(key(f) for f in forms if f))


def load_ecmdb(path: Path) -> dict[str, str]:
    """Main compound names first (a later entry with the same name wins), then IUPAC names."""
    records = [r for r in json.loads(path.read_text()) if (r.get("moldb_smiles") or "").strip()]
    by_name = {key(r["name"]): r["moldb_smiles"].strip() for r in records if (r.get("name") or "").strip()}
    for rec in records:
        for field in ("moldb_iupac", "moldb_traditional_iupac"):
            if (rec.get(field) or "").strip():
                by_name.setdefault(key(rec[field]), rec["moldb_smiles"].strip())
    return by_name


def load_db_csv_smiles() -> dict[str, str]:
    with DB_CSV_SMILES.open(newline="", encoding="utf-8") as f:
        out: dict[str, str] = {}
        for r in csv.DictReader(f):
            if r["smiles"]:
                out.setdefault(key(r["name"]), r["smiles"])
        return out


def lookup(name: str, ecmdb: dict, db_csv: dict, use_pubchem: bool, cache: dict) -> tuple[str, str]:
    if GENERIC.search(name):
        return "", "N_A"
    forms = candidates(name)
    for table, label in ((HARDCODED, "HARDCODED"), (ecmdb, "ECMDB"), (db_csv, "DB_CSV")):
        for form in forms:
            if form in table:
                return table[form], label
    if use_pubchem:
        if name not in cache:
            hit = pubchem_lookup.lookup(name, {}, timeout=20.0, pause=0.0)
            cache[name] = hit["smiles"]
            time.sleep(0.2)
        if cache[name]:
            return cache[name], "PUBCHEM"
    return "", "NOT_FOUND"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ecmdb", type=Path, required=True, help="ecmdb.json (ECMDB 2015 download)")
    ap.add_argument("--uniprot", type=Path, required=True, help="UniProt TSV.gz with catalytic activity")
    ap.add_argument("--pubchem", action="store_true", help="also query PubChem for names not found locally")
    ap.add_argument("--cache", type=Path, default=Path("pubchem_cache_ecoli.json"))
    ap.add_argument("--out", type=Path, default=Path("compound_smiles_lookup.csv"))
    args = ap.parse_args()

    tables = b.load_tables()
    names: set[str] = set()
    for row in b.uniprot_reactions(b.read_uniprot(args.uniprot)):
        names.update(row["lookup_names"])
        b.rename_to_ecmdb(row, tables["name_map"])
        names.update(row["stoichiometry"])

    ecmdb, db_csv = load_ecmdb(args.ecmdb), load_db_csv_smiles()
    cache = json.loads(args.cache.read_text()) if args.pubchem and args.cache.exists() else {}
    found = {n: lookup(n, ecmdb, db_csv, args.pubchem, cache) for n in sorted(names)}
    if args.pubchem:
        args.cache.write_text(json.dumps(cache, indent=1, sort_keys=True))

    with args.out.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["name", "smiles", "label"])
        w.writerows([n, s, lab] for n, (s, lab) in found.items())

    frozen_smiles, frozen_labels = tables["smiles"], tables["labels"]
    agree = collections.Counter()
    for n, (s, _) in found.items():
        if n not in frozen_smiles:
            continue
        frozen_label = frozen_labels.get(n, "")
        frozen_label = "manual fix 2026-11" if "manual fix" in frozen_label else frozen_label
        agree[(frozen_label, s == frozen_smiles[n])] += 1
    print(f"{len(found)} names -> {args.out}")
    print("Agreement with the frozen table (names that carry a SMILES cell), by frozen label:")
    for label in sorted({lab for lab, _ in agree}):
        same, total = agree[(label, True)], agree[(label, True)] + agree[(label, False)]
        print(f"  {label or '(no label)':28s} {same:5d} / {total:5d} same SMILES")


if __name__ == "__main__":
    main()
