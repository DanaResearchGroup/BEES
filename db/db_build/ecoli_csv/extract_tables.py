#!/usr/bin/env python3
"""Create the data tables used by build_ecoli_csv.py (run once; kept for transparency).

Inputs:
  --uniprot   the UniProt TSV.gz with catalytic activity, 2026-04-20 (see REBUILD_INSTRUCTIONS.md; not shipped)
  --ecmdb     ecmdb.json (ECMDB 2015 download, not shipped: https://ecmdb.ca/downloads),
              used only for the name map
  --released  the released db/ecoli.csv (BEES v1.0.0)

The name map is computed with the rules of the original 2026-05 script (exact name,
a short alias list, or name without bracketed parts; case-insensitive). Everything
else is read from the released file: per-name SMILES and labels, and whatever the
rules cannot explain (enzyme names, dropped/added rows, a few cells) as curation tables.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import re
from pathlib import Path

import build_ecoli_csv as b

ECMDB_ALIASES = {
    "atp": "adenosine triphosphate", "adp": "adenosine diphosphate", "amp": "adenosine monophosphate",
    "nad(+)": "nicotinamide adenine dinucleotide", "nadp(+)": "nicotinamide adenine dinucleotide phosphate",
    "nadph": "nicotinamide adenine dinucleotide phosphate - reduced",
    "nadh": "nicotinamide adenine dinucleotide - reduced", "h2o": "water", "co2": "carbon dioxide",
    "nh3": "ammonia", "nh4(+)": "ammonium", "phosphate": "phosphoric acid", "pi": "phosphoric acid",
    "ppi": "pyrophosphate", "coa": "coenzyme a", "o2": "oxygen", "h(+)": "hydron",
    "glutamate": "l-glutamic acid", "l-glutamate": "l-glutamic acid", "pyruvate": "pyruvic acid",
    "succinate": "succinic acid", "fumarate": "fumaric acid", "malate": "l-malic acid",
    "citrate": "citric acid", "isocitrate": "isocitric acid",
}


def write_table(name: str, header: list[str], rows: list[list]) -> None:
    with (b.DATA / name).open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)
    print(f"{name}: {len(rows)} rows")


def ecmdb_match(name: str, ecmdb: dict[str, str]) -> str | None:
    key = name.lower().strip()
    if key in ecmdb:
        return ecmdb[key]
    if key in ECMDB_ALIASES and ECMDB_ALIASES[key] in ecmdb:
        return ecmdb[ECMDB_ALIASES[key]]
    return ecmdb.get(re.sub(r"\(.*?\)", "", key).strip())


def most_common(values: dict) -> dict:
    return {k: c.most_common(1)[0][0] for k, c in values.items()}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--uniprot", type=Path, required=True)
    ap.add_argument("--ecmdb", type=Path, required=True)
    ap.add_argument("--released", type=Path, required=True)
    args = ap.parse_args()

    b.DATA.mkdir(exist_ok=True)
    with args.released.open(newline="", encoding="utf-8") as f:
        released = list(csv.DictReader(f))
    by_index = {r["index"]: r for r in released}

    rows = b.uniprot_reactions(b.read_uniprot(args.uniprot))
    built_index = {r["index"] for r in rows}

    # ECMDB name map
    ecmdb = {e["name"].lower().strip(): e["name"] for e in json.loads(args.ecmdb.read_text())
             if e.get("moldb_smiles") and e.get("name")}
    names = sorted({n for r in rows for n in r["stoichiometry"]})
    name_map = {n: m for n in names if (m := ecmdb_match(n, ecmdb))}
    write_table("ecmdb_name_map.csv", ["name", "ecmdb_name"], [[n, m] for n, m in name_map.items()])

    # Per-name SMILES, labels and ACP proxy provenance (majority over UniProt rows)
    smiles, labels, acp = (collections.defaultdict(collections.Counter) for _ in range(3))
    for r in released:
        if r["index"] not in built_index:
            continue
        meta = json.loads(r["meta"])
        for n, s in json.loads(r["compound_smiles"]).items():
            smiles[n][s] += 1
        for n, lab in meta.get("compound_smiles_source", {}).items():
            labels[n][lab] += 1
        for n, prov in meta.get("acp_ppant_proxy", {}).items():
            acp[n][json.dumps(prov, ensure_ascii=False)] += 1
    smiles, labels, acp = most_common(smiles), most_common(labels), most_common(acp)
    write_table("compound_smiles.csv", ["name", "smiles"], sorted(smiles.items()))
    write_table("smiles_source_labels.csv", ["name", "label", "acp_ppant_proxy"],
                [[n, lab, acp.get(n, "")] for n, lab in sorted(labels.items())])

    # Hand curation
    write_table("enzyme_names.csv", ["index", "enzyme"],
                [[r["index"], by_index[r["index"]]["enzyme"]] for r in rows
                 if r["index"] in by_index and by_index[r["index"]]["enzyme"] != r["enzyme"]])
    write_table("dropped_rows.csv", ["index", "Entry", "reaction"],
                [[r["index"], r["Entry"], r["reaction"]] for r in rows if r["index"] not in by_index])
    added = [[i] + [r[c] for c in b.COLUMNS]
             for i, r in enumerate(released, 1) if r["index"] not in built_index]
    write_table("added_rows.csv", ["row_number"] + b.COLUMNS, added)

    # Cells the per-name tables do not explain
    write_table("row_exceptions.csv", ["index", "field", "key", "value"], [])
    tables = b.load_tables()
    exceptions = []
    for row in rows:
        ref = by_index.get(row["index"])
        if ref is None:
            continue
        ecmdb_names = b.rename_to_ecmdb(row, tables["name_map"])
        b.attach_smiles(row, ecmdb_names, tables["smiles"], tables["labels"], tables["acp"])
        ref_meta, ref_smiles = json.loads(ref["meta"]), json.loads(ref["compound_smiles"])
        pairs = [("compound_smiles", row["compound_smiles"], ref_smiles),
                 ("compound_smiles_source", row["meta"]["compound_smiles_source"], ref_meta["compound_smiles_source"]),
                 ("meta", {k: v for k, v in row["meta"].items() if k != "compound_smiles_source"},
                  {k: v for k, v in ref_meta.items() if k != "compound_smiles_source"})]
        for field, got, want in pairs:
            for key in sorted(set(got) | set(want)):
                if got.get(key) != want.get(key):
                    value = json.dumps(want[key], ensure_ascii=False) if key in want else ""
                    exceptions.append([row["index"], field, key, value])
    write_table("row_exceptions.csv", ["index", "field", "key", "value"], exceptions)


if __name__ == "__main__":
    main()
