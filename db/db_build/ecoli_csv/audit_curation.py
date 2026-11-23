#!/usr/bin/env python3
"""Trace every hand-curation entry in data/ to the git commit that introduced it.

The curation tables were cut from the released ecoli.csv (extract_tables.py). This
script checks them independently against the history of db/ecoli.csv in a BEES clone:
for each entry it reports the first commit whose ecoli.csv already contains it.

Usage:
  python audit_curation.py --repo /path/to/BEES
"""

from __future__ import annotations

import argparse
import collections
import csv
import io
import json
import subprocess
from pathlib import Path

csv.field_size_limit(10**9)

DATA = Path(__file__).resolve().parent / "data"


def history(repo: Path) -> list[tuple[str, list[dict]]]:
    """(commit, rows of db/ecoli.csv) for every commit that changed the file, oldest first."""
    log = subprocess.run(["git", "-C", str(repo), "log", "--all", "--reverse", "--format=%h %ad",
                          "--date=short", "--", "db/ecoli.csv"], capture_output=True, text=True, check=True)
    out = []
    for line in log.stdout.split("\n"):
        if not line:
            continue
        commit = line.split()[0]
        text = subprocess.run(["git", "-C", str(repo), "show", f"{commit}:db/ecoli.csv"],
                              capture_output=True, text=True, check=True).stdout
        out.append((line, list(csv.DictReader(io.StringIO(text, newline="")))))
    return out


def read(name: str) -> list[dict]:
    with (DATA / name).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def cell(row: dict, field: str, key: str):
    meta = json.loads(row["meta"])
    target = {"compound_smiles": json.loads(row["compound_smiles"]),
              "compound_smiles_source": meta.get("compound_smiles_source", {}),
              "meta": meta}[field]
    return target.get(key)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo", type=Path, required=True)
    args = ap.parse_args()
    commits = history(args.repo)

    def first(test) -> str:
        return next((c for c, rows in commits if test(rows)), "NOT FOUND")

    checks = {
        "enzyme_names.csv": [lambda rows, e=e: any(r["index"] == e["index"] and r["enzyme"] == e["enzyme"] for r in rows)
                             for e in read("enzyme_names.csv")],
        "dropped_rows.csv": [lambda rows, e=e: not any(r["index"] == e["index"] and r["Entry"] == e["Entry"] for r in rows)
                             for e in read("dropped_rows.csv")],
        "added_rows.csv": [lambda rows, e=e: any(r["index"] == e["index"] and r["reaction"] == e["reaction"]
                                                 and r["enzyme"] == e["enzyme"]
                                                 and json.loads(r["compound_smiles"]) == json.loads(e["compound_smiles"])
                                                 for r in rows)
                           for e in read("added_rows.csv")],
        "row_exceptions.csv": [lambda rows, e=e: any(r["index"] == e["index"]
                                                     and cell(r, e["field"], e["key"]) == (json.loads(e["value"]) if e["value"] else None)
                                                     for r in rows)
                               for e in read("row_exceptions.csv")],
    }
    for table, tests in checks.items():
        counts = collections.Counter(first(t) for t in tests)
        print(f"{table} ({len(tests)} entries), first commit containing each entry:")
        for commit, n in counts.items():
            print(f"  {n:3d}  {commit}")


if __name__ == "__main__":
    main()
