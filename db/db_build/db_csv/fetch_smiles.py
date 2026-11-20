#!/usr/bin/env python3
"""Look up SMILES for every compound name in db.csv (how data/compound_smiles.csv was made).

This is a clean rewrite of the two PubChem lookup passes that produced the SMILES in
BEES db.csv (Feb 2026). It is NOT needed to build db.csv: build_db_csv.py uses the
frozen table data/compound_smiles.csv. PubChem changes over time, so a fresh run will
give a slightly different table; use --compare to see how close it is.

For each compound name:
  1. Clean the name (strip quotes, leading "> ", trailing "<"/">", extra spaces).
  2. PubChem name search (PUG REST) for the cleaned name, then for simple variants
     (alpha <-> Greek letter, "- " / " -" -> "-"). First hit wins.
Lookups are cached in a JSON file, so the run can be stopped and resumed.

A ChEMBL free-text fallback was also tried at the time but never merged into db.csv:
it returns wrong molecules for vague names (e.g. "'activated' tRNA" -> "[C]").

Usage:
  python fetch_smiles.py --db db.csv --out compound_smiles_refetched.csv \
      --compare data/compound_smiles.csv
  python fetch_smiles.py --db db.csv --limit 20     # quick test
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

csv.field_size_limit(10**9)

PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/SMILES/JSON"
GREEK = [("alpha", "\u03b1"), ("beta", "\u03b2"), ("gamma", "\u03b3"), ("delta", "\u03b4"),
         ("epsilon", "\u03b5"), ("lambda", "\u03bb"), ("omega", "\u03c9")]


def compound_names(db_path: Path) -> list[str]:
    """All species names in the BKMS rows of db.csv, sorted."""
    names: set[str] = set()
    with db_path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["meta"] and json.loads(row["meta"]).get("source") == "BKMS-react":
                names.update(json.loads(row["stoichiometry"]))
    return sorted(names)


def clean_name(name: str) -> str:
    name = name.strip().strip("'").strip('"')
    name = re.sub(r"^>\s*", "", name)
    name = re.sub(r"\s*[<>]$", "", name)
    return re.sub(r"\s+", " ", name).strip()


def name_variants(name: str) -> list[str]:
    variants = [name]
    for word, letter in GREEK:
        if word in name:
            variants.append(name.replace(word, letter))
        if letter in name:
            variants.append(name.replace(letter, word))
    variants.append(re.sub(r"\s+-", "-", re.sub(r"-\s+", "-", name)))
    return list(dict.fromkeys(v for v in variants if v))


def get_json(url: str, timeout: float, retries: int = 3) -> dict | None:
    """GET a JSON document; None on 404 or after repeated failures."""
    for attempt in range(retries + 1):
        try:
            with urlopen(Request(url, headers={"User-Agent": "BEES"}), timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except HTTPError as exc:
            if exc.code == 404:
                return None
        except (URLError, TimeoutError, json.JSONDecodeError):
            pass
        time.sleep(1.0 * (attempt + 1))
    return None


def pubchem(name: str, timeout: float) -> tuple[str, str] | None:
    data = get_json(PUBCHEM_URL.format(quote(name, safe="")), timeout)
    props = (data or {}).get("PropertyTable", {}).get("Properties") or []
    if props and props[0].get("SMILES"):
        return props[0]["SMILES"], f"https://pubchem.ncbi.nlm.nih.gov/compound/{props[0]['CID']}"
    return None


def lookup(name: str, cache: dict, timeout: float, pause: float) -> dict:
    """Return {'smiles', 'url'} for one compound name (empty smiles if not found)."""
    for variant in name_variants(clean_name(name)):
        if variant not in cache:
            hit = pubchem(variant, timeout)
            cache[variant] = list(hit) if hit else None
            time.sleep(pause)
        if cache[variant]:
            smiles, url = cache[variant]
            return {"smiles": smiles, "url": url}
    return {"smiles": "", "url": ""}


def compare(found: dict[str, dict], frozen_path: Path) -> None:
    with frozen_path.open(newline="", encoding="utf-8") as f:
        frozen = {r["name"]: r["smiles"] for r in csv.DictReader(f)}
    names = [n for n in found if n in frozen or found[n]["smiles"]]
    same = sum(found[n]["smiles"] == frozen.get(n) for n in names)
    print(f"Compared with {frozen_path.name}: {same}/{len(names)} names give the same SMILES")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", type=Path, required=True, help="db.csv built by build_db_csv.py")
    ap.add_argument("--out", type=Path, default=Path("compound_smiles_refetched.csv"))
    ap.add_argument("--cache", type=Path, default=Path("smiles_lookup_cache.json"))
    ap.add_argument("--compare", type=Path, help="frozen table to compare with")
    ap.add_argument("--limit", type=int, help="only the first N names (testing)")
    ap.add_argument("--timeout", type=float, default=20.0)
    ap.add_argument("--pause", type=float, default=0.2, help="seconds between web requests")
    args = ap.parse_args()

    names = compound_names(args.db)[: args.limit]
    cache = json.loads(args.cache.read_text()) if args.cache.exists() else {}

    found: dict[str, dict] = {}
    for i, name in enumerate(names, 1):
        found[name] = lookup(name, cache, args.timeout, args.pause)
        if i % 200 == 0 or i == len(names):
            args.cache.write_text(json.dumps(cache, indent=1, sort_keys=True))
            print(f"{i}/{len(names)} names looked up")

    with args.out.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["name", "smiles", "url"])
        for name, hit in found.items():
            if hit["smiles"]:
                w.writerow([name, hit["smiles"], hit["url"]])

    n_found = sum(bool(h["smiles"]) for h in found.values())
    print(f"SMILES found for {n_found}/{len(names)} names -> {args.out}")
    if args.compare:
        compare(found, args.compare)


if __name__ == "__main__":
    main()
