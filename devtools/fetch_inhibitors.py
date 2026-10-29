#!/usr/bin/env python3

"""Median Ki per (target, compound) from the ChEMBL REST API, for choosing Ki values.

Keeps binding assays (assay_type B) against the human target with an exact Ki ("="), and
drops every activity that carries a data_validity_comment. Records of one target/compound
pair are reduced to their median; the output is a CSV (values in nM and mM, with the ids of
the activities behind each median so it can be checked against ChEMBL).

    python devtools/fetch_inhibitors.py --target CHEMBL202 --compound CHEMBL34259 -o ki.csv

ChEMBL data are CC BY-SA 3.0 (https://www.ebi.ac.uk/chembl/); credit them where you use a value.
"""

import argparse
import csv
import statistics
import sys
from collections import defaultdict
from typing import Callable, Dict, Iterable, Iterator, List, Optional

BASE_URL = "https://www.ebi.ac.uk"
ENDPOINT = "/chembl/api/data/activity.json"
TO_NM = {"nM": 1.0, "uM": 1e3, "µM": 1e3, "mM": 1e6, "pM": 1e-3, "M": 1e9}
FIELDS = ["target_chembl_id", "molecule_chembl_id", "n_records", "median_ki_nM", "median_ki_mM",
          "min_ki_nM", "max_ki_nM", "activity_ids"]


def _http_get(url: str, params: Optional[dict] = None) -> dict:
    import requests

    resp = requests.get(url, params=params, timeout=60, headers={"Accept": "application/json"})
    resp.raise_for_status()
    return resp.json()


def iter_activities(target: str, compound: Optional[str] = None, get: Callable = _http_get,
                    organism: str = "Homo sapiens", page_size: int = 1000) -> Iterator[dict]:
    """Every Ki activity of ``target`` (optionally one compound) that passes the ChEMBL-side filters."""
    params = {
        "target_chembl_id": target, "standard_type": "Ki", "assay_type": "B",
        "target_organism": organism, "data_validity_comment__isnull": "true",
        "limit": page_size,
    }
    if compound:
        params["molecule_chembl_id"] = compound
    url, query = BASE_URL + ENDPOINT, params
    while url:
        payload = get(url, query)
        yield from payload.get("activities", [])
        nxt = (payload.get("page_meta") or {}).get("next")
        url, query = (BASE_URL + nxt if nxt and nxt.startswith("/") else nxt), None


def ki_nM(activity: dict) -> Optional[float]:
    """Ki in nM of one activity, or None when it must not be used."""
    if activity.get("data_validity_comment"):  # also guards against a server ignoring the filter
        return None
    if activity.get("standard_type") != "Ki" or activity.get("standard_relation") not in ("=", None):
        return None
    if activity.get("assay_type") not in ("B", None):
        return None
    try:
        value = float(activity["standard_value"])
    except (KeyError, TypeError, ValueError):
        return None
    factor = TO_NM.get(activity.get("standard_units"))
    if factor is None or value <= 0:
        return None
    return value * factor


def median_ki(activities: Iterable[dict]) -> List[dict]:
    """One row per (target, compound): median Ki over the usable activities."""
    groups: Dict[tuple, List[tuple]] = defaultdict(list)
    for act in activities:
        value = ki_nM(act)
        if value is not None:
            groups[(act.get("target_chembl_id"), act.get("molecule_chembl_id"))].append((value, act.get("activity_id")))
    rows = []
    for (target, molecule), recs in sorted(groups.items(), key=lambda kv: (str(kv[0][0]), str(kv[0][1]))):
        values = [v for v, _ in recs]
        med = statistics.median(values)
        rows.append({
            "target_chembl_id": target, "molecule_chembl_id": molecule, "n_records": len(values),
            "median_ki_nM": med, "median_ki_mM": med * 1e-6, "min_ki_nM": min(values), "max_ki_nM": max(values),
            "activity_ids": ";".join(str(a) for _, a in recs),
        })
    return rows


def write_csv(rows: List[dict], path: str) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in row.items()})


def fetch(targets: Iterable[str], compounds: Optional[Iterable[str]] = None, get: Callable = _http_get) -> List[dict]:
    activities: List[dict] = []
    for target in targets:
        for compound in (list(compounds) if compounds else [None]):
            activities.extend(iter_activities(target, compound, get=get))
    return median_ki(activities)


def main(argv=None) -> List[dict]:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--target", action="append", required=True, help="ChEMBL target id (repeatable)")
    p.add_argument("--compound", action="append", help="restrict to this ChEMBL molecule id (repeatable)")
    p.add_argument("-o", "--output", default="chembl_ki.csv")
    args = p.parse_args(argv)
    try:
        rows = fetch(args.target, args.compound)
    except Exception as exc:  # network blocked, API change, ...
        sys.exit(f"ChEMBL request failed: {exc}")
    if not rows:
        print("No usable Ki records.")
    write_csv(rows, args.output)
    for r in rows:
        print(f"{r['target_chembl_id']} {r['molecule_chembl_id']}: median Ki {r['median_ki_nM']:.4g} nM "
              f"({r['n_records']} records)")
    print(f"wrote {args.output}")
    return rows


if __name__ == "__main__":
    main()
