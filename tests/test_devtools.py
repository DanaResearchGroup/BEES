"""
Tests for the devtools scripts: ChEMBL Ki fetcher (mocked), release sweep, target scan.

To run: pytest -v tests/test_devtools.py
"""

import csv
import importlib.util
import math
import os

import numpy as np
import pytest

from bees.common import BEES_PATH

DEV = os.path.join(BEES_PATH, "devtools")


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(DEV, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    import sys
    sys.path.insert(0, DEV)
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.path.remove(DEV)
    return mod


fetch = _load("fetch_inhibitors")


def _act(i, value, units="nM", molecule="CHEMBL34259", target="CHEMBL202", **kw):
    base = {"activity_id": i, "target_chembl_id": target, "molecule_chembl_id": molecule, "standard_type": "Ki",
            "standard_relation": "=", "standard_value": str(value), "standard_units": units,
            "assay_type": "B", "data_validity_comment": None}
    base.update(kw)
    return base


class TestFetcher:
    def test_median_per_pair_with_unit_conversion(self):
        rows = fetch.median_ki([
            _act(1, 0.077), _act(2, 0.179), _act(3, 0.01),
            _act(4, 2.0, molecule="CHEMBL1"), _act(5, 0.004, units="uM", molecule="CHEMBL1"),
        ])
        mtx = next(r for r in rows if r["molecule_chembl_id"] == "CHEMBL34259")
        assert mtx["n_records"] == 3 and mtx["median_ki_nM"] == pytest.approx(0.077)
        assert mtx["median_ki_mM"] == pytest.approx(7.7e-8)
        assert mtx["activity_ids"] == "1;2;3"
        other = next(r for r in rows if r["molecule_chembl_id"] == "CHEMBL1")
        assert other["median_ki_nM"] == pytest.approx(3.0)  # median of 2 nM and 4 nM

    @pytest.mark.parametrize("bad", [
        {"data_validity_comment": "Potential transcription error"}, {"standard_relation": ">"},
        {"standard_type": "IC50"}, {"assay_type": "F"}, {"standard_units": "%"}, {"standard_value": None},
        {"standard_value": "0"},
    ])
    def test_unusable_records_are_dropped(self, bad):
        assert fetch.median_ki([_act(1, 5.0, **bad)]) == []

    def test_pagination_and_request_filters(self):
        calls = []
        pages = {
            None: {"activities": [_act(1, 1.0), _act(2, 3.0)], "page_meta": {"next": "/chembl/api/data/activity.json?offset=2"}},
            "next": {"activities": [_act(3, 2.0)], "page_meta": {"next": None}},
        }

        def get(url, params=None):
            calls.append((url, params))
            return pages[None] if params else pages["next"]

        rows = fetch.fetch(["CHEMBL202"], ["CHEMBL34259"], get=get)
        assert rows[0]["n_records"] == 3 and rows[0]["median_ki_nM"] == pytest.approx(2.0)
        first = calls[0][1]
        assert first["assay_type"] == "B" and first["target_organism"] == "Homo sapiens"
        assert first["standard_type"] == "Ki" and first["data_validity_comment__isnull"] == "true"
        assert first["molecule_chembl_id"] == "CHEMBL34259"
        assert calls[1][0].endswith("offset=2") and calls[1][1] is None

    def test_cli_writes_csv(self, tmp_path, monkeypatch):
        monkeypatch.setattr(fetch, "_http_get", lambda url, params=None: {"activities": [_act(1, 0.1)], "page_meta": {}})
        monkeypatch.setattr(fetch, "fetch", lambda t, c=None, get=None: fetch.median_ki([_act(1, 0.1)]))
        out = tmp_path / "ki.csv"
        fetch.main(["--target", "CHEMBL202", "-o", str(out)])
        row = next(csv.DictReader(open(out)))
        assert float(row["median_ki_nM"]) == pytest.approx(0.1) and row["n_records"] == "1"

    def test_network_failure_is_reported_not_raised_as_a_traceback(self, monkeypatch):
        def boom(url, params=None):
            raise ConnectionError("blocked")
        monkeypatch.setattr(fetch, "_http_get", boom)
        with pytest.raises(SystemExit, match="ChEMBL request failed"):
            fetch.main(["--target", "CHEMBL202"])


sweep = _load("release_sweep")
scan = _load("target_scan")
common = _load("sweep_common")

MTX_INPUT = os.path.join(BEES_PATH, "projects", "methotrexate_dhfr", "input.yml")


def test_arms_deliver_the_same_total_dose():
    arms = sweep.make_arms("drug", 0.002, 1.0, 10000.0)
    assert arms["control"] == []
    for name in ("free bolus", "liposome", "slow depot"):
        (sch,) = arms[name]
        assert float(sch.released(1e12)) + 0.0 == pytest.approx(0.002, rel=1e-6)
        assert sch.elimination_k == pytest.approx(math.log(2) / 10000.0)
    assert arms["liposome"][0].params["k_fast"] == pytest.approx(math.log(2) / (2 * 3600))
    assert arms["slow depot"][0].params["tau"] == 24 * 3600 and arms["slow depot"][0].params["beta"] == 0.6


def test_time_scale_multiplies_release_times_only():
    base = sweep.make_arms("drug", 0.002, 1.0, 100.0)
    fast = sweep.make_arms("drug", 0.002, 0.1, 100.0)
    assert fast["slow depot"][0].params["tau"] == pytest.approx(0.1 * base["slow depot"][0].params["tau"])
    assert fast["liposome"][0].params["k_slow"] == pytest.approx(10 * base["liposome"][0].params["k_slow"])
    assert fast["free bolus"][0].elimination_k == base["free bolus"][0].elimination_k


def test_summed_output_skips_labels_not_in_the_model():
    traj = {"a": np.array([1.0, 2.0]), "b": np.array([10.0, 20.0])}
    total, used, skipped = common.summed_output(traj, ["A", "b", "ghost"])
    assert list(total) == [11.0, 22.0] and used == ["A", "b"] and skipped == ["ghost"]
    with pytest.raises(ValueError, match="None of the output species"):
        common.summed_output(traj, ["ghost"])
    assert common.parse_labels(" a ; b;;c ") == ["a", "b", "c"]


def test_load_input_overrides(tmp_path):
    data = common.load_input(MTX_INPUT, str(tmp_path), end_time=1234.0, drug="Methotrexate", ki=5e-6)
    assert data["settings"]["end_time"] == 1234.0 and data["project_directory"].startswith(str(tmp_path))
    dhfr = next(e for e in data["enzymes"] if e["label"] == "DHFR")
    assert dhfr["inhibition"][0]["ki"] == 5e-6
    assert common.dosing_entry(data, "methotrexate")["profile"] == "biexponential"
    with pytest.raises(ValueError, match="no enzyme lists"):
        common.load_input(MTX_INPUT, str(tmp_path), drug="aspirin", ki=1.0)


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    out = tmp_path_factory.mktemp("sweep")
    data = common.load_input(MTX_INPUT, str(out), end_time=6 * 3600.0)
    bees_object, _, model = common.build_model(data)
    return bees_object, model


def test_release_sweep_on_the_methotrexate_model(built, tmp_path):
    bees_object, model = built
    runs, rows = sweep.run_sweep(model, bees_object, "methotrexate", ["dTMP", "ghost"], 0.002, 0.1, 10000.0, 6 * 3600.0)
    by = {r["arm"]: r for r in rows}
    assert by["control"]["output_lost_pct"] == 0.0 and by["control"]["peak_free_drug_uM"] == 0.0
    assert by["free bolus"]["peak_free_drug_uM"] == pytest.approx(2.0)
    # a depot spreading the same dose lowers the peak, and the drug costs output in every arm
    assert by["liposome"]["peak_free_drug_uM"] < by["free bolus"]["peak_free_drug_uM"] / 3
    assert all(by[a]["output_lost_pct"] > 1.0 for a in ("free bolus", "liposome", "slow depot"))
    assert all(by[a]["inhibition_auc_h"] > 0 for a in ("free bolus", "liposome", "slow depot"))
    csv_path, png_path = sweep.write_outputs(runs, rows, str(tmp_path), "test", "methotrexate", "dTMP")
    assert os.path.getsize(png_path) > 5000
    assert [r["arm"] for r in csv.DictReader(open(csv_path))] == list(sweep.ARMS)
    assert runs["control"]["skipped"] == ["ghost"]


def test_target_scan_ranks_enzymes_and_restores_kinetics(built, tmp_path):
    bees_object, model = built
    before = [(r.kinetics.kcat, r.kinetics.vmax) for r in model.core_reactions]
    rows, used, _ = scan.scan(model, bees_object, ["SHMT", "TS", "DHFR", "NoSuchEnzyme"], ["dTMP"], 0.1, [1, 3, 6], 6 * 3600.0)
    assert [(r.kinetics.kcat, r.kinetics.vmax) for r in model.core_reactions] == before
    by = {r["enzyme"]: r for r in rows}
    assert by["TS"]["lost_pct_6h"] > 50 and by["SHMT"]["lost_pct_6h"] > 50
    assert by["DHFR"]["lost_pct_6h"] < by["TS"]["lost_pct_6h"]  # DHFR is in excess
    assert by["NoSuchEnzyme"]["n_reactions"] == 0 and math.isnan(by["NoSuchEnzyme"]["lost_pct_6h"])
    assert rows[-1]["enzyme"] == "NoSuchEnzyme"
    csv_path, png_path = scan.write_outputs(rows, [1, 3, 6], str(tmp_path), "test", 0.1)
    assert os.path.getsize(png_path) > 3000 and os.path.exists(csv_path)
