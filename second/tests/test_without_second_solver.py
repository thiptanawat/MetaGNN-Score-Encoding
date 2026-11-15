"""Checks of results/second/without_second_solver.json against the deposited records it was computed from.

second/without_second_solver.py needs the second study root, which a checkout does not hold. These tests
need the checkout only. They recompute from the deposited records every quantity of the file that does
not read a measurement (the profiles removed, the range stages, the targets constant across profiles,
the fixed coordinates of the cap sweep) and compare the values the file gives for all 47 profiles with
the deposited outputs of the reserved run. The concordance of the 43 profiles reads the measurements
through the locked evaluator and is not recomputed here.
"""
import csv
import glob
import hashlib
import json
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
R2 = ROOT / "results" / "second"
W = json.loads((R2 / "without_second_solver.json").read_text())
FULL, LESS = W["all_profiles"], W["without_removed_profiles"]
PASS2 = json.loads((R2 / "reserved_v3" / "second_solver_arms.json").read_text())
TOL = 1e-5
ENCODINGS = ["magnitude_g0.50", "magnitude_g0.80", "magnitude_g1.00", "magnitude_g1.25", "magnitude_g2.00", "ordinal"]
CAPS = ["exact", "rel_1pct", "rel_5pct", "rel_20pct"]


def reserved_profiles():
    rows = csv.DictReader(open(R2 / "reserved_v3" / "predictions.tsv"), delimiter="\t")
    return sorted({r["context_id"] for r in rows if r["scenario"] == "primary" and r["arm"] in ENCODINGS})


def kept_profiles():
    return [c for c in reserved_profiles() if c not in W["profiles_removed"]]


def test_the_profiles_removed_are_those_of_the_second_solver():
    assert W["profiles_removed"] == sorted({x["context_id"] for x in PASS2["solved_in_pass2"]})
    assert sorted((x["context_id"], x["arm"]) for x in W["arms_completed_by_the_second_solver"]) == sorted((x["context_id"], x["arm"]) for x in PASS2["solved_in_pass2"])
    assert len(W["profiles_removed"]) == 4 and FULL["profiles"] == len(reserved_profiles()) == 47 and LESS["profiles"] == len(kept_profiles()) == 43


def test_the_file_names_the_plan_and_the_predictions_it_read():
    lock = json.loads((R2 / "protocol" / "PROTOCOL_LOCK.json").read_text())
    assert W["second_reconstruction_plan_sha256"] == lock["second_reconstruction_plan_sha256"]
    assert W["predictions_sha256"] == hashlib.sha256((R2 / "reserved_v3" / "predictions.tsv").read_bytes()).hexdigest()
    assert (ROOT / W["plan"]).is_file()
    assert all(W["all_profiles_reproduce_the_deposited_values"].values()) and len(W["all_profiles_reproduce_the_deposited_values"]) == 5


def test_the_values_for_all_profiles_are_the_deposited_ones():
    dep = json.loads((R2 / "reserved_v3" / "evaluation_primary.json").read_text())["primary_fixed_panel"]
    e = FULL["evaluation"]
    assert e["eligible_pairs"] == dep["observed_eligible_pairs"] and e["primary"] == dep["primary"]
    assert e["delta_grid_interval"] == dep["bootstrap"]["delta_grid_interval"] and e["delta_identity_interval"] == dep["bootstrap"]["delta_identity_interval"]
    assert (e["seed"], e["resamples"], e["valid_draws"]) == (dep["bootstrap"]["seed"], dep["bootstrap"]["resamples"], dep["bootstrap"]["valid_draws"])
    tr = json.loads((R2 / "reserved_v3" / "transform_primary.json").read_text())["primary_fixed_panel"]
    assert all(FULL["transform"][k] == tr[k] for k in FULL["transform"])
    summary = json.loads((R2 / "second_summary.json").read_text())
    for k in ("constant_under_every_encoding", "constant_under_every_encoding_and_at_uniform_value", "targets_moved_from_uniform_value_per_profile"):
        assert FULL["constant_targets"][k] == summary[k]
    sweep = json.loads((R2 / "closure" / "cap_sweep" / "summary.json").read_text())
    conc = json.loads((R2 / "closure" / "cap_sweep" / "cap_concordance.json").read_text())
    for arm in ("magnitude_g1.00", "ordinal"):
        for cap in CAPS:
            got = FULL["cap_sweep"][arm][cap]
            assert abs(got["mean_fixed_of_52"] - sweep[arm][cap]["mean_fixed_of_52"]) < 1e-12
            assert got["resolved_pairs"] == conc[arm][cap]["resolved_pairs"] and abs(got["resolved_share"] - conc[arm][cap]["resolved_share"]) < 1e-15


def test_the_same_evaluation_settings_for_the_kept_profiles():
    a, b = FULL["evaluation"], LESS["evaluation"]
    assert (b["seed"], b["resamples"]) == (a["seed"], a["resamples"]) and b["status"] == "complete_predictions"
    assert b["valid_draws"] + b["invalid_draws"] == b["resamples"] and b["origin_groups"] == LESS["origins"] < FULL["origins"]
    assert b["eligible_pairs"] < a["eligible_pairs"] and set(b["primary"]["absolute_macro_C"]) == set(ENCODINGS)
    for key in ("delta_grid", "delta_identity"):
        low, high = b[key + "_interval"]
        assert low <= b["primary"][key] <= high


def test_the_range_stages_of_the_kept_profiles():
    keep = set(kept_profiles()); rows = json.loads((R2 / "mechanism" / "reserved_range" / "summary.json").read_text())["summary_by_profile_arm"]
    for arm in ENCODINGS:
        mine = [r for r in rows if r["arm"] == arm and r["context"] in keep]
        got = LESS["range_stages"][arm]
        assert len(mine) == got["rows"] == 43
        assert sum(r["targets"] - r["B0_fixed"] for r in mine) == got["open"]
        assert sum(r["first_fixed_B1"] for r in mine) == got["first_fixed_by_cost_cap"] and sum(r["B0_fixed"] for r in mine) == got["fixed_by_network"]
        assert abs(got["share_first_fixed_by_cost_cap"] - got["first_fixed_by_cost_cap"] / got["open"]) < 1e-15


def test_the_constant_targets_of_the_kept_profiles():
    primary = sorted(r["exchange_id"] for r in csv.DictReader(open(ROOT / "second" / "target_map.tsv"), delimiter="\t") if r["primary"] == "True")
    assert len(primary) == LESS["constant_targets"]["targets"] == 52
    keep = kept_profiles(); point, uniform = {}, {}
    for r in csv.DictReader(open(R2 / "reserved_v3" / "predictions.tsv"), delimiter="\t"):
        if r["scenario"] != "primary" or r["status"] != "optimal" or r["exchange_id"] not in primary: continue
        if r["arm"] == "uniform_pfba": uniform[r["exchange_id"]] = float(r["point"])
        else: point[(r["context_id"], r["arm"], r["exchange_id"])] = float(r["point"])
    constant = [t for t in primary if all(max(point[(c, m, t)] for c in keep) - min(point[(c, m, t)] for c in keep) <= TOL for m in ENCODINGS)]
    at_uniform = [t for t in constant if all(abs(point[(c, m, t)] - uniform[t]) <= TOL for c in keep for m in ENCODINGS)]
    got = LESS["constant_targets"]
    assert len(constant) == got["constant_under_every_encoding"] and len(at_uniform) == got["constant_under_every_encoding_and_at_uniform_value"]
    for m in ENCODINGS:
        moved = [sum(abs(point[(c, m, t)] - uniform[t]) > TOL for t in primary) for c in keep]
        assert got["targets_moved_from_uniform_value_per_profile"][m] == dict(median=statistics.median(moved), min=min(moved), max=max(moved))


def test_the_fixed_coordinates_of_the_cap_sweep_for_the_kept_profiles():
    keep = set(kept_profiles()); jobs = [json.loads(Path(f).read_text()) for f in glob.glob(str(R2 / "closure" / "cap_sweep" / "jobs" / "*.json"))]
    for arm in ("magnitude_g1.00", "ordinal"):
        mine = [j for j in jobs if j["arm"] == arm and j["context"] in keep and j["status"] == "complete"]
        assert len(mine) == LESS["cap_sweep"][arm]["profiles"] == 43
        for cap in CAPS:
            got = LESS["cap_sweep"][arm][cap]
            assert abs(got["mean_fixed_of_52"] - sum(j["caps"][cap]["fixed"] for j in mine) / len(mine)) < 1e-12
            assert abs(got["resolved_share"] - got["resolved_pairs"] / LESS["evaluation"]["eligible_pairs"]) < 1e-15


def test_the_readings_of_the_plan():
    for part in (FULL, LESS):
        e, r = part["evaluation"], part["readings"]
        includes_zero = e["delta_grid_interval"][0] < 0 < e["delta_grid_interval"][1] and e["delta_identity_interval"][0] < 0 < e["delta_identity_interval"][1]
        assert r["encoding"] == includes_zero
        half_open = all(a["first_fixed_by_cost_cap"] > 0.5 * a["open"] for a in part["range_stages"].values())
        c = part["constant_targets"]
        assert r["attribution"] == (half_open and c["constant_under_every_encoding_and_at_uniform_value"] > 0.5 * c["constant_under_every_encoding"])
        k = part["cap_sweep"]
        assert r["allowance"] == all(k[a]["rel_1pct"]["resolved_share"] < k[a]["exact"]["resolved_share"] and k[a]["rel_1pct"]["mean_fixed_of_52"] < k[a]["exact"]["mean_fixed_of_52"]
                                     for a in ("magnitude_g1.00", "ordinal"))
