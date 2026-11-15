"""Known-answer checks of second/summarize.py on a tiny synthetic second root.

The root is built in a temporary directory: a contexts table with two partitions, an analysis plan
with three targets, two arm names and the point tolerance, and a predictions table whose points are
chosen by hand. No model, no solver and no measurement is involved.
"""
import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "summarize.py"
ARMS = ["magnitude_g1.00", "ordinal"]
TOL = 1e-5
UNIFORM = {"EX_same": 5.0, "EX_away": 5.0, "EX_vary": 1.0}
TEST_CONTEXTS = ["t1", "t2", "t3"]
# Points of the encodings in the reserved ("test") partition, by hand.
#   EX_same  constant at the uniform value; one profile differs from it by 4e-6, inside the tolerance
#   EX_away  constant at 3.0, away from the uniform value of 5.0
#   EX_vary  1, 2, 3 across the three profiles, so it varies and the first profile sits at its uniform value
POINTS = {
    "EX_same": {"t1": 5.0, "t2": 5.0 + 4e-6, "t3": 5.0},
    "EX_away": {"t1": 3.0, "t2": 3.0, "t3": 3.0},
    "EX_vary": {"t1": 1.0, "t2": 2.0, "t3": 3.0},
}


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def prediction_rows(drop=None, drop_uniform=None):
    rows = []
    for target, by_context in POINTS.items():
        for context, point in by_context.items():
            for arm in ARMS:
                rows.append(dict(context_id=context, arm=arm, scenario="primary", exchange_id=target, status="optimal", point=point))
    # A development profile, a second scenario and a failed arm must not enter the test-partition counts.
    for target in POINTS:
        for arm in ARMS:
            rows.append(dict(context_id="d1", arm=arm, scenario="primary", exchange_id=target, status="optimal", point=99.0))
            rows.append(dict(context_id="t1", arm=arm, scenario="half_serum", exchange_id=target, status="optimal", point=77.0))
            rows.append(dict(context_id="t2", arm=arm, scenario="primary", exchange_id=target, status="failed", point=55.0))
    for target, value in UNIFORM.items():
        rows.append(dict(context_id="", arm="uniform_pfba", scenario="primary", exchange_id=target, status="optimal", point=value))
    if drop is not None:
        rows = [r for r in rows if not (r["context_id"] == drop[0] and r["arm"] == drop[1] and r["exchange_id"] == drop[2]
                                        and r["scenario"] == "primary" and r["status"] == "optimal")]
    if drop_uniform is not None:
        rows = [r for r in rows if not (r["arm"] == "uniform_pfba" and r["exchange_id"] == drop_uniform)]
    return rows


def build_root(tmp_path, **drops):
    root = tmp_path / "second_root"
    (root / "manifests").mkdir(parents=True)
    (root / "protocol").mkdir()
    (root / "results").mkdir()
    contexts = [dict(context_id="d1", partition="development")] + [dict(context_id=c, partition="test") for c in TEST_CONTEXTS]
    write_tsv(root / "manifests/contexts.tsv", contexts, ["context_id", "partition"])
    plan = dict(point_tolerance=TOL, comparison_arms=ARMS, primary_targets=[dict(exchange_id=t) for t in UNIFORM])
    (root / "protocol/analysis_plan.json").write_text(json.dumps(plan))
    write_tsv(root / "results/predictions.tsv", prediction_rows(**drops),
              ["context_id", "arm", "scenario", "exchange_id", "status", "point"])
    return root


def run(root, partition="test"):
    return subprocess.run([sys.executable, str(SCRIPT), "--root", str(root), "--partition", partition,
                           "--predictions", "results/predictions.tsv", "--output", "results/summary.json"],
                          capture_output=True, text=True)


def test_counts_on_hand_built_cases(tmp_path):
    root = build_root(tmp_path)
    done = run(root)
    assert done.returncode == 0, done.stderr
    out = json.loads((root / "results/summary.json").read_text())
    assert out["profiles"] == 3 and out["targets"] == 3 and out["partition"] == "test"
    # EX_same and EX_away are constant, EX_vary is not, under both encodings
    assert out["constant_across_profiles"] == {arm: 2 for arm in ARMS}
    assert out["constant_under_every_encoding"] == 2
    # only EX_same sits at the uniform value in every profile (its 4e-6 offset is inside the tolerance)
    assert out["constant_and_at_uniform_value"] == {arm: 1 for arm in ARMS}
    assert out["constant_under_every_encoding_and_at_uniform_value"] == 1
    assert out["constant_under_every_encoding_not_at_uniform_value"] == ["EX_away"]
    assert out["varying_targets_under_every_encoding"] == ["EX_vary"]
    # moved from the uniform value: t1 has EX_away only, t2 and t3 have EX_away and EX_vary
    for arm in ARMS:
        assert out["targets_moved_from_uniform_value_per_profile"][arm] == {"median": 2, "min": 1, "max": 2}


def test_other_partition_scenario_and_failed_rows_are_ignored(tmp_path):
    root = build_root(tmp_path)
    assert run(root, partition="development").returncode == 0
    out = json.loads((root / "results/summary.json").read_text())
    # one profile (d1) with every point at 99: every target is constant, none is at the uniform value
    assert out["profiles"] == 1
    assert out["constant_under_every_encoding"] == 3
    assert out["constant_under_every_encoding_and_at_uniform_value"] == 0
    assert out["targets_moved_from_uniform_value_per_profile"]["ordinal"] == {"median": 3, "min": 3, "max": 3}


def test_a_difference_just_above_the_tolerance_is_not_constant(tmp_path, monkeypatch):
    monkeypatch.setitem(POINTS, "EX_same", {"t1": 5.0, "t2": 5.0 + 2e-5, "t3": 5.0})
    root = build_root(tmp_path)
    assert run(root).returncode == 0
    out = json.loads((root / "results/summary.json").read_text())
    assert out["constant_across_profiles"] == {arm: 1 for arm in ARMS}
    assert out["varying_targets_under_every_encoding"] == ["EX_same", "EX_vary"]


def test_a_missing_prediction_stops_with_an_error(tmp_path):
    root = build_root(tmp_path, drop=("t2", "ordinal", "EX_vary"))
    done = run(root)
    assert done.returncode != 0
    assert "predictions incomplete" in done.stderr
    assert "1 cells of the encodings missing" in done.stderr
    assert not (root / "results/summary.json").exists()


def test_a_missing_uniform_value_stops_with_an_error(tmp_path):
    root = build_root(tmp_path, drop_uniform="EX_away")
    done = run(root)
    assert done.returncode != 0
    assert "1 of the uniform cost" in done.stderr
    assert not (root / "results/summary.json").exists()
