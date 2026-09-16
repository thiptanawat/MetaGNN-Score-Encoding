"""Post-lock descriptive certificates from saved ranges; no CORE or LP access.

This presentation/validation script is deliberately outside the frozen scientific
source set. It does not change predictions, targets, thresholds or estimands.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def normalize_record(point, lo, hi, numerical_tolerance, sign_factor=1.0):
    """Apply the production midpoint convention before optional sign alignment."""
    point, lo, hi = float(point), float(lo), float(hi)
    if not all(math.isfinite(v) for v in (point, lo, hi)):
        raise ValueError("Nonfinite point or interval")
    if sign_factor not in (-1.0, 1.0):
        raise ValueError("Only explicit unit-magnitude orientation factors are supported")
    inverted = lo > hi
    if inverted:
        if lo - hi > numerical_tolerance:
            raise ValueError("Material range inversion")
        lo = hi = (lo + hi) / 2.0
    if point < lo - numerical_tolerance or point > hi + numerical_tolerance:
        raise ValueError("Material point-containment failure")
    if sign_factor == -1.0:
        point, lo, hi = -point, -hi, -lo
    return (point, lo, hi), inverted


def disjoint_witness(records, tolerance):
    """A separated pair of arm intervals certifies one coordinate changes."""
    high = max(records, key=lambda arm: records[arm][1])
    low = min(records, key=lambda arm: records[arm][2])
    gap = records[high][1] - records[low][2]
    if gap <= tolerance:
        return None
    return {"higher_arm": high, "lower_arm": low,
            "higher_lower_bound": records[high][1],
            "lower_upper_bound": records[low][2], "gap": gap}


def reversal_witness(first, second, tolerance):
    """Certify both possible context orderings under different arm caps."""
    positive = max(first, key=lambda arm: first[arm][1] - second[arm][2])
    negative = max(first, key=lambda arm: second[arm][1] - first[arm][2])
    positive_gap = first[positive][1] - second[positive][2]
    negative_gap = second[negative][1] - first[negative][2]
    if positive_gap <= tolerance or negative_gap <= tolerance:
        return None
    return {"first_above_second_arm": positive,
            "second_above_first_arm": negative,
            "first_above_second_gap": positive_gap,
            "second_above_first_gap": negative_gap,
            "minimum_gap": min(positive_gap, negative_gap)}


def distinct_origin_pairs(contexts, origins):
    return [(a, b) for n, a in enumerate(contexts) for b in contexts[n + 1:]
            if origins[a] != origins[b]]


def describe_panel(targets, contexts, origins, records, arms,
                   point_tolerance, interval_tolerance):
    pairs = distinct_origin_pairs(contexts, origins)
    counts = {"point_changed_cells": 0, "disjoint_interval_cells": 0,
              "disjoint_cells_without_point_change": 0,
              "point_strict_reversal_pairs": 0,
              "robust_interval_reversal_pairs": 0,
              "robust_reversals_without_point_reversal": 0}
    cell_witnesses, pair_witnesses, per_target = [], [], []
    for target in targets:
        metabolite, exchange = target["metabolite_id"], target["exchange_id"]
        row = {"metabolite_id": metabolite, "exchange_id": exchange,
               **{name: 0 for name in counts}}
        for context in contexts:
            values = {arm: records[(context, arm, exchange)] for arm in arms}
            points = [v[0] for v in values.values()]
            changed = max(points) - min(points) > point_tolerance
            witness = disjoint_witness(values, interval_tolerance)
            row["point_changed_cells"] += int(changed)
            row["disjoint_interval_cells"] += int(witness is not None)
            row["disjoint_cells_without_point_change"] += int(witness is not None and not changed)
            if witness:
                cell_witnesses.append({"context_id": context,
                                       "metabolite_id": metabolite,
                                       "exchange_id": exchange, **witness})
        for first_context, second_context in pairs:
            first = {arm: records[(first_context, arm, exchange)] for arm in arms}
            second = {arm: records[(second_context, arm, exchange)] for arm in arms}
            differences = [first[arm][0] - second[arm][0] for arm in arms]
            reversed_points = (max(differences) > point_tolerance and
                               min(differences) < -point_tolerance)
            witness = reversal_witness(first, second, interval_tolerance)
            row["point_strict_reversal_pairs"] += int(reversed_points)
            row["robust_interval_reversal_pairs"] += int(witness is not None)
            row["robust_reversals_without_point_reversal"] += int(witness is not None and not reversed_points)
            if witness:
                pair_witnesses.append({"first_context_id": first_context,
                                       "second_context_id": second_context,
                                       "metabolite_id": metabolite,
                                       "exchange_id": exchange, **witness})
        for name in counts:
            counts[name] += row[name]
        per_target.append(row)
    return {"contexts": len(contexts), "origins": len(set(origins.values())),
            "targets": len(targets), "context_target_cells": len(contexts) * len(targets),
            "distinct_origin_context_pairs": len(pairs),
            "distinct_origin_context_target_pairs": len(pairs) * len(targets),
            **counts, "per_target": per_target,
            "disjoint_cell_witnesses": cell_witnesses,
            "robust_reversal_witnesses": pair_witnesses}


def analyze(root, predictions, scenario, partition):
    root, predictions = Path(root).resolve(), Path(predictions).resolve()
    plan_path = root / "protocol/analysis_plan.json"
    contexts_path = root / "manifests/contexts.tsv"
    lock_path = root / "protocol/PROTOCOL_LOCK.json"
    plan = json.loads(plan_path.read_text())
    lock = json.loads(lock_path.read_text())
    if lock.get("status") != "locked" or not lock.get("reserved_evaluation_authorized"):
        raise ValueError("An authorized, completed protocol lock is required")
    if sha256(plan_path) != lock["analysis_plan_sha256"]:
        raise ValueError("Analysis plan differs from the locked plan")
    if sha256(contexts_path) != plan["input_sha256"]["manifests/contexts.tsv"]:
        raise ValueError("Origin metadata differs from the frozen analysis input")
    point_tolerance, interval_tolerance = plan["point_tolerance"], plan["interval_tolerance"]
    numerical_tolerance = plan["numerical_tolerance"]
    arms = [v for v in plan["comparison_arms"] if v.startswith("magnitude_g")]
    if len(arms) != 5 or len(set(arms)) != 5:
        raise ValueError("Expected the five frozen magnitude arms")
    with contexts_path.open() as handle:
        metadata = [r for r in csv.DictReader(handle, delimiter="\t") if r["partition"] == partition]
    contexts = sorted(row["context_id"] for row in metadata)
    origins = {row["context_id"]: row["origin_group"] for row in metadata}
    if not contexts or len(contexts) != len(origins):
        raise ValueError("Missing or duplicate partition context metadata")
    chemistry = plan["chemistry_targets"]
    exchanges = {r["exchange_id"]: r for r in chemistry}
    if len(exchanges) != len(chemistry):
        raise ValueError("Expected one chemical target per exchange")
    values, inversions = {}, []
    with predictions.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["scenario"] != scenario or row["arm"] not in arms:
                continue
            context, exchange = row["context_id"], row["exchange_id"]
            if context not in origins or exchange not in exchanges:
                raise ValueError("Unexpected context or exchange in the requested magnitude scenario")
            key = context, row["arm"], exchange
            if key in values:
                raise ValueError(f"Duplicate prediction: {key}")
            if row["status"] != "optimal":
                raise ValueError(f"Nonoptimal saved prediction: {key}")
            value, inverted = normalize_record(row["point"], row["range_lo"], row["range_hi"],
                                               numerical_tolerance,
                                               float(exchanges[exchange]["sign_factor"]))
            values[key] = value
            if inverted:
                inversions.append({"context_id": context, "arm": row["arm"], "exchange_id": exchange})
    expected = {(c, a, e) for c in contexts for a in arms for e in exchanges}
    if set(values) != expected:
        raise ValueError(f"Incomplete fixed comparison grid: missing {len(expected-set(values))} records")
    created_at = datetime.now(timezone.utc).isoformat()
    result = {
        "schema_version": 1, "analysis_status": "post_lock_post_primary_inspection_descriptive",
        "created_at_utc": created_at, "partition": partition, "scenario": scenario,
        "primary_endpoint_changed": False, "core_outcome_table_read": False,
        "new_lp_solves": 0, "intervals_reoptimized": False,
        "threshold_source": "exact locked analysis plan; no audit-specific tuning",
        "point_tolerance": point_tolerance, "interval_tolerance": interval_tolerance,
        "numerical_tolerance": numerical_tolerance, "arms": arms,
        "normalized_tolerance_level_inversion_records": inversions,
        "input_sha256": {"predictions.tsv": sha256(predictions),
                         "protocol/analysis_plan.json": sha256(plan_path),
                         "manifests/contexts.tsv": sha256(contexts_path),
                         "protocol/PROTOCOL_LOCK.json": sha256(lock_path)},
        "script_sha256": sha256(__file__),
        "input_paths": {"predictions": str(predictions), "root": str(root)},
        "interpretation": [
            "Positive gaps are numerically qualified coordinate-wise certificates under each arm's own saved near-optimal objective caps, not mathematically exact certificates.",
            "Disjoint same-context intervals show that the corresponding coordinate differs across at least two encodings throughout their reported ranges.",
            "A robust reversal requires opposite between-context orderings separated by more than the fixed tolerance under two different encodings.",
            "Overlap is inconclusive: it does not establish a common jointly feasible vector, absence of an encoding effect, or biological equality.",
            "Ranges are optimization outputs conditional on model, constraints and caps; they are not confidence intervals or in vivo flux bounds.",
            "Denominators include all possible different-origin context pairs, regardless of CORE separability; they are not independent sample sizes or accuracy denominators.",
            "This audit is supplementary and was specified after inspecting reserved primary findings; no CORE values, new LPs or modified primary estimands enter it."
        ]}
    for name, targets in [("primary_fixed_panel", plan["primary_targets"]), ("broader_chemistry", chemistry)]:
        result[name] = describe_panel(targets, contexts, origins, values, arms,
                                      point_tolerance, interval_tolerance)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--scenario", choices=["primary", "half_serum", "lower_task"], required=True)
    parser.add_argument("--partition", choices=["development", "test"], default="test")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.root, args.predictions, args.scenario, args.partition)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"output": str(args.output), "created_at_utc": result["created_at_utc"],
                      "primary": {k: v for k, v in result["primary_fixed_panel"].items()
                                  if not isinstance(v, list)}}, indent=2))


if __name__ == "__main__":
    main()
