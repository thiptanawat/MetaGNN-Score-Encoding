"""Fixed-panel, paired evaluation of constrained exchange predictions.

The command defaults to development data. Reserved evaluation requires the shared
strict JSON protocol lock, pinning the analysis plan, inputs and sources. No target selection uses
predictions, and no bootstrap draw silently changes its target panel.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


MAGNITUDE_ARMS = [f"magnitude_g{g:.2f}" for g in (0.5, 0.8, 1.0, 1.25, 2.0)]
COMPARISON_ARMS = MAGNITUDE_ARMS + ["ordinal"]
REFERENCE_ARM = "uniform_pfba"


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def finite(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def truthy(value):
    return str(value).strip().lower() in {"true", "1", "yes"}


def context_table(rows, partition):
    out = {}
    for row in rows:
        if row["partition"] != partition:
            continue
        key = row["context_id"]
        if key in out:
            raise ValueError(f"Duplicate context_id: {key}")
        if not row.get("origin_group"):
            raise ValueError(f"Missing biological origin: {key}")
        out[key] = row
    if not out:
        raise ValueError(f"No contexts in partition {partition!r}")
    return out


def mean_observations(rows, allowed_contexts):
    """Filter context before parsing values; replicates are averaged within context.

    A duplicate replicate identifier is an input error, not an extra observation.
    Empty/nonfinite measurements are missing outcomes and retain no numeric weight.
    """
    values, seen = defaultdict(list), set()
    for row in rows:
        context = row["context_id"]
        if context not in allowed_contexts:
            continue
        key = (context, row["metabolite_id"], row["replicate"])
        if key in seen:
            raise ValueError(f"Duplicate outcome replicate: {key}")
        seen.add(key)
        number = finite(row.get("value"))
        if number is not None:
            values[(row["metabolite_id"], context)].append(number)
    return {key: float(np.mean(val)) for key, val in values.items()}


@dataclass(frozen=True)
class Pair:
    metabolite_id: str
    exchange_id: str
    context_i: str
    context_j: str
    origin_i: str
    origin_j: str
    truth: int


def build_pairs(observations, contexts, targets, tolerances):
    """Equal-weight eligible cell-line pairs from distinct biological origins.

    Within-origin pairs are excluded prospectively. Origins containing several
    distinct cell lines contribute more context pairs; resampling keeps them whole.
    """
    pairs = []
    for target in targets:
        met = target["metabolite_id"]
        tolerance = finite(tolerances.get(met))
        if tolerance is None or tolerance < 0:
            raise ValueError(f"Invalid/missing assay tolerance for {met}")
        for ci, cj in itertools.combinations(sorted(contexts), 2):
            oi, oj = contexts[ci]["origin_group"], contexts[cj]["origin_group"]
            if oi == oj:
                continue
            a, b = observations.get((met, ci)), observations.get((met, cj))
            if a is None or b is None or abs(a - b) <= tolerance:
                continue
            pairs.append(Pair(met, target["exchange_id"], ci, cj, oi, oj,
                              1 if a > b else -1))
    return pairs


def point_score(a, b, truth, tolerance):
    if abs(a - b) <= tolerance:
        return 0.5
    return float((1 if a > b else -1) == truth)


def interval_order(a, b, tolerance):
    if a[0] > b[1] + tolerance:
        return 1
    if b[0] > a[1] + tolerance:
        return -1
    return 0


def prediction_index(rows, scenario):
    out = {}
    for row in rows:
        if row["scenario"] != scenario:
            continue
        key = (row["context_id"], row["arm"], row["exchange_id"])
        if key in out:
            raise ValueError(f"Duplicate prediction: {key}")
        out[key] = row
    if not out:
        raise ValueError(f"No predictions for scenario {scenario!r}")
    return out


def prediction_values(index, context, arm, target, numerical_tolerance):
    row = index.get((context, arm, target["exchange_id"]))
    if row is None and arm == REFERENCE_ARM:
        row = index.get(("__reference__", arm, target["exchange_id"]))
    if row is None or row.get("status") != "optimal":
        return None, None, "missing_or_failed_prediction"
    factor = finite(target.get("sign_factor", 1))
    if factor is None or factor == 0:
        raise ValueError(f"Invalid sign factor for {target['metabolite_id']}")
    point, lo, hi = (finite(row.get(key)) for key in ("point", "range_lo", "range_hi"))
    point = None if point is None else point * factor
    bounds, note = None, None
    if lo is not None and hi is not None:
        if lo > hi:
            if lo - hi > numerical_tolerance:
                note = "material_range_inversion"
            else:
                lo = hi = (lo + hi) / 2
                note = "tolerance_level_range_inversion"
        if note != "material_range_inversion":
            values = (lo * factor, hi * factor)
            bounds = (min(values), max(values))
    else:
        note = "missing_or_failed_range"
    if point is None:
        note = "missing_or_failed_point"
    return point, bounds, note


def score_predictions(index, pairs, targets, arms, point_tolerance=1e-6,
                      interval_tolerance=1e-6, numerical_tolerance=1e-6):
    target_by_id = {t["metabolite_id"]: t for t in targets}
    shape = (len(arms), len(pairs))
    point = np.full(shape, 0.5)
    interval = np.full(shape, 0.5)
    point_valid = np.zeros(shape, bool)
    range_valid = np.zeros(shape, bool)
    resolved = np.zeros(shape, bool)
    notes = defaultdict(int)
    cache = {}
    for ai, arm in enumerate(arms):
        for pi, pair in enumerate(pairs):
            vals = []
            for context in (pair.context_i, pair.context_j):
                key = (context, arm, pair.metabolite_id)
                if key not in cache:
                    cache[key] = prediction_values(index, context, arm,
                                                   target_by_id[pair.metabolite_id],
                                                   numerical_tolerance)
                    if cache[key][2]:
                        notes[cache[key][2]] += 1
                vals.append(cache[key])
            a, b = vals
            if a[0] is not None and b[0] is not None:
                point_valid[ai, pi] = True
                point[ai, pi] = point_score(a[0], b[0], pair.truth, point_tolerance)
            if a[1] is not None and b[1] is not None:
                range_valid[ai, pi] = True
                order = interval_order(a[1], b[1], interval_tolerance)
                resolved[ai, pi] = order != 0
                interval[ai, pi] = 0.5 if order == 0 else float(order == pair.truth)
    return {"point": point, "interval": interval, "point_valid": point_valid,
            "range_valid": range_valid, "resolved": resolved, "notes": dict(notes)}


def aggregate_fixed(scores, pairs, target_ids, weights=None):
    """Return a macro score only when EVERY fixed target has a denominator.

    No minimum-count threshold is applied to replicated/weighted observations.
    """
    target_ids = list(target_ids)
    if len(set(target_ids)) != len(target_ids):
        raise ValueError("Duplicate fixed-panel target")
    lookup = {met: i for i, met in enumerate(target_ids)}
    scores = np.asarray(scores, dtype=float)
    if scores.ndim == 1:
        scores = scores[None, :]
    if scores.shape[1] != len(pairs):
        raise ValueError("Score/pair shape mismatch")
    weights = np.ones(len(pairs)) if weights is None else np.asarray(weights, dtype=float)
    if weights.shape != (len(pairs),) or np.any(weights < 0) or not np.isfinite(weights).all():
        raise ValueError("Invalid pair weights")
    use = np.array([i for i, p in enumerate(pairs) if p.metabolite_id in lookup], dtype=int)
    indices = np.array([lookup[pairs[i].metabolite_id] for i in use], dtype=int)
    den = np.bincount(indices, weights=weights[use], minlength=len(target_ids))
    per = np.full((scores.shape[0], len(target_ids)), np.nan)
    for ai, values in enumerate(scores):
        sums = np.bincount(indices, weights=weights[use] * values[use], minlength=len(target_ids))
        np.divide(sums, den, out=per[ai], where=den > 0)
    valid = bool(len(target_ids) and np.all(den > 0) and np.isfinite(per).all())
    macro = per.mean(axis=1) if valid else np.full(scores.shape[0], np.nan)
    return {"valid": valid, "macro": macro, "per_target": per,
            "denominators": den, "missing_targets": [target_ids[i] for i in np.where(den == 0)[0]]}


def paired_group_bootstrap(scores, pairs, target_ids, arms, origins, *,
                           resamples=2000, seed=11, maximum_invalid_fraction=0.01):
    """Paired fixed-panel origin bootstrap, with no discarded-target repair.

    Intervals using the estimable draws are explicitly conditional. If invalidity
    exceeds the declared bound, only diagnostic quantiles are returned.
    """
    origins = sorted(set(origins))
    if len(origins) < 2 or resamples < 1:
        raise ValueError("At least two origins and one resample are required")
    if not 0 <= maximum_invalid_fraction < 1:
        raise ValueError("Invalid allowed bootstrap failure fraction")
    origin_idx = {g: i for i, g in enumerate(origins)}
    ii = np.array([origin_idx[p.origin_i] for p in pairs], dtype=int)
    jj = np.array([origin_idx[p.origin_j] for p in pairs], dtype=int)
    if np.any(ii == jj):
        raise ValueError("Bootstrap expects distinct-origin pairs; within-origin pairs must be excluded")
    oi = arms.index("ordinal")
    mi = [arms.index(arm) for arm in MAGNITUDE_ARMS]
    ident = arms.index("magnitude_g1.00")
    contrast_scores = np.vstack((scores[oi] - scores[mi].mean(axis=0),
                                scores[oi] - scores[ident]))
    rng = np.random.default_rng(seed)
    draws, invalid, missing_counts = [], 0, defaultdict(int)
    for _ in range(resamples):
        pick = rng.choice(len(origins), size=len(origins), replace=True)
        mult = np.bincount(pick, minlength=len(origins))
        weights = mult[ii] * mult[jj]
        result = aggregate_fixed(np.vstack((contrast_scores, scores)), pairs, target_ids, weights)
        if not result["valid"]:
            invalid += 1
            for met in result["missing_targets"]:
                missing_counts[met] += 1
            continue
        draws.append(result["macro"].tolist())
    fraction = invalid / resamples
    quantiles = None if not draws else np.percentile(np.asarray(draws), [2.5, 97.5], axis=0).T.tolist()
    reliable = quantiles is not None and fraction <= maximum_invalid_fraction
    return {"seed": seed, "resamples": resamples, "origin_groups": len(origins),
            "valid_draws": len(draws), "invalid_draws": invalid, "invalid_fraction": fraction,
            "missing_target_counts": dict(sorted(missing_counts.items())),
            "fixed_panel_size": len(target_ids), "maximum_invalid_fraction": maximum_invalid_fraction,
            "interval_status": "conditional_on_estimability" if reliable else "unreliable_invalid_draw_rate",
            "delta_grid_interval": quantiles[0] if reliable else None,
            "delta_identity_interval": quantiles[1] if reliable else None,
            "diagnostic_valid_draw_quantiles": None if quantiles is None else quantiles[:2],
            "absolute_C_intervals": None if not reliable else {arm: quantiles[2+i] for i, arm in enumerate(arms)},
            "interpretation": "Context-resampling uncertainty conditional on the frozen pipeline and fixed targets; not model-refitting or new-lineage uncertainty."}


def _float_or_none(number):
    return float(number) if np.isfinite(number) else None


def evaluate_panel(index, observations, contexts, targets, tolerances, plan, *,
                   resamples, seed, allow_drop_unobserved=False):
    pairs = build_pairs(observations, contexts, targets, tolerances)
    ids = [t["metabolite_id"] for t in targets]
    observed = {p.metabolite_id for p in pairs}
    excluded = [met for met in ids if met not in observed]
    # Primary membership remains fixed. The broad descriptive sensitivity reports
    # an explicitly labeled effective panel defined once from observed eligibility.
    evaluated_ids = [met for met in ids if met in observed] if allow_drop_unobserved else ids
    arms = COMPARISON_ARMS + [REFERENCE_ARM]
    scored = score_predictions(index, pairs, targets, arms,
                               plan["point_tolerance"], plan["interval_tolerance"],
                               plan["numerical_tolerance"])
    point = aggregate_fixed(scored["point"], pairs, evaluated_ids)
    interval = aggregate_fixed(scored["interval"], pairs, evaluated_ids)
    required = [arms.index(a) for a in COMPARISON_ARMS]
    complete = bool(point["valid"] and np.all(scored["point_valid"][required]))
    arm_rows, target_rows = [], []
    for ai, arm in enumerate(arms):
        resolved = scored["resolved"][ai]
        nresolved = int(resolved.sum())
        common_n = len(pairs)
        arm_rows.append({"arm": arm,
                         "macro_point_operational": _float_or_none(point["macro"][ai]),
                         "macro_interval_operational": _float_or_none(interval["macro"][ai]),
                         "eligible_pairs": common_n,
                         "failed_point_pairs": int((~scored["point_valid"][ai]).sum()),
                         "failed_range_pairs": int((~scored["range_valid"][ai]).sum()),
                         "resolved_interval_pairs": nresolved,
                         "resolved_accuracy_pooled": None if not nresolved else float(scored["interval"][ai, resolved].mean()),
                         "fixed_target_count": len(evaluated_ids)})
        for ti, met in enumerate(evaluated_ids):
            target_rows.append({"arm": arm, "metabolite_id": met,
                                "eligible_pairs": int(point["denominators"][ti]),
                                "point_C": _float_or_none(point["per_target"][ai, ti]),
                                "interval_C": _float_or_none(interval["per_target"][ai, ti])})
    # Operational failure-as-tie summaries are diagnostics only when a primary
    # prediction is missing/failed. They never silently become the primary result.
    primary = None
    bootstrap = None
    if complete:
        oi, mi, ii = arms.index("ordinal"), [arms.index(a) for a in MAGNITUDE_ARMS], arms.index("magnitude_g1.00")
        primary = {"delta_grid": float(point["macro"][oi] - point["macro"][mi].mean()),
                   "delta_identity": float(point["macro"][oi] - point["macro"][ii]),
                   "absolute_macro_C": {arms[i]: float(point["macro"][i]) for i in required}}
        if resamples:
            bootstrap = paired_group_bootstrap(scored["point"], pairs, evaluated_ids, arms,
                                                [r["origin_group"] for r in contexts.values()],
                                                resamples=resamples, seed=seed,
                                                maximum_invalid_fraction=plan["maximum_invalid_fraction"])
    common = np.all(scored["point_valid"][required], axis=0)
    common_result = aggregate_fixed(scored["point"], pairs, evaluated_ids, common.astype(float))
    return {"status": "complete_predictions" if complete else "primary_not_estimable_or_prediction_failure",
            "candidate_targets": len(ids), "contributing_target_panel": evaluated_ids,
            "zero_observed_pair_targets": excluded,
            "observed_eligible_pairs": len(pairs), "contexts": len(contexts),
            "origins": len({r["origin_group"] for r in contexts.values()}),
            "primary": primary, "bootstrap": bootstrap,
            "arms": arm_rows, "per_target": target_rows, "prediction_diagnostics": scored["notes"],
            "common_success_sensitivity": {"pairs": int(common.sum()),
                                           "fixed_panel_estimable": common_result["valid"],
                                           "absolute_macro_C": {arm: _float_or_none(common_result["macro"][i]) for i, arm in enumerate(arms)}}}


def check_test_lock(partition, plan_path, lock_file, root=None):
    if partition != "test":
        return
    root = Path(root) if root is not None else Path(plan_path).parent.parent
    if Path(plan_path).resolve() != (root / "protocol/analysis_plan.json").resolve():
        raise ValueError("Reserved evaluation must use the canonical locked analysis plan")
    from .verify_inputs import verify_protocol_lock
    verify_protocol_lock(root, lock_file)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--scenario", required=True)
    parser.add_argument("--partition", choices=["development", "test"], default="development")
    parser.add_argument("--lock-file", type=Path)
    parser.add_argument("--output", type=Path, required=True, help="JSON output path; TSV siblings are also written")
    args = parser.parse_args(argv)
    plan_path = args.plan or args.root / "protocol/analysis_plan.json"
    check_test_lock(args.partition, plan_path, args.lock_file, args.root)
    plan = json.loads(plan_path.read_text())
    # Source mutation after selection is never silently accepted.
    for relative, digest in plan["input_sha256"].items():
        if sha256(args.root / relative) != digest:
            raise ValueError(f"Input differs from selected analysis plan: {relative}")
    contexts = context_table(read_tsv(args.root / "manifests/contexts.tsv"), args.partition)
    observations = mean_observations(read_tsv(args.root / "data/outcomes_long.tsv"), contexts)
    index = prediction_index(read_tsv(args.predictions), args.scenario)
    tolerances = plan["assay_tolerances"]
    seed = plan["development_seed"] if args.partition == "development" else plan["test_seed"]
    primary = evaluate_panel(index, observations, contexts, plan["primary_targets"], tolerances, plan,
                             resamples=plan["bootstrap_resamples"], seed=seed)
    broad = evaluate_panel(index, observations, contexts, plan["chemistry_targets"], tolerances, plan,
                           resamples=0, seed=seed, allow_drop_unobserved=True)
    result = {"partition": args.partition, "scenario": args.scenario,
              "analysis_plan_sha256": sha256(plan_path), "prediction_sha256": sha256(args.predictions),
              "interpretation": "Development feasibility only" if args.partition == "development" else "Locked reserved evaluation",
              "primary_fixed_panel": primary, "broader_chemistry_descriptive_sensitivity": broad,
              "pair_policy": plan["pair_policy"]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    stem = args.output.with_suffix("")
    summary_rows = [dict(panel=name, **row) for name, value in (("primary", primary), ("broader_chemistry", broad)) for row in value["arms"]]
    target_rows = [dict(panel=name, **row) for name, value in (("primary", primary), ("broader_chemistry", broad)) for row in value["per_target"]]
    write_tsv(f"{stem}_summary.tsv", summary_rows, list(summary_rows[0]))
    write_tsv(f"{stem}_per_metabolite.tsv", target_rows, list(target_rows[0]) if target_rows else ["panel", "arm", "metabolite_id"])
    print(json.dumps({"output": str(args.output), "status": primary["status"], "primary": primary["primary"],
                      "bootstrap": primary["bootstrap"]}, indent=2, allow_nan=False))
    return 0 if primary["status"] == "complete_predictions" else 2


if __name__ == "__main__":
    raise SystemExit(main())
