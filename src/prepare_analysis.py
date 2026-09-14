"""Select a fixed, observed-only development estimability panel before predictions.

For each chemistry-eligible target, calculate exactly the probability that an
N-origin empirical bootstrap of the N development origins has no eligible pair.
Retain it if that risk is <= alpha / M_chemistry. A union bound then controls the
development full-panel invalid-draw risk by alpha. This is conditional on the
observed development pair graphs, not a promise about future test precision.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from .evaluation import (COMPARISON_ARMS, build_pairs, context_table, finite,
                         mean_observations, read_tsv, sha256, truthy, write_tsv)


def support_distribution(n):
    """Exact probability of each origin support under N draws from N origins.

    Every given support of size k has k! S(N,k) / N**N probability. Compute its
    numerator by inclusion-exclusion using integers. Exhaustive enumeration is
    bounded to 18 development origins, preventing an accidental exponential run.
    """
    if not 2 <= n <= 18:
        raise ValueError("Exact support enumeration requires 2..18 development origins; amend the method otherwise")
    supports = np.arange(1 << n, dtype=np.uint64)
    counts = np.array([int(s).bit_count() for s in supports])
    onto = [sum((-1) ** j * math.comb(k, j) * (k - j) ** n for j in range(k + 1)) for k in range(n + 1)]
    probabilities = np.array([onto[k] / n ** n for k in counts], dtype=float)
    if not np.isclose(probabilities.sum(), 1, atol=1e-12):
        raise ArithmeticError("Bootstrap support probabilities do not sum to one")
    return supports, probabilities


def estimability_risks(pairs, target_ids, origins):
    origins = sorted(set(origins))
    supports, probabilities = support_distribution(len(origins))
    origin_idx = {g: i for i, g in enumerate(origins)}
    grouped = {met: set() for met in target_ids}
    for pair in pairs:
        if pair.origin_i == pair.origin_j:
            raise ValueError("Within-origin pairs cannot define the selected estimand")
        grouped[pair.metabolite_id].add(tuple(sorted((origin_idx[pair.origin_i], origin_idx[pair.origin_j]))))
    risks, invalid_masks = {}, {}
    for met in target_ids:
        valid = np.zeros(len(supports), dtype=bool)
        for a, b in grouped[met]:
            valid |= ((supports & (1 << a)) != 0) & ((supports & (1 << b)) != 0)
        invalid_masks[met] = ~valid
        risks[met] = {"invalid_probability": float(probabilities[~valid].sum()),
                      "eligible_origin_edges": len(grouped[met]),
                      "origins_with_an_eligible_edge": len({x for edge in grouped[met] for x in edge})}
    return risks, invalid_masks, probabilities


def prepare(root, *, alpha=0.01, point_tolerance=1e-5, interval_tolerance=1e-5,
            numerical_tolerance=1e-6, bootstrap_resamples=2000):
    root = Path(root)
    if not 0 < alpha < 1:
        raise ValueError("Invalid maximum panel risk")
    if any(not math.isfinite(t) or t < 0 for t in (point_tolerance, interval_tolerance, numerical_tolerance)):
        raise ValueError("Numerical tolerances must be finite and nonnegative")
    if bootstrap_resamples < 1:
        raise ValueError("At least one bootstrap resample is required")
    paths = ["manifests/contexts.tsv", "manifests/metabolite_map.tsv",
             "data/outcomes_long.tsv", "data/assay_noise_development.tsv"]
    contexts = context_table(read_tsv(root / paths[0]), "development")
    observations = mean_observations(read_tsv(root / paths[2]), contexts)
    noise_rows = read_tsv(root / paths[3])
    tolerances = {}
    for row in noise_rows:
        met = row["metabolite_id"]
        if met in tolerances:
            raise ValueError(f"Duplicate noise rule: {met}")
        value = finite(row["tolerance"])
        if value is not None and value >= 0:
            tolerances[met] = value
    chemistry, exclusions, seen_met, seen_exchange = [], [], set(), set()
    for row in read_tsv(root / paths[1]):
        if row["status"] != "mapped" or not truthy(row["primary_eligible"]):
            continue
        met, exchange = row["metabolite_id"], row["exchange_id"]
        if met in seen_met or exchange in seen_exchange:
            raise ValueError(f"Primary map must be one-to-one: {met}, {exchange}")
        seen_met.add(met)
        seen_exchange.add(exchange)
        factor = finite(row.get("sign_factor", 1))
        if factor is None or factor == 0:
            raise ValueError(f"Invalid sign_factor for {met}")
        chemistry.append({"metabolite_id": met, "exchange_id": exchange,
                          "sign_factor": factor, "source_label": row.get("source_label", met)})
        if met not in tolerances:
            exclusions.append({"metabolite_id": met, "reason": "no_valid_development_assay_tolerance"})
    chemistry.sort(key=lambda t: t["metabolite_id"])
    if not chemistry:
        raise ValueError("No independently chemistry-eligible targets")
    valid_noise = [t for t in chemistry if t["metabolite_id"] in tolerances]
    pairs = build_pairs(observations, contexts, valid_noise, tolerances)
    origins = sorted({row["origin_group"] for row in contexts.values()})
    risks, invalid_masks, probabilities = estimability_risks(pairs, [t["metabolite_id"] for t in valid_noise], origins)
    cutoff = alpha / len(chemistry)
    selected = [t for t in valid_noise if risks[t["metabolite_id"]]["invalid_probability"] <= cutoff]
    if not selected:
        raise ValueError("No target meets the prespecified observed-only estimability rule; stop and amend transparently")
    joint_invalid = np.zeros(len(probabilities), bool)
    for target in selected:
        joint_invalid |= invalid_masks[target["metabolite_id"]]
    union_bound = sum(risks[t["metabolite_id"]]["invalid_probability"] for t in selected)
    exact_joint = float(probabilities[joint_invalid].sum())
    if exact_joint > union_bound + 1e-12 or union_bound > alpha + 1e-12:
        raise ArithmeticError("Panel invalidity bound failed")
    counts = {met: 0 for met in seen_met}
    for pair in pairs:
        counts[pair.metabolite_id] += 1
    selected_ids = {t["metabolite_id"] for t in selected}
    ledger = []
    for target in chemistry:
        met = target["metabolite_id"]
        risk = risks.get(met, {})
        ledger.append({**target, "eligible_development_context_pairs": counts[met],
                       **risk, "risk_cutoff": cutoff, "primary_selected": met in selected_ids,
                       "reason": "observed_only_estimability_pass" if met in selected_ids else
                                 ("no_valid_development_assay_tolerance" if met not in tolerances else "observed_pair_graph_not_estimable_at_prespecified_risk")})
    plan = {"version": 1, "status": "development_selection_unlocked",
            "selection_inputs": "Development observed values and assay tolerance; no predictions or reserved performance",
            "input_sha256": {path: sha256(root / path) for path in paths},
            "development_contexts": sorted(contexts), "development_origins": origins,
            "chemistry_candidate_count": len(chemistry),
            "chemistry_targets": valid_noise, "chemistry_without_noise": exclusions,
            "primary_targets": selected, "assay_tolerances": tolerances,
            "maximum_invalid_fraction": alpha, "per_target_invalid_probability_cutoff": cutoff,
            "development_exact_joint_invalid_probability": exact_joint,
            "development_union_bound_invalid_probability": union_bound,
            "selection_rule": "Keep targets with exact empirical N-origin bootstrap P(no eligible distinct-origin pair) <= alpha / M_chemistry",
            "risk_scope": "Exact conditional risk for the observed development graphs; not test-set invalidity or inferential coverage",
            "pair_policy": "Equal weight per eligible distinct-context pair from different biological origins; within-origin pairs excluded. Origins with multiple distinct cell lines contribute more pair rows; bootstrap carries each origin whole.",
            "point_tolerance": point_tolerance, "interval_tolerance": interval_tolerance,
            "numerical_tolerance": numerical_tolerance, "comparison_arms": COMPARISON_ARMS,
            "bootstrap_resamples": bootstrap_resamples, "development_seed": 11, "test_seed": 20260913,
            "interval_policy": "Fixed panel; invalid draws reported, never redraw/drop targets. Conditional interval suppressed when invalid fraction exceeds alpha.",
            "within_origin_amendment": "Primary comparisons now exclude within-origin pairs; original protocol left that rule unspecified. This was chosen before new prediction evaluation."}
    return plan, ledger, pairs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--alpha", type=float, default=0.01)
    parser.add_argument("--point-tolerance", type=float, default=1e-5)
    parser.add_argument("--interval-tolerance", type=float, default=1e-5)
    parser.add_argument("--numerical-tolerance", type=float, default=1e-6)
    parser.add_argument("--bootstrap-resamples", type=int, default=2000)
    args = parser.parse_args(argv)
    plan, ledger, pairs = prepare(args.root, alpha=args.alpha, point_tolerance=args.point_tolerance,
                                  interval_tolerance=args.interval_tolerance,
                                  numerical_tolerance=args.numerical_tolerance,
                                  bootstrap_resamples=args.bootstrap_resamples)
    output = args.root / "protocol/analysis_plan.json"
    if any((args.root / f"protocol/PROTOCOL_LOCK.{extension}").exists() for extension in ("md", "json")):
        raise ValueError("Refusing to overwrite analysis selection after a protocol lock exists")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(plan, indent=2, allow_nan=False) + "\n")
    fields = list(dict.fromkeys(k for row in ledger for k in row))
    write_tsv(args.root / "manifests/analysis_panel.tsv", ledger, fields)
    write_tsv(args.root / "manifests/development_observed_pairs.tsv", [asdict_pair(p) for p in pairs],
              ["metabolite_id", "exchange_id", "context_i", "context_j", "origin_i", "origin_j", "truth"])
    print(json.dumps({"plan": str(output), "sha256": sha256(output),
                      "chemistry_candidates": plan["chemistry_candidate_count"],
                      "primary_targets": len(plan["primary_targets"]),
                      "exact_development_invalid_probability": plan["development_exact_joint_invalid_probability"]}, indent=2))
    return 0


def asdict_pair(pair):
    return {key: getattr(pair, key) for key in Pair_fields}


Pair_fields = ["metabolite_id", "exchange_id", "context_i", "context_j", "origin_i", "origin_j", "truth"]


if __name__ == "__main__":
    raise SystemExit(main())
