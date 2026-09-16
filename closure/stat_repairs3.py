#!/usr/bin/env python3
"""Third round of supporting analyses on saved outputs. No LP is solved.

Repairs two points found in the 16 September 2026 external methods review of stat_repairs.py:
  1. The supporting common-subset bootstrap let a target drop out of the panel in draws where
     none of its pairs' origins were resampled (macro() skipped zero-weight targets). The locked
     primary analysis keeps the 52-target panel fixed and never met that case (0 of 2,000 draws);
     the six-target subset did. Here every draw with a zero-weight target is recorded as invalid
     and the interval is reported conditional on the valid draws, next to the changing-panel
     value for the record.
  2. The 3,113-pair subset was the intersection of pairs whose POINT predictions differ by more
     than the tie threshold in all six encodings ("common point-nontied"), not the intersection
     of pairs whose feasible RANGES are disjoint. Both subsets are computed and named.
Also records, per target, whether the midrank-minus-grid contrast can be nonzero (the 14-target
union) and the constant-target contribution (0.5 to every arm, 0 to any contrast).
"""
import sys, json, argparse, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from stat_repairs import load, build_pairs, scores, macro, ci, TIE, ARMS, MAG
from stat_repairs2 import read_pred_file


def origin_bootstrap_fixed_panel(S_by_arm, pairs, test, origin, B, seed):
    """Origin bootstrap that keeps the target panel fixed: a draw in which any target of the panel
    receives zero total pair weight is recorded as invalid instead of being scored on fewer targets."""
    groups = sorted(set(origin.values())); gi = {g: i for i, g in enumerate(groups)}
    oi = np.array([gi[origin[c]] for c in test])
    rng = np.random.default_rng(seed)
    draws = {a: [] for a in S_by_arm}; dgrid = []; dgrid_changing = []; invalid = 0; missing = {}
    for _ in range(B):
        mult = np.bincount(rng.choice(len(groups), len(groups), replace=True), minlength=len(groups))
        w = {m: mult[oi[arr[:, 0].astype(int)]] * mult[oi[arr[:, 1].astype(int)]] for m, arr in pairs.items()}
        zero = [m for m in pairs if w[m].sum() == 0]
        est_changing = {a: macro(S_by_arm[a], pairs, w) for a in S_by_arm}
        if "ordinal" in est_changing:
            dgrid_changing.append(est_changing["ordinal"] - np.mean([est_changing[a] for a in MAG]))
        if zero:
            invalid += 1
            for m in zero: missing[m] = missing.get(m, 0) + 1
            for a in draws: draws[a].append(float("nan"))
            dgrid.append(float("nan"))
            continue
        for a in draws: draws[a].append(est_changing[a])
        dgrid.append(est_changing["ordinal"] - np.mean([est_changing[a] for a in MAG]))
    return ({a: np.array(v) for a, v in draws.items()}, np.array(dgrid), np.array(dgrid_changing),
            {"invalid_draws": invalid, "draws": B, "targets_with_zero_weight": missing})


def subset_report(name, sub, P, test, ex_of, origin, B, seed, note):
    S = {arm: scores(sub, P, arm, test, ex_of) for arm in ARMS}
    est = {arm: macro(S[arm], sub) for arm in ARMS}
    dr, dg, dg_changing, acc = origin_bootstrap_fixed_panel(S, sub, test, origin, B, seed)
    return {"subset": name, "pairs": int(sum(len(v) for v in sub.values())), "targets": len(sub),
            "pairs_per_target": {m: int(len(v)) for m, v in sub.items()},
            "macro_C": est, "delta_grid": est["ordinal"] - float(np.mean([est[x] for x in MAG])),
            "delta_identity": est["ordinal"] - est["magnitude_g1.00"],
            "delta_grid_ci95_fixed_panel_conditional": ci(dg),
            "delta_grid_ci95_changing_panel_for_the_record": ci(dg_changing),
            "absolute_ci95_fixed_panel_conditional": {arm: ci(dr[arm]) for arm in ARMS},
            "bootstrap_accounting": acc, "seed": seed, "note": note}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", type=Path, required=True); ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--B", type=int, default=2000)
    a = ap.parse_args(); a.out.mkdir(parents=True, exist_ok=True); t0 = time.time()
    test, origin, targets, tol, obs, obs_mean, P = load(a.root, "primary")
    _, R = read_pred_file(a.root / "results/reserved_v3/predictions.tsv", "primary")
    ex_of = {m: ex for m, ex in targets}
    all_arms = ARMS + ["uniform_pfba"]
    if "uniform_pfba" in P:
        ref = {ex: v for (c, ex), v in P["uniform_pfba"].items() if c == "__reference__"}
        P["uniform_pfba"] = {(c, ex): v for c in test for ex, v in ref.items()}
    report = {}

    # ---- 0. reproduction, and the fixed-panel accounting of the primary analysis itself
    pairs = build_pairs(test, origin, targets, tol, obs_mean, P, all_arms)
    S = {arm: scores(pairs, P, arm, test, ex_of) for arm in all_arms}
    saved = json.load(open(a.root / "results/reserved_v3/independent_concordance_primary.json"))
    mine = {arm: macro(S[arm], pairs) for arm in all_arms}
    diffs = {arm: abs(mine[arm] - saved["absolute_C"][arm]) for arm in all_arms}
    npairs = int(sum(len(v) for v in pairs.values()))
    report["reproduction"] = {"eligible_pairs": npairs, "targets": len(pairs), "max_abs_diff": max(diffs.values())}
    print("reproduction: pairs %d | targets %d | max |diff| %.2e" % (npairs, len(pairs), max(diffs.values())), flush=True)
    if max(diffs.values()) > 1e-9 or npairs != saved["observed_pairs"]:
        print("ABORT"); sys.exit(2)
    S6 = {arm: S[arm] for arm in ARMS}
    dr, dg, dg_ch, acc = origin_bootstrap_fixed_panel(S6, pairs, test, origin, a.B, 20260913)
    report["primary_fixed_panel_accounting"] = {"delta_grid_ci95": ci(dg), "bootstrap_accounting": acc,
        "note": "seed of stat_repairs.py; the locked evaluator's own draws are in results/reserved_v3"}
    print("primary: invalid draws %d of %d | delta_grid CI %s" % (acc["invalid_draws"], a.B, ci(dg)), flush=True)

    # ---- 1. the two common subsets
    common_pt, common_iv = {}, {}
    for m, arr in pairs.items():
        ex = ex_of[m]
        keep_pt = np.ones(len(arr), bool); keep_iv = np.ones(len(arr), bool)
        for arm in ARMS:
            d = np.array([P[arm][(test[int(x)], ex)] - P[arm][(test[int(y)], ex)] for x, y in arr[:, :2]])
            keep_pt &= np.abs(d) > TIE
            lo_i = np.array([R[arm][(test[int(x)], ex)][0] for x in arr[:, 0]]); hi_i = np.array([R[arm][(test[int(x)], ex)][1] for x in arr[:, 0]])
            lo_j = np.array([R[arm][(test[int(y)], ex)][0] for y in arr[:, 1]]); hi_j = np.array([R[arm][(test[int(y)], ex)][1] for y in arr[:, 1]])
            keep_iv &= (lo_i > hi_j + TIE) | (lo_j > hi_i + TIE)
        if keep_pt.sum() >= 1: common_pt[m] = arr[keep_pt]
        if keep_iv.sum() >= 1: common_iv[m] = arr[keep_iv]
    report["common_point_nontied_subset"] = subset_report("common point-nontied", common_pt, P, test, ex_of, origin, a.B, 20260916,
        "pairs whose point predictions differ by more than the tie threshold in all six encodings; this is the 3,113-pair subset of stat_repairs.py under its correct name")
    report["common_interval_resolved_subset"] = subset_report("common interval-resolved", common_iv, P, test, ex_of, origin, a.B, 20260916,
        "pairs whose feasible ranges are disjoint by the tie threshold in all six encodings")
    for k in ("common_point_nontied_subset", "common_interval_resolved_subset"):
        r = report[k]
        print("%s: %d pairs over %d targets | delta_grid %+.4f | fixed-panel CI %s (invalid %d/%d) | changing-panel CI %s" % (
            r["subset"], r["pairs"], r["targets"], r["delta_grid"], [round(x, 5) for x in r["delta_grid_ci95_fixed_panel_conditional"]],
            r["bootstrap_accounting"]["invalid_draws"], a.B, [round(x, 5) for x in r["delta_grid_ci95_changing_panel_for_the_record"]]), flush=True)

    # ---- 2. constant targets and the contrast: per-target contributions
    per_target = {}
    for m, arr in pairs.items():
        est = {arm: float(np.mean(S[arm][m])) for arm in ARMS}
        const_mag = all(np.ptp([P[arm][(c, ex_of[m])] for c in test]) <= TIE for arm in MAG)
        const_all = const_mag and np.ptp([P["ordinal"][(c, ex_of[m])] for c in test]) <= TIE
        per_target[m] = {"pairs": int(len(arr)), "C": est, "delta_grid_contribution": (est["ordinal"] - float(np.mean([est[x] for x in MAG]))) / len(pairs),
                         "constant_in_all_magnitude_arms": bool(const_mag), "constant_in_all_six_arms": bool(const_all)}
    n_const_mag = sum(v["constant_in_all_magnitude_arms"] for v in per_target.values())
    n_const_all = sum(v["constant_in_all_six_arms"] for v in per_target.values())
    contrib_nonzero = [m for m, v in per_target.items() if abs(v["delta_grid_contribution"]) > 1e-12]
    report["contrast_arithmetic"] = {"targets": len(pairs), "constant_in_all_magnitude_arms": n_const_mag, "constant_in_all_six_arms": n_const_all,
        "targets_where_contrast_can_be_nonzero": len(pairs) - n_const_all,
        "targets_with_nonzero_contrast_contribution_observed": len(contrib_nonzero),
        "sum_of_contributions": float(sum(v["delta_grid_contribution"] for v in per_target.values())),
        "delta_grid_locked": saved["absolute_C"]["ordinal"] - float(np.mean([saved["absolute_C"][x] for x in MAG])),
        "note": "a target constant across profiles in every arm scores 0.5 in every arm and contributes exactly 0 to any contrast; the midrank-minus-grid contrast can be nonzero only on the targets that vary in at least one of the six arms",
        "per_target": per_target}
    print("contrast arithmetic: constant in all magnitude arms %d, in all six %d, contrast can be nonzero on %d, observed nonzero on %d, sum %.6f vs locked %.6f" % (
        n_const_mag, n_const_all, len(pairs) - n_const_all, len(contrib_nonzero), report["contrast_arithmetic"]["sum_of_contributions"], report["contrast_arithmetic"]["delta_grid_locked"]), flush=True)
    report["_seconds"] = round(time.time() - t0, 1)
    json.dump(report, open(a.out / "stat_repairs3.json", "w"), indent=1)
    print("FINISHED %.0fs" % report["_seconds"], flush=True)


if __name__ == "__main__":
    main()
