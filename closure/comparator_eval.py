#!/usr/bin/env python3
"""Score the two comparator formulations (bound-based encoding; RIPTiDe) with the study's own conventions.

Everything is computed on the locked pair set: the 46,447 assay-separable, different-origin pairs over the
52 fixed targets that the reimplemented evaluator reproduces to the last digit. For each comparator arm:
  point concordance (macro over targets; tie at 1e-5) with a 2,000-draw origin bootstrap interval;
  interval concordance with the study's interval rule (resolved when disjoint by the tie tolerance);
  across-profile constancy (targets whose point span across profiles is at most 1e-5; targets zero everywhere);
  exploration against the study's shared B0 width (spread of interval midpoints across profiles / B0 width)
  for the 49 targets the network left free, so the comparators sit on the same scale as Table 8.
A pair is dropped for a comparator only when that comparator has no prediction for one of its contexts
(RIPTiDe failed on one profile); the number of pairs actually scored is reported. No LP is solved here.
"""
import sys, csv, json, argparse
from collections import defaultdict
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from stat_repairs import load, build_pairs, TIE, MAG, ARMS, origin_bootstrap, macro, ci

EPS = 1e-5


def read_pred(path, scenario="primary"):
    P = defaultdict(dict); R = defaultdict(dict); status = defaultdict(dict)
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["status"] != "optimal" or r["point"] == "": continue
        if r.get("scenario", scenario) != scenario: continue
        if r["context_id"] == "__reference__": continue
        P[r["arm"]][(r["context_id"], r["exchange_id"])] = float(r["point"])
        R[r["arm"]][(r["context_id"], r["exchange_id"])] = [float(r["range_lo"]), float(r["range_hi"])]
    return P, R


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", type=Path, required=True); ap.add_argument("--pred", type=Path, required=True)
    ap.add_argument("--b0", type=Path, required=True, help="range_location_reserved_stages.tsv from the mechanism analysis (shared B0 per target)")
    ap.add_argument("--out", type=Path, required=True); ap.add_argument("--label", required=True); ap.add_argument("--B", type=int, default=2000)
    a = ap.parse_args(); a.out.mkdir(parents=True, exist_ok=True)
    test, origin, targets, tol, obs, obs_mean, P0 = load(a.root, "primary")
    ex_of = {m: ex for m, ex in targets}; primary_ex = [ex for m, ex in targets]
    ref = {ex: v for (c, ex), v in P0["uniform_pfba"].items() if c == "__reference__"}
    P0["uniform_pfba"] = {(c, ex): v for c in test for ex, v in ref.items()}
    pairs = build_pairs(test, origin, targets, tol, obs_mean, P0, ARMS + ["uniform_pfba"])
    assert sum(len(v) for v in pairs.values()) == 46447
    b0 = {}
    for r in csv.DictReader(open(a.b0), delimiter="\t"):
        if r["arm"] == "magnitude_g1.00": b0[r["exchange_id"]] = float(r["B0_width"])
    free = [ex for ex in primary_ex if b0.get(ex, 0.0) > EPS]
    P, R = read_pred(a.pred)
    out = {"label": a.label, "pairs_locked": 46447, "free_targets_study_B0": len(free), "arms": {}}
    S_by_arm = {}; pairs_by_arm = {}
    for arm in sorted(P):
        ctx_have = sorted({c for (c, ex) in P[arm]})
        # pairs scorable for this arm
        sub = {}
        for m, arr in pairs.items():
            ex = ex_of[m]
            keep = [row for row in arr if (test[int(row[0])], ex) in P[arm] and (test[int(row[1])], ex) in P[arm]]
            if keep: sub[m] = np.array(keep)
        n_pairs = int(sum(len(v) for v in sub.values()))
        S = {}; res = corr = 0; per_int = []; n_int_pairs = 0
        for m, arr in sub.items():
            ex = ex_of[m]
            pi = np.array([P[arm][(test[int(x)], ex)] for x in arr[:, 0]]); pj = np.array([P[arm][(test[int(y)], ex)] for y in arr[:, 1]])
            d = pi - pj; S[m] = np.where(np.abs(d) <= TIE, 0.5, (np.sign(d) == arr[:, 2]).astype(float))
            s = []
            for x, y, truth in arr:
                A = R[arm][(test[int(x)], ex)]; Bq = R[arm][(test[int(y)], ex)]
                if A[0] > Bq[1] + TIE: pr = 1
                elif Bq[0] > A[1] + TIE: pr = -1
                else: pr = 0
                if pr == 0: s.append(0.5)
                else:
                    v = 1.0 if pr == truth else 0.0; s.append(v); res += 1; corr += v
            per_int.append(float(np.mean(s))); n_int_pairs += len(s)
        S_by_arm[arm] = S; pairs_by_arm[arm] = sub
        # constancy and exploration over the primary panel
        const = zero = 0; expl = []; ident_free = 0; widths = []; per_target = {}
        for ex in primary_ex:
            pts = np.array([P[arm][(c, ex)] for c in ctx_have if (c, ex) in P[arm]])
            if len(pts) == 0: continue
            if pts.max() - pts.min() <= EPS: const += 1
            if np.all(np.abs(pts) <= EPS): zero += 1
            mids = np.array([0.5 * (R[arm][(c, ex)][0] + R[arm][(c, ex)][1]) for c in ctx_have if (c, ex) in R[arm]])
            widths.extend([R[arm][(c, ex)][1] - R[arm][(c, ex)][0] for c in ctx_have if (c, ex) in R[arm]])
            per_target[ex] = {"point_span": float(pts.max() - pts.min()), "midpoint_spread": float(mids.max() - mids.min()),
                              "mean_width": float(np.mean([R[arm][(c, ex)][1] - R[arm][(c, ex)][0] for c in ctx_have if (c, ex) in R[arm]]))}
            if ex in free:
                e = (mids.max() - mids.min()) / b0[ex]; expl.append(e); per_target[ex]["exploration_vs_study_B0"] = float(e)
                if mids.max() - mids.min() <= EPS: ident_free += 1
        widths = np.array(widths)
        out["arms"][arm] = {"profiles": len(ctx_have), "pairs_scored": n_pairs, "targets_scored": len(sub),
                            "point_macro_C": macro(S, sub), "interval_macro_C": float(np.mean(per_int)) if per_int else None,
                            "resolved_pairs": int(res), "resolved_share_of_scored": (res / n_int_pairs) if n_int_pairs else None,
                            "accuracy_among_resolved": (corr / res) if res else None,
                            "targets_constant_across_profiles_of_52": const, "targets_zero_everywhere_of_52": zero,
                            "free_targets_identical_across_profiles_of_%d" % len(free): ident_free,
                            "exploration_vs_study_B0_median": float(np.median(expl)) if expl else None,
                            "exploration_vs_study_B0_mean": float(np.mean(expl)) if expl else None,
                            "exploration_vs_study_B0_max": float(np.max(expl)) if expl else None,
                            "range_width_median": float(np.median(widths)) if len(widths) else None,
                            "range_width_mean": float(np.mean(widths)) if len(widths) else None,
                            "profile_target_ranges_above_eps": int((widths > EPS).sum()), "profile_target_ranges": int(len(widths)),
                            "per_target": per_target}
    # bootstrap intervals: absolute per arm on each arm's own pair set; grid contrast when all six arms exist on a common pair set
    boot = {}
    for arm in S_by_arm:
        d, _ = origin_bootstrap({arm: S_by_arm[arm]}, pairs_by_arm[arm], test, origin, a.B, 20260916)
        boot[arm] = ci(d[arm])
    out["absolute_ci95"] = boot
    if all(x in S_by_arm for x in ARMS):
        common = {m: pairs_by_arm["ordinal"][m] for m in pairs_by_arm["ordinal"]}
        d, dg = origin_bootstrap({x: S_by_arm[x] for x in ARMS}, common, test, origin, a.B, 20260916)
        out["delta_grid"] = out["arms"]["ordinal"]["point_macro_C"] - float(np.mean([out["arms"][x]["point_macro_C"] for x in MAG]))
        out["delta_grid_ci95"] = ci(dg)
    json.dump(out, open(a.out / ("comparator_%s.json" % a.label), "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
