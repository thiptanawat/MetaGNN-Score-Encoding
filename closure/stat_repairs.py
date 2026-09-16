#!/usr/bin/env python3
"""Statistical repairs on saved outputs. Nothing here solves an LP or reads a new prediction.

Order of operations, fixed so that no repair is trusted before the reimplementation is:
  0. Reproduce the locked evaluator's macro concordance for all seven arms to 1e-9. Abort otherwise.
  1. Common-resolved-subset contrast: the only comparison in the study with power, restricted
     to pairs every encoding resolves, with the same origin bootstrap.
  2. Assay tolerance recomputed from the reserved lines' own replicates, then eligibility and
     concordance recomputed under it. Tests the signed mechanism by which a transported
     development tolerance could drive the reserved concordance toward 0.5.
  3. Multiplicity: Bonferroni-adjusted bootstrap intervals for the 6 primary absolute estimates
     and the 2 contrasts, so the reader can see which "above 0.5" statements survive.
  4. Metabolite jackknife and a two-way (origin x target) bootstrap for the grid contrast.
  5. Bootstrap calibration by simulation under a known null (association destroyed by permuting
     context labels within target) and a known alternative (signal injected), reporting the
     empirical coverage of the origin bootstrap as implemented.
  6. Positive control for the exploration measure on synthetic intervals with known spread.
  7. Minimum detectable effect on the informative-target scale.
"""
import sys, csv, json, argparse, time
from collections import defaultdict
from pathlib import Path
import numpy as np

TIE = 1e-5
ARMS = ["magnitude_g0.50", "magnitude_g0.80", "magnitude_g1.00", "magnitude_g1.25", "magnitude_g2.00", "ordinal"]
MAG = [a for a in ARMS if a.startswith("magnitude")]


def load(root, scenario):
    ctx = {r["context_id"]: r for r in csv.DictReader(open(root / "manifests/contexts.tsv"), delimiter="\t")}
    test = sorted(c for c, r in ctx.items() if r["partition"] == "test")
    origin = {c: ctx[c]["origin_group"] for c in test}
    plan = json.loads((root / "protocol/analysis_plan.json").read_text())
    targets = [(t["metabolite_id"], t["exchange_id"]) for t in plan["primary_targets"]]
    tol = {m: float(v) for m, v in plan["assay_tolerances"].items()}
    obs = defaultdict(lambda: defaultdict(list))
    for r in csv.DictReader(open(root / "data/outcomes_long.tsv"), delimiter="\t"):
        if r["context_id"] in origin and r["value"] != "":
            obs[r["metabolite_id"]][r["context_id"]].append(float(r["value"]))
    obs_mean = {m: {c: float(np.mean(v)) for c, v in d.items()} for m, d in obs.items()}
    P = defaultdict(dict)
    for r in csv.DictReader(open(root / "results/reserved_v3/predictions.tsv"), delimiter="\t"):
        if r["scenario"] == scenario and r["status"] == "optimal":
            P[r["arm"]][(r["context_id"], r["exchange_id"])] = float(r["point"])
    return test, origin, targets, tol, obs, obs_mean, P


def build_pairs(test, origin, targets, tol, obs_mean, P, arms):
    """Eligible pair arrays per target: (i, j, truth) with different origins and measured gap > tol."""
    out = {}
    for m, ex in targets:
        rows = []
        o = obs_mean.get(m, {})
        for a in range(len(test)):
            for b in range(a + 1, len(test)):
                ci, cj = test[a], test[b]
                if origin[ci] == origin[cj] or ci not in o or cj not in o: continue
                if abs(o[ci] - o[cj]) <= tol[m]: continue
                if any((ci, ex) not in P[arm] or (cj, ex) not in P[arm] for arm in arms): continue
                rows.append((a, b, 1.0 if o[ci] > o[cj] else -1.0))
        if rows:
            out[m] = np.array(rows)
    return out


def scores(pairs, P, arm, test, ex_of):
    """Per-target vector of pair scores in {0, 0.5, 1} for one arm."""
    S = {}
    for m, arr in pairs.items():
        ex = ex_of[m]
        pi = np.array([P[arm][(test[int(a)], ex)] for a in arr[:, 0]])
        pj = np.array([P[arm][(test[int(b)], ex)] for b in arr[:, 1]])
        d = pi - pj
        s = np.where(np.abs(d) <= TIE, 0.5, (np.sign(d) == arr[:, 2]).astype(float))
        S[m] = s
    return S


def macro(S, pairs, w=None):
    """Macro concordance: equal weight per target, weighted mean within target."""
    vals = []
    for m, s in S.items():
        ww = np.ones(len(s)) if w is None else w[m]
        if ww.sum() > 0: vals.append(float((s * ww).sum() / ww.sum()))
    return float(np.mean(vals)) if vals else float("nan")


def origin_bootstrap(S_by_arm, pairs, test, origin, B, seed):
    """Resample origin groups; a pair's weight is the product of its two origins' multiplicities."""
    groups = sorted(set(origin.values())); gi = {g: i for i, g in enumerate(groups)}
    oi = np.array([gi[origin[c]] for c in test])
    rng = np.random.default_rng(seed)
    draws = {a: [] for a in S_by_arm}; dgrid = []
    for _ in range(B):
        mult = np.bincount(rng.choice(len(groups), len(groups), replace=True), minlength=len(groups))
        w = {m: mult[oi[arr[:, 0].astype(int)]] * mult[oi[arr[:, 1].astype(int)]] for m, arr in pairs.items()}
        est = {a: macro(S_by_arm[a], pairs, w) for a in S_by_arm}
        for a in est: draws[a].append(est[a])
        if "ordinal" in est and all(a in est for a in MAG):
            dgrid.append(est["ordinal"] - np.mean([est[a] for a in MAG]))
    return {a: np.array(v) for a, v in draws.items()}, np.array(dgrid)


def ci(x, alpha=0.05):
    x = x[np.isfinite(x)]
    return [float(np.percentile(x, 100 * alpha / 2)), float(np.percentile(x, 100 * (1 - alpha / 2)))]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", type=Path, required=True); ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--B", type=int, default=2000); ap.add_argument("--calib-datasets", type=int, default=100); ap.add_argument("--calib-B", type=int, default=500)
    a = ap.parse_args(); a.out.mkdir(parents=True, exist_ok=True); t0 = time.time()
    test, origin, targets, tol, obs, obs_mean, P = load(a.root, "primary")
    ex_of = {m: ex for m, ex in targets}
    all_arms = ARMS + ["uniform_pfba"]
    # the uniform arm is stored once under __reference__; broadcast it to every context
    if "uniform_pfba" in P:
        ref = {ex: v for (c, ex), v in P["uniform_pfba"].items() if c == "__reference__"}
        P["uniform_pfba"] = {(c, ex): v for c in test for ex, v in ref.items()}
    report = {}

    # ---- 0. reproduction of the locked evaluator
    pairs = build_pairs(test, origin, targets, tol, obs_mean, P, all_arms)
    S = {arm: scores(pairs, P, arm, test, ex_of) for arm in all_arms}
    saved = json.load(open(a.root / "results/reserved_v3/independent_concordance_primary.json"))
    mine = {arm: macro(S[arm], pairs) for arm in all_arms}
    npairs = int(sum(len(v) for v in pairs.values()))
    diffs = {arm: abs(mine[arm] - saved["absolute_C"][arm]) for arm in all_arms}
    report["reproduction"] = {"eligible_pairs_mine": npairs, "eligible_pairs_saved": saved["observed_pairs"],
                              "targets_mine": len(pairs), "targets_saved": saved["targets"],
                              "macro_C_mine": mine, "macro_C_saved": saved["absolute_C"], "max_abs_diff": max(diffs.values())}
    print("reproduction: pairs %d vs %d | targets %d vs %d | max |diff| %.2e" % (npairs, saved["observed_pairs"], len(pairs), saved["targets"], max(diffs.values())), flush=True)
    if max(diffs.values()) > 1e-9 or npairs != saved["observed_pairs"]:
        print("ABORT: reimplementation does not reproduce the locked evaluator"); json.dump(report, open(a.out / "stat_repairs.json", "w"), indent=1); sys.exit(2)

    # ---- 1. common-resolved-subset contrast
    common = {}
    for m, arr in pairs.items():
        keep = np.ones(len(arr), bool)
        for arm in ARMS:
            ex = ex_of[m]
            d = np.array([P[arm][(test[int(x)], ex)] - P[arm][(test[int(y)], ex)] for x, y in arr[:, :2]])
            keep &= np.abs(d) > TIE
        if keep.sum() >= 1: common[m] = arr[keep]
    S_c = {arm: scores(common, P, arm, test, ex_of) for arm in ARMS}
    est_c = {arm: macro(S_c[arm], common) for arm in ARMS}
    dr_c, dg_c = origin_bootstrap(S_c, common, test, origin, a.B, 20260916)
    report["common_resolved_subset"] = {"pairs": int(sum(len(v) for v in common.values())), "targets": len(common),
        "macro_C": est_c, "delta_grid": est_c["ordinal"] - float(np.mean([est_c[x] for x in MAG])), "delta_grid_ci95": ci(dg_c),
        "delta_identity": est_c["ordinal"] - est_c["magnitude_g1.00"],
        "absolute_ci95": {arm: ci(dr_c[arm]) for arm in ARMS},
        "note": "pairs resolved by all six encodings; subset chosen by predictions, not outcomes, but the resolution sets are encoding dependent so this is a supporting analysis"}
    print("common subset: %d pairs over %d targets | delta_grid %+.4f %s" % (report["common_resolved_subset"]["pairs"], len(common), report["common_resolved_subset"]["delta_grid"], ci(dg_c)), flush=True)

    # ---- 2. tolerance recomputed from reserved replicates
    tol_res = {}
    for m, _ in targets:
        d = [abs(v[0] - v[1]) for c, v in obs[m].items() if len(v) == 2]
        tol_res[m] = float(np.percentile(d, 90)) if d else tol[m]
    pairs_r = build_pairs(test, origin, targets, tol_res, obs_mean, P, all_arms)
    S_r = {arm: scores(pairs_r, P, arm, test, ex_of) for arm in all_arms}
    est_r = {arm: macro(S_r[arm], pairs_r) for arm in all_arms}
    dr_r, dg_r = origin_bootstrap({arm: S_r[arm] for arm in ARMS}, pairs_r, test, origin, a.B, 20260917)
    report["reserved_tolerance"] = {"tolerance_ratio_reserved_over_development": {m: tol_res[m] / tol[m] if tol[m] > 0 else None for m, _ in targets},
        "median_ratio": float(np.median([tol_res[m] / tol[m] for m, _ in targets if tol[m] > 0])),
        "eligible_pairs": int(sum(len(v) for v in pairs_r.values())), "eligible_pairs_development_tolerance": npairs,
        "targets": len(pairs_r), "macro_C": est_r, "delta_grid": est_r["ordinal"] - float(np.mean([est_r[x] for x in MAG])), "delta_grid_ci95": ci(dg_r),
        "absolute_ci95": {arm: ci(dr_r[arm]) for arm in ARMS}}
    print("reserved tolerance: median ratio %.2f | pairs %d -> %d | ordinal C %.4f | delta %+.4f %s" % (report["reserved_tolerance"]["median_ratio"], npairs, report["reserved_tolerance"]["eligible_pairs"], est_r["ordinal"], report["reserved_tolerance"]["delta_grid"], ci(dg_r)), flush=True)

    # ---- 3. multiplicity on the primary analysis
    dr, dg = origin_bootstrap({arm: S[arm] for arm in ARMS}, pairs, test, origin, a.B, 20260913)
    m_abs, m_con = 18, 6
    report["multiplicity"] = {"absolute_unadjusted_ci95": {arm: ci(dr[arm]) for arm in ARMS},
        "absolute_bonferroni_over_18": {arm: ci(dr[arm], 0.05 / m_abs) for arm in ARMS},
        "above_0.5_unadjusted": [arm for arm in ARMS if ci(dr[arm])[0] > 0.5],
        "above_0.5_bonferroni": [arm for arm in ARMS if ci(dr[arm], 0.05 / m_abs)[0] > 0.5],
        "delta_grid_unadjusted": ci(dg), "delta_grid_bonferroni_over_6": ci(dg, 0.05 / m_con),
        "note": "my bootstrap uses its own RNG, so unadjusted intervals are close to but not identical with the locked ones; the comparison of adjusted and unadjusted is within one bootstrap"}
    print("multiplicity: above 0.5 unadjusted %s | Bonferroni %s" % (report["multiplicity"]["above_0.5_unadjusted"], report["multiplicity"]["above_0.5_bonferroni"]), flush=True)

    # ---- 4. metabolite jackknife and two-way bootstrap
    mets = sorted(pairs); jk = []
    for m in mets:
        sub = {x: pairs[x] for x in mets if x != m}
        e = {arm: macro({x: S[arm][x] for x in sub}, sub) for arm in ARMS}
        jk.append(e["ordinal"] - float(np.mean([e[x] for x in MAG])))
    rng = np.random.default_rng(20260918); groups = sorted(set(origin.values())); gi = {g: i for i, g in enumerate(groups)}
    oi = np.array([gi[origin[c]] for c in test]); tw = []
    for _ in range(a.B):
        mult = np.bincount(rng.choice(len(groups), len(groups), replace=True), minlength=len(groups))
        pick = rng.choice(len(mets), len(mets), replace=True)
        w = {m: mult[oi[pairs[m][:, 0].astype(int)]] * mult[oi[pairs[m][:, 1].astype(int)]] for m in mets}
        e = {}
        for arm in ARMS:
            vals = [(S[arm][mets[i]] * w[mets[i]]).sum() / w[mets[i]].sum() for i in pick if w[mets[i]].sum() > 0]
            e[arm] = float(np.mean(vals))
        tw.append(e["ordinal"] - float(np.mean([e[x] for x in MAG])))
    report["target_resampling"] = {"jackknife_delta_grid_range": [float(min(jk)), float(max(jk))], "jackknife_n": len(jk),
        "two_way_origin_target_delta_grid_ci95": ci(np.array(tw)), "one_way_origin_delta_grid_ci95": ci(dg)}
    print("two-way bootstrap delta_grid %s vs one-way %s" % (ci(np.array(tw)), ci(dg)), flush=True)

    # ---- 5. bootstrap calibration by simulation
    def coverage(make_preds, truth_fn, label):
        hits_d, hits_c, widths = 0, 0, []
        for rep in range(a.calib_datasets):
            Pn = make_preds(np.random.default_rng(1000 + rep))
            Sn = {arm: scores(pairs, Pn, arm, test, ex_of) for arm in ARMS}
            drn, dgn = origin_bootstrap(Sn, pairs, test, origin, a.calib_B, 5000 + rep)
            lo, hi = ci(dgn); truth_d, truth_c = truth_fn(Pn)
            hits_d += lo <= truth_d <= hi; widths.append(hi - lo)
            lc, hc = ci(drn["ordinal"]); hits_c += lc <= truth_c <= hc
        return {"datasets": a.calib_datasets, "resamples": a.calib_B, "coverage_delta_grid": hits_d / a.calib_datasets,
                "coverage_absolute_ordinal": hits_c / a.calib_datasets, "mean_interval_width": float(np.mean(widths))}
    def null_preds(rng):
        Pn = {}
        for arm in ARMS:
            Pn[arm] = {}
            for m, ex in targets:
                vals = [P[arm][(c, ex)] for c in test]; perm = rng.permutation(len(test))
                for i, c in enumerate(test): Pn[arm][(c, ex)] = vals[perm[i]]
        return Pn
    def null_truth(Pn): return 0.0, 0.5
    report["calibration_null"] = coverage(null_preds, null_truth, "null")
    print("calibration null: coverage delta %.3f | absolute %.3f" % (report["calibration_null"]["coverage_delta_grid"], report["calibration_null"]["coverage_absolute_ordinal"]), flush=True)
    # known alternative: ordinal follows the outcome with noise, magnitude arms are permuted; truth by large Monte Carlo
    def alt_preds(rng, noise=1.0):
        Pn = null_preds(rng)
        for m, ex in targets:
            o = obs_mean.get(m, {}); sd = np.std([o[c] for c in test if c in o]) or 1.0
            for c in test:
                Pn["ordinal"][(c, ex)] = (o.get(c, 0.0) + rng.normal(0, noise * sd))
        return Pn
    truths = []
    for rep in range(200):
        Pn = alt_preds(np.random.default_rng(77000 + rep)); Sn = {arm: scores(pairs, Pn, arm, test, ex_of) for arm in ARMS}
        e = {arm: macro(Sn[arm], pairs) for arm in ARMS}; truths.append((e["ordinal"] - np.mean([e[x] for x in MAG]), e["ordinal"]))
    td, tc = float(np.mean([t[0] for t in truths])), float(np.mean([t[1] for t in truths]))
    report["calibration_alternative"] = coverage(alt_preds, lambda Pn: (td, tc), "alt"); report["calibration_alternative"]["true_delta_grid"] = td; report["calibration_alternative"]["true_ordinal_C"] = tc
    print("calibration alt: true delta %.4f | coverage delta %.3f | absolute %.3f" % (td, report["calibration_alternative"]["coverage_delta_grid"], report["calibration_alternative"]["coverage_absolute_ordinal"]), flush=True)

    # ---- 6. positive control for the exploration measure
    rng = np.random.default_rng(1); ex_pinned, ex_spread = [], []
    for _ in range(52):
        lo, hi = -10.0, 25.0; width = hi - lo
        mids_p = np.full(47, lo + 0.3 * width) + rng.normal(0, 1e-9, 47)
        mids_s = lo + rng.uniform(0, 1, 47) * width
        ex_pinned.append((mids_p.max() - mids_p.min()) / width); ex_spread.append((mids_s.max() - mids_s.min()) / width)
    report["exploration_positive_control"] = {"pinned_profiles_median_exploration": float(np.median(ex_pinned)),
        "uniformly_placed_profiles_median_exploration": float(np.median(ex_spread)),
        "expected_for_47_uniform_placements": float(1 - 2 / 48), "note": "synthetic intervals; the measure returns ~0 for pinned and ~0.96 for uniformly placed midpoints"}

    # ---- 7. minimum detectable effect on the informative targets
    n_const = sum(1 for m in pairs if all(np.all(S[arm][m] == 0.5) for arm in MAG))
    n_inf = len(pairs) - n_const; hw = (ci(dg)[1] - ci(dg)[0]) / 2
    report["minimum_detectable_effect"] = {"targets": len(pairs), "constant_targets_all_magnitude_arms": n_const, "informative_targets": n_inf,
        "delta_grid_halfwidth_macro": hw, "mde_on_informative_targets": hw * len(pairs) / max(n_inf, 1),
        "informative_mean_C_needed_for_observed": {arm: (mine[arm] * len(pairs) - 0.5 * n_const) / max(n_inf, 1) for arm in ARMS}}
    print("MDE: %d constant of %d | halfwidth %.4f | MDE informative %.4f" % (n_const, len(pairs), hw, report["minimum_detectable_effect"]["mde_on_informative_targets"]), flush=True)

    report["_seconds"] = round(time.time() - t0, 1)
    json.dump(report, open(a.out / "stat_repairs.json", "w"), indent=1, default=float)
    print("done %.0fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
