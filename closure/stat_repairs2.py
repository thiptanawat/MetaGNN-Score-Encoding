#!/usr/bin/env python3
"""Second round of supporting analyses on saved outputs. No LP is solved.

  1. Value-level comparison of the uniform cost with every weighted arm on the exchange coordinates:
     do the coordinates the readout fixes coincide in value with the uniform solution, or only in number?
  2. Minimum detectable effect restated: informative targets as the union over the arms entering the
     contrast, the half-width criterion (about 50 percent power) and the 80 percent power criterion.
  3. Simulated power curve: signal injected into the midrank arm on the informative targets at a series of
     true contrasts, with the origin bootstrap as implemented; power is the share of datasets whose interval
     excludes zero.
  4. Threshold sensitivity: the headline counts recomputed with the tie and resolution threshold at 1e-4 and
     1e-6 instead of 1e-5 (constant targets, fixed coordinates in the nested-set ledger, resolved pairs,
     macro concordance, grid contrast).
  5. Accuracy among resolved pairs with origin-bootstrap intervals, per encoding.
  6. Informative-subset concordance for any predictions file (targets whose predictions vary across profiles),
     used for the study arms and for the comparator formulations.
  7. The scale constant k0 recomputed from the reserved profiles' own evidence.
  8. Tolerance-level range inversions in the saved ranges: counts, magnitudes, and whether any inverted
     coordinate entered a resolved pair.
"""
import sys, csv, json, argparse, time
from collections import defaultdict
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from stat_repairs import load, build_pairs, scores, macro, origin_bootstrap, ci, ARMS, MAG, TIE


def read_pred_file(path, scenario="primary"):
    P = defaultdict(dict); R = defaultdict(dict)
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["status"] != "optimal" or r["point"] == "": continue
        if r.get("scenario", scenario) != scenario: continue
        P[r["arm"]][(r["context_id"], r["exchange_id"])] = float(r["point"])
        R[r["arm"]][(r["context_id"], r["exchange_id"])] = [float(r["range_lo"]), float(r["range_hi"])]
    return P, R


def macro_with_tie(pairs, P, arm, test, ex_of, tie):
    S = {}
    for m, arr in pairs.items():
        ex = ex_of[m]
        pi = np.array([P[arm][(test[int(a)], ex)] for a in arr[:, 0]]); pj = np.array([P[arm][(test[int(b)], ex)] for b in arr[:, 1]])
        d = pi - pj; S[m] = np.where(np.abs(d) <= tie, 0.5, (np.sign(d) == arr[:, 2]).astype(float))
    return S


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", type=Path, required=True); ap.add_argument("--ledger", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True); ap.add_argument("--B", type=int, default=2000)
    ap.add_argument("--power-datasets", type=int, default=100); ap.add_argument("--power-B", type=int, default=300)
    ap.add_argument("--comparators", nargs="*", default=[], help="label=path/to/predictions.tsv")
    a = ap.parse_args(); a.out.mkdir(parents=True, exist_ok=True); t0 = time.time(); rep = {}
    test, origin, targets, tol, obs, obs_mean, P = load(a.root, "primary")
    ex_of = {m: ex for m, ex in targets}; primary_ex = [ex for m, ex in targets]
    ref = {ex: v for (c, ex), v in P["uniform_pfba"].items() if c == "__reference__"}
    P["uniform_pfba"] = {(c, ex): v for c in test for ex, v in ref.items()}
    pairs = build_pairs(test, origin, targets, tol, obs_mean, P, ARMS + ["uniform_pfba"])
    assert sum(len(v) for v in pairs.values()) == 46447
    P_all, R_all = read_pred_file(a.root / "results/reserved_v3/predictions.tsv", "primary")
    panel96 = sorted({ex for (c, ex) in P_all["magnitude_g1.00"]})
    S_by_arm = {arm: scores(pairs, P, arm, test, ex_of) for arm in ARMS}

    # ---- 1. uniform versus weighted, value level
    val = {}
    for arm in ARMS:
        n_diff52, n_diff96, l1 = [], [], []
        for c in test:
            v52 = np.array([P[arm][(c, ex)] for ex in primary_ex]); u52 = np.array([ref[ex] for ex in primary_ex])
            v96 = np.array([P_all[arm][(c, ex)] for ex in panel96]); u96 = np.array([ref[ex] for ex in panel96])
            n_diff52.append(int((np.abs(v52 - u52) > TIE).sum())); n_diff96.append(int((np.abs(v96 - u96) > TIE).sum()))
            l1.append(2 * np.abs(v96 - u96).sum() / max(np.abs(v96).sum() + np.abs(u96).sum(), 1e-12))
        const_equal = const_total = 0; const_diff_targets = []
        for ex in primary_ex:
            pts = np.array([P[arm][(c, ex)] for c in test])
            if pts.max() - pts.min() <= TIE:
                const_total += 1
                if abs(pts.mean() - ref[ex]) <= TIE: const_equal += 1
                else: const_diff_targets.append(ex)
        val[arm] = {"targets_of_52_differing_from_uniform_median_over_profiles": float(np.median(n_diff52)),
                    "targets_of_52_differing_from_uniform_min_max": [int(min(n_diff52)), int(max(n_diff52))],
                    "exchanges_of_96_differing_from_uniform_median": float(np.median(n_diff96)),
                    "relative_l1_to_uniform_over_96_median": float(np.median(l1)), "relative_l1_to_uniform_over_96_max": float(np.max(l1)),
                    "constant_targets": const_total, "constant_targets_equal_to_uniform_value": const_equal,
                    "constant_targets_differing_from_uniform": const_diff_targets}
    rep["uniform_value_comparison"] = val
    print("1. uniform value comparison:", {k: (v["targets_of_52_differing_from_uniform_median_over_profiles"], v["constant_targets_equal_to_uniform_value"], v["constant_targets"]) for k, v in val.items()}, flush=True)

    # ---- 2. MDE restated
    const_by_arm = {arm: {ex for ex in primary_ex if np.ptp([P[arm][(c, ex)] for c in test]) <= TIE} for arm in ARMS}
    const_all_mag = set.intersection(*[const_by_arm[x] for x in MAG]); const_all_six = set.intersection(*[const_by_arm[x] for x in ARMS])
    informative_union = 52 - len(const_all_six); informative_mag = 52 - len(const_all_mag)
    dr, dg = origin_bootstrap(S_by_arm, pairs, test, origin, a.B, 20260916)
    hw = (ci(dg)[1] - ci(dg)[0]) / 2
    z = 1.959964; z80 = 0.841621
    rep["mde"] = {"halfwidth_macro": hw, "informative_targets_magnitude_arms": informative_mag, "informative_targets_any_arm": informative_union,
                  "constant_in_all_six_arms": len(const_all_six), "constant_in_all_magnitude_arms": len(const_all_mag),
                  "mde_halfwidth_on_9": hw * 52 / informative_mag, "mde_halfwidth_on_union": hw * 52 / informative_union,
                  "mde_80pct_power_on_9": hw * 52 / informative_mag * (z + z80) / z, "mde_80pct_power_on_union": hw * 52 / informative_union * (z + z80) / z,
                  "note": "half-width criterion corresponds to about 50 percent power; 80 percent power multiplies by (z_0.975 + z_0.80)/z_0.975 = 1.43"}
    print("2. MDE:", {k: (round(v, 4) if isinstance(v, float) else v) for k, v in rep["mde"].items() if k != "note"}, flush=True)

    # ---- 3. simulated power curve (signal injected into the midrank arm on the informative targets)
    inf_targets = [m for m, ex in targets if ex not in const_all_six]
    def null_preds(rng):
        Pn = {}
        for arm in ARMS:
            Pn[arm] = {}
            for m, ex in targets:
                vals = [P[arm][(c, ex)] for c in test]; perm = rng.permutation(len(test))
                for i, c in enumerate(test): Pn[arm][(c, ex)] = vals[perm[i]]
        return Pn
    def inject(rng, p):
        """With probability p per informative target, the midrank prediction follows the measured ordering (plus small noise)."""
        Pn = null_preds(rng)
        for m in inf_targets:
            ex = ex_of[m]
            if rng.random() < p:
                o = obs_mean.get(m, {}); sd = np.std([o[c] for c in test if c in o]) or 1.0
                for c in test: Pn["ordinal"][(c, ex)] = o.get(c, 0.0) + rng.normal(0, 0.25 * sd)
        return Pn
    def true_delta(p, n=60):
        vals = []
        for rep_i in range(n):
            Pn = inject(np.random.default_rng(90000 + rep_i), p); Sn = {arm: scores(pairs, Pn, arm, test, ex_of) for arm in ARMS}
            e = {arm: macro(Sn[arm], pairs) for arm in ARMS}; vals.append(e["ordinal"] - np.mean([e[x] for x in MAG]))
        return float(np.mean(vals))
    curve = []
    for target_inf in (0.02, 0.04, 0.06, 0.08, 0.10):
        target_macro = target_inf * len(inf_targets) / 52
        lo, hi = 0.0, 1.0
        for _ in range(10):                                   # bisection on the injection probability
            mid = 0.5 * (lo + hi); d = true_delta(mid, 30)
            if d < target_macro: lo = mid
            else: hi = mid
        p = 0.5 * (lo + hi); td = true_delta(p, 60)
        hits = 0
        for rep_i in range(a.power_datasets):
            Pn = inject(np.random.default_rng(7000 + rep_i), p); Sn = {arm: scores(pairs, Pn, arm, test, ex_of) for arm in ARMS}
            _, dgn = origin_bootstrap(Sn, pairs, test, origin, a.power_B, 8000 + rep_i)
            l, h = ci(dgn); hits += (l > 0) or (h < 0)
        curve.append({"target_contrast_informative_scale": target_inf, "target_contrast_macro": target_macro, "injection_probability": p,
                      "realized_true_contrast_macro": td, "realized_true_contrast_informative_scale": td * 52 / len(inf_targets),
                      "power": hits / a.power_datasets, "datasets": a.power_datasets, "resamples": a.power_B})
        print("3. power at informative contrast %.2f (macro %.4f, realized %.4f): %.2f" % (target_inf, target_macro, td, hits / a.power_datasets), flush=True)
    rep["power_curve"] = {"informative_targets": len(inf_targets), "curve": curve}

    # ---- 4. threshold sensitivity
    ledger = list(csv.DictReader(open(a.ledger), delimiter="\t"))
    thr_out = {}
    for thr in (1e-6, 1e-5, 1e-4):
        d = {}
        for arm in ARMS:
            const = sum(1 for ex in primary_ex if np.ptp([P[arm][(c, ex)] for c in test]) <= thr)
            Sx = macro_with_tie(pairs, P, arm, test, ex_of, thr); C = macro(Sx, pairs)
            res = 0
            for m, arr in pairs.items():
                ex = ex_of[m]
                for x, y, truth in arr:
                    A = R_all[arm][(test[int(x)], ex)]; Bq = R_all[arm][(test[int(y)], ex)]
                    if A[0] > Bq[1] + thr or Bq[0] > A[1] + thr: res += 1
            sub = [r for r in ledger if r["arm"] == arm]
            b0free = [r for r in sub if float(r["B0_width"]) > thr]
            first_b1 = sum(1 for r in b0free if float(r["B1_width"]) <= thr)
            d[arm] = {"constant_targets": const, "macro_C": C, "resolved_pairs": res, "free_coordinates": len(b0free), "first_fixed_B1": first_b1,
                      "share_first_fixed_B1": first_b1 / max(len(b0free), 1)}
        d["delta_grid"] = d["ordinal"]["macro_C"] - float(np.mean([d[x]["macro_C"] for x in MAG]))
        thr_out["%.0e" % thr] = d
        print("4. threshold %.0e: const %s | C ordinal %.4f | resolved ordinal %d | share B1 identity %.3f" % (thr, [d[x]["constant_targets"] for x in ARMS], d["ordinal"]["macro_C"], d["ordinal"]["resolved_pairs"], d["magnitude_g1.00"]["share_first_fixed_B1"]), flush=True)
    rep["threshold_sensitivity"] = thr_out

    # ---- 5. accuracy among resolved pairs with origin-bootstrap intervals
    groups = sorted(set(origin.values())); gi = {g: i for i, g in enumerate(groups)}; oi = np.array([gi[origin[c]] for c in test])
    acc = {}
    for arm in ARMS:
        rows = []                                              # (i, j, correct)
        for m, arr in pairs.items():
            ex = ex_of[m]
            for x, y, truth in arr:
                A = R_all[arm][(test[int(x)], ex)]; Bq = R_all[arm][(test[int(y)], ex)]
                pr = 1 if A[0] > Bq[1] + TIE else (-1 if Bq[0] > A[1] + TIE else 0)
                if pr: rows.append((int(x), int(y), 1.0 if pr == truth else 0.0))
        arr = np.array(rows); rng = np.random.default_rng(20260916); draws = []
        for _ in range(a.B):
            mult = np.bincount(rng.choice(len(groups), len(groups), replace=True), minlength=len(groups))
            w = mult[oi[arr[:, 0].astype(int)]] * mult[oi[arr[:, 1].astype(int)]]
            if w.sum() > 0: draws.append(float((arr[:, 2] * w).sum() / w.sum()))
        acc[arm] = {"resolved_pairs": len(rows), "accuracy": float(arr[:, 2].mean()), "ci95": ci(np.array(draws))}
    rep["accuracy_among_resolved"] = acc
    print("5. accuracy among resolved:", {k: (v["resolved_pairs"], round(v["accuracy"], 4), [round(x, 4) for x in v["ci95"]]) for k, v in acc.items()}, flush=True)

    # ---- 6. informative-subset concordance for the study arms and the comparators
    def informative_subset(Pf, arm, ctx_list):
        inf = [m for m, ex in targets if all((c, ex) in Pf[arm] for c in ctx_list) and np.ptp([Pf[arm][(c, ex)] for c in ctx_list]) > TIE]
        sub = {m: pairs[m] for m in inf if m in pairs}
        keep = {}
        for m, arr in sub.items():
            ex = ex_of[m]; rows = [row for row in arr if (test[int(row[0])], ex) in Pf[arm] and (test[int(row[1])], ex) in Pf[arm]]
            if rows: keep[m] = np.array(rows)
        S = scores(keep, Pf, arm, test, ex_of) if keep else {}
        return {"informative_targets": len(inf), "pairs": int(sum(len(v) for v in keep.values())), "macro_C_informative": macro(S, keep) if keep else None}
    info = {"study": {arm: informative_subset(P, arm, test) for arm in ARMS}}
    for spec in a.comparators:
        label, path = spec.split("=", 1); Pc, Rc = read_pred_file(path, "primary"); info[label] = {}
        for arm in sorted(Pc):
            ctx_list = sorted({c for (c, ex) in Pc[arm]}); info[label][arm] = informative_subset(Pc, arm, ctx_list)
    rep["informative_subset_concordance"] = info
    print("6. informative-subset C:", {k: {a2: (v2["informative_targets"], None if v2["macro_C_informative"] is None else round(v2["macro_C_informative"], 4)) for a2, v2 in v.items()} for k, v in info.items()}, flush=True)

    # ---- 7. k0 from the reserved profiles
    d = np.load(a.root / "data/scores_conservative.npz", allow_pickle=False); ctxs = d["contexts"].astype(str).tolist()
    A = d["A"]; res_cols = [ctxs.index(c) for c in test]
    dev_cols = [i for i, c in enumerate(ctxs) if c not in test]
    k_res = float(np.median(A[:, res_cols][A[:, res_cols] > 0])); k_dev = float(np.median(A[:, dev_cols][A[:, dev_cols] > 0])) if dev_cols else None
    k_locked = float(json.loads((a.root / "results/development_v3/run_summary_scalars.json").read_text())["k"])
    rep["scale_constant"] = {"k_locked": k_locked, "k_from_development_columns_here": k_dev, "k_from_reserved_profiles": k_res, "positive_values_reserved": int((A[:, res_cols] > 0).sum())}
    print("7. k0 locked %.3f, development recomputed %s, reserved %.3f" % (k_locked, k_dev, k_res), flush=True)

    # ---- 8. range inversions
    inv = {}
    for scen in ("primary", "half_serum", "lower_task"):
        Ps, Rs = read_pred_file(a.root / "results/reserved_v3/predictions.tsv", scen)
        n = 0; mags = []; in_resolved = 0
        for arm in Rs:
            for (c, ex), (lo, hi) in Rs[arm].items():
                if lo > hi: n += 1; mags.append(lo - hi)
        if scen == "primary":
            for arm in ARMS:
                for m, arr in pairs.items():
                    ex = ex_of[m]
                    for x, y, truth in arr:
                        A_ = Rs[arm][(test[int(x)], ex)]; B_ = Rs[arm][(test[int(y)], ex)]
                        if (A_[0] > A_[1] or B_[0] > B_[1]) and (A_[0] > B_[1] + TIE or B_[0] > A_[1] + TIE): in_resolved += 1
        inv[scen] = {"inverted_intervals_over_all_arms_and_96_exchanges": n, "max_inversion_magnitude": float(max(mags)) if mags else 0.0,
                     "inverted_coordinates_entering_a_resolved_primary_pair": in_resolved if scen == "primary" else None}
    rep["range_inversions"] = inv
    print("8. inversions:", inv, flush=True)
    rep["_seconds"] = round(time.time() - t0, 1)
    json.dump(rep, open(a.out / "stat_repairs2.json", "w"), indent=1)
    print("FINISHED %.0fs" % rep["_seconds"], flush=True)


if __name__ == "__main__":
    main()
