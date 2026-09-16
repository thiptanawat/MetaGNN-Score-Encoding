#!/usr/bin/env python3
"""Does loosening the cap trade pinning for information?

At the near-exact cap every coordinate is a point; at looser caps the ranges are wide and their
positions differ between profiles. Wide ranges can still resolve an ordering when they do not
overlap. This scores the cap-sweep ranges against the measured CORE ordering with the interval
rule, on the same 46,447 eligible pairs the locked evaluation used: a pair is resolved when the
two intervals are disjoint by the tie tolerance, scored 1 or 0 by direction, and unresolved
pairs receive 0.5. Coverage and conditional accuracy are reported together.

The pair set and the scoring conventions are the ones the reimplementation reproduced to the
last digit against the locked evaluator; only the predictions differ.
"""
import sys, json, argparse, glob
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from stat_repairs import load, build_pairs, TIE, MAG

CAPS = ["exact", "rel_1pct", "rel_5pct", "rel_20pct"]
import os
if os.environ.get("CAP_NAMES"): CAPS = os.environ["CAP_NAMES"].split(",")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--sweep", type=Path, required=True); ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(); a.out.mkdir(parents=True, exist_ok=True)
    test, origin, targets, tol, obs, obs_mean, P = load(a.root, "primary")
    ex_of = {m: ex for m, ex in targets}
    all_arms = ["magnitude_g0.50", "magnitude_g0.80", "magnitude_g1.00", "magnitude_g1.25", "magnitude_g2.00", "ordinal", "uniform_pfba"]
    ref = {ex: v for (c, ex), v in P["uniform_pfba"].items() if c == "__reference__"}
    P["uniform_pfba"] = {(c, ex): v for c in test for ex, v in ref.items()}
    pairs = build_pairs(test, origin, targets, tol, obs_mean, P, all_arms)
    assert sum(len(v) for v in pairs.values()) == 46447

    R = {}                                   # arm -> cap -> (context, ex) -> [lo, hi]
    for f in glob.glob(str(a.sweep / "jobs" / "*.json")):
        d = json.loads(Path(f).read_text())
        if d["status"] != "complete": continue
        for cap in CAPS:
            R.setdefault(d["arm"], {}).setdefault(cap, {})
            for ex, v in d["caps"][cap]["ranges"].items():
                R[d["arm"]][cap][(d["context"], ex)] = v
    if "uniform_pfba" in R:                   # context independent: broadcast the single job
        for cap in CAPS:
            one = {ex: v for (c, ex), v in R["uniform_pfba"][cap].items()}
            R["uniform_pfba"][cap] = {(c, ex): v for c in test for ex, v in one.items()}

    out = {}
    for arm in sorted(R):
        out[arm] = {}
        for cap in CAPS:
            per, res, corr = [], 0, 0
            for m, arr in pairs.items():
                ex = ex_of[m]; s = []
                for x, y, truth in arr:
                    A = R[arm][cap].get((test[int(x)], ex)); Bq = R[arm][cap].get((test[int(y)], ex))
                    if A is None or Bq is None: continue
                    if A[0] > Bq[1] + TIE: pr = 1
                    elif Bq[0] > A[1] + TIE: pr = -1
                    else: pr = 0
                    if pr == 0: s.append(0.5)
                    else:
                        v = 1.0 if pr == truth else 0.0; s.append(v); res += 1; corr += v
                if s: per.append(float(np.mean(s)))
            out[arm][cap] = {"macro_interval_C": float(np.mean(per)), "targets": len(per), "resolved_pairs": res,
                             "resolved_share": res / 46447, "accuracy_among_resolved": (corr / res) if res else None}
    if "ordinal" in out:
        for cap in CAPS:
            mags = [out[x][cap]["macro_interval_C"] for x in MAG if x in out]
            out.setdefault("_contrast", {})[cap] = {"delta_grid_interval": out["ordinal"][cap]["macro_interval_C"] - float(np.mean(mags)) if mags else None}
    json.dump(out, open(a.out / "cap_concordance.json", "w"), indent=1)
    print("%-16s %-10s %8s %10s %10s %10s" % ("arm", "cap", "macro C", "resolved", "share", "acc|res"))
    for arm in sorted(k for k in out if not k.startswith("_")):
        for cap in CAPS:
            v = out[arm][cap]
            print("%-16s %-10s %8.4f %10d %10.4f %10s" % (arm, cap, v["macro_interval_C"], v["resolved_pairs"], v["resolved_share"], ("%.4f" % v["accuracy_among_resolved"]) if v["accuracy_among_resolved"] is not None else "n/a"))


if __name__ == "__main__":
    main()
