#!/usr/bin/env python3
"""Summarize a nested-set ladder ledger (range_stage_ledger.tsv) into the two quantities the manuscript tables use:
the stage at which each profile-target coordinate first became fixed (Table 8 in the manuscript) and the position of the
surviving B2 intervals across profiles (Table 11): free targets, targets identical across all profiles, and exploration,
the spread of B2 midpoints across profiles divided by the B0 width. Run on the primary ledger it must reproduce the
published counts; run on the scenario ledgers it gives the sensitivity rows. No outcome is read."""
import sys, csv, json, argparse
from collections import defaultdict
from pathlib import Path
import numpy as np
EPS = 1e-5

ap = argparse.ArgumentParser(); ap.add_argument("--ledger", type=Path, required=True); ap.add_argument("--out", type=Path, required=True); ap.add_argument("--label", required=True)
a = ap.parse_args()
rows = list(csv.DictReader(open(a.ledger), delimiter="\t"))
arms = sorted({r["arm"] for r in rows}); ctxs = sorted({r["context"] for r in rows}); targets = sorted({r["exchange_id"] for r in rows})
out = {"label": a.label, "rows": len(rows), "profiles": len(ctxs), "targets": len(targets), "arms": {}}
for arm in arms:
    sub = [r for r in rows if r["arm"] == arm]
    stages = defaultdict(int)
    for r in sub: stages[r["first_fixed_stage"]] += 1
    b0w = {}; mids = defaultdict(list); b1w = defaultdict(list)
    for r in sub:
        b0w[r["exchange_id"]] = float(r["B0_width"])
        mids[r["exchange_id"]].append(0.5 * (float(r["B2_low"]) + float(r["B2_high"])))
        b1w[r["exchange_id"]].append(float(r["B1_width"]))
    free = [t for t in targets if b0w[t] > EPS]
    expl = {t: (max(mids[t]) - min(mids[t])) / b0w[t] for t in free}
    ident = sum(1 for t in free if max(mids[t]) - min(mids[t]) <= EPS)
    n_free_coords = len(sub) - stages["already_fixed_B0"]
    out["arms"][arm] = {"first_fixed": dict(stages), "free_coordinates": n_free_coords,
                        "share_of_free_first_fixed_B1": stages["first_fixed_B1"] / n_free_coords if n_free_coords else None,
                        "free_targets": len(free), "identical_across_profiles": ident,
                        "exploration_median": float(np.median(list(expl.values()))), "exploration_mean": float(np.mean(list(expl.values()))),
                        "mean_B0_width_free": float(np.mean([b0w[t] for t in free])),
                        "targets_fixed_by_B0": len(targets) - len(free)}
a.out.mkdir(parents=True, exist_ok=True)
json.dump(out, open(a.out / ("ladder_summary_%s.json" % a.label), "w"), indent=1)
for arm, v in out["arms"].items():
    print("%-16s B0 %4d  B1 %4d  B2 %4d  var %3d | share B1 %.3f | free %d ident %d | expl med %.2e mean %.4f" % (
        arm, v["first_fixed"].get("already_fixed_B0", 0), v["first_fixed"].get("first_fixed_B1", 0), v["first_fixed"].get("first_fixed_B2", 0),
        v["first_fixed"].get("still_variable_B2", 0), v["share_of_free_first_fixed_B1"], v["free_targets"], v["identical_across_profiles"], v["exploration_median"], v["exploration_mean"]))
