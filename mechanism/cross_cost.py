#!/usr/bin/env python3
"""Q2: do the encodings select different flux vectors, and does the difference survive the
projection onto the reported exchange coordinates?

If each encoding selects its own interior flux distribution but every one of them projects onto
almost the same exchange values, then the reported predictions are governed by the constraints
rather than by the transcript weights, and no choice of encoding can recover between-profile
information that the projection has already removed.

Two measurements per profile, both from saved vectors, no LP:

  cross-cost regret   score arm a's chosen vector against arm b's cost, relative to arm b's own
                      optimum. Small regret means the arms would nearly accept each other's
                      answers; large regret means they genuinely disagree about what is cheap.

  paired distance     relative L1 distance between two arms' full 10,600-reaction vectors, and
                      between their projections onto the 96 reported exchanges.

The self-cost reconstruction check runs first: scoring an arm against its own cost must return
its own recorded optimum, or the cost vectors are not the ones that produced these solutions.
"""
import csv, glob, hashlib, json, os, sys, time
from collections import defaultdict
import numpy as np
from mechanism import linear_cost, regret, relative_l1

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
ARMS = os.path.join(BASE, "cpu_study_2026-09-13", "results", "reserved_v3", "arms")
PANEL = os.path.join(BASE, "cpu_study_2026-09-13", "protocol", "analysis_plan.json")
SELF_TOL = 1e-6            # relative agreement required of the self-cost reconstruction
ACTIVE_TOL = 1e-6          # a reaction counts as differing above this absolute flux difference


def main():
    t0 = time.time()
    panel = [t["exchange_id"] for t in json.load(open(PANEL))["chemistry_targets"]]

    jobs = defaultdict(dict)                       # (scenario, context) -> arm -> record
    shared_uniform = {}                            # scenario -> record (context independent)
    ids = None
    files = sorted(glob.glob(os.path.join(ARMS, "*.json")))
    for i, f in enumerate(files):
        d = json.load(open(f))
        z = np.load(f[:-5] + ".npz")
        if ids is None:
            ids = [str(x) for x in z["reaction_ids"]]
            pidx = np.array([ids.index(e) for e in panel])
        rec = {"flux": z["flux"].astype(float), "cost": z["reaction_costs"].astype(float),
               "primary_optimum": float(d["primary_optimum"]), "primary_cap": float(d["primary_cap"]),
               "secondary_optimum": float(d["secondary_optimum"]),
               "secondary_cap": float(d["secondary_cap"]), "arm": d["arm"],
               "context": d["context_id"], "scenario": d["scenario"]}
        if d["arm"] == "uniform_pfba":
            shared_uniform[d["scenario"]] = rec
        else:
            jobs[(d["scenario"], d["context_id"])][d["arm"]] = rec
        if (i + 1) % 200 == 0:
            print("  loaded %d/%d  %.0fs" % (i + 1, len(files), time.time() - t0), flush=True)

    # ---- self-cost reconstruction, every arm, before anything is compared
    self_rows, worst_self = [], 0.0
    for key in list(jobs) + [("__uniform__", s) for s in shared_uniform]:
        pool = shared_uniform[key[1]] if key[0] == "__uniform__" else None
        recs = [pool] if pool else list(jobs[key].values())
        for r in recs:
            e = linear_cost(r["flux"], r["cost"])
            rel = (e - r["primary_optimum"]) / r["primary_optimum"]
            worst_self = max(worst_self, abs(rel))
            self_rows.append([r["scenario"], r["context"], r["arm"], r["primary_optimum"], e, rel])
    print("self-cost reconstruction: %d arms, worst relative difference %.3e"
          % (len(self_rows), worst_self), flush=True)
    if worst_self > SELF_TOL:
        print("ABORT: saved costs do not reproduce recorded optima"); sys.exit(2)

    # ---- pairwise regret and distance
    pair_rows, dist_rows = [], []
    for (scen, ctx), arms in sorted(jobs.items()):
        names = sorted(arms)
        for a in names:
            va = arms[a]["flux"]
            for b in names:
                rb = arms[b]
                cost_ab = linear_cost(va, rb["cost"])
                sec_ab = float(np.abs(va).sum())
                pair_rows.append([scen, ctx, a, b,
                                  regret(va, rb["cost"], rb["primary_optimum"]),
                                  int(cost_ab <= rb["primary_cap"]),
                                  int(sec_ab <= rb["secondary_cap"]),
                                  cost_ab, rb["primary_optimum"], rb["primary_cap"]])
            for b in names:
                if b <= a:
                    continue
                vb = arms[b]["flux"]
                diff = np.abs(va - vb)
                dist_rows.append([scen, ctx, a, b,
                                  relative_l1(va, vb),
                                  relative_l1(va[pidx], vb[pidx]),
                                  int((diff > ACTIVE_TOL).sum()),
                                  int((np.abs(va[pidx] - vb[pidx]) > ACTIVE_TOL).sum()),
                                  float(diff.max()), float(np.abs(va[pidx] - vb[pidx]).max())])

    with open("cross_cost_regret.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["scenario", "context", "vector_from_arm", "scored_against_arm", "regret",
                    "within_primary_cap", "within_secondary_cap", "cost_under_target",
                    "target_optimum", "target_cap"])
        w.writerows(pair_rows)
    with open("cross_cost_distance.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["scenario", "context", "arm_a", "arm_b", "relative_l1_full_vector",
                    "relative_l1_panel_exchanges", "reactions_differing", "panel_exchanges_differing",
                    "max_abs_difference_full", "max_abs_difference_panel"])
        w.writerows(dist_rows)
    with open("cross_cost_selfcheck.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["scenario", "context", "arm", "recorded_optimum", "recomputed", "relative"])
        w.writerows(self_rows)

    # ---- summaries
    summary = {"arms_loaded": len(files), "panel_exchanges": len(panel),
               "reactions_in_model": len(ids), "self_cost_worst_relative": worst_self,
               "active_tolerance": ACTIVE_TOL, "seconds": round(time.time() - t0, 1),
               "no_lp_solved": True, "no_outcome_read": True}
    for scen in sorted({r[0] for r in pair_rows}):
        off = [r for r in pair_rows if r[0] == scen and r[2] != r[3]]
        d = [r for r in dist_rows if r[0] == scen]
        summary[scen] = {
            "off_diagonal_comparisons": len(off),
            "regret_median": float(np.median([r[4] for r in off])),
            "regret_mean": float(np.mean([r[4] for r in off])),
            "regret_p95": float(np.percentile([r[4] for r in off], 95)),
            "regret_max": float(np.max([r[4] for r in off])),
            "share_within_the_other_arms_primary_cap": float(np.mean([r[5] for r in off])),
            "share_within_the_other_arms_secondary_cap": float(np.mean([r[6] for r in off])),
            "pairs_compared": len(d),
            "relative_l1_full_median": float(np.median([r[4] for r in d])),
            "relative_l1_full_mean": float(np.mean([r[4] for r in d])),
            "relative_l1_panel_median": float(np.median([r[5] for r in d])),
            "relative_l1_panel_mean": float(np.mean([r[5] for r in d])),
            "reactions_differing_median": float(np.median([r[6] for r in d])),
            "panel_exchanges_differing_median": float(np.median([r[7] for r in d])),
            "panel_exchanges_differing_mean": float(np.mean([r[7] for r in d])),
        }
    json.dump(summary, open("cross_cost.json", "w"), indent=1, default=float)
    print(json.dumps(summary, indent=1, default=float))


if __name__ == "__main__":
    main()
