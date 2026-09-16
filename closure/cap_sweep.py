#!/usr/bin/env python3
"""Cap-tolerance sweep: is the pinning of exchange coordinates a property of the near-exact cap?

The study caps the weighted-parsimony cost at its optimum plus 1e-7 absolute plus 1e-9 relative.
A cap that tight makes the admissible set close to the optimal face of the LP, where every
coordinate is fixed for generic costs. Published workflows operate far looser: RIPTiDe's default
objective fraction is 0.8 and GIMME uses a required-functionality threshold. This sweep re-ranges
the 52 reported targets under the B0 constraints plus the weighted cost capped at its optimum
times (1 + r) for r in {1e-9 (the study), 0.01, 0.05, 0.20}, for identity, midrank and uniform
costs across all 47 reserved profiles in the primary scenario.

No outcome is read. Each job builds the model once, solves the primary optimum once, and ranges
the panel under each cap in turn so that the four caps share an identical LP up to the cap.
"""
import sys
sys.dont_write_bytecode = True
import argparse, csv, json, time, traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import numpy as np

CAPS = [("exact", 1e-9, 1e-7), ("rel_1pct", 0.01, 0.0), ("rel_5pct", 0.05, 0.0), ("rel_20pct", 0.20, 0.0)]
FINE_CAPS = [("rel_1e-6", 1e-6, 0.0), ("rel_1e-5", 1e-5, 0.0), ("rel_1e-4", 1e-4, 0.0), ("rel_1e-3", 1e-3, 0.0)]
EPS = 1e-5
ARMS = ["magnitude_g1.00", "ordinal", "uniform_pfba"]


def worker(job):
    root = Path(job["root"]); out = Path(job["out"])
    sys.path.insert(0, str(root))
    from src.cpu_model import build_model, score_weights, set_absolute_cost
    rec = {"context": job["context"], "arm": job["arm"], "status": "running", "caps": {}}
    t0 = time.monotonic()
    try:
        config = json.loads((root / "configs/study.json").read_text())
        model, meta = build_model(root, config, "primary")
        d = np.load(root / "data/scores_conservative.npz", allow_pickle=False)
        ctxs = d["contexts"].astype(str).tolist()
        w = score_weights(model, d["rxn"].astype(str).tolist(), d["A"][:, ctxs.index(job["context"])],
                          job["k"], job["arm"], config["epsilon"])
        set_absolute_cost(model, w)
        opt = model.slim_optimize()
        if model.solver.status != "optimal":
            raise RuntimeError("primary solve %s" % model.solver.status)
        expr = model.objective.expression
        rec["primary_optimum"] = float(opt); rec["task_flux"] = meta["task_flux"]
        for name, rel, absr in job.get("caps", CAPS):
            cap = float(opt) + absr + abs(float(opt)) * rel
            con = model.problem.Constraint(expr, ub=cap, name="sweep_cap")
            model.add_cons_vars([con]); model.solver.update()
            ranges = {}
            for rid in job["panel"]:
                rx = model.reactions.get_by_id(rid)
                model.objective = model.problem.Objective(rx.flux_expression, direction="min")
                lo = model.slim_optimize()
                if model.solver.status != "optimal": raise RuntimeError("range_min %s %s" % (rid, model.solver.status))
                model.objective = model.problem.Objective(rx.flux_expression, direction="max")
                hi = model.slim_optimize()
                if model.solver.status != "optimal": raise RuntimeError("range_max %s %s" % (rid, model.solver.status))
                ranges[rid] = [float(lo), float(hi)]
            model.remove_cons_vars([con]); model.solver.update()
            rec["caps"][name] = {"cap": cap, "ranges": ranges,
                                 "fixed": sum(1 for v in ranges.values() if v[1] - v[0] <= EPS)}
        rec["status"] = "complete"
    except Exception as exc:
        rec["status"] = "failed"; rec["error"] = repr(exc); rec["traceback"] = traceback.format_exc()
    rec["seconds"] = time.monotonic() - t0
    (out / "jobs").mkdir(parents=True, exist_ok=True)
    (out / "jobs" / ("%s__%s.json" % (job["context"].replace(":", "_").replace("/", "_"), job["arm"]))
     ).write_text(json.dumps(rec, indent=1))
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True); ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=20); ap.add_argument("--fine", action="store_true", help="use the intermediate relative tolerances 1e-6 to 1e-3 instead of the coarse grid")
    a = ap.parse_args(); a.out.mkdir(parents=True, exist_ok=True)
    global CAPS
    if a.fine: CAPS = FINE_CAPS
    plan = json.loads((a.root / "protocol/analysis_plan.json").read_text())
    panel = sorted(t["exchange_id"] for t in plan["primary_targets"]); assert len(panel) == 52
    ctxs = sorted(r["context_id"] for r in csv.DictReader(open(a.root / "manifests/contexts.tsv"), delimiter="\t")
                  if r["partition"] == "test"); assert len(ctxs) == 47
    k = float(json.loads((a.root / "results/development_v3/run_summary_scalars.json").read_text())["k"])
    jobs = [{"root": str(a.root), "out": str(a.out), "context": c, "arm": arm, "k": k, "panel": panel, "caps": CAPS}
            for arm in ARMS for c in (ctxs if arm != "uniform_pfba" else ctxs[:1])]
    (a.out / "run_manifest.json").write_text(json.dumps({
        "started_at_utc": datetime.now(timezone.utc).isoformat(), "jobs": len(jobs), "arms": ARMS,
        "caps": CAPS, "profiles": ctxs, "targets": panel, "k": k, "outcomes_read": False}, indent=1))
    t0 = time.time(); res = []
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(worker, j) for j in jobs]):
            r = f.result(); res.append(r)
            if len(res) % 10 == 0 or r["status"] != "complete":
                print("%3d/%d %-14s %-16s %s %.0fs" % (len(res), len(jobs), r["context"], r["arm"], r["status"], time.time() - t0), flush=True)
    # summary: for each arm and cap, fixed count over 52 targets averaged over profiles, and the
    # across-profile location spread per target (the exploration numerator)
    summ = {}
    for arm in ARMS:
        rr = [r for r in res if r["arm"] == arm and r["status"] == "complete"]
        summ[arm] = {"profiles": len(rr)}
        for name, _, _ in CAPS:
            fixed = [r["caps"][name]["fixed"] for r in rr]
            spread = {}
            for rid in panel:
                mids = [0.5 * (r["caps"][name]["ranges"][rid][0] + r["caps"][name]["ranges"][rid][1]) for r in rr]
                spread[rid] = max(mids) - min(mids)
            widths = [r["caps"][name]["ranges"][rid][1] - r["caps"][name]["ranges"][rid][0] for r in rr for rid in panel]
            summ[arm][name] = {"mean_fixed_of_52": float(np.mean(fixed)), "min_fixed": int(min(fixed)), "max_fixed": int(max(fixed)),
                               "targets_constant_across_profiles": sum(1 for v in spread.values() if v <= EPS),
                               "median_width": float(np.median(widths)), "mean_width": float(np.mean(widths))}
    summ["_run"] = {"completed": sum(r["status"] == "complete" for r in res), "failed": sum(r["status"] != "complete" for r in res),
                    "seconds": round(time.time() - t0, 1), "finished_at_utc": datetime.now(timezone.utc).isoformat()}
    (a.out / "summary.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
