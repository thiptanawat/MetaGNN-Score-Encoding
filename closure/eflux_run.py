#!/usr/bin/env python3
"""Bound-based encoding: the falsification test the manuscript names.

In the study's objective formulation the transcript evidence enters as reaction costs, prices
internal reactions only, and leaves the reported exchanges unpriced. The manuscript states that
a formulation carrying the evidence in reaction bounds instead should, if the mechanism is
specific to objective costs, show substantially larger across-profile exploration of the
admissible set. This runs that formulation.

Bound-based encoding, following the E-Flux idea of Colijn et al. (2009) without claiming a
faithful reproduction of that method: for every supported internal reaction the magnitude of
each bound is scaled by the encoded score q in [0, 1), so ub = min(ub0, B q) and
lb = max(lb0, -B q), with B the reconstruction's default bound magnitude. A floor of B/1000 keeps
no reaction fully blocked, and is recorded as a declared deviation. Unsupported reactions keep
their default bounds, mirroring the conservative policy. The same six encodings map evidence to
q exactly as they map it to cost in the study. The task is 0.9 times the profile's own maximum
under its bounds, and the objective is uniform parsimony, so no transcript information enters
through the objective. The readout is the study's own select_flux with uniform costs, which
returns the point, the ranges under the parsimony cap, and the numerical audit. B0 for this
formulation, the profile's bounded admissible set with the task and no objective, is ranged
too, so exploration can be computed against it.

Predictions are written in the study's own format so the study's evaluator scores them
unchanged. No outcome is read here.
"""
import sys
sys.dont_write_bytecode = True
import argparse, csv, hashlib, json, time, traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import numpy as np

ARMS = ["magnitude_g0.50", "magnitude_g0.80", "magnitude_g1.00", "magnitude_g1.25", "magnitude_g2.00", "ordinal"]
B = 1000.0
FLOOR_FRACTION = 1e-3
EPS = 1e-5


def encoded_q(evidence, k, arm):
    from scipy.stats import rankdata
    z = np.asarray(evidence, float) / (np.asarray(evidence, float) + k)
    if arm == "ordinal":
        return (rankdata(z, method="average") - 0.5) / len(z)
    return z ** float(arm.removeprefix("magnitude_g"))


def worker(job):
    root = Path(job["root"]); out = Path(job["out"]); sys.path.insert(0, str(root))
    import cobra
    from src.cpu_model import build_model
    from src.cpu_readout import select_flux
    rec = {"run_id": job["run_id"], "context_id": job["context"], "arm": job["arm"], "scenario": "primary",
           "status": "running", "formulation": "bound_based"}
    t0 = time.monotonic()
    try:
        config = json.loads((root / "configs/study.json").read_text())
        model, meta = build_model(root, config, "primary")        # medium and provisional task
        d = np.load(root / "data/scores_conservative.npz", allow_pickle=False)
        rxns = d["rxn"].astype(str).tolist(); ctxs = d["contexts"].astype(str).tolist()
        q = encoded_q(d["A"][:, ctxs.index(job["context"])], job["k"], job["arm"])
        # release the provisional task, apply the bounds, then re-derive the task under them
        task = model.reactions.get_by_id(config["objective"])
        task.bounds = (0.0, 1000.0)
        changed = 0
        for rid, qq in zip(rxns, q):
            r = model.reactions.get_by_id(rid)
            if r.boundary: continue
            lim = max(job["B"] * float(qq), job["B"] * FLOOR_FRACTION)
            lo, hi = r.bounds
            r.bounds = (max(lo, -lim) if lo < 0 else lo, min(hi, lim) if hi > 0 else hi)
            changed += 1
        model.objective = config["objective"]
        gmax = model.slim_optimize()
        if model.solver.status != "optimal" or not np.isfinite(gmax) or gmax <= 0:
            raise RuntimeError("task maximum under bounds: %s %s" % (gmax, model.solver.status))
        if job.get("task_fixed") is not None:                # profile-independent task: the same absolute value in every profile
            g0 = float(job["task_fixed"])
            if g0 > float(gmax) + 1e-9:
                raise RuntimeError("fixed task %.6f exceeds this profile's maximum %.6f under its bounds" % (g0, gmax))
        else:
            g0 = float(gmax) * config["scenarios"]["primary"]["task_fraction"]
        task.bounds = (g0, g0)
        rec["model"] = {"g_max": float(gmax), "task_flux": g0, "bounded_reactions": changed, "B": job["B"], "task_rule": "fixed" if job.get("task_fixed") is not None else "0.9 x profile maximum"}
        # B0 for this formulation: bounds + task, no objective
        b0 = {}
        for rid in job["panel"]:
            rx = model.reactions.get_by_id(rid)
            model.objective = model.problem.Objective(rx.flux_expression, direction="min"); lo = model.slim_optimize()
            model.objective = model.problem.Objective(rx.flux_expression, direction="max"); hi = model.slim_optimize()
            b0[rid] = [float(lo), float(hi)]
        rec["b0_ranges"] = b0
        # uniform parsimony readout through the study's own selection routine
        uniform = {r.id: (0.0 if r.boundary else 1.0) for r in model.reactions}
        res = select_flux(model, uniform, job["panel"], config)
        flux = res.pop("flux"); ids = res.pop("reaction_ids")
        binding = 0; bounded = 0
        for rid, qq in zip(rxns, q):
            r = model.reactions.get_by_id(rid)
            if r.boundary: continue
            bounded += 1; v = float(flux[ids.index(rid)])
            if abs(abs(v) - abs(r.upper_bound if v > 0 else r.lower_bound)) <= 1e-6 and abs(v) > 1e-6: binding += 1
        rec["binding_bounds"] = binding; rec["bounded_reactions"] = bounded
        np.savez_compressed(out / "arms" / ("%s.npz" % job["token"]), flux=flux, reaction_ids=np.array(ids))
        rec.update(res); rec["status"] = "optimal"
    except Exception as exc:
        rec["status"] = "failed"; rec["error"] = repr(exc); rec["traceback"] = traceback.format_exc()
    rec["seconds"] = time.monotonic() - t0
    (out / "arms" / ("%s.json" % job["token"])).write_text(json.dumps(rec, indent=1, allow_nan=False))
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True); ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=20); ap.add_argument("--B", type=float, default=10.0)
    ap.add_argument("--profiles", type=int, default=0); ap.add_argument("--arms", default=",".join(ARMS))
    ap.add_argument("--task-fixed", type=float, default=None, help="profile-independent biomass flux applied to every profile (default: 0.9 x each profile's own maximum)")
    a = ap.parse_args(); (a.out / "arms").mkdir(parents=True, exist_ok=True)
    plan = json.loads((a.root / "protocol/analysis_plan.json").read_text())
    panel = sorted(t["exchange_id"] for t in plan["chemistry_targets"])      # the study reports all 96
    primary = sorted(t["exchange_id"] for t in plan["primary_targets"]); assert len(primary) == 52
    ctxs = sorted(r["context_id"] for r in csv.DictReader(open(a.root / "manifests/contexts.tsv"), delimiter="\t")
                  if r["partition"] == "test"); assert len(ctxs) == 47
    k = float(json.loads((a.root / "results/development_v3/run_summary_scalars.json").read_text())["k"])
    run_id = "bound-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    use_arms = a.arms.split(","); use_ctx = ctxs[:a.profiles] if a.profiles else ctxs
    jobs = [{"root": str(a.root), "out": str(a.out), "context": c, "arm": arm, "k": k, "panel": panel, "run_id": run_id, "B": a.B, "task_fixed": a.task_fixed,
             "token": hashlib.sha256(("%s|%s|bound" % (c, arm)).encode()).hexdigest()[:16]} for c in use_ctx for arm in use_arms]
    (a.out / "run_manifest.json").write_text(json.dumps({"run_id": run_id, "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "jobs": len(jobs), "arms": use_arms, "bound_magnitude": a.B, "floor_fraction": FLOOR_FRACTION, "task_rule": ("fixed at %.6f in every profile" % a.task_fixed) if a.task_fixed is not None else "0.9 x profile maximum under its bounds",
        "objective": "uniform parsimony; transcript information enters through bounds only", "k": k, "outcomes_read": False}, indent=1))
    t0 = time.time(); res = []
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(worker, j) for j in jobs]):
            r = f.result(); res.append(r)
            if len(res) % 20 == 0 or r["status"] != "optimal":
                print("%3d/%d %-14s %-16s %s %.0fs" % (len(res), len(jobs), r["context_id"], r["arm"], r["status"], time.time() - t0), flush=True)
                if r["status"] != "optimal": print("   ", r.get("error", "")[:200], flush=True)
    res.sort(key=lambda r: (r["context_id"], r["arm"]))
    with open(a.out / "predictions.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["run_id", "context_id", "arm", "exchange_id", "point", "range_lo", "range_hi", "scenario", "status"])
        for r in res:
            for rid in panel:
                pt = r.get("point", {}).get(rid, ""); rn = r.get("ranges", {}).get(rid, ["", ""])
                w.writerow([r["run_id"], r["context_id"], r["arm"], rid, pt, *rn, "primary", r["status"]])
    # exploration against the profile-specific B0
    summ = {"completed": sum(r["status"] == "optimal" for r in res), "failed": sum(r["status"] != "optimal" for r in res),
            "seconds": round(time.time() - t0, 1), "arms": {}}
    for arm in use_arms:
        rr = [r for r in res if r["arm"] == arm and r["status"] == "optimal"]
        if not rr: continue
        ex, const, fixed_b1, free = [], 0, 0, 0
        for rid in primary:
            b0w = float(np.mean([r["b0_ranges"][rid][1] - r["b0_ranges"][rid][0] for r in rr]))
            mids = [0.5 * (r["ranges"][rid][0] + r["ranges"][rid][1]) for r in rr]
            spread = max(mids) - min(mids)
            fixed_b1 += sum(1 for r in rr if r["ranges"][rid][1] - r["ranges"][rid][0] <= EPS)
            if b0w > EPS:
                free += 1; ex.append(spread / b0w)
                if spread <= EPS: const += 1
        summ["arms"][arm] = {"profiles": len(rr), "free_targets_mean_B0": free, "constant_across_profiles": const,
                             "exploration_median": float(np.median(ex)) if ex else None, "exploration_mean": float(np.mean(ex)) if ex else None,
                             "exploration_p90": float(np.percentile(ex, 90)) if ex else None,
                             "profile_target_fixed_under_parsimony_cap": fixed_b1, "of": len(rr) * len(primary),
                             "g_max_range": [min(r["model"]["g_max"] for r in rr), max(r["model"]["g_max"] for r in rr)],
                             "binding_bounds_mean": float(np.mean([r["binding_bounds"] for r in rr])), "bounded_reactions": rr[0]["bounded_reactions"]}
    (a.out / "summary.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
