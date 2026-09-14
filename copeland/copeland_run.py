#!/usr/bin/env python3
"""Run the frozen pipeline on the Copeland lung-fibroblast conditions.

The optimization interface is the study's own: the same GPR-supported reaction evidence, the
same score-to-cost map, the same sequential readout with pre-selection ranges and the same
numerical tolerances. Only the culture environment, the task value and the transcript source
differ, which is what makes this an independent evaluation rather than a re-run.

No measured glucose or lactate value is read here. Growth rate is a declared model input.
"""
import argparse, csv, hashlib, json, sys, time, traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import numpy as np

sys.dont_write_bytecode = True
TARGETS = ["EX_glc__D_e", "EX_lac__L_e"]
ARMS = ["magnitude_g0.50", "magnitude_g0.80", "magnitude_g1.00", "magnitude_g1.25",
        "magnitude_g2.00", "ordinal", "uniform_pfba"]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build_copeland_model(study_root, medium_module, config, task):
    """The pinned Recon3D under the MCDB 131 environment, with the task fixed as declared."""
    import cobra
    model = cobra.io.load_json_model(str(Path(study_root) / "data/raw/Recon3D.json"))
    model.solver = config["solver"]
    model.solver.configuration.tolerances.feasibility = config["solver_feasibility_tolerance"]
    model.solver.configuration._smcp.tol_dj = config["solver_optimality_tolerance"]
    model.solver.configuration.timeout = config["solver_timeout_seconds"]
    opened, missing, changes = medium_module.apply_medium(model)
    model.objective = config["objective"]
    gmax = model.slim_optimize()
    if model.solver.status != "optimal" or not np.isfinite(gmax) or gmax <= 0:
        raise RuntimeError("Invalid task maximum under the MCDB 131 medium: %s" % gmax)
    rx = model.reactions.get_by_id(config["objective"])
    if task["mode"] == "measured":
        value = float(task["mu"])
        if value > gmax:
            raise RuntimeError("Measured growth rate %.6f exceeds the model maximum %.6f; the "
                               "environment cannot support the observed growth" % (value, gmax))
    elif task["mode"] == "fraction":
        value = float(gmax) * float(task["fraction"])
    else:
        raise ValueError(task["mode"])
    rx.bounds = (value, value)
    return model, {"g_max": float(gmax), "task_flux": value, "objective": rx.id,
                   "objective_name": rx.name, "task_mode": task["mode"],
                   "opened": opened, "missing": missing, "bound_changes": changes}


def one_arm(job):
    root = Path(job["study_root"]); out = Path(job["output"])
    sys.path.insert(0, str(root)); sys.path.insert(0, str(Path(job["code"])))
    import copeland_medium
    from src.cpu_model import score_weights
    from src.cpu_readout import select_flux
    token = hashlib.sha256(("%s|%s|%s" % (job["context"], job["arm"], job["task_mode"]))
                           .encode()).hexdigest()[:16]
    dest = out / "arms" / ("%s.json" % token); dest.parent.mkdir(parents=True, exist_ok=True)
    rec = {"context_id": job["context"], "arm": job["arm"], "task_mode": job["task_mode"],
           "condition": job["condition"]}
    t0 = time.monotonic()
    try:
        d = np.load(job["scores"], allow_pickle=False)
        rxns = d["rxn"].astype(str).tolist(); ctxs = d["contexts"].astype(str).tolist()
        model, meta = build_copeland_model(root, copeland_medium, job["config"], job["task"])
        w = score_weights(model, rxns, d["A"][:, ctxs.index(job["context"])], job["k"],
                          job["arm"], job["config"]["epsilon"])
        res = select_flux(model, w, TARGETS, job["config"])
        np.savez_compressed(dest.with_suffix(".npz"), flux=res.pop("flux"),
                            reaction_ids=np.array(res.pop("reaction_ids")),
                            reaction_costs=np.array([w[r.id] for r in model.reactions]))
        rec.update(res); rec.update({"status": "optimal", "model": meta,
                                     "flux_file": "arms/%s.npz" % token,
                                     "flux_sha256": digest(dest.with_suffix(".npz"))})
    except Exception as exc:
        rec.update({"status": "failed", "error": repr(exc), "traceback": traceback.format_exc()})
    rec["seconds"] = time.monotonic() - t0
    dest.write_text(json.dumps(rec, indent=1, allow_nan=False) + "\n")
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--study-root", type=Path, required=True)
    ap.add_argument("--prepared", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=20)
    a = ap.parse_args()
    code = Path(__file__).resolve().parent
    a.out.mkdir(parents=True, exist_ok=True)
    config = json.loads((a.study_root / "configs/study.json").read_text())
    prep = json.loads((a.prepared / "prepare_summary.json").read_text())
    plan = json.loads((code / "copeland_plan.json").read_text())
    if plan["status"] != "frozen_before_any_prediction_was_computed":
        raise ValueError("The analysis plan is not frozen")

    samples = list(csv.DictReader(open(a.prepared / "sample_table.tsv"), delimiter="\t"))
    growth = list(csv.DictReader(open(a.prepared / "growth_rates.tsv"), delimiter="\t"))
    mu = {}
    for r in growth:
        mu.setdefault((r["oxygen"], r["treatment"]), []).append(float(r["mu"]))
    mu_mean = {k: float(np.mean(v)) for k, v in mu.items()}
    if len(mu_mean) != 4 or any(len(v) != 4 for v in mu.values()):
        raise ValueError("Expected four conditions with four replicates each")

    jobs = []
    for s in samples:
        cond = (s["oxygen"], s["treatment"])
        if cond not in mu_mean:
            raise ValueError("No measured growth rate for condition %s" % (cond,))
        for arm in ARMS:
            for mode, task in (("measured", {"mode": "measured", "mu": mu_mean[cond]}),
                               ("fraction", {"mode": "fraction",
                                             "fraction": config["scenarios"]["primary"]["task_fraction"]})):
                jobs.append({"study_root": str(a.study_root), "code": str(code),
                             "output": str(a.out), "context": s["sample"], "arm": arm,
                             "task_mode": mode, "task": task, "condition": list(cond),
                             "config": config, "k": prep["scale_constant_k"],
                             "scores": str(a.prepared / "scores_copeland.npz")})
    manifest = {"started_at_utc": datetime.now(timezone.utc).isoformat(),
                "jobs": len(jobs), "arms": ARMS, "targets": TARGETS,
                "condition_mean_growth_rate_per_hour": {"|".join(k): v for k, v in mu_mean.items()},
                "scale_constant_k": prep["scale_constant_k"],
                "plan_sha256": digest(code / "copeland_plan.json"),
                "medium_sha256": digest(code / "copeland_medium.py"),
                "runner_sha256": digest(Path(__file__).resolve()),
                "prepare_summary_sha256": digest(a.prepared / "prepare_summary.json"),
                "scores_sha256": digest(a.prepared / "scores_copeland.npz"),
                "target_outcomes_read": False}
    (a.out / "run_manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")

    t0 = time.time(); results = []
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futs = [pool.submit(one_arm, j) for j in jobs]
        for f in as_completed(futs):
            r = f.result(); results.append(r)
            if len(results) % 20 == 0 or r["status"] != "optimal":
                print("%3d/%d %-4s %-16s %-9s %s %.0fs" % (len(results), len(jobs), r["context_id"],
                      r["arm"], r["task_mode"], r["status"], time.time() - t0), flush=True)
                if r["status"] != "optimal":
                    print("   ", r.get("error", "")[:300], flush=True)

    rows = []
    for r in results:
        if r["status"] != "optimal":
            continue
        for rid in TARGETS:
            lo, hi = r["ranges"][rid]
            rows.append([r["context_id"], r["condition"][0], r["condition"][1], r["arm"],
                         r["task_mode"], rid, r["point"][rid], lo, hi,
                         r["primary_optimum"], r["model"]["task_flux"], r["model"]["g_max"]])
    with open(a.out / "predictions.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["context_id", "oxygen", "treatment", "arm", "task_mode", "exchange_id",
                    "point", "range_lo", "range_hi", "primary_optimum", "task_flux", "g_max"])
        w.writerows(rows)
    summary = {"completed": sum(r["status"] == "optimal" for r in results),
               "failed": sum(r["status"] != "optimal" for r in results),
               "jobs": len(jobs), "prediction_rows": len(rows),
               "seconds": round(time.time() - t0, 1),
               "failures": [{"context": r["context_id"], "arm": r["arm"],
                             "task_mode": r["task_mode"], "error": r.get("error")}
                            for r in results if r["status"] != "optimal"],
               "finished_at_utc": datetime.now(timezone.utc).isoformat()}
    (a.out / "run_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "failures"}, indent=1))
    if summary["failures"]:
        print("FAILURES:", json.dumps(summary["failures"][:5], indent=1)[:1500])


if __name__ == "__main__":
    main()
