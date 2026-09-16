#!/usr/bin/env python3
"""Third growth task for the independent evaluation: condition independent, at the mean measured growth rate.

The two declared tasks differ in two ways at once: the fixed-fraction task is identical across conditions but
sits at about 0.9 of the model maximum (order 3 units), whereas the measured task differs between conditions and
sits at 0.018 to 0.027 per hour. This task removes the operating-point difference: the biomass reaction is fixed
at the grand mean of the four condition-mean measured growth rates in every condition, so the task carries no
condition information and the transcripts are the only input that differs between conditions, at the same
operating point as the measured task. Everything else (medium, evidence, encodings, readout, tolerances) is the
frozen pipeline of copeland_run.py, whose functions are reused unchanged. No measured glucose or lactate value
is read here.
"""
import argparse, csv, json, sys, time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from copeland_run import one_arm, digest, ARMS, TARGETS


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--study-root", type=Path, required=True); ap.add_argument("--prepared", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True); ap.add_argument("--workers", type=int, default=12)
    a = ap.parse_args(); code = Path(__file__).resolve().parent; a.out.mkdir(parents=True, exist_ok=True)
    config = json.loads((a.study_root / "configs/study.json").read_text())
    prep = json.loads((a.prepared / "prepare_summary.json").read_text())
    samples = list(csv.DictReader(open(a.prepared / "sample_table.tsv"), delimiter="\t"))
    growth = list(csv.DictReader(open(a.prepared / "growth_rates.tsv"), delimiter="\t"))
    mu = {}
    for r in growth: mu.setdefault((r["oxygen"], r["treatment"]), []).append(float(r["mu"]))
    mu_mean = {k: float(np.mean(v)) for k, v in mu.items()}
    grand = float(np.mean(list(mu_mean.values())))
    jobs = []
    for s in samples:
        cond = (s["oxygen"], s["treatment"])
        for arm in ARMS:
            jobs.append({"study_root": str(a.study_root), "code": str(code), "output": str(a.out), "context": s["sample"], "arm": arm,
                         "task_mode": "mean_measured", "task": {"mode": "measured", "mu": grand}, "condition": list(cond),
                         "config": config, "k": prep["scale_constant_k"], "scores": str(a.prepared / "scores_copeland.npz")})
    (a.out / "run_manifest.json").write_text(json.dumps({"started_at_utc": datetime.now(timezone.utc).isoformat(), "jobs": len(jobs),
        "task": "biomass fixed at the grand mean of the four condition-mean measured growth rates, identical across conditions",
        "grand_mean_growth_rate_per_hour": grand, "condition_mean_growth_rate_per_hour": {"|".join(k): v for k, v in mu_mean.items()},
        "runner_sha256": digest(code / "copeland_run.py"), "target_outcomes_read": False}, indent=1) + "\n")
    t0 = time.time(); results = []
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(one_arm, j) for j in jobs]):
            r = f.result(); results.append(r)
            if len(results) % 14 == 0 or r["status"] != "optimal":
                print("%3d/%d %-4s %-16s %s %.0fs" % (len(results), len(jobs), r["context_id"], r["arm"], r["status"], time.time() - t0), flush=True)
    rows = []
    for r in results:
        if r["status"] != "optimal": continue
        for rid in TARGETS:
            lo, hi = r["ranges"][rid]
            rows.append([r["context_id"], r["condition"][0], r["condition"][1], r["arm"], r["task_mode"], rid, r["point"][rid], lo, hi,
                         r["primary_optimum"], r["model"]["task_flux"], r["model"]["g_max"]])
    with open(a.out / "predictions.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["context_id", "oxygen", "treatment", "arm", "task_mode", "exchange_id", "point", "range_lo", "range_hi", "primary_optimum", "task_flux", "g_max"])
        w.writerows(rows)
    summary = {"completed": sum(r["status"] == "optimal" for r in results), "failed": sum(r["status"] != "optimal" for r in results), "jobs": len(jobs),
               "seconds": round(time.time() - t0, 1), "grand_mean_growth_rate_per_hour": grand,
               "failures": [{"context": r["context_id"], "arm": r["arm"], "error": r.get("error")} for r in results if r["status"] != "optimal"]}
    (a.out / "run_summary.json").write_text(json.dumps(summary, indent=1) + "\n"); print(json.dumps(summary, indent=1)[:1500]); print("FINISHED", flush=True)


if __name__ == "__main__":
    main()
