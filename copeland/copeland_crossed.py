#!/usr/bin/env python3
"""Crossed control: which input actually moves the prediction, the transcripts or the task?

In the condition analysis each transcript profile is paired with its own condition's measured
growth rate, so the two inputs move together and their contributions cannot be separated. This
script breaks that pairing. Every one of the sixteen transcript profiles is run against every
one of the four condition growth rates, giving a full 16 by 4 factorial for each encoding.

If the prediction is governed by the transcripts, it should vary across profiles at a fixed
growth rate. If it is governed by the imposed task, it should vary across growth rates at a
fixed profile and barely at all across profiles. The two spreads are measured directly.

No measured glucose or lactate value is read here.
"""
import argparse, csv, json, sys, time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from copeland_run import ARMS, TARGETS, one_arm, digest


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
    samples = list(csv.DictReader(open(a.prepared / "sample_table.tsv"), delimiter="\t"))
    growth = list(csv.DictReader(open(a.prepared / "growth_rates.tsv"), delimiter="\t"))
    mu = {}
    for r in growth:
        mu.setdefault((r["oxygen"], r["treatment"]), []).append(float(r["mu"]))
    mu_mean = {k: float(np.mean(v)) for k, v in sorted(mu.items())}

    jobs = []
    for s in samples:
        for cond, value in mu_mean.items():
            for arm in ARMS:
                jobs.append({"study_root": str(a.study_root), "code": str(code),
                             "output": str(a.out), "context": s["sample"], "arm": arm,
                             "task_mode": "mu_%s_%s" % cond, "condition": list(cond),
                             "task": {"mode": "measured", "mu": value},
                             "config": config, "k": prep["scale_constant_k"],
                             "scores": str(a.prepared / "scores_copeland.npz")})
    (a.out / "run_manifest.json").write_text(json.dumps({
        "started_at_utc": datetime.now(timezone.utc).isoformat(), "jobs": len(jobs),
        "design": "16 transcript profiles crossed with 4 condition growth rates, all arms",
        "condition_mean_growth_rate_per_hour": {"|".join(k): v for k, v in mu_mean.items()},
        "target_outcomes_read": False,
        "crossed_runner_sha256": digest(Path(__file__).resolve())}, indent=1) + "\n")

    t0 = time.time(); results = []
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futs = [pool.submit(one_arm, j) for j in jobs]
        for f in as_completed(futs):
            r = f.result(); results.append(r)
            if len(results) % 60 == 0 or r["status"] != "optimal":
                print("%3d/%d %s %.0fs" % (len(results), len(jobs), r["status"],
                                           time.time() - t0), flush=True)

    rows = []
    for r in results:
        if r["status"] != "optimal":
            continue
        for rid in TARGETS:
            lo, hi = r["ranges"][rid]
            rows.append([r["context_id"], r["task_mode"], r["arm"], rid, r["point"][rid], lo, hi,
                         r["model"]["task_flux"]])
    with open(a.out / "crossed_predictions.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["context_id", "task_label", "arm", "exchange_id", "point",
                    "range_lo", "range_hi", "task_flux"])
        w.writerows(rows)

    # how much of the prediction moves with the task, and how much with the transcripts
    decomposition = {}
    for arm in sorted({r[2] for r in rows}):
        for rid in TARGETS:
            sub = [r for r in rows if r[2] == arm and r[3] == rid]
            if not sub:
                continue
            by_profile, by_task = {}, {}
            for c, t, _, _, p, _, _, _ in sub:
                by_profile.setdefault(c, []).append(p)     # one profile, four growth rates
                by_task.setdefault(t, []).append(p)        # one growth rate, sixteen profiles
            across_task = float(np.mean([max(v) - min(v) for v in by_profile.values()]))
            across_profile = float(np.mean([max(v) - min(v) for v in by_task.values()]))
            points = [r[4] for r in sub]
            total = float(max(points) - min(points))
            decomposition["%s|%s" % (arm, rid)] = {
                "mean_spread_across_growth_rates_at_fixed_transcript": across_task,
                "mean_spread_across_transcripts_at_fixed_growth_rate": across_profile,
                "total_spread": total,
                "ratio_transcript_to_task": (across_profile / across_task
                                             if across_task > 0 else None)}
    summary = {"completed": sum(r["status"] == "optimal" for r in results),
               "failed": sum(r["status"] != "optimal" for r in results), "jobs": len(jobs),
               "seconds": round(time.time() - t0, 1), "decomposition": decomposition,
               "finished_at_utc": datetime.now(timezone.utc).isoformat()}
    (a.out / "crossed_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print("\n%-16s %-14s %16s %16s %10s" % ("arm", "target", "spread by task",
                                            "spread by transcript", "ratio"))
    for k, v in decomposition.items():
        arm, rid = k.split("|")
        rr = v["ratio_transcript_to_task"]
        print("%-16s %-14s %16.6f %16.6f %10s" % (arm, rid,
              v["mean_spread_across_growth_rates_at_fixed_transcript"],
              v["mean_spread_across_transcripts_at_fixed_growth_rate"],
              ("%.4f" % rr) if rr is not None else "n/a"))


if __name__ == "__main__":
    main()
