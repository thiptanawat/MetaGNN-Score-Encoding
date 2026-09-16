#!/usr/bin/env python3
"""Two-stage readout at loosened tolerances: the complete study readout with the caps relaxed.

The single-stage sweep (cap_sweep.py) loosened the cost cap alone and ranged the targets under it, which leaves
open whether the study's two-stage readout, with the parsimony cap and the lexicographic point selection in
place, still pins or still returns an informative point at a loose tolerance. This runs the study's own
select_flux, unchanged, with the relative allowances of the two caps set to r for
  (a) the cost cap loosened to r and the parsimony cap kept at the study's near-exact value, and
  (b) both caps loosened to r,
for r of 1 and 20 percent, for the identity and midrank encodings on all 47 reserved profiles and for the
uniform cost once. Predictions and ranges are written in the study's format so comparator_eval.py scores the
point rule and the interval rule on the locked eligible pairs. No outcome is read.
"""
import sys
sys.dont_write_bytecode = True
import argparse, csv, hashlib, json, time, traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import numpy as np

EPS = 1e-5


def worker(job):
    root = Path(job["root"]); out = Path(job["out"]); sys.path.insert(0, str(root))
    import cobra
    from src.cpu_model import build_model, score_weights
    from src.cpu_readout import select_flux
    rec = {"run_id": job["run_id"], "context_id": job["context"], "arm": job["arm"], "scenario": "primary", "config": job["label"], "status": "running"}
    t0 = time.monotonic()
    try:
        config = json.loads((root / "configs/study.json").read_text())
        config = dict(config); config["cost_relative_allowance"] = job["r_cost"]; config["secondary_relative_allowance"] = job["r_parsimony"]
        model, meta = build_model(root, config, "primary")
        if job["arm"] == "uniform_pfba":
            w = {r.id: (0.0 if r.boundary else 1.0) for r in model.reactions}
        else:
            d = np.load(root / "data/scores_conservative.npz", allow_pickle=False)
            ctxs = d["contexts"].astype(str).tolist()
            w = score_weights(model, d["rxn"].astype(str).tolist(), d["A"][:, ctxs.index(job["context"])], job["k"], job["arm"], config["epsilon"])
        res = select_flux(model, w, job["panel"], config)
        res.pop("flux", None); res.pop("reaction_ids", None)
        rec.update(res); rec["status"] = "optimal"
    except Exception as exc:
        rec["status"] = "failed"; rec["error"] = repr(exc)[:400]; rec["traceback"] = traceback.format_exc()[-1500:]
    rec["seconds"] = time.monotonic() - t0
    (out / "arms").mkdir(parents=True, exist_ok=True)
    (out / "arms" / ("%s.json" % job["token"])).write_text(json.dumps(rec, indent=1, allow_nan=False))
    return rec


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", type=Path, required=True); ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=20); ap.add_argument("--r", default="0.01,0.20"); ap.add_argument("--arms", default="magnitude_g1.00,ordinal")
    ap.add_argument("--profiles", type=int, default=0)
    a = ap.parse_args(); a.out.mkdir(parents=True, exist_ok=True)
    plan = json.loads((a.root / "protocol/analysis_plan.json").read_text())
    panel = sorted(t["exchange_id"] for t in plan["chemistry_targets"])
    ctxs = sorted(r["context_id"] for r in csv.DictReader(open(a.root / "manifests/contexts.tsv"), delimiter="\t") if r["partition"] == "test")
    if a.profiles: ctxs = ctxs[:a.profiles]
    k = float(json.loads((a.root / "results/development_v3/run_summary_scalars.json").read_text())["k"])
    exact = json.loads((a.root / "configs/study.json").read_text())["secondary_relative_allowance"]
    run_id = "twostage-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    jobs = []
    for r in [float(x) for x in a.r.split(",")]:
        for label, rc, rp in (("cost_%g_parsimony_exact" % r, r, exact), ("both_%g" % r, r, r)):
            for arm in a.arms.split(",") + ["uniform_pfba"]:
                for c in (ctxs if arm != "uniform_pfba" else ["__reference__"]):
                    jobs.append({"root": str(a.root), "out": str(a.out / label), "context": c, "arm": arm, "k": k, "panel": panel, "run_id": run_id,
                                 "label": label, "r_cost": rc, "r_parsimony": rp, "token": hashlib.sha256(("%s|%s|%s" % (c, arm, label)).encode()).hexdigest()[:16]})
    (a.out / "run_manifest.json").write_text(json.dumps({"run_id": run_id, "jobs": len(jobs), "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "configs": sorted({j["label"] for j in jobs}), "outcomes_read": False, "readout": "src.cpu_readout.select_flux unchanged; only the two relative cap allowances differ"}, indent=1))
    t0 = time.time(); res = []
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(worker, j) for j in jobs]):
            r = f.result(); res.append(r)
            if len(res) % 20 == 0 or r["status"] != "optimal":
                print("%3d/%d %-14s %-16s %-28s %s %.0fs" % (len(res), len(jobs), r["context_id"], r["arm"], r["config"], r["status"], time.time() - t0), flush=True)
                if r["status"] != "optimal": print("   ", r.get("error", "")[:200], flush=True)
    for label in sorted({r["config"] for r in res}):
        rows = [r for r in res if r["config"] == label]; rows.sort(key=lambda r: (r["arm"], r["context_id"]))
        with open(a.out / label / "predictions.tsv", "w", newline="") as f:
            w = csv.writer(f, delimiter="\t"); w.writerow(["run_id", "context_id", "arm", "exchange_id", "point", "range_lo", "range_hi", "scenario", "status"])
            for r in rows:
                for rid in panel:
                    pt = r.get("point", {}).get(rid, ""); rn = r.get("ranges", {}).get(rid, ["", ""])
                    w.writerow([r["run_id"], r["context_id"], r["arm"], rid, pt, *rn, "primary", r["status"]])
        json.dump({"completed": sum(r["status"] == "optimal" for r in rows), "failed": sum(r["status"] != "optimal" for r in rows), "arms": sorted({r["arm"] for r in rows})},
                  open(a.out / label / "summary.json", "w"), indent=1)
    print(json.dumps({"completed": sum(r["status"] == "optimal" for r in res), "failed": sum(r["status"] != "optimal" for r in res), "seconds": round(time.time() - t0, 1)}), flush=True)
    print("FINISHED", flush=True)


if __name__ == "__main__":
    main()
