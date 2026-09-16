#!/usr/bin/env python3
"""Published comparator: RIPTiDe (Jenior et al. 2020) run in its native formulation.

RIPTiDe takes gene-level abundances, derives its own reaction weights through the model's
gene-protein-reaction rules, minimizes weighted flux subject to a minimum objective fraction
(default 0.8, far looser than the study's near-exact cap), prunes reactions that carry no flux,
and samples the remaining space. Its exchange prediction here is the sample median, with the
2.5th and 97.5th sample percentiles as the reported range. Nothing in the study's own objective,
cap or readout is reused, which is the point: it tests whether a published loose-cap sampling
workflow escapes the pinning that the study's readout imposes.

Inputs: the same medium, the same reconstruction and the same per-gene expression that fed the
study's reaction scores. The biomass objective is left free so RIPTiDe applies its own fraction.
No outcome is read.
"""
import sys
sys.dont_write_bytecode = True
import argparse, csv, hashlib, json, time, traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import numpy as np


def worker(job):
    root = Path(job["root"]); out = Path(job["out"]); sys.path.insert(0, str(root))
    import riptide, cobra, pandas as pd
    from src.cpu_model import build_model
    rec = {"run_id": job["run_id"], "context_id": job["context"], "arm": "riptide", "scenario": "primary", "status": "running"}
    t0 = time.monotonic()
    try:
        config = json.loads((root / "configs/study.json").read_text())
        model, meta = build_model(root, config, "primary")
        model.reactions.get_by_id(config["objective"]).bounds = (0.0, 1000.0)     # RIPTiDe applies its own fraction
        model.objective = config["objective"]
        d = np.load(root / "data/gene_expression.npz", allow_pickle=False)
        genes = d["genes"].astype(str).tolist(); ctxs = d["contexts"].astype(str).tolist()
        col = d["A"][:, ctxs.index(job["context"])]; expr = dict(zip(genes, col.tolist()))
        gmap = {r["model_gene_id"]: r["entrez_id"] for r in csv.DictReader(open(root / "manifests/model_gene_map.tsv"), delimiter="\t")}
        # RIPTiDe wants strictly positive abundances on the model's own gene identifiers
        floor = float(np.min(col[col > 0])) if (col > 0).any() else 1e-3
        transcriptome = {g.id: max(float(expr.get(gmap.get(g.id, ""), floor)), floor) for g in model.genes}
        r = riptide.contextualize(model, transcriptome=transcriptome, samples=job["samples"], fraction=job["fraction"], silent=True)
        fs = r.flux_samples
        pts, rng = {}, {}
        for rid in job["panel"]:
            if rid in fs.columns:
                v = fs[rid].to_numpy(float); pts[rid] = float(np.median(v)); rng[rid] = [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
            else:                                                     # pruned away: RIPTiDe found no flux through it
                pts[rid] = 0.0; rng[rid] = [0.0, 0.0]
        rec.update({"point": pts, "ranges": rng, "status": "optimal", "samples": int(fs.shape[0]), "reactions_retained": int(fs.shape[1]),
                    "pruned_exchanges": sum(1 for rid in job["panel"] if rid not in fs.columns),
                    "objective_fraction": job["fraction"], "riptide_seconds": getattr(r, "run_time", None)})
    except Exception as exc:
        rec.update({"status": "failed", "error": repr(exc)[:400], "traceback": traceback.format_exc()[-1500:]})
    rec["seconds"] = time.monotonic() - t0
    (out / "arms").mkdir(parents=True, exist_ok=True)
    (out / "arms" / ("%s.json" % hashlib.sha256(("%s|riptide" % job["context"]).encode()).hexdigest()[:16])).write_text(json.dumps(rec, indent=1, allow_nan=False))
    return rec


def run_single(a):
    """Run one context in this process and write its arm file (used inside an isolated subprocess)."""
    plan = json.loads((a.root / "protocol/analysis_plan.json").read_text())
    panel = sorted(t["exchange_id"] for t in plan["chemistry_targets"])
    job = {"root": str(a.root), "out": str(a.out), "context": a.single, "panel": panel, "run_id": a.run_id, "samples": a.samples, "fraction": a.fraction}
    r = worker(job); print(r["status"], flush=True); return 0 if r["status"] == "optimal" else 1


def arm_path(out, context):
    return Path(out) / "arms" / ("%s.json" % hashlib.sha256(("%s|riptide" % context).encode()).hexdigest()[:16])


def launch(job):
    """Each context runs in its own interpreter so a solver-level abort (a C assertion inside GLPK kills the
    process rather than raising) cannot take the other contexts down with it. Failures are recorded, not retried."""
    import subprocess
    t0 = time.time()
    cmd = [sys.executable, "-W", "ignore", __file__, "--root", job["root"], "--out", job["out"], "--single", job["context"], "--run-id", job["run_id"],
           "--samples", str(job["samples"]), "--fraction", str(job["fraction"])]
    p = subprocess.run(cmd, capture_output=True, text=True)
    ap_ = arm_path(job["out"], job["context"])
    if ap_.exists():
        rec = json.loads(ap_.read_text())
        if rec.get("run_id") == job["run_id"] or rec.get("status") == "optimal":
            return rec
    rec = {"run_id": job["run_id"], "context_id": job["context"], "arm": "riptide", "scenario": "primary", "status": "failed",
           "error": "subprocess exit %d: %s" % (p.returncode, (p.stderr or "").strip().splitlines()[-1][:300] if (p.stderr or "").strip() else "no stderr"),
           "traceback": (p.stderr or "")[-1500:], "seconds": time.time() - t0}
    ap_.parent.mkdir(parents=True, exist_ok=True); ap_.write_text(json.dumps(rec, indent=1, allow_nan=False))
    return rec


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", type=Path, required=True); ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=12); ap.add_argument("--samples", type=int, default=500); ap.add_argument("--fraction", type=float, default=0.8)
    ap.add_argument("--profiles", type=int, default=0); ap.add_argument("--single", default=None); ap.add_argument("--run-id", default=None)
    a = ap.parse_args(); a.out.mkdir(parents=True, exist_ok=True)
    if a.single:
        sys.exit(run_single(a))
    from concurrent.futures import ThreadPoolExecutor
    plan = json.loads((a.root / "protocol/analysis_plan.json").read_text())
    panel = sorted(t["exchange_id"] for t in plan["chemistry_targets"])
    ctxs = sorted(r["context_id"] for r in csv.DictReader(open(a.root / "manifests/contexts.tsv"), delimiter="\t") if r["partition"] == "test")
    if a.profiles: ctxs = ctxs[:a.profiles]
    manifest_path = a.out / "run_manifest.json"
    if manifest_path.exists() and a.run_id is None:                     # resume: keep the original run identifier
        run_id = json.loads(manifest_path.read_text())["run_id"]
    else:
        run_id = a.run_id or ("riptide-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"))
        manifest_path.write_text(json.dumps({"run_id": run_id, "jobs": len(ctxs), "samples": a.samples, "fraction": a.fraction,
            "started_at_utc": datetime.now(timezone.utc).isoformat(), "outcomes_read": False}, indent=1))
    done = {}
    for c in ctxs:
        p = arm_path(a.out, c)
        if p.exists():
            rec = json.loads(p.read_text())
            if rec.get("status") == "optimal": done[c] = rec
    todo = [c for c in ctxs if c not in done]
    print("resume: %d already optimal, %d to run, isolated subprocess per context" % (len(done), len(todo)), flush=True)
    jobs = [{"root": str(a.root), "out": str(a.out), "context": c, "panel": panel, "run_id": run_id, "samples": a.samples, "fraction": a.fraction} for c in todo]
    t0 = time.time(); res = list(done.values())
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(launch, j) for j in jobs]):
            r = f.result(); res.append(r); print("%2d/%d %-14s %s %.0fs" % (len(res), len(ctxs), r["context_id"], r["status"], time.time() - t0), flush=True)
            if r["status"] != "optimal": print("   ", str(r.get("error"))[:300], flush=True)
    res.sort(key=lambda r: r["context_id"])
    with open(a.out / "predictions.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(["run_id", "context_id", "arm", "exchange_id", "point", "range_lo", "range_hi", "scenario", "status"])
        for r in res:
            for rid in panel:
                pt = r.get("point", {}).get(rid, ""); rn = r.get("ranges", {}).get(rid, ["", ""]); w.writerow([r["run_id"], r["context_id"], "riptide", rid, pt, *rn, "primary", r["status"]])
    json.dump({"completed": sum(r["status"] == "optimal" for r in res), "failed": sum(r["status"] != "optimal" for r in res), "seconds": round(time.time() - t0, 1),
               "failed_contexts": sorted(r["context_id"] for r in res if r["status"] != "optimal"),
               "median_reactions_retained": float(np.median([r["reactions_retained"] for r in res if r["status"] == "optimal"])) if any(r["status"] == "optimal" for r in res) else None},
              open(a.out / "summary.json", "w"), indent=1)
    print("FINISHED", flush=True)


if __name__ == "__main__":
    main()
