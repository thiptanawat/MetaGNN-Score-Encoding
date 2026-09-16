#!/usr/bin/env python3
"""Cross-solver check: do the fixed/variable classifications and point values survive a change of LP solver?

The study used GLPK throughout, and its mechanism claims rest on distinctions at the 1e-7 to
1e-10 scale on degenerate LPs with near-exact caps, where solvers can land on different optimal
faces. This reruns the complete readout (primary optimum, caps, ranges, coordinate fixing and the
lexicographic point selection over the production list of 96 chemistry exchanges, in the same
order as src/run_study.py) for a representative subset of reserved arms with a second solver
and compares, per coordinate: the point value, the range width, the fixed classification at
1e-5, and the between-profile ordering signs among the subset. The first version of this script
(16 September 2026, before the external review) passed only the 52 primary targets to the
selection routine, which changes the lexicographic sequence; --panel primary reproduces that
variant and the default --panel chemistry is the production procedure. GLPK's exact rational-arithmetic variant is preferred when the
installed optlang exposes it, because it removes floating-point tolerance from the comparison;
otherwise the SciPy interface, which uses HiGHS, is used.

No outcome is read.
"""
import sys
sys.dont_write_bytecode = True
import argparse, csv, hashlib, json, time, traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np

EPS = 1e-5


def pick_solver():
    import optlang
    for name in ("glpk_exact", "scipy"):
        try:
            __import__("optlang.%s_interface" % name); return name
        except Exception: pass
    raise RuntimeError("no second solver available")


def worker(job):
    root = Path(job["root"]); sys.path.insert(0, str(root))
    import cobra
    from src.cpu_model import score_weights
    from src.cpu_readout import select_flux
    from src.medium_v2 import apply_medium
    rec = {"context": job["context"], "arm": job["arm"], "solver": job["solver"], "status": "running"}
    t0 = time.monotonic()
    try:
        config = json.loads((root / "configs/study.json").read_text())
        model = cobra.io.load_json_model(str(root / "data/raw/Recon3D.json"))
        model.solver = job["solver"]
        try: model.solver.configuration.tolerances.feasibility = config["solver_feasibility_tolerance"]
        except Exception: pass
        try: model.solver.configuration.timeout = 600
        except Exception: pass
        sc = config["scenarios"]["primary"]; apply_medium(model, {"assumed_serum": sc["serum_uptake"]})
        model.objective = config["objective"]; gmax = model.slim_optimize()
        task = model.reactions.get_by_id(config["objective"]); g0 = float(gmax) * sc["task_fraction"]; task.bounds = (g0, g0)
        d = np.load(root / "data/scores_conservative.npz", allow_pickle=False)
        ctxs = d["contexts"].astype(str).tolist()
        w = score_weights(model, d["rxn"].astype(str).tolist(), d["A"][:, ctxs.index(job["context"])], job["k"], job["arm"], config["epsilon"])
        res = select_flux(model, w, job["panel"], config)
        rec.update({"g_max": float(gmax), "point": res["point"], "ranges": res["ranges"], "primary_optimum": res["primary_optimum"], "status": "optimal"})
    except Exception as exc:
        rec.update({"status": "failed", "error": repr(exc)[:400], "traceback": traceback.format_exc()[-1500:]})
    rec["seconds"] = time.monotonic() - t0
    return rec


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", type=Path, required=True); ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--profiles", type=int, default=8); ap.add_argument("--workers", type=int, default=8); ap.add_argument("--solver", default=None)
    ap.add_argument("--panel", choices=["chemistry", "primary"], default="chemistry",
                    help="coordinate list handed to the selection routine: chemistry = the production list of 96 exchanges; primary = the 52-target variant of the first run")
    a = ap.parse_args(); a.out.mkdir(parents=True, exist_ok=True)
    solver = a.solver or pick_solver(); print("second solver:", solver, flush=True)
    plan = json.loads((a.root / "protocol/analysis_plan.json").read_text())
    primary = sorted(t["exchange_id"] for t in plan["primary_targets"])
    panel = sorted(t["exchange_id"] for t in plan["chemistry_targets"]) if a.panel == "chemistry" else primary
    print("selection list: %s (%d coordinates); comparison over the %d coordinates of that list" % (a.panel, len(panel), len(panel)), flush=True)
    ctxs = sorted(r["context_id"] for r in csv.DictReader(open(a.root / "manifests/contexts.tsv"), delimiter="\t") if r["partition"] == "test")
    rng = np.random.default_rng(20260916); sub = sorted(rng.choice(ctxs, a.profiles, replace=False).tolist())
    k = float(json.loads((a.root / "results/development_v3/run_summary_scalars.json").read_text())["k"])
    jobs = [{"root": str(a.root), "context": c, "arm": arm, "k": k, "panel": panel, "solver": solver} for c in sub for arm in ("magnitude_g1.00", "ordinal")]
    t0 = time.time(); res = []
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(worker, j) for j in jobs]):
            r = f.result(); res.append(r); print("%2d/%d %-14s %-16s %s %.0fs" % (len(res), len(jobs), r["context"], r["arm"], r["status"], time.time() - t0), flush=True)
            if r["status"] != "optimal": print("   ", r["error"][:300], flush=True)
    # compare against the saved GLPK arms
    rows = []; agree_fixed = agree_sign = total_fixed = total_sign = 0; pt_diffs = []
    saved = {}
    for r in res:
        if r["status"] != "optimal": continue
        tok = hashlib.sha256(("%s|%s|primary" % (r["context"], r["arm"])).encode()).hexdigest()[:16]
        s = json.loads((a.root / "results/reserved_v3/arms" / ("%s.json" % tok)).read_text()); saved[(r["context"], r["arm"])] = s
        for rid in panel:
            p2, p1 = r["point"][rid], s["point"][rid]; w2 = r["ranges"][rid][1] - r["ranges"][rid][0]; w1 = s["ranges"][rid][1] - s["ranges"][rid][0]
            f1, f2 = w1 <= EPS, w2 <= EPS; agree_fixed += (f1 == f2); total_fixed += 1; pt_diffs.append(abs(p2 - p1))
            rows.append([r["context"], r["arm"], rid, p1, p2, abs(p2 - p1), w1, w2, f1, f2, abs(r["primary_optimum"] - s["primary_optimum"]) / abs(s["primary_optimum"])])
    # between-profile ordering signs within the subset, per arm and target
    for arm in ("magnitude_g1.00", "ordinal"):
        rr = [r for r in res if r["status"] == "optimal" and r["arm"] == arm]
        for rid in panel:
            for i in range(len(rr)):
                for j in range(i + 1, len(rr)):
                    d2 = rr[i]["point"][rid] - rr[j]["point"][rid]
                    d1 = saved[(rr[i]["context"], arm)]["point"][rid] - saved[(rr[j]["context"], arm)]["point"][rid]
                    s1 = 0 if abs(d1) <= EPS else np.sign(d1); s2 = 0 if abs(d2) <= EPS else np.sign(d2)
                    agree_sign += (s1 == s2); total_sign += 1
    with open(a.out / "solver_check_rows.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(["context", "arm", "exchange_id", "point_glpk", "point_second", "abs_diff", "width_glpk", "width_second", "fixed_glpk", "fixed_second", "rel_diff_primary_optimum"]); w.writerows(rows)
    # the same summary restricted to the 52 primary targets, for comparison with the first run
    prim_rows = [row for row in rows if row[2] in set(primary)]
    prim = {"coordinates": len(prim_rows), "fixed_classification_agreement": (sum(row[8] == row[9] for row in prim_rows) / len(prim_rows)) if prim_rows else None,
            "point_abs_diff_max": max((row[5] for row in prim_rows), default=None), "points_within_1e-5": (sum(row[5] <= 1e-5 for row in prim_rows) / len(prim_rows)) if prim_rows else None}
    summ = {"second_solver": solver, "selection_list": a.panel, "selection_list_size": len(panel), "profiles": sub, "arms_run": len(jobs), "arms_optimal": sum(r["status"] == "optimal" for r in res),
            "coordinates": total_fixed, "fixed_classification_agreement": agree_fixed / max(total_fixed, 1),
            "primary_52_subset": prim,
            "point_abs_diff_median": float(np.median(pt_diffs)) if pt_diffs else None, "point_abs_diff_max": float(np.max(pt_diffs)) if pt_diffs else None,
            "points_within_1e-5": float(np.mean([d <= 1e-5 for d in pt_diffs])) if pt_diffs else None,
            "between_profile_sign_agreement": agree_sign / max(total_sign, 1), "sign_comparisons": total_sign,
            "seconds": round(time.time() - t0, 1), "failures": [{"context": r["context"], "arm": r["arm"], "error": r.get("error")} for r in res if r["status"] != "optimal"]}
    json.dump(summ, open(a.out / "solver_check.json", "w"), indent=1); print(json.dumps({k: v for k, v in summ.items() if k != "profiles"}, indent=1))


if __name__ == "__main__":
    main()
