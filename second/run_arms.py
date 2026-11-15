#!/usr/bin/env python3
"""Second reconstruction: run the arms of the main analysis in the second root.

    python3 second/run_arms.py --root SECOND_ROOT --partition development --output results/development_v3 --workers 12
    python3 second/run_arms.py --root SECOND_ROOT --partition test        --output results/reserved_v3    --workers 12
    python3 second/run_arms.py --root SECOND_ROOT --write-lock PLAN.md    # once, when the plan is fixed

Pass 1 calls src.run_study.main() of the root, unchanged: the same job list (every profile of the
partition under the five power encodings and the midrank encoding, and the uniform cost once), the
same worker, the solver of the main analysis (GLPK) at its tolerances.

Pass 2. On Human-GEM the lexicographic point selection of the readout returned an infeasible problem
under GLPK in 2 of the 67 development arms; the admissible ranges were computed in every arm. An arm
whose readout failed in pass 1 is therefore solved again by the same worker with the model built on
the second solver of the main analysis's cross-solver check (optlang's hybrid interface, which uses
HiGHS; closure/solver_check.py), at the same tolerances. Nothing else differs between the passes. Arms
that pass 1 solved are not touched (the worker reuses a saved arm whose inputs and flux file verify).
The arms solved in pass 2 are listed in <output>/second_solver_arms.json.

For the reserved partition src.run_study verifies a protocol lock. The lock of the main analysis pins
the inputs of the main root. The second root carries its own lock, written by --write-lock when the
plan of the second reconstruction is fixed: it pins the SHA-256 of that plan and of the inputs of this
root, and the reserved partition does not run without it.

No measured exchange is read here.
"""
import argparse, hashlib, json, sys, time
from pathlib import Path

LOCK_INPUTS = ("data/raw/Recon3D.json", "configs/study.json", "environment.lock", "data/scores_conservative.npz",
               "manifests/contexts.tsv", "manifests/metabolite_map.tsv", "protocol/analysis_plan.json",
               "src/cpu_model.py", "src/cpu_readout.py", "src/medium_v2.py", "src/run_study.py", "src/gpr.py", "src/evaluation.py")
PLAN_IN_ROOT = "protocol/SECOND_RECONSTRUCTION_PLAN.md"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_lock(root, plan):
    dst = root / PLAN_IN_ROOT; dst.write_bytes(Path(plan).read_bytes())
    lock = dict(status="locked", reserved_evaluation_authorized=True,
                scope="second reconstruction (Human-GEM); this lock replaces, in this root only, the lock of the main analysis",
                locked_at_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), plan=PLAN_IN_ROOT, plan_source=Path(plan).name,
                second_reconstruction_plan_sha256=sha(dst), input_sha256={p: sha(root / p) for p in LOCK_INPUTS})
    (root / "protocol/PROTOCOL_LOCK.json").write_text(json.dumps(lock, indent=1) + "\n")
    print(json.dumps({k: lock[k] for k in ("status", "locked_at_utc", "second_reconstruction_plan_sha256")}, indent=1))


def verify_second_lock(root, lock_path=None):
    root = Path(root); lock = json.loads((root / "protocol/PROTOCOL_LOCK.json").read_text())
    if lock.get("status") != "locked" or not lock.get("reserved_evaluation_authorized") or "second_reconstruction_plan_sha256" not in lock:
        raise RuntimeError("the reserved partition needs the lock of the second reconstruction (run_arms.py --write-lock)")
    if sha(root / PLAN_IN_ROOT) != lock["second_reconstruction_plan_sha256"]:
        raise RuntimeError("the plan in the root differs from the locked plan")
    for p, want in lock["input_sha256"].items():
        if sha(root / p) != want: raise RuntimeError(f"locked input differs: {p}")
    return lock


def build_model_second_solver(root, config, scenario):
    """src.cpu_model.build_model with the model on the second solver; the same steps in the same order."""
    import numpy as np, cobra
    from src.medium_v2 import apply_medium
    model = cobra.io.load_json_model(str(Path(root) / "data/raw/Recon3D.json"))
    model.solver = "hybrid"
    model.solver.configuration.tolerances.feasibility = config["solver_feasibility_tolerance"]
    model.solver.configuration.tolerances.optimality = config["solver_optimality_tolerance"]
    model.solver.configuration.timeout = 600
    sc = config["scenarios"][scenario]
    opened, missing, changes = apply_medium(model, {"assumed_serum": sc["serum_uptake"]})
    model.objective = config["objective"]
    gmax = model.slim_optimize()
    if model.solver.status != "optimal" or not np.isfinite(gmax) or gmax <= 0:
        raise RuntimeError(f"Invalid task maximum: {gmax}, {model.solver.status}")
    task = model.reactions.get_by_id(config["objective"])
    g0 = float(gmax) * sc["task_fraction"]
    task.bounds = (g0, g0)
    return model, {"g_max": float(gmax), "task_flux": g0, "objective": task.id, "objective_name": task.name, "scenario": scenario,
                   "opened": opened, "missing": missing, "bound_changes": changes, "solver": "hybrid (HiGHS), pass 2"}


def run(RS, argv):
    old = sys.argv; sys.argv = ["src.run_study"] + argv
    try:
        RS.main(); return 0
    except SystemExit as e:
        return int(e.code or 0)
    finally:
        sys.argv = old


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True); ap.add_argument("--partition", choices=["development", "test"])
    ap.add_argument("--output"); ap.add_argument("--workers", type=int, default=4); ap.add_argument("--write-lock", type=Path)
    a = ap.parse_args(); root = a.root.resolve()
    if a.write_lock: return write_lock(root, a.write_lock)
    sys.path.insert(0, str(root))
    import src.verify_inputs as V
    V.verify_protocol_lock = verify_second_lock                 # the reserved partition runs under the lock of the second reconstruction
    import src.run_study as RS
    assert Path(RS.ROOT).resolve() == root, "run from the second root's own src/"
    out = (root / a.output).resolve(); argv = ["--partition", a.partition, "--output", str(out), "--workers", str(a.workers)]
    if a.partition == "test": print("lock:", verify_second_lock(root)["second_reconstruction_plan_sha256"], flush=True)
    code1 = run(RS, argv)
    arms = sorted((out / "arms").glob("*.json"))
    failed = [json.loads(p.read_text()) for p in arms]; failed = [r for r in failed if r["status"] != "optimal"]
    record = dict(pass1_solver="glpk", pass1_exit=code1, arms=len(arms),
                  failed_in_pass1=[dict(context_id=r["context_id"], arm=r["arm"], error=(r.get("error") or "")[:200]) for r in failed])
    if failed:
        print(f"pass 2: {len(failed)} arm(s) on the second solver", flush=True)
        RS.build_model = build_model_second_solver
        code2 = run(RS, argv)
        again = [json.loads(p.read_text()) for p in arms]
        record.update(pass2_solver="hybrid (HiGHS)", pass2_exit=code2,
                      solved_in_pass2=[dict(context_id=r["context_id"], arm=r["arm"], seconds=r["seconds"]) for r in again
                                       if r["status"] == "optimal" and (r.get("model") or {}).get("solver", "").startswith("hybrid")],
                      still_failed=[dict(context_id=r["context_id"], arm=r["arm"], error=(r.get("error") or "")[:200]) for r in again if r["status"] != "optimal"])
    (out / "second_solver_arms.json").write_text(json.dumps(record, indent=1) + "\n")
    s = json.loads((out / "run_summary.json").read_text())
    scalars = {k: s[k] for k in ("k", "input_fingerprint", "run_id", "partition", "contexts", "completed", "failed", "seconds")}
    scalars["note"] = "Scalar extract of run_summary.json. The scale constant k is taken from development and is not recomputed on reserved data."
    (out / "run_summary_scalars.json").write_text(json.dumps(scalars, indent=1) + "\n")
    print(json.dumps({k: record.get(k) for k in ("arms", "pass1_exit", "pass2_exit")} | dict(failed_in_pass1=len(record["failed_in_pass1"]),
                     still_failed=len(record.get("still_failed", []))), indent=1))
    return 1 if record.get("still_failed") or (not failed and code1) else 0


if __name__ == "__main__":
    raise SystemExit(main())
