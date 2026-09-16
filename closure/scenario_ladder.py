#!/usr/bin/env python3
"""Reserved-cohort confirmation of the constraint-stage diagnosis.

The bounded development pilot found that the primary transcript-weighted cost cap, not the
network/medium/task constraints, is what fixes most exchange coordinates, and the location
analysis then found that the fixed coordinates coincide across profiles. The pilot's own stop
rule allows expansion once a concrete explanation needs confirming. This runs the same
computation on all 47 reserved profiles and all six encodings in the primary scenario.

The per-job worker is imported unchanged from the pilot so the arithmetic, the comparability
checks and the numerical audits are literally the same code. Only the job list differs.

Nothing here reads a CORE outcome. Saved B2 results are reused, never reoptimized.
"""
import sys
sys.dont_write_bytecode = True
import argparse, csv, hashlib, json, time
from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from nested_range_pilot_scen import (worker, nested_violation, stage_label, sha, read_tsv, dump,
                                EPS, TOL)

ARMS = ["magnitude_g0.50", "magnitude_g0.80", "magnitude_g1.00", "magnitude_g1.25",
        "magnitude_g2.00", "ordinal"]
SCENARIO = None
PARTITION = "test"


def main(root, out, workers, budget, scenario):
    global SCENARIO
    SCENARIO = scenario
    out.mkdir(parents=True, exist_ok=True)
    if (out / "run_manifest.json").exists():
        raise FileExistsError("Do not overwrite an existing run")
    config = json.loads((root / "configs/study.json").read_text())
    plan = json.loads((root / "protocol/analysis_plan.json").read_text())
    lock = json.loads((root / "protocol/PROTOCOL_LOCK.json").read_text())
    if lock.get("status") != "locked" or not lock.get("reserved_evaluation_authorized"):
        raise ValueError("Reserved evaluation is not authorized by the protocol lock")
    if plan["point_tolerance"] != EPS or plan["numerical_tolerance"] != TOL:
        raise ValueError("Tolerance mismatch against the frozen plan")

    panel = sorted(t["exchange_id"] for t in plan["primary_targets"])
    if len(panel) != 52 or len(set(panel)) != 52:
        raise ValueError("Unexpected panel")
    contexts = sorted(r["context_id"] for r in read_tsv(root / "manifests/contexts.tsv")
                      if r["partition"] == PARTITION)
    if len(contexts) != 47:
        raise ValueError("Unexpected reserved allocation: %d" % len(contexts))

    # the scale constant is frozen from development and must not be recomputed on reserved data
    k = float(json.loads((root / "results/development_v3/run_summary_scalars.json").read_text())["k"])

    sources = ["data/raw/Recon3D.json", "configs/study.json", "environment.lock",
               "data/scores_conservative.npz", "manifests/contexts.tsv",
               "manifests/metabolite_map.tsv", "protocol/analysis_plan.json",
               "protocol/PROTOCOL_LOCK.json", "src/cpu_model.py", "src/cpu_readout.py",
               "src/medium_v2.py", "src/run_study.py"]
    hashes = {p: sha(root / p) for p in sources}
    for name, want in lock["input_sha256"].items():
        if name in hashes and hashes[name] != want:
            raise ValueError("Locked input differs: %s" % name)

    selections = []
    for context in contexts:
        for arm in ARMS:
            token = hashlib.sha256(("%s|%s|%s" % (context, arm, SCENARIO)).encode()).hexdigest()[:16]
            saved_json = root / "results/reserved_v3/arms" / ("%s.json" % token)
            saved = json.loads(saved_json.read_text())
            saved_npz = root / "results/reserved_v3" / saved["flux_file"]
            if (saved["status"] != "optimal" or saved["context_id"] != context
                    or saved["arm"] != arm or saved["scenario"] != SCENARIO
                    or sha(saved_npz) != saved["flux_sha256"]):
                raise ValueError("Saved reserved record identity or hash failure: %s" % token)
            for stage in ("primary", "secondary"):
                prefix = "cost" if stage == "primary" else "secondary"
                expected = (saved["%s_optimum" % stage] + config["%s_absolute_allowance" % prefix]
                            + abs(saved["%s_optimum" % stage]) * config["%s_relative_allowance" % prefix])
                if abs(expected - saved["%s_cap" % stage]) > 1e-12:
                    raise ValueError("Saved cap formula mismatch: %s" % token)
            selections.append({"token": token, "context": context, "arm": arm, "stage": "B1",
                               "saved_json": str(saved_json), "saved_npz": str(saved_npz), "k": k, "scenario": SCENARIO})
    expected_jobs = len(contexts) * len(ARMS)
    if len(selections) != expected_jobs:
        raise ValueError("Job count mismatch")

    hashes["__driver__"] = sha(Path(__file__).resolve())
    hashes["__worker__"] = sha(Path(__file__).resolve().parent / "nested_range_pilot_scen.py")
    start = time.monotonic()
    deadline = start + budget
    common = {"root": str(root), "out": str(out), "panel": panel, "deadline": deadline}
    dump(out / "run_manifest.json", {
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "preflight_passed_before_LPs", "scenario": SCENARIO, "partition": PARTITION,
        "inputs_sha256": hashes, "expected_B0_jobs": 1, "expected_B1_jobs": expected_jobs,
        "workers": workers, "runtime_budget_seconds": budget, "arms": ARMS,
        "profiles": contexts, "targets": panel, "scale_constant_k": k,
        "outcomes_read": False, "original_files_mutated": False,
        "saved_B2_extrema_reoptimized": False,
        "python": sys.version, "python_executable": sys.executable})

    b0 = worker({**common, "token": "B0", "stage": "B0", "scenario": SCENARIO})
    print("B0", b0["status"], round(b0["seconds"], 2), flush=True)
    results = []
    if b0["status"] == "complete":
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(worker, {**common, **s}) for s in selections]
            for fut in as_completed(futures):
                r = fut.result(); results.append(r)
                print(len(results), r["context"], r["arm"], r["status"], round(r["seconds"], 2),
                      flush=True)

    rows, issues = [], []
    for r in results:
        if r["status"] != "complete":
            continue
        for rid in panel:
            z = b0["ranges"][rid]; one = r["ranges"][rid]
            two = r["saved_B2_ranges"][rid]; point = r["saved_B2_points"][rid]
            v01 = nested_violation(z, one); v12 = nested_violation(one, two)
            pv = max(0., z[0] - point, point - z[1], one[0] - point, point - one[1],
                     two[0] - point, point - two[1])
            if max(v01, v12, pv) > TOL:
                issues.append({"context": r["context"], "arm": r["arm"], "exchange_id": rid,
                               "B0_B1": v01, "B1_B2": v12, "point": pv})
            w0, w1, w2 = [x[1] - x[0] for x in (z, one, two)]
            rows.append({"context": r["context"], "arm": r["arm"], "exchange_id": rid,
                         "B0_low": z[0], "B0_high": z[1], "B1_low": one[0], "B1_high": one[1],
                         "B2_low": two[0], "B2_high": two[1],
                         "B0_width": w0, "B1_width": w1, "B2_width": w2,
                         "B0_fixed": w0 <= EPS, "B1_fixed": w1 <= EPS, "B2_fixed": w2 <= EPS,
                         "B0_to_B1_width_contraction": w0 - w1,
                         "B1_to_B2_width_contraction": w1 - w2,
                         "first_fixed_stage": stage_label(z, one, two),
                         "B0_B1_nesting_violation": v01, "B1_B2_nesting_violation": v12,
                         "saved_point_containment_violation": pv})
    if rows:
        with (out / "range_stage_ledger.tsv").open("w", newline="") as f:
            w = csv.DictWriter(f, list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)

    grouped = []
    for context in contexts:
        for arm in ARMS:
            sub = [r for r in rows if r["context"] == context and r["arm"] == arm]
            if not sub:
                continue
            e = {"context": context, "arm": arm, "targets": len(sub)}
            for st in ("B0", "B1", "B2"):
                e[st + "_fixed"] = sum(r[st + "_fixed"] for r in sub)
            for lab in ("already_fixed_B0", "first_fixed_B1", "first_fixed_B2", "still_variable_B2"):
                e[lab] = sum(r["first_fixed_stage"] == lab for r in sub)
            grouped.append(e)

    completed = sum(r["status"] == "complete" for r in results)
    unchanged = all(sha(root / p) == v for p, v in hashes.items() if not p.startswith("__"))
    res = {"status": ("complete_verified" if completed == expected_jobs and not issues and unchanged
                      else "partial_or_validation_failure"),
           "completed_at_utc": datetime.now(timezone.utc).isoformat(),
           "seconds": time.monotonic() - start, "scenario": SCENARIO, "partition": PARTITION,
           "B0_status": b0["status"], "B1_completed": completed,
           "B1_failed_or_incomplete": sum(r["status"] != "complete" for r in results),
           "B1_unstarted": expected_jobs - len(results), "reused_B2_records": completed,
           "target_context_arm_rows": len(rows), "validation_issues": issues,
           "all_input_hashes_unchanged": unchanged, "summary_by_profile_arm": grouped,
           "outcomes_read": False, "point_tolerance": EPS, "nesting_tolerance": TOL,
           "interpretation": ("Conditional within-profile coordinate restriction on the reserved "
                              "cohort. Location analysis is reported separately.")}
    dump(out / "summary.json", res)
    print(json.dumps({k2: v for k2, v in res.items() if k2 != "summary_by_profile_arm"}, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--workers", type=int, default=20)
    p.add_argument("--budget", type=int, default=5400)
    p.add_argument("--scenario", required=True, choices=["half_serum", "lower_task"])
    a = p.parse_args()
    main(a.root.resolve(), a.output.resolve(), a.workers, a.budget, a.scenario)
