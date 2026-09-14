#!/usr/bin/env python3
"""Control: does a uniform cost cap fix the reported exchanges as tightly as the transcript one?

The constraint-stage analysis attributed most coordinate fixing to the transcript-weighted cost
cap. That attribution is only meaningful if a cost cap carrying no transcript information does
less. The transcript-weighted objective prices internal reactions only, leaving every boundary
reaction at zero cost, so the reported exchanges are never priced directly and the attribution
could instead belong to the L1 cap itself at the tolerance used.

This runs the identical B0/B1/B2 ladder with the uniform-cost arm, whose weights carry no
profile information at all. The uniform arm is context independent under this formulation, so
one B1 calculation characterises it. Its fixed-coordinate count is then comparable with the
per-profile counts of the six transcript-weighted arms over the same 52 targets.

The per-job worker is imported unchanged from the pilot.
"""
import sys
sys.dont_write_bytecode = True
import argparse, csv, hashlib, json, time
from pathlib import Path
from datetime import datetime, timezone

sys.path.insert(0, str(Path(__file__).resolve().parent))
UNIFORM_JSON = None
from nested_range_pilot import worker, nested_violation, stage_label, sha, read_tsv, dump, EPS, TOL


def main(root, out, budget):
    global UNIFORM_JSON
    out.mkdir(parents=True, exist_ok=True)
    plan = json.loads((root / "protocol/analysis_plan.json").read_text())
    config = json.loads((root / "configs/study.json").read_text())
    panel = sorted(t["exchange_id"] for t in plan["primary_targets"])
    assert len(panel) == 52

    # the saved uniform arm for the primary scenario, named explicitly
    saved_json = Path(UNIFORM_JSON)
    saved = json.loads(saved_json.read_text())
    if saved["arm"] != "uniform_pfba" or saved["scenario"] != "primary":
        raise SystemExit("named record is not the primary uniform arm")
    saved_npz = root / "results/reserved_v3" / saved["flux_file"]
    if sha(saved_npz) != saved["flux_sha256"]:
        raise SystemExit("saved uniform flux hash mismatch")

    # any reserved context indexes the score matrix; the uniform arm ignores the evidence, and
    # the worker's own 1e-12 cost comparison against the saved vector proves that it did
    ctx = sorted(r["context_id"] for r in read_tsv(root / "manifests/contexts.tsv")
                 if r["partition"] == "test")[0]
    k = float(json.loads((root / "results/development_v3/run_summary_scalars.json").read_text())["k"])

    start = time.monotonic()
    common = {"root": str(root), "out": str(out), "panel": panel, "deadline": start + budget}
    b0 = worker({**common, "token": "B0", "stage": "B0"})
    print("B0", b0["status"], round(b0["seconds"], 2), flush=True)
    if b0["status"] != "complete":
        raise SystemExit("B0 failed: %s" % b0.get("error"))
    r = worker({**common, "token": "B1_uniform", "stage": "B1", "context": ctx,
                "arm": "uniform_pfba", "saved_json": str(saved_json),
                "saved_npz": str(saved_npz), "k": k})
    print("B1_uniform", r["status"], round(r["seconds"], 2), flush=True)
    if r["status"] != "complete":
        raise SystemExit("B1 uniform failed: %s" % r.get("error"))

    rows, issues = [], []
    for rid in panel:
        z = b0["ranges"][rid]; one = r["ranges"][rid]
        two = r["saved_B2_ranges"][rid]; point = r["saved_B2_points"][rid]
        v01 = nested_violation(z, one); v12 = nested_violation(one, two)
        pv = max(0., z[0] - point, point - z[1], one[0] - point, point - one[1],
                 two[0] - point, point - two[1])
        if max(v01, v12, pv) > TOL:
            issues.append({"exchange_id": rid, "B0_B1": v01, "B1_B2": v12, "point": pv})
        w0, w1, w2 = [x[1] - x[0] for x in (z, one, two)]
        rows.append({"arm": "uniform_pfba", "exchange_id": rid,
                     "B0_low": z[0], "B0_high": z[1], "B1_low": one[0], "B1_high": one[1],
                     "B2_low": two[0], "B2_high": two[1],
                     "B0_width": w0, "B1_width": w1, "B2_width": w2,
                     "B0_fixed": w0 <= EPS, "B1_fixed": w1 <= EPS, "B2_fixed": w2 <= EPS,
                     "first_fixed_stage": stage_label(z, one, two),
                     "B0_B1_nesting_violation": v01, "B1_B2_nesting_violation": v12,
                     "saved_point_containment_violation": pv})
    with (out / "range_stage_ledger_uniform.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)

    counts = {lab: sum(x["first_fixed_stage"] == lab for x in rows)
              for lab in ("already_fixed_B0", "first_fixed_B1", "first_fixed_B2",
                          "still_variable_B2")}
    free = len(rows) - counts["already_fixed_B0"]
    res = {"status": "complete_verified" if not issues else "validation_failure",
           "arm": "uniform_pfba", "scenario": "primary", "targets": len(rows),
           "context_independent": True, "counts": counts,
           "free_targets_after_B0": free,
           "share_of_free_first_fixed_by_B1": counts["first_fixed_B1"] / free if free else None,
           "validation_issues": issues, "seconds": time.monotonic() - start,
           "completed_at_utc": datetime.now(timezone.utc).isoformat(),
           "interpretation": ("A uniform cost cap carries no transcript information. Its "
                              "fixed-coordinate share is the reference against which the "
                              "transcript-weighted arms' share must be read.")}
    dump(out / "summary_uniform.json", res)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--budget", type=int, default=1200)
    p.add_argument("--uniform-json", type=Path, required=True)
    a = p.parse_args()
    UNIFORM_JSON = a.uniform_json
    main(a.root.resolve(), a.output.resolve(), a.budget)
