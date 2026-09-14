#!/usr/bin/env python3
"""Independent evaluation: does the reporting rule call the measured condition responses?

This is the first script in the Copeland analysis that reads a measured glucose or lactate
value. The plan, the medium, the preprocessing interface and the predictions were all fixed
before it ran.

For each declared contrast and target:
  * the measurement either resolves a direction, by complete separation of the four replicates
    of one condition from the four of the other, or it does not;
  * the reporting rule either returns a direction, by requiring every encoding in the
    uncertainty set and every cross-replicate profile pair to clear the threshold, or it
    abstains;
  * agreement is judged only where both resolve, and coverage is reported beside it, because a
    high conditional agreement at negligible coverage is not a useful result.

The context-independent uniform baseline is evaluated by the identical procedure. If it calls
the same responses, the response follows from the environment and task rather than from the
transcripts.
"""
import argparse, csv, hashlib, json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from copeland_rule import EPS, measured_direction, interval_direction, point_category, agreement

MEAS_NAME = {"EX_glc__D_e": "glucose", "EX_lac__L_e": "lactate"}


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def load_measured(compendium):
    import rdata
    df = rdata.conversion.convert(rdata.parser.parse_file(
        Path(compendium) / "data" / "fluxes.rda"))["fluxes"]
    df.columns = [str(c) for c in df.columns]
    keep = df[(df["cell_type"].astype(str) == "lf") &
              (df["experiment"].astype(str) == "05-bay")]
    out = defaultdict(list)
    for _, r in keep.iterrows():
        out[(str(r["metabolite"]), str(r["oxygen"]), str(r["treatment"]))].append(float(r["flux"]))
    return {k: sorted(v) for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", type=Path, required=True)
    ap.add_argument("--compendium", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    code = Path(__file__).resolve().parent
    plan = json.loads((code / "copeland_plan.json").read_text())
    sets = plan["uncertainty_set"]

    preds = list(csv.DictReader(open(a.run / "predictions.tsv"), delimiter="\t"))
    P = defaultdict(list)     # (task_mode, armset_member, oxygen, treatment, target) -> intervals
    for r in preds:
        P[(r["task_mode"], r["arm"], r["oxygen"], r["treatment"], r["exchange_id"])].append(
            (float(r["point"]), float(r["range_lo"]), float(r["range_hi"])))

    measured = load_measured(a.compendium)
    contrasts = (plan["contrasts"]["primary_drug_within_oxygen"] +
                 plan["contrasts"]["secondary_oxygen_within_drug"])
    families = {"primary_five_power": sets["primary"],
                "sensitivity_six_with_ordinal": sets["sensitivity"],
                "context_independent_baseline": ["uniform_pfba"]}

    rows, summary = [], {}
    for task_mode in sorted({r["task_mode"] for r in preds}):
        for fam, arms in families.items():
            counts = defaultdict(int)
            for c in contrasts:
                for target, mname in MEAS_NAME.items():
                    ma = measured.get((mname, c["a"]["oxygen"], c["a"]["treatment"]), [])
                    mb = measured.get((mname, c["b"]["oxygen"], c["b"]["treatment"]), [])
                    md = measured_direction(ma, mb)
                    ia = [(lo, hi) for arm in arms
                          for _, lo, hi in P[(task_mode, arm, c["a"]["oxygen"], c["a"]["treatment"], target)]]
                    ib = [(lo, hi) for arm in arms
                          for _, lo, hi in P[(task_mode, arm, c["b"]["oxygen"], c["b"]["treatment"], target)]]
                    pa = [p for arm in arms
                          for p, _, _ in P[(task_mode, arm, c["a"]["oxygen"], c["a"]["treatment"], target)]]
                    pb = [p for arm in arms
                          for p, _, _ in P[(task_mode, arm, c["b"]["oxygen"], c["b"]["treatment"], target)]]
                    idir = interval_direction(ia, ib)
                    cat, pdir = point_category(pa, pb)
                    ag_i = agreement(idir, md)
                    ag_p = agreement(pdir, md)
                    counts["contrasts"] += 1
                    counts["measurement_eligible"] += md != 0
                    counts["interval_reported"] += idir != 0
                    counts["point_reported"] += pdir != 0
                    if md != 0:
                        counts["interval_reported_on_eligible"] += idir != 0
                        counts["point_reported_on_eligible"] += pdir != 0
                    if ag_i is True: counts["interval_correct"] += 1
                    if ag_i is False: counts["interval_incorrect"] += 1
                    if ag_p is True: counts["point_correct"] += 1
                    if ag_p is False: counts["point_incorrect"] += 1
                    rows.append([task_mode, fam, c["id"], mname, target,
                                 c["a"]["oxygen"], c["a"]["treatment"],
                                 c["b"]["oxygen"], c["b"]["treatment"],
                                 md, round(float(np.mean(ma)), 3) if ma else "",
                                 round(float(np.mean(mb)), 3) if mb else "",
                                 idir, pdir, cat,
                                 "" if ag_i is None else int(ag_i),
                                 "" if ag_p is None else int(ag_p),
                                 len(ia), len(ib),
                                 round(min(lo for lo, _ in ia), 6) if ia else "",
                                 round(max(hi for _, hi in ia), 6) if ia else "",
                                 round(min(lo for lo, _ in ib), 6) if ib else "",
                                 round(max(hi for _, hi in ib), 6) if ib else ""])
            summary["%s|%s" % (task_mode, fam)] = dict(counts)

    with open(a.out / "copeland_decisions.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["task_mode", "encoding_family", "contrast", "metabolite", "exchange_id",
                    "a_oxygen", "a_treatment", "b_oxygen", "b_treatment",
                    "measured_direction", "a_measured_mean", "b_measured_mean",
                    "interval_direction", "point_direction", "point_category",
                    "interval_agrees", "point_agrees", "a_intervals", "b_intervals",
                    "a_min_lo", "a_max_hi", "b_min_lo", "b_max_hi"])
        w.writerows(rows)

    report = {"epsilon": EPS, "contrasts_declared": len(contrasts) * len(MEAS_NAME),
              "summary": summary,
              "measured_values": {"|".join(k): v for k, v in sorted(measured.items())},
              "inputs_sha256": {
                  "copeland_plan.json": digest(code / "copeland_plan.json"),
                  "copeland_rule.py": digest(code / "copeland_rule.py"),
                  "predictions.tsv": digest(a.run / "predictions.tsv"),
                  "run_manifest.json": digest(a.run / "run_manifest.json")},
              "interpretation_limits": plan["what_this_cannot_establish"]}
    json.dump(report, open(a.out / "copeland_evaluation.json", "w"), indent=1)

    print("%-10s %-32s %6s %6s %8s %8s %6s %6s" % ("task", "encoding family", "contr", "elig",
          "int_rep", "pt_rep", "i_ok", "i_bad"))
    for k, v in summary.items():
        t, f = k.split("|")
        print("%-10s %-32s %6d %6d %8d %8d %6d %6d" % (t, f, v["contrasts"],
              v["measurement_eligible"], v["interval_reported"], v["point_reported"],
              v.get("interval_correct", 0), v.get("interval_incorrect", 0)))
    print("\nEligible contrasts and what the rule did:")
    for r in rows:
        if r[9] != 0:
            print("  %-9s %-32s %-9s %-8s measured %+d | interval %+d | point %+d (%s)"
                  % (r[0], r[1], r[2], r[3], r[9], r[12], r[13], r[14]))


if __name__ == "__main__":
    main()
