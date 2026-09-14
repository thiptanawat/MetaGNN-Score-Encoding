#!/usr/bin/env python3
"""Check every numeric claim in the manuscript against the saved analysis outputs.

Each check states the claim, recomputes or re-reads the value from the artefact that produced
it, and compares. A claim that cannot be traced to an artefact is reported as UNCOVERED rather
than passed, so the gap is visible instead of implied.
"""
import csv, json, re, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
TEXT = (HERE / "Paper3_PeerJ_Manuscript.md").read_text()
SUPP = (HERE / "Supplemental_Article_S1.md").read_text()

results = []


def check(label, claim_in_text, expected, actual, tol=0.0):
    """Confirm the string appears in the manuscript and that the artefact value matches."""
    present = claim_in_text in TEXT or claim_in_text in SUPP
    if isinstance(expected, str):
        ok_value = str(actual) == expected
    else:
        ok_value = abs(float(actual) - float(expected)) <= tol
    results.append((label, present, ok_value, claim_in_text, expected, actual))


def load_json(rel):
    return json.load(open(BASE / rel))


def load_tsv(rel):
    return list(csv.DictReader(open(BASE / rel), delimiter="\t"))


# ---------------------------------------------------------------- encoding sensitivity
t = load_json("cpu_study_2026-09-13/results/reserved_v3/independent_concordance_primary.json")
check("midrank concordance, primary", "midrank concordance was 0.5054",
      0.5054, round(t["absolute_C"]["ordinal"], 4))
check("grid contrast", "-0.0013", -0.0013, round(t["delta_grid"], 4))
check("grid contrast interval low", "-0.0087", -0.0087, round(t["delta_grid_interval"][0], 4))
check("grid contrast interval high", "0.0069", 0.0069, round(t["delta_grid_interval"][1], 4))
check("eligible comparisons", "46,447", 46447, t["observed_pairs"])
check("origin groups", "44 origins", 44, t["origins"])

# ---------------------------------------------------------------- unanimity coverage
ag = load_json("strengthening_reference_checked_2026-09-14/computation/agreement_summary.json")
row = {(r["scenario"], r["arm_set"]): r for r in ag["summary"]}
p5 = row[("primary", "five_power")]
check("resolved share, five encodings", "5.73% of comparisons resolved",
      5.73, round(100 * p5["point_unanimous_nontied_coverage"], 2))
check("tied share, five encodings", "88.20% tied",
      88.20, round(100 * p5["unanimous_ties_coverage"], 2))
check("resolved count", "3,214 of these (5.73%)", 3214, p5["point_unanimous_nontied"])
check("tie count", "49,442 (88.20%)", 49442, p5["unanimous_ties"])
check("shared-encoding support", "3,204 (5.72%)", 3204, p5["matched_supported"])
check("independently varying support", "2,636 (4.70%)", 2636, p5["independent_supported"])
p6 = row[("primary", "six_including_ordinal")]
check("resolved count, six encodings", "1,443 (2.57%)", 1443, p6["point_unanimous_nontied"])

# ---------------------------------------------------------------- constraint-stage attribution
rr = load_json("mechanism_2026-09-14/reserved_range/summary.json")
check("reserved jobs completed", "282 calculations", 282, rr["B1_completed"])
check("reserved jobs failed", "every one of the 282 calculations completed", 0, rr["B1_failed_or_incomplete"])
check("reserved validation issues", "no validation issue", 0, len(rr["validation_issues"]))
check("reserved ledger rows", "14,664 ledger rows", 14664, rr["target_context_arm_rows"])
check("reserved hashes unchanged", "left every input hash unchanged", "True",
      str(rr["all_input_hashes_unchanged"]))

ledger = load_tsv("mechanism_2026-09-14/reserved_range/range_stage_ledger.tsv")
by_arm = defaultdict(lambda: defaultdict(int))
for r in ledger:
    by_arm[r["arm"]][r["first_fixed_stage"]] += 1
    by_arm[r["arm"]]["rows"] += 1
first_b1 = {a: v["first_fixed_B1"] for a, v in by_arm.items()}
share = {a: 100 * first_b1[a] / by_arm[a]["rows"] for a in by_arm}
check("rows per arm", "2,444 profile-target coordinates", 2444,
      min(v["rows"] for v in by_arm.values()))
check("first fixed at B1, minimum", "2,077", 2077, min(first_b1.values()))
check("first fixed at B1, maximum", "2,110", 2110, max(first_b1.values()))
free_share = {a: 100 * first_b1[a] / (by_arm[a]["rows"] - by_arm[a]["already_fixed_B0"])
              for a in by_arm}
check("free coordinates per encoding", "2,303 coordinates per encoding", 2303,
      min(by_arm[a]["rows"] - by_arm[a]["already_fixed_B0"] for a in by_arm))
check("share of free coordinates first fixed by the primary cost", "90 to 92 percent",
      True, 90.0 <= min(free_share.values()) and max(free_share.values()) <= 92.0)

loc = load_json("mechanism_2026-09-14/range_location.json")
res = loc["reserved_primary_stage_resolved"]
check("targets fixed by the network alone", "only 3 of 52 targets", 3,
      res["magnitude_g1.00"]["targets_fixed_by_B0_alone"])
check("free targets", "49 free targets", 49,
      res["magnitude_g1.00"]["targets_with_admissible_width_at_B0"])
consts = [v["of_those_constant_across_profiles_at_B2"] for a, v in res.items()
          if a.startswith("magnitude")]
check("constant across profiles, magnitude arms", "41 or 42 of 49", True,
      set(consts) <= {41, 42})
meds = [v["B2_exploration_median"] for a, v in res.items()]
check("median exploration order of magnitude", "4 parts in 10 billion", True,
      2.5e-10 <= min(meds) and max(meds) <= 4.5e-10)
check("mean B0 width", "35.0 canonical model units", 35.0,
      round(loc["development_primary"]["magnitude_g1.00"]["mean_B0_width_of_free_targets"], 1))

# ---------------------------------------------------------------- cross-encoding comparison
cc = load_json("mechanism_2026-09-14/cross_cost.json")
check("self-cost reconstruction", "1.5 parts in a billion", True,
      cc["self_cost_worst_relative"] < 2e-9)
check("arms rescored", "849 arms", 849, cc["arms_loaded"])
pr = cc["primary"]
check("median cross regret", "0.04 percent", 0.04, round(100 * pr["regret_median"], 2))
check("regret p95", "7.5 percent", 7.5, round(100 * pr["regret_p95"], 1))
check("regret max", "11.8 percent", 11.8, round(100 * pr["regret_max"], 1))
check("share outside the other cap", "97.4 percent", 97.4,
      round(100 * (1 - pr["share_within_the_other_arms_primary_cap"]), 1))
check("off-diagonal comparisons", "1,410 cross-encoding comparisons", 1410,
      pr["off_diagonal_comparisons"])
check("reactions differing", "median of 117", 117, pr["reactions_differing_median"])
check("full-vector distance", "median relative L1 distance of 0.054", 0.054,
      round(pr["relative_l1_full_median"], 3))
check("panel distance", "the median relative L1 distance was 0.016", 0.016,
      round(pr["relative_l1_panel_median"], 3))
check("panel exchanges differing", "median of 6 of the 96 exchanges", 6,
      pr["panel_exchanges_differing_median"])

# ---------------------------------------------------------------- independent evaluation
ce = load_json("copeland_2026-09-14/evaluation/copeland_evaluation.json")
s = ce["summary"]
check("declared contrasts", "8 declared contrasts", 8, s["measured|primary_five_power"]["contrasts"])
check("measurement-eligible", "the measurements resolved 2", 2,
      s["measured|primary_five_power"]["measurement_eligible"])
check("reported under the original task", "reported a direction for none of the 8", 0,
      s["fraction|primary_five_power"]["interval_reported"])
check("reported under measured growth", "reported a direction for all 8", 8,
      s["measured|primary_five_power"]["interval_reported"])
check("correct calls", "correct once and incorrect once", 1,
      s["measured|primary_five_power"]["interval_correct"])
check("incorrect calls", "correct once and incorrect once", 1,
      s["measured|primary_five_power"]["interval_incorrect"])
check("uniform reference coverage", "reported a direction for 4 of the 8", 4,
      s["measured|context_independent_baseline"]["interval_reported"])
mv = ce["measured_values"]
check("lactate 21 percent vehicle mean", "975.15", 975.15,
      round(float(np.mean(mv["lactate|21%|DMSO"])), 2))
check("lactate 21 percent inhibitor mean", "1,345.29", 1345.29,
      round(float(np.mean(mv["lactate|21%|BAY"])), 2))
check("lactate hypoxic inhibitor mean", "703.94", 703.94,
      round(float(np.mean(mv["lactate|0.5%|BAY"])), 2))
check("complete separation of the eligible contrast", "completely separated", True,
      min(mv["lactate|21%|BAY"]) > max(mv["lactate|21%|DMSO"]))

pp = load_json("copeland_2026-09-14/run/run_summary.json")
check("condition optimizations", "all 224 optimizations completed", 224, pp["completed"])
check("condition failures", "without failure", 0, pp["failed"])
prep = load_json("copeland_2026-09-14/prepared/prepare_summary.json")
check("supported reactions", "5,936 of 5,938", 5936, prep["reactions"]["supported_conservative"])
check("reactions with rules", "5,936 of 5,938", 5938, prep["reactions"]["internal_with_gpr"])
check("scale constant", "43.181", 43.181, round(prep["scale_constant_k"], 3))

cr = load_json("copeland_2026-09-14/crossed/crossed_summary.json")
check("crossed jobs", "All 448 optimizations completed", 448, cr["completed"])
d = cr["decomposition"]
check("task spread, lactate, identity", "0.2416", 0.2416,
      round(d["magnitude_g1.00|EX_lac__L_e"]["mean_spread_across_growth_rates_at_fixed_transcript"], 4))
check("transcript spread, lactate, identity", "0.0035", 0.0035,
      round(d["magnitude_g1.00|EX_lac__L_e"]["mean_spread_across_transcripts_at_fixed_growth_rate"], 4))
ratios = [v["ratio_transcript_to_task"] for k, v in d.items()
          if k.startswith("magnitude") and v["ratio_transcript_to_task"]]
check("magnitude ratio range", "between 1.1 and 2.9 percent", True,
      1.0 <= 100 * min(ratios) and 100 * max(ratios) <= 3.0)
ratios_o = [v["ratio_transcript_to_task"] for k, v in d.items()
            if k.startswith("ordinal") and v["ratio_transcript_to_task"]]
check("midrank ratio range", "between 4.6 and 5.2 percent", True,
      4.5 <= 100 * min(ratios_o) and 100 * max(ratios_o) <= 5.3)

# ---------------------------------------------------------------- endpoint reconstruction
ep = load_json("strengthening_reference_checked_2026-09-14/endpoint_audit/coverage_summary.json")
flagged = {v["source_cultures"]["screen_positive"] for v in ep["variants"].values()}
unresolved = {v["source_cultures"]["unresolved"] for v in ep["variants"].values()}
reproduced = {v["matches_Nilsson_84_positive_36_negative"] for v in ep["variants"].values()}
check("cultures flagged", "flagging 119 of 120 cultures", True, flagged == {119})
check("cultures unresolved", "with one unresolved", True, unresolved == {1})
check("published screen not reproduced", "did not recover its reported culture classification",
      True, reproduced == {False})
check("source cultures", "119 of 120 cultures", 120, ep["source_cultures"])

# ---------------------------------------------------------------- report
print("%-46s %-8s %-8s" % ("check", "in text", "value"))
print("-" * 66)
bad = 0
for label, present, ok, claim, exp, act in results:
    flag = "ok" if ok else "MISMATCH"
    tflag = "yes" if present else "NOT FOUND"
    if not ok or not present:
        bad += 1
    print("%-46s %-8s %-8s" % (label[:46], tflag, flag))
    if not ok:
        print("      claim %r expected %r got %r" % (claim, exp, act))
    elif not present:
        print("      expected string not located: %r" % claim)
print("-" * 66)
print("%d checks, %d problems" % (len(results), bad))

uncovered = [
    "Table 4 per-scenario sensitivity counts (340 and 383 of 2,444; 2,507 and 2,575 reversals)",
    "Table 5 per-encoding resolved-pair and width counts",
    "development-stage results reported in Supplementary Methods S9",
    "runtime figures quoted in the supplement",
]
print("\nUNCOVERED by this checker, stated so the gap is visible:")
for u in uncovered:
    print("  -", u)
sys.exit(1 if bad else 0)
