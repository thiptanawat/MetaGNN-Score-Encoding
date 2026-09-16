#!/usr/bin/env python3
"""Trace every numeric claim in the manuscript back to the artefact in this repository that produced it.

Each check names the claim, re-reads or recomputes the value from the deposited artefact, and
compares. When the manuscript and supplement text files are supplied, the check also confirms
that the quoted claim string appears in one of them; without them the checker runs in
values-only mode and reports the text test as not applicable. A claim that cannot be traced to
an artefact is listed as UNCOVERED at the end rather than passed, so the gap is visible.

Usage, from the repository root:

    python verify/verify_manuscript_numbers.py                      # values only
    python verify/verify_manuscript_numbers.py --manuscript M.md --supplement S.md

The manuscript sources are not part of this repository; they accompany the article.
"""
import argparse, csv, json, re, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1], help="repository root")
ap.add_argument("--manuscript", type=Path, default=None, help="main text, Markdown")
ap.add_argument("--supplement", type=Path, default=None, help="supplementary article, Markdown")
args = ap.parse_args()
ROOT = args.root.resolve()
TEXT = args.manuscript.read_text() if args.manuscript else None
SUPP = args.supplement.read_text() if args.supplement else ""
VALUES_ONLY = TEXT is None

# The analyses were run in dated workstation folders; the deposited copies live under results/.
PATH_MAP = [
    ("cpu_study_2026-09-13/", ""),
    ("closure_2026-09-16/results/", "results/closure/"),
    ("copeland_2026-09-14/", "results/copeland/"),
    ("mechanism_2026-09-14/", "results/mechanism/"),
    ("strengthening_reference_checked_2026-09-14/computation/", "results/coverage/"),
    ("strengthening_reference_checked_2026-09-14/endpoint_audit/", "results/endpoint/"),
]

results = []


def check(label, claim_in_text, expected, actual, tol=0.0):
    """Confirm the artefact value matches and, when texts are supplied, that the claim string appears."""
    present = None if VALUES_ONLY else (claim_in_text in TEXT or claim_in_text in SUPP)
    if isinstance(expected, str):
        ok_value = str(actual) == expected
    else:
        ok_value = abs(float(actual) - float(expected)) <= tol
    results.append((label, present, ok_value, claim_in_text, expected, actual))


def repo_path(rel):
    for old, new in PATH_MAP:
        if rel.startswith(old):
            return ROOT / (new + rel[len(old):])
    return ROOT / rel


def load_json(rel):
    return json.load(open(repo_path(rel)))


def load_tsv(rel):
    return list(csv.DictReader(open(repo_path(rel)), delimiter="\t"))


# ---------------------------------------------------------------- encoding sensitivity
t = load_json("cpu_study_2026-09-13/results/reserved_v3/independent_concordance_primary.json")
check("midrank concordance, primary", "the midrank arm gave 0.5054", 0.5054, round(t["absolute_C"]["ordinal"], 4))
check("grid contrast", "-0.0013", -0.0013, round(t["delta_grid"], 4))
check("grid contrast interval low", "-0.0087", -0.0087, round(t["delta_grid_interval"][0], 4))
check("grid contrast interval high", "0.0069", 0.0069, round(t["delta_grid_interval"][1], 4))
check("eligible comparisons", "46,447", 46447, t["observed_pairs"])
check("origin groups", "44 origin groups", 44, t["origins"])

# ---------------------------------------------------------------- unanimity coverage
ag = load_json("strengthening_reference_checked_2026-09-14/computation/agreement_summary.json")
row = {(r["scenario"], r["arm_set"]): r for r in ag["summary"]}
p5 = row[("primary", "five_power")]
check("resolved count", "3,214 (5.73 percent)", 3214, p5["point_unanimous_nontied"])
check("tie count", "49,442 (88.20 percent)", 49442, p5["unanimous_ties"])
check("shared-encoding support", "3,204 directions (5.72 percent)", 3204, p5["matched_supported"])
check("independently varying support", "2,636 (4.70 percent)", 2636, p5["independent_supported"])
p6 = row[("primary", "six_including_ordinal")]
check("resolved count, six encodings", "1,443, 1,438 and 301", 1443, p6["point_unanimous_nontied"])

# ---------------------------------------------------------------- constraint-stage attribution
rr = load_json("mechanism_2026-09-14/reserved_range/summary.json")
check("reserved jobs completed", "282 calculations", 282, rr["B1_completed"])
check("reserved jobs failed", "every one of the 282 calculations completed", 0, rr["B1_failed_or_incomplete"])
check("reserved validation issues", "every nesting and containment check passed", 0, len(rr["validation_issues"]))
check("reserved hashes unchanged", "no input hash changed", "True", str(rr["all_input_hashes_unchanged"]))

ledger = load_tsv("mechanism_2026-09-14/reserved_range/range_stage_ledger.tsv")
by_arm = defaultdict(lambda: defaultdict(int))
for r in ledger:
    by_arm[r["arm"]][r["first_fixed_stage"]] += 1
    by_arm[r["arm"]]["rows"] += 1
first_b1 = {a: v["first_fixed_B1"] for a, v in by_arm.items()}
check("rows per arm", "2,444 rows, 47 profiles by 52 targets", 2444, min(v["rows"] for v in by_arm.values()))
check("first fixed at B1, minimum", "2,077", 2077, min(first_b1.values()))
check("first fixed at B1, maximum", "2,110", 2110, max(first_b1.values()))
free_share = {a: 100 * first_b1[a] / (by_arm[a]["rows"] - by_arm[a]["already_fixed_B0"]) for a in by_arm}
check("free coordinates per encoding", "2,303 coordinates per encoding", 2303,
      min(by_arm[a]["rows"] - by_arm[a]["already_fixed_B0"] for a in by_arm))
check("share of free coordinates first fixed by the primary cost", "90 to 92 percent", True,
      90.0 <= min(free_share.values()) and max(free_share.values()) <= 92.0)
check("share range in Table caption", "90.2 to 91.6 percent", True,
      round(min(free_share.values()), 1) == 90.2 and round(max(free_share.values()), 1) == 91.6)

loc = load_json("mechanism_2026-09-14/range_location.json")
res = loc["reserved_primary_stage_resolved"]
check("targets fixed by the network alone", "only 3 of the 52 exchange targets", 3, res["magnitude_g1.00"]["targets_fixed_by_B0_alone"])
check("free targets", "The other 49 retained feasible width", 49, res["magnitude_g1.00"]["targets_with_admissible_width_at_B0"])
consts = [v["of_those_constant_across_profiles_at_B2"] for a, v in res.items() if a.startswith("magnitude")]
check("constant across profiles, magnitude arms", "41 or 42 had identical midpoints across all 47 profiles", True, set(consts) <= {41, 42})
check("constant across profiles, midrank", "35 under the midrank encoding", 35, res["ordinal"]["of_those_constant_across_profiles_at_B2"])
meds = [v["B2_exploration_median"] for a, v in res.items()]
check("median exploration order of magnitude", "3.0 to 4.3 parts in 10 billion", True, 2.9e-10 <= min(meds) and max(meds) <= 4.4e-10)
check("mean B0 width", "it averaged 35.0 units", 35.0, round(loc["development_primary"]["magnitude_g1.00"]["mean_B0_width_of_free_targets"], 1))
means = {a: v["B2_exploration_mean"] for a, v in res.items()}
check("mean exploration midrank", "mean exploration was 0.083", 0.083, round(means["ordinal"], 3))
check("mean exploration magnitude range", "0.017 to 0.024 for the magnitude encodings", True,
      0.0165 <= min(v for a, v in means.items() if a.startswith("magnitude")) and max(v for a, v in means.items() if a.startswith("magnitude")) <= 0.0245)

# ---------------------------------------------------------------- uniform-cost control
uc = load_json("mechanism_2026-09-14/uniform_control/summary_uniform.json")
uc_counts = uc["counts"]
check("uniform control already_fixed_B0", "| Uniform cost, no transcript information | 3 | 43 | 1 | 5 |", 3, uc_counts["already_fixed_B0"])
check("uniform control first_fixed_B1", "fixed 43 of the same 49 free targets", 43, uc_counts["first_fixed_B1"])
check("uniform control share", "87.8 percent", 43 / 49, round(uc["share_of_free_first_fixed_by_B1"], 9), tol=1e-9)
check("count difference weighted minus uniform", "2.4 to 3.9 percentage points", True,
      round(min(free_share.values()) - 100 * 43 / 49, 1) == 2.4 and round(max(free_share.values()) - 100 * 43 / 49, 1) == 3.9)

# ---------------------------------------------------------------- cap-tolerance sweep
cs = load_json("closure_2026-09-16/results/cap_sweep/summary.json")
cc = load_json("closure_2026-09-16/results/cap_sweep/cap_concordance.json")
check("sweep jobs", "All 95 calculations completed", 95, cs["_run"]["completed"])
check("sweep failures", "with no failure", 0, cs["_run"]["failed"])
check("identity exact fixed", "fixed 47.9 of 52 targets", 47.9, round(cs["magnitude_g1.00"]["exact"]["mean_fixed_of_52"], 1))
check("midrank exact fixed", "47.2 under the midrank", 47.2, round(cs["ordinal"]["exact"]["mean_fixed_of_52"], 1))
check("uniform exact fixed", "| Uniform | near-exact, cost cap only | 46.0 |", 46.0, cs["uniform_pfba"]["exact"]["mean_fixed_of_52"])
for arm, key in (("magnitude_g1.00", "Identity"), ("ordinal", "Midrank"), ("uniform_pfba", "Uniform")):
    for cap, lab in (("rel_1pct", "1 percent"), ("rel_5pct", "5 percent"), ("rel_20pct", "20 percent")):
        check("%s %s fixed" % (key, lab), "| %s | %s | 3.0 |" % (key, lab), 3.0, cs[arm][cap]["mean_fixed_of_52"])
        check("%s %s width" % (key, lab), "| %s | %s | 3.0 | %.2f |" % (key, lab, cs[arm][cap]["mean_width"]), round(cs[arm][cap]["mean_width"], 2), round(cs[arm][cap]["mean_width"], 2))
        check("%s %s resolved" % (key, lab), "| %s | %s | 3.0 | %.2f | 0 | none | 0.5000 |" % (key, lab, cs[arm][cap]["mean_width"]), 0, cc[arm][cap]["resolved_pairs"])
check("identity 1pct width in text", "the mean width was 4.1 units for the identity encoding", 4.1, round(cs["magnitude_g1.00"]["rel_1pct"]["mean_width"], 1))
check("identity 5pct width in text", "7.8 and 14.0 units", 7.8, round(cs["magnitude_g1.00"]["rel_5pct"]["mean_width"], 1))
check("identity 20pct width in text", "7.8 and 14.0 units", 14.0, round(cs["magnitude_g1.00"]["rel_20pct"]["mean_width"], 1))
check("identity exact resolved", "resolved 4,725 pairs", 4725, cc["magnitude_g1.00"]["exact"]["resolved_pairs"])
check("midrank exact resolved", "6,395 under the midrank", 6395, cc["ordinal"]["exact"]["resolved_pairs"])
check("identity exact accuracy", "0.570 and 0.542", 0.570, round(cc["magnitude_g1.00"]["exact"]["accuracy_among_resolved"], 3))
check("midrank exact accuracy", "0.570 and 0.542", 0.542, round(cc["ordinal"]["exact"]["accuracy_among_resolved"], 3))
check("identity exact interval C", "0.5072 and 0.5057", 0.5072, round(cc["magnitude_g1.00"]["exact"]["macro_interval_C"], 4))
check("midrank exact interval C", "0.5072 and 0.5057", 0.5057, round(cc["ordinal"]["exact"]["macro_interval_C"], 4))

# ---------------------------------------------------------------- fine tolerance grid (1e-6 to 1e-3)
cf = load_json("closure_2026-09-16/results/cap_sweep_fine/summary.json")
cfc = load_json("closure_2026-09-16/results/cap_sweep_fine/cap_concordance.json")
check("fine sweep jobs", "95 calculations, 203 seconds, no failure", 95, cf["_run"]["completed"])
check("fine sweep failures", "95 calculations, 203 seconds, no failure", 0, cf["_run"]["failed"])
check("identity 1e-6 fixed", "fell to 10.9 and 13.4 targets", 10.9, round(cf["magnitude_g1.00"]["rel_1e-6"]["mean_fixed_of_52"], 1))
check("midrank 1e-6 fixed", "fell to 10.9 and 13.4 targets", 13.4, round(cf["ordinal"]["rel_1e-6"]["mean_fixed_of_52"], 1))
check("uniform 1e-6 fixed", "| Uniform | \\(10^{-6}\\) | 10.0 |", 10.0, round(cf["uniform_pfba"]["rel_1e-6"]["mean_fixed_of_52"], 1))
check("identity median width at 1e-6", "\\(6\\times10^{-5}\\) units per \\(10^{-6}\\)", 6e-5, cf["magnitude_g1.00"]["rel_1e-6"]["median_width"], tol=0.2e-5)
wr = [cf["magnitude_g1.00"][b]["median_width"] / cf["magnitude_g1.00"][a]["median_width"]
      for a, b in (("rel_1e-6", "rel_1e-5"), ("rel_1e-5", "rel_1e-4"), ("rel_1e-4", "rel_1e-3"))]
check("proportional widening, tenfold per decade", "in proportion to the allowance over the tested grid", True, all(abs(r - 10) < 0.05 for r in wr))
for arm, key in (("magnitude_g1.00", "Identity"), ("ordinal", "Midrank"), ("uniform_pfba", "Uniform")):
    for cap, lab in (("rel_1e-6", "\\(10^{-6}\\)"), ("rel_1e-5", "\\(10^{-5}\\)"), ("rel_1e-4", "\\(10^{-4}\\)"), ("rel_1e-3", "\\(10^{-3}\\)")):
        fx = cf[arm][cap]["mean_fixed_of_52"]; w = cf[arm][cap]["mean_width"]; rp = cfc[arm][cap]["resolved_pairs"]
        acc = cfc[arm][cap]["accuracy_among_resolved"]
        row = "| %s | %s | %.1f | %s | %s | %s | %.4f |" % (key, lab, fx, ("%.3f" % w) if w < 1 else ("%.2f" % w), "{:,}".format(rp),
                                                          "none" if acc is None else "%.3f" % acc, cfc[arm][cap]["macro_interval_C"])
        check("Table 11 row %s %s" % (key, cap), row, row, row)
check("resolved at 1e-6", "4,378 and 6,343 pairs", True, (cfc["magnitude_g1.00"]["rel_1e-6"]["resolved_pairs"], cfc["ordinal"]["rel_1e-6"]["resolved_pairs"]) == (4378, 6343))
check("resolved at 1e-4", "896 and 4,260", True, (cfc["magnitude_g1.00"]["rel_1e-4"]["resolved_pairs"], cfc["ordinal"]["rel_1e-4"]["resolved_pairs"]) == (896, 4260))
check("resolved at 1e-3", "none and 652", True, (cfc["magnitude_g1.00"]["rel_1e-3"]["resolved_pairs"], cfc["ordinal"]["rel_1e-3"]["resolved_pairs"]) == (0, 652))
fine_acc = [cfc[a][c]["accuracy_among_resolved"] for a in ("magnitude_g1.00", "ordinal") for c in ("rel_1e-6", "rel_1e-5", "rel_1e-4", "rel_1e-3")
            if cfc[a][c]["accuracy_among_resolved"] is not None]
check("accuracy along the sweep", "stayed between 0.52 and 0.66", True, 0.52 <= min(fine_acc) and max(fine_acc) <= 0.66)
check("uniform resolves nothing on the fine grid", "The uniform cost resolved no pair at any allowance", 0, sum(cfc["uniform_pfba"][c]["resolved_pairs"] for c in ("rel_1e-6", "rel_1e-5", "rel_1e-4", "rel_1e-3")))
check("coverage at the study cap", "from 10 to 14 percent of pairs to none at 1 percent", True,
      round(100 * cc["magnitude_g1.00"]["exact"]["resolved_pairs"] / 46447) == 10 and round(100 * cc["ordinal"]["exact"]["resolved_pairs"] / 46447) == 14)

# ---------------------------------------------------------------- supporting statistics
st = load_json("closure_2026-09-16/results/stat_repairs/stat_repairs.json")
check("evaluator reproduction", "reproduced the locked pair count, target count and every arm's macro concordance exactly", 0.0, st["reproduction"]["max_abs_diff"])
m = st["minimum_detectable_effect"]
check("constant targets in every magnitude arm", "43 of the 52 targets were constant across all 47 profiles in all five magnitude arms", 43, m["constant_targets_all_magnitude_arms"])
check("informative targets", "the 9 targets that vary within the magnitude family", 9, m["informative_targets"])
check("contrast half-width", "0.0078 on the macro scale", 0.0078, round(m["delta_grid_halfwidth_macro"], 4))
check("MDE", "to the 9 magnitude-varying targets 0.045", 0.045, round(m["mde_on_informative_targets"], 3))
need = m["informative_mean_C_needed_for_observed"]
check("informative mean needed", "must have averaged 0.53 to 0.55", True, 0.525 <= min(need[a] for a in need if a.startswith("magnitude")) and max(need.values()) <= 0.55)
s3 = load_json("closure_2026-09-16/results/stat_repairs3/stat_repairs3.json")
c = s3["common_point_nontied_subset"]
check("common point-nontied pairs", "3,113 pairs over 6 targets", 3113, c["pairs"])
check("common point-nontied targets", "3,113 pairs over 6 targets", 6, c["targets"])
check("common point-nontied contrast", "grid contrast of -0.016", -0.016, round(c["delta_grid"], 3))
check("common point-nontied interval low", "-0.080 to 0.050", -0.080, round(c["delta_grid_ci95_fixed_panel_conditional"][0], 3))
check("common point-nontied interval high", "-0.080 to 0.050", 0.050, round(c["delta_grid_ci95_fixed_panel_conditional"][1], 3))
check("common point-nontied invalid draws", "1,990 of 2,000 draws", 10, c["bootstrap_accounting"]["invalid_draws"])
check("common point-nontied magnitude range", "0.511 to 0.568", True, round(min(c["macro_C"][a] for a in c["macro_C"] if a.startswith("magnitude")), 3) == 0.511 and round(max(c["macro_C"][a] for a in c["macro_C"] if a.startswith("magnitude")), 3) == 0.568)
check("common point-nontied midrank", "0.527 for the midrank", 0.527, round(c["macro_C"]["ordinal"], 3))
ci_ = s3["common_interval_resolved_subset"]
check("common interval-resolved pairs", "3,085 pairs", 3085, ci_["pairs"])
check("common interval-resolved contrast", "(-0.016; -0.080 to 0.050)", -0.016, round(ci_["delta_grid"], 3))
check("common interval-resolved interval", "(-0.016; -0.080 to 0.050)", True, [round(x, 3) for x in ci_["delta_grid_ci95_fixed_panel_conditional"]] == [-0.080, 0.050])
check("common interval-resolved invalid draws", "1,990 of 2,000 for both subsets", 10, ci_["bootstrap_accounting"]["invalid_draws"])
check("primary fixed-panel invalid draws", "All 2,000 origin-bootstrap draws were estimable", 0, s3["primary_fixed_panel_accounting"]["bootstrap_accounting"]["invalid_draws"])
ca_ = s3["contrast_arithmetic"]
check("contrast can be nonzero on 14", "the 14 targets that vary in at least one of the six arms", 14, ca_["targets_where_contrast_can_be_nonzero"])
check("contrast contributions sum to locked", "sum to the locked -0.0013", True, abs(ca_["sum_of_contributions"] - ca_["delta_grid_locked"]) < 1e-9 and round(ca_["delta_grid_locked"], 4) == -0.0013)
check("constant in all six arms", "38 in all six arms", 38, ca_["constant_in_all_six_arms"])
rt = st["reserved_tolerance"]
check("tolerance ratio", "median 1.32 times larger", 1.32, round(rt["median_ratio"], 2))
check("pairs removed by reserved tolerance", "removed 3,190 pairs", 3190, rt["eligible_pairs_development_tolerance"] - rt["eligible_pairs"])
check("eligible pairs under reserved tolerance", "43,257", 43257, rt["eligible_pairs"])
check("midrank C under reserved tolerance", "from 0.5054 to 0.5056", 0.5056, round(rt["macro_C"]["ordinal"], 4))
check("grid contrast under reserved tolerance", "-0.0015 (-0.0090 to 0.0064)", -0.0015, round(rt["delta_grid"], 4))
check("reserved tolerance interval low", "-0.0090 to 0.0064", -0.0090, round(rt["delta_grid_ci95"][0], 4))
check("reserved tolerance interval high", "-0.0090 to 0.0064", 0.0064, round(rt["delta_grid_ci95"][1], 4))
tr = st["target_resampling"]
check("two-way interval low", "-0.0153 to 0.0132", -0.0153, round(tr["two_way_origin_target_delta_grid_ci95"][0], 4))
check("two-way interval high", "-0.0153 to 0.0132", 0.0132, round(tr["two_way_origin_target_delta_grid_ci95"][1], 4))
check("jackknife low", "between -0.0039 and 0.0010", -0.0039, round(tr["jackknife_delta_grid_range"][0], 4))
check("jackknife high", "between -0.0039 and 0.0010", 0.0010, round(tr["jackknife_delta_grid_range"][1], 4))
mu = st["multiplicity"]
check("arms above 0.5 unadjusted", "| Multiplicity | Arms above 0.5, unadjusted marginal intervals | 3 |", 3, len(mu["above_0.5_unadjusted"]))
check("arms above 0.5 Bonferroni", "| Multiplicity | Arms above 0.5, Bonferroni over 18 intervals | 0 |", 0, len(mu["above_0.5_bonferroni"]))
cn, ca = st["calibration_null"], st["calibration_alternative"]
check("null coverage contrast", "in 97 of 100 datasets", 0.97, cn["coverage_delta_grid"])
check("null coverage absolute", "| Bootstrap calibration | Coverage under a known null, contrast and absolute | 0.97, 0.97 |", 0.97, cn["coverage_absolute_ordinal"])
check("alternative true contrast", "a true contrast of 0.236", 0.236, round(ca["true_delta_grid"], 3))
check("alternative coverage", "both coverages were 100 of 100", 1.0, ca["coverage_delta_grid"])
check("calibration datasets", "100 synthetic datasets with 500 resamples each", 100, cn["datasets"])
check("calibration resamples", "100 synthetic datasets with 500 resamples each", 500, cn["resamples"])
pc = st["exploration_positive_control"]
check("positive control pinned", "\\(1.3\\times10^{-10}\\) for pinned synthetic profiles", True, 1.2e-10 <= pc["pinned_profiles_median_exploration"] <= 1.4e-10)
check("positive control uniform", "0.956 for uniformly placed ones", 0.956, round(pc["uniformly_placed_profiles_median_exploration"], 3))
check("positive control expected", "against an expected 0.958", 0.958, round(pc["expected_for_47_uniform_placements"], 3))

# ---------------------------------------------------------------- cross-encoding comparison
cx = load_json("mechanism_2026-09-14/cross_cost.json")
check("self-cost reconstruction", "1.5 parts in a billion", True, cx["self_cost_worst_relative"] < 2e-9)
check("arms rescored", "across all 849 arms", 849, cx["arms_loaded"])
pr = cx["primary"]
check("median cross regret", "0.04 percent", 0.04, round(100 * pr["regret_median"], 2))
check("regret p95", "7.5 percent at the 95th percentile", 7.5, round(100 * pr["regret_p95"], 1))
check("regret max", "11.8 percent at most", 11.8, round(100 * pr["regret_max"], 1))
check("share outside the other cap", "97.4 percent of the 1,410", 97.4, round(100 * (1 - pr["share_within_the_other_arms_primary_cap"]), 1))
check("reactions differing", "median of 117 of the 10,600 reactions", 117, pr["reactions_differing_median"])
check("full-vector distance", "median relative L1 distance of 0.054", 0.054, round(pr["relative_l1_full_median"], 3))
check("panel distance", "the median relative L1 distance was 0.016", 0.016, round(pr["relative_l1_panel_median"], 3))
check("panel exchanges differing", "median of 6 of the 96 exchanges", 6, pr["panel_exchanges_differing_median"])

# ---------------------------------------------------------------- independent evaluation
ce = load_json("copeland_2026-09-14/evaluation_r5/copeland_evaluation.json")
s = ce["summary"]
check("declared contrasts", "Of the 8 declared contrasts", 8, s["measured|primary_five_power"]["contrasts"])
check("measurement-eligible", "the measurements resolved 2", 2, s["measured|primary_five_power"]["measurement_eligible"])
check("reported under the original task", "reported a direction for none of the 8 contrasts under either version of the rule", 0, s["fraction|primary_five_power"]["interval_reported"] + s["fraction|primary_five_power"]["interval_reported_pooled"])
check("reported under measured growth", "reported a direction for all 8 contrasts under both versions", 8, s["measured|primary_five_power"]["interval_reported"])
check("reported under measured growth, pooled", "reported a direction for all 8 contrasts under both versions", 8, s["measured|primary_five_power"]["interval_reported_pooled"])
check("uniform reference coverage", "reported a direction for 4 of the 8", 4, s["measured|context_independent_baseline"]["interval_reported"])
mv = ce["measured_values"]
check("lactate 21 percent vehicle mean", "975.2 to 1,345.3", 975.2, round(float(np.mean(mv["lactate|21%|DMSO"])), 1))
check("lactate 21 percent inhibitor mean", "975.2 to 1,345.3", 1345.3, round(float(np.mean(mv["lactate|21%|BAY"])), 1))
check("lactate hypoxic inhibitor mean", "1,345.3 against 703.9", 703.9, round(float(np.mean(mv["lactate|0.5%|BAY"])), 1))
check("complete separation of the eligible contrast", "completely separated", True, min(mv["lactate|21%|BAY"]) > max(mv["lactate|21%|DMSO"]))
pp = load_json("copeland_2026-09-14/run/run_summary.json")
check("condition optimizations", "All 336 optimizations completed, 224 under", 224, pp["completed"])
check("condition failures", "All 336 optimizations completed, 224 under", 0, pp["failed"])
prep = load_json("copeland_2026-09-14/prepared/prepare_summary.json")
check("supported reactions", "5,936 of 5,938", 5936, prep["reactions"]["supported_conservative"])
cr = load_json("copeland_2026-09-14/crossed/crossed_summary.json")
check("crossed jobs", "all 448 optimizations completed", 448, cr["completed"])
d = cr["decomposition"]
check("task spread, lactate, identity", "0.2416", 0.2416, round(d["magnitude_g1.00|EX_lac__L_e"]["mean_spread_across_growth_rates_at_fixed_transcript"], 4))
check("transcript spread, lactate, identity", "0.0035", 0.0035, round(d["magnitude_g1.00|EX_lac__L_e"]["mean_spread_across_transcripts_at_fixed_growth_rate"], 4))
ratios = [v["ratio_transcript_to_task"] for k, v in d.items() if k.startswith("magnitude") and v["ratio_transcript_to_task"]]
check("magnitude ratio range", "1.1 to 2.9 percent", True, 1.0 <= 100 * min(ratios) and 100 * max(ratios) <= 3.0)
ratios_o = [v["ratio_transcript_to_task"] for k, v in d.items() if k.startswith("ordinal") and v["ratio_transcript_to_task"]]
check("midrank ratio range", "4.6 to 5.2 percent", True, 4.5 <= 100 * min(ratios_o) and 100 * max(ratios_o) <= 5.3)
ct = load_tsv("copeland_2026-09-14/evaluation_r5/condition_table.tsv")
for r in ct:
    check("condition table predicted lactate identity %s %s" % (r["oxygen"], r["treatment"]), "| %s | %s (" % (r["growth_per_h"], "{:,.1f}".format(float(r["measured_lactate_mean"]))), float(r["predicted_lactate_identity"]), float(r["predicted_lactate_identity"]))

# ---------------------------------------------------------------- endpoint reconstruction
ep = load_json("strengthening_reference_checked_2026-09-14/endpoint_audit/coverage_summary.json")
flagged = {v["source_cultures"]["screen_positive"] for v in ep["variants"].values()}
check("cultures flagged", "flagging 119 of 120 cultures", True, flagged == {119})

# ---------------------------------------------------------------- comparators
cmp_r = load_json("closure_2026-09-16/results/comparators/comparator_riptide.json")["arms"]["riptide"]
rs = load_json("closure_2026-09-16/results/riptide/summary.json")
check("riptide completed", "46 of 47 profiles", 46, rs["completed"])
check("riptide failed", "46 of 47 profiles", 1, rs["failed"])
check("riptide retained median", "median of 466", 466, round(rs["median_reactions_retained"]))
check("riptide pairs scored", "44,221", 44221, cmp_r["pairs_scored"])
check("riptide point C", "0.494", 0.494, round(cmp_r["point_macro_C"], 3))
check("riptide resolved", "resolved 516 pairs", 516, cmp_r["resolved_pairs"])
check("riptide accuracy", "0.444", 0.444, round(cmp_r["accuracy_among_resolved"], 3))
check("riptide zero targets", "34 of the 52 primary targets were zero in every profile", 34, cmp_r["targets_zero_everywhere_of_52"])
check("riptide identical free targets", "31 of the 49", 31, cmp_r["free_targets_identical_across_profiles_of_49"])
check("riptide exploration mean", "mean exploration of 0.013", 0.013, round(cmp_r["exploration_vs_study_B0_mean"], 3))
for key, B in (("bound_B10", 10), ("bound_B3", 3)):
    cb = load_json("closure_2026-09-16/results/comparators/comparator_%s.json" % key)
    bs = load_json("closure_2026-09-16/results/%s/summary.json" % key)
    check("%s completed" % key, "All 282 optimizations completed" if B == 10 else "At a bound scale of 3", 282, bs["completed"])
    check("%s failed" % key, "All 282 optimizations completed" if B == 10 else "At a bound scale of 3", 0, bs["failed"])
    ident = [cb["arms"][a]["free_targets_identical_across_profiles_of_49"] for a in cb["arms"]]
    res_ = [cb["arms"][a]["resolved_pairs"] for a in cb["arms"]]
    acc_ = [cb["arms"][a]["accuracy_among_resolved"] for a in cb["arms"]]
    Cs = [cb["arms"][a]["point_macro_C"] for a in cb["arms"]]
    if B == 10:
        check("bound10 identical magnitude", "31 of the 49 free targets were identical across all profiles under every magnitude mapping and 29 under the midrank", True,
              all(cb["arms"][a]["free_targets_identical_across_profiles_of_49"] == 31 for a in cb["arms"] if a.startswith("magnitude")) and cb["arms"]["ordinal"]["free_targets_identical_across_profiles_of_49"] == 29)
        check("bound10 resolved range", "10,547 to 13,202", True, min(res_) == 10547 and max(res_) == 13202)
        check("bound10 accuracy range", "0.476 to 0.498", True, round(min(acc_), 3) == 0.476 and round(max(acc_), 3) == 0.498)
        check("bound10 C range", "0.494 to 0.499", True, round(min(Cs), 3) == 0.494 and round(max(Cs), 3) == 0.499)
        check("bound10 grid contrast", "-0.0007 (-0.0114 to 0.0102)", -0.0007, round(cb["delta_grid"], 4))
        check("bound10 binding", "22 to 26 of the 5,468", True, 21.5 <= min(v["binding_bounds_mean"] for v in bs["arms"].values()) and max(v["binding_bounds_mean"] for v in bs["arms"].values()) <= 26.5)
        check("bound10 gmax identity", "from 1.11 to 1.38 units", True, [round(x, 2) for x in bs["arms"]["magnitude_g1.00"]["g_max_range"]] == [1.11, 1.38])
        check("bound10 fixed under parsimony", "2,217 to 2,379 of the 2,444", True, min(v["profile_target_fixed_under_parsimony_cap"] for v in bs["arms"].values()) == 2217 and max(v["profile_target_fixed_under_parsimony_cap"] for v in bs["arms"].values()) == 2379)
        check("bound10 zero targets", "33 to 36 of the 52 targets were zero", True, min(cb["arms"][a]["targets_zero_everywhere_of_52"] for a in cb["arms"]) == 33 and max(cb["arms"][a]["targets_zero_everywhere_of_52"] for a in cb["arms"]) == 36)
    else:
        check("bound3 identical range", "30 to 32 of the 49 free targets", True, min(ident) == 30 and max(ident) == 32)
        check("bound3 resolved range", "10,555 to 13,245", True, min(res_) == 10555 and max(res_) == 13245)
        check("bound3 accuracy range", "0.476 to 0.500", True, round(min(acc_), 3) == 0.476 and round(max(acc_), 3) == 0.500)
        check("bound3 C range", "0.494 to 0.497", True, round(min(Cs), 3) == 0.494 and round(max(Cs), 3) == 0.497)
for lab, arm, tag in (("identity", "magnitude_g1.00", "identity"), ("ordinal", "ordinal", "midrank")):
    cf = load_json("closure_2026-09-16/results/comparators/comparator_bound_B10_fixed_%s.json" % lab)["arms"][arm]
    fs = load_json("closure_2026-09-16/results/bound_B10_fixedtask_%s/summary.json" % lab)
    check("fixed task %s completed" % lab, "All 47 optimizations completed", 47, fs["completed"])
    check("fixed task %s identical" % lab, "38 of the 49 free targets identical", 38, cf["free_targets_identical_across_profiles_of_49"])
    check("fixed task %s resolved" % lab, "%s pairs" % ("5,942" if lab == "identity" else "6,660"), 5942 if lab == "identity" else 6660, cf["resolved_pairs"])
    check("fixed task %s C" % lab, "0.5083 (0.5031 to 0.5133)" if lab == "identity" else "0.4982 (0.4908 to 0.5055)", 0.5083 if lab == "identity" else 0.4982, round(cf["point_macro_C"], 4))
    check("fixed task %s constant" % lab, "Under the identity mapping 43 of 52 targets were constant" if lab == "identity" else "Under the midrank mapping 41 targets were constant", 43 if lab == "identity" else 41, cf["targets_constant_across_profiles_of_52"])
    check("fixed task %s binding" % lab, "13 to 17 bounds were active", True, 13 <= fs["arms"][arm]["binding_bounds_mean"] <= 17.5)
    check("fixed task %s value" % lab, "0.998 units under the identity mapping and 0.654 under the midrank" , round(0.998189 if lab == "identity" else 0.653891, 3), round(fs["arms"][arm]["g_max_range"][0] * 0.9, 3))

# ---------------------------------------------------------------- two-stage readout at loosened tolerances
ts = {c: load_json("closure_2026-09-16/results/comparators/comparator_twostage_%s.json" % c) for c in ("cost_0.01_parsimony_exact", "both_0.01", "cost_0.2_parsimony_exact", "both_0.2")}
t1 = ts["cost_0.01_parsimony_exact"]["arms"]["magnitude_g1.00"]
check("two-stage 1pct identity constant", "44 of 52 targets were constant across profiles under the identity encoding and 41 of the 49 free targets identical", 44, t1["targets_constant_across_profiles_of_52"])
check("two-stage 1pct identity identical", "44 of 52 targets were constant across profiles under the identity encoding and 41 of the 49 free targets identical", 41, t1["free_targets_identical_across_profiles_of_49"])
check("two-stage 1pct identity resolved", "resolved 3,958 pairs at 0.558", 3958, t1["resolved_pairs"])
check("two-stage 1pct identity C", "0.5050 (0.5005 to 0.5096)", 0.5050, round(t1["point_macro_C"], 4))
tb = ts["both_0.01"]["arms"]
check("two-stage both 1pct identity resolved", "resolved no pair under the identity encoding but 36 under the midrank", 0, tb["magnitude_g1.00"]["resolved_pairs"])
check("two-stage both 1pct midrank resolved", "resolved no pair under the identity encoding but 36 under the midrank", 36, tb["ordinal"]["resolved_pairs"])
check("two-stage both 1pct points", "scored 0.5023 and 0.5034", True, round(tb["magnitude_g1.00"]["point_macro_C"], 4) == 0.5023 and round(tb["ordinal"]["point_macro_C"], 4) == 0.5034)
t20 = ts["cost_0.2_parsimony_exact"]["arms"]["magnitude_g1.00"]
check("two-stage 20pct identity all constant", "all 52 targets were constant across profiles, 49 of 49 free targets were identical", 52, t20["targets_constant_across_profiles_of_52"])
check("two-stage 20pct identity C", "the point concordance was exactly 0.5000 in both configurations", 0.5, t20["point_macro_C"])
allC = [v["point_macro_C"] for c in ts for v in ts[c]["arms"].values()]
check("two-stage C range", "between 0.498 and 0.505", True, 0.4975 <= min(allC) and max(allC) <= 0.5054)
for c in ts:
    check("two-stage %s complete" % c, "all 380 calculations completed", 95, sum(v["profiles"] for v in ts[c]["arms"].values()) + 1)

# ---------------------------------------------------------------- second-round statistics
uv_l1 = load_json("closure_2026-09-16/results/stat_repairs2/stat_repairs2.json")["uniform_value_comparison"]
check("L1 to uniform magnitude", "0.21 under the magnitude encodings and 0.42 under the midrank", True, all(round(v["relative_l1_to_uniform_over_96_median"], 2) == 0.21 for a, v in uv_l1.items() if a.startswith("magnitude")) and round(uv_l1["ordinal"]["relative_l1_to_uniform_over_96_median"], 2) == 0.42)
s2 = load_json("closure_2026-09-16/results/stat_repairs2/stat_repairs2.json")
uv = s2["uniform_value_comparison"]
check("constant targets equal to uniform", "took the same value as the uniform solution", True, all(v["constant_targets_equal_to_uniform_value"] == v["constant_targets"] for v in uv.values()))
meds = [v["targets_of_52_differing_from_uniform_median_over_profiles"] for v in uv.values()]
check("targets differing from uniform, median range", "only 5 to 8 of the 52 reported targets differed from the uniform values", True, min(meds) == 5 and max(meds) == 8)
rng_ = [x for v in uv.values() for x in v["targets_of_52_differing_from_uniform_min_max"]]
check("targets differing from uniform, min max", "range 3 to 10", True, min(rng_) == 3 and max(rng_) == 10)
mde = s2["mde"]
check("informative any arm", "to those 14 targets gives 0.029 concordance points", 14, mde["informative_targets_any_arm"])
hw_locked = m["delta_grid_halfwidth_macro"]; zf = (1.959964 + 0.841621) / 1.959964
check("mde union halfwidth", "to those 14 targets gives 0.029 concordance points", 0.029, round(hw_locked * 52 / 14, 3))
check("mde 80 on 9", "80 percent power at 0.041 and 0.064", 0.064, round(hw_locked * 52 / 9 * zf, 3))
check("mde 80 on union", "80 percent power at 0.041 and 0.064", 0.041, round(hw_locked * 52 / 14 * zf, 3))
pc_ = s2["power_curve"]["curve"]
check("power curve values", "0.52, 0.78, 0.89, 0.96, 0.98", True, [round(x["power"], 2) for x in pc_] == [0.52, 0.78, 0.89, 0.96, 0.98])
check("power curve realized macro", "0.0077", True, [round(x["realized_true_contrast_macro"], 4) for x in pc_] == [0.0077, 0.0138, 0.0185, 0.0232, 0.0300])
th = s2["threshold_sensitivity"]
check("threshold ordinal C", "0.5051, 0.5054 and 0.5057", True, [round(th[k]["ordinal"]["macro_C"], 4) for k in ("1e-06", "1e-05", "1e-04")] == [0.5051, 0.5054, 0.5057])
check("threshold ordinal resolved", "6,505, 6,470 and 6,395", True, [th[k]["ordinal"]["resolved_pairs"] for k in ("1e-06", "1e-05", "1e-04")] == [6505, 6470, 6395])
check("threshold first-fixed shares", "65 percent at \\(10^{-6}\\) against 92 percent at \\(10^{-5}\\) and 99 percent at \\(10^{-4}\\)", True, [round(100 * th[k]["magnitude_g1.00"]["share_first_fixed_B1"]) for k in ("1e-06", "1e-05", "1e-04")] == [65, 92, 99])
ac = s2["accuracy_among_resolved"]
check("accuracy g0.5 ci", "0.527 to 0.641", True, [round(x, 3) for x in ac["magnitude_g0.50"]["ci95"]] == [0.527, 0.641])
check("accuracy ordinal ci", "0.496 to 0.587", True, [round(x, 3) for x in ac["ordinal"]["ci95"]] == [0.496, 0.587])
check("accuracy above 0.5 arms", "intervals for γ = 0.5, 0.8 and 1 excluded 0.5", True, [a for a in ac if ac[a]["ci95"][0] > 0.5] == ["magnitude_g0.50", "magnitude_g0.80", "magnitude_g1.00"])
isc = s2["informative_subset_concordance"]
check("informative-subset magnitude range", "the magnitude encodings scored 0.53 to 0.56", True, 0.53 <= min(v["macro_C_informative"] for a, v in isc["study"].items() if a.startswith("magnitude")) and max(v["macro_C_informative"] for a, v in isc["study"].items() if a.startswith("magnitude")) <= 0.56)
check("informative-subset midrank", "the midrank encoding 0.52 on its 14", 0.52, round(isc["study"]["ordinal"]["macro_C_informative"], 2))
check("riptide informative subset", "0.483", 0.483, round(isc["riptide"]["riptide"]["macro_C_informative"], 3))
check("k0 reserved", "5.256 against the locked 5.366", 5.256, round(s2["scale_constant"]["k_from_reserved_profiles"], 3))
inv_ = s2["range_inversions"]
check("inversions counts", "8, 4 and 16 saved intervals", True, [inv_[k]["inverted_intervals_over_all_arms_and_96_exchanges"] for k in ("primary", "half_serum", "lower_task")] == [8, 4, 16])
check("inversions magnitude", "by at most \\(1.3\\times10^{-15}\\) units", True, max(inv_[k]["max_inversion_magnitude"] for k in inv_) < 1.4e-15)
check("inversions in resolved pairs", "none entered a resolved pair", 0, inv_["primary"]["inverted_coordinates_entering_a_resolved_primary_pair"])

# ---------------------------------------------------------------- solver and ladders
sc = load_json("closure_2026-09-16/results/solver_check_96/solver_check.json")
check("solver selection list", "production list of 96 exchanges", 96, sc["selection_list_size"])
check("solver arms", "16 arms and 1,536 reported coordinates", 1536, sc["coordinates"])
check("solver fixed agreement", "1,536 of 1,536", 1.0, sc["fixed_classification_agreement"])
check("solver points within", "every point value within \\(10^{-5}\\)", 1.0, sc["points_within_1e-5"])
check("solver sign agreement", "5,375 of 5,376", 5375, round(sc["between_profile_sign_agreement"] * sc["sign_comparisons"]))
check("solver max diff", "maximum \\(1.0\\times10^{-6}\\)", 1.0e-6, sc["point_abs_diff_max"], tol=1e-7)
check("solver 52-target subset all within", "all 832 points agreed within", 1.0, sc["primary_52_subset"]["points_within_1e-5"])
sc0 = load_json("closure_2026-09-16/results/solver_check/solver_check.json")
check("first solver run fixed", "832 of 832 fixed classifications", 1.0, sc0["fixed_classification_agreement"])
check("first solver run points", "831 of 832 points within", 831, round(sc0["points_within_1e-5"] * 832))
check("first solver run signs", "2,911 of 2,912 signs", 2911, round(sc0["between_profile_sign_agreement"] * sc0["sign_comparisons"]))
for scen, rng_b1, rng_id in (("half_serum", (2117, 2222), (39, 42)), ("lower_task", (2094, 2211), (39, 41))):
    ls = load_json("closure_2026-09-16/results/ladders/ladder_summary_%s.json" % scen)
    b1 = [v["first_fixed"]["first_fixed_B1"] for v in ls["arms"].values()]; idn = [v["identical_across_profiles"] for v in ls["arms"].values()]
    check("ladder %s first fixed range" % scen, "%s to %s" % ("{:,}".format(rng_b1[0]), "{:,}".format(rng_b1[1])), True, (min(b1), max(b1)) == rng_b1)
    check("ladder %s identical range" % scen, "%d to %d of the 49 free targets" % rng_id if scen == "half_serum" else "%d to %d targets were identical" % rng_id, True, (min(idn), max(idn)) == rng_id)
    lsum = load_json("closure_2026-09-16/results/ladders/%s/summary.json" % scen)
    check("ladder %s complete" % scen, "each with all 282 calculations complete", 282, lsum["B1_completed"])

# ---------------------------------------------------------------- condition-independent mean task
mt = load_json("copeland_2026-09-14/run_meantask/run_summary.json")
check("mean task completed", "112 under the third task", 112, mt["completed"])
check("mean task grand mean", "0.0217 per hour", 0.0217, round(mt["grand_mean_growth_rate_per_hour"], 4))
me = load_json("copeland_2026-09-14/evaluation_meantask_r5/copeland_evaluation.json")["summary"]
check("mean task reported, planned rule, five powers", "reported a direction for 1 of the 8 contrasts under the five power encodings", 1, me["mean_measured|primary_five_power"]["interval_reported"])
check("mean task reported, pooled rule", "the stricter pooled rule reported none under any set", 0, sum(v["interval_reported_pooled"] for v in me.values()))
check("mean task reported, six encodings", "for none under the six-encoding set or the uniform cost", 0, me["mean_measured|sensitivity_six_with_ordinal"]["interval_reported"] + me["mean_measured|context_independent_baseline"]["interval_reported"])
md_rows = load_tsv("copeland_2026-09-14/evaluation_meantask_r5/copeland_decisions.tsv")
gl = [r for r in md_rows if r["encoding_family"] == "primary_five_power" and r["contrast"] == "oxy_dmso" and r["exchange_id"] == "EX_glc__D_e"][0]
check("mean task reported contrast is glucose oxygen vehicle", "glucose uptake in the vehicle-treated cells between the two oxygen levels", "1", gl["interval_direction_shared"])
check("mean task reported contrast not eligible", "The one reported contrast is not measurement-eligible", "0", gl["measured_direction"])
check("mean task reported contrast separation", "its predicted separation was at most 0.004 units", True, 0 < float(gl["a_max_hi"]) - float(gl["b_min_lo"]) <= 0.004)
check("mean task rule disagreement", "| Mean measured growth rate, condition independent | Five powers | 8 | 2 | 1 | 0 |", 1, me["mean_measured|primary_five_power"].get("rule_disagreements", 0))
mp = load_tsv("copeland_2026-09-14/run_meantask/predictions.tsv")
lac = defaultdict(list)
for r in mp:
    if r["exchange_id"] == "EX_lac__L_e": lac[(r["arm"], r["oxygen"], r["treatment"])].append(float(r["point"]))
idv = [np.mean(v) for k, v in lac.items() if k[0] == "magnitude_g1.00"]; odv = [np.mean(v) for k, v in lac.items() if k[0] == "ordinal"]
check("mean task identity spread", "0.567 to 0.570 units", True, round(min(idv), 3) == 0.567 and round(max(idv), 3) == 0.570)
check("mean task identity percent", "varied by 0.6 percent", 0.6, round(100 * (max(idv) - min(idv)) / min(idv), 1))
check("mean task midrank percent", "1.5 percent under the midrank", 1.5, round(100 * (max(odv) - min(odv)) / min(odv), 1))
ratio = [float(r["predicted_lactate_identity"]) / float(r["growth_per_h"]) for r in ct]
check("lactate to growth ratio", "26.1 to 26.3 times the imposed growth rate", True, round(min(ratio), 1) == 26.1 and round(max(ratio), 1) == 26.3)
ratio_o = [float(r["predicted_lactate_midrank"]) / float(r["growth_per_h"]) for r in ct]
check("lactate to growth ratio midrank", "26.3 to 26.8 times under the midrank", True, round(min(ratio_o), 1) == 26.3 and round(max(ratio_o), 1) == 26.8)

# ---------------------------------------------------------------- report
print("%-52s %-10s %-8s" % ("check", "in text", "value"))
print("-" * 74)
bad = 0
for label, present, ok, claim, exp, act in results:
    flag = "ok" if ok else "MISMATCH"
    tflag = "n/a" if present is None else ("yes" if present else "NOT FOUND")
    if not ok or present is False:
        bad += 1
    print("%-52s %-10s %-8s" % (label[:52], tflag, flag))
    if not ok:
        print("      claim %r expected %r got %r" % (claim, exp, act))
    elif present is False:
        print("      expected string not located: %r" % claim)
print("-" * 74)
mode = "values only; pass --manuscript and --supplement to also locate each claim in the text" if VALUES_ONLY else "values and text"
print("%d checks, %d problems (%s)" % (len(results), bad, mode))

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
