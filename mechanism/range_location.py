#!/usr/bin/env python3
"""Q1: does the cost cap pin every profile to the SAME coordinate, or to different ones?

The constraint-stage pilot established that most within-profile range contraction happens when
the primary transcript-weighted cost cap is imposed. That is a statement about interval WIDTH.
It does not say where the surviving interval sits, so it cannot by itself explain why
predictions barely differ between profiles. This script measures the location.

Inputs, all frozen and hashed:
  bottleneck_pilot/range_stage_ledger.tsv   B0/B1/B2 intervals, 11 development profiles,
                                            identity and ordinal encodings, primary scenario
  cpu_study/results/reserved_v3/predictions.tsv   saved B2 intervals for the reserved cohort

No CORE outcome is read and no optimization is run.
"""
import csv, hashlib, json, os, sys
from collections import defaultdict
import numpy as np
from mechanism import WIDTH_TOL, width, location, is_fixed, spread, exploration

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
PILOT = os.path.join(BASE, "strengthening_reference_checked_2026-09-14", "bottleneck_pilot",
                     "range_stage_ledger.tsv")
RESERVED_LEDGER = os.path.join(HERE, "reserved_range", "range_stage_ledger.tsv")
RESERVED = os.path.join(BASE, "cpu_study_2026-09-13", "results", "reserved_v3", "predictions.tsv")
PANEL = os.path.join(BASE, "cpu_study_2026-09-13", "protocol", "analysis_plan.json")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_pilot(path=PILOT):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    by = defaultdict(dict)                       # (arm, exchange) -> context -> stage intervals
    for r in rows:
        by[(r["arm"], r["exchange_id"])][r["context"]] = {
            st: (float(r["%s_low" % st]), float(r["%s_high" % st])) for st in ("B0", "B1", "B2")}
    return rows, by


def analyse_pilot(by):
    """Per target and encoding, compare how wide each stage's interval is with how far the
    profiles move apart. The two are different quantities and the distinction is the point."""
    out, b0_mismatch = [], 0
    for (arm, ex), per_ctx in sorted(by.items()):
        ctxs = sorted(per_ctx)
        b0s = {per_ctx[c]["B0"] for c in ctxs}
        if len(b0s) != 1:                        # B0 must be profile independent by construction
            b0_mismatch += 1
        b0_lo, b0_hi = per_ctx[ctxs[0]]["B0"]
        rec = {"arm": arm, "exchange_id": ex, "profiles": len(ctxs),
               "B0_low": b0_lo, "B0_high": b0_hi, "B0_width": width(b0_lo, b0_hi),
               "B0_shared": len(b0s) == 1}
        for st in ("B1", "B2"):
            locs = [location(*per_ctx[c][st]) for c in ctxs]
            wids = [width(*per_ctx[c][st]) for c in ctxs]
            rec["%s_mean_width" % st] = float(np.mean(wids))
            rec["%s_max_width" % st] = float(np.max(wids))
            rec["%s_location_spread" % st] = spread(locs)
            rec["%s_exploration" % st] = exploration(locs, b0_lo, b0_hi)
            rec["%s_all_fixed" % st] = bool(all(w <= WIDTH_TOL for w in wids))
            rec["%s_constant_across_profiles" % st] = bool(spread(locs) <= WIDTH_TOL)
        out.append(rec)
    return out, b0_mismatch


def summarise_pilot(recs):
    s = {}
    for arm in sorted({r["arm"] for r in recs}):
        a = [r for r in recs if r["arm"] == arm]
        free = [r for r in a if r["B0_width"] > WIDTH_TOL]      # network left room to vary
        pinned = [r for r in free if r["B2_constant_across_profiles"]]
        moving = [r for r in free if not r["B2_constant_across_profiles"]]
        expl = [r["B2_exploration"] for r in free if np.isfinite(r["B2_exploration"])]
        s[arm] = {
            "targets": len(a),
            "targets_fixed_by_B0_alone": sum(1 for r in a if r["B0_width"] <= WIDTH_TOL),
            "targets_with_admissible_width_at_B0": len(free),
            "of_those_constant_across_profiles_at_B2": len(pinned),
            "of_those_varying_across_profiles_at_B2": len(moving),
            "B2_exploration_median": float(np.median(expl)) if expl else None,
            "B2_exploration_mean": float(np.mean(expl)) if expl else None,
            "B2_exploration_max": float(np.max(expl)) if expl else None,
            "B2_exploration_p90": float(np.percentile(expl, 90)) if expl else None,
            "mean_B0_width_of_free_targets": float(np.mean([r["B0_width"] for r in free])),
            "mean_B2_within_profile_width_of_free_targets": float(
                np.mean([r["B2_mean_width"] for r in free])),
            "mean_B2_location_spread_of_free_targets": float(
                np.mean([r["B2_location_spread"] for r in free])),
            "varying_targets": sorted(r["exchange_id"] for r in moving),
        }
    return s


def analyse_reserved():
    """The same location question on the reserved cohort, from saved B2 intervals."""
    panel = {t["exchange_id"] for t in json.load(open(PANEL))["chemistry_targets"]}
    by = defaultdict(lambda: defaultdict(dict))      # scenario -> (arm, ex) -> ctx -> (lo, hi)
    with open(RESERVED) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["exchange_id"] in panel and r["status"] == "optimal":
                by[r["scenario"]][(r["arm"], r["exchange_id"])][r["context_id"]] = (
                    float(r["range_lo"]), float(r["range_hi"]))
    out = {}
    for scen in sorted(by):
        per_arm = defaultdict(list)
        for (arm, ex), per_ctx in by[scen].items():
            locs = [location(*v) for v in per_ctx.values()]
            wids = [width(*v) for v in per_ctx.values()]
            per_arm[arm].append({
                "exchange_id": ex, "profiles": len(per_ctx),
                "location_spread": spread(locs),
                "mean_width": float(np.mean(wids)), "max_width": float(np.max(wids)),
                "constant_across_profiles": bool(spread(locs) <= WIDTH_TOL),
                "all_fixed_within_profile": bool(all(w <= WIDTH_TOL for w in wids))})
        out[scen] = {arm: {
            "targets": len(v),
            "constant_across_profiles": sum(1 for r in v if r["constant_across_profiles"]),
            "all_fixed_within_profile": sum(1 for r in v if r["all_fixed_within_profile"]),
            "narrow_within_but_varying_across": sum(
                1 for r in v if r["all_fixed_within_profile"] and not r["constant_across_profiles"]),
            "mean_within_profile_width": float(np.mean([r["mean_width"] for r in v])),
            "mean_across_profile_location_spread": float(np.mean([r["location_spread"] for r in v])),
            "median_across_profile_location_spread": float(
                np.median([r["location_spread"] for r in v])),
        } for arm, v in sorted(per_arm.items())}
        out["%s__per_target" % scen] = {a: v for a, v in per_arm.items()}
    return out


def main():
    rows, by = load_pilot()
    recs, mismatch = analyse_pilot(by)
    pilot_summary = summarise_pilot(recs)
    reserved = analyse_reserved()

    # the same stage-resolved location analysis on the 47 reserved profiles and six encodings
    res_summary, res_mismatch, res_rows, res_recs = None, None, 0, []
    if os.path.exists(RESERVED_LEDGER):
        res_rows_list, res_by = load_pilot(RESERVED_LEDGER)
        res_recs, res_mismatch = analyse_pilot(res_by)
        res_summary = summarise_pilot(res_recs)
        res_rows = len(res_rows_list)
        with open("range_location_reserved_stages.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(res_recs[0]), delimiter="\t")
            w.writeheader(); w.writerows(res_recs)

    with open("range_location_pilot.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(recs[0]), delimiter="\t")
        w.writeheader(); w.writerows(recs)

    report = {
        "question": "Do the contracted intervals sit at the same location in every profile?",
        "width_tolerance": WIDTH_TOL,
        "ledger_rows_read": len(rows),
        "B0_profile_dependence_violations": mismatch,
        "development_primary": pilot_summary,
        "reserved_primary_stage_resolved": res_summary,
        "reserved_ledger_rows_read": res_rows,
        "reserved_B0_profile_dependence_violations": res_mismatch,
        "reserved_by_scenario": {k: v for k, v in reserved.items() if not k.endswith("__per_target")},
        "inputs_sha256": {os.path.relpath(p, BASE): sha256(p) for p in (PILOT, RESERVED, PANEL)},
        "no_outcome_read": True, "no_lp_solved": True,
    }
    json.dump(report, open("range_location.json", "w"), indent=1, default=float)
    for scen in sorted(k for k in reserved if k.endswith("__per_target")):
        name = "range_location_reserved_%s.tsv" % scen.replace("__per_target", "")
        flat = [dict(arm=a, **r) for a, v in reserved[scen].items() for r in v]
        with open(name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(flat[0]), delimiter="\t")
            w.writeheader(); w.writerows(flat)

    print(json.dumps(report["development_primary"], indent=1, default=float)[:2600])
    print("\nRESERVED, primary scenario")
    for arm, v in report["reserved_by_scenario"]["primary"].items():
        print(" %-16s targets %d | constant across profiles %d | narrow within but varying across %d"
              " | mean within width %.3e | mean location spread %.4f"
              % (arm, v["targets"], v["constant_across_profiles"],
                 v["narrow_within_but_varying_across"], v["mean_within_profile_width"],
                 v["mean_across_profile_location_spread"]))
    print("\nB0 profile-dependence violations: development %s | reserved %s" % (mismatch, res_mismatch))
    if res_summary:
        print("\nRESERVED, stage-resolved (47 profiles, primary scenario)")
        print(" %-16s %7s %7s %7s %11s %11s %11s" % ("arm", "targets", "B0fix", "free",
              "constant", "expl_median", "expl_mean"))
        for arm, v in res_summary.items():
            print(" %-16s %7d %7d %7d %11d %11.3e %11.4f"
                  % (arm, v["targets"], v["targets_fixed_by_B0_alone"],
                     v["targets_with_admissible_width_at_B0"],
                     v["of_those_constant_across_profiles_at_B2"],
                     v["B2_exploration_median"], v["B2_exploration_mean"]))


if __name__ == "__main__":
    main()
