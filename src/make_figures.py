"""Plot saved analysis summaries only; no outcomes, model solves, or inference.

Presentation source is manifested separately from the scientific protocol lock.
Missing or non-reportable primary intervals cause a clear failure, not imputation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np

MAGNITUDE = [f"magnitude_g{g:.2f}" for g in (0.5, 0.8, 1.0, 1.25, 2.0)]
ARM_ORDER = ["uniform_pfba", *MAGNITUDE, "ordinal"]
SCENARIOS = ["primary", "half_serum", "lower_task"]
SCENARIO_LABELS = {"primary": "Primary: 90% maximum biomass", "half_serum": "Half assumed serum uptake",
                   "lower_task": "Lower task: 50% maximum biomass"}
ORANGE, DARK, GRAY, LIGHT = "#D75C13", "#693711", "#6B7280", "#D1D5DB"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def number(value, label):
    if isinstance(value, bool) or value is None:
        raise ValueError(f"Missing numeric summary: {label}")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"Nonfinite numeric summary: {label}")
    return value


def interval(value, label):
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise ValueError(f"Missing reportable interval: {label}")
    lo, hi = (number(v, label) for v in value)
    if lo > hi:
        raise ValueError(f"Reversed interval: {label}")
    return lo, hi


def arm_label(arm):
    if arm == "ordinal":
        return "Ordinal costs"
    if arm == "uniform_pfba":
        return "Uniform-cost reference"
    gamma = float(arm.removeprefix("magnitude_g"))
    return f"Magnitude, γ = {gamma:g}" + (" (identity)" if gamma == 1 else "")


def prepare_plot_data(evaluation, sensitivity, analysis_plan):
    """Validate matching, already-computed summaries; return display values.

    Simple descriptive count-to-percentage conversions are formatting only.
    No bootstrap quantile, concordance, prediction, or primary contrast is made here.
    """
    for key in ("scenario", "partition", "analysis_plan_sha256"):
        if not evaluation.get(key) or evaluation[key] != sensitivity.get(key):
            raise ValueError(f"Evaluation/sensitivity mismatch: {key}")
    if evaluation.get("prediction_sha256") != sensitivity.get("predictions_sha256"):
        raise ValueError("Evaluation/sensitivity prediction hashes differ")
    panel = evaluation["primary_fixed_panel"]
    if panel.get("status") != "complete_predictions" or panel.get("primary") is None:
        raise ValueError("Primary panel is not complete; no primary figure can be produced")
    bootstrap = panel.get("bootstrap") or {}
    if bootstrap.get("interval_status") != "conditional_on_estimability":
        raise ValueError("Bootstrap intervals are not reportable")
    primary = panel["primary"]
    rows = {r["arm"]: r for r in panel["arms"]}
    absolute_intervals = bootstrap.get("absolute_C_intervals") or {}
    absolute = []
    for arm in ARM_ORDER:
        row = rows.get(arm)
        if row is None or row.get("failed_point_pairs") != 0:
            if arm == "uniform_pfba":
                continue  # The optional reference must not become a failure-as-tie result.
            raise ValueError(f"Missing successful primary arm: {arm}")
        value = primary["absolute_macro_C"].get(arm)
        if arm == "uniform_pfba":
            value = row["macro_point_operational"]  # Valid only after zero-failure check above.
        value = number(value, f"absolute C {arm}")
        lo, hi = interval(absolute_intervals.get(arm), f"absolute C {arm}")
        if not (0 <= value <= 1 and 0 <= lo <= hi <= 1):
            raise ValueError(f"Concordance outside [0,1]: {arm}")
        absolute.append({"arm": arm, "C": value, "lo": lo, "hi": hi})
    contrasts = {}
    for contrast in ("delta_grid", "delta_identity"):
        lo, hi = interval(bootstrap.get(f"{contrast}_interval"), contrast)
        contrasts[contrast] = {"estimate": number(primary.get(contrast), contrast), "lo": lo, "hi": hi}
    desc = sensitivity["primary_fixed_panel"]
    point_tolerance = number(analysis_plan.get("point_tolerance"), "frozen prediction-tie tolerance")
    interval_tolerance = number(analysis_plan.get("interval_tolerance"), "frozen range-resolution tolerance")
    transform_tolerance = number(desc.get("prediction_tolerance"), "transform prediction tolerance")
    if point_tolerance <= 0 or interval_tolerance <= 0 or transform_tolerance != point_tolerance:
        raise ValueError("Transform threshold differs from the frozen analysis plan")
    ranges = sensitivity["range_diagnostics_primary"]
    if ranges.get("point_tolerance") != point_tolerance or ranges.get("interval_tolerance") != interval_tolerance:
        raise ValueError("Range diagnostics use different frozen thresholds")
    panel_size = len(panel["contributing_target_panel"])
    if desc["targets"] != panel_size or desc["contexts"] != panel["contexts"]:
        raise ValueError("Descriptive sensitivity uses a different context/target panel")
    all_pairs = number(desc["all_distinct_origin_context_target_pairs"], "all predicted pairs")
    cells = number(desc["context_target_cells"], "context-target cells")
    if all_pairs <= 0 or cells <= 0:
        raise ValueError("Descriptive denominators must be positive")
    range_diagnostics = []
    for arm in ARM_ORDER:
        row = ranges["per_arm"][arm]
        if row["context_target_cells"] != cells or row["total_targets"] != panel_size:
            raise ValueError("Within-context range diagnostics use a different panel")
        wide = number(row["within_context_width_above_tolerance_count"], "wide range count")
        constant = number(row["constant_prediction_targets"], "constant target count")
        zero = number(row["all_zero_prediction_targets"], "all-zero target count")
        if not (0 <= wide <= cells and 0 <= constant <= panel_size and 0 <= zero <= panel_size):
            raise ValueError("Invalid within-context range diagnostic count")
        range_diagnostics.append({"arm": arm, "wide_cells": int(wide),
                                  "constant_targets": int(constant), "zero_targets": int(zero)})
    changed = number(desc["changed_fraction"], "changed fraction")
    reversal = number(desc["strict_order_reversal_fraction"], "strict reversal fraction")
    tie = number(desc["tie_transition_without_strict_reversal_pairs"], "tie transition count") / all_pairs
    if not all(0 <= v <= 1 for v in (changed, reversal, tie)) or reversal + tie > 1 + 1e-12:
        raise ValueError("Invalid descriptive fractions")
    coverage = []
    for arm in ARM_ORDER:
        row = rows.get(arm)
        if row is None or row.get("failed_range_pairs") != 0:
            if arm == "uniform_pfba":
                continue
            raise ValueError(f"Range failure in primary arm: {arm}")
        den = number(row["eligible_pairs"], "observed eligible pairs")
        num = number(row["resolved_interval_pairs"], "resolved interval pairs")
        if den <= 0 or not 0 <= num <= den:
            raise ValueError("Invalid separated-range count")
        coverage.append({"arm": arm, "fraction": num / den, "resolved": int(num), "eligible": int(den)})
    return {"scenario": evaluation["scenario"], "partition": evaluation["partition"],
            "plan_sha256": evaluation["analysis_plan_sha256"], "absolute": absolute,
            "contrasts": contrasts, "coverage": coverage, "changed_fraction": changed,
            "reversal_fraction": reversal, "tie_transition_fraction": tie,
            "all_prediction_pairs": int(all_pairs), "context_target_cells": int(cells),
            "contexts": panel["contexts"], "origins": panel["origins"], "targets": panel_size,
            "point_tolerance": point_tolerance, "interval_tolerance": interval_tolerance,
            "range_diagnostics": range_diagnostics,
            "invalid_draws": bootstrap["invalid_draws"], "resamples": bootstrap["resamples"]}


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 12,
                         "axes.labelsize": 10, "figure.facecolor": "white", "axes.facecolor": "white",
                         "savefig.facecolor": "white", "axes.spines.top": False,
                         "axes.spines.right": False, "axes.edgecolor": "#9CA3AF",
                         "text.color": "#252525", "axes.labelcolor": "#252525",
                         "pdf.fonttype": 42, "ps.fonttype": 42})


def scope_text(data):
    stage = "Development feasibility" if data["partition"] == "development" else "Locked reserved evaluation"
    return f"{stage} · {data['contexts']} contexts / {data['origins']} origins · {data['targets']} fixed targets"


def tolerance_label(value):
    return f"{value:.0e}".replace("e-0", "e-").replace("e+0", "e+")


def finish(fig, out, stem, outputs):
    for suffix in ("png", "pdf"):
        path = out / f"{stem}.{suffix}"
        fig.savefig(path, dpi=300 if suffix == "png" else None)
        outputs[path.name] = sha256(path)
    plt.close(fig)


def draw_absolute(data, out, outputs):
    fig = plt.figure(figsize=(10.8, 5.6))
    full = fig.add_axes([.25, .29, .23, .49])
    detail = fig.add_axes([.57, .29, .39, .49], sharey=full)
    rows = data["absolute"]
    extent = max(.03, max(max(abs(r["lo"] - .5), abs(r["hi"] - .5)) for r in rows) * 1.15)
    for ax in (full, detail):
        for i, row in enumerate(rows):
            color = ORANGE if row["arm"] == "ordinal" else DARK if row["arm"] == "magnitude_g1.00" else GRAY
            ax.hlines(i, row["lo"], row["hi"], color=color, lw=2)
            ax.plot(row["C"], i, "o", color=color, ms=5)
        ax.set_ylim(len(rows) - .5, -.5)
        ax.axvline(.5, color=LIGHT, ls="--", zorder=0)
        ax.set_xlabel("Macro concordance C")
        ax.grid(axis="x", alpha=.16)
    full.set(yticks=range(len(rows)), yticklabels=[arm_label(r["arm"]) for r in rows], xlim=(0, 1), xticks=[0, .5, 1])
    detail.set_xlim(max(0, .5 - extent), min(1, .5 + extent))
    detail.tick_params(axis="y", left=False, labelleft=False)
    full.set_title("A  Full scale", loc="left", fontsize=11)
    detail.set_title("B  Detail around 0.5", loc="left", fontsize=11)
    fig.suptitle("Absolute agreement with extracellular measurements", x=.03, ha="left", y=.965, fontsize=14)
    fig.text(.03, .89, scope_text(data), fontsize=9)
    fig.text(.03, .13, "Dots: point estimates. Lines: marginal 95% origin-bootstrap intervals, conditional on fixed-panel estimability.", fontsize=8.5)
    fig.text(.03, .085, f"Absolute predicted differences ≤ {tolerance_label(data['point_tolerance'])} model units are ties. Invalid draws: {data['invalid_draws']}/{data['resamples']}.", fontsize=8.5)
    fig.text(.03, .04, "Pipeline and targets are fixed; no model-refitting or new-lineage uncertainty is represented.", fontsize=8.5)
    finish(fig, out, "figure1_absolute_concordance", outputs)


def draw_contrasts(data, out, outputs):
    fig, ax = plt.subplots(figsize=(9, 4.9))
    fig.subplots_adjust(left=.31, right=.96, bottom=.27, top=.78)
    values = []
    for i, d in enumerate(data):
        for key, offset, color, marker in (("delta_grid", -.13, ORANGE, "o"), ("delta_identity", .13, DARK, "s")):
            row = d["contrasts"][key]
            ax.hlines(i + offset, row["lo"], row["hi"], color=color, lw=2)
            ax.plot(row["estimate"], i + offset, marker, color=color, ms=6)
            values.extend((row["lo"], row["hi"], row["estimate"]))
    extent = max(.01, max(abs(v) for v in values) * 1.18)
    ax.set(xlim=(-extent, extent), yticks=range(len(data)), yticklabels=[SCENARIO_LABELS[d["scenario"]] for d in data])
    ax.set_ylim(len(data)-.55, -.55)
    ax.axvline(0, color=GRAY, lw=1, ls="--")
    ax.grid(axis="x", alpha=.16)
    ax.set_xlabel("Concordance difference (positive values favor ordinal costs)")
    handles = [plt.Line2D([], [], color=ORANGE, marker="o", label="Ordinal − mean of five magnitude arms"),
               plt.Line2D([], [], color=DARK, marker="s", label="Ordinal − identity magnitude arm")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(.025, .11), frameon=False, fontsize=9)
    title = "Paired contrasts across prespecified model settings" if len(data) > 1 else "Paired contrasts in the prespecified primary setting"
    fig.suptitle(title, x=.03, ha="left", y=.965, fontsize=14)
    fig.text(.03, .89, scope_text(data[0]), fontsize=9)
    fig.text(.03, .06, "Marginal 95% intervals: same origin resamples within each scenario; conditional on fixed-panel estimability.", fontsize=8.5)
    note = "Sensitivity settings are not independent replications. " if len(data) > 1 else ""
    fig.text(.03, .025, note + "No equivalence or noninferiority margin is asserted.", fontsize=8.5)
    finish(fig, out, "figure2_paired_contrasts", outputs)


def draw_sensitivity(data, out, outputs):
    primary = next(d for d in data if d["scenario"] == "primary")
    fig = plt.figure(figsize=(11.8, 8.3))
    axes = [fig.add_axes(bounds) for bounds in ((.09, .57, .33, .25), (.60, .57, .33, .25),
                                                (.245, .23, .28, .235), (.59, .23, .385, .235))]
    labels = [{"primary": "Primary", "half_serum": "Half serum", "lower_task": "Lower task"}[d["scenario"]] for d in data]
    x = np.arange(len(data))
    axes[0].bar(x, [d["changed_fraction"] for d in data], color=ORANGE, width=.6)
    axes[0].set(xticks=x, xticklabels=labels, ylim=(0, 1), title="A  Point sensitivity", ylabel="Context–target cells changed (%)")
    axes[1].bar(x, [d["reversal_fraction"] for d in data], color=ORANGE, width=.6, label="Strict reversal")
    axes[1].bar(x, [d["tie_transition_fraction"] for d in data], bottom=[d["reversal_fraction"] for d in data], color=LIGHT, width=.6, label="Tie change only")
    axes[1].set(xticks=x, xticklabels=labels, ylim=(0, 1), title="B  Order sensitivity", ylabel="All different-origin predicted pairs (%)")
    axes[1].legend(frameon=False, fontsize=8, loc="upper right")
    for ax in axes[:2]:
        ax.yaxis.set_major_formatter(PercentFormatter(1))
        ax.grid(axis="y", alpha=.16)
    rows = primary["coverage"]
    axes[2].barh(range(len(rows)), [r["fraction"] for r in rows], color=[ORANGE if r["arm"] == "ordinal" else GRAY for r in rows], height=.6)
    axes[2].set(yticks=range(len(rows)), yticklabels=[arm_label(r["arm"]) for r in rows], xlim=(0, 1), title="C  Between-context separation")
    axes[2].invert_yaxis()
    axes[2].xaxis.set_major_formatter(PercentFormatter(1))
    axes[2].set_xlabel("Observed-eligible pairs (%)")
    axes[2].grid(axis="x", alpha=.16)
    diagnostics = primary["range_diagnostics"]
    def short_label(arm):
        return {"ordinal": "Ordinal", "uniform_pfba": "Uniform"}.get(arm, "γ = " + arm.removeprefix("magnitude_g"))
    body = [[short_label(row["arm"]), str(row["wide_cells"]), str(row["constant_targets"]), str(row["zero_targets"])] for row in diagnostics]
    headings = ["Encoding", f"Wide / {primary['context_target_cells']}", f"Constant / {primary['targets']}", f"Zero / {primary['targets']}"]
    table = axes[3].table(cellText=body, colLabels=headings, cellLoc="center", bbox=(0, 0, 1, 1), colWidths=[.24, .23, .30, .23])
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    for (row, col), cell in table.get_celld().items():
        cell.set_facecolor("white")
        cell.set_edgecolor("#E5E7EB")
        if row == 0:
            cell.get_text().set_fontweight("bold")
        elif diagnostics[row-1]["arm"] == "ordinal":
            cell.get_text().set_color(ORANGE)
            cell.get_text().set_fontweight("bold")
    axes[3].set_title("D  Within-context width and constancy")
    axes[3].set_axis_off()
    fig.suptitle("Descriptive sensitivity, range separation and constancy", x=.03, ha="left", y=.965, fontsize=14)
    fig.text(.03, .89, scope_text(primary), fontsize=9)
    fig.text(.03, .145, "A–B: changes across five fixed-rank magnitude transforms; no observed-outcome filter or accuracy interpretation.", fontsize=8.5)
    fig.text(.03, .11, f"C–D: primary setting. C uses observed-eligible pairs with range gaps > {tolerance_label(primary['interval_tolerance'])}; these are not confidence intervals.", fontsize=8.5)
    fig.text(.03, .075, f"D: 'Wide' means within-context range width > {tolerance_label(primary['interval_tolerance'])}; 'Constant' means between-context point span ≤ {tolerance_label(primary['point_tolerance'])}.", fontsize=8.5)
    fig.text(.03, .04, f"'Zero' means all point magnitudes ≤ {tolerance_label(primary['point_tolerance'])}. Model units apply. Small range width does not certify exact uniqueness.", fontsize=8.5)
    finish(fig, out, "figure3_descriptive_sensitivity", outputs)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument("--results-dir", type=Path, required=True, help="Contains evaluation_SCENARIO.json and transform_SCENARIO.json")
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--scenarios", nargs="+", choices=SCENARIOS, default=SCENARIOS)
    args = ap.parse_args(argv)
    if "primary" not in args.scenarios or len(args.scenarios) != len(set(args.scenarios)):
        raise ValueError("Use each scenario once and include primary")
    ordered = [s for s in SCENARIOS if s in args.scenarios]
    plan_path = args.root / "protocol/analysis_plan.json"
    plan = json.loads(plan_path.read_text())
    plan_hash = sha256(plan_path)
    inputs, data = {str(plan_path.resolve()): plan_hash}, []
    for scenario in ordered:
        paths = [args.results_dir / f"{kind}_{scenario}.json" for kind in ("evaluation", "transform")]
        values = [json.loads(p.read_text()) for p in paths]
        if values[0].get("analysis_plan_sha256") != plan_hash:
            raise ValueError("Summary was not generated under the current frozen analysis plan")
        inputs.update({str(p.resolve()): sha256(p) for p in paths})
        d = prepare_plot_data(*values, plan)
        if d["scenario"] != scenario:
            raise ValueError("Scenario does not match summary filename")
        data.append(d)
    for key in ("partition", "plan_sha256", "contexts", "origins", "targets", "point_tolerance", "interval_tolerance"):
        if len({d[key] for d in data}) != 1:
            raise ValueError(f"Cross-scenario mismatch: {key}")
    if data[0]["partition"] == "test":
        from src.verify_inputs import verify_protocol_lock
        verify_protocol_lock(args.root)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    style()
    outputs = {}
    draw_absolute(data[0], args.output_dir, outputs)
    draw_contrasts(data, args.output_dir, outputs)
    draw_sensitivity(data, args.output_dir, outputs)
    manifest = {"source_sha256": {"src/make_figures.py": sha256(__file__)},
                "input_sha256": inputs, "output_sha256": outputs, "partition": data[0]["partition"],
                "scenarios": ordered, "analysis_plan_sha256": data[0]["plan_sha256"],
                "interpretation": "Plots saved summaries only; no statistical estimation or outcome construction.",
                "display_data": data}
    (args.output_dir / "figure_manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"output_directory": str(args.output_dir), "files": sorted(outputs)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
