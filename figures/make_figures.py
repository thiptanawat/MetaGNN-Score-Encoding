#!/usr/bin/env python3
"""
make_figures.py -- Figure 4 and Figure 5 for the PeerJ manuscript (paper3).

Figure 4  "Constraint-stage attribution and the position of restricted
           coordinates" (three panels: A stacked bar, B ECDF of widths,
           C exploration distribution).

Figure 5  "Independent condition-response evaluation" (three panels:
           A measured replicates, B predicted lactate vs measured growth
           rate, C crossed profile x growth-rate control).

This script reads ONLY the following inputs (paths resolved relative to
the paper3 base directory, which is taken to be the parent of the
directory this script lives in):

  mechanism_2026-09-14/reserved_range/range_stage_ledger.tsv
  mechanism_2026-09-14/range_location_reserved_stages.tsv
  copeland_2026-09-14/evaluation/copeland_evaluation.json
  copeland_2026-09-14/run/predictions.tsv
  copeland_2026-09-14/evaluation/copeland_decisions.tsv
  copeland_2026-09-14/crossed/crossed_predictions.tsv
  copeland_2026-09-14/crossed/crossed_summary.json

and writes:

  figures/Figure_4_Constraint_Stage.png
  figures/Figure_5_Independent_Evaluation.png
  figures/SHA256SUMS.txt

All numbers shown on the figures (counts, medians, means, spreads) are
computed from the data files at run time -- nothing is hard-coded.
"""
import os
import json
import hashlib

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import MaxNLocator

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(SCRIPT_DIR)  # .../paper3

MECH = os.path.join(BASE, "mechanism_2026-09-14")
COPE = os.path.join(BASE, "copeland_2026-09-14")
OUTDIR = os.path.join(SCRIPT_DIR, "figures")
os.makedirs(OUTDIR, exist_ok=True)

LEDGER_TSV = os.path.join(MECH, "reserved_range", "range_stage_ledger.tsv")
LOCATION_TSV = os.path.join(MECH, "range_location_reserved_stages.tsv")
EVAL_JSON = os.path.join(COPE, "evaluation", "copeland_evaluation.json")
PRED_TSV = os.path.join(COPE, "run", "predictions.tsv")
DECISIONS_TSV = os.path.join(COPE, "evaluation", "copeland_decisions.tsv")
CROSSED_TSV = os.path.join(COPE, "crossed", "crossed_predictions.tsv")
CROSSED_JSON = os.path.join(COPE, "crossed", "crossed_summary.json")

CLIP = 1e-12  # floor applied before log-transforming any width/exploration value

# ---------------------------------------------------------------------------
# Style: plain and neutral. No seaborn, no 3D, no chartjunk.
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 11,
    "axes.labelsize": 12.5,
    "xtick.labelsize": 10.5,
    "ytick.labelsize": 10.5,
    "legend.fontsize": 9.3,
    "axes.linewidth": 0.9,
    "axes.edgecolor": "black",
    "axes.facecolor": "white",
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "hatch.linewidth": 0.7,
})

ARMS6 = ["magnitude_g0.50", "magnitude_g0.80", "magnitude_g1.00",
         "magnitude_g1.25", "magnitude_g2.00", "ordinal"]
ARM_LABELS6 = {
    "magnitude_g0.50": "γ=0.50",
    "magnitude_g0.80": "γ=0.80",
    "magnitude_g1.00": "γ=1.00",
    "magnitude_g1.25": "γ=1.25",
    "magnitude_g2.00": "γ=2.00",
    "ordinal": "Ordinal",
}
ARMS7 = ARMS6 + ["uniform_pfba"]
ARM_LABELS7 = dict(ARM_LABELS6, uniform_pfba="Uniform pFBA")

# Distinct color + distinct marker shape per arm (redundant encoding so the
# arm identity in Figure 5B survives grayscale printing).
OKABE_ITO = ["#E69F00", "#56B4E9", "#009E73", "#F0E442",
             "#0072B2", "#D55E00", "#000000"]
ARM_MARKERS7 = {a: m for a, m in zip(ARMS7, ["o", "s", "^", "v", "D", "P", "X"])}
ARM_COLORS7 = {a: c for a, c in zip(ARMS7, OKABE_ITO)}


def panel_letter(ax, letter):
    ax.text(-0.16, 1.06, letter, transform=ax.transAxes,
            fontsize=17, fontweight="bold", va="bottom", ha="left")


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def png_size(path):
    """Read width/height straight out of the PNG IHDR chunk (no Pillow dependency)."""
    with open(path, "rb") as f:
        data = f.read(33)
    assert data[:8] == b"\x89PNG\r\n\x1a\n", f"not a PNG: {path}"
    w = int.from_bytes(data[16:20], "big")
    h = int.from_bytes(data[20:24], "big")
    return w, h


# ===========================================================================
# FIGURE 4 -- constraint-stage attribution
# ===========================================================================
def make_figure4(report):
    ledger = pd.read_csv(LEDGER_TSV, sep="\t")
    location = pd.read_csv(LOCATION_TSV, sep="\t")

    if len(ledger) != 14664:
        raise ValueError(f"range_stage_ledger.tsv: expected 14664 rows, got {len(ledger)}")

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(14.2, 4.9))

    # ---------------- Panel A: stacked bar of first_fixed_stage counts ----
    cats = ["already_fixed_B0", "first_fixed_B1", "first_fixed_B2", "still_variable_B2"]
    cat_labels = ["Fixed by B0", "First fixed by B1", "First fixed by B2", "Still variable at B2"]
    cat_colors = ["#3f3f3f", "#7d7d7d", "#b7b7b7", "#e9e9e9"]
    cat_hatches = [None, "///", "xxx", "..."]

    counts = (ledger.groupby(["arm", "first_fixed_stage"]).size()
              .unstack(fill_value=0).reindex(index=ARMS6, columns=cats, fill_value=0))
    totals = counts.sum(axis=1)
    for a, t in totals.items():
        if t != 2444:
            raise ValueError(f"arm {a}: expected 2444 profile-target coordinates, got {t}")

    x = np.arange(len(ARMS6))
    bottoms = np.zeros(len(ARMS6))
    b2_top = None
    for cat, color, hatch in zip(cats, cat_colors, cat_hatches):
        vals = counts[cat].values.astype(float)
        axA.bar(x, vals, bottom=bottoms, color=color, edgecolor="black",
                linewidth=0.7, hatch=hatch, width=0.64, zorder=3)
        if cat != "first_fixed_B2":
            # "first_fixed_B2" is skipped here -- its segment is thin and
            # hatched ("xxx"), so an inline count collides with the hatch
            # lines; it is called out above the bar instead (below).
            for xi, (v, b) in enumerate(zip(vals, bottoms)):
                if v >= 90:
                    txt_color = "white" if color in ("#3f3f3f", "#7d7d7d") else "black"
                    axA.text(xi, b + v / 2, f"{int(v):,}", ha="center", va="center",
                              fontsize=8.2, color=txt_color, zorder=4)
        bottoms += vals
        if cat == "first_fixed_B2":
            b2_top = bottoms.copy()

    # "First fixed by B2" counts: called out above the bar with a short
    # leader from the top of that (thin, hatched) segment -- outer row.
    for xi, v in enumerate(counts["first_fixed_B2"].values):
        axA.annotate(f"{int(v)}", xy=(xi, b2_top[xi]), xytext=(xi - 0.20, 2444 + 205),
                     ha="center", va="bottom", fontsize=8.2,
                     arrowprops=dict(arrowstyle="-", lw=0.6, color="black"))

    # "Still variable at B2" counts: kept close to the bar (inner row), and
    # fanned to the opposite side from the B2 leaders above so the two rows
    # of leader lines diverge instead of overlapping near the bar top.
    for xi, v in enumerate(counts["still_variable_B2"].values):
        axA.annotate(f"{int(v)}", xy=(xi, 2444), xytext=(xi + 0.20, 2444 + 95),
                     ha="center", va="bottom", fontsize=8.2,
                     arrowprops=dict(arrowstyle="-", lw=0.6, color="black"))

    axA.set_xticks(x)
    axA.set_xticklabels([ARM_LABELS6[a] for a in ARMS6], rotation=30, ha="right",
                        rotation_mode="anchor")
    axA.set_ylabel("profile–target coordinates")
    axA.set_xlabel("encoding arm")
    axA.set_ylim(0, 2444 + 360)
    axA.set_axisbelow(True)
    axA.yaxis.grid(True, linewidth=0.5, color="0.85", zorder=0)
    handles = [Patch(facecolor=c, edgecolor="black", hatch=h, label=l)
               for c, h, l in zip(cat_colors, cat_hatches, cat_labels)]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, -0.02), fontsize=9.6)
    panel_letter(axA, "A")

    # ---------------- Panel B: ECDF of widths at B0/B1/B2, arm = g1.00 -----
    sub = ledger[ledger.arm == "magnitude_g1.00"]
    if len(sub) != 2444:
        raise ValueError(f"magnitude_g1.00 subset: expected 2444 rows, got {len(sub)}")

    stage_specs = [
        ("B0_width", "B0: network + medium + growth task",
         dict(color="black", linestyle="-", linewidth=1.9)),
        ("B1_width", "B1: + transcript-weighted cost cap",
         dict(color="#4477AA", linestyle="--", linewidth=1.9)),
        ("B2_width", "B2: + parsimony (min total flux) cap",
         dict(color="#CC6677", linestyle=":", linewidth=2.6)),
    ]
    clip_counts = {}
    for col, label, sty in stage_specs:
        vals = sub[col].values.astype(float)
        clip_counts[col] = int((vals < CLIP).sum())
        logv = np.log10(np.clip(vals, CLIP, None))
        xs = np.sort(logv)
        ys = np.arange(1, len(xs) + 1) / len(xs)
        axB.step(xs, ys, where="post", label=label, **sty)

    axB.set_xlabel("log$_{10}$(interval width, model units)\nzeros clipped to $10^{-12}$")
    axB.set_ylabel("cumulative fraction of\nprofile–target coordinates")
    axB.set_ylim(-0.02, 1.05)
    axB.set_axisbelow(True)
    axB.grid(True, linewidth=0.5, color="0.88", zorder=0)
    axB.legend(loc="lower right", frameon=True, framealpha=0.9, fontsize=8.4)
    axB.text(0.02, 0.98, "arm: magnitude, γ=1.00", transform=axB.transAxes,
              fontsize=9, ha="left", va="top", style="italic", color="0.25")
    panel_letter(axB, "B")

    # ---------------- Panel C: exploration distribution per arm -----------
    free = location[location.B0_width > 1e-5].copy()
    if len(free) != 6 * 49:
        raise ValueError(f"free-target subset: expected {6*49} rows, got {len(free)}")
    n_clip_c = int((free["B2_exploration"] < CLIP).sum())
    free["explo_c"] = np.clip(free["B2_exploration"].values, CLIP, None)

    data_by_arm = [free.loc[free.arm == a, "explo_c"].values for a in ARMS6]
    positions = np.arange(1, len(ARMS6) + 1)
    axC.boxplot(data_by_arm, positions=positions, vert=False, widths=0.55,
                patch_artist=True, showfliers=False,
                boxprops=dict(facecolor="#d9d9d9", edgecolor="black", linewidth=0.9),
                medianprops=dict(color="black", linewidth=1.9),
                whiskerprops=dict(color="black", linewidth=0.9),
                capprops=dict(color="black", linewidth=0.9))
    rng = np.random.default_rng(0)
    for pos, vals in zip(positions, data_by_arm):
        jitter = rng.uniform(-0.17, 0.17, size=len(vals))
        axC.scatter(vals, pos + jitter, s=12, facecolor="black", edgecolor="none",
                    alpha=0.45, zorder=3)

    axC.set_xscale("log")
    axC.set_xlim(1e-11, 2.2)
    axC.axvline(1.0, color="black", linewidth=1.1, linestyle="--", zorder=1)
    axC.text(1.0, len(ARMS6) + 0.55, "profiles span full\nadmissible range",
              ha="center", va="bottom", fontsize=8.4)
    n_reserved_profiles = int(free["profiles"].iloc[0])
    if not (free["profiles"] == n_reserved_profiles).all():
        raise ValueError("range_location_reserved_stages.tsv: 'profiles' column is not "
                          "constant across free targets; cannot state a single reserved-"
                          "profile count on the panel.")
    axC.text(0.99, 0.90, f"{n_reserved_profiles} reserved profiles per arm", transform=axC.transAxes,
              fontsize=7.6, ha="right", va="top", style="italic", color="0.3")
    medians = free.groupby("arm")["explo_c"].median().reindex(ARMS6)
    med_text = "median B2 exploration\n" + "\n".join(
        f"{ARM_LABELS6[a]}: {medians[a]:.1e}" for a in ARMS6)
    axC.text(0.02, 0.04, med_text, transform=axC.transAxes, fontsize=7.6,
              ha="left", va="bottom", family="monospace",
              bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                        edgecolor="0.6", alpha=0.92))
    axC.set_yticks(positions)
    axC.set_yticklabels([ARM_LABELS6[a] for a in ARMS6])
    axC.set_ylim(0.3, len(ARMS6) + 1.3)
    axC.set_xlabel("B2 exploration (fraction of B0 width)\nzeros clipped to $10^{-12}$")
    axC.set_ylabel("encoding arm")
    axC.set_axisbelow(True)
    axC.xaxis.grid(True, which="major", linewidth=0.5, color="0.85", zorder=0)
    panel_letter(axC, "C")

    fig.tight_layout(rect=[0, 0.045, 1, 1], w_pad=2.6)
    out_path = os.path.join(OUTDIR, "Figure_4_Constraint_Stage.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    report["figure4"] = {
        "path": out_path,
        "panelA_category_totals_per_arm": {a: int(totals[a]) for a in ARMS6},
        "panelA_counts": {a: {c: int(counts.loc[a, c]) for c in cats} for a in ARMS6},
        "panelB_arm": "magnitude_g1.00",
        "panelB_n_rows": int(len(sub)),
        "panelB_n_clipped_zeros": clip_counts,
        "panelC_n_free_targets": int(len(free)),
        "panelC_n_clipped_zeros": n_clip_c,
        "panelC_median_B2_exploration_by_arm": {a: float(medians[a]) for a in ARMS6},
    }


# ===========================================================================
# FIGURE 5 -- independent condition-response evaluation
# ===========================================================================
CONDITION_ORDER = [("0.5%", "DMSO"), ("0.5%", "BAY"), ("21%", "DMSO"), ("21%", "BAY")]
# Within each oxygen level, DMSO (vehicle) always precedes BAY (treatment),
# so all four x-tick groups read in the same order left to right. Oxygen
# groups two adjacent conditions each (see CONDITION_ORDER): index pairs
# (0,1) share 0.5% O2, (2,3) share 21% O2. Panel A uses a two-tier
# x-axis -- short treatment name on the tick, oxygen level as a group label
# spanning its two ticks underneath -- because "0.5% O2 / DMSO" on every one
# of 8 closely spaced ticks does not fit without adjacent labels colliding.
OXYGEN_GROUPS = [(0, 1, "0.5% O$_2$"), (2, 3, "21% O$_2$")]
TREAT_MARKER = {"DMSO": "o", "BAY": "^"}
TREAT_COLOR = {"DMSO": "#4477AA", "BAY": "#CC6677"}


def _bracket(ax, x1, x2, y, h, text, fontsize=8.6):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="black",
            linewidth=1.0, clip_on=False, zorder=5)
    ax.text((x1 + x2) / 2.0, y + h * 1.08, text, ha="center", va="bottom",
              fontsize=fontsize, zorder=5)


def make_figure5(report):
    with open(EVAL_JSON) as f:
        ev = json.load(f)
    measured = ev["measured_values"]

    pred = pd.read_csv(PRED_TSV, sep="\t")
    decisions = pd.read_csv(DECISIONS_TSV, sep="\t")
    crossed = pd.read_csv(CROSSED_TSV, sep="\t")
    with open(CROSSED_JSON) as f:
        csumm = json.load(f)

    # Panel A carries 8 x-tick categories (vs. a plain scatter in B and C) and
    # needs more physical width for its labels not to touch; width_ratios
    # gives it that directly. (Scaling the abstract per-tick step alone
    # does not help: it enlarges the data range and the pixel spacing in the
    # same proportion, so the on-screen gap between labels barely changes.)
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(15.6, 5.1),
                                          gridspec_kw={"width_ratios": [1.55, 1, 1]})

    # ---------------- Panel A: measured replicate values -------------------
    metabolites = ["glucose", "lactate"]
    step = 1.0   # spacing between adjacent condition ticks
    gap = 3.3    # extra spacing between the glucose and lactate clusters (must be wide
                 # enough that the centered legend box does not reach either cluster --
                 # verified below, after layout, against the actual rendered geometry)
    xpos, tick_pos, tick_lab = {}, [], []
    xa = 0.0
    for met in metabolites:
        for ci, (ox, tr) in enumerate(CONDITION_ORDER):
            xpos[(met, ox, tr)] = xa
            tick_pos.append(xa)
            tick_lab.append(tr)  # oxygen level is added below as a group label
            xa += step
        xa += gap

    rng = np.random.default_rng(1)
    means = {}
    point_x, point_y = [], []  # every plotted replicate, for the legend-overlap check below
    for met in metabolites:
        for (ox, tr) in CONDITION_ORDER:
            key = f"{met}|{ox}|{tr}"
            vals = np.asarray(measured[key], dtype=float)
            x0 = xpos[(met, ox, tr)]
            jitter = rng.uniform(-0.10, 0.10, size=len(vals))
            filled = (ox == "21%")
            axA.scatter(x0 + jitter, vals, marker=TREAT_MARKER[tr],
                        facecolors=TREAT_COLOR[tr] if filled else "none",
                        edgecolors=TREAT_COLOR[tr], linewidths=1.4, s=46, zorder=3)
            point_x.extend(x0 + jitter)
            point_y.extend(vals)
            m = float(vals.mean())
            means[(met, ox, tr)] = m
            axA.plot([x0 - 0.24, x0 + 0.24], [m, m], color="black", linewidth=2.6, zorder=4)

    axA.axhline(0, color="0.35", linewidth=1.0, zorder=1)

    # Determine, from copeland_decisions.tsv, which contrasts are completely
    # separated (measured_direction == 1) rather than assuming which ones.
    dec_measured = decisions[decisions.task_mode == "measured"]
    one_family = dec_measured.encoding_family.iloc[0]
    dec_one = dec_measured[dec_measured.encoding_family == one_family]
    separated = dec_one[dec_one.measured_direction == 1]
    if len(separated) == 0:
        raise ValueError("No completely-separated contrasts found in copeland_decisions.tsv; "
                          "cannot annotate Figure 5 Panel A brackets as specified.")

    lac_max = max(np.asarray(measured[f"lactate|{ox}|{tr}"]).max() for ox, tr in CONDITION_ORDER)
    base_h = lac_max * 0.06
    base_y = lac_max * 1.08
    for i, (_, row) in enumerate(separated.sort_values("contrast").iterrows()):
        met = row["metabolite"]
        x1 = xpos[(met, row["a_oxygen"], row["a_treatment"])]
        x2 = xpos[(met, row["b_oxygen"], row["b_treatment"])]
        y = base_y + i * base_h * 2.6
        label = (f"{row['a_oxygen']} {row['a_treatment']} vs "
                 f"{row['b_oxygen']} {row['b_treatment']}: "
                 f"{row['a_measured_mean']:.0f} vs {row['b_measured_mean']:.0f}, "
                 f"fully separated")
        _bracket(axA, min(x1, x2), max(x1, x2), y, base_h, label)
    top_of_brackets = base_y + len(separated) * base_h * 2.6 + base_h * 2.2

    glc_vals = np.concatenate([np.asarray(measured[f"glucose|{ox}|{tr}"]) for ox, tr in CONDITION_ORDER])
    glc_center = float(np.mean([xpos[("glucose", ox, tr)] for ox, tr in CONDITION_ORDER]))
    lac_center = float(np.mean([xpos[("lactate", ox, tr)] for ox, tr in CONDITION_ORDER]))
    axA.text(glc_center, glc_vals.max() + abs(glc_vals.min()) * 0.06, "Glucose (uptake)",
              ha="center", va="bottom", fontsize=10.5, fontweight="bold")
    axA.text(lac_center, top_of_brackets, "Lactate (secretion)",
              ha="center", va="bottom", fontsize=10.5, fontweight="bold")

    y_bottom = glc_vals.min() * 1.12
    y_top = top_of_brackets + abs(glc_vals.min()) * 0.10 + 40
    axA.set_ylim(y_bottom, y_top)
    axA.set_xticks(tick_pos)
    axA.set_xticklabels(tick_lab, fontsize=9.3)
    axA.set_xlim(-0.6, max(tick_pos) + 0.6)
    axA.set_ylabel("measured exchange flux (fmol cell$^{-1}$ h$^{-1}$)\npositive = secretion, negative = uptake")

    # Oxygen-level group labels, one tier below the treatment tick labels,
    # each spanning its own pair of ticks (see OXYGEN_GROUPS above).
    xaxis_trans = axA.get_xaxis_transform()
    y_line, y_tick, y_text = -0.155, -0.135, -0.185
    for met in metabolites:
        for i0, i1, glabel in OXYGEN_GROUPS:
            gx0 = xpos[(met,) + CONDITION_ORDER[i0]]
            gx1 = xpos[(met,) + CONDITION_ORDER[i1]]
            axA.plot([gx0, gx1], [y_line, y_line], transform=xaxis_trans,
                      color="black", linewidth=0.8, clip_on=False)
            axA.plot([gx0, gx0], [y_tick, y_line], transform=xaxis_trans,
                      color="black", linewidth=0.8, clip_on=False)
            axA.plot([gx1, gx1], [y_tick, y_line], transform=xaxis_trans,
                      color="black", linewidth=0.8, clip_on=False)
            axA.text((gx0 + gx1) / 2.0, y_text, glabel, transform=xaxis_trans,
                      ha="center", va="top", fontsize=9.0, clip_on=False)

    # The vertical gap between the glucose and lactate clusters carries no
    # data at any y, so a centered legend cannot occlude any point or bracket.
    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="none",
               markeredgecolor="black", markersize=7, label="DMSO (vehicle)"),
        Line2D([0], [0], marker="^", linestyle="none", markerfacecolor="none",
               markeredgecolor="black", markersize=7, label="BAY (molidustat)"),
        Line2D([0], [0], marker="s", linestyle="none", markerfacecolor="0.25",
               markeredgecolor="0.25", markersize=7, label="21% O$_2$ (filled)"),
        Line2D([0], [0], marker="s", linestyle="none", markerfacecolor="none",
               markeredgecolor="0.25", markersize=7, label="0.5% O$_2$ (open)"),
        Line2D([0], [0], color="black", linewidth=2.6, label="group mean"),
    ]
    leg = axA.legend(handles=legend_handles, loc="center", ncol=1, fontsize=7.4,
                      frameon=True, framealpha=0.92, handletextpad=0.5, borderpad=0.6,
                      labelspacing=0.5)
    panel_letter(axA, "A")

    # ---------------- Panel B: predicted lactate vs measured growth rate ---
    m = pred[(pred.task_mode == "measured") & (pred.exchange_id == "EX_lac__L_e")].copy()
    if m.empty:
        raise ValueError("predictions.tsv: no measured/EX_lac__L_e rows found")
    grp = (m.groupby(["arm", "oxygen", "treatment"])
           .agg(task_flux=("task_flux", "mean"), point=("point", "mean"))
           .reset_index())

    flux_values = np.sort(grp["task_flux"].unique())
    spacing = np.min(np.diff(flux_values)) if len(flux_values) > 1 else 1.0
    n_arms = len(ARMS7)
    for i, arm in enumerate(ARMS7):
        garm = grp[grp.arm == arm].sort_values("task_flux")
        if garm.empty:
            continue
        dodge = (i - (n_arms - 1) / 2.0) * spacing * 0.06
        axB.plot(garm.task_flux + dodge, garm.point, linestyle="-", linewidth=1.1,
                  color=ARM_COLORS7[arm], alpha=0.9, zorder=2)
        axB.scatter(garm.task_flux + dodge, garm.point, marker=ARM_MARKERS7[arm],
                    s=48, facecolor=ARM_COLORS7[arm], edgecolor="black", linewidths=0.5,
                    label=ARM_LABELS7[arm], zorder=3)

    axB.set_xlabel("measured growth rate imposed on the task (h$^{-1}$)")
    axB.set_ylabel("predicted lactate exchange\n(canonical model units; positive = secretion)")
    axB.grid(True, linewidth=0.5, color="0.88", zorder=0)
    axB.set_axisbelow(True)
    # Extra headroom at upper-left for the 7-arm legend: growth rate and
    # predicted lactate move together, so low x / high y is always empty.
    b_ymin, b_ymax = grp["point"].min(), grp["point"].max()
    b_span = b_ymax - b_ymin
    axB.set_ylim(b_ymin - 0.06 * b_span, b_ymax + 0.42 * b_span)
    axB.legend(loc="upper left", ncol=1, fontsize=8.0, frameon=True, framealpha=0.9,
               handletextpad=0.4, borderpad=0.4, labelspacing=0.35)
    panel_letter(axB, "B")

    # ---------------- Panel C: crossed profile x growth-rate control -------
    key = "magnitude_g1.00|EX_lac__L_e"
    if key not in csumm.get("decomposition", {}):
        raise ValueError(f"crossed_summary.json: decomposition key '{key}' not found")
    decomp = csumm["decomposition"][key]

    cx = crossed[(crossed.arm == "magnitude_g1.00") & (crossed.exchange_id == "EX_lac__L_e")]
    if cx.empty:
        raise ValueError("crossed_predictions.tsv: no magnitude_g1.00/EX_lac__L_e rows found")
    n_profiles = cx.context_id.nunique()
    n_growth_rates = cx.task_flux.nunique()
    # alpha is deliberately low but linewidth/markersize deliberately visible:
    # the 16 lines are meant to read as a soft overlapping band here, not a
    # single crisp vector stroke -- the inset added below (after the
    # headroom and annotation box are set up) is what actually resolves the
    # 16 distinct values.
    for ctx, g in cx.groupby("context_id"):
        g = g.sort_values("task_flux")
        axC.plot(g.task_flux, g.point, color="#4477AA", alpha=0.22, linewidth=2.2,
                  marker="o", markersize=3.6, zorder=2)

    axC.set_xlabel("imposed growth rate (h$^{-1}$)")
    axC.set_ylabel("predicted lactate exchange\n(canonical model units; positive = secretion)")
    axC.grid(True, linewidth=0.5, color="0.88", zorder=0)
    axC.set_axisbelow(True)
    # Extra headroom at upper-left for the decomposition annotation box: at
    # the lowest imposed growth rate every profile's prediction is also at
    # its lowest, so low x / high y is empty.
    c_ymin, c_ymax = cx["point"].min(), cx["point"].max()
    c_span = c_ymax - c_ymin
    axC.set_ylim(c_ymin - 0.06 * c_span, c_ymax + 0.55 * c_span)

    spread_growth = decomp["mean_spread_across_growth_rates_at_fixed_transcript"]
    spread_profile = decomp["mean_spread_across_transcripts_at_fixed_growth_rate"]
    ratio = decomp["ratio_transcript_to_task"]
    note = (f"crossed design: {n_profiles} transcript profiles × {n_growth_rates} growth rates\n"
            f"arm γ=1.00, EX_lac__L_e\n"
            f"mean spread across growth rates: {spread_growth:.3f}\n"
            f"mean spread across profiles: {spread_profile:.4f}\n"
            f"ratio (profile / growth-rate): {ratio*100:.1f}%")
    axC.text(0.03, 0.97, note, transform=axC.transAxes, fontsize=8.3,
              ha="left", va="top",
              bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                        edgecolor="0.5", alpha=0.92), zorder=6)

    # ---- Inset: zoom on the single highest growth rate (largest across-
    # profile spread of the four) so the 16 nominally-coincident lines above
    # are resolved into 16 distinct values; the main panel keeps showing the
    # complete rise across all growth rates unchanged.
    flux_sorted = np.sort(cx["task_flux"].unique())
    zoom_flux, prev_flux = flux_sorted[-1], flux_sorted[-2]
    at_zoom = cx.loc[cx["task_flux"] == zoom_flux, "point"]
    zy_lo, zy_hi = float(at_zoom.min()), float(at_zoom.max())
    zy_pad = (zy_hi - zy_lo) * 0.35
    x_span = zoom_flux - prev_flux
    zx_lo, zx_hi = zoom_flux - 0.45 * x_span, zoom_flux + 0.10 * x_span

    rect = Rectangle((zx_lo, zy_lo - zy_pad), zx_hi - zx_lo,
                      (zy_hi + zy_pad) - (zy_lo - zy_pad),
                      fill=False, edgecolor="0.25", linewidth=0.9,
                      linestyle="--", zorder=4)
    axC.add_patch(rect)

    axins = axC.inset_axes([0.55, 0.08, 0.41, 0.40])
    for ctx, g in cx.groupby("context_id"):
        g = g.sort_values("task_flux")
        axins.plot(g.task_flux, g.point, color="#4477AA", alpha=0.6, linewidth=1.3,
                   marker="o", markersize=3.2, zorder=2)
    axins.set_xlim(zx_lo, zx_hi)
    axins.set_ylim(zy_lo - zy_pad, zy_hi + zy_pad)
    axins.set_xticks([])  # x-precision is not the point here; see caption
    axins.yaxis.set_major_locator(MaxNLocator(3))
    axins.ticklabel_format(axis="y", useOffset=False, style="plain")
    axins.tick_params(labelsize=6.2)
    axins.grid(True, linewidth=0.4, color="0.85", zorder=0)
    for spine in axins.spines.values():
        spine.set_edgecolor("0.25")
        spine.set_linewidth(0.9)
    axins.text(0.5, 1.15, f"zoom: all {n_profiles} profiles at {zoom_flux:.4f} h$^{{-1}}$",
               transform=axins.transAxes, ha="center", va="bottom", fontsize=6.4,
               style="italic",
               bbox=dict(boxstyle="round,pad=0.12", facecolor="white",
                         edgecolor="none", alpha=0.85))

    panel_letter(axC, "C")

    fig.tight_layout(w_pad=2.6)

    # Verify the Panel A tick labels actually fit without touching, checked
    # against the FINAL post-tight_layout geometry (matplotlib's own
    # renderer), not assumed from a font-size-to-width conversion.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    label_boxes = [t.get_window_extent(renderer) for t in axA.get_xticklabels()]
    for i in range(len(label_boxes) - 1):
        if label_boxes[i].x1 > label_boxes[i + 1].x0:
            raise ValueError(
                f"Figure 5 Panel A: x-tick labels {tick_lab[i]!r} and "
                f"{tick_lab[i + 1]!r} overlap by "
                f"{label_boxes[i].x1 - label_boxes[i + 1].x0:.1f}px; increase `step`.")

    # Verify the Panel A legend box (placed in the gap between the glucose
    # and lactate clusters) does not visually cover any plotted replicate --
    # checked in real pixel space so a future change in `gap`, the data, or
    # the legend content cannot silently reintroduce a hidden point.
    leg_box = leg.get_window_extent(renderer)
    pt_disp = axA.transData.transform(np.column_stack([point_x, point_y]))
    covered = ((pt_disp[:, 0] >= leg_box.x0) & (pt_disp[:, 0] <= leg_box.x1) &
               (pt_disp[:, 1] >= leg_box.y0) & (pt_disp[:, 1] <= leg_box.y1))
    if covered.any():
        raise ValueError(
            f"Figure 5 Panel A: the legend box covers {int(covered.sum())} plotted "
            f"replicate point(s); increase `gap` or reposition the legend.")

    out_path = os.path.join(OUTDIR, "Figure_5_Independent_Evaluation.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    report["figure5"] = {
        "path": out_path,
        "panelA_separated_contrasts": [
            {
                "metabolite": row["metabolite"],
                "a": f"{row['a_oxygen']} {row['a_treatment']}",
                "b": f"{row['b_oxygen']} {row['b_treatment']}",
                "a_mean": float(row["a_measured_mean"]),
                "b_mean": float(row["b_measured_mean"]),
            }
            for _, row in separated.sort_values("contrast").iterrows()
        ],
        "panelB_n_arms": int(grp.arm.nunique()),
        "panelB_growth_rates_h": [float(v) for v in flux_values],
        "panelB_mean_point_by_arm_condition": {
            f"{r.arm}|{r.oxygen}|{r.treatment}": float(r.point) for r in grp.itertuples()
        },
        "panelC_n_profiles": int(n_profiles),
        "panelC_decomposition": decomp,
    }


# ===========================================================================
def main():
    report = {}
    make_figure4(report)
    make_figure5(report)

    sums_path = os.path.join(OUTDIR, "SHA256SUMS.txt")
    lines = []
    for key in ("figure4", "figure5"):
        p = report[key]["path"]
        digest = sha256_of(p)
        w, h = png_size(p)
        size_bytes = os.path.getsize(p)
        lines.append(f"{digest}  {os.path.basename(p)}\n")
        report[key]["sha256"] = digest
        report[key]["pixel_size"] = [w, h]
        report[key]["file_size_bytes"] = size_bytes
    with open(sums_path, "w") as f:
        f.writelines(lines)

    print(json.dumps(report, indent=2, default=str))


if __name__ == "__main__":
    main()
