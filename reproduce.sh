#!/usr/bin/env bash
# Stages of the study, in order. Each stage writes into its own output directory and every
# stage after the first reads only what an earlier stage wrote.
#
#   bash reproduce.sh [tests|scores|study|evaluate|coverage|mechanism|closure|independent|figures|verify]
#
# Three third-party inputs are not redistributed here. SOURCES.md gives the retrieval route and
# the expected SHA-256 for each; place them where the paths below expect them before running.
# The solver is GLPK and no accelerator is required at any stage.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
export PYTHONPATH="$HERE"
STAGE="${1:-all}"
ROOT="${ROOT:-$HERE/work}"
say() { printf '\n== %s ==\n' "$1"; }

need() {
  [ -s "$1" ] || { echo "MISSING REQUIRED INPUT: $1" >&2; echo "  $2" >&2; exit 2; }
}

# ----------------------------------------------------------------- tests
if [ "$STAGE" = all ] || [ "$STAGE" = tests ]; then
say "known-answer checks, run before anything is interpreted"
python3 -m unittest discover -s mechanism -p 'test_*.py'
python3 -m unittest discover -s copeland  -p 'test_*.py'
python3 -m unittest discover -s coverage  -p 'test_*.py'
fi

# ----------------------------------------------------------------- inputs and scores
if [ "$STAGE" = all ] || [ "$STAGE" = scores ]; then
say "source verification and reaction evidence"
need "$ROOT/data/raw/Recon3D.json"        "see SOURCES.md"
need "$ROOT/data/raw/cellminer.zip"       "see SOURCES.md"
need "$ROOT/data/raw/core/NIHMS419088-supplement-Database_S1.xls" "see SOURCES.md"
python3 src/verify_inputs.py --root "$ROOT"
python3 src/prepare_data.py --root "$ROOT"
python3 src/prepare_analysis.py --root "$ROOT"
fi

# ----------------------------------------------------------------- the encoding study
if [ "$STAGE" = all ] || [ "$STAGE" = study ]; then
say "development and reserved arms"
python3 src/run_study.py --root "$ROOT" --partition development
python3 src/production_controls.py --root "$ROOT"
python3 src/run_study.py --root "$ROOT" --partition test
fi

if [ "$STAGE" = all ] || [ "$STAGE" = evaluate ]; then
say "concordance, contrasts and sensitivity"
python3 src/evaluation.py --root "$ROOT"
python3 src/transform_sensitivity.py --root "$ROOT"
python3 src/audit_results.py --root "$ROOT"
fi

# ----------------------------------------------------------------- diagnosis
if [ "$STAGE" = all ] || [ "$STAGE" = coverage ]; then
say "encoding-unanimity coverage"
python3 coverage/agreement_audit.py --root "$ROOT"
fi

if [ "$STAGE" = all ] || [ "$STAGE" = mechanism ]; then
say "constraint-stage attribution, interval position and cross-encoding cost"
python3 mechanism/reserved_range_stage.py --root "$ROOT" --output "$ROOT/reserved_range"
python3 mechanism/range_location.py
python3 mechanism/cross_cost.py
fi

# ----------------------------------------------------------------- closure analyses
# The tolerance sweeps, the supporting statistics, the two comparator formulations, the
# second-solver check and the sensitivity-scenario ladders. All read saved outputs of the stages
# above; the RIPTiDe run needs the riptide package and the solver check needs highspy.
if [ "$STAGE" = all ] || [ "$STAGE" = closure ]; then
say "optimality-tolerance sweeps and interval scoring"
python3 closure/cap_sweep.py --root "$ROOT" --out "$ROOT/closure/cap_sweep"
python3 closure/cap_sweep.py --root "$ROOT" --out "$ROOT/closure/cap_sweep_fine" --fine
python3 closure/cap_concordance.py --root "$ROOT" --sweep "$ROOT/closure/cap_sweep" --out "$ROOT/closure/cap_sweep"
CAP_NAMES=rel_1e-6,rel_1e-5,rel_1e-4,rel_1e-3 python3 closure/cap_concordance.py --root "$ROOT" \
        --sweep "$ROOT/closure/cap_sweep_fine" --out "$ROOT/closure/cap_sweep_fine"
python3 closure/cap_sweep_two_stage.py --root "$ROOT" --out "$ROOT/closure/two_stage"
say "bound-based interface, its fixed-task control, and RIPTiDe"
python3 closure/eflux_run.py --root "$ROOT" --out "$ROOT/closure/bound_B10" --B 10
python3 closure/eflux_run.py --root "$ROOT" --out "$ROOT/closure/bound_B3"  --B 3
# the fixed-task controls hold the growth demand at 0.9 of the minimum profile maximum of each
# mapping, the values recorded in the run manifests under results/closure/
python3 closure/eflux_run.py --root "$ROOT" --out "$ROOT/closure/bound_B10_fixedtask_identity" --B 10 \
        --arms magnitude_g1.00 --task-fixed 0.998189
python3 closure/eflux_run.py --root "$ROOT" --out "$ROOT/closure/bound_B10_fixedtask_ordinal"  --B 10 \
        --arms ordinal --task-fixed 0.653891
python3 closure/riptide_run.py --root "$ROOT" --out "$ROOT/closure/riptide"
# range_location_reserved_stages.tsv is written by mechanism/range_location.py into the working directory
for pair in "reserved_v3:study_cost" "closure/bound_B10:bound_B10" "closure/bound_B3:bound_B3" \
            "closure/bound_B10_fixedtask_identity:bound_B10_fixed_identity" \
            "closure/bound_B10_fixedtask_ordinal:bound_B10_fixed_ordinal" "closure/riptide:riptide"; do
run="${pair%%:*}"; label="${pair##*:}"
python3 closure/comparator_eval.py --root "$ROOT" --pred "$ROOT/results/$run/predictions.tsv" \
        --b0 "$HERE/range_location_reserved_stages.tsv" --out "$ROOT/closure/comparators" --label "$label"
done
say "supporting statistical analyses"
python3 closure/stat_repairs.py --root "$ROOT" --out "$ROOT/closure/stat_repairs"
python3 closure/stat_repairs2.py --root "$ROOT" --ledger "$ROOT/reserved_range/range_stage_ledger.tsv" \
        --out "$ROOT/closure/stat_repairs2" --comparators riptide="$ROOT/closure/riptide/predictions.tsv"
say "second solver and sensitivity-scenario ladders"
python3 closure/solver_check.py --root "$ROOT" --out "$ROOT/closure/solver_check"
for scen in half_serum lower_task; do
python3 closure/scenario_ladder.py --root "$ROOT" --output "$ROOT/closure/ladders/$scen" --scenario "$scen"
python3 closure/ladder_summary.py --ledger "$ROOT/closure/ladders/$scen/range_stage_ledger.tsv" \
        --out "$ROOT/closure/ladders/ladder_summary_$scen.json" --label "$scen"
done
fi

# ----------------------------------------------------------------- independent evaluation
if [ "$STAGE" = all ] || [ "$STAGE" = independent ]; then
say "independent condition-response evaluation"
python3 copeland/copeland_prepare.py  --study-root "$ROOT" \
        --compendium "$ROOT/copeland/compendium" --rnaseq "$ROOT/copeland/rnaseq" \
        --out "$ROOT/copeland/prepared"
python3 copeland/copeland_run.py      --study-root "$ROOT" --prepared "$ROOT/copeland/prepared" \
        --out "$ROOT/copeland/run"
python3 copeland/copeland_crossed.py  --study-root "$ROOT" --prepared "$ROOT/copeland/prepared" \
        --out "$ROOT/copeland/crossed"
python3 copeland/copeland_evaluate.py --run "$ROOT/copeland/run" \
        --compendium "$ROOT/copeland/compendium" --out "$ROOT/copeland/evaluation"
python3 copeland/copeland_meantask.py --study-root "$ROOT" --prepared "$ROOT/copeland/prepared" \
        --out "$ROOT/copeland/run_meantask"
python3 copeland/copeland_evaluate.py --run "$ROOT/copeland/run_meantask" \
        --compendium "$ROOT/copeland/compendium" --out "$ROOT/copeland/evaluation_meantask"
fi

# ----------------------------------------------------------------- outputs
if [ "$STAGE" = all ] || [ "$STAGE" = figures ]; then
say "figures"
python3 figures/make_figures.py
fi

if [ "$STAGE" = all ] || [ "$STAGE" = verify ]; then
say "trace every reported number back to its artefact"
python3 verify/verify_manuscript_numbers.py
fi

say "done"
