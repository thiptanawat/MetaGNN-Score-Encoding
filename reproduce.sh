#!/usr/bin/env bash
# Stages of the study, in order. Each stage writes into its own output directory and every
# stage after the first reads only what an earlier stage wrote.
#
#   bash reproduce.sh [tests|scores|study|evaluate|coverage|mechanism|independent|figures|verify]
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
