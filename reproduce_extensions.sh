#!/usr/bin/env bash
# The analyses added after the protocol lock, in dependency order. Every stage reads the locked
# study artefacts in this repository (which is laid out as the study root) and the outputs of the
# stages before it, and writes into its own directory under results/.
#
#   bash reproduce_extensions.sh STAGE [--workers N] [--out DIR]
#
#   STAGE   tests | coverage | mechanism | closure | independent | endpoint | figures | verify | all
#
# The locked study itself (data preparation, panel selection, controls, the encoding arms and their
# evaluation) is driven by reproduce.sh, which is the byte-identical driver named in
# protocol/PROTOCOL_LOCK.json; run it first if the per-arm records are to be regenerated rather
# than fetched (verify/fetch_release_assets.sh).
#
# Deposited outputs are never overwritten: by default every stage writes into a fresh directory
# named by --out (default: results_repeat/<timestamp>), and the deposited results/ tree stays as it
# is for comparison. Set OUT=results to write next to the deposited outputs; the scripts that keep
# run manifests refuse to overwrite an existing run.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
export PYTHONPATH="$HERE${PYTHONPATH:+:$PYTHONPATH}"
export STUDY_ROOT="$HERE"
PYTHON_BIN="${PYTHON:-$HERE/.venv/bin/python}"
STAGE="${1:-}"; shift || true
WORKERS=4
OUT="${OUT:-$HERE/results_repeat/$(date -u +%Y%m%dT%H%M%SZ)}"
while [ $# -gt 0 ]; do
  case "$1" in
    --workers) WORKERS="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    *) echo "unknown argument: $1" >&2; exit 2;;
  esac
done
case "$STAGE" in tests|coverage|mechanism|closure|independent|endpoint|figures|verify|all) ;;
  *) sed -n '2,20p' "$0"; exit 2;; esac
[ -x "$PYTHON_BIN" ] || { echo "Python is not executable: $PYTHON_BIN (create .venv from environment.lock or set PYTHON)" >&2; exit 2; }
case "$OUT" in /*) ;; *) OUT="$HERE/$OUT";; esac
mkdir -p "$OUT"
say() { printf '\n== %s ==\n' "$1"; }
need() { [ -s "$1" ] || { echo "MISSING REQUIRED INPUT: $1" >&2; echo "  $2" >&2; exit 2; }; }
ARMS="$HERE/results/reserved_v3/arms"
DEV_ARMS="$HERE/results/development_v3/arms"
DEV_SUMMARY="$HERE/results/development_v3/run_summary.json"
RL_STAGES="$HERE/results/mechanism/range_location_reserved_stages.tsv"   # written by mechanism/range_location.py

# ----------------------------------------------------------------- tests
if [ "$STAGE" = all ] || [ "$STAGE" = tests ]; then
say "known-answer checks of the post-lock analyses (the locked study's own suite is reproduce.sh --stage tests)"
"$PYTHON_BIN" -m unittest discover -s mechanism -p 'test_*.py'
"$PYTHON_BIN" -m unittest discover -s copeland  -p 'test_*.py'
"$PYTHON_BIN" -m unittest discover -s coverage  -p 'test_*.py'
"$PYTHON_BIN" -m pytest -q reports/validation_scripts/test_cross_encoding_ranges.py
fi

# ----------------------------------------------------------------- encoding-unanimity coverage
if [ "$STAGE" = all ] || [ "$STAGE" = coverage ]; then
say "encoding-unanimity coverage over the reserved predictions"
"$PYTHON_BIN" coverage/agreement_audit.py --root "$HERE" --output "$OUT/coverage"
fi

# ----------------------------------------------------------------- constraint-stage attribution
if [ "$STAGE" = all ] || [ "$STAGE" = mechanism ]; then
say "nested admissible sets: development pilot, reserved cohort, uniform-cost control"
need "$DEV_SUMMARY" "results/development_v3/run_summary.json: fetch development_v3_arms.tar.gz (verify/fetch_release_assets.sh) or rerun reproduce.sh --stage all --partition development"
need "$ARMS/12b0ae9b6ce61e27.json" "results/reserved_v3/arms: fetch reserved_v3_arms.tar.gz (verify/fetch_release_assets.sh) or rerun reproduce.sh --stage all --partition test"
"$PYTHON_BIN" mechanism/nested_range_pilot.py    --root "$HERE" --output "$OUT/mechanism/bottleneck_pilot"
"$PYTHON_BIN" mechanism/readback_check.py "$OUT/mechanism/bottleneck_pilot"
"$PYTHON_BIN" mechanism/reserved_range_stage.py --root "$HERE" --output "$OUT/mechanism/reserved_range" --workers "$WORKERS"
# the primary-scenario uniform-cost record is context independent; any reserved run holds exactly one
UNIFORM_JSON="$("$PYTHON_BIN" - "$ARMS" <<'PY'
import json, sys, glob
for f in sorted(glob.glob(sys.argv[1] + "/*.json")):
    d = json.load(open(f))
    if d.get("arm") == "uniform_pfba" and d.get("scenario") == "primary":
        print(f); break
PY
)"
"$PYTHON_BIN" mechanism/uniform_control.py --root "$HERE" --output "$OUT/mechanism/uniform_control" --uniform-json "$UNIFORM_JSON"
say "interval position and cross-encoding cost comparison (saved vectors, no LP)"
mkdir -p "$OUT/mechanism"
# both scripts write into the working directory and read the ledgers written above
( cd "$OUT/mechanism" && PILOT_LEDGER="$OUT/mechanism/bottleneck_pilot/range_stage_ledger.tsv" \
      RESERVED_LEDGER="$OUT/mechanism/reserved_range/range_stage_ledger.tsv" "$PYTHON_BIN" "$HERE/mechanism/range_location.py" )
( cd "$OUT/mechanism" && "$PYTHON_BIN" "$HERE/mechanism/cross_cost.py" )
fi

# ----------------------------------------------------------------- objective-allowance sweeps, comparators, statistics
if [ "$STAGE" = all ] || [ "$STAGE" = closure ]; then
say "objective-allowance sweeps (single cap, fine grid, complete two-cap readout) and their interval scoring"
need "$ARMS/12b0ae9b6ce61e27.json" "results/reserved_v3/arms: fetch reserved_v3_arms.tar.gz (verify/fetch_release_assets.sh)"
"$PYTHON_BIN" closure/cap_sweep.py --root "$HERE" --out "$OUT/closure/cap_sweep"      --workers "$WORKERS"
"$PYTHON_BIN" closure/cap_sweep.py --root "$HERE" --out "$OUT/closure/cap_sweep_fine" --workers "$WORKERS" --fine
"$PYTHON_BIN" closure/cap_concordance.py --root "$HERE" --sweep "$OUT/closure/cap_sweep" --out "$OUT/closure/cap_sweep"
CAP_NAMES=rel_1e-6,rel_1e-5,rel_1e-4,rel_1e-3 "$PYTHON_BIN" closure/cap_concordance.py --root "$HERE" \
        --sweep "$OUT/closure/cap_sweep_fine" --out "$OUT/closure/cap_sweep_fine"
"$PYTHON_BIN" closure/cap_sweep_two_stage.py --root "$HERE" --out "$OUT/closure/two_stage" --workers "$WORKERS"
say "bound-based interface, its common-task controls, and RIPTiDe (needs the riptide package)"
"$PYTHON_BIN" closure/eflux_run.py --root "$HERE" --out "$OUT/closure/bound_B10" --B 10 --workers "$WORKERS"
"$PYTHON_BIN" closure/eflux_run.py --root "$HERE" --out "$OUT/closure/bound_B3"  --B 3  --workers "$WORKERS"
# the common-task controls hold the growth demand at 0.9 of the smallest profile maximum of each
# mapping, the values recorded in results/closure/bound_B10_fixedtask_*/run_manifest.json
"$PYTHON_BIN" closure/eflux_run.py --root "$HERE" --out "$OUT/closure/bound_B10_fixedtask_identity" --B 10 \
        --arms magnitude_g1.00 --task-fixed 0.998189 --workers "$WORKERS"
"$PYTHON_BIN" closure/eflux_run.py --root "$HERE" --out "$OUT/closure/bound_B10_fixedtask_ordinal"  --B 10 \
        --arms ordinal --task-fixed 0.653891 --workers "$WORKERS"
"$PYTHON_BIN" closure/riptide_run.py --root "$HERE" --out "$OUT/closure/riptide" --samples 500 --fraction 0.8 --workers "$WORKERS"
say "comparator scoring against the shared network-only ranges"
B0="$RL_STAGES"; [ -s "$OUT/mechanism/range_location_reserved_stages.tsv" ] && B0="$OUT/mechanism/range_location_reserved_stages.tsv"
for pair in "results/reserved_v3:study_cost" "$OUT/closure/bound_B10:bound_B10" "$OUT/closure/bound_B3:bound_B3" \
            "$OUT/closure/bound_B10_fixedtask_identity:bound_B10_fixed_identity" \
            "$OUT/closure/bound_B10_fixedtask_ordinal:bound_B10_fixed_ordinal" "$OUT/closure/riptide:riptide" \
            "$OUT/closure/two_stage/both_0.01:twostage_both_0.01" "$OUT/closure/two_stage/both_0.2:twostage_both_0.2" \
            "$OUT/closure/two_stage/cost_0.01_parsimony_exact:twostage_cost_0.01_parsimony_exact" \
            "$OUT/closure/two_stage/cost_0.2_parsimony_exact:twostage_cost_0.2_parsimony_exact"; do
  run="${pair%%:*}"; label="${pair##*:}"
  case "$run" in /*) pred="$run/predictions.tsv";; *) pred="$HERE/$run/predictions.tsv";; esac
  "$PYTHON_BIN" closure/comparator_eval.py --root "$HERE" --pred "$pred" --b0 "$B0" --out "$OUT/closure/comparators" --label "$label"
done
say "supporting statistical analyses (saved outputs only, no LP)"
LEDGER="$HERE/results/mechanism/reserved_range/range_stage_ledger.tsv"
[ -s "$OUT/mechanism/reserved_range/range_stage_ledger.tsv" ] && LEDGER="$OUT/mechanism/reserved_range/range_stage_ledger.tsv"
"$PYTHON_BIN" closure/stat_repairs.py  --root "$HERE" --out "$OUT/closure/stat_repairs"
"$PYTHON_BIN" closure/stat_repairs2.py --root "$HERE" --ledger "$LEDGER" --out "$OUT/closure/stat_repairs2" \
        --comparators riptide="$OUT/closure/riptide/predictions.tsv"
"$PYTHON_BIN" closure/stat_repairs3.py --root "$HERE" --out "$OUT/closure/stat_repairs3"
say "second solver (needs highspy) over the production selection list, then over the 52-target list"
"$PYTHON_BIN" closure/solver_check.py --root "$HERE" --out "$OUT/closure/solver_check_96" --solver hybrid --panel chemistry --workers "$WORKERS"
"$PYTHON_BIN" closure/solver_check.py --root "$HERE" --out "$OUT/closure/solver_check"    --solver hybrid --panel primary   --workers "$WORKERS"
say "sensitivity-scenario ladders"
for scen in half_serum lower_task; do
  "$PYTHON_BIN" closure/scenario_ladder.py --root "$HERE" --output "$OUT/closure/ladders/$scen" --scenario "$scen" --workers "$WORKERS"
  "$PYTHON_BIN" closure/ladder_summary.py --ledger "$OUT/closure/ladders/$scen/range_stage_ledger.tsv" \
        --out "$OUT/closure/ladders/ladder_summary_$scen.json" --label "$scen"
done
"$PYTHON_BIN" closure/ladder_summary.py --ledger "$LEDGER" --out "$OUT/closure/ladders/ladder_summary_primary_selfcheck.json" --label primary
fi

# ----------------------------------------------------------------- independent condition-response evaluation
if [ "$STAGE" = all ] || [ "$STAGE" = independent ]; then
say "independent condition-response evaluation (needs the rdata package and data/raw/Recon3D.json)"
need "$HERE/data/raw/Recon3D.json" "see data/raw/README.md"
SRC="$HERE/work/copeland"
[ -d "$SRC/compendium/.git" ] && [ -d "$SRC/rnaseq/.git" ] || bash copeland/fetch_inputs.sh "$SRC"
"$PYTHON_BIN" copeland/copeland_prepare.py  --study-root "$HERE" --compendium "$SRC/compendium" --rnaseq "$SRC/rnaseq" \
        --out "$OUT/copeland/prepared"
"$PYTHON_BIN" copeland/copeland_run.py      --study-root "$HERE" --prepared "$OUT/copeland/prepared" --out "$OUT/copeland/run" --workers "$WORKERS"
"$PYTHON_BIN" copeland/copeland_crossed.py  --study-root "$HERE" --prepared "$OUT/copeland/prepared" --out "$OUT/copeland/crossed" --workers "$WORKERS"
"$PYTHON_BIN" copeland/copeland_evaluate.py --run "$OUT/copeland/run" --compendium "$SRC/compendium" --out "$OUT/copeland/evaluation"
"$PYTHON_BIN" copeland/copeland_meantask.py --study-root "$HERE" --prepared "$OUT/copeland/prepared" --out "$OUT/copeland/run_meantask" --workers "$WORKERS"
"$PYTHON_BIN" copeland/copeland_evaluate.py --run "$OUT/copeland/run_meantask" --compendium "$SRC/compendium" --out "$OUT/copeland/evaluation_meantask"
fi

# ----------------------------------------------------------------- endpoint reconstruction attempt
if [ "$STAGE" = all ] || [ "$STAGE" = endpoint ]; then
say "reconstruction attempt of the published depletion screen (reads the raw workbook; reported as unresolved)"
need "$HERE/data/raw/core/NIHMS419088-supplement-Database_S1.xls" "see data/raw/README.md"
ENDPOINT_OUT="$OUT/endpoint" "$PYTHON_BIN" endpoint/run_endpoint_audit.py
fi

# ----------------------------------------------------------------- figures
if [ "$STAGE" = all ] || [ "$STAGE" = figures ]; then
say "figures 1 to 3 from the reserved evaluation, figures 4 to 6 from the post-lock analyses"
"$PYTHON_BIN" -m src.make_figures --root "$HERE" --results-dir results/reserved_v3 --output-dir "$OUT/figures/reserved_v3"
FIGURE_OUT="$OUT/figures/manuscript" "$PYTHON_BIN" figures/make_figures.py
fi

# ----------------------------------------------------------------- verification of the deposited outputs
if [ "$STAGE" = all ] || [ "$STAGE" = verify ]; then
say "release integrity, lock reconciliation and the number checker over the deposited outputs"
"$PYTHON_BIN" verify/check_release.py --root "$HERE"
"$PYTHON_BIN" verify/verify_manuscript_numbers.py --root "$HERE"
fi

say "done ($STAGE); outputs under $OUT"
