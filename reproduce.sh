#!/usr/bin/env bash
# CPU study entry point. Default is read-only preflight, never reserved evaluation.
set -euo pipefail
ROOT="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
PYTHON_BIN="${PYTHON:-$ROOT/.venv/bin/python}"
STAGE=preflight
PARTITION=development
WORKERS=4
MAX_CONTEXTS=""
OUTPUT=""
RESUME=0
SCENARIOS=(primary half_serum lower_task)

usage() {
  cat <<'EOF'
Usage: bash reproduce.sh [--stage STAGE] [--partition development|test]
  --stage preflight|data|analysis|tests|controls|run|audit|evaluate|summary|all|smoke
  --workers N               CPU worker count (default 4)
  --scenarios CSV           default primary,half_serum,lower_task
  --output-dir PATH         default results/development or results/development_smoke
  --max-contexts N           development optimization smoke only; no cohort evaluation
  --resume                  reuse only hash-verified successful arm artifacts
Environment: PYTHON=/absolute/path/to/python (default study/.venv/bin/python)
No downloads, environment installation, protocol authorization or publication occurs.
EOF
}
fail() { printf '%s\n' "Error: $*" >&2; exit 2; }
need_value() { [[ $# -ge 2 && -n "$2" ]] || fail "Missing value for $1"; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    --stage) need_value "$@"; STAGE="$2"; shift 2;;
    --partition) need_value "$@"; PARTITION="$2"; shift 2;;
    --workers) need_value "$@"; WORKERS="$2"; shift 2;;
    --max-contexts) need_value "$@"; MAX_CONTEXTS="$2"; shift 2;;
    --output-dir) need_value "$@"; OUTPUT="$2"; shift 2;;
    --scenarios) need_value "$@"; IFS=',' read -r -a SCENARIOS <<< "$2"; shift 2;;
    --resume) RESUME=1; shift;;
    --help|-h) usage; exit 0;;
    *) fail "Unknown argument: $1";;
  esac
done
case "$STAGE" in preflight|data|analysis|tests|controls|run|audit|evaluate|summary|all|smoke) ;; *) fail "Unknown stage: $STAGE";; esac
case "$PARTITION" in development|test) ;; *) fail "Invalid partition: $PARTITION";; esac
[[ "$WORKERS" =~ ^[1-9][0-9]*$ ]] || fail "--workers must be positive"
if [[ "$STAGE" == smoke ]]; then STAGE=all; MAX_CONTEXTS=1; fi
if [[ -n "$MAX_CONTEXTS" ]]; then
  [[ "$MAX_CONTEXTS" =~ ^[1-9][0-9]*$ ]] || fail "--max-contexts must be positive"
  [[ "$PARTITION" == development ]] || fail "Subset smoke is permitted only for development"
  [[ "$STAGE" != evaluate ]] || fail "A subset smoke cannot evaluate the full cohort"
fi
[[ ${#SCENARIOS[@]} -gt 0 ]] || fail "At least one scenario is required"
for scenario in "${SCENARIOS[@]}"; do
  case "$scenario" in primary|half_serum|lower_task) ;; *) fail "Unknown scenario: $scenario";; esac
done
if [[ -z "$OUTPUT" ]]; then
  OUTPUT="$ROOT/results/$PARTITION"
  [[ -z "$MAX_CONTEXTS" ]] || OUTPUT="${OUTPUT}_smoke"
fi
[[ "$OUTPUT" == /* ]] || OUTPUT="$ROOT/$OUTPUT"
[[ -x "$PYTHON_BIN" ]] || fail "Python is not executable: $PYTHON_BIN. Create .venv and install environment.lock, or set PYTHON to an absolute executable path."
export PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONUNBUFFERED=1
export RUN_ID="${RUN_ID:-cpu-$(date -u +%Y%m%dT%H%M%SZ)}"
verify() { "$PYTHON_BIN" -m src.verify_inputs --root "$ROOT" --check "$@"; }

# All execution stages start with exact raw inputs and pinned dependencies.
verify preflight
if [[ "$PARTITION" == test ]]; then
  verify lock
  case "$STAGE" in data|analysis) fail "Reserved mode cannot rebuild locked data or analysis";; esac
fi
[[ "$STAGE" != preflight ]] || exit 0
mkdir -p "$ROOT/logs" "$OUTPUT"
LOG="$ROOT/logs/${RUN_ID}_${PARTITION}_${STAGE}.log"
exec > >(tee -a "$LOG") 2>&1
printf 'Stage=%s partition=%s Python=%s output=%s\n' "$STAGE" "$PARTITION" "$PYTHON_BIN" "$OUTPUT"

run_data() {
  if [[ -e "$ROOT/protocol/PROTOCOL_LOCK.json" || -e "$ROOT/protocol/PROTOCOL_LOCK.md" ]]; then
    [[ "$STAGE" == all ]] || fail "Data rebuild is forbidden after the protocol lock"
    verify prepared
  else
    "$PYTHON_BIN" -m src.prepare_data --root "$ROOT"
    verify seal-data
  fi
}
run_analysis() {
  verify prepared
  if [[ -e "$ROOT/protocol/PROTOCOL_LOCK.json" || -e "$ROOT/protocol/PROTOCOL_LOCK.md" ]]; then
    [[ "$STAGE" == all ]] || fail "Analysis selection is forbidden after the protocol lock"
    verify analysis
  else
    "$PYTHON_BIN" -m src.prepare_analysis --root "$ROOT"
    verify analysis
  fi
}
run_tests() { "$PYTHON_BIN" -m pytest -q tests; }
run_controls() {
  verify analysis
  [[ -f "$ROOT/src/production_controls.py" ]] || fail "Production controls are not yet implemented; this stage cannot be certified"
  "$PYTHON_BIN" -m src.production_controls --root "$ROOT"
  verify controls
}
run_optimization() {
  verify analysis
  verify controls
  [[ "$PARTITION" != test ]] || verify lock
  if [[ -d "$OUTPUT/arms" ]]; then
    [[ "$RESUME" == 1 ]] || fail "Output contains arm artifacts; choose a new --output-dir or explicitly --resume"
  fi
  if [[ "$RESUME" == 1 ]]; then verify cache --output "$OUTPUT"; fi
  local command=("$PYTHON_BIN" -m src.run_study --partition "$PARTITION" --scenarios "${SCENARIOS[@]}" \
    --workers "$WORKERS" --output "$OUTPUT")
  [[ -z "$MAX_CONTEXTS" ]] || command+=(--max-contexts "$MAX_CONTEXTS")
  "${command[@]}"
  verify complete --output "$OUTPUT"
  run_audit
}
run_audit() {
  verify complete --output "$OUTPUT"
  local temporary="$OUTPUT/independent_audit.json.tmp"
  if "$PYTHON_BIN" -m src.audit_results --root "$ROOT" --run-directory "$OUTPUT" > "$temporary"; then
    mv "$temporary" "$OUTPUT/independent_audit.json"
    verify audit --output "$OUTPUT"
  else
    mv "$temporary" "$OUTPUT/independent_audit.json"
    fail "Independent saved-vector audit failed; inspect $OUTPUT/independent_audit.json"
  fi
}
run_evaluation() {
  if [[ -n "$MAX_CONTEXTS" ]]; then
    printf '%s\n' 'Subset smoke completed; cohort evaluation is intentionally absent.'
    return
  fi
  verify analysis
  verify complete --output "$OUTPUT"
  verify audit --output "$OUTPUT"
  if [[ "$PARTITION" == test ]]; then
    verify lock
  fi
  for scenario in "${SCENARIOS[@]}"; do
    local command=("$PYTHON_BIN" -m src.evaluation --root "$ROOT" --partition "$PARTITION" \
      --predictions "$OUTPUT/predictions.tsv" --scenario "$scenario" \
      --output "$OUTPUT/evaluation_${scenario}.json")
    [[ "$PARTITION" != test ]] || command+=(--lock-file "$ROOT/protocol/PROTOCOL_LOCK.json")
    "${command[@]}"
    "$PYTHON_BIN" -m src.transform_sensitivity --root "$ROOT" --partition "$PARTITION" \
      --predictions "$OUTPUT/predictions.tsv" --scenario "$scenario" \
      --output "$OUTPUT/transform_${scenario}.json"
  done
}
run_summary() {
  local arguments=(summary --output "$OUTPUT" --partition "$PARTITION" --scenarios "${SCENARIOS[@]}")
  [[ -z "$MAX_CONTEXTS" ]] || arguments+=(--smoke)
  verify "${arguments[@]}"
}
if [[ "$STAGE" == all ]]; then
  run_data; run_analysis; run_tests; run_controls; run_optimization; run_evaluation; run_summary
else
  case "$STAGE" in
    data) run_data;; analysis) run_analysis;; tests) run_tests;; controls) run_controls;;
    run) run_optimization;; audit) run_audit;; evaluate) run_evaluation;; summary) run_summary;;
  esac
fi
printf 'Completed requested stage %s. Log: %s\n' "$STAGE" "$LOG"
