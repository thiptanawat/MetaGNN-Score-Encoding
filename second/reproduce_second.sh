#!/usr/bin/env bash
# The analysis on a second reconstruction (Human-GEM 2.1.1), stage by stage, in the order of the plan
# (second/SECOND_RECONSTRUCTION_PLAN_2026-11-14.md). Each stage prints what it does. Nothing is deleted
# or overwritten: a stage stops if the output it would write already exists.
#
#   bash second/reproduce_second.sh STAGE [--python PATH] [--humangem DIR] [--root DIR] [--main DIR]
#                                         [--workers N] [--build-outputs DIR] [--resume]
#
#   STAGE   build        write the second root from the main repository and the Human-GEM checkout
#           development  run the 11 development profiles (6 encodings plus the uniform cost) in the second root
#           pilot        the development pilot of the range stages (needs the development arms). Its ledger
#                        is read by the interval-position step of evaluate and by nothing else. It stops at
#                        the time budget fixed in mechanism/nested_range_pilot.py (480 seconds), so the number
#                        of calculations it completes depends on the machine (8 of 22 in the deposited run)
#           lock         write the lock of the second root from the plan; do this before the reserved run
#           reserved     run the 47 reserved profiles (needs the lock)
#           ranges       range stages, cap sweep and solver check (needs the reserved arms)
#           evaluate     the only stage that reads measurements: evaluation, transform sensitivity,
#                        summary of constant targets, interval position, cap concordance (needs the
#                        reserved arms, the ranges and the pilot ledger)
#
# Settings, from the environment or the options (an option wins):
#   PYTHON      interpreter with the packages of environment.lock, plus highspy and osqp
#               (environment_extensions.txt); default MAIN/.venv/bin/python
#   HUMANGEM    checkout of SysBioChalmers/Human-GEM at tag v2.1.1 (commit 6812be2a...); needed by build
#   ROOT        the second root; default MAIN/work/second_root (work/ is not part of the deposited tree)
#   MAIN        the main study root, here the repository; default: the parent of this script's directory
#   WORKERS     CPU workers for the arm runs and range stages; default 4
#   BUILD_OUT   where build writes exchange_map.tsv, target_map.tsv, medium_additions.tsv and
#               build_report.json; default ROOT.build_outputs (the copies in second/ are not touched)
#   RESUME=1    allow development and reserved to continue in an output directory that holds arms
set -euo pipefail
HERE="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
STAGE="${1:-}"; [ $# -eq 0 ] || shift
MAIN="${MAIN:-$(dirname "$HERE")}"
PYTHON_BIN="${PYTHON:-}"
HUMANGEM="${HUMANGEM:-}"
ROOT="${ROOT:-}"
WORKERS="${WORKERS:-4}"
BUILD_OUT="${BUILD_OUT:-}"
RESUME="${RESUME:-0}"
PLAN="$HERE/SECOND_RECONSTRUCTION_PLAN_2026-11-14.md"
PLAN_SHA256="0de8cbf4a737829b4c005b33122568d1abb342166162917f5d5481d4b7cd9a81"
HUMANGEM_COMMIT="6812be2a5c1504a89ea9e84d770786a84f5f3ce2"

usage() { sed -n '2,31p' "$0"; }
fail() { printf 'Error: %s\n' "$*" >&2; exit 2; }
say() { printf '\n== %s ==\n' "$*"; }
while [ $# -gt 0 ]; do
  case "$1" in
    --python) [ $# -ge 2 ] || fail "missing value for $1"; PYTHON_BIN="$2"; shift 2;;
    --humangem) [ $# -ge 2 ] || fail "missing value for $1"; HUMANGEM="$2"; shift 2;;
    --root) [ $# -ge 2 ] || fail "missing value for $1"; ROOT="$2"; shift 2;;
    --main) [ $# -ge 2 ] || fail "missing value for $1"; MAIN="$2"; shift 2;;
    --workers) [ $# -ge 2 ] || fail "missing value for $1"; WORKERS="$2"; shift 2;;
    --build-outputs) [ $# -ge 2 ] || fail "missing value for $1"; BUILD_OUT="$2"; shift 2;;
    --resume) RESUME=1; shift;;
    --help|-h) usage; exit 0;;
    *) fail "unknown argument: $1";;
  esac
done
case "$STAGE" in build|development|lock|reserved|ranges|pilot|evaluate) ;; *) usage; exit 2;; esac
[[ "$WORKERS" =~ ^[1-9][0-9]*$ ]] || fail "--workers must be a positive integer"
MAIN="$(CDPATH= cd -- "$MAIN" && pwd)"
[ -n "$PYTHON_BIN" ] || PYTHON_BIN="$MAIN/.venv/bin/python"
[ -x "$PYTHON_BIN" ] || fail "Python is not executable: $PYTHON_BIN (set PYTHON or --python to an absolute path)"
[ -n "$ROOT" ] || ROOT="$MAIN/work/second_root"
case "$ROOT" in /*) ;; *) ROOT="$PWD/$ROOT";; esac
[ -n "$BUILD_OUT" ] || BUILD_OUT="$ROOT.build_outputs"
case "$BUILD_OUT" in /*) ;; *) BUILD_OUT="$PWD/$BUILD_OUT";; esac

sha() { "$PYTHON_BIN" -c 'import hashlib,sys;print(hashlib.sha256(open(sys.argv[1],"rb").read()).hexdigest())' "$1"; }
check_plan() {
  [ -f "$PLAN" ] || fail "plan not found: $PLAN"
  [ "$(sha "$PLAN")" = "$PLAN_SHA256" ] || fail "the plan differs from the fixed one (expected SHA-256 $PLAN_SHA256)"
  echo "plan SHA-256 verified: $PLAN_SHA256"
}
need_modules() { "$PYTHON_BIN" -c "import importlib.util,sys
missing=[m for m in sys.argv[1:] if importlib.util.find_spec(m) is None]
sys.exit('missing Python packages: '+', '.join(missing)) if missing else None" "$@" || fail "install them as environment.lock and environment_extensions.txt describe"; }
need_root() { [ -f "$ROOT/data/raw/Recon3D.json" ] && [ -f "$ROOT/src/run_study.py" ] || fail "no second root at $ROOT (run the build stage first)"; }
need_absent() { [ ! -e "$1" ] || fail "$1 already exists; choose a new location, nothing is overwritten"; }
is_locked() {
  "$PYTHON_BIN" - "$ROOT/protocol/PROTOCOL_LOCK.json" <<'PY'
import json, sys
lock = json.load(open(sys.argv[1]))
sys.exit(0 if lock.get("status") == "locked" and "second_reconstruction_plan_sha256" in lock else 1)
PY
}
need_lock() { is_locked || fail "the second root has no lock of the second reconstruction; run the lock stage first"; }
arms_guard() {   # a directory that holds arms is continued only on request
  if [ -d "$1/arms" ] && [ "$RESUME" != 1 ]; then fail "$1/arms exists; set RESUME=1 (or --resume) to continue it, or use a new root"; fi
}
export PYTHONUNBUFFERED=1
printf 'stage=%s root=%s main=%s python=%s workers=%s\n' "$STAGE" "$ROOT" "$MAIN" "$PYTHON_BIN" "$WORKERS"

case "$STAGE" in
build)
  say "build: the second root from the main repository and Human-GEM 2.1.1"
  [ -n "$HUMANGEM" ] || fail "HUMANGEM (or --humangem) must name a Human-GEM checkout at tag v2.1.1"
  [ -f "$HUMANGEM/model/Human-GEM.xml" ] || fail "no model/Human-GEM.xml under $HUMANGEM"
  need_modules cobra numpy scipy
  need_absent "$ROOT"; need_absent "$BUILD_OUT"
  check_plan
  if [ -d "$HUMANGEM/.git" ] || [ -f "$HUMANGEM/.git" ]; then
    got="$(git -C "$HUMANGEM" rev-parse HEAD)"
    [ "$got" = "$HUMANGEM_COMMIT" ] || fail "Human-GEM is at $got, expected $HUMANGEM_COMMIT (tag v2.1.1)"
    echo "Human-GEM commit verified: $got"
  else
    echo "note: $HUMANGEM is not a git checkout; the file digests below are the only check"
  fi
  "$PYTHON_BIN" - "$HUMANGEM" "$HERE/build_report.json" <<'PY'
import hashlib, json, pathlib, sys
hg, rep = pathlib.Path(sys.argv[1]), json.load(open(sys.argv[2]))["humangem"]["files"]
bad = [f for f, h in rep.items() if hashlib.sha256((hg / f).read_bytes()).hexdigest() != h]
sys.exit("Human-GEM files differ from build_report.json: " + ", ".join(bad)) if bad else print("Human-GEM input digests match build_report.json")
PY
  mkdir -p "$BUILD_OUT"
  env -u PYTHONPATH "$PYTHON_BIN" "$HERE/build_root.py" --main "$MAIN" --humangem "$HUMANGEM" --out "$ROOT" --second "$BUILD_OUT"
  say "build: comparison of the four build outputs with the copies deposited in second/"
  for f in exchange_map.tsv target_map.tsv medium_additions.tsv build_report.json; do
    if cmp -s "$BUILD_OUT/$f" "$HERE/$f"; then echo "identical: $f"; else echo "DIFFERS:   $f (compare $BUILD_OUT/$f with $HERE/$f)"; fi
  done
  ;;
development)
  say "development: 11 profiles, five power encodings, the midrank encoding, and the uniform cost"
  need_root; need_modules cobra numpy scipy highspy osqp
  arms_guard "$ROOT/results/development_v3"
  export PYTHONPATH="$ROOT"
  "$PYTHON_BIN" "$HERE/run_arms.py" --root "$ROOT" --partition development --output results/development_v3 --workers "$WORKERS"
  ;;
lock)
  say "lock: pin the plan and the inputs of the second root (once, before any reserved profile is run)"
  need_root; check_plan
  need_absent "$ROOT/protocol/SECOND_RECONSTRUCTION_PLAN.md"
  if is_locked; then fail "the second root is already locked"; fi
  "$PYTHON_BIN" "$HERE/run_arms.py" --root "$ROOT" --write-lock "$PLAN"
  ;;
reserved)
  say "reserved: 47 profiles under the lock of the second reconstruction; no measurement is read"
  need_root; need_modules cobra numpy scipy highspy osqp; need_lock; check_plan
  arms_guard "$ROOT/results/reserved_v3"
  export PYTHONPATH="$ROOT"
  "$PYTHON_BIN" "$HERE/run_arms.py" --root "$ROOT" --partition test --output results/reserved_v3 --workers "$WORKERS"
  ;;
pilot)
  say "pilot: range stages on the 11 development profiles; no measurement is read"
  need_root; need_modules cobra numpy scipy
  need_absent "$ROOT/results/mechanism/bottleneck_pilot/range_stage_ledger.tsv"
  bash "$HERE/run_ranges.sh" "$ROOT" "$PYTHON_BIN" "$WORKERS" pilot
  ;;
ranges)
  say "ranges: nested range stages, uniform control, cap sweep and solver check on the reserved arms; no measurement is read"
  need_root; need_modules cobra numpy scipy highspy osqp; need_lock
  [ -s "$ROOT/results/reserved_v3/predictions.tsv" ] || fail "no reserved predictions; run the reserved stage first"
  need_absent "$ROOT/results/mechanism/reserved_range"
  need_absent "$ROOT/results/closure/cap_sweep"
  need_absent "$ROOT/results/closure/solver_check_96"
  bash "$HERE/run_ranges.sh" "$ROOT" "$PYTHON_BIN" "$WORKERS" reserved
  ;;
evaluate)
  say "evaluate: the stage that reads the measured exchanges (after the reserved arms and ranges are complete)"
  need_root; need_modules numpy scipy; need_lock
  [ -s "$ROOT/results/reserved_v3/predictions.tsv" ] || fail "no reserved predictions; run the reserved stage first"
  compgen -G "$ROOT/results/closure/cap_sweep/jobs/*.json" >/dev/null || fail "no cap sweep; run the ranges stage first"
  PILOT="$ROOT/results/mechanism/bottleneck_pilot/range_stage_ledger.tsv"
  [ -s "$PILOT" ] || fail "no pilot ledger at $PILOT; run the pilot stage first (the interval-position step reads it)"
  for out in results/reserved_v3/evaluation_primary.json results/reserved_v3/transform_primary.json results/second_summary.json \
             results/closure/cap_sweep/cap_concordance.json results/mechanism/range_location.json; do need_absent "$ROOT/$out"; done
  say "evaluate 1/6: place the two measurement files of the main root, checked against the digests in the locked analysis plan"
  for f in data/outcomes_long.tsv data/assay_noise_development.tsv; do
    [ -f "$MAIN/$f" ] || fail "main root lacks $f"
    if [ -e "$ROOT/$f" ]; then echo "present: $f"; else cp -p "$MAIN/$f" "$ROOT/$f"; echo "copied:  $f"; fi
  done
  "$PYTHON_BIN" - "$ROOT" <<'PY'
import hashlib, json, pathlib, sys
root = pathlib.Path(sys.argv[1]); plan = json.load(open(root / "protocol/analysis_plan.json"))["input_sha256"]
bad = [p for p, h in plan.items() if hashlib.sha256((root / p).read_bytes()).hexdigest() != h]
sys.exit("differs from the locked analysis plan: " + ", ".join(bad)) if bad else print("inputs match protocol/analysis_plan.json")
PY
  export PYTHONPATH="$ROOT"
  PRED="$ROOT/results/reserved_v3/predictions.tsv"
  say "evaluate 2/6: concordance of every encoding and the two contrasts (src.evaluation, primary scenario, reserved partition)"
  "$PYTHON_BIN" "$HERE/in_second_root.py" "$ROOT" src.evaluation --root "$ROOT" --partition test --predictions "$PRED" \
      --scenario primary --output "$ROOT/results/reserved_v3/evaluation_primary.json"
  say "evaluate 3/6: predictions changed by the encoding and orderings reversed (src.transform_sensitivity)"
  "$PYTHON_BIN" "$HERE/in_second_root.py" "$ROOT" src.transform_sensitivity --root "$ROOT" --partition test --predictions "$PRED" \
      --scenario primary --output "$ROOT/results/reserved_v3/transform_primary.json"
  say "evaluate 4/6: targets constant across profiles and their relation to the uniform solution (second/summarize.py)"
  "$PYTHON_BIN" "$HERE/summarize.py" --root "$ROOT" --partition test --predictions results/reserved_v3/predictions.tsv --output results/second_summary.json
  say "evaluate 5/6: position of the surviving interval (mechanism/range_location.py; reads the pilot ledger)"
  ( cd "$ROOT/results/mechanism" && STUDY_ROOT="$ROOT" PILOT_LEDGER="$PILOT" \
      RESERVED_LEDGER="$ROOT/results/mechanism/reserved_range/range_stage_ledger.tsv" "$PYTHON_BIN" "$ROOT/mechanism/range_location.py" )
  say "evaluate 6/6: pairs certified by the interval rule under each cap (closure/cap_concordance.py)"
  "$PYTHON_BIN" "$ROOT/closure/cap_concordance.py" --root "$ROOT" --sweep "$ROOT/results/closure/cap_sweep" --out "$ROOT/results/closure/cap_sweep"
  ;;
esac
say "done ($STAGE)"
