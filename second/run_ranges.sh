#!/bin/bash
# Second reconstruction: the range stages and the solver check of the main analysis, run in the second
# root with the scripts of the main analysis (the copies in the root). No measurement is read here.
#
#   bash second/run_ranges.sh SECOND_ROOT PYTHON WORKERS pilot      # development pilot of the range stages
#   bash second/run_ranges.sh SECOND_ROOT PYTHON WORKERS reserved   # after the reserved arms (run_arms.py --partition test)
set -euo pipefail
ROOT=$(cd "$1" && pwd); PY=$2; W=${3:-12}; STAGE=${4:-reserved}
OUT=$ROOT/results; M=$ROOT/mechanism; C=$ROOT/closure
export PYTHONPATH=$ROOT
if [ "$STAGE" = pilot ]; then
  P=$OUT/mechanism/bottleneck_pilot; mkdir -p "$P"
  [ -f "$P/AMENDMENT.md" ] || printf '%s\n' "Pilot of the range stages on the second reconstruction (Human-GEM): the script of the main analysis, unchanged, on the 11 development profiles of the second root. It reads no measurement." > "$P/AMENDMENT.md"
  "$PY" "$M/nested_range_pilot.py" --root "$ROOT" --output "$P"
  exit 0
fi
"$PY" "$M/reserved_range_stage.py" --root "$ROOT" --output "$OUT/mechanism/reserved_range" --workers "$W" --budget 14400
U=$("$PY" - "$OUT/reserved_v3/arms" <<'PY'
import json, sys, glob
for f in sorted(glob.glob(sys.argv[1] + "/*.json")):
    r = json.load(open(f))
    if r.get("arm") == "uniform_pfba" and r.get("scenario") == "primary": print(f); break
PY
)
"$PY" "$M/uniform_control.py" --root "$ROOT" --output "$OUT/mechanism/uniform_control" --uniform-json "$U"
"$PY" "$C/cap_sweep.py" --root "$ROOT" --out "$OUT/closure/cap_sweep" --workers "$W"
"$PY" "$C/solver_check.py" --root "$ROOT" --out "$OUT/closure/solver_check_96" --solver hybrid --panel chemistry --workers 8
echo RANGES_DONE
