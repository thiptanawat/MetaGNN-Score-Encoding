#!/usr/bin/env bash
# Download the per-arm records from the GitHub release that matches this checkout, verify their
# digests against RELEASE_MANIFEST.md, and unpack them: those of the locked study into
# results/reserved_v3/arms and results/development_v3/arms (with the two run_summary.json files),
# those of the second reconstruction into results/second/reserved_v3 and
# results/second/development_v3 (arms/ and run_summary.json in each).
#
#   bash verify/fetch_release_assets.sh [reserved|development|second|all]     default: all
#
# The archives are large (160 megabytes together). They are outputs, not inputs: reproduce.sh
# --stage all --partition development|test regenerates the first two from the raw inputs, and
# reports/fresh_arm_reproduction_v3/comparison.json records that a fresh arm matched byte for byte;
# second/reproduce_second.sh regenerates the third from a Human-GEM checkout.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"
WHICH="${1:-all}"
TAG="$(sed -n 's/^Release tag: *`\([^`]*\)`.*/\1/p' RELEASE_MANIFEST.md | head -1)"
BASE_URL="https://github.com/thiptanawat/MetaGNN-Score-Encoding/releases/download/$TAG"
digest() { if command -v shasum >/dev/null 2>&1; then shasum -a 256 "$1" | cut -d' ' -f1; else sha256sum "$1" | cut -d' ' -f1; fi; }
expected() { grep -E "^\| \`$1\` " RELEASE_MANIFEST.md | head -1 | awk -F'|' '{gsub(/[` ]/,"",$4); print $4}'; }
fetch() {
  local name="$1" want
  want="$(expected "$name")"
  [ -n "$want" ] || { echo "RELEASE_MANIFEST.md has no digest for $name" >&2; exit 2; }
  mkdir -p work/assets
  if [ ! -s "work/assets/$name" ]; then
    echo "downloading $BASE_URL/$name"
    curl -L --fail --retry 3 -o "work/assets/$name" "$BASE_URL/$name"
  fi
  local got; got="$(digest "work/assets/$name")"
  [ "$got" = "$want" ] || { echo "$name: digest $got differs from the recorded $want" >&2; exit 3; }
  echo "$name: digest verified"
}
get() {
  local name="$1" part="$2"
  fetch "$name"
  tar -xzf "work/assets/$name" -C "results/$part"
  echo "unpacked into results/$part ($(ls "results/$part/arms" | wc -l | tr -d ' ') files under arms/)"
}
# The archive of the second reconstruction carries its paths from the repository root.
get_second() {
  local name="second_reconstruction_arms.tar.gz" part
  fetch "$name"
  tar -xzf "work/assets/$name" -C "$HERE"
  for part in reserved_v3 development_v3; do
    echo "unpacked into results/second/$part ($(ls "results/second/$part/arms" | wc -l | tr -d ' ') files under arms/)"
  done
}
case "$WHICH" in
  reserved)    get reserved_v3_arms.tar.gz reserved_v3;;
  development) get development_v3_arms.tar.gz development_v3;;
  second)      get_second;;
  all)         get reserved_v3_arms.tar.gz reserved_v3; get development_v3_arms.tar.gz development_v3; get_second;;
  *) echo "usage: $0 [reserved|development|second|all]" >&2; exit 2;;
esac
python3 verify/check_release.py --root "$HERE" || true
