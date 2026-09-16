#!/usr/bin/env bash
# Retrieve the two public repositories that the independent evaluation reads, at the commits the
# study used, and refuse to continue if either checkout does not resolve to the expected commit.
#
#   bash copeland/fetch_inputs.sh [DESTINATION]     default: work/copeland
#
# copeland_prepare.py reads data/growth_rates.rda from the compendium and data/lf_hyp_bay_rnaseq.rda
# from the expression package and records `git rev-parse HEAD` of each; it does not clone or check
# out anything itself. This script is the acquisition step that SOURCES.md describes.
set -euo pipefail
DEST="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/work/copeland}"
COMPENDIUM_URL="https://github.com/oldhamlab/Copeland.2023.hypoxia.flux.git"
COMPENDIUM_COMMIT="354bd99459367f373c7b8a641bca4667bf96c779"
RNASEQ_URL="https://github.com/wmoldham/rnaseq.lf.hypoxia.molidustat.git"
RNASEQ_COMMIT="1b6563bac8ba50019e90f7e6c6ed7d2f49a90ea5"
mkdir -p "$DEST"
fetch() {
  local name="$1" url="$2" commit="$3" dir="$DEST/$1"
  if [ ! -d "$dir/.git" ]; then
    git clone --quiet "$url" "$dir"
  fi
  git -C "$dir" fetch --quiet origin
  git -C "$dir" checkout --quiet --detach "$commit"
  local head
  head="$(git -C "$dir" rev-parse HEAD)"
  if [ "$head" != "$commit" ]; then
    echo "$name: checked-out commit $head is not the expected $commit" >&2
    exit 3
  fi
  printf '%s  %s  %s\n' "$name" "$head" "$dir"
}
fetch compendium "$COMPENDIUM_URL" "$COMPENDIUM_COMMIT"
fetch rnaseq     "$RNASEQ_URL"     "$RNASEQ_COMMIT"
# the two files the preparation step reads; the expression package's digest is the one recorded in
# results/copeland/prepared/prepare_summary.json (rnaseq_rda_sha256), the growth-rate table's the
# digest of that file at the pinned commit
RNASEQ_RDA_SHA256="4088ad11271e5a6def001018b8911fec8e0d7868462a1f390b86f9a146e4901f"
GROWTH_RDA_SHA256="d051b317077f40e6e5d7047a07a32ceb6adeb87c07cd31e903949122ca4757b6"   # data/growth_rates.rda at the pinned commit
digest() { if command -v shasum >/dev/null 2>&1; then shasum -a 256 "$1" | cut -d' ' -f1; else sha256sum "$1" | cut -d' ' -f1; fi; }
for f in "$DEST/compendium/data/growth_rates.rda" "$DEST/rnaseq/data/lf_hyp_bay_rnaseq.rda"; do
  [ -s "$f" ] || { echo "missing input file: $f" >&2; exit 3; }
  printf '%s  %s\n' "$(digest "$f")" "$f"
done
got="$(digest "$DEST/rnaseq/data/lf_hyp_bay_rnaseq.rda")"
[ "$got" = "$RNASEQ_RDA_SHA256" ] || { echo "lf_hyp_bay_rnaseq.rda digest $got differs from the recorded $RNASEQ_RDA_SHA256" >&2; exit 3; }
got="$(digest "$DEST/compendium/data/growth_rates.rda")"
[ "$got" = "$GROWTH_RDA_SHA256" ] || { echo "growth_rates.rda digest $got differs from the recorded $GROWTH_RDA_SHA256" >&2; exit 3; }
echo FETCH_OK
