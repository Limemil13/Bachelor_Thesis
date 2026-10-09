#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  printf 'Usage: %s <BGN|DCN|EPYC|FMOD|LUM|OGN|PRELP>\n' "$0" >&2
  exit 2
fi

GENE="$1"
case "${GENE}" in
  BGN|DCN|EPYC|FMOD|LUM|OGN|PRELP) ;;
  *)
    printf 'Unsupported compact-panel gene: %s\n' "${GENE}" >&2
    exit 2
    ;;
esac

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
IQTREE_BIN="${IQTREE_BIN:-iqtree2}"

ALIGNMENT="${PROJECT_ROOT}/analyses/protein_analysis/alignments/canonical/${GENE}_canonical_aligned.faa"
TREE_DIR="${PROJECT_ROOT}/analyses/phylogenetics/compact_panel/trees/${GENE}"
TREE_PREFIX="${TREE_DIR}/${GENE}"

if [[ ! -s "${ALIGNMENT}" ]]; then
  # Compact trees reuse the validated canonical MSA for the selected gene.
  printf 'Missing canonical alignment: %s\n' "${ALIGNMENT}" >&2
  exit 1
fi

mkdir -p "${TREE_DIR}"
# IQ-TREE chooses the model and calculates ultrafast-bootstrap/aLRT support.
"${IQTREE_BIN}" \
  -s "${ALIGNMENT}" \
  -m MFP \
  -B 1000 \
  -alrt 1000 \
  -T AUTO \
  -redo \
  --prefix "${TREE_PREFIX}"

printf 'Compact canonical tree complete: %s\n' "${TREE_PREFIX}.treefile"
