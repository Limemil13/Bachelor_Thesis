#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
cd "${PROJECT_ROOT}"

GENE="${1:-LUM}"
TREE="analyses/phylogenetics/per_gene_trees/${GENE}/final_tree/${GENE}_large_final.treefile"
IQTREE_REPORT="analyses/phylogenetics/per_gene_trees/${GENE}/final_tree/${GENE}_large_final.iqtree"
REPORT="analyses/phylogenetics/per_gene_trees/${GENE}/candidates/${GENE}_large_filter_report.tsv"
OUT="analyses/phylogenetics/per_gene_trees/itol/${GENE}_large_final_itol"
ITOL_TREE="${OUT}.tree"

if [ ! -s "${TREE}" ]; then
  echo "ERROR: completed tree is missing: ${TREE}" >&2
  exit 1
fi

if [ ! -s "${IQTREE_REPORT}" ]; then
  echo "ERROR: IQ-TREE is incomplete; final report is missing: ${IQTREE_REPORT}" >&2
  exit 1
fi

if [ ! -s "${REPORT}" ]; then
  echo "ERROR: filter report is missing: ${REPORT}" >&2
  exit 1
fi

mkdir -p "$(dirname "${OUT}")"
tr '|' '_' < "${TREE}" > "${ITOL_TREE}"

python3 analyses/phylogenetics/per_gene_trees/scripts/make_large_itol_annotations.py \
  --tree "${ITOL_TREE}" \
  --gene "${GENE}" \
  --out-prefix "${OUT}"

python3 analyses/phylogenetics/per_gene_trees/scripts/make_large_taxgroup_annotations.py \
  --tree "${ITOL_TREE}" \
  --gene "${GENE}" \
  --filter-report "${REPORT}" \
  --out-prefix "${OUT}"

echo "DONE: ${GENE} large-tree iTOL bundle under $(dirname "${OUT}")"
