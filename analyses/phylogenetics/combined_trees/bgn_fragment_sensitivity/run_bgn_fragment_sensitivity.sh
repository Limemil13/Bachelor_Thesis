#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
cd "${PROJECT_ROOT}"

BASE="analyses/phylogenetics/combined_trees/bgn_fragment_sensitivity"
INPUT_DIR="${BASE}/inputs"
RESULT_DIR="${BASE}/results"
SOURCE="analyses/protein_analysis/candidates/slrp_candidates_canonical.faa"
CANONICAL_TREE="analyses/phylogenetics/combined_trees/seven_gene_tree/SLRP_selected_working_genes.treefile"
PREFIX="${RESULT_DIR}/SLRP_partial_BGN_excluded"

PYTHON="${PYTHON_BIN:-python3}"
MAFFT="${MAFFT_BIN:-mafft}"
TRIMAL="${TRIMAL_BIN:-trimal}"
IQTREE="${IQTREE_BIN:-iqtree2}"
SEED="${IQTREE_SEED:-797597}"

mkdir -p "${INPUT_DIR}" "${RESULT_DIR}"

"${PYTHON}" "${BASE}/build_fragment_excluded_input.py" \
  --input "${SOURCE}" \
  --output "${INPUT_DIR}/slrp_candidates_partial_BGN_excluded.faa" \
  --excluded-table "${RESULT_DIR}/excluded_partial_BGN_records.tsv" \
  --count-table "${RESULT_DIR}/input_gene_counts.tsv"

"${MAFFT}" --quiet --auto \
  "${INPUT_DIR}/slrp_candidates_partial_BGN_excluded.faa" \
  > "${RESULT_DIR}/SLRP_partial_BGN_excluded_aligned.faa"

"${TRIMAL}" \
  -in "${RESULT_DIR}/SLRP_partial_BGN_excluded_aligned.faa" \
  -out "${RESULT_DIR}/SLRP_partial_BGN_excluded_aligned_trimmed.faa" \
  -automated1

"${IQTREE}" \
  -s "${RESULT_DIR}/SLRP_partial_BGN_excluded_aligned_trimmed.faa" \
  -m MFP \
  -B 1000 \
  -alrt 1000 \
  -T AUTO \
  -seed "${SEED}" \
  -redo \
  --prefix "${PREFIX}"

"${PYTHON}" "${BASE}/compare_bgn_trees.py" \
  --canonical-tree "${CANONICAL_TREE}" \
  --sensitivity-tree "${PREFIX}.treefile" \
  --out-dir "${RESULT_DIR}"

{
  printf 'field\tvalue\n'
  printf 'canonical_input\t%s\n' "${SOURCE}"
  printf 'canonical_tree\t%s\n' "${CANONICAL_TREE}"
  printf 'sensitivity_input_count\t101\n'
  printf 'excluded_accessions\tXP_038638277.1,XP_048475905.1\n'
  printf 'iqtree_seed\t%s\n' "${SEED}"
  printf 'mafft_version\t%s\n' "$("${MAFFT}" --version 2>&1 | head -n 1)"
  printf 'trimal_version\t%s\n' "$("${TRIMAL}" --version 2>&1 | awk 'NF {print; exit}')"
  printf 'iqtree_version\t%s\n' "$("${IQTREE}" --version 2>&1 | head -n 1)"
} > "${RESULT_DIR}/run_manifest.tsv"

echo "BGN fragment-exclusion sensitivity analysis complete."
