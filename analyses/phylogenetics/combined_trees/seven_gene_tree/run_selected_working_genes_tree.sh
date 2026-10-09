#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
cd "${PROJECT_ROOT}"

BASE="analyses/phylogenetics/combined_trees/seven_gene_tree"
SOURCE="analyses/protein_analysis/candidates/slrp_candidates_canonical.faa"
TRIMAL="${TRIMAL_BIN:-trimal}"

cp "${SOURCE}" "${BASE}/SLRP_selected_working_genes.faa"

# Align the combined canonical set and remove unreliable alignment columns.
mafft --quiet --auto "${BASE}/SLRP_selected_working_genes.faa" \
  > "${BASE}/SLRP_selected_working_genes_aligned.faa"

"${TRIMAL}" \
  -in "${BASE}/SLRP_selected_working_genes_aligned.faa" \
  -out "${BASE}/SLRP_selected_working_genes_aligned_trimmed.faa" \
  -automated1

# Infer the maximum-likelihood tree with model selection and two support measures.
iqtree2 \
  -s "${BASE}/SLRP_selected_working_genes_aligned_trimmed.faa" \
  -m MFP \
  -B 1000 \
  -alrt 1000 \
  -T AUTO \
  -redo \
  --prefix "${BASE}/SLRP_selected_working_genes"

echo "Selected-working-gene tree complete."
