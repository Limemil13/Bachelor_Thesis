#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
cd "${PROJECT_ROOT}"

BASE="analyses/phylogenetics/combined_trees/six_gene_tree"
SOURCE="analyses/protein_analysis/candidates/slrp_candidates_canonical.faa"
TRIMAL="${TRIMAL_BIN:-trimal}"
PYTHON="${PYTHON_BIN:-python3}"

cp "${SOURCE}" "${BASE}/SLRP_selected_working_genes.faa"

mafft --quiet --auto "${BASE}/SLRP_selected_working_genes.faa" \
  > "${BASE}/SLRP_selected_working_genes_aligned.faa"

"${TRIMAL}" \
  -in "${BASE}/SLRP_selected_working_genes_aligned.faa" \
  -out "${BASE}/SLRP_selected_working_genes_aligned_trimmed.faa" \
  -automated1

iqtree2 \
  -s "${BASE}/SLRP_selected_working_genes_aligned_trimmed.faa" \
  -m MFP \
  -B 1000 \
  -alrt 1000 \
  -T AUTO \
  -redo \
  --prefix "${BASE}/SLRP_selected_working_genes"

cp "${BASE}/SLRP_selected_working_genes.treefile" \
   "${BASE}/SLRP_selected_working_genes_itol.tree"
sed -i 's/|/_/g' "${BASE}/SLRP_selected_working_genes_itol.tree"

"${PYTHON}" "${BASE}/make_selected_gene_itol_annotations.py" \
  --tree "${BASE}/SLRP_selected_working_genes_itol.tree" \
  --out-prefix "${BASE}/SLRP_selected_working_genes_itol"

echo "Selected-working-gene tree complete."
