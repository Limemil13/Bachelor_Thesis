#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
cd "${PROJECT_ROOT}"

if [ "$#" -gt 0 ]; then
  GENES=("$@")
else
  GENES=(BGN FMOD OGN EPYC PRELP LUM)
fi

for gene in "${GENES[@]}"; do
  echo "======================================"
  echo "Large tree for $gene"
  echo "======================================"

  BASE="analyses/phylogenetics/per_gene_trees/${gene}"

  mkdir -p "${BASE}/ncbi_download"
  mkdir -p "${BASE}/candidates"
  mkdir -p "${BASE}/alignments"
  mkdir -p "${BASE}/trimmed_alignments"
  mkdir -p "${BASE}/final_tree"
  mkdir -p "${BASE}/logs"

  ZIP="${BASE}/ncbi_download/${gene}_orthologs.zip"
  UNZIP_DIR="${BASE}/ncbi_download/${gene}_orthologs"

  echo "=== Downloading NCBI orthologs for $gene ==="

  conda run -n protspace_env datasets download gene symbol "${gene}" \
    --taxon human \
    --ortholog all \
    --filename "${ZIP}"

  echo "=== Unzipping $gene ==="

  rm -rf "${UNZIP_DIR}"
  unzip -o "${ZIP}" -d "${UNZIP_DIR}" > "${BASE}/logs/unzip.log"

  RAW="${UNZIP_DIR}/ncbi_dataset/data/protein.faa"

  if [ ! -f "${RAW}" ]; then
    echo "ERROR: protein.faa not found for $gene"
    exit 1
  fi

  cp "${RAW}" "${BASE}/candidates/${gene}_large_candidates_raw.faa"

  echo "=== Filtering $gene ==="

  python analyses/phylogenetics/per_gene_trees/scripts/filter_large_slrp_fasta.py \
    --input "${BASE}/candidates/${gene}_large_candidates_raw.faa" \
    --output "${BASE}/candidates/${gene}_large_candidates_filtered.faa" \
    --report "${BASE}/candidates/${gene}_large_filter_report.tsv" \
    --gene "${gene}" \
    --one-per-locus

  echo "=== MAFFT $gene ==="

  mafft --auto "${BASE}/candidates/${gene}_large_candidates_filtered.faa" \
    > "${BASE}/alignments/${gene}_large_filtered_aligned.faa"

  echo "=== trimAl $gene ==="

  trimal \
    -in "${BASE}/alignments/${gene}_large_filtered_aligned.faa" \
    -out "${BASE}/trimmed_alignments/${gene}_large_filtered_aligned_trimmed.faa" \
    -automated1

  echo "=== IQ-TREE $gene ==="

  iqtree2 \
    -s "${BASE}/trimmed_alignments/${gene}_large_filtered_aligned_trimmed.faa" \
    -m MFP \
    -B 1000 \
    -alrt 1000 \
    -T AUTO \
    --prefix "${BASE}/final_tree/${gene}_large_final"

  echo "DONE: $gene"
done
