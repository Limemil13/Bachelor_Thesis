#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
IQTREE_BIN="${IQTREE_BIN:-iqtree2}"

if [[ $# -gt 0 ]]; then
  GENES=("$@")
else
  GENES=(BGN FMOD PRELP EPYC LUM OGN)
fi

cd "${PROJECT_ROOT}"

for gene in "${GENES[@]}"; do
  case "${gene}" in
    BGN|EPYC|FMOD|LUM|OGN|PRELP) ;;
    *) printf 'Unsupported gene: %s\n' "${gene}" >&2; exit 2 ;;
  esac

  BASE="analyses/phylogenetics/per_gene_trees/${gene}"
  RAW="${BASE}/candidates/${gene}_large_candidates_raw.faa"
  FILTERED="${BASE}/candidates/${gene}_large_candidates_locus_filtered.faa"
  REPORT="${BASE}/candidates/${gene}_large_locus_filter_report.tsv"
  ALIGNED="${BASE}/alignments/${gene}_large_locus_filtered_aligned.faa"
  TRIMMED="${BASE}/trimmed_alignments/${gene}_large_locus_filtered_aligned_trimmed.faa"
  FINAL_DIR="${BASE}/corrected_final_tree"
  PREFIX="${FINAL_DIR}/${gene}_large_corrected"

  if [[ ! -s "${RAW}" ]]; then
    printf 'Missing raw NCBI protein FASTA: %s\n' "${RAW}" >&2
    exit 1
  fi

  mkdir -p "${FINAL_DIR}"
  "${PYTHON_BIN}" analyses/phylogenetics/per_gene_trees/scripts/filter_large_slrp_fasta.py \
    --input "${RAW}" \
    --output "${FILTERED}" \
    --report "${REPORT}" \
    --gene "${gene}" \
    --one-per-locus

  mafft --auto "${FILTERED}" > "${ALIGNED}"
  trimal -in "${ALIGNED}" -out "${TRIMMED}" -automated1
  "${IQTREE_BIN}" \
    -s "${TRIMMED}" \
    -m MFP \
    -B 1000 \
    -alrt 1000 \
    -T AUTO \
    --prefix "${PREFIX}"
done
