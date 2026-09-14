#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSPACE="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
SYNVOY_ROOT="${1:-${SYNVOY_ROOT:-}}"
if [[ -z "${SYNVOY_ROOT}" ]]; then
  printf 'Set SYNVOY_ROOT or pass the SynVoy checkout as argument 1.\n' >&2
  exit 2
fi
RESULTS="${SYNVOY_ROOT}/results/gene_structure"
PARTS="${RESULTS}/lum_by_species"
mkdir -p "${PARTS}"

cd "${SYNVOY_ROOT}"
python3 "${WORKSPACE}/analyses/overview/scripts/extract_gene_structure_streaming.py" \
  --gff-list "${RESULTS}/main_annotation_gffs.txt" \
  --gene LUM \
  --out "${RESULTS}/LUM_gene_structure_5species.tsv"

python3 "${WORKSPACE}/analyses/overview/scripts/merge_lum_gene_structure.py" \
  --existing "${RESULTS}/main_slrp_gene_structure_5species.tsv" \
  --lum "${RESULTS}/LUM_gene_structure_5species.tsv" \
  --output "${RESULTS}/main_slrp_gene_structure_5species.tsv"

python3 "${WORKSPACE}/analyses/overview/scripts/select_gene_structure_representatives.py" \
  --input "${RESULTS}/main_slrp_gene_structure_5species.tsv" \
  --out-dir "${RESULTS}"
python3 scripts/plot_gene_structure_5species.py

cat "${RESULTS}/main_slrp_gene_structure_5species_summary_clean.tsv"
