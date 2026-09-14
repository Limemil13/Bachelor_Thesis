#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
MEME_BIN="${MEME_BIN:-meme}"
MAST_BIN="${MAST_BIN:-mast}"
STREME_BIN="${STREME_BIN:-streme}"
BASE="$REPO/analyses/protein_analysis/motifs"
INPUTS="$BASE/inputs"
RAW="${MOTIF_RAW_DIR:-$BASE/raw_outputs_7gene_20260914}"
SEED="${MOTIF_SEED:-20260826}"

mkdir -p "$INPUTS" "$RAW"
"$PYTHON_BIN" "$BASE/scripts/prepare_protein_motif_inputs.py" \
  --fasta "$REPO/analyses/protein_analysis/candidates/slrp_candidates_canonical.faa" \
  --protein-table "$REPO/analyses/protein_analysis/tables/protein_conservation_domain_msa_signalp_summary.tsv" \
  --output-dir "$INPUTS"

"$MEME_BIN" "$INPUTS/slrp_canonical_103_mature.faa" -protein \
  -oc "$RAW/meme_all_slrp" -mod zoops -nmotifs 12 -minw 6 -maxw 40 \
  -evt 0.05 -seed "$SEED" -maxsize 1000000 -nostatus
"$MAST_BIN" -oc "$RAW/mast_all_slrp" -comp -ev 100 \
  "$RAW/meme_all_slrp/meme.txt" "$INPUTS/slrp_canonical_103_mature.faa"

for gene in BGN DCN FMOD PRELP EPYC LUM OGN; do
  "$MEME_BIN" "$INPUTS/per_gene/${gene}_mature.faa" -protein \
    -oc "$RAW/meme_${gene}" -mod zoops -nmotifs 12 -minw 6 -maxw 40 \
    -evt 0.05 -seed "$SEED" -maxsize 1000000 -nostatus
  "$MAST_BIN" -oc "$RAW/mast_${gene}_motifs_all_proteins" -comp -ev 1000000 \
    "$RAW/meme_${gene}/meme.txt" "$INPUTS/slrp_canonical_103_mature.faa"
  "$STREME_BIN" --protein \
    --p "$INPUTS/per_gene/${gene}_mature.faa" \
    --n "$INPUTS/per_gene/${gene}_other_genes_control.faa" \
    --minw 6 --maxw 30 --thresh 0.05 --patience 5 --seed "$SEED" \
    --oc "$RAW/streme_${gene}_vs_other_genes"
done

echo "Protein motif analysis complete: $BASE"
