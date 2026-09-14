#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
SYNVOY_ROOT="${SYNVOY_ROOT:?Set SYNVOY_ROOT to the SynVoy checkout}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"
FASTA_SHUFFLE_BIN="${FASTA_SHUFFLE_BIN:-fasta-shuffle-letters}"
STREME_BIN="${STREME_BIN:-streme}"
TOMTOM_BIN="${TOMTOM_BIN:-tomtom}"
AME_BIN="${AME_BIN:-ame}"
FIMO_BIN="${FIMO_BIN:-fimo}"
BASE="$REPO/analyses/regulatory"
INPUTS="$BASE/inputs"
RAW="$BASE/raw_outputs_7gene_dcn_20260902"
DB="$BASE/databases"
JASPAR="$DB/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt"
TARGETED="$DB/JASPAR2026_cartilage_growth_plate_targeted.meme"
SEED=20260826

mkdir -p "$INPUTS" "$RAW" "$DB"

"$PYTHON_BIN" "$BASE/scripts/extract_promoters.py" \
  --representatives "$REPO/analyses/gene_structure/tables/gene_structure_representative_5species.tsv" \
  --genome-dir "$SYNVOY_ROOT/pro_panel/genomes/fna" \
  --samtools "$SAMTOOLS_BIN" \
  --output-dir "$INPUTS"

if [[ ! -s "$JASPAR" ]]; then
  curl -L --fail --retry 4 \
    'https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt' \
    -o "$JASPAR"
fi

"$PYTHON_BIN" "$BASE/scripts/subset_meme_database.py" --input "$JASPAR" --output "$TARGETED"

PROMOTERS="$INPUTS/slrp_proximal_promoters_5species.fna"
SHUFFLED="$INPUTS/slrp_proximal_promoters_dinucleotide_shuffled.fna"
"$FASTA_SHUFFLE_BIN" -dna -kmer 2 -seed "$SEED" "$PROMOTERS" "$SHUFFLED"

"$STREME_BIN" --dna --p "$PROMOTERS" --n "$SHUFFLED" \
  --minw 6 --maxw 15 --thresh 0.05 --patience 5 --seed "$SEED" \
  --oc "$RAW/streme_all_vs_shuffled"
"$TOMTOM_BIN" -oc "$RAW/tomtom_all_vs_jaspar" -thresh 0.10 \
  "$RAW/streme_all_vs_shuffled/streme.txt" "$JASPAR"

for gene in BGN DCN FMOD PRELP EPYC LUM OGN; do
  "$STREME_BIN" --dna \
    --p "$INPUTS/per_gene/proximal/${gene}_positive.fna" \
    --n "$INPUTS/per_gene/proximal/${gene}_other_genes_control.fna" \
    --minw 6 --maxw 15 --thresh 0.05 --patience 5 --seed "$SEED" \
    --oc "$RAW/streme_${gene}_vs_other_genes"
  "$TOMTOM_BIN" -oc "$RAW/tomtom_${gene}_vs_jaspar" -thresh 0.10 \
    "$RAW/streme_${gene}_vs_other_genes/streme.txt" "$JASPAR"
done

mkdir -p "$RAW/ame_all_jaspar_vs_shuffled"
# AME 5.5.8 can fail while rendering its HTML report when long FASTA labels
# contain locus punctuation. Text mode preserves the complete statistical table
# and avoids that report-only failure.
"$AME_BIN" --control --shuffle-- --kmer 2 --seed "$SEED" \
  --method fisher --scoring totalhits --evalue-report-threshold 10 --text \
  "$PROMOTERS" "$JASPAR" > "$RAW/ame_all_jaspar_vs_shuffled/ame.tsv"

mkdir -p "$RAW/ame_targeted_vs_shuffled"
"$AME_BIN" --control --shuffle-- --kmer 2 --seed "$SEED" \
  --method fisher --scoring totalhits --evalue-report-threshold 1000000 --text \
  "$PROMOTERS" "$TARGETED" > "$RAW/ame_targeted_vs_shuffled/ame.tsv"

"$FIMO_BIN" --qv-thresh --thresh 0.05 --oc "$RAW/fimo_targeted" \
  "$TARGETED" "$PROMOTERS"

# Preserve an explicitly exploratory, uncorrected scan for cross-species
# recurrence summaries. It must not be described as statistically significant.
"$FIMO_BIN" --thresh 1e-4 --oc "$RAW/fimo_targeted_p1e4_exploratory" \
  "$TARGETED" "$PROMOTERS"

echo "Regulatory motif analysis complete: $BASE"
