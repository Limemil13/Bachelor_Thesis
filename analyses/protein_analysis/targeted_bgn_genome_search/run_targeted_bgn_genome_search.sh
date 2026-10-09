#!/usr/bin/env bash
set -euo pipefail

# Targeted genome-level search for BGN in the two unresolved species.
# Whole-genome TBLASTN finds candidate regions; Miniprot then tests the top windows
# for a spliced protein model while keeping memory use manageable.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="${REPO_ROOT:-$(cd "$SCRIPT_DIR/../../.." && pwd)}"
: "${SYNVOY_ROOT:?Set SYNVOY_ROOT to the SynVoy installation directory}"
SYNVOY="$SYNVOY_ROOT"
WORK="${BGN_SEARCH_WORK:-$REPO/tmp/targeted_bgn_genome_search}"
OUT="$REPO/analyses/protein_analysis/targeted_bgn_genome_search/results"
MINIPROT="${MINIPROT_BIN:-miniprot}"
WINDOW_SCRIPT="$REPO/analyses/protein_analysis/targeted_bgn_genome_search/select_candidate_windows.py"

BGN_PANEL="$REPO/analyses/protein_analysis/candidates/canonical_by_gene/BGN_canonical.faa"
HUMAN_SLRP="$REPO/analyses/protein_analysis/diagnostics/BGN/chicken_annotation_conflict/human_slrp_reference_plus_chicken_bgn_locus.faa"
OPOSSUM_GENOME="${OPOSSUM_GENOME:-$SYNVOY/pro_panel/targets_15/fna/opossum.fna}"
ELEPHANT_GENOME="${ELEPHANT_GENOME:-$SYNVOY/pro_panel/targets_15/fna/elephant_shark.fna}"

mkdir -p "$WORK/queries" "$WORK/db" "$WORK/windows" "$OUT"
# PID and stage files support monitoring and safe resumption.
printf '%s\n' "$$" > "$OUT/run.pid"
touch "$OUT/run_status.tsv"
if [[ ! -s "$OUT/run_status.tsv" ]]; then
  printf 'stage\tspecies\tstatus\ttime\n' > "$OUT/run_status.tsv"
fi

record_status() {
  printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$(date --iso-8601=seconds)" \
    >> "$OUT/run_status.tsv"
}

# Prepare fixed BGN and SLRP-control query sets once.
if [[ ! -s "$WORK/queries/human_slrp_reference.faa" ]]; then
  awk '/^>/{keep=($0 !~ /^>CHICKEN_BGN_LOCUS/)} keep' "$HUMAN_SLRP" \
    > "$WORK/queries/human_slrp_reference.faa"
fi

if [[ ! -s "$WORK/queries/opossum_bgn_queries.faa" ]]; then
  # Use several mammalian BGN queries for opossum.
  awk '
    /^>/{keep=($0 ~ /^>(human|mouse|rattus_norvegicus|canis_lupus_familiaris|bos_taurus)\|BGN\|/)}
    keep
  ' "$BGN_PANEL" > "$WORK/queries/opossum_bgn_queries.faa"
fi

if [[ ! -s "$WORK/queries/elephant_shark_bgn_queries.faa" ]]; then
  # Use a broad vertebrate BGN query set for elephant shark.
  awk '
    /^>/{keep=($0 ~ /^>(human|anolis_carolinensis|frog|zebrafish|lepisosteus_oculatus|latimeria_chalumnae|scyliorhinus_canicula|whale_shark)\|BGN\|/)}
    keep
  ' "$BGN_PANEL" > "$WORK/queries/elephant_shark_bgn_queries.faa"
fi

build_database() {
  local species="$1"
  local genome="$2"
  local database="$WORK/db/$species"

  # Reuse a completed nucleotide BLAST database.
  if [[ -e "${database}.ndb" || -e "${database}.nin" ]]; then
    record_status makeblastdb "$species" SKIPPED_EXISTING
    return
  fi

  record_status makeblastdb "$species" STARTED
  makeblastdb -in "$genome" -dbtype nucl -parse_seqids -out "$database" \
    > "$OUT/${species}.makeblastdb.log" 2>&1
  record_status makeblastdb "$species" COMPLETED
}

run_tblastn() {
  local species="$1"
  local queries="$2"
  local database="$WORK/db/$species"
  local result="$OUT/${species}.bgn_lineage_queries.tblastn.tsv"

  # Preserve a completed TBLASTN table when resuming.
  if [[ -s "$result" ]]; then
    record_status tblastn "$species" SKIPPED_EXISTING
    return
  fi

  record_status tblastn "$species" STARTED
  tblastn -query "$queries" -db "$database" -task tblastn \
    -evalue 1e-5 -seg yes -soft_masking true -max_target_seqs 150 \
    -max_hsps 30 -num_threads 6 \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovhsp stitle' \
    -out "$result"
  record_status tblastn "$species" COMPLETED
}

extract_windows() {
  local species="$1"
  local database="$WORK/db/$species"
  local hits="$OUT/${species}.bgn_lineage_queries.tblastn.tsv"
  local table="$OUT/${species}.candidate_windows.tsv"
  local fasta="$WORK/windows/${species}.candidate_windows.fna"

  # Cluster HSPs and add flanking sequence for Miniprot.
  python3 "$WINDOW_SCRIPT" "$hits" "$table" --top 25 --padding 100000 \
    --cluster-gap 100000

  : > "$fasta"
  tail -n +2 "$table" | while IFS=$'\t' read -r rank seqid start end strand score query_count hsp_count; do
    printf '>%s__%s_%s__rank%s__strand%s\n' "$seqid" "$start" "$end" "$rank" "$strand" \
      >> "$fasta"
    blastdbcmd -db "$database" -entry "$seqid" -range "${start}-${end}" \
      -outfmt '%s' >> "$fasta"
  done

  # Always include the expected opossum syntenic block as a control window.
  if [[ "$species" == "opossum" ]]; then
    printf '>NC_077235.2__15000000_17000000__expected_BGN_synteny\n' >> "$fasta"
    blastdbcmd -db "$database" -entry NC_077235.2 -range 15000000-17000000 \
      -outfmt '%s' >> "$fasta"
  fi
}

run_miniprot_windows() {
  local species="$1"
  local queries="$2"
  local prefix="$3"
  local windows="$WORK/windows/${species}.candidate_windows.fna"
  local bgn_result="$OUT/${species}.bgn_candidate_windows.miniprot.gff"
  local control_result="$OUT/${species}.slrp_control_candidate_windows.miniprot.gff"

  # Align BGN queries to the selected windows.
  if [[ ! -s "$bgn_result" ]]; then
    record_status miniprot_bgn_windows "$species" STARTED
    "$MINIPROT" -t 4 -j 2 -N 50 --outn=50 --outs=0.40 --outc=0.20 \
      --gff --aln --trans -P "$prefix" "$windows" "$queries" > "$bgn_result"
    record_status miniprot_bgn_windows "$species" COMPLETED
  else
    record_status miniprot_bgn_windows "$species" SKIPPED_EXISTING
  fi

  # Align the broader SLRP reference set as a paralog control.
  if [[ ! -s "$control_result" ]]; then
    record_status miniprot_slrp_control "$species" STARTED
    "$MINIPROT" -t 4 -j 2 -N 30 --outn=30 --outs=0.60 --outc=0.25 \
      --gff --aln --trans -P "${prefix}S" "$windows" \
      "$WORK/queries/human_slrp_reference.faa" > "$control_result"
    record_status miniprot_slrp_control "$species" COMPLETED
  else
    record_status miniprot_slrp_control "$species" SKIPPED_EXISTING
  fi
}

for species in opossum elephant_shark; do
  # Run species sequentially to limit resource use.
  if [[ "$species" == "opossum" ]]; then
    genome="$OPOSSUM_GENOME"
    queries="$WORK/queries/opossum_bgn_queries.faa"
    prefix=OPO
  else
    genome="$ELEPHANT_GENOME"
    queries="$WORK/queries/elephant_shark_bgn_queries.faa"
    prefix=ESH
  fi

  build_database "$species" "$genome"
  run_tblastn "$species" "$queries"
  extract_windows "$species"
  run_miniprot_windows "$species" "$queries" "$prefix"
done

record_status all both COMPLETED
