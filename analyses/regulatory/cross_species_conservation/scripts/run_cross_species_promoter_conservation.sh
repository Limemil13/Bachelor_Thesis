#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PROJECT_ROOT="$(cd "$ANALYSIS_DIR/../../.." && pwd)"
SOURCE_DIR="$PROJECT_ROOT/analyses/regulatory/inputs/per_gene"
RESULTS_DIR="$ANALYSIS_DIR/results"
STATUS_FILE="$RESULTS_DIR/run_status.tsv"
PID_FILE="$RESULTS_DIR/run.pid"
GENES=(BGN DCN FMOD PRELP EPYC LUM OGN)
WINDOWS=(core proximal)

mkdir -p "$RESULTS_DIR"
if [[ -s "$PID_FILE" ]]; then
  existing_pid="$(tr -d '[:space:]' < "$PID_FILE")"
  if [[ "$existing_pid" =~ ^[0-9]+$ ]] && kill -0 "$existing_pid" 2>/dev/null; then
    echo "A promoter-conservation controller is already running with PID $existing_pid." >&2
    exit 3
  fi
fi
printf '%s\n' "$$" > "$PID_FILE"

if [[ ! -s "$STATUS_FILE" ]]; then
  printf 'timestamp\twindow\tgene\tstage\tstate\tmessage\n' > "$STATUS_FILE"
fi

record() {
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(date --iso-8601=seconds)" "$1" "$2" "$3" "$4" "$5" >> "$STATUS_FILE"
}

on_error() {
  local exit_code=$?
  record controller all controller failed "exit code $exit_code; completed alignments retained"
  exit "$exit_code"
}
trap on_error ERR

for tool in mafft python3; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "Required tool not found on PATH: $tool" >&2
    exit 4
  fi
done

record controller all preparation running started
for window in "${WINDOWS[@]}"; do
  for gene in "${GENES[@]}"; do
    input="$SOURCE_DIR/$window/${gene}_positive.fna"
    output_dir="$RESULTS_DIR/alignments/$window"
    output="$output_dir/${gene}.aln.fna"
    log="$output_dir/${gene}.mafft.log"
    mkdir -p "$output_dir"
    if [[ ! -s "$input" ]]; then
      echo "Missing promoter FASTA: $input" >&2
      exit 5
    fi
    if [[ -s "$output" ]]; then
      record "$window" "$gene" mafft cached existing_alignment
      continue
    fi
    record "$window" "$gene" mafft running started
    mafft --quiet --localpair --maxiterate 1000 --thread 1 "$input" > "$output" 2> "$log"
    record "$window" "$gene" mafft complete finished
  done
done
record controller all preparation complete alignments_finished

record controller all summary running started
python3 "$SCRIPT_DIR/summarize_promoter_conservation.py" --analysis-dir "$ANALYSIS_DIR"
record controller all summary complete finished
printf 'complete\n' > "$RESULTS_DIR/run.complete"
record controller all controller complete all_stages_complete
