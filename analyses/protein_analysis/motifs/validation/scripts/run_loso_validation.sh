#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PROJECT_ROOT="$(cd "$ANALYSIS_DIR/../../../.." && pwd)"
INPUT_DIR="$ANALYSIS_DIR/inputs"
RESULTS_DIR="$ANALYSIS_DIR/results"
STATUS_FILE="$RESULTS_DIR/run_status.tsv"
PID_FILE="$RESULTS_DIR/run.pid"
JOBS="${MOTIF_JOBS:-4}"
SEED="${MOTIF_SEED:-20260826}"
GENES=(BGN DCN FMOD PRELP EPYC LUM OGN)
MOTIF_ENV_BIN="${MOTIF_ENV_BIN:-$HOME/.local/share/mamba/envs/thesis-motif/bin}"
if [[ -d "$MOTIF_ENV_BIN" ]]; then
  export PATH="$MOTIF_ENV_BIN:$PATH"
fi
for required_tool in python3 meme mast; do
  if ! command -v "$required_tool" >/dev/null 2>&1; then
    echo "Required tool not found on PATH: $required_tool" >&2
    exit 4
  fi
done
mkdir -p "$INPUT_DIR" "$RESULTS_DIR"

if [[ -s "$PID_FILE" ]]; then
  existing_pid="$(tr -d '[:space:]' < "$PID_FILE")"
  if [[ "$existing_pid" =~ ^[0-9]+$ ]] && kill -0 "$existing_pid" 2>/dev/null; then
    echo "A protein-pattern validation controller is already running with PID $existing_pid." >&2
    exit 3
  fi
fi
printf '%s\n' "$$" > "$PID_FILE"

if [[ ! -s "$STATUS_FILE" ]]; then
  printf 'timestamp\tfold\tgene\tstage\tstate\telapsed_seconds\tmessage\n' > "$STATUS_FILE"
fi

record() {
  local fold="$1" gene="$2" stage="$3" state="$4" elapsed="$5" message="$6"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(date --iso-8601=seconds)" "$fold" "$gene" "$stage" "$state" "$elapsed" "$message" \
    >> "$STATUS_FILE"
}

on_error() {
  local exit_code=$?
  record controller all controller failed "" "exit code $exit_code; completed fold outputs retained"
  exit "$exit_code"
}
trap on_error ERR

run_model() {
  local fold="$1" gene="$2" train="$3" test="$4" output="$5"
  mkdir -p "$output"
  if [[ ! -s "$output/meme.complete" ]]; then
    local started
    started="$(date +%s)"
    record "$fold" "$gene" meme running "" started
    meme "$train" -protein -oc "$output/meme" -mod zoops -nmotifs 12 \
      -minw 6 -maxw 40 -evt 0.05 -seed "$SEED" -maxsize 1000000 -nostatus
    local elapsed=$(( $(date +%s) - started ))
    printf 'complete\n' > "$output/meme.complete"
    record "$fold" "$gene" meme complete "$elapsed" finished
  fi
  if [[ ! -s "$output/mast.complete" ]]; then
    local started
    started="$(date +%s)"
    record "$fold" "$gene" mast running "" started
    mast -oc "$output/mast" -comp -ev 1000000 "$output/meme/meme.txt" "$test"
    local elapsed=$(( $(date +%s) - started ))
    printf 'complete\n' > "$output/mast.complete"
    record "$fold" "$gene" mast complete "$elapsed" finished
  fi
}

record controller all preparation running "" started
python3 "$SCRIPT_DIR/prepare_loso_validation.py" \
  --mature-fasta "$PROJECT_ROOT/analyses/protein_analysis/motifs/inputs/slrp_canonical_103_mature.faa" \
  --output-dir "$INPUT_DIR"
record controller all preparation complete "" finished

mapfile -t FOLDS < <(find "$INPUT_DIR/folds" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)
for fold in "${FOLDS[@]}"; do
  test="$INPUT_DIR/folds/$fold/test.faa"
  record "$fold" all fold running "" started
  active=0
  pids=()
  for gene in "${GENES[@]}"; do
    run_model \
      "$fold" "$gene" \
      "$INPUT_DIR/folds/$fold/train/$gene.faa" \
      "$test" \
      "$RESULTS_DIR/folds/$fold/$gene" &
    pids+=("$!")
    active=$((active + 1))
    if (( active >= JOBS )); then
      wait "${pids[0]}"
      pids=("${pids[@]:1}")
      active=$((active - 1))
    fi
  done
  for pid in "${pids[@]}"; do
    wait "$pid"
  done
  printf 'complete\n' > "$RESULTS_DIR/folds/$fold/fold.complete"
  record "$fold" all fold complete "" finished
done

record final_models all final_models running "" started
active=0
pids=()
for gene in "${GENES[@]}"; do
  run_model \
    final_models "$gene" \
    "$INPUT_DIR/final_models/train/$gene.faa" \
    "$INPUT_DIR/challenge_controls.faa" \
    "$RESULTS_DIR/final_models/$gene" &
  pids+=("$!")
  active=$((active + 1))
  if (( active >= JOBS )); then
    wait "${pids[0]}"
    pids=("${pids[@]:1}")
    active=$((active - 1))
  fi
done
for pid in "${pids[@]}"; do
  wait "$pid"
done
record final_models all final_models complete "" finished

record controller all summary running "" started
python3 "$SCRIPT_DIR/summarize_loso_validation.py" --analysis-dir "$ANALYSIS_DIR"
record controller all summary complete "" finished
printf 'complete\n' > "$RESULTS_DIR/run.complete"
record controller all controller complete "" "all leave-one-species-out stages complete"
