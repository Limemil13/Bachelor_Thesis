#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-class1}"
if [[ "$MODE" != "class1" && "$MODE" != "class2" && "$MODE" != "all" ]]; then
  echo "Usage: $0 [class1|class2|all]" >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PROJECT_ROOT="$(cd "$ANALYSIS_DIR/../../.." && pwd)"
INPUT_DIR="$ANALYSIS_DIR/inputs"
RESULTS_ROOT="$ANALYSIS_DIR/results"
STATUS_FILE="$RESULTS_ROOT/run_status.tsv"
PID_FILE="$RESULTS_ROOT/run.pid"
THREADS="${IQTREE_THREADS:-4}"
mkdir -p "$INPUT_DIR" "$RESULTS_ROOT"

if [[ -s "$PID_FILE" ]]; then
  existing_pid="$(tr -d '[:space:]' < "$PID_FILE")"
  if [[ "$existing_pid" =~ ^[0-9]+$ ]] && kill -0 "$existing_pid" 2>/dev/null; then
    echo "A targeted-extension controller is already running with PID $existing_pid." >&2
    exit 3
  fi
fi
printf '%s\n' "$$" > "$PID_FILE"

if [[ ! -s "$STATUS_FILE" ]]; then
  printf 'timestamp\tanalysis\tstage\tstate\telapsed_seconds\tmessage\n' > "$STATUS_FILE"
fi

current_analysis="controller"
current_stage="initialization"
record() {
  local state="$1"
  local elapsed="$2"
  local message="$3"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(date --iso-8601=seconds)" "$current_analysis" "$current_stage" "$state" "$elapsed" "$message" \
    >> "$STATUS_FILE"
}

on_error() {
  local exit_code=$?
  record "failed" "" "exit code $exit_code; completed outputs were retained"
  exit "$exit_code"
}
trap on_error ERR

run_timed() {
  local analysis="$1"
  local stage="$2"
  shift 2
  current_analysis="$analysis"
  current_stage="$stage"
  local started
  started="$(date +%s)"
  record "running" "" "started"
  "$@"
  local elapsed=$(( $(date +%s) - started ))
  record "complete" "$elapsed" "finished"
}

prepare_inputs() {
  python3 "$SCRIPT_DIR/build_targeted_datasets.py" \
    --project-root "$PROJECT_ROOT" \
    --analysis "$MODE" \
    --output-dir "$INPUT_DIR"
}

run_tree() {
  local analysis="$1"
  local result_dir="$RESULTS_ROOT/$analysis"
  local prefix="$result_dir/${analysis}_targeted"
  local dataset="$INPUT_DIR/${analysis}_targeted_dataset.faa"
  mkdir -p "$result_dir"

  if [[ -s "$result_dir/run.complete" ]]; then
    current_analysis="$analysis"
    current_stage="all"
    record "skipped" "0" "existing completed result retained"
    return
  fi

  if [[ ! -s "$result_dir/mafft.complete" ]]; then
    run_timed "$analysis" "mafft" bash -c \
      'mafft --genafpair --maxiterate 1000 "$1" > "$2" 2> "$3"' \
      _ "$dataset" "$result_dir/aligned.faa" "$result_dir/mafft.log"
    printf 'complete\n' > "$result_dir/mafft.complete"
  fi

  if [[ ! -s "$result_dir/trimming.complete" ]]; then
    run_timed "$analysis" "gap_trimming" python3 \
      "$PROJECT_ROOT/analyses/phylogenetics/combined_trees/rooted_class_i_tree/scripts/trim_gap_alignment.py" \
      "$result_dir/aligned.faa" "$result_dir/aligned_trimmed.faa" \
      --maximum-gap-fraction 0.50
    printf 'complete\n' > "$result_dir/trimming.complete"
  fi

  if [[ ! -s "$result_dir/iqtree.complete" ]]; then
    if [[ -s "$prefix.treefile" && -s "$prefix.iqtree" ]] && grep -q 'Total wall-clock time used' "$prefix.iqtree"; then
      current_analysis="$analysis"
      current_stage="iqtree"
      record "complete" "" "complete IQ-TREE output already present"
    else
      run_timed "$analysis" "iqtree" iqtree2 \
        -s "$result_dir/aligned_trimmed.faa" \
        -m MFP -B 1000 -alrt 1000 -seed 12345 -T "$THREADS" \
        --prefix "$prefix"
    fi
    printf 'complete\n' > "$result_dir/iqtree.complete"
  fi

  run_timed "$analysis" "tree_summary" python3 "$SCRIPT_DIR/summarize_targeted_tree.py" \
    --analysis "$analysis" --tree "$prefix.treefile" --output-dir "$result_dir"

  run_timed "$analysis" "figure" python3 "$SCRIPT_DIR/render_targeted_tree.py" \
    --analysis "$analysis" --nodes "$result_dir/tree_nodes.tsv" \
    --edges "$result_dir/tree_edges.tsv" --output "$result_dir/${analysis}_targeted_tree.pdf"

  printf 'complete\n' > "$result_dir/run.complete"
  current_analysis="$analysis"
  current_stage="all"
  record "complete" "" "all requested stages complete"
}

run_timed "controller" "dataset_preparation" prepare_inputs
if [[ "$MODE" == "class1" || "$MODE" == "all" ]]; then
  run_tree "class1"
fi
if [[ "$MODE" == "class2" || "$MODE" == "all" ]]; then
  run_tree "class2"
fi
current_analysis="controller"
current_stage="all"
record "complete" "" "requested mode $MODE complete"
