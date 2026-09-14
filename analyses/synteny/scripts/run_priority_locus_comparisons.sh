#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 SYNVOY_ROOT THESIS_ROOT" >&2
  exit 2
fi

synvoy_root=$1
thesis_root=$2
comparison_script="$thesis_root/analyses/synteny/scripts/compare_gff_neighborhoods.py"
human_gff="$synvoy_root/pro_panel/genomes/gff/human.gff"
output_root="$thesis_root/analyses/synteny/diagnostics/priority_locus_reviews"

run_comparison() {
  local gene=$1
  local species=$2
  local human_selector=$3
  local target_selector=$4
  local target_gff="$synvoy_root/pro_panel/targets_15/gff/${species}.gff"
  local output_dir="$output_root/${gene}_${species}"

  python3 "$comparison_script" \
    --locus "human_${gene}" "$human_gff" "$human_selector" \
    --locus "${species}_${gene}" "$target_gff" "$target_selector" \
    --flank-count 15 \
    --output-dir "$output_dir"
}

run_comparison BGN catshark BGN@NC_000023.11 LOC119955784@NC_052166.1

run_comparison DCN coelacanth DCN@NC_000012.12 DCN@NC_088145.1
run_comparison DCN spotted_gar DCN@NC_000012.12 dcn@NC_090702.1
run_comparison DCN zebrafish DCN@NC_000012.12 dcn@NC_133179.1

run_comparison EPYC coelacanth EPYC@NC_000012.12 EPYC@NC_088145.1
run_comparison EPYC elephant_shark EPYC@NC_000012.12 LOC103179840@NW_024704746.1
run_comparison EPYC zebrafish EPYC@NC_000012.12 epyc@NC_133179.1

run_comparison FMOD catshark FMOD@NC_000001.11 LOC119978081@NC_052160.1
run_comparison FMOD zebrafish FMOD@NC_000001.11 fmodb@NC_133183.1

run_comparison LUM catshark LUM@NC_000012.12 lum@NC_052165.1
run_comparison LUM coelacanth LUM@NC_000012.12 LUM@NC_088145.1
run_comparison LUM spotted_gar LUM@NC_000012.12 lum@NC_090702.1

run_comparison PRELP catshark PRELP@NC_000001.11 prelp@NC_052160.1
run_comparison PRELP zebrafish PRELP@NC_000001.11 prelp@NC_133186.1
