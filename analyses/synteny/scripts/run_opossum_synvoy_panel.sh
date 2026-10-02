#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
SYNVOY_ROOT="${SYNVOY_ROOT:?Set SYNVOY_ROOT to the SynVoy checkout}"
RUN_TAG="${RUN_TAG:-matched_GCF_027887165.2_20260922}"
SKIP_INPUT_VALIDATION="${SKIP_INPUT_VALIDATION:-false}"

# SynVoy checks for Conda before launching Nextflow. Set CONDA_BIN when a
# non-interactive shell does not inherit the Conda installation directory.
if [[ -n "${CONDA_BIN:-}" ]]; then
  export PATH="${CONDA_BIN}:${PATH}"
elif [[ -n "${CONDA_EXE:-}" ]]; then
  export PATH="$(dirname "${CONDA_EXE}"):${PATH}"
fi

declare -A QUERY=(
  [BGN]="bgn_human.faa"
  [DCN]="dcn_human.faa"
  [FMOD]="fmod_human.faa"
  [PRELP]="prelp_human_clean.faa"
  [EPYC]="epyc_human.faa"
  [LUM]="lum_human_clean.faa"
  [OGN]="ogn_human_clean.faa"
)

declare -A TOKENS=(
  [BGN]="BGN,biglycan"
  [DCN]="DCN,decorin"
  [FMOD]="FMOD,fibromodulin"
  [PRELP]="PRELP,prolargin"
  [EPYC]="EPYC,epiphycan,DSPG3"
  [LUM]="LUM,lumican"
  [OGN]="OGN,osteoglycin,mimecan"
)

GENES=("$@")
if ((${#GENES[@]} == 0)); then
  GENES=(BGN DCN FMOD PRELP EPYC LUM OGN)
fi

cd "$SYNVOY_ROOT"
if [[ "$SKIP_INPUT_VALIDATION" != "true" ]]; then
  python3 "${PROJECT_ROOT}/analyses/synteny/scripts/validate_synvoy_target_gffs.py" \
    --fna-dir "$SYNVOY_ROOT/pro_panel/targets_15/fna" \
    --gff-dir "$SYNVOY_ROOT/pro_panel/targets_15/gff" \
    --expected 14
fi

mkdir -p results/rerun_launcher_logs
for gene in "${GENES[@]}"; do
  gene="${gene^^}"
  if [[ -z "${QUERY[$gene]:-}" ]]; then
    printf 'Unsupported gene: %s\n' "$gene" >&2
    exit 2
  fi

  lower_gene="${gene,,}"
  outdir="results/${lower_gene}_human_opossum_${RUN_TAG}"
  log="results/rerun_launcher_logs/${lower_gene}_opossum_${RUN_TAG}.log"
  nextflow_name="${lower_gene}_opossum_${RUN_TAG//./_}"
  if [[ -e "$outdir" ]]; then
    printf 'Refusing to overwrite existing output: %s\n' "$outdir" >&2
    exit 3
  fi

  cmd=(
    ./run_synvoy.sh
    -name "$nextflow_name"
    -profile standard
    -c stage_genomes_timeout.config
    --mode pro
    --query "pro_panel/queries/${QUERY[$gene]}"
    --home_genome pro_panel/genomes/fna/human.fna
    --home_gff pro_panel/genomes/gff/human.gff
    --home_species "Homo sapiens"
    --target_genomes pro_panel/targets_15/fna/opossum.fna
    --target_gffs pro_panel/targets_15/gff/opossum.gff
    --home_goi_gene "$gene"
    --goi_family_tokens "${TOKENS[$gene]}"
    --strict_goi_family false
    --auto_apply_preset false
    --mmseqs_split_memory_limit 8G
    --iterative_search_cpus 1
    --outdir "$outdir"
  )

  printf 'Launching %s opossum-only rerun -> %s\n' "$gene" "$outdir"
  printf 'Command:'
  printf ' %q' "${cmd[@]}"
  printf '\n'
  "${cmd[@]}" 2>&1 | tee "$log"
done
