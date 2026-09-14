#!/usr/bin/env bash
set -euo pipefail

# Sequential launcher for the annotation-aware SynVoy reruns prepared on
# 2026-09-02. Run from WSL; outputs use new names and never replace old results.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
SYNVOY_ROOT="${SYNVOY_ROOT:?Set SYNVOY_ROOT to the SynVoy checkout}"
CONDA_SH="${CONDA_SH:-}"
RUN_TAG="dev_20260903"
DRY_RUN=false
GENES=()

usage() {
  cat <<'EOF'
Usage: run_updated_synvoy_panel.sh [--dry-run] [--run-tag TAG] [GENE ...]

The default genes are BGN FMOD PRELP EPYC DCN. Each run is written to a new
results/<gene>_human_15species_<TAG> directory; existing output is never
overwritten.
EOF
}

while (($#)); do
  case "$1" in
    --dry-run)
      DRY_RUN=true
      shift
      ;;
    --run-tag)
      if (($# < 2)); then
        printf '%s\n' 'Missing value after --run-tag.' >&2
        usage >&2
        exit 2
      fi
      RUN_TAG="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      while (($#)); do
        GENES+=("$1")
        shift
      done
      ;;
    -*)
      printf 'Unknown option: %s\n' "$1" >&2
      usage >&2
      exit 2
      ;;
    *)
      GENES+=("$1")
      shift
      ;;
  esac
done

if ((${#GENES[@]} == 0)); then
  # These runs predate the August adjudication/coverage/ownership update.
  GENES=(BGN FMOD PRELP EPYC DCN)
fi

if [[ ! "$RUN_TAG" =~ ^[A-Za-z0-9._-]+$ ]]; then
  printf 'Unsafe run tag: %s\n' "$RUN_TAG" >&2
  exit 2
fi

declare -A QUERY=(
  [BGN]="bgn_human.faa"
  [FMOD]="fmod_human.faa"
  [PRELP]="prelp_human_clean.faa"
  [EPYC]="epyc_human.faa"
  [DCN]="dcn_human.faa"
  [OGN]="ogn_human_clean.faa"
  [LUM]="lum_human_clean.faa"
)

declare -A TOKENS=(
  [BGN]="BGN,biglycan"
  [FMOD]="FMOD,fibromodulin"
  [PRELP]="PRELP,prolargin"
  [EPYC]="EPYC,epiphycan,DSPG3"
  [DCN]="DCN,decorin"
  [OGN]="OGN,osteoglycin,mimecan"
  [LUM]="LUM,lumican"
)

if [[ -n "${CONDA_SH}" ]]; then
  source "${CONDA_SH}"
fi
cd "$SYNVOY_ROOT"

python3 "${PROJECT_ROOT}/analyses/synteny/scripts/validate_synvoy_target_gffs.py" \
  --fna-dir "$SYNVOY_ROOT/pro_panel/targets_15/fna" \
  --gff-dir "$SYNVOY_ROOT/pro_panel/targets_15/gff" \
  --expected 14

mkdir -p results/rerun_launcher_logs
PROVENANCE="results/rerun_launcher_logs/${RUN_TAG}_provenance.tsv"
if [[ ! -e "$PROVENANCE" ]]; then
  printf 'date\tbranch\tcommit\tgenes\n' > "$PROVENANCE"
fi
printf '%s\t%s\t%s\t%s\n' \
  "$(date --iso-8601=seconds)" \
  "$(git branch --show-current)" \
  "$(git rev-parse HEAD)" \
  "${GENES[*]}" >> "$PROVENANCE"

for gene in "${GENES[@]}"; do
  gene="${gene^^}"
  if [[ -z "${QUERY[$gene]:-}" ]]; then
    printf 'Unsupported gene: %s\n' "$gene" >&2
    exit 2
  fi

  lower_gene="${gene,,}"
  outdir="results/${lower_gene}_human_15species_${RUN_TAG}"
  log="results/rerun_launcher_logs/${lower_gene}_${RUN_TAG}.log"
  if [[ -e "$outdir" ]]; then
    printf 'Refusing to overwrite existing output: %s\n' "$outdir" >&2
    printf 'Resume that run explicitly or choose a new RUN_TAG.\n' >&2
    exit 3
  fi

  cmd=(
    ./run_synvoy.sh
    -name "${lower_gene}_15species_${RUN_TAG}"
    -profile standard
    -c stage_genomes_timeout.config
    --mode pro
    --query "pro_panel/queries/${QUERY[$gene]}"
    --home_genome pro_panel/genomes/fna/human.fna
    --home_gff pro_panel/genomes/gff/human.gff
    --home_species "Homo sapiens"
    --target_genomes "pro_panel/targets_15/fna/*.fna"
    --target_gffs "pro_panel/targets_15/gff/*.gff"
    --home_goi_gene "$gene"
    --goi_family_tokens "${TOKENS[$gene]}"
    --strict_goi_family false
    --auto_apply_preset false
    --mmseqs_split_memory_limit 8G
    --iterative_search_cpus 1
    --outdir "$outdir"
  )

  printf 'Launching %s -> %s\n' "$gene" "$outdir"
  printf 'Command:'
  printf ' %q' "${cmd[@]}"
  printf '\n'
  if [[ "$DRY_RUN" == false ]]; then
    "${cmd[@]}" 2>&1 | tee "$log"
  fi
done
