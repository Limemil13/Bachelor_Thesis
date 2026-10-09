#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
cd "${PROJECT_ROOT}"

BASE="analyses/phylogenetics/combined_trees/slrp_family_reference_tree"
DOWNLOADS="${BASE}/downloads"

mkdir -p "${DOWNLOADS}"

# Human SLRP reference/context genes.
GENES=(
  ASPN BGN DCN ECM2
  FMOD LUM KERA OMD PRELP
  OGN EPYC OPTC
  CHAD CHADL NYX TSKU PODN PODNL1
  LINGO1
)

echo "=== Download missing human reference proteins from NCBI Datasets ==="

for gene in "${GENES[@]}"; do
  # Prefer existing local queries; otherwise download the human gene package.
  if [ -f "data/queries/${gene}.faa" ]; then
    echo "Local query exists for ${gene}, skipping download."
    continue
  fi

  mkdir -p "${DOWNLOADS}/${gene}"

  ZIP="${DOWNLOADS}/${gene}/${gene}_human.zip"
  OUTDIR="${DOWNLOADS}/${gene}/${gene}_human"

  if [ -f "${OUTDIR}/ncbi_dataset/data/protein.faa" ]; then
    echo "Already downloaded ${gene}."
    continue
  fi

  echo "Downloading ${gene}..."

  if conda run -n protspace_env datasets download gene symbol "${gene}" \
      --taxon "Homo sapiens" \
      --filename "${ZIP}"; then

    rm -rf "${OUTDIR}"
    unzip -o "${ZIP}" -d "${OUTDIR}" > "${DOWNLOADS}/${gene}/unzip.log"

  else
    echo "WARNING: download failed for ${gene}, continuing."
  fi
done

echo "=== Build selected reference FASTA ==="

# Select one plausible human reference protein for every context gene.
python "${BASE}/scripts/build_slrp_reference_context_fasta.py"

echo "=== Check selected proteins ==="

column -t -s $'\t' "${BASE}/SLRP_reference_context_selected.tsv" | less -S

echo "=== MAFFT E-INS-i alignment ==="

# Use a consistency-based alignment suitable for related but divergent proteins.
mafft --genafpair --maxiterate 1000 \
  "${BASE}/SLRP_reference_context.faa" \
  > "${BASE}/SLRP_reference_context_aligned.faa"

echo "=== trimAl gap filtering, remove columns with >95% gaps ==="

trimal \
  -in "${BASE}/SLRP_reference_context_aligned.faa" \
  -out "${BASE}/SLRP_reference_context_aligned_trimmed.faa" \
  -gt 0.05

echo "=== IQ-TREE ==="

# Infer the context tree with automatic model selection and branch support.
iqtree2 \
  -s "${BASE}/SLRP_reference_context_aligned_trimmed.faa" \
  -m MFP \
  -B 1000 \
  -alrt 1000 \
  -T AUTO \
  --prefix "${BASE}/SLRP_reference_context"

echo "DONE."
echo "Tree: ${BASE}/SLRP_reference_context.treefile"
