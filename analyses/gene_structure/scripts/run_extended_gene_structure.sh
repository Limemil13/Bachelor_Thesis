#!/usr/bin/env bash
set -euo pipefail

# Rebuild the extended five-species gene-structure package from the same
# SynVoy/NCBI genome annotations used by the thesis. Override these variables
# when the local installation differs.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../../.." && pwd)}"
SYNVOY_ROOT="${SYNVOY_ROOT:?Set SYNVOY_ROOT to the SynVoy checkout}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
PFAM_HMM="${PFAM_HMM:-${PROJECT_ROOT}/data/domains/Pfam-A.hmm}"

cd "${PROJECT_ROOT}"

"${PYTHON_BIN}" analyses/gene_structure/scripts/build_extended_gene_structure.py \
  --synvoy-root "${SYNVOY_ROOT}"

mkdir -p analyses/gene_structure/extended/domains
hmmscan --cpu "${HMMER_CPU:-7}" \
  --domtblout analyses/gene_structure/extended/domains/representative_proteins.pfam.domtblout \
  -o analyses/gene_structure/extended/domains/representative_proteins.hmmscan.log \
  "${PFAM_HMM}" \
  analyses/gene_structure/extended/sequences/representative_proteins_5species.faa

"${PYTHON_BIN}" analyses/gene_structure/scripts/map_domains_to_exons.py \
  --domtblout analyses/gene_structure/extended/domains/representative_proteins.pfam.domtblout

"${PYTHON_BIN}" analyses/evolutionary_rates/scripts/analyze_codon_constraint.py
