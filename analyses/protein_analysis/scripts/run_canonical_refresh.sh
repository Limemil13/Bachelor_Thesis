#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
PFAM_DB="${PFAM_DB:-data/domains/Pfam-A.hmm}"
HMMER_CPU="${HMMER_CPU:-7}"

cd "${PROJECT_ROOT}"

"${PYTHON_BIN}" analyses/protein_analysis/scripts/canonical_dataset.py

DOMAIN_DIR="analyses/protein_analysis/domains/domain_scan_canonical"
mkdir -p "${DOMAIN_DIR}"
hmmscan \
  -o "${DOMAIN_DIR}/domain_candidates_canonical.hmmscan.log" \
  --domtblout "${DOMAIN_DIR}/domain_candidates_canonical.pfam.domtblout" \
  --cpu "${HMMER_CPU}" \
  "${PFAM_DB}" \
  analyses/protein_analysis/candidates/slrp_candidates_canonical.faa

ALIGNMENT_DIR="analyses/protein_analysis/alignments/canonical"
mkdir -p "${ALIGNMENT_DIR}"
mafft --quiet --auto analyses/protein_analysis/candidates/canonical_by_gene/BGN_canonical.faa \
  > "${ALIGNMENT_DIR}/BGN_canonical_aligned.faa"
mafft --quiet --auto analyses/protein_analysis/candidates/canonical_by_gene/DCN_canonical.faa \
  > "${ALIGNMENT_DIR}/DCN_canonical_aligned.faa"
mafft --quiet --auto analyses/protein_analysis/candidates/canonical_by_gene/EPYC_canonical.faa \
  > "${ALIGNMENT_DIR}/EPYC_canonical_aligned.faa"
mafft --quiet --auto analyses/protein_analysis/candidates/canonical_by_gene/FMOD_canonical.faa \
  > "${ALIGNMENT_DIR}/FMOD_canonical_aligned.faa"
mafft --quiet --auto analyses/protein_analysis/candidates/canonical_by_gene/OGN_canonical.faa \
  > "${ALIGNMENT_DIR}/OGN_canonical_aligned.faa"
mafft --quiet --auto analyses/protein_analysis/candidates/canonical_by_gene/PRELP_canonical.faa \
  > "${ALIGNMENT_DIR}/PRELP_canonical_aligned.faa"
mafft --quiet --auto analyses/protein_analysis/candidates/canonical_by_gene/LUM_canonical.faa \
  > "${ALIGNMENT_DIR}/LUM_canonical_aligned.faa"

"${PYTHON_BIN}" analyses/protein_analysis/scripts/parse_domain_scan.py
"${PYTHON_BIN}" analyses/protein_analysis/scripts/msa_qc.py
"${PYTHON_BIN}" analyses/protein_analysis/scripts/build_msa_conservation_outputs.py
"${PYTHON_BIN}" analyses/protein_analysis/signal_peptides/scripts/prepare_signalp_input.py
"${PYTHON_BIN}" analyses/protein_analysis/signal_peptides/scripts/merge_signalp_summaries.py \
  --existing analyses/protein_analysis/signal_peptides/tables/baseline_77/signalp_summary_clean_77.tsv \
  --lum analyses/protein_analysis/signal_peptides/tables/signalp_lum_slow_summary_clean.tsv \
  --supplemental analyses/protein_analysis/signal_peptides/tables/signalp_pending_17_slow_summary_clean.tsv \
  --manifest analyses/protein_analysis/candidates/canonical_candidate_manifest.tsv \
  --out-dir analyses/protein_analysis/signal_peptides/tables
"${PYTHON_BIN}" analyses/protein_analysis/scripts/build_protein_conservation_summary.py
"${PYTHON_BIN}" analyses/protein_analysis/scripts/make_sequence_logos.py

printf 'Canonical refresh complete. SignalP and ProtSpace external reruns remain documented separately.\n'
