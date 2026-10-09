#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
PFAM_DB="${PFAM_DB:-data/domains/Pfam-A.hmm}"
HMMER_CPU="${HMMER_CPU:-7}"

cd "${PROJECT_ROOT}"

# Rebuild the canonical manifest and synchronized per-gene/combined FASTAs.
"${PYTHON_BIN}" analyses/protein_analysis/scripts/canonical_dataset.py

DOMAIN_DIR="analyses/protein_analysis/domains/domain_scan_canonical"
# Scan every canonical protein against the pinned Pfam HMM database.
mkdir -p "${DOMAIN_DIR}"
hmmscan \
  -o "${DOMAIN_DIR}/domain_candidates_canonical.hmmscan.log" \
  --domtblout "${DOMAIN_DIR}/domain_candidates_canonical.pfam.domtblout" \
  --cpu "${HMMER_CPU}" \
  "${PFAM_DB}" \
  analyses/protein_analysis/candidates/slrp_candidates_canonical.faa

ALIGNMENT_DIR="analyses/protein_analysis/alignments/canonical"
# Align each gene separately so conservation and gap metrics compare ortholog sets.
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

# Parse the retained tool outputs, prepare the SignalP submission, then rebuild
# summaries and logos. SignalP itself is an external run; its final 103-row
# summary is retained in the repository.
"${PYTHON_BIN}" analyses/protein_analysis/scripts/parse_domain_scan.py
"${PYTHON_BIN}" analyses/protein_analysis/scripts/msa_qc.py
"${PYTHON_BIN}" analyses/protein_analysis/scripts/build_msa_conservation_outputs.py
"${PYTHON_BIN}" analyses/protein_analysis/signal_peptides/scripts/prepare_signalp_input.py
"${PYTHON_BIN}" analyses/protein_analysis/scripts/build_protein_conservation_summary.py
"${PYTHON_BIN}" analyses/protein_analysis/scripts/make_sequence_logos.py

printf 'Canonical refresh complete. SignalP and ProtSpace external reruns remain documented separately.\n'
