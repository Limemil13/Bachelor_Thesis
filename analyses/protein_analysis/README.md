# Protein analysis

This directory contains the current **103-protein, seven-gene** comparative
protein package (BGN, DCN, EPYC, FMOD, LUM, OGN, and PRELP).

- `candidates/`: curated per-gene FASTAs, canonical concatenated FASTA, and manifest.
- `domains/`: Pfam/HMMER outputs.
- `alignments/canonical/`: MAFFT MSAs, QC, consensus, identity, and conservation.
- `signal_peptides/`: SignalP inputs, retained raw outputs, mappings, and merged calls.
- `protspace/`: 107-sequence input (103 SLRPs plus four controls), embeddings, QC, and plots.
- `motifs/`: MEME/STREME/MAST analysis of the seven-gene mature-protein set.
- `tables/`: protein- and gene-level integrated summaries.
- `figures/`: protein overview and sequence logos.
- `manual_review/`: alignment-review packet and recorded decisions.
- `diagnostics/`: targeted annotation investigations; not canonical output.
- `archive/`: superseded historical runs retained for provenance.
- `scripts/`: reproducible builders and refresh scripts.

The protein-level source of truth is
`tables/protein_conservation_domain_msa_signalp_summary.tsv`. It currently has
103 rows, and SignalP coverage is complete (100 positive and three negative).
The completed 17-protein upload and raw result bundle are retained under
`signal_peptides/input/` and `signal_peptides/raw_outputs/`. See
`CANONICAL_DATASET.md` for accessions, curation rules,
metrics, and interpretation boundaries. Trees are under `../phylogenetics/`.
