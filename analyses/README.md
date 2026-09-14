# Thesis analysis directory map

This directory contains the current bachelor-thesis analyses. Use the files in
`overview/` to navigate; do not choose a result merely because it has a
shorter filename or an older timestamp.

## Authoritative sections

- `overview/`: one-row-per-gene evidence, decisions, methods,
  reproducibility checks, and result locations.
- `gene_structure/`: the local 40-row representative-transcript analysis,
  coding-block tables, provenance, and five-species figures.
- `protein_analysis/`: compact canonical protein workflow. The one editable
  protein-level source of truth is
  `protein_analysis/tables/protein_conservation_domain_msa_signalp_summary.tsv`.
  The supplementary mature-protein MEME/STREME/MAST module is under
  `protein_analysis/motifs/`.
- `phylogenetics/`: every phylogenetic analysis, divided by purpose:
  `compact_panel/`, `combined_trees/`, and `per_gene_trees/`.
- `expression/`: expression inputs, scripts, tables, and figures.
- `synteny/`: parsed synteny evidence and prioritized manual review.
- `regulatory/`: bounded five-species promoter extraction and exploratory
  STREME/Tomtom/AME/FIMO analysis.
- `phenotypes/`: pinned MGI/MP and HPO evidence summaries; downloaded raw files
  are gitignored.
- `evolutionary_rates/`: supplementary protein-guided codon alignments and
  human-reference NG86 dN/dS estimates with explicit saturation flags.

## Interpretation entry points

- `overview/README.md`: short navigation page and directory map.
- `overview/RESULTS_INTERPRETATION_GUIDE.md`: what every result and figure means.
- `overview/GENE_FUNCTION_REFERENCE.md`: detailed, cited biology and explicit
  unknowns for BGN, FMOD, PRELP, EPYC, LUM, OGN, and context gene OMD.
- `overview/SIX_GENE_PROFILES.md`: concise bridge from gene function to this
  project's expression and comparative results.
- `overview/WORKFLOW_OVERVIEW.md`: what was done, in biological and computational order.
- `overview/THESIS_CONCLUSION.md`: defensible overall conclusion and claims to avoid.
- `gene_structure/README.md`: where the gene-structure analysis lives and how to read it.

## Provenance rules

- Historical analysis snapshots have been moved to the local-only top-level
  `archive/analysis_history/` directory.
- Directories named `diagnostics` contain QC evidence, not final thesis figures.
- Tool-native IQ-TREE, SignalP, HMMER, ProtSpace, and SynVoy files remain beside
  their workflows for reproducibility.
- Raw databases stay under top-level `data/`; they are not copied here.
- Use `overview/AUTHORITATIVE_RESULT_LOCATIONS.md` rather than historical path
  dumps or filenames found in dated logs.
