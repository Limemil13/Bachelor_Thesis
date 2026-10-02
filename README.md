# SLRP bachelor-thesis and synteny workspace

This repository contains the SLRP bachelor-thesis analyses together with the
Python synteny pipeline from which the work developed. The authoritative thesis
results are under `analyses/`; older SFTP/HRNR outputs and RNA heatmap scripts
are preserved under `archive/legacy_sftp_hrnr_pipeline/` and are not the source
of truth for the current seven-gene comparative panel.

For a quick project map, start with `analyses/overview/README.md`. It links the
authoritative result locations, interpretation guide, gene profiles, methods,
and workflow documentation.

## Current thesis layout

- `analyses/overview/`: panel decision, reproducibility status, methods,
  integrated tables, and the authoritative result map.
- `analyses/protein_analysis/`: 103 canonical proteins, domains, MSA,
  conservation, SignalP, ProtSpace, mature-protein motifs, and the manual-review
  packet.
- `analyses/phylogenetics/`: all tree work, separated into compact per-gene
  trees, combined/family trees, and large per-gene NCBI trees.
- `analyses/evolutionary_rates/`: supplementary protein-guided codon
  alignments and bounded pairwise dN/dS estimates.
- `analyses/expression/`: NCBI GEO GSE305415 P0, Bgee, and direct mouse/rat
  growth-plate expression analyses.
- `analyses/gene_structure/`: representative five-species coding-gene
  structure tables, provenance, and figures.
- `analyses/synteny/`: parsed SynVoy evidence, complete 98-row P1/P2/P3 review,
  exact-GFF locus diagnostics, and the rebuilt outcome figure. Raw SynVoy
  reports remain in the WSL SynVoy project; the current opossum locus layer is
  excluded until its mismatched FASTA/GFF pair is rerun.
- `analyses/regulatory/`: bounded five-species promoter extraction and
  exploratory STREME/Tomtom/AME/FIMO results.
- `analyses/phenotypes/`: current MGI/MP and pinned HPO phenotype summaries;
  raw downloadable database files remain gitignored.
- `data/`: large raw databases, genomes, queries, and Pfam assets; do not treat
  this as a results folder.
- `runs/`: retained compact reciprocal-BLAST provenance used to build candidate
  panels.
- `archive/`: historical/diagnostic outputs, the legacy SFTP/HRNR project, and
  dated cleanup records. This directory remains local and is not pushed.
- `thesis/`: active, Overleaf-ready LaTeX thesis source. Start with
  `thesis/README.md`; the official university front-page template is already
  integrated in `main.tex`.

The large `data/`, `runs/`, `archive/`, virtual-environment, and temporary
directories are local working material. GitHub contains code, small curated
inputs, authoritative tables/figures, documentation, and thesis text—not
downloadable genomes, databases, private papers, or workflow caches.

Start with `analyses/overview/RESULTS_INTERPRETATION_GUIDE.md`,
`analyses/overview/WORKFLOW_OVERVIEW.md`, and
`analyses/overview/AUTHORITATIVE_RESULT_LOCATIONS.md`.

## Legacy synteny pipeline

The Python pipeline identifies a focal gene and compares gene neighborhoods
across vertebrate genomes. It was originally used for an SFTP-focused project;
those old outputs are preserved separately from the current SLRP thesis.

### Overview

The pipeline performs:

- Identification of a focal gene in a reference genome
- Extraction of neighboring genes within a defined genomic window
- Mapping of these neighbors to target genomes using BLAST (BLASTP, TBLASTN, BLASTN fallback)
- Optional reciprocal validation of hits
- Scoring of conserved gene neighborhoods
- Generation of comparative tables and visualizations

The approach allows robust comparison of local gene order (synteny) across species.

### Features

- Flexible CLI interface (Typer-based)
- Supports multiple genomes and all-vs-all comparisons
- Automatic handling of missing annotations via sequence-based fallback
- Reciprocal BLAST validation for increased reliability
- Output as JSON, TSV, and visualization plots

### Local setup

Install the Python package and development tools from the repository root:

```bash
python -m pip install -e ".[dev]"
```

The mapping commands also require the relevant external executables on
`PATH`, including NCBI BLAST+ and, for genome downloads, NCBI Datasets.

Before committing changes, run:

```bash
ruff format --check analyses synteny_pipeline
ruff check analyses synteny_pipeline
python analyses/overview/verify_thesis_analysis_outputs.py
```
