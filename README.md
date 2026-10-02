# Comparative analysis of vertebrate SLRP genes

This repository contains the code and processed results used to examine the
conservation of selected small leucine-rich proteoglycan (SLRP) genes across a
broad vertebrate species panel. The main comparative panel contains BGN, DCN,
EPYC, FMOD, LUM, OGN, and PRELP; OMD is included in selected contextual
analyses.

The project combines expression evidence, protein-domain and sequence
conservation, phylogenetic reconstruction, microsynteny, gene structure,
coding-sequence constraint, promoter motifs, and phenotype annotations. Each
evidence layer is retained separately so that uncertain gene models and
conflicting assignments remain visible.

## Repository structure

- `analyses/`: method-specific scripts, curated inputs, tables, figures, and QC
  records.
- `analyses/overview/`: integrated evidence tables, the workflow summary, and a
  map of authoritative outputs.
- `synteny_pipeline/`: Python package for extracting and comparing annotated
  gene neighbourhoods.
- `pyproject.toml`: Python package metadata and development dependencies.

Start with [`analyses/overview/README.md`](analyses/overview/README.md) for the
short result map. Exact input choices, parameters, implementations, and output
locations are recorded in
[`analyses/overview/METHODS_ALGORITHM_NOTES.md`](analyses/overview/METHODS_ALGORITHM_NOTES.md).

## Reproducibility

Large genomes, public database snapshots, raw sequencing data, private papers,
workflow caches, and the manuscript are not stored in this repository. Their
derived, compact results and provenance records are retained where licensing
and file size permit.

Install the Python package and development dependencies from the repository
root:

```bash
python -m pip install -e ".[dev]"
```

Run the repository checks with:

```bash
ruff format --check analyses synteny_pipeline
ruff check analyses synteny_pipeline
python analyses/overview/verify_thesis_analysis_outputs.py
```

External command-line programs used by individual workflows include NCBI
BLAST+, NCBI Datasets, MAFFT, trimAl, IQ-TREE, HMMER, SignalP, MEME Suite, and
SynVoy. Tool versions and analysis-specific settings are documented beside the
corresponding scripts.

## Current limitation

The synchronized FASTA/GFF input audit failed for the opossum SynVoy target.
Opossum protein, alignment, and tree evidence remains usable, but its synteny
rows are excluded until that target is rerun with a matched assembly pair.
