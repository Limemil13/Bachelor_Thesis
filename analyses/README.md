# Analysis directory

Each subdirectory contains the scripts and compact outputs for one evidence
layer. Raw downloads and large external databases are excluded from version
control; retained provenance files identify the public sources and processed
inputs used in the analyses.

## Sections

- `expression/`: juvenile skeletal-lineage, growth-plate, and Bgee expression
  summaries.
- `protein_analysis/`: candidate proteins, reciprocal searches, Pfam domains,
  MSAs, SignalP, sequence logos, ProtSpace, and protein motifs.
- `phylogenetics/`: compact per-gene trees, the combined panel tree, reference
  family tree, and iTOL annotations.
- `synteny/`: SynVoy summaries, locus audits, review decisions, and diagnostic
  neighbourhood comparisons.
- `gene_structure/`: representative transcript, CDS-exon, splice-phase, and
  domain-exon analyses.
- `evolutionary_rates/`: protein-guided codon alignments and pairwise coding
  constraint estimates.
- `regulatory/`: exploratory promoter extraction and motif analysis.
- `phenotypes/`: filtered MGI and HPO associations.
- `literature/`: structured reference evidence without copyrighted PDFs.
- `overview/`: integrated tables, workflow documentation, and output map.

## Entry points

- `overview/AUTHORITATIVE_RESULT_LOCATIONS.md`: locations of final tables and
  figures.
- `overview/METHODS_ALGORITHM_NOTES.md`: inputs, settings, scripts, outputs, and
  interpretation limits for each method.
- `overview/WORKFLOW_OVERVIEW.md`: order and purpose of the analysis steps.
- `overview/THESIS_REPRODUCIBILITY_STATUS.md`: current completion and QC status.
- `overview/verify_thesis_analysis_outputs.py`: consistency checks across the
  retained analysis package.

Generated evidence tables should be rebuilt from their source tables and
scripts rather than edited directly.
