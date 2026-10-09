# Comparative analysis of vertebrate SLRP genes

Here you can find the main scripts and processed result files used for my
bachelor thesis on small leucine-rich proteoglycan genes. The comparative panel
consists of BGN, DCN, EPYC, FMOD, LUM, OGN and PRELP. OMD was included in
selected analyses as a context gene.

The analyses cover expression, protein sequences and domains, phylogeny,
microsynteny, gene structure, coding-sequence constraint, promoter motifs and
phenotype annotations. Suspicious candidates are kept in the evidence tables so
that exclusions and uncertain assignments can be traced.

## Repository structure

- `analyses/expression/`: mouse expression analyses.
- `analyses/protein_analysis/`: candidate curation, sequence, domain, SignalP,
  motif and protein-space analyses.
- `analyses/phylogenetics/`: compact gene trees and focused reference trees.
- `analyses/synteny/`: SynVoy summaries and manual locus checks.
- `analyses/gene_structure/`: coding structure and domain-exon comparisons.
- `analyses/evolutionary_rates/`: exploratory pairwise coding constraint.
- `analyses/regulatory/`: exploratory promoter-sequence analyses.
- `analyses/phenotypes/`: curated mouse and human phenotype summaries.
- `analyses/overview/`: final cross-analysis evidence tables and summary figures.

The main outputs and the role of each analysis are listed in
[`analyses/overview/README.md`](analyses/overview/README.md).
