# Comparative analysis of vertebrate SLRP genes

This repository contains the main scripts and processed results used for my
bachelor thesis on small leucine-rich proteoglycan genes. The comparative panel
consists of BGN, DCN, EPYC, FMOD, LUM, OGN and PRELP. OMD was included in
selected analyses as a context gene.

The analyses cover expression, protein sequences and domains, phylogeny,
microsynteny, gene structure, coding-sequence constraint, promoter motifs and
phenotype annotations. Suspicious candidates are kept in the evidence tables so
that exclusions and uncertain assignments can be traced.

The thesis presents the main findings. More detailed tables and
analysis-specific outputs are kept here so that the candidate decisions and
reported numbers can be checked.

## Repository structure

- `analyses/expression/`: mouse expression analyses.
- `analyses/protein_analysis/`: candidate curation, sequence, domain, SignalP,
  motif and protein-space analyses.
- `analyses/phylogenetics/`: compact gene trees, the combined seven-gene tree,
  and the expanded BGN--DCN and FMOD--LUM family trees.
- `analyses/synteny/`: SynVoy summaries and manual locus checks.
- `analyses/gene_structure/`: coding structure and domain-exon comparisons.
- `analyses/evolutionary_rates/`: exploratory pairwise coding constraint.
- `analyses/regulatory/`: exploratory promoter-sequence analyses.
- `analyses/phenotypes/`: curated mouse and human phenotype summaries.
- `analyses/overview/`: final cross-analysis evidence tables and summary figures.

The main outputs and the role of each analysis are listed in
[`analyses/overview/README.md`](analyses/overview/README.md).

## Running the final checks

The repository contains processed results as well as the scripts used to make
them. From the repository root, the overall integrity check can be run with:

```bash
python analyses/overview/verify_thesis_analysis_outputs.py
```

The three final extensions were run in Ubuntu/WSL. Their reproducible entry
points are:

```bash
bash analyses/phylogenetics/targeted_extensions/scripts/run_targeted_extensions.sh all
bash analyses/protein_analysis/motifs/validation/scripts/run_loso_validation.sh
bash analyses/regulatory/cross_species_conservation/scripts/run_cross_species_promoter_conservation.sh
```

The phylogenetic wrapper records the exact IQ-TREE command and settings. The
repository does not include program binaries, restart files, logs or other
temporary outputs.

## Main tool versions

The final focused analyses recorded MAFFT 7.505, IQ-TREE 2.0.7 and MEME Suite
5.5.8. Signal peptides were predicted with SignalP 6.0. The final Python
scripts were checked with Python 3.13, Biopython 1.87, pandas 3.0.1, NumPy
2.4.2 and Matplotlib 3.10.8. Analysis-specific parameters are kept in the
corresponding scripts and IQ-TREE reports.
