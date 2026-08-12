# Bachelor Thesis: Expression-guided evolutionary characterization of growth-plate-associated SLRPs

This repository contains scripts, curated result tables, and summary figures for a bachelor thesis on small leucine-rich proteoglycans (SLRPs) with relevance to cartilage and growth plate biology.

## Project idea

The thesis combines mouse growth plate RNA-seq expression analysis with comparative genomics to identify and characterize SLRP genes associated with cartilage and growth plate biology.

Main workflow:

- Mouse growth plate RNA-seq
- SLRP expression screening
- Candidate gene prioritization
- BLAST / reciprocal BLAST / length filtering
- Domain and protein architecture analysis
- Multiple sequence alignment and phylogeny
- Synteny analysis with SynVoy
- Final orthology/confidence summary

## Main candidate genes

Current main candidate genes:

- EPYC
- FMOD
- BGN
- OGN
- PRELP

Additional candidate / comparison genes:

- OMD
- DCN
- LUM
- ASPN

## Repository structure

- `data/accessions/`: accession lists
- `data/metadata/`: gene panels and metadata
- `data/processed/`: small curated input files
- `scripts/rnaseq/`: RNA-seq expression scripts
- `scripts/phylogeny/`: alignment, tree, and iTOL helper scripts
- `results/rnaseq/`: expression tables
- `results/expression_figures/`: expression heatmaps
- `results/domain_tables/`: domain and MSA quality tables
- `results/phylogeny/`: alignments, trimmed alignments, and iTOL files
- `results/protspace/`: protein embedding inputs and summaries

## Current status

Completed or partly completed:

- Mouse growth plate RNA-seq quantification with Salmon
- Growth plate vs comparison tissue expression table
- Expression heatmaps
- Candidate prioritization based on growth plate expression
- BLASTP / reciprocal BLASTP / length filtering for selected SLRPs
- Pfam/HMMER domain analysis
- MAFFT alignments and preliminary phylogenetic trees
- iTOL label/group files
- ProtSpace protein embedding analysis
- SynVoy synteny analysis for several selected genes

Main tasks still in progress:

- Complete/debug SynVoy for OGN
- Consistently trim and rerun final phylogenies
- Add signal peptide analysis
- Optionally add MEME motif analysis
- Integrate expression, domains, signal peptides, MSA, phylogeny, and synteny into a final confidence table

## Data availability

Large raw data files are not stored in this repository. RNA-seq SRA accessions and gene metadata are listed under `data/accessions/` and `data/metadata/`.

Raw sequencing data, transcriptomes, proteomes, genomes, Salmon indices, Pfam databases, and SynVoy work directories should be downloaded or generated locally.
