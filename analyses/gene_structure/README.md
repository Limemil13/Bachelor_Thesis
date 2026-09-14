# Gene-structure analysis

Last audited: 2026-09-02

This directory is the local, authoritative home of the representative
gene-structure analysis. It contains the restored SynVoy/NCBI GFF3-derived
tables, raw/provenance files, a cross-gene overview, and full-transcript figures.

## Where to find the result

- Representative 40 gene-species rows:
  `tables/gene_structure_representative_5species.tsv`
- Per-gene summary: `tables/gene_structure_summary_5species.tsv`
- Individual coding blocks: `tables/gene_structure_cds_blocks_5species.tsv`
- Cross-gene variability summary:
  `tables/gene_structure_extended_summary.tsv`
- Structural-outlier table:
  `tables/gene_structure_structural_outliers_5species.tsv`
- Overview figure: `figures/gene_structure_conservation_overview.png`
- Per-gene figures: `figures/<GENE>_cds_gene_structure_5species.png`
- Original extraction and legacy files: `raw/`
- GFF provenance lists: `provenance/`

## Extended full-transcript package

The higher-resolution extension is under `extended/`. It adds the parts that
the original CDS-block screen intentionally did not test:

- full representative exon, CDS, derived 5' UTR, and derived 3' UTR features:
  `extended/tables/full_transcript_features_5species.tsv`
- reconstructed-CDS translation and completeness QC for all 40 loci:
  `extended/tables/full_gene_structure_qc.tsv`
- coding-exon boundaries, intron phase, and GFF phase checks:
  `extended/tables/coding_exon_splice_phase_5species.tsv`
- one row per homologous splice junction:
  `extended/tables/splice_junction_conservation_summary.tsv`
- alternative-transcript sensitivity:
  `extended/tables/isoform_sensitivity_5species.tsv`
- one-row-per-gene synthesis:
  `extended/tables/extended_gene_structure_gene_summary.tsv`
- complete transcript plot:
  `extended/figures/full_gene_structure_5species.png`
- splice-boundary and intron-phase plot:
  `extended/figures/coding_exon_splice_phase_map.png`
- Pfam-to-exon map and its tables:
  `extended/figures/human_domain_exon_architecture.png` and
  `extended/tables/representative_domain_exon_gene_summary.tsv`

The reproducible entry point is
`scripts/run_extended_gene_structure.sh`. It reconstructs CDS/protein sequences
from indexed genome FASTAs, reruns the local Pfam scan, maps significant
LRR-family hits to coding exons, and refreshes the supplementary pairwise
coding-constraint analysis.

## What was done

Gene, transcript, exon, and CDS parent relationships were streamed from NCBI
GFF3 annotations. For each of BGN, DCN, EPYC, FMOD, LUM, OGN, OMD, and PRELP in
human, mouse, cow, chicken, and zebrafish, one representative protein-coding
transcript was selected. The analysis records strand, transcript accession,
CDS block count and length, estimated protein length, and genomic locus span.

The merged 86-row raw extraction included alternative transcripts.
Representative selection reduced this to 40 comparable gene-species
observations. The
structural-outlier screen tested changes in CDS exon count and conspicuous
protein-length differences; it found no outliers.

The extension returned to the exact selected NCBI transcripts and genome
assemblies. Every reconstructed CDS was divisible by three, began with a
methionine, ended with a terminal stop, and lacked internal stops. All 170
coding blocks passed the annotation-phase continuity checks. The 26 comparable
splice junctions retained the same intron phase in all five species. Exact
boundary positions were closest for BGN (all six junction ranges at 3 aa) and
OMD (2 aa); the larger 5--18-aa ranges in the other genes indicate modest local
boundary drift while the reading-frame relationship remains conserved.

DCN has eight representative exons and seven CDS blocks in all five species.
All six homologous splice junctions conserve intron phase, although their
coding-boundary positions vary by 17 aa across the species panel. This supports
a stable frame-compatible gene organization while showing more local boundary
drift than BGN.

Eleven of the 40 loci have at least one annotated coding isoform that changes
protein length and/or CDS-exon count. This does not invalidate the selected
representatives: it identifies the species for which conclusions depend most on
transcript choice. In particular, the compact selected zebrafish EPYC model is
more compatible with the conserved protein than its very long alternative.

The representative-protein Pfam rescan recovered 294 significant LRR-family
repeat/domain hits and at least one hit in every translation. Domain models
cross a coding-exon boundary in at least one representative of every gene
except OMD. This links the conserved protein architecture to the stable coding
framework without implying that each exon is an independent functional module.

## How to read the figures

The per-gene plots are transcript-oriented. Rectangles represent coding
segments and connecting lines represent intronic distance. Negative-strand genes
are flipped into transcript direction so that species can be compared directly.
The x-axis is a within-gene relative locus scale; box/line distances should not
be compared as absolute coordinates across different gene figures.

The original overview has two distinct messages:

1. CDS exon counts are invariant within each tested gene: BGN and DCN have 7,
   EPYC and OGN have 6, and FMOD/LUM/OMD/PRELP have 2.
2. Protein length changes little, but genomic span changes substantially. Span
   fold ranges are 2.02 for BGN, 2.97 for EPYC, 2.13 for FMOD, 4.15 for LUM,
   3.65 for OGN, 2.24 for OMD, and 2.04 for PRELP.

All original seven-gene protein-length coefficients of variation are below 2%;
DCN likewise has a narrow representative range (355--374 aa). The result
therefore supports conserved coding organization and protein size, with lineage
variation concentrated mostly in introns. OMD is shown as a class/context gene;
DCN is now included in the completed comparative package following the
2026-09-02 panel audit.

The full-transcript plot adds UTRs but should be read conservatively: UTR length
depends heavily on transcript annotation and is much less stable than the CDS.
The splice-phase map is the stronger comparative result because it shows that
matching exon counts generally represent homologous frame-compatible coding
junctions, not only the same number of boxes.

## Important limitations

- The main comparison remains a representative-transcript screen. The isoform
  sensitivity table inventories annotated coding alternatives but does not test
  their expression or biological dominance.
- Human DCN has several short annotated coding isoforms. The analysis uses the
  complete 360-aa representative; the alternative-transcript table records this
  dependence explicitly.
- Five well-spaced reference species provide a strong overview but not the same
  taxonomic sampling as the compact protein panel.
- GFF3 models, especially predicted `XM_` transcripts, can reflect annotation
  rather than biology.
- Chicken BGN is an unresolved locus/product conflict: the locus structure can
  be described, but it does not validate the disputed ASPN-like protein as BGN.
- UTR differences are annotation- and transcript-choice-sensitive and should not
  be presented as demonstrated regulatory evolution.
- Conserved splice phase is strong structural evidence but does not by itself
  prove orthology or unchanged protein function.

## Is more gene-structure analysis needed?

The three highest-value extensions---splice phase, isoform sensitivity, and
Pfam/domain-to-exon mapping---are now complete. No additional gene-structure
analysis is required for the bachelor thesis. Expanding to more species would
be justified only for a specific evolutionary-origin question and would need
equivalent annotation quality.

Promoter, enhancer, methylation, or broad UTR analyses would be new regulatory
projects and are not required to answer the current comparative coding-structure
question.
