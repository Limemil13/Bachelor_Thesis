# Bachelor thesis reproducibility status

Last audited: 2026-10-01

## Current status

The local comparative package has been regenerated around a clean **103-protein,
seven-gene panel**: BGN 13, DCN 15, EPYC 15, FMOD 14, LUM 16, OGN 14, and
PRELP 16. OMD remains a context gene. The integrity checker passes for the
protein, expression, phylogeny, structure, evolutionary-rate, promoter,
phenotype, and current SynVoy outputs.

SignalP is now complete for all 103 current proteins, including the final 17
sequences (15 DCN plus two replacement BGN models). The updated-tool DCN
SynVoy run completed at `results/dcn_human_15species_dev_20260903` and is parsed
into the current evidence package. Updated BGN also completed at
`results/bgn_human_15species_dev_20260905` and is integrated. The updated FMOD
report at `results/fmod_human_15species_dev_20260906` has been compared with
the historical report, manually reviewed, and integrated into all 14 FMOD rows
of the 98-row table.
PRELP started at `results/prelp_human_15species_dev_20260914` and has a
resumable iterative-search checkpoint after wave 4 of 7, but no final report.
An updated-tool EPYC run has not been launched. No SynVoy process was active
during this audit. Historical PRELP and EPYC reports remain the sources for
their rows; their updated-tool reruns are optional sensitivity analyses.
The staged procedure is documented in
`analyses/synteny/SYNVOY_UPDATED_RERUN_RUNBOOK.md`.

## Authoritative protein summary

The one-row-per-protein conservation/domain/MSA/SignalP table belongs at:

`analyses/protein_analysis/tables/protein_conservation_domain_msa_signalp_summary.tsv`

It contains 103 rows and currently has no pending SignalP calls. The
one-row-per-gene aggregate is adjacent in
`protein_conservation_domain_msa_signalp_gene_summary.tsv`.

## Verified local layers

### Protein, domains, MSA, and secretion

- Canonical manifest: 103 unique gene/species/accession records.
- Pfam: 693 significant domain hits and 103 per-protein summaries; all 103 are
  SLRP-like by the current domain rules except one OGN record that remains a
  known manual/automated caveat.
- MAFFT: seven canonical alignments with conservation profiles, identity
  matrices, consensuses, and sequence logos.
- Integrated QC: 84 Supported, 13 Watch, and 6 Needs inspection.
- SignalP calls: 100 positive and three negative. The negative records are
  spotted-gar OGN, fragmentary whale-shark BGN, and N-terminally incomplete
  opossum DCN. The final 17-sequence input and raw job bundle are retained for
  reproducibility.

### Protein embeddings and motifs

- ProtSpace input: 103 canonical SLRPs plus four controls; embeddings are
  107 x 1,024 and the current dated bundle is retained.
- ProtSpace canonical QC: 85 Supported, 1 Watch, 17 Needs inspection. These are
  exploratory outlier flags, not rejection decisions.
- MEME/STREME/MAST: 116 MEME motifs, 35 enriched motifs, 721 protein/model
  comparisons, and 103/103 proteins with their own gene model ranked best.
  The 2026-09-14 rerun used all final SignalP calls: 100 positive signal
  peptides were removed and the three SignalP-negative proteins were untrimmed.
  Because gene models were learned from the same small curated sets and no
  independent hold-out was available, this is supportive classification rather
  than a stand-alone test of orthology or biochemical function.

### Phylogeny

- Seven compact per-gene trees and a combined 103-tip tree were rebuilt with
  ModelFinder and 1,000 SH-aLRT/1,000 ultrafast-bootstrap replicates.
- The combined alignment has 172 trimmed columns and uses LG+I+G4.
- DCN, EPYC, FMOD, OGN, and PRELP are unrooted-monophyletic.
- LUM has a 15/16 pure main split; lamprey remains outside and tentative.
- BGN is not one pure split after inclusion of DCN. A direct class-I diagnostic
  nevertheless finds a closer BGN neighbour for every BGN tip and a closer DCN
  neighbour for every DCN tip. Deep/partial BGN models remain manual-review
  priorities.
- The curated compact trees are sufficient for the central bachelor-thesis
  result. Corrected hundreds-of-sequences trees are optional supplementary work.

### Gene structure and coding constraint

- Five representative species x eight genes = 40 transcripts.
- All 40 reconstructed CDSs pass start/stop/divisibility/internal-stop QC.
- Twenty-six comparable splice junctions conserve intron phase.
- All 40 translations have LRR-family Pfam support; 294 significant hits were
  mapped to coding exons.
- Thirty-two human-target NG86 comparisons produce 25 finite estimates, all
  below one. Deep saturated or undefined comparisons remain explicitly limited.

### Expression

- Previously generated Salmon quantifications of public NCBI GEO GSE305415
  WT samples (BioProject PRJNA1305489) were re-extracted and matched the
  curated values; no blind Salmon rerun is needed.
- P0 Prx1-lineage TPMs are a juvenile skeletal-lineage screen, not isolated
  growth-plate tissue and not a tissue-specificity test.
- GSE114919 mouse/rat microdissected growth-plate zones provide the direct
  anatomical validation. Bgee is qualitative and coverage-limited.
- Expression supports candidate prioritization but does not establish orthology
  or gene function by itself.

### Promoters and phenotypes

- Seven-gene promoter analysis: 35 loci, 70 strand-aware windows, 40 de-novo
  motifs. No targeted cartilage/growth-plate TF motif survived panel-wide
  correction. Keep this exploratory result supplementary and do not rerank genes
  from it.
- MGI/HPO were refreshed for all seven comparative genes from pinned local
  releases. DCN adds 24 single-gene MP terms (4 skeletal/cartilage/joint) and
  10 HPO phenotype terms. Counts reflect curation depth, not effect size.

### Synteny

- BGN, DCN, EPYC, FMOD, LUM, OGN, and PRELP reports are parsed into 98
  gene-target rows: 26 P1, 68 P2, and four P3 review rows. Every row now has a
  recorded decision (62 accepted, 20 tentative, 12 ambiguous, four rejected),
  and all pending queues are empty.
- The DCN updated-tool run passed QC for 14/14 genomes and retained 10 HIGH and
  one MEDIUM post-filter candidates. Amphioxus, spotted gar, and zebrafish have
  no retained post-filter DCN candidate; detailed review recovered canonical
  coelacanth, spotted-gar, and zebrafish loci and left amphioxus unresolved.
- Updated BGN passed QC for 14/14 genomes and retained 10 HIGH and one MEDIUM
  annotation. Updated FMOD passed QC for 14/14 genomes and retained 10 HIGH and
  five MEDIUM calls; it is the authoritative FMOD source in the review table.
  PRELP has a partial resumable run; updated EPYC has not been launched.
- SynVoy is locus-prioritization/microsynteny support, not an automatic
  one-to-one orthology truth table. The completed review confirms 67 direct
  canonical-accession/locus overlaps and documents all exceptions. All seven
  opossum rows remain technically ambiguous because the paired FASTA/GFF
  sequence identifiers mismatch; rerun is required before citing opossum
  synteny.

## Reproducible builders

- Protein refresh: `analyses/protein_analysis/scripts/run_canonical_refresh.sh`
- Protein motifs: `analyses/protein_analysis/motifs/scripts/run_protein_motif_analysis.sh`
- Compact tree builders/review: `analyses/phylogenetics/compact_panel/scripts/`
- Integrated tables/figures: `analyses/overview/scripts/`
- Integrity checker: `analyses/overview/verify_thesis_analysis_outputs.py`

Run the integrity checker after any future SignalP or SynVoy-table update.
Generated source-of-truth tables should be
rebuilt, not manually overwritten.
