# Comparative promoter-motif analysis

## Thesis question and scope

This exploratory module asks whether the seven selected SLRP loci contain
recurrent promoter-sequence patterns or known transcription-factor motifs that
could generate testable regulatory hypotheses. It does **not** test enhancer
activity, transcription-factor occupancy, chromatin accessibility, or causal
control of growth-plate expression.

The analysis uses one representative transcript for each of BGN, DCN, FMOD,
PRELP, EPYC, LUM, and OGN in human, mouse, cow, chicken, and zebrafish (35 loci).
The
same transcripts and NCBI genome assemblies used by the representative
gene-structure analysis were retained. Two strand-aware windows were extracted
relative to the annotated transcription start site:

- core promoter: -500 to +100 bp;
- proximal promoter: -2,000 to +200 bp.

All 70 sequences had their expected length, none was clipped at a contig edge,
and none contained ambiguous bases. Proximal-window GC content ranged from
32.48% to 59.06%, which is important because composition can influence motif
discovery and scanning.

## Analyses

1. STREME 5.5.8 discovered motifs in all 35 proximal promoters relative to
   dinucleotide-shuffled controls. Separate gene-level searches compared the
   five promoters of one gene against the 30 promoters of the other genes.
2. Tomtom compared discovered motifs with the 1,019 nonredundant vertebrate
   CORE profiles in JASPAR 2026.
3. AME tested JASPAR motif enrichment relative to dinucleotide-shuffled
   promoters. A second 40-profile panel focused on SOX5/6/9, RUNX2/3, SMAD,
   HIF1A, GLI, AP-1, NF-kB, WNT/TCF, MEF2C, STAT3, CREB1, ATF4, FOXA2, and
   NFATC1 motifs.
4. FIMO scanned the 35 promoters with the targeted panel. The strict scan used
   a q-value threshold of 0.05. A separate p <= 1e-4 scan was retained only to
   summarize exploratory recurrence across species.

## Main result

The strict targeted FIMO scan returned zero q <= 0.05 sites. No targeted
cartilage/growth-plate motif reached AME panel-corrected E < 0.05. HIF1A was the
top targeted AME result (raw p = 7.49e-5; adjusted p = 0.0638), but its
panel-corrected E-value was 2.55 and it is therefore not a significant
panel-wide enrichment.

Twelve all-JASPAR profiles reached AME E < 0.05. Several corresponded to A-rich,
T-rich, or repetitive patterns, and the global STREME motifs were likewise
dominated by low-complexity sequence. These results are best interpreted as
promoter-composition patterns rather than evidence for a shared growth-plate
regulatory mechanism.

The gene-level STREME comparisons had only five positive sequences each and
could not create holdout sets. Their five reported motifs per gene are
training-only discoveries. Tomtom returned reference-profile matches for two
global motifs, one EPYC motif, and two LUM motifs at q <= 0.10, but none provides
validated regulatory evidence.

The exploratory FIMO recurrence heatmap can be used in supplementary material
to identify hypotheses for later ChIP-seq/ATAC-seq or reporter testing. Its
uncorrected motif occurrences must not be described as transcription-factor
binding sites.

## Authoritative outputs

- `tables/promoter_sequence_qc.tsv`
- `tables/promoter_denovo_motif_summary.tsv`
- `tables/promoter_ame_all_jaspar_significant.tsv`
- `tables/promoter_ame_targeted_tf_summary.tsv`
- `tables/promoter_targeted_tf_strict_scan_summary.tsv`
- `tables/promoter_targeted_tf_exploratory_recurrence.tsv`
- `figures/promoter_proximal_gc_heatmap.png`
- `figures/promoter_targeted_tf_ame_enrichment.png`
- `figures/promoter_targeted_tf_species_recurrence.png`

The seven-gene raw MEME Suite reports are under
`raw_outputs_7gene_dcn_20260902/`; the earlier six-gene run remains in
`raw_outputs/` for provenance. Extracted sequences and metadata are under
`inputs/`. The analysis is reproduced by
`scripts/run_regulatory_motif_analysis.sh` followed by
`scripts/summarize_regulatory_motifs.py`.

## Recommended thesis placement

Keep this module as a short exploratory Results paragraph and a supplementary
table/figure. The defensible conclusion is that promoter motif analysis did not
identify a statistically supported shared cartilage-regulatory signature in
this small five-species-per-gene design. This negative result does not weaken
the protein-conservation story and should not be used to rerank the seven genes.
