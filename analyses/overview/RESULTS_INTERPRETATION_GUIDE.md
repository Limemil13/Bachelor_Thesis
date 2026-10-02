
# Results and figure interpretation guide

Last audited: 2026-09-03

This is the entry point for understanding the results. It separates biological
results from quality-control observations and records the claim that each figure
can support. The analyses did not all use the same species set: the compact
protein panel contains 103 proteins (13-16 per gene), gene structure uses five
reference species, completed SynVoy covers 14 target genomes for seven genes, and
direct growth-plate expression covers mouse and rat. Therefore, do not write
that every result was obtained in "the same 15 species."

## Integrated result

The seven selected growth-plate-relevant SLRPs show strong conservation of their
secreted leucine-rich-repeat protein core, coding-exon organization, and
gene-level phylogenetic identity across representative vertebrates. This is
evidence of evolutionary constraint, not an absence of a result. The variation
is concentrated in expression context, intron length, lineage-specific
duplication/paralogy, and uncertain deep-lineage annotations.

The evidence does not justify the stronger claim that growth-plate function is
conserved in every sampled species: expression was not measured in matched
growth plates across the protein panel. Expression identifies a biologically
relevant candidate set; sequence, structure, trees, and synteny then test its
evolutionary conservation and annotation quality.

Open the integrated dashboard first:

- Figure: `figures/six_gene_evidence_dashboard.png` (editable SVG beside it)
- Values: `tables/six_gene_evidence_dashboard_values.tsv`
- Interpretation: BGN has the strongest expression anchor, while DCN supplies
  the necessary class-I paralog control. FMOD and PRELP
  are also strong anchors. EPYC and LUM remain justified main genes, with
  specific sequence/synteny caveats. OGN is the most variable and therefore the
  most informative exploratory gene, not a failed result.

## Expression

### Main comparative-panel expression figure

- Figure: `../expression/figures/six_gene_expression_evidence.png`
- Direct-data table:
  `../expression/gse114919/tables/gse114919_slrp_tibia_cross_condition_summary.tsv`
- GSE305415 WT P0 data table (stored Salmon quantifications):
  `../expression/tables/mouse_expression_descriptive_metadata_corrected.tsv`

The GSE305415 WT P0 Prx1-lineage dataset detects all seven candidates. Its mean
TPM ranking is BGN, OGN, FMOD, DCN, LUM, PRELP, EPYC. These are cells isolated from
digested juvenile knee joints, not laser-captured growth-plate zones. This
panel therefore supports expression in a relevant juvenile skeletal lineage,
not tissue specificity.

The independent GSE114919 panel is the anatomically direct validation: it uses
laser-captured proliferative and hypertrophic zones from juvenile mouse and rat
growth plates. BGN, FMOD, and PRELP are the strongest and most consistent genes;
EPYC and LUM are mid-ranked; OGN is present in mouse but weak or partly absent in
rat. Values must be interpreted within each experiment and species rather than
pooled with the local TPMs.

### Age contrasts and Bgee overview

- Figure: `../expression/figures/expression_age_and_bgee_overview.png`
- Age contrasts:
  `../expression/gse114919/tables/gse114919_slrp_tibia_age_contrasts.tsv`
- Bgee summaries: `../expression/tables/bgee_slrp_relevant_expression_summary.tsv`

Panel A shows four-week minus one-week changes on the authors' processed
expression scale. BGN, FMOD, and PRELP increase in both zones of mouse and rat.
EPYC decreases strongly in rat; LUM decreases in rat PZ and is almost unchanged
in rat HZ; OGN is mixed. These are descriptive contrasts, not a formal
differential-expression model, so they indicate possible developmental
regulation without establishing significance.

Panel B counts matching curated Bgee calls. Mouse has much better coverage than
zebrafish or chicken. A zero is a lack of a matching database call, not evidence
that a gene is biologically absent or unexpressed.

### Legacy descriptive TPM heatmaps

- `../expression/figures/mouse_expression_descriptive_raw_TPM.png`
- `../expression/figures/mouse_expression_descriptive_log10_TPM.png`

The raw panel emphasizes high-abundance genes; the log panel makes lower values
visible. Both mix P0 Prx1-lineage samples with unrelated aging-atlas tissues,
ages, and studies. Keep them as diagnostic/descriptive figures or in an
appendix. Do not call them differential expression or a clean tissue-specificity
comparison.

### Biological decision

Juvenile RNA is appropriate for a growth-plate thesis because the growth plate
is developmentally active in young animals. The limitation is anatomical and
cell-state specificity, not youth. Candidate selection should therefore combine
direct growth-plate expression, known skeletal biology, protein identity,
domain architecture, phylogeny, and synteny; RNA abundance alone is not a sound
ortholog or thesis-candidate criterion.

## Protein conservation, domains, MSA, and SignalP

### Source-of-truth tables

- Protein-level integrated table:
  `../protein_analysis/tables/protein_conservation_domain_msa_signalp_summary.tsv`
- Gene-level integrated table:
  `../protein_analysis/tables/protein_conservation_domain_msa_signalp_gene_summary.tsv`
- MSA summary: `../protein_analysis/tables/msa_conservation_gene_summary.tsv`

The curated panel contains BGN 13, DCN 15, EPYC 15, FMOD 14, LUM 16, OGN 14,
and PRELP 16 proteins. Chicken `XP_414298.2` is excluded from protein analyses because the
BGN locus annotation conflicts with an ASPN-like product and phylogenetic
placement. The lamprey class-II sequence is retained only once, under LUM.

### Protein-conservation overview

- Figure: `../protein_analysis/figures/protein_conservation_overview.png`

Mean within-gene MSA identities are BGN 72.14%, DCN 71.64%, PRELP 66.20%, FMOD
64.96%, EPYC 62.17%, LUM 60.21%, and OGN 58.04%. OGN remains the most sequence-variable
panel member, while every gene retains a substantial conserved core. Removing
the amphioxus SLRP-like sequence raised OGN alignment occupancy to 89.20% and
its fully conserved ungapped columns from 43 to 79, showing that most of the
former gap-rich structure was driven by the questionable deep-lineage record.

The percentages are summaries across related species, not independent replicate
measurements. They support constraint on homologous protein structure but do not
measure tissue specificity or phenotype.

### Conservation profiles and sequence logos

- Profiles:
  `../protein_analysis/alignments/canonical/conservation/<GENE>_conservation_profile.png`
- Logos: `../protein_analysis/figures/sequence_logos/<GENE>_sequence_logo.png`

The seven profiles show where occupancy and residue conservation rise or fall
along each alignment. The logos show residue frequencies with stack height
weighted by information content and non-gap occupancy. BGN has a dense,
continuous conserved core. OGN has more gap-heavy and low-occupancy segments,
consistent with its lower mean identity and deep-lineage uncertainty. The other
four genes are intermediate and retain clear repeated conserved blocks expected
for LRR-containing SLRPs. Logo coordinates are MAFFT columns, not human residue
numbers; exact domain/residue claims require mapping to a reference sequence.

### Pfam domains

All retained proteins have significant SLRP/LRR-family Pfam evidence except one
OGN sequence (13/14 positive). The weak/negative OGN domain case is zebrafish.
The excluded amphioxus sequence has strong SLRP-like LRR evidence but poor
OGN-specific orthology evidence, illustrating that Pfam supports family
membership but cannot by itself resolve close paralogs.

### Signal peptides

Retained slow-sequential SignalP results exist for all 103 current proteins:
100 are positive and three are negative. The negative records are spotted-gar
OGN `XP_006631165.2`, fragmentary whale-shark BGN `XP_048475905.1`, and the
N-terminally incomplete opossum DCN `XP_001363160.3`. The spotted-gar call should trigger a
manual N-terminal completeness check; it is not proof that the biological
protein is non-secreted. The numerous PNGs under
`../protein_analysis/signal_peptides/output/better_run/` and `slower_run/` are
tool-native per-sequence diagnostic plots from historical runs. Use the merged
canonical table, not the presence of one old plot, as the final result. The
historical plots also contain the now-excluded chicken BGN/ASPN record and the
old duplicated lamprey record.

### Five-sequence review queue

The five-protein packet is under `../protein_analysis/manual_review/`. OGN
opossum and zebrafish were retained, spotted gar was retained as tentative, and
amphioxus was excluded from the confident OGN set and preserved only as
uncertain SLRP-like provenance. EPYC whale shark was retained as tentative:
its conserved middle-to-C-terminal region supports EPYC, while earlier
mismatches, reduced coverage, and an alternative reciprocal hit remain visible
as caveats.

## Supplementary motif and phenotype evidence

- Protein MEME/MAST assigns all 103 candidates most strongly to their current
  gene's motif model. This supports the curated assignments, including tentative
  whale-shark EPYC and spotted-gar OGN, but largely recapitulates the conserved
  LRR/cysteine-rich sequence blocks already visible in the MSA.
- MGI phenotype annotations give direct skeletal/cartilage/joint terms for BGN,
  DCN, FMOD, and EPYC, broader ECM phenotypes for LUM and OGN, and no skeletal-category
  PRELP term in the current mouse release. HPO supplies direct human disease
  annotations for BGN and DCN. Counts reflect database depth, not biological rank.
- Comparative promoter analysis has a negative corrected result: no targeted
  motif reaches AME panel-corrected E < 0.05 and strict FIMO reports no q <= 0.05
  site. The exploratory recurrence figure is hypothesis-generating only.

## ProtSpace

- Figure: `../protein_analysis/protspace/figures/protspace_canonical_projections.png`
- Quantitative summary:
  `../protein_analysis/protspace/tables/protspace_canonical_gene_summary.tsv`
- Per-protein review:
  `../protein_analysis/protspace/tables/protspace_canonical_embedding_review.tsv`

PCA, UMAP, and t-SNE show the same embeddings in different two-dimensional
projections. Visible proximity is exploratory because each projection distorts
some distances. The more defensible leave-one-out centroid results assign
BGN 10/13, DCN 15/15, EPYC 14/15, FMOD 12/14, LUM 15/16, OGN 14/14, and PRELP 16/16 to
their own gene centroid. Exceptions include deep/partial BGN models tending
toward DCN, whale-shark EPYC toward OGN, catshark/whale-shark FMOD toward LUM,
and lamprey LUM toward FMOD. These patterns agree with the focused manual-review questions but
do not independently prove orthology.

## Phylogenies

### Compact per-gene trees

- Figure: `../phylogenetics/figures/compact_per_gene_trees.png`
- Tool outputs: `../phylogenetics/compact_panel/trees/<GENE>/`

These are the most interpretable gene-level trees because they use the curated
canonical sequences. BGN, DCN, EPYC, FMOD, LUM, OGN, and PRELP contain 13, 15,
15, 14, 16, 14, and 16 tips, respectively. The displayed Newick root is arbitrary;
interpret supported splits and candidate placement, not a left-to-right
ancestral direction. The FMOD tree was corrected to remove the lamprey sequence
that had been duplicated between FMOD and LUM.

### Combined seven-gene tree

- Figure: `../phylogenetics/figures/six_gene_combined_tree.png`
- Review table:
  `../phylogenetics/combined_trees/six_gene_tree/compact_tree_gene_clade_review.tsv`
- Flagged-neighbor table:
  `../phylogenetics/combined_trees/six_gene_tree/compact_tree_flagged_tip_neighbors.tsv`

DCN, EPYC, FMOD, OGN, and PRELP form pure unrooted gene splits. The largest LUM
split contains 15/16 LUM sequences, with the lamprey LUM-like sequence among
FMOD-like proteins. The revised vertebrate OGN set forms a 14/14 pure split
(97.7/100 label). BGN does not form one pure split, but all BGN and DCN tips have
a closer same-gene neighbour in a focused class-I analysis. Thus, most groups
are strongly separated while the BGN and deep class-II patterns remain qualified. This is a stronger and more
interesting result than a blanket statement that "everything is conserved."

### Phylogenetic result overview and large trees

- Figure: `../phylogenetics/figures/phylogenetic_result_overview.png`
- Inventory: `tables/phylogeny_result_inventory.tsv`

The downloaded NCBI datasets contain hundreds to more than one thousand records
per gene. The existing large trees retain 521-781 tips per gene, but their old
filter parsed the last bracketed NCBI header field as a species name. Fields such
as `isoform=X1` could therefore collapse unrelated loci. Those tree files are
retained for exploration and provenance only. Corrected locus-filtered candidate
sets are ready (BGN 639, FMOD 546, PRELP 848, EPYC 632, LUM 838, OGN 855), but
the six long IQ-TREE reruns are pending. Do not cite topology or tip counts from
the legacy large trees as final evidence. The curated compact trees already
support the central thesis conclusion.

## Species-lineage comparison

- Figure: `figures/species_lineage_evidence_summary.png`
- Table: `tables/species_lineage_evidence_summary.tsv`

Panel A summarizes median identity to the human reference across the canonical
proteins available for each species. The broad gradient from lamprey/sharks to
mammals is compatible with evolutionary divergence, but it is not monotonic:
spotted gar (68.1%) exceeds zebrafish (61.2%) and coelacanth (59.8%). Different
genes, evolutionary rates, and annotation quality contribute to these values,
and lamprey has only three retained genes.

Panel B summarizes the completed synteny decisions. Spotted gar and coelacanth
each have five accepted loci, while every shark has three accepted loci plus
tentative or ambiguous records. This supports recognizable early
jawed-vertebrate loci while preserving annotation caveats. Opossum is the
critical technical control: its proteins are relatively close to human, but
all synteny rows are unusable because its FASTA/GFF versions mismatch.
Amphioxus defines an uncertainty boundary, not demonstrated absence.

The useful species-level finding is that biological divergence, duplication,
paralog confusion, model quality, and input compatibility can be distinguished
only by combining evidence layers. Do not rank species by one median or infer a
linear evolutionary progression from the display order.

## Gene structure

- Overview: `../gene_structure/figures/gene_structure_conservation_overview.png`
- Full transcript structures:
  `../gene_structure/extended/figures/full_gene_structure_5species.png`
- Coding-exon/splice-phase map:
  `../gene_structure/extended/figures/coding_exon_splice_phase_map.png`
- Human Pfam/domain-to-exon map:
  `../gene_structure/extended/figures/human_domain_exon_architecture.png`
- Per-gene figures:
  `../gene_structure/figures/<GENE>_cds_gene_structure_5species.png`
- Main table:
  `../gene_structure/tables/gene_structure_representative_5species.tsv`
- Full explanation: `../gene_structure/README.md`

Across human, mouse, cow, chicken, and zebrafish, the representative coding
structure is identical in exon count within every tested gene: BGN 7; EPYC and
OGN 6; FMOD, LUM, OMD, and PRELP 2. Protein lengths are highly stable (all
coefficients of variation below 2%), whereas genomic spans vary about two- to
four-fold. This means that coding organization and protein size are strongly
conserved while most locus-size variation lies in introns.

Each per-gene figure is transcript-oriented: boxes are CDS segments and lines
are introns. Their x-axis is relative genomic position, not a cross-gene common
scale. OMD is included as a class/context gene but is not part of the main seven.
Chicken BGN remains a locus/product annotation conflict; a conserved BGN-locus
structure does not prove that the disputed translated model is a valid BGN
protein.

The extension reconstructs all 40 CDSs from the corresponding genome assembly.
Every translation is complete by the applied QC checks, all coding blocks
pass GFF phase continuity, and all 26 comparable coding junctions retain their
intron phase across the five species. This is stronger than invariant exon
count: it shows a homologous frame-compatible coding pattern even where a
boundary shifts by several residues. Eleven loci have annotated coding alternatives
that could alter exon count or protein length, so the selected representative
must remain explicit. UTR lengths are descriptive and annotation-sensitive.

## Supplementary coding constraint

- Figure: `../evolutionary_rates/figures/pairwise_dn_ds_human_reference.png`
- Pairwise table: `../evolutionary_rates/tables/pairwise_dn_ds_human_reference.tsv`
- Interpretation and limits: `../evolutionary_rates/README.md`

Twenty-five finite human-reference NG86 estimates are below one, which is
compatible with predominant purifying constraint. Six human--zebrafish
comparisons are unavailable because dS is non-positive or undefined, and three
human--chicken comparisons have high-dS saturation warnings. Use this as a
supplementary nucleotide-level consistency check, not as a formal test of
selection on every lineage or residue.

## Synteny

- Figure: `../synteny/figures/synvoy_species_confidence_heatmap.png`
- Evidence table: `../synteny/tables/synvoy_gene_species_evidence.tsv`
- Completed decisions: `../synteny/diagnostics/priority_locus_reviews/all_synvoy_review_decisions.tsv`
- Status and interpretation: `../synteny/SYNTENY_ANALYSIS_STATUS.md`

The upper heatmap summarizes the best retained SynVoy confidence per target
genome, not the number of verified orthologs. High-confidence best calls occur
in 9/14 BGN, 9/14 EPYC, 11/14 FMOD, 10/14 LUM, 13/14 OGN, and 11/14 PRELP target
genomes. The lower heatmap shows the completed integrated locus decision.
Across all 98 rows, 62 are accepted, 20 tentative, 12 ambiguous, and four
rejected. The figure therefore separates automated recovery from the final
multi-evidence interpretation.

All seven panel-wide SynVoy runs and all 26 P1, 68 P2, and four P3 reviews are
complete. Sixty-seven canonical protein accessions directly overlap the best
post-ownership locus. Detailed exceptions show that alternative retained loci
can outperform the summary's highest-ranked fragment, as in spotted-gar FMOD,
and that deep-lineage hits can belong to unrelated genes, as in amphioxus EPYC
and FMOD. All opossum synteny rows are excluded pending a rerun with a matched
FASTA/GFF pair. Failure to recover a call can result from annotation, assembly,
naming, contig, or search-window limitations and must not be written as certain
gene loss without independent evidence.

## What belongs in the thesis

Recommended main figures:

1. `figures/thesis_workflow_overview.png`
2. `../expression/figures/six_gene_expression_evidence.png`
3. `../protein_analysis/figures/protein_conservation_overview.png`
4. `../phylogenetics/figures/six_gene_combined_tree.png`
5. `../gene_structure/figures/gene_structure_conservation_overview.png`
6. `figures/six_gene_evidence_dashboard.png` as the integrated summary

Recommended supplementary/diagnostic material includes the six sequence logos,
six per-gene trees, seven per-gene structure plots, ProtSpace projections,
SynVoy heatmap and review tables, Bgee/age overview, and individual SignalP
diagnostics. The two mixed-study TPM heatmaps and all legacy large-tree
topologies should not be central result figures.
