# Thesis-ready figure captions

Last updated: 2026-09-23

These captions describe the current figures and explicitly preserve the scope
and limitations needed for defensible interpretation. Replace `Figure X` with
the final thesis numbering.

Frozen copies of the selected main and supplementary figures are integrated
under `thesis/figures/`; LaTeX assigns the final numbering automatically.

## RNA expression figure

**Figure X. Expression evidence used to prioritize the seven-gene SLRP panel.**
(A) Mean expression of BGN, DCN, FMOD, PRELP, EPYC, LUM, and OGN
in three samples of postnatal-day-0-derived Prx1-lineage chondrocyte material from NCBI GEO
GSE305415 (BioProject PRJNA1305489). Bars show mean
TPM and error bars show standard deviation; the y-axis is logarithmic. These
samples represent a relevant juvenile skeletal lineage but are not anatomically
isolated growth-plate zones. Deposited series- and sample-level records conflict
about culture before RNA extraction, so no exact cell state is assigned. (B) Mean within-condition rank in processed
microdissected tibial growth-plate data from mouse and rat (GSE114919), averaged
across the examined age-by-zone conditions; rank 1 denotes the highest of the
nine screened SLRPs. Values are interpreted within each dataset and species and
were not pooled with the TPM measurements. Expression supports biological
prioritization but does not establish function or orthology.

Source figure: `analyses/expression/figures/six_gene_expression_evidence.png`

**Supplementary Figure X. Replicate-level expression in the GSE114919
growth-plate dataset.** Points show the five deposited biological replicates
for each available mouse or rat age--bone--zone condition; black bars denote
condition means. Values are displayed on the authors' published normalized-
value scale without logarithmic transformation. The panel is descriptive:
absolute values are not compared between species and no inferential test is
implied.

Source figure:
`analyses/expression/gse114919/figures/gse114919_slrp_replicates.png`

## Combined phylogenetic tree

**Figure X. Unrooted maximum-likelihood phylogeny of the 103-protein curated
SLRP panel.** The trimmed MAFFT alignment contained BGN (13 proteins), DCN
(15), EPYC (15), FMOD (14), OGN (14), PRELP (16), and LUM (16). IQ-TREE selected
LG+I+G4 using the Bayesian information criterion and estimated branch support
with 1,000 SH-aLRT and 1,000 ultrafast-bootstrap replicates. Tip colours denote
gene assignment. DCN, EPYC, FMOD, OGN, and PRELP form pure gene-specific
unrooted splits; 15 of 16 LUM sequences form the main LUM split, while lamprey
`XP_075930353.1` is retained as a tentative LUM-like deep-lineage assignment.
BGN does not form one pure split in this seven-gene tree, although every BGN and
DCN tip has a closer same-gene tip in a focused class-I diagnostic.
The displayed orientation is arbitrary and does not imply evolutionary
direction.

Source figure: `analyses/phylogenetics/figures/six_gene_combined_tree.png`

## Per-gene phylogenies

**Figure X. Maximum-likelihood phylogenies of the seven curated gene-specific
protein sets.** BGN, DCN, EPYC, FMOD, OGN, PRELP, and LUM contain 13, 15, 15,
14, 14, 16, and 16 nonredundant proteins, respectively. Alignments were produced with
MAFFT; IQ-TREE ModelFinder selected the substitution model independently for
each gene, and support was estimated with 1,000 SH-aLRT and 1,000
ultrafast-bootstrap replicates. Branch lengths represent substitutions per site.
Trees are unrooted; topology and support are used to assess within-gene sequence
placement, not to infer ancestral direction. The amphioxus SLRP-like sequence is
absent from the final OGN tree following manual exclusion from the confident
ortholog set.

Source figure: `analyses/phylogenetics/figures/compact_per_gene_trees.png`

## Species-lineage evidence summary

**Figure X. Species-lineage summary of curated SLRP evidence.** (A) Median
amino-acid identity to the human reference across the canonical proteins
retained for each species. Labels report the number of proteins; colours
identify lineage groups. The values are descriptive because some species
contribute fewer than seven genes. (B) Accepted, tentative, ambiguous, and
rejected decisions from the completed seven-gene synteny review. Lamprey and
rat were not SynVoy targets, and human was the home-locus reference. Opossum
ambiguity reflects incompatible FASTA/GFF versions. The amphioxus protein and
synteny searches used different *Branchiostoma* species, so the figure does not
imply gene absence.

Source figure: `analyses/overview/figures/species_lineage_evidence_summary.png`

## SynVoy overview

**Figure X. SynVoy locus-level evidence and completed review decisions for
seven genes across 14 non-human target genomes.** (A) Best retained gene-of-interest
confidence per gene and target genome after ownership/paralog filtering: H,
HIGH; M, MEDIUM; dash, no retained HIGH or MEDIUM candidate. A missing retained
candidate is not equivalent to confirmed gene loss. (B) Completed integrated
locus decision after exact-GFF accession matching and detailed exception
review: A, accepted; T, tentative; question mark, ambiguous; R, rejected. The
98 gene-by-target rows comprise 62 accepted, 20 tentative, 12 ambiguous, and
four rejected decisions. All opossum rows are ambiguous because the target
FASTA and GFF3 sequence identifiers were incompatible; independently supported
opossum proteins are not invalidated by this locus-level exclusion.

Source figure: `analyses/synteny/figures/synvoy_species_confidence_heatmap.png`

## Gene-specific SynVoy plot template

**Figure X. Conserved genomic neighbourhood around [GENE] in human and [TARGET
SPECIES].** Genes are plotted in genomic order and arrow direction represents
transcriptional orientation; homologous or retained gene-of-interest annotations
are linked/coloured according to the SynVoy output. The [GENE] candidate is
interpreted together with neighbouring-gene conservation, coordinate ownership,
protein reciprocal-hit evidence, and phylogenetic placement. Local rearrangement,
assembly fragmentation, or annotation failure can reduce visible neighbourhood
conservation and is not by itself evidence of gene loss.

Use this only after filling in the actual gene, species, plot-specific symbols,
and manual review result. Do not use one generic caption for all SynVoy panels
without checking the plotted legend.

## Gene-structure overview

**Figure X. Conservation of representative protein-coding gene structure across
five vertebrate reference annotations.** One representative transcript was
selected for each gene in human, mouse, cow, chicken, and zebrafish. (A) Number
of CDS-containing exons. Exon number was constant within each gene: BGN and DCN, seven;
EPYC and OGN, six; FMOD, PRELP, LUM, and OMD, two. (B) Estimated protein length
derived from the combined CDS. (C) Genomic gene span on a logarithmic scale.
Protein length and coding-exon count were stable, whereas locus span varied more
strongly because of intron-length differences. Results depend on transcript and
annotation choice and therefore support gene-model plausibility rather than
proving orthology alone. The chicken BGN transcript passes structural QC but
encodes an ASPN-like product, and the zebrafish BGN row represents *bgna* rather
than the canonical-panel *bgnb*.

Source figure: `analyses/gene_structure/figures/gene_structure_conservation_overview.png`

## Per-gene visual structure figures

**Figure X. Representative CDS organization of [GENE] in human, mouse, cow,
chicken, and zebrafish.** Filled boxes represent CDS segments and connecting
lines represent intervening genomic sequence, ordered in transcript orientation.
Coordinates are relative within each transcript model and the x-axis is not a
common genomic scale across species. The figure compares coding-exon number and
relative organization; differences in total span largely reflect intron length.
The selected transcript accession and predicted-model status should be reported
in the accompanying table.

Source figures: `analyses/gene_structure/figures/<GENE>_cds_gene_structure_5species.png`

## Full representative transcript structures

**Supplementary Figure X. Full representative transcript organization across
five vertebrates.** Exons are shown in transcript orientation for BGN, DCN,
FMOD, PRELP, EPYC, LUM, OGN, and context gene OMD in human, mouse, cow, chicken, and
zebrafish. Grey denotes complete exon extent, blue and orange denote derived
5-prime and 3-prime UTR segments, and gene-specific colour denotes CDS. The
x-axis is normalized within each transcript so that exon order can be compared;
the labelled genomic span is the appropriate measure of absolute locus size.
Coding structure is much more stable than intron spacing or annotated UTR
length. UTR differences are sensitive to transcript-end annotation and are not
interpreted as demonstrated regulatory divergence.

Source figure: `analyses/gene_structure/extended/figures/full_gene_structure_5species.png`

## Coding-exon boundaries and splice phase

**Supplementary Figure X. Conservation of coding-exon boundaries and intron
phase across five vertebrates.** Horizontal bars show the translated length of
the selected representative transcript, alternating shades delimit coding
exons, and symbols mark intervening intron phases 0, 1, or 2. All 26 homologous
junctions retain the same intron phase in human, mouse, cow, chicken, and
zebrafish. Exact junction position can shift modestly while preserving the
reading-frame relationship. The figure supports conserved coding organization
but does not establish tissue-specific isoform use or orthology by itself.

Source figure: `analyses/gene_structure/extended/figures/coding_exon_splice_phase_map.png`

## Domain-to-exon architecture

**Supplementary Figure X. Human LRR-family Pfam hits projected onto coding
exons.** Alternating gene-specific shades delimit the coding exons of each human
representative protein, and purple bars show significant LRR-family Pfam hits
($i$-value $\leq10^{-3}$). The map links the conserved SLRP/LRR protein core to
the underlying coding structure and shows whether domain models cross splice
boundaries. Overlap does not mean that each exon is an independent functional
module, and Pfam family support is not a gene-specific orthology test.

Source figure: `analyses/gene_structure/extended/figures/human_domain_exon_architecture.png`

## Pairwise coding-sequence constraint

**Supplementary Figure X. Human-reference pairwise coding-sequence constraint
for eight SLRP genes (seven comparative genes plus OMD context).** Protein-guided codon alignments were used to estimate
whole-sequence $d_N/d_S$ with the Nei--Gojobori method for human versus mouse,
cow, chicken, and zebrafish. All 25 finite estimates are below one, consistent
with predominant purifying constraint. The ASPN-like chicken BGN-locus product
was excluded; undefined/non-positive $d_S$ values are not plotted, and EPYC and
LUM chicken comparisons are retained with saturation caution because $d_S>2$.
The analysis is descriptive and does not test individual branches or sites.

Source figure: `analyses/evolutionary_rates/figures/pairwise_dn_ds_human_reference.png`

## Protein/domain/MSA/SignalP overview

**Figure X. Protein conservation and quality-control summary for the 103-protein
canonical SLRP panel.** (A) Mean within-gene pairwise identity in the MAFFT
alignments. (B) Mean alignment-column occupancy and modal residue conservation
at high-occupancy columns. (C) Numbers of proteins classified as Supported,
Watch, or Needs inspection by the integrated rules;
these categories are quality-control flags, not statistical confidence
intervals. (D) Proportion of proteins with significant SLRP/LRR-family Pfam
support and a SignalP-positive N terminus. All 103 proteins have retained
SignalP calls (100 positive and three negative). Conserved
secretion and LRR architecture support secreted SLRP identity, but do not alone
prove gene-specific orthology or conserved biological function.

Source figure: `analyses/protein_analysis/figures/protein_conservation_overview.png`

## Mature-protein motif-model matrix

**Supplementary Figure X. Compatibility of the 103 canonical SLRP proteins with
gene-specific mature-protein motif models.** SignalP-positive signal peptides
were removed from all 100 positive proteins using the final SignalP calls; the
three SignalP-negative proteins remained untrimmed. MEME models were learned
separately for BGN, DCN, FMOD, PRELP, EPYC, LUM, and OGN and then scanned
against all proteins with MAST.
Cells show the median per-protein sum of -log10 motif-hit p-values for the
indicated assigned gene and tested motif model; the value is a ranking and
visualization score rather than a calibrated new statistical test. The diagonal
maximum for every gene, and the best own-gene model for all 103 individual
proteins, support gene-associated sequence conservation. Because the models
were learned from the same small curated sets and no independent hold-out was
available, this is not an independent accuracy estimate. Motifs should not be
interpreted as validated biochemical binding sites.

Source figure: `analyses/protein_analysis/motifs/figures/protein_gene_motif_model_matrix.png`

## Promoter sequence QC

**Supplementary Figure X. GC composition of proximal promoters from the
five-species comparative panel.** Strand-aware -2,000/+200-bp windows were
extracted relative to one representative transcription start site for each of
seven genes in human, mouse, cow, chicken, and zebrafish. All 35 proximal windows had their
expected length, none was clipped, and none contained ambiguous bases. GC
content ranged from 32.48% to 59.06%, illustrating composition differences that
can influence de-novo motif discovery and sequence scanning.

Source figure: `analyses/regulatory/figures/promoter_proximal_gc_heatmap.png`

## Targeted promoter-motif enrichment

**Supplementary Figure X. AME enrichment results for a targeted
cartilage/osteogenesis-related transcription-factor motif panel.** Bars show
−log10 of the AME E-value for the 40-profile panel. HIF1A was the strongest
targeted result (E = 2.55), but no targeted motif reached E < 0.05. The figure
therefore does not demonstrate shared transcription-factor regulation of the
seven genes.

Source figure: `analyses/regulatory/figures/promoter_targeted_tf_ame_enrichment.png`

## Exploratory promoter recurrence

**Supplementary Figure X. Exploratory cross-species recurrence of targeted
motif occurrences in proximal promoters.** Cells report the number of five
species in which at least one FIMO occurrence was found at the uncorrected
threshold p <= 1e-4 for each gene and targeted transcription-factor family. The
strict q <= 0.05 scan returned no sites. These recurrences are hypothesis-
generating sequence matches and must not be described as transcription-factor
binding sites or conserved regulatory mechanisms.

Source figure: `analyses/regulatory/figures/promoter_targeted_tf_species_recurrence.png`

## Curated mouse phenotype evidence

**Supplementary Figure X. Breadth of curated single-gene mouse phenotype
annotations for the seven SLRPs.** Counts are unique Mammalian Phenotype terms
assigned to focal single-gene genotypes in the MGI files downloaded on 26 August
2026 and grouped through ontology ancestors into overlapping broad categories.
BGN, DCN, FMOD, and EPYC include direct skeletal/cartilage/joint annotations, whereas
LUM and OGN are represented more strongly by other ECM-associated categories.
Counts measure database annotation depth, not phenotype severity, causal effect
size, or cross-gene importance.

Source figure: `analyses/phenotypes/figures/mgi_single_gene_phenotype_heatmap.png`
