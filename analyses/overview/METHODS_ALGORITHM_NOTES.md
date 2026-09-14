# Methods and algorithm notes

These notes describe what each main computational step contributes. They are a
practical companion to `../../thesis/chapters/methods.tex`: the LaTeX file is the
thesis draft, while this file records exactly how each method was applied and
where its reproducible implementation and output live.

## Exact implementation map

| Method | How it was used in this project | Main settings | Reproducible implementation | Authoritative output |
|---|---|---|---|---|
| Juvenile RNA-seq screen | Existing Salmon `quant.sf` files for three P0 mouse Prx1-lineage knee-joint samples were re-extracted. Transcript TPMs were summed by gene; mean, sample SD, coefficient of variation, and rank were calculated across WT1--WT3. The heterogeneous PRJNA936435 samples were retained only as a metadata/descriptive audit. | No new Salmon quantification was run because the preserved quantifications had 81.9--92.0% mapping and reproduced the curated values. Do not call this differential expression or isolated growth-plate tissue. | `../expression/build_corrected_mouse_expression.py` and `../expression/sample_metadata_corrected.tsv` | `../expression/tables/mouse_expression_descriptive_metadata_corrected.tsv` |
| Direct growth-plate expression | Deposited GSE114919 mouse and rat matrices were subset to the focal SLRPs and summarized within species, age, bone, and proliferative/hypertrophic zone. Tibial one- and four-week groups contain five replicates each. | Retain the authors' processed scale; do not pool it with Salmon TPMs or compare absolute mouse and rat values. | `../expression/gse114919/summarize_gse114919_growth_plate.py` | `../expression/gse114919/tables/gse114919_slrp_condition_summary.tsv` and `gse114919_slrp_tibia_cross_condition_summary.tsv` |
| Bgee expression context | Bgee 16 calls for mouse, chicken, and zebrafish were retrieved and filtered to relevant anatomy/stage terms. | Qualitative presence/absence-of-evidence only; unmatched calls are missing evidence, not biological absence. | `../expression/fetch_bgee_slrp_expression.py` and `summarize_bgee_slrp_expression.py` | `../expression/tables/bgee_slrp_relevant_expression_summary.tsv` |
| Candidate discovery and reciprocal BLASTP | Human proteins were queries against local NCBI/RefSeq proteomes; candidate proteins were searched back against the human SLRP reference panel. Annotation, query coverage, length, duplicate accessions, and reciprocal assignment were combined. | BLASTP 2.12.0+, E-value `1e-5`; initial median length-ratio screen 0.7--1.3. Exceptions were reviewed rather than rejected solely by length. | `../protein_analysis/scripts/build_reciprocal_blastp_panel.py` and gene-specific curation scripts in the same directory | `../protein_analysis/candidates/canonical_candidate_manifest.tsv` and `tables/final_candidate_level_confidence.tsv` |
| Pfam domain scan | All 103 canonical proteins were scanned against Pfam-A with HMMER `hmmscan`; retained hits were collapsed into protein- and gene-level SLRP/LRR summaries. | HMMER 3.4; independent domain E-value `<=1e-3`. Domain support establishes SLRP-family architecture, not a unique paralog assignment. | `../protein_analysis/scripts/run_canonical_refresh.sh` and `../protein_analysis/scripts/tables.py` | `../protein_analysis/tables/protein_conservation_domain_msa_signalp_summary.tsv` |
| Signal peptide prediction | Canonical proteins were submitted to SignalP 6.0 and merged by stable sequence identifier. | Eukaryote, slow-sequential mode. All 103 canonical proteins have retained calls: 100 positive and three negative. The final 17-sequence job was `6AA7E4B90014723A58E1A86B`. | `../protein_analysis/signal_peptides/scripts/prepare_signalp_input.py` and `merge_signalp_summaries.py` | `../protein_analysis/signal_peptides/tables/signalp_summary_clean.tsv` and the integrated protein table |
| Multiple-sequence alignment and conservation | One canonical protein alignment per gene was generated. Pairwise identity used only columns containing residues in both sequences; column occupancy, modal-residue conservation, majority consensus, identity matrices, and occupancy-weighted sequence logos were then calculated. | MAFFT 7.505 `--auto`; untrimmed alignments retained for manual review. trimAl 1.5 `-automated1` was used for the combined tree input. | `../protein_analysis/scripts/run_canonical_refresh.sh` and `build_msa_conservation_outputs.py` | `../protein_analysis/alignments/canonical/`, `../protein_analysis/tables/msa_conservation_gene_summary.tsv`, and `../protein_analysis/figures/sequence_logos/` |
| Maximum-likelihood phylogeny and iTOL | Seven compact per-gene trees and one combined 103-tip tree were inferred from the canonical alignments. Tip labels and gene/taxonomic colour files were generated for iTOL. | IQ-TREE 2.0.7, `-m MFP`, `-alrt 1000`, `-bb 1000`, automatic threads. Trees are interpreted as unrooted unless an explicit outgroup is stated. | `../phylogenetics/compact_panel/scripts/` and `../phylogenetics/combined_trees/six_gene_tree/run_selected_working_genes_tree.sh` | `../phylogenetics/compact_panel/trees/` and `../phylogenetics/combined_trees/six_gene_tree/SLRP_selected_working_genes.treefile` |
| ProtSpace/ProtT5 | The 103 canonical proteins plus four controls were embedded with ProtT5 and projected by PCA, UMAP, and t-SNE. A leave-one-out nearest-centroid test was calculated in the original embedding space. | Projection plots are exploratory; the quantitative centroid check is the primary embedding QC. | `../protein_analysis/protspace/scripts/` | `../protein_analysis/protspace/tables/protspace_canonical_embedding_review.tsv` and `protspace_canonical_gene_summary.tsv` |
| SynVoy microsynteny | Each human focal locus was compared with 14 target genomes and their GFF3 files. Completed reports were parsed after candidate ownership/paralog filtering into one gene-by-species evidence table and P1/P2/P3 review strata. | Human home locus; protein mode; 10 flanking genes. Updated BGN/DCN runs used MMseqs sensitivity 9.5, non-strict gene-family tokens, automatic generic presets disabled, 8 GB MMseqs split-memory limit, and one iterative-search CPU. Older completed reports retain their own run provenance. | `../synteny/scripts/build_synvoy_review.py`, `run_updated_synvoy_panel.sh`, and `SYNVOY_MANUAL_REVIEW_GUIDE.md` | `../synteny/tables/synvoy_gene_species_evidence.tsv` and the three completed-review tables |
| Manual protein and locus review | AliView/Jalview review checked termini, cysteines, LRR-core continuity, unique indels, and agreement with tree/synteny/SignalP. For all 98 locus rows, the canonical protein accession was located in the exact GFF and compared with the post-ownership SynVoy coordinate. Twenty-five P1 exceptions used explicit feature overlap and oriented flanking-gene comparison; P2/P3 used reproducible accession--locus confirmation. | Default detailed window: 15 genes per side (40 for catshark BGN). Final terms are `accepted`, `tentative`, `rejected`, or `ambiguous`. FASTA/GFF sequence IDs were audited; all opossum locus calls are excluded pending matched-input rerun. | `../synteny/scripts/audit_synvoy_canonical_loci.py`, `finalize_all_synvoy_reviews.py`, `compare_gff_neighborhoods.py`, and `../synteny/SYNTENY_MANUAL_REVIEW_GUIDE.md` | `../synteny/diagnostics/priority_locus_reviews/all_synvoy_review_decisions.tsv`, `canonical_locus_audit.tsv`, and `tables/final_candidate_level_confidence.tsv` |
| Elephant-shark compound-model diagnostic | The nine annotated CDS segments of `XP_007897806.2` were reconstructed from the exact SynVoy FASTA/GFF3 pair, translated, compared with eight human SLRP references by local BLASTP, and scanned with Pfam. The two SLRP-like halves were evaluated separately from locus-level synteny. | Complete 2055-bp CDS; 684 aa; terminal stop and no internal stop. BLASTP: N half DCN-best, C half FMOD-best. HMMER 3.4/Pfam-A detects LRRs and two LRRNT caps. Hydrophobic stretches were inspected directly and are not formal SignalP calls. | `../synteny/scripts/extract_gff_protein.py` and `compare_protein_to_reference_queries.py` | `../synteny/diagnostics/priority_locus_reviews/FMOD_elephant_shark/XP_007897806.2_diagnostic_summary.tsv` and raw outputs in the same directory |
| Representative gene structure | Gene--transcript--exon--CDS relationships were parsed for human, mouse, cow, chicken, and zebrafish. One representative transcript per gene/species was selected; CDSs were reconstructed and checked for reading-frame completeness. UTRs, intron phases, alternative coding transcripts, and Pfam-to-exon overlap were summarized. | Eight genes including context gene OMD; 40 representative transcripts. Checks: length divisible by three, start methionine, terminal stop, no internal stop, and GFF phase continuity. | `../gene_structure/scripts/build_extended_gene_structure.py`, `map_domains_to_exons.py`, and `run_extended_gene_structure.sh` | `../gene_structure/extended/tables/full_gene_structure_qc.tsv` and `../gene_structure/extended/figures/full_gene_structure_5species.svg` |
| Pairwise coding constraint | Representative proteins were aligned and corresponding CDSs were back-translated. Human was compared separately with mouse, cow, chicken, and zebrafish. | Pairwise Nei--Gojobori (Biopython); gap/invalid codon pairs omitted; non-positive/undefined dS reported unavailable; dS `>=2` flagged as saturation-prone. Not a branch/site selection test. | `../evolutionary_rates/scripts/analyze_codon_constraint.py` | `../evolutionary_rates/tables/pairwise_dn_ds_human_reference.tsv` |
| Mature-protein motifs | Signal peptides were removed at the predicted cleavage site from all 100 final SignalP-positive proteins; the three SignalP-negative proteins remained untrimmed. Per-gene MEME models were scanned across the complete panel with MAST, and STREME contrasted each gene with the other genes. | MEME 5.5.8 protein ZOOPS; widths 6--40 aa; E-value stop 0.05. No independent holdout, so motif scores are supporting classification rather than functional-site evidence. | `../protein_analysis/motifs/scripts/run_protein_motif_analysis.sh` | `../protein_analysis/motifs/tables/` |
| Promoter motifs | Strand-aware core and proximal promoters were extracted from the five representative assemblies, checked, and analysed with STREME, Tomtom, AME, and FIMO. | Core -500/+100 bp; proximal -2000/+200 bp. STREME widths 6--15 bp, dinucleotide-shuffled controls, seed 20260826. Tomtom: JASPAR 2026 CORE vertebrates. AME: Fisher/total-hit. Strict FIMO q `<=0.05`; p `<=1e-4` retained only as exploratory. | `../regulatory/scripts/` | `../regulatory/tables/` and `../regulatory/figures/` |
| Phenotype associations | Pinned MGI genotype/MP and HPO gene--phenotype/disease files were filtered for the seven genes; MP parent categories were propagated and single-gene versus compound genotypes separated. | MGI downloads 2026-08-26; MP ontology 2026-07-22; HPO v2026-06-23. Counts measure curation depth, not effect size. | `../phenotypes/scripts/` | `../phenotypes/tables/mgi_phenotype_evidence_summary.tsv` and `hpo_human_phenotype_evidence_summary.tsv` |
| Evidence integration | Expression, literature, protein, tree, synteny, gene-structure, manual-review, motif, and phenotype evidence were joined without converting them into one unsupported numerical truth score. | Each evidence layer retains its own limitations; conflicts remain visible in notes and confidence decisions. | `scripts/build_final_evidence_tables.py`, `scripts/build_results_synthesis.py`, and `verify_thesis_analysis_outputs.py` | `tables/final_gene_level_evidence.tsv` and `tables/final_candidate_level_confidence.tsv` |

The current LaTeX Methods draft already converts this implementation into
continuous thesis prose. Before final submission, update the SignalP counts and
any opossum rerun result once from the authoritative tables; do not manually
synchronize several copies during active analysis.

## Candidate ortholog selection: BLASTP and reciprocal checking

BLASTP finds local amino-acid similarities using word matches that seed gapped
extensions. Scores are evaluated against an empirical substitution matrix and
converted into E-values. For the compact panel, the best candidate in each local
species proteome was checked by searching it back against the human reference
proteome. Agreement with the original human query is useful orthology evidence,
but it is not a proof: recent paralogs, incomplete proteomes, and deep-lineage
duplications can still mislead reciprocal best-hit logic. That is why domain,
tree, annotation, and later synteny evidence are combined with it.

## Protein domains: HMMER hmmscan against Pfam

Pfam families are represented by profile hidden Markov models. A profile HMM
captures position-specific residue preferences and insertion/deletion behavior,
so it is more sensitive to remote homology than a single pairwise sequence
comparison. `hmmscan` scores each protein against the library, and the parsed
table retains statistically significant domain matches. SLRP-like LRR/domain
architecture supports family membership; it does not by itself distinguish
closely related SLRP paralogs.

## Multiple sequence alignment: MAFFT

MAFFT aligns homologous proteins using fast Fourier transform-derived similarity
information followed by progressive and iterative refinement, depending on the
selected strategy. The alignment establishes homologous columns for conservation
statistics and phylogenetic inference. Long terminal insertions, fragments, and
incorrect gene models can create gap-heavy or poorly aligned regions, which are
why sequence-level MSA QC and manual review are retained.

## Alignment trimming: trimAl `automated1`

trimAl removes columns judged unreliable using gap and similarity properties.
The `automated1` heuristic selects a trimming strategy from features of the
input alignment. Trimming can reduce noise in a protein phylogeny, but it also
removes data; both untrimmed MSA QC and the exact trimmed tree input are kept.

## Conservation summaries

Pairwise identity is the fraction of aligned comparable residues that match.
Column conservation records residue agreement at each MSA position, while the
majority consensus reports the most frequent residue at each retained column.
The amino-acid logos display residue frequencies scaled by Shannon information
content; the stack height is also multiplied by non-gap occupancy so sparsely
represented alignment columns do not look maximally conserved. Logo coordinates
are MAFFT alignment columns, not human-protein residue numbers.
These measures summarize conservation; they are not independent observations
when species share evolutionary history, and a conserved protein is not
automatically specific to the growth plate.

## Maximum-likelihood phylogeny: IQ-TREE

IQ-TREE searches for the topology and branch lengths that maximize the
likelihood of the observed alignment under an amino-acid substitution model.
ModelFinder compares candidate models with an information criterion. Ultrafast
bootstrap approximates resampling support efficiently, while SH-aLRT evaluates
local branch support by comparing nearby topologies. Support describes stability
under the analysis assumptions, not proof of orthology. The combined SLRP tree
is unrooted until a biologically justified root/outgroup is applied; the left-to-
right display order in a Newick viewer must not be interpreted as ancestry.

For large per-gene trees, the completed Newick is converted to iTOL-safe tip
identifiers and paired with display labels plus a broad taxonomic color strip by
`analyses/phylogenetics/per_gene_trees/scripts/finalize_large_gene_itol.sh`.

## Signal peptides: SignalP 6.0

SignalP 6.0 uses a protein language-model representation and a sequence-labeling
model to classify N-terminal targeting peptides and estimate cleavage sites.
The eukaryotic slow-sequential run was used here. A positive classical signal
peptide is biologically expected for secreted extracellular-matrix SLRPs. A
negative or uncertain result is a QC flag, not automatic evidence that the gene
is wrong, because incomplete N termini and annotation errors can cause failure.

## Protein-space analysis: ProtT5 embeddings and ProtSpace

ProtT5 converts each amino-acid sequence into a high-dimensional numerical
representation learned from large protein-sequence corpora. ProtSpace compares
the per-protein embeddings and projects them with PCA, UMAP, and t-SNE. PCA is a
linear variance-preserving projection; UMAP builds a weighted neighbor graph;
t-SNE emphasizes local neighborhoods. Apparent 2-D clusters depend on projection
parameters and should be treated as supporting pattern evidence. The quantitative
leave-one-out centroid check is preferable to judging only the plot.

## Gene structure from NCBI GFF3

The streaming extractor links gene, transcript, CDS, and exon features by their
GFF3 parent relationships, then selects one representative protein-coding
transcript per gene/species. It records CDS exon count, genomic span, strand,
and protein length. Comparable exon/CDS structure can support annotation and
evolutionary consistency, but transcript choice and assembly quality matter;
predicted `XM_` models are therefore prioritized for manual review.

The extended extractor subtracts CDS intervals from exons to derive annotated
5' and 3' UTR segments, reconstructs CDS sequence from the indexed genome, and
checks translation completeness. GFF `phase` records how many bases of a codon
are carried over at the start of a CDS block; comparing this with cumulative CDS
length tests reading-frame continuity. Intron phase describes whether a splice
junction falls between codons (0) or after the first (1) or second (2) base of a
codon. Conserved intron phase is stronger structural evidence than equal exon
count, although it still depends on transcript annotation.

Pfam HMM hits were projected from amino-acid to cumulative-CDS coordinates and
intersected with coding exons. Crossing a splice boundary means the statistical
domain/repeat model is encoded by more than one exon; it does not demonstrate
that an exon is an independent protein module.

## Pairwise coding constraint (NG86)

Protein-guided codon alignments preserve codon triplets while using amino-acid
alignment to place gaps. The Nei--Gojobori method estimates nonsynonymous and
synonymous substitutions per corresponding site and reports their ratio.
Values below one are compatible with purifying constraint, but pairwise whole-
sequence ratios average across sites and branches. When synonymous changes are
saturated, dS can be unstable, non-positive, or undefined and the ratio should
not be interpreted. The thesis analysis therefore leaves those deep comparisons
unavailable rather than changing methods selectively.

## SynVoy synteny

SynVoy compares conserved genomic neighborhoods around the human query gene in
the selected species. Conserved neighbors provide evidence that a candidate lies
in an orthologous genomic context, especially when protein paralogs are close.
Failure to display an SLRP is not automatically absence: assembly gaps, gene-name
differences, annotation quality, contig boundaries, and search-window choices can
all disrupt the neighborhood. The JSON review tables prioritized questionable
calls instead of requiring equal manual effort for every species. The completed
review then mapped each canonical protein accession back to the exact target
GFF, tested overlap with the post-ownership candidate, and escalated coordinate
conflicts to explicit feature and flanking-gene comparison. Whole blocks were
oriented relative to the focal gene before order was assessed, so a reversed
chromosomal segment could still count as conserved. Exact gene-symbol matching
was kept conservative; `LOC` labels and lineage-specific suffixes were not
silently treated as exact matches. A separate sequence-ID audit exposed the
opossum FASTA/GFF mismatch and prevented technically invalid scores from being
interpreted biologically.

## RNA-expression evidence

The local P0 measurements are descriptive TPM values from Prx1-lineage cells
isolated around the juvenile knee; they establish expression in a relevant
chondroprogenitor lineage but are not anatomically isolated growth-plate zones.
The mouse/rat GSE114919 laser-capture data provide the more direct proliferative
and hypertrophic growth-plate evidence. Bgee contributes curated expression
calls. Values from these different experiments are ranked or interpreted within
their own datasets and are not pooled into a single quantitative expression
scale.

Juvenile RNA data make biological sense for a growth-plate thesis because the
growth plate is an active developmental structure. The main limitation is tissue
and cell-state specificity, not the young age itself. Gene selection should
therefore combine expression with direct skeletal biology and evolutionary/QC
evidence rather than use expression alone.
