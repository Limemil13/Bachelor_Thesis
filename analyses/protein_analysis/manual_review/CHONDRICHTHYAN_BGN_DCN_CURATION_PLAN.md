# Chondrichthyan BGN/DCN locus-curation plan

Last updated: 2026-09-03

## Purpose and hypothesis

This is a bounded discovery extension, not a replacement for the seven-gene
comparative analysis. It tests the hypothesis that BGN and DCN were already
distinct in early jawed vertebrates but incomplete or incorrect genome
annotations obscure their assignment in cartilaginous fishes.

The focal species are:

- elephant shark/chimaera (*Callorhinchus milii*);
- catshark (*Scyliorhinus canicula*);
- whale shark (*Rhincodon typus*).

Human and spotted gar provide annotated BGN/DCN locus and protein anchors.
ASPN is the close class-I paralog control. Lamprey can provide deep-lineage
context but should not be treated as the unquestioned tree root.

## Why these loci are worth resolving

- Historical cartilaginous-fish records assigned to BGN were found to be
  DCN/DCN-like during the current panel audit.
- The replacement catshark BGN `XP_038638277.1` is only 268 aa and covers about
  68% of the human query.
- The replacement whale-shark BGN `XP_048475905.1` is only 161 aa and covers
  about 39% of the human query.
- Elephant-shark `NP_001279248.1`, catshark `XP_038636759.1`, and whale-shark
  `XP_048463588.1` support DCN rather than BGN.
- Older SLRP-family literature did not recover BGN confidently in the available
  shark annotations. That is an annotation/search observation, not proof of
  true absence.

## Existing protein inputs

- BGN panel:
  `analyses/protein_analysis/candidates/canonical_by_gene/BGN_canonical.faa`
- DCN panel:
  `analyses/protein_analysis/candidates/canonical_by_gene/DCN_canonical.faa`
- Human ASPN/class-I context:
  `analyses/phylogenetics/combined_trees/slrp_family_reference_tree/SLRP_reference_context.faa`
- Current BGN/DCN diagnostic:
  `analyses/phylogenetics/compact_panel/tables/class_I_BGN_DCN_tree_discrimination.tsv`
- Candidate decisions:
  `analyses/overview/tables/final_candidate_level_confidence.tsv`

The updated BGN and DCN SynVoy results should be added as locus inputs after
their reports exist. Do not infer absence from an empty annotation track.

## UGENE/manual reconstruction procedure

For each candidate locus:

1. Import the genomic interval and its GFF/GTF annotation into UGENE.
2. Search the interval separately with human and spotted-gar BGN, human and
   spotted-gar DCN, and human ASPN proteins.
3. Search again with individual reference coding exons when a full-protein
   search misses a short or highly divergent exon.
4. Require candidate exon matches to be collinear and on one genomic strand.
5. Inspect exon boundaries, splice-site dinucleotides, reading frame, start and
   stop codons, assembly gaps, and premature stops.
6. Assemble a provisional CDS only from explicitly recorded genomic
   coordinates. Translate it and retain both nucleotide and amino-acid forms.
7. Compare the translation with BGN, DCN, and ASPN alignments. Pay particular
   attention to continuous LRR-core alignment, conserved cysteines, terminal
   completeness, and large candidate-specific insertions/deletions.
8. Re-run Pfam/HMMER and SignalP for any revised translation.
9. Add the translation to the focused BGN-DCN-ASPN phylogeny. Do not classify it
   solely from the nearest BLAST hit or from a two-dimensional embedding.
10. Confirm the locus using conserved neighboring genes and distinguish a gene
    model problem from an assembly/scaffold boundary.

## Evidence table to complete

Record one row per focal locus with these fields:

```text
species
assembly
scaffold_or_chromosome
candidate_gene
candidate_accession_or_new_model_id
strand
genomic_start
genomic_end
coding_exon_count
start_codon_status
stop_codon_status
splice_site_status
frameshift_or_internal_stop
assembly_gap_status
protein_length_aa
Pfam_status
SignalP_status
BGN_DCN_ASPN_tree_placement
flanking_gene_support
final_decision
evidence_note
```

Allowed final decisions are:

- `complete orthologous gene model`;
- `partial orthologous gene model`;
- `probable pseudogene/remnant`;
- `assembly-limited`;
- `wrong paralog/annotation`;
- `unresolved`.

## Thesis-facing outputs

If at least one locus can be resolved, create:

1. a BGN/DCN chondrichthyan locus-evidence table;
2. one BGN and one DCN microsynteny map;
3. a focused BGN-DCN-ASPN tree;
4. a compact exon-model diagram for reconstructed or partial candidates;
5. a supplementary UGENE screenshot showing one representative reconstruction.

The discovery claim must match the evidence. Appropriate wording is
"homology-guided reconstruction of a previously incomplete candidate model" or
"resolution of a paralog annotation conflict." Do not write "first discovery,"
"gene loss," or "functional ortholog" without a literature audit and the
necessary locus and functional evidence.

## Stop rule

Do not expand this into an unrestricted survey of more fish species. Stop after
the three focal cartilaginous-fish genomes unless a specific assembly gap can be
resolved by one clearly better assembly. The core thesis remains complete even
if every focal locus is classified as assembly-limited or unresolved.
