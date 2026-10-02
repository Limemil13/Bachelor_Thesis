# Thesis workflow overview

Last updated: 2026-09-07

## Biological question

Which SLRP genes are credible growth-plate/cartilage ECM candidates, and how
well are their orthologous proteins and loci conserved across a lineage-spanning
vertebrate panel?

The final full comparative panel is BGN, DCN, FMOD, PRELP, EPYC, LUM, and OGN.
OMD is retained as cartilage--bone-interface context. OGN is the most exploratory
main member.

## Workflow and role of each layer

1. **Expression-guided prioritization.** GSE305415-derived P0 Prx1-lineage
   TPMs, re-extracted from stored Salmon quantifications, provide a juvenile
   skeletal-lineage screen. GSE114919 mouse/rat microdissected zones
   provide direct growth-plate support, and Bgee adds qualitative context.
   Expression prioritizes genes but does not prove function or orthology.

2. **Literature/function evidence.** Primary papers establish known ECM,
   collagen, cartilage, bone, developmental, and signalling roles and state what
   remains unknown. Literature provides the biological rationale against which
   conservation results are interpreted.

3. **Candidate discovery and reciprocal checking.** Human sequences anchor
   searches. Length, annotation, reciprocal human SLRP hit, and locus evidence
   distinguish likely orthologs from close paralogs and incomplete models.

4. **Protein architecture.** Pfam/HMMER tests for the LRR-rich SLRP scaffold;
   SignalP tests for a classical N-terminal secretion signal. These are
   complementary family/QC evidence, not gene-specific proof.

5. **MSA and conservation.** MAFFT alignments support global comparison.
   Pairwise identity, column occupancy, consensus, cysteine/LRR continuity, and
   sequence logos describe conserved cores and expose truncations or unique
   insertions.

6. **Protein-space and motif checks.** ProtT5/ProtSpace provides an
   alignment-independent outlier view. MEME/STREME/MAST tests gene-associated
   sequence patterns. Both are secondary diagnostics; neither establishes a
   binding site or replaces phylogeny.

7. **Phylogeny.** Curated IQ-TREE per-gene and combined trees test whether
   candidate labels form coherent unrooted groups. The 103-tip tree is the
   central result. A focused BGN/DCN nearest-neighbour diagnostic qualifies the
   non-monophyletic BGN pattern. Large trees are optional supplementary analyses.

8. **Gene structure and coding constraint.** NCBI GFF3 gene/transcript/exon/CDS
   relationships support representative model completeness, conserved coding-
   exon organization, intron phase, and LRR-domain placement. Protein-guided
   NG86 provides a bounded whole-CDS constraint check. These layers support
   stable models but do not demonstrate tissue-specific isoform use.

9. **Synteny.** SynVoy identifies/prioritizes candidate loci and conserved
   neighborhoods. Manual review verifies gene models, flanking genes,
   order/orientation, protein evidence, and tree placement. Apparent absence can
   reflect annotation or assembly failure and is not automatically gene loss.

10. **Supplementary regulatory and phenotype context.** Five-species promoter
    scans test exploratory motif recurrence; no targeted motif survived global
    correction. Pinned MGI/HPO releases provide phenotype context but annotation
    counts are not effect sizes.

11. **Evidence integration.** The final gene and candidate tables combine the
    layers without calculating a biologically unjustified single score. Final
    decisions use concordance plus explicit caveats.

## Current quantitative snapshot

- 103 canonical proteins: BGN 13, DCN 15, EPYC 15, FMOD 14, LUM 16, OGN 14,
  PRELP 16.
- 693 significant Pfam hits; seven canonical MSAs/logos.
- 103 SignalP calls: 100 positive, three negative; none pending.
- ProtSpace: 103 SLRPs plus four controls.
- Protein motifs: all 103 rank their own gene model best, with no independent
  hold-out set.
- Compact combined tree: 103 tips; five of seven genes unrooted-monophyletic;
  LUM 15/16; BGN qualified by the class-I diagnostic.
- Gene structure: 40 transcripts, eight genes x five species; all reconstructed
  CDSs complete; 26 comparable junctions conserve intron phase.
- Promoters: 35 loci/70 windows; no corrected targeted TF-motif result.
- Phenotypes: seven-gene MGI/HPO summaries, including DCN.
- SynVoy: 98 rows for seven genes, including completed updated-tool DCN and BGN
  runs. The FMOD refresh was interrupted before its final report; PRELP and
  EPYC refreshes have not been launched.

## Evolutionary scope

The sampled species cover major vertebrate lineages rather than every species.
Living tunicates, amphioxus, and lamprey are not “the species where SLRPs
originated”; they are extant lineages useful for bracketing early chordate and
vertebrate history. Literature supports a tunicate SLRP-like precursor and
early vertebrate expansion. Deep candidates are therefore expected to be harder
to assign gene-specifically. See `SLRP_EVOLUTIONARY_ORIGIN_AND_SPECIES.md`.

An **outgroup** is selected to orient a tree based on known evolutionary
relationships. An **outlier** is an unusual observation. The thesis does not
need a deliberately bad sequence; real, documented edge cases already test the
workflow without forcing an artificial outlier.

## Source order

1. Curated FASTAs and canonical manifest.
2. Pfam, SignalP, MSA/conservation, ProtSpace, and motifs.
3. Compact trees and clade diagnostics.
4. Gene structure and coding constraint.
5. Expression and literature.
6. SynVoy evidence plus manual decisions.
7. Final gene/candidate evidence tables and thesis figures.

Exact locations are in `AUTHORITATIVE_RESULT_LOCATIONS.md`; exact tool logic and
limitations are in `METHODS_ALGORITHM_NOTES.md`.
