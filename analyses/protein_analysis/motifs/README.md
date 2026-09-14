# Protein motif analysis

Last rebuilt: 2026-09-14

## Scope

This supplementary module asks whether the 103 curated proteins retain
gene-associated sequence patterns beyond broad SLRP/LRR Pfam annotation and
whether any candidate behaves as a motif-model outlier. It is an unaligned
sequence-pattern check, not a biochemical motif-function assay.

The current panel is BGN 13, DCN 15, FMOD 14, PRELP 16, EPYC 15, LUM 16, and
OGN 14. SignalP-positive N termini were removed at the recorded cleavage site.
The three SignalP-negative proteins (whale-shark BGN, opossum DCN, and
spotted-gar OGN) remained untrimmed. Exact processing is recorded in
`inputs/protein_motif_input_metadata.tsv`; thus all final SignalP calls were
incorporated before the current MEME/MAST/STREME run.

## Analyses

1. MEME 5.5.8 was run in protein/ZOOPS mode with widths 6--40 aa and an
   E-value stopping threshold of 0.05 on the combined set and on each gene.
2. Each gene-specific model was scanned against all 103 proteins with MAST.
   `sum(-log10(p))` over significant hits is used only as a ranking/plotting
   score, not a newly calibrated statistical test.
3. STREME compared each gene with the other six genes and returned five
   patterns per comparison. The small groups did not permit a separate hold-out
   set, so these are training-set enrichment patterns.

## Result

- 116 MEME motifs: 37 combined-panel and 79 gene-level motifs.
- 35 STREME gene-enriched patterns.
- 721 protein/model score rows.
- All 103 proteins ranked the model learned from their assigned gene first.

This supports gene-associated conservation and supplies no additional motif-
model outlier. It does not erase independent partial-model, reciprocal-hit,
tree, synteny, or SignalP concerns. Because the same curated proteins were used
to learn and evaluate the models and no independent hold-out set was available,
103/103 is not an independent classification-accuracy estimate. The motifs are
mostly conserved LRR/cysteine-rich sequence blocks and must not be labelled as
validated collagen-binding, receptor-binding, or signalling sites.

## Authoritative outputs

- `tables/protein_meme_motif_summary.tsv`
- `tables/protein_streme_gene_enriched_motifs.tsv`
- `tables/protein_candidate_motif_model_scores.tsv`
- `tables/protein_candidate_motif_model_assignments.tsv`
- `tables/protein_global_motif_gene_coverage.tsv`
- `figures/protein_gene_motif_model_matrix.png`
- `figures/protein_candidate_motif_model_margins.png`
- Tool-native dated reports: `raw_outputs_7gene_20260914/`
- Previous summary snapshot: `archive/summary_20260903/`

The Ghostscript warnings from MEME affected only its optional EPS-to-PNG
conversion. XML/text reports and the project-generated summary figures are
complete.

## Reproduction and thesis placement

Use `scripts/prepare_protein_motif_inputs.py`,
`scripts/run_protein_motif_analysis.sh`, and
`scripts/summarize_protein_motifs.py`. Keep the result in a short supplementary
paragraph/figure; Pfam, MSA, and phylogeny remain the primary interpretable
protein/orthology evidence.
