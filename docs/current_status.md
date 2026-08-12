# Current thesis status

## Core thesis direction

Expression-guided evolutionary characterization of growth-plate/cartilage-associated SLRP genes.

The thesis starts with mouse growth plate RNA-seq expression analysis and uses this to prioritize SLRP genes for downstream evolutionary analysis.

## RNA-seq status

Completed:

- Mouse growth plate RNA-seq samples quantified with Salmon.
- Comparison tissues included: skin, heart, kidney, and brain.
- Gene-level TPM table generated.
- Expression heatmaps generated using raw TPM and log10(TPM + 1).
- Cartilage markers confirmed growth plate identity.
- SLRPs separated into growth-plate-enriched and broader ECM-associated groups.

Main expression-supported candidate genes:

- EPYC
- FMOD
- BGN
- OGN
- OMD
- PRELP

## Orthology / phylogeny / domain status

Completed or partly completed for:

- BGN
- FMOD
- OGN
- EPYC
- PRELP

Existing analysis layers:

- BLASTP candidate discovery
- reciprocal BLASTP validation
- length filtering
- Pfam/HMMER domain analysis
- MAFFT alignments
- MSA quality summaries
- preliminary IQ-TREE trees
- iTOL tree annotation files

OGN requires additional care because of alignment gaps and possible sequence outliers.

## Synteny status

SynVoy has been run for several prioritized genes, including FMOD, EPYC, BGN, and PRELP. OGN remains to be rerun/debugged.

Synteny is mainly used as orthology support. Conserved genomic neighborhood supports that candidate genes are true orthologs.

## Protein architecture / motif status

Already available:

- Pfam/HMMER LRR-domain analysis.

Planned:

- Signal peptide prediction to test whether candidates retain secreted ECM protein architecture.
- MEME motif analysis as optional support for conserved motifs.
- Targeted gene structure checks for suspicious/outlier candidates.

## Next steps

1. Add signal peptide analysis.
2. Update domain/protein architecture summary table.
3. Finish OGN SynVoy.
4. Trim and rerun final gene trees consistently.
5. Add MEME motif analysis if useful.
6. Build final evidence/confidence table with expression, length, domain status, signal peptide, MSA quality, tree placement, synteny support, and final confidence.
