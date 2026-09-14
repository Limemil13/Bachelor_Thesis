# Supplementary coding-sequence constraint analysis

Last audited: 2026-08-29

This module asks a different question from the protein MSA: whether the coding
sequences of the five reference vertebrates show an excess of amino-acid-
changing substitutions. It is a deliberately bounded supplementary analysis,
not a branch-site selection test.

## Method

For each of BGN, FMOD, PRELP, EPYC, LUM, OGN, and context gene OMD, the five
representative proteins used in the extended gene-structure analysis were
aligned with MAFFT. Their reconstructed CDSs were then back-translated through
that protein alignment. Gap-containing and invalid codon pairs were excluded,
and human was compared separately with mouse, cow, chicken, and zebrafish.
Pairwise dN and dS were estimated with the Nei--Gojobori 1986 (NG86) method.

The known chicken BGN-locus/ASPN-like product was excluded from estimation; it
is not a valid BGN ortholog comparison. Estimates with undefined/non-positive
dS are reported as unavailable, and dS >= 2 is flagged because synonymous-site
saturation can destabilize the ratio.

## Results

- Pairwise table: `tables/pairwise_dn_ds_human_reference.tsv`
- Gene summary: `tables/codon_constraint_gene_summary.tsv`
- Protein and codon alignments: `alignments/`
- Figure: `figures/pairwise_dn_ds_human_reference.png`
- Script: `scripts/analyze_codon_constraint.py`

There are 28 possible human--target rows. One is excluded for the chicken BGN
orthology conflict. Of the remaining rows, 21 produce a finite ratio and every
one is below 1. Human--zebrafish dS is non-interpretable for six of seven genes;
PRELP is the exception. Chicken EPYC and LUM have finite ratios but dS above 2
and are flagged for saturation caution. The cleanest interpretation is that the
available pairwise estimates are compatible with purifying constraint on the
representative coding sequences.

## Limits and thesis wording

Do not write that this analysis proves purifying selection in every lineage or
at every residue. Pairwise NG86 averages across the whole sequence, does not
model branches or individual sites, and performs poorly at deep divergence once
synonymous substitutions saturate. A defensible sentence is:

> All interpretable human--target NG86 estimates were below one, consistent
> with predominant purifying constraint on the representative SLRP coding
> sequences; deep vertebrate comparisons were often unavailable because dS was
> saturated or undefined.

This layer is useful because it complements protein identity, domain structure,
phylogeny, and splice conservation with a nucleotide-level measure, while its
limitations remain visible in the table.
