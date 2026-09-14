# Chicken BGN/ASPN neighborhood diagnostic

This comparison tests whether the chicken locus labelled `BGN` is more
consistent with the conserved BGN or ASPN genomic neighborhood. It uses the
same NCBI GFF3 assemblies supplied to SynVoy and considers the ten nearest
protein-coding genes on each side of each focal gene.

## Reproduce

Set `SYNVOY_ROOT` to the local SynVoy checkout, then run from the thesis
repository root:

```bash
python analyses/synteny/scripts/compare_gff_neighborhoods.py \
  --locus chicken_BGN "$SYNVOY_ROOT/pro_panel/targets_15/gff/chicken.gff" BGN@NC_052543.1 \
  --locus chicken_DCN "$SYNVOY_ROOT/pro_panel/targets_15/gff/chicken.gff" DCN@NC_052532.1 \
  --locus anole_BGN "$SYNVOY_ROOT/pro_panel/targets_15/gff/anole.gff" BGN@NC_085842.1 \
  --locus human_BGN "$SYNVOY_ROOT/pro_panel/genomes/gff/human.gff" BGN@NC_000023.11 \
  --locus human_ASPN "$SYNVOY_ROOT/pro_panel/genomes/gff/human.gff" ASPN@NC_000009.12 \
  --flank-count 10 \
  --output-dir analyses/synteny/diagnostics/chicken_bgn_asporin_neighborhood
```

The `GENE@SEQID` form fixes the comparison to the primary chromosome when the
GFF3 also contains alternate-locus records.

## Result

- Chicken `BGN` versus human `ASPN`: six exact shared flanking genes (`CENPP`,
  `ECM2`, `IPPK`, `NOL8`, `OGN`, and `OMD`); five of six retain a consistent
  order after orienting the neighborhoods by the focal-gene strand.
- Chicken `BGN` versus human `BGN`: no exact shared flanking genes in the
  10-gene windows.
- Chicken `BGN` versus anole `BGN`: no exact shared flanking genes.
- Anole `BGN` versus human `BGN`: ten exact shared flanking genes, all ten in
  consistent order. This positive control shows that the method detects the
  conserved BGN neighborhood across amniotes.
- Chicken `BGN` versus chicken `DCN`: no shared flanking genes, confirming that
  the SynVoy DCN cross-hit is a separate locus.

Together with the ASPN reciprocal best hit and focused protein tree, the
neighborhood evidence strongly supports interpreting chicken `XP_414298.2` as
an ASPN-like protein at the conserved ASPN locus, despite the GFF3 gene symbol
`BGN`. The record remains inappropriate for the canonical BGN protein panel.

Exact gene-symbol matching is intentionally conservative. Failure to share a
symbol is not by itself evidence of gene loss because annotation gaps, renamed
orthologs, and assembly breaks can hide otherwise conserved neighborhoods.
