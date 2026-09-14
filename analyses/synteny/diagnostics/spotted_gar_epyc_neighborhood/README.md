# Spotted-gar EPYC neighborhood review

This review tests the annotated spotted-gar `epyc` locus directly because the
EPYC SynVoy run retained no HIGH/MEDIUM spotted-gar GOI call. It uses the exact
NCBI GFF3 assemblies supplied to SynVoy and compares the 15 nearest
protein-coding genes on each side of human and spotted-gar EPYC.

## Reproduce

Set `SYNVOY_ROOT` to the local SynVoy checkout and run from the thesis
repository root:

```bash
python analyses/synteny/scripts/compare_gff_neighborhoods.py \
  --locus human_EPYC "$SYNVOY_ROOT/pro_panel/genomes/gff/human.gff" EPYC@NC_000012.12 \
  --locus spotted_gar_EPYC "$SYNVOY_ROOT/pro_panel/targets_15/gff/spotted_gar.gff" epyc@NC_090702.1 \
  --flank-count 15 \
  --output-dir analyses/synteny/diagnostics/spotted_gar_epyc_neighborhood
```

## Result

- Human EPYC is at `NC_000012.12:90963682-91004972`; spotted-gar `epyc` is at
  `NC_090702.1:16496431-16523089` and encodes `XP_006633678.1`.
- Eleven exact non-focal gene symbols are shared in the 30-gene windows:
  `CEP290`, `DCN`, `DUSP6`, `KERA`, `LUM`, `MGAT4C`, `NTS`, `POC1B`, `RASSF9`,
  `RLIG1`, and `TMTC3`.
- All 11 shared genes retain a consistent order after orienting each block by
  the EPYC strand (`order_fraction = 1.000`). The opposite genomic direction is
  a whole-block inversion/orientation difference, not a loss of synteny.
- The closest SLRP cluster is conserved as human
  `EPYC-KERA-LUM-DCN` and spotted-gar `dcn-lum-kera-epyc` in ascending genomic
  coordinates. These become the same order when read relative to the focal
  gene's transcriptional orientation.
- `ATP2B1/atp2b1a` and `KITLG/kitlga` also occupy corresponding positions but
  are not counted by the conservative exact-symbol comparison.

The canonical spotted-gar locus therefore has strong manual synteny support.
Together with the reciprocal EPYC hit, EPYC-clade placement, intact LRR-domain
architecture, and strong SignalP result, this supports accepting
`XP_006633678.1` as spotted-gar EPYC. SynVoy's missing retained call is a
candidate-detection failure; its two provisional MEDIUM hits were correctly
removed as KERA- and LUM-related paralogs.
