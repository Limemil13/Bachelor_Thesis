# Mouse RNA expression interpretation

The focal RNA-seq samples are three wild-type samples from GSE305415
(PRJNA1305489): P0 FACS-sorted Prx1-lineage chondroprogenitors isolated from
digested knee joints. They should not be described as bulk whole-growth-plate
samples. WT1 and WT2 are paired-end; WT3 is single-end.

Data provenance: GSE305415 is a public NCBI Gene Expression Omnibus (GEO)
study; its raw sequencing reads are available through NCBI SRA under
BioProject PRJNA1305489. The focal WT runs used here are SRR34978216,
SRR34978215, and SRR34978214. This project re-extracted expression values from
previously generated Salmon `quant.sf` files stored on the computer. "Local"
describes those working files, not the origin of the samples or sequencing data.
The original Salmon index-building/quantification command was not preserved.
Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305415

The four files previously labeled as skin, heart, kidney, and brain are from
the aging atlas PRJNA936435. Their verified samples are skin at 6 months,
heart at 30 months, kidney at 30 months, and skin at 15 months. There is no
brain sample: `SRR23517588` was mislabeled locally and is `Skin_15mo_2`.

Because these comparators differ in study, age, tissue/cell composition,
library context, and replication, they are retained only as a descriptive
panel. They must not be used for formal differential expression or as a clean
estimate of growth-plate tissue specificity. Candidate prioritization should
use reproducible expression in the three focal samples, cartilage-marker
co-expression, prior biological relevance, and independent orthology/protein
evidence. A matched developmental cartilage or growth-plate dataset is needed
for a defensible tissue-enrichment comparison.

## Independent growth-plate validation added

GSE114919 provides a directly relevant validation design: laser-capture
microdissected proliferative (PZ) and hypertrophic (HZ) zones from mouse and
rat growth plates, with tibia sampled at 1 and 4 weeks and five replicates per
tibial condition. The published processed workbooks and all derived tables are
under `gse114919/`.

The authors' normalized values are retained on their original scale. They are
not pooled with the GSE305415-derived Salmon TPMs, and absolute values are not
compared between species. Cross-species interpretation uses detection and
within-condition rank.

Across tibial conditions, BGN, FMOD, and PRELP are the strongest and most
reproducible main-panel genes. EPYC is consistently mid-ranked. OGN is present
in all mouse groups but weak/partly absent in rat. LUM and DCN are consistently
expressed comparators. OMD and ASPN have zero values in all tested rat tibial
conditions.

The mouse processed matrix contains 29 rather than 30 sample columns because
1-week phalanx PZ replicate 3 is absent. All tibial conditions remain n=5, so
the missing phalanx column does not affect the thesis validation.

Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114919

## Bgee support

The Bgee 16 workflow under this folder retrieves and preserves raw JSON for
nine SLRPs in mouse, chicken, and zebrafish. Its output is useful as curated
qualitative support, but the data types, stages, and anatomical coverage are
not matched across species. Absence of a Bgee call is treated as missing
database evidence, not as biological absence.

## One-Embedding decision

One-Embedding was evaluated and is not an RNA-expression method. It is a
protein-language-model embedding compression/transport codec. Its current
residue-conservation utility is explicitly heuristic and has no trained probe,
so it is not used as conservation evidence. The thesis should use the existing
MSA-derived conservation profiles instead.

Repository reviewed: https://github.com/jcoludar/One-Embedding
