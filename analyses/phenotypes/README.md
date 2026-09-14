# Curated phenotype association analysis

## Thesis question and scope

This module asks whether curated mouse mutant phenotypes and human
disease-phenotype annotations independently support skeletal, cartilage, joint,
growth, or extracellular-matrix relevance for the seven selected SLRPs. It is a
database evidence summary, not a cross-species phenotype-evolution analysis.

## Data sources

- Mouse Genome Informatics `MGI_GenePheno.rpt`, downloaded 2026-08-26.
- Mammalian Phenotype ontology release 2026-07-22
  (`MPheno_OBO.ontology`).
- Human Phenotype Ontology pinned release v2026-06-23,
  `genes_to_phenotype.txt` and `genes_to_disease.txt`.

The MGI parser retained all records whose genotype composition contained Bgn,
Dcn, Fmod, Prelp, Epyc, Lum, or Ogn. Single-gene and compound-genotype records are
separated. Broad categories were assigned from the term and its ontology
ancestors; categories may overlap. Counts represent annotation depth and must
not be interpreted as comparable effect sizes.

## Main mouse result

The analysis recovered 118 MGI records and 105 unique gene/MP-term/evidence-scope
combinations. All recovered focal records were single-gene genotypes in this
release. Unique single-gene MP terms were:

- BGN: 19 terms, including 14 skeletal/cartilage/joint terms;
- DCN: 24 terms, including 4 skeletal/cartilage/joint terms;
- FMOD: 8 terms, including 4 skeletal/cartilage/joint terms;
- PRELP: 12 terms, with no term assigned to the skeletal/cartilage/joint
  category in this release;
- EPYC: 3 terms, including short femur and osteoarthritis;
- LUM: 31 terms, dominated by corneal and skin/connective-tissue phenotypes,
  with one tendon phenotype in the skeletal category;
- OGN: 8 terms, dominated by corneal, skin/collagen, cardiovascular, and
  metabolic annotations, with no skeletal-category term recovered.

The phenotype layer therefore strengthens the direct skeletal interpretation
of BGN, DCN, FMOD, and EPYC. It supports broader ECM/connective-tissue roles for LUM
and OGN. PRELP's absence of a skeletal MP annotation does not negate the direct
growth-plate expression and experimental literature used elsewhere in the
thesis.

## Human result

The pinned HPO release contained 89 focal phenotype rows: 79 for BGN and 10 for
DCN. HPO mapped BGN to Mendelian OMIM:300106 and OMIM:300989 and DCN to
Mendelian OMIM:610048. No HPO gene-phenotype row was recovered for FMOD, PRELP,
EPYC, LUM, or OGN. This is absence of a curated
human disease annotation in that release, not evidence that those genes lack
human biological functions.

## Authoritative outputs

- `tables/mgi_focal_gene_phenotype_records.tsv`
- `tables/mgi_focal_gene_phenotype_terms.tsv`
- `tables/mgi_phenotype_category_summary.tsv`
- `tables/mgi_phenotype_evidence_summary.tsv`
- `tables/hpo_human_focal_gene_phenotypes.tsv`
- `tables/hpo_human_focal_gene_diseases.tsv`
- `tables/hpo_human_phenotype_evidence_summary.tsv`
- `figures/mgi_single_gene_phenotype_heatmap.png`

The tables and figure are reproduced with
`scripts/build_mgi_phenotype_summary.py` and
`scripts/summarize_hpo_focal_genes.py`.

## Recommended thesis placement

Use the heatmap as a supplementary figure or compact Results figure and place
the exact terms in a supplementary table. In the Discussion, use the phenotype
data as an independent biological evidence layer. Do not rank genes by the
number of database annotations, because heavily studied genes and phenotype
domains accumulate more records.
