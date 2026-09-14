# Completed priority-locus review

Last updated: 2026-09-10

This directory contains the reproducible curation layer for all 98 SynVoy
gene-by-species combinations. The 25 P1 rows received detailed locus review;
the remaining 69 P2 rows and four P3 controls received an exact-GFF canonical
accession/locus confirmation. All generated pending queues are now empty. The
complete decisions are stored in `all_synvoy_review_decisions.tsv` and merged
into `../../tables/synvoy_gene_species_evidence.tsv`.

## Detailed P1 outcome

| Decision | Rows | Interpretation |
|---|---:|---|
| accepted | 7 | The canonical locus, flanking block and protein/tree evidence agree. |
| tentative | 8 | The ortholog or co-ortholog is plausible, but one important weakness remains. |
| ambiguous | 6 | The available assembly/annotation cannot resolve a gene-specific locus. |
| rejected | 4 | The SynVoy interval belongs to another gene, or no gene-specific evidence supports it. |

The accepted calls are EPYC in coelacanth, elephant shark and spotted gar;
LUM in catshark, coelacanth and spotted gar; and PRELP in zebrafish.

The rejected calls are chicken BGN (`XP_414298.2`, ASPN-like), amphioxus EPYC
(RUN-domain protein-like locus), amphioxus FMOD (BRD3-like locus) and amphioxus
OGN (all raw calls reassigned to other genes). Rejection applies to the tested
candidate assignment; it is not a general claim that the entire species lacks
every homologous SLRP.

## Full-panel outcome

Across all 98 gene-by-species rows, the completed review records 62 accepted,
20 tentative, 12 ambiguous and four rejected calls. Sixty-seven rows were
confirmed by direct overlap between the canonical protein accession and the
best post-ownership SynVoy locus. The remaining rows were resolved by detailed
review, except that all seven opossum gene rows remain ambiguous pending a
clean target rerun.

The rapid canonical-locus audit exposed one additional hidden conflict:
spotted-gar FMOD. SynVoy's top HIGH interval spans lumican-like/PRELP models,
but a second retained HIGH rescue overlaps canonical `fmoda`
(`XP_006628521.2`). Eight human flanking genes are shared and six preserve
order, so the canonical `fmoda` locus is accepted while the summary's top
coordinate is not treated as FMOD.

Elephant-shark `LOC103182549` (`XP_007897806.2`) is retained only as tentative
FMOD-like provenance. Its locus has partial broad-block support, but the GFF
reconstructs a complete 684-aa, nine-CDS product that is probably a compound
annotation. BLASTP assigns residues 1--330 most strongly to DCN and residues
376--684 most strongly to FMOD. Pfam detects LRRs across both halves and two
LRR N-terminal caps, while residues 330--347 form a second hydrophobic,
signal-peptide-like segment. The locus remains useful FMOD evidence, but the
compound protein is excluded from the confident FMOD sequence panel.

## Method

The review used the exact FASTA/GFF target pair supplied to SynVoy. For each
canonical locus, the existing `compare_gff_neighborhoods.py` script extracted
gene windows from the human and target GFFs, oriented both windows relative to
the gene of interest, and counted exact shared flanking-gene symbols and their
relative order. The default comparison used 15 genes on each side. Catshark
BGN was broadened to 40 genes on each side because the local annotation
contains many `LOC` records.

Exact-symbol matching is deliberately conservative. A low count can reflect
symbol changes, lineage-specific duplicates, `LOC` annotations or assembly
fragmentation. Therefore, the final decision also used:

- whether the canonical RefSeq protein maps to the candidate locus;
- reciprocal protein similarity and query coverage;
- gene-specific phylogenetic placement;
- domain and SignalP support;
- the identity of any gene actually overlapping the SynVoy interval; and
- known duplicate structure, especially for teleost genes.

The scripts and exact intermediate tables are retained so every decision can
be regenerated rather than depending on screenshots.

## Main findings

1. The vertebrate class-III block is strongly conserved. Coelacanth,
   spotted-gar and zebrafish comparisons preserve combinations of
   `EPYC-KERA-LUM-DCN` plus distal anchors such as `BTG1`, `CEP290`, `TMTC3`,
   `DUSP6`, `POC1B` and `ATP2B1`.
2. SynVoy's retained fragment is not automatically the canonical ortholog.
   Canonical EPYC, LUM and DCN loci were sometimes missed or a paralogous
   fragment received the provisional GOI label. The final call therefore uses
   the post-ownership report plus independent locus/protein evidence.
3. Zebrafish contains two plausible FMOD co-orthologs. `fmodb`
   (`XP_073767056.1`) is the better full-length human-FMOD match (50.9% identity,
   99% query coverage), while `fmoda` (`NP_001025243.1`) has higher local
   identity (56.6%, 80% query coverage) and retains the ancestral
   FMOD-PRELP neighbourhood. Both must be discussed as duplicate evidence.
4. Deep-lineage Branchiostoma calls are frequently not gene-specific. The
   reviewed EPYC and FMOD intervals overlap unrelated annotated genes, LUM is
   intergenic, and all OGN raw candidates were reassigned. The manually reviewed
   OGN sequence `XP_066265713.1` is from *Branchiostoma lanceolatum*, whereas
   the SynVoy target is *B. floridae*; the two must not be treated as one locus.
5. Chicken `XP_414298.2` has an annotation conflict: the GFF label is BGN, but
   its product, focused tree and six locus neighbours support ASPN. It remains
   excluded from the canonical BGN panel.

## Critical input warning: opossum

The opossum SynVoy target is not technically interpretable. The GFF uses
chromosome sequence IDs including `NC_077235.2`, while the paired FASTA uses
older IDs including `NC_077235.1`. Only 4 of 13 annotated sequence IDs match
exactly; six differ by version and three annotated scaffolds are absent from
the FASTA. All opossum SynVoy scores must be excluded until the target is rerun
with a synchronized FASTA/GFF pair.

This warning concerns the SynVoy locus layer only. Independently supported
opossum proteins can remain in protein/MSA/tree analyses, but synteny must not
be claimed for them from the current run.

The audit is in `fasta_gff_seqid_audit_summary.tsv` and
`fasta_gff_seqid_audit_mismatches.tsv`. The validator was strengthened so this
type of mismatch is detected before future runs.

## Key files

| File or directory | Purpose |
|---|---|
| `manual_review_decisions.tsv` | Curated decisions and full evidence notes for all 25 P1 rows. |
| `all_synvoy_review_decisions.tsv` | Complete decision source for all 98 P1/P2/P3 rows. |
| `canonical_locus_audit.tsv` | Exact-GFF canonical-accession versus SynVoy-coordinate audit for all rows. |
| `synvoy_raw_candidates.tsv` | Coordinates and ownership outcomes of all raw HIGH/MEDIUM candidates. |
| `fasta_gff_seqid_audit_summary.tsv` | Cross-target FASTA/GFF identifier compatibility check. |
| `<GENE>_<species>/pairwise_comparison.tsv` | Reproducible human-target anchor and order metrics. |
| `<GENE>_<species>/neighborhoods.tsv` | Extracted oriented GFF neighbourhoods. |
| `*_features.tsv` / `*_neighbors.tsv` | Exact feature overlap and local context for suspicious fragments. |
| `FMOD_zebrafish_blastp.tsv` | Human FMOD comparison supporting the two zebrafish co-orthologs. |
| `BGN_opossum_gap_tblastn.tsv` | Focused no-hit result on the older FASTA interval; not decisive because the inputs mismatch. |

## Rebuild

Run the reusable locus comparisons from the repository root with:

```bash
bash analyses/synteny/scripts/run_priority_locus_comparisons.sh
```

Rebuild the master review tables with:

```bash
SYNVOY_ROOT=/path/to/SynVoy
python3 analyses/synteny/scripts/build_synvoy_review.py \
  --synvoy-results "$SYNVOY_ROOT/results" \
  --conservation-table analyses/protein_analysis/tables/protein_conservation_domain_msa_signalp_summary.tsv \
  --protspace-table analyses/protein_analysis/protspace/tables/protspace_canonical_embedding_qc.tsv \
  --input-audit analyses/synteny/diagnostics/priority_locus_reviews/fasta_gff_seqid_audit_summary.tsv \
  --manual-decisions analyses/synteny/diagnostics/priority_locus_reviews/all_synvoy_review_decisions.tsv \
  --output-dir analyses/synteny/tables
```

After rebuilding, all three queue files should contain zero data rows;
`synvoy_manual_review_completed.tsv`, `synvoy_batch_confirmation_completed.tsv`
and `synvoy_spot_check_completed.tsv` should contain 25, 69 and four rows,
respectively.
