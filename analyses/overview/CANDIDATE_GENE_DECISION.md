# Candidate-gene decision

Last updated: 2026-09-14

## 2026-09-02 audit update

DCN is promoted from comparator to the complete comparative package. Its
15-species protein panel, Pfam scan, MSA/conservation outputs, compact and
combined trees, five-species gene-structure analysis, coding constraint,
ProtSpace, protein motifs, promoter analysis, and phenotype summaries are
complete. SignalP is complete for all 15 DCN proteins (14 positive; the
N-terminally incomplete opossum fragment is negative). Its updated SynVoy run
and detailed locus review are integrated.

The same audit found four historical BGN entries that are DCN or DCN-like
records. The proposed 13-protein BGN panel was approved and the clean
103-protein downstream package has been rebuilt. Full accessions, replacements,
evidence, and affected outputs are recorded in
`history/DCN_PROMOTION_AND_BGN_AUDIT_2026-09-02.md`.

## Recommendation

Use **BGN, DCN, FMOD, PRELP, EPYC, OGN, and LUM** as the expanded seven-gene
panel, but do not describe them as the seven most highly or specifically
expressed SLRPs. Describe them as an **expression-supported, class-spanning panel
of juvenile cartilage/growth-plate-associated SLRPs** selected using several
evidence types:

- reproducible expression in three P0 Prx1-lineage chondroprogenitor samples
  from public NCBI GEO GSE305415 (re-extracted from stored Salmon files);
- independent, zone-resolved mouse and rat growth-plate expression;
- prior cartilage, growth-plate, or skeletal evidence;
- representation of SLRP classes I (BGN), II (FMOD, PRELP, and LUM), and III
  (EPYC and OGN); and
- feasibility/completeness of the existing cross-vertebrate protein,
  phylogeny, and synteny work.

DCN began as an established expression/biology comparator, but the 2026-09-02
audit showed that it merits the complete comparative package and is required to
control BGN/DCN paralog assignments. Keep ASPN and OMD as secondary context.

Suggested Methods wording:

> Seven SLRP genes (BGN, DCN, FMOD, PRELP, EPYC, OGN, and LUM) were prioritized as an
> expression-supported, class-spanning panel associated with juvenile
> cartilage and growth-plate biology. Selection integrated reproducible
> expression in P0 Prx1-lineage chondroprogenitors, independent rodent
> growth-plate evidence, prior skeletal biology, and suitability for
> cross-vertebrate orthology and protein analysis; it was not based solely on
> expression magnitude.

## Why this is defensible

The independent GSE114919 analysis uses laser-capture-microdissected
proliferative and hypertrophic growth-plate zones from mouse and rat, at 1 and
4 weeks, with five tibial replicates per condition. BGN, FMOD, and PRELP are
the strongest and most reproducible members across these conditions. EPYC is
consistently mid-ranked. These data directly repair the tissue/age mismatch in
the old adult-tissue comparison.

BGN also has an in-vivo skeletal phenotype, FMOD is directly localized in the
growth plate, PRELP is a validated proliferative-zone marker, and EPYC is
expressed in epiphyseal/growth-plate cartilage and has a mild knockout
phenotype. LUM is reproducibly expressed in both rodent growth-plate datasets,
has a conserved two-CDS-exon structure, and completed the cleanest protein QC
profile in the panel (16/16 domain- and SignalP-supported proteins). Source
papers are linked in
`analyses/overview/candidate_gene_evidence.tsv`.

## The important exception: OGN

OGN is the least certain of the six as a cross-species growth-plate gene:

- it is the second most abundant of the nine SLRPs in the local P0 samples;
- it is detected in all four mouse tibial GSE114919 conditions, but ranks
  below BGN, FMOD, PRELP, LUM, EPYC, and DCN;
- it is low and not detected in every rat replicate;
- it has the largest number of protein-QC exceptions in the completed
  canonical panel; and
- its strongest literature rationale is bone biology and informative SLRP
  class-III evolution rather than a clean growth-plate-specific knockout
  phenotype.

Therefore OGN should be retained as an **exploratory evolutionary/class-III
candidate**, not used as evidence that the panel is simply the top
growth-plate-expression set.

LUM has now been added to the full protein, domain, MSA, SignalP, ProtSpace,
gene-structure, compact-tree, large-tree, and iTOL workflows. If the supervisor
requires a strictly growth-plate-focused final figure set, omit OGN from that
specific headline comparison and retain it as an exploratory class-III
evolutionary case; do not remove its completed analyses from the thesis record.

## Chicken BGN annotation caveat

BGN remains the class-I anchor, but chicken `XP_414298.2` is not part of the
canonical BGN protein panel. Although its current transcript is assigned to a
BGN locus, NCBI names the product asporin, reciprocal protein comparison selects
human ASPN, and a focused SLRP-reference tree places it with ASPN at 99.3
SH-aLRT/100 UFBoot. Report this as a chicken locus/product annotation conflict
and resolve it with neighborhood/synteny evidence if a species-complete BGN
claim is needed. It does not weaken the biological rationale for the retained
13-protein BGN panel, but the partial shark replacements remain tentative.

## Why OMD should not be added now

OMD is biologically relevant to bone and subchondral homeostasis, but it is not
a stronger fit for this thesis than the expanded seven-gene panel:

- local P0 expression is reproducible but lowest among the nine tested SLRPs;
- mouse growth-plate values are relatively low;
- the rat GSE114919 growth-plate matrix reports zero values in all tested
  conditions;
- recent Omd-null work found altered transverse cortical-bone size without an
  effect on longitudinal bone growth; and
- OMD currently lacks the completed 15-species protein, SignalP, ProtSpace,
  tree, and SynVoy layers available for the main panel.

Use OMD in the introduction/discussion or a compact comparator table. Do not
add it to gene structure alone and then imply that an eight-gene pipeline is
complete.

## Does juvenile RNA make sense for a growth-plate thesis?

Yes. Growth plates are active developmental structures, so juvenile data are
more relevant than adult tissue for a question about longitudinal bone growth.
The limitation is not the young age; it is anatomical specificity.

The local P0 samples are FACS-sorted Prx1-lineage chondroprogenitors from
digested knee joints, not isolated growth-plate zones. They support the claim
that the genes are expressed in the relevant juvenile lineage, but not a claim
of growth-plate zone specificity or tissue enrichment. GSE114919 supplies the
missing direct growth-plate/zone validation. Bgee remains qualitative support
because its species, stages, and data types are not matched consistently.

## Source-of-truth table

The integrated nine-gene evidence table is:

`analyses/overview/candidate_gene_evidence.tsv`

The ranking is intentionally transparent rather than collapsed into an
arbitrary single numerical score.
