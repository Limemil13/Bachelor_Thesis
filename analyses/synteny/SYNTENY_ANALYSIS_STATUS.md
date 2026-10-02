# SLRP synteny analysis status

Last audited: 2026-10-02

## Decision

SynVoy reports are complete and parsed for BGN, DCN, EPYC, FMOD, LUM, OGN, and
PRELP in the 98-row master table. The updated-tool DCN run completed on
4 September 2026 at `results/dcn_human_15species_dev_20260903`, and the updated
BGN run completed on 6 September 2026 at
`results/bgn_human_15species_dev_20260905`. The updated FMOD report at
`results/fmod_human_15species_dev_20260906` has also completed, been compared
with the historical report, and been integrated into all 14 FMOD rows. The
PRELP refresh has a resumable checkpoint after iterative-search wave 4 of 7 but
no final report, and an updated-tool EPYC run has not been launched. Historical
PRELP and EPYC reports therefore remain the sources for those rows; these two
refreshes are optional sensitivity analyses rather than requirements for the
thesis. OGN and LUM were produced with the newer workflow and are not part of
this refresh batch.

All 98 current gene-by-target rows have now been reviewed. Historical calls
remain useful for comparison, and each future updated report should be compared
with its older counterpart before replacing corresponding rows.

The DCN report contains 10 HIGH and one MEDIUM retained post-filter candidate.
Amphioxus, spotted gar, and zebrafish have no retained post-filter DCN candidate;
coelacanth has the MEDIUM call. Detailed review recovered canonical ordered
coelacanth, spotted-gar, and zebrafish loci and left amphioxus unresolved; none
of these automated misses is evidence that DCN is absent.

The updated BGN report contains 10 HIGH and one MEDIUM retained annotation,
with only nine self-consistency flags compared with 36 in the historical run.
Best retained calls are HIGH in nine genomes. Amphioxus, catshark, chicken,
elephant shark, and opossum have no retained post-filter BGN candidate. Review
retained catshark tentatively, rejected the chicken ASPN-like record, and left
amphioxus, elephant shark, and opossum ambiguous rather than inferring losses.

The LUM run completed successfully on 1 September 2026 after 14 h 12 min 58 s.
All 14 target genomes passed staging/QC, and all 166 workflow tasks completed.
The final report is stored locally at:

`$SYNVOY_ROOT/results/lum_human_15species/synvoy_report.json`

The run used the human LUM protein and locus, the same 14 target genomes as the
other panel genes, LUM/lumican family tokens, non-strict family matching,
automatic generic presets disabled, an 8 GB MMseqs split-memory limit, and one
iterative-search CPU. These settings are recorded in the run's
`intermediate/locate_gene/effective_params.json` and Nextflow log.

## Final automated summary

| Gene | HIGH retained candidates | MEDIUM retained candidates | Best-HIGH genomes | P1 detailed rows | P2 rapid rows | P3 control rows |
|---|---:|---:|---:|---:|---:|---:|
| BGN | 10 | 1 | 9 | 5 | 6 | 3 |
| DCN | 10 | 1 | 10 | 4 | 10 | 0 |
| EPYC | 9 | 15 | 9 | 5 | 8 | 1 |
| FMOD | 10 | 5 | 10 | 4 | 10 | 0 |
| LUM | 10 | 39 | 10 | 4 | 10 | 0 |
| OGN | 13 | 22 | 13 | 1 | 13 | 0 |
| PRELP | 11 | 11 | 11 | 3 | 11 | 0 |
| **Total** | **73** | **94** | **72** | **26** | **68** | **4** |

HIGH/MEDIUM values are retained annotations/candidates, not counts of verified
one-to-one orthologs. A genome can contain more than one retained candidate.
The best-HIGH column counts genomes whose top retained call is HIGH.

LUM has a best HIGH call in 10 genomes and a best MEDIUM call in amphioxus,
catshark, coelacanth, and spotted gar. Most LUM HIGH calls are rescue-derived
and/or accompanied by multiple candidates or self-consistency flags. This makes
the run useful for locating and prioritizing loci, but not sufficient by itself
to claim 14 confirmed LUM orthologs. Completed review accepted the canonical
catshark, coelacanth, and spotted-gar LUM loci and left amphioxus ambiguous.

## Completed review summary

| Decision | Gene-by-target rows |
|---|---:|
| accepted | 62 |
| tentative | 20 |
| ambiguous | 12 |
| rejected | 4 |

Sixty-seven canonical protein accessions overlap the best post-ownership
SynVoy candidate in the exact target GFF. Detailed review resolved the remaining
P1 exceptions and two additional FMOD cases. The rejected rows are chicken BGN
(ASPN-like product) and amphioxus EPYC, FMOD, and OGN assignments. Spotted-gar
FMOD was accepted at an alternative retained `fmoda` locus rather than the
summary's top interval. Elephant-shark FMOD remains tentative because its
fibromodulin-like 684-aa model contains two SLRP-like halves: a DCN/BGN-like
N-terminal half and an FMOD-best C-terminal half. The locus is retained, but
the likely compound/fused protein model is excluded from the confident panel.

All seven opossum rows are technically ambiguous: the local GFF and FASTA use
incompatible sequence-accession versions. No current opossum SynVoy score may
be cited until the target is rerun with a synchronized FASTA/GFF pair.

## What remains

1. Replace the opossum target with a synchronized FASTA/GFF pair and rerun its
   affected comparisons before using opossum synteny in the thesis.
2. The outstanding SignalP calls and focused BGN/DCN protein review are now
   complete. The elephant-shark FMOD-like protein remains triaged as a likely
   compound/fused annotation; resolving it requires a revised model or
   transcript evidence rather than another routine domain run.
3. PRELP and EPYC updated-tool reruns may be completed as sensitivity analyses,
   but the documented historical reports are sufficient for the present thesis.
4. Rebuild the evidence tables and repeat locus review only for rows changed by
   a new report. Do not interpret a missing displayed gene as gene loss.

No additional broad synteny method is required for the bachelor thesis. If a
specific locus remains unresolved after NCBI/GFF review, a targeted protein-to-
genome search or a better assembly is justified for that locus only.

## Interpretation boundary

SynVoy is qualitative microsynteny and locus-prioritization evidence. Final
orthology claims must combine locus annotation and flanking genes with reciprocal
protein evidence, SLRP/LRR architecture, SignalP, and phylogenetic placement.
