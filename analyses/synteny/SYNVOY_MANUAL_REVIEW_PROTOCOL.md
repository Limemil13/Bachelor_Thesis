# SynVoy manual-review protocol

Last updated: 2026-09-10

Status: complete for the current reports. All 98 P1/P2/P3 rows have a recorded
decision; the three pending queues contain zero rows.

## What was automated

Seven final 14-target SynVoy reports were parsed from:

- `bgn_human_15species_dev_20260905`
- `dcn_human_15species_dev_20260903`
- `epyc_human_15species`
- `fmod_human_15species`
- `lum_human_15species`
- `ogn_human_15species_standard_1408`
- `prelp_human_15species_fresh`

The parser validates retained HIGH and MEDIUM counts against each report's final
post-ownership summary. Records labeled `paralog_not_goi` or marked
`coverage_demoted` are excluded. Low-confidence family matches remain ambiguity
evidence and are not promoted to orthologs.

Outputs:

- `tables/synvoy_gene_species_evidence.tsv`: all 98 gene-by-species rows.
- `tables/synvoy_manual_review_completed.tsv`: 25 completed P1 exceptions.
- `tables/synvoy_batch_confirmation_completed.tsv`: 69 completed P2 confirmations.
- `tables/synvoy_spot_check_completed.tsv`: four completed P3 controls.
- `diagnostics/priority_locus_reviews/all_synvoy_review_decisions.tsv`: complete 98-row decision source.
- `tables/synvoy_run_summary.tsv`: seven-run totals and report provenance.

## Review tiers

### P1: detailed review (26 rows; complete)

These rows have no retained HIGH/MEDIUM candidate or only a MEDIUM best call.
They are the loci most likely to change a biological conclusion.

For each row:

1. Open the exact target assembly and `best_start`--`best_end` interval in NCBI
   Genome Data Viewer or inspect the corresponding GFF3.
2. Determine whether the interval overlaps the expected gene, a named paralog,
   an uncharacterized/partial model, or no annotated model.
3. Record the confirmed symbol and protein accession when available.
4. Compare local order and orientation with the SynVoy SVG/HTML and the human
   flanking-gene block.
5. Check length, SLRP/LRR domains, SignalP, reciprocal protein support, and the
   gene-specific tree.
6. Enter `accepted`, `tentative`, `rejected`, or `ambiguous` with a concise
   reason in the master table.

No post-ownership HIGH/MEDIUM candidate was retained for BGN amphioxus,
catshark, chicken, elephant shark, or opossum; DCN amphioxus, spotted gar, or
zebrafish; EPYC spotted gar; FMOD zebrafish; OGN amphioxus; and PRELP amphioxus
or zebrafish. MEDIUM-only P1 rows are DCN coelacanth; EPYC amphioxus,
coelacanth, elephant shark, or zebrafish; FMOD amphioxus or catshark; LUM
amphioxus, catshark, coelacanth, or spotted gar; and PRELP catshark.

### P2: rapid confirmation (68 rows; complete)

These rows have a HIGH best call plus at least one automated exception, such as
multiple candidates, a rescue/partial model, a consistency flag, a deep-lineage
assembly, or a protein-QC flag. The canonical accession was located in the exact
target GFF and compared with the best post-ownership SynVoy coordinate. One
hidden conflict, spotted-gar FMOD, was escalated and resolved using the
alternative retained locus and a full flanking-gene comparison.

### P3: spot check (4 rows; complete)

BGN dog, mouse, and zebrafish plus EPYC mouse are grade-A positive controls:
each has one complete HIGH call, expected ownership, no automated exception,
and supported protein evidence. All four canonical accessions overlap the
corresponding SynVoy locus.

## Completed outcome

Across all 98 rows, 62 are accepted, 20 tentative, 12 ambiguous and four
rejected. All seven opossum rows are ambiguous because the target FASTA and GFF
sequence IDs are incompatible; they require a matched-input rerun. This is a
synteny limitation, not a reason to remove independently supported opossum
proteins from the MSA/tree analyses.

## What not to claim

SynVoy is supporting microsynteny evidence, not a one-ortholog-per-species truth
table. A retained HIGH call can still be a rescue or partial model, and a missing
call can reflect annotation, assembly, naming, contig, or window limitations.
Final orthology claims must integrate annotation, protein/domain evidence, and
the sequence tree. Re-running unchanged inputs is not a substitute for review.
