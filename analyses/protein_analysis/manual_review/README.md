# Manual sequence review: priority proteins

The packet in this directory limits manual work to the five proteins with the
integrated status `Needs inspection`:

- EPYC whale shark `XP_048463897.1`
- OGN opossum `XP_056660002.1`
- OGN zebrafish `NP_001013588.1`
- OGN spotted gar `XP_006631165.2`
- OGN amphioxus `XP_066265713.1`

## Files

- `manual_sequence_review_5.faa`: the five unaligned protein sequences.
- `manual_sequence_review_5.tsv`: all protein/domain/MSA/SignalP QC fields.
- `EPYC_manual_review_alignment_subset.faa`: whale-shark EPYC plus human and
  its two closest same-gene compact-tree neighbors.
- `OGN_manual_review_alignment_subset.faa`: the four OGN targets plus human,
  mouse, and the closest OGN compact-tree tip to amphioxus.
- `manual_sequence_review_5_tree_neighbors.tsv`: five closest compact-tree tips
  for every target.
- `manual_sequence_review_5_synteny_context.tsv`: the matching SynVoy evidence,
  candidate counts, P1/P2 tier, model status, and review reasons for each target.
- `manual_sequence_review_decisions.tsv`: fill in your final decision and one
  sentence of evidence for each target; automated starting interpretations are
  provided but are not substitutes for the visual review.
- `focused_bgn_dcn_manual_review.tsv`: reproducible alignment, SignalP, domain,
  tree, motif, and synteny evidence for the four later BGN/DCN review targets.

## Inspect in this order

1. Open each alignment subset in AliView/Jalview and enable residue coloring.
2. Check that the target has an intact N terminus and does not begin or end in a
   conspicuously truncated, gap-rich region.
3. Check the cysteine pattern and the central LRR-rich region against the human
   and nearest-species comparators. An isolated long insertion or deletion is a
   gene-model warning.
4. Check whether low identity/query coverage is distributed across the protein
   (plausible divergence) or confined to a suspicious terminal/additional block.
5. For spotted-gar OGN, inspect the N terminus especially carefully because it is
   the only sequence without a SignalP-positive result.
6. For amphioxus OGN, treat orthology as unconfirmed even if the sequence looks
   SLRP-like: it falls outside the 14-tip OGN split and most nearest tree tips are
   other SLRP genes. This candidate needs the strongest annotation/synteny check.
7. Use the synteny-context table as corroborating evidence: the four P2 targets
   have a retained HIGH candidate but still need model/sequence confirmation;
   amphioxus OGN is P1 and has no retained HIGH/MEDIUM GOI candidate.

Record one of: `retain`, `retain as tentative`, `replace`, or `exclude`, plus one
sentence of evidence. Do not edit the source-of-truth summary directly; record
the decision first, then update the candidate manifest and rerun the canonical
refresh only if a sequence is replaced or excluded.

## Review status on 2026-08-25

Four OGN decisions have been recorded from manual AliView/Jalview inspection:
opossum and zebrafish are retained, spotted gar is retained as tentative because
SignalP support is absent, and amphioxus is excluded from the confident OGN
ortholog set and retained only as an uncertain SLRP-like provenance candidate.

Whale-shark EPYC `XP_048463897.1` was retained as tentative after manual
inspection. Its N terminus begins similarly to the comparison sequences and the
middle-to-C-terminal region is strongly conserved, but numerous mismatches in
the earlier portion, reduced query coverage, and the alternative reciprocal hit
prevent an unqualified confident call.

The applied visual rule did not require exact terminal amino-acid identity.
Candidates were assessed for intact termini, conserved cysteine positions,
continuous LRR-rich alignment, absence of large target-specific indels, and
agreement with tree, synteny, domain, and SignalP evidence.

## Focused BGN/DCN review completed on 2026-09-14

The two partial BGN models and two exceptional DCN models were reviewed after
the missing SignalP run was completed. Catshark BGN is retained as a tentative
partial ortholog because its secretion signal and N-terminal/LRR structure are
supported but its C-terminal region is missing. Whale-shark BGN is retained as
a tentative locus-associated fragment and is excluded from claims about a
complete secreted BGN protein. Opossum DCN is retained as a high-identity
C-terminal DCN fragment, but it cannot support secretion or complete-domain
claims and its current synteny result remains unusable. Catshark DCN is retained
as a confident ortholog: full-span MSA, domain, motif, tree, and locus evidence
agree, while the predicted 45--46 SignalP cleavage site is reported as an
N-terminal model caveat rather than hidden.
