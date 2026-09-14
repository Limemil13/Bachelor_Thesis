# Comparative SLRP panel: current status

Last updated: 2026-09-07

The historical filename is retained for existing links. The current full panel
has **seven genes**: BGN, DCN, FMOD, PRELP, EPYC, LUM, and OGN. OMD is a context
gene. Candidate selection integrates juvenile-lineage expression, direct
mouse/rat growth-plate expression, literature, SLRP-class coverage, orthology,
and technical completeness; it is not a rank by RNA abundance alone.

## Current decision

- Principal genes: BGN, DCN, FMOD, PRELP, EPYC, and LUM.
- Exploratory main gene: OGN.
- Context only: OMD.

The local protein package is a verified 103-protein set. The former 15-protein
BGN set was decontaminated to 13 after removal of DCN/DCN-like assignments, and
DCN was promoted as a full 15-protein member and paralog control.

## Important qualifications

- BGN does not form one pure split in the combined seven-gene tree, although a
  focused BGN/DCN diagnostic correctly separates all tips by nearest same-gene
  neighbour. Catshark and whale-shark BGN are partial/tentative.
- Whale-shark EPYC is retained as tentative.
- Spotted-gar OGN is retained as tentative and is SignalP negative.
- Amphioxus OGN is excluded from the confident set and retained as uncertain
  SLRP-like provenance.
- Lamprey LUM is outside the 15/16 main LUM split and remains tentative.
- Seventeen SignalP calls (15 DCN and two BGN replacements) are pending.
- SynVoy is complete for all seven panel genes, including updated-tool DCN and
  BGN runs. The optional updated-tool FMOD refresh was interrupted before its
  final report; PRELP and EPYC refreshes have not been launched. Historical
  completed reports remain the current evidence for those three genes.

## Does juvenile RNA make sense?

Yes. A growth plate is a developmental cartilage, so juvenile data are
biologically appropriate. The limitation is anatomical specificity. Local P0
Prx1-lineage samples are a relevant screen, whereas GSE114919 microdissected
zones provide direct growth-plate support. The two scales must not be pooled or
described as one differential-expression experiment.

## Does the evolutionary comparison make sense?

Yes. It tests conservation, paralog separation, incomplete models, and locus
continuity. The informative outcome is a conserved vertebrate ECM framework
plus specific exceptions—not a claim that all orthologs have identical
growth-plate function.

For exact counts and remaining work, use `THESIS_REPRODUCIBILITY_STATUS.md` and
`THESIS_NEXT_STEPS.md`. For result paths, use
`AUTHORITATIVE_RESULT_LOCATIONS.md`.
