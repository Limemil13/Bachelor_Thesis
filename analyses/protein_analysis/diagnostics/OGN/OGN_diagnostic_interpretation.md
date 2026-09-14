# OGN diagnostic interpretation

## Summary

OGN showed warnings in the MSA quality check and SignalP analysis. The issue is not global for all OGN candidates; most OGN proteins have normal length and a predicted signal peptide.

## Dog OGN

Candidate:
OGN_canis_lupus_familiaris_XP_038383338.1

Evidence:
- BLASTP identity: 90.60%
- Query coverage: 100%
- Reciprocal best hit: human NP_148935.1 mimecan / OGN
- Protein header: XP_038383338.1 mimecan [Canis lupus familiaris]
- Full protein length: 361 aa
- SignalP full sequence: OTHER / no signal peptide
- A normal OGN-like start motif, MKTLQST, occurs at position 65.
- Trimming the sequence to this motif gives a 297 aa protein.
- SignalP on the trimmed sequence predicts a strong Sec/SPI signal peptide.

Interpretation:
The dog OGN candidate is likely the correct OGN/mimecan ortholog, but the protein model probably has an incorrect N-terminal extension or wrong annotated start site. This explains the missing SignalP prediction and the abnormal MSA gap behavior.

Recommended handling:
Keep dog OGN as an ortholog candidate, but mark the original protein model as N-terminally misannotated. For diagnostic/corrected alignment, use the trimmed version starting at MKTLQST.

## Spotted gar OGN

Candidate:
OGN_lepisosteus_oculatus_XP_006631165.2

Evidence:
- BLASTP identity: 69.86%
- Query coverage: around 73–74%
- Reciprocal best hit: human NP_148935.1 mimecan / OGN
- Protein header: XP_006631165.2 osteoglycin, paralog a [Lepisosteus oculatus]
- Protein length: 295 aa
- SignalP prediction: OTHER / no signal peptide

Interpretation:
The spotted gar candidate has BLAST and reciprocal BLAST support and is annotated as osteoglycin paralog a, but lacks a SignalP-predicted signal peptide. Unlike dog, this is not caused by a long N-terminal extension. It may represent a divergent N-terminus, uncertain cleavage site, annotation issue, or a paralog-specific difference.

Recommended handling:
Keep as candidate but mark as SignalP-negative / needs validation. Check whether an alternative OGN-like candidate exists in spotted gar and later validate by synteny.
