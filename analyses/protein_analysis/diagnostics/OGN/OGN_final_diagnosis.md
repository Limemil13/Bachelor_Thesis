# OGN final diagnostic interpretation

## Main conclusion

OGN is not globally problematic. Most OGN candidates have normal protein length, expected SLRP/LRR domain architecture, and predicted N-terminal signal peptides.

The main issues are candidate-specific:
- dog OGN has an incorrect N-terminal extension
- spotted gar OGN lacks a SignalP-predicted signal peptide despite good BLAST support

## Dog OGN

Candidate:
OGN_canis_lupus_familiaris_XP_038383338.1

Evidence:
- Protein header: XP_038383338.1 mimecan [Canis lupus familiaris]
- BLASTP identity to human OGN: 90.60%
- Query coverage: 100%
- Reciprocal best hit: human NP_148935.1 mimecan / OGN
- Original protein length: 361 aa
- Normal OGN-like start motif MKTLQST begins after 64 aa
- Trimmed protein length from MKTLQST: 297 aa
- SignalP full protein: no signal peptide
- SignalP trimmed protein: strong Sec/SPI signal peptide, cleavage site 19-20

Interpretation:
The dog OGN candidate is very likely the correct OGN/mimecan ortholog, but the predicted protein has an incorrect N-terminal extension or wrong annotated start site. This explains both the SignalP failure and the abnormal alignment behavior.

Handling:
Keep dog OGN as candidate, but mark original protein model as N-terminally misannotated. For diagnostic/corrected alignment, the trimmed version starting at MKTLQST can be used.

## Spotted gar OGN

Candidate:
OGN_lepisosteus_oculatus_XP_006631165.2

Evidence:
- Protein header: XP_006631165.2 osteoglycin, paralog a [Lepisosteus oculatus]
- BLASTP identity to human OGN: 69.86%
- Query coverage: around 73%
- E-value: 2.24e-109
- Bitscore: 320
- Reciprocal best hit: human NP_148935.1 mimecan / OGN
- SignalP: no signal peptide predicted
- Protein length: 295 aa

Interpretation:
The spotted gar candidate is likely the best available OGN/osteoglycin candidate in the local proteome, but its N-terminal signal peptide is not detected by SignalP. This may reflect an unusual/divergent N-terminus, uncertain cleavage site, or annotation issue.

Handling:
Keep the spotted gar candidate, but mark it as SignalP-negative and requiring synteny validation.
