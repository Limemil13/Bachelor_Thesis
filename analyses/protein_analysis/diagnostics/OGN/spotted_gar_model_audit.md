# Spotted-gar OGN model audit

The exact spotted-gar annotation used by SynVoy contains one transcript at the
`ogna` locus:

- assembly: `GCF_040954835.1`
- locus: `NC_090699.1:61481801-61491048`
- transcript: `XM_006631102.3`
- protein: `XP_006631165.2`
- model source: NCBI Gnomon
- coding structure: six CDS blocks

No alternative annotated transcript or alternative N terminus is present in
this GFF. The protein begins
`MSKLKPLIFSLILVLWVILAAVTAF`; residues 8--25 form a strongly hydrophobic
N-terminal segment compatible with a signal peptide or signal-anchor-like
region. SignalP 6 nevertheless classified the sequence as `OTHER` with an SP
score of 0.167 and did not assign a cleavage site.

The sequence is retained as tentative OGN. Its same-locus HIGH SynVoy call,
reciprocal protein match, pure OGN-tree placement, continuous LRR-rich
alignment, and hydrophobic N terminus support the assignment. The absence of an
annotated alternative model and the negative formal SignalP call prevent an
upgrade to a confident secretion claim. A second validated secretion/topology
predictor or experimental N-terminal evidence would be needed to resolve this
specific uncertainty.
