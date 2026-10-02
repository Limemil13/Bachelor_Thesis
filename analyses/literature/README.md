# SLRP literature evidence

This directory stores citation metadata and short, original evidence notes. It
does **not** contain the source PDFs.

The private PDFs remain in a local literature library outside this repository
and must not be committed to a public repository. A local, gitignored inventory
records their filenames and review status.

- `slrp_literature_evidence.tsv`: selected findings that are directly relevant
  to the gene profiles, Results/Discussion drafts, and biological interpretation.
- `PAPER_EVIDENCE_AND_THESIS_USE.md`: human-readable overview of every paper in
  the thesis bibliography plus the additional local lumican review, including
  the supported claim, current thesis location, source check, and main caveat.
- `paper_evidence_and_thesis_use.tsv`: sortable detailed version of the same
  overview, including the current citation context.
- `scripts/build_paper_use_map.py`: rebuilds both overview files after chapter
  text or citation locations change.
- `../overview/SIX_GENE_PROFILES.md`: integration of these findings with the
  computational results.
- `../../thesis/bibliography.bib`: thesis-ready citation records.

The supplied Perplexity material was used only as a literature map. Statements
were checked against the local PDFs or linked open-access primary articles. No
verbatim abstract text is stored here.
