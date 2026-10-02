# Integrated result overview

This directory joins the method-specific analyses without collapsing them into
a single numerical score. Conflicts and limitations remain recorded in the
candidate- and gene-level tables.

## Main files

- `AUTHORITATIVE_RESULT_LOCATIONS.md`: final output map.
- `METHODS_ALGORITHM_NOTES.md`: implementation and parameter record.
- `WORKFLOW_OVERVIEW.md`: biological and computational workflow.
- `THESIS_REPRODUCIBILITY_STATUS.md`: current completion status and known
  limitations.
- `tables/final_gene_level_evidence.tsv`: one row per focal gene.
- `tables/final_candidate_level_confidence.tsv`: one row per reviewed candidate
  or diagnostic record.
- `tables/results_interpretation_summary.tsv`: compact cross-analysis result
  summary.
- `verify_thesis_analysis_outputs.py`: consistency checks for the retained
  tables, figures, and scripts.

The filenames containing `six_gene` are retained for compatibility with older
builders. Their current contents use the seven-gene comparative panel: BGN,
DCN, EPYC, FMOD, LUM, OGN, and PRELP.
