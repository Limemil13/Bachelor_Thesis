# Thesis overview: start here

This directory is the control centre for the bachelor-thesis project. The
analysis-specific inputs and outputs remain in the neighbouring analysis
directories; this folder contains the cross-analysis interpretation and final
decision material.

## Read these first

1. `AUTHORITATIVE_RESULT_LOCATIONS.md` — the exact table and figure to use for
   each result.
2. `RESULTS_INTERPRETATION_GUIDE.md` — how to read the results without
   overstating them.
3. `GENE_FUNCTION_REFERENCE.md` and `SIX_GENE_PROFILES.md` — cited gene biology
   and its connection to this thesis.
4. `METHODS_ALGORITHM_NOTES.md` — exact per-method inputs, settings, scripts,
   outputs, algorithms, and interpretation limits for the Methods chapter.
5. `WORKFLOW_OVERVIEW.md` — the analysis sequence and links between evidence
   layers.

## Directory contents

- `tables/`: final gene-level, candidate-level, and result-synthesis tables.
- `figures/`: thesis-wide overview and evidence-summary figures.
- `scripts/`: builders for the final tables and figures.
- The few filenames containing `six_gene` are historical names;
their current contents use the seven-gene comparative panel (BGN, DCN, EPYC,
FMOD, LUM, OGN, and PRELP). They are retained for now because analysis builders
and thesis links already depend on those paths.

## Editing rule

Do not hand-edit generated evidence tables merely to improve a conclusion.
Update the appropriate source table, run its builder, and then run
`verify_thesis_analysis_outputs.py`.
