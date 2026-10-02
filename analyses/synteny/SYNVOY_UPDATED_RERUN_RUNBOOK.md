# Updated SynVoy rerun runbook

Last updated: 2026-10-01

## Scope

This batch refreshes BGN, FMOD, PRELP, and EPYC with the annotation-aware
SynVoy code. DCN, BGN, and FMOD are complete and integrated. The FMOD report
at `results/fmod_human_15species_dev_20260906` is the authoritative source for
all 14 FMOD rows. PRELP was started at
`results/prelp_human_15species_dev_20260914`; its existing iterative-search
checkpoint records four of seven waves, but there is no final report. An
updated EPYC run has not been started. No SynVoy process was active on
2026-09-21. OGN and LUM are not first-priority reruns.

Do not use an old `--run-tag dev_20260906 PRELP` launch command: the actual
PRELP output is dated `dev_20260914` and should be resumed only after checking
its exact Nextflow invocation and cache. No SynVoy run is authorized by this
status note alone.

The launcher runs one new gene at a time, validates all 14 target genome/GFF
pairs before starting, and refuses to overwrite an existing result directory.

Before using the launcher, define the two checkout locations for the current
machine:

```bash
export PROJECT_ROOT=/path/to/Bachelor_Thesis
export SYNVOY_ROOT=/path/to/SynVoy
# Optional when Conda is not already initialized:
export CONDA_SH=/path/to/conda.sh
```

## 1. Check the existing work before any new run

FMOD already has a final report and has been parsed, compared with the
historical FMOD rows, manually reviewed, and integrated. PRELP has an existing
dated output directory and a partial Nextflow cache. Confirm its latest
checkpoint, input FASTA/GFF pairing, available disk space, and original
launcher invocation before using Nextflow `-resume`. Do not start PRELP with
a new tag simply because its final report is absent.

EPYC is the only fresh updated-tool full-gene run still to start, if that
sensitivity analysis remains in scope. Use a new date-specific tag and the
launcher dry-run first; the dry-run should validate all 14 target genome/GFF
pairs. Run only when explicitly requested, one gene at a time. Keep logs and
dated result directories; do not overwrite historical reports or resumable
work.

The separate opossum matched-input issue affects all seven gene reviews.
Resolve and rerun those comparisons before presenting opossum synteny as
confirmed. This is distinct from an updated full-gene PRELP or EPYC run.

## 4. Check completion

```bash
test -s "$SYNVOY_ROOT/results/fmod_human_15species_dev_20260906/synvoy_report.json"
test -s "$SYNVOY_ROOT/results/prelp_human_15species_dev_20260914/synvoy_report.json"
```

The first check currently succeeds; the second currently fails because the
PRELP run is incomplete. An updated EPYC output directory does not yet exist.
For any later EPYC launch, record its actual dated tag before checking a path.

The expected final directories are:

```text
results/bgn_human_15species_dev_20260905/
results/fmod_human_15species_dev_20260906/
results/prelp_human_15species_dev_20260914/  # partial; no report yet
results/epyc_human_15species_<future_tag>/  # not launched
results/dcn_human_15species_dev_20260903/
```

Do not overwrite or delete historical result directories. A non-empty
`synvoy_report.json` is the minimum completion indicator; the report still
needs to be parsed and manually reviewed before an orthology decision.

## 5. If the batch fails

Copy the final 80-120 lines of both the panel console log and the current
per-gene log. Also report the gene that was running and whether its output
directory exists. Do not immediately delete `.nextflow`, `work`, or the dated
result directory: those may be needed for a safe resume. A new tag is for a
genuinely clean rerun, not a substitute for diagnosing a resumable failure.

## 6. After a refresh report exists

The dated FMOD report has completed this sequence. Apply it to PRELP and EPYC
only if their optional updated reports are completed. Compare with the
current 98-row table, regenerate P1/P2/P3 review queues, and then update the
final evidence tables and thesis figures. Automated HIGH/MEDIUM calls remain
candidate-level synteny evidence, not final one-to-one orthology assignments.
