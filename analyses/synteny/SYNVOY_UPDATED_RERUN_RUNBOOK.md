# Updated SynVoy rerun runbook

Last updated: 2026-09-07

## Scope

This batch refreshes BGN, FMOD, PRELP, and EPYC with the annotation-aware
SynVoy code. DCN and BGN have completed and are integrated. The FMOD refresh at
`results/fmod_human_15species_dev_20260906` was interrupted by `SIGHUP` during
iterative-search wave 6 of 7 and has no final report. Preserve that result
directory and its Nextflow work for a later explicit resume; do not launch a
duplicate FMOD run. PRELP and EPYC have not been launched. OGN and LUM are
deliberately excluded because their completed runs are newer.

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

## 1. Open Ubuntu/WSL and run the dry check

```bash
cd "$PROJECT_ROOT"
bash analyses/synteny/scripts/run_updated_synvoy_panel.sh \
  --dry-run \
  --run-tag dev_20260906 \
  PRELP EPYC
```

The dry check must report that all 14 target genome/GFF pairs pass and print
three `./run_synvoy.sh` commands. It does not run Nextflow. If it reports an
existing output directory, inspect that directory instead of deleting it and
choose a new dated tag only when it is not a resumable run.

## 2. Start the next overnight run

DCN completed after about 16 hours and BGN after about 18 hours. No saved result
is a symlink into `work/`. Clear completed work that is no longer needed for
resume and check disk space before continuing. The remaining order is FMOD,
then PRELP and EPYC. Each invocation appends its provenance and refuses to
overwrite an existing gene result.

First create the log directory because shell redirection happens before the
launcher itself starts:

```bash
mkdir -p "$SYNVOY_ROOT/results/rerun_launcher_logs"
```

The former FMOD launch command must not be run again because its dated output
directory now exists. Resume that exact run only when the SynVoy work is
explicitly restarted. After FMOD has a checked final report, use the launcher
for PRELP and then EPYC.

Example for a new PRELP run:

```bash
cd "$PROJECT_ROOT"
nohup bash analyses/synteny/scripts/run_updated_synvoy_panel.sh \
  --run-tag dev_20260906 \
  PRELP \
  > "$SYNVOY_ROOT/results/rerun_launcher_logs/prelp_dev_20260906.console.log" 2>&1 &
echo $!
```

Record the printed process ID. Closing the terminal should not stop a `nohup`
run. Run only one gene at a time while Windows disk space remains tight.

## 3. Monitor without repeatedly restarting it

```bash
tail -f "$SYNVOY_ROOT/results/rerun_launcher_logs/prelp_dev_20260906.console.log"
```

Press `Ctrl+C` to stop following the log; that does not stop SynVoy. Check the
process separately with:

```bash
ps -fp YOUR_PROCESS_ID
```

Per-gene logs are written under:

```text
$SYNVOY_ROOT/results/rerun_launcher_logs/
```

## 4. Check completion

```bash
for gene in fmod prelp epyc; do
  report="$SYNVOY_ROOT/results/${gene}_human_15species_dev_20260906/synvoy_report.json"
  if [[ -s "$report" ]]; then
    echo "COMPLETE  $gene  $report"
  else
    echo "MISSING   $gene  $report"
  fi
done
```

The expected final directories are:

```text
results/bgn_human_15species_dev_20260905/
results/fmod_human_15species_dev_20260906/
results/prelp_human_15species_dev_20260906/
results/epyc_human_15species_dev_20260906/
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

## 6. After all three remaining refresh reports exist

Return to the project and parse the dated reports, compare them with the
current 98-row table, regenerate P1/P2/P3 review queues, and then update the
final evidence tables and thesis figures. Automated HIGH/MEDIUM calls remain
candidate-level synteny evidence, not final one-to-one orthology assignments.
