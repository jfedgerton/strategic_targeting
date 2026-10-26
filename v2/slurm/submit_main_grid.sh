#!/bin/bash
# Submit the main grid as two dependent arrays of 02_main_grid.sbatch (PI decision 2026-10-25): the M1 cells first, then the
# M2 cells with --dependency=afterok on the M1 array, so every M2 task finds the M1 fit of its cell in v2/output and starts
# its chains at that fit's posterior means (02_fit_model.R --init_from default; a missing M1 fit file falls back to random inits,
# an unusable one stops the task).
# The array indices are those of 02_main_grid.sbatch (0-95: lag slowest, then model, then sample): in each lag block of 24,
# indices 0-11 are the M1 models and 12-23 the M2 models. One afterok on the whole M1 array is used instead of per-cell
# dependencies (simpler; the M2 array waits until every M1 task has finished successfully). Audit fix 2 (2026-10-26): the M2
# array is submitted with --kill-on-invalid-dep=yes, so one failed M1 task cancels the M2 array instead of leaving it pending
# forever (resubmit the failed M1 task and the M2 array by hand: M1_IDX=<its index> M2_IDX=<its M2 indices> bash this script).
#
# Usage (from the repo root):
#   bash v2/slurm/submit_main_grid.sh [extra sbatch options for BOTH arrays, e.g. --export=ALL,DRY_RUN=TRUE --account=open]
# Environment overrides:
#   M1_IDX / M2_IDX   the array index lists (defaults = the full grid). One lag alone, e.g. lag 3: M1_IDX=24-35 M2_IDX=36-47
#                     (lag 1: 0-11 / 12-23; lag 5: 48-59 / 60-71; lag 10: 72-83 / 84-95).
#   M1_OPTS / M2_OPTS sbatch options for one array only, e.g. M2_OPTS="--time=10-00:00:00" (the M2 cells need about twice the
#                     M1 walltime). Word-split by the shell; no spaces inside one option.
#   M1_JOB            the job id of an M1 array submitted EARLIER: no M1 array is submitted, the M2 array depends on that job
#                     (afterok:$M1_JOB, kill-on-invalid-dep=yes). For a second M2 array on the same M1 fits (the K = 1 array below).
#   M2_SWITCHES       spec.R switches for the M2 array only, exported to its tasks as SWITCHES (02_main_grid.sbatch appends them to
#                     the 02_fit_model.R call), e.g. M2_SWITCHES="--ame_K=1". The M1 array never sees them (an M1 cell refuses a K).
#                     Needs the default environment export of sbatch (no --export, or --export=ALL,...); --export=NONE drops it.
# --array in "$@", M1_OPTS or M2_OPTS is refused (it would override both index lists; use M1_IDX / M2_IDX).
# Examples:
#   bash v2/slurm/submit_main_grid.sh                                                   # the 96-cell grid
#   M1_IDX=24-35 M2_IDX=36-47 bash v2/slurm/submit_main_grid.sh                         # lag 3 only
#   M2_OPTS="--time=10-00:00:00" bash v2/slurm/submit_main_grid.sh                      # longer walltime for the M2 array
#   bash v2/slurm/submit_main_grid.sh --account=open --export=ALL,DRY_RUN=TRUE --time=00:30:00 --mem=16G --cpus-per-task=1   # dry run (open)
# AME K = 1 beside the main K = 2 at lag 3 (PI request 2026-10-26; README "AME K = 1 at lag 3"): the lag-3 M1 array once, then TWO M2 arrays
# on it, K = 2 (the main grid) and K = 1 (cells <cell>_K1, which never collide with the K = 2 names and start from the same M1 fits):
#   M1_IDX=24-35 M2_IDX=36-47 bash v2/slurm/submit_main_grid.sh                         # (a) + (b): M1 array <id>, M2 K = 2 array afterok:<id>
#   M1_JOB=<id> M2_IDX=36-47 M2_SWITCHES="--ame_K=1" bash v2/slurm/submit_main_grid.sh  # (c): M2 K = 1 array afterok:<id>, no new M1 array
set -eo pipefail
[ -f v2/slurm/02_main_grid.sbatch ] || { echo "run from the repo root" >&2; exit 1; }
M1_IDX="${M1_IDX:-0-11,24-35,48-59,72-83}"
M2_IDX="${M2_IDX:-12-23,36-47,60-71,84-95}"
for a in "$@" ${M1_OPTS:-} ${M2_OPTS:-}; do
  case "$a" in --array|--array=*|-a|-a*) echo "ERROR: do not pass --array; set M1_IDX / M2_IDX instead (got '$a')" >&2; exit 1 ;; esac
done
if [ -n "${M2_SWITCHES:-}" ]; then
  case "$M2_SWITCHES " in *"--ame_K=2 "*) echo "ERROR: M2_SWITCHES='$M2_SWITCHES' is the main K = 2 grid; run it without M2_SWITCHES" >&2; exit 1 ;; esac
fi
if [ -n "${M1_JOB:-}" ]; then
  [[ "$M1_JOB" =~ ^[0-9]+$ ]] || { echo "ERROR: M1_JOB must be a Slurm job id (got '$M1_JOB')" >&2; exit 1; }
  m1=$M1_JOB
  echo "M1 array NOT submitted: the M2 array depends on the existing job $m1 (M1_JOB)"
else
  # sbatch --parsable prints "<jobid>" or "<jobid>;<cluster>"; keep the job id only
  m1_raw=$(sbatch --parsable --array="$M1_IDX" "$@" ${M1_OPTS:-} v2/slurm/02_main_grid.sbatch)
  m1=${m1_raw%%;*}
  [[ "$m1" =~ ^[0-9]+$ ]] || { echo "ERROR: sbatch did not return a job id for the M1 array (got '$m1_raw')" >&2; exit 1; }
  echo "M1 array submitted: job $m1 (indices $M1_IDX)"
fi
# the M2 array inherits SWITCHES through sbatch's environment export (default ALL); M2_SWITCHES is never given to the M1 array
m2_raw=$(SWITCHES="${M2_SWITCHES:-}" sbatch --parsable --array="$M2_IDX" --dependency="afterok:$m1" --kill-on-invalid-dep=yes "$@" ${M2_OPTS:-} v2/slurm/02_main_grid.sbatch)
m2=${m2_raw%%;*}
[[ "$m2" =~ ^[0-9]+$ ]] || { echo "ERROR: sbatch did not return a job id for the M2 array (got '$m2_raw'); the M1 array $m1 is running" >&2; exit 1; }
echo "M2 array submitted: job $m2 (indices $M2_IDX${M2_SWITCHES:+, switches '$M2_SWITCHES'}), afterok:$m1, kill-on-invalid-dep=yes"
