#!/bin/bash
# Submit 09_joiner_m2c.sbatch (PI-approved 2026-11-04) as three arrays, from the repo root. Nothing is submitted unless this script is run.
#   (a) tasks 0-1    M1c --joiner_set=opposing, fmid / mid, long settings          -> job A
#   (b) tasks 2-5    M2c K1 / K2, joiner_set all (init: the EXISTING long M1c fits) -> no dependency
#   (c) tasks 6-9    M2c K1 / K2, joiner_set opposing (init: the fits of (a))      -> --dependency=afterok:A --kill-on-invalid-dep=yes
# Environment passed through: TAG, LONG_ITER, LONG_WARMUP, LONG_ADAPT_DELTA, CHECK_LONG, M2_ITER, M2_WARMUP, M2_ADAPT_DELTA, THREADS_PER_CHAIN.
# Extra arguments go to every sbatch call (e.g. --time=...). A_JOB=<id> reuses an (a) array submitted earlier instead of submitting a new one.
# Dry check first (open account): sbatch --account=open --export=ALL,DRY_RUN=TRUE --time=00:30:00 --mem=16G --cpus-per-task=1 v2/slurm/09_joiner_m2c.sbatch
set -euo pipefail
[ -f v2/slurm/09_joiner_m2c.sbatch ] || { echo "run from the repo root" >&2; exit 1; }
for a in "$@"; do case $a in --array*|--dependency*) echo "ERROR: $a is set by this script" >&2; exit 1 ;; esac; done
if [ -z "${A_JOB:-}" ]; then
  A_JOB=$(sbatch --parsable --array=0-1 "$@" v2/slurm/09_joiner_m2c.sbatch)
  echo "(a) M1c joiner_set=opposing array: job $A_JOB (tasks 0-1)"
else
  echo "(a) reusing job $A_JOB"
fi
B_JOB=$(sbatch --parsable --array=2-5 "$@" v2/slurm/09_joiner_m2c.sbatch)
echo "(b) M2c joiner_set=all array: job $B_JOB (tasks 2-5; inits from the existing long M1c fits)"
C_JOB=$(sbatch --parsable --array=6-9 --dependency=afterok:$A_JOB --kill-on-invalid-dep=yes "$@" v2/slurm/09_joiner_m2c.sbatch)
echo "(c) M2c joiner_set=opposing array: job $C_JOB (tasks 6-9; afterok:$A_JOB)"
