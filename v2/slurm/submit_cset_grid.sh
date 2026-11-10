#!/bin/bash
# Submit 11_cset_grid.sbatch (the choice-set grid, PI decision 2026-11-08; README section 13) as two arrays, from the repo root. Nothing is submitted unless
# this script is run.
#   (a) the M1c tasks (default M1_IDX=2-3: dymid all / polrel) on the open account, as the sbatch header says                                -> job A
#   (b) the M2c K = 2 tasks (default M2_IDX=6-7, inits from the M1c fits of (a)) on jfe4_cr_default / basic with the run-7 resources,
#       --dependency=afterok:A --kill-on-invalid-dep=yes (one dependency on the whole M1c array, as in submit_joiner_m2c.sh)
# The pairs all / polrel cells (tasks 0-1 and 4-5) exist from run 7 and are skipped by default. To rerun all eight: M1_IDX=0-3 M2_IDX=4-7.
# The MID-history choice set (mid3) left the grid on 2026-11-10 (audit P1-1; --variant=mid3 is refused).
# Environment passed through: TAG, LONG_ITER, LONG_WARMUP, LONG_ADAPT_DELTA, M2_ITER, M2_WARMUP, M2_ADAPT_DELTA, THREADS_PER_CHAIN.
# M2_OPTS overrides the M2c array's account / partition / time (default "--account=jfe4_cr_default --partition=basic --time=7-00:00:00").
# Extra arguments go to both sbatch calls. A_JOB=<id> reuses an (a) array submitted earlier instead of submitting a new one.
# Check the allocation before submitting (b): get_balance. Dry check first (open account):
#   sbatch --export=ALL,DRY_RUN=TRUE --array=0-7 --time=00:30:00 --mem=16G --cpus-per-task=1 v2/slurm/11_cset_grid.sbatch
set -euo pipefail
[ -f v2/slurm/11_cset_grid.sbatch ] || { echo "run from the repo root" >&2; exit 1; }
for a in "$@"; do case $a in --array*|--dependency*) echo "ERROR: $a is set by this script" >&2; exit 1 ;; esac; done
M1_IDX=${M1_IDX:-2-3}; M2_IDX=${M2_IDX:-6-7}
M2_OPTS=${M2_OPTS:---account=jfe4_cr_default --partition=basic --time=7-00:00:00}
if [ -z "${A_JOB:-}" ]; then
  A_JOB=$(sbatch --parsable --array=$M1_IDX --export=ALL "$@" v2/slurm/11_cset_grid.sbatch); A_JOB=${A_JOB%%;*}
  echo "(a) M1c array (open account): job $A_JOB (tasks $M1_IDX)"
else
  echo "(a) reusing job $A_JOB"
fi
B_JOB=$(sbatch --parsable --array=$M2_IDX --export=ALL $M2_OPTS --dependency=afterok:$A_JOB --kill-on-invalid-dep=yes "$@" v2/slurm/11_cset_grid.sbatch); B_JOB=${B_JOB%%;*}
echo "(b) M2c K2 array ($M2_OPTS): job $B_JOB (tasks $M2_IDX; afterok:$A_JOB)"
