#!/bin/bash
# Submit run 1 (13_opp_weights.sbatch; README section 14) as three chained jobs, from the repo root. Nothing is submitted unless this script is run.
#   (a) the auxiliary model (STEP=aux; open account, one task, AUX_OPTS)                                                          -> job A
#   (b) the M1c cells (STEP=fit; array FIT_IDX, default 1-7,9-15,17-23,25-31 = the 28 cells with lambda 0.5 or 1; the four lambda = 0 cells 0, 8, 16, 24 are
#       bit-identical to the existing _tlong fits and are left out unless FIT_IDX names them, e.g. FIT_IDX=0-31; open account / basic partition / requeue, 48 h, the M1 resources of 12_run0.sbatch, as the
#       sbatch header says), --dependency=afterok:A --kill-on-invalid-dep=yes                                                      -> job B
#   (c) the mixture over the five auxiliary draws at lambda = 1 (STEP=combine; open account, one task, MIX_OPTS), --dependency=afterok:B (the whole array)
# Environment passed through: TAG, LONG_ITER, LONG_WARMUP, LONG_ADAPT_DELTA, THREADS_PER_CHAIN, OPP_OUTCOME (fit / combine; default mid), OPP_OUTCOMES (aux).
# FIT_OPTS adds sbatch options to the fit array only (e.g. "--account=jfe4_cr_default --partition=basic --time=4-00:00:00"; check get_balance first).
# A_JOB=<id> reuses an auxiliary job submitted earlier. A partial FIT_IDX that leaves out a lambda = 1 cell makes (c) fail on the missing fits: set COMBINE=FALSE.
# Extra arguments go to all three sbatch calls. Dry check of the fit cells once the weights exist: see the header of 13_opp_weights.sbatch.
set -euo pipefail
[ -f v2/slurm/13_opp_weights.sbatch ] || { echo "run from the repo root" >&2; exit 1; }
for a in "$@"; do case $a in --array*|--dependency*) echo "ERROR: $a is set by this script" >&2; exit 1 ;; esac; done
FIT_IDX=${FIT_IDX:-1-7,9-15,17-23,25-31}
AUX_OPTS=${AUX_OPTS:---time=06:00:00 --mem=32G --cpus-per-task=1}
MIX_OPTS=${MIX_OPTS:---time=06:00:00 --mem=64G --cpus-per-task=1}
if [ -z "${A_JOB:-}" ]; then
  A_JOB=$(sbatch --parsable --job-name=v2_oppw_aux --array=0 --export=ALL,STEP=aux $AUX_OPTS "$@" v2/slurm/13_opp_weights.sbatch); A_JOB=${A_JOB%%;*}
  echo "(a) auxiliary model (open account): job $A_JOB"
else
  echo "(a) reusing job $A_JOB"
fi
B_JOB=$(sbatch --parsable --array=$FIT_IDX --export=ALL,STEP=fit ${FIT_OPTS:-} --dependency=afterok:$A_JOB --kill-on-invalid-dep=yes "$@" v2/slurm/13_opp_weights.sbatch); B_JOB=${B_JOB%%;*}
echo "(b) M1c cells: job $B_JOB (tasks $FIT_IDX; afterok:$A_JOB)"
if [ "${COMBINE:-TRUE}" = "TRUE" ]; then
  C_JOB=$(sbatch --parsable --job-name=v2_oppw_mix --array=0 --export=ALL,STEP=combine $MIX_OPTS --dependency=afterok:$B_JOB --kill-on-invalid-dep=yes "$@" v2/slurm/13_opp_weights.sbatch); C_JOB=${C_JOB%%;*}
  echo "(c) mixture at lambda = 1: job $C_JOB (afterok:$B_JOB)"
fi
