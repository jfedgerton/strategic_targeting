#!/bin/bash
# Submit 09_joiner_m2c.sbatch (PI decisions 2026-11-04; --entry=later, --roles=$ROLES) as two arrays, from the repo root. Nothing is submitted unless
# this script is run. ROLES defaults to bcj, the MAIN coding (PI decision 2026-11-04, especially for M2c); ROLES=bc is the optional appendix coding.
#   (a) tasks 0-7    M1c all, M1cr all, M1c joiner_set=opposing, M1b (fmid / mid), long settings      -> job A
#                    (SKIP_M1B=TRUE: tasks 0-5 only; the M1b cells carry no roles, so submit them with the FIRST coding only)
#       + tasks 16-17 M1c --variant=polrel (fmid / mid), long settings
#   (b) tasks 8-15   M2c K1 / K2 x joiner_set all / opposing (fmid / mid; inits from the M1c fits of (a)) -> --dependency=afterok:A --kill-on-invalid-dep=yes
#       + tasks 18-21 M2c K1 / K2 --variant=polrel (inits from tasks 16-17)
# EVENTS=pairs (default: every cell --entry=later, v2/data_entry, suffix _ent; unchanged) | dymid (every cell --events=dymid, v2/data_dymid, suffix _dy;
# the interacting-dyads foundation, PI-approved 2026-11-04): same cells, same ROLES default (bcj).
# Environment passed through: EVENTS, ROLES, TAG, LONG_ITER, LONG_WARMUP, LONG_ADAPT_DELTA, CHECK_LONG, M2_ITER, M2_WARMUP, M2_ADAPT_DELTA, THREADS_PER_CHAIN.
# Extra arguments go to both sbatch calls (e.g. --time=...). A_JOB=<id> reuses an (a) array submitted earlier instead of submitting a new one.
# Example: bash v2/slurm/submit_joiner_m2c.sh (bcj, main); optionally later ROLES=bc SKIP_M1B=TRUE bash v2/slurm/submit_joiner_m2c.sh (appendix)
# Dry check first (open account): sbatch --account=open --export=ALL,ROLES=bcj,DRY_RUN=TRUE --time=00:30:00 --mem=16G --cpus-per-task=1 v2/slurm/09_joiner_m2c.sbatch
set -euo pipefail
[ -f v2/slurm/09_joiner_m2c.sbatch ] || { echo "run from the repo root" >&2; exit 1; }
ROLES=${ROLES:-bcj}      # bcj = the main coding (PI decision 2026-11-04); ROLES=bc = optional appendix coding
case $ROLES in bc|bcj) export ROLES ;; *) echo "ROLES must be bcj (default) or bc" >&2; exit 1 ;; esac
EVENTS=${EVENTS:-pairs}  # pairs = the entry build (default, unchanged); dymid = the interacting-dyads build
case $EVENTS in pairs|dymid) export EVENTS ;; *) echo "EVENTS must be pairs (default) or dymid" >&2; exit 1 ;; esac
for a in "$@"; do case $a in --array*|--dependency*) echo "ERROR: $a is set by this script" >&2; exit 1 ;; esac; done
if [ -z "${A_JOB:-}" ]; then
  IDX=0-7,16-17; [ "${SKIP_M1B:-FALSE}" = "TRUE" ] && IDX=0-5,16-17
  A_JOB=$(sbatch --parsable --array=$IDX --export=ALL "$@" v2/slurm/09_joiner_m2c.sbatch)
  echo "(a) M1 array (roles $ROLES, events $EVENTS): job $A_JOB (tasks $IDX)"
else
  echo "(a) reusing job $A_JOB"
fi
B_JOB=$(sbatch --parsable --array=8-15,18-21 --export=ALL --dependency=afterok:$A_JOB --kill-on-invalid-dep=yes "$@" v2/slurm/09_joiner_m2c.sbatch)
echo "(b) M2c array (roles $ROLES, events $EVENTS): job $B_JOB (tasks 8-15,18-21; afterok:$A_JOB)"
