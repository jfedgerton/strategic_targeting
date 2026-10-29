#!/bin/bash
# Submit the permutation test of one cell (08_perm_refit.sbatch; 2026-10-29, PI request: the calibrated test on the contiguity + distance
# main cell before the sweep). Three modes, through the environment:
#   MODE=rescore   one job on the open account (reliance diagnostic, no refit):      --account=open --time=24:00:00
#   MODE=refit     the refit ARRAY 1..N_PERM on jfe4_cr_default / basic (one refit of a and b per task), THEN the combine job on the open
#                  account with --dependency=afterok:<array> --kill-on-invalid-dep=yes (a failed task cancels the combine instead of leaving
#                  it pending; re-run the failed index with sbatch --array=<p> and then MODE=combine)
#   MODE=combine   the combine job alone (every part 1..N_PERM must already exist)
# Cell: SAMPLE (fmid) LAG (3) FAMILY (m1) SWITCHES ("--distance=capital"; SWITCHES="" = the main-spec cell) N_PERM (20); CHAINS / ITER / WARMUP
# override the refit sampler settings (default 4 / 1500 / 750, the sequential refit mode's). Extra sbatch options in "$@" go to EVERY job
# (e.g. --time=2-00:00:00 for the array once the per-task time is known: the open jobs take it too, so prefer ARRAY_OPTS / COMBINE_OPTS /
# RESCORE_OPTS for one job only). --array is refused (the index range is 1..N_PERM). THROTTLE=k limits the concurrent array tasks (%k).
# Usage (from the repo root):
#   MODE=rescore bash v2/slurm/submit_perm_refit.sh
#   MODE=refit bash v2/slurm/submit_perm_refit.sh
#   MODE=refit SAMPLE=mid bash v2/slurm/submit_perm_refit.sh
#   MODE=refit ARRAY_OPTS="--time=2-00:00:00" THROTTLE=10 bash v2/slurm/submit_perm_refit.sh
#   MODE=combine bash v2/slurm/submit_perm_refit.sh
set -eo pipefail
[ -f v2/slurm/08_perm_refit.sbatch ] || { echo "run from the repo root" >&2; exit 1; }
M=${MODE:?set MODE=<rescore|refit|combine>}
NP=${N_PERM:-20}
[[ "$NP" =~ ^[0-9]+$ ]] && [ "$NP" -ge 1 ] || { echo "ERROR: N_PERM must be a positive integer (got '$NP')" >&2; exit 1; }
for a in "$@" ${ARRAY_OPTS:-} ${COMBINE_OPTS:-} ${RESCORE_OPTS:-}; do
  case "$a" in --array|--array=*|-a|-a*) echo "ERROR: do not pass --array; the index range is 1..N_PERM (got '$a')" >&2; exit 1 ;; esac
done
jobid() { local raw=$1; local id=${raw%%;*}; [[ "$id" =~ ^[0-9]+$ ]] || { echo "ERROR: sbatch did not return a job id (got '$raw')" >&2; exit 1; }; echo "$id"; }
export SAMPLE=${SAMPLE:-fmid} LAG=${LAG:-3} FAMILY=${FAMILY:-m1} N_PERM=$NP
export SWITCHES=${SWITCHES-"--distance=capital"}
cell_desc="$SAMPLE lag $LAG $FAMILY [${SWITCHES:-main spec}]"
case $M in
  rescore)
    raw=$(MODE=rescore sbatch --parsable --export=ALL --account=open --time=24:00:00 "$@" ${RESCORE_OPTS:-} v2/slurm/08_perm_refit.sbatch)
    j=$(jobid "$raw")
    echo "perm rescore submitted: job $j ($cell_desc, n_perm $NP; open account, reliance diagnostic)" ;;
  refit)
    raw=$(MODE=refit sbatch --parsable --export=ALL --array=1-${NP}${THROTTLE:+%$THROTTLE} "$@" ${ARRAY_OPTS:-} v2/slurm/08_perm_refit.sbatch)
    arr=$(jobid "$raw")
    echo "perm refit array submitted: job $arr (indices 1-$NP, one refit of a and b per task; $cell_desc; jfe4_cr_default / basic)"
    raw=$(MODE=combine sbatch --parsable --export=ALL --account=open --time=04:00:00 --dependency="afterok:$arr" --kill-on-invalid-dep=yes "$@" ${COMBINE_OPTS:-} v2/slurm/08_perm_refit.sbatch)
    c=$(jobid "$raw")
    echo "perm combine submitted: job $c (afterok:$arr, kill-on-invalid-dep=yes; writes measurement_<cell>_permutation_refit.csv + _draws.csv; open account)" ;;
  combine)
    raw=$(MODE=combine sbatch --parsable --export=ALL --account=open --time=04:00:00 "$@" ${COMBINE_OPTS:-} v2/slurm/08_perm_refit.sbatch)
    c=$(jobid "$raw")
    echo "perm combine submitted: job $c ($cell_desc, parts 1..$NP; open account)" ;;
  *) echo "unknown MODE $M (rescore | refit | combine)" >&2; exit 1 ;;
esac
