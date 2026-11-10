#!/bin/bash
#SBATCH --job-name=v2_cset_grid
#SBATCH --account=open
#SBATCH --partition=basic   # the open partition no longer exists; --account=open bills the free READ credits
#SBATCH --requeue
#SBATCH --time=48:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32   # = chains (4) x threads_per_chain (8)
#SBATCH --mem=128G
#SBATCH --array=2-3
#SBATCH --output=%x_%A_%a.out
#SBATCH --open-mode=append
# THE CHOICE-SET GRID (PI decisions 2026-11-08 and 2026-11-10; README section 13): choice sets {all, polrel} x fatal outcomes {fmid on the pairs events = all
# fatal (--entry=later, v2/data_entry, suffix _ent), fmid on the interacting dyads = jointly fatal (--events=dymid, v2/data_dymid, suffix _dy)} x models
# {M1c, then M2c K = 2 initialised from the matching M1c}; roles bcj, --comparison=directed, lag 3, joiner_set all, tag $TAG (default long -> _tlong).
# The lag-3 MID-history choice set (mid3) was DROPPED from the grid on 2026-11-10 (audit P1-1; --variant=mid3 is refused by spec.R).
#   idx  events  variant  model  needs (init fit in v2/output; M2c only)                         status
#    0   pairs   all      M1c    -                                                              EXISTS (run 7: 09_joiner_m2c.sbatch task 0); not in the default arrays
#    1   pairs   polrel   M1c    -                                                              EXISTS (run 7: task 16); not in the default arrays
#    2   dymid   all      M1c    -                                                              new
#    3   dymid   polrel   M1c    -                                                              new
#    4   pairs   all      M2c    fit_fmid_lag3_all_pr_M1c_cdir_bcj_ent_t<TAG>.rds     (task 0)   EXISTS (run 7: task 10); not in the default arrays
#    5   pairs   polrel   M2c    fit_fmid_lag3_polrel_pr_M1c_cdir_bcj_ent_t<TAG>.rds  (task 1)   EXISTS (run 7: task 20); not in the default arrays
#    6   dymid   all      M2c    fit_fmid_lag3_all_pr_M1c_cdir_bcj_dy_t<TAG>.rds      (task 2)   new
#    7   dymid   polrel   M2c    fit_fmid_lag3_polrel_pr_M1c_cdir_bcj_dy_t<TAG>.rds   (task 3)   new
# The four pairs cells were fitted in run 7 (the full-set fits take many hours to days), so they are NOT rerun by default: the default array of this file
# and of submit_cset_grid.sh is 2-3 (M1c) and 6-7 (M2c). They stay in the table so that a rerun is one --array away (0-1, then 4-5).
# Accounts: the M1c tasks run on the free open account (the header above: --account=open --partition=basic --requeue, the convention of every open-account
# job of this repository, 48 h; a preempted task is requeued and starts again from the beginning, the log is appended). The M2c tasks need the paid
# allocation with the run-7 resources: submit_cset_grid.sh passes --account=jfe4_cr_default --partition=basic --time=7-00:00:00 (same 32 cpus / 128G) on
# the command line, which overrides the header, with --dependency=afterok:<M1c array> --kill-on-invalid-dep=yes. Running an M2c index with a plain
# `sbatch` of this file would put it on open: refused below unless ALLOW_M2C_OPEN=TRUE or the job's account is not open.
# 02_fit_model.R REFUSES to sample an M2c cell whose M1c fit is missing (no random-init fallback) and refuses an M1c fit of another cell.
# Sampler: M1c LONG_ITER / LONG_WARMUP / LONG_ADAPT_DELTA (defaults 6000 / 3000 / 0.9, the long role-fit settings of 09_joiner_m2c.sbatch); M2c the M2
# defaults of spec.R (4000 / 2000 / 0.95) unless M2_ITER / M2_WARMUP / M2_ADAPT_DELTA are set.
# DRY RUN (no sampling, no init fit needed): sbatch --export=ALL,DRY_RUN=TRUE --array=0-7 --time=00:30:00 --mem=16G --cpus-per-task=1 v2/slurm/11_cset_grid.sbatch
# Needs both builds (01_build_data.R --entry=later and --events=dymid) and the compiled m1_context / m2_context programs. Nothing is submitted by this file.
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
i=${SLURM_ARRAY_TASK_ID:?run as an array}
case $i in
   0) ev=pairs; vr=all;    m=M1c ;;
   1) ev=pairs; vr=polrel; m=M1c ;;
   2) ev=dymid; vr=all;    m=M1c ;;
   3) ev=dymid; vr=polrel; m=M1c ;;
   4) ev=pairs; vr=all;    m=M2c ;;
   5) ev=pairs; vr=polrel; m=M2c ;;
   6) ev=dymid; vr=all;    m=M2c ;;
   7) ev=dymid; vr=polrel; m=M2c ;;
  *) echo "array index $i out of range 0-7" >&2; exit 1 ;;
esac
s=fmid; ROLES=bcj; K=2
case $ev in pairs) DATA_ARG="--entry=later"; DS=_ent ;; dymid) DATA_ARG="--events=dymid"; DS=_dy ;; esac
DRY=${DRY_RUN:-FALSE}
THR=${THREADS_PER_CHAIN:-8}
TAG=${TAG:-long}
[[ "$TAG" =~ ^[A-Za-z0-9]+$ ]] || { echo "TAG must be letters and digits (got '$TAG')" >&2; exit 1; }
BASE="--sample=$s --lag=3 --model=$m --comparison=directed $DATA_ARG --variant=$vr --roles=$ROLES --joiner_set=all --tag=$TAG --dry_run=$DRY --threads_per_chain=$THR"
if [ "$m" = "M1c" ]; then
  SAMP="--iter=${LONG_ITER:-6000} --warmup=${LONG_WARMUP:-3000} --adapt_delta=${LONG_ADAPT_DELTA:-0.9}"
else
  SAMP="--ame_K=$K${M2_ITER:+ --iter=$M2_ITER}${M2_WARMUP:+ --warmup=$M2_WARMUP}${M2_ADAPT_DELTA:+ --adapt_delta=$M2_ADAPT_DELTA}"
  if [ "$DRY" != "TRUE" ]; then
    if [ "${SLURM_JOB_ACCOUNT:-open}" = "open" ] && [ "${ALLOW_M2C_OPEN:-FALSE}" != "TRUE" ]; then
      echo "ERROR: task $i is an M2c cell on the open account; submit it through v2/slurm/submit_cset_grid.sh (jfe4_cr_default / basic) or set ALLOW_M2C_OPEN=TRUE" >&2; exit 1
    fi
    want="v2/output/fit_${s}_lag3_${vr}_pr_M1c_cdir_${ROLES}${DS}_t${TAG}.rds"
    [ -f "$want" ] || { echo "ERROR: task $i ($s M2c K$K variant=$vr events=$ev) needs $want (see the table above)" >&2; exit 1; }
  fi
fi
echo "task $i -> sample=$s model=$m roles=$ROLES events=$ev ($DATA_ARG) variant=$vr K=$([ "$m" = M2c ] && echo $K || echo -) tag=$TAG; $SAMP; threads_per_chain=$THR cpus=${SLURM_CPUS_PER_TASK:-?} account=${SLURM_JOB_ACCOUNT:-?} dry_run=$DRY"
Rscript v2/R/02_fit_model.R $BASE $SAMP
