#!/bin/bash
#SBATCH --job-name=v2_run0
#SBATCH --account=open
#SBATCH --partition=basic   # the open partition no longer exists; --account=open bills the free READ credits
#SBATCH --requeue
#SBATCH --time=48:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32   # = chains (4) x threads_per_chain (8)
#SBATCH --mem=128G
#SBATCH --array=0-19
#SBATCH --output=%x_%A_%a.out
#SBATCH --open-mode=append
# RUN 0 (built 2026-11-10; README section 13): functional-form and comparison cells, M1 only, on the FULL choice set. Every cell: --comparison=directed,
# lag 3, tag $TAG (default long -> _tlong), the long role-fit sampler settings; the M1c cells --roles=bcj. Events: pairs = --entry=later (v2/data_entry,
# suffix _ent), dymid = --events=dymid (v2/data_dymid, suffix _dy; fmid there = jointly fatal). Index = 4 x cell + 2 x events + sample.
#   idx    cell  switch                 model  events  sample   cell name (pairs fmid shown)
#   0-3    0a    --reach=capdist        M1c    pairs fmid, pairs mid, dymid fmid, dymid mid     fmid_lag3_all_pr_M1c_cdir_reach_bcj_ent_tlong
#   4-7    0b    --subsample=prtarget   M1c    (same order)                                     fmid_lag3_all_pr_M1c_cdir_bcj_sprt_ent_tlong
#   8-11   0c    --subsample=bilateral  M1b    (same order)                                     fmid_lag3_all_pr_M1b_cdir_sbil_ent_tlong
#   12-15  0d    --subsample=majinit    M1c    (same order)                                     fmid_lag3_all_pr_M1c_cdir_bcj_smaj_ent_tlong
#   16-19  0d    --subsample=mininit    M1c    (same order)                                     fmid_lag3_all_pr_M1c_cdir_bcj_smin_ent_tlong
# What each isolates: 0a the distance functional form (squared log distance, distance x initiator CINC level); 0b the events the polrel set drops, with
# the set itself left full (same events, roles and D8 weights as the polrel cell); 0c bilateral events alone (pooled M1b: one role); 0d major-power and minor-power initiators separately (BOTH stages).
# No dependencies (M1 only). Free open account; a preempted task is requeued and starts again from the beginning (log appended).
# DRY RUN: sbatch --export=ALL,DRY_RUN=TRUE --time=00:30:00 --mem=16G --cpus-per-task=1 v2/slurm/12_run0.sbatch
# Needs both builds (01_build_data.R --entry=later and --events=dymid) and the compiled m1_context / m1_latent_additive programs. Nothing is submitted by this file.
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
i=${SLURM_ARRAY_TASK_ID:?run as an array}
[[ "$i" =~ ^[0-9]+$ ]] && [ "$i" -le 19 ] || { echo "array index $i out of range 0-19" >&2; exit 1; }
case $((i / 4)) in
  0) m=M1c; SW="--reach=capdist" ;;
  1) m=M1c; SW="--subsample=prtarget" ;;
  2) m=M1b; SW="--subsample=bilateral" ;;
  3) m=M1c; SW="--subsample=majinit" ;;
  4) m=M1c; SW="--subsample=mininit" ;;
esac
case $(((i % 4) / 2)) in 0) ev=pairs; DATA_ARG="--entry=later" ;; 1) ev=dymid; DATA_ARG="--events=dymid" ;; esac
case $((i % 2)) in 0) s=fmid ;; 1) s=mid ;; esac
DRY=${DRY_RUN:-FALSE}
THR=${THREADS_PER_CHAIN:-8}
TAG=${TAG:-long}
[[ "$TAG" =~ ^[A-Za-z0-9]+$ ]] || { echo "TAG must be letters and digits (got '$TAG')" >&2; exit 1; }
if [ "$m" = "M1b" ]; then ROLE_ARGS=""; else ROLE_ARGS="--roles=bcj"; fi
SAMP="--iter=${LONG_ITER:-6000} --warmup=${LONG_WARMUP:-3000} --adapt_delta=${LONG_ADAPT_DELTA:-0.9}"
echo "task $i -> sample=$s model=$m events=$ev ($DATA_ARG) $SW tag=$TAG; $SAMP; threads_per_chain=$THR cpus=${SLURM_CPUS_PER_TASK:-?} dry_run=$DRY"
Rscript v2/R/02_fit_model.R --sample=$s --lag=3 --model=$m --comparison=directed $DATA_ARG $ROLE_ARGS $SW --tag=$TAG --dry_run=$DRY --threads_per_chain=$THR $SAMP
