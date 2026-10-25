#!/bin/bash
#SBATCH --job-name=v2_measure
#SBATCH --account=open
#SBATCH --partition=basic   # the open partition no longer exists; --account=open bills the free READ credits
#SBATCH --requeue
#SBATCH --time=47:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-5
#SBATCH --output=%x_%A_%a.out
#SBATCH --open-mode=append   # a requeued task appends to its log instead of truncating it (2026-10-25 C7)
# Refit-based measurement tests on the MAIN cell (sample mid, lag 5, variant all = the full all-states choice set, PageRank;
# 2026-10-25 C6: formerly polrel). Needs the production a / b fits of that cell from 08_main_grid.sbatch (fit_mid_lag5_all_pr_M1a.rds,
# ..._M1b.rds, ..._M2a_K2.rds, ..._M2b_K2.rds); run_measurement_real.R refuses pilots, stale fits and fits that do not score the
# assembled cell. Task layout (one task per family for the permutation and inflation tests, one task per cutoff and family for the
# temporal holdout, whose constants are frozen on the training years):
#   0 = m1 permutation + inflation      1 = m2 permutation + inflation
#   2 = m1 holdout cutoff 1990          3 = m2 holdout cutoff 1990
#   4 = m1 holdout cutoff 2000          5 = m2 holdout cutoff 2000
# perm_mode=rescore is the fast RELIANCE DIAGNOSTIC (not a calibrated p-value); perm_mode=refit gives the calibrated permutation
# p-value at the cost of 2 x n_perm refits (override with --export=ALL,PERM_MODE=refit,N_PERM=20).
# Slurm runs a spooled copy of this script, so locate the repo via the submit dir (submit from the repo root).
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
FAMS=(m1 m2); i=${SLURM_ARRAY_TASK_ID:?run as an array}; fam=${FAMS[$(( i % 2 ))]}
PM=${PERM_MODE:-rescore}; NP=${N_PERM:-50}
COMMON="--sample=mid --lag=5 --variant=all --centrality=pr --family=$fam --K=2 --r0=0.5 --chains=4 --iter=1500 --warmup=750"
case $i in
  0|1) echo "task $i: family $fam permutation ($PM, n_perm $NP) + inflation"
       Rscript v2/R/run_measurement_real.R $COMMON --tests=permutation,inflation --n_perm=$NP --perm_mode=$PM --deltas=0.25,0.5,1.0 ;;
  2|3) echo "task $i: family $fam temporal holdout, cutoff 1990"
       Rscript v2/R/run_measurement_real.R $COMMON --tests=holdout --cutoff=1990 ;;
  4|5) echo "task $i: family $fam temporal holdout, cutoff 2000"
       Rscript v2/R/run_measurement_real.R $COMMON --tests=holdout --cutoff=2000 ;;
  *) echo "array index $i out of range 0-5" >&2; exit 1 ;;
esac
