#!/bin/bash
#SBATCH --job-name=v2_measure
#SBATCH --account=open
#SBATCH --partition=open
#SBATCH --requeue
#SBATCH --time=47:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-3
#SBATCH --output=%x_%A_%a.out
# Refit-based measurement tests (temporal holdout, permutation, measurement-error inflation)
# on the main cell, for M1 and M2 and for two holdout cutoffs. Needs the a/b fits of the cell.
# Slurm runs a spooled copy of this script, so locate the repo via the submit dir (submit from the repo root).
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
FAMS=(m1 m2); CUTS=(1990 2000)
i=$SLURM_ARRAY_TASK_ID; fam=${FAMS[$(( i % 2 ))]}; cut=${CUTS[$(( i / 2 ))]}
Rscript v2/R/run_measurement_real.R --sample=mid --lag=5 --variant=polrel --family=$fam --K=2 \
  --cutoff=$cut --n_perm=50 --perm_mode=rescore --deltas=0.25,0.5,1.0 --r0=0.5 --chains=4 --iter=1500 --warmup=750
