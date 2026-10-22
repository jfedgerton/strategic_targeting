#!/bin/bash
#SBATCH --job-name=v2_bench
#SBATCH --account=open
#SBATCH --partition=open
#SBATCH --requeue
#SBATCH --time=24:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-3
#SBATCH --output=%x_%A_%a.out
# Benchmarks (directed-dyad logit, old gap regression) for each lag on the main sample.
# Slurm runs a spooled copy of this script, so locate the repo via the submit dir (submit from the repo root).
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
LAGS=(1 3 5 10); l=${LAGS[$SLURM_ARRAY_TASK_ID]}
Rscript v2/R/benchmarks.R --sample=mid --lag=$l --variant=polrel --centrality=pr --chains=4 --iter=2000 --warmup=1000
