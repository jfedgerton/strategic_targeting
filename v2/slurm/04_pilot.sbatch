#!/bin/bash
#SBATCH --job-name=v2_pilot
#SBATCH --account=open
#SBATCH --partition=open
#SBATCH --requeue
#SBATCH --time=06:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=%x_%j.out
# Timing pilot: 1 chain x 60 iterations on the main cell with the full choice set and with the
# polrel set. Read "sampling finished in X min" from the logs and scale: per-iteration time x
# (iter 2000) gives the per-chain wall time; chains run in parallel on 4 cores.
# Slurm runs a spooled copy of this script, so locate the repo via the submit dir (submit from the repo root).
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
Rscript v2/R/02_fit_model.R --sample=mid --lag=5 --variant=polrel --model=M2b --pilot=TRUE
Rscript v2/R/02_fit_model.R --sample=mid --lag=5 --variant=all    --model=M2b --pilot=TRUE
Rscript v2/R/02_fit_model.R --sample=mid --lag=5 --variant=all    --model=M2b --pilot=TRUE --sampled_m=50
