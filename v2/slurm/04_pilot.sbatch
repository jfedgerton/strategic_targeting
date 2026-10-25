#!/bin/bash
#SBATCH --job-name=v2_pilot
#SBATCH --account=open
#SBATCH --partition=basic   # the open partition no longer exists; --account=open bills the free READ credits
#SBATCH --requeue
#SBATCH --time=12:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8    # = chains (1) x threads_per_chain (8)
#SBATCH --mem=32G                                  # basic bills 4 GB per core: 8 cores x 4 GB
#SBATCH --array=1-4
#SBATCH --output=%x_%A_%a.out
#SBATCH --open-mode=append   # a requeued task appends to its log instead of truncating it (2026-10-25 C7)
# Timing pilots: 1 chain x 60 iterations each, run in parallel as an array, with THREADS_PER_CHAIN
# (default 8) reduce_sum threads; compare with --export=ALL,THREADS_PER_CHAIN=1 --cpus-per-task=1 --mem=4G
# to measure the within-chain speed-up on Roar hardware.
#   1 = all states, M1b    2 = all states, M2b    3 = polrel, M2b    4 = all states, M2b, 50 sampled alternatives
# Read "sampling finished in X min" from each log; per-iteration time x 2000 gives the per-chain wall time.
# Slurm runs a spooled copy of this script, so locate the repo via the submit dir (submit from the repo root).
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
THR=${THREADS_PER_CHAIN:-8}
echo "pilot task ${SLURM_ARRAY_TASK_ID} threads_per_chain=$THR cpus=${SLURM_CPUS_PER_TASK:-?}"
case "${SLURM_ARRAY_TASK_ID:?run as an array}" in
  1) Rscript v2/R/02_fit_model.R --sample=mid --lag=5 --variant=all    --model=M1b --pilot=TRUE --threads_per_chain=$THR ;;
  2) Rscript v2/R/02_fit_model.R --sample=mid --lag=5 --variant=all    --model=M2b --pilot=TRUE --threads_per_chain=$THR ;;
  3) Rscript v2/R/02_fit_model.R --sample=mid --lag=5 --variant=polrel --model=M2b --pilot=TRUE --threads_per_chain=$THR ;;
  4) Rscript v2/R/02_fit_model.R --sample=mid --lag=5 --variant=all    --model=M2b --pilot=TRUE --sampled_m=50 --threads_per_chain=$THR ;;
  *) echo "unknown array index ${SLURM_ARRAY_TASK_ID}" >&2; exit 1 ;;
esac
