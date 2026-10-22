#!/bin/bash
#SBATCH --job-name=v2_pilot
#SBATCH --account=open
#SBATCH --partition=basic   # the open partition no longer exists; --account=open bills the free READ credits
#SBATCH --requeue
#SBATCH --time=12:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem=16G           # basic bills 4 GB per core, so 16 GB on 1 core is billed as 4 cores
#SBATCH --array=1-4
#SBATCH --output=%x_%A_%a.out
# Timing pilots: 1 chain x 60 iterations each, run in parallel as an array.
#   1 = all states, M1b    2 = all states, M2b    3 = polrel, M2b    4 = all states, M2b, 50 sampled alternatives
# Read "sampling finished in X min" from each log; per-iteration time x 2000 gives the per-chain wall time.
# Slurm runs a spooled copy of this script, so locate the repo via the submit dir (submit from the repo root).
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
case "${SLURM_ARRAY_TASK_ID:?run as an array}" in
  1) Rscript v2/R/02_fit_model.R --sample=mid --lag=5 --variant=all    --model=M1b --pilot=TRUE ;;
  2) Rscript v2/R/02_fit_model.R --sample=mid --lag=5 --variant=all    --model=M2b --pilot=TRUE ;;
  3) Rscript v2/R/02_fit_model.R --sample=mid --lag=5 --variant=polrel --model=M2b --pilot=TRUE ;;
  4) Rscript v2/R/02_fit_model.R --sample=mid --lag=5 --variant=all    --model=M2b --pilot=TRUE --sampled_m=50 ;;
  *) echo "unknown array index ${SLURM_ARRAY_TASK_ID}" >&2; exit 1 ;;
esac
