#!/bin/bash
#SBATCH --job-name=v2_build_data
#SBATCH --account=open
#SBATCH --partition=basic   # the open partition no longer exists; --account=open bills the free READ credits
#SBATCH --requeue
#SBATCH --time=02:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=%x_%j.out
# Builds every v2 data file (v2/data/*.rds + manifest.csv). ~3 min locally; 2 h is generous.
# Slurm runs a spooled copy of this script, so locate the repo via the submit dir (submit from the repo root).
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
set -eo pipefail   # audit 2026-10-23: never build data if model compilation failed (old binary + new data would run silently)
Rscript v2/R/compile_models.R
Rscript v2/R/01_build_data.R
