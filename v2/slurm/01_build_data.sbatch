#!/bin/bash
#SBATCH --job-name=v2_build_data
#SBATCH --account=open
#SBATCH --partition=open
#SBATCH --requeue
#SBATCH --time=02:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=%x_%j.out
# Builds every v2 data file (v2/data/*.rds + manifest.csv). ~3 min locally; 2 h is generous.
source "$(dirname "$0")/00_env.sh"
cd "$V2_ROOT/.." || exit 1
Rscript v2/R/compile_models.R
Rscript v2/R/01_build_data.R
