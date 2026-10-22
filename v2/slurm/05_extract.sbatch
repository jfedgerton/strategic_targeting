#!/bin/bash
#SBATCH --job-name=v2_extract
#SBATCH --account=open
#SBATCH --partition=open
#SBATCH --requeue
#SBATCH --time=04:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=%x_%j.out
# Typology + AME extraction for every saved fit, then LOO / misclassification for each a-b pair.
# Slurm runs a spooled copy of this script, so locate the repo via the submit dir (submit from the repo root).
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
for f in v2/output/fit_*.rds; do Rscript v2/R/03_extract.R --fit="$f" --r0=0.5; done
for fa in v2/output/fit_*_M1a*.rds v2/output/fit_*_M2a*.rds; do
  fb="${fa/M1a/M1b}"; fb="${fb/M2a/M2b}"
  [ -f "$fb" ] && Rscript v2/R/measurement_tests.R --fit_a="$fa" --fit_b="$fb" --r0=0.5
done
