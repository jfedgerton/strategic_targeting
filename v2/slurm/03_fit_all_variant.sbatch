#!/bin/bash
#SBATCH --job-name=v2_fit_all
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --time=120:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-15
#SBATCH --output=%x_%A_%a.out
# Full choice set (variant=all, ~1.7M alternative rows per lag) for the a/b models only: expected > 48 h, so paid account.
# Superseded for the main results by 08_main_grid.sbatch (all eight models x four lags); kept for ad hoc reruns.
# Cells: 4 models x 4 lags for the main sample (mid, PageRank). Run 04_pilot.sbatch FIRST to get the
# per-iteration time on Roar hardware and adjust --time (or switch to --sampled_m=50, which cuts the
# alternative rows ~5x and is consistent under uniform sampling; see README).
# Slurm runs a spooled copy of this script, so locate the repo via the submit dir (submit from the repo root).
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
LAGS=(1 3 5 10); MODELS=(M1a M1b M2a M2b)
i=$SLURM_ARRAY_TASK_ID
l=${LAGS[$(( i % 4 ))]}; i=$(( i / 4 )); m=${MODELS[$(( i % 4 ))]}
echo "cell: sample=mid lag=$l model=$m centrality=pr variant=all"
Rscript v2/R/02_fit_model.R --sample=mid --lag=$l --variant=all --model=$m --centrality=pr \
  --K=2 --chains=4 --iter=2000 --warmup=1000 --adapt_delta=0.9 --max_treedepth=10
