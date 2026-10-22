#!/bin/bash
#SBATCH --job-name=v2_fit
#SBATCH --account=open
#SBATCH --partition=open
#SBATCH --requeue
#SBATCH --time=47:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-95
#SBATCH --output=%x_%A_%a.out
# Job array over sample x lag x model (x variant) for the OPEN partition (< 48 h).
# Default cells here use the politically relevant choice set (variant=polrel, ~220k alternative
# rows) which fits comfortably in 48 h at 4 chains x 2000 iterations. The "all" variant
# (~1.7M rows) goes through 03_fit_all_variant.sbatch on jfe4_cr_default/basic.
#
# Index layout: 3 samples x 4 lags x 4 models x 2 centralities = 96 cells
# Slurm runs a spooled copy of this script, so locate the repo via the submit dir (submit from the repo root).
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
SAMPLES=(mid fmid war); LAGS=(1 3 5 10); MODELS=(M1a M1b M2a M2b); CENTS=(pr eig)
i=$SLURM_ARRAY_TASK_ID
s=${SAMPLES[$(( i % 3 ))]}; i=$(( i / 3 ))
l=${LAGS[$(( i % 4 ))]};    i=$(( i / 4 ))
m=${MODELS[$(( i % 4 ))]};  i=$(( i / 4 ))
c=${CENTS[$(( i % 2 ))]}
echo "cell: sample=$s lag=$l model=$m centrality=$c variant=polrel"
Rscript v2/R/02_fit_model.R --sample=$s --lag=$l --variant=polrel --model=$m --centrality=$c \
  --K=2 --chains=4 --iter=2000 --warmup=1000 --adapt_delta=0.9 --max_treedepth=10
