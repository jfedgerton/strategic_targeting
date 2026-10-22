#!/bin/bash
#SBATCH --job-name=v2_main_extract
#SBATCH --account=open
#SBATCH --partition=basic
#SBATCH --requeue
#SBATCH --time=08:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --array=0-11
#SBATCH --output=%x_%A_%a.out
# After 08_main_grid.sbatch: per (sample, lag) cell group, run extraction + typology (incl. the appendix
# 9-cell classification, reach and focus outputs) for each of the eight fits, then the cross-model
# LOO comparison (material vs relational vs independent vs joint, and M1 vs M2) and the a-vs-b
# misclassification. 12 tasks = 3 samples x 4 lags.
# Submit from the repo root: sbatch v2/slurm/09_main_extract.sbatch
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
i=${SLURM_ARRAY_TASK_ID:?run as an array}
SAMPLES=(mid fmid war); LAGS=(1 3 5 10)
s=${SAMPLES[$(( i % 3 ))]}; l=${LAGS[$(( i / 3 ))]}
echo "task $i -> sample=$s lag=$l"
for m in M1m M1r M1a M1b M2m M2r M2a M2b; do
  f="v2/output/fit_${s}_lag${l}_all_pr_${m}.rds"; [ -f "$f" ] || f="v2/output/fit_${s}_lag${l}_all_pr_${m}_K2.rds"
  if [ -f "$f" ]; then Rscript v2/R/03_extract.R --fit="$f" --r0=0.5; else echo "missing fit for $m: $f" >&2; fi
done
Rscript v2/R/compare_models.R --sample=$s --lag=$l --variant=all --centrality=pr --K=2 --r0=0.5
