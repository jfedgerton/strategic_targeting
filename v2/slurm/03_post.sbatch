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
# POST-PROCESSING, after 02_main_grid.sbatch: per (sample, lag) cell group, 03_extract.R (typology, mu_state, mu_event, the
# appendix 9-cell classification, reach and focus; AME alignment) for each of the eight fits, then 04_compare.R (LOO, joint vs
# independent, M1 vs M2, a-vs-b misclassification). 12 tasks = 3 samples x 4 lags.
# Submit from the repo root: sbatch v2/slurm/03_post.sbatch
# Variant groups: --export=ALL,SWITCHES="--ww2=drop" runs the same for the _noww2 fits (any switches of spec.R).
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
i=${SLURM_ARRAY_TASK_ID:?run as an array}
SAMPLES=(mid fmid war); LAGS=(1 3 5 10)
s=${SAMPLES[$(( i % 3 ))]}; l=${LAGS[$(( i / 3 ))]}
echo "task $i -> sample=$s lag=$l"
SW=${SWITCHES:-}
# the cell suffix of the switches (empty for the main specification), from spec.R itself
SUF=$(Rscript -e 'source("v2/R/spec.R"); cat(spec_suffix(parse_spec(commandArgs(TRUE))))' $SW)
for m in M1m M1r M1a M1b M2m M2r M2a M2b; do
  case $m in M2*) k="_K2" ;; *) k="" ;; esac
  f="v2/output/fit_${s}_lag${l}_all_pr_${m}${k}${SUF}.rds"
  if [ -f "$f" ]; then Rscript v2/R/03_extract.R --fit="$f" --r0=0.5; else echo "missing fit for $m: $f" >&2; fi
done
Rscript v2/R/04_compare.R --sample=$s --lag=$l --r0=0.5 $SW
