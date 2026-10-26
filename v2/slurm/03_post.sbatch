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
# Variant groups: --export=ALL,SWITCHES="--ww2=drop" runs the same for the _noww2 fits (any switches of spec.R; the fit names come
# from cell_name() in spec.R, so the K suffix and the switch suffixes are never hand-written here).
# Audit fix 4 (2026-10-26): the task FAILS (exit 1) unless all eight fits of the group are present, so 04_compare.R never writes
# compare_* CSVs of a partial group unnoticed; --export=ALL,ALLOW_PARTIAL=TRUE runs on the fits that exist (04_compare.R skips the
# missing ones and says so). A SWITCHES value spec.R rejects, a failed extraction or a failed comparison also exit 1 (set -e).
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
set -eo pipefail
i=${SLURM_ARRAY_TASK_ID:?run as an array}
SAMPLES=(mid fmid war); LAGS=(1 3 5 10)
[ "$i" -ge 0 ] && [ "$i" -le 11 ] || { echo "array index $i out of range 0-11" >&2; exit 1; }
s=${SAMPLES[$(( i % 3 ))]}; l=${LAGS[$(( i / 3 ))]}
SW=${SWITCHES:-}; PARTIAL=${ALLOW_PARTIAL:-FALSE}
echo "task $i -> sample=$s lag=$l switches='$SW' allow_partial=$PARTIAL"
# the eight cell names of the group from spec.R itself (parse_spec() stops on a bad switch value, which stops this task)
CELLS=$(Rscript -e 'source("v2/R/spec.R"); a <- commandArgs(TRUE); spec <- parse_spec(a); s <- get_arg(a, "sample", ""); l <- as.integer(get_arg(a, "lag", 0)); writeLines(vapply(MODELS, function(m) cell_name(s, l, m, spec), ""))' --sample=$s --lag=$l $SW) \
  || { echo "ERROR: spec.R rejected SWITCHES='$SW' (or the sample / lag); nothing run" >&2; exit 1; }
n_cells=$(wc -l <<< "$CELLS"); [ "$n_cells" -eq 8 ] || { echo "ERROR: expected 8 cell names, got $n_cells" >&2; exit 1; }
found=0; missing=""
for c in $CELLS; do
  f="v2/output/fit_${c}.rds"
  if [ -f "$f" ]; then found=$(( found + 1 )); else missing="${missing:+$missing }$f"; fi
done
echo "fits present: $found of 8${missing:+; missing: $missing}"
if [ "$found" -lt 8 ]; then
  if [ "$PARTIAL" != "TRUE" ]; then echo "ERROR: $(( 8 - found )) fit(s) of group ${s}_lag${l} missing; set --export=ALL,ALLOW_PARTIAL=TRUE to post-process the $found present" >&2; exit 1; fi
  [ "$found" -gt 0 ] || { echo "ERROR: no fit of the group present; nothing to do" >&2; exit 1; }
  echo "WARNING: ALLOW_PARTIAL=TRUE; post-processing $found of 8 fits (the compare_* tables cover these only)"
fi
for c in $CELLS; do
  f="v2/output/fit_${c}.rds"
  [ -f "$f" ] || continue
  Rscript v2/R/03_extract.R --fit="$f" --r0=0.5
done
Rscript v2/R/04_compare.R --sample=$s --lag=$l --r0=0.5 $SW
echo "task $i done: group ${s}_lag${l}${SW:+ ($SW)} with $found fits"
