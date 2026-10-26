#!/bin/bash
# Submit the main grid as two dependent arrays of 02_main_grid.sbatch (PI decision 2026-10-25): the 48 M1 cells first, then the
# 48 M2 cells with --dependency=afterok on the M1 array, so every M2 task finds the M1 fit of its cell in v2/output and starts
# its chains at that fit's posterior means (02_fit_model.R --init_from default; a missing M1 fit falls back to random inits).
# The array indices are those of 02_main_grid.sbatch (0-95: lag slowest, then model, then sample): in each lag block of 24,
# indices 0-11 are the M1 models and 12-23 the M2 models. One afterok on the whole M1 array is used instead of per-cell
# dependencies (simpler; the M2 array waits until every M1 cell has finished successfully).
# Usage (from the repo root): bash v2/slurm/submit_main_grid.sh [extra sbatch options, e.g. --export=ALL,DRY_RUN=TRUE --time=00:30:00]
set -eo pipefail
[ -f v2/slurm/02_main_grid.sbatch ] || { echo "run from the repo root" >&2; exit 1; }
M1_IDX="0-11,24-35,48-59,72-83"
M2_IDX="12-23,36-47,60-71,84-95"
m1=$(sbatch --parsable --array=$M1_IDX "$@" v2/slurm/02_main_grid.sbatch)
echo "M1 array submitted: job $m1 (indices $M1_IDX)"
m2=$(sbatch --parsable --array=$M2_IDX --dependency=afterok:$m1 "$@" v2/slurm/02_main_grid.sbatch)
echo "M2 array submitted: job $m2 (indices $M2_IDX), afterok:$m1"
