#!/bin/bash
#SBATCH --job-name=v2_main
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --requeue
#SBATCH --time=7-00:00:00    # PLACEHOLDER: set from a Roar pilot (02_fit_model.R --pilot=TRUE, README "Quick start") before the real run
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32   # = chains (4) x threads_per_chain (8)
#SBATCH --mem=128G                                 # basic bills 4 GB per core: 32 cores x 4 GB
#SBATCH --array=0-95
#SBATCH --output=%x_%A_%a.out
#SBATCH --open-mode=append   # a requeued task appends to its log instead of truncating it (2026-10-25 C7)
# MAIN-RESULTS GRID: sample {mid, fmid, war} x model {M1m, M1r, M1a, M1b, M2m, M2r, M2a, M2b} x lag {1, 3, 5, 10}
# = 96 cells of the MAIN specification (spec.R defaults: full all-states choice set, PageRank, K = 2 for M2, the sampler
# settings of SAMPLER per family), seed 123, within-chain parallelism via reduce_sum with THREADS_PER_CHAIN threads per
# chain (default 8 -> 32 cpus; override with --export=ALL,THREADS_PER_CHAIN=k together with --cpus-per-task=4k and --mem=<16k>G).
#
# Real run      : bash v2/slurm/submit_main_grid.sh
#                 = two arrays of this script: the 48 M1 cells (indices 0-11,24-35,48-59,72-83), then the 48 M2 cells
#                 (12-23,36-47,60-71,84-95) with --dependency=afterok on the M1 array, so each M2 task starts its chains at the
#                 posterior means of its cell's M1 fit (02_fit_model.R --init_from default; PI decision 2026-10-25). A plain
#                 `sbatch v2/slurm/02_main_grid.sbatch` runs all 96 at once; an M2 task whose M1 fit is not there yet falls
#                 back to random inits (logged loudly, recorded in the fit).
# Roar DRY RUN  : bash v2/slurm/submit_main_grid.sh --export=ALL,DRY_RUN=TRUE --time=00:30:00 --mem=16G --cpus-per-task=1
#                 (every task assembles and validates its Stan data with --dry_run=TRUE and exits; no sampling.
#                  Command-line sbatch options override the directives above, so the dry run is cheap.)
# The M2 cells run 4000 iterations (2000 warmup, adapt_delta 0.95; SAMPLER in spec.R) and need about twice the M1 walltime.
# Variants are NOT run from this script (their cells would carry a switch suffix); use 04_robustness.sbatch.
# Submit from the repo root; nothing in this file submits anything by itself.
# Slurm runs a spooled copy of this script, so locate the repo via the submit dir.
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
i=${SLURM_ARRAY_TASK_ID:?run as an array}
# ---- explicit case table: array index -> (sample, model, lag); lag slowest, then model, then sample
case $i in
   0) s=mid; m=M1m; l=1 ;;
   1) s=fmid; m=M1m; l=1 ;;
   2) s=war; m=M1m; l=1 ;;
   3) s=mid; m=M1r; l=1 ;;
   4) s=fmid; m=M1r; l=1 ;;
   5) s=war; m=M1r; l=1 ;;
   6) s=mid; m=M1a; l=1 ;;
   7) s=fmid; m=M1a; l=1 ;;
   8) s=war; m=M1a; l=1 ;;
   9) s=mid; m=M1b; l=1 ;;
  10) s=fmid; m=M1b; l=1 ;;
  11) s=war; m=M1b; l=1 ;;
  12) s=mid; m=M2m; l=1 ;;
  13) s=fmid; m=M2m; l=1 ;;
  14) s=war; m=M2m; l=1 ;;
  15) s=mid; m=M2r; l=1 ;;
  16) s=fmid; m=M2r; l=1 ;;
  17) s=war; m=M2r; l=1 ;;
  18) s=mid; m=M2a; l=1 ;;
  19) s=fmid; m=M2a; l=1 ;;
  20) s=war; m=M2a; l=1 ;;
  21) s=mid; m=M2b; l=1 ;;
  22) s=fmid; m=M2b; l=1 ;;
  23) s=war; m=M2b; l=1 ;;
  24) s=mid; m=M1m; l=3 ;;
  25) s=fmid; m=M1m; l=3 ;;
  26) s=war; m=M1m; l=3 ;;
  27) s=mid; m=M1r; l=3 ;;
  28) s=fmid; m=M1r; l=3 ;;
  29) s=war; m=M1r; l=3 ;;
  30) s=mid; m=M1a; l=3 ;;
  31) s=fmid; m=M1a; l=3 ;;
  32) s=war; m=M1a; l=3 ;;
  33) s=mid; m=M1b; l=3 ;;
  34) s=fmid; m=M1b; l=3 ;;
  35) s=war; m=M1b; l=3 ;;
  36) s=mid; m=M2m; l=3 ;;
  37) s=fmid; m=M2m; l=3 ;;
  38) s=war; m=M2m; l=3 ;;
  39) s=mid; m=M2r; l=3 ;;
  40) s=fmid; m=M2r; l=3 ;;
  41) s=war; m=M2r; l=3 ;;
  42) s=mid; m=M2a; l=3 ;;
  43) s=fmid; m=M2a; l=3 ;;
  44) s=war; m=M2a; l=3 ;;
  45) s=mid; m=M2b; l=3 ;;
  46) s=fmid; m=M2b; l=3 ;;
  47) s=war; m=M2b; l=3 ;;
  48) s=mid; m=M1m; l=5 ;;
  49) s=fmid; m=M1m; l=5 ;;
  50) s=war; m=M1m; l=5 ;;
  51) s=mid; m=M1r; l=5 ;;
  52) s=fmid; m=M1r; l=5 ;;
  53) s=war; m=M1r; l=5 ;;
  54) s=mid; m=M1a; l=5 ;;
  55) s=fmid; m=M1a; l=5 ;;
  56) s=war; m=M1a; l=5 ;;
  57) s=mid; m=M1b; l=5 ;;
  58) s=fmid; m=M1b; l=5 ;;
  59) s=war; m=M1b; l=5 ;;
  60) s=mid; m=M2m; l=5 ;;
  61) s=fmid; m=M2m; l=5 ;;
  62) s=war; m=M2m; l=5 ;;
  63) s=mid; m=M2r; l=5 ;;
  64) s=fmid; m=M2r; l=5 ;;
  65) s=war; m=M2r; l=5 ;;
  66) s=mid; m=M2a; l=5 ;;
  67) s=fmid; m=M2a; l=5 ;;
  68) s=war; m=M2a; l=5 ;;
  69) s=mid; m=M2b; l=5 ;;
  70) s=fmid; m=M2b; l=5 ;;
  71) s=war; m=M2b; l=5 ;;
  72) s=mid; m=M1m; l=10 ;;
  73) s=fmid; m=M1m; l=10 ;;
  74) s=war; m=M1m; l=10 ;;
  75) s=mid; m=M1r; l=10 ;;
  76) s=fmid; m=M1r; l=10 ;;
  77) s=war; m=M1r; l=10 ;;
  78) s=mid; m=M1a; l=10 ;;
  79) s=fmid; m=M1a; l=10 ;;
  80) s=war; m=M1a; l=10 ;;
  81) s=mid; m=M1b; l=10 ;;
  82) s=fmid; m=M1b; l=10 ;;
  83) s=war; m=M1b; l=10 ;;
  84) s=mid; m=M2m; l=10 ;;
  85) s=fmid; m=M2m; l=10 ;;
  86) s=war; m=M2m; l=10 ;;
  87) s=mid; m=M2r; l=10 ;;
  88) s=fmid; m=M2r; l=10 ;;
  89) s=war; m=M2r; l=10 ;;
  90) s=mid; m=M2a; l=10 ;;
  91) s=fmid; m=M2a; l=10 ;;
  92) s=war; m=M2a; l=10 ;;
  93) s=mid; m=M2b; l=10 ;;
  94) s=fmid; m=M2b; l=10 ;;
  95) s=war; m=M2b; l=10 ;;
  *) echo "array index $i out of range 0-95" >&2; exit 1 ;;
esac
DRY=${DRY_RUN:-FALSE}
THR=${THREADS_PER_CHAIN:-8}
echo "task $i -> sample=$s model=$m lag=$l (main specification; sampler defaults of spec.R) threads_per_chain=$THR cpus=${SLURM_CPUS_PER_TASK:-?} dry_run=$DRY"
Rscript v2/R/02_fit_model.R --sample=$s --lag=$l --model=$m --dry_run=$DRY --threads_per_chain=$THR
