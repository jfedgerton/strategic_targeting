#!/bin/bash
#SBATCH --job-name=v2_robust
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --requeue
#SBATCH --time=7-00:00:00    # PLACEHOLDER like 02_main_grid.sbatch: set from the main-grid timings
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32   # = chains (4) x threads_per_chain (8)
#SBATCH --mem=128G                                 # basic bills 4 GB per core
#SBATCH --array=0-5
#SBATCH --output=%x_%A_%a.out
#SBATCH --open-mode=append
# ROBUSTNESS CHECKS, one preset per submission (05_robustness.R; presets in spec.R, ROBUSTNESS_PRESETS):
#   refit presets  noww2 | sideA | noweights | fullcontrols | winsor : each array task refits ONE cell of the cell list with the
#                  preset's switch flipped (the fit carries the switch suffix, so main-grid outputs are never overwritten).
#                  Default cell list = ROBUSTNESS_CELLS in spec.R (lag 3, M1b and M2b, three samples; 6 tasks); override with
#                  --export=ALL,PRESET=noww2,CELLS=mid:3:M1b,war:5:M2b and --array=0-<n-1>.
#   context        the context split (dual typology stage a): M1b, lag 3, mid and fmid, for context=bilateral (_bi), context=coalition
#                  (_co) and coalition without the top-3 coalition disputes (_co_drop3); default tie convention; 6 tasks (array 0-5:
#                  0-1 bilateral, 2-3 coalition, 4-5 coalition_drop3; cell order mid, fmid). Post: 06_context_compare.R per sample.
#   perm_error     the permutation and measurement-error tests on the main cell (mid, lag 5): task 0 = family m1, task 1 = m2;
#                  submit with --array=0-1 --account=open --cpus-per-task=4 --mem=32G --time=47:00:00. Needs the a / b fits of
#                  02_main_grid.sbatch. PERM_MODE=refit, N_PERM=20 give the calibrated p-value instead of the rescore diagnostic.
# Examples (submit from the repo root):
#   sbatch --export=ALL,PRESET=noww2 v2/slurm/04_robustness.sbatch
#   sbatch --export=ALL,PRESET=sideA v2/slurm/04_robustness.sbatch          # needs v2/data_sideA (01_build.sbatch with DIRECTION=sideA)
#   sbatch --export=ALL,PRESET=context v2/slurm/04_robustness.sbatch        # 6 tasks (the default array 0-5), then Rscript v2/R/06_context_compare.R --sample=mid --lag=3
#   sbatch --export=ALL,PRESET=noww2,DRY_RUN=TRUE --time=00:30:00 --mem=16G --cpus-per-task=1 v2/slurm/04_robustness.sbatch
#   sbatch --export=ALL,PRESET=perm_error --array=0-1 --account=open --cpus-per-task=4 --mem=32G --time=47:00:00 v2/slurm/04_robustness.sbatch
# Post-processing of a refit preset: 03_post.sbatch with --export=ALL,SWITCHES="--ww2=drop" (the preset's switch).
# M2 cells start from the M1 fit of the same cell and switches when it exists in v2/output (02_fit_model.R --init_from default), so
# submit the M1 cells of a preset first (e.g. CELLS=mid:3:M1b,fmid:3:M1b,war:3:M1b --array=0-2) and the M2 cells with
# --dependency=afterok:<that job>; an M2 cell without its M1 fit falls back to random inits (logged loudly, recorded in the fit).
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
i=${SLURM_ARRAY_TASK_ID:?run as an array}
P=${PRESET:?set --export=ALL,PRESET=<noww2|sideA|noweights|fullcontrols|winsor|context|perm_error>}
THR=${THREADS_PER_CHAIN:-8}; DRY=${DRY_RUN:-FALSE}
case $P in
  perm_error)
    FAMS=(m1 m2); fam=${FAMS[$(( i % 2 ))]}
    echo "task $i: perm_error, family $fam (perm_mode ${PERM_MODE:-rescore}, n_perm ${N_PERM:-50})"
    Rscript v2/R/05_robustness.R --preset=perm_error --family=$fam --tests=permutation,inflation \
      --n_perm=${N_PERM:-50} --perm_mode=${PERM_MODE:-rescore} --deltas=0.25,0.5,1.0 --r0=0.5 --chains=4 --iter=1500 --warmup=750 ;;
  noww2|sideA|noweights|fullcontrols|winsor)
    echo "task $i: preset $P, cell index $i of ${CELLS:-the default list} threads_per_chain=$THR dry_run=$DRY"
    Rscript v2/R/05_robustness.R --preset=$P --cell_index=$i ${CELLS:+--cells=$CELLS} --dry_run=$DRY --threads_per_chain=$THR ;;
  context)
    echo "task $i: preset context (variant x cell index $i of 6: 0-1 bilateral, 2-3 coalition, 4-5 coalition_drop3; mid, fmid) threads_per_chain=$THR dry_run=$DRY"
    Rscript v2/R/05_robustness.R --preset=context --cell_index=$i --dry_run=$DRY --threads_per_chain=$THR ;;
  *) echo "unknown PRESET $P" >&2; exit 1 ;;
esac
