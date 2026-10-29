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
#                  Default cell list = ROBUSTNESS_CELLS in spec.R (lag 3, M1b and M2b, three samples; 6 tasks). Submit a refit
#                  preset with submit_robustness.sh (below): it runs the M1 cells first and the M2 cells afterok, so every M2 cell
#                  starts from the M1 fit of its cell and switches (02_fit_model.R --init_from default). A plain sbatch of this file
#                  starts all six cells at once and its M2 cells fall back to random inits (an M1 fit that exists but is unusable
#                  stops the task). CELLS=sample:lag:model,... overrides the list (through the ENVIRONMENT, not --export=ALL,CELLS=...,
#                  because the commas of the list are --export separators): CELLS=mid:3:M1b,war:5:M1b sbatch --export=ALL ... --array=0-1
#   context        the context split (dual typology stage a): M1b, lag 3, mid and fmid, for context=bilateral (_bi), context=coalition
#                  (_co) and coalition without the top-3 coalition disputes (_co_drop3); default tie convention; 6 tasks (array 0-5:
#                  0-1 bilateral, 2-3 coalition, 4-5 coalition_drop3; cell order mid, fmid). M1b only: no dependency needed.
#                  Post: 06_context_compare.R per sample.
#   context_b      dual typology stage b (2026-10-27): the M1c model (stan/m1_context.stan: group-specific ideal-point shifts delta_lead /
#                  delta_join and tolerance scales sigma_g over bilateral / coalition leader / joiner events, shared rho, multi_mid_jan in W)
#                  at lag 3 on mid and fmid, every switch at its default; 2 tasks (array 0-1: 0 mid, 1 fmid; cells <sample>_lag3_all_pr_M1c).
#                  M1 family only (M1 sampler settings), no dependency, no war. Post: 03_extract.R --fit=v2/output/fit_<sample>_lag3_all_pr_M1c.rds
#   major_ix       the PI-approved stage-2 switches of 2026-10-27: M1b, lag 3, fmid and mid, for major=joint (_major), interact=prefixed (_ix)
#                  and both (_major_ix); 6 tasks (array 0-5: 0-1 major, 2-3 interact, 4-5 both; cell order fmid, mid). M1b only: no
#                  dependency needed. Post: 03_post.sbatch is for whole groups; read the fits with 03_extract.R --fit=<file> per cell.
#   dist           PI request 2026-10-28 (distance and distance + major powers; contiguity stays in): M1b, lag 3, fmid and mid, for
#                  distance=capital (_dist) and distance=capital + major=joint (_dist_major); 4 tasks (array 0-3: 0-1 distance, 2-3
#                  distance_major; cell order fmid, mid). M1b only: no dependency needed. Submit with --array=0-3 (submit_robustness.sh does).
#   dist_nocontig  PI request 2026-10-28, same day (the dist cells WITHOUT contiguity; distance replaces it): M1b, lag 3, fmid and mid, for
#                  distance=capital + contiguity=none (_dist_nocontig) and + major=joint (_dist_nocontig_major); 4 tasks with the SAME index
#                  mapping as dist (array 0-3: 0-1 distance, 2-3 distance_major; cell order fmid, mid), so the two compare task by task.
#                  M1b only: no dependency needed. Submit with --array=0-3 (submit_robustness.sh does).
#   perm_error     the permutation, measurement-error and temporal-holdout tests on the main cell (ROBUSTNESS_MAIN_CELL: mid, lag 5),
#                  the task layout of the former 06_measurement.sbatch (restored 2026-10-26, audit fix 6):
#                    0 = m1 permutation + inflation      1 = m2 permutation + inflation
#                    2 = m1 holdout cutoff 1990          3 = m2 holdout cutoff 1990
#                    4 = m1 holdout cutoff 2000          5 = m2 holdout cutoff 2000
#                  submit with --account=open --cpus-per-task=4 --mem=32G --time=47:00:00 (the default array 0-5).
#                  NEEDS the production a / b fits of BOTH families of the mid lag-5 cell from the main grid (fit_mid_lag5_all_pr_M1a,
#                  _M1b, _M2a_K2, _M2b_K2): a lag-3-first main run (M1_IDX=24-35 M2_IDX=36-47) does NOT produce them; run the lag-5
#                  block first (M1_IDX=48-59 M2_IDX=60-71 bash v2/slurm/submit_main_grid.sh) or the tasks stop at "a / b fits missing".
#                  PERM_MODE=refit, N_PERM=20 give the calibrated p-value instead of the rescore diagnostic; refit mode refits 2 x N_PERM
#                  models per family and will NOT finish within the 48 h, preemptible open allocation: submit it on basic
#                  (--account=jfe4_cr_default --partition=basic, several days of walltime). Since 2026-10-29 the permutation test of ONE
#                  cell (any sample / lag, --distance=capital allowed) runs through 08_perm_refit.sbatch / submit_perm_refit.sh instead:
#                  the refit mode as an array of one permutation per task plus a combine job, the rescore mode as one open-account job.
# Examples (submit from the repo root):
#   PRESET=noww2 bash v2/slurm/submit_robustness.sh                                          # M1 cells, then M2 cells afterok
#   PRESET=sideA bash v2/slurm/submit_robustness.sh                                          # needs v2/data_sideA (01_build.sbatch with DIRECTION=sideA)
#   PRESET=context bash v2/slurm/submit_robustness.sh                                        # 6 tasks, then Rscript v2/R/06_context_compare.R --sample=mid --lag=3
#   PRESET=context_b bash v2/slurm/submit_robustness.sh                                      # 2 M1c tasks (mid, fmid; lag 3), then 03_extract.R on each fit
#   PRESET=context_b bash v2/slurm/submit_robustness.sh --account=open --export=ALL,DRY_RUN=TRUE --time=00:30:00 --mem=16G --cpus-per-task=1   # assembly check
#   PRESET=major_ix bash v2/slurm/submit_robustness.sh                                       # 6 M1b tasks (fmid, mid x _major, _ix, _major_ix)
#   PRESET=dist bash v2/slurm/submit_robustness.sh                                           # 4 M1b tasks (fmid, mid x _dist, _dist_major)
#   PRESET=dist_nocontig bash v2/slurm/submit_robustness.sh                                  # 4 M1b tasks (fmid, mid x _dist_nocontig, _dist_nocontig_major)
#   PRESET=noww2 bash v2/slurm/submit_robustness.sh --account=open --export=ALL,DRY_RUN=TRUE --time=00:30:00 --mem=16G --cpus-per-task=1
#   sbatch --export=ALL,PRESET=perm_error --account=open --cpus-per-task=4 --mem=32G --time=47:00:00 v2/slurm/04_robustness.sbatch
#   sbatch --export=ALL,PRESET=perm_error,PERM_MODE=refit,N_PERM=20 --cpus-per-task=4 --mem=32G --time=7-00:00:00 v2/slurm/04_robustness.sbatch   # basic
# Post-processing of a refit preset: 03_post.sbatch with --export=ALL,SWITCHES="--ww2=drop" (the preset's switch).
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
i=${SLURM_ARRAY_TASK_ID:?run as an array}
P=${PRESET:?set --export=ALL,PRESET=<noww2|sideA|noweights|fullcontrols|winsor|context|context_b|major_ix|dist|dist_nocontig|perm_error>}
THR=${THREADS_PER_CHAIN:-8}; DRY=${DRY_RUN:-FALSE}
case $P in
  perm_error)
    FAMS=(m1 m2); fam=${FAMS[$(( i % 2 ))]}
    PM=${PERM_MODE:-rescore}; NP=${N_PERM:-50}
    COMMON="--preset=perm_error --family=$fam --r0=0.5 --chains=4 --iter=1500 --warmup=750"
    case $i in
      0|1) echo "task $i: perm_error, family $fam permutation ($PM, n_perm $NP) + inflation"
           Rscript v2/R/05_robustness.R $COMMON --tests=permutation,inflation --n_perm=$NP --perm_mode=$PM --deltas=0.25,0.5,1.0 ;;
      2|3) echo "task $i: perm_error, family $fam temporal holdout, cutoff 1990"
           Rscript v2/R/05_robustness.R $COMMON --tests=holdout --cutoff=1990 ;;
      4|5) echo "task $i: perm_error, family $fam temporal holdout, cutoff 2000"
           Rscript v2/R/05_robustness.R $COMMON --tests=holdout --cutoff=2000 ;;
      *) echo "perm_error: array index $i out of range 0-5" >&2; exit 1 ;;
    esac ;;
  noww2|sideA|noweights|fullcontrols|winsor)
    echo "task $i: preset $P, cell index $i of ${CELLS:-the default list} threads_per_chain=$THR dry_run=$DRY"
    Rscript v2/R/05_robustness.R --preset=$P --cell_index=$i ${CELLS:+--cells=$CELLS} --dry_run=$DRY --threads_per_chain=$THR ;;
  context)
    echo "task $i: preset context (variant x cell index $i of 6: 0-1 bilateral, 2-3 coalition, 4-5 coalition_drop3; mid, fmid) threads_per_chain=$THR dry_run=$DRY"
    Rscript v2/R/05_robustness.R --preset=context --cell_index=$i --dry_run=$DRY --threads_per_chain=$THR ;;
  context_b)
    [ "$i" -le 1 ] || { echo "context_b: array index $i out of range 0-1 (0 mid, 1 fmid)" >&2; exit 1; }
    echo "task $i: preset context_b (the M1c model, lag 3; cell index $i of 2: 0 mid, 1 fmid) threads_per_chain=$THR dry_run=$DRY"
    Rscript v2/R/05_robustness.R --preset=context_b --cell_index=$i --dry_run=$DRY --threads_per_chain=$THR ;;
  major_ix)
    echo "task $i: preset major_ix (variant x cell index $i of 6: 0-1 major, 2-3 interact, 4-5 both; fmid, mid) threads_per_chain=$THR dry_run=$DRY"
    Rscript v2/R/05_robustness.R --preset=major_ix --cell_index=$i --dry_run=$DRY --threads_per_chain=$THR ;;
  dist)
    [ "$i" -le 3 ] || { echo "dist: array index $i out of range 0-3 (0-1 distance, 2-3 distance_major; fmid, mid)" >&2; exit 1; }
    echo "task $i: preset dist (variant x cell index $i of 4: 0-1 distance, 2-3 distance_major; fmid, mid) threads_per_chain=$THR dry_run=$DRY"
    Rscript v2/R/05_robustness.R --preset=dist --cell_index=$i --dry_run=$DRY --threads_per_chain=$THR ;;
  dist_nocontig)
    [ "$i" -le 3 ] || { echo "dist_nocontig: array index $i out of range 0-3 (0-1 distance, 2-3 distance_major; fmid, mid)" >&2; exit 1; }
    echo "task $i: preset dist_nocontig (variant x cell index $i of 4: 0-1 distance, 2-3 distance_major; fmid, mid; contiguity=none) threads_per_chain=$THR dry_run=$DRY"
    Rscript v2/R/05_robustness.R --preset=dist_nocontig --cell_index=$i --dry_run=$DRY --threads_per_chain=$THR ;;
  *) echo "unknown PRESET $P" >&2; exit 1 ;;
esac
