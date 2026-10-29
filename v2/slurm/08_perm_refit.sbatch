#!/bin/bash
#SBATCH --job-name=v2_perm
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --requeue
#SBATCH --time=4-00:00:00    # one refit task = 2 serial-chain refits (a, b); see "Resources" below; override from the command line
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4   # = chains (4) x 1 thread per chain: the refits of permutation_test() run serial chains (grainsize = N1)
#SBATCH --mem=32G                                 # the assembled frames (N2 = 1.72 M rows) + the two production fits + two refits; basic bills 4 GB per core
#SBATCH --output=%x_%A_%a.out
#SBATCH --open-mode=append
# PERMUTATION TEST OF ONE CELL (05_robustness.R --preset=perm_error, --tests=permutation), 2026-10-29: the PI's request to run the
# calibrated permutation test on the NEW main cell (fatal MIDs, lag 3, M1 family, contiguity + capital distance: fit_fmid_lag3_all_pr_M1a_dist /
# _M1b_dist), with the refit mode split over a Slurm array instead of one multi-day job. Submit with submit_perm_refit.sh (below); this file
# is the task body. MODE (required, through the environment):
#   rescore   ONE task: the reliance diagnostic (the original fits rescored on N_PERM permuted data sets; no refit, no p-value).
#             Fast enough for the open account: --account=open --time=24:00:00 (submit_perm_refit.sh does that).
#   refit     an ARRAY 1..N_PERM: task p refits the a and b models ONCE each on the data with the relational gap shuffled with seed
#             V2_SEED + p (chains 4, iter 1500, warmup 750, the settings of the sequential refit mode; CHAINS / ITER / WARMUP override
#             them for every task) and writes output/measurement_<cell>_permutation_refit_parts/perm_<pppp>.csv (stage-2 LOO elpd of a and b).
#   combine   ONE task, afterok on the array: reads the N_PERM parts (refuses a missing, duplicated or extra index, or parts with other
#             settings), recomputes the observed statistic from the a / b fits, writes measurement_<cell>_permutation_refit.csv and
#             _draws.csv with p = (1 + exceedances) / (1 + N_PERM), exactly the sequential output. Open account, 4 h.
# Cell (environment, defaults = the PI's cell): SAMPLE=fmid LAG=3 FAMILY=m1 SWITCHES="--distance=capital" N_PERM=20. SWITCHES is the spec.R
# switch list of the cell (the cell name carries its suffix: fmid_lag3_all_pr_M1_dist); the driver refuses what the permutation test does not take
# (contiguity=none, major, interact, winsor, pairs != all). SWITCHES="" runs the main-spec cell of SAMPLE / LAG.
# Resources of a refit task (why 4 cores, 32 GB, 4 days): permutation_test() samples each refit with rstan::sampling(chains = 4) and
# mc.cores = min(4, cores), i.e. four chains in parallel, each SERIAL (build_stan_data() leaves grainsize = N1, so reduce_sum runs one slice
# per chain whatever STAN_NUM_THREADS says), so a task uses 4 cores and no more. Memory: the assembled fmid lag-3 frames (1.72 M alternative
# rows), the two production fits (4000 draws x 14,597 stage-1 log-lik each) and the two refits (3000 draws each) fit in 32 GB with room.
# Walltime: one refit is 1500 iterations on serial chains against a main-grid M1 fit of 2000 iterations with 8 threads per chain; take the
# "sampling finished in" line of v2/logs/02_fit_fmid_lag3_all_pr_M1b_dist_*.log (T_main), then one refit ~ 0.75 x T_main x (the reduce_sum
# speed-up, 4-6 on 8 threads) ~ 3-4.5 x T_main and one task (a + b, plus two LOOs) ~ 6-9 x T_main. The 4-day directive is a ceiling; pass
# --time=<h> on the command line once T_main is known. The array runs N_PERM tasks concurrently as the queue allows (--array=1-N_PERM%k limits it).
# Examples (from the repo root):
#   MODE=rescore bash v2/slurm/submit_perm_refit.sh                                  # fmid lag 3 m1 _dist, N_PERM (20) permutations, open account
#   MODE=rescore N_PERM=50 bash v2/slurm/submit_perm_refit.sh                        # the 50 of the old perm_error default
#   MODE=refit bash v2/slurm/submit_perm_refit.sh                                    # array 1-20 on basic, then the combine job afterok (open)
#   MODE=refit SAMPLE=mid bash v2/slurm/submit_perm_refit.sh                         # the same on mid
#   MODE=refit N_PERM=20 bash v2/slurm/submit_perm_refit.sh --time=2-00:00:00        # a shorter walltime for the array once T_main is known
#   MODE=combine bash v2/slurm/submit_perm_refit.sh                                  # re-run the combiner alone (every part 1..N_PERM must exist)
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
M=${MODE:?set MODE=<rescore|refit|combine> (through submit_perm_refit.sh)}
SAMPLE=${SAMPLE:-fmid}; LAG=${LAG:-3}; FAMILY=${FAMILY:-m1}; NP=${N_PERM:-20}
SW=${SWITCHES-"--distance=capital"}          # SWITCHES="" (set and empty) = the main-spec cell; unset = the distance cell
CH=${CHAINS:-4}; IT=${ITER:-1500}; WU=${WARMUP:-750}
COMMON="--preset=perm_error --family=$FAMILY --sample=$SAMPLE --lag=$LAG --r0=0.5 --chains=$CH --iter=$IT --warmup=$WU --n_perm=$NP $SW"
case $M in
  rescore)
    echo "perm rescore: $SAMPLE lag $LAG $FAMILY [$SW], n_perm $NP (reliance diagnostic, no refit)"
    Rscript v2/R/05_robustness.R $COMMON --tests=permutation --perm_mode=rescore ;;
  refit)
    p=${SLURM_ARRAY_TASK_ID:?MODE=refit runs as an array 1..N_PERM (task p = permutation p)}
    [ "$p" -ge 1 ] && [ "$p" -le "$NP" ] || { echo "refit: array index $p outside 1..$NP (N_PERM)" >&2; exit 1; }
    echo "perm refit task $p of $NP: $SAMPLE lag $LAG $FAMILY [$SW], seed 123 + $p, chains $CH iter $IT warmup $WU, 4 serial chains"
    Rscript v2/R/05_robustness.R $COMMON --tests=permutation --perm_mode=refit --perm_index=$p ;;
  combine)
    echo "perm combine: $SAMPLE lag $LAG $FAMILY [$SW], parts 1..$NP -> measurement_<cell>_permutation_refit.csv"
    Rscript v2/R/05_robustness.R $COMMON --tests=permutation_combine --perm_mode=refit ;;
  *) echo "unknown MODE $M (rescore | refit | combine)" >&2; exit 1 ;;
esac
