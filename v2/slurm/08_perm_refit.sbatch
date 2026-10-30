#!/bin/bash
#SBATCH --job-name=v2_perm
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --requeue
#SBATCH --time=2-00:00:00    # one refit task = 2 refits (a, b) with 8 threads per chain, ~15-25 h expected; see "Resources" below; override from the command line
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32   # = chains (4) x threads_per_chain (8): the refits run reduce_sum over state-years (2026-10-30)
#SBATCH --mem=128G                                 # basic bills 4 GB per core: 32 cores x 4 GB (the frames + the two production fits + two refits need far less)
#SBATCH --output=%x_%A_%a.out
#SBATCH --open-mode=append
# PERMUTATION TEST OF ONE CELL (05_robustness.R --preset=perm_error, --tests=permutation), 2026-10-29: the PI's request to run the
# calibrated permutation test on the NEW main cell (fatal MIDs, lag 3, M1 family, contiguity + capital distance: fit_fmid_lag3_all_pr_M1a_dist /
# _M1b_dist), with the refit mode split over a Slurm array instead of one multi-day job. Submit with submit_perm_refit.sh (below); this file
# is the task body. MODE (required, through the environment):
#   rescore   ONE task: the reliance diagnostic (the original fits rescored on N_PERM permuted data sets; no refit, no p-value).
#             Fast enough for the open account: --account=open --time=24:00:00 --cpus-per-task=4 --mem=32G (submit_perm_refit.sh does that).
#   refit     an ARRAY 1..N_PERM: task p refits the a and b models ONCE each on the data with the relational gap shuffled with seed
#             V2_SEED + p (chains 4, iter 1500, warmup 750, the settings of the sequential refit mode; CHAINS / ITER / WARMUP override
#             them for every task), each chain with THREADS_PER_CHAIN (8) reduce_sum threads (2026-10-30; --threads_per_chain=k of the
#             driver: rstan_options(threads_per_chain = k) and grainsize = choose_grainsize(N1, k), as 02_fit_model.R), and writes
#             output/measurement_<cell>_permutation_refit_parts/perm_<pppp>.csv (stage-2 LOO elpd of a and b, with the threads and grainsize).
#   combine   ONE task, afterok on the array: reads the N_PERM parts (refuses a missing, duplicated or extra index, or parts with other
#             settings, the threads / grainsize among them, so it is passed the SAME THREADS_PER_CHAIN), recomputes the observed statistic
#             from the a / b fits, writes measurement_<cell>_permutation_refit.csv and _draws.csv with p = (1 + exceedances) / (1 + N_PERM),
#             exactly the sequential output. Open account, 4 h, 4 cores (it samples nothing).
# Cell (environment, defaults = the PI's cell): SAMPLE=fmid LAG=3 FAMILY=m1 SWITCHES="--distance=capital" N_PERM=20. SWITCHES is the spec.R
# switch list of the cell (the cell name carries its suffix: fmid_lag3_all_pr_M1_dist); SWITCHES="--distance=capital --major=joint" is the
# distance + major-power cell (fmid_lag3_all_pr_M1_dist_major, fits fit_fmid_lag3_all_pr_M1a_dist_major / _M1b_dist_major; 2026-10-30); the
# driver refuses what the permutation test does not take (contiguity=none, major without the distance, interact, winsor, pairs != all).
# SWITCHES="" runs the main-spec cell of SAMPLE / LAG.
# Resources of a refit task (why 32 cores, 128 GB, 2 days): each refit samples rstan::sampling(chains = 4) with mc.cores = 4 and
# THREADS_PER_CHAIN (8) reduce_sum threads per chain, so a task uses chains x threads = 32 cores (submit_perm_refit.sh passes
# --cpus-per-task and --mem for CHAINS x THREADS_PER_CHAIN when they differ from the defaults). Memory: basic bills 4 GB per core, so 32 cores
# are 128 GB whatever the job needs; the assembled fmid lag-3 frames (1.72 M alternative rows), the two production fits (4000 draws x 14,597
# stage-1 log-lik each) and the two refits (3000 draws each) fit in 32 GB. Walltime: one refit is 1500 iterations with 8 threads per chain
# against a main-grid M1 fit of 2000 iterations with the same threads; take the "sampling finished in" line of
# v2/logs/02_fit_fmid_lag3_all_pr_M1b_dist_*.log (T_main), then one refit ~ 0.75 x T_main and one task (a + b, plus two LOOs) ~ 1.5-2 x T_main,
# i.e. about 15-25 h for the lag-3 fmid and mid cells; the 2-day directive leaves headroom. Pass --time=<d-hh:mm:ss> on the command line to
# change it. The array runs N_PERM tasks concurrently as the queue allows (--array=1-N_PERM%k limits it).
# Examples (from the repo root):
#   MODE=rescore bash v2/slurm/submit_perm_refit.sh                                  # fmid lag 3 m1 _dist, N_PERM (20) permutations, open account
#   MODE=rescore N_PERM=50 bash v2/slurm/submit_perm_refit.sh                        # the 50 of the old perm_error default
#   MODE=refit bash v2/slurm/submit_perm_refit.sh                                    # array 1-20 on basic (32 cores, 2 days per task), then the combine job afterok (open)
#   MODE=refit SAMPLE=mid bash v2/slurm/submit_perm_refit.sh                         # the same on mid
#   MODE=refit SWITCHES="--distance=capital --major=joint" bash v2/slurm/submit_perm_refit.sh   # the distance + major cell (fits ..._M1a_dist_major / _M1b_dist_major)
#   MODE=refit THREADS_PER_CHAIN=4 bash v2/slurm/submit_perm_refit.sh                # 4 threads per chain: 16 cores, 64 GB per task
#   MODE=refit N_PERM=20 bash v2/slurm/submit_perm_refit.sh --time=1-00:00:00        # a shorter walltime once T_main is known
#   MODE=combine bash v2/slurm/submit_perm_refit.sh                                  # re-run the combiner alone (every part 1..N_PERM must exist; same THREADS_PER_CHAIN)
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
M=${MODE:?set MODE=<rescore|refit|combine> (through submit_perm_refit.sh)}
SAMPLE=${SAMPLE:-fmid}; LAG=${LAG:-3}; FAMILY=${FAMILY:-m1}; NP=${N_PERM:-20}
SW=${SWITCHES-"--distance=capital"}          # SWITCHES="" (set and empty) = the main-spec cell; unset = the distance cell
CH=${CHAINS:-4}; IT=${ITER:-1500}; WU=${WARMUP:-750}
THR=${THREADS_PER_CHAIN:-8}                  # reduce_sum threads per refit chain (refit and combine: the parts carry it; rescore samples nothing)
[[ "$THR" =~ ^[0-9]+$ ]] && [ "$THR" -ge 1 ] || { echo "THREADS_PER_CHAIN must be a positive integer (got '$THR')" >&2; exit 1; }
COMMON="--preset=perm_error --family=$FAMILY --sample=$SAMPLE --lag=$LAG --r0=0.5 --chains=$CH --iter=$IT --warmup=$WU --n_perm=$NP $SW"
case $M in
  rescore)
    echo "perm rescore: $SAMPLE lag $LAG $FAMILY [$SW], n_perm $NP (reliance diagnostic, no refit)"
    Rscript v2/R/05_robustness.R $COMMON --tests=permutation --perm_mode=rescore ;;
  refit)
    p=${SLURM_ARRAY_TASK_ID:?MODE=refit runs as an array 1..N_PERM (task p = permutation p)}
    [ "$p" -ge 1 ] && [ "$p" -le "$NP" ] || { echo "refit: array index $p outside 1..$NP (N_PERM)" >&2; exit 1; }
    echo "perm refit task $p of $NP: $SAMPLE lag $LAG $FAMILY [$SW], seed 123 + $p, chains $CH iter $IT warmup $WU, $THR threads per chain (cpus=${SLURM_CPUS_PER_TASK:-?})"
    Rscript v2/R/05_robustness.R $COMMON --tests=permutation --perm_mode=refit --perm_index=$p --threads_per_chain=$THR ;;
  combine)
    echo "perm combine: $SAMPLE lag $LAG $FAMILY [$SW], parts 1..$NP (threads per chain $THR) -> measurement_<cell>_permutation_refit.csv"
    Rscript v2/R/05_robustness.R $COMMON --tests=permutation_combine --perm_mode=refit --threads_per_chain=$THR ;;
  *) echo "unknown MODE $M (rescore | refit | combine)" >&2; exit 1 ;;
esac
