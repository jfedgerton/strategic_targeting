#!/bin/bash
#SBATCH --job-name=v2_sim
#SBATCH --account=open
#SBATCH --partition=basic   # the open partition no longer exists; --account=open bills the free READ credits
#SBATCH --requeue
#SBATCH --time=24:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-4
#SBATCH --output=%x_%A_%a.out
#SBATCH --open-mode=append
# Simulation rerun under the CURRENT model (2026-10-25, B5): mu0, gap scaling, in_iv / in_s2 masks, minimal controls with the power
# levels and log system size, event weights, r0 = 0.5, x_ref = initiator mean; plus the persistent-position variant (random-walk
# step SD 0.06 instead of 0.3). sim/SIM_REPORT.md is marked superseded until this has run and write_sim_report.R is rerun.
#   0 = run_sim.R (M1a, M1b, M2a, M2b on M1b and M2b truths)       1 = run_sim.R --persistent (persistent positions)
#   2 = run_sim_1d.R                                                3 = run_measurement_tests.R
#   4 = run_benchmarks_sim.R (brms; separate task because of its compile memory)
# Nothing in this file submits anything by itself; run from the repo root with sbatch v2/slurm/10_simulations.sbatch (open account).
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
i=${SLURM_ARRAY_TASK_ID:?run as an array}
case $i in
  0) Rscript v2/sim/run_sim.R 4 2000 1000 ;;
  1) Rscript v2/sim/run_sim.R 4 2000 1000 "" persistent ;;
  2) Rscript v2/sim/run_sim_1d.R 4 2000 1000 ;;
  3) Rscript v2/sim/run_measurement_tests.R 4 1500 750 ;;
  4) Rscript v2/sim/run_benchmarks_sim.R 4 2000 1000 ;;
  *) echo "array index $i out of range 0-4" >&2; exit 1 ;;
esac
