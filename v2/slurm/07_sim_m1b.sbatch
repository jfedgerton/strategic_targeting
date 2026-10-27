#!/bin/bash
#SBATCH --job-name=v2_sim_m1b
#SBATCH --account=open
#SBATCH --partition=basic   # the open partition no longer exists; --account=open bills the free READ credits
#SBATCH --requeue
#SBATCH --time=12:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=0-19
#SBATCH --output=%x_%A_%a.out
#SBATCH --open-mode=append
# M1b simulation recovery (PI request 2026-10-27, after the M1c recovery bias): 20 replicates of sim/sim_context_b.R --model=M1b on the
# observed fmid lag-3 M1b design (the main-specification cell: W without multi_mid_jan, no groups), replicate r = SLURM_ARRAY_TASK_ID with
# seed 123 + r, y1 / event slots / targets all simulated from the M1b model at the M1c bilateral truth, 4 chains x 2000 iterations
# (1000 warmup). Each task writes v2/sim/results/m1b_recovery/m1b_recovery_r<r>.csv and m1b_design_r<r>.csv (never the M1c files);
# afterwards `Rscript v2/sim/summarise_context_b.R --model=M1b` writes the coverage table v2/sim/results/m1b_coverage.csv (needs >= 20
# replicates unless --min_reps is lowered). Requires the data (01_build_data.R) and the compiled programs (compile_models.R) on the
# submitting checkout. Nothing in this file submits anything by itself; run from the repo root with: sbatch v2/slurm/07_sim_m1b.sbatch
# Subsample sizes: --n_init / --n_noninit default to 250 / 250 (sim_context_b.R); N_INIT / N_NONINIT in the environment override them.
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
r=${SLURM_ARRAY_TASK_ID:?run as an array}
Rscript v2/sim/sim_context_b.R --model=M1b --rep="$r" --n_init="${N_INIT:-250}" --n_noninit="${N_NONINIT:-250}" --chains=4 --iter=2000 --warmup=1000
