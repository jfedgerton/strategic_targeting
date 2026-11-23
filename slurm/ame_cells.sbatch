#!/bin/bash
#SBATCH --job-name=ame_cells
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --requeue
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32   # 4 chains x 8 threads per chain
#SBATCH --mem=128G                                 # basic bills 4 GB per core
#SBATCH --array=0-2
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --open-mode=append
## AME cells 1-3 (docs/ame_design.md rev 2; PI 2026-11-23): 0 = main (joint, K = 0), 1 = cut posterior, 2 = k1 (joint, K = 1).
## Fit-only: reads the data bundle written on the Mac by `Rscript R/ame/ame_cell.R --prep_only=1 --bundle=...`. Submit from the repo root:
##   BUNDLE=ame/ame_bundle_fmid_w1_book.rds sbatch slurm/ame_cells.sbatch
module purge; module load gcc/12.2.0; module load r/4.2.1-gcc-8.5.0
export LD_LIBRARY_PATH="$(dirname "$(gcc -print-file-name=libstdc++.so)"):${LD_LIBRARY_PATH}"
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/x86_64-pc-linux-gnu-library/4.2}"
export R_MAKEVARS_USER="$PWD/env/Makevars_roar"; export PKG_CPPFLAGS="${PKG_CPPFLAGS:-} -DSTAN_THREADS"; export OMP_NUM_THREADS=1
TH=8; export STAN_NUM_THREADS=$TH
BUNDLE=${BUNDLE:?set BUNDLE=path to the AME data bundle}; [ -s "$BUNDLE" ] || { echo "missing $BUNDLE" >&2; exit 1; }
CELLS=(main cut k1); C=${CELLS[$SLURM_ARRAY_TASK_ID]}
OUT=${OUT:-ame/out_$(basename "$BUNDLE" .rds)}; mkdir -p "$OUT" logs
echo "host $(hostname); R $(R --version | head -1); cell $C; bundle $BUNDLE; out $OUT; commit $(git rev-parse --short HEAD)"
## one compile: task 0 compiles (rstan auto_write caches stan/ame_joint.rds); the others wait for a cache newer than the .stan file
if [ "$SLURM_ARRAY_TASK_ID" != "0" ]; then for i in $(seq 1 240); do [ stan/ame_joint.rds -nt stan/ame_joint.stan ] && break; sleep 30; done; fi
Rscript R/ame/ame_cell.R --cell=$C --bundle="$BUNDLE" --out="$OUT" --chains=4 --warmup=${WARMUP:-1000} --iter=${ITER:-1000} --threads=$TH --seed=123
