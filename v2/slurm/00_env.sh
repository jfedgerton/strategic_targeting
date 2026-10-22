#!/bin/bash
# Shared environment for v2 jobs on Roar. Sourced by every sbatch script.
# Stan on Roar needs a C++17 toolchain and the matching libstdc++ at runtime.
module purge
module load gcc/12.2.0
# Roar has no R 4.3.x module with rstan installed; jfe4's working Stan setup (used by the old pipeline)
# is r/4.2.1-gcc-8.5.0 + ~/R/x86_64-pc-linux-gnu-library/4.2 (rstan 2.32.7, StanHeaders 2.32.10). Verified 2026-10-22.
module load r/4.2.1-gcc-8.5.0
# libstdc++ fix: the gcc module's runtime library must be found before the system one
export LD_LIBRARY_PATH="$(dirname "$(gcc -print-file-name=libstdc++.so)"):${LD_LIBRARY_PATH}"
# V2_ROOT defaults to <submit dir>/v2, so jobs must be submitted from the repo root.
export V2_ROOT="${V2_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}/v2}"
if [ ! -f "$V2_ROOT/R/01_build_data.R" ]; then
  echo "ERROR: V2_ROOT=$V2_ROOT does not contain R/01_build_data.R; submit from the repo root" >&2
  exit 1
fi
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/x86_64-pc-linux-gnu-library/4.2}"
export OMP_NUM_THREADS=1
mkdir -p "$V2_ROOT/logs" "$V2_ROOT/output"
echo "host $(hostname); gcc $(gcc -dumpversion); R $(R --version | head -1); V2_ROOT=$V2_ROOT"
