#!/bin/bash
# Shared environment for v2 jobs on Roar. Sourced by every sbatch script.
# Stan on Roar needs a C++17 toolchain and the matching libstdc++ at runtime.
module purge
module load gcc/12.2.0
module load r/4.3.3 2>/dev/null || module load r
# libstdc++ fix: the gcc module's runtime library must be found before the system one
export LD_LIBRARY_PATH="$(dirname "$(gcc -print-file-name=libstdc++.so)"):${LD_LIBRARY_PATH}"
export V2_ROOT="${V2_ROOT:-/storage/group/LiberalArts/default/jfe4_collab/strategic_targeting/v2}"
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/x86_64-pc-linux-gnu-library/4.3}"
export OMP_NUM_THREADS=1
mkdir -p "$V2_ROOT/logs" "$V2_ROOT/output"
echo "host $(hostname); gcc $(gcc -dumpversion); R $(R --version | head -1); V2_ROOT=$V2_ROOT"
