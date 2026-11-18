#!/usr/bin/env bash
## ============================================================================
## run_rel_sweep.sh -- relational-measure sweep, end to end, on a local machine (PI 2026-11-18: run on the laptop)
## Run from the repo root:  bash v2/R/run_rel_sweep.sh <netdata_dir> <cores>
##   netdata_dir : folder holding dipl_dyad.csv, destaRaw.rda, dca_dyad.csv
##                 (Dropbox cooperative_regimes/data/; dca_dyad.csv and dipl_dyad.csv are the same files used for IO-9933)
##   cores       : parallel workers (forked; macOS/Linux)
## Needs: R >= 4.3 with dplyr, tidyr, peacesciencer (+ download_extdata() once), countrycode, fixest, igraph, igoR, survival;
##        the entry data build in v2/data_entry (01_build_data.R --entry=later) already present.
## Stages (each run writes to a fresh v2/output/rel_sweep/run_<timestamp>/):
##   1 layers  2 long frames (5 type schemes)  3 positions (core W1-W4: 768 specs with all PTAs + 384 PTA-containing specs in force only;
##   dipl arm W1-W4 x 768; IGO arm 6) = 7,686 positions -> x 5 schemes x 2 year scales = 76,860 spec-runs (153,720 coxph fits, M0 + M1)
##   4 plausibility (no models)  5 sweep (every position x every scheme)
## ============================================================================
set -euo pipefail
NET="$1"; CORES="${2:-4}"
## GPT review P2-7: every run writes to a FRESH directory and the grid is an explicit file list, never a glob of a reused folder
RUN=$(date +%Y%m%d_%H%M%S)
OUT=v2/output/rel_sweep/run_$RUN; mkdir -p "$OUT/positions" "$OUT/face" "$OUT/sweep" v2/logs
LOG=v2/logs/rel_sweep_$RUN.log
echo "commit $(git rev-parse HEAD)$(git diff --quiet HEAD -- v2/R || echo '-dirty')  cores $CORES  netdata $NET  run dir $OUT" | tee "$LOG"
Rscript -e 'cat(R.version.string, "\n"); for (p in c("peacesciencer","fixest","igraph","igoR","survival","countrycode")) cat(p, as.character(packageVersion(p)), "\n")' | tee -a "$LOG"

## 1. layer panel
Rscript v2/R/20_build_layers.R "$NET" "$OUT/layers_dy.rds" 2>&1 | tee -a "$LOG"

## 2. long frames; every scheme also yields the pooled model (M0)
for T in bc bcj bcjs bcd s22; do
  Rscript v2/R/10_clogit_roles.R --comparison=directed --roles=bcj --entry=later --types=$T --rel_scale=year --boot=0 --dump_long=TRUE 2>&1 | tail -2 | tee -a "$LOG"
  mv v2/output/rel_sweep/long_${T}_relyear_fmid.rds "$OUT/"   # the dump script writes to v2/output/rel_sweep; move it into this run
done

## 3. positions (PI 2026-11-18: W1-W4, clipped residuals; PTAs both ways; dipl arm at pta=all only)
for PT in all inforce; do
  Rscript v2/R/21_rel_measures.R "$OUT/layers_dy.rds" "$OUT/positions" --arm=core --pta=$PT --cores="$CORES" 2>&1 | tee -a "$LOG"
done
Rscript v2/R/21_rel_measures.R "$OUT/layers_dy.rds" "$OUT/positions" --arm=dipl --pta=all --cores="$CORES" 2>&1 | tee -a "$LOG"
Rscript v2/R/21_rel_measures.R "$OUT/layers_dy.rds" "$OUT/positions" --arm=igo 2>&1 | tee -a "$LOG"

## the grid = exactly these files (refuse anything else in the folder)
POS=( "$OUT"/positions/positions_core_W{1raw,2contig,3dist,4fe}_ptaall.rds "$OUT"/positions/positions_core_W{1raw,2contig,3dist,4fe}_ptainforce.rds
      "$OUT"/positions/positions_dipl_W{1raw,2contig,3dist,4fe}_ptaall.rds "$OUT"/positions/positions_igo_arm.rds )
for F in "${POS[@]}"; do [ -f "$F" ] || { echo "missing expected positions file $F" | tee -a "$LOG"; exit 1; }; done
N_FOUND=$(ls "$OUT"/positions/positions_*.rds | wc -l); [ "$N_FOUND" -eq "${#POS[@]}" ] || { echo "unexpected positions files in $OUT/positions" | tee -a "$LOG"; exit 1; }
(cd "$OUT" && shasum -a 256 layers_dy.rds positions/*.rds long_*.rds > manifest_sha256.txt)

## 4. plausibility (read this before the sweep output)
Rscript v2/R/23_rel_face.R "$OUT/long_bc_relyear_fmid.rds" "$OUT/face" "${POS[@]}" 2>&1 | tee -a "$LOG"

## 5. sweep: every position x every scheme x both year-scale conventions (t as 10_clogit_roles.R; t-1)
for YS in t tm1; do
  for T in bc bcj bcjs bcd s22; do
    for F in "${POS[@]}"; do
      Rscript v2/R/22_rel_sweep.R "$OUT/long_${T}_relyear_fmid.rds" "$F" "$OUT/sweep" --yscale=$YS --cores="$CORES" 2>&1 | grep -E "CENSUS|ERROR|validation|Error" | tee -a "$LOG"
    done
  done
done
echo "done $(date)" | tee -a "$LOG"
