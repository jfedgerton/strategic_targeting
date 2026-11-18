#!/usr/bin/env bash
## ============================================================================
## run_rel_sweep.sh -- relational-measure sweep, end to end, on a local machine (PI 2026-11-18: run on the laptop)
## Run from the repo root:  bash v2/R/run_rel_sweep.sh <netdata_dir> <cores>
##   netdata_dir : folder holding dipl_dyad.csv, destaRaw.rda, dca_dyad.csv
##                 (Dropbox cooperative_regimes/data/; dca_dyad.csv and dipl_dyad.csv are the same files used for IO-9933)
##   cores       : parallel workers (forked; macOS/Linux)
## Needs: R >= 4.3 with dplyr, tidyr, peacesciencer (+ download_extdata() once), countrycode, fixest, igraph, igoR, survival;
##        the entry data build in v2/data_entry (01_build_data.R --entry=later) already present.
## Stages (each writes to v2/output/rel_sweep/; rerunning a stage overwrites it):
##   1 layers  2 long frames (5 type schemes)  3 positions (core 5 x 768, dipl arm 5 x 768, IGO arm 4)
##   4 plausibility (no models)  5 sweep (every position x every scheme)
## ============================================================================
set -euo pipefail
NET="$1"; CORES="${2:-4}"
OUT=v2/output/rel_sweep; mkdir -p "$OUT/positions" "$OUT/face" "$OUT/sweep" v2/logs
LOG=v2/logs/rel_sweep_$(date +%Y%m%d_%H%M%S).log
echo "commit $(git rev-parse HEAD)  cores $CORES  netdata $NET" | tee "$LOG"
Rscript -e 'cat(R.version.string, "\n"); for (p in c("peacesciencer","fixest","igraph","igoR","survival","countrycode")) cat(p, as.character(packageVersion(p)), "\n")' | tee -a "$LOG"

## 1. layer panel
Rscript v2/R/20_build_layers.R "$NET" "$OUT/layers_dy.rds" 2>&1 | tee -a "$LOG"

## 2. long frames; every scheme also yields the pooled model (M0)
for T in bc bcj bcjs bcd s22; do
  Rscript v2/R/10_clogit_roles.R --comparison=directed --roles=bcj --entry=later --types=$T --rel_scale=year --boot=0 --dump_long=TRUE 2>&1 | tail -2 | tee -a "$LOG"
done

## 3. positions
Rscript v2/R/21_rel_measures.R "$OUT/layers_dy.rds" "$OUT/positions" --arm=core --cores="$CORES" 2>&1 | tee -a "$LOG"
Rscript v2/R/21_rel_measures.R "$OUT/layers_dy.rds" "$OUT/positions" --arm=dipl --cores="$CORES" 2>&1 | tee -a "$LOG"
Rscript v2/R/21_rel_measures.R "$OUT/layers_dy.rds" "$OUT/positions" --arm=igo 2>&1 | tee -a "$LOG"

## 4. plausibility (read this before the sweep output)
Rscript v2/R/23_rel_face.R "$OUT/long_bc_relyear_fmid.rds" "$OUT/face" "$OUT"/positions/positions_*.rds 2>&1 | tee -a "$LOG"

## 5. sweep
for T in bc bcj bcjs bcd s22; do
  for F in "$OUT"/positions/positions_*.rds; do
    Rscript v2/R/22_rel_sweep.R "$OUT/long_${T}_relyear_fmid.rds" "$F" "$OUT/sweep" --cores="$CORES" 2>&1 | tail -3 | tee -a "$LOG"
  done
done
echo "done $(date)" | tee -a "$LOG"
