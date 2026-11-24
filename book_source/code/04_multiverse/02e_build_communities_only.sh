#!/bin/bash
# Measure + communities only (no outcome models), for validation-first builds (Jared 2026-11-23).
# usage: 02e_build_communities_only.sh LAYER_SET WEIGHT [GEO]
set -u
R=/opt/renv/bin/Rscript; BOOK=${BOOK:?}; MV=${MV:?}
S=$1; W=$2; GEO=${3:-dist}; GS=""; [ "$GEO" != "dist" ] && GS="_$GEO"
if [[ "+$S+" == *"+igo+"* ]]; then I=${IGO_MODE:-dyadic}; else I=none; fi; IS=""; [ "$I" = "hyper" ] && IS="hyper"; [ "$I" = "mia" ] && IS="mia"; [ "$I" = "mia2" ] && IS="mia2"; [ "$I" = "mia2nd" ] && IS="mia2nd"; [ "$I" = "var" ] && IS="iv${IGO_VAR:?}"; [ -n "${IGO_WT:-}" ] && [ "$IGO_WT" != "1" ] && IS="${IS}w${IGO_WT/./}"; [ -n "${MIA2_CAP:-}" ] && [ "$MIA2_CAP" != "10" ] && IS="${IS}cap${MIA2_CAP}"; [ "${MIA2_NOREG:-0}" = "1" ] && IS="${IS}noreg"
PV=${PTA_VARIANT:-}; PS=""; [ -n "$PV" ] && PS="_pv$PV"
TAG=c9_${S//+/.}_${W}${IS}${GS}${PS}; ROOT=$MV/$TAG
mkdir -p $ROOT/data/processed/variants $ROOT/results/main $ROOT/logs
for f in layers_dyad_year.rds state_year_covariates.rds; do ln -sf $BOOK/data/processed/$f $ROOT/data/processed/$f; done
ln -sf main_measure_allyigo.rds $ROOT/data/processed/main_measure.rds
env BOOK_ROOT=$BOOK OUT_ROOT=$ROOT LAYERS=custom IGO=$I WEIGHT=$W LAYER_SET=$S GEO=$GEO $R $BOOK/code/04_multiverse/01_build_measure.R > $ROOT/logs/measure.log 2>&1 || { echo "measure failed $TAG"; tail -5 $ROOT/logs/measure.log; exit 1; }
env NEG_NULL=none EDGES_FILE=$ROOT/data/processed/main_measure_edges_allyigo.csv.gz COMM_FILE=$ROOT/data/processed/variants/communities_glm_allyigo.csv SUMM_FILE=$ROOT/results/main/communities_summary.csv python3 $BOOK/code/02_network_universes/02_leiden_communities.py > $ROOT/logs/leiden.log 2>&1 || { echo "leiden failed $TAG"; exit 1; }
echo "built $TAG"
