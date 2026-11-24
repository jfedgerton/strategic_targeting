#!/bin/bash
###############################################################################
# run_all.sh   (Jared 2026-11-24)
# Full replication pipeline for the book's current network, communities, and H1-H3 models,
# from data/raw only. Seed 123 throughout.
#
# Stages (run all by default; STAGES="1 2 3" runs a subset, in order):
#   1  data       layers_dyad_year.rds from the Tremors panel; prebuilt inputs copied in; IGO and PTA side files
#   2  measures   four builds of the JOP five-layer measure (IGO count, def, off, DCA, old JOP PTA):
#                   W  = weighted   sry + contig12, DCA carried to 2014     (communities, H1)
#                   Wt = weighted   sry + contig12, truncated at 2010       (H2 comparison column)
#                   U  = unweighted raw, DCA carried to 2014                (H3)
#                   Ut = unweighted raw, truncated at 2010                  (H2)
#   3  communities multislice Leiden on W tie_pos (resolution 1, interslice weight 1, zero-tie state-years singletons)
#   4  frames     H1 conflict frame (on W) and diffusion institutions frames (on U and W)
#   5  models     H1 (LPM + logit, with and without lagged fatal MID), H2 churn table, H3 diffusion hazard (+ no-shock/stable samples)
#   6  validate   historical-bloc benchmarks for the communities
#
# Requirements: R 4.5 with data.table, fixest, igraph, haven, countrycode, dplyr (here /opt/mm/envs/r);
#               Python 3 with igraph, leidenalg 0.12, pandas. Set RSCRIPT / PY to override.
# Outputs: data/processed/ (side files), multiverse/<tag>/ (measures, frames), results/ (tables), logs/
###############################################################################
set -euo pipefail
BOOK=$(cd "$(dirname "$0")" && pwd); export BOOK
RSCRIPT=${RSCRIPT:-/opt/mm/envs/r/bin/Rscript}; PY=${PY:-python3}
STAGES=${STAGES:-"1 2 3 4 5 6"}
MV=$BOOK/multiverse; P=$BOOK/data/processed; RES=$BOOK/results; LOG=$BOOK/logs
mkdir -p $P/variants $MV $RES $LOG $BOOK/results/main
cd $BOOK
SET=igo+def+off+dca+pta
W=c9_igo.def.off.dca.pta_sry_contig12; U=c9_igo.def.off.dca.pta_raw
COMM=$P/communities_multislice_w1.csv
has() { [[ " $STAGES " == *" $1 "* ]]; }
step() { local n=$1; shift; echo "$(date +%T) start $n"; if "$@" > $LOG/$n.log 2>&1; then echo "$(date +%T) done  $n"; else echo "FAILED $n (see logs/$n.log)"; tail -20 $LOG/$n.log; exit 1; fi; }
check() { [ -s "$1" ] || { echo "MISSING expected output $1"; exit 1; }; }

# per-tag root: symlinks to the shared inputs the chapter scripts expect
mkroot() { local root=$MV/$1; mkdir -p $root/data/processed/variants $root/results/main $root/logs
  for f in layers_dyad_year.rds state_year_covariates.rds; do ln -sf $P/$f $root/data/processed/$f; done
  ln -sf main_measure_allyigo.rds $root/data/processed/main_measure.rds; }
build() { # tag weight geo dca_mode
  local tag=$1 w=$2 geo=$3 dm=$4; mkroot $tag
  step measure_$tag env BOOK_ROOT=$BOOK OUT_ROOT=$MV/$tag LAYERS=custom IGO=dyadic WEIGHT=$w LAYER_SET=$SET GEO=$geo DCA_MODE=$dm \
       $RSCRIPT $BOOK/code/04_multiverse/01_build_measure.R
  check $MV/$tag/data/processed/main_measure_edges_allyigo.csv.gz; }

# ---- 1 data ---------------------------------------------------------------
if has 1; then
  for f in pta_layers.rds atop_un_layers.rds state_year_covariates.rds; do cp $BOOK/data/raw/prebuilt/$f $P/$f; done
  step 01_import_layers env BOOK_ROOT=$BOOK $RSCRIPT code/01_data/01_import_layers.R; check $P/layers_dyad_year.rds
  step 06_igo_variants      env BOOK=$BOOK $RSCRIPT code/01_data/06_igo_variants.R;          check $P/igo_variants.rds
  step 00_igo_hyperedge     env BOOK_ROOT=$BOOK $RSCRIPT code/04_multiverse/00_igo_hyperedge_projection.R; check $P/igo_hyper_size.rds
  step 04_pta_desta_costly  env BOOK=$BOOK $RSCRIPT code/01_data/04_pta_desta_costly.R;     check $P/pta_costly.rds
  step 04b_pta_variants     env BOOK=$BOOK $RSCRIPT code/01_data/04b_pta_variants.R;        check $P/pta_variants.rds
  step 05_igo_mia           env BOOK=$BOOK $RSCRIPT code/01_data/05_igo_mia_backfill.R;     check $P/igo_mia.rds
  step 05b_igo_mia2         env BOOK=$BOOK $RSCRIPT code/01_data/05b_igo_mia_fixed.R;       check $P/igo_mia2.rds
  step 05b_igo_mia2_cap999  env BOOK=$BOOK CAP=999 $RSCRIPT code/01_data/05b_igo_mia_fixed.R
  step 05b_igo_mia2_noreg   env BOOK=$BOOK CAP=999 NOREG=1 $RSCRIPT code/01_data/05b_igo_mia_fixed.R
fi

# ---- 2 measures -----------------------------------------------------------
if has 2; then
  build $W             sry contig12 carry
  build ${W}__trunc2010 sry contig12 trunc2010
  build $U             raw dist     carry       # GEO is unused when WEIGHT=raw
  build ${U}__trunc2010 raw dist    trunc2010
fi

# ---- 3 communities --------------------------------------------------------
if has 3; then
  step communities_multislice env GAMMA=1 OMEGA=1 EDGES_FILE=$MV/$W/data/processed/main_measure_edges_allyigo.csv.gz OUT_DIR=$MV/communities \
       $PY code/02_network_universes/02c_leiden_temporal.py
  cp $MV/communities/communities_temporal_g1_w1.csv $COMM; check $COMM
  for t in $W $U; do cp $COMM $MV/$t/data/processed/variants/communities_glm_allyigo.csv; done
fi

# ---- 4 frames -------------------------------------------------------------
if has 4; then
  step frame_conflict env ROOT=$MV/$W $RSCRIPT code/05_robustness/00_conflict_frame.R; check $MV/$W/data/processed/robust_conflict_frame.rds
  for t in $U $W; do
    step frame_institutions_$t bash -c "cd $MV/$t && BOOK_ROOT=$MV/$t MEASURE=allyigo $RSCRIPT $BOOK/code/03_chapters/03_orders_spread_institutions.R"
    check $MV/$t/data/processed/state_year_institutions_allyigo.rds
  done
fi

# ---- 5 models -------------------------------------------------------------
if has 5; then
  for L in 0 1; do
    step h1_lag$L env BOOK=$BOOK FRAME=$MV/$W/data/processed/robust_conflict_frame.rds COMM_FILE=$COMM LAG=$L OUT=$RES/h1_multislice_lag$L.csv \
         $RSCRIPT code/03_chapters/02h_h1_multislice.R
  done
  step h2_churn env MV=$MV OUT=$RES LAGALL=1 COMM_FILE=$COMM SPECS="Unweighted=$U,Weighted=$W" OTAG=multislice \
       $RSCRIPT code/03_chapters/01c_churn_tremors_rev_table.R
  for t in $U $W; do mkdir -p $RES/h3_$t
    step h3_$t env MV=$MV OUT=$RES/h3_$t TAG=$t $RSCRIPT code/03_chapters/03l_diffusion_hazard.R
  done
  for s in noshock stable both; do mkdir -p $RES/h3_${U}_$s
    step h3_${U}_$s env MV=$MV OUT=$RES/h3_${U}_$s TAG=$U SAMPLE=$s COMM_FILE=$COMM $RSCRIPT code/03_chapters/03o_diffusion_no_shocks.R
  done
fi

# ---- 6 validate -----------------------------------------------------------
if has 6; then
  mkdir -p $RES/validation
  step validate env MV=$MV TAGS="multislice_w1=$W" OUT=$RES/validation LAYERS_FILE=$P/layers_dyad_year.rds \
       $RSCRIPT code/99_checks/validate_communities.R
fi
echo "$(date +%T) ALL DONE. Tables in $RES"
