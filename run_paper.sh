#!/usr/bin/env bash
## run_paper.sh -- the paper pipeline on the main measure (PI 2026-11-24): networks, main table, Table 1 and figures, the network-measure grid
## (306 measures x fatal MIDs and all MIDs), bootstrapped grid cells, grid summaries and the AME data bundle. Run after run_all.sh (TYPES=init)
## has finished in the same output folder (it supplies layers_dy.rds, the stage-1 dumps and a positions key).
## Steps (each resumable; a finished step is skipped on rerun):
##   1. networks: book PTA/DCA layers; MAIN measure (unweighted book measure, BOOKraw); weighted book measure (BOOKc12, appendix);
##      grid layers and grid positions (306 measures, R/lib/grid.R)
##   2. stage-2 builds, every window x {fmid, mid}, on the main measure, each accepted against its stage-1 init dump (in $OUT/paper/builds,
##      separate from run_stage2.sh builds so the two runners never overwrite each other; audit 2026-11-24)
##   3. main table: two-stage bootstrap of the main cell (B_MAIN draws), fmid and mid at the main window; stage 1 / stage 2 / joint table;
##      marginal-effect plots with 95% bootstrap bands (B_CELL draws, fmid)
##   4. Table 1 and draft figures (main measure, main window)
##   5. grid sweep: 306 measures x main window x {fmid, mid} (joint two-stage FE point fits)
##   6. bootstrapped grid cells (R/grid_pick_boot_cells.R; B_CELL; the main-window main cells come from step 3)
##   7. grid summaries (distributions, decision regressions, specification curves)
##   8. AME data bundle for the joint AME (main measure, main window, fmid); the fits run on Roar (slurm/ame_cells.sbatch)
## Usage: OUT=output_init ./run_paper.sh            (optional: P_STEPS="1 2 3 4 5 6 7 8", B_MAIN=2000, B_CELL=500)
set -euo pipefail
cd "$(dirname "$0")"; source config.sh
O="${OUT:-output_init}"; S2=$O/stage2; PP=$O/paper; mkdir -p $S2/{edges,logs} $PP/{builds,main_table,marginal_effects,table1,sweep,boot,summary,ame}
STEPS="${P_STEPS:-1 2 3 4 5 6 7 8}"; B_MAIN="${B_MAIN:-2000}"; B_CELL="${B_CELL:-500}"
LOG=$PP/run_$(date +%Y%m%d_%H%M%S).log
say() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }
die() { say "STOP: $*"; exit 1; }
runR() { local PAT=$1; shift; local T; T=$(mktemp); local rc=0; Rscript "$@" > "$T" 2>&1 || rc=$?; cat "$T" >> "$LOG"; grep -E "$PAT" "$T" | cut -c1-250 || true; rm -f "$T"; return $rc; }
has() { [[ " $STEPS " == *" $1 "* ]]; }
h() { shasum -a 256 "$1" | cut -c1-16; }
fresh() { [ -s "$1.stamp" ] && [ "$(cat "$1.stamp")" = "$2" ]; }
stamp() { echo "$2" > "$1.stamp.tmp" && mv "$1.stamp.tmp" "$1.stamp"; }
wtag() { [ "$1" -gt 0 ] && echo "_w${1}d" || echo ""; }
GAPS=""; gap() { local rc=$1 NM=$2; [ $rc -eq 0 ] && return 0; [ $rc -eq 3 ] && { say "WARNING: $NM finished with ERROR cells (rerun to retry)"; GAPS="$GAPS $NM"; return 0; }; die "$NM failed (exit $rc)"; }
[ -e $O/.created_by_run_all ] || die "$O is not a run_all.sh output folder"
[ "${S2_FORCE:-0}" = "1" ] || ! pgrep -f "run_all.sh|R/07_fit_sweep.R" > /dev/null || die "run_all.sh looks still running (S2_FORCE=1 overrides)"
[ -s $O/layers_dy.rds ] || die "missing $O/layers_dy.rds"
[ -e data_entry ] || die "data_entry link missing (made by run_all.sh)"
[ -s $NETDATA/dca_dyad.csv ] || die "missing $NETDATA/dca_dyad.csv"
POSD=$O/positions; KEYF=$POSD/positions_core_W2contig_ptaall.rds; [ -s $KEYF ] || die "missing positions key $KEYF (from run_all.sh)"
MAINF=$POSD/positions_core_${MAIN_SPEC%%|*}.rds
GMAIN='Graw_count_book|igo+def+off+dca+pta|pr'   # the grid's copy of the main measure (equal to MAIN_SPEC, checked by R/grid_build.R)
say "paper run | commit $(git rev-parse --short HEAD) | out $PP | main $MAIN_SPEC | cores $CORES | windows $WINDOWS | main window $WINDOW_MAIN | steps $STEPS | B $B_MAIN/$B_CELL"

## ---- 1. networks ----
if has 1; then
  BL=$S2/edges/book_layers.rds
  S="$(h R/book_layers.R) $(h $O/layers_dy.rds) $(h data_book/desta_list_of_treaties_02_02.txt) $(h data_book/desta23_dyadic_withdrawal.csv) $(h $NETDATA/dca_dyad.csv)"
  { [ -s $BL ] && fresh $BL "$S"; } || { runR "BOOK_LAYERS_DONE|linked|DCA|Error" R/book_layers.R $O/layers_dy.rds data_book/desta_list_of_treaties_02_02.txt data_book/desta23_dyadic_withdrawal.csv $NETDATA/dca_dyad.csv $BL || die "book layers failed"; stamp $BL "$S"; }
  BRP=$POSD/positions_core_BOOKraw_ptadesta.rds; BRE=$S2/edges/edges_book_raw.rds; S="$(h R/book_measure_raw.R) $(h $O/layers_dy.rds) $(h $KEYF) $(h $BL)"
  { [ -s $BRP ] && [ -s $BRE ] && fresh $BRP "$S"; } || { runR "BOOK_RAW_DONE|pooled|Error" R/book_measure_raw.R $O/layers_dy.rds $BL $KEYF $BRP $BRE || die "main measure failed"; stamp $BRP "$S"; }
  BKP=$POSD/positions_core_BOOKc12_ptadesta.rds; BKE=$S2/edges/edges_book_c12.rds; S="$(h R/book_measure_c12.R) $(h $O/layers_dy.rds) $(h $KEYF) $(h $BL)"
  { [ -s $BKP ] && [ -s $BKE ] && fresh $BKP "$S"; } || { runR "BOOK_MEASURE_DONE|PTA layer|Error" R/book_measure_c12.R $O/layers_dy.rds $KEYF $BKP $BKE $BL || die "weighted book measure failed"; stamp $BKP "$S"; }
  GL=$S2/edges/grid_layers.rds; S="$(h R/grid_layers.R) $(h $O/layers_dy.rds) $(h $BL) $(cat data_book/MD5SUMS | shasum -a 256 | cut -c1-16)"
  { [ -s $GL ] && fresh $GL "$S"; } || { runR "GRID_LAYERS_DONE|IGO count|PTA depth|Error" R/grid_layers.R $O/layers_dy.rds $BL data_book $GL --cores=$CORES || die "grid layers failed"; stamp $GL "$S"; }
  GS=$POSD/grid_specs.csv; S="$(h R/grid_build.R) $(h R/lib/grid.R) $(h $GL) $(h $KEYF) $(h $BRP)"
  { [ -s $GS ] && fresh $GS "$S"; } || { runR "GRID_BUILD_DONE|specs in|vs BOOKraw|constant|Error" R/grid_build.R $GL $KEYF $BRP $POSD --cores=$CORES || die "grid positions failed"; stamp $GS "$S"; }
  rm -f $POSD/positions_core_BOOKc12_ptaany.rds
  say "1 networks done"
fi
[ -s $MAINF ] || die "missing $MAINF (step 1)"
GRIDPOS=""; if has 5 || has 6; then for B in raw dca1 ally sryc12 srydist syc12 sydist; do F=$POSD/positions_core_G$B.rds; [ -s $F ] || die "missing grid positions $F (step 1)"; GRIDPOS="$GRIDPOS${GRIDPOS:+,}$F"; done   # exactly the 7 blocks (audit 2026-11-24)
  [ "$(($(wc -l < $POSD/grid_specs.csv) - 1))" -eq 306 ] || die "grid_specs.csv does not list 306 measures"; fi

## ---- 2. builds ----
if has 2; then
  for W in $WINDOWS; do for SM in fmid mid; do
    BD=$PP/builds/w${W}_${SM}; DP=$O/dumps/long_init_bon$(wtag $W)_relyear_${SM}.rds; [ -s $DP ] || die "missing stage-1 dump $DP"
    BST="$(h R/s2_01_build_riskset.R) $(h R/lib/s2.R) $(h R/lib/s2_cell.R) $(h R/lib/sry.R) $(h $S2/edges/edges_book_raw.rds) $(h $DP) $(h $MAINF) $(h $O/layers_dy.rds) $(h data_entry/manifest.csv) $MAIN_SPEC $W $SM"
    [ -s $BD/s2_meta.rds ] && grep -q "written (release build)" $BD/s2_build.log 2>/dev/null && fresh $BD/.build "$BST" && continue
    rm -rf $BD; FD=$([ $SM = fmid ] && echo stage1 || echo mid)
    runR "acceptance|FAILED|written|Error" R/s2_01_build_riskset.R --data_entry=data_entry --layers=$O/layers_dy.rds --positions=$MAINF --poscol="$MAIN_SPEC" \
      --out=$BD --window=$W --fatal_def=$FD --dump_bcjs=$DP --edges_dir=$S2/edges || die "build w$W $SM failed"
    stamp $BD/.build "$BST"
  done; done
  say "2 builds done"
fi

## ---- 3. main table (main measure, main window, fmid and mid) ----
if has 3; then
  MC=$PP/main_table/cells.csv
  printf 'file,spec_id,window,sample,rule,n_boot,cell_id\n%s,%s,%s,fmid,main,%s,m_fmid\n%s,%s,%s,mid,main,%s,m_mid\n' "$(basename $MAINF)" "$MAIN_SPEC" $WINDOW_MAIN $B_MAIN "$(basename $MAINF)" "$MAIN_SPEC" $WINDOW_MAIN $B_MAIN > $MC
  runR "draws ok|census|INVALID|SKIPPED|S2_BOOT_DONE|S2_BOOT_GAPS|Error" R/s2_07_boot_cells.R --cells=$MC --posdir=$POSD --builds_root=$PP/builds --dumps_dir=$O/dumps \
    --edges_dir=$S2/edges --out=$PP/main_table/boot --cores=$CORES && rc=0 || rc=$?; gap $rc main_boot
  for SM in fmid mid; do
    runR "MAIN_TABLE_DONE|Error" R/main_table.R --build=$PP/builds/w${WINDOW_MAIN}_${SM} --dump=$O/dumps/long_init_bon$(wtag $WINDOW_MAIN)_relyear_${SM}.rds --positions=$MAINF \
      --spec="$MAIN_SPEC" --edges_dir=$S2/edges --cell_dir=$PP/main_table/boot/m_$SM --out=$PP/main_table/$SM --cores=$CORES || die "main table $SM failed"
  done
  runR "MARGINAL_EFFECTS_DONE|bootstrap:|Error" R/marginal_effects.R --build=$PP/builds/w${WINDOW_MAIN}_fmid --dump=$O/dumps/long_init_bon$(wtag $WINDOW_MAIN)_relyear_fmid.rds \
    --positions=$MAINF --spec="$MAIN_SPEC" --edges_dir=$S2/edges --n_boot=$B_CELL --out=$PP/marginal_effects --cores=$CORES || die "marginal effects failed"
  say "3 main table and marginal effects done"
fi

## ---- 4. Table 1 and figures ----
if has 4; then
  runR "TABLE1_FIGURES_DONE|Error" R/table1_figures.R --positions=$MAINF --spec="$MAIN_SPEC" --edges_dir=$S2/edges \
    --build_fmid=$PP/builds/w${WINDOW_MAIN}_fmid --dump_fmid=$O/dumps/long_init_bon$(wtag $WINDOW_MAIN)_relyear_fmid.rds \
    --build_mid=$PP/builds/w${WINDOW_MAIN}_mid --dump_mid=$O/dumps/long_init_bon$(wtag $WINDOW_MAIN)_relyear_mid.rds --out=$PP/table1 || die "Table 1 failed"
  say "4 Table 1 and figures done"
fi

## ---- 5. grid sweep ----
if has 5; then
  for SM in fmid mid; do
    runR "chunk|measures|census|S2_SWEEP_DONE|S2_SWEEP_GAPS|Error" R/s2_05_sweep.R --build=$PP/builds/w${WINDOW_MAIN}_${SM} --dump=$O/dumps/long_init_bon$(wtag $WINDOW_MAIN)_relyear_${SM}.rds \
      --edges_dir=$S2/edges --out=$PP/sweep/w${WINDOW_MAIN}_${SM} --positions=$GRIDPOS --cores=$CORES && rc=0 || rc=$?; gap $rc sweep_$SM
  done
  say "5 grid sweep done"
fi

## ---- 6. bootstrapped grid cells ----
if has 6; then
  CST="$(h R/grid_pick_boot_cells.R) $(h $POSD/grid_specs.csv) $MAIN_SPEC $WINDOW_MAIN $WINDOWS $B_CELL"   # regenerate the list when its inputs change (audit 2026-11-24)
  { [ -s $PP/boot/cells.csv ] && fresh $PP/boot/cells.csv "$CST"; } || { runR "wrote|Error" R/grid_pick_boot_cells.R --posdir=$POSD --main_spec="$MAIN_SPEC" --main_window=$WINDOW_MAIN --windows="$WINDOWS" \
    --b_main=$B_MAIN --b_cell=$B_CELL --out=$PP/boot/cells.csv || die "cell list failed"; stamp $PP/boot/cells.csv "$CST"; }
  runR "draws ok|census|INVALID|SKIPPED|S2_BOOT_DONE|S2_BOOT_GAPS|Error" R/s2_07_boot_cells.R --cells=$PP/boot/cells.csv --posdir=$POSD --builds_root=$PP/builds --dumps_dir=$O/dumps \
    --edges_dir=$S2/edges --out=$PP/boot --cores=$CORES && rc=0 || rc=$?; gap $rc grid_boot
  say "6 bootstrapped grid cells done"
fi

## ---- 7. grid summaries ----
if has 7; then
  runR "GRID_SUMMARY_DONE|Error" R/grid_summary.R --sweeps=$PP/sweep/w${WINDOW_MAIN}_fmid,$PP/sweep/w${WINDOW_MAIN}_mid --main_spec="$GMAIN" --out=$PP/summary || die "grid summary failed"
  say "7 grid summaries done"
fi

## ---- 8. AME data bundle ----
if has 8; then
  AB=$PP/ame/ame_bundle_w${WINDOW_MAIN}_fmid_main.rds
  runR "bundle|Error" R/ame/ame_cell.R --prep_only=1 --bundle=$AB --build=$PP/builds/w${WINDOW_MAIN}_fmid --dump=$O/dumps/long_init_bon$(wtag $WINDOW_MAIN)_relyear_fmid.rds \
    --positions=$MAINF --spec="$MAIN_SPEC" --edges_dir=$S2/edges --out=$PP/ame || die "AME bundle failed"
  say "8 AME bundle: $AB (copy to Roar; BUNDLE=... sbatch slurm/ame_cells.sbatch)"
fi
[ -z "$GAPS" ] && say "PAPER_DONE" || { say "PAPER_DONE_WITH_GAPS:$GAPS (rerun the same command to retry ERROR cells)"; exit 3; }
