#!/usr/bin/env bash
## run_stage2.sh -- the stage-2 appendix run (PI 2026-11-22), after the stage-1 run (run_all.sh, TYPES=init) has FINISHED in the same output folder.
## Do not start (or git pull) while run_all.sh is running: bash reads its script as it goes.
## Steps (each resumable; a finished step is skipped on rerun):
##   1. per-layer edge files for every arm (R/02_build_positions.R --edges_out), for measure-specific ties (option A)
##   2. stage-2 builds, every window x {fmid, mid}, each accepted against its stage-1 init dump
##   3. layer 1a: parallel multiverse, ALL measures x main window x fmid (joint two-stage FE point fits)
##   4. layer 2:  two-stage bootstrap of the rule-based cell list (main cell B_MAIN draws, others B_CELL)
##   5. layer 1b: parallel multiverse, PLAUSIBLE measures x main window x mid
##   6. reference-model robustness checks on the main build (R/s2_02_fit_reference.R, all checks)
## Usage: OUT=output_init ./run_stage2.sh            (optional: S2_STEPS="1 2 3 4 5 6", B_MAIN=2000, B_CELL=500)
set -euo pipefail
cd "$(dirname "$0")"; source config.sh
O="${OUT:-output_init}"; S2=$O/stage2; mkdir -p $S2/{edges,builds,sweep,boot,reference,logs}
STEPS="${S2_STEPS:-1 2 3 4 5 6}"; B_MAIN="${B_MAIN:-2000}"; B_CELL="${B_CELL:-500}"
LOG=$S2/logs/run_$(date +%Y%m%d_%H%M%S).log
say() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }
die() { say "STOP: $*"; exit 1; }
runR() { local PAT=$1; shift; local T; T=$(mktemp); local rc=0; Rscript "$@" > "$T" 2>&1 || rc=$?; cat "$T" >> "$LOG"; grep -E "$PAT" "$T" | cut -c1-250 || true; rm -f "$T"; return $rc; }
has() { [[ " $STEPS " == *" $1 "* ]]; }
wtag() { [ "$1" -gt 0 ] && echo "_w${1}d" || echo ""; }
[ -e $O/.created_by_run_all ] || die "$O is not a run_all.sh output folder"
## audit P2: pgrep also matches a parent shell that chains the two scripts; S2_FORCE=1 skips this check when you know stage 1 is done
[ "${S2_FORCE:-0}" = "1" ] || ! pgrep -f "run_all.sh|R/07_fit_sweep.R" > /dev/null || die "run_all.sh (or a stage-1 fit) looks still running; start stage 2 after it finishes (S2_FORCE=1 overrides)"
[ -s $O/layers_dy.rds ] || die "missing $O/layers_dy.rds"
[ -e data_entry ] || die "data_entry link missing (made by run_all.sh)"
POSD=$O/positions; MAINF=$POSD/positions_core_W2contig_ptaall.rds
[ "$TYPES" = "init" ] || say "note: config TYPES=$TYPES; stage 2 always uses the init dumps"
say "stage-2 appendix run | commit $(git rev-parse --short HEAD) | out $S2 | cores $CORES | windows $WINDOWS | main window $WINDOW_MAIN | steps $STEPS | B $B_MAIN/$B_CELL"

## ---- 1. edges (mirrors the arms of run_all.sh) ----
if has 1; then
  ed() { local NM=$1 FILES=$2; shift 2; local ok=1 f; for f in $FILES; do [ -s $S2/edges/$f.rds ] || ok=0; done; [ $ok = 1 ] && return 0
    runR "Error|edges" R/02_build_positions.R $O/layers_dy.rds $S2/edges "$@" --edges_out=$S2/edges || die "edges $NM failed"
    for f in $FILES; do [ -s $S2/edges/$f.rds ] || die "edges $NM did not write $f.rds"; done; }
  ed core_all "edges_core_W1raw_ptaall edges_core_W2contig_ptaall edges_core_W3dist_ptaall edges_core_W4fe_ptaall" --arm=core --pta=all --cores=$CORES
  ed core_inforce "edges_core_W1raw_ptainforce edges_core_W2contig_ptainforce edges_core_W3dist_ptainforce edges_core_W4fe_ptainforce" --arm=core --pta=inforce --cores=$CORES
  ed core_any "edges_core_W1raw_ptaany edges_core_W2contig_ptaany edges_core_W3dist_ptaany edges_core_W4fe_ptaany" --arm=core --pta=any --cores=$CORES
  for P in all inforce any; do
    ed core_sry_$P "edges_core_W6srybook_pta$P edges_core_W7sryyear_pta$P edges_core_W6srybookmm_pta$P edges_core_W7sryyearmm_pta$P" --arm=core --pta=$P --weightings=W6srybook,W7sryyear,W6srybookmm,W7sryyearmm --cores=$CORES
  done
  ed dipl "edges_dipl_W1raw_ptaall edges_dipl_W2contig_ptaall edges_dipl_W3dist_ptaall edges_dipl_W4fe_ptaall" --arm=dipl --pta=all --cores=$CORES
  ed igo "edges_igo_arm" --arm=igo
  say "1 edges done: $(ls $S2/edges | wc -l | tr -d ' ') files"
fi

## ---- 2. builds ----
if has 2; then
  for W in $WINDOWS; do for SM in fmid mid; do
    BD=$S2/builds/w${W}_${SM}; [ -s $BD/s2_meta.rds ] && grep -q "written (release build)" $BD/s2_build.log 2>/dev/null && continue
    DP=$O/dumps/long_init_bon$(wtag $W)_relyear_${SM}.rds; [ -s $DP ] || die "missing stage-1 dump $DP"
    FD=$([ $SM = fmid ] && echo stage1 || echo mid)
    runR "acceptance|FAILED|written|Error" R/s2_01_build_riskset.R --data_entry=data_entry --layers=$O/layers_dy.rds --positions=$MAINF --poscol="$MAIN_SPEC" \
      --out=$BD --window=$W --fatal_def=$FD --dump_bcjs=$DP || die "build w$W $SM failed"
  done; done
  say "2 builds done"
fi

ALLPOS=$(ls "$POSD"/positions_*.rds 2>/dev/null | grep -v positions_cmp | paste -sd, - || true)
{ has 3 || has 5; } && [ -z "$ALLPOS" ] && die "no positions files in $POSD"
MAINB=$S2/builds/w${WINDOW_MAIN}_fmid
## ---- 3. layer 1a: all measures x main window x fmid ----
if has 3; then
  runR "chunk|measures|S2_SWEEP_DONE|Error" R/s2_05_sweep.R --build=$MAINB --dump=$O/dumps/long_init_bon$(wtag $WINDOW_MAIN)_relyear_fmid.rds \
    --edges_dir=$S2/edges --out=$S2/sweep/w${WINDOW_MAIN}_fmid_all --positions=$ALLPOS --cores=$CORES || die "sweep (all, fmid) failed"
  say "3 layer 1a done"
fi

## ---- 4. layer 2: bootstrapped cells ----
if has 4; then
  [ -s $S2/boot/cells.csv ] || runR "wrote|Error" R/s2_06_pick_boot_cells.R --posdir=$POSD --main_spec="$MAIN_SPEC" --main_window=$WINDOW_MAIN --windows="$WINDOWS" \
    --b_main=$B_MAIN --b_cell=$B_CELL --out=$S2/boot/cells.csv || die "cell list failed"
  runR "draws ok|S2_BOOT_DONE|Error" R/s2_07_boot_cells.R --cells=$S2/boot/cells.csv --posdir=$POSD --builds_root=$S2/builds --dumps_dir=$O/dumps \
    --edges_dir=$S2/edges --out=$S2/boot --cores=$CORES || die "bootstrap cells failed"
  say "4 layer 2 done"
fi

## ---- 5. layer 1b: plausible measures x main window x mid ----
if has 5; then
  runR "chunk|measures|S2_SWEEP_DONE|Error" R/s2_05_sweep.R --build=$S2/builds/w${WINDOW_MAIN}_mid --dump=$O/dumps/long_init_bon$(wtag $WINDOW_MAIN)_relyear_mid.rds \
    --edges_dir=$S2/edges --out=$S2/sweep/w${WINDOW_MAIN}_mid_plausible --positions=$ALLPOS --validation=$O/validation/measure_validation.csv --set=plausible --cores=$CORES || die "sweep (plausible, mid) failed"
  say "5 layer 1b done"
fi

## ---- 6. reference-model robustness (main build) ----
if has 6; then
  [ -e $S2/reference/.done ] || { runR "S2_FIT_DONE|FAILED|Error" R/s2_02_fit_reference.R --build=$MAINB --out=$S2/reference --n_perm=1000 --n_null=1000 --n_boot=$B_MAIN --n_noise=20 --cores=$CORES --checks=all || die "reference fitter failed"; touch $S2/reference/.done; }
  say "6 reference robustness done"
fi
say "STAGE2_DONE"
