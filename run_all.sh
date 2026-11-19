#!/usr/bin/env bash
## ============================================================================
## run_all.sh -- the full multiverse loop (2026-11-19). Run from the repo root:   bash run_all.sh
## Every step writes into ./output and is SKIPPED when its output already exists, so a stopped run can be restarted with the same command.
## Steps (see README.md for what each one answers):
##   1 positions       layers and all 7,686 relational positions (or reuse, config.sh)
##   2 validation      outcome-blind plausibility + power-capture metrics for every measure
##   3 measure sets    Phase A (main + 200), no-ATOP, main only, 25 random plausible
##   4 dumps           event frames: three types x windows 0-6 x {fatal, all MID onsets}; four types x windows; unit weights
##   5 multiverse      every measure x every window (fatal; all-MID step A); step B (failed bargains) on every measure
##   6 robustness      every robustness check on the 201 Phase A measures at the main window
##   7 influence       leave-one-dispute-out and leave-one-initiator-out on the main measure, then drop the top 5 / top 10 disputes
##   8 resampling      pigeonhole bootstrap (main measure at every window; 25 random measures), permutation tests
##   9 event study     relational trajectory of chosen targets around the dispute
##  10 summary         appendix tables
## ============================================================================
set -euo pipefail
cd "$(dirname "$0")"; source config.sh
## TEST=1 bash run_all.sh : mechanics check only (windows 0 and 1, the Phase A measures as the "multiverse", tiny B); never use for results
if [ "${TEST:-0}" = "1" ]; then WINDOWS="0 1"; B_BOOT=20; B_PERM=20; N_RANDOM_BOOT=2; MERROR_K="0.5"; fi
O=output; mkdir -p $O/{positions,validation,sets,dumps,sweep,resample,influence,eventstudy,summary,logs}
LOG=$O/logs/run_$(date +%Y%m%d_%H%M%S).log
say() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }
done_or() { [ -s "$1" ] && { say "skip (exists): $1"; return 0; } || return 1; }   # done_or <output> || <command>
say "commit $(git rev-parse HEAD 2>/dev/null)$(git diff --quiet HEAD 2>/dev/null || echo '-dirty') | cores $CORES | types $TYPES | main window $WINDOW_MAIN"
Rscript -e 'cat(R.version.string, "\n"); for (p in c("survival","data.table","dplyr","fixest","igraph","peacesciencer","countrycode","sandwich","ggplot2")) cat(p, as.character(packageVersion(p)), "\n")' | tee -a "$LOG"

## ---- 0. inputs ----
[ -e data_entry ] || ln -s "$DATA_ENTRY" data_entry
for f in events.rds participants_timing.rds sy_positions.rds manifest.csv; do [ -s data_entry/$f ] || { say "missing data_entry/$f"; exit 1; }; done

## ---- 1. positions ----
if [ "$REUSE_POSITIONS" = "1" ]; then
  for F in "$REUSE_CORE_IGO"/positions_core_*.rds "$REUSE_CORE_IGO"/positions_igo_arm.rds "$REUSE_DIPL"/positions_dipl_*.rds; do
    [ -e $O/positions/$(basename "$F") ] || ln -s "$F" $O/positions/; done
  done_or $O/layers_dy.rds || Rscript R/01_build_layers.R "$NETDATA" $O/layers_dy.rds 2>&1 | tee -a "$LOG"
else
  done_or $O/layers_dy.rds || Rscript R/01_build_layers.R "$NETDATA" $O/layers_dy.rds 2>&1 | tee -a "$LOG"
  for PT in all inforce; do ls $O/positions/positions_core_*_pta$PT.rds >/dev/null 2>&1 || Rscript R/02_build_positions.R $O/layers_dy.rds $O/positions --arm=core --pta=$PT --cores=$CORES 2>&1 | tee -a "$LOG"; done
  ls $O/positions/positions_dipl_*.rds >/dev/null 2>&1 || Rscript R/02_build_positions.R $O/layers_dy.rds $O/positions --arm=dipl --pta=all --cores=$CORES 2>&1 | tee -a "$LOG"
  done_or $O/positions/positions_igo_arm.rds || Rscript R/02_build_positions.R $O/layers_dy.rds $O/positions --arm=igo 2>&1 | tee -a "$LOG"
fi
POS=( $O/positions/positions_*.rds ); [ ${#POS[@]} -eq 13 ] || { say "expected 13 positions files, found ${#POS[@]}"; exit 1; }
(cd $O && shasum -a 256 positions/*.rds layers_dy.rds > positions_sha256.txt) 2>/dev/null || true

## ---- 2. validation ----
echo "$MAIN_SPEC" > $O/sets/main_id.txt
VPOS=( "${POS[@]}" ); [ "${TEST:-0}" = "1" ] && VPOS=( $O/positions/positions_core_W2contig_ptaall.rds )
done_or $O/validation/measure_validation.csv || Rscript R/03_validate_measures.R $O/validation "${VPOS[@]}" --top=$O/sets/main_id.txt --cores=$CORES 2>&1 | tee -a "$LOG"

## ---- 3. measure sets ----
done_or $O/sets/phaseA.rds || Rscript R/04_subset_positions.R $O/positions $O/validation/measure_validation.csv $O/sets/phaseA.rds --set=ids --ids=inputs/phaseA_specs.txt 2>&1 | tee -a "$LOG"
done_or $O/sets/main.rds || Rscript R/04_subset_positions.R $O/positions $O/validation/measure_validation.csv $O/sets/main.rds --set=ids --ids=$O/sets/main_id.txt 2>&1 | tee -a "$LOG"
done_or $O/sets/noatop.rds || Rscript R/04_subset_positions.R $O/positions $O/validation/measure_validation.csv $O/sets/noatop.rds --set=noatop 2>&1 | tee -a "$LOG"
[ -s $O/sets/random_ids.txt ] || Rscript -e 'set.seed(123); v <- read.csv(commandArgs(TRUE)[1]); p <- setdiff(v$spec_id[v$plausible], commandArgs(TRUE)[3]); writeLines(sort(sample(p, as.integer(commandArgs(TRUE)[4]))), commandArgs(TRUE)[2])' \
  $O/validation/measure_validation.csv $O/sets/random_ids.txt "$MAIN_SPEC" $N_RANDOM_BOOT
done_or $O/sets/random.rds || Rscript R/04_subset_positions.R $O/positions $O/validation/measure_validation.csv $O/sets/random.rds --set=ids --ids=$O/sets/random_ids.txt 2>&1 | tee -a "$LOG"

## ---- 4. event frames ----
dump() {  # dump <types> <window> <sample> [extra flags]   -> output/dumps/long_<types>_bon[_w<W>d][_uw]_relyear_<sample>.rds
  local T=$1 W=$2 SM=$3; shift 3; local WF=""; local WT=""; [ "$W" -gt 0 ] && { WF="--coorig_months=$W --coorig_bound=definite"; WT="_w${W}d"; }
  local UW=""; [[ " $* " == *"--unit_weights=TRUE"* ]] && UW="_uw"
  local F=$O/dumps/long_${T}_bon${WT}_relyear${UW}_${SM}.rds
  done_or $F && return 0
  Rscript R/05_dump_events.R --comparison=directed --roles=bcj --entry=later --types=$T --rel_scale=year --bilateral=onset $WF --sample=$SM --boot=0 --dump_long=TRUE "$@" 2>&1 \
    | grep -E "P1-2|check:|promotions undone|types $T:|fatal-MID|dump_long: wrote|Error" | tee -a "$LOG"
  mv output/rel_sweep/$(basename $F) $F
}
for W in $WINDOWS; do dump $TYPES $W fmid; dump $TYPES $W mid; dump bcjs $W fmid; done
dump $TYPES $WINDOW_MAIN fmid --unit_weights=TRUE
wtag() { [ "$1" -gt 0 ] && echo "_w${1}d" || echo ""; }
MAIN_F=$O/dumps/long_${TYPES}_bon$(wtag $WINDOW_MAIN)_relyear_fmid.rds
MAIN_MID=$O/dumps/long_${TYPES}_bon$(wtag $WINDOW_MAIN)_relyear_mid.rds

## ---- sweep helper: one call = one dump x one positions file; skipped when its output exists ----
sweep() {  # sweep <outdir> <dump> <positions> [options]
  local OD=$1 DP=$2 PS=$3; shift 3; mkdir -p $OD
  local MARK=$OD/.done_$(echo "$(basename $DP) $(basename $PS) $*" | shasum | cut -c1-16)
  [ -e $MARK ] && return 0
  Rscript R/07_fit_sweep.R $DP $PS $OD --bilateral=onset --cores=$CORES "$@" 2>&1 | grep -E "CENSUS|ERROR|Error|dropped|restricted|weight kept" | tee -a "$LOG"
  touch $MARK
}

## ---- 5. multiverse: every measure x every window ----
MV=( "${POS[@]}" ); [ "${TEST:-0}" = "1" ] && MV=( $O/sets/main.rds )
for W in $WINDOWS; do
  for F in "${MV[@]}"; do
    sweep $O/sweep/multiverse_fmid_w$W $O/dumps/long_${TYPES}_bon$(wtag $W)_relyear_fmid.rds $F
    sweep $O/sweep/multiverse_mid_w$W  $O/dumps/long_${TYPES}_bon$(wtag $W)_relyear_mid.rds $F
  done
done
for F in "${MV[@]}"; do sweep $O/sweep/multiverse_stepB_w$WINDOW_MAIN $MAIN_MID $F --stepB=TRUE; done
PA=$O/sets/phaseA.rds; [ "${TEST:-0}" = "1" ] && PA=$O/sets/main.rds
for W in $WINDOWS; do sweep $O/sweep/fourtypes_w$W $O/dumps/long_bcjs_bon$(wtag $W)_relyear_fmid.rds $PA; done

## ---- 6. robustness on the 201 Phase A measures (main window) ----
R6=$O/sweep/robustness
done_or $O/dumps/var_aggr.rds || Rscript R/06_make_variants.R $MAIN_F $O/dumps/var_aggr.rds --variant=aggr 2>&1 | tee -a "$LOG"
done_or $O/dumps/var_expdef.rds || Rscript R/06_make_variants.R $MAIN_F $O/dumps/var_expdef.rds --variant=expdef --anchor=onset --layers=$O/layers_dy.rds 2>&1 | tee -a "$LOG"
done_or $O/dumps/var_expdefentry.rds || Rscript R/06_make_variants.R $MAIN_F $O/dumps/var_expdefentry.rds --variant=expdef --anchor=entry --layers=$O/layers_dy.rds 2>&1 | tee -a "$LOG"
for L in hard reviewed; do done_or $O/dumps/var_ds$L.rds || Rscript R/06_make_variants.R $MAIN_F $O/dumps/var_ds$L.rds --variant=dropside --level=$L 2>&1 | tee -a "$LOG"; done
sweep $R6/base $MAIN_F $PA
sweep $R6/aggr $O/dumps/var_aggr.rds $PA
sweep $R6/expdef $O/dumps/var_expdef.rds $PA
sweep $R6/expdefentry $O/dumps/var_expdefentry.rds $PA
sweep $R6/dropside_hard $O/dumps/var_dshard.rds $PA
sweep $R6/dropside_reviewed $O/dumps/var_dsreviewed.rds $PA
for L in 5 10; do
  sweep $R6/lag${L}_yearscale $MAIN_F $PA --pos_lag=$L
  sweep $R6/samp${L}_yearscale $MAIN_F $PA --sample_lag=$L
  sweep $R6/lag${L}_pooled $MAIN_F $PA --pos_lag=$L --scale=pooled
  sweep $R6/samp${L}_pooled $MAIN_F $PA --sample_lag=$L --scale=pooled
done
sweep $R6/pooledscale $MAIN_F $PA --scale=pooled
sweep $R6/decomp10_linear $MAIN_F $PA --decomp=10
sweep $R6/decomp10_full $MAIN_F $PA --decomp=10 --decomp_terms=full
sweep $R6/decomp5_full $MAIN_F $PA --decomp=5 --decomp_terms=full
sweep $R6/yscale_tm1 $MAIN_F $PA --yscale=tm1
sweep $R6/unit_weights $O/dumps/long_${TYPES}_bon$(wtag $WINDOW_MAIN)_relyear_uw_fmid.rds $PA
sweep $R6/lead_t $MAIN_F $PA --pos_lag=0
sweep $R6/lead_t1 $MAIN_F $PA --pos_lag=-1
sweep $R6/lead_t1_sample_tm1 $MAIN_F $PA --sample_lag=-1
for K in $MERROR_K; do
  sweep $R6/merror_def_$K $MAIN_F $PA --merror=def --merror_k=$K --layers=$O/layers_dy.rds
  sweep $R6/merror_contig_$K $MAIN_F $PA --merror=contig --merror_k=$K --layers=$O/layers_dy.rds
  sweep $R6/merror_random_$K $MAIN_F $PA --merror=random --merror_k=$K
done
sweep $R6/noatop_base $MAIN_F $O/sets/noatop.rds
sweep $R6/noatop_expdef $O/dumps/var_expdef.rds $O/sets/noatop.rds
sweep $R6/stepB $MAIN_MID $PA --stepB=TRUE
sweep $R6/stepA $MAIN_MID $PA

## ---- 7. influence and top-dispute drops ----
done_or $O/influence/disp_loo.csv || Rscript R/09_influence.R $MAIN_F $O/sets/main.rds "$MAIN_SPEC" $O/influence/disp --unit=disp --cores=$CORES 2>&1 | tee -a "$LOG"
done_or $O/influence/init_loo.csv || Rscript R/09_influence.R $MAIN_F $O/sets/main.rds "$MAIN_SPEC" $O/influence/init --unit=init --cores=$CORES 2>&1 | tee -a "$LOG"
for N in 5 10; do
  IDS=$(head -$N $O/influence/disp_top.txt | paste -sd, -)
  done_or $O/dumps/var_droptop$N.rds || Rscript R/06_make_variants.R $MAIN_F $O/dumps/var_droptop$N.rds --variant=drop --disputes=$IDS 2>&1 | tee -a "$LOG"
  sweep $R6/drop_top$N $O/dumps/var_droptop$N.rds $PA
done

## ---- 8. resampling ----
rs_() { local NM=$1 DP=$2 PS=$3 SP=$4; shift 4; local F=$O/resample/$NM.csv; done_or $F || Rscript R/08_resample.R $DP $PS "$SP" $F "$@" --cores=$CORES 2>&1 | tail -4 | tee -a "$LOG"; }
for W in $WINDOWS; do rs_ boot_main_fmid_w$W $O/dumps/long_${TYPES}_bon$(wtag $W)_relyear_fmid.rds $O/sets/main.rds "$MAIN_SPEC" --mode=boot --B=$B_BOOT; done
rs_ boot_main_mid_w$WINDOW_MAIN $MAIN_MID $O/sets/main.rds "$MAIN_SPEC" --mode=boot --B=$B_BOOT
rs_ perm_dR_main $MAIN_F $O/sets/main.rds "$MAIN_SPEC" --mode=perm_dR --B=$B_PERM
rs_ perm_type_main $MAIN_F $O/sets/main.rds "$MAIN_SPEC" --mode=perm_type --B=$B_PERM
i=0; while read -r SP; do i=$((i+1)); rs_ boot_random_$(printf %02d $i) $MAIN_F $O/sets/random.rds "$SP" --mode=boot --B=$B_BOOT; done < $O/sets/random_ids.txt

## ---- 9. event study ----
done_or $O/eventstudy/es_main.csv || Rscript R/10_event_study.R $MAIN_F $O/sets/main.rds "$MAIN_SPEC" $O/eventstudy/es_main --B=$B_BOOT 2>&1 | tail -12 | tee -a "$LOG"
i=0; while read -r SP; do i=$((i+1)); done_or $O/eventstudy/es_random_$(printf %02d $i).csv || Rscript R/10_event_study.R $MAIN_F $O/sets/random.rds "$SP" $O/eventstudy/es_random_$(printf %02d $i) --B=200 2>&1 | tail -2 | tee -a "$LOG"; done < $O/sets/random_ids.txt

## ---- 10. summary ----
Rscript R/11_summarize.R $O/sweep $O/validation/measure_validation.csv "$MAIN_SPEC" $O/summary 2>&1 | tee -a "$LOG"
say "LOOP_DONE"
