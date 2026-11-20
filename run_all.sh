#!/usr/bin/env bash
## ============================================================================
## run_all.sh -- the full multiverse loop (2026-11-19; restart and provenance rules rebuilt after the loop review P1-1..P1-3)
##   bash run_all.sh            production run into ./output
##   TEST=1 bash run_all.sh     mechanics check into ./output_test (2 windows, main measure only, tiny B) -- never mixes with ./output
## Provenance (P1-3): at start the run computes a FINGERPRINT of the code commit, config.sh, the data_entry manifest and the sha256 of every
## positions file. An output folder belongs to exactly one fingerprint: if ./output was built from anything else the run stops (start a new
## folder with OUT=... or remove the old one). Within one fingerprint a finished step is skipped; a step is "finished" only when its output was
## written atomically and, for sweeps, contains no ERROR rows. Sweeps with ERROR rows are recorded as INCOMPLETE and the run ends non-zero.
## Production refuses a dirty working tree (uncommitted code).
## ============================================================================
set -euo pipefail
cd "$(dirname "$0")"; source config.sh
TEST="${TEST:-0}"
if [ "$TEST" = "1" ]; then O="${OUT:-output_test}"; WINDOWS="0 1"; B_BOOT=20; B_PERM=20; N_RANDOM_BOOT=2; MERROR_K="0.5"; N_NOISE_DRAWS=1; else O="${OUT:-output}"; fi
mkdir -p $O/{positions,validation,sets,dumps,sweep,resample,influence,eventstudy,summary,logs,status} output/rel_sweep
LOG=$O/logs/run_$(date +%Y%m%d_%H%M%S).log
say() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }
die() { say "STOP: $*"; exit 1; }

## ---- 0. preflight (P1-1) ----
for W in $WINDOWS; do [[ "$W" =~ ^[0-6]$ ]] || die "window $W outside 0-6"; done
[[ " $WINDOWS " == *" $WINDOW_MAIN "* ]] || die "main window $WINDOW_MAIN is not in WINDOWS"
Rscript -e 'p <- c("survival","data.table","dplyr","igraph","peacesciencer","countrycode","sandwich","ggplot2","parallel"); m <- p[!vapply(p, requireNamespace, TRUE, quietly = TRUE)]; if (length(m)) stop("missing R packages: ", paste(m, collapse = ", "))'
if [ "$TEST" != "1" ] && ! git diff --quiet HEAD 2>/dev/null; then die "uncommitted changes in the repo; commit them first (production runs are tied to a commit)"; fi
if [ -e data_entry ] && [ "$(readlink data_entry)" != "$DATA_ENTRY" ]; then die "data_entry points to $(readlink data_entry), config says $DATA_ENTRY"; fi
[ -e data_entry ] || ln -s "$DATA_ENTRY" data_entry
for f in events.rds participants_timing.rds sy_positions.rds manifest.csv; do [ -s data_entry/$f ] || die "missing data_entry/$f"; done

## ---- 1. positions (link or build), then fingerprint ----
EXPECT="positions_core_W1raw_ptaall positions_core_W1raw_ptainforce positions_core_W2contig_ptaall positions_core_W2contig_ptainforce positions_core_W3dist_ptaall positions_core_W3dist_ptainforce positions_core_W4fe_ptaall positions_core_W4fe_ptainforce positions_dipl_W1raw_ptaall positions_dipl_W2contig_ptaall positions_dipl_W3dist_ptaall positions_dipl_W4fe_ptaall positions_igo_arm"
link_pos() { local F=$1 D=$O/positions/$(basename "$1"); [ -s "$F" ] || die "positions source missing: $F"
  if [ -L "$D" ]; then [ "$(readlink "$D")" = "$F" ] || die "$D points to $(readlink "$D"), config says $F"; else ln -s "$F" "$D"; fi; }
if [ "$REUSE_POSITIONS" = "1" ]; then
  for F in "$REUSE_CORE_IGO"/positions_core_*.rds "$REUSE_CORE_IGO"/positions_igo_arm.rds "$REUSE_DIPL"/positions_dipl_*.rds; do link_pos "$F"; done
fi
[ -s $O/layers_dy.rds ] || { Rscript R/01_build_layers.R "$NETDATA" $O/layers_dy.rds.tmp.rds 2>&1 | tee -a "$LOG"; mv $O/layers_dy.rds.tmp.rds $O/layers_dy.rds; }
if [ "$REUSE_POSITIONS" != "1" ]; then
  for PT in all inforce; do for W in W1raw W2contig W3dist W4fe; do [ -s $O/positions/positions_core_${W}_pta$PT.rds ] || NEED_CORE_$PT=1; done; done
  [ "${NEED_CORE_all:-0}" = 1 ] && Rscript R/02_build_positions.R $O/layers_dy.rds $O/positions --arm=core --pta=all --cores=$CORES 2>&1 | tee -a "$LOG"
  [ "${NEED_CORE_inforce:-0}" = 1 ] && Rscript R/02_build_positions.R $O/layers_dy.rds $O/positions --arm=core --pta=inforce --cores=$CORES 2>&1 | tee -a "$LOG"
  ls $O/positions/positions_dipl_W4fe_ptaall.rds >/dev/null 2>&1 || Rscript R/02_build_positions.R $O/layers_dy.rds $O/positions --arm=dipl --pta=all --cores=$CORES 2>&1 | tee -a "$LOG"
  [ -s $O/positions/positions_igo_arm.rds ] || Rscript R/02_build_positions.R $O/layers_dy.rds $O/positions --arm=igo 2>&1 | tee -a "$LOG"
fi
for e in $EXPECT; do [ -s $O/positions/$e.rds ] || die "positions file missing: $e.rds"; done
POS=(); for e in $EXPECT; do POS+=( $O/positions/$e.rds ); done
SHA=$( (cd $O/positions && for e in $EXPECT; do shasum -a 256 -- "$(readlink -f $e.rds)" | cut -d' ' -f1; done) | shasum -a 256 | cut -c1-16 )
FP=$( { git rev-parse HEAD; cat config.sh; echo "TEST=$TEST WINDOWS=$WINDOWS B=$B_BOOT/$B_PERM N=$N_RANDOM_BOOT K=$MERROR_K"; shasum -a 256 data_entry/manifest.csv; echo "$SHA"; } | shasum -a 256 | cut -c1-16 )
if [ -s $O/FINGERPRINT ]; then [ "$(cat $O/FINGERPRINT)" = "$FP" ] || die "$O was built from different code/config/inputs (fingerprint $(cat $O/FINGERPRINT), now $FP); use OUT=<new folder> or remove $O"
else echo "$FP" > $O/FINGERPRINT; fi
say "fingerprint $FP | commit $(git rev-parse --short HEAD) | positions sha $SHA | cores $CORES | types $TYPES | main window $WINDOW_MAIN | test $TEST"
Rscript -e 'cat(R.version.string, "\n"); for (p in c("survival","data.table","dplyr","igraph","peacesciencer","countrycode","sandwich","ggplot2")) cat(p, as.character(packageVersion(p)), "\n")' | tee -a "$LOG"
MAN=$O/status/manifest_sweeps.tsv; [ -s $MAN ] || printf "group\tpath\tn_requested\n" > $MAN
INC=$O/status/incomplete.txt; : > $INC
nspec() { Rscript -e 'cat(ncol(readRDS(commandArgs(TRUE)[1])$R))' "$1"; }

## ---- 2. validation (outcome-blind) ----
echo "$MAIN_SPEC" > $O/sets/main_id.txt
VPOS=( "${POS[@]}" ); [ "$TEST" = "1" ] && VPOS=( $O/positions/positions_core_W2contig_ptaall.rds )
[ -s $O/validation/measure_validation.csv ] || Rscript R/03_validate_measures.R $O/validation "${VPOS[@]}" --top=$O/sets/main_id.txt --cores=$CORES 2>&1 | tee -a "$LOG"
[ -s $O/validation/measure_validation.csv ] || die "validation failed"

## ---- 3. measure sets ----
sub_set() { [ -s $O/sets/$1.rds ] || Rscript R/04_subset_positions.R $O/positions $O/validation/measure_validation.csv $O/sets/$1.rds "${@:2}" 2>&1 | tee -a "$LOG"; [ -s $O/sets/$1.rds ] || die "set $1 failed"; }
sub_set phaseA --set=ids --ids=inputs/phaseA_specs.txt
sub_set main --set=ids --ids=$O/sets/main_id.txt
sub_set noatop --set=noatop
[ -s $O/sets/random_ids.txt ] || Rscript -e 'set.seed(123); v <- read.csv(commandArgs(TRUE)[1]); p <- setdiff(v$spec_id[v$plausible], commandArgs(TRUE)[3]); writeLines(sort(sample(p, as.integer(commandArgs(TRUE)[4]))), commandArgs(TRUE)[2])' \
  $O/validation/measure_validation.csv $O/sets/random_ids.txt "$MAIN_SPEC" $N_RANDOM_BOOT
sub_set random --set=ids --ids=$O/sets/random_ids.txt

## ---- 4. event frames ----
wtag() { [ "$1" -gt 0 ] && echo "_w${1}d" || echo ""; }
dump() {  # dump <types> <window> <sample> [--unit_weights=TRUE]
  local T=$1 W=$2 SM=$3; shift 3; local WF=""; [ "$W" -gt 0 ] && WF="--coorig_months=$W --coorig_bound=definite"
  local UW=""; [[ " $* " == *"--unit_weights=TRUE"* ]] && UW="_uw"
  local F=$O/dumps/long_${T}_bon$(wtag $W)_relyear${UW}_${SM}.rds
  [ -s $F ] && return 0
  Rscript R/05_dump_events.R --comparison=directed --roles=bcj --entry=later --types=$T --rel_scale=year --bilateral=onset $WF --sample=$SM --boot=0 --dump_long=TRUE "$@" 2>&1 \
    | grep -E "participant table|P1-2|check:|note:|promotions undone|types $T:|fatal|dump_long: wrote|Error" | tee -a "$LOG" || true
  [ -s output/rel_sweep/$(basename $F) ] || die "dump failed: $F"; mv output/rel_sweep/$(basename $F) $F
}
for W in $WINDOWS; do dump $TYPES $W fmid; dump $TYPES $W mid; dump bcjs $W fmid; done
dump $TYPES $WINDOW_MAIN fmid --unit_weights=TRUE
MAIN_F=$O/dumps/long_${TYPES}_bon$(wtag $WINDOW_MAIN)_relyear_fmid.rds
MAIN_MID=$O/dumps/long_${TYPES}_bon$(wtag $WINDOW_MAIN)_relyear_mid.rds

## ---- sweep: one dump x one positions file. Finished = an atomically written CSV with no ERROR rows; recorded in the manifest ----
sweep() {  # sweep <group> <dump> <positions> [options]
  local G=$1 DP=$2 PS=$3; shift 3; local OD=$O/sweep/$G; mkdir -p $OD
  local KEY=$(echo "$G $(basename $DP) $(basename $PS) $*" | shasum | cut -c1-16); local ST=$O/status/sweep_$KEY
  if [ -s $ST ] && [ -s "$(cat $ST)" ]; then return 0; fi
  local OUTL; OUTL=$(Rscript R/07_fit_sweep.R $DP $PS $OD --bilateral=onset --cores=$CORES "$@" 2>&1 | tee -a "$LOG" | grep -E "^WROTE |CENSUS|Error" || true)
  echo "$OUTL" | grep -E "CENSUS|Error" | cut -c1-250
  local F; F=$(echo "$OUTL" | sed -n 's/^WROTE //p' | tail -1)
  if [ -z "$F" ] || [ ! -s "$F" ]; then echo "$G $(basename $PS) $*: no output" >> $INC; return 0; fi
  if grep -q '"ERROR"' "$F"; then echo "$G $(basename $PS) $*: ERROR rows in $F" >> $INC; return 0; fi
  printf "%s\t%s\t%s\n" "$G" "$F" "$(nspec $PS)" >> $MAN; echo "$F" > $ST
}

## ---- 5. multiverse ----
MV=( "${POS[@]}" ); PA=$O/sets/phaseA.rds; [ "$TEST" = "1" ] && { MV=( $O/sets/main.rds ); PA=$O/sets/main.rds; }
for W in $WINDOWS; do for F in "${MV[@]}"; do
  sweep multiverse_fmid_w$W $O/dumps/long_${TYPES}_bon$(wtag $W)_relyear_fmid.rds $F
  sweep multiverse_mid_w$W  $O/dumps/long_${TYPES}_bon$(wtag $W)_relyear_mid.rds $F
done; done
for F in "${MV[@]}"; do sweep multiverse_stepB_w$WINDOW_MAIN $MAIN_MID $F --stepB=TRUE; done      # originator decisions at onset, dispute fatality
for W in $WINDOWS; do sweep fourtypes_w$W $O/dumps/long_bcjs_bon$(wtag $W)_relyear_fmid.rds $PA; done

## ---- 6. robustness on the Phase A measures (main window) ----
var() { local NM=$1; shift; [ -s $O/dumps/var_$NM.rds ] || Rscript R/06_make_variants.R $MAIN_F $O/dumps/var_$NM.rds "$@" 2>&1 | tee -a "$LOG" | grep -E "Error|dropside|drop:|expdef:|aggr:" || true; [ -s $O/dumps/var_$NM.rds ] || die "variant $NM failed"; }
var aggr --variant=aggr
var expdef --variant=expdef --anchor=onset --layers=$O/layers_dy.rds
var expdefentry --variant=expdef --anchor=entry --layers=$O/layers_dy.rds
for L in hard reviewed hard_disp reviewed_disp; do var ds_$L --variant=dropside --level=$L; done
sweep rob_base $MAIN_F $PA
sweep rob_aggr $O/dumps/var_aggr.rds $PA
sweep rob_expdef $O/dumps/var_expdef.rds $PA
sweep rob_expdefentry $O/dumps/var_expdefentry.rds $PA
for L in hard reviewed hard_disp reviewed_disp; do sweep rob_dropside_$L $O/dumps/var_ds_$L.rds $PA; done
for L in 5 10; do
  sweep rob_lag${L}_yearscale $MAIN_F $PA --pos_lag=$L
  sweep rob_samp${L}_yearscale $MAIN_F $PA --sample_lag=$L
  sweep rob_lag${L}_pooled $MAIN_F $PA --pos_lag=$L --scale=pooled
  sweep rob_samp${L}_pooled $MAIN_F $PA --sample_lag=$L --scale=pooled
done
sweep rob_pooledscale $MAIN_F $PA --scale=pooled
sweep rob_decomp10_linear $MAIN_F $PA --decomp=10
sweep rob_decomp10_full $MAIN_F $PA --decomp=10 --decomp_terms=full
sweep rob_decomp5_full $MAIN_F $PA --decomp=5 --decomp_terms=full
sweep rob_yscale_tm1 $MAIN_F $PA --yscale=tm1
sweep rob_unit_weights $O/dumps/long_${TYPES}_bon$(wtag $WINDOW_MAIN)_relyear_uw_fmid.rds $PA
sweep rob_lead_t $MAIN_F $PA --pos_lag=0
sweep rob_lead_t1 $MAIN_F $PA --pos_lag=-1
sweep rob_lead_t1_sample_tm1 $MAIN_F $PA --sample_lag=-1
for K in $MERROR_K; do
  for S in 1 -1; do
    sweep rob_merror_def_${S}_$K $MAIN_F $PA --merror=def --merror_k=$K --merror_sign=$S --layers=$O/layers_dy.rds
    sweep rob_merror_contig_${S}_$K $MAIN_F $PA --merror=contig --merror_k=$K --merror_sign=$S --layers=$O/layers_dy.rds
  done
  for D in $(seq 1 $N_NOISE_DRAWS); do sweep rob_merror_random_${K}_d$D $MAIN_F $PA --merror=random --merror_k=$K --merror_draw=$D; done
done
sweep rob_noatop_base $MAIN_F $O/sets/noatop.rds
sweep rob_noatop_expdef $O/dumps/var_expdef.rds $O/sets/noatop.rds
sweep rob_stepA $MAIN_MID $PA
sweep rob_stepB_orig_disp $MAIN_MID $PA --stepB=TRUE
sweep rob_stepB_all_pair $MAIN_MID $PA --stepB=TRUE --stepB_unit=all --stepB_outcome=pair      # old membership version, comparison only

## ---- 7. influence; influence-chosen and history-chosen dispute drops ----
inf() { [ -s $O/influence/$1_loo.csv ] && [ -s $O/influence/$1_top.txt ] && return 0
  Rscript R/09_influence.R $MAIN_F $O/sets/main.rds "$MAIN_SPEC" $O/influence/$1 --unit=$1 --cores=$CORES 2>&1 | tee -a "$LOG" | grep -E "LOO|refits|Error" || true
  [ -s $O/influence/$1_loo.csv ] && [ -s $O/influence/$1_top.txt ] || die "influence ($1) failed"; }
inf disp; inf init
for N in 5 10; do var droptop$N --variant=drop --disputes=$(head -$N $O/influence/disp_top.txt | paste -sd, -); sweep rob_drop_top$N $O/dumps/var_droptop$N.rds $PA; done
for X in un punitive; do   # history-chosen lists (docs/historical_exclusions.md); run only once the list file exists (approved)
  if [ -s inputs/exclusions_$X.txt ]; then var hist_$X --variant=drop --disputes=$(grep -v '^#' inputs/exclusions_$X.txt | cut -d' ' -f1 | paste -sd, -); sweep rob_drop_hist_$X $O/dumps/var_hist_$X.rds $PA
  else say "skip: inputs/exclusions_$X.txt not present (list awaiting approval)"; echo "history drop $X: list not approved" >> $INC; fi
done

## ---- 8. resampling ----
rs() { local NM=$1 DP=$2 PS=$3 SP=$4; shift 4; local F=$O/resample/$NM.csv; [ -s $F ] && return 0
  Rscript R/08_resample.R $DP $PS "$SP" $F "$@" --cores=$CORES 2>&1 | tee -a "$LOG" | grep -E "usable|Error" || true; [ -s $F ] || echo "resample $NM failed" >> $INC; }
for W in $WINDOWS; do rs boot_main_fmid_w$W $O/dumps/long_${TYPES}_bon$(wtag $W)_relyear_fmid.rds $O/sets/main.rds "$MAIN_SPEC" --mode=boot --B=$B_BOOT; done
rs boot_main_mid_w$WINDOW_MAIN $MAIN_MID $O/sets/main.rds "$MAIN_SPEC" --mode=boot --B=$B_BOOT
rs diag_perm_dR_main $MAIN_F $O/sets/main.rds "$MAIN_SPEC" --mode=perm_dR --B=$B_PERM
rs diag_perm_type_main $MAIN_F $O/sets/main.rds "$MAIN_SPEC" --mode=perm_type --B=$B_PERM
i=0; while read -r SP; do i=$((i+1)); rs boot_random_$(printf %02d $i) $MAIN_F $O/sets/random.rds "$SP" --mode=boot --B=$B_BOOT; done < $O/sets/random_ids.txt

## ---- 9. event study: both alignments; balanced and pre-only ----
es() { local NM=$1 PS=$2 SP=$3; shift 3; [ -s $O/eventstudy/$NM.csv ] && return 0
  Rscript R/10_event_study.R $MAIN_F $PS "$SP" $O/eventstudy/$NM "$@" 2>&1 | tee -a "$LOG" | grep -E "kept|align|Error" || true; [ -s $O/eventstudy/$NM.csv ] || echo "event study $NM failed" >> $INC; }
for A in entry onset; do
  es es_main_${A} $O/sets/main.rds "$MAIN_SPEC" --align=$A --B=$B_BOOT
  es es_main_${A}_preonly $O/sets/main.rds "$MAIN_SPEC" --align=$A --post=0 --B=$B_BOOT
  es es_main_${A}_unbalanced $O/sets/main.rds "$MAIN_SPEC" --align=$A --balanced=FALSE --B=$B_BOOT
  i=0; while read -r SP; do i=$((i+1)); es es_random_$(printf %02d $i)_${A} $O/sets/random.rds "$SP" --align=$A --B=200; done < $O/sets/random_ids.txt
done

## ---- 10. summary (manifest-listed sweeps only) ----
Rscript R/11_summarize.R $MAN $O/validation/measure_validation.csv "$MAIN_SPEC" $O/summary 2>&1 | tee -a "$LOG"
if [ -s $INC ]; then say "LOOP_FINISHED_WITH_GAPS: $(wc -l < $INC) items in $INC"; cat $INC | tee -a "$LOG"; exit 2; fi
say "LOOP_DONE"
