#!/usr/bin/env bash
## ============================================================================
## run_all.sh -- the full multiverse loop (2026-11-19; restart and provenance rules rebuilt after the loop review P1-1..P1-3)
##   bash run_all.sh            production run into ./output
##   TEST=1 bash run_all.sh     mechanics check into ./output_test (2 windows, main measure only, tiny B) -- never mixes with ./output
## Provenance (P1-3): at start the run computes a FINGERPRINT of the code (all R scripts and this runner), config.sh, the data_entry manifest and the sha256 of every
## positions file. An output folder belongs to exactly one fingerprint: if ./output was built from anything else the run stops (start a new
## folder with OUT=... or remove the old one). Within one fingerprint a finished step is skipped; a step is "finished" only when its output was
## written atomically and, for sweeps, contains no ERROR rows. Sweeps with ERROR rows are recorded as INCOMPLETE and the run ends non-zero.
## Production refuses a dirty working tree (uncommitted code).
## Delta review 2026-11-20: the fingerprint also covers the layer inputs and file, R and package versions; artifacts built from supplied inputs
## (measure-id lists, exclusion lists, variants) carry a stamp of those inputs and are rebuilt when it differs; a sweep is keyed by the hashes of
## its dump and positions files; every R step must exit 0 (outputs are never accepted from existence alone); dumps are written straight to this
## run's folder (temp file + rename); the sweep manifest is rebuilt from validated per-task records; a non-empty folder this runner did not
## create is refused.
## ============================================================================
set -euo pipefail
cd "$(dirname "$0")"; source config.sh
## GPT stage-1 review S1: TYPES is environment-overridable, so only the schemes this runner is built for are accepted, and the resolved value enters the fingerprint
case "$TYPES" in bcd|init) ;; *) echo "TYPES must be bcd or init (got '$TYPES')" >&2; exit 1 ;; esac
TEST="${TEST:-0}"
if [ "$TEST" = "1" ]; then O="${OUT:-output_test}"; WINDOWS="0 1"; B_BOOT=20; B_PERM=20; N_RANDOM_BOOT=2; MERROR_K="0.5"; N_NOISE_DRAWS=1; else O="${OUT:-output}"; fi
## delta review N2: refuse a non-empty output folder this runner did not create (it cannot be tied to a fingerprint)
if [ -d "$O" ] && [ ! -e "$O/.created_by_run_all" ] && [ -n "$(find "$O" -type f -print -quit)" ]; then echo "STOP: $O is not empty and was not created by run_all.sh; use OUT=<new folder>"; exit 1; fi
mkdir -p $O/{positions,validation,sets/ids,dumps,sweep,resample,influence,eventstudy,summary,logs,status}; touch $O/.created_by_run_all
LOG=$O/logs/run_$(date +%Y%m%d_%H%M%S).log
say() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }
die() { say "STOP: $*"; exit 1; }
## run an R step: full output to the log, selected lines to the screen; returns the R exit status (delta review N3/N5: never inferred from files)
runR() { local PAT=$1; shift; local T; T=$(mktemp); local rc=0; Rscript "$@" > "$T" 2>&1 || rc=$?; cat "$T" >> "$LOG"; grep -E "$PAT" "$T" | cut -c1-250 || true; rm -f "$T"; return $rc; }
## content hashes, memoized per invocation
## content hashes (GPT review E1): memoized only WITHIN this invocation, in a fresh temporary folder (portable to macOS bash 3.2, no associative
## arrays). Nothing is reused across restarts, so a replaced file is always re-hashed by content. Within the invocation an entry is keyed by
## path + size + mtime AND explicitly invalidated (hinv) whenever this runner writes the file, so a rewrite in the same second cannot reuse it.
HC=$(mktemp -d "${TMPDIR:-/tmp}/run_all_hcache.XXXXXX"); trap 'rm -rf "$HC"' EXIT
hk() { printf '%s' "$(readlink -f "$1")" | shasum | cut -c1-16; }
h() { local f sz mt k c; f=$(readlink -f "$1"); sz=$(wc -c < "$f" | tr -d ' '); mt=$(date -r "$f" +%s); k=$(hk "$f")
  if [ -s "$HC/$k" ]; then read -r c < "$HC/$k"; if [ "${c% *}" = "$sz $mt" ]; then echo "${c##* }"; return 0; fi; fi
  c=$(shasum -a 256 "$f" | cut -c1-16); echo "$sz $mt $c" > "$HC/$k.tmp" && mv "$HC/$k.tmp" "$HC/$k"; echo "$c"; }
hm() { h "$1" > /dev/null; }
hinv() { local f; for f in "$@"; do rm -f "$HC/$(hk "$f")"; done; }
## stamps: an artifact built from supplied inputs is reused only if its .stamp matches the inputs now (delta review N2)
fresh() { [ -s "$1" ] && [ -s "$1.stamp" ] && [ "$(cat "$1.stamp")" = "$2" ]; }
stamp() { hinv "$1"; echo "$2" > "$1.stamp.tmp" && mv "$1.stamp.tmp" "$1.stamp"; }   # every stamped artifact was just (re)written

## ---- 0. preflight (P1-1) ----
for W in $WINDOWS; do [[ "$W" =~ ^[0-6]$ ]] || die "window $W outside 0-6"; done
[[ " $WINDOWS " == *" $WINDOW_MAIN "* ]] || die "main window $WINDOW_MAIN is not in WINDOWS"
PKGS="survival data.table dplyr igraph peacesciencer countrycode sandwich ggplot2 parallel fixest"; [ "$REUSE_POSITIONS" != "1" ] && PKGS="$PKGS igoR"   # N1: rebuild needs igoR; the new arms (2026-11-22) always need fixest
Rscript -e 'p <- strsplit(commandArgs(TRUE)[1], " ")[[1]]; m <- p[!vapply(p, requireNamespace, TRUE, quietly = TRUE)]; if (length(m)) stop("missing R packages: ", paste(m, collapse = ", "))' "$PKGS"
ENVQ=$(Rscript -e 'cat(R.version.string, "|"); for (p in strsplit(commandArgs(TRUE)[1], " ")[[1]]) cat(p, as.character(packageVersion(p)), "|")' "$PKGS")
if [ "$TEST" != "1" ] && ! git diff --quiet HEAD 2>/dev/null; then die "uncommitted changes in the repo; commit them first (production runs are tied to a commit)"; fi
if [ -e data_entry ] && [ "$(readlink data_entry)" != "$DATA_ENTRY" ]; then die "data_entry points to $(readlink data_entry), config says $DATA_ENTRY"; fi
[ -e data_entry ] || ln -s "$DATA_ENTRY" data_entry
for f in events.rds participants_timing.rds sy_positions.rds manifest.csv; do [ -s data_entry/$f ] || die "missing data_entry/$f"; done

## ---- 1. positions (link or build), then fingerprint ----
EXPECT="positions_core_W1raw_ptaall positions_core_W1raw_ptainforce positions_core_W2contig_ptaall positions_core_W2contig_ptainforce positions_core_W3dist_ptaall positions_core_W3dist_ptainforce positions_core_W4fe_ptaall positions_core_W4fe_ptainforce positions_dipl_W1raw_ptaall positions_dipl_W2contig_ptaall positions_dipl_W3dist_ptaall positions_dipl_W4fe_ptaall positions_igo_arm"
## PI 2026-11-22: book-aligned weightings (W6srybook, W7sryyear) under all three PTA codings, and unweighted PTA (ptaany) for W1-W4.
## These arms are always built in this run folder (never taken from the reuse folders, which predate them).
NEW_ANY="positions_core_W1raw_ptaany positions_core_W2contig_ptaany positions_core_W3dist_ptaany positions_core_W4fe_ptaany"
NEW_SRY_ALL="positions_core_W6srybook_ptaall positions_core_W7sryyear_ptaall"
NEW_SRY_INF="positions_core_W6srybook_ptainforce positions_core_W7sryyear_ptainforce"
NEW_SRY_ANY="positions_core_W6srybook_ptaany positions_core_W7sryyear_ptaany"
EXPECT="$EXPECT $NEW_ANY $NEW_SRY_ALL $NEW_SRY_INF $NEW_SRY_ANY"
link_pos() { local F=$1 D=$O/positions/$(basename "$1"); [ -s "$F" ] || die "positions source missing: $F"
  if [ -L "$D" ]; then [ "$(readlink "$D")" = "$F" ] || die "$D points to $(readlink "$D"), config says $F"; else ln -s "$F" "$D"; fi; }
if [ "$REUSE_POSITIONS" = "1" ]; then
  for F in "$REUSE_CORE_IGO"/positions_core_W{1raw,2contig,3dist,4fe}_pta{all,inforce}.rds "$REUSE_CORE_IGO"/positions_igo_arm.rds "$REUSE_DIPL"/positions_dipl_W{1raw,2contig,3dist,4fe}_ptaall.rds; do link_pos "$F"; done
fi
NETSHA=$( (for f in dipl_dyad.csv destaRaw.rda dca_dyad.csv; do [ -s "$NETDATA/$f" ] || die "missing $NETDATA/$f"; h "$NETDATA/$f"; done; h R/01_build_layers.R) | shasum -a 256 | cut -c1-16 )
## GPT review E1: an existing run folder records its network-input hash; a different input stops the run BEFORE layers are rebuilt in it
if [ -s $O/NETSHA ] && [ "$(cat $O/NETSHA)" != "$NETSHA" ]; then die "network inputs changed since $O was built (was $(cat $O/NETSHA), now $NETSHA); use OUT=<new folder>"; fi
if ! fresh $O/layers_dy.rds "$NETSHA"; then rm -f $O/layers_dy.rds; runR "Error|wrote" R/01_build_layers.R "$NETDATA" $O/layers_dy.rds.tmp.rds || die "layers build failed"; mv $O/layers_dy.rds.tmp.rds $O/layers_dy.rds; stamp $O/layers_dy.rds "$NETSHA"; fi
## GPT review D5: each arm is stamped with the layers file and builder hashes; an arm whose stamp differs (or is missing) is rebuilt, so an
## interrupted first build cannot mix arms made from different layer versions
AST="$(h $O/layers_dy.rds) $(h R/02_build_positions.R) $ENVQ"
arm() { local NM=$1 FILES=$2; shift 2; local ok=1 f; for f in $FILES; do [ -s $O/positions/$f.rds ] || ok=0; done
    if [ $ok = 1 ] && [ "$(cat $O/positions/.arm_$NM.stamp 2>/dev/null)" = "$AST" ]; then return 0; fi
    for f in $FILES; do rm -f $O/positions/$f.rds; done
    runR "Error|wrote" R/02_build_positions.R $O/layers_dy.rds $O/positions "$@" || die "positions arm $NM failed"
    for f in $FILES; do [ -s $O/positions/$f.rds ] || die "positions arm $NM did not write $f.rds"; hinv $O/positions/$f.rds; done
    echo "$AST" > $O/positions/.arm_$NM.stamp.tmp && mv $O/positions/.arm_$NM.stamp.tmp $O/positions/.arm_$NM.stamp; }
if [ "$REUSE_POSITIONS" != "1" ]; then
  arm core_all "positions_core_W1raw_ptaall positions_core_W2contig_ptaall positions_core_W3dist_ptaall positions_core_W4fe_ptaall" --arm=core --pta=all --cores=$CORES
  arm core_inforce "positions_core_W1raw_ptainforce positions_core_W2contig_ptainforce positions_core_W3dist_ptainforce positions_core_W4fe_ptainforce" --arm=core --pta=inforce --cores=$CORES
  arm dipl "positions_dipl_W1raw_ptaall positions_dipl_W2contig_ptaall positions_dipl_W3dist_ptaall positions_dipl_W4fe_ptaall" --arm=dipl --pta=all --cores=$CORES
  arm igo "positions_igo_arm" --arm=igo
fi
arm core_any "$NEW_ANY" --arm=core --pta=any --cores=$CORES
arm core_sry_all "$NEW_SRY_ALL" --arm=core --pta=all --weightings=W6srybook,W7sryyear --cores=$CORES
arm core_sry_inforce "$NEW_SRY_INF" --arm=core --pta=inforce --weightings=W6srybook,W7sryyear --cores=$CORES
arm core_sry_any "$NEW_SRY_ANY" --arm=core --pta=any --weightings=W6srybook,W7sryyear --cores=$CORES
for e in $EXPECT; do [ -s $O/positions/$e.rds ] || die "positions file missing: $e.rds"; done
POS=(); for e in $EXPECT; do POS+=( $O/positions/$e.rds ); done
SHA=$( (for e in $EXPECT; do h $O/positions/$e.rds; done) | shasum -a 256 | cut -c1-16 )
ENV=$ENVQ
## code identity = the content of the scripts and this runner (doc-only commits do not invalidate outputs); the commit is logged
FP=$( { cat run_all.sh R/*.R R/lib/*.R | shasum -a 256; cat config.sh; echo "TEST=$TEST TYPES=$TYPES WINDOWS=$WINDOWS B=$B_BOOT/$B_PERM N=$N_RANDOM_BOOT K=$MERROR_K"; shasum -a 256 data_entry/manifest.csv; echo "$SHA"; echo "layers $(h $O/layers_dy.rds) net $NETSHA"; echo "$ENV"; } | shasum -a 256 | cut -c1-16 )
if [ -s $O/FINGERPRINT ]; then [ "$(cat $O/FINGERPRINT)" = "$FP" ] || die "$O was built from different code/config/inputs (fingerprint $(cat $O/FINGERPRINT), now $FP); use OUT=<new folder> or remove $O"
else echo "$FP" > $O/FINGERPRINT; echo "$NETSHA" > $O/NETSHA; fi
say "fingerprint $FP | commit $(git rev-parse --short HEAD) | positions sha $SHA | cores $CORES | types $TYPES | main window $WINDOW_MAIN | test $TEST"
say "environment: $ENV"
MAN=$O/status/manifest_sweeps.tsv
INC=$O/status/incomplete.txt; : > $INC
CUR=$O/status/current_tasks.txt; : > $CUR

## ---- 2. validation (outcome-blind) ----
echo "$MAIN_SPEC" > $O/sets/main_id.txt
VPOS=( "${POS[@]}" ); [ "$TEST" = "1" ] && VPOS=( $O/positions/positions_core_W2contig_ptaall.rds )
[ -s $O/validation/measure_validation.csv ] || { runR "Error|wrote|plausible" R/03_validate_measures.R $O/validation "${VPOS[@]}" --top=$O/sets/main_id.txt --cores=$CORES || die "validation failed"; }
[ -s $O/validation/measure_validation.csv ] || die "validation failed"

## ---- 3. measure sets ----
## a set is rebuilt when its stamp (hash of the validation file, of the id list it was asked for, and its arguments) differs (delta review N2)
sub_set() { local NM=$1 ST="$(h $O/validation/measure_validation.csv) ${*:2}"; local A; for A in "${@:2}"; do [[ "$A" == --ids=* ]] && ST="$ST $(h "${A#--ids=}")"; done
  fresh $O/sets/$NM.rds "$ST" && return 0; rm -f $O/sets/$NM.rds
  runR "Error|wrote|set " R/04_subset_positions.R $O/positions $O/validation/measure_validation.csv $O/sets/$NM.rds "${@:2}" || die "set $NM failed"; stamp $O/sets/$NM.rds "$ST"; }
sub_set phaseA --set=ids --ids=inputs/phaseA_specs.txt
sub_set main --set=ids --ids=$O/sets/main_id.txt
sub_set noatop --set=noatop
RST="$(h $O/validation/measure_validation.csv) $MAIN_SPEC $N_RANDOM_BOOT"
fresh $O/sets/random_ids.txt "$RST" || { Rscript -e 'set.seed(123); v <- read.csv(commandArgs(TRUE)[1]); p <- setdiff(v$spec_id[v$plausible], commandArgs(TRUE)[3]); writeLines(sort(sample(p, as.integer(commandArgs(TRUE)[4]))), commandArgs(TRUE)[2])' \
  $O/validation/measure_validation.csv $O/sets/random_ids.txt "$MAIN_SPEC" $N_RANDOM_BOOT; stamp $O/sets/random_ids.txt "$RST"; }
sub_set random --set=ids --ids=$O/sets/random_ids.txt
## expected spec ids of every positions file and set (N6: the summarizer checks the analysed ids against these exactly); each id file is stamped with its source's hash
## GPT review D4: the helper fails (non-zero) unless the ids file exists with the CURRENT stamp; callers check it
ids_of() { local F=$1 I=$O/sets/ids/$(basename "$1" .rds).txt S; S=$(h "$F")
  if ! fresh "$I" "$S"; then rm -f "$I" "$I.stamp"; Rscript -e 'writeLines(colnames(readRDS(commandArgs(TRUE)[1])$R), commandArgs(TRUE)[2])' "$F" "$I.tmp" || return 1
    mv "$I.tmp" "$I" || return 1; stamp "$I" "$S" || return 1; fi
  fresh "$I" "$S" || return 1; echo "$I"; }
for F in "${POS[@]}" $O/sets/{phaseA,main,noatop,random}.rds; do ids_of "$F" > /dev/null || die "expected-ids file for $F could not be built"; done
cmp -s <(grep -v '^$' inputs/phaseA_specs.txt) $O/sets/ids/phaseA.txt || die "sets/phaseA.rds does not hold exactly the ids in inputs/phaseA_specs.txt"
cmp -s $O/sets/random_ids.txt $O/sets/ids/random.txt || die "sets/random.rds does not hold exactly the ids in random_ids.txt"
[ "$(cat $O/sets/ids/main.txt)" = "$MAIN_SPEC" ] || die "sets/main.rds does not hold exactly the main measure"

## ---- 4. event frames ----
wtag() { [ "$1" -gt 0 ] && echo "_w${1}d" || echo ""; }
dump() {  # dump <types> <window> <sample> [--unit_weights=TRUE]
  ## PI 2026-11-21: TYPES=init (stage 1: bilateral initiators vs coalition originators) is built from the four-type frame by R/06 --variant=initiators
  if [ "$1" = "init" ]; then local W=$2 SM=$3 UWT=""; [ "${4:-}" = "--unit_weights=TRUE" ] && UWT="_uw"
    shift; dump bcjs "$@"; local SRC=$O/dumps/long_bcjs_bon$(wtag $W)_relyear${UWT}_$SM.rds F=$O/dumps/long_init_bon$(wtag $W)_relyear${UWT}_$SM.rds
    [ -s $F ] && return 0
    runR "initiators|Error" R/06_make_variants.R $SRC $F.tmp.rds --variant=initiators || { rm -f $F.tmp.rds; die "initiator frame failed: $F"; }
    Rscript -e 'a <- commandArgs(TRUE); d <- readRDS(a[1]); stopifnot(d$types == "init", identical(d$roles, c("bil", "coal")), all(d$D$type %in% 0:1), d$sample == a[2], as.integer(d$coorig_months) == as.integer(a[3]))' \
      $F.tmp.rds $SM $W || { rm -f $F.tmp.rds; die "initiator frame $F is not the requested frame"; }
    mv $F.tmp.rds $F; hinv $F; return 0; fi
  local T=$1 W=$2 SM=$3; shift 3; local WF=""; [ "$W" -gt 0 ] && WF="--coorig_months=$W --coorig_bound=definite"
  local UW=""; [[ " $* " == *"--unit_weights=TRUE"* ]] && UW="_uw"
  local F=$O/dumps/long_${T}_bon$(wtag $W)_relyear${UW}_${SM}.rds
  [ -s $F ] && return 0
  ## delta review N3: written straight to this run's folder (temp + rename inside R/05); accepted only on exit 0 and a readable frame of the requested kind
  runR "participant table|P1-2|N7|check:|note:|promotions undone|types $T:|fatal|dump_long: wrote|Error" R/05_dump_events.R --comparison=directed --roles=bcj --entry=later --types=$T \
    --rel_scale=year --bilateral=onset $WF --sample=$SM --boot=0 --dump_long=TRUE --dump_out=$F "$@" || { rm -f $F; die "dump failed: $F"; }
  Rscript -e 'a <- commandArgs(TRUE); d <- readRDS(a[1]); stopifnot(is.data.frame(d$D), nrow(d$D) > 0, d$types == a[2], d$sample == a[3], as.integer(d$coorig_months) == as.integer(a[4]), d$bilateral_rule == "onset")' \
    $F $T $SM $W || { rm -f $F; die "dump $F is not the requested frame"; }
  hinv $F
}
for W in $WINDOWS; do dump $TYPES $W fmid; dump $TYPES $W mid; dump bcjs $W fmid; done
dump $TYPES $WINDOW_MAIN fmid --unit_weights=TRUE
MAIN_F=$O/dumps/long_${TYPES}_bon$(wtag $WINDOW_MAIN)_relyear_fmid.rds
MAIN_MID=$O/dumps/long_${TYPES}_bon$(wtag $WINDOW_MAIN)_relyear_mid.rds

## ---- sweep: one dump x one positions file ----
## Task record (delta review N5): status/sweep_<key> = "group <tab> csv <tab> n ids <tab> ids file", written atomically after the CSV is validated;
## the key includes the HASHES of the dump and positions files (N2), so a rebuilt variant or set reruns its sweeps. A record is trusted only if
## its CSV still exists without ERROR rows. The manifest is rebuilt from the records before the summary (never appended).
sweep() {  # sweep <group> <dump> <positions> [options]
  local G=$1 DP=$2 PS=$3; shift 3; local OD=$O/sweep/$G; mkdir -p $OD
  local KEY; KEY=$(echo "$G $(h $DP) $(h $PS) $*" | shasum | cut -c1-16); local ST=$O/status/sweep_$KEY
  echo "$KEY" >> $CUR    # GPT review D1: this invocation's planned tasks; only these can enter the manifest
  if [ -s $ST ]; then local F0; F0=$(cut -f2 $ST); [ -s "$F0" ] && ! grep -q '"ERROR"' "$F0" && return 0; fi
  local T; T=$(mktemp); local rc=0
  Rscript R/07_fit_sweep.R $DP $PS $OD --bilateral=onset --cores=$CORES "$@" > "$T" 2>&1 || rc=$?
  cat "$T" >> "$LOG"; grep -E "CENSUS|Error" "$T" | cut -c1-250 || true
  local F; F=$(sed -n 's/^WROTE //p' "$T" | tail -1); rm -f "$T"
  if [ "$rc" != 0 ] || [ -z "$F" ] || [ ! -s "$F" ]; then echo "$G $(basename $PS) $*: failed (exit $rc)" >> $INC; return 0; fi
  if grep -q '"ERROR"' "$F"; then echo "$G $(basename $PS) $*: ERROR rows in $F" >> $INC; return 0; fi
  local I; I=$(ids_of $PS) || { echo "$G $(basename $PS) $*: expected-ids file failed" >> $INC; return 0; }
  ## GPT review D3: the single atomic record is the only source of truth; an older record pointing to the same CSV is superseded and removed
  local R0; for R0 in $O/status/sweep_*; do [ -e "$R0" ] || continue; case "$R0" in *.tmp) continue;; esac
    [ "$R0" != "$ST" ] && [ "$(cut -f2 "$R0")" = "$F" ] && rm -f "$R0"; done
  printf "%s\t%s\t%s\t%s\n" "$G" "$F" "$(wc -l < $I | tr -d ' ')" "$I" > $ST.tmp && mv $ST.tmp $ST
}

## ---- 5. multiverse ----
MV=( "${POS[@]}" ); PA=$O/sets/phaseA.rds; [ "$TEST" = "1" ] && { MV=( $O/sets/main.rds ); PA=$O/sets/main.rds; }
for W in $WINDOWS; do for F in "${MV[@]}"; do
  sweep multiverse_fmid_w$W $O/dumps/long_${TYPES}_bon$(wtag $W)_relyear_fmid.rds $F
  sweep multiverse_mid_w$W  $O/dumps/long_${TYPES}_bon$(wtag $W)_relyear_mid.rds $F
done; done
for F in "${MV[@]}"; do sweep multiverse_stepB_w$WINDOW_MAIN $MAIN_MID $F --stepB=TRUE; done      # originator decisions at onset, dispute fatality
## intensity samples (PI 2026-11-20; full sweep): jointly fatal and dispute fatality level >= 2..6, every measure x every window, built from each
## window's fatal-MID frame; fits that fail are recorded INVALID per model and the run continues
varin() { local IN=$1 NM=$2; shift 2; hm $IN; local F=$O/dumps/var_$NM.rds ST; ST="$(h $IN) $*"; fresh $F "$ST" && return 0; rm -f $F
  runR "Error|sample " R/06_make_variants.R $IN $F "$@" || die "variant $NM failed"; [ -s $F ] || die "variant $NM wrote nothing"; stamp $F "$ST"; }
for W in $WINDOWS; do for S in $INTENSITY; do
  varin $O/dumps/long_${TYPES}_bon$(wtag $W)_relyear_fmid.rds int_${S}_w$W --variant=sample --keep=$S
  for F in "${MV[@]}"; do sweep multiverse_${S}_w$W $O/dumps/var_int_${S}_w$W.rds $F; done
done; done
for W in $WINDOWS; do sweep fourtypes_w$W $O/dumps/long_bcjs_bon$(wtag $W)_relyear_fmid.rds $PA; done

## ---- 6. robustness on the Phase A measures (main window) ----
## variants are stamped with the hash of their input dump and their arguments (incl. any dispute list), so a changed list rebuilds them (N2)
var() { local NM=$1; shift; hm $MAIN_F; local F=$O/dumps/var_$NM.rds ST; ST="$(h $MAIN_F) $*"; fresh $F "$ST" && return 0; rm -f $F
  runR "Error|dropside|drop:|expdef:|aggr" R/06_make_variants.R $MAIN_F $F "$@" || die "variant $NM failed"; [ -s $F ] || die "variant $NM wrote nothing"; stamp $F "$ST"; }
var aggr --variant=aggr --window=complete          # delta review P2-11: complete, GML-reconciled history windows only
var aggr_base --variant=aggr --window=complete_base  # the same strata without the control (comparator)
var aggr_censored --variant=aggr --window=censored   # every stratum; truncated counts after 2010 (lower-bound appendix sensitivity)
var expdef --variant=expdef --anchor=onset --layers=$O/layers_dy.rds
var expdefentry --variant=expdef --anchor=entry --layers=$O/layers_dy.rds
## GPT stage-1 review S2: kept under init too -- coalition originators include initiating-side joiners promoted by the window, whose side
## came from the case-reviewed map (side_case tags survive the init filter); R06 logs how many strata each level removes (0 is reported, not skipped)
DSL="hard reviewed hard_disp reviewed_disp"
for L in $DSL; do var ds_$L --variant=dropside --level=$L; done
sweep rob_base $MAIN_F $PA
sweep rob_aggr $O/dumps/var_aggr.rds $PA
sweep rob_aggr_base $O/dumps/var_aggr_base.rds $PA
sweep rob_aggr_censored $O/dumps/var_aggr_censored.rds $PA
sweep rob_expdef $O/dumps/var_expdef.rds $PA
sweep rob_expdefentry $O/dumps/var_expdefentry.rds $PA
for L in $DSL; do sweep rob_dropside_$L $O/dumps/var_ds_$L.rds $PA; done
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
  local IQ="chooser-bil:u_R_less,defender-chooser:u_R_less"; [ "$TYPES" = "init" ] && IQ="coal-bil:u_R_less"
  runR "LOO|refits|Error" R/09_influence.R $MAIN_F $O/sets/main.rds "$MAIN_SPEC" $O/influence/$1 --unit=$1 --q=$IQ --cores=$CORES || die "influence ($1) failed"
  [ -s $O/influence/$1_loo.csv ] && [ -s $O/influence/$1_top.txt ] || die "influence ($1) wrote no output"; }
inf disp; inf init
for N in 5 10; do var droptop$N --variant=drop --disputes=$(head -$N $O/influence/disp_top.txt | paste -sd, -); sweep rob_drop_top$N $O/dumps/var_droptop$N.rds $PA; done
for X in un strict strict_noborder expanded noruleA anycoder; do   # history-chosen lists (docs/historical_exclusions.md; PI approved 2026-11-20); run only if the list file exists
  if [ -s inputs/exclusions_$X.txt ]; then var hist_$X --variant=drop --disputes=$(grep -v '^#' inputs/exclusions_$X.txt | cut -d' ' -f1 | paste -sd, -); sweep rob_drop_hist_$X $O/dumps/var_hist_$X.rds $PA
  else say "skip: inputs/exclusions_$X.txt not present (list awaiting approval)"; echo "history drop $X: list not approved" >> $INC; fi
done

## ---- 8. resampling ----
rs() { local NM=$1 DP=$2 PS=$3 SP=$4; shift 4; local F=$O/resample/$NM.csv; [ -s $F ] && return 0
  runR "usable|Error" R/08_resample.R $DP $PS "$SP" $F "$@" --cores=$CORES || { rm -f $F; echo "resample $NM failed" >> $INC; return 0; }; [ -s $F ] || echo "resample $NM wrote nothing" >> $INC; }
for W in $WINDOWS; do rs boot_main_fmid_w$W $O/dumps/long_${TYPES}_bon$(wtag $W)_relyear_fmid.rds $O/sets/main.rds "$MAIN_SPEC" --mode=boot --B=$B_BOOT; done
rs boot_main_mid_w$WINDOW_MAIN $MAIN_MID $O/sets/main.rds "$MAIN_SPEC" --mode=boot --B=$B_BOOT
rs diag_perm_dR_main $MAIN_F $O/sets/main.rds "$MAIN_SPEC" --mode=perm_dR --B=$B_PERM
rs diag_perm_type_main $MAIN_F $O/sets/main.rds "$MAIN_SPEC" --mode=perm_type --B=$B_PERM
i=0; while read -r SP; do i=$((i+1)); rs boot_random_$(printf %02d $i) $MAIN_F $O/sets/random.rds "$SP" --mode=boot --B=$B_BOOT; done < $O/sets/random_ids.txt

## ---- 9. event study: both alignments; balanced and pre-only ----
es() { local NM=$1 PS=$2 SP=$3; shift 3; [ -s $O/eventstudy/$NM.csv ] && return 0
  runR "kept|align|Error" R/10_event_study.R $MAIN_F $PS "$SP" $O/eventstudy/$NM "$@" || { rm -f $O/eventstudy/$NM.csv; echo "event study $NM failed" >> $INC; return 0; }; [ -s $O/eventstudy/$NM.csv ] || echo "event study $NM wrote nothing" >> $INC; }
for A in entry onset; do
  es es_main_${A} $O/sets/main.rds "$MAIN_SPEC" --align=$A --B=$B_BOOT
  es es_main_${A}_preonly $O/sets/main.rds "$MAIN_SPEC" --align=$A --post=0 --B=$B_BOOT
  es es_main_${A}_unbalanced $O/sets/main.rds "$MAIN_SPEC" --align=$A --balanced=FALSE --B=$B_BOOT
  i=0; while read -r SP; do i=$((i+1)); es es_random_$(printf %02d $i)_${A} $O/sets/random.rds "$SP" --align=$A --B=200; done < $O/sets/random_ids.txt
done
es es_main_onset_fallback $O/sets/main.rds "$MAIN_SPEC" --align=onset --onset_fallback=TRUE --B=$B_BOOT   # missing/conflicting starts at the event year (flagged sensitivity)

## ---- 10. summary: manifest rebuilt from the validated task records (N5), then the summarizer reads only those ----
{ printf "group\tpath\tn_requested\tids_file\n"
  while read -r K; do ST=$O/status/sweep_$K                          # D1: planned tasks of THIS invocation only (superseded records never qualify)
    [ -s "$ST" ] || continue; F0=$(cut -f2 "$ST"); if [ -s "$F0" ] && ! grep -q '"ERROR"' "$F0"; then cat "$ST"; fi
  done < <(sort -u $CUR); } > $MAN.tmp
mv $MAN.tmp $MAN
runR "summarized|Error" R/11_summarize.R $MAN $O/validation/measure_validation.csv "$MAIN_SPEC" $O/summary || { echo "summary failed" >> $INC; }
if [ -s $INC ]; then say "LOOP_FINISHED_WITH_GAPS: $(wc -l < $INC) items in $INC"; cat $INC | tee -a "$LOG"; exit 2; fi
say "LOOP_DONE"
