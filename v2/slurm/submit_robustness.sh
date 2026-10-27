#!/bin/bash
# Submit one robustness preset of 04_robustness.sbatch with the M1 cells BEFORE the M2 cells (audit fix 3, 2026-10-26; mirrors
# submit_main_grid.sh). A plain `sbatch --export=ALL,PRESET=noww2 v2/slurm/04_robustness.sbatch` starts all six cells of the
# preset at once, so its M2 cells never find their M1 fits and fall back to random inits. This script splits the preset's cell
# list by family: the M1 cells as one array (CELLS = the M1 cells, --array=0-<n1-1>), then the M2 cells as a second array with
# --dependency=afterok on the first and --kill-on-invalid-dep=yes (a failed M1 task cancels the M2 array instead of leaving it
# pending). Each M2 task then starts from the M1 fit of its cell and switches (02_fit_model.R --init_from default).
#
# Usage (from the repo root):
#   PRESET=<noww2|sideA|noweights|fullcontrols|winsor|context|major_ix> [CELLS=sample:lag:model,...] bash v2/slurm/submit_robustness.sh [extra sbatch options]
# Environment:
#   PRESET            the preset (required). The refit presets split their cell list by family; `context` and `major_ix` are M1b only
#                     and are submitted as ONE array each (indices 0-5 of CONTEXT_PRESET / MAJOR_IX_PRESET, no dependency); `perm_error` is refused here (it needs
#                     the main-grid a / b fits, not a cell list; submit it directly, see 04_robustness.sbatch).
#   CELLS             the cell list of a refit preset, default = ROBUSTNESS_CELLS of spec.R (the list below must match it; tested).
#   M1_OPTS / M2_OPTS sbatch options for one array only (e.g. M2_OPTS="--time=10-00:00:00").
# CELLS and PRESET reach the job through the environment (sbatch --export=ALL), not through --export=ALL,CELLS=..., because the
# cell list contains commas, which --export treats as separators. --array in "$@" / M1_OPTS / M2_OPTS is refused.
# Examples:
#   PRESET=noww2 bash v2/slurm/submit_robustness.sh
#   PRESET=sideA bash v2/slurm/submit_robustness.sh                                   # needs v2/data_sideA
#   PRESET=winsor CELLS=mid:3:M1b,mid:3:M2b bash v2/slurm/submit_robustness.sh         # one sample
#   PRESET=context bash v2/slurm/submit_robustness.sh                                 # 6 M1b tasks, then Rscript v2/R/06_context_compare.R
#   PRESET=major_ix bash v2/slurm/submit_robustness.sh                                # 6 M1b tasks: fmid, mid x _major, _ix, _major_ix
#   PRESET=noww2 bash v2/slurm/submit_robustness.sh --account=open --export=ALL,DRY_RUN=TRUE --time=00:30:00 --mem=16G --cpus-per-task=1
set -eo pipefail
[ -f v2/slurm/04_robustness.sbatch ] || { echo "run from the repo root" >&2; exit 1; }
P=${PRESET:?set PRESET=<noww2|sideA|noweights|fullcontrols|winsor|context|major_ix>}
DEFAULT_CELLS="mid:3:M1b,fmid:3:M1b,war:3:M1b,mid:3:M2b,fmid:3:M2b,war:3:M2b"   # = ROBUSTNESS_CELLS in spec.R
for a in "$@" ${M1_OPTS:-} ${M2_OPTS:-}; do
  case "$a" in --array|--array=*|-a|-a*) echo "ERROR: do not pass --array; the index ranges follow the cell list (got '$a')" >&2; exit 1 ;; esac
done
jobid() { local raw=$1; local id=${raw%%;*}; [[ "$id" =~ ^[0-9]+$ ]] || { echo "ERROR: sbatch did not return a job id (got '$raw')" >&2; exit 1; }; echo "$id"; }
case $P in
  perm_error)
    echo "ERROR: perm_error is not a cell-list preset; submit it directly (see 04_robustness.sbatch):" >&2
    echo "  sbatch --export=ALL,PRESET=perm_error --account=open --cpus-per-task=4 --mem=32G --time=47:00:00 v2/slurm/04_robustness.sbatch" >&2
    exit 1 ;;
  context)
    [ -z "${CELLS:-}" ] || { echo "ERROR: the context preset fixes its cells (CONTEXT_PRESET in spec.R); unset CELLS" >&2; exit 1; }
    raw=$(PRESET=context sbatch --parsable --export=ALL --array=0-5 "$@" ${M1_OPTS:-} v2/slurm/04_robustness.sbatch)
    j=$(jobid "$raw")
    echo "context array submitted: job $j (indices 0-5: 0-1 bilateral, 2-3 coalition, 4-5 coalition_drop3; mid, fmid; M1b only, no dependency)" ;;
  major_ix)
    [ -z "${CELLS:-}" ] || { echo "ERROR: the major_ix preset fixes its cells (MAJOR_IX_PRESET in spec.R); unset CELLS" >&2; exit 1; }
    raw=$(PRESET=major_ix sbatch --parsable --export=ALL --array=0-5 "$@" ${M1_OPTS:-} v2/slurm/04_robustness.sbatch)
    j=$(jobid "$raw")
    echo "major_ix array submitted: job $j (indices 0-5: 0-1 major, 2-3 interact, 4-5 both; fmid, mid; M1b only, no dependency)" ;;
  noww2|sideA|noweights|fullcontrols|winsor)
    CELLS=${CELLS:-$DEFAULT_CELLS}
    m1_cells=""; m2_cells=""
    IFS=',' read -r -a cells <<< "$CELLS"
    for c in "${cells[@]}"; do
      model=${c##*:}
      case "$model" in
        M1[mrab]) m1_cells="${m1_cells:+$m1_cells,}$c" ;;
        M2[mrab]) m2_cells="${m2_cells:+$m2_cells,}$c" ;;
        *) echo "ERROR: cell '$c' is not sample:lag:model with a model of M1m|M1r|M1a|M1b|M2m|M2r|M2a|M2b" >&2; exit 1 ;;
      esac
    done
    n1=$(( $(tr -cd ',' <<< "$m1_cells" | wc -c) + 1 )); [ -n "$m1_cells" ] || n1=0
    n2=$(( $(tr -cd ',' <<< "$m2_cells" | wc -c) + 1 )); [ -n "$m2_cells" ] || n2=0
    m1=""
    if [ "$n1" -gt 0 ]; then
      raw=$(PRESET="$P" CELLS="$m1_cells" sbatch --parsable --export=ALL --array=0-$((n1 - 1)) "$@" ${M1_OPTS:-} v2/slurm/04_robustness.sbatch)
      m1=$(jobid "$raw")
      echo "preset $P, M1 array submitted: job $m1 (CELLS=$m1_cells, indices 0-$((n1 - 1)))"
    fi
    if [ "$n2" -gt 0 ]; then
      if [ -n "$m1" ]; then dep=(--dependency="afterok:$m1" --kill-on-invalid-dep=yes); else dep=(); echo "WARNING: no M1 cell in CELLS; the M2 cells start now and fall back to random inits unless their M1 fits are already in v2/output" >&2; fi
      raw=$(PRESET="$P" CELLS="$m2_cells" sbatch --parsable --export=ALL --array=0-$((n2 - 1)) "${dep[@]}" "$@" ${M2_OPTS:-} v2/slurm/04_robustness.sbatch)
      m2=$(jobid "$raw")
      echo "preset $P, M2 array submitted: job $m2 (CELLS=$m2_cells, indices 0-$((n2 - 1)))${m1:+, afterok:$m1, kill-on-invalid-dep=yes}"
    fi ;;
  *) echo "unknown PRESET $P" >&2; exit 1 ;;
esac
