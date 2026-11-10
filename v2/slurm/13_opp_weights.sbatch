#!/bin/bash
#SBATCH --job-name=v2_oppw
#SBATCH --account=open
#SBATCH --partition=open
#SBATCH --requeue
#SBATCH --time=48:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32   # = chains (4) x threads_per_chain (8)
#SBATCH --mem=128G
#SBATCH --array=0-31
#SBATCH --output=%x_%A_%a.out
#SBATCH --open-mode=append
# RUN 1 (built 2026-11-10; README section 14): the MID-propensity OPPORTUNITY-WEIGHT SENSITIVITY (not a consideration-set model). Three steps, chosen by STEP;
# submit them chained with v2/slurm/submit_opp_weights.sh (aux -> fit array afterok -> combine afterok). Nothing is submitted by this file.
#   STEP=aux      one job (submit with --array=0): 08_opportunity_weights.R for both foundations (--entry=later -> data_entry, --events=dymid -> data_dymid) and
#                 the outcomes in OPP_OUTCOMES (default "mid host3"); writes v2/output/oppw_<foundation>_<outcome>_lag3.rds and the diagnostics. Minutes per file.
#   STEP=fit      (default) the M1c cells: --comparison=directed --roles=bcj, lag 3, full choice set, tag $TAG (default long -> _tlong) with the long role-fit
#                 sampler settings of run 7 / run 0. Index = 8 x (2 x events + sample) + v:
#                   events 0 pairs (--entry=later, _ent; fmid = all-fatal), 1 dymid (--events=dymid, _dy; fmid = jointly fatal); sample 0 fmid, 1 mid
#                   v = 0 lambda 0 (_opw0; the current model, zero offset)   1 lambda 0.5 (_opw05)   2 lambda 1 at the auxiliary MLE (_opw1)
#                       3..7 lambda 1 at the auxiliary draws 1..5 (_opw1_d1 .. _opw1_d5)
#                   e.g. 0 fmid_lag3_all_pr_M1c_cdir_bcj_ent_opw0_tlong, 13 mid_lag3_all_pr_M1c_cdir_bcj_ent_opw1_d3_tlong, 18 fmid_lag3_all_pr_M1c_cdir_bcj_dy_opw1_tlong
#                 OPP_OUTCOME=host3 runs the same grid on the hostility >= 3 weights (suffix _oh3 after the weight); default mid.
#   STEP=combine  one job (submit with --array=0): 09_opp_combine.R for the four lambda = 1 cells (the equally weighted mixture over the five draws).
# The fit cells use the open account like run 0 (same data size and sampler settings as its full-set M1c cells); a preempted task is requeued and starts
# again from the beginning (log appended). Whether a long M1c fit ends inside 48 h on the open partition has NOT been measured here: check a run-0 log
# ("sampling finished in") before submitting, and move the array to the allocation (--account=jfe4_cr_default --partition=basic) if it does not.
# DRY RUN of the fit cells (needs the weights; open account): sbatch --export=ALL,DRY_RUN=TRUE --time=00:30:00 --mem=16G --cpus-per-task=1 v2/slurm/13_opp_weights.sbatch
# Needs both builds (01_build_data.R --entry=later and --events=dymid) and the compiled m1_context_offset program (compile_models.R).
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
STEP=${STEP:-fit}
TAG=${TAG:-long}
[[ "$TAG" =~ ^[A-Za-z0-9]+$ ]] || { echo "TAG must be letters and digits (got '$TAG')" >&2; exit 1; }
OUTC=${OPP_OUTCOME:-mid}
case $OUTC in mid|host3) ;; *) echo "OPP_OUTCOME must be mid or host3 (got '$OUTC')" >&2; exit 1 ;; esac
case $STEP in
  aux)
    for o in ${OPP_OUTCOMES:-mid host3}; do
      case $o in mid|host3) ;; *) echo "OPP_OUTCOMES holds '$o' (mid | host3)" >&2; exit 1 ;; esac
      for DATA_ARG in "--entry=later" "--events=dymid"; do
        echo "aux: 08_opportunity_weights.R $DATA_ARG --opp_outcome=$o"
        Rscript v2/R/08_opportunity_weights.R $DATA_ARG --opp_outcome=$o || exit 1
      done
    done
    ;;
  fit)
    i=${SLURM_ARRAY_TASK_ID:?run as an array}
    [[ "$i" =~ ^[0-9]+$ ]] && [ "$i" -le 31 ] || { echo "array index $i out of range 0-31" >&2; exit 1; }
    case $(((i / 8) / 2)) in 0) ev=pairs; DATA_ARG="--entry=later" ;; 1) ev=dymid; DATA_ARG="--events=dymid" ;; esac
    case $(((i / 8) % 2)) in 0) s=fmid ;; 1) s=mid ;; esac
    case $((i % 8)) in
      0) OPP="--opp_weight=0" ;;
      1) OPP="--opp_weight=0.5" ;;
      2) OPP="--opp_weight=1" ;;
      *) OPP="--opp_weight=1 --opp_draw=$((i % 8 - 2))" ;;
    esac
    [ "$OUTC" = "mid" ] || OPP="$OPP --opp_outcome=$OUTC"
    DRY=${DRY_RUN:-FALSE}
    THR=${THREADS_PER_CHAIN:-8}
    SAMP="--iter=${LONG_ITER:-6000} --warmup=${LONG_WARMUP:-3000} --adapt_delta=${LONG_ADAPT_DELTA:-0.9}"
    echo "task $i -> sample=$s model=M1c events=$ev ($DATA_ARG) $OPP tag=$TAG; $SAMP; threads_per_chain=$THR cpus=${SLURM_CPUS_PER_TASK:-?} dry_run=$DRY"
    Rscript v2/R/02_fit_model.R --sample=$s --lag=3 --model=M1c --comparison=directed $DATA_ARG --roles=bcj $OPP --tag=$TAG --dry_run=$DRY --threads_per_chain=$THR $SAMP
    ;;
  combine)
    OPP="--opp_weight=1"; [ "$OUTC" = "mid" ] || OPP="$OPP --opp_outcome=$OUTC"
    for DATA_ARG in "--entry=later" "--events=dymid"; do
      for s in fmid mid; do
        echo "combine: sample=$s $DATA_ARG $OPP tag=$TAG"
        Rscript v2/R/09_opp_combine.R --sample=$s --comparison=directed $DATA_ARG --roles=bcj $OPP --tag=$TAG || exit 1
      done
    done
    ;;
  *) echo "STEP must be aux, fit or combine (got '$STEP')" >&2; exit 1 ;;
esac
