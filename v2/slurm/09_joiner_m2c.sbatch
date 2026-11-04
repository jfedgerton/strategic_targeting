#!/bin/bash
#SBATCH --job-name=v2_jset_m2c
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --requeue
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32   # = chains (4) x threads_per_chain (8)
#SBATCH --mem=128G                                 # basic bills 4 GB per core: 32 cores x 4 GB
#SBATCH --array=0-21
#SBATCH --output=%x_%A_%a.out
#SBATCH --open-mode=append
# ENTRY DIRECTION, ROLE CODINGS bc / bcj, JOINER CHOICE SETS AND THE AME ROLE MODEL M2c (PI decisions 2026-11-04). Every cell: --entry=later
# (v2/data_entry), --comparison=directed, lag 3; the role cells take --roles=$ROLES, ROLES=bcj by default (the MAIN coding, PI decision 2026-11-04) or ROLES=bc (optional appendix coding; one submission per coding).
#   idx  sample  model  K  joiner_set  needs (init fit in v2/output; M2c only)
#    0   fmid    M1c    -  all         -   (baseline and M2c init)
#    1   mid     M1c    -  all         -
#    2   fmid    M1cr   -  all         -   (baseline)
#    3   mid     M1cr   -  all         -
#    4   fmid    M1c    -  opposing    -
#    5   mid     M1c    -  opposing    -
#    6   fmid    M1b    -  -           -   (the pooled baseline on the entry events; no roles, so its cell is the same for both codings: run it ONCE)
#    7   mid     M1b    -  -           -
#    8   fmid    M2c    1  all         fit_fmid_lag3_all_pr_M1c_cdir_<ROLES>_ent_t<TAG>.rds        (task 0)
#    9   mid     M2c    1  all         fit_mid_lag3_all_pr_M1c_cdir_<ROLES>_ent_t<TAG>.rds         (task 1)
#   10   fmid    M2c    2  all         fit_fmid_lag3_all_pr_M1c_cdir_<ROLES>_ent_t<TAG>.rds        (task 0)
#   11   mid     M2c    2  all         fit_mid_lag3_all_pr_M1c_cdir_<ROLES>_ent_t<TAG>.rds         (task 1)
#   12   fmid    M2c    1  opposing    fit_fmid_lag3_all_pr_M1c_cdir_<ROLES>_jset_ent_t<TAG>.rds   (task 4)
#   13   mid     M2c    1  opposing    fit_mid_lag3_all_pr_M1c_cdir_<ROLES>_jset_ent_t<TAG>.rds    (task 5)
#   14   fmid    M2c    2  opposing    fit_fmid_lag3_all_pr_M1c_cdir_<ROLES>_jset_ent_t<TAG>.rds   (task 4)
#   15   mid     M2c    2  opposing    fit_mid_lag3_all_pr_M1c_cdir_<ROLES>_jset_ent_t<TAG>.rds    (task 5)
#  politically relevant choice set (--variant=polrel; cell names <sample>_lag3_polrel_pr_...; joiner_set all):
#   16   fmid    M1c    -  all         -   (polrel baseline and M2c init)
#   17   mid     M1c    -  all         -
#   18   fmid    M2c    1  all         fit_fmid_lag3_polrel_pr_M1c_cdir_<ROLES>_ent_t<TAG>.rds     (task 16)
#   19   mid     M2c    1  all         fit_mid_lag3_polrel_pr_M1c_cdir_<ROLES>_ent_t<TAG>.rds      (task 17)
#   20   fmid    M2c    2  all         fit_fmid_lag3_polrel_pr_M1c_cdir_<ROLES>_ent_t<TAG>.rds     (task 16)
#   21   mid     M2c    2  all         fit_mid_lag3_polrel_pr_M1c_cdir_<ROLES>_ent_t<TAG>.rds      (task 17)
# Dependencies: tasks 8-15 need the M1c fits of tasks 0, 1, 4, 5 and tasks 18-21 those of tasks 16-17, of the SAME ROLES -> submit_joiner_m2c.sh
# submits the M1 tasks (0-7,16-17, or 0-5,16-17 with SKIP_M1B=TRUE for the second coding) as one array and the M2c tasks 8-15,18-21 with --dependency=afterok:<that array> --kill-on-invalid-dep=yes (one
# dependency on the whole M1 array, as in submit_main_grid.sh: a failed M1cr / M1b task also cancels the M2c array; resubmit with A_JOB=<id>).
# 02_fit_model.R REFUSES to sample an M2c cell whose M1c fit is missing (no random-init fallback) and refuses an M1c fit of another roles /
# joiner_set / entry / tag / cell.
# Sampler settings:
#   M1c / M1cr / M1b (tasks 0-7): the long role-fit settings, LONG_ITER / LONG_WARMUP / LONG_ADAPT_DELTA / TAG (defaults 6000 / 3000 / 0.9 / long).
#     These defaults are the README's documented long-refit example (--iter=6000 --warmup=3000 --tag=long), NOT read from the existing long role fits,
#     which are not in the repository: before sampling, each M1c / M1cr task reads the existing aprime long fit of the same model and sample on the
#     keep data (fit_<s>_lag3_all_pr_<M1c|M1cr>_cdir_aprime_t<TAG>.rds, copy it into v2/output) and STOPS unless its recorded iter / warmup /
#     adapt_delta equal these; set CHECK_LONG=FALSE to skip that check. M1b uses the same long settings (no reference fit to check).
#   M2c (tasks 8-15): --tag=<TAG> (so the init lookup finds the tagged M1c fit), sampler M2_ITER / M2_WARMUP / M2_ADAPT_DELTA, default empty =
#     the M2 defaults of spec.R (4000 / 2000 / 0.95).
# DRY RUN (open account, no sampling, no init fit needed):
#   sbatch --account=open --export=ALL,ROLES=bcj,DRY_RUN=TRUE --time=00:30:00 --mem=16G --cpus-per-task=1 v2/slurm/09_joiner_m2c.sbatch
# Post: Rscript v2/R/03_extract.R --fit=v2/output/fit_<cell>.rds for each fit.
# Needs the entry build (Rscript v2/R/01_build_data.R --entry=later -> v2/data_entry with its manifest) and the recompiled m1_context / m1_context_rho
# / m2_context programs (compile_models.R). NOTE: the edited m1_context / m1_context_rho sources change their hashes, so 03_extract.R refuses M1c /
# M1cr fits sampled with the earlier sources (assert_fit_stan_current).
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
i=${SLURM_ARRAY_TASK_ID:?run as an array}
vr=all
case $i in
   0) s=fmid; m=M1c;  K=;  js=all ;;
   1) s=mid;  m=M1c;  K=;  js=all ;;
   2) s=fmid; m=M1cr; K=;  js=all ;;
   3) s=mid;  m=M1cr; K=;  js=all ;;
   4) s=fmid; m=M1c;  K=;  js=opposing ;;
   5) s=mid;  m=M1c;  K=;  js=opposing ;;
   6) s=fmid; m=M1b;  K=;  js=all ;;
   7) s=mid;  m=M1b;  K=;  js=all ;;
   8) s=fmid; m=M2c;  K=1; js=all ;;
   9) s=mid;  m=M2c;  K=1; js=all ;;
  10) s=fmid; m=M2c;  K=2; js=all ;;
  11) s=mid;  m=M2c;  K=2; js=all ;;
  12) s=fmid; m=M2c;  K=1; js=opposing ;;
  13) s=mid;  m=M2c;  K=1; js=opposing ;;
  14) s=fmid; m=M2c;  K=2; js=opposing ;;
  15) s=mid;  m=M2c;  K=2; js=opposing ;;
  16) s=fmid; m=M1c;  K=;  js=all; vr=polrel ;;
  17) s=mid;  m=M1c;  K=;  js=all; vr=polrel ;;
  18) s=fmid; m=M2c;  K=1; js=all; vr=polrel ;;
  19) s=mid;  m=M2c;  K=1; js=all; vr=polrel ;;
  20) s=fmid; m=M2c;  K=2; js=all; vr=polrel ;;
  21) s=mid;  m=M2c;  K=2; js=all; vr=polrel ;;
  *) echo "array index $i out of range 0-21" >&2; exit 1 ;;
esac
ROLES=${ROLES:-bcj}      # PI decision 2026-11-04: bcj is the main coding (especially for M2c); ROLES=bc is the optional appendix coding
case $ROLES in bc|bcj) ;; *) echo "ROLES must be bc or bcj (got '$ROLES')" >&2; exit 1 ;; esac
DRY=${DRY_RUN:-FALSE}
THR=${THREADS_PER_CHAIN:-8}
TAG=${TAG:-long}
[[ "$TAG" =~ ^[A-Za-z0-9]+$ ]] || { echo "TAG must be letters and digits (got '$TAG')" >&2; exit 1; }
if [ "$m" = "M1b" ]; then ROLE_ARGS=""; else ROLE_ARGS="--roles=$ROLES --joiner_set=$js"; fi
BASE="--sample=$s --lag=3 --model=$m --comparison=directed --entry=later --variant=$vr $ROLE_ARGS --tag=$TAG --dry_run=$DRY --threads_per_chain=$THR"
if [ "$m" != "M2c" ]; then
  IT=${LONG_ITER:-6000}; WU=${LONG_WARMUP:-3000}; AD=${LONG_ADAPT_DELTA:-0.9}
  SAMP="--iter=$IT --warmup=$WU --adapt_delta=$AD"
  if [ "$DRY" != "TRUE" ] && [ "${CHECK_LONG:-TRUE}" = "TRUE" ] && [ "$m" != "M1b" ]; then
    ref="v2/output/fit_${s}_lag3_all_pr_${m}_cdir_aprime_t${TAG}.rds"
    [ -f "$ref" ] || { echo "ERROR: $ref (the existing long aprime fit whose settings this fit must match) is not in v2/output; copy it there or set CHECK_LONG=FALSE" >&2; exit 1; }
    Rscript -e "a <- readRDS('$ref')\$args; ok <- identical(as.integer(a\$iter), as.integer($IT)) && identical(as.integer(a\$warmup), as.integer($WU)) && isTRUE(all.equal(as.numeric(a\$adapt_delta), $AD)); cat(sprintf('existing long aprime fit %s: iter %s warmup %s adapt_delta %s; this task: $IT / $WU / $AD -> %s\n', '$ref', a\$iter, a\$warmup, a\$adapt_delta, if (ok) 'match' else 'MISMATCH')); if (!ok) quit(status = 1)" \
      || { echo "ERROR: LONG_ITER / LONG_WARMUP / LONG_ADAPT_DELTA differ from the existing long aprime fit; export the recorded values" >&2; exit 1; }
  fi
else
  SAMP="--ame_K=$K${M2_ITER:+ --iter=$M2_ITER}${M2_WARMUP:+ --warmup=$M2_WARMUP}${M2_ADAPT_DELTA:+ --adapt_delta=$M2_ADAPT_DELTA}"
  if [ "$DRY" != "TRUE" ]; then
    want="v2/output/fit_${s}_lag3_${vr}_pr_M1c_cdir_${ROLES}$([ "$js" = opposing ] && echo _jset)_ent_t${TAG}.rds"
    [ -f "$want" ] || { echo "ERROR: task $i ($s M2c K$K joiner_set=$js) needs $want (see the table above)" >&2; exit 1; }
  fi
fi
echo "task $i -> sample=$s model=$m roles=$ROLES (M1b: none) entry=later variant=$vr K=${K:--} joiner_set=$js tag=$TAG; $SAMP; threads_per_chain=$THR cpus=${SLURM_CPUS_PER_TASK:-?} dry_run=$DRY"
Rscript v2/R/02_fit_model.R $BASE $SAMP
