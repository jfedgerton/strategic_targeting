#!/bin/bash
#SBATCH --job-name=v2_jset_m2c
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --requeue
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32   # = chains (4) x threads_per_chain (8)
#SBATCH --mem=128G                                 # basic bills 4 GB per core: 32 cores x 4 GB
#SBATCH --array=0-9
#SBATCH --output=%x_%A_%a.out
#SBATCH --open-mode=append
# JOINER CHOICE SETS AND THE AME ROLE MODEL M2c (PI-approved 2026-11-04). Directed comparison, aprime roles, lag 3, fmid and mid. 10 tasks:
#   idx  sample  model  K  joiner_set  needs (init fit in v2/output, M2c only)
#    0   fmid    M1c    -  opposing    -
#    1   mid     M1c    -  opposing    -
#    2   fmid    M2c    1  all         fit_fmid_lag3_all_pr_M1c_cdir_aprime_t<TAG>.rds        (the EXISTING long M1c fit; copy it into v2/output)
#    3   mid     M2c    1  all         fit_mid_lag3_all_pr_M1c_cdir_aprime_t<TAG>.rds         (existing)
#    4   fmid    M2c    2  all         fit_fmid_lag3_all_pr_M1c_cdir_aprime_t<TAG>.rds        (existing)
#    5   mid     M2c    2  all         fit_mid_lag3_all_pr_M1c_cdir_aprime_t<TAG>.rds         (existing)
#    6   fmid    M2c    1  opposing    fit_fmid_lag3_all_pr_M1c_cdir_aprime_jset_t<TAG>.rds   (written by task 0)
#    7   mid     M2c    1  opposing    fit_mid_lag3_all_pr_M1c_cdir_aprime_jset_t<TAG>.rds    (written by task 1)
#    8   fmid    M2c    2  opposing    fit_fmid_lag3_all_pr_M1c_cdir_aprime_jset_t<TAG>.rds   (task 0)
#    9   mid     M2c    2  opposing    fit_mid_lag3_all_pr_M1c_cdir_aprime_jset_t<TAG>.rds    (task 1)
# Dependencies: tasks 6-9 need tasks 0-1 (afterok on that array, kill-on-invalid-dep); tasks 2-5 need no new job, only the existing long M1c
# fits (_cdir_aprime_t<TAG>) in v2/output. 02_fit_model.R REFUSES to sample an M2c cell whose M1c fit is missing (no random-init fallback), and
# refuses an M1c fit of another joiner_set / tag / cell. Submit with submit_joiner_m2c.sh (three sbatch calls; nothing here submits by itself).
# Sampler settings:
#   M1c (tasks 0-1): the long role-fit settings, LONG_ITER / LONG_WARMUP / LONG_ADAPT_DELTA / TAG (defaults 6000 / 3000 / 0.9 / long). These defaults
#     are the README's documented long-refit example (--iter=6000 --warmup=3000 --tag=long), NOT read from the existing M1c long fits, which are not
#     in the repository: before sampling, each M1c task reads the existing fit_<s>_lag3_all_pr_M1c_cdir_aprime_t<TAG>.rds and STOPS unless its
#     recorded iter / warmup / adapt_delta equal these (so the jset fit is comparable to it); set CHECK_LONG=FALSE to skip that check.
#   M2c (tasks 2-9): --tag=<TAG> (so the init lookup finds the tagged M1c fit; the cell names carry _t<TAG>), sampler M2_ITER / M2_WARMUP /
#     M2_ADAPT_DELTA, default empty = the M2 defaults of spec.R (4000 / 2000 / 0.95).
# DRY RUN (open account, no sampling, no init fit needed):
#   sbatch --account=open --export=ALL,DRY_RUN=TRUE --time=00:30:00 --mem=16G --cpus-per-task=1 v2/slurm/09_joiner_m2c.sbatch
# Post: Rscript v2/R/03_extract.R --fit=v2/output/fit_<cell>.rds for each of the 10 fits.
# Needs the data built at or after the 2026-11-04 commit (joiner_sets_aprime.rds in v2/data/manifest.csv) and the recompiled m1_context /
# m1_context_rho / m2_context programs (compile_models.R). NOTE: the edited m1_context source changes its hash, so 03_extract.R refuses M1c fits
# sampled with the earlier source (assert_fit_stan_current); the M2c init from such a fit is allowed (the init path checks data and cell, not the hash).
source "${SLURM_SUBMIT_DIR:?submit from the repo root}/v2/slurm/00_env.sh" || exit 1
cd "$V2_ROOT/.." || exit 1
i=${SLURM_ARRAY_TASK_ID:?run as an array}
case $i in
  0) s=fmid; m=M1c; K=;  js=opposing ;;
  1) s=mid;  m=M1c; K=;  js=opposing ;;
  2) s=fmid; m=M2c; K=1; js=all ;;
  3) s=mid;  m=M2c; K=1; js=all ;;
  4) s=fmid; m=M2c; K=2; js=all ;;
  5) s=mid;  m=M2c; K=2; js=all ;;
  6) s=fmid; m=M2c; K=1; js=opposing ;;
  7) s=mid;  m=M2c; K=1; js=opposing ;;
  8) s=fmid; m=M2c; K=2; js=opposing ;;
  9) s=mid;  m=M2c; K=2; js=opposing ;;
  *) echo "array index $i out of range 0-9" >&2; exit 1 ;;
esac
DRY=${DRY_RUN:-FALSE}
THR=${THREADS_PER_CHAIN:-8}
TAG=${TAG:-long}
[[ "$TAG" =~ ^[A-Za-z0-9]+$ ]] || { echo "TAG must be letters and digits (got '$TAG')" >&2; exit 1; }
BASE="--sample=$s --lag=3 --model=$m --comparison=directed --roles=aprime --joiner_set=$js --tag=$TAG --dry_run=$DRY --threads_per_chain=$THR"
if [ "$m" = "M1c" ]; then
  IT=${LONG_ITER:-6000}; WU=${LONG_WARMUP:-3000}; AD=${LONG_ADAPT_DELTA:-0.9}
  SAMP="--iter=$IT --warmup=$WU --adapt_delta=$AD"
  if [ "$DRY" != "TRUE" ] && [ "${CHECK_LONG:-TRUE}" = "TRUE" ]; then
    ref="v2/output/fit_${s}_lag3_all_pr_M1c_cdir_aprime_t${TAG}.rds"
    [ -f "$ref" ] || { echo "ERROR: $ref (the existing long M1c fit whose settings this jset fit must match) is not in v2/output; copy it there or set CHECK_LONG=FALSE" >&2; exit 1; }
    Rscript -e "a <- readRDS('$ref')\$args; ok <- identical(as.integer(a\$iter), as.integer($IT)) && identical(as.integer(a\$warmup), as.integer($WU)) && isTRUE(all.equal(as.numeric(a\$adapt_delta), $AD)); cat(sprintf('existing long M1c fit %s: iter %s warmup %s adapt_delta %s; this task: $IT / $WU / $AD -> %s\n', '$ref', a\$iter, a\$warmup, a\$adapt_delta, if (ok) 'match' else 'MISMATCH')); if (!ok) quit(status = 1)" \
      || { echo "ERROR: LONG_ITER / LONG_WARMUP / LONG_ADAPT_DELTA differ from the existing long M1c fit; export the recorded values" >&2; exit 1; }
  fi
else
  SAMP="--ame_K=$K${M2_ITER:+ --iter=$M2_ITER}${M2_WARMUP:+ --warmup=$M2_WARMUP}${M2_ADAPT_DELTA:+ --adapt_delta=$M2_ADAPT_DELTA}"
  if [ "$DRY" != "TRUE" ]; then
    want="v2/output/fit_${s}_lag3_all_pr_M1c_cdir_aprime$([ "$js" = opposing ] && echo _jset)_t${TAG}.rds"
    [ -f "$want" ] || { echo "ERROR: task $i ($s M2c K$K joiner_set=$js) needs $want (see the table above)" >&2; exit 1; }
  fi
fi
echo "task $i -> sample=$s model=$m K=${K:--} joiner_set=$js tag=$TAG; $SAMP; threads_per_chain=$THR cpus=${SLURM_CPUS_PER_TASK:-?} dry_run=$DRY"
Rscript v2/R/02_fit_model.R $BASE $SAMP
