## ============================================================================
## 06_context_compare.R -- compare the bilateral and coalition context fits of one cell (dual typology stage a, 2026-10-25)
##
## Usage: Rscript v2/R/06_context_compare.R --sample=mid --lag=3 [--model=M1b] [--coal_ties=inclusive|strict] [other switches of spec.R]
##        [--r0=0.5] [--min_decisions=3] [--allow_stale=FALSE]
## Loads fit_<cell>_bi.rds and fit_<cell>_co.rds (the cell name of the sample, lag, model and the OTHER switches; --coal_ties=strict
## adds _strict to both) and, when present, fit_<cell>_co_drop3.rds (the coalition context refitted without the 3 disputes contributing
## the most coalition events). Refuses pilots and stale fits (unless --allow_stale), and fits of data that no longer match the manifest.
## Writes to v2/output/:
##   context_compare_<group>.csv         the long statistics table (block context / leverage / role; see context_compare.R)
##   context_compare_<group>_crosstab.csv MAP-type cross-tabulations
##   context_compare_<group>_states.csv  per-state weighted decisions, MAP types and mu means per context
##   context_compare_<group>.md          a short summary
## where <group> = the cell name without the context suffix: <sample>_lag<h>_all_pr_<model>[_K<K>][_orig|_join][_strict][_noww2]...
## (the suffixes in cell_name() order, so fmid_lag3_all_pr_M1b_strict_noww2 goes with the fits _co_strict_noww2 and _bi_strict_noww2).
## ============================================================================
suppressPackageStartupMessages({ library(rstan); library(dplyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
R_DIR <- dirname(this_file)
source(file.path(R_DIR, "spec.R")); source(file.path(R_DIR, "helpers_data.R")); source(file.path(R_DIR, "helpers_model.R"))
source(file.path(R_DIR, "typology.R")); source(file.path(R_DIR, "context_compare.R"))
Sys.setenv(V2_ROOT = normalizePath(file.path(R_DIR, "..")))
args <- commandArgs(trailingOnly = TRUE)
refuse_archived_options(args, "06_context_compare.R")    # audit fix 14
SAMPLE <- get_arg(args, "sample", "mid"); LAG <- as.integer(get_arg(args, "lag", 3)); MODEL <- get_arg(args, "model", "M1b")
R0 <- get_arg(args, "r0", R0); MIN_DEC <- get_arg(args, "min_decisions", 3); ALLOW_STALE <- get_arg(args, "allow_stale", FALSE)
stopifnot(SAMPLE %in% ALL_SAMPLES, LAG %in% LAGS, MODEL %in% MODELS, model_D_mode(MODEL) == 3L, is.finite(R0), R0 > 0, MIN_DEC > 0)
if (length(grep("^--(context|coal_drop)=", args))) stop("06_context_compare.R sets --context and --coal_drop itself; pass --coal_ties and the other switches only")
TIES <- get_arg(args, "coal_ties", "inclusive"); stopifnot(TIES %in% SPEC_CHOICES$coal_ties)
base_args <- args[!grepl("^--(sample|lag|model|r0|min_decisions|allow_stale|coal_ties)=", args)]
spec_bi <- parse_spec(c(base_args, "--context=bilateral", paste0("--coal_ties=", TIES)))
spec_co <- parse_spec(c(base_args, "--context=coalition", paste0("--coal_ties=", TIES)))
spec_d3 <- parse_spec(c(base_args, "--context=coalition", paste0("--coal_ties=", TIES), "--coal_drop=top3"))
## the output group name = the bilateral cell name without its context suffix, so the suffix order is cell_name()'s (pairs, [context],
## coal_ties, coal_drop, ww2, controls, weights, changes, direction; audit fix 23, 2026-10-26: before, _strict was appended last)
group <- sub(SPEC_SUFFIX$context[["bilateral"]], "", cell_name(SAMPLE, LAG, MODEL, spec_bi), fixed = TRUE)
stopifnot(identical(sub(SPEC_SUFFIX$context[["coalition"]], "", cell_name(SAMPLE, LAG, MODEL, spec_co), fixed = TRUE), group))
log_open(paste0("06_context_compare_", group))
out_dir <- file.path(v2_root(), "output"); data_dir <- spec_data_dir(spec_bi, v2_root()); manifest_path <- file.path(data_dir, "manifest.csv")
logf("group %s; ties %s; r0 %.2f; min weighted decisions %g; data %s", group, TIES, R0, MIN_DEC, basename(data_dir))
load_fit <- function(spec, required) {
  f <- file.path(out_dir, paste0("fit_", cell_name(SAMPLE, LAG, MODEL, spec), ".rds"))
  if (!file.exists(f) && ALLOW_STALE && file.exists(sub("\\.rds$", "_STALE.rds", f))) f <- sub("\\.rds$", "_STALE.rds", f)
  if (!file.exists(f)) { if (required) stop("missing fit: ", basename(f)); logf("optional fit missing: %s", basename(f)); return(NULL) }
  obj <- read_fit_or_null(f)
  if (is.null(obj)) { if (required) stop("unreadable fit: ", basename(f)); return(NULL) }
  assert_production_fit(obj, allow_stale = ALLOW_STALE); assert_fit_current(obj, manifest_path); assert_fit_stan_current(obj)   # audit fix 5
  stopifnot(identical(obj$spec$context, spec$context), identical(obj$spec$coal_ties, spec$coal_ties), identical(obj$spec$coal_drop, spec$coal_drop))
  logf("loaded %s: E = %d, N1 = %d, sum(y1) = %d, weights %.1f", basename(f), obj$dat$E, obj$dat$N1, sum(obj$dat$y1), sum(obj$dat$ev_w))
  obj
}
fits <- list(bilateral = load_fit(spec_bi, TRUE), coalition = load_fit(spec_co, TRUE), coalition_drop3 = load_fit(spec_d3, FALSE))
## the joiner flag per event comes from the build's events table (manifest-checked)
events <- assert_matches_manifest(file.path(data_dir, "events.rds"), manifest_path)
cc <- context_compare(fits, events |> filter(!escalation) |> select(initiator, target, year, is_joiner_pair), r0 = R0, min_decisions = MIN_DEC)
write.csv(cc$stats, file.path(out_dir, paste0("context_compare_", group, ".csv")), row.names = FALSE)
write.csv(cc$crosstab, file.path(out_dir, paste0("context_compare_", group, "_crosstab.csv")), row.names = FALSE)
write.csv(cc$states, file.path(out_dir, paste0("context_compare_", group, "_states.csv")), row.names = FALSE)
md <- context_compare_md(cc, sprintf("Context split, %s (%s ties)", group, TIES))
if (!is.null(fits$coalition_drop3) && !is.null(fits$coalition_drop3$args$context_top3))
  md <- c(md, "", sprintf("Disputes dropped in the leverage fit: %s.", paste(sprintf("%d (%d coalition events)", fits$coalition_drop3$args$context_top3$dispnum, fits$coalition_drop3$args$context_top3$n_events), collapse = ", ")))
writeLines(md, file.path(out_dir, paste0("context_compare_", group, ".md")))
print(cc$stats |> select(block, comparison, axis, statistic, estimate, lo, hi, n_states))
logf("06_context_compare done for %s", group)
