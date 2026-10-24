## ============================================================================
## compare_models.R -- cross-model comparison for one (sample, lag) cell group of the main grid
##
## Usage: Rscript v2/R/compare_models.R --sample=mid --lag=5 [--variant=all] [--centrality=pr] [--K=2] [--r0=0.5]
## Loads the eight fits fit_<sample>_lag<h>_<variant>_<cent>_{M1m,M1r,M1a,M1b,M2m,M2r,M2a,M2b}[_K<K>].rds
## (missing fits are skipped and listed) and writes to v2/output/:
##   compare_loo_<group>.csv        per model: LOO ELPD (stage 2 choice, stage 1 initiation), SE, Pareto-k > 0.7 counts
##   compare_contrasts_<group>.csv  planned contrasts with elpd_diff and se_diff from loo_compare:
##                                    within each family: r vs m, a vs m, a vs r, b vs a (material vs relational vs
##                                    independent vs joint); across families: M2 vs M1 for each of m, r, a, b
##   compare_misclass_<group>.csv   a vs b MAP-type disagreement (initiating states), within M1 and within M2
## All models score the SAME events and choice sets (same assemble_cell() output), so LOO is like with like.
## ============================================================================
suppressPackageStartupMessages({ library(rstan); library(loo); library(dplyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(this_file), "00_helpers.R")); source(file.path(dirname(this_file), "typology.R"))
source(file.path(dirname(this_file), "measurement_tests.R"))
Sys.setenv(V2_ROOT = normalizePath(file.path(dirname(this_file), "..")))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) { h <- grep(paste0("^--", name, "="), args, value = TRUE); if (!length(h)) return(default)
  v <- sub(paste0("^--", name, "="), "", h[1]); if (is.numeric(default)) as.numeric(v) else v }
SAMPLE <- get_arg("sample", "mid"); LAG <- as.integer(get_arg("lag", 5)); VARIANT <- get_arg("variant", "all")
CENT <- get_arg("centrality", "pr"); K <- as.integer(get_arg("K", 2)); R0 <- get_arg("r0", 0.5)
group <- sprintf("%s_lag%d_%s_%s", SAMPLE, LAG, VARIANT, CENT)
log_open(paste0("compare_", group))
out_dir <- file.path(v2_root(), "output")
MODELS <- c("M1m", "M1r", "M1a", "M1b", "M2m", "M2r", "M2a", "M2b")
fits <- list(); loos <- list(); rows <- list()
for (m in MODELS) {
  f <- file.path(out_dir, sprintf("fit_%s_%s%s.rds", group, m, if (grepl("^M2", m)) sprintf("_K%d", K) else ""))
  if (!file.exists(f)) { logf("missing fit: %s", basename(f)); next }
  obj <- readRDS(f); assert_fit_current(obj, file.path(v2_root(), "data", "manifest.csv")); fits[[m]] <- obj
  l2 <- loo_from_fit(obj$fit, "log_lik_stage2"); l1 <- loo_from_fit(obj$fit, "log_lik_stage1")
  loos[[m]] <- list(l2 = l2, l1 = l1)
  rows[[m]] <- data.frame(model = m, D_mode = obj$dat$D_mode, joint = obj$dat$joint, E = obj$dat$E, N1 = obj$dat$N1,
                          elpd_stage2 = l2$estimates["elpd_loo", 1], se_stage2 = l2$estimates["elpd_loo", 2], bad_k_stage2 = sum(l2$diagnostics$pareto_k > 0.7),
                          elpd_stage1 = l1$estimates["elpd_loo", 1], se_stage1 = l1$estimates["elpd_loo", 2], bad_k_stage1 = sum(l1$diagnostics$pareto_k > 0.7))
  logf("%s: elpd2 = %.1f (%.1f), elpd1 = %.1f (%.1f)", m, rows[[m]]$elpd_stage2, rows[[m]]$se_stage2, rows[[m]]$elpd_stage1, rows[[m]]$se_stage1)
}
stopifnot(length(fits) > 0)
## every loaded fit must score the same events (like with like)
Es <- vapply(fits, function(o) o$dat$E, numeric(1)); stopifnot(length(unique(Es)) == 1)
evs <- lapply(fits, function(o) paste(o$dat$ev_sy, o$dat$alt_target[o$dat$ev_row])); stopifnot(all(vapply(evs, identical, logical(1), evs[[1]])))
tab <- bind_rows(rows) |> mutate(delta_elpd_stage2_vs_best = elpd_stage2 - max(elpd_stage2), delta_elpd_stage1_vs_best = elpd_stage1 - max(elpd_stage1))
write.csv(tab, file.path(out_dir, paste0("compare_loo_", group, ".csv")), row.names = FALSE)
## planned contrasts
contrast <- function(a, b, label) {
  if (!all(c(a, b) %in% names(loos))) return(NULL)
  c2 <- loo::loo_compare(loos[[a]]$l2, loos[[b]]$l2); c1 <- loo::loo_compare(loos[[a]]$l1, loos[[b]]$l1)
  data.frame(contrast = label, model_a = a, model_b = b,
             stage2_elpd_diff_b_minus_a = loos[[b]]$l2$estimates["elpd_loo", 1] - loos[[a]]$l2$estimates["elpd_loo", 1], stage2_se_diff = c2[2, "se_diff"],
             stage1_elpd_diff_b_minus_a = loos[[b]]$l1$estimates["elpd_loo", 1] - loos[[a]]$l1$estimates["elpd_loo", 1], stage1_se_diff = c1[2, "se_diff"])
}
cons <- list()
for (fam in c("M1", "M2")) {
  cons[[length(cons) + 1]] <- contrast(paste0(fam, "m"), paste0(fam, "r"), paste0(fam, ": relational vs material"))
  cons[[length(cons) + 1]] <- contrast(paste0(fam, "m"), paste0(fam, "a"), paste0(fam, ": independent (both) vs material"))
  cons[[length(cons) + 1]] <- contrast(paste0(fam, "r"), paste0(fam, "a"), paste0(fam, ": independent (both) vs relational"))
  cons[[length(cons) + 1]] <- contrast(paste0(fam, "a"), paste0(fam, "b"), paste0(fam, ": joint vs independent"))
}
for (sp in c("m", "r", "a", "b")) cons[[length(cons) + 1]] <- contrast(paste0("M1", sp), paste0("M2", sp), paste0(sp, ": M2 (AME) vs M1"))
con_tab <- bind_rows(cons); print(con_tab)
write.csv(con_tab, file.path(out_dir, paste0("compare_contrasts_", group, ".csv")), row.names = FALSE)
## misclassification a vs b within each family
mis <- list()
for (fam in c("M1", "M2")) {
  a <- paste0(fam, "a"); b <- paste0(fam, "b")
  if (!all(c(a, b) %in% names(fits))) next
  init_states <- which(tabulate(fits[[a]]$dat$sy_state[fits[[a]]$dat$ev_sy], nbins = fits[[a]]$dat$S) > 0)
  mc <- misclassification_ab(fits[[a]]$fit, fits[[b]]$fit, r0 = R0, only_states = init_states)
  mis[[fam]] <- data.frame(family = fam, r0 = R0, n_initiating_states = length(init_states), disagreement_rate = mc$rate, mean_abs_prob_shift = mc$mean_abs_prob_shift)
}
if (length(mis)) write.csv(bind_rows(mis), file.path(out_dir, paste0("compare_misclass_", group, ".csv")), row.names = FALSE)
logf("compare_models done for %s", group)
