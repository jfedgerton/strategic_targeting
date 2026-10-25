## ============================================================================
## run_measurement_real.R -- the refit-based measurement tests on ONE real-data cell
## (temporal holdout, permutation test, systematic measurement error). LOO and
## misclassification between two saved fits are done by measurement_tests.R's CLI.
##
## Usage: Rscript v2/R/run_measurement_real.R --sample=mid --lag=5 --variant=all --family=m1 \
##          [--tests=holdout,permutation,inflation] [--centrality=pr] [--K=2] [--cutoff=1990] [--n_perm=50] \
##          [--perm_mode=rescore|refit] [--deltas=0.25,0.5,1.0] [--r0=0.5] [--chains=4] [--iter=1500] [--warmup=750] \
##          [--event_weights=dispute|none] [--direction=cowwar|sideA] [--allow_stale=FALSE]
##   tests     : which tests to run (any subset; 06_measurement.sbatch runs permutation + inflation in one task per family
##               and the holdout in one task per cutoff and family)
##   deltas    : measurement-error inflation of the relational LEVEL (log PageRank) of state-years with recent conflict, in log
##               points; the gaps, the initiator and target changes and the level column are rebuilt from the perturbed levels (B8)
##   --sampled_m is NOT supported here and is rejected (C5): the holdout / permutation / inflation refits would need the sampled
##               sets re-drawn consistently for every refit, which is not implemented
## Needs the a / b fits of the cell from 02_fit_model.R (production fits of the current data that score exactly the assembled
## cell; checked with assert_production_fit(), assert_fit_current() and assert_fit_matches_cell(), C1 / C4).
## Writes v2/output/measurement_<cell>_*.csv.
## ============================================================================
suppressPackageStartupMessages({ library(rstan); library(loo); library(dplyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(this_file), "00_helpers.R")); source(file.path(dirname(this_file), "typology.R"))
source(file.path(dirname(this_file), "measurement_tests.R")); source(file.path(dirname(this_file), "assemble_cell.R"))
Sys.setenv(V2_ROOT = normalizePath(file.path(dirname(this_file), "..")))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) { h <- grep(paste0("^--", name, "="), args, value = TRUE); if (!length(h)) return(default)
  v <- sub(paste0("^--", name, "="), "", h[1]); if (is.logical(default)) as.logical(v) else if (is.numeric(default)) as.numeric(v) else v }
SAMPLE <- get_arg("sample", "mid"); LAG <- as.integer(get_arg("lag", 5)); VARIANT <- get_arg("variant", "all"); CENT <- get_arg("centrality", "pr")
FAMILY <- get_arg("family", "m1"); K <- as.integer(get_arg("K", 2)); CUTOFF <- as.integer(get_arg("cutoff", 1990))
N_PERM <- as.integer(get_arg("n_perm", 50)); PERM_MODE <- get_arg("perm_mode", "rescore")
DELTAS <- as.numeric(strsplit(get_arg("deltas", "0.25,0.5,1.0"), ",")[[1]]); R0 <- get_arg("r0", 0.5)
CHAINS <- as.integer(get_arg("chains", 4)); ITER <- as.integer(get_arg("iter", 1500)); WARMUP <- as.integer(get_arg("warmup", 750))
TESTS <- strsplit(get_arg("tests", "holdout,permutation,inflation"), ",")[[1]]
EVENT_WEIGHTS <- get_arg("event_weights", "dispute"); DIRECTION <- get_arg("direction", "cowwar"); ALLOW_STALE <- get_arg("allow_stale", FALSE)
## C5: --sampled_m is rejected explicitly instead of being propagated half-way
if (length(grep("^--sampled_m=", args)) > 0 && as.integer(sub("^--sampled_m=", "", grep("^--sampled_m=", args, value = TRUE)[1])) != 0L)
  stop("--sampled_m is not supported by run_measurement_real.R: the refit-based tests need full choice sets (the sampled IV rows and per-event sets would have to be re-drawn for every holdout / permutation / inflation refit, which is not implemented)")
stopifnot(all(TESTS %in% c("holdout", "permutation", "inflation")), FAMILY %in% c("m1", "m2"), PERM_MODE %in% c("rescore", "refit"),
          EVENT_WEIGHTS %in% c("dispute", "none"), DIRECTION %in% c("cowwar", "sideA"))
cell <- sprintf("%s_lag%d_%s_%s_%s%s%s", SAMPLE, LAG, VARIANT, CENT, toupper(FAMILY), if (EVENT_WEIGHTS == "none") "_w1" else "", if (DIRECTION == "sideA") "_sideA" else "")
log_open(paste0("measurement_real_", cell)); set.seed(V2_SEED)
rstan_options(auto_write = FALSE); options(mc.cores = min(CHAINS, parallel::detectCores()))
out_dir <- file.path(v2_root(), "output"); data_dir <- file.path(v2_root(), if (DIRECTION == "sideA") "data_sideA" else "data")
logf("cell %s; tests %s; data %s", cell, paste(TESTS, collapse = ","), basename(data_dir))
## compiled models, checked against the checked-out sources (C3)
models <- list(m1 = readRDS(file.path(v2_root(), "stan", "build", "m1_latent_additive.rds")),
               m2 = readRDS(file.path(v2_root(), "stan", "build", "m2_ame.rds")))
for (nm in names(models)) assert_compiled_current(models[[nm]], file.path(v2_root(), "stan", paste0(c(m1 = "m1_latent_additive", m2 = "m2_ame")[nm], ".stan")),
                                                  file.path(v2_root(), "stan", "build", paste0(c(m1 = "m1_latent_additive", m2 = "m2_ame")[nm], ".sha256")))
A <- assemble_cell(SAMPLE, LAG, VARIANT, CENT, TRUE, 0L, EVENT_WEIGHTS = EVENT_WEIGHTS, data_dir = data_dir)
## saved fits follow 02_fit_model.R naming: fit_<sample>_lag<h>_<variant>_<cent>_<M1a|M1b>[_K2][_w1][_sideA].rds
suffix <- paste0(if (FAMILY == "m2") sprintf("_K%d", K) else "", if (EVENT_WEIGHTS == "none") "_w1" else "", if (DIRECTION == "sideA") "_sideA" else "")
fa_path <- file.path(out_dir, sprintf("fit_%s_lag%d_%s_%s_%sa%s.rds", SAMPLE, LAG, VARIANT, CENT, toupper(FAMILY), suffix))
fb_path <- file.path(out_dir, sprintf("fit_%s_lag%d_%s_%s_%sb%s.rds", SAMPLE, LAG, VARIANT, CENT, toupper(FAMILY), suffix))
if (ALLOW_STALE) { if (!file.exists(fa_path)) fa_path <- sub("\\.rds$", "_STALE.rds", fa_path); if (!file.exists(fb_path)) fb_path <- sub("\\.rds$", "_STALE.rds", fb_path) }
stopifnot(file.exists(fa_path), file.exists(fb_path))
obj_a <- readRDS(fa_path); obj_b <- readRDS(fb_path)
## C1 / C2 / C4: production fits of the CURRENT data that score exactly this assembled cell
dat_chk <- build_stan_data(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 0L, x_ref = A$x_ref, K = if (FAMILY == "m2") K else NULL, allow_orphan_initiations = TRUE)
for (o in list(obj_a, obj_b)) {
  assert_production_fit(o, allow_stale = ALLOW_STALE)
  assert_fit_current(o, file.path(data_dir, "manifest.csv"))
  assert_fit_matches_cell(o, dat_chk, A$W_cols, A$X_cols, A$Z_cols)
  stopifnot(o$dat$D_mode == 3L)
}
stopifnot(obj_a$dat$joint == 0L, obj_b$dat$joint == 1L)
logf("fits %s and %s verified: production, current data, same E = %d, same events, same W / X / Z", basename(fa_path), basename(fb_path), dat_chk$E)
fit_a <- obj_a$fit; fit_b <- obj_b$fit
## positions (relational LEVELS) and the gap constant for the holdout freezing (B9) and the level perturbation (B8)
pos <- assert_matches_manifest(file.path(data_dir, "sy_positions.rds"), file.path(data_dir, "manifest.csv"))
levels <- pos |> transmute(state_idx = match(ccode, A$state_codes), year = as.integer(year), M, R = R_pr) |> filter(!is.na(state_idx))
scale_tbl <- read.csv(file.path(data_dir, "change_scale.csv"), stringsAsFactors = FALSE)
GAP_R <- scale_tbl$pooled_sd[is.na(scale_tbl$h) & scale_tbl$variable == "gap_R_pr"]; GAP_M <- scale_tbl$pooled_sd[is.na(scale_tbl$h) & scale_tbl$variable == "gap_M"]
stopifnot(length(GAP_R) == 1, length(GAP_M) == 1, is.finite(GAP_R), is.finite(GAP_M))
## ---- temporal holdout (constants frozen on the training years, B9) ---------------------------------------------------------------
if ("holdout" %in% TESTS) {
  freeze <- list(positions = levels, gap_raw = c(dM = "dM_logratio", dR = "dR_pr_logratio"), levels = list(dM = "M_lvl", dR = "R_lvl"),
                 change_raw = c(dM_i = "dM_raw_h", dR_i = "dRpr_raw_h"), target_change = c(dM_i = "dM_j_h", dR_i = "dRpr_j_h"),
                 ctrl_cent = setdiff(union(w_controls(), x_controls()), W_LEVELS), z_cent = intersect(A$Z_cols, Z_CENTRED))
  th <- temporal_holdout(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, cutoff = CUTOFF, family = FAMILY, K = K,
                         stan_models = models, chains = CHAINS, iter = ITER, warmup = WARMUP, x_ref = A$x_ref, freeze = freeze)
  write.csv(th$summary, file.path(out_dir, sprintf("measurement_%s_temporal_holdout_%d.csv", cell, CUTOFF)), row.names = FALSE)
  logf("temporal holdout written (cutoff %d)", CUTOFF)
}
## ---- permutation (B7: refit mode = calibrated p-value; rescore mode = reliance diagnostic) ----------------------------------------
if ("permutation" %in% TESTS) {
  pt <- permutation_test(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, fit_a, fit_b, n_perm = N_PERM, mode = PERM_MODE,
                         family = FAMILY, K = K, stan_models = models, chains = CHAINS, iter = ITER, warmup = WARMUP, x_ref = A$x_ref)
  write.csv(pt$summary, file.path(out_dir, sprintf("measurement_%s_permutation_%s.csv", cell, PERM_MODE)), row.names = FALSE)
  write.csv(data.frame(perm = seq_along(pt$permuted_advantage), adv = pt$permuted_advantage, a = pt$permuted_a, b = pt$permuted_b),
            file.path(out_dir, sprintf("measurement_%s_permutation_%s_draws.csv", cell, PERM_MODE)), row.names = FALSE)
  logf("permutation (%s) written: %s", PERM_MODE, pt$summary$interpretation)
}
## ---- systematic measurement error (B8): inflate the relational LEVEL of state-years with recent conflict ---------------------------
## recent conflict at the POSITION-year: the state initiated in that year (its measured centrality that year is inflated)
if ("inflation" %in% TESTS) {
  recent <- A$sy |> transmute(state_idx, year, recent = y1)
  cls_b <- classify_types(mu_state_draws(fit_b), r0 = R0)
  me <- measurement_error_check(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, recent = recent, deltas = DELTAS,
                                baseline_cls = cls_b, levels = levels, h = LAG, gap_scale_R = GAP_R, family = FAMILY, joint = 1L, K = K, stan_models = models,
                                chains = CHAINS, iter = ITER, warmup = WARMUP, r0 = R0, x_ref = A$x_ref)
  write.csv(me, file.path(out_dir, sprintf("measurement_%s_error_inflation.csv", cell)), row.names = FALSE)
  logf("measurement-error inflation written (deltas %s log points of R)", paste(DELTAS, collapse = ","))
}
logf("run_measurement_real done")
