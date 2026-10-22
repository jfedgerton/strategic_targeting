## ============================================================================
## run_measurement_real.R -- the refit-based measurement tests on ONE real-data cell
## (temporal holdout, permutation test, systematic measurement error). LOO and
## misclassification between two saved fits are done by measurement_tests.R's CLI.
##
## Usage: Rscript v2/R/run_measurement_real.R --sample=mid --lag=5 --variant=polrel --family=m1 \
##          [--centrality=pr] [--K=2] [--cutoff=1990] [--n_perm=50] [--perm_mode=rescore|refit] \
##          [--deltas=0.5,1.0] [--r0=0.5] [--chains=4] [--iter=1500] [--warmup=750] [--sampled_m=0]
## Needs fit_a / fit_b for the cell from 02_fit_model.R (for the permutation rescore and the
## baseline classification). Writes v2/output/measurement_<cell>_*.csv.
## NOT run in this phase (no real-data model runs); exercised on simulated data by
## v2/sim/run_measurement_tests.R.
## ============================================================================
suppressPackageStartupMessages({ library(rstan); library(loo); library(dplyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(this_file), "00_helpers.R")); source(file.path(dirname(this_file), "typology.R"))
source(file.path(dirname(this_file), "measurement_tests.R")); source(file.path(dirname(this_file), "assemble_cell.R"))
Sys.setenv(V2_ROOT = normalizePath(file.path(dirname(this_file), "..")))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) { h <- grep(paste0("^--", name, "="), args, value = TRUE); if (!length(h)) return(default)
  v <- sub(paste0("^--", name, "="), "", h[1]); if (is.logical(default)) as.logical(v) else if (is.numeric(default)) as.numeric(v) else v }
SAMPLE <- get_arg("sample", "mid"); LAG <- as.integer(get_arg("lag", 5)); VARIANT <- get_arg("variant", "polrel"); CENT <- get_arg("centrality", "pr")
FAMILY <- get_arg("family", "m1"); K <- as.integer(get_arg("K", 2)); CUTOFF <- as.integer(get_arg("cutoff", 1990))
N_PERM <- as.integer(get_arg("n_perm", 50)); PERM_MODE <- get_arg("perm_mode", "rescore")
DELTAS <- as.numeric(strsplit(get_arg("deltas", "0.5,1.0"), ",")[[1]]); R0 <- get_arg("r0", 0.5)
CHAINS <- as.integer(get_arg("chains", 4)); ITER <- as.integer(get_arg("iter", 1500)); WARMUP <- as.integer(get_arg("warmup", 750))
SAMPLED_M <- as.integer(get_arg("sampled_m", 0))
cell <- sprintf("%s_lag%d_%s_%s_%s", SAMPLE, LAG, VARIANT, CENT, toupper(FAMILY))
log_open(paste0("measurement_real_", cell)); set.seed(V2_SEED)
rstan_options(auto_write = FALSE); options(mc.cores = min(CHAINS, parallel::detectCores()))
out_dir <- file.path(v2_root(), "output")
models <- list(m1 = readRDS(file.path(v2_root(), "stan", "build", "m1_latent_additive.rds")),
               m2 = readRDS(file.path(v2_root(), "stan", "build", "m2_ame.rds")))
A <- assemble_cell(SAMPLE, LAG, VARIANT, CENT, TRUE, SAMPLED_M)
## saved fits follow 02_fit_model.R naming: fit_<sample>_lag<h>_<variant>_<cent>_<M1a|M1b>[_K2].rds
fa_path <- file.path(out_dir, sprintf("fit_%s_lag%d_%s_%s_%sa%s.rds", SAMPLE, LAG, VARIANT, CENT, toupper(FAMILY), if (FAMILY == "m2") sprintf("_K%d", K) else ""))
fb_path <- sub("a(_K[0-9]+)?\\.rds$", "b\\1.rds", fa_path)
stopifnot(file.exists(fa_path), file.exists(fb_path))
fit_a <- readRDS(fa_path)$fit; fit_b <- readRDS(fb_path)$fit
## temporal holdout
th <- temporal_holdout(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, cutoff = CUTOFF, family = FAMILY, K = K,
                       stan_models = models, chains = CHAINS, iter = ITER, warmup = WARMUP, x_ref = A$x_ref)
write.csv(th$summary, file.path(out_dir, sprintf("measurement_%s_temporal_holdout_%d.csv", cell, CUTOFF)), row.names = FALSE)
## permutation
pt <- permutation_test(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, fit_a, fit_b, n_perm = N_PERM, mode = PERM_MODE,
                       family = FAMILY, K = K, stan_models = models, chains = CHAINS, iter = ITER, warmup = WARMUP, x_ref = A$x_ref)
write.csv(pt$summary, file.path(out_dir, sprintf("measurement_%s_permutation_%s.csv", cell, PERM_MODE)), row.names = FALSE)
## systematic measurement error: recent conflict = state initiated in year t-1 (positions are at t-1)
recent <- A$sy |> transmute(state_idx, year = year + 1L, recent = y1)
cls_b <- classify_types(mu_state_draws(fit_b), r0 = R0)
me <- measurement_error_check(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, recent = recent, deltas = DELTAS,
                              baseline_cls = cls_b, family = FAMILY, joint = 1L, K = K, stan_models = models,
                              chains = CHAINS, iter = ITER, warmup = WARMUP, r0 = R0, x_ref = A$x_ref, dR_in_X = "dR_i")
write.csv(me, file.path(out_dir, sprintf("measurement_%s_error_inflation.csv", cell)), row.names = FALSE)
logf("run_measurement_real done")
