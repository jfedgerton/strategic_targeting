## ============================================================================
## run_measurement_tests.R -- exercise measurement_tests.R and benchmarks.R on SIMULATED data
##
## Usage: Rscript v2/sim/run_measurement_tests.R [chains] [iter] [warmup]
## Uses a smaller simulated world (S = 30, Tt = 12, truth = M1b; simulate_current(), the CURRENT model layout since 2026-10-25, B5)
## so that the refits fit in a 2-core container. Benchmarks: run_benchmarks_sim.R. Results -> v2/sim/results/measurement_*.csv ;
## summarised in SIM_REPORT.md. The holdout freezes its preprocessing constants on the training years (B9), the permutation test runs
## in rescore mode (a reliance diagnostic; refit mode gives the calibrated p-value), and the measurement-error check perturbs the
## relational LEVELS and rebuilds the gaps, changes and level column (B8).
## ============================================================================
suppressPackageStartupMessages({ library(rstan); library(loo); library(dplyr); library(brms) })
args <- commandArgs(trailingOnly = TRUE)
CHAINS <- if (length(args) >= 1) as.integer(args[1]) else 2L
ITER   <- if (length(args) >= 2) as.integer(args[2]) else 500L
WARMUP <- if (length(args) >= 3) as.integer(args[3]) else 250L
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
sim_dir <- normalizePath(dirname(this_file)); Sys.setenv(V2_ROOT = normalizePath(file.path(sim_dir, "..")))
source(file.path(sim_dir, "..", "R", "00_helpers.R")); source(file.path(sim_dir, "..", "R", "typology.R"))
source(file.path(sim_dir, "..", "R", "measurement_tests.R")); source(file.path(sim_dir, "..", "R", "benchmarks.R"))
source(file.path(sim_dir, "simulate.R"))
log_open("run_measurement_tests"); set.seed(V2_SEED)
rstan_options(auto_write = FALSE); options(mc.cores = min(CHAINS, parallel::detectCores()))
res_dir <- file.path(sim_dir, "results"); dir.create(res_dir, showWarnings = FALSE)
models <- list(m1 = readRDS(file.path(v2_root(), "stan", "build", "m1_latent_additive.rds")),
               m2 = readRDS(file.path(v2_root(), "stan", "build", "m2_ame.rds")))
R0 <- 0.5
sim <- simulate_current("M1b", S = 30L, Tt = 12L, seed = V2_SEED)
W_cols <- sim$cols$W; X_cols <- sim$cols$X; Z_cols <- sim$cols$Z; x_ref <- sim$x_ref
logf("small world: S=%d N1=%d N2=%d E=%d (weight-1/2 events %d; IV-only rows %d, stage-2-only rows %d)", sim$meta$S, sim$meta$N1, sim$meta$N2, sim$meta$E, 2 * sim$meta$n_coalition, sim$meta$n_iv_only, sim$meta$n_s2_only)
sy <- sim$sy; alts <- sim$alts; ev <- sim$ev
excl <- c("z_re", "z_b", "V_raw", "U_free", "L_Om_joint", "L_Om", "re")

## ---- baseline fits a and b (M1) ------------------------------------------------------------------------------
fit_ab <- function(joint) {
  dat <- build_stan_data(sy, alts, ev, W_cols, X_cols, Z_cols, joint = joint, x_ref = x_ref)
  list(fit = sampling(models$m1, data = dat, chains = CHAINS, iter = ITER, warmup = WARMUP, seed = V2_SEED, refresh = 0,
                      control = list(adapt_delta = 0.9), pars = excl, include = FALSE), dat = dat)
}
t0 <- Sys.time(); fa <- fit_ab(0L); fb <- fit_ab(1L); fit_a <- fa$fit; fit_b <- fb$fit
logf("baseline a and b fitted in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

## ---- 1. LOO a vs b (stage 1 and 2 separately) -----------------------------------------------------------------
lc <- loo_compare_ab(fit_a, fit_b, fa$dat, fb$dat); print(lc)                   # stage 2, stage 1 conditional, joint state-year (B1)
write.csv(lc, file.path(res_dir, "measurement_loo_ab.csv"), row.names = FALSE)

## ---- 2. misclassification a vs b ------------------------------------------------------------------------------
mc <- misclassification_ab(fit_a, fit_b, r0 = R0)
cls_a <- classify_types(mu_state_draws(fit_a), r0 = R0); cls_b <- classify_types(mu_state_draws(fit_b), r0 = R0)
mis <- data.frame(rate_a_vs_b = mc$rate, mean_abs_prob_shift = mc$mean_abs_prob_shift,
                  accuracy_a = mean(cls_a$map_type == sim$types), accuracy_b = mean(cls_b$map_type == sim$types))
print(mis); write.csv(mis, file.path(res_dir, "measurement_misclassification.csv"), row.names = FALSE)
write.csv(as.data.frame(mc$table), file.path(res_dir, "measurement_misclassification_table.csv"), row.names = FALSE)

## ---- 3. temporal holdout --------------------------------------------------------------------------------------
cutoff <- sort(unique(sy$year))[9]              # train on the first 9 usable years, test on the last 3
t0 <- Sys.time()
freeze <- list(positions = sim$levels, gap_raw = c(dM = "dM_logratio", dR = "dR_pr_logratio"), levels = list(dM = "M_lvl", dR = "R_lvl"),
               change_raw = c(dM_i = "dM_raw_h", dR_i = "dRpr_raw_h"), target_change = c(dM_i = "dM_j_h", dR_i = "dRpr_j_h"),
               ctrl_cent = c("ongoing_tm1", "log_prior_mids", "dem_tm1", "log_n_system_tm1"), z_cent = c("log_dy_prior_mids", "dM_j_h", "dRpr_j_h"))
th <- temporal_holdout(sy, alts, ev, W_cols, X_cols, Z_cols, cutoff = cutoff, family = "m1", stan_models = models,
                       chains = CHAINS, iter = ITER, warmup = WARMUP, n_draws = 100, x_ref = x_ref, freeze = freeze)
logf("temporal holdout fitted in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
print(th$summary); write.csv(th$summary, file.path(res_dir, "measurement_temporal_holdout.csv"), row.names = FALSE)

## ---- 4. permutation test (rescore mode; refit mode is available but too slow here) -----------------------------
pt <- permutation_test(sy, alts, ev, W_cols, X_cols, Z_cols, fit_a, fit_b, n_perm = 20, mode = "rescore", n_draws = 100, x_ref = x_ref)
print(pt$summary); write.csv(pt$summary, file.path(res_dir, "measurement_permutation.csv"), row.names = FALSE)
write.csv(data.frame(perm = seq_along(pt$permuted_advantage), adv = pt$permuted_advantage, a = pt$permuted_a, b = pt$permuted_b),
          file.path(res_dir, "measurement_permutation_draws.csv"), row.names = FALSE)

## ---- 5. systematic measurement error ---------------------------------------------------------------------------
## "recent conflict" at the POSITION-year: the state initiated in that year, so its measured relational level that year is inflated
## (deltas in log points of R, like the real-data run); the gaps, changes and level column are rebuilt from the perturbed levels (B8)
recent <- sy |> transmute(state_idx, year, recent = y1)
t0 <- Sys.time()
me <- measurement_error_check(sy, alts, ev, W_cols, X_cols, Z_cols, recent = recent, deltas = c(0.25, 0.5, 1.0),
                              baseline_cls = cls_b, levels = sim$levels, h = 1L, gap_scale_R = sim$truth$gap_R,
                              family = "m1", joint = 1L, stan_models = models,
                              chains = CHAINS, iter = ITER, warmup = WARMUP, r0 = R0, x_ref = x_ref)
logf("measurement-error refits in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
me$baseline_accuracy <- mean(cls_b$map_type == sim$types)
print(me); write.csv(me, file.path(res_dir, "measurement_error_inflation.csv"), row.names = FALSE)

## Benchmarks on the simulated frames are in run_benchmarks_sim.R (separate R session: the brms compile
## needs ~4 GB and was OOM-killed here alongside the six fits).
logf("run_measurement_tests finished")
