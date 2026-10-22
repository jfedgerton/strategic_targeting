## ============================================================================
## run_measurement_tests.R -- exercise measurement_tests.R and benchmarks.R on SIMULATED data
##
## Usage: Rscript v2/sim/run_measurement_tests.R [chains] [iter] [warmup]
## Uses a smaller simulated world (S = 30, Tt = 12, truth = M1b) so that the six refits fit in
## a 2-core container. Results -> v2/sim/results/measurement_*.csv ; summarised in SIM_REPORT.md.
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
W_cols <- c("W1", "W2"); X_cols <- paste0("X", 1:4); Z_cols <- c("Z1", "Z2"); R0 <- 0.8

sim <- simulate_targeting("M1b", S = 30L, Tt = 12L, seed = V2_SEED)
logf("small world: S=%d N1=%d N2=%d E=%d", sim$meta$S, sim$meta$N1, sim$meta$N2, sim$meta$E)
sy <- sim$sy; alts <- sim$alts; ev <- sim$ev
excl <- c("z_re", "z_b", "V_raw", "U_free", "L_Om_joint", "L_Om", "re")

## ---- baseline fits a and b (M1) ------------------------------------------------------------------------------
fit_ab <- function(joint) {
  dat <- build_stan_data(sy, alts, ev, W_cols, X_cols, Z_cols, joint = joint, x_ref = rep(0, 4))
  sampling(models$m1, data = dat, chains = CHAINS, iter = ITER, warmup = WARMUP, seed = V2_SEED, refresh = 0,
           control = list(adapt_delta = 0.9), pars = excl, include = FALSE)
}
t0 <- Sys.time(); fit_a <- fit_ab(0L); fit_b <- fit_ab(1L)
logf("baseline a and b fitted in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

## ---- 1. LOO a vs b (stage 1 and 2 separately) -----------------------------------------------------------------
lc <- loo_compare_ab(fit_a, fit_b); print(lc)
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
th <- temporal_holdout(sy, alts, ev, W_cols, X_cols, Z_cols, cutoff = cutoff, family = "m1", stan_models = models,
                       chains = CHAINS, iter = ITER, warmup = WARMUP, n_draws = 100)
logf("temporal holdout fitted in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
print(th$summary); write.csv(th$summary, file.path(res_dir, "measurement_temporal_holdout.csv"), row.names = FALSE)

## ---- 4. permutation test (rescore mode; refit mode is available but too slow here) -----------------------------
pt <- permutation_test(sy, alts, ev, W_cols, X_cols, Z_cols, fit_a, fit_b, n_perm = 20, mode = "rescore", n_draws = 100)
print(pt$summary); write.csv(pt$summary, file.path(res_dir, "measurement_permutation.csv"), row.names = FALSE)
write.csv(data.frame(perm = seq_along(pt$permuted_advantage), adv = pt$permuted_advantage, a = pt$permuted_a, b = pt$permuted_b),
          file.path(res_dir, "measurement_permutation_draws.csv"), row.names = FALSE)

## ---- 5. systematic measurement error ---------------------------------------------------------------------------
## "recent conflict" in the sim = the state initiated in the previous year (y1 at t-1), positions are at t-1
recent <- sy |> mutate(sy_row = row_number()) |> transmute(state_idx, year = year + 1L, recent = y1)
t0 <- Sys.time()
me <- measurement_error_check(sy, alts, ev, W_cols, X_cols, Z_cols, recent = recent, deltas = c(0.5, 1.0),
                              baseline_cls = cls_b, family = "m1", joint = 1L, stan_models = models,
                              chains = CHAINS, iter = ITER, warmup = WARMUP, r0 = R0, dR_in_X = "X2")
logf("measurement-error refits in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
me$baseline_accuracy <- mean(cls_b$map_type == sim$types)
print(me); write.csv(me, file.path(res_dir, "measurement_error_inflation.csv"), row.names = FALSE)

## ---- 6. benchmarks on the simulated frames (assembly + short brms fits) ---------------------------------------------
ev_keys <- data.frame(state_idx = alts$state_idx[ev$alt_row], target_idx = alts$target_idx[ev$alt_row], year = alts$year[ev$alt_row])
d1 <- assemble_dyad_logit(alts, ev_keys, id_cols = c("state_idx", "target_idx", "year"), z_cols = Z_cols)
stopifnot(sum(d1$y) == nrow(ev))
d1$off <- 0
t0 <- Sys.time()
b1 <- brm(y ~ cinc_ratio + dR + Z1 + Z2 + offset(off) + (1 | state_idx) + (1 | target_idx), data = d1, family = bernoulli(),
          prior = c(set_prior("normal(0, 1)", class = "b"), set_prior("student_t(3, 0, 1.5)", class = "sd")),
          chains = CHAINS, cores = CHAINS, iter = ITER, warmup = WARMUP, seed = V2_SEED, refresh = 0)
d2 <- assemble_gap_regression(sy, alts, ev, X_cols)
b2 <- brm(gap_formula(X_cols), data = d2, prior = set_prior("normal(0, 1)", class = "b"),
          chains = CHAINS, cores = CHAINS, iter = ITER, warmup = WARMUP, seed = V2_SEED, refresh = 0)
logf("benchmarks fitted in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
fe1 <- as.data.frame(fixef(b1)); fe1$model <- "dyad_logit"; fe1$term <- rownames(fe1)
fe2 <- as.data.frame(fixef(b2)); fe2$model <- "gap_mvbind"; fe2$term <- rownames(fe2)
write.csv(bind_rows(fe1, fe2), file.path(res_dir, "benchmarks_sim.csv"), row.names = FALSE)
## truth for the gap regression's dM ~ X1 slope is NOT beta[1,1]: the realised target position is a
## softmax draw around mu, so the slope recovers beta only up to shrinkage; recorded for the report.
logf("run_measurement_tests finished")
