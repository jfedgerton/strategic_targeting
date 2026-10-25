## ============================================================================
## run_benchmarks_sim.R -- exercise benchmarks.R (directed-dyad logit, old gap regression) on the
## simulated frames of the 30-state world. Separate from run_measurement_tests.R because the brms
## compile needs several GB. Usage: Rscript v2/sim/run_benchmarks_sim.R [chains] [iter] [warmup]
## ============================================================================
suppressPackageStartupMessages({ library(brms); library(dplyr) })
args <- commandArgs(trailingOnly = TRUE)
CHAINS <- if (length(args) >= 1) as.integer(args[1]) else 2L
ITER   <- if (length(args) >= 2) as.integer(args[2]) else 500L
WARMUP <- if (length(args) >= 3) as.integer(args[3]) else 250L
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
sim_dir <- normalizePath(dirname(this_file)); Sys.setenv(V2_ROOT = normalizePath(file.path(sim_dir, "..")))
source(file.path(sim_dir, "..", "R", "00_helpers.R")); source(file.path(sim_dir, "..", "R", "benchmarks.R")); source(file.path(sim_dir, "simulate.R"))
log_open("run_benchmarks_sim"); set.seed(V2_SEED)
res_dir <- file.path(sim_dir, "results")
sim <- simulate_current("M1b", S = 30L, Tt = 12L, seed = V2_SEED)          # the current-model world (2026-10-25, B5)
sy <- sim$sy; alts <- sim$alts |> filter(in_s2 == 1L); ev <- sim$ev; X_cols <- sim$cols$X; Z_cols <- sim$cols$Z
ev <- ev |> mutate(alt_row = match(paste(sim$alts$sy_row[alt_row], sim$alts$target_idx[alt_row]), paste(alts$sy_row, alts$target_idx)))
stopifnot(!anyNA(ev$alt_row))
## ---- 6. benchmarks on the simulated frames (assembly + short brms fits) ---------------------------------------------
ev_keys <- data.frame(state_idx = alts$state_idx[ev$alt_row], target_idx = alts$target_idx[ev$alt_row], year = alts$year[ev$alt_row])
## cinc_ratio = plogis(raw log material gap) (B10); simulate_current() carries dM_logratio like the real alternatives files
d1 <- assemble_dyad_logit(alts, ev_keys, id_cols = c("state_idx", "target_idx", "year"), z_cols = Z_cols, logratio_col = "dM_logratio")
stopifnot(sum(d1$y) == nrow(ev))
d1$off <- 0
t0 <- Sys.time()
b1 <- brm(bf(as.formula(paste("y ~ cinc_ratio + dR +", paste(Z_cols, collapse = " + "), "+ offset(off) + (1 | state_idx) + (1 | target_idx)"))), data = d1, family = bernoulli(),
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
logf("run_benchmarks_sim finished")
