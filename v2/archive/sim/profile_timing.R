## ARCHIVED 2026-10-25: one-off gradient profile and thread timing (2026-10-22); not part of the pipeline.
## ============================================================================
## profile_timing.R -- where the gradient time goes, and the reduce_sum speed-up
##
## Usage: Rscript v2/sim/profile_timing.R [n_grad_evals] [iter_timing]
## 1. Gradient-time breakdown on the 40-state / 25-year simulated cell (N1 = 1000, N2 = 39,000), M1b and M2b:
##    log_prob + gradient timed for (a) the full model, (b) the same data with E = 0 (stage 1 + IV + priors
##    only), (c) prior_only = 1 (priors + transformed parameters only). Differences give the stage-2 event
##    terms, stage 1 + IV (dominated by the utility computation over all alternative rows), and priors.
##    rstan cannot read Stan profile() sections (that needs CmdStan), so this is the readable breakdown.
## 2. Sampling time for 1 chain x iter_timing iterations with threads_per_chain = 1, 2, 4 (seed 123; the
##    container has 2 cores, so 4 threads is oversubscribed and reported as such).
## Writes v2/sim/results/profile_breakdown.csv and threads_timing.csv.
## ============================================================================
suppressPackageStartupMessages({ library(rstan); library(dplyr) })
args <- commandArgs(trailingOnly = TRUE)
N_EVAL <- if (length(args) >= 1) as.integer(args[1]) else 20L
ITER_T <- if (length(args) >= 2) as.integer(args[2]) else 40L
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
sim_dir <- normalizePath(dirname(this_file)); Sys.setenv(V2_ROOT = normalizePath(file.path(sim_dir, "..", "..")))
source(file.path(sim_dir, "..", "..", "R", "spec.R")); source(file.path(sim_dir, "..", "..", "R", "helpers_data.R")); source(file.path(sim_dir, "..", "..", "R", "helpers_model.R")); source(file.path(sim_dir, "..", "..", "sim", "simulate.R"))
source(file.path(sim_dir, "..", "..", "tests", "testthat", "helper-fixed-draw.R"))
log_open("profile_timing"); set.seed(V2_SEED); rstan_options(auto_write = FALSE)
res_dir <- file.path(sim_dir, "results")
models <- list(m1_latent_additive = readRDS(file.path(v2_root(), "stan", "build", "m1_latent_additive.rds")),
               m2_ame = readRDS(file.path(v2_root(), "stan", "build", "m2_ame.rds")))
cores <- parallel::detectCores()
logf("container cores: %d", cores)

grad_ms <- function(sm, dat, init, n = N_EVAL) {
  rstan_options(threads_per_chain = 1)
  fit <- sampling(sm, data = dat, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(init), seed = 123, refresh = 0)
  up <- unconstrain_pars(fit, complete_init(fit, init))
  invisible(log_prob(fit, up, gradient = TRUE))                      # warm-up call
  t0 <- Sys.time(); for (i in seq_len(n)) invisible(log_prob(fit, up, gradient = TRUE))
  1000 * as.numeric(difftime(Sys.time(), t0, units = "secs")) / n
}
rows <- list(); trows <- list()
for (fam in c("M1b", "M2b")) {
  mn <- if (fam == "M1b") "m1_latent_additive" else "m2_ame"; with_uv <- fam == "M2b"
  sim <- simulate_targeting(fam, seed = V2_SEED)
  dat <- build_stan_data(sim$sy, sim$alts, sim$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"), joint = 1L, x_ref = rep(0, 4),
                         K = if (with_uv) 2L else NULL)
  init <- stan_inits(dat, make_pars(dat, with_uv))$init
  ## (a) full, (b) no events, (c) priors only
  dat_noev <- dat; dat_noev$E <- 0L
  for (nm in c("ev_sy", "ev_row", "ev_start", "ev_end")) dat_noev[[nm]] <- integer(0)
  dat_prior <- dat; dat_prior$prior_only <- 1L
  t_full <- grad_ms(models[[mn]], dat, init); t_noev <- grad_ms(models[[mn]], dat_noev, init); t_prior <- grad_ms(models[[mn]], dat_prior, init)
  rows[[fam]] <- data.frame(model = fam, N1 = dat$N1, N2 = dat$N2, E = dat$E, ms_full = t_full,
                            ms_stage2_events = t_full - t_noev, ms_stage1_and_iv_utilities = t_noev - t_prior, ms_priors_and_transforms = t_prior,
                            share_stage2_events = (t_full - t_noev) / t_full, share_stage1_iv_utilities = (t_noev - t_prior) / t_full,
                            share_priors = t_prior / t_full)
  logf("%s gradient: full %.1f ms | events %.1f | stage1+IV utilities %.1f | priors %.1f", fam, t_full, t_full - t_noev, t_noev - t_prior, t_prior)
  ## threads timing: 1 chain x ITER_T iterations, identical seed and data except grainsize
  for (thr in c(1L, 2L, 4L)) {
    rstan_options(threads_per_chain = thr)
    d <- dat; d$grainsize <- choose_grainsize(dat$N1, thr)
    t0 <- Sys.time()
    f <- sampling(models[[mn]], data = d, chains = 1, iter = ITER_T, warmup = ITER_T %/% 2, seed = V2_SEED, refresh = 0,
                  control = list(adapt_delta = 0.9, max_treedepth = 10))
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    sp <- get_sampler_params(f, inc_warmup = TRUE)[[1]]
    trows[[paste(fam, thr)]] <- data.frame(model = fam, threads_per_chain = thr, grainsize = d$grainsize, cores_available = cores,
                                           oversubscribed = thr > cores, iterations = ITER_T, elapsed_s = el,
                                           s_per_iteration = el / ITER_T, leapfrogs_total = sum(sp[, "n_leapfrog__"]),
                                           ms_per_gradient = 1000 * el / sum(sp[, "n_leapfrog__"]))
    logf("%s threads %d: %.1f s for %d iterations (%.2f ms per gradient over %d leapfrogs)", fam, thr, el, ITER_T,
         1000 * el / sum(sp[, "n_leapfrog__"]), sum(sp[, "n_leapfrog__"]))
  }
  rstan_options(threads_per_chain = 1)
}
write.csv(bind_rows(rows), file.path(res_dir, "profile_breakdown.csv"), row.names = FALSE)
write.csv(bind_rows(trows), file.path(res_dir, "threads_timing.csv"), row.names = FALSE)
logf("profile_timing done")
