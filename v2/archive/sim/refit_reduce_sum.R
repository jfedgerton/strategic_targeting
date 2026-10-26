## ARCHIVED 2026-10-25: one-off reduce_sum refit (2026-10-22); not part of the pipeline.
## refit_reduce_sum.R -- refit M1b on the M1b simulated data with the reduce_sum program and
## threads_per_chain = 2 (chains sequential), then compare with the earlier serial fit.
## Usage: Rscript v2/sim/refit_reduce_sum.R [chains] [iter] [warmup] [threads]
suppressPackageStartupMessages({ library(rstan); library(dplyr) })
args <- commandArgs(trailingOnly = TRUE)
CHAINS <- if (length(args) >= 1) as.integer(args[1]) else 2L; ITER <- if (length(args) >= 2) as.integer(args[2]) else 600L
WARMUP <- if (length(args) >= 3) as.integer(args[3]) else 300L; THREADS <- if (length(args) >= 4) as.integer(args[4]) else 2L
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
sim_dir <- normalizePath(dirname(this_file)); Sys.setenv(V2_ROOT = normalizePath(file.path(sim_dir, "..", "..")))
source(file.path(sim_dir, "..", "..", "R", "spec.R")); source(file.path(sim_dir, "..", "..", "R", "helpers_data.R")); source(file.path(sim_dir, "..", "..", "R", "helpers_model.R")); source(file.path(sim_dir, "..", "..", "sim", "simulate.R"))
log_open("refit_reduce_sum"); rstan_options(auto_write = FALSE); rstan_options(threads_per_chain = THREADS); options(mc.cores = 1)
sim <- simulate_targeting("M1b", seed = V2_SEED)
dat <- build_stan_data(sim$sy, sim$alts, sim$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"), joint = 1L, x_ref = rep(0, 4),
                       grainsize = choose_grainsize(nrow(sim$sy), THREADS))
logf("threads = %d, grainsize = %d, chains sequential", THREADS, dat$grainsize)
sm <- readRDS(file.path(v2_root(), "stan", "build", "m1_latent_additive.rds"))
t0 <- Sys.time()
fit <- sampling(sm, data = dat, chains = CHAINS, iter = ITER, warmup = WARMUP, seed = V2_SEED, refresh = 0,
                control = list(adapt_delta = 0.9, max_treedepth = 10), pars = c("z_re", "z_b", "V_raw", "U_free", "L_Om_joint", "L_Om", "re"), include = FALSE)
el <- as.numeric(difftime(Sys.time(), t0, units = "mins")); logf("done in %.1f min", el)
saveRDS(list(fit = fit, dat = dat, truth_family = "M1b", elapsed_min = el, threads = THREADS), file.path(sim_dir, "fits", "M1b_M1b_reduce_sum.rds"))
