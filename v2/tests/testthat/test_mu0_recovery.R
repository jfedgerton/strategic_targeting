## The mu intercept mu0 (population-average ideal gap, added 2026-10-23) must be recovered from data simulated
## with a known non-zero mu0. Simulator truth: mu_n = mu0 + X_n beta + a_i with the planted type offsets a_i
## (type centres average to zero by construction, so the identified population mean at X = 0 is mu0 up to the
## realised mean of a_i, about 0.35 / sqrt(S)). One M1b fit (S = 30, 10 years, 2 chains x 800), about 3 min on 2 cores.
## Requires the compiled M1 model (Rscript v2/R/compile_models.R).
suppressPackageStartupMessages({ library(testthat); library(rstan) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
## simulate_targeting(), build_stan_data(), skip_if_no_build() come from helper-fixed-draw.R

test_that("M1b recovers a known non-zero mu0 within its 90% interval", {
  sm <- skip_if_no_build("m1_latent_additive")
  mu0_true <- c(0.8, -0.5)
  S <- 30L
  sim <- simulate_targeting("M1b", S = S, Tt = 10L, seed = 123L, target_rate = 0.5, mu0 = mu0_true)
  expect_equal(sim$truth$mu0, mu0_true)
  expect_equal(sim$truth$mu_state_ref0, rep(1, S) %o% mu0_true + sim$truth$a_i)
  dat <- build_stan_data(sim$sy, sim$alts, sim$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"),
                         joint = 1L, x_ref = rep(0, 4))
  expect_equal(dat$prior_mu0_sd, 2)
  rstan_options(threads_per_chain = 1)
  fit <- sampling(sm, data = dat, chains = 2, iter = 800, warmup = 400, seed = 123, refresh = 0,
                  cores = min(2L, parallel::detectCores()), control = list(adapt_delta = 0.9),
                  pars = c("mu0", "beta", "sigma_x", "sd_re", "tau", "mu_state"))
  post <- rstan::extract(fit, pars = "mu0")$mu0           # [draws, 2]
  expect_equal(ncol(post), 2L)
  lo <- apply(post, 2, quantile, 0.05); hi <- apply(post, 2, quantile, 0.95); pm <- colMeans(post)
  a_bar <- colMeans(sim$truth$a_i)                         # realised offset of the planted a_i from zero
  cat(sprintf("\nmu0 recovery: truth = (%.2f, %.2f); posterior mean = (%.3f, %.3f); 90%% interval = [%.3f, %.3f] x [%.3f, %.3f]; mean(a_i) = (%.3f, %.3f)\n",
              mu0_true[1], mu0_true[2], pm[1], pm[2], lo[1], hi[1], lo[2], hi[2], a_bar[1], a_bar[2]))
  rh <- summary(fit, pars = "mu0")$summary[, "Rhat"]
  expect_true(all(rh < 1.1))
  expect_true(all(mu0_true >= lo & mu0_true <= hi))        # 90% interval covers the truth on both axes
  expect_true(all(sign(pm) == sign(mu0_true)))             # and the intercept is not shrunk to zero
  expect_gt(lo[1], 0)                                      # the larger intercept excludes zero
  ## mu_state at x_ref = 0 is mu0 + a_i: its posterior mean tracks the simulator's mu_state_ref0
  ms <- apply(rstan::extract(fit, pars = "mu_state")$mu_state, c(2, 3), mean)
  expect_gt(cor(as.numeric(ms), as.numeric(sim$truth$mu_state_ref0)), 0.8)
})
