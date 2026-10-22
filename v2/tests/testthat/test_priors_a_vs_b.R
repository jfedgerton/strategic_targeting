## The "a" (independent) and "b" (joint) models must put the SAME marginal prior on the correlations
## they share, corr(c, a1) and corr(c, a2): LKJ(2) on a 3x3 correlation matrix implies
## (r + 1)/2 ~ Beta(2.5, 2.5), Var(r) = 1/6. The "a" model uses that Beta prior directly with
## corr(a1, a2) = 0 and PD enforced. Checked analytically in R and by prior-only Stan runs.
suppressPackageStartupMessages({ library(testthat); library(rstan) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "00_helpers.R")); source(file.path(root, "sim", "simulate.R"))

test_that("Beta(2.5, 2.5) matches the LKJ(2) 3x3 marginal; the PD truncation barely moves the variance", {
  set.seed(123)
  ## LKJ(eta) marginal for a d x d matrix: (r+1)/2 ~ Beta(a, a), a = eta - 1 + d/2 = 2.5 for eta = 2, d = 3
  r_lkj <- 2 * rbeta(2e5, 2.5, 2.5) - 1
  expect_equal(var(r_lkj), 1 / 6, tolerance = 0.02)
  ## "a" prior: independent Beta marginals on r1, r2 restricted to r1^2 + r2^2 < 1
  r1 <- 2 * rbeta(4e5, 2.5, 2.5) - 1; r2 <- 2 * rbeta(4e5, 2.5, 2.5) - 1
  keep <- r1^2 + r2^2 < 1
  expect_gt(mean(keep), 0.9)                                    # truncation removes < 10% of the mass
  expect_equal(var(r1[keep]), 1 / 6, tolerance = 0.06)          # about 1/6 (slightly below, from the truncation)
  expect_equal(var(r2[keep]), var(r1[keep]), tolerance = 0.02)  # symmetric across the two dimensions
})

test_that("prior-only Stan runs: shared correlations have matching marginal variance (about 1/6) in a and b", {
  p <- file.path(root, "stan", "build", "m1_latent_additive.rds")
  if (!file.exists(p)) skip("compiled model missing")
  sm <- readRDS(p)
  sim <- simulate_targeting("M1b", S = 6L, Tt = 3L, seed = 123L, target_rate = 0.6)
  vars <- list()
  for (joint in c(0L, 1L)) {
    dat <- build_stan_data(sim$sy, sim$alts, sim$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"),
                           joint = joint, x_ref = rep(0, 4), prior_only = 1L)
    fit <- sampling(sm, data = dat, chains = 2, iter = 3000, warmup = 1000, seed = 123, refresh = 0,
                    pars = c("r_c_a1", "r_c_a2", "R_Om"))
    d <- rstan::extract(fit, pars = c("r_c_a1", "r_c_a2", "R_Om"))
    vars[[joint + 1]] <- c(v1 = var(d$r_c_a1), v2 = var(d$r_c_a2), m1 = mean(d$r_c_a1), m2 = mean(d$r_c_a2),
                           v_a1a2 = var(d$R_Om[, 2, 3]))
    expect_equal(as.numeric(d$R_Om[, 1, 2]), as.numeric(d$r_c_a1), tolerance = 1e-8)   # exported corr = the prior's r
  }
  for (k in 1:2) {
    expect_equal(unname(vars[[k]]["v1"]), 1 / 6, tolerance = 0.12)
    expect_equal(unname(vars[[k]]["v2"]), 1 / 6, tolerance = 0.12)
    expect_lt(abs(vars[[k]]["m1"]), 0.05); expect_lt(abs(vars[[k]]["m2"]), 0.05)
  }
  expect_equal(unname(vars[[1]]["v1"]), unname(vars[[2]]["v1"]), tolerance = 0.12)
  expect_equal(unname(vars[[1]]["v_a1a2"]), 0)                   # corr(a1, a2) fixed at 0 in "a"
  expect_equal(unname(vars[[2]]["v_a1a2"]), 1 / 6, tolerance = 0.12)
})
