## The "a" (independent) and "b" (joint) models must put the SAME marginal prior on the correlations they share, corr(c, a1) and
## corr(c, a2): LKJ(2) on a 3x3 correlation matrix implies (r + 1)/2 ~ Beta(2.5, 2.5), Var(r) = 1/6. Since the review fix of
## 2026-10-25 (B4) the "a" model uses the LKJ(2) density restricted to corr(a1, a2) = 0, p(r1, r2) propto 1 - r1^2 - r2^2 on the
## unit disk, whose marginals are exactly Beta(2.5, 2.5); the previous code (two Beta marginals truncated to the disk) gave
## Var about 0.158. The "b" model's Sigma correlation rho_x gets the same marginal explicitly. Checked analytically in R and by
## prior-only Stan runs.
suppressPackageStartupMessages({ library(testthat); library(rstan) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "spec.R")); source(file.path(root, "R", "helpers_data.R")); source(file.path(root, "R", "helpers_model.R")); source(file.path(root, "sim", "simulate.R"))

test_that("the disk density 1 - r1^2 - r2^2 has exact Beta(2.5, 2.5) marginals; the old truncated-Beta scheme does not", {
  set.seed(123)
  ## LKJ(eta) marginal for a d x d matrix: (r+1)/2 ~ Beta(a, a), a = eta - 1 + d/2 = 2.5 for eta = 2, d = 3
  r_lkj <- 2 * rbeta(4e5, 2.5, 2.5) - 1
  expect_equal(var(r_lkj), 1 / 6, tolerance = 0.01)
  ## the marginal of r1 under p propto (1 - r1^2 - r2^2) on the disk is (4/3)(1 - r1^2)^(3/2), i.e. Beta(2.5, 2.5) on (r1 + 1)/2
  f <- function(r1) sapply(r1, function(a) integrate(function(b) 1 - a^2 - b^2, -sqrt(1 - a^2), sqrt(1 - a^2))$value)
  g <- seq(-0.95, 0.95, by = 0.05)
  expect_equal(f(g) / (pi / 2), dbeta((g + 1) / 2, 2.5, 2.5) / 2, tolerance = 1e-8)   # disk integral = pi/2; Beta density in r has Jacobian 1/2
  ## rejection sampling from the disk density: Var(r1) = Var(r2) = 1/6
  n <- 2e6; u <- runif(n, -1, 1); v <- runif(n, -1, 1); keep <- runif(n) < (1 - u^2 - v^2) & u^2 + v^2 < 1
  expect_equal(var(u[keep]), 1 / 6, tolerance = 0.01); expect_equal(var(v[keep]), 1 / 6, tolerance = 0.01)
  ## the parameterisation of the Stan program: r1 and q = r2 / sqrt(1 - r1^2) are independent, (r1 + 1)/2 ~ Beta(2.5, 2.5), (q + 1)/2 ~ Beta(2, 2)
  q <- v[keep] / sqrt(1 - u[keep]^2)
  expect_lt(abs(cor(u[keep], q)), 0.01); expect_equal(var(q), 1 / 5, tolerance = 0.01)      # Beta(2, 2) on (q+1)/2: Var(q) = 4 * 1/20 = 1/5
  ## the OLD prior (independent Beta(2.5, 2.5) marginals restricted to r1^2 + r2^2 < 1) is not that density: Var about 0.158
  r1 <- 2 * rbeta(4e5, 2.5, 2.5) - 1; r2 <- 2 * rbeta(4e5, 2.5, 2.5) - 1
  keep_old <- r1^2 + r2^2 < 1
  expect_lt(var(r1[keep_old]), 1 / 6 - 0.005); expect_equal(var(r1[keep_old]), 0.158, tolerance = 0.02)
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
    fit <- sampling(sm, data = dat, chains = 4, iter = 4000, warmup = 1000, seed = 123, refresh = 0,
                    pars = c("r_c_a1", "r_c_a2", "R_Om", "rho"))
    d <- rstan::extract(fit, pars = c("r_c_a1", "r_c_a2", "R_Om", "rho"))
    vars[[joint + 1]] <- c(v1 = var(d$r_c_a1), v2 = var(d$r_c_a2), m1 = mean(d$r_c_a1), m2 = mean(d$r_c_a2),
                           v_a1a2 = var(d$R_Om[, 2, 3]), v_rho = var(d$rho),
                           ks1 = unname(suppressWarnings(ks.test(d$r_c_a1, function(x) pbeta((x + 1) / 2, 2.5, 2.5)))$statistic),
                           ks2 = unname(suppressWarnings(ks.test(d$r_c_a2, function(x) pbeta((x + 1) / 2, 2.5, 2.5)))$statistic))
    expect_equal(as.numeric(d$R_Om[, 1, 2]), as.numeric(d$r_c_a1), tolerance = 1e-8)   # exported corr = the prior's r
  }
  for (k in 1:2) {
    expect_equal(unname(vars[[k]]["v1"]), 1 / 6, tolerance = 0.06)     # 0.158 (the old "a" prior) would fail here
    expect_equal(unname(vars[[k]]["v2"]), 1 / 6, tolerance = 0.06)
    expect_lt(abs(vars[[k]]["m1"]), 0.04); expect_lt(abs(vars[[k]]["m2"]), 0.04)
    expect_lt(unname(vars[[k]]["ks1"]), 0.05); expect_lt(unname(vars[[k]]["ks2"]), 0.05)   # whole marginal, not just the variance
  }
  expect_equal(unname(vars[[1]]["v1"]), unname(vars[[2]]["v1"]), tolerance = 0.08)
  expect_equal(unname(vars[[1]]["v_a1a2"]), 0)                   # corr(a1, a2) fixed at 0 in "a"
  expect_equal(unname(vars[[2]]["v_a1a2"]), 1 / 6, tolerance = 0.08)
  expect_equal(unname(vars[[1]]["v_rho"]), 0)                    # Sigma diagonal in "a"
  expect_equal(unname(vars[[2]]["v_rho"]), 1 / 6, tolerance = 0.08)   # rho_x: explicit Beta(2.5, 2.5) marginal (was uniform, Var 1/3)
})
