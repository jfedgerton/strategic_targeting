## The vectorised posterior-predictive scorer used by the holdout and permutation tests must agree
## with the brute-force stage-2 likelihood, and permute_dR must be a within-choice-set permutation.
suppressPackageStartupMessages({ library(testthat); library(dplyr) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "00_helpers.R")); source(file.path(root, "R", "typology.R"))
source(file.path(root, "R", "measurement_tests.R")); source(file.path(root, "sim", "simulate.R"))

sim <- simulate_targeting("M2b", S = 7L, Tt = 4L, K = 2L, seed = 123L, target_rate = 0.6)
dat <- build_stan_data(sim$sy, sim$alts, sim$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"), joint = 1L, x_ref = rep(0, 4), K = 2L)
set.seed(5)
S <- dat$S
pars <- list(mu0 = c(0.6, -0.3), beta = matrix(rnorm(8, 0, 0.3), 4, 2), a_i = matrix(rnorm(S * 2), S, 2), b_j = rnorm(S), gamma = c(0.4, -0.2),
             sigma_x = c(0.9, 1.1), rho = 0.3, U = matrix(rnorm(S * 2), S, 2), V = matrix(rnorm(S * 2), S, 2))
## two identical "draws" so the posterior mean equals the single-draw probability (mu0 non-zero on purpose)
draws <- list(mu0 = rbind(pars$mu0, pars$mu0),
              beta = array(rep(pars$beta, each = 2), c(2, 4, 2)), a_i = array(rep(pars$a_i, each = 2), c(2, S, 2)),
              b_j = rbind(pars$b_j, pars$b_j), sigma_x = rbind(pars$sigma_x, pars$sigma_x), rho = c(0.3, 0.3),
              gamma = rbind(pars$gamma, pars$gamma), U = array(rep(pars$U, each = 2), c(2, S, 2)), V = array(rep(pars$V, each = 2), c(2, S, 2)))

test_that("stage2_scores matches brute-force log-likelihood and ranks", {
  sc <- stage2_scores(dat, draws)
  bf <- stage2_loglik_R(dat, pars)
  expect_equal(sc$log_p_chosen, bf$loglik, tolerance = 1e-8)
  expect_equal(nrow(sc), dat$E)
  ## rank of chosen target from brute-force utilities
  rk <- sapply(seq_len(dat$E), function(e) { rows <- dat$cs_start[dat$ev_sy[e]]:dat$cs_end[dat$ev_sy[e]]; rank(-bf$util[rows])[match(dat$ev_row[e], rows)] })
  expect_equal(sc$mean_rank, as.numeric(rk))
  expect_true(all(sc$set_size == S - 1))
  ## both reimplementations fail loudly without mu0
  expect_error(stage2_scores(dat, draws[setdiff(names(draws), "mu0")]), "mu0")
  expect_error(stage2_loglik_R(dat, pars[setdiff(names(pars), "mu0")]), "mu0")
})

test_that("permute_dR shuffles within each choice set and keeps everything else", {
  p <- permute_dR(sim$alts, seed = 1)
  expect_equal(p$dM, sim$alts$dM); expect_equal(p$target_idx, sim$alts$target_idx)
  for (n in unique(sim$alts$sy_row)) expect_equal(sort(p$dR[p$sy_row == n]), sort(sim$alts$dR[sim$alts$sy_row == n]))
  expect_false(isTRUE(all.equal(p$dR, sim$alts$dR)))
})

test_that("permute_dR shuffles only within the stage-2 rows when IV-only rows are present", {
  a <- sim$alts
  set.seed(7)
  a$in_s2 <- 1L
  ## mark two non-chosen rows per state-year as IV-only (in_s2 = 0)
  for (n in unique(a$sy_row)) { r <- which(a$sy_row == n); a$in_s2[r[seq_len(min(2, length(r) - 2))]] <- 0L }
  a$in_iv <- 1L
  p <- permute_dR(a, seed = 1)
  expect_equal(p$dR[a$in_s2 == 0L], a$dR[a$in_s2 == 0L])
  expect_equal(p$dM, a$dM); expect_equal(p$in_s2, a$in_s2); expect_equal(p$sy_row, a$sy_row)
  for (n in unique(a$sy_row)) {
    k <- a$sy_row == n & a$in_s2 == 1L
    expect_equal(sort(p$dR[k]), sort(a$dR[k]))
  }
  expect_false(isTRUE(all.equal(p$dR[a$in_s2 == 1L], a$dR[a$in_s2 == 1L])))
})

test_that("misclassification helper and holdout split keep one state index space", {
  mu <- array(rnorm(200 * S * 2), c(200, S, 2))
  ca <- classify_types(mu, r0 = 0.5); expect_equal(nrow(ca), S)
  ## build_stan_data with an S override larger than the data
  d2 <- build_stan_data(sim$sy, sim$alts, sim$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"), joint = 0L, x_ref = rep(0, 4), S = S + 3)
  expect_equal(d2$S, S + 3)
  expect_error(build_stan_data(sim$sy, sim$alts, sim$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"), joint = 0L, S = 2))
})

test_that("extract_draws handles rstan's 1-d arrays for scalars and subsets every parameter", {
  p <- file.path(root, "stan", "build", "m1_latent_additive.rds")
  if (!file.exists(p)) skip("compiled model missing")
  suppressPackageStartupMessages(library(rstan))
  d7 <- build_stan_data(sim$sy, sim$alts, sim$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"), joint = 1L, x_ref = rep(0, 4))
  f <- sampling(readRDS(p), data = d7, chains = 1, iter = 30, warmup = 10, seed = 123, refresh = 0)
  dr <- extract_draws(f, n_use = 5)
  expect_equal(length(dr$rho), 5); expect_equal(dim(dr$beta), c(5, 4, 2)); expect_equal(dim(dr$sigma_x), c(5, 2))
  expect_equal(dim(dr$mu0), c(5, 2))
  sc <- stage2_scores(d7, dr)
  expect_equal(nrow(sc), d7$E); expect_true(all(is.finite(sc$log_p_chosen)))
})
