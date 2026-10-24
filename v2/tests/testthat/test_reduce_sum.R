## reduce_sum (within-chain parallelism) must not change the model: for a fixed parameter draw the
## log density, its gradient and the generated log_lik must agree between the serial path
## (grainsize >= N1, one slice) and sliced evaluation with threads_per_chain = 2 and 4, for M1 and M2
## and every D_mode (m, r, a, b). Requires compiled models.
suppressPackageStartupMessages({ library(testthat); library(rstan) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
## make_pars() and stan_inits() come from helper-fixed-draw.R

## complete_init() comes from helper-fixed-draw.R
eval_fixed <- function(sm, dat, init, threads) {
  rstan_options(threads_per_chain = threads)
  fit <- sampling(sm, data = dat, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(init), seed = 123, refresh = 0)
  up <- unconstrain_pars(fit, complete_init(fit, init))
  lp <- log_prob(fit, up, adjust_transform = TRUE, gradient = TRUE)
  list(lp = as.numeric(lp), grad = as.numeric(attr(lp, "gradient")),
       ll2 = as.numeric(rstan::extract(fit, "log_lik_stage2")$log_lik_stage2[1, ]),
       ll1 = as.numeric(rstan::extract(fit, "log_lik_stage1")$log_lik_stage1[1, ]))
}

check_threads <- function(model_name, family, joint, with_uv, D_mode, K = 2L) {
  p <- file.path(root, "stan", "build", paste0(model_name, ".rds"))
  if (!file.exists(p)) skip(paste("compiled model missing:", p))
  sm <- readRDS(p)
  sim <- simulate_targeting(family, S = 7L, Tt = 5L, K = K, seed = 123L, target_rate = 0.6)
  base <- build_stan_data(sim$sy, sim$alts, sim$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"),
                          joint = joint, x_ref = rep(0, 4), K = if (with_uv) K else NULL, D_mode = D_mode)
  pars <- make_pars(base, with_uv); init <- stan_inits(base, pars)$init
  serial <- eval_fixed(sm, base, init, threads = 1)                    # grainsize = N1: one slice
  expect_true(is.finite(serial$lp)); expect_equal(length(serial$ll2), base$E)
  for (threads in c(2L, 4L)) {
    dat <- base; dat$grainsize <- choose_grainsize(base$N1, threads)  # several slices
    expect_lt(dat$grainsize, base$N1)
    par <- eval_fixed(sm, dat, init, threads = threads)
    expect_equal(par$lp, serial$lp, tolerance = 1e-10, info = sprintf("%s D_mode %d threads %d lp", model_name, D_mode, threads))
    expect_equal(par$grad, serial$grad, tolerance = 1e-8, info = sprintf("%s D_mode %d threads %d gradient", model_name, D_mode, threads))
    expect_equal(par$ll2, serial$ll2, tolerance = 1e-12)             # generated quantities are serial: identical
    expect_equal(par$ll1, serial$ll1, tolerance = 1e-12)
  }
  ## grainsize 1 (one state-year per slice) is the finest partition
  dat <- base; dat$grainsize <- 1L
  fine <- eval_fixed(sm, dat, init, threads = 2L)
  expect_equal(fine$lp, serial$lp, tolerance = 1e-10); expect_equal(fine$grad, serial$grad, tolerance = 1e-8)
  rstan_options(threads_per_chain = 1)
  invisible(TRUE)
}

for (fam in c("M1", "M2")) for (code in c("m", "r", "a", "b")) {
  mn <- if (fam == "M1") "m1_latent_additive" else "m2_ame"
  D_mode <- switch(code, m = 1L, r = 2L, 3L); joint <- as.integer(code == "b")
  test_that(sprintf("reduce_sum: %s%s identical lp / gradient / log_lik with threads 2 and 4", fam, code),
            check_threads(mn, if (fam == "M1") "M1b" else "M2b", joint, fam == "M2", D_mode))
}
test_that("reduce_sum: M2 with K = 1 and K = 3", {
  check_threads("m2_ame", "M2b", 1L, TRUE, 3L, K = 1L)
  check_threads("m2_ame", "M2b", 1L, TRUE, 3L, K = 3L)
})
test_that("reduce_sum: sampled-alternatives mode (per-event sets) also agrees across slicings", {
  p <- file.path(root, "stan", "build", "m1_latent_additive.rds"); if (!file.exists(p)) skip("compiled model missing")
  sm <- readRDS(p)
  sim <- simulate_targeting("M1b", S = 8L, Tt = 3L, seed = 123L, target_rate = 0.6)
  sa <- sample_alternatives(sim$alts, sim$ev, m = 3)
  base <- build_stan_data(sim$sy, sa$alts, sa$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"), joint = 1L, sampled = 1L,
                          n_full = sa$n_full, x_ref = rep(0, 4))
  pars <- make_pars(base, FALSE); init <- stan_inits(base, pars)$init
  serial <- eval_fixed(sm, base, init, 1L); dat <- base; dat$grainsize <- 2L; par <- eval_fixed(sm, dat, init, 2L)
  expect_equal(par$lp, serial$lp, tolerance = 1e-10); expect_equal(par$grad, serial$grad, tolerance = 1e-8)
  expect_equal(par$ll2, stage2_loglik_R(base, stan_inits(base, pars)$pars)$loglik, tolerance = 1e-6)
  rstan_options(threads_per_chain = 1)
})
test_that("reduce_sum: masked choice sets (cs_mixed path) agree between the serial and the sliced evaluation, M1 and M2", {
  for (mn in c("m1_latent_additive", "m2_ame")) {
    p <- file.path(root, "stan", "build", paste0(mn, ".rds")); if (!file.exists(p)) skip("compiled model missing")
    sm <- readRDS(p); with_uv <- mn == "m2_ame"
    sim <- simulate_targeting(if (with_uv) "M2b" else "M1b", S = 7L, Tt = 5L, K = 2L, seed = 123L, target_rate = 0.6)
    M <- mask_sets(sim)
    base <- build_stan_data(M$sy, M$alts, M$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"), joint = 1L, x_ref = rep(0, 4), K = if (with_uv) 2L else NULL)
    expect_true(any(base$in_iv == 0L) && any(base$in_s2 == 0L))
    pars <- make_pars(base, with_uv); init <- stan_inits(base, pars)$init
    serial <- eval_fixed(sm, base, init, 1L)
    expect_true(is.finite(serial$lp))
    for (gs in c(1L, 3L)) { dat <- base; dat$grainsize <- gs; par <- eval_fixed(sm, dat, init, 2L)
      expect_equal(par$lp, serial$lp, tolerance = 1e-10); expect_equal(par$grad, serial$grad, tolerance = 1e-8)
      expect_equal(par$ll2, serial$ll2, tolerance = 1e-12); expect_equal(par$ll1, serial$ll1, tolerance = 1e-12) }
    expect_equal(serial$ll2, stage2_loglik_R(base, stan_inits(base, pars)$pars)$loglik, tolerance = 1e-6)
    rstan_options(threads_per_chain = 1)
  }
})
