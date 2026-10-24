## Stage-2 likelihood in Stan must match a brute-force R computation on a tiny example.
## Uses algorithm = "Fixed_param" so generated quantities are evaluated at supplied inits.
## Requires compiled models in v2/stan/build (Rscript v2/R/compile_models.R).
suppressPackageStartupMessages({ library(testthat); library(rstan) })
## tiny_data(), make_pars(), stan_inits(), skip_if_no_build() come from helper-fixed-draw.R

check_model <- function(model_name, family, joint, with_uv, K = 2L, D_mode = 3L) {
  sm <- skip_if_no_build(model_name)
  dat <- tiny_data(family, joint, K = if (with_uv) K else NULL, D_mode = D_mode)
  p <- make_pars(dat, with_uv)
  st <- stan_inits(dat, p)
  fit <- sampling(sm, data = dat, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0,
                  init = list(st$init), seed = 123, refresh = 0)
  ll_stan <- as.numeric(rstan::extract(fit, "log_lik_stage2")$log_lik_stage2[1, ])
  ll_R <- stage2_loglik_R(dat, st$pars)$loglik
  expect_equal(length(ll_stan), dat$E)
  expect_equal(ll_stan, ll_R, tolerance = 1e-6)
  ## stage-1 log-lik: eta1 = alpha0 + W alpha + c_i + tau * lse
  lse <- stage2_loglik_R(dat, st$pars)$lse
  eta1 <- st$init$alpha0 + dat$W %*% st$init$alpha + st$pars$c_i[dat$sy_state] + st$init$tau * lse
  ll1_R <- dbinom(dat$y1, 1, plogis(as.numeric(eta1)), log = TRUE)
  ll1_stan <- as.numeric(rstan::extract(fit, "log_lik_stage1")$log_lik_stage1[1, ])
  expect_equal(ll1_stan, ll1_R, tolerance = 1e-6)
  ## mu_state at x_ref = 0 equals mu0 + a_i (mu0 is non-zero in make_pars, so a missing intercept would fail here)
  mu <- rstan::extract(fit, "mu_state")$mu_state[1, , ]
  D <- if (D_mode == 3L) 2L else 1L
  mu_expected <- rep(1, dat$S) %o% st$pars$mu0 + matrix(st$pars$a_i, ncol = D)
  expect_true(all(st$pars$mu0 != 0))
  expect_equal(unname(as.numeric(mu)), unname(as.numeric(mu_expected)), tolerance = 1e-6)
  ## mu_event equals mu0 + X[ev_sy] beta + a_i
  mue <- matrix(rstan::extract(fit, "mu_event")$mu_event[1, , ], ncol = D)
  mue_expected <- rep(1, dat$E) %o% st$pars$mu0 + dat$X[dat$ev_sy, , drop = FALSE] %*% st$pars$beta +
    matrix(st$pars$a_i, ncol = D)[dat$sy_state[dat$ev_sy], , drop = FALSE]
  expect_equal(unname(as.numeric(mue)), unname(as.numeric(mue_expected)), tolerance = 1e-6)
  ## Sigma assembled correctly
  Sg <- rstan::extract(fit, "Sigma")$Sigma[1, , ]
  if (D_mode == 3L) {
    expect_equal(Sg[1, 2], st$pars$rho * prod(st$pars$sigma_x), tolerance = 1e-8)
    if (joint == 0) expect_equal(Sg[1, 2], 0)
  } else {
    expect_equal(as.numeric(Sg), st$pars$sigma_x[1]^2, tolerance = 1e-8)
  }
  invisible(TRUE)
}

test_that("M1 joint: Stan stage-2 log-lik equals brute-force R", check_model("m1_latent_additive", "M1b", 1L, FALSE))
test_that("M1 independent: Stan stage-2 log-lik equals brute-force R", check_model("m1_latent_additive", "M1b", 0L, FALSE))
test_that("M2 joint (K=2): Stan stage-2 log-lik equals brute-force R", check_model("m2_ame", "M2b", 1L, TRUE))
test_that("M2 independent (K=2): Stan stage-2 log-lik equals brute-force R", check_model("m2_ame", "M2b", 0L, TRUE))
test_that("M1 material-only (D_mode 1): Stan stage-2 log-lik equals brute-force R", check_model("m1_latent_additive", "M1b", 0L, FALSE, D_mode = 1L))
test_that("M1 relational-only (D_mode 2): Stan stage-2 log-lik equals brute-force R", check_model("m1_latent_additive", "M1b", 0L, FALSE, D_mode = 2L))
test_that("M2 material-only (D_mode 1, K=2): Stan stage-2 log-lik equals brute-force R", check_model("m2_ame", "M2b", 0L, TRUE, D_mode = 1L))
test_that("M2 relational-only (D_mode 2, K=2): Stan stage-2 log-lik equals brute-force R", check_model("m2_ame", "M2b", 0L, TRUE, D_mode = 2L))
test_that("M2 joint (K=1): Stan stage-2 log-lik equals brute-force R", check_model("m2_ame", "M2b", 1L, TRUE, K = 1L))
test_that("M2 joint (K=3): Stan stage-2 log-lik equals brute-force R", check_model("m2_ame", "M2b", 1L, TRUE, K = 3L))

test_that("softmax probabilities sum to one within each choice set (R reference)", {
  dat <- tiny_data("M1b", 1L)
  p <- make_pars(dat, FALSE)
  out <- stage2_loglik_R(dat, p)
  for (n in seq_len(dat$N1)) {
    rows <- dat$cs_start[n]:dat$cs_end[n]
    expect_equal(sum(exp(out$util[rows] - out$lse[n])), 1, tolerance = 1e-10)
  }
})

test_that("sampled alternatives: IV rows are outcome-blind, each event has its own set, Stan log-lik matches R", {
  sim <- simulate_targeting("M1b", S = 8L, Tt = 3L, seed = 123L, target_rate = 0.6)
  m <- 3
  sa <- sample_alternatives(sim$alts, sim$ev, m = m)
  n_iv <- sa$n_iv
  iv <- sa$alts[seq_len(n_iv), ]
  expect_true(all(table(iv$sy_row) == m))                                         # exactly m IV rows per state-year
  ## chosen targets are NOT forced into the IV rows: inclusion rate close to m / n_full, far below 1
  chosen_key <- paste(sim$alts$sy_row[sim$ev$alt_row], sim$alts$target_idx[sim$ev$alt_row])
  inc <- mean(chosen_key %in% paste(iv$sy_row, iv$target_idx))
  expect_lt(inc, 0.9); expect_gt(inc, 0.05)
  ## event blocks: chosen target first, then m others from the same state-year, independent per event
  expect_true(all(sa$alts$sy_row[sa$ev$alt_row] == sa$ev$sy_row))
  expect_true(all(sa$alts$target_idx[sa$ev$alt_row] == sim$alts$target_idx[sim$ev$alt_row]))
  expect_true(all(sa$ev$ev_end - sa$ev$ev_start + 1 == m + 1))
  expect_true(all(sa$ev$alt_row == sa$ev$ev_start))
  expect_true(all(sa$ev$ev_start > n_iv))
  dat <- build_stan_data(sim$sy, sa$alts, sa$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"),
                         joint = 1L, sampled = 1L, n_full = sa$n_full, x_ref = rep(0, 4))
  expect_equal(dat$sampled, 1L); expect_true(all(dat$n_full == 7)); expect_equal(dat$N2, n_iv + (m + 1) * dat$E)
  ## Stan (Fixed_param) stage-2 log-lik on the sampled data equals the brute-force R computation
  sm <- skip_if_no_build("m1_latent_additive")
  p <- make_pars(dat, FALSE); st <- stan_inits(dat, p)
  fit <- sampling(sm, data = dat, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(st$init), seed = 123, refresh = 0)
  ll_stan <- as.numeric(rstan::extract(fit, "log_lik_stage2")$log_lik_stage2[1, ])
  bf <- stage2_loglik_R(dat, st$pars)
  expect_equal(ll_stan, bf$loglik, tolerance = 1e-6)
  ## and the IV correction log(n_full / n_sub) is applied in stage 1
  eta1 <- st$init$alpha0 + dat$W %*% st$init$alpha + st$pars$c_i[dat$sy_state] + st$init$tau * (bf$lse + log(dat$n_full / m))
  ll1_stan <- as.numeric(rstan::extract(fit, "log_lik_stage1")$log_lik_stage1[1, ])
  expect_equal(ll1_stan, dbinom(dat$y1, 1, plogis(as.numeric(eta1)), log = TRUE), tolerance = 1e-6)
})

## ---- a and b likelihoods unchanged versus the pre-D_mode programs (commit 82a3b24) -----------------------
## The reference models are compiled from `git show 82a3b24:v2/stan/<file>` into V2_REF_BUILD
## (set the env var to a directory holding m1_latent_additive_82a3b24.rds and m2_ame_82a3b24.rds).
ref_dir <- Sys.getenv("V2_REF_BUILD", unset = "")
check_reference <- function(model_name, family, joint, with_uv) {
  if (ref_dir == "" || !file.exists(file.path(ref_dir, paste0(model_name, "_82a3b24.rds")))) skip("reference build not available (set V2_REF_BUILD)")
  sm_new <- skip_if_no_build(model_name); sm_ref <- readRDS(file.path(ref_dir, paste0(model_name, "_82a3b24.rds")))
  dat <- tiny_data(family, joint, K = if (with_uv) 2L else NULL)
  p <- make_pars(dat, with_uv); st <- stan_inits(dat, p)
  f_new <- sampling(sm_new, data = dat, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(st$init), seed = 123, refresh = 0)
  f_ref <- sampling(sm_ref, data = dat, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(st$init), seed = 123, refresh = 0)
  for (par in c("log_lik_stage2", "log_lik_stage1", "mu_state", "Sigma", "Omega"))
    expect_equal(as.numeric(rstan::extract(f_new, par)[[1]]), as.numeric(rstan::extract(f_ref, par)[[1]]), tolerance = 1e-12, info = par)
}
test_that("reference: M1a unchanged", check_reference("m1_latent_additive", "M1b", 0L, FALSE))
test_that("reference: M1b unchanged", check_reference("m1_latent_additive", "M1b", 1L, FALSE))
test_that("reference: M2a unchanged", check_reference("m2_ame", "M2b", 0L, TRUE))
test_that("reference: M2b unchanged", check_reference("m2_ame", "M2b", 1L, TRUE))
