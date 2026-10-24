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

## ---- two choice sets per state-year (Codex review P1-1, 2026-10-24) ----------------------------------------------------------
## IV[n] runs over the in_iv rows, the stage-2 softmax over the in_s2 rows of the same block. Stan must match the brute-force R
## computation with masked sets, and must reproduce the flag-free program exactly when every row is in both sets.
check_masked <- function(model_name, family, joint, with_uv, K = 2L, D_mode = 3L) {
  sm <- skip_if_no_build(model_name)
  sim <- simulate_targeting(family, S = 7L, Tt = 5L, K = K, seed = 123L, target_rate = 0.6)
  M <- mask_sets(sim)
  Kk <- if (with_uv) K else NULL
  dat <- build_stan_data(M$sy, M$alts, M$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"), joint = joint, x_ref = rep(0, 4), K = Kk, D_mode = D_mode)
  expect_true(any(dat$in_iv == 0L)); expect_true(any(dat$in_s2 == 0L)); expect_true(all(dat$in_s2[dat$ev_row] == 1L))
  p <- make_pars(dat, with_uv); st <- stan_inits(dat, p)
  fit <- sampling(sm, data = dat, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(st$init), seed = 123, refresh = 0)
  bf <- stage2_loglik_R(dat, st$pars)
  expect_equal(as.numeric(rstan::extract(fit, "log_lik_stage2")$log_lik_stage2[1, ]), bf$loglik, tolerance = 1e-6)
  eta1 <- st$init$alpha0 + dat$W %*% st$init$alpha + st$pars$c_i[dat$sy_state] + st$init$tau * bf$lse
  expect_equal(as.numeric(rstan::extract(fit, "log_lik_stage1")$log_lik_stage1[1, ]), dbinom(dat$y1, 1, plogis(as.numeric(eta1)), log = TRUE), tolerance = 1e-6)
  ## the masked likelihood differs from the whole-block one where the sets differ, and only there
  d_all <- dat; d_all$in_iv[] <- 1L; d_all$in_s2[] <- 1L
  bf_all <- stage2_loglik_R(d_all, st$pars)
  blk_mixed <- sapply(seq_len(dat$N1), function(n) { r <- dat$cs_start[n]:dat$cs_end[n]; any(dat$in_iv[r] == 0L | dat$in_s2[r] == 0L) })
  iv_mixed <- sapply(seq_len(dat$N1), function(n) { r <- dat$cs_start[n]:dat$cs_end[n]; any(dat$in_iv[r] == 0L) })
  s2_mixed <- sapply(seq_len(dat$N1), function(n) { r <- dat$cs_start[n]:dat$cs_end[n]; any(dat$in_s2[r] == 0L) })
  expect_true(all(bf$lse[!iv_mixed] == bf_all$lse[!iv_mixed])); expect_true(all(bf$lse[iv_mixed] < bf_all$lse[iv_mixed]))   # fewer rows in the IV sum
  expect_true(all(bf$loglik[!s2_mixed[dat$ev_sy]] == bf_all$loglik[!s2_mixed[dat$ev_sy]])); expect_true(all(bf$loglik[s2_mixed[dat$ev_sy]] > bf_all$loglik[s2_mixed[dat$ev_sy]]))
  ## probabilities sum to one over the STAGE-2 rows of each initiating state-year
  for (n in unique(dat$ev_sy)) { r <- dat$cs_start[n]:dat$cs_end[n]; r2 <- r[dat$in_s2[r] == 1L]
    e <- which(dat$ev_sy == n)[1]; lse2 <- bf$util[dat$ev_row[e]] - bf$loglik[e]
    expect_equal(sum(exp(bf$util[r2] - lse2)), 1, tolerance = 1e-10) }
  invisible(dat)
}
test_that("masked sets, M1 joint: Stan equals brute-force R", check_masked("m1_latent_additive", "M1b", 1L, FALSE))
test_that("masked sets, M1 independent: Stan equals brute-force R", check_masked("m1_latent_additive", "M1b", 0L, FALSE))
test_that("masked sets, M1 material-only: Stan equals brute-force R", check_masked("m1_latent_additive", "M1b", 0L, FALSE, D_mode = 1L))
test_that("masked sets, M2 joint (K=2): Stan equals brute-force R", check_masked("m2_ame", "M2b", 1L, TRUE))
test_that("masked sets, M2 relational-only (K=2): Stan equals brute-force R", check_masked("m2_ame", "M2b", 0L, TRUE, D_mode = 2L))

test_that("masked sets: Stan rejects an empty IV set, a chosen row outside the stage-2 set, and unresolved flags with sampled = 1", {
  sm <- skip_if_no_build("m1_latent_additive")
  sim <- simulate_targeting("M1b", S = 6L, Tt = 4L, seed = 123L, target_rate = 0.6)
  dat <- build_stan_data(sim$sy, sim$alts, sim$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"), joint = 1L, x_ref = rep(0, 4))
  p <- make_pars(dat, FALSE); st <- stan_inits(dat, p)
  ## rstan does not throw on a reject() in transformed data: it returns a stanfit in error mode (mode 2) with no draws
  rejected <- function(d) { f <- tryCatch(suppressMessages(suppressWarnings(sampling(sm, data = d, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(st$init), seed = 123, refresh = 0))),
                                          error = function(e) NULL); is.null(f) || f@mode == 2L }
  expect_false(rejected(dat))
  d1 <- dat; d1$in_iv[d1$cs_start[1]:d1$cs_end[1]] <- 0L                                                           # empty IV set
  expect_true(rejected(d1))
  d2 <- dat; d2$in_s2[d2$ev_row[1]] <- 0L                                                                           # chosen row outside the stage-2 set
  expect_true(rejected(d2))
  d3 <- dat; n <- d3$ev_sy[1]; r <- d3$cs_start[n]:d3$cs_end[n]; d3$in_s2[setdiff(r, d3$ev_row[d3$ev_sy == n])] <- 0L   # chosen-only stage-2 set
  expect_true(rejected(d3))
  d4 <- dat; d4$in_iv[1] <- 0L; d4$in_s2[1] <- 0L                                                                   # a row in neither set
  expect_true(rejected(d4))
  sa <- sample_alternatives(sim$alts, sim$ev, m = 2)
  ds <- build_stan_data(sim$sy, sa$alts, sa$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"), joint = 1L, sampled = 1L, n_full = sa$n_full, x_ref = rep(0, 4))
  expect_false(rejected(ds))
  ds$in_s2[1] <- 0L                                                                                                 # unresolved flag with sampled = 1
  expect_true(rejected(ds))
})

## ---- equivalence with the pre-review programs (commit 3a9e7d0, whole block = both sets) --------------------------------------
## The reference models are compiled from `git show 3a9e7d0:v2/stan/<file>` into V2_REF_BUILD as <model>_3a9e7d0.rds.
## (a) With every row in both sets the new program must give the SAME log density, gradient and generated quantities.
## (b) With masked sets the new program's stage 1 must equal the old program run on the IV rows alone (no events), its stage 2 the
##     old program run on the stage-2 rows alone, and its log density the old IV-only density plus the new stage-2 log-liks.
check_reference_3a9e7d0 <- function(model_name, family, joint, with_uv, K = 2L, D_mode = 3L) {
  ref <- file.path(ref_dir, paste0(model_name, "_3a9e7d0.rds"))
  if (ref_dir == "" || !file.exists(ref)) skip("reference build 3a9e7d0 not available (set V2_REF_BUILD)")
  sm_new <- skip_if_no_build(model_name); sm_ref <- readRDS(ref)
  expect_false(grepl("in_iv", sm_ref@model_code, fixed = TRUE)); expect_true(grepl("array[N2] int<lower=0, upper=1> in_iv;", sm_new@model_code, fixed = TRUE))
  Kk <- if (with_uv) K else NULL
  cols <- list(c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"))
  ## (a) identical sets
  sim <- simulate_targeting(family, S = 7L, Tt = 5L, K = K, seed = 123L, target_rate = 0.6)
  dat <- build_stan_data(sim$sy, sim$alts, sim$ev, cols[[1]], cols[[2]], cols[[3]], joint = joint, x_ref = rep(0, 4), K = Kk, D_mode = D_mode)
  expect_true(all(dat$in_iv == 1L & dat$in_s2 == 1L))
  p <- make_pars(dat, with_uv); init <- stan_inits(dat, p)$init
  new <- eval_fixed_draw(sm_new, dat, init); old <- eval_fixed_draw(sm_ref, dat, init)   # rstan ignores the extra in_iv / in_s2 data for the old program
  expect_equal(new$lp, old$lp, tolerance = 1e-12); expect_equal(new$grad, old$grad, tolerance = 1e-10)
  expect_equal(new$ll1, old$ll1, tolerance = 1e-12); expect_equal(new$ll2, old$ll2, tolerance = 1e-12); expect_equal(new$mu_state, old$mu_state, tolerance = 1e-12)
  ## (b) masked sets
  M <- mask_sets(sim)
  d_union <- build_stan_data(M$sy, M$alts, M$ev, cols[[1]], cols[[2]], cols[[3]], joint = joint, x_ref = rep(0, 4), K = Kk, anchor = dat$anchor, D_mode = D_mode)
  d_iv <- build_stan_data(M$sy, M$iv_only, M$ev_none, cols[[1]], cols[[2]], cols[[3]], joint = joint, x_ref = rep(0, 4), K = Kk, anchor = dat$anchor, D_mode = D_mode,
                          allow_orphan_initiations = TRUE, S = dat$S)
  d_s2 <- build_stan_data(M$sy, M$s2_only, M$ev_s2, cols[[1]], cols[[2]], cols[[3]], joint = joint, x_ref = rep(0, 4), K = Kk, anchor = dat$anchor, D_mode = D_mode, S = dat$S)
  expect_equal(d_iv$E, 0L); expect_equal(d_s2$E, d_union$E); expect_lt(d_iv$N2, d_union$N2); expect_lt(d_s2$N2, d_union$N2)
  nu <- eval_fixed_draw(sm_new, d_union, init); o_iv <- eval_fixed_draw(sm_ref, d_iv, init); o_s2 <- eval_fixed_draw(sm_ref, d_s2, init)
  expect_equal(nu$ll1, o_iv$ll1, tolerance = 1e-12)                                       # stage 1 = old program on the IV rows
  expect_equal(nu$ll2, o_s2$ll2, tolerance = 1e-12)                                       # stage 2 = old program on the stage-2 rows
  expect_equal(nu$lp, o_iv$lp + sum(nu$ll2), tolerance = 1e-10)                            # density = priors + stage 1 (IV rows) + stage 2
  ## and the whole-block evaluation of the same union data differs (the masks matter)
  expect_false(isTRUE(all.equal(nu$ll1, eval_fixed_draw(sm_ref, d_union, init)$ll1, tolerance = 1e-8)))
  invisible(TRUE)
}
test_that("reference 3a9e7d0: M1a identical when the sets coincide; masked sets decompose into the old program", check_reference_3a9e7d0("m1_latent_additive", "M1b", 0L, FALSE))
test_that("reference 3a9e7d0: M1b", check_reference_3a9e7d0("m1_latent_additive", "M1b", 1L, FALSE))
test_that("reference 3a9e7d0: M1m", check_reference_3a9e7d0("m1_latent_additive", "M1b", 0L, FALSE, D_mode = 1L))
test_that("reference 3a9e7d0: M2b (K=2)", check_reference_3a9e7d0("m2_ame", "M2b", 1L, TRUE))
test_that("reference 3a9e7d0: M2a (K=2)", check_reference_3a9e7d0("m2_ame", "M2b", 0L, TRUE))
test_that("reference 3a9e7d0: M2r (K=2)", check_reference_3a9e7d0("m2_ame", "M2b", 0L, TRUE, D_mode = 2L))
