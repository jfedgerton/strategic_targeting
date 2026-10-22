## Stage-2 likelihood in Stan must match a brute-force R computation on a tiny example.
## Uses algorithm = "Fixed_param" so generated quantities are evaluated at supplied inits.
## Requires compiled models in v2/stan/build (Rscript v2/R/compile_models.R).
suppressPackageStartupMessages({ library(testthat); library(rstan) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "00_helpers.R"))
source(file.path(root, "sim", "simulate.R"))

skip_if_no_build <- function(m) {
  p <- file.path(root, "stan", "build", paste0(m, ".rds"))
  if (!file.exists(p)) skip(paste("compiled model missing:", p))
  readRDS(p)
}

tiny_data <- function(family, joint, K = NULL) {
  set.seed(123)
  sim <- simulate_targeting(family, S = 6L, Tt = 4L, K = if (is.null(K)) 2L else K, seed = 123L, target_rate = 0.6)
  build_stan_data(sim$sy, sim$alts, sim$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"),
                  joint = joint, x_ref = rep(0, 4), K = K)
}

make_pars <- function(dat, with_uv) {
  set.seed(7)
  S <- dat$S
  rho <- if (dat$joint == 1) 0.4 else 0
  list(beta = matrix(rnorm(dat$P * 2, 0, 0.3), dat$P, 2),
       a_i = matrix(rnorm(S * 2), S, 2), c_i = rnorm(S, 0, 0.5),
       b_j = { z <- rnorm(S); 0.6 * (z - mean(z)) }, z_b = NULL,
       gamma = c(0.5, -0.3), sigma_x = c(0.9, 1.2), rho = rho,
       U = if (with_uv) matrix(rnorm(S * dat$K), S, dat$K) else NULL,
       V = if (with_uv) { v <- matrix(rnorm(S * dat$K), S, dat$K); sweep(v, 2, colMeans(v)) } else NULL)
}

## Translate "natural" parameters into the Stan parameterisation (non-centred etc.)
stan_inits <- function(dat, p) {
  S <- dat$S
  sd_re <- c(0.5, 1, 1)
  if (dat$joint == 1) {
    L <- t(chol(matrix(c(1, 0.2, 0.1, 0.2, 1, 0.3, 0.1, 0.3, 1), 3)))
  } else {
    r1 <- 0.2; q <- 0.1; r2 <- q * sqrt(1 - r1^2)
    L <- t(chol(matrix(c(1, r1, r2, r1, 1, 0, r2, 0, 1), 3)))
  }
  re <- cbind(p$c_i, p$a_i)                             # S x 3 desired
  z_re <- solve(diag(sd_re) %*% L) %*% t(re)            # 3 x S so that (diag(sd) L z)' = re
  z_b <- p$b_j / 0.6                                    # mean-zero already, so b_j = 0.6 * (z_b - mean) = b_j
  init <- list(alpha0 = -0.5, alpha = c(0.3, -0.2), tau = 0.4, beta = p$beta,
               sigma_x = p$sigma_x, sd_re = sd_re, z_re = z_re, sigma_b = 0.6, z_b = z_b, gamma = p$gamma)
  if (dat$joint == 1) { init$rho_x <- array(p$rho, 1); init$L_Om_joint <- array(L, c(1, 3, 3)) }
  else { init$r_ca1 <- array(0.2, 1); init$q_ca2 <- array(0.1, 1) }
  if (!is.null(p$U)) {
    K <- dat$K; anchor <- dat$anchor
    ## force the anchor block of the supplied U to be lower-triangular with positive diagonal
    for (k in seq_len(K)) { for (l in seq_len(K)) { if (l > k) p$U[anchor[k], l] <- 0 }; p$U[anchor[k], k] <- abs(p$U[anchor[k], k]) + 0.1 }
    init$u_diag <- array(sapply(seq_len(K), function(k) p$U[anchor[k], k]), dim = K)
    init$u_lower <- array(if (K > 1) unlist(lapply(2:K, function(k) p$U[anchor[k], 1:(k - 1)])) else numeric(0),
                          dim = K * (K - 1) / 2)
    ## U columns sum to zero: the last non-anchor row is determined by the others
    nonanchor <- setdiff(seq_len(S), anchor)
    last <- nonanchor[length(nonanchor)]
    p$U[last, ] <- -colSums(p$U[-last, , drop = FALSE])
    init$U_free <- p$U[nonanchor[-length(nonanchor)], , drop = FALSE]
    init$V_raw <- p$V
  }
  list(init = init, pars = p)
}

check_model <- function(model_name, family, joint, with_uv, K = 2L) {
  sm <- skip_if_no_build(model_name)
  dat <- tiny_data(family, joint, K = if (with_uv) K else NULL)
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
  ## mu_state at x_ref = 0 equals a_i
  mu <- rstan::extract(fit, "mu_state")$mu_state[1, , ]
  expect_equal(unname(mu), unname(st$pars$a_i), tolerance = 1e-6)
  ## Sigma assembled correctly
  Sg <- rstan::extract(fit, "Sigma")$Sigma[1, , ]
  expect_equal(Sg[1, 2], st$pars$rho * prod(st$pars$sigma_x), tolerance = 1e-8)
  if (joint == 0) expect_equal(Sg[1, 2], 0)
  invisible(TRUE)
}

test_that("M1 joint: Stan stage-2 log-lik equals brute-force R", check_model("m1_latent_additive", "M1b", 1L, FALSE))
test_that("M1 independent: Stan stage-2 log-lik equals brute-force R", check_model("m1_latent_additive", "M1b", 0L, FALSE))
test_that("M2 joint (K=2): Stan stage-2 log-lik equals brute-force R", check_model("m2_ame", "M2b", 1L, TRUE))
test_that("M2 independent (K=2): Stan stage-2 log-lik equals brute-force R", check_model("m2_ame", "M2b", 0L, TRUE))
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
