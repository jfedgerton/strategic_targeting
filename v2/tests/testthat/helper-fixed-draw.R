## helper-fixed-draw.R -- shared by test_stage2_likelihood.R and test_reduce_sum.R (testthat sources helper-*.R first)
## Builds tiny Stan data sets and a fixed "natural" parameter draw translated into Stan inits.
suppressPackageStartupMessages(library(rstan))
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "00_helpers.R")); source(file.path(root, "sim", "simulate.R"))

skip_if_no_build <- function(m) {
  p <- file.path(root, "stan", "build", paste0(m, ".rds"))
  if (!file.exists(p)) skip(paste("compiled model missing:", p))
  readRDS(p)
}

tiny_data <- function(family, joint, K = NULL, D_mode = 3L) {
  set.seed(123)
  sim <- simulate_targeting(family, S = 6L, Tt = 4L, K = if (is.null(K)) 2L else K, seed = 123L, target_rate = 0.6)
  build_stan_data(sim$sy, sim$alts, sim$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"),
                  joint = joint, x_ref = rep(0, 4), K = K, D_mode = D_mode)
}

make_pars <- function(dat, with_uv) {
  set.seed(7)
  S <- dat$S
  D <- if (dat$D_mode == 3L) 2L else 1L
  rho <- if (dat$joint == 1 && D == 2L) 0.4 else 0
  ## mu0 is deliberately non-zero so that any mu reimplementation that omits the intercept fails the likelihood tests
  list(mu0 = c(0.7, -0.4)[seq_len(D)],
       beta = matrix(rnorm(dat$P * D, 0, 0.3), dat$P, D),
       a_i = matrix(rnorm(S * D), S, D), c_i = rnorm(S, 0, 0.5),
       b_j = { z <- rnorm(S); 0.6 * (z - mean(z)) }, z_b = NULL,
       gamma = c(0.5, -0.3), sigma_x = c(0.9, 1.2)[seq_len(D)], rho = rho,
       U = if (with_uv) matrix(rnorm(S * dat$K), S, dat$K) else NULL,
       V = if (with_uv) { v <- matrix(rnorm(S * dat$K), S, dat$K); sweep(v, 2, colMeans(v)) } else NULL)
}

## Translate "natural" parameters into the Stan parameterisation (non-centred etc.)
stan_inits <- function(dat, p) {
  S <- dat$S
  D <- if (dat$D_mode == 3L) 2L else 1L
  sd_re <- c(0.5, 1, 1)[seq_len(1 + D)]
  if (D == 1L) {
    L <- t(chol(matrix(c(1, 0.2, 0.2, 1), 2)))
  } else if (dat$joint == 1) {
    L <- t(chol(matrix(c(1, 0.2, 0.1, 0.2, 1, 0.3, 0.1, 0.3, 1), 3)))
  } else {
    r1 <- 0.2; q <- 0.1; r2 <- q * sqrt(1 - r1^2)
    L <- t(chol(matrix(c(1, r1, r2, r1, 1, 0, r2, 0, 1), 3)))
  }
  re <- cbind(p$c_i, p$a_i)                             # S x (1 + D) desired
  z_re <- solve(diag(sd_re) %*% L) %*% t(re)            # (1 + D) x S so that (diag(sd) L z)' = re
  z_b <- p$b_j / 0.6                                    # mean-zero already, so b_j = 0.6 * (z_b - mean) = b_j
  init <- list(alpha0 = -0.5, alpha = c(0.3, -0.2), tau = 0.4, mu0 = array(p$mu0, dim = D), beta = p$beta,
               sigma_x = array(p$sigma_x, dim = D), sd_re = sd_re, z_re = z_re, sigma_b = 0.6, z_b = z_b, gamma = p$gamma)
  if (D == 1L) { init$r_ca_1d <- array(0.2, 1) }
  else if (dat$joint == 1) { init$rho_x <- array(p$rho, 1); init$L_Om_joint <- array(L, c(1, 3, 3)) }
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


## unconstrain_pars() needs every parameter, including the zero-size ones (e.g. r_ca1 when joint = 1);
## sampling() accepts a partial init list, so complete it from the fit's declared dimensions
PARAM_BLOCK <- c("alpha0", "alpha", "tau", "mu0", "beta", "sigma_x", "rho_x", "sd_re", "L_Om_joint", "r_ca1", "q_ca2", "r_ca_1d",
                 "z_re", "sigma_b", "z_b", "gamma", "u_diag", "u_lower", "U_free", "V_raw")   # the parameters {} blocks of M1 / M2
complete_init <- function(fit, init) {
  for (nm in intersect(PARAM_BLOCK, fit@model_pars)) if (is.null(init[[nm]])) {
    d <- fit@par_dims[[nm]]
    if (length(d) == 0 || prod(d) != 0) stop("missing non-empty init: ", nm)
    init[[nm]] <- array(numeric(0), dim = d)
  }
  init
}

## ---- two choice sets per state-year (Codex review P1-1, 2026-10-24) ------------------------------------------------------
## mask_sets(): from a simulated cell, mark one non-chosen row of every other initiating state-year as IV only (in_s2 = 0)
## and one non-chosen row of the remaining initiating state-years as stage-2 only (in_iv = 0); non-initiating state-years keep
## identical sets, as in the real build. Returns the frames with flags plus the IV-only and stage-2-only sub-frames that the
## pre-review program (whole block = both sets) needs to reproduce stage 1 and stage 2 separately.
mask_sets <- function(sim) {
  alts <- sim$alts |> mutate(in_iv = 1L, in_s2 = 1L, .row = row_number())
  ev <- sim$ev
  init_sy <- sort(unique(ev$sy_row))
  stopifnot(length(init_sy) >= 2)
  for (k in seq_along(init_sy)) {
    n <- init_sy[k]
    rows <- alts$.row[alts$sy_row == n]
    free <- setdiff(rows, ev$alt_row[ev$sy_row == n])
    stopifnot(length(free) >= 2)                      # keeps >= 2 stage-2 rows and >= 1 IV row
    if (k %% 2 == 1) alts$in_s2[free[1]] <- 0L else alts$in_iv[free[1]] <- 0L
  }
  stopifnot(any(alts$in_iv == 0L), any(alts$in_s2 == 0L))
  iv_only <- alts |> filter(in_iv == 1L) |> select(-in_iv, -in_s2, -.row)
  s2_rows <- which(alts$in_s2 == 1L)
  s2_only <- alts[s2_rows, ] |> select(-in_iv, -in_s2, -.row)
  ev_s2 <- ev |> mutate(alt_row = match(alt_row, s2_rows))
  stopifnot(!anyNA(ev_s2$alt_row), all(s2_only$sy_row[ev_s2$alt_row] == ev_s2$sy_row))
  list(sy = sim$sy, alts = alts |> select(-.row), ev = ev, iv_only = iv_only, s2_only = s2_only, ev_s2 = ev_s2,
       ev_none = ev[0, , drop = FALSE])
}
## log density (with Jacobian), gradient and generated log-liks at a fixed init (Fixed_param run, then log_prob)
eval_fixed_draw <- function(sm, dat, init, threads = 1L) {
  rstan_options(threads_per_chain = threads)
  fit <- sampling(sm, data = dat, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(init), seed = 123, refresh = 0)
  up <- unconstrain_pars(fit, complete_init(fit, init))
  lp <- log_prob(fit, up, adjust_transform = TRUE, gradient = TRUE)
  out <- list(lp = as.numeric(lp), grad = as.numeric(attr(lp, "gradient")),
              ll1 = as.numeric(rstan::extract(fit, "log_lik_stage1")$log_lik_stage1[1, ]),
              mu_state = as.numeric(rstan::extract(fit, "mu_state")$mu_state[1, , ]))
  out$ll2 <- if (dat$E > 0) as.numeric(rstan::extract(fit, "log_lik_stage2")$log_lik_stage2[1, ]) else numeric(0)
  rstan_options(threads_per_chain = 1)
  out
}
