## G as data in m1_context.stan (2026-10-31) and the role-specific-rho program m1_context_rho.stan (model M1cr, PI request 2026-10-31).
## (a) G = 3 reproduces the pre-2026-10-31 program: the source of commit 5153f8e is compiled here and its log density, gradient and
##     generated log-likelihoods are compared with the current program at fixed parameter points of a toy cell;
## (b) m1_context_rho with rho_g = (rho, rho, rho) is m1_context plus (G - 1) Beta(2.5, 2.5) prior terms, and with distinct rho_g its
##     log-likelihoods equal the brute-force R likelihood (IV with the reference group's rho_g[1]);
## (c) spec routing of M1cr.
suppressPackageStartupMessages({ library(testthat); library(dplyr); library(rstan) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "typology.R")); source(file.path(root, "R", "context_model.R")); source(file.path(root, "sim", "simulate.R"))

toy_ctx <- function(seed = 123L) {
  set.seed(seed)
  sim <- simulate_targeting("M1b", S = 6L, Tt = 6L, K = 2L, seed = seed, target_rate = 0.6)
  dat <- build_stan_data(sim$sy, sim$alts, sim$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"), joint = 1L, x_ref = rep(0.1, 4), D_mode = 3L)
  dat$G <- 3L; dat$ev_grp <- as.integer((seq_len(dat$E) - 1L) %% 3L)
  dat
}
## one random parameter point, as named lists for the old (delta_lead / delta_join) and the new (delta) parameterisation
rand_point <- function(dat, seed) {
  set.seed(seed); S <- dat$S
  L <- t(chol(matrix(c(1, 0.2, 0.1, 0.2, 1, 0.3, 0.1, 0.3, 1), 3)))
  com <- list(alpha0 = rnorm(1), alpha = rnorm(dat$P1, 0, 0.5), tau = rnorm(1, 0, 0.5), mu0 = rnorm(2), beta = matrix(rnorm(dat$P * 2, 0, 0.3), dat$P, 2),
              sigma_g = matrix(exp(rnorm(6, 0, 0.3)), 3, 2), rho_x = runif(1, -0.8, 0.8), sd_re = exp(rnorm(3, 0, 0.3)), L_Om = L,
              z_re = matrix(rnorm(3 * S), 3, S), sigma_b = exp(rnorm(1, 0, 0.3)), z_b = rnorm(S), gamma = rnorm(dat$Pz, 0, 0.5))
  d <- rbind(rnorm(2, 0, 0.7), rnorm(2, 0, 0.7))
  list(old = c(com, list(delta_lead = d[1, ], delta_join = d[2, ])), new = c(com, list(delta = d)))
}
gq <- function(fit, par) as.numeric(rstan::extract(fit, par)[[par]])          # one draw: every element in storage order

test_that("m1_context with G = 3 gives the log density, gradient and generated quantities of the pre-2026-10-31 program (commit 5153f8e)", {
  smc <- skip_if_no_build("m1_context")
  old_src <- tryCatch(suppressWarnings(system2("git", c("-C", root, "show", "5153f8e:v2/stan/m1_context.stan"), stdout = TRUE, stderr = FALSE)), error = function(e) NULL)
  if (is.null(old_src) || !length(old_src) || !is.null(attr(old_src, "status"))) skip("git history with commit 5153f8e not available")
  expect_true(any(grepl("row_vector[2] delta_lead;", old_src, fixed = TRUE))); expect_true(any(grepl("int<lower=1> P;", old_src, fixed = TRUE)))
  sm_old <- stan_model(model_code = paste(old_src, collapse = "\n"), model_name = "m1_context_5153f8e", auto_write = FALSE)
  dat <- toy_ctx(); dat_old <- dat; dat_old$G <- NULL
  rstan_options(threads_per_chain = 1)
  p1 <- rand_point(dat, 1L)
  f_old <- sampling(sm_old, data = dat_old, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(p1$old), seed = 123, refresh = 0)
  f_new <- sampling(smc, data = dat, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(p1$new), seed = 123, refresh = 0)
  expect_equal(f_old@mode, 0L); expect_equal(f_new@mode, 0L)
  for (s in 1:6) {
    p <- rand_point(dat, 100L + s)
    u_old <- unconstrain_pars(f_old, p$old); u_new <- unconstrain_pars(f_new, p$new)
    expect_identical(u_old, u_new)                                                 # same unconstrained layout (delta[1] = delta_lead, delta[2] = delta_join)
    for (adj in c(TRUE, FALSE)) expect_identical(log_prob(f_old, u_old, adjust_transform = adj), log_prob(f_new, u_new, adjust_transform = adj))   # bit-identical: same terms, same order
    expect_equal(grad_log_prob(f_old, u_old), grad_log_prob(f_new, u_new), tolerance = 1e-12, check.attributes = FALSE)
  }
  ## generated quantities at the first point: the log-liks and ideal points are identical, mu0_g carries mu0_lead / mu0_join
  for (par in c("log_lik_stage1", "log_lik_stage2", "mu_state", "mu_event", "Sigma", "Sigma_g", "Omega")) expect_identical(gq(f_old, par), gq(f_new, par), info = par)
  m0g <- rstan::extract(f_new, "mu0_g")$mu0_g[1, , ]
  expect_identical(m0g[1, ], as.numeric(p1$new$mu0)); expect_identical(as.numeric(m0g[2, ]), gq(f_old, "mu0_lead")); expect_identical(as.numeric(m0g[3, ]), gq(f_old, "mu0_join"))
  ## the same with no X at all would need the new program only (P = 0 is refused by the old one); the threads do not change the new density
  dat2 <- dat; dat2$grainsize <- 2L; rstan_options(threads_per_chain = 2)
  f2 <- sampling(smc, data = dat2, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(p1$new), seed = 123, refresh = 0)
  expect_equal(log_prob(f2, unconstrain_pars(f2, p1$new)), log_prob(f_new, unconstrain_pars(f_new, p1$new)), tolerance = 1e-10)
  rstan_options(threads_per_chain = 1)
})

test_that("m1_context with G = 2 and G = 4: the brute-force R likelihood matches the generated log-liks (the group count is data)", {
  smc <- skip_if_no_build("m1_context")
  rstan_options(threads_per_chain = 1)
  for (G in c(2L, 4L)) {
    dat <- toy_ctx(); dat$G <- G; dat$ev_grp <- as.integer((seq_len(dat$E) - 1L) %% G)
    set.seed(G); S <- dat$S; L <- diag(3)
    d <- matrix(rnorm(2 * (G - 1), 0, 0.6), G - 1, 2); sg <- matrix(exp(rnorm(2 * G, 0, 0.3)), G, 2)
    a_c <- matrix(rnorm(3 * S, 0, 0.5), 3, S); z_b <- rnorm(S)
    init <- list(alpha0 = -0.3, alpha = c(0.2, -0.1), tau = 0.3, mu0 = c(0.4, -0.2), delta = d, beta = matrix(0.1, dat$P, 2), sigma_g = sg, rho_x = 0.3,
                 sd_re = c(1, 1, 1), L_Om = L, z_re = a_c, sigma_b = 0.7, z_b = z_b, gamma = c(0.4, -0.2))
    fit <- sampling(smc, data = dat, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(init), seed = 123, refresh = 0)
    expect_equal(fit@mode, 0L); expect_equal(fit@par_dims$delta, c(G - 1L, 2L)); expect_equal(fit@par_dims$mu0_g, c(G, 2L))
    R <- stage2_loglik_R_context(dat, list(mu0 = init$mu0, delta = d, beta = init$beta, a_i = t(a_c[2:3, ]), b_j = 0.7 * (z_b - mean(z_b)), gamma = init$gamma, sigma_g = sg, rho = 0.3))
    expect_equal(gq(fit, "log_lik_stage2"), R$loglik, tolerance = 1e-8, info = paste("G", G))
    eta <- init$alpha0 + as.numeric(dat$W %*% init$alpha) + a_c[1, dat$sy_state] + init$tau * R$lse
    expect_equal(gq(fit, "log_lik_stage1"), dbinom(dat$y1, 1, plogis(eta), log = TRUE), tolerance = 1e-8, info = paste("G", G))
  }
})

test_that("m1_context_rho: equal rho_g is m1_context plus (G - 1) prior terms; distinct rho_g match the brute-force likelihood (IV with rho_g[1])", {
  smc <- skip_if_no_build("m1_context"); smr <- skip_if_no_build("m1_context_rho")
  src <- readLines(file.path(root, "stan", "m1_context_rho.stan")); src_c <- readLines(file.path(root, "stan", "m1_context.stan"))
  expect_true(any(grepl("vector<lower=-1, upper=1>[G] rho_g;", src, fixed = TRUE))); expect_false(any(grepl("rho_x", src[grep("^functions", src):length(src)], fixed = TRUE)))
  expect_true(any(grepl("for (g in 1:G) target += beta_lpdf((rho_g[g] + 1) / 2 | 2.5, 2.5);", src, fixed = TRUE)))
  ## the programs differ only in the correlation: substituting rho_x back for rho_g gives m1_context.stan's body line for line
  body <- function(x) x[grep("^functions", x):length(x)]
  back <- body(src)
  back <- sub("vector<lower=-1, upper=1>[G] rho_g;          // the correlation of Sigma_g per group: rho_g[g + 1] for group g (1 = the reference)", "real<lower=-1, upper=1> rho_x;               // the shared correlation of Sigma_g", back, fixed = TRUE)
  back <- sub("for (g in 1:G) target += beta_lpdf((rho_g[g] + 1) / 2 | 2.5, 2.5);", "target += beta_lpdf((rho_x + 1) / 2 | 2.5, 2.5);", back, fixed = TRUE)
  back <- sub("vector rho, vector b_j", "real rho, vector b_j", back, fixed = TRUE)
  back <- gsub("rho[1]", "rho", back, fixed = TRUE); back <- gsub("rho[g]", "rho", back, fixed = TRUE)
  back <- gsub("rho_g[ev_grp[ev] + 1]", "rho_x", back, fixed = TRUE); back <- gsub("rho_g[1]", "rho_x", back, fixed = TRUE); back <- gsub("rho_g[g]", "rho_x", back, fixed = TRUE)
  back <- gsub("rho_g", "rho_x", back, fixed = TRUE); back <- gsub("m1_context_rho:", "m1_context:", back, fixed = TRUE)
  expect_identical(back, body(src_c))
  dat <- toy_ctx(); rstan_options(threads_per_chain = 1)
  p <- rand_point(dat, 7L)$new
  pr <- p; pr$rho_x <- NULL; pr$rho_g <- rep(p$rho_x, 3)
  fc <- sampling(smc, data = dat, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(p), seed = 123, refresh = 0)
  fr <- sampling(smr, data = dat, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(pr), seed = 123, refresh = 0)
  expect_equal(fr@mode, 0L); expect_equal(fr@par_dims$rho_g, 3L)
  for (par in c("log_lik_stage1", "log_lik_stage2", "mu_state", "mu_event", "Sigma_g", "mu0_g")) expect_equal(gq(fr, par), gq(fc, par), tolerance = 1e-12, info = par)
  lp_c <- log_prob(fc, unconstrain_pars(fc, p), adjust_transform = FALSE); lp_r <- log_prob(fr, unconstrain_pars(fr, pr), adjust_transform = FALSE)
  expect_equal(lp_r - lp_c, 2 * dbeta((p$rho_x + 1) / 2, 2.5, 2.5, log = TRUE), tolerance = 1e-8)
  ## distinct rho_g: the per-group correlation reaches each group's utility, the IV uses rho_g[1]
  pr2 <- pr; pr2$rho_g <- c(0.5, -0.6, 0.2)
  fr2 <- sampling(smr, data = dat, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(pr2), seed = 123, refresh = 0)
  re <- t(diag(p$sd_re) %*% p$L_Om %*% p$z_re)
  R <- stage2_loglik_R_context(dat, list(mu0 = p$mu0, delta = p$delta, beta = p$beta, a_i = re[, 2:3], b_j = p$sigma_b * (p$z_b - mean(p$z_b)), gamma = p$gamma, sigma_g = p$sigma_g, rho = pr2$rho_g))
  expect_equal(gq(fr2, "log_lik_stage2"), R$loglik, tolerance = 1e-8)
  eta <- p$alpha0 + as.numeric(dat$W %*% p$alpha) + re[dat$sy_state, 1] + p$tau * R$lse
  expect_equal(gq(fr2, "log_lik_stage1"), dbinom(dat$y1, 1, plogis(eta), log = TRUE), tolerance = 1e-8)
  Sg <- rstan::extract(fr2, "Sigma_g")$Sigma_g[1, , , ]
  for (g in 1:3) expect_equal(Sg[g, 1, 2], pr2$rho_g[g] * prod(p$sigma_g[g, ]), tolerance = 1e-12)
  ## changing a non-reference rho_g leaves stage 1 untouched; changing rho_g[1] moves it
  pr3 <- pr2; pr3$rho_g[3] <- -0.1
  fr3 <- sampling(smr, data = dat, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(pr3), seed = 123, refresh = 0)
  expect_identical(gq(fr3, "log_lik_stage1"), gq(fr2, "log_lik_stage1")); expect_false(isTRUE(all.equal(gq(fr3, "log_lik_stage2"), gq(fr2, "log_lik_stage2"))))
  pr4 <- pr2; pr4$rho_g[1] <- 0.1
  fr4 <- sampling(smr, data = dat, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(pr4), seed = 123, refresh = 0)
  expect_false(isTRUE(all.equal(gq(fr4, "log_lik_stage1"), gq(fr2, "log_lik_stage1"))))
})

test_that("spec: M1cr routes to m1_context_rho, is a joint both-axes context model, and is accepted by the conventional comparison", {
  expect_true("M1cr" %in% CONTEXT_MODELS); expect_true("M1cr" %in% ALL_MODELS); expect_false("M1cr" %in% MODELS)
  expect_equal(stan_program("M1cr"), "m1_context_rho"); expect_equal(stan_program("M1c"), "m1_context"); expect_true("m1_context_rho" %in% STAN_PROGRAMS)
  expect_equal(model_joint("M1cr"), 1L); expect_equal(model_D_mode("M1cr"), 3L); expect_equal(model_family("M1cr"), "M1")
  for (m in MODELS) expect_equal(model_joint(m), as.integer(m %in% c("M1b", "M2b")), info = m)   # the old regex's values for every main-grid model
  expect_equal(cell_name("fmid", 3L, "M1cr", parse_spec(character(0))), "fmid_lag3_all_pr_M1cr")
  expect_equal(cell_name("fmid", 3L, "M1cr", parse_spec("--comparison=conventional")), "fmid_lag3_all_pr_M1cr_conv")
  expect_error(cell_name("fmid", 3L, "M1cr", parse_spec("--rho_prior=flat")), "prior switches apply to the M1 m / r / a / b programs only")
})
