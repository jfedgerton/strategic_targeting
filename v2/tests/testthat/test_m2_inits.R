## M2 initial values from the matching M1 fit (PI decision 2026-10-25): tiny synthetic M1 and M2 fits.
## Requires the compiled models (Rscript v2/R/compile_models.R); simulate_targeting(), build_stan_data(), skip_if_no_build() come
## from helper-fixed-draw.R.
suppressPackageStartupMessages({ library(testthat); library(rstan) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)

EXCL <- c("z_re", "z_b", "V_raw", "U_free", "L_Om_joint", "L_Om", "re")
tiny_pair <- function(joint, D_mode = 3L, K = 2L) {
  set.seed(123)
  sim <- simulate_targeting("M1b", S = 6L, Tt = 4L, K = K, seed = 123L, target_rate = 0.6)
  d1 <- build_stan_data(sim$sy, sim$alts, sim$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"), joint = joint, x_ref = rep(0, 4), D_mode = D_mode)
  d2 <- build_stan_data(sim$sy, sim$alts, sim$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"), joint = joint, x_ref = rep(0, 4), D_mode = D_mode, K = K)
  list(d1 = d1, d2 = d2)
}
tiny_m1_fit <- function(sm, dat) {
  rstan_options(threads_per_chain = 1)
  fit <- sampling(sm, data = dat, chains = 1, iter = 60, warmup = 30, seed = 123, refresh = 0, pars = EXCL, include = FALSE)
  list(fit = fit, dat = dat, cell = "tiny_M1", W_cols = c("W1", "W2"), X_cols = paste0("X", 1:4), Z_cols = c("Z1", "Z2"),
       args = list(pilot = FALSE, chains = 1L, iter = 60L, warmup = 30L))
}
## the M2 parameters with a non-zero size for this (D, joint): everything the init list must name
m2_param_dims <- function(fit2) { d <- fit2@par_dims[fit2@model_pars]; d[vapply(d, function(x) length(x) == 0 || prod(x) > 0, logical(1))] }

for (case in list(list(joint = 1L, D_mode = 3L, label = "M2b from M1b"), list(joint = 0L, D_mode = 3L, label = "M2a from M1a"), list(joint = 0L, D_mode = 1L, label = "M2m from M1m"))) {
  test_that(sprintf("%s: the init list names every M2 parameter with the right dimensions, copies the M1 means and starts the AME factors near 0", case$label), {
    sm1 <- skip_if_no_build("m1_latent_additive"); sm2 <- skip_if_no_build("m2_ame")
    tp <- tiny_pair(case$joint, case$D_mode); o1 <- tiny_m1_fit(sm1, tp$d1)
    D <- if (case$D_mode == 3L) 2L else 1L; S <- tp$d2$S; K <- tp$d2$K
    ## (a) no jitter: the init reproduces the M1 posterior means exactly
    in0 <- m1_inits_for_m2(o1, tp$d2, chains = 1, seed = 123, jitter_sd = 0)
    expect_length(in0, 1L); expect_equal(attr(in0, "n_copied"), length(in0[[1]]) - 4L)        # everything but the four AME blocks is copied
    ex <- rstan::extract(o1$fit, pars = c("alpha0", "alpha", "tau", "mu0", "beta", "gamma", "sigma_x", "sd_re", "sigma_b", "c_i", "a_i", "b_j", "R_Om"))
    expect_equal(in0[[1]]$alpha0, mean(ex$alpha0)); expect_equal(as.numeric(in0[[1]]$alpha), colMeans(ex$alpha)); expect_equal(as.numeric(in0[[1]]$gamma), colMeans(ex$gamma))
    expect_equal(as.numeric(in0[[1]]$mu0), as.numeric(apply(ex$mu0, 2, mean))); expect_equal(in0[[1]]$beta, matrix(apply(ex$beta, c(2, 3), mean), tp$d2$P, D))
    expect_equal(as.numeric(in0[[1]]$sd_re), colMeans(ex$sd_re)); expect_equal(in0[[1]]$sigma_b, mean(ex$sigma_b))
    expect_true(all(vapply(in0[[1]][c("alpha", "mu0", "sigma_x", "sd_re", "gamma", "z_b", "u_diag")], function(v) !is.null(dim(v)), logical(1))))   # 1-d arrays, so length-1 vectors keep their dimension
    ## a Fixed_param draw of the M2 program AT the init returns the init: the transformed (c_i, a_i, b_j) equal the M1 posterior means
    f2 <- sampling(sm2, data = tp$d2, init = in0, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, seed = 1, refresh = 0)
    expect_equal(f2@mode, 0L)                                                                              # the inits were accepted
    e2 <- rstan::extract(f2, pars = c("c_i", "a_i", "b_j", "R_Om", "U", "V", "u_diag", "alpha", "mu0"), permuted = FALSE)
    g <- function(nm) { v <- e2[1, 1, grepl(paste0("^", nm, "\\["), dimnames(e2)[[3]])]; as.numeric(v) }
    expect_equal(g("c_i"), colMeans(ex$c_i), tolerance = 1e-6)
    expect_equal(g("a_i"), as.numeric(apply(ex$a_i, c(2, 3), mean)), tolerance = 1e-6)
    expect_equal(g("b_j"), colMeans(ex$b_j), tolerance = 1e-6)
    R_m1 <- apply(ex$R_Om, c(2, 3), mean); if (case$joint == 0L && D == 2L) R_m1[2, 3] <- R_m1[3, 2] <- 0
    expect_equal(g("R_Om"), as.numeric(R_m1), tolerance = 1e-6)
    expect_equal(g("alpha"), colMeans(ex$alpha), tolerance = 1e-8); expect_equal(g("mu0"), as.numeric(apply(ex$mu0, 2, mean)), tolerance = 1e-8)
    ## (b) names and dimensions against the M2 program itself (every non-empty parameter is named, with its dims)
    need <- m2_param_dims(f2); need <- need[!names(need) %in% c("L_Om", "rho", "r_c_a1", "r_c_a2", "re", "c_i", "a_i", "b_j", "U", "V")]   # transformed parameters
    need <- need[!grepl("^(log_lik|mu_state|mu_event|Sigma|Omega|R_Om|lp__)", names(need))]                                                 # generated quantities
    expect_setequal(names(in0[[1]]), names(need))
    for (nm in names(need)) {
      v <- in0[[1]][[nm]]; dv <- if (is.null(dim(v))) length(v) else dim(v)
      expect_equal(as.integer(dv), as.integer(if (length(need[[nm]]) == 0) 1L else need[[nm]]), info = nm)
    }
    ## (c) the AME factors start near 0 and u_diag is positive
    expect_true(all(in0[[1]]$u_diag > 0 & in0[[1]]$u_diag <= 0.1)); expect_true(all(abs(in0[[1]]$U_free) <= 0.1)); expect_true(all(abs(in0[[1]]$V_raw) <= 0.1))
    expect_true(all(abs(g("V")) <= 0.2)); expect_true(all(abs(g("u_diag")) <= 0.1))
    expect_equal(dim(in0[[1]]$U_free), c(S - K - 1L, K)); expect_equal(dim(in0[[1]]$V_raw), c(S, K)); expect_length(in0[[1]]$u_lower, K * (K - 1L) / 2L)
    ## (d) jittered chains differ from each other and from the means, by a small amount; the jitter is seeded per chain
    inj <- m1_inits_for_m2(o1, tp$d2, chains = 3, seed = 123, jitter_sd = 0.05)
    expect_length(inj, 3L)
    expect_false(identical(inj[[1]]$alpha, inj[[2]]$alpha)); expect_false(identical(inj[[1]]$z_re, inj[[3]]$z_re)); expect_false(identical(inj[[1]]$V_raw, inj[[2]]$V_raw))
    expect_lt(max(abs(inj[[1]]$alpha - in0[[1]]$alpha)), 0.3); expect_lt(max(abs(inj[[2]]$gamma - in0[[1]]$gamma)), 0.3)
    expect_true(all(inj[[1]]$sd_re > 0), all(inj[[2]]$sigma_x > 0))
    expect_identical(m1_inits_for_m2(o1, tp$d2, chains = 2, seed = 123)[[2]], inj[[2]])                       # reproducible
    ## the jittered inits are accepted by the M2 program as well
    f3 <- sampling(sm2, data = tp$d2, init = inj[1:2], algorithm = "Fixed_param", chains = 2, iter = 1, warmup = 0, seed = 1, refresh = 0)
    expect_equal(f3@mode, 0L); expect_equal(dim(rstan::extract(f3, pars = "U", permuted = FALSE))[2], 2L)
  })
}

test_that("dimension mismatches between the M1 fit and the M2 data are errors, and the M1 fit of a cell is located by name", {
  sm1 <- skip_if_no_build("m1_latent_additive")
  tp <- tiny_pair(1L, 3L); o1 <- tiny_m1_fit(sm1, tp$d1)
  d_bad <- tp$d2; d_bad$P1 <- 3L; expect_error(m1_inits_for_m2(o1, d_bad, chains = 1), "P1 differs")
  d_bad <- tp$d2; d_bad$joint <- 0L; expect_error(m1_inits_for_m2(o1, d_bad, chains = 1), "joint differs")
  d_bad <- tp$d2; d_bad$K <- NULL; expect_error(m1_inits_for_m2(o1, d_bad, chains = 1))
  td <- tempfile("out"); dir.create(td); spec <- parse_spec(character(0))
  expect_true(is.na(locate_m1_fit("mid", 3L, "M2b", spec, td)))
  file.create(file.path(td, "fit_mid_lag3_all_pr_M1b.rds"))
  expect_equal(basename(locate_m1_fit("mid", 3L, "M2b", spec, td)), "fit_mid_lag3_all_pr_M1b.rds")
  expect_true(is.na(locate_m1_fit("mid", 3L, "M2b", parse_spec("--ww2=drop"), td)))                # the switches must match
  file.create(file.path(td, "fit_mid_lag3_all_pr_M1a_noww2.rds"))
  expect_equal(basename(locate_m1_fit("mid", 3L, "M2a", parse_spec("--ww2=drop"), td)), "fit_mid_lag3_all_pr_M1a_noww2.rds")
  expect_error(locate_m1_fit("mid", 3L, "M1b", spec, td))
})

test_that("audit fix 17: a binding correlation clamp in the M2 init construction is logged", {
  sm1 <- skip_if_no_build("m1_latent_additive")
  tp <- tiny_pair(1L, 3L); o1 <- tiny_m1_fit(sm1, tp$d1)
  expect_silent_log <- function(expr) expect_false(any(grepl("clamp binds", capture.output(expr))))
  expect_silent_log(m1_inits_for_m2(o1, tp$d2, chains = 2, seed = 123, jitter_sd = 0.05))                 # the default jitter never binds here
  out <- capture.output(res <- m1_inits_for_m2(o1, tp$d2, chains = 2, seed = 123, jitter_sd = 4))          # N(0, 0.8) on the correlations: the clamp binds (seed 123)
  expect_true(any(grepl("WARNING: m1_inits_for_m2 chain [12]: the .* init clamp binds", out)))
  expect_true(any(grepl("correlation init\\(s\\) clamped to \\+/-0.95 over 2 chains", out)))
  expect_true(attr(res, "n_clamped") >= length(grep("clamp binds", out)))                                  # one line per binding call, >= 1 value each
  expect_equal(attr(m1_inits_for_m2(o1, tp$d2, chains = 1, seed = 123, jitter_sd = 0), "n_clamped"), 0L)
  expect_true(all(abs(res[[1]]$rho_x) <= 0.95))
})
