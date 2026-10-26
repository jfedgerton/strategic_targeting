## The vectorised posterior-predictive scorer used by the holdout and permutation tests must agree
## with the brute-force stage-2 likelihood, and permute_dR must be a within-choice-set permutation.
suppressPackageStartupMessages({ library(testthat); library(dplyr) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "spec.R")); source(file.path(root, "R", "helpers_data.R")); source(file.path(root, "R", "helpers_model.R")); source(file.path(root, "R", "typology.R"))
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
  ## B6: with identical draws the per-draw hit indicators and the posterior-mean-probability ranks agree with the single-draw rank
  expect_equal(sc$hit1_pdraw, as.numeric(rk == 1)); expect_equal(sc$hit5_pdraw, as.numeric(rk <= 5))
  expect_equal(sc$rank_pmean, as.numeric(rk)); expect_equal(sc$hit1_pmean, rk == 1)
  expect_equal(sc$w, rep(1, dat$E))                                        # no weights in the simulator's ev: 1
  pars <<- pars
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

test_that("B6: per-draw hit indicators differ from thresholding the mean rank; the two labelled variants are both reported", {
  ## two draws that disagree: in draw 1 the chosen target ranks first, in draw 2 it ranks third (mean rank 2 -> the old hit1 = mean rank <= 1.5
  ## was FALSE and hit5 TRUE; the per-draw indicator gives P(top-1) = 0.5)
  dr2 <- draws
  dr2$b_j[2, ] <- dr2$b_j[2, ] * 0                                                # change draw 2 only
  dr2$b_j[2, dat$alt_target[dat$ev_row[1]]] <- -30                                 # the first event's chosen target becomes very unattractive in draw 2
  sc <- stage2_scores(dat, dr2)
  bf <- stage2_loglik_R(dat, pars)                                                 # draw 1 = the original parameters
  rows1 <- dat$cs_start[dat$ev_sy[1]]:dat$cs_end[dat$ev_sy[1]]; r1 <- rank(-bf$util[rows1])[match(dat$ev_row[1], rows1)]
  expect_equal(sc$hit1_pdraw[1], 0.5 * (r1 == 1))                                # draw 1 decides; draw 2 ranks the chosen target last
  expect_false(sc$mean_rank[1] <= 1.5)                                            # the OLD rule (mean rank <= 1.5) would say "no hit" even when draw 1 hit
  expect_equal(sc$hit5_pdraw[1], 0.5 * (r1 <= 5))
  expect_true(all(sc$hit1_pdraw >= 0 & sc$hit1_pdraw <= 1)); expect_true(all(sc$hit5_pdraw >= sc$hit1_pdraw))
  expect_true(all(sc$rank_pmean >= 1)); expect_true(is.logical(sc$hit1_pmean))
  ## the posterior-mean-probability rank of the first event is worse than or equal to draw 1's rank 1 and better than the collapsed draw-2 rank
  expect_true(sc$rank_pmean[1] >= 1)
  ## with event weights the scorer carries them and the weighted log score differs from the plain one
  dw <- dat; dw$ev_w <- seq(0.2, 1, length.out = dat$E)
  scw <- stage2_scores(dw, draws)
  expect_equal(scw$w, dw$ev_w); expect_equal(scw$log_p_chosen, stage2_scores(dat, draws)$log_p_chosen)
  expect_false(isTRUE(all.equal(sum(scw$w * scw$log_p_chosen), sum(scw$log_p_chosen))))
})

test_that("B3: stage2_scores runs for every model code (1-D codes keep beta and a_i as 3-d arrays in rstan)", {
  p1 <- file.path(root, "stan", "build", "m1_latent_additive.rds"); p2 <- file.path(root, "stan", "build", "m2_ame.rds")
  if (!file.exists(p1) || !file.exists(p2)) skip("compiled models missing")
  suppressPackageStartupMessages(library(rstan))
  sms <- list(M1 = readRDS(p1), M2 = readRDS(p2))
  codes <- c("M1m", "M1r", "M1a", "M1b", "M2m", "M2r", "M2a", "M2b")
  for (code in codes) {
    fam <- substr(code, 1, 2); sp <- substr(code, 3, 3)
    D_mode <- switch(sp, m = 1L, r = 2L, a = 3L, b = 3L); joint <- as.integer(sp == "b")
    d <- build_stan_data(sim$sy, sim$alts, sim$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"), joint = joint, x_ref = rep(0, 4),
                         K = if (fam == "M2") 2L else NULL, D_mode = D_mode)
    f <- sampling(sms[[fam]], data = d, chains = 1, iter = 24, warmup = 12, seed = 123, refresh = 0,
                  pars = c("z_re", "z_b", "V_raw", "U_free", "L_Om_joint", "L_Om", "re"), include = FALSE)
    dr <- extract_draws(f, n_use = 6)
    if (D_mode != 3L) { expect_equal(length(dim(dr$beta)), 3L); expect_equal(dim(dr$beta)[3], 1L); expect_equal(dim(dr$a_i)[3], 1L) }   # the 3-d arrays that crashed the old scorer
    sc <- stage2_scores(d, dr)
    expect_equal(nrow(sc), d$E, info = code); expect_true(all(is.finite(sc$log_p_chosen)), info = code)
    expect_true(all(c("hit1_pdraw", "hit5_pdraw", "hit1_pmean", "hit5_pmean", "rank_pmean", "w") %in% names(sc)))
    ## the scorer's log score at a single draw equals the brute-force likelihood at that draw
    one <- lapply(dr, function(a) { nd <- length(dim(a)); if (nd <= 1) a[1] else if (nd == 2) a[1, , drop = FALSE] else a[1, , , drop = FALSE] })
    one$rho <- dr$rho[1]
    D <- if (D_mode == 3L) 2L else 1L
    pars <- list(mu0 = as.numeric(one$mu0), beta = matrix(one$beta, d$P, D), a_i = matrix(one$a_i, d$S, D), b_j = as.numeric(one$b_j),
                 gamma = as.numeric(one$gamma), sigma_x = as.numeric(one$sigma_x), rho = one$rho,
                 U = if (fam == "M2") matrix(one$U, d$S) else NULL, V = if (fam == "M2") matrix(one$V, d$S) else NULL)
    expect_equal(stage2_scores(d, one)$log_p_chosen, stage2_loglik_R(d, pars)$loglik, tolerance = 1e-8, info = code)
    ## the joint state-year log-lik (B1) sums to stage 1 + the weighted stage 2 per draw, and loo runs on it
    llj <- loglik_joint_sy(f, d); expect_equal(dim(llj)[3], d$N1)
    expect_equal(as.numeric(apply(llj, c(1, 2), sum)), as.numeric(apply(loo::extract_log_lik(f, "log_lik_stage1", merge_chains = FALSE), c(1, 2), sum) +
                                                                  apply(loo::extract_log_lik(f, "log_lik_stage2", merge_chains = FALSE), c(1, 2), sum)), tolerance = 1e-8)
    if (code == "M1b") { lj <- suppressWarnings(loo_joint_sy(f, d)); expect_true(is.finite(lj$estimates["elpd_loo", 1]))
      lc <- suppressWarnings(loo_compare_ab(f, f, d, d)); expect_equal(lc$stage, c("stage2_choice", "stage1_conditional", "joint_stateyear")); expect_equal(lc$elpd_diff_b_minus_a, c(0, 0, 0)) }
  }
})

test_that("B1: the joint state-year log-lik uses the event weights", {
  ll1 <- array(rnorm(2 * 3 * 4), c(2, 3, 4)); ll2 <- array(rnorm(2 * 3 * 3), c(2, 3, 3))
  d <- list(N1 = 4L, E = 3L, ev_sy = c(1L, 1L, 3L), ev_w = c(0.5, 0.5, 1))
  out <- loglik_joint_sy(NULL, d, ll1 = ll1, ll2 = ll2)
  expect_equal(out[, , 1], ll1[, , 1] + 0.5 * ll2[, , 1] + 0.5 * ll2[, , 2]); expect_equal(out[, , 3], ll1[, , 3] + ll2[, , 3]); expect_equal(out[, , 2], ll1[, , 2])
  out1 <- loglik_joint_sy(NULL, d, weighted = FALSE, ll1 = ll1, ll2 = ll2); expect_equal(out1[, , 1], ll1[, , 1] + ll2[, , 1] + ll2[, , 2])
  expect_error(loglik_joint_sy(NULL, list(N1 = 5L, E = 3L, ev_sy = c(1L, 1L, 3L)), ll1 = ll1, ll2 = ll2))
})

test_that("B8: perturb_relational_levels() with delta = 0 reproduces the frames; a positive delta moves gaps, changes and the level together", {
  sc <- simulate_current("M1b", S = 12L, Tt = 6L, seed = 3L)
  lv <- sc$levels; flag <- sc$sy |> transmute(state_idx, year, flag = y1)
  p0 <- perturb_relational_levels(sc$sy, sc$alts, lv, flag, delta = 0, h = 1L, gap_scale_R = sc$truth$gap_R)
  expect_equal(p0$alts$dR, sc$alts$dR, tolerance = 1e-10); expect_equal(p0$sy$dR_i, sc$sy$dR_i, tolerance = 1e-10)
  expect_equal(p0$sy$R_lvl, sc$sy$R_lvl, tolerance = 1e-10); expect_equal(p0$alts$dRpr_j_h, sc$alts$dRpr_j_h, tolerance = 1e-10)
  expect_equal(p0$alts$dM, sc$alts$dM); expect_equal(p0$sy$dM_i, sc$sy$dM_i)                      # the material side is untouched
  p1 <- perturb_relational_levels(sc$sy, sc$alts, lv, flag, delta = 0.5, h = 1L, gap_scale_R = sc$truth$gap_R)
  ## a flagged initiator level at t-1 raises its gaps to every target by delta / gap constant exactly (targets unflagged at t-1)
  fl <- lv |> inner_join(flag, by = c("state_idx", "year")) |> filter(flag == 1)
  n <- which(paste(sc$sy$state_idx, sc$sy$year - 1L) %in% paste(fl$state_idx, fl$year))[1]
  rows <- which(sc$alts$sy_row == n); tg_flag <- paste(sc$alts$target_idx[rows], sc$alts$year[rows] - 1L) %in% paste(fl$state_idx, fl$year)
  expect_equal(p1$alts$dR[rows][!tg_flag] - sc$alts$dR[rows][!tg_flag], rep(0.5 / sc$truth$gap_R, sum(!tg_flag)), tolerance = 1e-10)
  expect_equal(p1$alts$dR[rows][tg_flag] - sc$alts$dR[rows][tg_flag], rep(0, sum(tg_flag)), tolerance = 1e-10)   # both ends inflated: no change
  ## the change is rebuilt with BOTH window endpoints, re-centred and re-scaled: the pooled SD moves and the year means stay 0
  expect_false(isTRUE(all.equal(p1$sd_change, sd(sc$sy$dRpr_raw_h - ave(sc$sy$dRpr_raw_h, sc$sy$year)))))
  expect_lt(max(abs(tapply(p1$sy$dR_i, p1$sy$year, mean))), 1e-10); expect_equal(sd(p1$sy$dR_i), 1, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(p1$sy$R_lvl, sc$sy$R_lvl))); expect_lt(max(abs(tapply(p1$sy$R_lvl * sc$truth$gap_R, p1$sy$year, mean))), 1e-10)
  expect_false(isTRUE(all.equal(p1$alts$dRpr_j_h, sc$alts$dRpr_j_h))); expect_lt(abs(mean(p1$alts$dRpr_j_h[p1$alts$in_iv == 1L])), 1e-10)
  expect_equal(p1$n_flagged, nrow(fl))
  expect_error(perturb_relational_levels(sc$sy, sc$alts, lv |> filter(year > 2), flag, 0.5, 1L, sc$truth$gap_R), "missing")
})

test_that("B9: freezing the preprocessing constants on all rows reproduces the frames; on a cutoff it changes them", {
  sc <- simulate_current("M1b", S = 12L, Tt = 8L, seed = 4L)
  spec <- list(positions = sc$levels, gap_raw = c(dM = "dM_logratio", dR = "dR_pr_logratio"), levels = list(dM = "M_lvl", dR = "R_lvl"),
               change_raw = c(dM_i = "dM_raw_h", dR_i = "dRpr_raw_h"), target_change = c(dM_i = "dM_j_h", dR_i = "dRpr_j_h"),
               ctrl_cent = c("ongoing_tm1", "log_prior_mids", "dem_tm1", "log_n_system_tm1"), z_cent = c("log_dy_prior_mids", "dM_j_h", "dRpr_j_h"))
  all_rows <- rep(TRUE, nrow(sc$sy))
  fz <- freeze_preprocessing(sc$sy, sc$alts, all_rows, spec, sc$cols$X)
  for (cc in c("dM_i", "dR_i", "M_lvl", "R_lvl", "ongoing_tm1", "log_n_system_tm1")) expect_equal(fz$sy[[cc]], sc$sy[[cc]], tolerance = 1e-8, info = cc)
  for (cc in c("dM", "dR", "dM_j_h", "dRpr_j_h", "log_dy_prior_mids")) expect_equal(fz$alts[[cc]], sc$alts[[cc]], tolerance = 1e-8, info = cc)
  expect_equal(unname(fz$x_ref), unname(sc$x_ref), tolerance = 1e-10)
  expect_equal(unname(fz$consts$gap_new["dM"]), sc$truth$gap_M, tolerance = 1e-8); expect_equal(fz$consts$sd_new_dR_i, sc$truth$sd_change_R, tolerance = 1e-8)
  ## frozen on the first half of the years: the constants come from the training years only and every row is rescaled by them
  cut <- median(sc$sy$year); tr <- sc$sy$year <= cut
  fc <- freeze_preprocessing(sc$sy, sc$alts, tr, spec, sc$cols$X)
  pos_tr <- sc$levels |> filter(year >= min(sc$sy$year) - 1L, year <= cut - 1L)
  expect_equal(unname(fc$consts$gap_new["dM"]), gap_scale(pos_tr$M, pos_tr$year)); expect_false(isTRUE(all.equal(fc$consts$gap_new["dM"], fc$consts$gap_old["dM"])))
  expect_equal(fc$alts$dM, sc$alts$dM_logratio / fc$consts$gap_new["dM"], tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(sd(fc$sy$dR_i[tr]), 1, tolerance = 1e-10); expect_false(isTRUE(all.equal(sd(fc$sy$dR_i), 1)))
  expect_lt(abs(mean(fc$sy$ongoing_tm1[tr])), 1e-10); expect_lt(abs(mean(fc$alts$dM_j_h[fc$alts$in_iv == 1L & fc$alts$sy_row %in% which(tr)])), 1e-10)
  expect_equal(unname(fc$x_ref), unname(colMeans(as.matrix(fc$sy[tr & fc$sy$y1 == 1L, sc$cols$X]))))
  ## the levels are rescaled by old / new gap constant (their within-year centring is split-invariant)
  expect_equal(fc$sy$M_lvl, sc$sy$M_lvl * fc$consts$gap_old["dM"] / fc$consts$gap_new["dM"], tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("B7: the permutation summary labels its modes and the refit p-value formula is (1 + exceedances) / (1 + n_perm)", {
  src <- readLines(file.path(root, "R", "measurement_tests.R"))
  expect_true(any(grepl("(1 + n_exceed) / (1 + n_perm)", src, fixed = TRUE)))
  expect_true(any(grepl("RELIANCE DIAGNOSTIC", src, fixed = TRUE)))
  expect_true(any(grepl('loo_from_fit(fit_a, "log_lik_stage2")$estimates["elpd_loo", 1]', src, fixed = TRUE)))   # observed statistic by the SAME procedure in refit mode
})
