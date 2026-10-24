## ============================================================================
## measurement_tests.R -- the PI's planned measurement tests for the a (independent) vs
## b (joint) comparison. Functions first; a small CLI at the bottom for fitted real-data cells.
##
##   loo_compare_ab(fit_a, fit_b)                 LOO, stage 1 and stage 2 kept separate
##   misclassification_ab(fit_a, fit_b, r0)       share of states whose MAP type differs
##   temporal_holdout(...)                        train <= cutoff year, score later events
##                                                (hit rate top-1 / top-5, mean rank, log score)
##   permutation_test(...)                        shuffle dR across alternatives within each
##                                                choice set; rescore (default) or refit; compare
##                                                the joint model's advantage with the observed one
##   measurement_error_check(...)                 inflate R for states with recent conflict by
##                                                delta, refit, track how classification changes
##
## All refits use seed 123. Every build_stan_data() call passes allow_orphan_initiations = TRUE,
## matching 02_fit_model.R: a state-year can have y1 = 1 with no usable stage-2 event when its
## only target lacks positions at t-1 (108 such rows in mid / lag 5 / polrel / cutoff 1990). Functions take the sy / alts / ev frames used by build_stan_data()
## so they run unchanged on simulated data (v2/sim/run_measurement_tests.R) and on the real
## frames assembled by 02_fit_model.R.
## ============================================================================
suppressPackageStartupMessages({ library(rstan); library(loo); library(dplyr) })
.mt_dir <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), error = function(e) ".")
if (!exists("build_stan_data")) source(file.path(.mt_dir, "00_helpers.R"))
if (!exists("classify_types")) source(file.path(.mt_dir, "typology.R"))

## ---- 1. LOO a vs b ----------------------------------------------------------------------------------
loo_from_fit <- function(fit, par) {
  ll <- loo::extract_log_lik(fit, par, merge_chains = FALSE)
  loo::loo(ll, r_eff = loo::relative_eff(exp(ll)))
}
loo_compare_ab <- function(fit_a, fit_b) {
  l2a <- loo_from_fit(fit_a, "log_lik_stage2"); l2b <- loo_from_fit(fit_b, "log_lik_stage2")
  l1a <- loo_from_fit(fit_a, "log_lik_stage1"); l1b <- loo_from_fit(fit_b, "log_lik_stage1")
  c2 <- loo::loo_compare(l2a, l2b); c1 <- loo::loo_compare(l1a, l1b)
  data.frame(stage = c("stage2_choice", "stage1_initiation"),
             elpd_a = c(l2a$estimates["elpd_loo", 1], l1a$estimates["elpd_loo", 1]),
             elpd_b = c(l2b$estimates["elpd_loo", 1], l1b$estimates["elpd_loo", 1]),
             elpd_diff_b_minus_a = c(l2b$estimates["elpd_loo", 1] - l2a$estimates["elpd_loo", 1],
                                     l1b$estimates["elpd_loo", 1] - l1a$estimates["elpd_loo", 1]),
             se_diff = c(c2[2, "se_diff"], c1[2, "se_diff"]),
             n_pareto_k_gt_0.7_a = c(sum(l2a$diagnostics$pareto_k > 0.7), sum(l1a$diagnostics$pareto_k > 0.7)),
             n_pareto_k_gt_0.7_b = c(sum(l2b$diagnostics$pareto_k > 0.7), sum(l1b$diagnostics$pareto_k > 0.7)))
}

## ---- 2. misclassification a vs b ------------------------------------------------------------------------
misclassification_ab <- function(fit_a, fit_b, r0, only_states = NULL) {
  ca <- classify_types(mu_state_draws(fit_a), r0 = r0); cb <- classify_types(mu_state_draws(fit_b), r0 = r0)
  stopifnot(nrow(ca) == nrow(cb))
  keep <- if (is.null(only_states)) rep(TRUE, nrow(ca)) else seq_len(nrow(ca)) %in% only_states
  tab <- table(a = factor(ca$map_type[keep], TYPE_LEVELS), b = factor(cb$map_type[keep], TYPE_LEVELS))
  list(rate = mean(ca$map_type[keep] != cb$map_type[keep]), table = tab,
       mean_abs_prob_shift = mean(abs(as.matrix(ca[keep, paste0("p_", TYPE_LEVELS)]) - as.matrix(cb[keep, paste0("p_", TYPE_LEVELS)]))))
}

## ---- vectorised posterior predictive scoring of stage 2 ------------------------------------------------------
## draws: output of extract_draws(fit); dat: Stan data list (possibly new / permuted data with the same S)
extract_draws <- function(fit, n_use = 100, seed = V2_SEED) {
  pars <- c("mu0", "beta", "a_i", "b_j", "sigma_x", "rho", "gamma")
  has_uv <- "U" %in% fit@model_pars
  if (has_uv) pars <- c(pars, "U", "V")
  dr <- rstan::extract(fit, pars = pars)
  n <- length(dr$rho); set.seed(seed); idx <- sort(sample.int(n, min(n_use, n)))
  ## rstan returns scalars as 1-d arrays (dim = n), vectors as [n, k], matrices as [n, r, c]
  lapply(dr, function(a) { nd <- length(dim(a))
    if (nd <= 1) as.numeric(a)[idx] else if (nd == 2) a[idx, , drop = FALSE] else a[idx, , , drop = FALSE] })
}
stage2_scores <- function(dat, draws) {
  ## owner state-year of every row: IV rows first, then (if sampled) the per-event blocks
  alt_sy <- rep(seq_len(dat$N1), times = dat$cs_end - dat$cs_start + 1)
  ev_block <- rep(0L, length(alt_sy))                       # 0 = IV row, e = row of event e's own set
  if (dat$sampled == 1) {
    alt_sy <- c(alt_sy, rep(dat$ev_sy, times = dat$ev_end - dat$ev_start + 1))
    ev_block <- c(ev_block, rep(seq_len(dat$E), times = dat$ev_end - dat$ev_start + 1))
  }
  stopifnot(length(alt_sy) == dat$N2)
  ## group id for the softmax: state-year for IV rows, N1 + e for event blocks (full sets: the state-year)
  grp <- ifelse(ev_block > 0L, dat$N1 + ev_block, alt_sy)
  ev_grp <- if (dat$sampled == 1) dat$N1 + seq_len(dat$E) else dat$ev_sy
  ## choice-set flags (review P1-1): with full sets the stage-2 softmax runs over the state-year's in_s2 rows only; the
  ## IV-only rows (in_s2 = 0) get utility -Inf, so they carry no probability mass and rank last. With sampled = 1 every
  ## row has in_s2 = 1 (asserted by build_stan_data()).
  if (is.null(dat$in_s2)) stop("stage2_scores: dat$in_s2 is missing (build the data with build_stan_data())")
  stopifnot(length(dat$in_s2) == dat$N2, all(dat$in_s2[dat$ev_row] == 1L))
  excl <- dat$sampled == 0 & dat$in_s2 == 0L
  set_size <- if (dat$sampled == 1) dat$ev_end - dat$ev_start + 1 else as.numeric(tabulate(alt_sy[dat$in_s2 == 1L], nbins = dat$N1))[dat$ev_sy]
  stopifnot(all(set_size >= 2))
  n_d <- length(draws$rho); E <- dat$E
  p_sum <- numeric(E); rank_sum <- numeric(E)
  zg0 <- if (dat$Pz > 0) dat$Z else NULL
  D <- if (dat$D_mode == 3L) 2L else 1L
  if (is.null(draws$mu0)) stop("stage2_scores: draws$mu0 is missing (mu = mu0 + X beta + a_i)")
  mu0_mat <- matrix(draws$mu0, ncol = D)                     # [n_d, D]; rstan gives [n, 1] for row_vector[1]
  stopifnot(nrow(mu0_mat) == n_d)
  for (d in seq_len(n_d)) {
    beta <- matrix(if (D == 2L) draws$beta[d, , ] else draws$beta[d, ], dat$P, D)
    a_d <- matrix(if (D == 2L) draws$a_i[d, , ] else draws$a_i[d, ], dat$S, D)
    Mu <- rep(1, dat$N1) %o% mu0_mat[d, ] + dat$X %*% beta + a_d[dat$sy_state, , drop = FALSE]
    sx <- matrix(draws$sigma_x, ncol = D)[d, ]; rho <- draws$rho[d]
    if (D == 2L) {
      d1 <- (dat$x[, 1] - Mu[alt_sy, 1]) / sx[1]; d2 <- (dat$x[, 2] - Mu[alt_sy, 2]) / sx[2]
      util <- -0.5 * (d1^2 - 2 * rho * d1 * d2 + d2^2) / (1 - rho^2)
    } else {
      d1 <- (dat$x[, dat$D_mode] - Mu[alt_sy, 1]) / sx[1]
      util <- -0.5 * d1^2
    }
    util <- util + draws$b_j[d, dat$alt_target]
    if (!is.null(zg0)) util <- util + as.numeric(zg0 %*% draws$gamma[d, ])
    if (!is.null(draws$U)) {
      U <- matrix(draws$U[d, , ], dat$S); V <- matrix(draws$V[d, , ], dat$S)
      util <- util + rowSums(U[dat$sy_state[alt_sy], , drop = FALSE] * V[dat$alt_target, , drop = FALSE])
    }
    util[excl] <- -Inf
    mx <- ave(util, grp, FUN = max)
    lse_tab <- log(as.numeric(rowsum(exp(util - mx), grp))) + as.numeric(tapply(mx, grp, function(v) v[1]))
    names(lse_tab) <- sort(unique(grp))
    logp <- util - lse_tab[as.character(grp)]
    p_sum <- p_sum + exp(logp[dat$ev_row])
    ## rank of the chosen target within its own choice set (1 = best) for this draw
    rk <- ave(-util, grp, FUN = rank)
    rank_sum <- rank_sum + rk[dat$ev_row]
  }
  p <- p_sum / n_d
  data.frame(event = seq_len(E), p_chosen = p, log_p_chosen = log(p), mean_rank = rank_sum / n_d,
             set_size = set_size,
             hit1 = (rank_sum / n_d) <= 1.5, hit5 = (rank_sum / n_d) <= 5.5)
}

## ---- 3. temporal holdout ---------------------------------------------------------------------------------------------
## Train on sy$year <= cutoff, score events with year > cutoff. Fits both a and b for the family.
temporal_holdout <- function(sy, alts, ev, W_cols, X_cols, Z_cols, cutoff, family = c("m1", "m2"), K = 2L,
                             stan_models, chains = 2, iter = 600, warmup = 300, n_draws = 100, x_ref = NULL) {
  family <- match.arg(family)
  stopifnot("year" %in% names(sy), any(sy$year <= cutoff), any(sy$year > cutoff))
  S_all <- max(c(sy$state_idx, alts$target_idx))
  split_part <- function(rows_keep) {
    sy_p <- sy[rows_keep, , drop = FALSE]
    map <- match(seq_len(nrow(sy)), which(rows_keep))      # old sy_row -> new
    alts_p <- alts |> filter(sy_row %in% which(rows_keep)) |> mutate(sy_row = map[sy_row]) |> arrange(sy_row)
    alts_p$.new_row <- seq_len(nrow(alts_p))
    old_rows <- which(alts$sy_row %in% which(rows_keep))
    ev_p <- ev |> filter(sy_row %in% which(rows_keep)) |>
      mutate(alt_row = match(alt_row, old_rows), sy_row = map[sy_row])
    stopifnot(!anyNA(ev_p$alt_row), all(alts_p$sy_row[ev_p$alt_row] == ev_p$sy_row))
    list(sy = sy_p, alts = alts_p |> select(-.new_row), ev = ev_p)
  }
  tr <- split_part(sy$year <= cutoff); te <- split_part(sy$year > cutoff)
  if (is.null(x_ref)) x_ref <- rep(0, length(X_cols))
  out <- list()
  for (joint in c(0L, 1L)) {
    dat_tr <- build_stan_data(tr$sy, tr$alts, tr$ev, W_cols, X_cols, Z_cols, joint = joint, x_ref = x_ref,
                              K = if (family == "m2") K else NULL, S = S_all, allow_orphan_initiations = TRUE)
    fit <- sampling(stan_models[[family]], data = dat_tr, chains = chains, iter = iter, warmup = warmup, seed = V2_SEED,
                    refresh = 0, control = list(adapt_delta = 0.9), pars = c("z_re", "z_b", "V_raw", "U_free", "L_Om_joint", "L_Om", "re"), include = FALSE)
    dat_te <- build_stan_data(te$sy, te$alts, te$ev, W_cols, X_cols, Z_cols, joint = joint, x_ref = x_ref,
                              K = if (family == "m2") K else NULL, S = S_all, anchor = dat_tr$anchor, allow_orphan_initiations = TRUE)
    sc <- stage2_scores(dat_te, extract_draws(fit, n_draws))
    out[[if (joint == 1) "b" else "a"]] <- list(fit = fit, scores = sc,
      summary = data.frame(spec = if (joint == 1) "b_joint" else "a_independent", cutoff = cutoff,
                           n_train_events = dat_tr$E, n_test_events = dat_te$E,
                           hit_rate_top1 = mean(sc$hit1), hit_rate_top5 = mean(sc$hit5),
                           mean_rank = mean(sc$mean_rank), mean_set_size = mean(sc$set_size),
                           mean_log_score = mean(sc$log_p_chosen),
                           baseline_log_score_uniform = mean(-log(sc$set_size))))
  }
  out$summary <- bind_rows(out$a$summary, out$b$summary)
  out
}

## ---- 4. permutation test --------------------------------------------------------------------------------------------
## Shuffle dR across the alternatives of each choice set (dM, Z, targets fixed), so dR carries no
## information about which target was chosen. mode = "rescore": rescore the ORIGINAL fits on the
## permuted data (fast); mode = "refit": refit a and b on each permuted data set (slow).
permute_dR <- function(alts, seed) {
  set.seed(seed)
  alts |> group_by(sy_row) |> mutate(dR = if (n() > 1) sample(dR) else dR) |> ungroup()
}
permutation_test <- function(sy, alts, ev, W_cols, X_cols, Z_cols, fit_a, fit_b, n_perm = 20, mode = c("rescore", "refit"),
                             family = "m1", K = 2L, stan_models = NULL, chains = 2, iter = 500, warmup = 250, n_draws = 100, x_ref = NULL) {
  mode <- match.arg(mode)
  if (is.null(x_ref)) x_ref <- rep(0, length(X_cols))
  Kk <- if (family == "m2") K else NULL
  dat_a <- build_stan_data(sy, alts, ev, W_cols, X_cols, Z_cols, joint = 0L, x_ref = x_ref, K = Kk, allow_orphan_initiations = TRUE)
  dat_b <- build_stan_data(sy, alts, ev, W_cols, X_cols, Z_cols, joint = 1L, x_ref = x_ref, K = Kk, anchor = dat_a$anchor, allow_orphan_initiations = TRUE)
  dr_a <- extract_draws(fit_a, n_draws); dr_b <- extract_draws(fit_b, n_draws)
  obs_a <- sum(stage2_scores(dat_a, dr_a)$log_p_chosen); obs_b <- sum(stage2_scores(dat_b, dr_b)$log_p_chosen)
  observed_adv <- obs_b - obs_a
  perm_adv <- numeric(n_perm); perm_a <- numeric(n_perm); perm_b <- numeric(n_perm)
  for (p in seq_len(n_perm)) {
    alts_p <- permute_dR(alts, seed = V2_SEED + p)
    da <- build_stan_data(sy, alts_p, ev, W_cols, X_cols, Z_cols, joint = 0L, x_ref = x_ref, K = Kk, anchor = dat_a$anchor, allow_orphan_initiations = TRUE)
    db <- build_stan_data(sy, alts_p, ev, W_cols, X_cols, Z_cols, joint = 1L, x_ref = x_ref, K = Kk, anchor = dat_a$anchor, allow_orphan_initiations = TRUE)
    if (mode == "refit") {
      fa <- sampling(stan_models[[family]], data = da, chains = chains, iter = iter, warmup = warmup, seed = V2_SEED, refresh = 0,
                     pars = c("z_re", "z_b", "V_raw", "U_free", "L_Om_joint", "L_Om", "re"), include = FALSE)
      fb <- sampling(stan_models[[family]], data = db, chains = chains, iter = iter, warmup = warmup, seed = V2_SEED, refresh = 0,
                     pars = c("z_re", "z_b", "V_raw", "U_free", "L_Om_joint", "L_Om", "re"), include = FALSE)
      pa <- loo_from_fit(fa, "log_lik_stage2")$estimates["elpd_loo", 1]; pb <- loo_from_fit(fb, "log_lik_stage2")$estimates["elpd_loo", 1]
    } else {
      pa <- sum(stage2_scores(da, dr_a)$log_p_chosen); pb <- sum(stage2_scores(db, dr_b)$log_p_chosen)
    }
    perm_a[p] <- pa; perm_b[p] <- pb; perm_adv[p] <- pb - pa
  }
  list(observed_a = obs_a, observed_b = obs_b, observed_advantage_b_minus_a = observed_adv,
       permuted_a = perm_a, permuted_b = perm_b, permuted_advantage = perm_adv,
       p_value = mean(perm_adv >= observed_adv),
       summary = data.frame(mode = mode, n_perm = n_perm, observed_log_score_a = obs_a, observed_log_score_b = obs_b,
                            observed_advantage = observed_adv, mean_permuted_advantage = mean(perm_adv),
                            mean_permuted_log_score_a = mean(perm_a), mean_permuted_log_score_b = mean(perm_b),
                            p_value_adv_ge_observed = mean(perm_adv >= observed_adv)))
}

## ---- 5. systematic measurement error ---------------------------------------------------------------------------------
## recent: data.frame(state_idx, year, recent = 0/1) at the state-year level (positions are at t-1, so
## callers pass the flag for the state-year whose R enters the alternatives). R is inflated by delta
## (SD units) for flagged states, dR = R_i - R_j is rebuilt, the model is refit, and MAP types compared
## with the baseline classification. Mimics a confounding process in which conflict raises measured
## centrality (e.g. conflict-related IGO activity) without changing true position.
measurement_error_check <- function(sy, alts, ev, W_cols, X_cols, Z_cols, recent, deltas, baseline_cls, family = "m1", joint = 1L,
                                    K = 2L, stan_models, chains = 2, iter = 600, warmup = 300, r0 = 0.5, x_ref = NULL,
                                    dR_in_X = "X2") {
  stopifnot(all(c("state_idx", "year", "recent") %in% names(recent)), "state_idx" %in% names(alts), "year" %in% names(alts))
  if (is.null(x_ref)) x_ref <- rep(0, length(X_cols))
  Kk <- if (family == "m2") K else NULL
  rows <- list()
  for (delta in deltas) {
    a2 <- alts |>
      left_join(recent |> transmute(state_idx, year, rec_i = recent), by = c("state_idx", "year")) |>
      left_join(recent |> transmute(target_idx = state_idx, year, rec_j = recent), by = c("target_idx", "year")) |>
      mutate(rec_i = coalesce(rec_i, 0), rec_j = coalesce(rec_j, 0), dR = dR + delta * (rec_i - rec_j)) |> select(-rec_i, -rec_j)
    sy2 <- sy
    if (!is.null(dR_in_X) && dR_in_X %in% names(sy2)) {
      ## the initiator's own relational change also shifts if its measured R is inflated this year but not h years ago
      sy2 <- sy2 |> left_join(recent |> transmute(state_idx, year, rec_i = recent), by = c("state_idx", "year")) |>
        mutate(rec_i = coalesce(rec_i, 0)); sy2[[dR_in_X]] <- sy2[[dR_in_X]] + delta * sy2$rec_i; sy2$rec_i <- NULL
    }
    dat <- build_stan_data(sy2, a2, ev, W_cols, X_cols, Z_cols, joint = joint, x_ref = x_ref, K = Kk, allow_orphan_initiations = TRUE)
    fit <- sampling(stan_models[[family]], data = dat, chains = chains, iter = iter, warmup = warmup, seed = V2_SEED, refresh = 0,
                    control = list(adapt_delta = 0.9), pars = c("z_re", "z_b", "V_raw", "U_free", "L_Om_joint", "L_Om", "re"), include = FALSE)
    cls <- classify_types(mu_state_draws(fit), r0 = r0)
    stopifnot(nrow(cls) == nrow(baseline_cls))
    rows[[length(rows) + 1]] <- data.frame(delta = delta, share_states_reclassified = mean(cls$map_type != baseline_cls$map_type),
                                           mean_shift_mu2 = mean(cls$mu2_mean - baseline_cls$mu2_mean),
                                           mean_abs_shift_mu1 = mean(abs(cls$mu1_mean - baseline_cls$mu1_mean)),
                                           n_bully = sum(cls$map_type == "bully"), n_challenger = sum(cls$map_type == "challenger"),
                                           n_peer = sum(cls$map_type == "peer"))
  }
  bind_rows(rows)
}

## ---- CLI for fitted real-data cells (LOO and misclassification only; the refit-based tests
## are run through v2/slurm/measurement.sbatch which calls run_measurement_real.R) -----------------------------
if (sys.nframe() == 0 && length(grep("^--fit_a=", commandArgs(TRUE))) == 1) {
  args <- commandArgs(TRUE)
  fa <- readRDS(sub("^--fit_a=", "", grep("^--fit_a=", args, value = TRUE)))
  fb <- readRDS(sub("^--fit_b=", "", grep("^--fit_b=", args, value = TRUE)))
  r0 <- as.numeric(sub("^--r0=", "", c(grep("^--r0=", args, value = TRUE), "--r0=0.5")[1]))
  out_dir <- dirname(sub("^--fit_b=", "", grep("^--fit_b=", args, value = TRUE)))
  tag <- paste0(fa$cell, "_VS_", fb$cell)
  lc <- loo_compare_ab(fa$fit, fb$fit); print(lc)
  write.csv(lc, file.path(out_dir, paste0("loo_ab_", tag, ".csv")), row.names = FALSE)
  init_states <- which(tabulate(fa$dat$sy_state[fa$dat$ev_sy], nbins = fa$dat$S) > 0)
  mc <- misclassification_ab(fa$fit, fb$fit, r0 = r0, only_states = init_states)
  cat(sprintf("misclassification rate a vs b (initiating states, r0 = %.2f): %.3f\n", r0, mc$rate)); print(mc$table)
  write.csv(as.data.frame(mc$table), file.path(out_dir, paste0("misclass_ab_", tag, ".csv")), row.names = FALSE)
}
