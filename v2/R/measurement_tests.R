## ============================================================================
## measurement_tests.R -- the PI's planned measurement tests for the a (independent) vs
## b (joint) comparison. Functions first; a small CLI at the bottom for fitted real-data cells.
##
##   loo_from_fit(fit, par)                       LOO of one generated log-lik vector
##   loglik_joint_sy(fit, dat) / loo_joint_sy()   JOINT state-year LOO (review 2026-10-25, B1): pointwise log-lik of
##                                                state-year n = log_lik_stage1[n] + sum over its events of
##                                                ev_w[e] * log_lik_stage2[e], so leaving a state-year out removes its
##                                                initiation AND its target choices. The stage-1-only LOO is the
##                                                CONDITIONAL one (the state-year's events stay in) and is labelled so.
##   loo_compare_ab(fit_a, fit_b, dat_a, dat_b)   LOO, stage 2 (choice), stage 1 (conditional) and joint state-year
##   misclassification_ab(fit_a, fit_b, r0)       share of states whose MAP type differs
##   stage2_scores(dat, draws)                    posterior predictive scores per event (works for every D_mode; B3)
##   temporal_holdout(...)                        train <= cutoff year, score later events; the pooled preprocessing
##                                                constants are FROZEN on the training years (B9)
##   permutation_test(...)                        shuffle dR within each stage-2 set; refit mode scores observed and
##                                                permuted data with the SAME procedure and gives a calibrated p-value;
##                                                rescore mode is a reliance diagnostic (B7)
##   perturb_relational_levels(...)               inflate the relational LEVEL of flagged state-years and rebuild the
##                                                gaps, the initiator / target changes and the level column (B8)
##   measurement_error_check(...)                 refit on the perturbed data, track how classification changes
##
## All refits use seed 123. Every build_stan_data() call passes allow_orphan_initiations = TRUE,
## matching 02_fit_model.R: a state-year can have y1 = 1 with no usable stage-2 event when its
## only target lacks positions at t-1. Functions take the sy / alts / ev frames used by build_stan_data()
## so they run unchanged on simulated data (v2/sim/run_sim.R --task=measurement) and on the real
## frames assembled by assemble_cell().
## ============================================================================
suppressPackageStartupMessages({ library(rstan); library(loo); library(dplyr) })
.mt_dir <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), error = function(e) ".")
if (!exists("build_stan_data")) { source(file.path(.mt_dir, "spec.R")); source(file.path(.mt_dir, "helpers_data.R")); source(file.path(.mt_dir, "helpers_model.R")) }
if (!exists("classify_types")) source(file.path(.mt_dir, "typology.R"))
EXCL_PARS <- c("z_re", "z_b", "V_raw", "U_free", "L_Om_joint", "L_Om", "re")

## ---- 1. LOO ----------------------------------------------------------------------------------------
loo_from_fit <- function(fit, par) {
  ll <- loo::extract_log_lik(fit, par, merge_chains = FALSE)
  loo::loo(ll, r_eff = loo::relative_eff(exp(ll)))
}
## joint state-year pointwise log-lik [iter, chain, N1] (B1). weighted = TRUE uses the stage-2 event weights dat$ev_w (D8), i.e. the
## state-year's contribution to the composite log density; with weights all 1 (or weighted = FALSE) it is the plain sum.
loglik_joint_sy <- function(fit, dat, weighted = TRUE, ll1 = NULL, ll2 = NULL) {
  if (is.null(ll1)) ll1 <- loo::extract_log_lik(fit, "log_lik_stage1", merge_chains = FALSE)       # [iter, chain, N1]; passed directly by the unit test
  stopifnot(length(dim(ll1)) == 3, dim(ll1)[3] == dat$N1)
  out <- ll1
  if (dat$E > 0) {
    if (is.null(ll2)) ll2 <- loo::extract_log_lik(fit, "log_lik_stage2", merge_chains = FALSE)
    stopifnot(length(dim(ll2)) == 3, dim(ll2)[3] == dat$E, length(dat$ev_sy) == dat$E)
    w <- if (weighted && !is.null(dat$ev_w)) as.numeric(dat$ev_w) else rep(1, dat$E)
    stopifnot(length(w) == dat$E, all(is.finite(w)))
    for (e in seq_len(dat$E)) out[, , dat$ev_sy[e]] <- out[, , dat$ev_sy[e]] + w[e] * ll2[, , e]
    ## the sum over state-years of the joint pointwise log-lik equals stage 1 + the weighted stage 2, draw by draw
    tot1 <- apply(ll1, c(1, 2), sum); tot2 <- apply(sweep(ll2, 3, w, "*"), c(1, 2), sum)
    stopifnot(max(abs(apply(out, c(1, 2), sum) - (tot1 + tot2))) < 1e-8)
  }
  out
}
loo_joint_sy <- function(fit, dat, weighted = TRUE) {
  ll <- loglik_joint_sy(fit, dat, weighted)
  loo::loo(ll, r_eff = loo::relative_eff(exp(ll)))
}
## a vs b with three rows: stage2_choice (stage-2 terms only), stage1_conditional (the Bernoulli term alone, the state-year's
## events staying in) and joint_stateyear (stage 1 + the state-year's weighted stage-2 terms). dat_a / dat_b are the fits' Stan
## data lists (needed for the event -> state-year map); when NULL only the two conditional rows are returned.
loo_compare_ab <- function(fit_a, fit_b, dat_a = NULL, dat_b = NULL) {
  l2a <- loo_from_fit(fit_a, "log_lik_stage2"); l2b <- loo_from_fit(fit_b, "log_lik_stage2")
  l1a <- loo_from_fit(fit_a, "log_lik_stage1"); l1b <- loo_from_fit(fit_b, "log_lik_stage1")
  c2 <- loo::loo_compare(l2a, l2b); c1 <- loo::loo_compare(l1a, l1b)
  row <- function(stage, la, lb, cmp) data.frame(stage = stage, elpd_a = la$estimates["elpd_loo", 1], elpd_b = lb$estimates["elpd_loo", 1],
                                                elpd_diff_b_minus_a = lb$estimates["elpd_loo", 1] - la$estimates["elpd_loo", 1], se_diff = cmp[2, "se_diff"],
                                                n_pareto_k_gt_0.7_a = sum(la$diagnostics$pareto_k > 0.7), n_pareto_k_gt_0.7_b = sum(lb$diagnostics$pareto_k > 0.7))
  out <- bind_rows(row("stage2_choice", l2a, l2b, c2), row("stage1_conditional", l1a, l1b, c1))
  if (!is.null(dat_a) && !is.null(dat_b)) {
    stopifnot(identical(dat_a$ev_sy, dat_b$ev_sy), dat_a$N1 == dat_b$N1)
    lja <- loo_joint_sy(fit_a, dat_a); ljb <- loo_joint_sy(fit_b, dat_b)
    out <- bind_rows(out, row("joint_stateyear", lja, ljb, loo::loo_compare(lja, ljb)))
  }
  out
}

## ---- 2. misclassification a vs b ------------------------------------------------------------------------
misclassification_ab <- function(fit_a, fit_b, r0, only_states = NULL) {
  ca <- classify_types(mu_state_draws(fit_a), r0 = r0); cb <- classify_types(mu_state_draws(fit_b), r0 = r0)
  stopifnot(nrow(ca) == nrow(cb))
  keep <- if (is.null(only_states)) rep(TRUE, nrow(ca)) else seq_len(nrow(ca)) %in% only_states
  lev <- if ("p_bully" %in% names(ca)) TYPE_LEVELS else TYPE_LEVELS_1D
  tab <- table(a = factor(ca$map_type[keep], lev), b = factor(cb$map_type[keep], lev))
  list(rate = mean(ca$map_type[keep] != cb$map_type[keep]), table = tab,
       mean_abs_prob_shift = mean(abs(as.matrix(ca[keep, paste0("p_", lev)]) - as.matrix(cb[keep, paste0("p_", lev)]))))
}

## ---- vectorised posterior predictive scoring of stage 2 ------------------------------------------------------
## draws: output of extract_draws(fit); dat: Stan data list (possibly new / permuted data with the same S)
extract_draws <- function(fit, n_use = 100, seed = V2_SEED) {
  pars <- c("mu0", "beta", "a_i", "b_j", "sigma_x", "rho", "gamma")
  has_uv <- "U" %in% fit@model_pars
  if (has_uv) pars <- c(pars, "U", "V")
  dr <- rstan::extract(fit, pars = pars)
  n <- length(dr$rho); set.seed(seed); idx <- sort(sample.int(n, min(n_use, n)))
  ## rstan returns scalars as 1-d arrays (dim = n), vectors as [n, k], matrices as [n, r, c]; a matrix[P, 1] stays 3-d [n, P, 1]
  lapply(dr, function(a) { nd <- length(dim(a))
    if (nd <= 1) as.numeric(a)[idx] else if (nd == 2) a[idx, , drop = FALSE] else a[idx, , , drop = FALSE] })
}
## one draw of a parameter as an nrow x ncol matrix, whatever rank rstan gave the array (B3: matrix[P, 1] and matrix[S, 1] are 3-d)
draw_matrix <- function(a, d, nrow, ncol) {
  nd <- length(dim(a))
  v <- if (nd == 3) a[d, , ] else if (nd == 2) a[d, ] else a[d]
  stopifnot(length(v) == nrow * ncol)
  matrix(v, nrow, ncol)
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
  ## choice-set flags (review P1-1): with full sets the stage-2 softmax runs over the state-year's in_s2 rows only; the
  ## IV-only rows (in_s2 = 0) get utility -Inf, so they carry no probability mass and rank last. With sampled = 1 every
  ## row has in_s2 = 1 (asserted by build_stan_data()).
  if (is.null(dat$in_s2)) stop("stage2_scores: dat$in_s2 is missing (build the data with build_stan_data())")
  stopifnot(length(dat$in_s2) == dat$N2, all(dat$in_s2[dat$ev_row] == 1L))
  excl <- dat$sampled == 0 & dat$in_s2 == 0L
  set_size <- if (dat$sampled == 1) dat$ev_end - dat$ev_start + 1 else as.numeric(tabulate(alt_sy[dat$in_s2 == 1L], nbins = dat$N1))[dat$ev_sy]
  stopifnot(all(set_size >= 2))
  ev_w <- if (is.null(dat$ev_w)) rep(1, dat$E) else as.numeric(dat$ev_w)
  n_d <- length(draws$rho); E <- dat$E
  p_sum <- numeric(E); rank_sum <- numeric(E); hit1_sum <- numeric(E); hit5_sum <- numeric(E); p_row_sum <- numeric(dat$N2)
  zg0 <- if (dat$Pz > 0) dat$Z else NULL
  D <- if (dat$D_mode == 3L) 2L else 1L
  if (is.null(draws$mu0)) stop("stage2_scores: draws$mu0 is missing (mu = mu0 + X beta + a_i)")
  mu0_mat <- matrix(draws$mu0, ncol = D)                     # [n_d, D]; rstan gives [n, 1] for row_vector[1]
  stopifnot(nrow(mu0_mat) == n_d)
  for (d in seq_len(n_d)) {
    beta <- draw_matrix(draws$beta, d, dat$P, D)
    a_d <- draw_matrix(draws$a_i, d, dat$S, D)
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
      U <- draw_matrix(draws$U, d, dat$S, dim(draws$U)[3]); V <- draw_matrix(draws$V, d, dat$S, dim(draws$V)[3])
      util <- util + rowSums(U[dat$sy_state[alt_sy], , drop = FALSE] * V[dat$alt_target, , drop = FALSE])
    }
    util[excl] <- -Inf
    mx <- ave(util, grp, FUN = max)
    lse_tab <- log(as.numeric(rowsum(exp(util - mx), grp))) + as.numeric(tapply(mx, grp, function(v) v[1]))
    names(lse_tab) <- sort(unique(grp))
    logp <- util - lse_tab[as.character(grp)]
    p_sum <- p_sum + exp(logp[dat$ev_row])
    p_row_sum <- p_row_sum + exp(logp)
    ## rank of the chosen target within its own choice set (1 = best) for THIS draw, and the per-draw hit indicators (B6)
    rk <- ave(-util, grp, FUN = function(v) rank(v, ties.method = "min"))
    rank_sum <- rank_sum + rk[dat$ev_row]
    hit1_sum <- hit1_sum + (rk[dat$ev_row] <= 1); hit5_sum <- hit5_sum + (rk[dat$ev_row] <= 5)
  }
  p <- p_sum / n_d
  ## rank of the chosen target by the POSTERIOR-MEAN probability of every row in its set (B6, the second labelled variant)
  p_row_mean <- p_row_sum / n_d
  rk_pmean <- ave(-p_row_mean, grp, FUN = function(v) rank(v, ties.method = "min"))[dat$ev_row]
  data.frame(event = seq_len(E), p_chosen = p, log_p_chosen = log(p), mean_rank = rank_sum / n_d,
             set_size = set_size, w = ev_w,
             hit1_pdraw = hit1_sum / n_d, hit5_pdraw = hit5_sum / n_d,     # posterior probability that the chosen target ranks first / top 5 (per-draw indicators averaged)
             rank_pmean = rk_pmean, hit1_pmean = rk_pmean <= 1, hit5_pmean = rk_pmean <= 5)   # rank of the posterior-mean probability
}

## ---- 3. temporal holdout ---------------------------------------------------------------------------------------------
## Pooled preprocessing constants frozen on the training years (review 2026-10-25, B9). The assembled frames carry quantities
## scaled or centred over ALL years: the gaps (dM, dR = raw log ratio / gap constant), the initiator changes (dM_i, dR_i =
## year-centred raw change / pooled SD), the target changes in Z (same table), the power levels (M_lvl, R_lvl = within-year
## centred level / gap constant), the grand-centred W / X controls, the IV-row-centred continuous Z columns and x_ref.
## freeze_preprocessing() recomputes every constant from the TRAINING rows only and applies it to train and test alike:
##   gap constants   : gap_scale() on the training positions (spec$positions: state_idx / ccode, year, M, R); the gaps are rebuilt
##                     from the raw log ratios spec$gap_raw = c(dM = "<col>", dR = "<col>") in alts; the levels are rescaled by
##                     old / new gap constant
##   change SDs      : SD of the year-centred raw change over the training sy rows (spec$change_raw = c(dM_i = "<col>", dR_i = "<col>"));
##                     dM_i / dR_i are rebuilt and the target changes spec$target_change = c(dM_i = "dM_j_h", dR_i = "dRpr_j_h") are
##                     rescaled by old / new SD (their within-year centring does not depend on the split)
##   centring        : spec$ctrl_cent (grand-centred W / X controls) over the training sy rows; spec$z_cent over the training IV rows
##   x_ref           : mean X over the initiating training state-years
## The year-centring of the changes and levels is within year, so it needs no freezing. With every row in training the output
## equals the input (asserted by the unit test).
freeze_preprocessing <- function(sy, alts, train_sy_rows, spec, X_cols) {
  stopifnot(is.logical(train_sy_rows), length(train_sy_rows) == nrow(sy), any(train_sy_rows), "sy_row" %in% names(alts))
  tr_alt <- alts$sy_row %in% which(train_sy_rows)
  consts <- list()
  ## the constant the frames were built with, recovered as the (constant) ratio of raw to scaled values
  ratio <- function(raw, scaled) { k <- abs(scaled) > 1e-8; r <- raw[k] / scaled[k]; if (sd(r) > 1e-6 * abs(median(r))) stop("freeze_preprocessing: raw / scaled ratio is not one constant"); median(r) }
  ## gaps
  if (!is.null(spec$gap_raw)) {
    stopifnot(all(spec$gap_raw %in% names(alts)), !is.null(spec$positions), all(c("year", "M", "R") %in% names(spec$positions)))
    pos_tr <- spec$positions |> filter(year >= min(sy$year) - 1L, year <= max(sy$year[train_sy_rows]) - 1L)   # the position years the training state-years use (t - 1)
    stopifnot(nrow(pos_tr) > 2)
    g_new <- c(dM = gap_scale(pos_tr$M, pos_tr$year), dR = gap_scale(pos_tr$R, pos_tr$year))
    g_old <- c(dM = ratio(alts[[spec$gap_raw["dM"]]], alts$dM), dR = ratio(alts[[spec$gap_raw["dR"]]], alts$dR))   # the constants the frames were built with
    stopifnot(all(is.finite(g_old)), all(g_old > 0))
    alts$dM <- alts[[spec$gap_raw["dM"]]] / g_new["dM"]; alts$dR <- alts[[spec$gap_raw["dR"]]] / g_new["dR"]
    if (!is.null(spec$levels)) for (k in names(spec$levels)) { cc <- spec$levels[[k]]; if (cc %in% names(sy)) sy[[cc]] <- sy[[cc]] * g_old[k] / g_new[k] }
    consts$gap_old <- g_old; consts$gap_new <- g_new
  }
  ## changes
  if (!is.null(spec$change_raw)) {
    stopifnot(all(spec$change_raw %in% names(sy)), all(names(spec$change_raw) %in% names(sy)))
    for (k in names(spec$change_raw)) {
      raw <- sy[[spec$change_raw[k]]]
      cen <- raw - ave(raw, sy$year)                                                        # within-year centring (split-invariant)
      sd_new <- sd(cen[train_sy_rows]); sd_old <- ratio(cen, sy[[k]])
      stopifnot(is.finite(sd_new), sd_new > 0, is.finite(sd_old), sd_old > 0)
      sy[[k]] <- cen / sd_new
      tc <- spec$target_change[k]
      if (!is.null(tc) && !is.na(tc) && tc %in% names(alts)) alts[[tc]] <- alts[[tc]] * sd_old / sd_new
      consts[[paste0("sd_old_", k)]] <- sd_old; consts[[paste0("sd_new_", k)]] <- sd_new
    }
  }
  ## centring of the W / X controls over the training state-years
  for (cc in intersect(spec$ctrl_cent, names(sy))) sy[[cc]] <- sy[[cc]] - mean(sy[[cc]][train_sy_rows])
  ## centring of the continuous Z columns over the training IV rows
  iv <- if ("in_iv" %in% names(alts)) alts$in_iv == 1L else rep(TRUE, nrow(alts))
  for (cc in intersect(spec$z_cent, names(alts))) alts[[cc]] <- alts[[cc]] - mean(alts[[cc]][tr_alt & iv])
  x_ref <- colMeans(as.matrix(sy[train_sy_rows & sy$y1 == 1L, X_cols, drop = FALSE]))
  stopifnot(all(is.finite(x_ref)))
  list(sy = sy, alts = alts, x_ref = x_ref, consts = consts)
}
## Train on sy$year <= cutoff, score events with year > cutoff. Fits both a and b for the family. freeze = the spec for
## freeze_preprocessing() (NULL = use the frames as given, the pre-B9 behaviour, flagged in the summary).
temporal_holdout <- function(sy, alts, ev, W_cols, X_cols, Z_cols, cutoff, family = c("m1", "m2"), K = 2L,
                             stan_models, chains = 2, iter = 600, warmup = 300, n_draws = 100, x_ref = NULL, freeze = NULL) {
  family <- match.arg(family)
  stopifnot("year" %in% names(sy), any(sy$year <= cutoff), any(sy$year > cutoff))
  S_all <- max(c(sy$state_idx, alts$target_idx))
  frozen <- !is.null(freeze)
  if (frozen) {
    fz <- freeze_preprocessing(sy, alts, sy$year <= cutoff, freeze, X_cols)
    sy <- fz$sy; alts <- fz$alts; x_ref <- fz$x_ref
    logf("temporal holdout: preprocessing constants frozen on years <= %d: %s", cutoff, paste(sprintf("%s=%.4f", names(unlist(fz$consts)), unlist(fz$consts)), collapse = ", "))
  }
  split_part <- function(rows_keep) {
    sy_p <- sy[rows_keep, , drop = FALSE]
    map <- match(seq_len(nrow(sy)), which(rows_keep))      # old sy_row -> new
    alts_p <- alts |> filter(sy_row %in% which(rows_keep)) |> mutate(sy_row = map[sy_row]) |> arrange(sy_row)
    old_rows <- which(alts$sy_row %in% which(rows_keep))
    ev_p <- ev |> filter(sy_row %in% which(rows_keep)) |>
      mutate(alt_row = match(alt_row, old_rows), sy_row = map[sy_row])
    stopifnot(!anyNA(ev_p$alt_row), all(alts_p$sy_row[ev_p$alt_row] == ev_p$sy_row))
    list(sy = sy_p, alts = alts_p, ev = ev_p)
  }
  tr <- split_part(sy$year <= cutoff); te <- split_part(sy$year > cutoff)
  if (is.null(x_ref)) x_ref <- rep(0, length(X_cols))
  out <- list()
  for (joint in c(0L, 1L)) {
    dat_tr <- build_stan_data(tr$sy, tr$alts, tr$ev, W_cols, X_cols, Z_cols, joint = joint, x_ref = x_ref,
                              K = if (family == "m2") K else NULL, S = S_all, allow_orphan_initiations = TRUE)
    fit <- sampling(stan_models[[family]], data = dat_tr, chains = chains, iter = iter, warmup = warmup, seed = V2_SEED,
                    refresh = 0, control = list(adapt_delta = 0.9), pars = EXCL_PARS, include = FALSE)
    dat_te <- build_stan_data(te$sy, te$alts, te$ev, W_cols, X_cols, Z_cols, joint = joint, x_ref = x_ref,
                              K = if (family == "m2") K else NULL, S = S_all, anchor = dat_tr$anchor, allow_orphan_initiations = TRUE)
    sc <- stage2_scores(dat_te, extract_draws(fit, n_draws))
    out[[if (joint == 1) "b" else "a"]] <- list(fit = fit, scores = sc,
      summary = data.frame(spec = if (joint == 1) "b_joint" else "a_independent", cutoff = cutoff, preprocessing_frozen_on_training = frozen,
                           n_train_events = dat_tr$E, n_test_events = dat_te$E,
                           hit_rate_top1_pdraw = mean(sc$hit1_pdraw), hit_rate_top5_pdraw = mean(sc$hit5_pdraw),
                           hit_rate_top1_pmean = mean(sc$hit1_pmean), hit_rate_top5_pmean = mean(sc$hit5_pmean),
                           mean_rank = mean(sc$mean_rank), mean_set_size = mean(sc$set_size),
                           mean_log_score = mean(sc$log_p_chosen), weighted_mean_log_score = sum(sc$w * sc$log_p_chosen) / sum(sc$w),
                           baseline_log_score_uniform = mean(-log(sc$set_size))))
  }
  out$summary <- bind_rows(out$a$summary, out$b$summary)
  out
}

## ---- 4. permutation test --------------------------------------------------------------------------------------------
## Shuffle dR across the alternatives of each STAGE-2 choice set (dM, Z, targets fixed), so dR carries no
## information about which target was chosen. Since review P1-1 a state-year's row block is the union of its
## stage-1 IV set and its stage-2 set (flags in_iv / in_s2). The shuffle runs only over the in_s2 = 1 rows, so
## the null distribution uses exactly the targets the stage-2 likelihood scores; IV-only rows (in_s2 = 0)
## keep their dR. Without an in_s2 column every row is treated as in_s2 = 1 (the pre-review behaviour).
## Modes (review 2026-10-25, B7). "refit": observed AND permuted data are scored by the SAME procedure, the stage-2 LOO elpd
## of a and b refitted on each data set (the observed statistic comes from fit_a / fit_b, which ARE the refits on the observed
## data); p-value = (1 + #{permuted advantage >= observed}) / (1 + n_perm), a calibrated permutation p-value. "rescore": the
## ORIGINAL fits are rescored on the permuted data (no refit); its exceedance share is a RELIANCE DIAGNOSTIC (how much the fitted
## b model's weighted stage-2 log score depends on the actual dR pattern), not a calibrated p-value, and is labelled so.
permute_dR <- function(alts, seed) {
  set.seed(seed)
  s2 <- if ("in_s2" %in% names(alts)) as.integer(alts$in_s2) else rep(1L, nrow(alts))
  stopifnot(!anyNA(s2), all(s2 %in% 0:1))
  out <- alts |> mutate(.s2 = s2, .row = row_number()) |> group_by(sy_row, .s2) |>
    mutate(dR = if (.s2[1] == 1L && n() > 1) sample(dR) else dR) |> ungroup() |> arrange(.row)
  stopifnot(nrow(out) == nrow(alts), identical(out$sy_row, alts$sy_row),
            all(out$dR[s2 == 0L] == alts$dR[s2 == 0L]))
  out |> select(-.s2, -.row)
}
permutation_test <- function(sy, alts, ev, W_cols, X_cols, Z_cols, fit_a, fit_b, n_perm = 20, mode = c("rescore", "refit"),
                             family = "m1", K = 2L, stan_models = NULL, chains = 2, iter = 500, warmup = 250, n_draws = 100, x_ref = NULL) {
  mode <- match.arg(mode)
  if (is.null(x_ref)) x_ref <- rep(0, length(X_cols))
  Kk <- if (family == "m2") K else NULL
  dat_a <- build_stan_data(sy, alts, ev, W_cols, X_cols, Z_cols, joint = 0L, x_ref = x_ref, K = Kk, allow_orphan_initiations = TRUE)
  dat_b <- build_stan_data(sy, alts, ev, W_cols, X_cols, Z_cols, joint = 1L, x_ref = x_ref, K = Kk, anchor = dat_a$anchor, allow_orphan_initiations = TRUE)
  wscore <- function(sc) sum(sc$w * sc$log_p_chosen)                      # the weighted stage-2 log score (composite likelihood weights)
  if (mode == "refit") {
    stopifnot(!is.null(stan_models))
    ## the observed statistic by the SAME procedure as the permuted ones: stage-2 LOO elpd of the fits on the observed data
    obs_a <- loo_from_fit(fit_a, "log_lik_stage2")$estimates["elpd_loo", 1]; obs_b <- loo_from_fit(fit_b, "log_lik_stage2")$estimates["elpd_loo", 1]
    statistic <- "stage2_loo_elpd_of_refits"
  } else {
    dr_a <- extract_draws(fit_a, n_draws); dr_b <- extract_draws(fit_b, n_draws)
    obs_a <- wscore(stage2_scores(dat_a, dr_a)); obs_b <- wscore(stage2_scores(dat_b, dr_b))
    statistic <- "weighted_stage2_log_score_rescored_with_the_original_fits"
  }
  observed_adv <- obs_b - obs_a
  perm_adv <- numeric(n_perm); perm_a <- numeric(n_perm); perm_b <- numeric(n_perm)
  for (p in seq_len(n_perm)) {
    alts_p <- permute_dR(alts, seed = V2_SEED + p)
    da <- build_stan_data(sy, alts_p, ev, W_cols, X_cols, Z_cols, joint = 0L, x_ref = x_ref, K = Kk, anchor = dat_a$anchor, allow_orphan_initiations = TRUE)
    db <- build_stan_data(sy, alts_p, ev, W_cols, X_cols, Z_cols, joint = 1L, x_ref = x_ref, K = Kk, anchor = dat_a$anchor, allow_orphan_initiations = TRUE)
    if (mode == "refit") {
      fa <- sampling(stan_models[[family]], data = da, chains = chains, iter = iter, warmup = warmup, seed = V2_SEED, refresh = 0, pars = EXCL_PARS, include = FALSE)
      fb <- sampling(stan_models[[family]], data = db, chains = chains, iter = iter, warmup = warmup, seed = V2_SEED, refresh = 0, pars = EXCL_PARS, include = FALSE)
      pa <- loo_from_fit(fa, "log_lik_stage2")$estimates["elpd_loo", 1]; pb <- loo_from_fit(fb, "log_lik_stage2")$estimates["elpd_loo", 1]
    } else {
      pa <- wscore(stage2_scores(da, dr_a)); pb <- wscore(stage2_scores(db, dr_b))
    }
    perm_a[p] <- pa; perm_b[p] <- pb; perm_adv[p] <- pb - pa
  }
  n_exceed <- sum(perm_adv >= observed_adv)
  p_value <- if (mode == "refit") (1 + n_exceed) / (1 + n_perm) else NA_real_
  list(observed_a = obs_a, observed_b = obs_b, observed_advantage_b_minus_a = observed_adv,
       permuted_a = perm_a, permuted_b = perm_b, permuted_advantage = perm_adv,
       p_value = p_value, mode = mode,
       summary = data.frame(mode = mode, statistic = statistic,
                            interpretation = if (mode == "refit") "calibrated permutation p-value: (1 + exceedances) / (1 + n_perm)" else "RELIANCE DIAGNOSTIC (original fits rescored on permuted dR), NOT a calibrated p-value",
                            n_perm = n_perm, observed_a = obs_a, observed_b = obs_b,
                            observed_advantage = observed_adv, mean_permuted_advantage = mean(perm_adv),
                            mean_permuted_a = mean(perm_a), mean_permuted_b = mean(perm_b),
                            n_exceedances = n_exceed, exceedance_share = n_exceed / n_perm, p_value = p_value))
}

## ---- 5. systematic measurement error ---------------------------------------------------------------------------------
## perturb_relational_levels(sy, alts, levels, flag, delta, h, gap_scale_R, cols) (review 2026-10-25, B8). The measured relational
## LEVEL R of every flagged (state, position-year) is inflated by delta log points (flag: data.frame(state_idx, year, flag) at the
## position-year; levels: data.frame(state_idx, year, R) of the raw relational levels for every state and position-year the frames
## use). Everything that depends on the levels is then REBUILT the way the build builds it:
##   gap        alts[[cols$dR]]      = (R'_i(t-1) - R'_j(t-1)) / gap_scale_R
##   change     sy[[cols$dR_i]]      = (year-centred (R'_i(t-1) - R'_i(t-1-h))) / pooled SD of the year-centred perturbed change over sy
##   target chg alts[[cols$dR_j]]    = the same transformed change of the target, re-centred over the IV rows like assemble_cell()
##   level      sy[[cols$R_lvl]]     = (R'_i(t-1) - mean over the sy states of that year) / gap_scale_R   (when present)
## so BOTH window endpoints, the centring and the scaling move with the perturbation instead of one delta added to quantities on
## different scales. delta = 0 returns the inputs up to numerical noise when the frames were built by the same rules (asserted by
## the unit test; the real frames are reproduced to 1e-8 because assemble_cell() keeps the raw columns).
perturb_relational_levels <- function(sy, alts, levels, flag, delta, h, gap_scale_R,
                                      cols = list(dR = "dR", dR_i = "dR_i", dR_j = "dRpr_j_h", R_lvl = "R_lvl")) {
  stopifnot(all(c("state_idx", "year", "R") %in% names(levels)), all(c("state_idx", "year", "flag") %in% names(flag)),
            all(c("state_idx", "year", "y1") %in% names(sy)), all(c("sy_row", "target_idx") %in% names(alts)),
            is.finite(delta), h >= 1, is.finite(gap_scale_R), gap_scale_R > 0, all(flag$flag %in% c(0, 1)))
  assert_no_dup_keys(levels, c("state_idx", "year"), "relational levels"); assert_no_dup_keys(flag, c("state_idx", "year"), "level flags")
  lv <- levels |> left_join(flag, by = c("state_idx", "year")) |> mutate(flag = coalesce(flag, 0), R2 = R + delta * flag)
  getR <- function(idx, yr) { v <- lv$R2[match(paste(idx, yr), paste(lv$state_idx, lv$year))]; if (anyNA(v)) stop("perturb_relational_levels: a needed (state, year) level is missing"); v }
  ## initiator levels at both window endpoints, and the change
  Ri_tm1 <- getR(sy$state_idx, sy$year - 1L); Ri_w0 <- getR(sy$state_idx, sy$year - 1L - h)
  raw <- Ri_tm1 - Ri_w0
  cen <- raw - ave(raw, sy$year); sd_new <- sd(cen)
  stopifnot(is.finite(sd_new), sd_new > 0)
  sy2 <- sy; sy2[[cols$dR_i]] <- cen / sd_new
  if (!is.null(cols$R_lvl) && cols$R_lvl %in% names(sy2)) sy2[[cols$R_lvl]] <- (Ri_tm1 - ave(Ri_tm1, sy$year)) / gap_scale_R
  ## gaps and the target change
  Rj_tm1 <- getR(alts$target_idx, alts$year - 1L)
  alts2 <- alts; alts2[[cols$dR]] <- (Ri_tm1[alts$sy_row] - Rj_tm1) / gap_scale_R
  if (!is.null(cols$dR_j) && cols$dR_j %in% names(alts2)) {
    ## the target's transformed change from the same year-centring and SD; a target without a level at t-1-h keeps 0 (dch_j_missing)
    has_w0 <- paste(alts$target_idx, alts$year - 1L - h) %in% paste(lv$state_idx, lv$year)
    tc <- rep(0, nrow(alts))
    raw_j <- lv$R2[match(paste(alts$target_idx[has_w0], alts$year[has_w0] - 1L), paste(lv$state_idx, lv$year))] -
             lv$R2[match(paste(alts$target_idx[has_w0], alts$year[has_w0] - 1L - h), paste(lv$state_idx, lv$year))]
    ymean <- tapply(raw, sy$year, mean)                                             # the year means of the initiator table (same table in the build)
    tc[has_w0] <- (raw_j - as.numeric(ymean[as.character(alts$year[has_w0])])) / sd_new
    if (anyNA(tc)) stop("perturb_relational_levels: a target change needs a year with no state-year in sy")
    iv <- if ("in_iv" %in% names(alts2)) alts2$in_iv == 1L else rep(TRUE, nrow(alts2))
    alts2[[cols$dR_j]] <- tc - mean(tc[iv])
  }
  stopifnot(all(is.finite(sy2[[cols$dR_i]])), all(is.finite(alts2[[cols$dR]])))
  attr(alts2, "sd_change") <- sd_new
  list(sy = sy2, alts = alts2, sd_change = sd_new, n_flagged = sum(lv$flag == 1))
}
## recent: data.frame(state_idx, year, recent = 0/1) at the POSITION-year level: the measured relational level of (state, year) is
## inflated when the state had recent conflict in that year (the callers pass y1 of (state, year)). Mimics a confounding process in
## which conflict raises measured centrality (e.g. conflict-related IGO activity) without changing true position. deltas are in
## log points of R (the level scale); every derived quantity is rebuilt by perturb_relational_levels().
measurement_error_check <- function(sy, alts, ev, W_cols, X_cols, Z_cols, recent, deltas, baseline_cls, levels, h, gap_scale_R,
                                    family = "m1", joint = 1L, K = 2L, stan_models, chains = 2, iter = 600, warmup = 300, r0 = 0.5, x_ref = NULL,
                                    cols = list(dR = "dR", dR_i = "dR_i", dR_j = "dRpr_j_h", R_lvl = "R_lvl")) {
  stopifnot(all(c("state_idx", "year", "recent") %in% names(recent)), "state_idx" %in% names(alts), "year" %in% names(alts))
  if (is.null(x_ref)) x_ref <- rep(0, length(X_cols))
  Kk <- if (family == "m2") K else NULL
  flag <- recent |> transmute(state_idx, year, flag = recent)
  rows <- list()
  for (delta in deltas) {
    pt <- perturb_relational_levels(sy, alts, levels, flag, delta, h, gap_scale_R, cols)
    dat <- build_stan_data(pt$sy, pt$alts, ev, W_cols, X_cols, Z_cols, joint = joint, x_ref = x_ref, K = Kk, allow_orphan_initiations = TRUE)
    fit <- sampling(stan_models[[family]], data = dat, chains = chains, iter = iter, warmup = warmup, seed = V2_SEED, refresh = 0,
                    control = list(adapt_delta = 0.9), pars = EXCL_PARS, include = FALSE)
    cls <- classify_types(mu_state_draws(fit), r0 = r0)
    stopifnot(nrow(cls) == nrow(baseline_cls))
    rows[[length(rows) + 1]] <- data.frame(delta = delta, n_flagged_levels = pt$n_flagged, sd_change_rebuilt = pt$sd_change,
                                           share_states_reclassified = mean(cls$map_type != baseline_cls$map_type),
                                           mean_shift_mu2 = mean(cls$mu2_mean - baseline_cls$mu2_mean),
                                           mean_abs_shift_mu1 = mean(abs(cls$mu1_mean - baseline_cls$mu1_mean)),
                                           n_bully = sum(cls$map_type == "bully"), n_challenger = sum(cls$map_type == "challenger"),
                                           n_peer = sum(cls$map_type == "peer"))
  }
  bind_rows(rows)
}

## ---- fit / cell consistency (review 2026-10-25, C4) ------------------------------------------------------------------------------
## assert_fit_matches_cell(obj, dat): the saved fit must score exactly the assembled cell (E, N1, S, Z / W / X columns, D_mode,
## event identities = (state-year, chosen target) pairs, event weights)
assert_fit_matches_cell <- function(obj, dat, W_cols, X_cols, Z_cols, what = "fit") {
  d <- obj$dat
  if (d$E != dat$E) stop(sprintf("%s %s has E = %d, the assembled cell %d", what, obj$cell, d$E, dat$E))
  if (d$N1 != dat$N1 || d$S != dat$S) stop(sprintf("%s %s: N1 / S %d / %d differ from the cell's %d / %d", what, obj$cell, d$N1, d$S, dat$N1, dat$S))
  if (!identical(obj$Z_cols, Z_cols)) stop(sprintf("%s %s: Z columns (%s) differ from the cell's (%s)", what, obj$cell, paste(obj$Z_cols, collapse = ","), paste(Z_cols, collapse = ",")))
  if (!identical(obj$W_cols, W_cols) || !identical(obj$X_cols, X_cols)) stop(sprintf("%s %s: W / X columns differ from the cell's", what, obj$cell))
  if (!identical(as.integer(d$ev_sy), as.integer(dat$ev_sy)) || !identical(as.integer(d$alt_target[d$ev_row]), as.integer(dat$alt_target[dat$ev_row])))
    stop(sprintf("%s %s: the events (state-year, chosen target) differ from the assembled cell", what, obj$cell))
  if (!is.null(d$ev_w) && !is.null(dat$ev_w) && max(abs(d$ev_w - dat$ev_w)) > 1e-10) stop(sprintf("%s %s: the event weights differ from the assembled cell", what, obj$cell))
  invisible(TRUE)
}

## ---- CLI for fitted real-data cells (LOO and misclassification only; the refit-based tests
## are run through v2/slurm/04_robustness.sbatch, preset perm_error, which calls 05_robustness.R) -------------------------
## Usage: Rscript v2/R/measurement_tests.R --fit_a=<fit_*a.rds> --fit_b=<fit_*b.rds> [--r0=0.5] [--allow_stale=FALSE]
if (sys.nframe() == 0 && length(grep("^--fit_a=", commandArgs(TRUE))) == 1) {
  args <- commandArgs(TRUE)
  allow_stale <- as.logical(sub("^--allow_stale=", "", c(grep("^--allow_stale=", args, value = TRUE), "--allow_stale=FALSE")[1]))
  fa <- readRDS(sub("^--fit_a=", "", grep("^--fit_a=", args, value = TRUE)))
  fb <- readRDS(sub("^--fit_b=", "", grep("^--fit_b=", args, value = TRUE)))
  r0 <- as.numeric(sub("^--r0=", "", c(grep("^--r0=", args, value = TRUE), "--r0=0.5")[1]))
  out_dir <- dirname(sub("^--fit_b=", "", grep("^--fit_b=", args, value = TRUE)))
  ## C1 / C2 / C4: production fits only, current data, and both fits must score the same cell
  for (o in list(fa, fb)) { assert_production_fit(o, allow_stale = allow_stale); assert_fit_current(o, file.path(v2_root(), if (is.null(o$data_dir)) "data" else o$data_dir, "manifest.csv")) }
  assert_fit_matches_cell(fa, fb$dat, fb$W_cols, fb$X_cols, fb$Z_cols, what = "fit_a vs fit_b:")
  stopifnot(fa$dat$D_mode == fb$dat$D_mode)
  tag <- paste0(fa$cell, "_VS_", fb$cell)
  lc <- loo_compare_ab(fa$fit, fb$fit, fa$dat, fb$dat); print(lc)
  write.csv(lc, file.path(out_dir, paste0("loo_ab_", tag, ".csv")), row.names = FALSE)
  init_states <- which(tabulate(fa$dat$sy_state[fa$dat$ev_sy], nbins = fa$dat$S) > 0)
  mc <- misclassification_ab(fa$fit, fb$fit, r0 = r0, only_states = init_states)
  cat(sprintf("misclassification rate a vs b (initiating states, r0 = %.2f): %.3f\n", r0, mc$rate)); print(mc$table)
  write.csv(as.data.frame(mc$table), file.path(out_dir, paste0("misclass_ab_", tag, ".csv")), row.names = FALSE)
}
