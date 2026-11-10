## ============================================================================
## opp_weights.R -- functions only. Run 1 (built 2026-11-10; README section 14): the MID-propensity OPPORTUNITY WEIGHT q_ij and its use as a
## fixed stage-2 / stage-1 offset lambda * ln q_ij in the M1c model. A SENSITIVITY specification, not a consideration-set model: q_ij is called an
## opportunity weight everywhere (never a "consideration probability"; design synthesis of 2026-11-08, S1 / S2).
##
## Sourced by 08_opportunity_weights.R (the auxiliary model), 02_fit_model.R (the offset and the per-cell diagnostics), 09_opp_combine.R (the mixture
## over the auxiliary draws) and tests/testthat/test_opp_weights.R. Needs spec.R, helpers_data.R (date_keys, date_cmp, multilateral_jan_flags,
## mid3_history, purrr_free_seq, assert_no_dup_keys) and dplyr / tidyr.
##
## PI DEFAULTS TAKEN HERE (flagged for sign-off in the README): estimand = target selection (Choice A), so NO target-strength variable enters the
## auxiliary predictors A; outcome = Dyadic MID onsets (hostility >= 3 as a second version); the history predictor is an INDICATOR of an onset in
## [t - 3, t - 1]; A4 uses the disputes that are multilateral at 1 January of t (>= 3 participants present, at least one on each side); expanding
## window with the burn-in through 1865 fitted in sample; decade effects with the rule of oppw_decade_levels(); five full-vector draws.
## ============================================================================
OPPW_LAG              <- 3L       # the lag of the cells the weights are built for (the directed comparison is lag 3 only)
OPPW_BURN_IN          <- 1865L    # years <= OPPW_BURN_IN are fitted IN SAMPLE (one fit on all of them); every later year t is predicted from years < t
OPPW_N_DRAWS          <- 5L       # coherent draws of the full auxiliary coefficient vector (lambda = 1 only)
OPPW_HIST_WINDOW      <- 3L       # MID history: an onset in [t - 3, t - 1] (never year t)
OPPW_DECADE_MIN_YEARS <- 3L       # a decade gets its own effect once the training window holds this many of its years AND an onset in it
OPPW_LOW_Q            <- 0.05     # "low q": below this quantile of q over ALL alternative rows of the cell (pooled)
OPPW_OUTCOMES         <- c(mid = "y_mid", host3 = "y_host3")
## the auxiliary predictors A (Choice A: opportunity only). All dated t - 1 or 1 January of t.
OPPW_A_COLS      <- c("contig_land", "contig_sea", "log_capdist", "log_capdist_x_init_M", "dymid_hist3", "a4_joinable")
OPPW_A_BASELINE  <- c("contig_land", "contig_sea", "log_capdist")          # the contiguity + distance baseline of the diagnostics
## columns of the alternatives / panel files that measure the TARGET's strength (or a gap built from it): never read by oppw_frame(), asserted absent from A
OPPW_TARGET_STRENGTH <- c("M_j", "Rpr_j", "Reig_j", "major2", "major_count", "polrel", "dM", "dR", "dR_pr", "dR_eig", "dM_logratio", "dR_pr_logratio", "dM_j_h", "dRpr_j_h",
                          "tgt_major", "tgt_major_x_init_major", "one_major", "both_major")
## the ONLY columns oppw_frame() takes from the alternatives file and the panel
OPPW_ALTS_COLS  <- c("ccode1", "ccode2", "year", "in_iv", "in_s2", "contig_land", "contig_any", "log_capdist")
OPPW_PANEL_COLS <- c("ccode", "year", "M_lvl")
stopifnot(!any(OPPW_TARGET_STRENGTH %in% c(OPPW_A_COLS, OPPW_ALTS_COLS, OPPW_PANEL_COLS)))

## oppw_key(ccode1, ccode2, year): one exact numeric key per directed dyad-year (COW codes < 1000)
oppw_key <- function(ccode1, ccode2, year) {
  stopifnot(!anyNA(ccode1), !anyNA(ccode2), !anyNA(year), all(ccode1 > 0 & ccode1 < 1000), all(ccode2 > 0 & ccode2 < 1000), all(year > 0 & year < 10000))
  (as.numeric(year) * 1000 + as.numeric(ccode1)) * 1000 + as.numeric(ccode2)
}
## oppw_region(ccode): the COW region of a state (the standard code ranges), for the calibration tables only
oppw_region <- function(ccode) {
  stopifnot(!anyNA(ccode), all(ccode >= 2 & ccode <= 990))
  as.character(cut(ccode, breaks = c(1, 199, 399, 599, 699, 990), labels = c("americas", "europe", "africa", "mideast", "asia_oceania")))
}

## ---- A4: the joinable-dispute indicator at 1 January of t --------------------------------------------------------------------------------------------
## oppw_jan_presence(psd, dd, sy): the participants PRESENT at 1 January of t in the disputes that are multilateral at that date, and the states "already in"
## each such dispute. The disputes and the presence rule are those of multilateral_jan_flags() (helpers_data.R; the M1c stage-1 control multi_mid_jan): a
## spell is present when it entered in a year before t (or is day-dated on or before 1 January of t) and has not exited; a dispute is active with at least
## one present participant on each side and multilateral with two or more on a side, which is the same as ACTIVE WITH >= 3 PRESENT PARTICIPANTS. A state is
## "already in" a dispute when ANY of its spells entered before that date (a former participant is a party). Nothing dated in t enters. The presence is
## recomputed here (the helper does not return it) and checked against the helper's dispute table: same dispute-years, same number present.
## psd = participant_spell_dates(dd, dispute_date_bounds(cow_mid_disps, dd), part = gml_part), dd = COW's directed dispute-years (the source of
## multi_mid_jan under BOTH event foundations); sy = (ccode, year) state-years for the helper's flags (returned for the consistency check).
oppw_jan_presence <- function(psd, dd, sy) {
  mj <- multilateral_jan_flags(sy, psd, dd)
  multi <- attr(mj, "disputes")
  stopifnot(all(c("dispnum", "year", "n_side_a", "n_side_b") %in% names(multi)), all(multi$n_side_a >= 1L), all(multi$n_side_b >= 1L), all(multi$n_side_a + multi$n_side_b >= 3L))
  sp <- psd |> select(dispnum, ccode, spell_first, spell_last, entry_res, entry_y, entry_mkey, entry_dkey, exit_res, exit_y, exit_mkey, exit_dkey)
  jan_dkey <- function(t) date_keys(rep(1L, length(t)), t, rep(1L, length(t)), rep(1L, length(t)))$dkey
  entered_by_jan <- function(entry_res, entry_y, entry_dkey, t) entry_y < t | (entry_res == 1L & entry_dkey <= jan_dkey(t))
  ex <- sp |> mutate(year = purrr_free_seq(spell_first, spell_last)) |> tidyr::unnest(year) |> mutate(year = as.integer(year))
  ref <- date_keys(rep(1L, nrow(ex)), ex$year, rep(1L, nrow(ex)), rep(1L, nrow(ex)))
  cx <- date_cmp(list(res = ex$exit_res, y = ex$exit_y, mkey = ex$exit_mkey, dkey = ex$exit_dkey), ref)
  ex$present <- entered_by_jan(ex$entry_res, ex$entry_y, ex$entry_dkey, ex$year) & !(cx$cmp < 0L)
  present <- ex |> filter(present) |> distinct(dispnum, ccode, year) |> semi_join(multi, by = c("dispnum", "year"))
  ## the recomputed presence reproduces the helper's counts in every multilateral dispute-year (a side switcher is one participant in both)
  chk <- present |> count(dispnum, year, name = "n_present") |> full_join(multi, by = c("dispnum", "year"))
  if (anyNA(chk$n_present) || anyNA(chk$n_side_a) || any(chk$n_present != chk$n_side_a + chk$n_side_b))
    stop("oppw_jan_presence: the recomputed January presence differs from multilateral_jan_flags() in ", sum(is.na(chk$n_present) | is.na(chk$n_side_a) | chk$n_present != chk$n_side_a + chk$n_side_b), " dispute-years")
  stopifnot(all(chk$n_present >= 3L))
  already_in <- multi |> select(dispnum, year) |> inner_join(sp |> select(dispnum, ccode, entry_res, entry_y, entry_dkey), by = "dispnum", relationship = "many-to-many") |>
    filter(entered_by_jan(entry_res, entry_y, entry_dkey, year)) |> distinct(dispnum, ccode, year)
  stopifnot(nrow(anti_join(present, already_in, by = c("dispnum", "ccode", "year"))) == 0L)        # every present participant is already in
  list(present = present, already_in = already_in, multi = multi, flags = mj |> select(ccode, year, multi_mid_jan))
}
## oppw_a4(keys, jp): A4_ij = 1 when j is PRESENT (either side) at 1 January of t in a dispute with >= 3 present participants to which i is NOT a party
## (i has no spell in it entered before that date). keys = data.frame(ccode1 = i, ccode2 = j, year = t); jp = oppw_jan_presence(). Returns 0/1 per row.
oppw_a4 <- function(keys, jp) {
  stopifnot(all(c("ccode1", "ccode2", "year") %in% names(keys)), all(keys$ccode1 != keys$ccode2))
  n_j <- jp$present |> count(ccode2 = ccode, year, name = "n_j")
  both <- jp$present |> transmute(dispnum, ccode2 = ccode, year) |>
    inner_join(jp$already_in |> transmute(dispnum, ccode1 = ccode, year), by = c("dispnum", "year"), relationship = "many-to-many") |> count(ccode1, ccode2, year, name = "n_both")
  out <- keys |> select(ccode1, ccode2, year) |> left_join(n_j, by = c("ccode2", "year")) |> left_join(both, by = c("ccode1", "ccode2", "year")) |>
    mutate(n_j = coalesce(n_j, 0L), n_both = coalesce(n_both, 0L))
  stopifnot(nrow(out) == nrow(keys), all(out$n_both <= out$n_j))
  a4 <- as.integer(out$n_j - out$n_both > 0L)
  ## consistency with the stage-1 control: a state with a joinable target has a multilateral dispute it had not entered (multi_mid_jan = 1)
  fl <- keys |> select(ccode = ccode1, year) |> left_join(jp$flags, by = c("ccode", "year"))
  if (any(a4 == 1L & (is.na(fl$multi_mid_jan) | fl$multi_mid_jan != 1L))) stop("oppw_a4: a joinable-dispute row whose initiator has multi_mid_jan = 0")
  a4
}

## ---- the auxiliary frame --------------------------------------------------------------------------------------------------------------------------
## oppw_frame(alts, panel, events_dymid, jp): one row per directed dyad-year (i, j, t) of the cell foundation's lag-3 alternatives file (every row, IV
## or stage 2), sorted by year, with the auxiliary outcomes and predictors.
##   alts          the foundation's alts_lag3_all.rds; ONLY the columns OPPW_ALTS_COLS are read (no target-strength column can reach A)
##   panel         its stage1_panel_lag3.rds; ONLY ccode, year, M_lvl (the initiator's log CINC level at t - 1, within-year centred, as in W)
##   events_dymid  events.rds of the DYMID build (v2/data_dymid) whatever the cell's foundation: one row per interacting dyad and onset year, every
##                 hostility level, escalation rows included (a new dyadic onset of a pair already fighting elsewhere), every sample flag
##   jp            oppw_jan_presence()
## Outcomes (undirected): y_mid = a DYMID onset between i and j in year t; y_host3 = one whose dyadic hostility max(hostlev1, hostlev2) is >= 3.
## Predictors: contig_land (conttype 1), contig_sea (conttype 2..5 = contig_any - contig_land; the build keeps no finer conttype), log_capdist (log1p km
## at t - 1), M_lvl_i, dymid_hist3 (INDICATOR of a DYMID onset between i and j in [t - 3, t - 1], mid3_history(); all hostility levels in both versions),
## a4_joinable. The product log_capdist x M_lvl_i is formed per training window (oppw_design()). Rows without a capital distance are dropped and counted
## (the directed assembly drops them too). in_iv marks the rows the model is FITTED on (the January stage-1 IV set, already-fighting at January removed).
oppw_frame <- function(alts, panel, events_dymid, jp, hist_window = OPPW_HIST_WINDOW) {
  stopifnot(all(OPPW_ALTS_COLS %in% names(alts)), all(OPPW_PANEL_COLS %in% names(panel)), all(c("initiator", "target", "year", "hostlev1", "hostlev2") %in% names(events_dymid)))
  fr <- as.data.frame(alts)[, OPPW_ALTS_COLS]              # the only alternatives columns read
  pn <- as.data.frame(panel)[, OPPW_PANEL_COLS]            # the only panel columns read
  for (cc in c("ccode1", "ccode2", "year", "in_iv", "in_s2")) fr[[cc]] <- as.integer(fr[[cc]])
  assert_no_dup_keys(fr, c("ccode1", "ccode2", "year"), "oppw_frame: alternatives"); assert_no_dup_keys(pn, c("ccode", "year"), "oppw_frame: panel")
  stopifnot(all(fr$in_iv %in% 0:1), all(fr$in_s2 %in% 0:1), all(fr$in_iv == 1L | fr$in_s2 == 1L), !anyNA(fr$contig_land), !anyNA(fr$contig_any),
            all(fr$contig_land %in% 0:1), all(fr$contig_any %in% 0:1), all(fr$contig_land <= fr$contig_any))
  n_nodist <- sum(is.na(fr$log_capdist)); fr <- fr[!is.na(fr$log_capdist), , drop = FALSE]
  fr$contig_sea <- as.integer(fr$contig_any - fr$contig_land); fr$contig_any <- NULL
  m <- match(oppw_key(fr$ccode1, 1L, fr$year), oppw_key(pn$ccode, 1L, pn$year))
  if (anyNA(m)) stop(sprintf("oppw_frame: %d alternative rows whose initiator state-year is not in the stage-1 panel", sum(is.na(m))))
  fr$M_lvl_i <- pn$M_lvl[m]
  ## outcomes: undirected DYMID onsets in year t
  ev <- as.data.frame(events_dymid)[, c("initiator", "target", "year", "hostlev1", "hostlev2")]
  stopifnot(!anyNA(ev$initiator), !anyNA(ev$target), !anyNA(ev$year), !anyNA(ev$hostlev1), !anyNA(ev$hostlev2), all(ev$initiator != ev$target))
  ev$a <- pmin(ev$initiator, ev$target); ev$b <- pmax(ev$initiator, ev$target); ev$host3 <- pmax(ev$hostlev1, ev$hostlev2) >= 3L
  k_fr <- oppw_key(pmin(fr$ccode1, fr$ccode2), pmax(fr$ccode1, fr$ccode2), fr$year)
  fr$y_mid   <- as.integer(k_fr %in% oppw_key(ev$a, ev$b, ev$year))
  fr$y_host3 <- as.integer(k_fr %in% oppw_key(ev$a[ev$host3], ev$b[ev$host3], ev$year[ev$host3]))
  stopifnot(all(fr$y_host3 <= fr$y_mid))
  ## history: an onset in [t - window, t - 1], never t (mid3_history() asserts it)
  h <- mid3_history(events_dymid, hist_window)
  fr$dymid_hist3 <- as.integer(oppw_key(fr$ccode1, fr$ccode2, fr$year) %in% oppw_key(h$ccode1, h$ccode2, h$year))
  fr$a4_joinable <- oppw_a4(fr, jp)
  fr$decade <- 10L * (fr$year %/% 10L)
  fr <- fr[order(fr$year, fr$ccode1, fr$ccode2), , drop = FALSE]; rownames(fr) <- NULL
  stopifnot(!anyNA(fr), all(is.finite(fr$log_capdist)), all(is.finite(fr$M_lvl_i)), !any(OPPW_TARGET_STRENGTH %in% names(fr)))
  attr(fr, "n_nodist") <- n_nodist
  fr
}

## ---- one window: design and fit ----------------------------------------------------------------------------------------------------------------------
## oppw_decade_levels(train_year, train_y): the decades that get their own effect in a training window, in order. A decade qualifies when the window
## holds at least OPPW_DECADE_MIN_YEARS of its years and at least one onset in it; a decade that does not (the first years of a new decade, a decade without
## an onset: its effect would be unidentified or minus infinity) is MERGED INTO THE LATEST EARLIER QUALIFYING DECADE (the earliest decades into the first
## qualifying one). Returns list(levels = qualifying decades, map = function(decade) -> level), where the map also carries a decade LATER than every
## training decade (the prediction year at the start of a new decade) to the last level.
oppw_decade_levels <- function(train_year, train_y, min_years = OPPW_DECADE_MIN_YEARS) {
  dec <- 10L * (train_year %/% 10L)
  tab <- data.frame(decade = sort(unique(dec)))
  tab$n_years <- vapply(tab$decade, function(d) length(unique(train_year[dec == d])), integer(1))
  tab$n_y <- vapply(tab$decade, function(d) sum(train_y[dec == d]), numeric(1))
  ok <- tab$decade[tab$n_years >= min_years & tab$n_y > 0]
  if (length(ok) == 0L) stop("oppw_decade_levels: no decade of the training window has ", min_years, " years and an onset")
  map <- function(decade) { i <- findInterval(decade, ok); ok[pmax(i, 1L)] }
  list(levels = ok, map = map, table = tab)
}
## oppw_design(fr, a_cols, dl, centre): the model matrix of the rows `fr` for one window: intercept, the A columns standardised with the TRAINING window's
## means and SDs (centre$mean / centre$sd; the product log_capdist_x_init_M = (log_capdist - the training mean of log_capdist) x M_lvl_i is formed first),
## then one 0/1 column per decade level after the first. centre = NULL computes the constants from `fr` itself (the training call) and returns them.
oppw_design <- function(fr, a_cols, dl, centre = NULL) {
  stopifnot(all(a_cols %in% OPPW_A_COLS), !any(OPPW_TARGET_STRENGTH %in% a_cols))
  ld_mean <- if (is.null(centre)) mean(fr$log_capdist) else centre$ld_mean
  raw <- vapply(a_cols, function(cc) if (cc == "log_capdist_x_init_M") (fr$log_capdist - ld_mean) * fr$M_lvl_i else as.numeric(fr[[cc]]), numeric(nrow(fr)))
  if (is.null(centre)) {
    centre <- list(ld_mean = ld_mean, mean = colMeans(raw), sd = apply(raw, 2, stats::sd))
    if (any(!is.finite(centre$sd)) || any(centre$sd <= 0)) stop("oppw_design: constant predictor in the training window: ", paste(a_cols[!(centre$sd > 0)], collapse = ", "))
  }
  X <- sweep(sweep(raw, 2, centre$mean, "-"), 2, centre$sd, "/")
  lev <- dl$map(fr$decade)
  D <- if (length(dl$levels) > 1L) vapply(dl$levels[-1], function(d) as.numeric(lev == d), numeric(nrow(fr))) else matrix(0, nrow(fr), 0L)
  X <- cbind(1, X, D); colnames(X) <- c("(Intercept)", a_cols, if (length(dl$levels) > 1L) paste0("decade_", dl$levels[-1]))
  list(X = X, centre = centre)
}
## oppw_fit(X, y, start): the logit MLE and its asymptotic covariance (X' W X)^-1. Newton / IRLS on the normal equations (crossprod), which is what
## stats::glm.fit(family = binomial()) iterates through a QR decomposition; several times faster on the 1.7 million rows x 26 columns of a late window.
## The same maximum: oppw_expanding() checks the first and the last window against stats::glm.fit (oppw_fit_glm()) and stops on a difference, and
## test_opp_weights.R checks it on toy data. Step halving when the log likelihood falls; stops unless converged, finite and without separation.
oppw_loglik <- function(eta, y) sum(ifelse(y == 1, stats::plogis(eta, log.p = TRUE), stats::plogis(-eta, log.p = TRUE)))
oppw_fit <- function(X, y, start = NULL, maxit = 100L, tol = 1e-10) {
  stopifnot(is.matrix(X), nrow(X) == length(y), all(y %in% 0:1), sum(y) > 0, sum(y) < length(y), all(X[, 1] == 1))
  b <- if (is.null(start)) c(stats::qlogis(mean(y)), numeric(ncol(X) - 1L)) else start
  stopifnot(length(b) == ncol(X), all(is.finite(b)))
  eta <- drop(X %*% b); ll <- oppw_loglik(eta, y); conv <- FALSE; it <- 0L
  while (it < maxit) {
    it <- it + 1L
    mu <- stats::plogis(eta); w <- mu * (1 - mu)
    H <- crossprod(X * sqrt(w)); g <- crossprod(X, y - mu)
    step <- drop(tryCatch(chol2inv(chol(H)) %*% g, error = function(e) stop("oppw_fit: the information matrix is not positive definite (collinear design or separation)")))
    half <- 0L
    repeat {
      b_new <- b + step; eta_new <- drop(X %*% b_new); ll_new <- oppw_loglik(eta_new, y)
      if (is.finite(ll_new) && ll_new >= ll - 1e-10 * abs(ll)) break
      step <- step / 2; half <- half + 1L
      if (half > 30L) stop("oppw_fit: step halving failed")
    }
    done <- max(abs(b_new - b)) < tol * (1 + max(abs(b_new)))
    b <- b_new; eta <- eta_new; ll <- ll_new
    if (done) { conv <- TRUE; break }
  }
  if (!conv) stop("oppw_fit: the auxiliary logit did not converge in ", maxit, " iterations")
  if (max(abs(b)) > 30) stop("oppw_fit: a coefficient beyond 30 in absolute value (separation)")
  mu <- stats::plogis(eta); V <- chol2inv(chol(crossprod(X * sqrt(mu * (1 - mu)))))
  stopifnot(all(is.finite(V)), all(diag(V) > 0))
  list(coef = unname(b), V = V, n = length(y), n_y = sum(y), iter = it, loglik = ll)
}
## oppw_fit_glm(X, y): the same fit by stats::glm.fit (the reference implementation; coefficients and covariance)
oppw_fit_glm <- function(X, y) {
  f <- suppressWarnings(stats::glm.fit(X, y, family = stats::binomial(), control = stats::glm.control(epsilon = 1e-12, maxit = 100L)))
  if (!isTRUE(f$converged) || isTRUE(f$boundary) || anyNA(f$coefficients) || f$rank < ncol(X)) stop("oppw_fit_glm: glm.fit did not converge to an interior full-rank solution")
  R <- qr.R(f$qr)[seq_len(ncol(X)), seq_len(ncol(X)), drop = FALSE]
  V <- chol2inv(R); V <- V[order(f$qr$pivot), order(f$qr$pivot), drop = FALSE]
  list(coef = unname(f$coefficients), V = V)
}

## ---- the expanding window ----------------------------------------------------------------------------------------------------------------------------
## oppw_windows(years, burn_in): the fits. Window 1 trains on the years <= burn_in and predicts them (IN SAMPLE) and the year burn_in + 1 (held out);
## every later year t has its own window trained on the years < t. Returns data.frame(window, train_max, pred_from, pred_to).
oppw_windows <- function(years, burn_in = OPPW_BURN_IN) {
  years <- sort(unique(as.integer(years)))
  stopifnot(length(years) >= 2L, min(years) <= burn_in, max(years) > burn_in)
  later <- years[years > burn_in + 1L]
  w <- data.frame(window = seq_len(1L + length(later)), train_max = c(burn_in, later - 1L), pred_from = c(min(years), later), pred_to = c(burn_in + 1L, later))
  stopifnot(all(w$train_max[-1] < w$pred_from[-1]), w$train_max[1] == burn_in)
  w
}
## oppw_expanding(fr, outcome, a_cols, burn_in, n_draws, seed): the expanding-window auxiliary logit.
## For every window: the model is fitted on the IV rows (in_iv = 1) of the training years ONLY (the training subset is cut before anything is computed
## from it: the centring constants, the decade levels and the coefficients see no row of a later year), then every row (IV or stage-2 only) of the
## window's prediction years gets ln q = log inv_logit(x' b), computed on the log scale (never -Inf).
## Uncertainty: n_draws standard-normal vectors z_d (one per draw, length = the largest coefficient count, set.seed(seed)) are drawn ONCE; in every window
## the d-th coefficient vector is b + L z_d[1:p] with L the lower Cholesky factor of that window's asymptotic covariance. Each draw is therefore a draw of the
## FULL coefficient vector of every window, with common random numbers across windows, and gives one coherent ln q column for all rows (never per-row draws).
## Returns list(lq = matrix nrow(fr) x (1 + n_draws), columns d0 (MLE), d1..; windows = the window table with fit statistics; coefs = long table of the
## MLE coefficients per window; held_out = logical per row (year > burn_in); z = the draw vectors). glm_check = TRUE refits the first and the last window
## with stats::glm.fit and stops unless coefficients and covariance agree.
oppw_expanding <- function(fr, outcome = "y_mid", a_cols = OPPW_A_COLS, burn_in = OPPW_BURN_IN, n_draws = OPPW_N_DRAWS, seed = V2_SEED, verbose = FALSE, glm_check = TRUE) {
  stopifnot(outcome %in% OPPW_OUTCOMES, all(c(outcome, "year", "in_iv", "decade", "log_capdist", "M_lvl_i", setdiff(a_cols, "log_capdist_x_init_M")) %in% names(fr)),
            !is.unsorted(fr$year), n_draws >= 0L, !any(OPPW_TARGET_STRENGTH %in% a_cols))
  w <- oppw_windows(fr$year, burn_in)
  p_max <- 1L + length(a_cols) + length(unique(fr$decade))
  set.seed(seed); z <- matrix(stats::rnorm(n_draws * p_max), nrow = n_draws, ncol = p_max)       # row d = the d-th draw's standard-normal vector
  lq <- matrix(NA_real_, nrow(fr), 1L + n_draws, dimnames = list(NULL, paste0("d", 0:n_draws)))
  w$n_train <- NA_integer_; w$n_y <- NA_real_; w$p <- NA_integer_; w$n_decades <- NA_integer_; w$train_year_max <- NA_integer_; w$iter <- NA_integer_; w$glm_max_diff <- NA_real_
  coefs <- vector("list", nrow(w)); start <- NULL; prev_names <- NULL
  for (i in seq_len(nrow(w))) {
    tr <- fr[fr$year <= w$train_max[i] & fr$in_iv == 1L, , drop = FALSE]                           # the training subset: nothing below sees another row
    stopifnot(nrow(tr) > 0L, max(tr$year) <= w$train_max[i], i == 1L || max(tr$year) < w$pred_from[i])
    dl <- oppw_decade_levels(tr$year, tr[[outcome]])
    dz <- oppw_design(tr, a_cols, dl)
    st <- if (!is.null(start) && all(prev_names %in% colnames(dz$X))) { s <- stats::setNames(numeric(ncol(dz$X)), colnames(dz$X)); s[prev_names] <- start; unname(s) } else NULL
    f <- oppw_fit(dz$X, tr[[outcome]], start = st)
    if (glm_check && (i == 1L || i == nrow(w))) {                                                   # the first and the last window against stats::glm.fit
      fg <- oppw_fit_glm(dz$X, tr[[outcome]])
      ## coefficients to 1e-6; the covariance to 1e-4 relative (glm.fit forms it from the weights of its last-but-one iterate)
      if (max(abs(fg$coef - f$coef)) > 1e-6 || max(abs(fg$V - f$V)) > 1e-4 * max(abs(fg$V)))
        stop(sprintf("oppw_expanding: window %d differs from stats::glm.fit (max coefficient difference %.3g, max relative covariance difference %.3g)", i, max(abs(fg$coef - f$coef)), max(abs(fg$V - f$V)) / max(abs(fg$V))))
      w$glm_max_diff[i] <- max(abs(fg$coef - f$coef))
    }
    start <- f$coef; prev_names <- colnames(dz$X)
    pr <- which(fr$year >= w$pred_from[i] & fr$year <= w$pred_to[i])
    Xp <- oppw_design(fr[pr, , drop = FALSE], a_cols, dl, centre = dz$centre)$X
    stopifnot(identical(colnames(Xp), colnames(dz$X)), all(is.na(lq[pr, 1L])))
    p <- length(f$coef); L <- t(chol(f$V))
    B <- cbind(f$coef, if (n_draws > 0L) f$coef + L %*% t(z[, seq_len(p), drop = FALSE]))          # p x (1 + n_draws): MLE, then the full-vector draws
    lq[pr, ] <- stats::plogis(Xp %*% B, log.p = TRUE)
    w$n_train[i] <- nrow(tr); w$n_y[i] <- f$n_y; w$p[i] <- p; w$n_decades[i] <- length(dl$levels); w$train_year_max[i] <- max(tr$year); w$iter[i] <- f$iter
    coefs[[i]] <- data.frame(window = i, train_max = w$train_max[i], term = colnames(dz$X), coef = f$coef, se = sqrt(diag(f$V)))
    if (verbose && (i %% 20L == 1L || i == nrow(w))) logf("  window %d / %d: years <= %d (%d IV rows, %.0f onsets, %d coefficients) -> predicts %d..%d", i, nrow(w), w$train_max[i], nrow(tr), f$n_y, p, w$pred_from[i], w$pred_to[i])
    rm(tr, dz, Xp)
  }
  stopifnot(!anyNA(lq), all(is.finite(lq)), all(lq <= 0), all(w$train_year_max[-1] < w$pred_from[-1]))
  list(lq = lq, windows = w, coefs = do.call(rbind, coefs), held_out = fr$year > burn_in, z = z, a_cols = a_cols, outcome = outcome, burn_in = burn_in, seed = seed)
}

## ---- auxiliary diagnostics (no AUC: uninformative at a 0.2% positive rate) -----------------------------------------------------------------------------
## oppw_logloss(y, lq): mean negative log likelihood of the 0/1 outcome under q = exp(lq)
oppw_logloss <- function(y, lq) { stopifnot(length(y) == length(lq), all(lq <= 0)); -mean(ifelse(y == 1L, lq, log1p(-exp(lq)))) }
## oppw_pr(y, lq): precision-recall summary: prevalence, average precision (the step-wise area under the PR curve) and the precision at recall 0.25 / 0.5 / 0.75
oppw_pr <- function(y, lq) {
  stopifnot(length(y) == length(lq), sum(y) > 0)
  o <- order(lq, decreasing = TRUE); ys <- y[o]; tp <- cumsum(ys); prec <- tp / seq_along(ys); rec <- tp / sum(ys)
  at <- function(r) prec[which(rec >= r)[1]]
  data.frame(n = length(y), n_y = sum(y), prevalence = mean(y), avg_precision = sum(prec[ys == 1L]) / sum(ys), prec_at_rec25 = at(0.25), prec_at_rec50 = at(0.5), prec_at_rec75 = at(0.75))
}
## oppw_calibration(y, lq, group, what): one row per group: rows, onsets, observed rate, mean q, observed / predicted, log loss
oppw_calibration <- function(y, lq, group, what) {
  stopifnot(length(y) == length(lq), length(group) == length(y), !anyNA(group))
  g <- split(seq_along(y), group)
  out <- do.call(rbind, lapply(names(g), function(nm) { i <- g[[nm]]
    data.frame(by = what, group = nm, n = length(i), n_y = sum(y[i]), obs_rate = mean(y[i]), mean_q = mean(exp(lq[i])), obs_over_pred = mean(y[i]) / mean(exp(lq[i])), logloss = oppw_logloss(y[i], lq[i])) }))
  rownames(out) <- NULL; out
}

## ---- the weights file and the stage-2 offset ----------------------------------------------------------------------------------------------------------
## spec_opp_*(spec): the three run-1 switches ("none" / "0" / "mid" for a spec saved before 2026-11-10, which has none of the fields)
## oppw_foundation(spec): the name of the cell foundation's data directory (data_entry | data_dymid | data), the key of the weights file
oppw_foundation <- function(spec) basename(spec_data_dir(spec, "."))
## oppw_path(spec, out_dir): the weights file of a spec, output/oppw_<foundation>_<outcome>_lag3.rds (written by 08_opportunity_weights.R)
oppw_path <- function(spec, out_dir) file.path(out_dir, sprintf("oppw_%s_%s_lag%d.rds", oppw_foundation(spec), spec_opp_outcome(spec), OPPW_LAG))
## oppw_offset(alts, W, lambda, draw, manifest): ln q of every assembled alternative row and the offset lambda * ln q.
##   alts     the assembled alternatives (ccode1, ccode2, year), in the row order of the Stan data
##   W        the weights object of 08_opportunity_weights.R: keys (ccode1, ccode2, year), lq (rows x d0..d5), input (file, sha256 of the foundation files)
##   lambda   0 | 0.5 | 1;  draw 0 (the auxiliary MLE) | 1..5 (a full-vector draw; lambda = 1 only)
##   manifest the cell's data manifest: the weights must have been built from the SAME alternatives and panel files (sha256), or the run stops
## lambda = 0 returns an exact zero vector (the Stan data and the log density of the cell without the switch).
oppw_offset <- function(alts, W, lambda, draw = 0L, manifest = NULL) {
  stopifnot(length(lambda) == 1L, lambda %in% c(0, 0.5, 1), length(draw) == 1L, draw %in% 0:OPPW_N_DRAWS, draw == 0L || lambda == 1,
            all(c("ccode1", "ccode2", "year") %in% names(alts)), all(c("keys", "lq", "input") %in% names(W)), nrow(W$keys) == nrow(W$lq), paste0("d", draw) %in% colnames(W$lq))
  if (!is.null(manifest)) {
    mm <- merge(W$input[W$input$role == "foundation", c("file", "sha256")], manifest[, c("file", "sha256")], by = "file", suffixes = c("_w", "_now"))
    if (nrow(mm) != sum(W$input$role == "foundation") || any(mm$sha256_w != mm$sha256_now))
      stop("oppw_offset: the opportunity weights were built from other data files than this cell's (sha256 of ", paste(W$input$file[W$input$role == "foundation"], collapse = ", "), "); rerun 08_opportunity_weights.R")
  }
  m <- match(oppw_key(alts$ccode1, alts$ccode2, alts$year), oppw_key(W$keys$ccode1, W$keys$ccode2, W$keys$year))
  if (anyNA(m)) stop(sprintf("oppw_offset: %d alternative rows of the cell have no opportunity weight", sum(is.na(m))))
  lq <- unname(W$lq[m, paste0("d", draw)])
  stopifnot(length(lq) == nrow(alts), all(is.finite(lq)), all(lq <= 0))
  list(lq = lq, off = if (lambda == 0) numeric(length(lq)) else lambda * lq)
}
## oppw_neff(q): c(sum_q, n_eff) with n_eff = (sum q)^2 / sum q^2 (the effective number of alternatives under the weights; between 1 and length(q))
oppw_neff <- function(q) { stopifnot(length(q) >= 1L, all(is.finite(q)), all(q > 0)); c(sum_q = sum(q), n_eff = sum(q)^2 / sum(q^2)) }
## oppw_cell_stats(dat, lq, group_labels, low_prob): the per-cell diagnostics, from the Stan data and ln q of its rows.
##   threshold  the low_prob quantile of q over ALL alternative rows of the cell (pooled over state-years; "low q" = strictly below it)
##   by_role    one row per stage-2 role, then "all_events": usable events (count, composite-weight sum), the events whose CHOSEN target has low q (count and
##              share of events; weight sum and share of the composite weight), the median q of the chosen targets; and the medians over the role's distinct
##              initiating state-years of sum q, N_eff and the set size, all over the state-year's IV rows (the rows the IV sums over). A state-year with
##              events of two roles counts in both. "non_initiators" (y1 = 0) and "all_state_years" carry the state-year columns only.
oppw_cell_stats <- function(dat, lq, group_labels, low_prob = OPPW_LOW_Q) {
  stopifnot(length(lq) == dat$N2, all(is.finite(lq)), all(lq <= 0), length(group_labels) == dat$G, length(dat$ev_grp) == dat$E)
  q <- exp(lq); thr <- unname(stats::quantile(q, low_prob, names = FALSE))
  sy_of_row <- rep(seq_len(dat$N1), times = dat$cs_end - dat$cs_start + 1L)
  stopifnot(length(sy_of_row) == dat$N2)
  iv <- dat$in_iv == 1L
  sq <- as.numeric(tapply(q[iv], factor(sy_of_row[iv], levels = seq_len(dat$N1)), sum)); sq2 <- as.numeric(tapply(q[iv]^2, factor(sy_of_row[iv], levels = seq_len(dat$N1)), sum))
  n_iv <- tabulate(sy_of_row[iv], nbins = dat$N1); neff <- sq^2 / sq2
  stopifnot(!anyNA(sq), all(n_iv >= 1L), all(neff >= 1 - 1e-9), all(neff <= n_iv + 1e-9))
  q_ch <- q[dat$ev_row]; low <- q_ch < thr
  sy_row <- function(label, idx) data.frame(role = label, n_sy = length(idx), sum_q_median = stats::median(sq[idx]), n_eff_median = stats::median(neff[idx]), n_iv_median = stats::median(n_iv[idx]))
  ev_row <- function(label, e) data.frame(role = label, n_events = length(e), w_events = sum(dat$ev_w[e]), n_low = sum(low[e]), share_low_n = mean(low[e]),
                                          w_low = sum(dat$ev_w[e][low[e]]), share_low_w = sum(dat$ev_w[e][low[e]]) / sum(dat$ev_w[e]), q_chosen_median = stats::median(q_ch[e]))
  rows <- lapply(seq_len(dat$G), function(g) { e <- which(dat$ev_grp == g - 1L); cbind(ev_row(group_labels[g], e), sy_row(group_labels[g], sort(unique(dat$ev_sy[e])))[, -1]) })
  rows[[length(rows) + 1L]] <- cbind(ev_row("all_events", seq_len(dat$E)), sy_row("all_events", sort(unique(dat$ev_sy)))[, -1])
  blank <- rows[[1]][1, 2:8]; blank[] <- NA
  for (r in list(sy_row("non_initiators", which(dat$y1 == 0L)), sy_row("all_state_years", seq_len(dat$N1)))) rows[[length(rows) + 1L]] <- cbind(role = r$role, blank, r[, -1])
  by_role <- do.call(rbind, rows); rownames(by_role) <- NULL
  list(low_prob = low_prob, threshold = thr, by_role = by_role, q_summary = stats::quantile(q, c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)))
}

## ---- the mixture over the auxiliary draws (lambda = 1) -------------------------------------------------------------------------------------------------
## oppw_mixture(draw_list, point): the equally weighted mixture of the draw-conditional posteriors. draw_list = one matrix of posterior draws (iterations x
## parameters, same columns, SAME number of rows: equal weights by construction) per auxiliary draw; point = the matrix of the fit at the auxiliary MLE.
## Returns one row per parameter: the point fit's mean and SD, the mixture's mean, SD and 5 / 50 / 95% quantiles, the SD of the conditional means across the
## auxiliary draws (the between-draw part), the ratio mixture SD / point SD and exceeds_10pct = ratio > 1.10.
oppw_mixture <- function(draw_list, point, tol = 1.10) {
  stopifnot(is.list(draw_list), length(draw_list) >= 2L, is.matrix(point), all(vapply(draw_list, is.matrix, logical(1))))
  cn <- colnames(point); stopifnot(!is.null(cn), all(vapply(draw_list, function(m) identical(colnames(m), cn), logical(1))))
  nr <- vapply(draw_list, nrow, integer(1))
  if (length(unique(nr)) != 1L) stop("oppw_mixture: the conditional fits hold different numbers of posterior draws (", paste(nr, collapse = ", "), "); the mixture would not be equally weighted")
  mix <- do.call(rbind, draw_list)
  cm <- vapply(draw_list, colMeans, numeric(length(cn)))
  if (length(cn) == 1L) cm <- matrix(cm, nrow = 1L)
  out <- data.frame(parameter = cn, point_mean = colMeans(point), point_sd = apply(point, 2, stats::sd), mix_mean = colMeans(mix), mix_sd = apply(mix, 2, stats::sd),
                    mix_q05 = apply(mix, 2, stats::quantile, 0.05), mix_q50 = apply(mix, 2, stats::quantile, 0.5), mix_q95 = apply(mix, 2, stats::quantile, 0.95),
                    between_draw_sd = apply(cm, 1, stats::sd), n_aux_draws = length(draw_list), n_post_draws_each = nr[1])
  out$sd_ratio <- out$mix_sd / out$point_sd; out$exceeds_10pct <- out$sd_ratio > tol
  rownames(out) <- NULL
  list(table = out, mix = mix)
}
