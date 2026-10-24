## ============================================================================
## 00_helpers.R -- small shared helpers for the v2 pipeline
## Sourced by every script. Keep this file small and concrete.
## ============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

V2_SEED <- 123L

## ---- paths -----------------------------------------------------------------
v2_root <- function() {
  root <- Sys.getenv("V2_ROOT", unset = "")
  if (root == "") {
    ## fall back: directory that contains this file's parent (v2/)
    args <- commandArgs(trailingOnly = FALSE)
    fa <- sub("^--file=", "", args[grepl("^--file=", args)])
    if (length(fa) == 1) root <- normalizePath(file.path(dirname(fa), ".."))
    else root <- normalizePath("v2", mustWork = FALSE)
  }
  root
}

## ---- logging -----------------------------------------------------------------
.v2_log_file <- NULL
log_open <- function(name) {
  dir.create(file.path(v2_root(), "logs"), showWarnings = FALSE, recursive = TRUE)
  .v2_log_file <<- file.path(v2_root(), "logs",
                             paste0(name, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
  cat("", file = .v2_log_file)
  logf("log opened: %s", .v2_log_file)
  invisible(.v2_log_file)
}
logf <- function(fmt, ...) {
  msg <- sprintf(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", fmt), ...)
  cat(msg, "\n", sep = "")
  if (!is.null(.v2_log_file)) cat(msg, "\n", sep = "", file = .v2_log_file, append = TRUE)
  invisible(msg)
}

## ---- assertions ------------------------------------------------------------------
assert_no_dup_keys <- function(df, keys, what = "data") {
  n_dup <- nrow(df) - nrow(distinct(df[, keys, drop = FALSE]))
  if (n_dup != 0) stop(sprintf("%s: %d duplicate rows on keys (%s)", what, n_dup, paste(keys, collapse = ", ")))
  invisible(TRUE)
}
assert_no_na <- function(df, cols, what = "data") {
  bad <- cols[vapply(cols, function(cc) anyNA(df[[cc]]), logical(1))]
  if (length(bad) > 0) stop(sprintf("%s: NA present in model columns: %s", what, paste(bad, collapse = ", ")))
  invisible(TRUE)
}
assert_finite <- function(df, cols, what = "data") {
  bad <- cols[vapply(cols, function(cc) !all(is.finite(df[[cc]])), logical(1))]
  if (length(bad) > 0) stop(sprintf("%s: non-finite values (NA / NaN / Inf) in columns: %s", what, paste(bad, collapse = ", ")))
  invisible(TRUE)
}
assert_standardized <- function(df, col, by = "year", tol = 1e-6, what = "data") {
  chk <- df |> group_by(across(all_of(by))) |>
    summarise(m = mean(.data[[col]]), s = sd(.data[[col]]), n = n(), .groups = "drop") |>
    filter(n > 1)
  if (any(abs(chk$m) > tol) || any(abs(chk$s - 1) > 1e-3))
    stop(sprintf("%s: %s is not standardized within %s (max |mean| = %.2e, sd range [%.4f, %.4f])",
                 what, col, by, max(abs(chk$m)), min(chk$s), max(chk$s)))
  invisible(TRUE)
}

## ---- power changes (PI decision 2026-10-23; see README "Positions" and "Changes") -----------
## pos  : state-year levels with columns ccode, year, M (log CINC share), R_pr (log IGO PageRank)
## at_t : the (ccode, year) state-years at t for which a change is wanted (the state-year panel)
## h    : window length. The change at (ccode, t) is level(t-1) - level(t-1-h), by CALENDAR joins.
## Returns list(tbl, scale). tbl has one row per (ccode, year) in at_t with levels at both dates:
##   dM_raw_h, dRpr_raw_h : the raw change in log units (log ratio of the state's own share across the window)
##   dM_h, dRpr_h         : the raw change minus the mean raw change of all states with a change in that year,
##                          divided by ONE pooled constant per lag and variable (the SD of the year-centred change
##                          over all rows of the table). Zero means "changed like the average state that year";
##                          the pooled SD is 1 by construction.
## scale is a data.frame(h, variable, pooled_sd, n_rows, n_years) with the two constants.
## The same table serves the initiator change in stage 1 (W, X) and the target change in stage 2 (Z).
build_change_table <- function(pos, at_t, h) {
  stopifnot(all(c("ccode", "year", "M", "R_pr") %in% names(pos)), all(c("ccode", "year") %in% names(at_t)),
            length(h) == 1, h >= 1)
  h <- as.integer(h)
  assert_no_dup_keys(pos, c("ccode", "year"), "change table: positions")
  assert_finite(pos, c("M", "R_pr"), "change table: positions")
  tbl <- at_t |> distinct(ccode, year) |>
    inner_join(pos |> transmute(ccode, year = year + 1L, .M_tm1 = M, .R_tm1 = R_pr), by = c("ccode", "year")) |>
    inner_join(pos |> transmute(ccode, year = year + 1L + h, .M_w0 = M, .R_w0 = R_pr), by = c("ccode", "year")) |>
    mutate(dM_raw_h = .M_tm1 - .M_w0, dRpr_raw_h = .R_tm1 - .R_w0) |>
    group_by(year) |>
    mutate(.dM_c = dM_raw_h - mean(dM_raw_h), .dR_c = dRpr_raw_h - mean(dRpr_raw_h)) |> ungroup()
  if (nrow(tbl) < 3) stop(sprintf("change table h=%d: only %d state-years have levels at t-1 and t-1-h", h, nrow(tbl)))
  sd_M <- sd(tbl$.dM_c); sd_R <- sd(tbl$.dR_c)
  if (!is.finite(sd_M) || sd_M <= 0 || !is.finite(sd_R) || sd_R <= 0)
    stop(sprintf("change table h=%d: pooled SD of the year-centred change is not positive (M %.3g, R_pr %.3g)", h, sd_M, sd_R))
  tbl <- tbl |> mutate(dM_h = .dM_c / sd_M, dRpr_h = .dR_c / sd_R, h = h) |>
    select(ccode, year, h, dM_raw_h, dRpr_raw_h, dM_h, dRpr_h) |> arrange(year, ccode)
  assert_no_dup_keys(tbl, c("ccode", "year"), sprintf("change table h=%d", h))
  assert_finite(tbl, c("dM_raw_h", "dRpr_raw_h", "dM_h", "dRpr_h"), sprintf("change table h=%d", h))
  ## invariants: within-year mean 0 (tolerance), pooled SD exactly 1 after scaling
  yr <- tbl |> group_by(year) |> summarise(mM = mean(dM_h), mR = mean(dRpr_h), .groups = "drop")
  if (any(abs(yr$mM) > 1e-8) || any(abs(yr$mR) > 1e-8))
    stop(sprintf("change table h=%d: within-year mean of the centred change is not 0 (max |mean| M %.2e, R_pr %.2e)", h, max(abs(yr$mM)), max(abs(yr$mR))))
  if (abs(sd(tbl$dM_h) - 1) > 1e-8 || abs(sd(tbl$dRpr_h) - 1) > 1e-8)
    stop(sprintf("change table h=%d: pooled SD after scaling is not 1 (M %.6f, R_pr %.6f)", h, sd(tbl$dM_h), sd(tbl$dRpr_h)))
  list(tbl = tbl,
       scale = data.frame(h = h, variable = c("dM_h", "dRpr_h"), pooled_sd = c(sd_M, sd_R),
                          n_rows = nrow(tbl), n_years = n_distinct(tbl$year)))
}

## ---- fatal-MID coding (PI decision 2026-10-23; see README "Samples") -----------------
## fatal_cow : COW coding (TRUE / FALSE / NA when every COW fatality code is missing)
## fatal_mic : MIC 1.0 confrontation fatalmin >= 1 (NA when the dispute is not in MIC)
## war       : COW dispute hostility level 5
## fatal_gml : GML MID 2.2.1 coding of the same directed pair (NA when missing or not in GML)
## Main fmid = war (fatal by definition) OR COW fatal OR (COW missing AND MIC fatal)
##             OR (COW missing AND MIC missing AND GML fatal).
## Each source only fills cases the earlier ones leave missing; none overrides a coded value.
## Returns logical, never NA.
code_fatal_main <- function(fatal_cow, fatal_mic, war, fatal_gml) {
  n <- length(fatal_cow)
  stopifnot(length(fatal_mic) == n, length(war) == n, length(fatal_gml) == n)
  (war %in% TRUE) | (fatal_cow %in% TRUE) | (is.na(fatal_cow) & fatal_mic %in% TRUE) |
    (is.na(fatal_cow) & is.na(fatal_mic) & fatal_gml %in% TRUE)
}
MIC_CONF_SHA256 <- "41cef02d0dcb23f91658b2367864b018a517e13d5b0278394dc12775fa5a5ee0"   # mic-conf-1.0.csv

## ---- manifest / digests ------------------------------------------------------------
file_digest <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
write_manifest <- function(entries, path) {
  ## entries: list of data.frames with file, nrow, ncol, sha256, note
  man <- bind_rows(entries) |> mutate(built = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  write.csv(man, path, row.names = FALSE)
  invisible(man)
}
assert_matches_manifest <- function(path, manifest_path) {
  man <- read.csv(manifest_path, stringsAsFactors = FALSE)
  row <- man[man$file == basename(path), , drop = FALSE]
  if (nrow(row) != 1) stop("file not in manifest: ", basename(path))
  obj <- readRDS(path)
  if (nrow(obj) != row$nrow) stop(sprintf("%s: nrow %d != manifest %d", basename(path), nrow(obj), row$nrow))
  if (file_digest(path) != row$sha256) stop(sprintf("%s: sha256 digest does not match manifest", basename(path)))
  obj
}

## ---- Stan data assembly ----------------------------------------------------------------
## Builds the Stan data list from three data frames:
##   sy:   state-years, one row per (state_idx, year); columns y1, state_idx, W cols, X cols
##   alts: alternative rows, one row per (sy_row, target_idx); columns sy_row, target_idx, dM, dR, Z cols
##   ev:   events, columns sy_row, alt_row (index into alts after ordering by sy_row)
## alts must be sorted by sy_row; every sy row must own >= 1 alternative.
## prior_mu0_sd: sd of the normal(0, .) prior on the mu intercept mu0 (the population-average ideal gap, in the
## units of x; since 2026-10-23 x is a log ratio of CINC shares / PageRank, no longer within-year SD units);
## 2 is the README value and was set when x was in SD units (flagged for PI review, not changed here).
build_stan_data <- function(sy, alts, ev, W_cols, X_cols, Z_cols, joint, sampled = 0L,
                            n_full = NULL, x_ref = NULL, K = NULL, anchor = NULL,
                            prior_beta_sd = 1, prior_mu0_sd = 2, allow_orphan_initiations = FALSE, S = NULL, prior_only = 0L,
                            D_mode = 3L, grainsize = NULL) {
  stopifnot(is.data.frame(sy), is.data.frame(alts), is.data.frame(ev))
  stopifnot(is.numeric(prior_mu0_sd), length(prior_mu0_sd) == 1, is.finite(prior_mu0_sd), prior_mu0_sd > 0)
  sampled <- as.integer(sampled)
  if (sampled == 1) {
    ## alts = IV rows (first n_iv, sorted by sy_row) followed by per-event blocks; ev carries ev_start/ev_end
    stopifnot(all(c("ev_start", "ev_end") %in% names(ev)), !is.null(n_full))
    n_iv <- min(ev$ev_start) - 1L
    iv_rows <- seq_len(n_iv)
    stopifnot(all(diff(alts$sy_row[iv_rows]) >= 0), all(ev$ev_start > n_iv), all(ev$ev_end <= nrow(alts)),
              all(ev$alt_row >= ev$ev_start & ev$alt_row <= ev$ev_end))
  } else {
    n_iv <- nrow(alts); iv_rows <- seq_len(n_iv)
    stopifnot(all(diff(alts$sy_row) >= 0))                     # sorted by state-year
  }
  stopifnot(setequal(unique(alts$sy_row[iv_rows]), seq_len(nrow(sy))))  # every state-year has an IV set
  stopifnot(all(sy$y1 %in% c(0L, 1L)))
  stopifnot(all(ev$sy_row %in% seq_len(nrow(sy))), all(ev$alt_row %in% seq_len(nrow(alts))))
  stopifnot(all(alts$sy_row[ev$alt_row] == ev$sy_row))         # chosen target belongs to the event's state-year
  stopifnot(all(sy$y1[ev$sy_row] == 1L))                        # every event sits in an initiating state-year
  n_orphan <- sum(sy$y1[setdiff(seq_len(nrow(sy)), ev$sy_row)] == 1L)
  ## y1 == 1 without a stage-2 event happens in real data when the only target lacks positions;
  ## the simulator never produces it, so it is an error unless explicitly allowed.
  if (n_orphan > 0 && !allow_orphan_initiations) stop(n_orphan, " state-years have y1 = 1 but no event")
  S_data <- max(c(sy$state_idx, alts$target_idx))
  if (is.null(S)) S <- S_data else stopifnot(S >= S_data)   # override keeps one index space across subsets
  cs <- alts[iv_rows, , drop = FALSE] |> mutate(row = row_number()) |> group_by(sy_row) |>
    summarise(s = min(row), e = max(row), n = n(), .groups = "drop") |> arrange(sy_row)
  stopifnot(all(cs$e - cs$s + 1 == cs$n))
  if (is.null(n_full)) n_full <- cs$n
  stopifnot(length(n_full) == nrow(sy), all(n_full >= cs$n))
  if (sampled == 1) { ev_start <- as.integer(ev$ev_start); ev_end <- as.integer(ev$ev_end) }
  else { ev_start <- as.integer(cs$s[ev$sy_row]); ev_end <- as.integer(cs$e[ev$sy_row]) }
  stopifnot(all(ev$alt_row >= ev_start & ev$alt_row <= ev_end))  # chosen target inside the event's own set
  if (is.null(x_ref)) x_ref <- colMeans(as.matrix(sy[, X_cols, drop = FALSE]))
  W <- as.matrix(sy[, W_cols, drop = FALSE]); X <- as.matrix(sy[, X_cols, drop = FALSE])
  Z <- as.matrix(alts[, Z_cols, drop = FALSE])
  assert_no_na(sy, c("y1", "state_idx", W_cols, X_cols), "stage-1 panel")
  assert_no_na(alts, c("sy_row", "target_idx", "dM", "dR", Z_cols), "alternatives")
  dat <- list(S = S, N1 = nrow(sy), N2 = nrow(alts), E = nrow(ev),
              P1 = ncol(W), P = ncol(X), Pz = ncol(Z),
              y1 = as.integer(sy$y1), sy_state = as.integer(sy$state_idx),
              W = W, X = X,
              cs_start = as.integer(cs$s), cs_end = as.integer(cs$e), n_full = as.integer(n_full),
              alt_target = as.integer(alts$target_idx), x = cbind(alts$dM, alts$dR), Z = Z,
              ev_sy = as.integer(ev$sy_row), ev_row = as.integer(ev$alt_row), ev_start = ev_start, ev_end = ev_end,
              joint = as.integer(joint), sampled = sampled,
              x_ref = as.numeric(x_ref), prior_beta_sd = prior_beta_sd, prior_mu0_sd = prior_mu0_sd,
              prior_only = as.integer(prior_only), D_mode = as.integer(D_mode))
  ## reduce_sum_static slice size over state-years: NULL -> one slice (exact serial sum, threads = 1)
  dat$grainsize <- as.integer(if (is.null(grainsize)) dat$N1 else max(1L, grainsize))
  stopifnot(dat$D_mode %in% 1:3, dat$grainsize >= 1)
  dimnames(dat$W) <- NULL; dimnames(dat$X) <- NULL; dimnames(dat$Z) <- NULL
  if (!is.null(K)) {
    stopifnot(K >= 1, K <= S - 1)
    if (is.null(anchor)) anchor <- choose_anchors(sy, ev, K, S)
    stopifnot(length(anchor) == K, !anyDuplicated(anchor), all(anchor %in% seq_len(S)))
    dat$K <- as.integer(K); dat$anchor <- array(as.integer(anchor), dim = K)   # array so K = 1 is a 1-d array for Stan
  }
  dat
}

## Anchor states for AME identification: the K states with the most initiation events.
## Ties broken by state index (deterministic). NOT string-sorted ccodes (audit F6).
choose_anchors <- function(sy, ev, K, S = max(sy$state_idx)) {
  cnt <- tabulate(sy$state_idx[ev$sy_row], nbins = S)
  ord <- order(-cnt, seq_along(cnt))
  ord[seq_len(K)]
}

## ---- sampled alternatives --------------------------------------------------------------
## Two outcome-separated samples (fixes the Codex review items 3 and 4):
##   IV rows    : for EVERY state-year, m alternatives drawn uniformly without replacement from
##                its full set, ignoring outcomes (chosen targets are NOT forced in, so the IV
##                estimate has the same distribution for initiating and non-initiating rows).
##   event rows : EACH event gets its own set = its chosen target + m alternatives drawn
##                uniformly from the remaining alternatives, independently per event (so two
##                events of one state-year never share an outcome-dependent set).
## Returns list(alts = IV rows then event rows (sorted), ev with alt_row / ev_start / ev_end into
## the event block, n_full per state-year, n_iv = number of IV rows).
## Full choice sets (m = 0 / sampled = 0) are the main specification; this is an option only.
sample_alternatives <- function(alts, ev, m, seed = V2_SEED) {
  stopifnot(m >= 1, all(diff(alts$sy_row) >= 0))
  set.seed(seed)
  alts <- alts |> mutate(.row = row_number())
  ## (i) IV rows: uniform m per state-year, outcome-blind
  iv <- alts |> group_by(sy_row) |> mutate(.u = runif(n()), .rank = rank(.u, ties.method = "first")) |>
    filter(.rank <= m) |> ungroup() |> arrange(sy_row, .row) |> select(-.u, -.rank)
  n_full <- alts |> count(sy_row) |> arrange(sy_row) |> pull(n)
  stopifnot(length(n_full) == length(unique(alts$sy_row)))
  ## (ii) event rows: per event, chosen + m others from the same state-year's full set
  ev_blocks <- vector("list", nrow(ev))
  for (e in seq_len(nrow(ev))) {
    rows_n <- alts$.row[alts$sy_row == ev$sy_row[e]]
    others <- setdiff(rows_n, ev$alt_row[e])
    pick <- if (length(others) > m) sample(others, m) else others
    blk <- alts[c(ev$alt_row[e], sort(pick)), , drop = FALSE]
    blk$.event <- e
    ev_blocks[[e]] <- blk
  }
  evrows <- bind_rows(ev_blocks)
  n_iv <- nrow(iv)
  out_alts <- bind_rows(iv |> mutate(.event = 0L), evrows)
  ev2 <- ev
  ev2$ev_start <- n_iv + match(seq_len(nrow(ev)), evrows$.event)                 # first row of each block (the chosen target)
  ev2$ev_end   <- n_iv + nrow(evrows) + 1L - match(seq_len(nrow(ev)), rev(evrows$.event))
  ev2$alt_row  <- ev2$ev_start                                                     # chosen target is the block's first row
  stopifnot(all(out_alts$sy_row[ev2$alt_row] == ev2$sy_row),
            all(out_alts$.row[ev2$alt_row] == alts$.row[ev$alt_row]),            # same chosen target as before
            all(ev2$ev_end >= ev2$ev_start), all(ev2$ev_start > n_iv))
  list(alts = out_alts |> select(-.row, -.event), ev = ev2, n_full = n_full, n_iv = n_iv)
}

## ---- brute-force stage-2 likelihood in R (used by tests and the Fixed_param check) -----
## pars: list(mu0 (D), beta (P x D), a_i (S x D), b_j (S), gamma (Pz), sigma_x (D), rho, U (S x K) or NULL, V)
## D = 2 when dat$D_mode == 3; D = 1 uses column dat$D_mode of x (1 = dM, 2 = dR).
## mu = mu0 + X beta + a_i; mu0 is required (no silent zero default) so a caller that forgets it fails here.
stage2_loglik_R <- function(dat, pars) {
  D <- if (dat$D_mode == 3L) 2L else 1L
  xcols <- if (D == 2L) 1:2 else dat$D_mode
  if (is.null(pars$mu0)) stop("stage2_loglik_R: pars$mu0 (length D intercept of mu) is missing")
  stopifnot(length(pars$mu0) == D, all(is.finite(pars$mu0)))
  Mu <- rep(1, dat$N1) %o% as.numeric(pars$mu0) +
    dat$X %*% matrix(pars$beta, ncol = D) + matrix(pars$a_i, ncol = D)[dat$sy_state, , drop = FALSE]
  zg <- if (dat$Pz > 0) as.numeric(dat$Z %*% pars$gamma) else rep(0, dat$N2)
  util <- numeric(dat$N2); lse <- numeric(dat$N1)
  Sinv <- if (D == 2L) solve(matrix(c(pars$sigma_x[1]^2, pars$rho * prod(pars$sigma_x),
                                      pars$rho * prod(pars$sigma_x), pars$sigma_x[2]^2), 2))
          else matrix(1 / pars$sigma_x[1]^2, 1, 1)
  one_util <- function(r, n) {
    d <- dat$x[r, xcols] - Mu[n, ]
    u <- -0.5 * as.numeric(t(d) %*% Sinv %*% d) + pars$b_j[dat$alt_target[r]] + zg[r]
    if (!is.null(pars$U)) u <- u + sum(pars$U[dat$sy_state[n], ] * pars$V[dat$alt_target[r], ])
    u
  }
  for (n in seq_len(dat$N1)) {
    rows <- dat$cs_start[n]:dat$cs_end[n]
    for (r in rows) util[r] <- one_util(r, n)
    lse[n] <- log(sum(exp(util[rows])))
  }
  loglik <- numeric(dat$E)
  for (e in seq_len(dat$E)) {
    rows <- dat$ev_start[e]:dat$ev_end[e]
    if (dat$sampled == 1) for (r in rows) util[r] <- one_util(r, dat$ev_sy[e])
    loglik[e] <- util[dat$ev_row[e]] - log(sum(exp(util[rows])))
  }
  list(loglik = loglik, lse = lse, util = util)
}

## ---- Procrustes alignment of U/V draws across chains ---------------------------------
## U_draws: array [iter, chain, S, K]. Rotates each chain's draws (U R, V R) to align the
## chain-mean U with the reference chain's mean U via orthogonal Procrustes (R = W V' from
## the SVD of U_c' U_ref). u_i'v_j is invariant to R, so the likelihood is unchanged.
procrustes_align <- function(U_draws, V_draws, ref_chain = 1) {
  stopifnot(length(dim(U_draws)) == 4, all(dim(U_draws) == dim(V_draws)))
  n_iter <- dim(U_draws)[1]; n_chain <- dim(U_draws)[2]; K <- dim(U_draws)[4]
  U_ref <- apply(U_draws[, ref_chain, , , drop = FALSE], c(3, 4), mean)
  rotations <- vector("list", n_chain)
  for (ch in seq_len(n_chain)) {
    U_c <- apply(U_draws[, ch, , , drop = FALSE], c(3, 4), mean)
    sv <- svd(t(U_c) %*% U_ref)
    R <- sv$u %*% t(sv$v)                       # K x K orthogonal
    rotations[[ch]] <- R
    for (it in seq_len(n_iter)) {
      U_draws[it, ch, , ] <- matrix(U_draws[it, ch, , ], ncol = K) %*% R
      V_draws[it, ch, , ] <- matrix(V_draws[it, ch, , ], ncol = K) %*% R
    }
  }
  list(U = U_draws, V = V_draws, rotations = rotations)
}

## grainsize for reduce_sum_static: with 1 thread one slice (serial); otherwise about 4 slices per thread
choose_grainsize <- function(N1, threads) if (threads <= 1) N1 else max(1L, as.integer(ceiling(N1 / (4 * threads))))
