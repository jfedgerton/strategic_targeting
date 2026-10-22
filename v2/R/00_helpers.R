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
assert_standardized <- function(df, col, by = "year", tol = 1e-6, what = "data") {
  chk <- df |> group_by(across(all_of(by))) |>
    summarise(m = mean(.data[[col]]), s = sd(.data[[col]]), n = n(), .groups = "drop") |>
    filter(n > 1)
  if (any(abs(chk$m) > tol) || any(abs(chk$s - 1) > 1e-3))
    stop(sprintf("%s: %s is not standardized within %s (max |mean| = %.2e, sd range [%.4f, %.4f])",
                 what, col, by, max(abs(chk$m)), min(chk$s), max(chk$s)))
  invisible(TRUE)
}

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
build_stan_data <- function(sy, alts, ev, W_cols, X_cols, Z_cols, joint, sampled = 0L,
                            n_full = NULL, x_ref = NULL, K = NULL, anchor = NULL,
                            prior_beta_sd = 1, allow_orphan_initiations = FALSE, S = NULL) {
  stopifnot(is.data.frame(sy), is.data.frame(alts), is.data.frame(ev))
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
              x_ref = as.numeric(x_ref), prior_beta_sd = prior_beta_sd)
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
## pars: list(beta (P x 2), a_i (S x 2), b_j (S), gamma (Pz), sigma_x (2), rho, U (S x K) or NULL, V)
stage2_loglik_R <- function(dat, pars) {
  Mu <- dat$X %*% pars$beta + pars$a_i[dat$sy_state, , drop = FALSE]
  zg <- if (dat$Pz > 0) as.numeric(dat$Z %*% pars$gamma) else rep(0, dat$N2)
  util <- numeric(dat$N2); lse <- numeric(dat$N1)
  Sinv <- solve(matrix(c(pars$sigma_x[1]^2, pars$rho * prod(pars$sigma_x),
                         pars$rho * prod(pars$sigma_x), pars$sigma_x[2]^2), 2))
  one_util <- function(r, n) {
    d <- dat$x[r, ] - Mu[n, ]
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
