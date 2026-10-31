## ============================================================================
## helpers_model.R -- model-side helper functions of the v2 pipeline (functions only, no side effects)
## Sourced after spec.R and helpers_data.R: fit provenance (assert_fit_current), parameter labels, fit safety (pilots, stale
## fits, compiled-model provenance), the Stan data assembly build_stan_data(), AME anchors, sampled alternatives (archived
## option, kept for the tests), the brute-force stage-2 likelihood, Procrustes alignment and the reduce_sum grainsize.
## ============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

## ---- fit provenance (audit 2026-10-23): a fit must have been made from the CURRENT data files ----------------
## obj$manifest holds the manifest rows (file, sha256) of the data the fit used; compare with the current manifest.
assert_fit_current <- function(obj, manifest_path) {
  stopifnot(!is.null(obj$manifest), nrow(obj$manifest) > 0, file.exists(manifest_path))
  cur <- read.csv(manifest_path, stringsAsFactors = FALSE)
  m <- merge(obj$manifest[, c("file", "sha256")], cur[, c("file", "sha256")], by = "file", all.x = TRUE, suffixes = c("_fit", "_now"))
  bad <- m$file[is.na(m$sha256_now) | m$sha256_fit != m$sha256_now]
  if (length(bad)) stop(sprintf("fit %s was made from data that no longer match the current build (%s); refit or rebuild", obj$cell, paste(bad, collapse = ", ")))
  invisible(TRUE)
}

## ---- parameter labels for the summary CSVs (review 2026-10-25, B2) ----------------------------------------------------------
## axis_labels(D_mode): the targeting axes of a cell. label_parameters(pars, W_cols, X_cols, Z_cols, D_mode, K) returns a data.frame
## (parameter, axis, covariate) for rstan summary row names: alpha[k] -> W_cols[k]; beta[p,d] -> X_cols[p] on axis d; gamma[k] ->
## Z_cols[k]; mu0[d], sigma_x[d], mu_state[s,d], mu_event[e,d] -> axis d; sd_re[k], Omega / R_Om indices -> the random-effect order
## c, a_<axis1>[, a_<axis2>]; Sigma[i,j] -> the axis pair; u_diag[k] -> latent dimension k. Anything else gets NA labels.
axis_labels <- function(D_mode) { stopifnot(D_mode %in% 1:3); if (D_mode == 3L) c("material", "relational") else c("material", "relational")[D_mode] }
label_parameters <- function(pars, W_cols, X_cols, Z_cols, D_mode, groups = STAGE2_GROUPS) {
  ax <- axis_labels(D_mode); re <- c("c", paste0("a_", ax))
  base <- sub("\\[.*$", "", pars)
  has <- grepl("^[A-Za-z_0-9]+\\[[0-9]+(,[0-9]+)?\\]$", pars)
  i1 <- ifelse(has, suppressWarnings(as.integer(sub("^[^\\[]+\\[([0-9]+).*$", "\\1", pars))), NA_integer_)
  i2 <- ifelse(grepl(",", pars, fixed = TRUE) & has, suppressWarnings(as.integer(sub("^[^,]+,([0-9]+)\\]$", "\\1", pars))), NA_integer_)
  pick <- function(v, i) ifelse(is.na(i) | i > length(v), NA_character_, v[pmin(pmax(i, 1L), length(v))])
  axis <- rep(NA_character_, length(pars)); cov <- rep(NA_character_, length(pars))
  w <- base == "alpha"; cov[w] <- pick(W_cols, i1[w])
  w <- base == "beta"; cov[w] <- pick(X_cols, i1[w]); axis[w] <- pick(ax, i2[w])
  w <- base == "gamma"; cov[w] <- pick(Z_cols, i1[w])
  w <- base %in% c("mu0", "sigma_x"); axis[w] <- pick(ax, i1[w])
  w <- base %in% c("mu_state", "mu_event"); axis[w] <- pick(ax, i2[w])
  w <- base == "sd_re"; cov[w] <- pick(re, i1[w])
  w <- base %in% c("Omega", "R_Om"); cov[w] <- paste(pick(re, i1[w]), pick(re, i2[w]), sep = " x ")
  w <- base == "Sigma"; axis[w] <- paste(pick(ax, i1[w]), pick(ax, i2[w]), sep = " x ")
  w <- base == "u_diag"; cov[w] <- paste0("latent_dim_", i1[w])
  ## M1c / M1cr (stage b; G groups as data since 2026-10-31): delta[g, d] -> group g (0-based, so label groups[g + 1]), axis d;
  ## mu0_g[k, d] and sigma_g[k, d] -> group label groups[k], axis d; rho_g[k] (M1cr) -> group label groups[k]; Sigma_g[k, i, j] is not summarised
  w <- base == "delta"; axis[w] <- pick(ax, i2[w]); cov[w] <- pick(groups, i1[w] + 1L)
  w <- base %in% c("mu0_g", "sigma_g"); axis[w] <- pick(ax, i2[w]); cov[w] <- pick(groups, i1[w])
  w <- base == "rho_g"; axis[w] <- paste(ax, collapse = " x "); cov[w] <- pick(groups, i1[w])
  w <- base == "rho_x"; axis[w] <- paste(ax, collapse = " x ")
  w <- base == "tau"; cov[w] <- "IV"
  w <- base == "alpha0"; cov[w] <- "intercept"
  out <- data.frame(parameter = pars, axis = axis, covariate = cov, stringsAsFactors = FALSE)
  stopifnot(nrow(out) == length(pars), !anyNA(out$covariate[base == "alpha"]), !anyNA(out$covariate[base == "beta"]), !anyNA(out$covariate[base == "gamma"]))
  out
}

## ---- fit safety (review 2026-10-25, C1 / C2): pilots and stale fits are refused by the production tools -----------------
## A fit object records its sampling configuration in obj$args (pilot flag, chains, iter, warmup) since 2026-10-25; fits saved
## before that are recognised as pilots by the pilot settings themselves (1 chain x 60 iterations). obj$stale is non-NULL when the
## manifest changed while the fit sampled (saved as fit_<cell>_STALE.rds); such a fit passes only with allow_stale = TRUE.
assert_production_fit <- function(obj, allow_stale = FALSE, what = "fit") {
  a <- obj$args
  if (is.null(a)) stop(sprintf("%s %s has no saved sampling configuration (args); refusing to use it", what, obj$cell))
  is_pilot <- isTRUE(a$pilot) || (is.null(a$pilot) && !is.null(a$iter) && !is.null(a$chains) && a$iter <= 60 && a$chains == 1)
  if (is_pilot) stop(sprintf("%s %s is a PILOT fit (chains %s, iter %s); production extraction refuses pilots", what, obj$cell, a$chains, a$iter))
  if (!is.null(obj$stale)) {
    if (!allow_stale) stop(sprintf("%s %s is STALE (the data manifest changed while it sampled: %s); pass --allow_stale to use it anyway", what, obj$cell, paste(obj$stale$changed_files, collapse = ", ")))
    logf("WARNING: using STALE fit %s (manifest changed while sampling: %s) because --allow_stale was given", obj$cell, paste(obj$stale$changed_files, collapse = ", "))
  }
  ## the stanfit must hold every requested chain with its full set of post-warmup draws (review 2026-10-26, P1; see assert_fit_chains)
  if (!is.null(obj$fit)) assert_fit_chains(obj$fit, chains = a$chains, iter = a$iter, warmup = a$warmup, what = what, cell = obj$cell)
  invisible(TRUE)
}
## ---- actual chains versus the requested configuration (review 2026-10-26, P1) -------------------------------------------------------
## With parallel sampling rstan drops a chain that fails (e.g. a rejected init) and returns the surviving chains with fit@mode == 0, so a
## job could save a "production" fit of 3 chains, record the 4 it requested and exit 0. fit_chain_summary(fit) reads the ACTUAL chains and
## post-warmup draws per chain from the stanfit (fit@sim$n_save - fit@sim$warmup2); assert_fit_chains() requires them to equal the requested
## chains and (iter - warmup) / thin, for each chain. 02_fit_model.R calls it before saving; assert_production_fit() repeats it downstream.
fit_chain_summary <- function(fit) {
  stopifnot(inherits(fit, "stanfit"))
  if (fit@mode != 0L) return(list(chains = 0L, draws = integer(0), thin = NA_integer_))
  n_save <- as.integer(fit@sim$n_save); warm <- as.integer(fit@sim$warmup2)
  stopifnot(length(n_save) == fit@sim$chains, length(warm) == length(n_save))
  list(chains = as.integer(fit@sim$chains), draws = n_save - warm, thin = as.integer(fit@sim$thin))
}
assert_fit_chains <- function(fit, chains, iter, warmup, what = "fit", cell = "") {
  if (is.null(chains) || is.null(iter) || is.null(warmup)) stop(sprintf("%s %s records no chains / iter / warmup; cannot check the sampled chains", what, cell))
  cs <- fit_chain_summary(fit)
  if (cs$chains != chains) stop(sprintf("%s %s holds %d chain(s) but %d were requested: a chain failed (rstan drops a failed parallel chain and returns the rest with mode 0); not a production fit",
                                        what, cell, cs$chains, as.integer(chains)))
  expect_draws <- as.integer(ceiling((iter - warmup) / cs$thin))
  if (length(cs$draws) != chains || any(cs$draws != expect_draws))
    stop(sprintf("%s %s: post-warmup draws per chain %s, expected %d from each of %d chains (iter %d, warmup %d, thin %d)", what, cell,
                 paste(cs$draws, collapse = "/"), expect_draws, as.integer(chains), as.integer(iter), as.integer(warmup), cs$thin))
  invisible(TRUE)
}
## read_fit_or_null(path): readRDS wrapped so an unreadable / truncated fit is reported and treated as missing (C2)
read_fit_or_null <- function(path) {
  obj <- tryCatch(readRDS(path), error = function(e) { logf("UNREADABLE fit %s (%s); treated as missing", basename(path), conditionMessage(e)); NULL })
  if (!is.null(obj) && (is.null(obj$fit) || is.null(obj$dat))) { logf("fit file %s lacks fit / dat; treated as missing", basename(path)); return(NULL) }
  obj
}
## ---- fit / cell consistency (review 2026-10-25, C4; moved here from measurement_tests.R 2026-10-26, audit fix 1) --------------
## assert_fit_matches_cell(obj, dat): the saved fit must score exactly the assembled cell (E, N1, S, Z / W / X columns, D_mode,
## event identities = (state-year, chosen target) pairs, event weights). 02_fit_model.R uses it for the M1 fit an M2 fit starts from.
assert_fit_matches_cell <- function(obj, dat, W_cols, X_cols, Z_cols, what = "fit") {
  d <- obj$dat
  if (is.null(d) || is.null(d$E)) stop(sprintf("%s %s has no saved Stan data (dat); refusing to use it", what, obj$cell))
  if (d$E != dat$E) stop(sprintf("%s %s has E = %d, the assembled cell %d", what, obj$cell, d$E, dat$E))
  if (d$N1 != dat$N1 || d$S != dat$S) stop(sprintf("%s %s: N1 / S %d / %d differ from the cell's %d / %d", what, obj$cell, d$N1, d$S, dat$N1, dat$S))
  if (!identical(obj$Z_cols, Z_cols)) stop(sprintf("%s %s: Z columns (%s) differ from the cell's (%s)", what, obj$cell, paste(obj$Z_cols, collapse = ","), paste(Z_cols, collapse = ",")))
  if (!identical(obj$W_cols, W_cols) || !identical(obj$X_cols, X_cols)) stop(sprintf("%s %s: W / X columns differ from the cell's", what, obj$cell))
  if (!identical(as.integer(d$ev_sy), as.integer(dat$ev_sy)) || !identical(as.integer(d$alt_target[d$ev_row]), as.integer(dat$alt_target[dat$ev_row])))
    stop(sprintf("%s %s: the events (state-year, chosen target) differ from the assembled cell", what, obj$cell))
  if (!is.null(d$ev_w) && !is.null(dat$ev_w) && max(abs(d$ev_w - dat$ev_w)) > 1e-10) stop(sprintf("%s %s: the event weights differ from the assembled cell", what, obj$cell))
  invisible(TRUE)
}

## ---- fit / program consistency (audit fix 5, 2026-10-26): a fit must come from the CURRENT Stan program -----------------------------
## 02_fit_model.R saves the sha256 of the .stan source it sampled with in obj$args$stan_source_sha256 (C3). The post-processing
## tools (03_extract.R, 04_compare.R, 05_robustness.R perm_error, 06_context_compare.R) compare it with the sidecar hash
## stan/build/<family>.sha256 written by compile_models.R AND with the hash of the checked-out source stan/<family>.stan (review
## 2026-10-26, P2: the sidecar alone does not notice an edited source), so fits of an earlier program can never be mixed with fits of
## the current one in the same output directory. A fit without the recorded hash (saved before C3) is refused as well.
assert_fit_stan_current <- function(obj, root = v2_root(), what = "fit") {
  sha <- obj$args$stan_source_sha256
  if (is.null(sha) || !is.character(sha) || length(sha) != 1 || !grepl("^[0-9a-f]{64}$", sha))
    stop(sprintf("%s %s records no Stan source hash (args$stan_source_sha256; saved before 2026-10-25 C3); refit it with the current program", what, obj$cell))
  if (is.null(obj$args$model)) stop(sprintf("%s %s records no model code (args$model)", what, obj$cell))
  ## m1_latent_additive | m2_ame | m1_context (M1c, stage b) | m1_latent_prior (an M1 fit with a prior switch off its default, 2026-10-30).
  ## The program is derived from the model code AND the saved switches (a fit saved before the prior switches existed has no such fields and
  ## is the main program's); a fit that recorded which program it sampled (args$stan_program, since 2026-10-27) must agree, so a fit of one
  ## program can never be checked against the other's sidecar
  family <- stan_program(obj$args$model, obj$spec)
  if (!is.null(obj$args$stan_program) && !identical(obj$args$stan_program, family))
    stop(sprintf("%s %s records the program %s but its model code / switches (%s) imply %s; the fit is not what its name says", what, obj$cell, obj$args$stan_program, obj$args$model, family))
  sidecar <- file.path(root, "stan", "build", paste0(family, ".sha256"))
  if (!file.exists(sidecar)) stop("compiled-model sidecar hash missing: ", sidecar, "; rerun v2/R/compile_models.R")
  side <- trimws(readLines(sidecar, warn = FALSE)[1])
  src <- file.path(root, "stan", paste0(family, ".stan"))
  if (!file.exists(src)) stop("Stan source missing: ", src)
  src_sha <- stan_source_sha256(src)
  if (!identical(side, src_sha))
    stop(sprintf("sidecar hash %s does not match the source hash %s of %s: the checked-out source changed after compilation; rerun v2/R/compile_models.R and refit", side, src_sha, basename(src)))
  if (!identical(side, sha))
    stop(sprintf("%s %s was sampled with %s source %s, the compiled program is %s (stan/build/%s.sha256): fits of different programs in one output directory; deploy into a fresh clone or empty v2/output and refit",
                 what, obj$cell, family, substr(sha, 1, 12), substr(side, 1, 12), family))
  invisible(TRUE)
}

## ---- compiled-model provenance (review 2026-10-25, C3) --------------------------------------------------------------------
## compile_models.R writes <model>.sha256 next to <model>.rds with the sha256 of the .stan SOURCE it compiled. Before sampling,
## assert_compiled_current() checks that the compiled program's model_code is the checked-out source (text compared after
## normalising line endings and trailing whitespace) and that the sidecar hash equals the source hash, and returns both hashes
## so the fit can store them. This replaces the two-marker string check (mu0 / in_iv) of 2026-10-23 and 2026-10-24.
stan_source_sha256 <- function(path) digest::digest(normalise_stan_text(readChar(path, file.info(path)$size, useBytes = TRUE)), algo = "sha256", serialize = FALSE)
normalise_stan_text <- function(txt) { txt <- gsub("\r\n", "\n", txt); txt <- gsub("[ \t]+\n", "\n", txt); sub("[[:space:]]+$", "", txt) }
assert_compiled_current <- function(sm, src_path, sidecar_path = NULL) {
  stopifnot(inherits(sm, "stanmodel"), file.exists(src_path))
  src_txt <- normalise_stan_text(readChar(src_path, file.info(src_path)$size, useBytes = TRUE))
  code_txt <- normalise_stan_text(paste(sm@model_code, collapse = "\n"))
  if (!identical(src_txt, code_txt))
    stop(sprintf("compiled %s does not match the checked-out source %s (model_code differs): stale binary; rerun v2/R/compile_models.R", sm@model_name, basename(src_path)))
  src_sha <- digest::digest(src_txt, algo = "sha256", serialize = FALSE)
  code_sha <- digest::digest(code_txt, algo = "sha256", serialize = FALSE)
  if (!is.null(sidecar_path)) {
    if (!file.exists(sidecar_path)) stop("compiled-model sidecar hash missing: ", sidecar_path, "; rerun v2/R/compile_models.R")
    side <- trimws(readLines(sidecar_path, warn = FALSE)[1])
    if (!identical(side, src_sha)) stop(sprintf("sidecar hash %s does not match the source hash %s of %s; rerun v2/R/compile_models.R", side, src_sha, basename(src_path)))
  }
  list(stan_source_sha256 = src_sha, compiled_model_code_sha256 = code_sha)
}

## ---- Stan data assembly ----------------------------------------------------------------
## Builds the Stan data list from three data frames:
##   sy:   state-years, one row per (state_idx, year); columns y1, state_idx, W cols, X cols
##   alts: alternative rows, one row per (sy_row, target_idx); columns sy_row, target_idx, dM, dR, Z cols
##   ev:   events, columns sy_row, alt_row (index into alts after ordering by sy_row)
## alts must be sorted by sy_row; every sy row must own >= 1 alternative.
## Choice-set flags (Codex review P1-1, PI choice "A'", 2026-10-24). alts may carry two 0/1 columns, in_iv (the row is
## in the state-year's stage-1 IV set, defined at January of t) and in_s2 (the row is in the stage-2 choice set, defined
## at the state-year's reference date). A state-year's row block is the UNION of both sets; the Stan programs take
## IV[n] over the in_iv rows and the stage-2 softmax over the in_s2 rows. Absent columns mean both sets are the whole
## block (the simulator, the tests, every pre-review caller). Every state-year needs >= 1 in_iv row; every event's
## chosen row needs in_s2 = 1 and its stage-2 set >= 2 rows (a chosen-only set has no information); with sampled = 1
## the sets are resolved by sample_alternatives() and every row carries in_iv = in_s2 = 1.
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
  if (!"in_iv" %in% names(alts)) alts$in_iv <- 1L
  if (!"in_s2" %in% names(alts)) alts$in_s2 <- 1L
  alts$in_iv <- as.integer(alts$in_iv); alts$in_s2 <- as.integer(alts$in_s2)
  stopifnot(!anyNA(alts$in_iv), !anyNA(alts$in_s2), all(alts$in_iv %in% 0:1), all(alts$in_s2 %in% 0:1), all(alts$in_iv == 1L | alts$in_s2 == 1L))
  if (sampled == 1) {
    ## alts = IV rows (first n_iv, sorted by sy_row) followed by per-event blocks; ev carries ev_start/ev_end.
    ## sample_alternatives() draws the IV rows from the in_iv set and the event blocks from the in_s2 set, so the flags
    ## are resolved: every remaining row is in its own set.
    stopifnot(all(c("ev_start", "ev_end") %in% names(ev)), !is.null(n_full))
    if (!all(alts$in_iv == 1L & alts$in_s2 == 1L)) stop("build_stan_data: with sampled = 1 every row must have in_iv = in_s2 = 1 (use sample_alternatives())")
    n_iv <- min(ev$ev_start) - 1L
    iv_rows <- seq_len(n_iv)
    stopifnot(all(diff(alts$sy_row[iv_rows]) >= 0), all(ev$ev_start > n_iv), all(ev$ev_end <= nrow(alts)),
              all(ev$alt_row >= ev$ev_start & ev$alt_row <= ev$ev_end))
    if (any(ev$ev_end - ev$ev_start + 1L < 2L)) stop("build_stan_data: an event's sampled set holds only the chosen target")
  } else {
    n_iv <- nrow(alts); iv_rows <- seq_len(n_iv)
    stopifnot(all(diff(alts$sy_row) >= 0))                     # sorted by state-year
  }
  stopifnot(setequal(unique(alts$sy_row[iv_rows]), seq_len(nrow(sy))))  # every state-year has a row block
  ## every state-year has a NON-EMPTY IV set (rows with in_iv = 1)
  n_iv_sy <- tabulate(alts$sy_row[iv_rows][alts$in_iv[iv_rows] == 1L], nbins = nrow(sy))
  if (any(n_iv_sy == 0L)) stop(sprintf("build_stan_data: %d state-years have an EMPTY stage-1 IV set (rows exist but none has in_iv = 1)", sum(n_iv_sy == 0L)))
  stopifnot(all(sy$y1 %in% c(0L, 1L)))
  stopifnot(all(ev$sy_row %in% seq_len(nrow(sy))), all(ev$alt_row %in% seq_len(nrow(alts))))
  stopifnot(all(alts$sy_row[ev$alt_row] == ev$sy_row))         # chosen target belongs to the event's state-year
  stopifnot(all(sy$y1[ev$sy_row] == 1L))                        # every event sits in an initiating state-year
  ## the chosen row is in the stage-2 set, and that set is not the chosen row alone
  if (any(alts$in_s2[ev$alt_row] != 1L)) stop(sprintf("build_stan_data: %d events have a chosen row outside the stage-2 set (in_s2 = 0)", sum(alts$in_s2[ev$alt_row] != 1L)))
  if (sampled == 0) {
    n_s2_sy <- tabulate(alts$sy_row[alts$in_s2 == 1L], nbins = nrow(sy))
    bad_s2 <- unique(ev$sy_row[n_s2_sy[ev$sy_row] < 2L])
    if (length(bad_s2)) stop(sprintf("build_stan_data: %d initiating state-years have an EMPTY or CHOSEN-ONLY stage-2 set (< 2 rows with in_s2 = 1)", length(bad_s2)))
  }
  n_orphan <- sum(sy$y1[setdiff(seq_len(nrow(sy)), ev$sy_row)] == 1L)
  ## y1 == 1 without a stage-2 event happens in real data when the only target lacks positions;
  ## the simulator never produces it, so it is an error unless explicitly allowed.
  if (n_orphan > 0 && !allow_orphan_initiations) stop(n_orphan, " state-years have y1 = 1 but no event")
  S_data <- max(c(sy$state_idx, alts$target_idx))
  if (is.null(S)) S <- S_data else stopifnot(S >= S_data)   # override keeps one index space across subsets
  cs <- alts[iv_rows, , drop = FALSE] |> mutate(row = row_number()) |> group_by(sy_row) |>
    summarise(s = min(row), e = max(row), n = n(), n_iv = sum(in_iv), .groups = "drop") |> arrange(sy_row)
  stopifnot(all(cs$e - cs$s + 1 == cs$n), all(cs$n_iv >= 1L))
  if (is.null(n_full)) n_full <- cs$n_iv                        # full IV-set size (the sampling correction's numerator)
  stopifnot(length(n_full) == nrow(sy), all(n_full >= cs$n_iv))
  if (sampled == 1) { ev_start <- as.integer(ev$ev_start); ev_end <- as.integer(ev$ev_end) }
  else { ev_start <- as.integer(cs$s[ev$sy_row]); ev_end <- as.integer(cs$e[ev$sy_row]) }
  stopifnot(all(ev$alt_row >= ev_start & ev$alt_row <= ev_end))  # chosen target inside the event's own block
  if (is.null(x_ref)) x_ref <- colMeans(as.matrix(sy[, X_cols, drop = FALSE]))
  ## stage-2 event weights (PI decision 2026-10-25, D8): ev$w when present (assemble_cell(), 1 / the initiator's targets in the
  ## same dispute and year), else 1 for every event (the simulator and older callers). They multiply the stage-2 log-lik terms
  ## in the Stan programs; log_lik_stage2 stays unweighted.
  ev_w <- if ("w" %in% names(ev)) as.numeric(ev$w) else rep(1, nrow(ev))
  stopifnot(length(ev_w) == nrow(ev), all(is.finite(ev_w)), all(ev_w > 0), all(ev_w <= 1))
  W <- as.matrix(sy[, W_cols, drop = FALSE]); X <- as.matrix(sy[, X_cols, drop = FALSE])
  if (length(X_cols) == 0L) X <- matrix(numeric(0), nrow = nrow(sy), ncol = 0L)   # no X (conventional comparison, 2026-10-31): a numeric N1 x 0 matrix, not the logical one as.matrix() returns for an empty selection
  Z <- as.matrix(alts[, Z_cols, drop = FALSE])
  assert_no_na(sy, c("y1", "state_idx", W_cols, X_cols), "stage-1 panel")
  assert_no_na(alts, c("sy_row", "target_idx", "dM", "dR", Z_cols), "alternatives")
  dat <- list(S = S, N1 = nrow(sy), N2 = nrow(alts), E = nrow(ev),
              P1 = ncol(W), P = ncol(X), Pz = ncol(Z),
              y1 = as.integer(sy$y1), sy_state = as.integer(sy$state_idx),
              W = W, X = X,
              cs_start = as.integer(cs$s), cs_end = as.integer(cs$e), n_full = as.integer(n_full),
              alt_target = as.integer(alts$target_idx), x = cbind(alts$dM, alts$dR), Z = Z,
              in_iv = as.integer(alts$in_iv), in_s2 = as.integer(alts$in_s2),
              ev_sy = as.integer(ev$sy_row), ev_row = as.integer(ev$alt_row), ev_start = ev_start, ev_end = ev_end, ev_w = ev_w,
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
## Choice-set flags (Codex review P1-1): the IV rows are drawn from the state-year's in_iv rows (n_full counts them),
## the event blocks from its in_s2 rows; absent flags mean the whole block. The output rows carry in_iv = in_s2 = 1
## because every row is then in the set it was drawn for.
## Full choice sets (m = 0 / sampled = 0) are the main specification; this is an option only.
sample_alternatives <- function(alts, ev, m, seed = V2_SEED) {
  stopifnot(m >= 1, all(diff(alts$sy_row) >= 0))
  if (!"in_iv" %in% names(alts)) alts$in_iv <- 1L
  if (!"in_s2" %in% names(alts)) alts$in_s2 <- 1L
  stopifnot(all(alts$in_iv %in% 0:1), all(alts$in_s2 %in% 0:1), all(alts$in_iv == 1L | alts$in_s2 == 1L))
  set.seed(seed)
  alts <- alts |> mutate(.row = row_number())
  ## (i) IV rows: uniform m per state-year from the IV set, outcome-blind
  iv <- alts |> filter(in_iv == 1L) |> group_by(sy_row) |> mutate(.u = runif(n()), .rank = rank(.u, ties.method = "first")) |>
    filter(.rank <= m) |> ungroup() |> arrange(sy_row, .row) |> select(-.u, -.rank)
  n_full <- alts |> filter(in_iv == 1L) |> count(sy_row) |> arrange(sy_row) |> pull(n)
  stopifnot(length(n_full) == length(unique(alts$sy_row)), setequal(unique(iv$sy_row), unique(alts$sy_row)))
  ## (ii) event rows: per event, chosen + m others from the same state-year's stage-2 set
  if (any(alts$in_s2[ev$alt_row] != 1L)) stop("sample_alternatives: a chosen row is not in the stage-2 set (in_s2 = 0)")
  ev_blocks <- vector("list", nrow(ev))
  for (e in seq_len(nrow(ev))) {
    rows_n <- alts$.row[alts$sy_row == ev$sy_row[e] & alts$in_s2 == 1L]
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
  out_alts$in_iv <- 1L; out_alts$in_s2 <- 1L                                       # resolved: every row is in the set it was drawn for
  list(alts = out_alts |> select(-.row, -.event), ev = ev2, n_full = n_full, n_iv = n_iv)
}

## ---- brute-force stage-2 likelihood in R (used by tests and the Fixed_param check) -----
## pars: list(mu0 (D), beta (P x D), a_i (S x D), b_j (S), gamma (Pz), sigma_x (D), rho, U (S x K) or NULL, V)
## D = 2 when dat$D_mode == 3; D = 1 uses column dat$D_mode of x (1 = dM, 2 = dR).
## mu = mu0 + X beta + a_i; mu0 is required (no silent zero default) so a caller that forgets it fails here.
## lse[n] is over the state-year's in_iv rows, the event softmax over its in_s2 rows (sampled = 0) or the event's own
## block (sampled = 1); dat$in_iv / dat$in_s2 are required (build_stan_data() always writes them).
stage2_loglik_R <- function(dat, pars) {
  D <- if (dat$D_mode == 3L) 2L else 1L
  xcols <- if (D == 2L) 1:2 else dat$D_mode
  if (is.null(pars$mu0)) stop("stage2_loglik_R: pars$mu0 (length D intercept of mu) is missing")
  stopifnot(length(pars$mu0) == D, all(is.finite(pars$mu0)))
  if (is.null(dat$in_iv) || is.null(dat$in_s2)) stop("stage2_loglik_R: dat$in_iv / dat$in_s2 are missing (build the data with build_stan_data())")
  stopifnot(length(dat$in_iv) == dat$N2, length(dat$in_s2) == dat$N2)
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
    rows_iv <- rows[dat$in_iv[rows] == 1L]
    stopifnot(length(rows_iv) >= 1)
    lse[n] <- log(sum(exp(util[rows_iv])))
  }
  loglik <- numeric(dat$E)
  for (e in seq_len(dat$E)) {
    rows <- dat$ev_start[e]:dat$ev_end[e]
    if (dat$sampled == 1) for (r in rows) util[r] <- one_util(r, dat$ev_sy[e])
    else rows <- rows[dat$in_s2[rows] == 1L]
    stopifnot(dat$ev_row[e] %in% rows)
    loglik[e] <- util[dat$ev_row[e]] - log(sum(exp(util[rows])))
  }
  ## stage-2 weights (D8): the composite log-likelihood weights every event term by dat$ev_w (1 when absent)
  ev_w <- if (is.null(dat$ev_w)) rep(1, dat$E) else dat$ev_w
  stopifnot(length(ev_w) == dat$E)
  list(loglik = loglik, lse = lse, util = util, ev_w = ev_w, loglik_weighted = loglik * ev_w, lp_stage2 = sum(loglik * ev_w))
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

## ---- M2 initial values from the matching M1 fit (PI decision 2026-10-25) -------------------------------------------------------------
## m1_inits_for_m2(obj_m1, dat_m2, chains, seed, jitter_sd, uv_range): one init list per chain for an M2 (AME) fit, built from the
## posterior MEANS of the M1 fit of the same cell (same sample, lag, switches; M2x from M1x). Every parameter the two programs share
## with identical dimensions is copied: alpha0, alpha, tau, mu0, beta, sigma_x, sd_re, sigma_b, gamma, the correlation parameters of
## the model's variant (rho_x and L_Om_joint in the b models, r_ca1 / q_ca2 in the a models, r_ca_1d in the 1-D models) and the
## random-effect standard-normals z_re / z_b. The M1 fits do not save z_re / z_b / L_Om_joint (excluded at sampling), so they are
## reconstructed from the saved transformed quantities: re = (c_i, a_i) = (diag(sd_re) L_Om z_re)' gives z_re = (diag(sd_re) L_Om)^-1 re',
## with L_Om the lower Cholesky factor of the posterior-mean correlation R_Om; b_j = sigma_b (z_b - mean z_b) gives z_b = b_j / sigma_b.
## The AME factors (u_diag, u_lower, U_free, V_raw) start near 0: uniform on uv_range (u_diag, which must be positive, on the upper half
## of the range), drawn with seed + chain. Each chain is jittered around the means (seed + chain): N(0, jitter_sd) added to the
## unconstrained parameters, exp(N(0, jitter_sd)) multiplying the positive scales, N(0, jitter_sd / 5) added to the correlations with the
## result kept inside (-0.95, 0.95) and positive definite, so the chains differ and R-hat stays meaningful.
## Dimensions are asserted against dat_m2 (S, N1, P1, P, Pz, D_mode, joint, K); a mismatch is an error, never a silent partial init.
m1_inits_for_m2 <- function(obj_m1, dat_m2, chains, seed = V2_SEED, jitter_sd = 0.05, uv_range = c(-0.1, 0.1)) {
  stopifnot(!is.null(obj_m1$fit), !is.null(obj_m1$dat), !is.null(dat_m2$K), chains >= 1, is.finite(jitter_sd), jitter_sd >= 0,
            length(uv_range) == 2, uv_range[1] < uv_range[2], uv_range[2] > 0)
  d1 <- obj_m1$dat
  for (nm in c("S", "N1", "P1", "P", "Pz", "D_mode", "joint")) if (!identical(as.integer(d1[[nm]]), as.integer(dat_m2[[nm]])))
    stop(sprintf("m1_inits_for_m2: %s differs between the M1 fit (%s) and the M2 data (%s)", nm, d1[[nm]], dat_m2[[nm]]))
  S <- as.integer(dat_m2$S); K <- as.integer(dat_m2$K); D <- if (dat_m2$D_mode == 3L) 2L else 1L; Dre <- D + 1L; joint <- as.integer(dat_m2$joint)
  stopifnot(K >= 1, K <= S - 1)
  fit <- obj_m1$fit
  pars <- c("alpha0", "alpha", "tau", "mu0", "beta", "sigma_x", "sd_re", "sigma_b", "gamma", "c_i", "a_i", "b_j", "R_Om",
            if (D == 2L && joint == 1L) "rho_x", if (D == 2L && joint == 0L) c("r_ca1", "q_ca2"), if (D == 1L) "r_ca_1d")
  stopifnot(all(pars %in% fit@model_pars))
  ex <- rstan::extract(fit, pars = pars, permuted = TRUE)
  pm <- function(nm) { a <- ex[[nm]]; if (is.null(dim(a)) || length(dim(a)) == 1L) mean(a) else apply(a, seq_along(dim(a))[-1], mean) }
  m <- list(alpha0 = pm("alpha0"), alpha = as.numeric(pm("alpha")), tau = pm("tau"), mu0 = as.numeric(pm("mu0")),
            beta = matrix(pm("beta"), nrow = dat_m2$P, ncol = D), sigma_x = as.numeric(pm("sigma_x")), sd_re = as.numeric(pm("sd_re")),
            sigma_b = pm("sigma_b"), gamma = as.numeric(pm("gamma")), c_i = as.numeric(pm("c_i")), a_i = matrix(pm("a_i"), nrow = S, ncol = D),
            b_j = as.numeric(pm("b_j")), R_Om = matrix(pm("R_Om"), Dre, Dre))
  stopifnot(length(m$alpha) == dat_m2$P1, length(m$mu0) == D, length(m$gamma) == dat_m2$Pz, length(m$sd_re) == Dre, length(m$sigma_x) == D,
            length(m$c_i) == S, length(m$b_j) == S, all(m$sigma_x > 0), all(m$sd_re > 0), m$sigma_b > 0, all(is.finite(unlist(m))))
  if (D == 2L && joint == 1L) m$rho <- mean(ex$rho_x)
  if (D == 2L && joint == 0L) { m$r1 <- mean(ex$r_ca1); m$q2 <- mean(ex$q_ca2) }
  if (D == 1L) m$r1d <- mean(ex$r_ca_1d)
  ## the correlation clamps (jittered corr(c, a), rho_x, the q of the a models) are expected never to bind at jitter_sd = 0.05; when
  ## one does, the M1 posterior mean sits near the boundary and the init was moved: say so in the log (audit fix 17, 2026-10-26)
  n_clamped <- 0L
  clamp <- function(x, lim = 0.95, what = "correlation") {
    out <- pmin(pmax(x, -lim), lim)
    hit <- which(abs(out - x) > 0)
    if (length(hit)) { n_clamped <<- n_clamped + length(hit)
      logf("WARNING: m1_inits_for_m2 chain %d: the %s init clamp binds (%s -> %s)", ch, what, paste(sprintf("%.4f", x[hit]), collapse = ", "), paste(sprintf("%.4f", out[hit]), collapse = ", ")) }
    out
  }
  ch <- 0L
  lower_chol <- function(R) { R <- (R + t(R)) / 2; diag(R) <- 1; stopifnot(min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) > 1e-8); t(chol(R)) }
  n_copied <- length(c("alpha0", "alpha", "tau", "mu0", "beta", "sigma_x", "sd_re", "sigma_b", "gamma", "z_re", "z_b",
                       if (D == 2L && joint == 1L) c("rho_x", "L_Om_joint"), if (D == 2L && joint == 0L) c("r_ca1", "q_ca2"), if (D == 1L) "r_ca_1d"))
  inits <- vector("list", chains)
  for (ch in seq_len(chains)) {
    set.seed(seed + ch)
    jit <- function(x) x + rnorm(length(x), 0, jitter_sd)
    pos <- function(x) x * exp(rnorm(length(x), 0, jitter_sd))
    cj  <- function(x) clamp(x + rnorm(length(x), 0, jitter_sd / 5), what = "corr(c, a) / corr(a1, a2)")
    ## correlation of (c, a): jittered off-diagonals, kept PD; the a models fix corr(a1, a2) = 0 by construction
    R <- m$R_Om; R[upper.tri(R)] <- cj(R[upper.tri(R)]); R[lower.tri(R)] <- t(R)[lower.tri(R)]
    if (D == 2L && joint == 0L) R[2, 3] <- R[3, 2] <- 0
    if (min(eigen((R + t(R)) / 2, symmetric = TRUE, only.values = TRUE)$values) <= 1e-6) R <- m$R_Om      # jitter broke PD: keep the means
    L <- lower_chol(R)
    sd_re <- pos(m$sd_re); sigma_b <- pos(m$sigma_b)
    re <- cbind(m$c_i, m$a_i)                                                     # S x Dre, the posterior means of (c_i, a_i)
    z_re <- solve(diag(sd_re, Dre) %*% L, t(re))                                  # Dre x S: re' = diag(sd_re) L z_re
    z_re <- matrix(jit(as.numeric(z_re)), Dre, S)
    z_b <- jit(m$b_j / sigma_b)
    ## vectors as 1-d arrays, so a length-1 vector (mu0, sigma_x when D = 1) keeps its declared dimension in Stan
    init <- list(alpha0 = jit(m$alpha0), alpha = as.array(jit(m$alpha)), tau = jit(m$tau), mu0 = as.array(jit(m$mu0)),
                 beta = matrix(jit(as.numeric(m$beta)), dat_m2$P, D), sigma_x = as.array(pos(m$sigma_x)), sd_re = as.array(sd_re),
                 sigma_b = sigma_b, gamma = as.array(jit(m$gamma)), z_re = z_re, z_b = as.array(z_b))
    if (D == 2L && joint == 1L) { init$rho_x <- as.array(clamp(m$rho + rnorm(1, 0, jitter_sd / 5), what = "rho_x")); init$L_Om_joint <- array(L, dim = c(1L, Dre, Dre)) }
    if (D == 2L && joint == 0L) { init$r_ca1 <- as.array(R[1, 2]); q <- R[1, 3] / sqrt(1 - R[1, 2]^2); init$q_ca2 <- as.array(clamp(q, what = "q_ca2 = r2 / sqrt(1 - r1^2)")) }
    if (D == 1L) init$r_ca_1d <- as.array(R[1, 2])
    ## AME factors near 0 (u_diag must be positive: the upper half of the range)
    init$u_diag <- as.array(runif(K, uv_range[2] / 2, uv_range[2]))
    init$u_lower <- array(runif(K * (K - 1L) / 2L, uv_range[1], uv_range[2]), dim = K * (K - 1L) / 2L)
    init$U_free <- matrix(runif((S - K - 1L) * K, uv_range[1], uv_range[2]), S - K - 1L, K)
    init$V_raw <- matrix(runif(S * K, uv_range[1], uv_range[2]), S, K)
    stopifnot(all(is.finite(unlist(init))), all(init$sigma_x > 0), all(init$sd_re > 0), init$sigma_b > 0, all(abs(init$u_diag) <= uv_range[2]),
              all(abs(init$U_free) <= max(abs(uv_range))), all(abs(init$V_raw) <= max(abs(uv_range))))
    inits[[ch]] <- init
  }
  if (n_clamped > 0L) logf("WARNING: m1_inits_for_m2: %d correlation init(s) clamped to +/-0.95 over %d chains (see the lines above)", n_clamped, chains)
  attr(inits, "means") <- m
  attr(inits, "n_copied") <- n_copied
  attr(inits, "n_clamped") <- n_clamped
  inits
}
## m1_cell_for_m2(sample, lag, model_m2, spec): the cell name of the matching M1 fit (M2x -> M1x, same switches, K dropped: m1_spec(), so
## every AME K of a cell starts from the SAME M1 fit). locate_m1_fit(...): its file in out_dir, or NA when no such file exists.
m1_cell_for_m2 <- function(sample, lag, model_m2, spec) {
  stopifnot(model_family(model_m2) == "M2")
  cell_name(sample, lag, sub("^M2", "M1", model_m2), m1_spec(spec))
}
locate_m1_fit <- function(sample, lag, model_m2, spec, out_dir) {
  f <- file.path(out_dir, paste0("fit_", m1_cell_for_m2(sample, lag, model_m2, spec), ".rds"))
  if (file.exists(f)) f else NA_character_
}

## ---- usable cores (review 2026-10-26, minor): parallel::detectCores() can return NA (seen on a reviewer's Mac) ---------------------------
## n_cores(detected, slurm): the detected count when it is a positive number; otherwise SLURM_CPUS_PER_TASK when set, else 1, with a log line
n_cores <- function(detected = parallel::detectCores(), slurm = Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")) {
  n <- suppressWarnings(as.integer(detected))
  if (length(n) == 1 && !is.na(n) && n >= 1L) return(n)
  s <- suppressWarnings(as.integer(slurm))
  out <- if (length(s) == 1 && !is.na(s) && s >= 1L) s else 1L
  logf("WARNING: parallel::detectCores() returned %s; using %d core(s) (%s)", if (length(detected) == 1) as.character(detected) else "an invalid value", out,
       if (out == 1L && (length(s) != 1 || is.na(s) || s < 1L)) "fallback 1" else "SLURM_CPUS_PER_TASK")
  out
}
