## R/lib/sry.R -- book-aligned "sry" tie score for one layer (PI 2026-11-22), used by R/02_build_positions.R and tests/test_sry.R.
## GPT review 2026-11-22 (docs/REVIEW_sry.md) fixes:
##   SRY1 numerical tolerance: raw residuals below tol are zero; a layer-year whose residual scale is below tol is applicable with score 0;
##        the combined signed score is zeroed below COMBO_TOL before its positive part is normalized (sry_combo).
##   SRY2 fit validity: any warning raised by fixest stops the arm; GLM convergence (convStatus) and finite fitted values are required.
##   SRY3 omissions: required predictors must be finite on applicable rows; every row fixest omits must belong to a fixed-effect group it
##        reports as removed (perfect fit), or (Poisson) be a zero outcome removed by its separation check; anything else stops.
##        A layer whose applicable outcome is constant (incl. all perfectly predicted) gets score 0 without a fit.
suppressMessages(library(fixest))
SRY_COMBO_TOL <- 1e-10

## v: outcome per dyad-year row (raw counts for pois, 0/1 for logit, numeric for ols); ok: applicable rows; X: data.frame with a, b, year,
## lcap, contig_land per dyad-year row; model: "W6srybook" | "W7sryyear"; fam: "logit" | "pois" | "ols"; yr_split: row indices by year.
## Returns the signed score (NA where not applicable) with attr "info".
sry_layer <- function(v, ok, X, model, fam, yr_split, glm_iter = 25L) {
  stopifnot(length(v) == nrow(X), length(ok) == nrow(X), model %in% c("W6srybook", "W7sryyear"), fam %in% c("logit", "pois", "ols"))
  ok <- ok & !is.na(v); n1 <- sum(ok)
  x <- rep(NA_real_, length(v)); if (n1 == 0L) { attr(x, "info") <- "no applicable rows"; return(x) }
  if (fam == "logit") stopifnot(all(v[ok] %in% c(0, 1)))
  if (fam == "pois")  stopifnot(all(v[ok] >= 0), all(abs(v[ok] - round(v[ok])) < 1e-8))
  need <- if (model == "W6srybook") c("lcap") else c("lcap", "contig_land")
  for (cc in need) if (any(!is.finite(X[[cc]][ok]))) stop(sprintf("sry: non-finite %s on %d applicable rows", cc, sum(!is.finite(X[[cc]][ok]))))   # SRY3
  d <- X[ok, c("a", "b", "year", need)]; d$v <- v[ok]
  dd <- rbind(d, transform(d, a = d$b, b = d$a))                  # doubled: original orientation first, reverse second
  fml <- if (model == "W6srybook") v ~ lcap | a + b + year else v ~ lcap + contig_land | year
  scale_v <- if (fam == "logit") 1 else max(1, max(abs(dd$v)))
  tol <- 1e-8 * scale_v
  ## a layer whose every applicable row is perfectly predicted (logit all 0 / all 1; pois all 0) has no fit: score 0 everywhere applicable
  if (length(unique(dd$v)) == 1L) {                              # constant outcome (any family): fitted exactly by the intercept, no ties
    x[ok] <- 0; attr(x, "info") <- sprintf("n %d: constant outcome (%g), score 0", n1, dd$v[1]); return(x) }
  warns <- character(0)
  m <- withCallingHandlers(
    switch(fam, ols = feols(fml, data = dd, notes = FALSE, warn = TRUE),
                logit = feglm(fml, data = dd, family = binomial("logit"), notes = FALSE, warn = TRUE, glm.iter = glm_iter),
                pois = fepois(fml, data = dd, notes = FALSE, warn = TRUE, glm.iter = glm_iter)),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (length(warns)) stop(sprintf("sry %s/%s: fixest warning(s): %s", model, fam, paste(unique(warns), collapse = " | ")))   # SRY2
  if (fam != "ols" && !isTRUE(m$convStatus)) stop(sprintf("sry %s/%s: GLM did not converge (%s iterations)", model, fam, m$iterations))
  idx <- obs(m); fv <- fitted(m)
  stopifnot(!anyDuplicated(idx), all(idx >= 1L & idx <= nrow(dd)), length(fv) == length(idx), all(is.finite(fv)))
  ex <- rep(NA_real_, nrow(dd)); ex[idx] <- fv
  om <- which(is.na(ex)); n_fe <- 0L; n_sep <- 0L
  if (length(om)) {                                               # SRY3: reconcile every omitted row
    rem <- m$fixef_removed; in_rem <- rep(FALSE, length(om))
    if (!is.null(rem)) for (k in names(rem)) if (length(rem[[k]])) in_rem <- in_rem | as.character(dd[[k]][om]) %in% as.character(rem[[k]])
    n_fe <- sum(in_rem)
    sep <- !in_rem & fam == "pois" & dd$v[om] == 0                # Poisson separation removals: zero outcomes fitted at the boundary
    n_sep <- sum(sep)
    if (any(!in_rem & !sep)) stop(sprintf("sry %s/%s: %d omitted rows are not explained by removed FE groups or Poisson separation",
                                          model, fam, sum(!in_rem & !sep)))
    ex[om] <- dd$v[om]                                            # perfectly predicted: residual 0
  }
  r2 <- dd$v - ex
  asym <- max(abs(r2[seq_len(n1)] - r2[n1 + seq_len(n1)])); if (asym > 1e-4) stop(sprintf("sry: orientations disagree (%.2g)", asym))
  r <- (r2[seq_len(n1)] + r2[n1 + seq_len(n1)]) / 2
  r[abs(r) < tol] <- 0                                            # SRY1: numerical residuals are zero
  x[ok] <- r
  n_flat <- 0L
  for (ix in yr_split) { j <- ix[ok[ix]]; if (length(j)) { mx <- max(abs(x[j]))
    if (mx < tol) { x[j] <- 0; n_flat <- n_flat + 1L } else x[j] <- x[j] / mx } }   # SRY1: below-tolerance layer-years stay applicable, score 0
  attr(x, "info") <- sprintf("n %d, removed FE groups %d rows, Poisson separation %d rows, flat layer-years %d, share > 0 %.3f, GLM iterations %s",
                             n1, n_fe, n_sep, n_flat, mean(x[ok] > 0), if (fam == "ols") "-" else m$iterations)
  x
}

## combo: mean over the layers that apply to each dyad-year (signed), tiny cancellations zeroed, positive part rescaled to [0,1] within year
sry_combo <- function(M, yr_split) {
  k <- rowSums(!is.na(M)); sc <- ifelse(k > 0, rowSums(M, na.rm = TRUE) / pmax(k, 1), 0)
  sc[abs(sc) < SRY_COMBO_TOL] <- 0                                # SRY1
  S <- pmax(sc, 0); E <- numeric(length(S)); for (ix in yr_split) { m <- max(S[ix]); if (m > 0) E[ix] <- S[ix] / m }
  E
}
