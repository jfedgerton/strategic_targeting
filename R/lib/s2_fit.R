## R/lib/s2_fit.R -- stage-2 estimation functions shared by R/s2_02_fit_reference.R, R/s2_03_joint_fe.R and R/s2_05_sweep.R.
## fit_mnl(): four-category alternative-specific logit (A, B, O = originator entry; stay out = reference) with row-specific availability
## (GPT P1-1, P2-3), optional row weights (bootstrap multiplicities), Newton-Raphson with an ANALYTIC Hessian (rev 5; BFGS fallback), and
## the stage-1-style two-way covariance (candidate x dispute; covariance-level PSD repair, GPT P1-2). Requires W_REF in the calling
## environment. Returns coefficients on the original scale with the intercept at the covariate means, plus the centering constants, so
## mnl_predict() reproduces fitted probabilities from raw covariates (GPT P2-9).
S2_CATS_FIT <- c(A = "A", B = "B", O = "orig")
## A: at W = 0 available in every month (s2_mask: no window), so every row; at W >= 1 only if the dispute runs past month W (audit 2026-11-24)
mnl_avail <- function(X, alts) cbind(A = if (W_REF == 0L) rep(TRUE, nrow(X)) else X$H > W_REF, B = rep(TRUE, nrow(X)), O = rep(W_REF >= 1L, nrow(X)) & X$likely == 1)[, alts, drop = FALSE]

fit_mnl <- function(X, design, label = "", vc = TRUE, w = NULL, start = NULL) {
  if (is.null(w)) w <- rep(1, nrow(X)); stopifnot(length(w) == nrow(X), all(w >= 0))
  if (W_REF == 0L && "O" %in% names(design)) { stopifnot(!any(X$y == "orig")); design$O <- NULL }   # W = 0: no originator alternative (s2_mask); O quantities come out NA (audit 2026-11-24)
  drops <- list()
  for (al in names(design)) { ev <- X[y == S2_CATS_FIT[[al]] & w > 0]; keep <- vapply(design[[al]], function(v) length(unique(ev[[v]])) > 1L, TRUE)
    drops[[al]] <- design[[al]][!keep]; design[[al]] <- design[[al]][keep] }
  allv <- unique(unlist(design)); stopifnot(all(allv %in% names(X)), !anyNA(X[, ..allv]))
  m <- vapply(allv, function(v) sum(w * X[[v]]) / sum(w), 0); s <- vapply(allv, function(v) sqrt(sum(w * (X[[v]] - sum(w * X[[v]]) / sum(w))^2) / sum(w)), 0)
  if (any(s < 1e-12)) stop(sprintf("%s: constant regressor(s) %s", label, paste(allv[s < 1e-12], collapse = " ")))
  Z <- lapply(design, function(vs) cbind(1, sweep(sweep(as.matrix(X[, ..vs]), 2, m[vs]), 2, s[vs], "/")))
  alts <- names(design)
  Y <- sapply(alts, function(al) as.numeric(X$y == S2_CATS_FIT[[al]])); if (!is.matrix(Y)) Y <- matrix(Y, ncol = length(alts), dimnames = list(NULL, alts))
  AV <- mnl_avail(X, alts)
  if (any(Y == 1 & !AV)) stop(sprintf("%s: an observed outcome is unavailable under the one-shot mask (%d rows)", label, sum(Y == 1 & !AV)))
  k <- vapply(Z, ncol, 0L); idx <- split(seq_len(sum(k)), factor(rep(alts, k), levels = alts))
  util <- function(p) sapply(alts, function(al) drop(Z[[al]] %*% p[idx[[al]]]))
  probs <- function(U) { U[!AV] <- -Inf; mx <- pmax(0, apply(U, 1, max)); E <- exp(U - mx); den <- exp(-mx) + rowSums(E); list(P = E / den, lse = mx + log(den)) }
  ## Y * U uses finite utilities: Y is 0 wherever an alternative is unavailable, so no 0 x -Inf arises
  nll <- function(p) { U <- util(p); pr <- probs(U); -sum(w * (rowSums(Y * U) - pr$lse)) }
  gr <- function(p, P = probs(util(p))$P) -unlist(lapply(alts, function(al) crossprod(Z[[al]], w * (Y[, al] - P[, al]))))
  hess <- function(p, P = probs(util(p))$P) {   # analytic: H_ab = sum_i w_i z_ia z_ib' p_ia (1[a = b] - p_ib)
    H <- matrix(0, sum(k), sum(k))
    for (a in alts) for (b in alts) { wt <- w * P[, a] * ((a == b) - P[, b]); H[idx[[a]], idx[[b]]] <- crossprod(Z[[a]] * wt, Z[[b]]) }
    (H + t(H)) / 2 }
  p <- if (!is.null(start) && length(start) == sum(k)) start else unlist(lapply(alts, function(al) c(log(max(sum(w * Y[, al]), 0.5) / sum(w)), rep(0, k[[al]] - 1L))))
  f0 <- nll(p); conv <- FALSE
  for (it in 1:100) {                            # Newton-Raphson with step halving
    P <- probs(util(p))$P; g <- gr(p, P); H <- hess(p, P)
    step <- tryCatch(solve(H, g), error = function(e) NULL); if (is.null(step)) break
    t <- 1; repeat { pn <- p - t * step; fn <- nll(pn); if (is.finite(fn) && fn <= f0 + 1e-10) break; t <- t / 2; if (t < 1e-8) break }
    if (t < 1e-8) break
    p <- pn; dec <- f0 - fn; f0 <- fn
    if (max(abs(gr(p))) < 1e-6 || dec < 1e-12) { conv <- TRUE; break } }
  gmax <- max(abs(gr(p)))
  if (!conv || gmax > 1e-3) {                     # fallback: BFGS from the current point
    o <- optim(p, nll, function(q) gr(q), method = "BFGS", control = list(maxit = 10000, reltol = 1e-14)); p <- o$par; f0 <- o$value; gmax <- max(abs(gr(p))) }
  if (gmax > 1e-3) stop(sprintf("%s: no convergence (max |grad| %.2e)", label, gmax))
  nm <- unlist(lapply(alts, function(al) paste0(al, ":", c("(int@mean)", design[[al]]))))
  sc <- unlist(lapply(alts, function(al) c(1, s[design[[al]]])))
  b <- setNames(p / sc, nm)
  V <- NULL; cond <- NA_real_
  if (vc) {
    P <- probs(util(p))$P; H <- hess(p, P); ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values; cond <- max(ev) / max(min(ev), 1e-300)
    Hi <- solve(H)
    Sc <- do.call(cbind, lapply(alts, function(al) Z[[al]] * (w * (Y[, al] - P[, al]))))   # per-row score contributions
    cl <- function(g) { S <- rowsum(Sc, g, reorder = FALSE); crossprod(S) }
    M2 <- cl(X$ccode) + cl(X$dispnum) - cl(paste(X$ccode, X$dispnum))
    V <- Hi %*% M2 %*% Hi; V <- V / outer(sc, sc); V <- (V + t(V)) / 2
    e <- eigen(V, symmetric = TRUE); n_neg <- sum(e$values < 0); neg_mass <- sum(abs(e$values[e$values < 0])) / sum(abs(e$values))
    if (n_neg) V <- e$vectors %*% diag(pmax(e$values, 0)) %*% t(e$vectors)
    dimnames(V) <- list(nm, nm); attr(V, "n_neg") <- n_neg; attr(V, "neg_share") <- neg_mass; attr(V, "eig_raw") <- e$values
  }
  list(b = b, V = V, ll = -f0, drops = drops, n = nrow(X), n_eff = sum(w), gmax = gmax, design = design, center = m, scale = s, cond = cond, par_scaled = p)
}

## predicted probabilities from raw covariates (GPT P2-9): U_a = alpha_mean + sum_j beta_j (x_j - m_j); unavailable alternatives get 0
mnl_predict <- function(f, X, b = f$b) {
  alts <- names(f$design); AV <- mnl_avail(X, alts)
  U <- sapply(alts, function(al) { vs <- f$design[[al]]; bb <- b[paste0(al, ":", vs)]
    b[[paste0(al, ":(int@mean)")]] + if (length(vs)) drop(sweep(as.matrix(X[, ..vs]), 2, f$center[vs]) %*% bb) else 0 })
  if (!is.matrix(U)) U <- matrix(U, ncol = length(alts), dimnames = list(NULL, alts))
  U[!AV] <- -Inf; mx <- pmax(0, apply(U, 1, max)); E <- exp(U - mx); E / (exp(-mx) + rowSums(E))
}

## surface contrasts (stage-1 style): utility change for a gap of +1 / -1 on one axis, the other at 0; linear in b, so SE by the delta method
SURF_TERMS <- c(A = "sA_", B = "sB_", O = "sA_")
surf_contrasts <- function(f) {
  out <- list()
  for (al in names(SURF_TERMS)) { pre <- SURF_TERMS[[al]]
    for (ax in c("M", "R")) for (sg in c(1, -1)) {
      a <- setNames(numeric(length(f$b)), names(f$b))
      tt <- paste0(al, ":", pre, c(paste0("d", ax), paste0("d", ax, "2")))
      if (!all(tt %in% names(a))) { out[[length(out) + 1L]] <- data.table(alt = al, axis = ax, gap = sg, est = NA_real_, se = NA_real_, status = "term dropped"); next }   # GPT P2-3: a dropped term is not an estimated zero
      a[tt[1]] <- sg; a[tt[2]] <- 1
      est <- sum(a * f$b); se <- if (!is.null(f$V)) sqrt(max(drop(t(a) %*% f$V %*% a), 0)) else NA_real_
      out[[length(out) + 1L]] <- data.table(alt = al, axis = ax, gap = sg, est = est, se = se, status = "ok") } }
  rbindlist(out)
}
wald <- function(f, terms) { terms <- intersect(terms, names(f$b)); if (!length(terms)) return(c(chi2 = NA, df = 0, p = NA))
  d <- f$b[terms]; S <- f$V[terms, terms, drop = FALSE]; w <- tryCatch(as.numeric(t(d) %*% solve(S, d)), error = function(e) NA_real_)
  c(chi2 = w, df = length(terms), p = if (is.na(w)) NA_real_ else pchisq(w, length(terms), lower.tail = FALSE)) }
surf_block <- unlist(lapply(names(SURF_TERMS), function(al) paste0(al, ":", SURF_TERMS[[al]], c("dM", "dM2", "dR", "dR2", "dMR"))))

## surface terms from the four gaps (shared by the fitter, sweep and checks)
add_surface <- function(X) { X[, `:=`(sA_dM = dM_kB, sA_dM2 = dM_kB^2, sA_dR = dR_kB, sA_dR2 = dR_kB^2, sA_dMR = dM_kB * dR_kB,
                                      sB_dM = dM_kA, sB_dM2 = dM_kA^2, sB_dR = dR_kA, sB_dR2 = dR_kA^2, sB_dMR = dM_kA * dR_kA)]; X }

## predicted-probability contrasts (memo M5 / GPT P2-9): average predicted probability of each alternative when the joiner's gap to the side
## it fights moves from the 25th to the 75th percentile of that gap (observed values elsewhere; surface terms rebuilt). A and O fight B (gaps
## dM_kB / dR_kB); B fights A (dM_kA / dR_kA). Uncertainty: parametric simulation from the two-way covariance (n_sim draws, seed fixed by
## the caller); reported as the point contrast and the 2.5 / 97.5 percentiles of the simulated contrasts.
prob_contrasts <- function(f, X, n_sim = 500L) {
  out <- list(); stopifnot(!is.null(f$V))
  bs <- MASS::mvrnorm(n_sim, f$b, f$V)
  for (al in names(f$design)) for (ax in c("M", "R")) {
    g <- if (al == "B") paste0("d", ax, "_kA") else paste0("d", ax, "_kB")
    q <- quantile(X[[g]], c(.25, .75))
    X1 <- copy(X); X1[, (g) := q[1]]; add_surface(X1); X2 <- copy(X); X2[, (g) := q[2]]; add_surface(X2)
    pt <- mean(mnl_predict(f, X2)[, al]) - mean(mnl_predict(f, X1)[, al])
    sims <- apply(bs, 1, function(bb) { bb <- setNames(bb, names(f$b)); mean(mnl_predict(f, X2, bb)[, al]) - mean(mnl_predict(f, X1, bb)[, al]) })
    out[[length(out) + 1L]] <- data.table(alt = al, gap = g, q25 = q[1], q75 = q[2], base_p = mean(mnl_predict(f, X)[, al]), contrast = pt,
                                          lo = quantile(sims, .025), hi = quantile(sims, .975)) }
  rbindlist(out)
}

## two-way pigeonhole bootstrap weights (candidate state x dispute, as stage 1): each draw resamples candidate states and disputes independently
## with replacement; a row's weight is the product of its candidate's and its dispute's multiplicities
pigeonhole_w <- function(X) {
  cs <- sort(unique(X$ccode)); ds <- sort(unique(X$dispnum))
  mc <- tabulate(sample.int(length(cs), length(cs), replace = TRUE), length(cs)); md <- tabulate(sample.int(length(ds), length(ds), replace = TRUE), length(ds))
  mc[match(X$ccode, cs)] * md[match(X$dispnum, ds)]
}
