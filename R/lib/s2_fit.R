## R/lib/s2_fit.R -- stage-2 estimation functions shared by R/s2_02_fit_reference.R and R/s2_03_joint_fe.R (moved verbatim from the fitter).
## fit_mnl(): four-category alternative-specific logit with row-specific availability (GPT P1-1) and stage-1-style two-way covariance repair
## (GPT P1-2). Requires W_REF in the calling environment (the reference window). Also returns the centering constants (GPT P2-9).
## ---- model -----------------------------------------------------------------------------------------------------------------------------
## design: list(A, B, O) of column names; returns coefficients on the original scale, two-way vcov, log-likelihood, dropped terms
fit_mnl <- function(X, design, label = "", vc = TRUE) {
  cats <- c(A = "A", B = "B", O = "orig")
  drops <- list()
  for (al in names(design)) { ev <- X[y == cats[[al]]]; keep <- vapply(design[[al]], function(v) length(unique(ev[[v]])) > 1L, TRUE)
    drops[[al]] <- design[[al]][!keep]; design[[al]] <- design[[al]][keep] }
  allv <- unique(unlist(design)); stopifnot(all(allv %in% names(X)), !anyNA(X[, ..allv]))
  m <- vapply(allv, function(v) mean(X[[v]]), 0); s <- vapply(allv, function(v) sd(X[[v]]), 0)
  if (any(s < 1e-12)) stop(sprintf("%s: constant regressor(s) %s", label, paste(allv[s < 1e-12], collapse = " ")))
  Z <- lapply(design, function(vs) cbind(1, sweep(sweep(as.matrix(X[, ..vs]), 2, m[vs]), 2, s[vs], "/")))
  Y <- sapply(names(design), function(al) as.numeric(X$y == cats[[al]]))
  ## GPT P1-1: row-specific availability (memo N1, one-shot): A only when the horizon reaches past the window (H > W); originator entry only
  ## when W >= 1 (month 0 always lies inside the window); B always. Unavailable alternatives get probability exactly 0.
  AV <- cbind(A = X$H > W_REF, B = rep(TRUE, nrow(X)), O = rep(W_REF >= 1L, nrow(X)))[, names(design), drop = FALSE]
  if (any(Y == 1 & !AV)) stop(sprintf("%s: an observed outcome is unavailable under the one-shot mask (%d rows)", label, sum(Y == 1 & !AV)))
  k <- vapply(Z, ncol, 0L); idx <- split(seq_len(sum(k)), rep(names(k), k))
  util <- function(p) sapply(names(Z), function(al) drop(Z[[al]] %*% p[idx[[al]]]))
  probs <- function(U) { U[!AV] <- -Inf; mx <- pmax(0, apply(U, 1, max)); E <- exp(U - mx); den <- exp(-mx) + rowSums(E); list(P = E / den, lse = mx + log(den)) }
  ## Y * U below uses the finite utilities: Y is 0 wherever an alternative is unavailable, so no 0 x -Inf arises
  nll <- function(p) { U <- util(p); pr <- probs(U); -sum(rowSums(Y * U) - pr$lse) }
  gr <- function(p) { pr <- probs(util(p))$P; -unlist(lapply(names(Z), function(al) crossprod(Z[[al]], Y[, al] - pr[, al]))) }
  st <- unlist(lapply(names(Z), function(al) c(log(mean(Y[, al])), rep(0, k[[al]] - 1L))))
  o <- optim(st, nll, gr, method = "BFGS", control = list(maxit = 10000, reltol = 1e-14))
  gmax <- max(abs(gr(o$par))); if (o$convergence != 0 || gmax > 1e-3) stop(sprintf("%s: no convergence (code %d, max |grad| %.2e)", label, o$convergence, gmax))
  nm <- unlist(lapply(names(design), function(al) paste0(al, ":", c("(int@mean)", design[[al]]))))
  sc <- unlist(lapply(names(design), function(al) c(1, s[design[[al]]])))
  b <- setNames(o$par / sc, nm)
  V <- NULL
  if (vc) {
    H <- optimHess(o$par, nll, gr)                                         # observed information on the scaled parameters
    pr <- probs(util(o$par))$P
    Sc <- do.call(cbind, lapply(names(Z), function(al) Z[[al]] * (Y[, al] - pr[, al])))   # per-row score contributions
    cl <- function(g) { S <- rowsum(Sc, g, reorder = FALSE); crossprod(S) }
    Hi <- solve(H)
    M2 <- cl(X$ccode) + cl(X$dispnum) - cl(paste(X$ccode, X$dispnum))
    ## GPT P1-2: as stage 1 (R/07 vc_2w), build the raw two-way covariance in the reported parameterization FIRST, then clip its negative
    ## eigenvalues; the repair is recorded
    V <- Hi %*% M2 %*% Hi; V <- V / outer(sc, sc); V <- (V + t(V)) / 2
    e <- eigen(V, symmetric = TRUE); n_neg <- sum(e$values < 0); neg_mass <- sum(abs(e$values[e$values < 0])) / sum(abs(e$values))
    if (n_neg) V <- e$vectors %*% diag(pmax(e$values, 0)) %*% t(e$vectors)
    dimnames(V) <- list(nm, nm); attr(V, "n_neg") <- n_neg; attr(V, "neg_share") <- neg_mass; attr(V, "eig_raw") <- e$values
  }
  list(b = b, V = V, ll = -o$value, drops = drops, n = nrow(X), gmax = gmax, design = design, center = m, scale = s)
}
## surface contrasts (stage-1 style): utility change for a gap of +1 / -1 on one axis, the other at 0; linear in b, so SE by the delta method
SURF_TERMS <- c(A = "sA_", B = "sB_", O = "sA_")
surf_contrasts <- function(f) {
  out <- list()
  for (al in names(SURF_TERMS)) { pre <- SURF_TERMS[[al]]
    for (ax in c("M", "R")) for (sg in c(1, -1)) {
      a <- setNames(numeric(length(f$b)), names(f$b))
      for (t in c(paste0("d", ax), paste0("d", ax, "2"))) { nmx <- paste0(al, ":", pre, t); if (nmx %in% names(a)) a[nmx] <- if (endsWith(t, "2")) 1 else sg }
      est <- sum(a * f$b); se <- if (!is.null(f$V)) sqrt(max(drop(t(a) %*% f$V %*% a), 0)) else NA_real_
      out[[length(out) + 1L]] <- data.table(alt = al, axis = ax, gap = sg, est = est, se = se) } }
  rbindlist(out)
}
wald <- function(f, terms) { terms <- intersect(terms, names(f$b)); if (!length(terms)) return(c(chi2 = NA, df = 0, p = NA))
  d <- f$b[terms]; S <- f$V[terms, terms, drop = FALSE]; w <- tryCatch(as.numeric(t(d) %*% solve(S, d)), error = function(e) NA_real_)
  c(chi2 = w, df = length(terms), p = if (is.na(w)) NA_real_ else pchisq(w, length(terms), lower.tail = FALSE)) }
surf_block <- unlist(lapply(names(SURF_TERMS), function(al) paste0(al, ":", SURF_TERMS[[al]], c("dM", "dM2", "dR", "dR2", "dMR"))))
