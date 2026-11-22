## R/s2_02_fit_reference.R -- stage-2 reference one-shot model and its opposing-side robustness checks (memo section 13, surface amendment).
## Four-category alternative-specific logit (A = join initiating side, B = join target side, O = originator entry; stay out = reference) on the
## reference estimation cohort (in_ref), with the frozen design S2_REF_DESIGN (stage-1 surface toward the opposing side + mirrored stage-1
## controls + recruitment terms). Inference: two-way cluster-robust (candidate state x dispute; CGM, negative eigenvalues set to 0).
## Reported quantities: coefficients and the stage-1-style surface contrasts per alternative (utility change for a +1 / -1 gap on each axis),
## each with a two-way SE.
## Robustness (PI 2026-11-22): A1-A4 opposing-side variants (40), design = linear, B1 permutation of candidates' gap vectors within dispute,
## B2 permutation of side labels among third-party entrants within dispute, B3 alliance- / contiguity-correlated measurement error in positions
## (both signs, k = 0.25 / 0.5 / 1) and random noise.
## Separation rule (fixed before fitting): a term is dropped from an alternative's equation when every observed outcome of that alternative
## has the same value of the term (the MLE does not exist); drops are logged.
## Usage: Rscript R/s2_02_fit_reference.R --build=DIR --out=DIR [--n_perm=100] [--checks=all|none]
suppressPackageStartupMessages(library(data.table))
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file)), "lib", "s2.R"))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
fl <- sub("^--([a-z_]+)=.*$", "\\1", args); if (any(!fl %in% c("build", "out", "n_perm", "checks")) || anyDuplicated(fl)) stop("unknown or duplicated option")
BUILD <- opt("build"); OUT <- opt("out"); N_PERM <- as.integer(opt("n_perm", "100")); CHECKS <- opt("checks", "all")
stopifnot(!is.null(BUILD), !is.null(OUT), CHECKS %in% c("all", "none"), N_PERM >= 0)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); set.seed(123)
LOG <- file.path(OUT, "s2_fit.log"); cat("", file = LOG)
logf <- function(...) { s <- sprintf(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOG, append = TRUE) }
GIT <- tryCatch(system2("git", c("-C", dirname(normalizePath(this_file)), "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
logf("stage-2 reference fit %s; build %s; git %s; R %s", format(Sys.time()), BUILD, GIT, getRversion())

## ---- data ----------------------------------------------------------------------------------------------------------------------------
META <- readRDS(file.path(BUILD, "s2_meta.rds")); stopifnot(META$mode == "release", META$fatal_def == "stage1", META$window == 1L)
ONE <- as.data.table(readRDS(file.path(BUILD, "s2_oneshot_w1.rds")))
VARS <- as.data.table(readRDS(file.path(BUILD, "s2_surface_variants.rds"))); SI <- readRDS(file.path(BUILD, "s2_surface_inputs.rds"))
X0 <- ONE[in_ref == TRUE]
for (v in names(X0)) if (is.logical(X0[[v]]) && v %in% unlist(c(S2_REF_DESIGN, S2_REF_DESIGN_LINEAR))) X0[, (v) := as.numeric(get(v))]
logf("reference cohort %d rows, %d disputes, %d candidate states; outcomes %s", nrow(X0), uniqueN(X0$dispnum), uniqueN(X0$ccode), paste(sprintf("%s %d", names(table(X0$y)), table(X0$y)), collapse = ", "))
add_surface <- function(X) { X[, `:=`(sA_dM = dM_kB, sA_dM2 = dM_kB^2, sA_dR = dR_kB, sA_dR2 = dR_kB^2, sA_dMR = dM_kB * dR_kB,
                                      sB_dM = dM_kA, sB_dM2 = dM_kA^2, sB_dR = dR_kA, sB_dR2 = dR_kA^2, sB_dMR = dM_kA * dR_kA)]; X }

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
  k <- vapply(Z, ncol, 0L); idx <- split(seq_len(sum(k)), rep(names(k), k))
  util <- function(p) sapply(names(Z), function(al) drop(Z[[al]] %*% p[idx[[al]]]))
  probs <- function(U) { mx <- pmax(0, apply(U, 1, max)); E <- exp(U - mx); den <- exp(-mx) + rowSums(E); list(P = E / den, lse = mx + log(den)) }
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
    e <- eigen(M2, symmetric = TRUE); n_neg <- sum(e$values < 0); M2 <- e$vectors %*% diag(pmax(e$values, 0)) %*% t(e$vectors)
    V <- Hi %*% M2 %*% Hi; V <- V / outer(sc, sc); dimnames(V) <- list(nm, nm); attr(V, "n_neg") <- n_neg
  }
  list(b = b, V = V, ll = -o$value, drops = drops, n = nrow(X), gmax = gmax, design = design)
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

## ---- 1. reference fit --------------------------------------------------------------------------------------------------------------
X0 <- add_surface(X0)
F0 <- fit_mnl(X0, S2_REF_DESIGN, "reference")
for (al in names(F0$drops)) if (length(F0$drops[[al]])) logf("separation rule: dropped from %s: %s", al, paste(F0$drops[[al]], collapse = ", "))
se0 <- sqrt(diag(F0$V))
COEF <- data.table(term = names(F0$b), est = F0$b, se = se0, z = F0$b / se0, p = 2 * pnorm(-abs(F0$b / se0)))
fwrite(COEF, file.path(OUT, "reference_coefs.csv"))
SC0 <- surf_contrasts(F0); SC0[, `:=`(lo = est - 1.96 * se, hi = est + 1.96 * se)]
fwrite(SC0, file.path(OUT, "reference_surface_contrasts.csv"))
wb <- wald(F0, surf_block)
logf("reference fit: n %d, logLik %.2f, max |grad| %.1e, CGM negative eigenvalues %d; surface block Wald chi2 %.2f on %d df, p = %.4f", F0$n, F0$ll, F0$gmax, attr(F0$V, "n_neg"), wb["chi2"], wb["df"], wb["p"])
logf("surface contrasts (utility change; gap = joiner minus the side it fights; 95%% two-way intervals):")
for (r in seq_len(nrow(SC0))) logf("  %s  %s %+d: %7.3f  [%7.3f, %7.3f]", SC0$alt[r], SC0$axis[r], SC0$gap[r], SC0$est[r], SC0$lo[r], SC0$hi[r])
logf("coefficients (two-way SE):"); for (r in seq_len(nrow(COEF))) logf("  %-22s %8.3f  (%.3f)  p %.3f", COEF$term[r], COEF$est[r], COEF$se[r], COEF$p[r])
if (CHECKS == "none") { logf("S2_FIT_DONE (checks skipped)"); quit(status = 0) }

## ---- 2. A1-A4 opposing-side variants and the linear design ----------------------------------------------------------------------------
keys <- sub("^dR_kB__", "", grep("^dR_kB__", names(VARS), value = TRUE))
VR <- list()
for (key in keys) {
  Xv <- merge(X0[, !c("dM_kA", "dR_kA", "dM_kB", "dR_kB", paste0("sA_", c("dM", "dM2", "dR", "dR2", "dMR")), paste0("sB_", c("dM", "dM2", "dR", "dR2", "dMR"))), with = FALSE],
              VARS[, c("dispnum", "ccode", paste0(c("dM_kA", "dR_kA", "dM_kB", "dR_kB"), "__", key)), with = FALSE], by = c("dispnum", "ccode"))
  setnames(Xv, paste0(c("dM_kA", "dR_kA", "dM_kB", "dR_kB"), "__", key), c("dM_kA", "dR_kA", "dM_kB", "dR_kB"))
  n_all <- nrow(Xv); Xv <- Xv[is.finite(dM_kA) & is.finite(dR_kA) & is.finite(dM_kB) & is.finite(dR_kB)]; Xv <- add_surface(Xv)
  f <- tryCatch(fit_mnl(Xv, S2_REF_DESIGN, key), error = function(e) { logf("variant %s failed: %s", key, conditionMessage(e)); NULL })
  if (!is.null(f)) VR[[key]] <- surf_contrasts(f)[, `:=`(variant = key, n = f$n, dropped_rows = n_all - f$n)]
}
FL <- tryCatch(fit_mnl(X0, S2_REF_DESIGN_LINEAR, "linear"), error = function(e) { logf("linear design failed: %s", conditionMessage(e)); NULL })
VRt <- rbindlist(VR); fwrite(VRt, file.path(OUT, "variants_surface_contrasts.csv"))
VRt[, `:=`(lo = est - 1.96 * se, hi = est + 1.96 * se)]
sumv <- VRt[, .(variants = .N, same_sign = mean(sign(est) == sign(SC0$est[SC0$alt == alt[1] & SC0$axis == axis[1] & SC0$gap == gap[1]])), excl0 = mean(lo > 0 | hi < 0),
                min = min(est), max = max(est)), by = .(alt, axis, gap)]
logf("A1-A4: %d of %d variants fitted; per contrast, share with the reference sign / share whose 95%% interval excludes 0 / range:", length(VR), length(keys))
for (r in seq_len(nrow(sumv))) logf("  %s %s %+d: same sign %.2f, excludes 0 %.2f, range [%.2f, %.2f]", sumv$alt[r], sumv$axis[r], sumv$gap[r], sumv$same_sign[r], sumv$excl0[r], sumv$min[r], sumv$max[r])
if (!is.null(FL)) { cl <- data.table(term = names(FL$b), est = FL$b, se = sqrt(diag(FL$V))); fwrite(cl, file.path(OUT, "linear_design_coefs.csv"))
  logf("design = linear: logLik %.2f (surface %.2f); terms %s", FL$ll, F0$ll, paste(sprintf("%s %.2f (%.2f)", cl$term, cl$est, cl$se)[grepl("cent|ratio", cl$term)], collapse = "; ")) }

## ---- 3. B1 / B2 permutations (statistic = likelihood-ratio of the 15-term surface block; no vcov needed) ---------------------------
NULL_DESIGN <- lapply(S2_REF_DESIGN, function(v) setdiff(v, c(S2_SURF_A, S2_SURF_B)))
lr_stat <- function(X, ll_null = NULL) {
  f <- tryCatch(fit_mnl(X, S2_REF_DESIGN, "perm", vc = FALSE), error = function(e) NULL); if (is.null(f)) return(NA_real_)
  if (is.null(ll_null)) { g <- tryCatch(fit_mnl(X, NULL_DESIGN, "perm-null", vc = FALSE), error = function(e) NULL); if (is.null(g)) return(NA_real_); ll_null <- g$ll }
  2 * (f$ll - ll_null) }
LL0 <- fit_mnl(X0, NULL_DESIGN, "null", vc = FALSE)$ll; LR_OBS <- 2 * (F0$ll - LL0)
logf("surface block likelihood ratio (15 df): %.2f", LR_OBS)
P1 <- P2 <- rep(NA_real_, N_PERM)
gcols <- c("dM_kA", "dR_kA", "dM_kB", "dR_kB")
for (b in seq_len(N_PERM)) {
  Xp <- copy(X0); Xp[, ri := sample.int(.N), by = dispnum]                 # B1: shuffle candidates' gap vectors within dispute
  ord <- Xp[, .I[ri], by = dispnum]$V1
  Xp[, (gcols) := X0[ord, ..gcols]]; Xp <- add_surface(Xp); P1[b] <- lr_stat(Xp, LL0)    # the null model has no surface: unchanged
  Xq <- copy(X0); ent <- Xq[y %in% c("A", "B"), which = TRUE]                # B2: shuffle A/B labels among third-party entrants within dispute
  Xq[ent, y := sample(y), by = dispnum]; P2[b] <- lr_stat(Xq)
  if (b %% 20 == 0) logf("  permutations %d / %d", b, N_PERM)
}
fwrite(data.table(perm = seq_len(N_PERM), B1_within_dispute_gaps = P1, B2_side_labels = P2), file.path(OUT, "permutations.csv"))
logf("B1 (gap vectors shuffled within dispute): observed LR %.2f; permutation p = %.3f (%d valid; permuted LR median %.2f)", LR_OBS, (1 + sum(P1 >= LR_OBS, na.rm = TRUE)) / (1 + sum(!is.na(P1))), sum(!is.na(P1)), median(P1, na.rm = TRUE))
logf("B2 (side labels shuffled among entrants within dispute): observed LR %.2f; permutation p = %.3f (%d valid; permuted LR median %.2f)", LR_OBS, (1 + sum(P2 >= LR_OBS, na.rm = TRUE)) / (1 + sum(!is.na(P2))), sum(!is.na(P2)), median(P2, na.rm = TRUE))

## ---- 4. B3 measurement error mimicking confounding (positions perturbed, surface rebuilt) ---------------------------------------------
PV <- as.data.table(SI$PV); CONF <- as.data.table(SI$conf)
PV <- merge(PV, CONF, by = c("ccode", "year"), all.x = TRUE); PV[is.na(n_def), n_def := 0L]; PV[is.na(n_contig), n_contig := 0L]
zyr <- function(x, yr) ave(x, yr, FUN = function(v) if (sd(v) > 0) (v - mean(v)) / sd(v) else 0 * v)
PV[, `:=`(z_def = zyr(log1p(n_def), year), z_contig = zyr(log1p(n_contig), year), z_rand = rnorm(.N))]
MER <- list()
for (cf in c("def", "contig", "rand")) for (sg in if (cf == "rand") 1 else c(1, -1)) for (k in c(0.25, 0.5, 1)) {
  Pp <- copy(PV)[!is.na(R)]; Pp[, R := R + sg * k * ave(R, year, FUN = sd) * get(paste0("z_", cf))]
  sv <- s2_surface(SI$rows, SI$mem_roster, Pp, SI$sM, SI$sR, "mean")
  Xm <- merge(X0[, !c(gcols, paste0("sA_", c("dM", "dM2", "dR", "dR2", "dMR")), paste0("sB_", c("dM", "dM2", "dR", "dR2", "dMR"))), with = FALSE], sv[, c("dispnum", "ccode", gcols), with = FALSE], by = c("dispnum", "ccode"))
  Xm <- add_surface(Xm[is.finite(dR_kA) & is.finite(dR_kB)])
  rk <- mean(vapply(split(seq_len(nrow(Pp)), Pp$year), function(ix) { o <- PV[!is.na(R)][ix]; if (length(ix) > 3) suppressWarnings(cor(o$R, Pp$R[ix], method = "spearman")) else NA }, 0), na.rm = TRUE)
  f <- tryCatch(fit_mnl(Xm, S2_REF_DESIGN, paste("merror", cf, sg, k)), error = function(e) NULL)
  if (!is.null(f)) MER[[length(MER) + 1L]] <- surf_contrasts(f)[, `:=`(confounder = cf, sign = sg, k = k, rank_corr = rk)]
}
MERt <- rbindlist(MER); fwrite(MERt, file.path(OUT, "merror_surface_contrasts.csv"))
logf("B3: %d perturbed fits (alliance / contiguity confounders, both signs, k = 0.25 / 0.5 / 1; random noise); R contrasts by fit (A R+1, B R+1, O R+1):", length(MER))
for (g in split(MERt, by = c("confounder", "sign", "k"))) logf("  %-6s sign %+d k %.2f (rank corr %.2f): A %.2f, B %.2f, O %.2f", g$confounder[1], g$sign[1], g$k[1], g$rank_corr[1],
  g[alt == "A" & axis == "R" & gap == 1, est], g[alt == "B" & axis == "R" & gap == 1, est], g[alt == "O" & axis == "R" & gap == 1, est])
logf("S2_FIT_DONE")
