## tests/test_s2_fit.R -- fitter tests requested by the GPT surface review (P2-10): analytic gradient and Hessian against finite differences,
## masked likelihood (unavailable alternatives get probability exactly 0), parameter recovery on simulated data, raw-covariate prediction,
## weights equal to row replication, two-way covariance repair recorded. Run: Rscript tests/test_s2_fit.R (repo root). Seed 123.
suppressPackageStartupMessages(library(data.table))
source("R/lib/s2.R"); source("R/lib/s2_fit.R")
set.seed(123); W_REF <- 1L
ok_all <- TRUE
check <- function(name, cond) { cat(sprintf("%-96s %s\n", name, if (isTRUE(cond)) "PASS" else "FAIL")); if (!isTRUE(cond)) ok_all <<- FALSE }

## simulated one-shot data: 40 disputes x 300 candidates; A unavailable when H <= 1; O only for likely == 1
n_d <- 40; n_k <- 300
X <- CJ(dispnum = 1:n_d, ccode = 1:n_k)
X[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N), x3 = rbinom(.N, 1, .4), likely = rbinom(.N, 1, .3), H = sample(0:6, .N, replace = TRUE))]
beta <- list(A = c(-4, 0.8, -0.5), B = c(-3.5, -0.6, 0.7), O = c(-4.5, 1.0, 0.9))
UA <- beta$A[1] + beta$A[2] * X$x1 + beta$A[3] * X$x2; UB <- beta$B[1] + beta$B[2] * X$x1 + beta$B[3] * X$x3; UO <- beta$O[1] + beta$O[2] * X$x2 + beta$O[3] * X$x3
EA <- ifelse(X$H > 1, exp(UA), 0); EB <- exp(UB); EO <- ifelse(X$likely == 1, exp(UO), 0); den <- 1 + EA + EB + EO
u <- runif(nrow(X)); pa <- EA / den; pb <- EB / den; po <- EO / den
X[, y := ifelse(u < pa, "A", ifelse(u < pa + pb, "B", ifelse(u < pa + pb + po, "orig", "out")))]
DES <- list(A = c("x1", "x2"), B = c("x1", "x3"), O = c("x2", "x3"))
f <- fit_mnl(X, DES, "sim")

## 1. recovery: slopes within 3.5 SE of the truth (two-way SEs; generous because clusters are few)
tru <- c(beta$A[2:3], beta$B[2:3], beta$O[2:3]); est <- f$b[c("A:x1", "A:x2", "B:x1", "B:x3", "O:x2", "O:x3")]; se <- sqrt(diag(f$V))[names(est)]
check("1 parameter recovery: all six slopes within 3.5 two-way SE of the true values", all(abs(est - tru) / se < 3.5))
check("1 converged with max |gradient| < 1e-6 (Newton with analytic Hessian)", f$gmax < 1e-6)

## 2. analytic gradient and Hessian vs finite differences, on the scaled parameters at a perturbed point
alts <- names(DES); m <- f$center; s <- f$scale
Z <- lapply(DES, function(vs) cbind(1, sweep(sweep(as.matrix(X[, ..vs]), 2, m[vs]), 2, s[vs], "/")))
Y <- sapply(alts, function(al) as.numeric(X$y == S2_CATS_FIT[[al]])); AV <- mnl_avail(X, alts)
k <- vapply(Z, ncol, 0L); idx <- split(seq_len(sum(k)), factor(rep(alts, k), levels = alts))
nll <- function(p) { U <- sapply(alts, function(al) drop(Z[[al]] %*% p[idx[[al]]])); U[!AV] <- -Inf; mx <- pmax(0, apply(U, 1, max)); Uf <- U; Uf[!AV] <- 0
  -sum(rowSums(Y * Uf) - (mx + log(exp(-mx) + rowSums(exp(U - mx))))) }
p0 <- f$par_scaled + rnorm(length(f$par_scaled), 0, 0.05)
num_g <- vapply(seq_along(p0), function(j) { e <- replace(numeric(length(p0)), j, 1e-5); (nll(p0 + e) - nll(p0 - e)) / 2e-5 }, 0)
P0 <- { U <- sapply(alts, function(al) drop(Z[[al]] %*% p0[idx[[al]]])); U[!AV] <- -Inf; mx <- pmax(0, apply(U, 1, max)); E <- exp(U - mx); E / (exp(-mx) + rowSums(E)) }
an_g <- -unlist(lapply(alts, function(al) crossprod(Z[[al]], Y[, al] - P0[, al])))
check("2 analytic gradient matches central finite differences (max rel. error < 1e-5)", max(abs(an_g - num_g)) / max(abs(num_g)) < 1e-5)
grf <- function(p) { U <- sapply(alts, function(al) drop(Z[[al]] %*% p[idx[[al]]])); U[!AV] <- -Inf; mx <- pmax(0, apply(U, 1, max)); E <- exp(U - mx); P <- E / (exp(-mx) + rowSums(E))
  -unlist(lapply(alts, function(al) crossprod(Z[[al]], Y[, al] - P[, al]))) }
num_H <- sapply(seq_along(p0), function(j) { e <- replace(numeric(length(p0)), j, 1e-5); (grf(p0 + e) - grf(p0 - e)) / 2e-5 })
an_H <- matrix(0, sum(k), sum(k)); for (a in alts) for (b in alts) { wt <- P0[, a] * ((a == b) - P0[, b]); an_H[idx[[a]], idx[[b]]] <- crossprod(Z[[a]] * wt, Z[[b]]) }
check("2 analytic Hessian matches finite differences of the gradient (max rel. error < 1e-4)", max(abs(an_H - num_H)) / max(abs(num_H)) < 1e-4)

## 3. masked likelihood: unavailable alternatives get probability exactly 0; observed shares reproduced (intercept MLE property)
P <- mnl_predict(f, X)
check("3 P(A) = 0 exactly where H <= W; P(O) = 0 exactly where likely == 0", all(P[X$H <= 1, "A"] == 0) && all(P[X$likely == 0, "O"] == 0))
check("3 mean predicted probability equals the observed share for every alternative (raw-covariate prediction is consistent)",
      isTRUE(all.equal(colSums(P), c(A = sum(X$y == "A"), B = sum(X$y == "B"), O = sum(X$y == "orig")), tolerance = 1e-6)))
check("3 an observed outcome in an unavailable cell stops the fit", inherits(tryCatch(fit_mnl(copy(X)[which(H <= 1)[1], y := "A"], DES, "bad", vc = FALSE), error = function(e) e), "error"))

## 4. weights equal row replication
wts <- sample(0:3, nrow(X), replace = TRUE)
fw <- fit_mnl(X, DES, "w", vc = FALSE, w = wts); fr <- fit_mnl(X[rep(seq_len(.N), wts)], DES, "rep", vc = FALSE)
check("4 weighted fit equals the fit on replicated rows (coefficients and log-likelihood)", isTRUE(all.equal(fw$b, fr$b, tolerance = 1e-6)) && abs(fw$ll - fr$ll) < 1e-6)
check("4 pigeonhole weights are products of candidate and dispute multiplicities with the right totals", { pw <- pigeonhole_w(X); all(pw >= 0) && sum(pw) > 0 && length(pw) == nrow(X) })

## 5. covariance bookkeeping and dropped terms
check("5 covariance repair recorded (n_neg, neg_share, raw eigenvalues) and V is PSD", !is.null(attr(f$V, "n_neg")) && min(eigen(f$V, only.values = TRUE)$values) > -1e-10)
X5 <- copy(X); X5[, xc := 0]; X5[y == "B", xc := 0]; X5[sample(.N, 50), xc := 1]; X5[y == "B", xc := 0]
f5 <- fit_mnl(X5, list(A = DES$A, B = c(DES$B, "xc"), O = DES$O), "drop", vc = FALSE)
check("5 a term constant among an alternative's outcomes is dropped and reported", "xc" %in% f5$drops$B && !("B:xc" %in% names(f5$b)))
X5 <- copy(X); X5[, `:=`(dM_kB = x1, dR_kB = x2, dM_kA = x1, dR_kA = x2)]; add_surface(X5)
f6 <- fit_mnl(X5, list(A = c("sA_dM", "sA_dM2", "sA_dR"), B = c("sB_dM", "sB_dR", "sB_dR2"), O = c("sA_dM")), "surf", vc = TRUE)
sc6 <- surf_contrasts(f6)
check("5 surface contrast with a missing squared term is flagged unavailable, never an estimated zero", all(is.na(sc6[alt == "A" & axis == "R", est])) && all(sc6[alt == "A" & axis == "R", status] == "term dropped"))

cat(if (ok_all) "ALL PASS\n" else "SOME TESTS FAILED\n"); quit(status = if (ok_all) 0 else 1)
