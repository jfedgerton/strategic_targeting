## tests/test_sry.R -- acceptance tests for R/lib/sry.R requested by the GPT review of 2026-11-22 (docs/REVIEW_sry.md).
## Run: Rscript tests/test_sry.R   (from the repo root). Synthetic data only; seed 123.
suppressMessages({ library(fixest) })
source("R/lib/sry.R")
set.seed(123)
ok_all <- TRUE
check <- function(name, cond) { cat(sprintf("%-78s %s\n", name, if (isTRUE(cond)) "PASS" else "FAIL")); if (!isTRUE(cond)) ok_all <<- FALSE }
expect_error <- function(expr) inherits(tryCatch({ force(expr); NULL }, error = function(e) e), "error")

## synthetic dyad-years: 12 states, years 1..4, all dyads
st <- 1:12; yrs <- 1:4
X <- do.call(rbind, lapply(yrs, function(y) { p <- t(combn(st, 2)); data.frame(a = p[, 1], b = p[, 2], year = y) }))
X$lcap <- log1p(runif(nrow(X), 100, 5000)); X$contig_land <- rbinom(nrow(X), 1, 0.2)
ys <- split(seq_len(nrow(X)), X$year); ok <- rep(TRUE, nrow(X))

## 1. exactly fitted layer (GPT SRY1 counterexample): the model reproduces the data up to floating-point error, so every raw residual is
##    numerical noise; without a tolerance the layer-year max division would turn that noise into weight-1 ties
for (mdl in c("W6srybook", "W7sryyear")) {
  v <- 0.3 * X$year + 0.7 * X$lcap + (if (mdl == "W7sryyear") 0.2 * X$contig_land else 0)   # exact under that model's design
  s <- sry_layer(v, ok, X, mdl, "ols", ys)
  check(sprintf("1 %s: exact OLS fit gives score 0 everywhere (applicable, not NA)", mdl), all(s == 0) && !anyNA(s))
}
vc <- rep(4, nrow(X))                                                 # constant nonzero counts over the whole layer
s <- sry_layer(vc, ok, X, "W7sryyear", "pois", ys); check("1b constant nonzero counts give score 0, no error", all(s == 0) && !anyNA(s))
vc2 <- rpois(nrow(X), 3); vc2[X$year == 3] <- 4                       # Poisson: one constant year among varying years
s <- sry_layer(vc2, ok, X, "W7sryyear", "pois", ys); check("1b' Poisson fit with a constant year runs, scores in [-1, 1]", all(abs(s) <= 1) && !anyNA(s))
vy <- rnorm(nrow(X)); vy[X$year == 2] <- 0.3 * 2 + 0.7 * X$lcap[X$year == 2]   # one year exact, others noisy: real residuals elsewhere
s <- sry_layer(vy, ok, X, "W7sryyear", "ols", ys)
check("1c noisy years keep real residuals (not zeroed by the tolerance)", sd(s[X$year != 2]) > 0)

## 2. perfectly predicted FE groups: state 1 never has a tie (logit) -> its rows score 0 after removal; others vary
vb <- rbinom(nrow(X), 1, 0.3); vb[X$a == 1 | X$b == 1] <- 0
s <- sry_layer(vb, ok, X, "W6srybook", "logit", ys)
check("2 mixed: removed FE group (state 1, all 0) scores 0", all(s[X$a == 1 | X$b == 1] == 0))
check("2 mixed: remaining rows vary and are in [-1, 1]", sd(s[!(X$a == 1 | X$b == 1)]) > 0 && all(abs(s) <= 1))
## 2b. wholly perfectly predicted layer: all zeros / all ones
check("2b all-0 logit layer: score 0, no error", all(sry_layer(rep(0, nrow(X)), ok, X, "W6srybook", "logit", ys) == 0))
check("2b all-1 logit layer: score 0, no error", all(sry_layer(rep(1, nrow(X)), ok, X, "W6srybook", "logit", ys) == 0))
check("2b all-0 pois layer: score 0, no error",  all(sry_layer(rep(0, nrow(X)), ok, X, "W7sryyear", "pois", ys) == 0))

## 3. missing required predictors -> stop (W7 needs contiguity; W6 does not)
Xm <- X; Xm$contig_land[5] <- NA
check("3 missing contiguity stops W7", expect_error(sry_layer(vb, ok, Xm, "W7sryyear", "logit", ys)))
check("3 missing contiguity is irrelevant to W6 (no error)", !expect_error(sry_layer(vb, ok, Xm, "W6srybook", "logit", ys)))
Xd <- X; Xd$lcap[7] <- NA
check("3 missing distance stops W6", expect_error(sry_layer(vb, ok, Xd, "W6srybook", "logit", ys)))
check("3 missing distance on a non-applicable row is ignored", { ok2 <- ok; ok2[7] <- FALSE; !expect_error(sry_layer(vb, ok2, Xd, "W6srybook", "logit", ys)) })

## 4. deliberately insufficient GLM iterations -> stop
vb2 <- rbinom(nrow(X), 1, plogis(-1 + 2 * X$contig_land))
check("4 glm.iter = 1 stops (logit)", expect_error(sry_layer(vb2, ok, X, "W7sryyear", "logit", ys, glm_iter = 1L)))
check("4 glm.iter = 1 stops (pois)",  expect_error(sry_layer(vc2, ok, X, "W7sryyear", "pois", ys, glm_iter = 1L)))
check("4 default iterations succeed", !expect_error(sry_layer(vb2, ok, X, "W7sryyear", "logit", ys)))

## 5. signed-layer cancellation -> zero, not a tie; denominators count applicable layers only
M <- cbind(c(0.3, 1e-13, 0.5, NA), c(-0.3, 0, -0.1, 0.4))
yy <- list(1:4)
E <- sry_combo(M, yy)
check("5 exact cancellation (0.3, -0.3) -> 0", E[1] == 0)
check("5 tiny residual (1e-13, 0) -> 0", E[2] == 0)
check("5 NA layer excluded from the mean (row 4 = 0.4, the year max -> 1)", E[4] == 1)
check("5 row 3 = mean(0.5, -0.1) = 0.2 -> 0.2 / 0.4 = 0.5", abs(E[3] - 0.5) < 1e-12)

## 6. row / orientation alignment: a planted positive tie is recovered for the right dyad
vb3 <- rbinom(nrow(X), 1, 0.2); tgt <- which(X$a == 3 & X$b == 7)
vb3[tgt] <- 1; vb3[X$a %in% c(3, 7) & !(seq_len(nrow(X)) %in% tgt)] <- 0
s <- sry_layer(vb3, ok, X, "W6srybook", "logit", ys)
check("6 planted 3-7 tie scores > 0 in every year", all(s[tgt] > 0))

cat(if (ok_all) "ALL PASS\n" else "SOME TESTS FAILED\n"); quit(status = if (ok_all) 0 else 1)
