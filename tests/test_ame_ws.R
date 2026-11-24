## tests/test_ame_ws.R -- the Williams-Savitsky rescaling (R/ame/ame_ws.R) must reproduce the cluster-robust sandwich of a weighted logit.
## "Posterior" draws ~ N(MLE, model-based covariance) stand in for a weighted pseudo-posterior; after rescaling, the draws' covariance must
## match sandwich::vcovCL (HC0, cluster adjustment G / (G - 1)). Seed 123.
suppressPackageStartupMessages({ library(sandwich); library(MASS) })
source("R/ame/ame_ws.R"); set.seed(123)
G <- 150; n_g <- sample(2:8, G, replace = TRUE); g <- rep(seq_len(G), n_g); n <- length(g)
x1 <- rnorm(n); x2 <- rbinom(n, 1, .4); re <- rnorm(G, 0, 1)[g]            # cluster effect makes rows within a cluster correlated
y <- rbinom(n, 1, plogis(-.5 + .8 * x1 - .6 * x2 + re)); w <- runif(n, .3, 1.5)
m <- suppressWarnings(glm(y ~ x1 + x2, family = binomial, weights = w))
X <- model.matrix(m); b <- coef(m)
llc <- function(p) { eta <- drop(X %*% p); l <- w * (y * eta - log1p(exp(eta))); rowsum(l, g)[, 1] }
S <- ws_jacobian(llc, b)
dr <- mvrnorm(20000, b, vcov(m)); colnames(dr) <- names(b)
r <- ws_rescale(dr, S)
Vcl <- vcovCL(m, cluster = g, type = "HC0", cadjust = TRUE)
cat("sandwich SEs (WS analytic):", round(sqrt(diag(r$Sigma_sand)), 4), "\n")
cat("vcovCL SEs               :", round(sqrt(diag(Vcl)), 4), "\n")
cat("adjusted draws SEs       :", round(apply(r$draws, 2, sd), 4), "\n")
cat("posterior SEs (unadj.)   :", round(sqrt(diag(vcov(m))), 4), "\n")
stopifnot(max(abs(sqrt(diag(r$Sigma_sand)) / sqrt(diag(Vcl)) - 1)) < 0.02,          # analytic sandwich matches vcovCL (2%)
          max(abs(apply(r$draws, 2, sd) / sqrt(diag(Vcl)) - 1)) < 0.03,             # rescaled draws carry it (Monte Carlo 3%)
          max(abs(colMeans(r$draws) - colMeans(dr))) < 1e-10)                        # mean unchanged
cat("TEST_AME_WS_OK\n")
