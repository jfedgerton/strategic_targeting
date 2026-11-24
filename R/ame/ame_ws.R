## R/ame/ame_ws.R -- Williams-Savitsky (2021, Bayesian Analysis) adjustment of a weighted pseudo-posterior (PI 2026-11-24, design section 6).
## The weighted likelihood makes the AME posterior a pseudo-posterior whose spread ignores dispute clustering. The adjustment keeps the
## posterior mean and rescales the draws of the FIXED parameters so their covariance equals the cluster-robust sandwich:
##   Sigma_post = cov(draws) (approximates H^-1);  J = G/(G-1) * sum_g s_g s_g'  (s_g = dispute g's score of the weighted log-likelihood,
##   prior excluded, at the posterior mean);  Sigma_sand = Sigma_post J Sigma_post;  draws_adj = mean + (draws - mean) U_post^-1 U_sand,
##   with U the upper Cholesky factors. Random effects (b, eta, del, eps, u, v) are held at their posterior means when the scores are taken.
## Functions only; the script is R/ame/ame_ws_adjust.R; the test is tests/test_ame_ws.R.

## numerical Jacobian of a vector-valued function (central differences): rows = clusters, columns = parameters
ws_jacobian <- function(f, p, h = 1e-5) {
  f0 <- f(p); Jm <- matrix(0, length(f0), length(p))
  for (k in seq_along(p)) { e <- rep(0, length(p)); hk <- h * max(1, abs(p[k])); e[k] <- hk; Jm[, k] <- (f(p + e) - f(p - e)) / (2 * hk) }
  Jm
}
## the rescaling itself; draws: S x P matrix of fixed-parameter draws; scores: G x P cluster scores at the posterior mean
ws_rescale <- function(draws, scores) {
  stopifnot(is.matrix(draws), is.matrix(scores), ncol(draws) == ncol(scores), nrow(scores) > ncol(scores))
  G <- nrow(scores); mu <- colMeans(draws); Sp <- cov(draws)
  Jc <- crossprod(scores) * G / (G - 1)
  Ss <- Sp %*% Jc %*% Sp; Ss <- (Ss + t(Ss)) / 2
  Up <- tryCatch(chol(Sp), error = function(e) stop("posterior covariance of the fixed parameters is not positive definite: ", conditionMessage(e)))
  Us <- tryCatch(chol(Ss), error = function(e) stop("sandwich covariance is not positive definite: ", conditionMessage(e)))
  A <- backsolve(Up, diag(ncol(draws))) %*% Us                       # U_post^-1 U_sand
  adj <- sweep(sweep(draws, 2, mu) %*% A, 2, mu, "+"); colnames(adj) <- colnames(draws)
  list(draws = adj, Sigma_post = Sp, Sigma_sand = Ss, ratio_sd = sqrt(diag(Ss) / diag(Sp)))
}

## per-dispute weighted log-likelihood of the joint AME (stan/ame_joint.stan) as a function of the fixed parameters
##   p = c(beta1, as.vector(alpha) [3 x NE_ERA, column-major], betaA, betaB, betaO, theta); latent = posterior means of b, eta, del, bbar,
##   eps (ND x 3), and for K = 1 u, v, lam1, lam2. Returns a named vector over the clusters (dispute ids).
ame_unpack <- function(p, d) { i <- 0; nx <- function(n) { r <- p[i + seq_len(n)]; i <<- i + n; r }
  list(beta1 = nx(d$P1), alpha = matrix(nx(3 * d$NE_ERA), 3, d$NE_ERA), betaA = nx(d$PA), betaB = nx(d$PB), betaO = nx(d$PO), theta = nx(3)) }
ame_ll_by_cluster <- function(p, d, lat, ev_disp, disp_ids) {
  q <- ame_unpack(p, d)
  ## stage 1 (only if fitted)
  ll1 <- numeric(0); c1 <- integer(0)
  if (d$use_s1) {
    ev <- rep(seq_len(d$E), d$e_end - d$e_start + 1L)
    u <- drop(d$X1 %*% q$beta1) + lat$b[d$tgt1] + lat$eta[cbind(d$tgt1, d$era1)] + d$coal1 * lat$del[d$tgt1]
    if (d$K) u <- u + lat$lam1 * lat$u[d$ini1] * lat$v[d$tgt1]
    m <- ave(u, ev, FUN = max); lse <- log(rowsum(exp(u - m), ev, reorder = TRUE))[, 1] + m[d$e_start]
    ll1 <- d$w_e * (u[d$e_chosen] - lse); c1 <- ev_disp }
  ## stage 2 (only if fitted)
  ll2 <- numeric(0); c2 <- integer(0)
  if (d$use_s2) {
    dd <- d$d2; ed <- d$era_d[dd]; la <- lb <- 0
    if (d$K) { vA <- drop(d$W_oppA %*% lat$v); vB <- drop(d$W_oppB %*% lat$v); la <- lat$lam2 * lat$u[d$k2] * vA[dd]; lb <- lat$lam2 * lat$u[d$k2] * vB[dd] }
    UA <- drop(d$XA %*% q$betaA) + q$alpha[cbind(1, ed)] + q$theta[1] * lat$bbar[dd] + lat$eps[cbind(dd, 1)] + la
    UB <- drop(d$XB %*% q$betaB) + q$alpha[cbind(2, ed)] + q$theta[2] * lat$bbar[dd] + lat$eps[cbind(dd, 2)] + lb
    UO <- drop(d$XO %*% q$betaO) + q$alpha[cbind(3, ed)] + q$theta[3] * lat$bbar[dd] + lat$eps[cbind(dd, 3)] + la
    Um <- cbind(0, ifelse(d$av2[, 1] == 1, UA, -Inf), ifelse(d$av2[, 2] == 1, UB, -Inf), ifelse(d$av2[, 3] == 1, UO, -Inf))
    mx <- apply(Um, 1, max); lse <- mx + log(rowSums(exp(Um - mx)))
    ll2 <- Um[cbind(seq_len(nrow(Um)), d$y2)] - lse; c2 <- disp_ids[dd] }
  s <- rowsum(c(ll1, ll2), c(c1, c2)); setNames(s[, 1], rownames(s))
}
