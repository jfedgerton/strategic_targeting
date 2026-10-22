## ============================================================================
## simulate.R -- generate synthetic data from the v2 two-stage spatial choice model
##
## simulate_targeting(family = "M1b" | "M2b", S = 40, Tt = 25, K = 2, seed = 123)
## returns list(sy, alts, ev, truth, types)
##   sy    : state-year panel (y1, state_idx, year, W1, W2, X1..X4)
##   alts  : flattened choice sets, one row per (state-year, alternative)
##   ev    : one event per initiating state-year (sy_row, alt_row)
##   truth : every generating parameter (for recovery checks)
##   types : planted type per state (bully / peer / challenger / hybrid_pm / hybrid_mp)
##
## Planted types are offsets a_i (type centre + small noise) in the (dM, dR) plane,
## so at the reference X = 0 the initiator's mu_i sits in the planted region:
##   bully      (+,+) far        centre ( 1.6,  1.6)
##   peer       near origin      centre ( 0.0,  0.0)
##   challenger (-,-)            centre (-1.6, -1.6)
##   hybrid_pm  (+,-)            centre ( 1.6, -1.6)
##   hybrid_mp  (-,+)            centre (-1.6,  1.6)
## c_i is generated as a linear function of a_i plus noise, so Omega has non-zero
## correlations. Because a_i is a type mixture (not MVN), the "true Omega" reported is
## the realised sample covariance of (c_i, a_i1, a_i2) across states.
## ============================================================================
simulate_targeting <- function(family = c("M1b", "M2b"), S = 40L, Tt = 25L, K = 2L,
                               seed = 123L, target_rate = 0.45) {
  family <- match.arg(family)
  set.seed(seed)
  ## ---- types ---------------------------------------------------------------------
  type_levels <- c("bully", "peer", "challenger", "hybrid_pm", "hybrid_mp")
  ## type shares (bully .25, peer .30, challenger .25, hybrids .10 each); S = 40 gives 10/12/10/4/4
  shares <- c(bully = 0.25, peer = 0.30, challenger = 0.25, hybrid_pm = 0.10, hybrid_mp = 0.10)
  n_type <- floor(shares * S)
  n_type["peer"] <- n_type["peer"] + (S - sum(n_type))
  stopifnot(sum(n_type) == S, all(n_type >= 0))
  types <- rep(names(n_type), n_type)
  types <- sample(types)                               # random assignment to state indices
  centres <- rbind(bully = c(1.6, 1.6), peer = c(0, 0), challenger = c(-1.6, -1.6),
                   hybrid_pm = c(1.6, -1.6), hybrid_mp = c(-1.6, 1.6))
  a_i <- centres[types, ] + matrix(rnorm(S * 2, 0, 0.35), S, 2)
  ## c_i correlated with a_i: aggressive (bully-like) states initiate more
  c_i <- 0.45 * (a_i[, 1] + a_i[, 2]) / 2 + rnorm(S, 0, 0.5)
  re <- cbind(c = c_i, a1 = a_i[, 1], a2 = a_i[, 2])
  Omega_true <- cov(re)
  ## ---- positions M, R: random walks standardised within year ----------------------
  T_all <- Tt + 2                                       # two burn-in years for the change window
  M <- matrix(0, S, T_all); R <- matrix(0, S, T_all)
  M[, 1] <- rnorm(S); R[, 1] <- 0.6 * M[, 1] + rnorm(S, 0, 0.8)
  for (t in 2:T_all) {
    M[, t] <- M[, t - 1] + rnorm(S, 0, 0.3)
    R[, t] <- R[, t - 1] + rnorm(S, 0, 0.3)
  }
  M <- apply(M, 2, function(v) as.numeric(scale(v)))
  R <- apply(R, 2, function(v) as.numeric(scale(v)))
  ## ---- structural parameters -----------------------------------------------------------
  beta <- rbind(c(0.5, 0.3),      # X1 = dM_i (own material change, h = 1)
                c(0.2, 0.6),      # X2 = dR_i (own relational change)
                c(-0.3, 0.1),     # X3 = control (continuous)
                c(0.1, -0.2))     # X4 = control (continuous)
  sigma_x <- c(0.8, 0.8); rho <- 0.5
  Sigma <- matrix(c(sigma_x[1]^2, rho * prod(sigma_x), rho * prod(sigma_x), sigma_x[2]^2), 2)
  Sinv <- solve(Sigma)
  gamma <- c(0.8, -0.5)          # Z1 = contiguity dummy, Z2 = alliance dummy
  sigma_b <- 0.5
  b_j <- rnorm(S, 0, sigma_b); b_j <- b_j - mean(b_j)
  alpha <- c(0.6, -0.3)          # W1 = major power dummy, W2 = continuous control
  tau <- 0.5
  U <- NULL; V <- NULL
  if (family == "M2b") {
    U <- matrix(rnorm(S * K, 0, 0.7), S, K)
    V <- matrix(rnorm(S * K, 0, 0.7), S, K)
    V <- sweep(V, 2, colMeans(V))
  }
  ## ---- time-invariant dyadic controls (symmetric) ---------------------------------------
  contig <- matrix(rbinom(S * S, 1, 0.15), S, S); contig[lower.tri(contig)] <- t(contig)[lower.tri(contig)]
  ally   <- matrix(rbinom(S * S, 1, 0.20), S, S); ally[lower.tri(ally)] <- t(ally)[lower.tri(ally)]
  major  <- as.integer(seq_len(S) %in% order(-rowMeans(M))[1:6])
  ## ---- state-year panel ---------------------------------------------------------------------
  years <- 3:T_all                                       # usable years (need t-1 and t-2)
  sy <- expand.grid(state_idx = seq_len(S), year = years) |> arrange(year, state_idx)
  sy$X1 <- M[cbind(sy$state_idx, sy$year - 1)] - M[cbind(sy$state_idx, sy$year - 2)]
  sy$X2 <- R[cbind(sy$state_idx, sy$year - 1)] - R[cbind(sy$state_idx, sy$year - 2)]
  sy$X3 <- rnorm(nrow(sy)); sy$X4 <- rnorm(nrow(sy))
  sy$W1 <- major[sy$state_idx]; sy$W2 <- rnorm(nrow(sy))
  N1 <- nrow(sy)
  Xm <- as.matrix(sy[, c("X1", "X2", "X3", "X4")])
  Mu <- Xm %*% beta + a_i[sy$state_idx, ]
  ## ---- choice sets and systematic utility ------------------------------------------------
  alts <- sy |> mutate(sy_row = row_number()) |> select(sy_row, state_idx, year) |>
    tidyr::crossing(target_idx = seq_len(S)) |> filter(target_idx != state_idx) |>
    arrange(sy_row, target_idx)
  alts$dM <- M[cbind(alts$state_idx, alts$year - 1)] - M[cbind(alts$target_idx, alts$year - 1)]
  alts$dR <- R[cbind(alts$state_idx, alts$year - 1)] - R[cbind(alts$target_idx, alts$year - 1)]
  alts$Z1 <- contig[cbind(alts$state_idx, alts$target_idx)]
  alts$Z2 <- ally[cbind(alts$state_idx, alts$target_idx)]
  d1 <- alts$dM - Mu[alts$sy_row, 1]; d2 <- alts$dR - Mu[alts$sy_row, 2]
  quad <- Sinv[1, 1] * d1^2 + 2 * Sinv[1, 2] * d1 * d2 + Sinv[2, 2] * d2^2
  util <- -0.5 * quad + b_j[alts$target_idx] + alts$Z1 * gamma[1] + alts$Z2 * gamma[2]
  if (family == "M2b") util <- util + rowSums(U[alts$state_idx, , drop = FALSE] * V[alts$target_idx, , drop = FALSE])
  alts$util <- util
  lse <- tapply(util, alts$sy_row, function(u) log(sum(exp(u))))
  stopifnot(length(lse) == N1)
  IV <- as.numeric(lse)
  ## ---- stage 1: calibrate alpha0 to the target initiation rate, then draw ------------------
  lin <- alpha[1] * sy$W1 + alpha[2] * sy$W2 + c_i[sy$state_idx] + tau * IV
  f <- function(a0) mean(plogis(a0 + lin)) - target_rate
  alpha0 <- uniroot(f, c(-30, 30))$root
  sy$y1 <- rbinom(N1, 1, plogis(alpha0 + lin))
  sy$IV_true <- IV
  ## ---- stage 2: one target per initiating state-year ------------------------------------------
  ev_rows <- which(sy$y1 == 1L)
  alt_index <- split(seq_len(nrow(alts)), alts$sy_row)
  chosen <- vapply(ev_rows, function(n) {
    rows <- alt_index[[n]]
    p <- exp(util[rows] - lse[n])
    rows[sample.int(length(rows), 1, prob = p)]
  }, integer(1))
  ev <- data.frame(sy_row = ev_rows, alt_row = chosen)
  stopifnot(all(alts$sy_row[ev$alt_row] == ev$sy_row))
  truth <- list(alpha0 = alpha0, alpha = alpha, tau = tau, beta = beta, sigma_x = sigma_x, rho = rho,
                Sigma = Sigma, gamma = gamma, sigma_b = sigma_b, b_j = b_j, c_i = c_i, a_i = a_i,
                Omega = Omega_true, sd_re = sqrt(diag(Omega_true)), R_Om = cov2cor(Omega_true),
                U = U, V = V, UV = if (!is.null(U)) U %*% t(V) else NULL, family = family, K = K,
                mu_state_ref0 = a_i)       # mu_i at X = 0 is exactly a_i
  list(sy = sy, alts = alts, ev = ev, truth = truth, types = types,
       meta = list(S = S, Tt = Tt, N1 = N1, N2 = nrow(alts), E = nrow(ev), seed = seed))
}
