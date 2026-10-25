## ============================================================================
## simulate.R -- generate synthetic data from the v2 two-stage spatial choice model
##
## simulate_targeting(family = "M1b" | "M2b", S = 40, Tt = 25, K = 2, seed = 123, mu0 = c(0, 0))
## returns list(sy, alts, ev, truth, types)
##   mu0   : length-2 intercept of mu (population-average ideal gap), mu_i = mu0 + X beta + a_i; default 0 keeps
##           the planted types centred as described below, so mu_state at X = 0 is mu0 + a_i (truth$mu_state_ref0)
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
## family "M1m": one-dimensional truth -- the stage-2 utility uses only dM with scalar mu_i = X beta[,1] + a_i1
## and scalar variance sigma_x[1]^2 (beta[,2], a_i2, rho, dR play no role); used by run_sim_1d.R.
simulate_targeting <- function(family = c("M1b", "M2b", "M1m"), S = 40L, Tt = 25L, K = 2L,
                               seed = 123L, target_rate = 0.45, mu0 = c(0, 0)) {
  family <- match.arg(family)
  stopifnot(is.numeric(mu0), length(mu0) == 2, all(is.finite(mu0)))
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
  Mu <- rep(1, N1) %o% mu0 + Xm %*% beta + a_i[sy$state_idx, ]       # mu_n = mu0 + X_n beta + a_i
  ## ---- choice sets and systematic utility ------------------------------------------------
  alts <- sy |> mutate(sy_row = row_number()) |> select(sy_row, state_idx, year) |>
    tidyr::crossing(target_idx = seq_len(S)) |> filter(target_idx != state_idx) |>
    arrange(sy_row, target_idx)
  alts$dM <- M[cbind(alts$state_idx, alts$year - 1)] - M[cbind(alts$target_idx, alts$year - 1)]
  alts$dR <- R[cbind(alts$state_idx, alts$year - 1)] - R[cbind(alts$target_idx, alts$year - 1)]
  alts$Z1 <- contig[cbind(alts$state_idx, alts$target_idx)]
  alts$Z2 <- ally[cbind(alts$state_idx, alts$target_idx)]
  d1 <- alts$dM - Mu[alts$sy_row, 1]; d2 <- alts$dR - Mu[alts$sy_row, 2]
  quad <- if (family == "M1m") d1^2 / sigma_x[1]^2 else Sinv[1, 1] * d1^2 + 2 * Sinv[1, 2] * d1 * d2 + Sinv[2, 2] * d2^2
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
  mu_state_ref0 <- rep(1, S) %o% mu0 + a_i                        # mu_i at X = 0 is mu0 + a_i
  truth <- list(alpha0 = alpha0, alpha = alpha, tau = tau, mu0 = mu0, beta = beta, sigma_x = sigma_x, rho = rho,
                Sigma = Sigma, gamma = gamma, sigma_b = sigma_b, b_j = b_j, c_i = c_i, a_i = a_i,
                Omega = Omega_true, sd_re = sqrt(diag(Omega_true)), R_Om = cov2cor(Omega_true),
                U = U, V = V, UV = if (!is.null(U)) U %*% t(V) else NULL, family = family, K = K,
                mu_state_ref0 = mu_state_ref0)
  ## one-axis truth labels for the 1-D family (material axis): weaker / peer / stronger with band 0.8
  types_1d <- ifelse(mu_state_ref0[, 1] > 0.8, "weaker", ifelse(mu_state_ref0[, 1] < -0.8, "stronger", "peer"))
  list(sy = sy, alts = alts, ev = ev, truth = truth, types = types, types_1d = types_1d,
       meta = list(S = S, Tt = Tt, N1 = N1, N2 = nrow(alts), E = nrow(ev), seed = seed))
}

## ============================================================================
## simulate_current() -- the CURRENT model (review 2026-10-25, B5): the simulated world has the layout of the real cells
##   positions   : M, R random walks in log-share units (step SD `step_sd`; 0.3 = the original world, 0.06 = the persistent-position
##                 variant), NOT standardised within year; the stage-2 gaps are raw gaps divided by one constant per axis
##                 (gap_scale(), as the build does), so one unit matches the real data's scaled gap axes
##   levels (D2) : M_lvl, R_lvl = within-year centred level at t-1 / gap constant
##   changes     : dM_i, dR_i = (raw change over [t-2, t-1] minus its year mean) / pooled SD (h = 1), the target's dM_j_h, dRpr_j_h
##                 from the same table (every state has levels, so dch_j_missing is constant and absent)
##   W (minimal) : dM_i, dR_i, M_lvl, R_lvl, ongoing_tm1 (0/1), log_prior_mids, dem_tm1 (0/1), log_n_system_tm1 (D3; a centred
##                 year trend here, the system has a fixed size)                                         P1 = 8
##   X           : dM_i, dR_i                                                                              P  = 2
##   Z           : contig_land, joint_dem, log_dy_prior_mids, dy_mid_tm1, dM_j_h, dRpr_j_h                 Pz = 6
##   mu          : mu0 + X beta + a_i with a NON-ZERO mu0 (default (0.4, -0.3)); planted type centres as before
##   choice sets : a share `masked_share` of state-years has one stage-2-only row (in_iv = 0; the IV is summed over the in_iv rows)
##                 and, when the state-year initiates, one IV-only row (in_s2 = 0); non-initiators keep identical sets (their
##                 stage-2-only row leaves the block, as in the real build)
##   events      : one target per initiating state-year, plus for a share `coalition_share` of them a second target drawn from the
##                 same softmax (a "dispute" against two Side-B members): both events get weight 1/2 (D8), the rest weight 1
##   x_ref       : mean X over the initiating state-years; r0 = 0.5 (the real-data defaults); truth types are the deterministic
##                 classification of the TRUE mu_i at x_ref with that r0
## family "M1m": stage-2 utility on dM only (scalar mu_i = mu0[1] + X beta[, 1] + a_i1). Returns list(sy, alts, ev, truth, types,
## types_1d, levels, cols, x_ref, meta).
## ============================================================================
simulate_current <- function(family = c("M1b", "M2b", "M1m"), S = 40L, Tt = 25L, K = 2L, seed = 123L, target_rate = 0.45,
                             mu0 = c(0.4, -0.3), step_sd = 0.3, masked_share = 0.15, coalition_share = 0.2, r0 = 0.5) {
  family <- match.arg(family)
  stopifnot(is.numeric(mu0), length(mu0) == 2, all(is.finite(mu0)), step_sd > 0, masked_share >= 0, masked_share <= 1, coalition_share >= 0, coalition_share <= 1)
  set.seed(seed)
  ## ---- types and random effects (as in simulate_targeting) ------------------------------------------------------------------
  shares <- c(bully = 0.25, peer = 0.30, challenger = 0.25, hybrid_pm = 0.10, hybrid_mp = 0.10)
  n_type <- floor(shares * S); n_type["peer"] <- n_type["peer"] + (S - sum(n_type))
  planted <- sample(rep(names(n_type), n_type))
  centres <- rbind(bully = c(1.6, 1.6), peer = c(0, 0), challenger = c(-1.6, -1.6), hybrid_pm = c(1.6, -1.6), hybrid_mp = c(-1.6, 1.6))
  a_i <- centres[planted, ] + matrix(rnorm(S * 2, 0, 0.35), S, 2)
  c_i <- 0.45 * (a_i[, 1] + a_i[, 2]) / 2 + rnorm(S, 0, 0.5)
  Omega_true <- cov(cbind(c = c_i, a1 = a_i[, 1], a2 = a_i[, 2]))
  ## ---- positions: random walks in log-share units; NO within-year standardisation ------------------------------------------
  T_all <- Tt + 2
  M <- matrix(0, S, T_all); R <- matrix(0, S, T_all)
  M[, 1] <- rnorm(S, 0, 2.3); R[, 1] <- 0.6 * M[, 1] / 2.3 * 0.7 + rnorm(S, 0, 0.5)        # spreads like log CINC share (2.3) and log PageRank (0.7)
  for (t in 2:T_all) { M[, t] <- M[, t - 1] + rnorm(S, 0, step_sd); R[, t] <- R[, t - 1] + rnorm(S, 0, step_sd * 0.5) }
  years <- 3:T_all
  pos <- data.frame(state_idx = rep(seq_len(S), T_all), year = rep(seq_len(T_all), each = S), M = as.numeric(M), R = as.numeric(R))
  ## gap constants over the position years the state-years use (t - 1 for t in years, i.e. 2 .. T_all - 1); the real build uses every
  ## position year, which differs only by the unused first and last years
  pos_used <- pos |> filter(year >= min(years) - 1L, year <= max(years) - 1L)
  gap_M <- gap_scale(pos_used$M, pos_used$year); gap_R <- gap_scale(pos_used$R, pos_used$year)
  lvl <- pos |> group_by(year) |> mutate(M_lvl = (M - mean(M)) / gap_M, R_lvl = (R - mean(R)) / gap_R) |> ungroup()
  ## ---- state-year panel ----------------------------------------------------------------------------------------------------
  sy <- expand.grid(state_idx = seq_len(S), year = years) |> arrange(year, state_idx)
  N1 <- nrow(sy)
  raw_M <- M[cbind(sy$state_idx, sy$year - 1)] - M[cbind(sy$state_idx, sy$year - 2)]
  raw_R <- R[cbind(sy$state_idx, sy$year - 1)] - R[cbind(sy$state_idx, sy$year - 2)]
  cen_M <- raw_M - ave(raw_M, sy$year); cen_R <- raw_R - ave(raw_R, sy$year)
  sd_M <- sd(cen_M); sd_R <- sd(cen_R)
  sy$dM_raw_h <- raw_M; sy$dRpr_raw_h <- raw_R
  sy$dM_i <- cen_M / sd_M; sy$dR_i <- cen_R / sd_R
  sy$M_tm1 <- M[cbind(sy$state_idx, sy$year - 1)]; sy$Rpr_tm1 <- R[cbind(sy$state_idx, sy$year - 1)]
  sy$M_lvl <- lvl$M_lvl[match(paste(sy$state_idx, sy$year - 1), paste(lvl$state_idx, lvl$year))]
  sy$R_lvl <- lvl$R_lvl[match(paste(sy$state_idx, sy$year - 1), paste(lvl$state_idx, lvl$year))]
  dem <- rbinom(S, 1, 0.4)
  sy$ongoing_tm1 <- rbinom(N1, 1, 0.3); sy$dem_tm1 <- dem[sy$state_idx]; sy$log_prior_mids <- log1p(rpois(N1, 1.2))
  sy$log_n_system_tm1 <- log(30 + sy$year)                                               # a year trend (D3 column); centred below like the real controls
  ctrl <- c("ongoing_tm1", "log_prior_mids", "dem_tm1", "log_n_system_tm1")
  for (cc in ctrl) sy[[cc]] <- sy[[cc]] - mean(sy[[cc]])
  W_cols <- c("dM_i", "dR_i", "M_lvl", "R_lvl", ctrl); X_cols <- c("dM_i", "dR_i")
  Z_cols <- c("contig_land", "joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "dM_j_h", "dRpr_j_h")
  ## ---- structural parameters -------------------------------------------------------------------------------------------------
  beta <- rbind(c(0.5, 0.3), c(0.2, 0.6))                              # X = dM_i, dR_i
  sigma_x <- c(0.8, 0.8); rho <- 0.5
  Sigma <- matrix(c(sigma_x[1]^2, rho * prod(sigma_x), rho * prod(sigma_x), sigma_x[2]^2), 2); Sinv <- solve(Sigma)
  gamma <- c(0.8, 0.3, 0.4, -0.5, -0.2, 0.2)                           # Z
  sigma_b <- 0.5; b_j <- rnorm(S, 0, sigma_b); b_j <- b_j - mean(b_j)
  alpha <- c(0.3, -0.2, 0.4, 0.2, 0.5, 0.3, -0.3, 0.3)                 # W (dM_i, dR_i, M_lvl, R_lvl, ongoing, log_prior_mids, dem, log_n_system)
  tau <- 0.5
  U <- NULL; V <- NULL
  if (family == "M2b") { U <- matrix(rnorm(S * K, 0, 0.7), S, K); V <- matrix(rnorm(S * K, 0, 0.7), S, K); V <- sweep(V, 2, colMeans(V)) }
  Xm <- as.matrix(sy[, X_cols]); Mu <- rep(1, N1) %o% mu0 + Xm %*% beta + a_i[sy$state_idx, ]
  ## ---- alternatives -------------------------------------------------------------------------------------------------------------
  contig <- matrix(rbinom(S * S, 1, 0.15), S, S); contig[lower.tri(contig)] <- t(contig)[lower.tri(contig)]
  alts <- sy |> mutate(sy_row = row_number()) |> select(sy_row, state_idx, year) |> tidyr::crossing(target_idx = seq_len(S)) |>
    filter(target_idx != state_idx) |> arrange(sy_row, target_idx)
  alts$dM_logratio <- M[cbind(alts$state_idx, alts$year - 1)] - M[cbind(alts$target_idx, alts$year - 1)]
  alts$dR_pr_logratio <- R[cbind(alts$state_idx, alts$year - 1)] - R[cbind(alts$target_idx, alts$year - 1)]
  alts$dM <- alts$dM_logratio / gap_M; alts$dR <- alts$dR_pr_logratio / gap_R
  alts$contig_land <- contig[cbind(alts$state_idx, alts$target_idx)]
  alts$joint_dem <- dem[alts$state_idx] * dem[alts$target_idx]
  alts$log_dy_prior_mids <- log1p(rpois(nrow(alts), 0.3)); alts$dy_mid_tm1 <- rbinom(nrow(alts), 1, 0.05)
  alts$dM_j_h <- sy$dM_i[match(paste(alts$target_idx, alts$year), paste(sy$state_idx, sy$year))]
  alts$dRpr_j_h <- sy$dR_i[match(paste(alts$target_idx, alts$year), paste(sy$state_idx, sy$year))]
  stopifnot(!anyNA(alts$dM_j_h), !anyNA(alts$dRpr_j_h))
  ## masks (in_iv / in_s2): a share of state-years gets one stage-2-only row now; the IV-only row is set for the initiators below
  alts$in_iv <- 1L; alts$in_s2 <- 1L
  masked_sy <- which(runif(N1) < masked_share)
  s2only_row <- integer(0); ivonly_cand <- integer(0)
  for (n in masked_sy) { rows <- which(alts$sy_row == n); pick <- sample(rows, 2); alts$in_iv[pick[1]] <- 0L; s2only_row <- c(s2only_row, pick[1]); ivonly_cand <- c(ivonly_cand, pick[2]) }
  ## the continuous Z columns are centred over the IV rows (as assemble_cell() does)
  for (cc in c("log_dy_prior_mids", "dM_j_h", "dRpr_j_h")) alts[[cc]] <- alts[[cc]] - mean(alts[[cc]][alts$in_iv == 1L])
  d1 <- alts$dM - Mu[alts$sy_row, 1]; d2 <- alts$dR - Mu[alts$sy_row, 2]
  quad <- if (family == "M1m") d1^2 / sigma_x[1]^2 else Sinv[1, 1] * d1^2 + 2 * Sinv[1, 2] * d1 * d2 + Sinv[2, 2] * d2^2
  util <- -0.5 * quad + b_j[alts$target_idx] + as.numeric(as.matrix(alts[, Z_cols]) %*% gamma)
  if (family == "M2b") util <- util + rowSums(U[alts$state_idx, , drop = FALSE] * V[alts$target_idx, , drop = FALSE])
  alts$util <- util
  u_iv <- ifelse(alts$in_iv == 1L, util, -Inf)
  IV <- as.numeric(tapply(u_iv, alts$sy_row, function(u) log(sum(exp(u)))))
  stopifnot(length(IV) == N1, all(is.finite(IV)))
  ## ---- stage 1 ---------------------------------------------------------------------------------------------------------------
  lin <- as.numeric(as.matrix(sy[, W_cols]) %*% alpha) + c_i[sy$state_idx] + tau * IV
  alpha0 <- uniroot(function(a0) mean(plogis(a0 + lin)) - target_rate, c(-30, 30))$root
  sy$y1 <- rbinom(N1, 1, plogis(alpha0 + lin)); sy$IV_true <- IV
  ## IV-only rows for the masked INITIATORS; the stage-2-only rows of masked non-initiators leave the block (both sets identical)
  init_masked <- masked_sy[sy$y1[masked_sy] == 1L]
  alts$in_s2[ivonly_cand[match(init_masked, masked_sy)]] <- 0L
  drop_rows <- s2only_row[sy$y1[masked_sy] == 0L]
  if (length(drop_rows)) alts <- alts[-drop_rows, ]
  alts <- alts |> arrange(sy_row, target_idx)
  stopifnot(all(alts$in_iv == 1L | alts$in_s2 == 1L))
  ## ---- stage 2: one target per initiating state-year, a second one for a share of them (coalition, weight 1/2 each) --------------
  ev_rows <- which(sy$y1 == 1L)
  alt_index <- split(seq_len(nrow(alts)), alts$sy_row)
  ev_list <- list()
  for (n in ev_rows) {
    rows <- alt_index[[as.character(n)]]; rows <- rows[alts$in_s2[rows] == 1L]
    stopifnot(length(rows) >= 3)
    p <- exp(alts$util[rows] - max(alts$util[rows])); first <- rows[sample.int(length(rows), 1, prob = p)]
    if (runif(1) < coalition_share) {
      rest <- setdiff(rows, first); p2 <- exp(alts$util[rest] - max(alts$util[rest])); second <- rest[sample.int(length(rest), 1, prob = p2)]
      ev_list[[length(ev_list) + 1]] <- data.frame(sy_row = n, alt_row = c(first, second), w = 0.5)
    } else ev_list[[length(ev_list) + 1]] <- data.frame(sy_row = n, alt_row = first, w = 1)
  }
  ev <- bind_rows(ev_list)
  stopifnot(all(alts$sy_row[ev$alt_row] == ev$sy_row), all(alts$in_s2[ev$alt_row] == 1L), !anyDuplicated(paste(ev$sy_row, ev$alt_row)))
  ## ---- truth, x_ref, types ------------------------------------------------------------------------------------------------------
  x_ref <- colMeans(Xm[sy$y1 == 1L, , drop = FALSE])
  mu_state_ref <- rep(1, S) %o% mu0 + rep(1, S) %o% as.numeric(x_ref %*% beta) + a_i
  r <- sqrt(rowSums(mu_state_ref^2))
  types <- ifelse(r < r0, "peer", ifelse(mu_state_ref[, 1] > 0 & mu_state_ref[, 2] > 0, "bully", ifelse(mu_state_ref[, 1] < 0 & mu_state_ref[, 2] < 0, "challenger",
                  ifelse(mu_state_ref[, 1] > 0, "hybrid_pm", "hybrid_mp"))))
  types_1d <- ifelse(mu_state_ref[, 1] > r0, "weaker", ifelse(mu_state_ref[, 1] < -r0, "stronger", "peer"))
  truth <- list(alpha0 = alpha0, alpha = alpha, tau = tau, mu0 = mu0, beta = beta, sigma_x = sigma_x, rho = rho, Sigma = Sigma, gamma = gamma,
                sigma_b = sigma_b, b_j = b_j, c_i = c_i, a_i = a_i, Omega = Omega_true, sd_re = sqrt(diag(Omega_true)), R_Om = cov2cor(Omega_true),
                U = U, V = V, UV = if (!is.null(U)) U %*% t(V) else NULL, family = family, K = K, mu_state_ref = mu_state_ref, planted = planted,
                gap_M = gap_M, gap_R = gap_R, sd_change_M = sd_M, sd_change_R = sd_R, step_sd = step_sd, r0 = r0)
  list(sy = sy, alts = alts, ev = ev, truth = truth, types = types, types_1d = types_1d, levels = pos, x_ref = x_ref,
       cols = list(W = W_cols, X = X_cols, Z = Z_cols),
       meta = list(S = S, Tt = Tt, N1 = N1, N2 = nrow(alts), E = nrow(ev), seed = seed, n_masked_sy = length(masked_sy), n_coalition = sum(ev$w < 1) / 2,
                   n_iv_only = sum(alts$in_s2 == 0L), n_s2_only = sum(alts$in_iv == 0L)))
}
