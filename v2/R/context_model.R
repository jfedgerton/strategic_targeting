## ============================================================================
## context_model.R -- R-side functions of the M1c / M1cr models (dual typology stage b, PI-approved spec 2026-10-27; G groups as data and the
## role-specific-rho model M1cr since 2026-10-31)
## Sourced after spec.R, helpers_data.R, helpers_model.R and typology.R. Functions only, no side effects.
##   stage2_loglik_R_context(dat, pars)  brute-force log-likelihood of the M1c / M1cr programs (both stages' pieces) for the tests
##   context_model_summary(fit, groups)  the stage-b quantities of one M1c / M1cr fit as a long table: delta[g], sigma_g, rho (shared rho_x or
##                                       per-group rho_g), the group population ideal points mu0_g, Sigma_g per group
##   context_group_typology(fit, ...)    the typology per group: the existing classifier on mu_state + delta[g] (reference group: delta = 0)
## Group g (0-based ev_grp) has the label groups[g + 1]; groups defaults to STAGE2_GROUPS (G = 3: bilateral, leader, joiner).
## ============================================================================

## ---- brute-force M1c / M1cr likelihood in R (mirrors stage2_loglik_R of helpers_model.R) ----------------------------------------------
## pars: list(mu0 (2), delta ((G - 1) x 2 matrix: row g = the shift of group g), beta (P x 2), a_i (S x 2), b_j (S), gamma (Pz), sigma_g (G x 2
## matrix: row g + 1 = group g), rho (a scalar: the shared rho of M1c; or a length-G vector: rho_g of M1cr, element g + 1 = group g)).
## dat: build_stan_data() output plus G and ev_grp. Returns per-event unweighted log-lik, the reference-group IV (lse) per state-year, the
## weighted stage-2 sum and the event ideal points.
stage2_loglik_R_context <- function(dat, pars) {
  G <- dat$G
  stopifnot(dat$D_mode == 3L, dat$joint == 1L, dat$sampled == 0L, !is.null(G), G >= 2L, !is.null(dat$ev_grp), length(dat$ev_grp) == dat$E, all(dat$ev_grp %in% 0:(G - 1L)))
  for (nm in c("mu0", "delta", "beta", "a_i", "b_j", "gamma", "sigma_g", "rho")) if (is.null(pars[[nm]])) stop("stage2_loglik_R_context: pars$", nm, " is missing")
  delta <- matrix(pars$delta, ncol = 2)
  stopifnot(length(pars$mu0) == 2, nrow(delta) == G - 1L, identical(dim(pars$sigma_g), c(as.integer(G), 2L)), all(pars$sigma_g > 0),
            length(pars$rho) %in% c(1L, G), all(abs(pars$rho) < 1))
  rho <- rep(as.numeric(pars$rho), length.out = G)
  Mu <- rep(1, dat$N1) %o% as.numeric(pars$mu0) + dat$X %*% matrix(pars$beta, ncol = 2) + matrix(pars$a_i, ncol = 2)[dat$sy_state, , drop = FALSE]
  zg <- if (dat$Pz > 0) as.numeric(dat$Z %*% pars$gamma) else rep(0, dat$N2)
  shift <- rbind(c(0, 0), delta)
  Sinv <- lapply(seq_len(G), function(k) solve(matrix(c(pars$sigma_g[k, 1]^2, rho[k] * prod(pars$sigma_g[k, ]), rho[k] * prod(pars$sigma_g[k, ]), pars$sigma_g[k, 2]^2), 2)))
  util_g <- function(rows, n, g) {
    mu_e <- Mu[n, ] + shift[g + 1L, ]
    vapply(rows, function(r) { d <- dat$x[r, ] - mu_e; -0.5 * as.numeric(t(d) %*% Sinv[[g + 1L]] %*% d) + pars$b_j[dat$alt_target[r]] + zg[r] }, numeric(1))
  }
  lse <- numeric(dat$N1)
  for (n in seq_len(dat$N1)) {
    rows <- dat$cs_start[n]:dat$cs_end[n]
    rows_iv <- rows[dat$in_iv[rows] == 1L]
    stopifnot(length(rows_iv) >= 1)
    lse[n] <- log(sum(exp(util_g(rows_iv, n, 0L))))
  }
  loglik <- numeric(dat$E); mu_event <- matrix(NA_real_, dat$E, 2)
  for (e in seq_len(dat$E)) {
    n <- dat$ev_sy[e]; g <- dat$ev_grp[e]
    rows <- dat$cs_start[n]:dat$cs_end[n]; rows <- rows[dat$in_s2[rows] == 1L]
    stopifnot(dat$ev_row[e] %in% rows)
    u <- util_g(rows, n, g)
    loglik[e] <- u[match(dat$ev_row[e], rows)] - log(sum(exp(u)))
    mu_event[e, ] <- Mu[n, ] + shift[g + 1L, ]
  }
  ev_w <- if (is.null(dat$ev_w)) rep(1, dat$E) else dat$ev_w
  list(loglik = loglik, lse = lse, ev_w = ev_w, lp_stage2 = sum(loglik * ev_w), mu_event = mu_event, Mu = Mu)
}

## ---- stage-b quantities of one M1c / M1cr fit -------------------------------------------------------------------------------------------
## One long table (quantity, group, axis, mean, sd, lo, hi, Rhat, n_eff) over: mu0 (reference group), delta[g] (shift of group g), mu0_g (the
## group population ideal points), sigma_g per group and axis, the correlation (rho_x shared, or rho_g per group for M1cr), and the G
## Sigma_g (variances and the covariance). lo / hi are the 5 / 95% quantiles (a 90% interval, as elsewhere in the outputs).
context_model_summary <- function(fit, groups = STAGE2_GROUPS) {
  G <- length(groups)
  rho_par <- if ("rho_g" %in% fit@model_pars) "rho_g" else "rho_x"
  stopifnot(inherits(fit, "stanfit"), all(c("mu0", "delta", "sigma_g", rho_par, "mu0_g", "Sigma_g") %in% fit@model_pars),
            identical(as.integer(fit@par_dims$sigma_g), c(as.integer(G), 2L)), identical(as.integer(fit@par_dims$delta), c(as.integer(G) - 1L, 2L)))
  ax <- c("material", "relational")
  sm <- summary(fit, pars = c("mu0", "delta", "mu0_g", "sigma_g", rho_par, "Sigma_g"), probs = c(0.05, 0.95))$summary
  rn <- rownames(sm)
  base <- sub("\\[.*$", "", rn)
  i1 <- suppressWarnings(as.integer(sub("^[^\\[]+\\[([0-9]+).*$", "\\1", rn)))
  i2 <- suppressWarnings(as.integer(sub("^[^\\[]+\\[[0-9]+,([0-9]+).*$", "\\1", rn)))
  i3 <- suppressWarnings(as.integer(sub("^[^\\[]+\\[[0-9]+,[0-9]+,([0-9]+)\\]$", "\\1", rn)))
  grp <- rep(NA_character_, length(rn)); axis <- rep(NA_character_, length(rn)); qty <- base
  w <- base == "mu0"; grp[w] <- groups[1]; axis[w] <- ax[i1[w]]; qty[w] <- "mu0 (population ideal point)"
  w <- base == "mu0_g"; grp[w] <- groups[i1[w]]; axis[w] <- ax[i2[w]]; qty[w] <- "mu0_g (group population ideal point)"
  w <- base == "delta"; grp[w] <- groups[i1[w] + 1L]; axis[w] <- ax[i2[w]]; qty[w] <- sprintf("delta (shift vs %s)", groups[1])
  w <- base == "sigma_g"; grp[w] <- groups[i1[w]]; axis[w] <- ax[i2[w]]; qty[w] <- "sigma_g (tolerance scale)"
  w <- base == "rho_x"; grp[w] <- "shared"; axis[w] <- "material x relational"; qty[w] <- "rho (shared tolerance correlation)"
  w <- base == "rho_g"; grp[w] <- groups[i1[w]]; axis[w] <- "material x relational"; qty[w] <- "rho_g (group tolerance correlation)"
  w <- base == "Sigma_g"; grp[w] <- groups[i1[w]]; axis[w] <- paste(ax[i2[w]], ax[i3[w]], sep = " x "); qty[w] <- "Sigma_g"
  out <- data.frame(quantity = qty, group = grp, axis = axis, mean = sm[, "mean"], sd = sm[, "sd"], lo = sm[, "5%"], hi = sm[, "95%"],
                    Rhat = sm[, "Rhat"], n_eff = sm[, "n_eff"], parameter = rn, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  stopifnot(!anyNA(out$group), !anyNA(out$axis), nrow(out) == 2 + 2 * (G - 1) + 2 * G + 2 * G + (if (rho_par == "rho_g") G else 1) + 4 * G)
  ## the redundant Sigma_g[g, 2, 1] rows are dropped (symmetric)
  out <- out[!(base == "Sigma_g" & i2 == 2L & i3 == 1L), ]
  rownames(out) <- NULL
  out
}

## ---- typology per group: the existing classifier on mu_state + delta[g] -------------------------------------------------------------------
## mu_state (generated quantity) is the REFERENCE-group ideal point at x_ref; the ideal point of group g >= 1 of every state is mu_state +
## delta[g], added draw by draw. n_events_by_group: an S x G matrix of the state's usable events per group (columns in group order), so each
## table flags the states that never acted in that group (their classification in the group rests on a_i from the other groups and on delta[g]).
## Returns a named list (one per group label) of classify_types() tables with the extra columns group, n_events_in_group, n_events_as_initiator.
context_group_typology <- function(fit, state_labels, r0, n_events_by_group, groups = STAGE2_GROUPS) {
  G <- length(groups)
  stopifnot(inherits(fit, "stanfit"), is.matrix(n_events_by_group), ncol(n_events_by_group) == G, nrow(n_events_by_group) == length(state_labels))
  mu <- mu_state_draws(fit)                                       # [draws, S, 2]
  stopifnot(dim(mu)[3] == 2L, dim(mu)[2] == length(state_labels))
  dl <- rstan::extract(fit, pars = "delta")$delta                 # [draws, G - 1, 2]
  stopifnot(identical(dim(dl), c(dim(mu)[1], as.integer(G) - 1L, 2L)))
  out <- vector("list", G); names(out) <- groups
  for (g in seq_len(G)) {
    mu_g <- mu
    if (g > 1L) for (d in 1:2) mu_g[, , d] <- mu[, , d] + dl[, g - 1L, d]      # the same draw's shift added to every state
    cls <- classify_types(mu_g, r0 = r0, state_labels = state_labels)
    cls$group <- groups[g]
    cls$n_events_in_group <- as.integer(n_events_by_group[, g])
    cls$n_events_as_initiator <- as.integer(rowSums(n_events_by_group))
    cls$r0 <- r0
    out[[g]] <- cls
  }
  ## the reference table equals the plain classification of mu_state (shift 0), asserted
  stopifnot(identical(out[[1]]$map_type, classify_types(mu, r0 = r0, state_labels = state_labels)$map_type))
  out
}
