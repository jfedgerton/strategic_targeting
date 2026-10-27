## ============================================================================
## context_model.R -- R-side functions of the M1c model (dual typology stage b, PI-approved spec 2026-10-27)
## Sourced after spec.R, helpers_data.R, helpers_model.R and typology.R. Functions only, no side effects.
##   stage2_loglik_R_context(dat, pars)  brute-force log-likelihood of the M1c program (both stages' pieces) for the tests
##   context_model_summary(fit)          the stage-b quantities of one M1c fit as a long table: delta_lead, delta_join, sigma_g, rho,
##                                       the group population ideal points mu0 / mu0 + delta_lead / mu0 + delta_join, Sigma_g per group
##   context_group_typology(fit, ...)    the typology per group: the existing classifier on mu_state + delta_g (bilateral: delta = 0)
## ============================================================================

## ---- brute-force M1c likelihood in R (mirrors stage2_loglik_R of helpers_model.R) -----------------------------------------------------
## pars: list(mu0 (2), delta_lead (2), delta_join (2), beta (P x 2), a_i (S x 2), b_j (S), gamma (Pz), sigma_g (3 x 2 matrix: rows bilateral,
## leader, joiner), rho). dat: build_stan_data() output plus ev_grp. Returns per-event unweighted log-lik, the bilateral IV (lse) per
## state-year, the weighted stage-2 sum and the event ideal points.
stage2_loglik_R_context <- function(dat, pars) {
  stopifnot(dat$D_mode == 3L, dat$joint == 1L, dat$sampled == 0L, !is.null(dat$ev_grp), length(dat$ev_grp) == dat$E, all(dat$ev_grp %in% 0:2))
  for (nm in c("mu0", "delta_lead", "delta_join", "beta", "a_i", "b_j", "gamma", "sigma_g", "rho")) if (is.null(pars[[nm]])) stop("stage2_loglik_R_context: pars$", nm, " is missing")
  stopifnot(length(pars$mu0) == 2, length(pars$delta_lead) == 2, length(pars$delta_join) == 2, identical(dim(pars$sigma_g), c(3L, 2L)), all(pars$sigma_g > 0), abs(pars$rho) < 1)
  Mu <- rep(1, dat$N1) %o% as.numeric(pars$mu0) + dat$X %*% matrix(pars$beta, ncol = 2) + matrix(pars$a_i, ncol = 2)[dat$sy_state, , drop = FALSE]
  zg <- if (dat$Pz > 0) as.numeric(dat$Z %*% pars$gamma) else rep(0, dat$N2)
  shift <- rbind(c(0, 0), as.numeric(pars$delta_lead), as.numeric(pars$delta_join))
  Sinv <- lapply(1:3, function(g) solve(matrix(c(pars$sigma_g[g, 1]^2, pars$rho * prod(pars$sigma_g[g, ]), pars$rho * prod(pars$sigma_g[g, ]), pars$sigma_g[g, 2]^2), 2)))
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

## ---- stage-b quantities of one M1c fit ------------------------------------------------------------------------------------------------
## One long table (quantity, group, axis, mean, sd, lo, hi, Rhat, n_eff) over: mu0 (bilateral), delta_lead, delta_join, mu0_lead, mu0_join
## (the group population ideal points), sigma_g per group and axis, rho_x, and the three Sigma_g (variances and the covariance). lo / hi are
## the 5 / 95% quantiles (a 90% interval, as elsewhere in the outputs).
context_model_summary <- function(fit) {
  stopifnot(inherits(fit, "stanfit"), all(c("mu0", "delta_lead", "delta_join", "sigma_g", "rho_x", "mu0_lead", "mu0_join", "Sigma_g") %in% fit@model_pars))
  ax <- c("material", "relational")
  sm <- summary(fit, pars = c("mu0", "delta_lead", "delta_join", "mu0_lead", "mu0_join", "sigma_g", "rho_x", "Sigma_g"), probs = c(0.05, 0.95))$summary
  rn <- rownames(sm)
  base <- sub("\\[.*$", "", rn)
  i1 <- suppressWarnings(as.integer(sub("^[^\\[]+\\[([0-9]+).*$", "\\1", rn)))
  i2 <- suppressWarnings(as.integer(sub("^[^\\[]+\\[[0-9]+,([0-9]+).*$", "\\1", rn)))
  i3 <- suppressWarnings(as.integer(sub("^[^\\[]+\\[[0-9]+,[0-9]+,([0-9]+)\\]$", "\\1", rn)))
  grp <- rep(NA_character_, length(rn)); axis <- rep(NA_character_, length(rn)); qty <- base
  w <- base == "mu0"; grp[w] <- "bilateral"; axis[w] <- ax[i1[w]]; qty[w] <- "mu0 (population ideal point)"
  w <- base == "mu0_lead"; grp[w] <- "leader"; axis[w] <- ax[i1[w]]; qty[w] <- "mu0 + delta_lead (population ideal point)"
  w <- base == "mu0_join"; grp[w] <- "joiner"; axis[w] <- ax[i1[w]]; qty[w] <- "mu0 + delta_join (population ideal point)"
  w <- base == "delta_lead"; grp[w] <- "leader"; axis[w] <- ax[i1[w]]; qty[w] <- "delta_lead (shift vs bilateral)"
  w <- base == "delta_join"; grp[w] <- "joiner"; axis[w] <- ax[i1[w]]; qty[w] <- "delta_join (shift vs bilateral)"
  w <- base == "sigma_g"; grp[w] <- STAGE2_GROUPS[i1[w]]; axis[w] <- ax[i2[w]]; qty[w] <- "sigma_g (tolerance scale)"
  w <- base == "rho_x"; grp[w] <- "shared"; axis[w] <- "material x relational"; qty[w] <- "rho (shared tolerance correlation)"
  w <- base == "Sigma_g"; grp[w] <- STAGE2_GROUPS[i1[w]]; axis[w] <- paste(ax[i2[w]], ax[i3[w]], sep = " x "); qty[w] <- "Sigma_g"
  out <- data.frame(quantity = qty, group = grp, axis = axis, mean = sm[, "mean"], sd = sm[, "sd"], lo = sm[, "5%"], hi = sm[, "95%"],
                    Rhat = sm[, "Rhat"], n_eff = sm[, "n_eff"], parameter = rn, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  stopifnot(!anyNA(out$group), !anyNA(out$axis), nrow(out) == 2 + 2 + 2 + 2 + 2 + 6 + 1 + 12)
  ## the redundant Sigma_g[g, 2, 1] rows are dropped (symmetric)
  out <- out[!(base == "Sigma_g" & i2 == 2L & i3 == 1L), ]
  rownames(out) <- NULL
  out
}

## ---- typology per group: the existing classifier on mu_state + delta_g -------------------------------------------------------------------
## mu_state (generated quantity) is the BILATERAL ideal point at x_ref; the leader / joiner ideal point of every state is mu_state + delta_g,
## added draw by draw. n_events_by_group: an S x 3 matrix of the state's usable events per group (columns bilateral, leader, joiner), so each
## table flags the states that never acted in that group (their classification in the group rests on a_i from the other groups and on delta_g).
## Returns a named list (bilateral, leader, joiner) of classify_types() tables with the extra columns group, n_events_in_group, n_events_as_initiator.
context_group_typology <- function(fit, state_labels, r0, n_events_by_group) {
  stopifnot(inherits(fit, "stanfit"), is.matrix(n_events_by_group), ncol(n_events_by_group) == 3L, nrow(n_events_by_group) == length(state_labels))
  mu <- mu_state_draws(fit)                                       # [draws, S, 2]
  stopifnot(dim(mu)[3] == 2L, dim(mu)[2] == length(state_labels))
  dl <- rstan::extract(fit, pars = "delta_lead")$delta_lead; dj <- rstan::extract(fit, pars = "delta_join")$delta_join   # [draws, 2]
  stopifnot(nrow(dl) == dim(mu)[1], ncol(dl) == 2L, nrow(dj) == dim(mu)[1], ncol(dj) == 2L)
  shifts <- list(bilateral = matrix(0, dim(mu)[1], 2), leader = dl, joiner = dj)
  out <- vector("list", 3); names(out) <- STAGE2_GROUPS
  for (g in seq_along(STAGE2_GROUPS)) {
    mu_g <- mu
    for (d in 1:2) mu_g[, , d] <- mu[, , d] + shifts[[g]][, d]      # the same draw's shift added to every state
    cls <- classify_types(mu_g, r0 = r0, state_labels = state_labels)
    cls$group <- STAGE2_GROUPS[g]
    cls$n_events_in_group <- as.integer(n_events_by_group[, g])
    cls$n_events_as_initiator <- as.integer(rowSums(n_events_by_group))
    cls$r0 <- r0
    out[[g]] <- cls
  }
  ## the bilateral table equals the plain classification of mu_state (shift 0), asserted
  stopifnot(identical(out$bilateral$map_type, classify_types(mu, r0 = r0, state_labels = state_labels)$map_type))
  out
}
