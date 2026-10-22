## ============================================================================
## typology.R -- classify initiators from posterior mu_i draws into the typology
##
## Regions in the (dM, dR) plane of the initiator's preferred target position mu_i:
##   peer        : r = ||mu_i|| < r0 (within radius r0 of the origin)
##   bully       : r >= r0, mu_1 > 0, mu_2 > 0      (+,+)
##   challenger  : r >= r0, mu_1 < 0, mu_2 < 0      (-,-)
##   hybrid_pm   : r >= r0, mu_1 > 0, mu_2 < 0      (+,-)  materially stronger, relationally weaker
##   hybrid_mp   : r >= r0, mu_1 < 0, mu_2 > 0      (-,+)
## r0 is an argument with no hidden default inside the classifier; callers must pass it
## (the sim and real-data scripts document their choice). Both dM and dR are differences
## of within-year standardised scores, so r0 is in SD units.
##
## Posterior type probabilities are the share of draws falling in each region.
## theta = atan2(mu_2, mu_1) is summarised by its circular mean; r by mean and quantiles.
##
## Every function checks array shapes (audit F1: the old polar_coords.r silently
## collapsed a 1 x N x K slice to a scalar).
## ============================================================================
TYPE_LEVELS <- c("bully", "peer", "challenger", "hybrid_pm", "hybrid_mp")
## One-dimensional models (D_mode 1 or 2): mu_i is a scalar on the single axis.
##   weaker   : mu >  r0   (prefers targets weaker than itself on that axis, x = own - target > 0)
##   peer     : |mu| <= r0
##   stronger : mu < -r0
TYPE_LEVELS_1D <- c("weaker", "peer", "stronger")

## mu_draws: numeric array [n_draws, S, 2]
## allow_identical_entities: per-event draws legitimately repeat when several events share one
## state-year (mu_event), so 03_extract.R passes TRUE there; for mu_state the default FALSE catches
## the F1 failure mode (one vector recycled to every entity).
classify_types <- function(mu_draws, r0, state_labels = NULL, allow_identical_entities = FALSE) {
  stopifnot(is.numeric(r0), length(r0) == 1, r0 > 0)
  stopifnot(is.array(mu_draws), is.numeric(mu_draws), length(dim(mu_draws)) == 3, dim(mu_draws)[3] %in% 1:2)
  if (dim(mu_draws)[3] == 1) return(classify_types_1d(mu_draws, r0, state_labels, allow_identical_entities))
  n_draws <- dim(mu_draws)[1]; S <- dim(mu_draws)[2]
  stopifnot(n_draws > 1, S >= 1, !anyNA(mu_draws))
  if (is.null(state_labels)) state_labels <- as.character(seq_len(S))
  stopifnot(length(state_labels) == S, !anyDuplicated(state_labels))    # unique entity ids
  ## F1 guard (shape / identity, not a statistic): the draws must not be one recycled vector.
  ## A recycled scalar or vector makes every entity's draws identical to entity 1's.
  if (S > 1 && !allow_identical_entities) {
    all_same <- all(vapply(2:S, function(s) isTRUE(all.equal(mu_draws[, s, ], mu_draws[, 1, ])), logical(1)))
    if (all_same) stop("classify_types: every entity has identical draws; the input looks recycled (audit F1)")
  }
  m1 <- mu_draws[, , 1, drop = TRUE]; m2 <- mu_draws[, , 2, drop = TRUE]
  if (S == 1) { m1 <- matrix(m1, ncol = 1); m2 <- matrix(m2, ncol = 1) }
  stopifnot(identical(dim(m1), c(n_draws, S)), identical(dim(m2), c(n_draws, S)))
  r <- sqrt(m1^2 + m2^2)
  theta <- atan2(m2, m1)
  region <- matrix(NA_character_, n_draws, S)
  region[r < r0] <- "peer"
  region[r >= r0 & m1 > 0 & m2 > 0] <- "bully"
  region[r >= r0 & m1 < 0 & m2 < 0] <- "challenger"
  region[r >= r0 & m1 > 0 & m2 <= 0] <- "hybrid_pm"
  region[r >= r0 & m1 <= 0 & m2 > 0] <- "hybrid_mp"
  region[r >= r0 & m1 <= 0 & m2 <= 0 & is.na(region)] <- "challenger"   # axis cases
  stopifnot(!anyNA(region))
  probs <- t(sapply(seq_len(S), function(s) {
    tab <- table(factor(region[, s], levels = TYPE_LEVELS)) / n_draws
    as.numeric(tab)
  }))
  if (S == 1) probs <- matrix(probs, nrow = 1)
  colnames(probs) <- paste0("p_", TYPE_LEVELS)
  stopifnot(identical(dim(probs), c(S, 5L)), all(abs(rowSums(probs) - 1) < 1e-8))
  circ_mean <- function(th) atan2(mean(sin(th)), mean(cos(th)))
  out <- data.frame(state = state_labels,
                    mu1_mean = colMeans(m1), mu2_mean = colMeans(m2),
                    mu1_lo = apply(m1, 2, quantile, 0.05), mu1_hi = apply(m1, 2, quantile, 0.95),
                    mu2_lo = apply(m2, 2, quantile, 0.05), mu2_hi = apply(m2, 2, quantile, 0.95),
                    r_mean = colMeans(r), r_lo = apply(r, 2, quantile, 0.05), r_hi = apply(r, 2, quantile, 0.95),
                    theta_circmean = apply(theta, 2, circ_mean),
                    theta_lo = apply(theta, 2, quantile, 0.05), theta_hi = apply(theta, 2, quantile, 0.95),
                    probs, stringsAsFactors = FALSE)
  out$map_type <- TYPE_LEVELS[apply(probs, 1, which.max)]
  out$map_prob <- apply(probs, 1, max)
  rownames(out) <- NULL
  stopifnot(nrow(out) == S, identical(dim(r), c(n_draws, S)))            # output rows = entities, r per entity x draw
  out
}

## Extract mu_state draws from an rstan fit as [n_draws, S, 2]
mu_state_draws <- function(fit) {
  arr <- rstan::extract(fit, pars = "mu_state")$mu_state
  if (length(dim(arr)) == 2) arr <- array(arr, c(dim(arr), 1L))     # D = 1: rstan drops the unit dimension
  stopifnot(length(dim(arr)) == 3, dim(arr)[3] %in% 1:2)
  arr
}

## Per-event classification (mu at the event's own X): same shape contract.
mu_event_draws <- function(fit) {
  arr <- rstan::extract(fit, pars = "mu_event")$mu_event
  if (length(dim(arr)) == 2) arr <- array(arr, c(dim(arr), 1L))
  stopifnot(length(dim(arr)) == 3, dim(arr)[3] %in% 1:2)
  arr
}

## Summarise the typology table: counts and mean posterior probability per MAP type
typology_summary <- function(cls) {
  stopifnot(all(c("map_type", "map_prob") %in% names(cls)))
  cls |> dplyr::count(map_type, name = "n_states") |>
    dplyr::left_join(cls |> dplyr::group_by(map_type) |> dplyr::summarise(mean_map_prob = mean(map_prob), .groups = "drop"),
                     by = "map_type")
}

## ---- one-dimensional classification: mu_draws is [n_draws, S, 1] ---------------------------
classify_types_1d <- function(mu_draws, r0, state_labels = NULL, allow_identical_entities = FALSE) {
  stopifnot(is.array(mu_draws), length(dim(mu_draws)) == 3, dim(mu_draws)[3] == 1, !anyNA(mu_draws))
  n_draws <- dim(mu_draws)[1]; S <- dim(mu_draws)[2]
  stopifnot(n_draws > 1, S >= 1)
  if (is.null(state_labels)) state_labels <- as.character(seq_len(S))
  stopifnot(length(state_labels) == S, !anyDuplicated(state_labels))
  m <- matrix(mu_draws[, , 1], n_draws, S)
  if (S > 1 && !allow_identical_entities) {
    all_same <- all(vapply(2:S, function(s) isTRUE(all.equal(m[, s], m[, 1])), logical(1)))
    if (all_same) stop("classify_types_1d: every entity has identical draws; the input looks recycled (audit F1)")
  }
  region <- matrix("peer", n_draws, S)
  region[m > r0] <- "weaker"; region[m < -r0] <- "stronger"
  probs <- t(sapply(seq_len(S), function(s) as.numeric(table(factor(region[, s], levels = TYPE_LEVELS_1D)) / n_draws)))
  if (S == 1) probs <- matrix(probs, nrow = 1)
  colnames(probs) <- paste0("p_", TYPE_LEVELS_1D)
  stopifnot(identical(dim(probs), c(S, 3L)), all(abs(rowSums(probs) - 1) < 1e-8))
  out <- data.frame(state = state_labels, mu_mean = colMeans(m),
                    mu_lo = apply(m, 2, quantile, 0.05), mu_hi = apply(m, 2, quantile, 0.95),
                    r_mean = colMeans(abs(m)), probs, stringsAsFactors = FALSE)
  out$map_type <- TYPE_LEVELS_1D[apply(probs, 1, which.max)]
  out$map_prob <- apply(probs, 1, max)
  rownames(out) <- NULL
  stopifnot(nrow(out) == S)
  out
}
