## ============================================================================
## context_compare.R -- functions of the context split comparison (dual typology stage a, PI-approved spec 2026-10-25)
## Sourced by 06_context_compare.R (and the tests). Inputs are saved fit objects of 02_fit_model.R for the SAME sample, lag, model and
## other switches, differing in the `context` switch: bilateral (_bi), coalition (_co) and, optionally, coalition without the top-3
## coalition disputes (_co_drop3). Nothing here fits anything or adds a parameter: every statistic is read off the saved draws and data.
##
## context_compare(fits, joiner_lookup, r0, min_decisions = 3) returns
##   $stats      long table (block, comparison, axis, statistic, estimate, lo, hi, n_states, note):
##               block "context": coalition vs bilateral, per axis (material, relational), over the initiating states with at least
##                 min_decisions WEIGHTED stage-2 decisions (sum of ev_w) in BOTH fits:
##                 cor_pmean        correlation of the posterior-mean mu_state across the two contexts (point estimate, no interval)
##                 cor_draw         the same correlation computed per posterior draw, draws paired by index across the two fits
##                                  (estimate = posterior mean, lo / hi = 5 / 95% quantiles)
##                 mean_shift       mean over those states of mu_state(coalition) - mu_state(bilateral), per draw (mean, 5 / 95%)
##                 type_agreement   share of those states with the same MAP type (classify_types(), r0) in both fits
##               block "leverage": coalition_drop3 vs coalition, the same four statistics (NULL when the drop3 fit is absent)
##               block "role": within each context, raw WEIGHTED mean chosen-target gap (the scaled x of the chosen alternative) for
##                 originator-pair and joiner-pair events, per axis; and, within the coalition context, across the states with at least
##                 min_decisions weighted decisions in BOTH roles, the correlation of the per-state weighted mean chosen gaps and the mean
##                 difference joiner - originator. mu_state is one a_i per state and is NOT split by role, so these are data descriptives,
##                 not posterior quantities (no interval; the note says so).
##   $crosstab   MAP-type cross-tabulations (block, type_a, type_b, n) for the context and leverage blocks
##   $states     the per-state table behind the context block (state, weighted decisions and MAP type per context, mu means)
## ============================================================================
suppressPackageStartupMessages({ library(dplyr) })

## weighted stage-2 decisions per state (sum of the event weights of the state's events)
weighted_decisions <- function(dat) {
  stopifnot(all(c("S", "ev_sy", "sy_state", "ev_w") %in% names(dat)), length(dat$ev_w) == length(dat$ev_sy))
  w <- rep(0, dat$S); idx <- dat$sy_state[dat$ev_sy]
  for (k in seq_along(idx)) w[idx[k]] <- w[idx[k]] + dat$ev_w[k]
  stopifnot(length(w) == dat$S, all(is.finite(w)), abs(sum(w) - sum(dat$ev_w)) < 1e-8)
  w
}
## per-event table of a fit: initiator / target ccodes, year, weight, the chosen target's scaled gap on both axes
fit_events <- function(obj) {
  d <- obj$dat
  stopifnot(d$D_mode == 3L, !is.null(obj$state_codes), length(obj$state_codes) == d$S, !is.null(obj$sy_keys), nrow(obj$sy_keys) == d$N1,
            ncol(d$x) == 2L, nrow(d$x) == d$N2)
  data.frame(initiator = obj$state_codes[d$sy_state[d$ev_sy]], target = obj$state_codes[d$alt_target[d$ev_row]], year = obj$sy_keys$year[d$ev_sy],
             w = d$ev_w, x_material = d$x[d$ev_row, 1], x_relational = d$x[d$ev_row, 2])
}
## the four statistics of fit b against fit a over the states `R` (indices), draws paired by index
pair_block <- function(mu_a, mu_b, R, r0, labels, block, comparison) {
  stopifnot(identical(dim(mu_a), dim(mu_b)), length(dim(mu_a)) == 3, dim(mu_a)[3] == 2L, length(R) >= 3)
  n_draws <- dim(mu_a)[1]
  ca <- classify_types(mu_a, r0 = r0, state_labels = labels); cb <- classify_types(mu_b, r0 = r0, state_labels = labels)
  rows <- list(); ct <- list()
  for (d in 1:2) {
    ax <- c("material", "relational")[d]
    A <- mu_a[, R, d, drop = TRUE]; B <- mu_b[, R, d, drop = TRUE]        # draws x states
    stopifnot(identical(dim(A), c(n_draws, length(R))))
    cor_pm <- cor(colMeans(A), colMeans(B))
    cor_dr <- vapply(seq_len(n_draws), function(k) suppressWarnings(cor(A[k, ], B[k, ])), numeric(1))
    shift <- rowMeans(B - A)
    stopifnot(is.finite(cor_pm), all(is.finite(shift)))
    rows[[length(rows) + 1]] <- data.frame(block = block, comparison = comparison, axis = ax,
      statistic = c("cor_pmean", "cor_draw", "mean_shift"),
      estimate = c(cor_pm, mean(cor_dr, na.rm = TRUE), mean(shift)),
      lo = c(NA, quantile(cor_dr, 0.05, na.rm = TRUE), quantile(shift, 0.05)), hi = c(NA, quantile(cor_dr, 0.95, na.rm = TRUE), quantile(shift, 0.95)),
      n_states = length(R), note = c("correlation of posterior means (point estimate)", sprintf("per-draw correlation, %d draws paired by index; mean and 5 / 95%% quantiles", n_draws),
                                     "mean over the states of mu_b - mu_a per draw; mean and 5 / 95% quantiles"), stringsAsFactors = FALSE)
  }
  agree <- ca$map_type[R] == cb$map_type[R]
  rows[[length(rows) + 1]] <- data.frame(block = block, comparison = comparison, axis = "both", statistic = "type_agreement", estimate = mean(agree), lo = NA, hi = NA,
                                         n_states = length(R), note = sprintf("share of states with the same MAP type at r0 = %.2f (%d of %d)", r0, sum(agree), length(R)), stringsAsFactors = FALSE)
  tab <- table(factor(ca$map_type[R], levels = TYPE_LEVELS), factor(cb$map_type[R], levels = TYPE_LEVELS))
  ct <- as.data.frame(tab, stringsAsFactors = FALSE); names(ct) <- c("type_a", "type_b", "n"); ct$block <- block; ct$comparison <- comparison
  stopifnot(sum(ct$n) == length(R))
  list(stats = bind_rows(rows), crosstab = ct[, c("block", "comparison", "type_a", "type_b", "n")], cls_a = ca, cls_b = cb)
}
context_compare <- function(fits, joiner_lookup, r0, min_decisions = 3) {
  stopifnot(is.list(fits), all(c("bilateral", "coalition") %in% names(fits)), is.numeric(r0), r0 > 0, min_decisions > 0,
            all(c("initiator", "target", "year", "is_joiner_pair") %in% names(joiner_lookup)))
  bi <- fits$bilateral; co <- fits$coalition; d3 <- fits$coalition_drop3
  for (o in c(list(bi, co), if (!is.null(d3)) list(d3))) stopifnot(!is.null(o$fit), !is.null(o$dat), o$dat$D_mode == 3L, identical(o$state_codes, bi$state_codes))
  ## stage 1 is unchanged across contexts: the same state-years and the same y1
  stopifnot(bi$dat$N1 == co$dat$N1, identical(bi$dat$y1, co$dat$y1), identical(bi$dat$sy_state, co$dat$sy_state))
  if (!is.null(d3)) stopifnot(d3$dat$N1 == co$dat$N1, identical(d3$dat$y1, co$dat$y1))
  labels <- as.character(bi$state_codes)
  mu_bi <- mu_state_draws(bi$fit); mu_co <- mu_state_draws(co$fit)
  if (!identical(dim(mu_bi), dim(mu_co))) stop(sprintf("the bilateral and coalition fits have different mu_state draw dimensions (%s vs %s); draws are paired by index, so the fits need the same chains x iterations", paste(dim(mu_bi), collapse = "x"), paste(dim(mu_co), collapse = "x")))
  w_bi <- weighted_decisions(bi$dat); w_co <- weighted_decisions(co$dat)
  R <- which(w_bi >= min_decisions & w_co >= min_decisions)
  if (length(R) < 3) stop(sprintf("only %d states have at least %g weighted decisions in both contexts", length(R), min_decisions))
  pb <- pair_block(mu_bi, mu_co, R, r0, labels, "context", "coalition_vs_bilateral")
  stats <- pb$stats; crosstab <- pb$crosstab
  states <- data.frame(state = labels, w_bilateral = w_bi, w_coalition = w_co, in_comparison = seq_along(labels) %in% R,
                       type_bilateral = pb$cls_a$map_type, type_coalition = pb$cls_b$map_type,
                       mu1_bilateral = pb$cls_a$mu1_mean, mu2_bilateral = pb$cls_a$mu2_mean, mu1_coalition = pb$cls_b$mu1_mean, mu2_coalition = pb$cls_b$mu2_mean, stringsAsFactors = FALSE)
  ## leverage: coalition without the top-3 coalition disputes against the full coalition fit
  if (!is.null(d3)) {
    mu_d3 <- mu_state_draws(d3$fit)
    if (!identical(dim(mu_d3), dim(mu_co))) stop("the coalition and coalition_drop3 fits have different mu_state draw dimensions")
    w_d3 <- weighted_decisions(d3$dat)
    R3 <- which(w_co >= min_decisions & w_d3 >= min_decisions)
    if (length(R3) < 3) stop("fewer than 3 states have enough weighted decisions in both the coalition and the coalition_drop3 fits")
    pl <- pair_block(mu_co, mu_d3, R3, r0, labels, "leverage", "coalition_drop3_vs_coalition")
    stats <- bind_rows(stats, pl$stats); crosstab <- bind_rows(crosstab, pl$crosstab)
    states$w_coalition_drop3 <- w_d3; states$type_coalition_drop3 <- pl$cls_b$map_type
  }
  ## role split: raw weighted mean chosen-target gaps (mu_state is not role-specific; said in the note)
  role_note <- "DATA DESCRIPTIVE, not a posterior quantity: mu_state is one a_i per state and the fitted model does not split it by role; weighted mean of the chosen target's scaled gap"
  role_rows <- list()
  per_state_role <- NULL
  for (ctx in c("bilateral", "coalition")) {
    fe <- fit_events(fits[[ctx]]) |> left_join(joiner_lookup |> distinct(initiator, target, year, is_joiner_pair), by = c("initiator", "target", "year"))
    if (anyNA(fe$is_joiner_pair)) stop(sprintf("%d events of the %s fit have no is_joiner_pair in the events table", sum(is.na(fe$is_joiner_pair)), ctx))
    for (role in c("originator", "joiner")) {
      sub <- fe[fe$is_joiner_pair == (role == "joiner"), ]
      for (d in 1:2) {
        ax <- c("material", "relational")[d]; xv <- sub[[paste0("x_", ax)]]
        role_rows[[length(role_rows) + 1]] <- data.frame(block = "role", comparison = paste(ctx, role, sep = "_"), axis = ax, statistic = "weighted_mean_chosen_gap",
          estimate = if (nrow(sub) > 0) sum(sub$w * xv) / sum(sub$w) else NA_real_, lo = NA, hi = NA, n_states = n_distinct(sub$initiator),
          note = sprintf("%s; %d events, weights sum %.1f", role_note, nrow(sub), sum(sub$w)), stringsAsFactors = FALSE)
      }
    }
    if (ctx == "coalition") {
      ps <- fe |> group_by(initiator, is_joiner_pair) |> summarise(w = sum(w), m1 = sum(w * x_material) / sum(w), m2 = sum(w * x_relational) / sum(w), .groups = "drop")
      wide <- ps |> filter(!is_joiner_pair) |> select(initiator, w_orig = w, m1_orig = m1, m2_orig = m2) |>
        inner_join(ps |> filter(is_joiner_pair) |> select(initiator, w_join = w, m1_join = m1, m2_join = m2), by = "initiator") |>
        filter(w_orig >= min_decisions, w_join >= min_decisions)
      per_state_role <- wide
      for (d in 1:2) {
        ax <- c("material", "relational")[d]; a <- wide[[paste0("m", d, "_orig")]]; b <- wide[[paste0("m", d, "_join")]]
        role_rows[[length(role_rows) + 1]] <- data.frame(block = "role", comparison = "coalition_joiner_vs_originator_by_state", axis = ax,
          statistic = c("cor_state_means", "mean_shift"), estimate = c(if (nrow(wide) >= 3) cor(a, b) else NA_real_, if (nrow(wide) >= 1) mean(b - a) else NA_real_), lo = NA, hi = NA,
          n_states = nrow(wide), note = paste0(role_note, sprintf("; states with at least %g weighted decisions in BOTH roles within the coalition context", min_decisions)), stringsAsFactors = FALSE)
      }
    }
  }
  stats <- bind_rows(stats, bind_rows(role_rows))
  rownames(stats) <- NULL
  stopifnot(!anyNA(stats$statistic), all(stats$block %in% c("context", "leverage", "role")))
  list(stats = stats, crosstab = crosstab, states = states, per_state_role = per_state_role, min_decisions = min_decisions, r0 = r0)
}
## a short markdown summary of context_compare() output
context_compare_md <- function(cc, title) {
  st <- cc$stats
  f <- function(block, comparison, ax, stat) { r <- st[st$block == block & st$comparison == comparison & st$axis == ax & st$statistic == stat, ]; stopifnot(nrow(r) == 1); r }
  line <- function(block, comparison, ax) {
    a <- f(block, comparison, ax, "cor_pmean"); b <- f(block, comparison, ax, "cor_draw"); m <- f(block, comparison, ax, "mean_shift")
    sprintf("| %s | %.3f | %.3f [%.3f, %.3f] | %.3f [%.3f, %.3f] | %d |", ax, a$estimate, b$estimate, b$lo, b$hi, m$estimate, m$lo, m$hi, a$n_states)
  }
  out <- c(sprintf("# %s", title), "",
           sprintf("States with at least %g weighted stage-2 decisions in both fits of a comparison; r0 = %.2f; draws paired by index. Stage 1 is identical across contexts (asserted); the fitted model has one mu_state per state, so the role split reports raw weighted chosen-target gaps, not posterior quantities.", cc$min_decisions, cc$r0), "",
           "## Coalition vs bilateral (mu_state, coalition minus bilateral)", "",
           "| axis | cor of posterior means | per-draw cor, mean [5%, 95%] | mean shift [5%, 95%] | states |", "|---|---|---|---|---|",
           line("context", "coalition_vs_bilateral", "material"), line("context", "coalition_vs_bilateral", "relational"), "",
           sprintf("MAP type agreement: %.1f%% (%s)", 100 * f("context", "coalition_vs_bilateral", "both", "type_agreement")$estimate, f("context", "coalition_vs_bilateral", "both", "type_agreement")$note), "",
           "Cross-tab (rows bilateral, columns coalition):", "", crosstab_md(cc$crosstab[cc$crosstab$block == "context", ]), "")
  if (any(st$block == "leverage")) {
    out <- c(out, "## Leverage: coalition without the top-3 coalition disputes vs coalition (drop3 minus coalition)", "",
             "| axis | cor of posterior means | per-draw cor, mean [5%, 95%] | mean shift [5%, 95%] | states |", "|---|---|---|---|---|",
             line("leverage", "coalition_drop3_vs_coalition", "material"), line("leverage", "coalition_drop3_vs_coalition", "relational"), "",
             sprintf("MAP type agreement: %.1f%% (%s)", 100 * f("leverage", "coalition_drop3_vs_coalition", "both", "type_agreement")$estimate, f("leverage", "coalition_drop3_vs_coalition", "both", "type_agreement")$note), "",
             crosstab_md(cc$crosstab[cc$crosstab$block == "leverage", ]), "")
  } else out <- c(out, "## Leverage", "", "No coalition_drop3 fit was given.", "")
  rr <- st[st$block == "role" & st$statistic == "weighted_mean_chosen_gap", ]
  out <- c(out, "## Originators vs joiners (raw weighted mean chosen-target gap, scaled units; data descriptive)", "",
           "| context | role | material | relational | events (weights) |", "|---|---|---|---|---|")
  for (cmp in unique(rr$comparison)) {
    a <- rr[rr$comparison == cmp & rr$axis == "material", ]; b <- rr[rr$comparison == cmp & rr$axis == "relational", ]
    out <- c(out, sprintf("| %s | %s | %.3f | %.3f | %s |", sub("_.*$", "", cmp), sub("^.*_", "", cmp), a$estimate, b$estimate, sub("^.*; ", "", a$note)))
  }
  rs <- st[st$block == "role" & st$comparison == "coalition_joiner_vs_originator_by_state", ]
  out <- c(out, "", sprintf("Within the coalition context, per state (joiner minus originator; %d states with enough decisions in both roles): material cor %.3f, shift %.3f; relational cor %.3f, shift %.3f.",
                            rs$n_states[1], rs$estimate[rs$axis == "material" & rs$statistic == "cor_state_means"], rs$estimate[rs$axis == "material" & rs$statistic == "mean_shift"],
                            rs$estimate[rs$axis == "relational" & rs$statistic == "cor_state_means"], rs$estimate[rs$axis == "relational" & rs$statistic == "mean_shift"]))
  out
}
crosstab_md <- function(ct) {
  lv <- TYPE_LEVELS
  m <- matrix(0L, length(lv), length(lv), dimnames = list(lv, lv))
  for (k in seq_len(nrow(ct))) m[ct$type_a[k], ct$type_b[k]] <- ct$n[k]
  c(paste0("| | ", paste(lv, collapse = " | "), " |"), paste0("|---|", paste(rep("---", length(lv)), collapse = "|"), "|"),
    vapply(lv, function(a) paste0("| ", a, " | ", paste(m[a, ], collapse = " | "), " |"), character(1)))
}
