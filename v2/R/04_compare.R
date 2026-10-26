## ============================================================================
## 04_compare.R -- cross-model comparison for one (sample, lag) cell group: LOO, joint vs independent, M1 vs M2
##
## Usage: Rscript v2/R/04_compare.R --sample=mid --lag=5 [switches of spec.R] [--r0=0.5] [--allow_stale=FALSE]
## Loads the eight fits fit_<sample>_lag<h>_all_pr_{M1m,M1r,M1a,M1b,M2m,M2r,M2a,M2b}[_K<K>]<switch suffix>.rds (the switches name
## the variant group, e.g. --ww2=drop compares the _noww2 fits among themselves)
## (missing or UNREADABLE fits are skipped and listed; pilot fits are refused; STALE fits only with --allow_stale) and writes:
##   compare_loo_<group>.csv        per model: LOO ELPD of (i) stage 2 (target choice), (ii) stage 1 CONDITIONAL (the Bernoulli
##                                  term alone, the state-year's events staying in) and (iii) the JOINT state-year LOO
##                                  (log_lik_stage1[n] + the weighted sum of the state-year's log_lik_stage2; review 2026-10-25,
##                                  B1), with SEs and Pareto-k > 0.7 counts
##   compare_contrasts_<group>.csv  planned contrasts with elpd_diff and se_diff from loo_compare for all three LOOs:
##                                  within each family: r vs m, a vs m, a vs r, b vs a; across families: M2 vs M1 for m, r, a, b
##   compare_misclass_<group>.csv   a vs b MAP-type disagreement (initiating states), within M1 and within M2
##   compare_m1m2_<group>.csv       M1x vs M2x (x = m, r, a, b; PI question "do the effects survive the M2 latent factors?"):
##                                  posterior means and 95% intervals of tau, mu0, beta, gamma, Sigma[1,2] and corr(a_1, a_2)
##                                  side by side, sign agreement, whether each interval excludes 0, plus (per x) the correlation of
##                                  the mu_state posterior means between M1x and M2x and the type-agreement rate at r0
## All models score the SAME events and choice sets (same assemble_cell() output), asserted below, so LOO is like with like.
## ============================================================================
suppressPackageStartupMessages({ library(rstan); library(loo); library(dplyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
R_DIR <- dirname(this_file)
source(file.path(R_DIR, "spec.R")); source(file.path(R_DIR, "helpers_data.R")); source(file.path(R_DIR, "helpers_model.R"))
source(file.path(R_DIR, "typology.R")); source(file.path(R_DIR, "measurement_tests.R"))
Sys.setenv(V2_ROOT = normalizePath(file.path(R_DIR, "..")))
args <- commandArgs(trailingOnly = TRUE)
SAMPLE <- get_arg(args, "sample", "mid"); LAG <- as.integer(get_arg(args, "lag", 5))
stopifnot(SAMPLE %in% ALL_SAMPLES, LAG %in% LAGS)
spec <- parse_spec(args); R0 <- get_arg(args, "r0", R0); ALLOW_STALE <- get_arg(args, "allow_stale", FALSE)
refuse_archived_options(args, "04_compare.R")            # the same archived flags 02_fit_model.R refuses (audit fix 14)
group <- compare_group_name(SAMPLE, LAG, spec)
log_open(paste0("compare_", group))
logf("group %s; spec: %s", group, spec_string(spec))
out_dir <- file.path(v2_root(), "output")
manifest_path <- file.path(spec_data_dir(spec, v2_root()), "manifest.csv")
fits <- list(); loos <- list(); rows <- list()
for (m in MODELS) {
  ## the M1 fits carry no K (cell_spec(): m1_spec(spec) for the M1 family), so a non-default --ame_K compares the K=3 M2 fits with the
  ## one set of M1 fits (review 2026-10-26, P2)
  f <- file.path(out_dir, paste0("fit_", cell_name(SAMPLE, LAG, m, cell_spec(m, spec)), ".rds"))
  if (!file.exists(f) && ALLOW_STALE && file.exists(sub("\\.rds$", "_STALE.rds", f))) f <- sub("\\.rds$", "_STALE.rds", f)
  if (!file.exists(f)) { logf("missing fit: %s", basename(f)); next }
  obj <- read_fit_or_null(f)                       # C2: an unreadable fit is logged and treated as missing
  if (is.null(obj)) next
  assert_production_fit(obj, allow_stale = ALLOW_STALE)        # C1 / C2
  assert_fit_current(obj, manifest_path); assert_fit_stan_current(obj)   # current data AND current Stan program (audit fix 5)
  fits[[m]] <- obj
  l2 <- loo_from_fit(obj$fit, "log_lik_stage2"); l1 <- loo_from_fit(obj$fit, "log_lik_stage1"); lj <- loo_joint_sy(obj$fit, obj$dat)
  loos[[m]] <- list(l2 = l2, l1 = l1, lj = lj)
  rows[[m]] <- data.frame(model = m, D_mode = obj$dat$D_mode, joint = obj$dat$joint, E = obj$dat$E, N1 = obj$dat$N1, sum_event_weights = sum(obj$dat$ev_w),
                          elpd_stage2 = l2$estimates["elpd_loo", 1], se_stage2 = l2$estimates["elpd_loo", 2], bad_k_stage2 = sum(l2$diagnostics$pareto_k > 0.7),
                          elpd_stage1_conditional = l1$estimates["elpd_loo", 1], se_stage1_conditional = l1$estimates["elpd_loo", 2], bad_k_stage1_conditional = sum(l1$diagnostics$pareto_k > 0.7),
                          elpd_joint_stateyear = lj$estimates["elpd_loo", 1], se_joint_stateyear = lj$estimates["elpd_loo", 2], bad_k_joint_stateyear = sum(lj$diagnostics$pareto_k > 0.7))
  logf("%s: elpd2 = %.1f (%.1f), elpd1 conditional = %.1f (%.1f), joint state-year = %.1f (%.1f)", m, rows[[m]]$elpd_stage2, rows[[m]]$se_stage2,
       rows[[m]]$elpd_stage1_conditional, rows[[m]]$se_stage1_conditional, rows[[m]]$elpd_joint_stateyear, rows[[m]]$se_joint_stateyear)
}
stopifnot(length(fits) > 0)
## every loaded fit must score the same events (like with like), with the same weights
Es <- vapply(fits, function(o) o$dat$E, numeric(1)); stopifnot(length(unique(Es)) == 1)
evs <- lapply(fits, function(o) paste(o$dat$ev_sy, o$dat$alt_target[o$dat$ev_row])); stopifnot(all(vapply(evs, identical, logical(1), evs[[1]])))
ws <- lapply(fits, function(o) round(o$dat$ev_w, 10)); stopifnot(all(vapply(ws, identical, logical(1), ws[[1]])))
tab <- bind_rows(rows) |> mutate(delta_elpd_stage2_vs_best = elpd_stage2 - max(elpd_stage2),
                                 delta_elpd_stage1_conditional_vs_best = elpd_stage1_conditional - max(elpd_stage1_conditional),
                                 delta_elpd_joint_stateyear_vs_best = elpd_joint_stateyear - max(elpd_joint_stateyear))
write.csv(tab, file.path(out_dir, paste0("compare_loo_", group, ".csv")), row.names = FALSE)
## planned contrasts
contrast <- function(a, b, label) {
  if (!all(c(a, b) %in% names(loos))) return(NULL)
  c2 <- loo::loo_compare(loos[[a]]$l2, loos[[b]]$l2); c1 <- loo::loo_compare(loos[[a]]$l1, loos[[b]]$l1); cj <- loo::loo_compare(loos[[a]]$lj, loos[[b]]$lj)
  data.frame(contrast = label, model_a = a, model_b = b,
             stage2_elpd_diff_b_minus_a = loos[[b]]$l2$estimates["elpd_loo", 1] - loos[[a]]$l2$estimates["elpd_loo", 1], stage2_se_diff = c2[2, "se_diff"],
             stage1_conditional_elpd_diff_b_minus_a = loos[[b]]$l1$estimates["elpd_loo", 1] - loos[[a]]$l1$estimates["elpd_loo", 1], stage1_conditional_se_diff = c1[2, "se_diff"],
             joint_stateyear_elpd_diff_b_minus_a = loos[[b]]$lj$estimates["elpd_loo", 1] - loos[[a]]$lj$estimates["elpd_loo", 1], joint_stateyear_se_diff = cj[2, "se_diff"])
}
cons <- list()
for (fam in c("M1", "M2")) {
  cons[[length(cons) + 1]] <- contrast(paste0(fam, "m"), paste0(fam, "r"), paste0(fam, ": relational vs material"))
  cons[[length(cons) + 1]] <- contrast(paste0(fam, "m"), paste0(fam, "a"), paste0(fam, ": independent (both) vs material"))
  cons[[length(cons) + 1]] <- contrast(paste0(fam, "r"), paste0(fam, "a"), paste0(fam, ": independent (both) vs relational"))
  cons[[length(cons) + 1]] <- contrast(paste0(fam, "a"), paste0(fam, "b"), paste0(fam, ": joint vs independent"))
}
for (sp in c("m", "r", "a", "b")) cons[[length(cons) + 1]] <- contrast(paste0("M1", sp), paste0("M2", sp), paste0(sp, ": M2 (AME) vs M1"))
con_tab <- bind_rows(cons); print(con_tab)
write.csv(con_tab, file.path(out_dir, paste0("compare_contrasts_", group, ".csv")), row.names = FALSE)
## misclassification a vs b within each family
mis <- list()
for (fam in c("M1", "M2")) {
  a <- paste0(fam, "a"); b <- paste0(fam, "b")
  if (!all(c(a, b) %in% names(fits))) next
  init_states <- which(tabulate(fits[[a]]$dat$sy_state[fits[[a]]$dat$ev_sy], nbins = fits[[a]]$dat$S) > 0)
  mc <- misclassification_ab(fits[[a]]$fit, fits[[b]]$fit, r0 = R0, only_states = init_states)
  mis[[fam]] <- data.frame(family = fam, r0 = R0, n_initiating_states = length(init_states), disagreement_rate = mc$rate, mean_abs_prob_shift = mc$mean_abs_prob_shift)
}
if (length(mis)) write.csv(bind_rows(mis), file.path(out_dir, paste0("compare_misclass_", group, ".csv")), row.names = FALSE)

## ---- M1x vs M2x: do the effects survive the M2 latent factors? (PI question, 2026-10-25) --------------------------------------------
## For each x the structural parameters of M1x and M2x side by side: posterior mean and 95% interval, sign agreement of the means,
## whether each interval excludes 0; then the correlation of the mu_state posterior means (per axis, over the initiating states) and the
## type-agreement rate at r0 (MAP types of the initiating states). The M2 program is untouched; this only reads the two fits.
m1m2 <- list(); m1m2_state <- list()
for (sp in c("m", "r", "a", "b")) {
  a <- paste0("M1", sp); b <- paste0("M2", sp)
  if (!all(c(a, b) %in% names(fits))) next
  oa <- fits[[a]]; ob <- fits[[b]]
  stopifnot(oa$dat$D_mode == ob$dat$D_mode, identical(oa$W_cols, ob$W_cols), identical(oa$X_cols, ob$X_cols), identical(oa$Z_cols, ob$Z_cols))
  D_mode <- oa$dat$D_mode; D <- if (D_mode == 3L) 2L else 1L
  pars <- c("tau", "mu0", "beta", "gamma", if (D == 2L) "Sigma", if (D == 2L) "R_Om")
  summ <- function(o) { s <- summary(o$fit, pars = pars, probs = c(0.025, 0.975))$summary
    data.frame(parameter = rownames(s), mean = s[, "mean"], lo = s[, "2.5%"], hi = s[, "97.5%"], stringsAsFactors = FALSE) }
  sa <- summ(oa); sb <- summ(ob)
  keep <- sa$parameter %in% c("tau") | grepl("^(mu0|beta|gamma)\\[", sa$parameter) | sa$parameter == "Sigma[1,2]" | sa$parameter == "R_Om[2,3]"
  sa <- sa[keep, ]; sb <- sb[match(sa$parameter, sb$parameter), ]
  stopifnot(identical(sa$parameter, sb$parameter), !anyNA(sb$mean))
  lab <- label_parameters(sa$parameter, oa$W_cols, oa$X_cols, oa$Z_cols, D_mode)
  lab$covariate[sa$parameter == "R_Om[2,3]"] <- "corr(a_1, a_2)"
  m1m2[[sp]] <- data.frame(spec = sp, model_m1 = a, model_m2 = b, lab,
                           m1_mean = sa$mean, m1_lo95 = sa$lo, m1_hi95 = sa$hi, m2_mean = sb$mean, m2_lo95 = sb$lo, m2_hi95 = sb$hi,
                           sign_agree = sign(sa$mean) == sign(sb$mean),
                           m1_excludes_zero = sa$lo > 0 | sa$hi < 0, m2_excludes_zero = sb$lo > 0 | sb$hi < 0,
                           note = ifelse(sa$parameter %in% c("Sigma[1,2]", "R_Om[2,3]") & sp == "a", "fixed at 0 in the a models", ""), stringsAsFactors = FALSE)
  ## mu_state posterior means and MAP types over the initiating states
  init_states <- which(tabulate(oa$dat$sy_state[oa$dat$ev_sy], nbins = oa$dat$S) > 0)
  mu_a <- mu_state_draws(oa$fit); mu_b <- mu_state_draws(ob$fit)
  stopifnot(dim(mu_a)[2] == dim(mu_b)[2], dim(mu_a)[3] == dim(mu_b)[3])
  ca <- classify_types(mu_a, r0 = R0); cb <- classify_types(mu_b, r0 = R0)
  cors <- vapply(seq_len(D), function(d) cor(apply(mu_a[, init_states, d, drop = FALSE], 2, mean), apply(mu_b[, init_states, d, drop = FALSE], 2, mean)), numeric(1))
  m1m2_state[[sp]] <- data.frame(spec = sp, model_m1 = a, model_m2 = b, r0 = R0, n_initiating_states = length(init_states),
                                 cor_mu_state_axis1 = cors[1], cor_mu_state_axis2 = if (D == 2L) cors[2] else NA_real_,
                                 axis1 = axis_labels(D_mode)[1], axis2 = if (D == 2L) axis_labels(D_mode)[2] else NA_character_,
                                 type_agreement_rate = mean(ca$map_type[init_states] == cb$map_type[init_states]), stringsAsFactors = FALSE)
  logf("M1%s vs M2%s: sign agreement on %d of %d parameters; cor(mu_state means) %s; type agreement %.3f", sp, sp, sum(m1m2[[sp]]$sign_agree), nrow(m1m2[[sp]]),
       paste(sprintf("%.3f", cors), collapse = " / "), m1m2_state[[sp]]$type_agreement_rate)
}
if (length(m1m2)) {
  write.csv(bind_rows(m1m2), file.path(out_dir, paste0("compare_m1m2_", group, ".csv")), row.names = FALSE)
  write.csv(bind_rows(m1m2_state), file.path(out_dir, paste0("compare_m1m2_states_", group, ".csv")), row.names = FALSE)
}
logf("04_compare done for %s", group)
