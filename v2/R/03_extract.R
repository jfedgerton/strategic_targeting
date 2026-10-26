## ============================================================================
## 03_extract.R -- post-process one fit: typology from posterior mu_i, AME alignment
##
## Usage: Rscript v2/R/03_extract.R --fit=v2/output/fit_<cell>.rds [--r0=0.5] [--band_M=r0] [--band_R=r0] [--allow_stale=FALSE]
## (r0 defaults to R0 in spec.R)
## Pilot fits (fit_*_pilot.rds, or a saved pilot configuration) are refused; a STALE fit (manifest changed while sampling) needs
## --allow_stale=TRUE (2026-10-25 C1 / C2). The manifest the fit is checked against is the one in the fit's own data directory
## (v2/data, or v2/data_sideA for --direction=sideA fits). The typology / mu_state CSVs of one-dimensional cells carry an `axis`
## column (material / relational; B2).
## Writes next to the fit:
##   typology_<cell>.csv   per-state posterior type probabilities, MAP type, r / theta summaries
##   mu_state_<cell>.csv   posterior mean and 90% interval of mu_i at x_ref
##   mu_event_<cell>.csv   per-event mu (initiator at its own X) with MAP type
##   ame_U_<cell>.csv / ame_V_<cell>.csv   (M2 only) Procrustes-aligned posterior means
##   ame_pairs_<cell>.csv  (M2 only) initiator x target u_i'v_j posterior mean, strongest 200 pairs
##   typology9_<cell>.csv  (appendix, D = 2) 9-cell classification with per-axis peer bands band_M / band_R
##   reach_<cell>.csv      (appendix) posterior mean / 90% interval of r_i = ||mu_i||, |mu_i1|, |mu_i2|
##   focus_sigma_<cell>.csv, focus_states_<cell>.csv  (appendix, descriptive) pooled Sigma summary and per-state
##                         empirical dispersion of chosen-target positions around mean mu_event (n < 5 flagged)
## r0 (peer radius, in scaled gap units: log ratio divided by the per-axis constant in change_scale.csv,
## so 0.5 is about a 5-fold CINC ratio on the material axis) is a documented argument; R0 = 0.5 in spec.R is the v2
## default and a PI decision to confirm.
## ============================================================================
suppressPackageStartupMessages({ library(rstan); library(dplyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
R_DIR <- dirname(this_file)
source(file.path(R_DIR, "spec.R")); source(file.path(R_DIR, "helpers_data.R")); source(file.path(R_DIR, "helpers_model.R"))
source(file.path(R_DIR, "typology.R"))
Sys.setenv(V2_ROOT = normalizePath(file.path(R_DIR, "..")))
args <- commandArgs(trailingOnly = TRUE)
fit_path <- sub("^--fit=", "", grep("^--fit=", args, value = TRUE))
r0 <- get_arg(args, "r0", R0)
allow_stale <- get_arg(args, "allow_stale", FALSE)
stopifnot(length(fit_path) == 1, file.exists(fit_path), is.finite(r0), r0 > 0, !is.na(allow_stale))
obj <- readRDS(fit_path); fit <- obj$fit; cell <- obj$cell
assert_production_fit(obj, allow_stale = allow_stale)
assert_fit_current(obj, file.path(v2_root(), if (is.null(obj$data_dir)) "data" else obj$data_dir, "manifest.csv"))
assert_fit_stan_current(obj)            # the fit's Stan source hash must be the compiled program's (stan/build/<family>.sha256; audit fix 5)
out_dir <- dirname(fit_path)
log_open(paste0("03_extract_", cell))
logf("fit %s; r0 = %.2f", cell, r0)

## ---- typology from mu_state (at x_ref) -----------------------------------------------------
mu <- mu_state_draws(fit)
labels <- if (!is.null(obj$state_codes)) as.character(obj$state_codes) else as.character(seq_len(dim(mu)[2]))
stopifnot(dim(mu)[2] == length(labels))
cls <- classify_types(mu, r0 = r0, state_labels = labels)
## states that never initiate in the sample have a_i identified only by the prior and Omega; flag them
n_init <- tabulate(obj$dat$sy_state[obj$dat$ev_sy], nbins = obj$dat$S)
cls$n_events_as_initiator <- n_init
cls$r0 <- r0
if (dim(mu)[3] == 1) cls$axis <- axis_labels(obj$dat$D_mode)     # B2: a one-dimensional cell names its axis
write.csv(cls, file.path(out_dir, paste0("typology_", cell, ".csv")), row.names = FALSE)
mu_cols <- intersect(c("mu1_mean", "mu1_lo", "mu1_hi", "mu2_mean", "mu2_lo", "mu2_hi", "mu_mean", "mu_lo", "mu_hi"), names(cls))
write.csv(cls |> select(state, all_of(mu_cols), any_of("axis"), n_events_as_initiator),
          file.path(out_dir, paste0("mu_state_", cell, ".csv")), row.names = FALSE)
logf("targeting dimensions D = %d (%s)", dim(mu)[3], if (dim(mu)[3] == 1) "one-axis typology: weaker / peer / stronger" else "bully / peer / challenger / hybrids")
print(typology_summary(cls |> filter(n_events_as_initiator > 0)))
logf("typology written; %d states with >= 1 initiation", sum(n_init > 0))

## ---- per-event mu (initiator at its own X) -------------------------------------------------------
mue <- mu_event_draws(fit)
ev_cls <- classify_types(mue, r0 = r0, state_labels = as.character(seq_len(dim(mue)[2])), allow_identical_entities = TRUE)
ev_tab <- data.frame(event = seq_len(dim(mue)[2]), sy_row = obj$dat$ev_sy,
                     initiator = labels[obj$dat$sy_state[obj$dat$ev_sy]],
                     target = labels[obj$dat$alt_target[obj$dat$ev_row]],
                     year = if (!is.null(obj$sy_keys)) obj$sy_keys$year[obj$dat$ev_sy] else NA,
                     ev_cls |> select(-state))
write.csv(ev_tab, file.path(out_dir, paste0("mu_event_", cell, ".csv")), row.names = FALSE)

## ---- APPENDIX: 9-cell classification (D = 2 only), reach, focus -----------------------------------------
band <- r0                                   # per-axis bands default to r0; override with --band_M / --band_R
band_M <- as.numeric(sub("^--band_M=", "", c(grep("^--band_M=", args, value = TRUE), paste0("--band_M=", band))[1]))
band_R <- as.numeric(sub("^--band_R=", "", c(grep("^--band_R=", args, value = TRUE), paste0("--band_R=", band))[1]))
if (dim(mu)[3] == 2) {
  c9 <- classify_cells9(mu, band_M = band_M, band_R = band_R, state_labels = labels)
  c9$n_events_as_initiator <- n_init
  write.csv(c9, file.path(out_dir, paste0("typology9_", cell, ".csv")), row.names = FALSE)
  logf("9-cell classification written (band_M = %.2f, band_R = %.2f)", band_M, band_R)
} else logf("D = 1: the single-axis 3-level classification in typology_%s.csv is the only typology output", cell)
rch <- reach_summary(mu, state_labels = labels); rch$n_events_as_initiator <- n_init
write.csv(rch, file.path(out_dir, paste0("reach_", cell, ".csv")), row.names = FALSE)
sig <- rstan::extract(fit, pars = "Sigma")$Sigma
if (length(dim(sig)) == 1) sig <- array(sig, c(length(sig), 1, 1))
write.csv(focus_sigma(sig, D_mode = obj$dat$D_mode), file.path(out_dir, paste0("focus_sigma_", cell, ".csv")), row.names = FALSE)
D <- dim(mu)[3]; xcols <- if (D == 2) 1:2 else obj$dat$D_mode
x_chosen <- obj$dat$x[obj$dat$ev_row, xcols, drop = FALSE]
mu_ev_mean <- apply(mue, c(2, 3), mean)
fs <- focus_states(x_chosen, mu_ev_mean, obj$dat$sy_state[obj$dat$ev_sy], labels, D_mode = obj$dat$D_mode)
write.csv(fs, file.path(out_dir, paste0("focus_states_", cell, ".csv")), row.names = FALSE)
logf("reach and focus (descriptive; pooled Sigma) written; %d initiating states have < 5 events", sum(fs$unreliable_lt_min_events))

## ---- AME: Procrustes alignment across chains, then pool -------------------------------------------------
if (!is.null(obj$dat$K)) {
  S <- obj$dat$S; K <- obj$dat$K
  Ud <- rstan::extract(fit, pars = "U", permuted = FALSE); Vd <- rstan::extract(fit, pars = "V", permuted = FALSE)
  U4 <- array(Ud, c(dim(Ud)[1], dim(Ud)[2], S, K)); V4 <- array(Vd, c(dim(Vd)[1], dim(Vd)[2], S, K))
  ## identification check before exporting u'v: U and V columns are centred in the Stan program
  ## (translation freedom U[,k] + h, b_j - h V[j,k] removed), so every draw's column means must be ~0
  stopifnot(max(abs(apply(U4, c(1, 2, 4), mean))) < 1e-8, max(abs(apply(V4, c(1, 2, 4), mean))) < 1e-8)
  al <- procrustes_align(U4, V4, ref_chain = 1)
  for (ch in seq_along(al$rotations)) logf("chain %d Procrustes rotation max |R - I| = %.3f", ch, max(abs(al$rotations[[ch]] - diag(K))))
  U_mean <- apply(al$U, c(3, 4), mean); V_mean <- apply(al$V, c(3, 4), mean)
  write.csv(data.frame(state = labels, U_mean, anchor = seq_len(S) %in% obj$dat$anchor), file.path(out_dir, paste0("ame_U_", cell, ".csv")), row.names = FALSE)
  write.csv(data.frame(state = labels, V_mean), file.path(out_dir, paste0("ame_V_", cell, ".csv")), row.names = FALSE)
  ## u_i'v_j posterior mean (rotation invariant), strongest pairs among initiators with events
  n_it <- dim(al$U)[1]; n_ch <- dim(al$U)[2]
  UV <- matrix(0, S, S)
  for (ch in seq_len(n_ch)) for (it in seq_len(n_it)) UV <- UV + matrix(al$U[it, ch, , ], S, K) %*% t(matrix(al$V[it, ch, , ], S, K))
  UV <- UV / (n_it * n_ch)
  pairs <- expand.grid(i = seq_len(S), j = seq_len(S)) |> filter(i != j, n_init[i] > 0) |>
    mutate(initiator = labels[i], target = labels[j], uv = UV[cbind(i, j)]) |> arrange(desc(abs(uv))) |> slice_head(n = 200)
  write.csv(pairs |> select(initiator, target, uv), file.path(out_dir, paste0("ame_pairs_", cell, ".csv")), row.names = FALSE)
  logf("AME extraction written (K = %d)", K)
}
logf("03_extract done")
