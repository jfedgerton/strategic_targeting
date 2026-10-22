## ============================================================================
## 03_extract.R -- post-process one fit: typology from posterior mu_i, AME alignment
##
## Usage: Rscript v2/R/03_extract.R --fit=v2/output/fit_<cell>.rds [--r0=0.5]
## Writes next to the fit:
##   typology_<cell>.csv   per-state posterior type probabilities, MAP type, r / theta summaries
##   mu_state_<cell>.csv   posterior mean and 90% interval of mu_i at x_ref
##   mu_event_<cell>.csv   per-event mu (initiator at its own X) with MAP type
##   ame_U_<cell>.csv / ame_V_<cell>.csv   (M2 only) Procrustes-aligned posterior means
##   ame_pairs_<cell>.csv  (M2 only) initiator x target u_i'v_j posterior mean, strongest 200 pairs
## r0 (peer radius, SD units of the within-year standardised positions) is a documented
## argument; 0.5 is the v2 default and a PI decision to confirm.
## ============================================================================
suppressPackageStartupMessages({ library(rstan); library(dplyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(this_file), "00_helpers.R"))
source(file.path(dirname(this_file), "typology.R"))
Sys.setenv(V2_ROOT = normalizePath(file.path(dirname(this_file), "..")))
args <- commandArgs(trailingOnly = TRUE)
fit_path <- sub("^--fit=", "", grep("^--fit=", args, value = TRUE))
r0 <- as.numeric(sub("^--r0=", "", c(grep("^--r0=", args, value = TRUE), "--r0=0.5")[1]))
stopifnot(length(fit_path) == 1, file.exists(fit_path), is.finite(r0), r0 > 0)
obj <- readRDS(fit_path); fit <- obj$fit; cell <- obj$cell
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
write.csv(cls, file.path(out_dir, paste0("typology_", cell, ".csv")), row.names = FALSE)
mu_cols <- intersect(c("mu1_mean", "mu1_lo", "mu1_hi", "mu2_mean", "mu2_lo", "mu2_hi", "mu_mean", "mu_lo", "mu_hi"), names(cls))
write.csv(cls |> select(state, all_of(mu_cols), n_events_as_initiator),
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
