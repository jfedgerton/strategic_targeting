## ============================================================================
## sim_context_b.R -- simulation recovery of the M1c model (dual typology stage b, PI-approved spec 2026-10-27) on the OBSERVED design,
##                    ONE replicate per call (audit fix 5, 2026-10-27: a clean recovery test; the first version kept the observed y1)
##
## Usage: Rscript v2/sim/sim_context_b.R [--rep=0] [--n_init=250] [--n_noninit=250] [--chains=2] [--iter=1000] [--warmup=500]
##        [--sample=fmid] [--lag=3] [--pilot=FALSE] [--smoke=FALSE] [--model=M1c|M1b]
##   --model=M1b  (PI request 2026-10-27, after the M1c recovery bias) runs the SAME recovery test on the main-specification model M1b
##                (stan/m1_latent_additive.stan, joint "b" variant): the design is the real <sample> lag-<lag> M1b cell
##                (assemble_cell_spec(model = "M1b"): W without multi_mid_jan, P1 = 8, no groups), the subsample keeps --n_init observed
##                initiating state-years at random (the leader-first rule has no groups to apply to) plus --n_noninit non-initiators, the
##                truth is the M1c BILATERAL truth trimmed to the M1b dimensions (mu0 = (-1.0, 0.3), sigma_x = (1.2, 0.9), rho = 0.4, the
##                same beta, gamma, alpha, tau, sd_re, R_Om, sigma_b), y1 and every event's target come from the single softmax over the
##                full in_s2 set, the slot patterns are the observed (weight) multisets of the kept initiators, and the fit uses the
##                compiled m1_latent_additive program with the M1 sampler settings. Outputs never collide with the M1c files:
##                v2/sim/results/m1b_recovery[_smoke]/m1b_recovery_r<r>.csv and m1b_design_r<r>.csv, v2/sim/fits/m1b_recovery[_smoke]/m1b_sim_r<r>.rds;
##                `Rscript v2/sim/summarise_context_b.R --model=M1b` aggregates them; slurm/07_sim_m1b.sbatch runs 20 replicates. Recovery
##                blocks: stage_2_core (mu0, sigma_x, rho_x, Sigma; the M1b counterpart of the M1c stage_b block), stage_2 (beta, gamma),
##                stage_1 (alpha0, alpha, tau), re_scales (sd_re, sigma_b), mu_state. The default --model=M1c is unchanged.
##   --rep=r      replicate index; the seed is V2_SEED + r (123 + r) and every output carries _r<r>. slurm/06_sim_context_b.sbatch runs
##                r = 0..19 as an array; sim/summarise_context_b.R turns the replicates into a coverage table.
##   --smoke=TRUE writes to v2/sim/results/context_b_smoke/ instead of context_b/, so a short local check never mixes with the Roar
##                replicates (the summariser refuses replicates run with different settings).
##
## Design (NOT simulated) = the real <sample> lag-<lag> M1c cell (assemble_cell_spec(model = "M1c") on the manifest-checked data): the
## state-years, their real choice sets (sizes, positions x, dyadic controls Z, stage-1 controls W, in_iv / in_s2 flags), S = 217 and the
## empirical distribution of EVENT-SLOT PATTERNS (the multiset of (group, D8 weight) of an initiating state-year's usable events). The
## cell is SUBSAMPLED to fit a 2-core / 8 GB container: the observed initiating state-years with a coalition-leader event are kept first
## (all of them when they fit in --n_init, else a random draw), then random other observed initiating state-years up to --n_init in
## total, plus --n_noninit random observed non-initiating state-years. The observed y1 is used ONLY to choose this subsample and the
## target initiation rate; it never enters the simulated data.
## SIMULATED, from the M1c model at KNOWN parameters (TRUTH below), following simulate_targeting() (sim/simulate.R):
##   random effects  (c_i, a_i) ~ MVN(0, Omega) per state, b_j ~ sigma_b * centred N(0, 1)
##   stage 1         IV_n = log-sum-exp of the BILATERAL utility over the in_iv rows; eta_n = alpha0 + W_n alpha + c_i + tau IV_n with
##                   alpha0 calibrated (uniroot) so that the mean initiation probability equals the kept subsample's observed initiation
##                   rate; y1_n ~ Bernoulli(inv_logit(eta_n))
##   event slots     every simulated initiator draws one slot pattern from the pool of the kept observed initiators' patterns (with
##                   replacement; a pattern with more slots than the state-year's stage-2 rows minus one is redrawn)
##   stage 2         each slot's target is drawn independently from the softmax of its GROUP's utility (mu0 + delta_g, Sigma_g) over the
##                   state-year's FULL in_s2 set (with replacement, matching the likelihood, which scores every event against the complete set)
## Recovery targets: EVERY parameter the fit reports (mu0, delta_lead, delta_join, mu0 + delta_g, sigma_g (3 x 2), rho, beta, gamma,
## alpha0, alpha, tau, sd_re (3), sigma_b) plus mu_state (bilateral ideal points of the states with events: correlation and coverage).
## Output: v2/sim/results/context_b/context_b_recovery_r<r>.csv (parameter, truth, posterior mean, sd, 5 / 95%, covered, bias, z, Rhat,
##         n_eff, rep, seed, run settings), context_b_design_r<r>.csv (the realised design counts), v2/sim/fits/context_b[_smoke]/context_b_sim_r<r>.rds
##         (gitignored), a log in v2/logs. The fit uses the compiled m1_context program (compile_models.R) with the default M1 adapt_delta
##         / max_treedepth. Smoke check (a 2-core container, ~1 h): --smoke=TRUE --n_init=160 --n_noninit=140 --chains=2 --iter=500 --warmup=250.
## ============================================================================
suppressPackageStartupMessages({ library(rstan); library(dplyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
sim_dir <- normalizePath(dirname(this_file)); R_DIR <- file.path(sim_dir, "..", "R")
Sys.setenv(V2_ROOT = normalizePath(file.path(sim_dir, "..")))
source(file.path(R_DIR, "spec.R")); source(file.path(R_DIR, "helpers_data.R")); source(file.path(R_DIR, "helpers_model.R"))
source(file.path(R_DIR, "assemble_cell.R")); source(file.path(R_DIR, "typology.R")); source(file.path(R_DIR, "context_model.R"))
args <- commandArgs(trailingOnly = TRUE)
REP <- as.integer(get_arg(args, "rep", 0)); SEED <- V2_SEED + REP
N_INIT <- as.integer(get_arg(args, "n_init", 250)); N_NONINIT <- as.integer(get_arg(args, "n_noninit", 250))
CHAINS <- as.integer(get_arg(args, "chains", 2)); ITER <- as.integer(get_arg(args, "iter", 1000)); WARMUP <- as.integer(get_arg(args, "warmup", 500))
SAMPLE <- get_arg(args, "sample", "fmid"); LAG <- as.integer(get_arg(args, "lag", 3)); PILOT <- get_arg(args, "pilot", FALSE); SMOKE <- get_arg(args, "smoke", FALSE)
MODEL <- get_arg(args, "model", "M1c")
if (length(grep("^--seed=", args))) stop("--seed is not an option: the seed is V2_SEED + --rep")
if (!MODEL %in% c("M1c", "M1b")) stop("--model must be M1c (default) or M1b")
stopifnot(REP >= 0, N_INIT >= 10, N_NONINIT >= 0, CHAINS >= 1, ITER > WARMUP, SAMPLE %in% ALL_SAMPLES, LAG %in% LAGS)
M1C <- MODEL == "M1c"
log_open(if (M1C) sprintf("sim_context_b_r%d", REP) else sprintf("sim_m1b_r%d", REP)); set.seed(SEED)
rstan_options(auto_write = FALSE); options(mc.cores = min(CHAINS, n_cores()))
if (M1C) {
  fit_dir <- file.path(sim_dir, "fits", if (SMOKE) "context_b_smoke" else "context_b"); dir.create(fit_dir, showWarnings = FALSE, recursive = TRUE); res_dir <- file.path(sim_dir, "results", if (SMOKE) "context_b_smoke" else "context_b"); dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
  prefix <- "context_b"; fit_stem <- "context_b_sim"
} else {
  fit_dir <- file.path(sim_dir, "fits", if (SMOKE) "m1b_recovery_smoke" else "m1b_recovery"); dir.create(fit_dir, showWarnings = FALSE, recursive = TRUE); res_dir <- file.path(sim_dir, "results", if (SMOKE) "m1b_recovery_smoke" else "m1b_recovery"); dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
  prefix <- "m1b"; fit_stem <- "m1b_sim"
}
tag <- sprintf("r%d", REP)
logf("%s simulation recovery on the observed %s lag-%d design, replicate %d (seed %d): n_init %d, n_noninit %d, chains %d, iter %d, warmup %d%s", MODEL, SAMPLE, LAG, REP, SEED, N_INIT, N_NONINIT, CHAINS, ITER, WARMUP, if (SMOKE) sprintf(" (SMOKE: results go to results/%s/)", basename(res_dir)) else "")

## ---- the observed design -----------------------------------------------------------------------------------------------------------
A <- assemble_cell_spec(SAMPLE, LAG, parse_spec(character(0)), model = MODEL)   # M1c appends multi_mid_jan to W; M1b is the main-grid assembly
sy <- A$sy; alts <- A$alts; ev <- A$ev
stopifnot("grp" %in% names(ev), all(ev$grp %in% 0:2))                           # the groups are computed for every model (assemble_cell.R)
if (M1C) stopifnot(tail(A$W_cols, 1) == "multi_mid_jan") else stopifnot(!"multi_mid_jan" %in% A$W_cols, identical(A$W_cols, c("dM_i", "dR_i", w_controls())))
S <- length(A$state_codes)
init_sy <- sort(unique(ev$sy_row)); lead_sy <- sort(unique(ev$sy_row[ev$grp == 1L]))
if (M1C) {
  keep_lead <- if (length(lead_sy) <= N_INIT) lead_sy else sort(sample(lead_sy, N_INIT))
  keep_init <- sort(c(keep_lead, sample(setdiff(init_sy, keep_lead), N_INIT - length(keep_lead))))
} else {
  ## M1b has no groups, so the leader-first rule does not apply: N_INIT observed initiating state-years at random
  stopifnot(length(init_sy) >= N_INIT)
  keep_lead <- integer(0); keep_init <- sort(sample(init_sy, N_INIT))
}
non_init <- setdiff(seq_len(nrow(sy)), init_sy)
keep_non <- sort(sample(non_init, min(N_NONINIT, length(non_init))))
keep <- sort(c(keep_init, keep_non))
logf("design: %d observed initiating state-years with usable events (%d with a leader event, %d kept), %d initiating kept; %d of %d non-initiating state-years kept; N1 = %d",
     length(init_sy), length(lead_sy), length(keep_lead), length(keep_init), length(keep_non), length(non_init), length(keep))
## subset the frames: sy rows, their alternative blocks, their events (the slot patterns); re-index sy_row / alt_row
new_sy <- match(seq_len(nrow(sy)), keep)                     # old sy_row -> new (NA when dropped)
sy_s <- sy[keep, ]
alts_s <- alts |> mutate(.old_row = row_number()) |> filter(sy_row %in% keep) |> mutate(sy_row = new_sy[sy_row])
ev_s <- ev |> filter(sy_row %in% keep) |> mutate(sy_row = new_sy[sy_row], alt_row = match(alt_row, alts_s$.old_row))
stopifnot(!anyNA(ev_s$alt_row), !anyNA(ev_s$sy_row), all(alts_s$sy_row[ev_s$alt_row] == ev_s$sy_row), all(diff(alts_s$sy_row) >= 0))
alts_s <- alts_s |> select(-.old_row)
N1 <- nrow(sy_s)
blk <- split(seq_len(nrow(alts_s)), alts_s$sy_row); stopifnot(length(blk) == N1, all(names(blk) == as.character(seq_len(N1))))
n_s2 <- vapply(blk, function(r) sum(alts_s$in_s2[r]), 1L); n_iv <- vapply(blk, function(r) sum(alts_s$in_iv[r]), 1L)
stopifnot(all(n_iv >= 1L), all(n_s2 >= 2L))                  # every kept state-year can carry at least one event
## the slot-pattern pool: the kept observed initiators' (group, weight) multisets
pool <- lapply(split(ev_s[, c("grp", "w")], ev_s$sy_row), function(d) d[order(d$grp, d$w), ])
stopifnot(length(pool) == length(keep_init))
pool_size <- vapply(pool, nrow, 1L)
obs_rate <- length(keep_init) / N1
if (M1C) { logf("slot-pattern pool: %d patterns (%d distinct), %d slots (bilateral %d, leader %d, joiner %d), slots per pattern median %d max %d; target initiation rate %.3f (the kept subsample's observed rate)",
     length(pool), n_distinct(vapply(pool, function(d) paste(d$grp, round(d$w, 4), collapse = " "), "")), sum(pool_size), sum(ev_s$grp == 0L), sum(ev_s$grp == 1L), sum(ev_s$grp == 2L),
     as.integer(median(pool_size)), max(pool_size), obs_rate)
} else { logf("slot-pattern pool (M1b: weight multisets, the group labels play no part): %d patterns (%d distinct), %d slots, slots per pattern median %d max %d; target initiation rate %.3f (the kept subsample's observed rate)",
     length(pool), n_distinct(vapply(pool, function(d) paste(round(d$w, 4), collapse = " "), "")), sum(pool_size), as.integer(median(pool_size)), max(pool_size), obs_rate) }
x_ref <- colMeans(as.matrix(sy_s[sy_s$y1 == 1L, A$X_cols, drop = FALSE]))     # the observed-initiator mean of X, kept as the design's x_ref
dat0 <- build_stan_data(sy_s, alts_s, ev_s, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, x_ref = x_ref, allow_orphan_initiations = TRUE, D_mode = 3L, S = S)
stopifnot(dat0$N1 == N1, dat0$P == 2L)

## ---- TRUTH (scaled gap units; one unit = a 26-fold CINC ratio / a 1.7-fold PageRank ratio) ------------------------------------------
set.seed(SEED + 1L)
TRUTH <- list(mu0 = c(-1.0, 0.3), delta_lead = c(-1.2, 0.8), delta_join = c(0.9, -0.5),
              sigma_g = rbind(bilateral = c(1.2, 0.9), leader = c(2.2, 1.4), joiner = c(1.5, 2.0)), rho = 0.4,
              beta = matrix(c(0.3, -0.2, 0.1, 0.4), 2, 2), gamma = c(1.0, -0.3, 0.4, 0.8, 0.1, -0.1, 0.2)[seq_len(dat0$Pz)],
              alpha = c(0.2, -0.1, 0.5, 0.3, 0.4, 0.3, -0.2, 0.1, 0.3)[seq_len(dat0$P1)], tau = 0.5,
              sd_re = c(0.9, 0.7, 0.6), R_Om = matrix(c(1, 0.3, 0.1, 0.3, 1, 0.2, 0.1, 0.2, 1), 3), sigma_b = 0.6)
stopifnot(length(TRUTH$gamma) == dat0$Pz, length(TRUTH$alpha) == dat0$P1)
if (!M1C) {
  ## M1b: one ideal point and one Sigma = the M1c bilateral truth; no shifts, no group scales (alpha has P1 = 8 entries: no multi_mid_jan)
  stopifnot(dat0$P1 == 8L)
  TRUTH$delta_lead <- NULL; TRUTH$delta_join <- NULL; TRUTH$sigma_x <- unname(TRUTH$sigma_g["bilateral", ]); TRUTH$sigma_g <- NULL
  TRUTH$Sigma <- matrix(c(TRUTH$sigma_x[1]^2, TRUTH$rho * prod(TRUTH$sigma_x), TRUTH$rho * prod(TRUTH$sigma_x), TRUTH$sigma_x[2]^2), 2)
}
L_true <- t(chol(TRUTH$R_Om))
z_re <- matrix(rnorm(3 * S), 3, S)
re <- t(diag(TRUTH$sd_re) %*% L_true %*% z_re)                 # S x 3: (c_i, a_i1, a_i2)
TRUTH$c_i <- re[, 1]; TRUTH$a_i <- re[, 2:3]
z_b <- rnorm(S); TRUTH$b_j <- TRUTH$sigma_b * (z_b - mean(z_b))
TRUTH$mu_state <- rep(1, S) %o% TRUTH$mu0 + rep(1, S) %o% as.numeric(x_ref %*% TRUTH$beta) + TRUTH$a_i

## ---- utilities at the truth ---------------------------------------------------------------------------------------------------------
Mu <- rep(1, N1) %o% TRUTH$mu0 + dat0$X %*% TRUTH$beta + TRUTH$a_i[dat0$sy_state, , drop = FALSE]
zg <- as.numeric(dat0$Z %*% TRUTH$gamma)
if (M1C) {
  shift <- rbind(c(0, 0), TRUTH$delta_lead, TRUTH$delta_join)
  Sinv <- lapply(1:3, function(g) solve(matrix(c(TRUTH$sigma_g[g, 1]^2, TRUTH$rho * prod(TRUTH$sigma_g[g, ]), TRUTH$rho * prod(TRUTH$sigma_g[g, ]), TRUTH$sigma_g[g, 2]^2), 2)))
} else {
  shift <- rbind(c(0, 0)); Sinv <- list(solve(TRUTH$Sigma))                   # M1b: every event uses group 0 (the single utility)
}
util_rows <- function(rows, n, g) {
  d <- sweep(dat0$x[rows, , drop = FALSE], 2, Mu[n, ] + shift[g + 1L, ])
  -0.5 * rowSums((d %*% Sinv[[g + 1L]]) * d) + TRUTH$b_j[dat0$alt_target[rows]] + zg[rows]
}
## ---- stage 1: the inclusive value, alpha0 calibrated to the observed rate, y1 drawn (simulate_targeting() logic) ----------------------
IV <- vapply(seq_len(N1), function(n) { rows <- dat0$cs_start[n]:dat0$cs_end[n]; rows <- rows[dat0$in_iv[rows] == 1L]; u <- util_rows(rows, n, 0L); m <- max(u); m + log(sum(exp(u - m))) }, numeric(1))
stopifnot(all(is.finite(IV)))
lin <- as.numeric(dat0$W %*% TRUTH$alpha) + TRUTH$c_i[dat0$sy_state] + TRUTH$tau * IV
TRUTH$alpha0 <- uniroot(function(a0) mean(plogis(a0 + lin)) - obs_rate, c(-40, 40), tol = 1e-10)$root
set.seed(SEED + 2L)
y1_sim <- rbinom(N1, 1L, plogis(TRUTH$alpha0 + lin))
if (sum(y1_sim) < 5L) stop("fewer than 5 simulated initiators; raise --n_init")
logf("stage 1 simulated: IV mean %.3f SD %.3f; alpha0 calibrated to %.3f; %d of %d state-years initiate (rate %.3f; observed initiators among them %d)",
     mean(IV), sd(IV), TRUTH$alpha0, sum(y1_sim), N1, mean(y1_sim), sum(y1_sim == 1L & sy_s$y1 == 1L))
sy_sim <- sy_s; sy_sim$y1 <- as.integer(y1_sim)
## ---- event slots and stage 2: a pattern per simulated initiator; each event drawn independently from its group softmax over the FULL stage-2 set
## (with replacement, as the likelihood scores every event against the complete set; Codex review 2026-10-27, P2)
set.seed(SEED + 3L)
ev_list <- vector("list", N1); n_redraw <- 0L
for (n in which(y1_sim == 1L)) {
  rows_all <- dat0$cs_start[n]:dat0$cs_end[n]; rows_all <- rows_all[dat0$in_s2[rows_all] == 1L]
  stopifnot(length(rows_all) >= 2L)
  pat <- pool[[sample.int(length(pool), 1L)]]
  picks <- integer(nrow(pat))
  for (k in seq_len(nrow(pat))) {
    u <- util_rows(rows_all, n, if (M1C) pat$grp[k] else 0L); p <- exp(u - max(u)); p <- p / sum(p)   # M1b: the single softmax for every slot
    picks[k] <- rows_all[sample.int(length(rows_all), 1L, prob = p)]   # full set every slot: no removal of earlier picks
  }
  ev_list[[n]] <- data.frame(sy_row = n, alt_row = picks, grp = pat$grp, w = pat$w)
}
ev_sim <- bind_rows(ev_list)
stopifnot(nrow(ev_sim) >= 1L, all(dat0$in_s2[ev_sim$alt_row] == 1L), all(alts_s$sy_row[ev_sim$alt_row] == ev_sim$sy_row),
          all(sy_sim$y1[ev_sim$sy_row] == 1L), setequal(unique(ev_sim$sy_row), which(y1_sim == 1L)), all(ev_sim$w > 0), all(ev_sim$w <= 1))
n_dup_events <- sum(duplicated(paste(ev_sim$sy_row, ev_sim$alt_row)))   # repeated (state-year, target) pairs: allowed under independent draws, logged
logf("stage 2 simulated: %d events; %d repeat an earlier (state-year, target) pair", nrow(ev_sim), n_dup_events)
if (M1C) for (g in 0:2) if (sum(ev_sim$grp == g) == 0L) stop("no simulated event in group ", STAGE2_GROUPS[g + 1L], "; raise --n_init")
if (!M1C) ev_sim$grp <- NULL                                  # M1b: the observed group label of a slot is not part of the simulated data
dat_sim <- build_stan_data(sy_sim, alts_s, ev_sim, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, x_ref = x_ref, allow_orphan_initiations = FALSE, D_mode = 3L, S = S)
if (M1C) dat_sim$ev_grp <- as.integer(ev_sim$grp)
stopifnot(identical(dat_sim$W, dat0$W), identical(dat_sim$Z, dat0$Z), identical(dat_sim$x, dat0$x), identical(dat_sim$in_iv, dat0$in_iv), identical(dat_sim$in_s2, dat0$in_s2),
          sum(dat_sim$y1) == sum(y1_sim), dat_sim$E == nrow(ev_sim))
design <- if (M1C) { data.frame(group = STAGE2_GROUPS, n_events = as.integer(tabulate(dat_sim$ev_grp + 1L, nbins = 3L)),
                     w_sum = vapply(0:2, function(g) sum(dat_sim$ev_w[dat_sim$ev_grp == g]), numeric(1)),
                     n_states = vapply(0:2, function(g) n_distinct(dat_sim$sy_state[dat_sim$ev_sy[dat_sim$ev_grp == g]]), 1L),
                     n_state_years = vapply(0:2, function(g) n_distinct(dat_sim$ev_sy[dat_sim$ev_grp == g]), 1L)) } else { data.frame(group = "all", n_events = dat_sim$E, w_sum = sum(dat_sim$ev_w), n_states = n_distinct(dat_sim$sy_state[dat_sim$ev_sy]), n_state_years = n_distinct(dat_sim$ev_sy)) }
logf("Stan data (simulated y1 and events on the observed design): S=%d N1=%d N2=%d E=%d P1=%d P=%d Pz=%d; sum(y1)=%d; groups: %s; choice sets median %d (min %d, max %d); pattern redraws %d",
     dat_sim$S, dat_sim$N1, dat_sim$N2, dat_sim$E, dat_sim$P1, dat_sim$P, dat_sim$Pz, sum(dat_sim$y1), paste(sprintf("%s %d events (w %.1f, %d states)", design$group, design$n_events, design$w_sum, design$n_states), collapse = "; "),
     as.integer(median(dat_sim$cs_end - dat_sim$cs_start + 1L)), min(dat_sim$cs_end - dat_sim$cs_start + 1L), max(dat_sim$cs_end - dat_sim$cs_start + 1L), n_redraw)
## the simulated targets must favour the group's ideal point: the mean chosen gap per group against the truth
chosen <- dat_sim$x[dat_sim$ev_row, , drop = FALSE]
if (M1C) { for (g in 0:2) logf("simulated chosen gaps, group %s: mean (%.2f, %.2f) over %d events; truth population ideal point (%.2f, %.2f)", STAGE2_GROUPS[g + 1L],
                    mean(chosen[dat_sim$ev_grp == g, 1]), mean(chosen[dat_sim$ev_grp == g, 2]), sum(dat_sim$ev_grp == g), TRUTH$mu0[1] + shift[g + 1L, 1], TRUTH$mu0[2] + shift[g + 1L, 2]) } else { logf("simulated chosen gaps: mean (%.2f, %.2f) over %d events; truth population ideal point (%.2f, %.2f)", mean(chosen[, 1]), mean(chosen[, 2]), dat_sim$E, TRUTH$mu0[1], TRUTH$mu0[2]) }
## the brute-force likelihood at the truth is finite (sanity)
R_true <- if (M1C) { stage2_loglik_R_context(dat_sim, list(mu0 = TRUTH$mu0, delta_lead = TRUTH$delta_lead, delta_join = TRUTH$delta_join, beta = TRUTH$beta, a_i = TRUTH$a_i, b_j = TRUTH$b_j, gamma = TRUTH$gamma, sigma_g = unname(TRUTH$sigma_g), rho = TRUTH$rho)) } else { stage2_loglik_R(dat_sim, list(mu0 = TRUTH$mu0, beta = TRUTH$beta, a_i = TRUTH$a_i, b_j = TRUTH$b_j, gamma = TRUTH$gamma, sigma_x = TRUTH$sigma_x, rho = TRUTH$rho)) }
stopifnot(is.finite(R_true$lp_stage2))
logf("stage-2 weighted log-likelihood at the truth: %.1f (mean per event %.3f; uniform choice would give %.3f)", R_true$lp_stage2, mean(R_true$loglik), mean(-log(dat_sim$cs_end - dat_sim$cs_start + 1)))

## ---- fit ---------------------------------------------------------------------------------------------------------------------------
family <- stan_program(MODEL)                                  # m1_context (M1c) | m1_latent_additive (M1b, the main grid's program)
sm <- readRDS(file.path(v2_root(), "stan", "build", paste0(family, ".rds")))
assert_compiled_current(sm, file.path(v2_root(), "stan", paste0(family, ".stan")), file.path(v2_root(), "stan", "build", paste0(family, ".sha256")))
if (PILOT) { CHAINS <- 1L; ITER <- 40L; WARMUP <- 20L; logf("PILOT: 1 chain x 40 iterations (timing only)") }
dat_sim$grainsize <- dat_sim$N1
rstan_options(threads_per_chain = 1)
t0 <- Sys.time()
fit <- sampling(sm, data = dat_sim, chains = CHAINS, iter = ITER, warmup = WARMUP, seed = SEED, refresh = max(1L, ITER %/% 10L),
                control = list(adapt_delta = SAMPLER$M1$adapt_delta, max_treedepth = SAMPLER$M1$max_treedepth),
                pars = c("z_re", "z_b", "L_Om", "re", "log_lik_stage1", "log_lik_stage2", "mu_event", if (!M1C) "L_Om_joint"), include = FALSE)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
if (fit@mode != 0L) stop("sampling did not run")
assert_fit_chains(fit, CHAINS, ITER, WARMUP)
sp <- get_sampler_params(fit, inc_warmup = FALSE)
n_div <- sum(sapply(sp, function(m) sum(m[, "divergent__"]))); n_td <- sum(sapply(sp, function(m) sum(m[, "treedepth__"] >= SAMPLER$M1$max_treedepth)))
logf("sampling finished in %.1f min; divergences %d; max-treedepth hits %d", elapsed, n_div, n_td)
if (PILOT) { logf("PILOT done; no recovery table"); quit(save = "no", status = 0) }

## ---- recovery table ---------------------------------------------------------------------------------------------------------------
## M1b (2026-10-27): sigma_x[1], sigma_x[2], rho_x[1] as the program names them, plus the generated Sigma (Sigma[1,2] is the PI's question;
## the mirror entry Sigma[2,1] is dropped) in the block stage_2_core
pars <- if (M1C) { c("mu0", "delta_lead", "delta_join", "mu0_lead", "mu0_join", "sigma_g", "rho_x", "beta", "gamma", "tau", "alpha0", "alpha", "sd_re", "sigma_b") } else { c("mu0", "sigma_x", "rho_x", "Sigma", "beta", "gamma", "tau", "alpha0", "alpha", "sd_re", "sigma_b") }
sm_tab <- summary(fit, pars = pars, probs = c(0.05, 0.95))$summary
if (!M1C) { stopifnot(all(c("sigma_x[1]", "sigma_x[2]", "rho_x[1]", "Sigma[1,2]", "Sigma[2,1]") %in% rownames(sm_tab)), !any(grepl("^(delta_|sigma_g|mu0_lead|mu0_join)", rownames(sm_tab)))); sm_tab <- sm_tab[rownames(sm_tab) != "Sigma[2,1]", , drop = FALSE] }
truth_of <- function(nm) {
  base <- sub("\\[.*$", "", nm); i1 <- suppressWarnings(as.integer(sub("^[^\\[]+\\[([0-9]+).*$", "\\1", nm))); i2 <- suppressWarnings(as.integer(sub("^[^\\[]+\\[[0-9]+,([0-9]+)\\]$", "\\1", nm)))
  switch(base, mu0 = TRUTH$mu0[i1], delta_lead = TRUTH$delta_lead[i1], delta_join = TRUTH$delta_join[i1], mu0_lead = TRUTH$mu0[i1] + TRUTH$delta_lead[i1],
         mu0_join = TRUTH$mu0[i1] + TRUTH$delta_join[i1], sigma_g = TRUTH$sigma_g[i1, i2], sigma_x = TRUTH$sigma_x[i1], Sigma = TRUTH$Sigma[i1, i2], rho_x = TRUTH$rho, beta = TRUTH$beta[i1, i2], gamma = TRUTH$gamma[i1],
         tau = TRUTH$tau, alpha0 = TRUTH$alpha0, alpha = TRUTH$alpha[i1], sd_re = TRUTH$sd_re[i1], sigma_b = TRUTH$sigma_b, NA_real_)
}
rec <- data.frame(parameter = rownames(sm_tab), truth = vapply(rownames(sm_tab), truth_of, numeric(1)), mean = sm_tab[, "mean"], sd = sm_tab[, "sd"],
                  lo = sm_tab[, "5%"], hi = sm_tab[, "95%"], Rhat = sm_tab[, "Rhat"], n_eff = sm_tab[, "n_eff"], stringsAsFactors = FALSE)
stopifnot(!anyNA(rec$truth))
rec$covered <- rec$truth >= rec$lo & rec$truth <= rec$hi
rec$bias <- rec$mean - rec$truth
rec$z <- rec$bias / rec$sd
rec$block <- if (M1C) { case_when(grepl("^(mu0|delta_lead|delta_join|mu0_lead|mu0_join|sigma_g|rho_x)", rec$parameter) ~ "stage_b", grepl("^(beta|gamma)", rec$parameter) ~ "stage_2",
                       grepl("^(tau|alpha0|alpha)", rec$parameter) ~ "stage_1", grepl("^(sd_re|sigma_b)", rec$parameter) ~ "re_scales", TRUE ~ "other") } else { case_when(grepl("^(mu0|sigma_x|rho_x|Sigma)", rec$parameter) ~ "stage_2_core", grepl("^(beta|gamma)", rec$parameter) ~ "stage_2",
                       grepl("^(tau|alpha0|alpha)", rec$parameter) ~ "stage_1", grepl("^(sd_re|sigma_b)", rec$parameter) ~ "re_scales", TRUE ~ "other") }
stopifnot(all(rec$block != "other"))
core_block <- if (M1C) "stage_b" else "stage_2_core"
rec$recovery_target <- TRUE                                   # y1, the slots and the targets are all simulated from the model: every parameter is a target
rownames(rec) <- NULL
## mu_state (bilateral ideal points) of the states with events: correlation of the posterior mean with the truth, coverage per state
ms <- rstan::extract(fit, pars = "mu_state")$mu_state
with_ev <- sort(unique(dat_sim$sy_state[dat_sim$ev_sy]))
ms_mean <- apply(ms, c(2, 3), mean); ms_lo <- apply(ms, c(2, 3), quantile, 0.05); ms_hi <- apply(ms, c(2, 3), quantile, 0.95)
cov_ms <- mean(TRUTH$mu_state[with_ev, ] >= ms_lo[with_ev, ] & TRUTH$mu_state[with_ev, ] <= ms_hi[with_ev, ])
cor_ms <- c(cor(ms_mean[with_ev, 1], TRUTH$mu_state[with_ev, 1]), cor(ms_mean[with_ev, 2], TRUTH$mu_state[with_ev, 2]))
rec <- bind_rows(rec, data.frame(parameter = c("mu_state[material] cor(post mean, truth) over states with events", "mu_state[relational] cor(post mean, truth) over states with events", "mu_state 90% coverage over states with events x axes"),
                                 truth = c(1, 1, 0.9), mean = c(cor_ms, cov_ms), sd = NA, lo = NA, hi = NA, Rhat = NA, n_eff = NA, covered = NA, bias = NA, z = NA, block = "mu_state", recovery_target = TRUE))
rec$rep <- REP; rec$seed <- SEED; rec$n_init <- N_INIT; rec$n_noninit <- N_NONINIT; rec$chains <- CHAINS; rec$iter <- ITER; rec$warmup <- WARMUP
rec$elapsed_min <- elapsed; rec$divergences <- n_div; rec$treedepth_hits <- n_td; rec$sample <- SAMPLE; rec$lag <- LAG
rec_file <- file.path(res_dir, sprintf("%s_recovery_%s.csv", prefix, tag)); dsn_file <- file.path(res_dir, sprintf("%s_design_%s.csv", prefix, tag))
if (M1C) stopifnot(rec_file == file.path(res_dir, sprintf("context_b_recovery_%s.csv", tag)), dsn_file == file.path(res_dir, sprintf("context_b_design_%s.csv", tag)))
write.csv(rec, rec_file, row.names = FALSE)
write.csv(design |> mutate(rep = REP, seed = SEED, N1 = dat_sim$N1, N2 = dat_sim$N2, E = dat_sim$E, S = dat_sim$S, sum_y1 = sum(dat_sim$y1), alpha0_true = TRUTH$alpha0,
                           n_states_with_events = length(with_ev), n_pattern_redraws = n_redraw, n_dup_events = n_dup_events),
          dsn_file, row.names = FALSE)
fit_file <- if (M1C) file.path(fit_dir, sprintf("context_b_sim_%s.rds", tag)) else file.path(fit_dir, sprintf("m1b_sim_%s.rds", tag))
saveRDS(list(fit = fit, dat = dat_sim, truth = TRUTH, design = design, keep = keep, args = list(rep = REP, seed = SEED, n_init = N_INIT, n_noninit = N_NONINIT, chains = CHAINS, iter = ITER, warmup = WARMUP, sample = SAMPLE, lag = LAG, model = MODEL)),
        fit_file)   # smoke fits never overwrite production fits (Codex review 2026-10-27, P2); M1b fits have their own directory
tg <- rec |> filter(!is.na(covered))
logf("recovery over the %d parameters: 90%% coverage %.2f, mean |bias| %.3f, max |z| %.2f, max Rhat %.3f, min n_eff %.0f", nrow(tg), mean(tg$covered), mean(abs(tg$bias)), max(abs(tg$z)), max(tg$Rhat), min(tg$n_eff))
for (b in c(core_block, "stage_2", "stage_1", "re_scales")) { tb <- tg |> filter(block == b); logf("  block %-12s %2d parameters: coverage %.2f, mean |bias| %.3f, max |z| %.2f", b, nrow(tb), mean(tb$covered), mean(abs(tb$bias)), max(abs(tb$z))) }
key <- rec |> filter(block %in% c(core_block, "stage_1", "re_scales"), !grepl("^alpha\\[", parameter))
for (r in seq_len(nrow(key))) logf("  %-14s truth %6.2f  mean %6.2f  [%6.2f, %6.2f]  covered %-5s  bias %6.2f  Rhat %.3f", key$parameter[r], key$truth[r], key$mean[r], key$lo[r], key$hi[r], key$covered[r], key$bias[r], key$Rhat[r])
logf("mu_state recovery (states with events, n = %d): cor material %.3f, relational %.3f; 90%% coverage %.2f", length(with_ev), cor_ms[1], cor_ms[2], cov_ms)
logf("written: %s, %s, %s", rec_file, dsn_file, fit_file)
