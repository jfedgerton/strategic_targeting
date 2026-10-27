## ============================================================================
## sim_context_b.R -- simulation recovery of the M1c model (dual typology stage b, PI-approved spec 2026-10-27) on the OBSERVED design
##
## Usage: Rscript v2/sim/sim_context_b.R [--n_init=250] [--n_noninit=250] [--chains=2] [--iter=1000] [--warmup=500] [--seed=123]
##        [--sample=fmid] [--lag=3] [--pilot=FALSE]
##
## Design = the real fmid lag-3 M1c cell (assemble_cell_spec(model = "M1c") on the manifest-checked data): the real state-years, their
## real choice sets (sizes, positions x, dyadic controls Z, W, in_iv / in_s2 flags), the real event slots (which state-year, how many
## events, in which stage-b group, with which D8 weight) and the real stage-1 outcome y1. The cell is SUBSAMPLED to fit a 2-core / 8 GB
## container: every initiating state-year with a coalition-leader event is kept, plus a random draw of the other initiating
## state-years up to --n_init in total, plus --n_noninit random non-initiating state-years (the stage-1 rows). The state index space
## (S = 217) is kept. What is SIMULATED: the chosen target of every event slot, drawn from the M1c stage-2 softmax of its group over the
## state-year's in_s2 rows (without replacement within a state-year) at KNOWN parameters (TRUTH below). What is NOT simulated: y1 and
## the event slots (the observed design; the stage-1 parameters alpha0, alpha, tau, c_i are fitted to the real y1 through the simulated
## inclusive value and are reported but are not recovery targets). Recovery targets: mu0, delta_lead, delta_join, sigma_g (3 x 2), rho,
## beta, gamma, and mu_state (bilateral ideal points of the states with events).
## Output: v2/sim/results/context_b_recovery.csv (parameter, truth, posterior mean, sd, 5 / 95%, covered, bias, Rhat, n_eff),
##         v2/sim/results/context_b_recovery_design.csv (the realised design counts), v2/sim/fits/context_b_sim.rds, a log in v2/logs.
## Seed 123 (audit F11). The fit uses the compiled m1_context program (compile_models.R) with the default M1 adapt_delta / max_treedepth.
## ============================================================================
suppressPackageStartupMessages({ library(rstan); library(dplyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
sim_dir <- normalizePath(dirname(this_file)); R_DIR <- file.path(sim_dir, "..", "R")
Sys.setenv(V2_ROOT = normalizePath(file.path(sim_dir, "..")))
source(file.path(R_DIR, "spec.R")); source(file.path(R_DIR, "helpers_data.R")); source(file.path(R_DIR, "helpers_model.R"))
source(file.path(R_DIR, "assemble_cell.R")); source(file.path(R_DIR, "typology.R")); source(file.path(R_DIR, "context_model.R"))
args <- commandArgs(trailingOnly = TRUE)
N_INIT <- as.integer(get_arg(args, "n_init", 250)); N_NONINIT <- as.integer(get_arg(args, "n_noninit", 250))
CHAINS <- as.integer(get_arg(args, "chains", 2)); ITER <- as.integer(get_arg(args, "iter", 1000)); WARMUP <- as.integer(get_arg(args, "warmup", 500))
SEED <- as.integer(get_arg(args, "seed", V2_SEED)); SAMPLE <- get_arg(args, "sample", "fmid"); LAG <- as.integer(get_arg(args, "lag", 3)); PILOT <- get_arg(args, "pilot", FALSE)
stopifnot(N_INIT >= 10, N_NONINIT >= 0, CHAINS >= 1, ITER > WARMUP, SAMPLE %in% ALL_SAMPLES, LAG %in% LAGS)
log_open("sim_context_b"); set.seed(SEED)
rstan_options(auto_write = FALSE); options(mc.cores = min(CHAINS, n_cores()))
dir.create(file.path(sim_dir, "fits"), showWarnings = FALSE); res_dir <- file.path(sim_dir, "results"); dir.create(res_dir, showWarnings = FALSE)
logf("M1c simulation recovery on the observed %s lag-%d design: n_init %d, n_noninit %d, chains %d, iter %d, warmup %d, seed %d", SAMPLE, LAG, N_INIT, N_NONINIT, CHAINS, ITER, WARMUP, SEED)

## ---- the observed design -----------------------------------------------------------------------------------------------------------
A <- assemble_cell_spec(SAMPLE, LAG, parse_spec(character(0)), model = "M1c")
sy <- A$sy; alts <- A$alts; ev <- A$ev
stopifnot("grp" %in% names(ev), tail(A$W_cols, 1) == "multi_mid_jan")
S <- length(A$state_codes)
init_sy <- sort(unique(ev$sy_row)); lead_sy <- sort(unique(ev$sy_row[ev$grp == 1L]))
stopifnot(length(lead_sy) <= N_INIT)
keep_init <- sort(c(lead_sy, sample(setdiff(init_sy, lead_sy), N_INIT - length(lead_sy))))
non_init <- setdiff(seq_len(nrow(sy)), init_sy)
keep_non <- sort(sample(non_init, min(N_NONINIT, length(non_init))))
keep <- sort(c(keep_init, keep_non))
logf("design: %d initiating state-years with usable events (%d with a leader event, all kept), %d kept; %d of %d non-initiating state-years kept; N1 = %d",
     length(init_sy), length(lead_sy), length(keep_init), length(keep_non), length(non_init), length(keep))
## subset the frames: sy rows, their alternative blocks, their events; re-index sy_row / alt_row
new_sy <- match(seq_len(nrow(sy)), keep)                     # old sy_row -> new (NA when dropped)
sy_s <- sy[keep, ]
alts_s <- alts |> mutate(.old_row = row_number()) |> filter(sy_row %in% keep) |> mutate(sy_row = new_sy[sy_row])
ev_s <- ev |> filter(sy_row %in% keep) |> mutate(sy_row = new_sy[sy_row], alt_row = match(alt_row, alts_s$.old_row))
stopifnot(!anyNA(ev_s$alt_row), !anyNA(ev_s$sy_row), all(alts_s$sy_row[ev_s$alt_row] == ev_s$sy_row), all(diff(alts_s$sy_row) >= 0))
alts_s <- alts_s |> select(-.old_row)
## the kept state-years' y1 stays observed; state-years kept without events must have y1 = 0 unless they are orphans of the full cell
stopifnot(all(sy_s$y1[ev_s$sy_row] == 1L), all(sy_s$y1[-unique(ev_s$sy_row)] %in% 0:1))
x_ref <- colMeans(as.matrix(sy_s[sy_s$y1 == 1L, A$X_cols, drop = FALSE]))
dat <- build_stan_data(sy_s, alts_s, ev_s, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, x_ref = x_ref, allow_orphan_initiations = TRUE, D_mode = 3L, S = S)
dat$ev_grp <- as.integer(ev_s$grp)
design <- data.frame(group = STAGE2_GROUPS, n_events = as.integer(tabulate(dat$ev_grp + 1L, nbins = 3L)),
                     w_sum = vapply(0:2, function(g) sum(dat$ev_w[dat$ev_grp == g]), numeric(1)),
                     n_states = vapply(0:2, function(g) n_distinct(dat$sy_state[dat$ev_sy[dat$ev_grp == g]]), 1L),
                     n_state_years = vapply(0:2, function(g) n_distinct(dat$ev_sy[dat$ev_grp == g]), 1L))
logf("Stan data: S=%d N1=%d N2=%d E=%d P1=%d P=%d Pz=%d; sum(y1)=%d; groups: %s; choice sets median %d (min %d, max %d)",
     dat$S, dat$N1, dat$N2, dat$E, dat$P1, dat$P, dat$Pz, sum(dat$y1), paste(sprintf("%s %d events (w %.1f, %d states)", design$group, design$n_events, design$w_sum, design$n_states), collapse = "; "),
     as.integer(median(dat$cs_end - dat$cs_start + 1L)), min(dat$cs_end - dat$cs_start + 1L), max(dat$cs_end - dat$cs_start + 1L))
stopifnot(all(design$n_events > 0))

## ---- TRUTH (scaled gap units; one unit = a 26-fold CINC ratio / a 1.7-fold PageRank ratio) ------------------------------------------
set.seed(SEED + 1L)
TRUTH <- list(mu0 = c(-1.0, 0.3), delta_lead = c(-1.2, 0.8), delta_join = c(0.9, -0.5),
              sigma_g = rbind(bilateral = c(1.2, 0.9), leader = c(2.2, 1.4), joiner = c(1.5, 2.0)), rho = 0.4,
              beta = matrix(c(0.3, -0.2, 0.1, 0.4), 2, 2), gamma = c(1.0, -0.3, 0.4, 0.8, 0.1, -0.1, 0.2)[seq_len(dat$Pz)],
              alpha0 = -2.5, alpha = c(0.2, -0.1, 0.5, 0.3, 0.4, 0.3, -0.2, 0.1, 0.3)[seq_len(dat$P1)], tau = 0.5,
              sd_re = c(0.9, 0.7, 0.6), R_Om = matrix(c(1, 0.3, 0.1, 0.3, 1, 0.2, 0.1, 0.2, 1), 3), sigma_b = 0.6)
stopifnot(length(TRUTH$gamma) == dat$Pz, length(TRUTH$alpha) == dat$P1, dat$P == 2L)
L_true <- t(chol(TRUTH$R_Om))
z_re <- matrix(rnorm(3 * S), 3, S)
re <- t(diag(TRUTH$sd_re) %*% L_true %*% z_re)                 # S x 3: (c_i, a_i1, a_i2)
TRUTH$c_i <- re[, 1]; TRUTH$a_i <- re[, 2:3]
z_b <- rnorm(S); TRUTH$b_j <- TRUTH$sigma_b * (z_b - mean(z_b))
TRUTH$mu_state <- rep(1, S) %o% TRUTH$mu0 + rep(1, S) %o% as.numeric(x_ref %*% TRUTH$beta) + TRUTH$a_i

## ---- simulate the chosen targets from the M1c stage-2 model (the slots, groups, weights and y1 are the design) ------------------------
Mu <- rep(1, dat$N1) %o% TRUTH$mu0 + dat$X %*% TRUTH$beta + TRUTH$a_i[dat$sy_state, , drop = FALSE]
zg <- as.numeric(dat$Z %*% TRUTH$gamma)
shift <- rbind(c(0, 0), TRUTH$delta_lead, TRUTH$delta_join)
Sinv <- lapply(1:3, function(g) solve(matrix(c(TRUTH$sigma_g[g, 1]^2, TRUTH$rho * prod(TRUTH$sigma_g[g, ]), TRUTH$rho * prod(TRUTH$sigma_g[g, ]), TRUTH$sigma_g[g, 2]^2), 2)))
util_rows <- function(rows, n, g) {
  d <- sweep(dat$x[rows, , drop = FALSE], 2, Mu[n, ] + shift[g + 1L, ])
  -0.5 * rowSums((d %*% Sinv[[g + 1L]]) * d) + TRUTH$b_j[dat$alt_target[rows]] + zg[rows]
}
set.seed(SEED + 2L)
ev_sim <- ev_s
taken <- vector("list", dat$N1)
for (e in seq_len(dat$E)) {
  n <- dat$ev_sy[e]; g <- dat$ev_grp[e]
  rows <- dat$cs_start[n]:dat$cs_end[n]; rows <- rows[dat$in_s2[rows] == 1L]
  rows <- setdiff(rows, taken[[n]])                             # without replacement within the state-year
  stopifnot(length(rows) >= 1)
  u <- util_rows(rows, n, g); p <- exp(u - max(u)); p <- p / sum(p)
  pick <- rows[sample.int(length(rows), 1L, prob = p)]
  taken[[n]] <- c(taken[[n]], pick)
  ev_sim$alt_row[e] <- pick
}
stopifnot(all(dat$in_s2[ev_sim$alt_row] == 1L), all(alts_s$sy_row[ev_sim$alt_row] == ev_sim$sy_row), !anyDuplicated(paste(ev_sim$sy_row, ev_sim$alt_row)))
dat_sim <- build_stan_data(sy_s, alts_s, ev_sim, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, x_ref = x_ref, allow_orphan_initiations = TRUE, D_mode = 3L, S = S)
dat_sim$ev_grp <- as.integer(ev_sim$grp)
stopifnot(identical(dat_sim$ev_sy, dat$ev_sy), identical(dat_sim$ev_grp, dat$ev_grp), identical(dat_sim$ev_w, dat$ev_w), identical(dat_sim$y1, dat$y1))
## the simulated targets must favour the group's ideal point: the mean chosen gap per group against the truth
chosen <- dat_sim$x[dat_sim$ev_row, , drop = FALSE]
for (g in 0:2) logf("simulated chosen gaps, group %s: mean (%.2f, %.2f) over %d events; truth population ideal point (%.2f, %.2f)", STAGE2_GROUPS[g + 1L],
                    mean(chosen[dat_sim$ev_grp == g, 1]), mean(chosen[dat_sim$ev_grp == g, 2]), sum(dat_sim$ev_grp == g), TRUTH$mu0[1] + shift[g + 1L, 1], TRUTH$mu0[2] + shift[g + 1L, 2])
## the brute-force likelihood at the truth is finite (sanity)
R_true <- stage2_loglik_R_context(dat_sim, list(mu0 = TRUTH$mu0, delta_lead = TRUTH$delta_lead, delta_join = TRUTH$delta_join, beta = TRUTH$beta, a_i = TRUTH$a_i, b_j = TRUTH$b_j, gamma = TRUTH$gamma, sigma_g = unname(TRUTH$sigma_g), rho = TRUTH$rho))
stopifnot(is.finite(R_true$lp_stage2))
logf("stage-2 weighted log-likelihood at the truth: %.1f (mean per event %.3f; uniform choice would give %.3f)", R_true$lp_stage2, mean(R_true$loglik), mean(-log(dat_sim$cs_end - dat_sim$cs_start + 1)))

## ---- fit ---------------------------------------------------------------------------------------------------------------------------
sm <- readRDS(file.path(v2_root(), "stan", "build", "m1_context.rds"))
assert_compiled_current(sm, file.path(v2_root(), "stan", "m1_context.stan"), file.path(v2_root(), "stan", "build", "m1_context.sha256"))
if (PILOT) { CHAINS <- 1L; ITER <- 40L; WARMUP <- 20L; logf("PILOT: 1 chain x 40 iterations (timing only)") }
dat_sim$grainsize <- dat_sim$N1
rstan_options(threads_per_chain = 1)
t0 <- Sys.time()
fit <- sampling(sm, data = dat_sim, chains = CHAINS, iter = ITER, warmup = WARMUP, seed = SEED, refresh = max(1L, ITER %/% 10L),
                control = list(adapt_delta = SAMPLER$M1$adapt_delta, max_treedepth = SAMPLER$M1$max_treedepth),
                pars = c("z_re", "z_b", "L_Om", "re", "log_lik_stage1", "log_lik_stage2", "mu_event"), include = FALSE)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
if (fit@mode != 0L) stop("sampling did not run")
sp <- get_sampler_params(fit, inc_warmup = FALSE)
logf("sampling finished in %.1f min; divergences %d; max-treedepth hits %d", elapsed, sum(sapply(sp, function(m) sum(m[, "divergent__"]))), sum(sapply(sp, function(m) sum(m[, "treedepth__"] >= SAMPLER$M1$max_treedepth))))
if (PILOT) { logf("PILOT done; no recovery table"); quit(save = "no", status = 0) }

## ---- recovery table ---------------------------------------------------------------------------------------------------------------
pars <- c("mu0", "delta_lead", "delta_join", "mu0_lead", "mu0_join", "sigma_g", "rho_x", "beta", "gamma", "tau", "alpha0", "alpha", "sd_re", "sigma_b")
sm_tab <- summary(fit, pars = pars, probs = c(0.05, 0.95))$summary
truth_of <- function(nm) {
  base <- sub("\\[.*$", "", nm); i1 <- suppressWarnings(as.integer(sub("^[^\\[]+\\[([0-9]+).*$", "\\1", nm))); i2 <- suppressWarnings(as.integer(sub("^[^\\[]+\\[[0-9]+,([0-9]+)\\]$", "\\1", nm)))
  switch(base, mu0 = TRUTH$mu0[i1], delta_lead = TRUTH$delta_lead[i1], delta_join = TRUTH$delta_join[i1], mu0_lead = TRUTH$mu0[i1] + TRUTH$delta_lead[i1],
         mu0_join = TRUTH$mu0[i1] + TRUTH$delta_join[i1], sigma_g = TRUTH$sigma_g[i1, i2], rho_x = TRUTH$rho, beta = TRUTH$beta[i1, i2], gamma = TRUTH$gamma[i1],
         tau = TRUTH$tau, alpha0 = TRUTH$alpha0, alpha = TRUTH$alpha[i1], sd_re = TRUTH$sd_re[i1], sigma_b = TRUTH$sigma_b, NA_real_)
}
rec <- data.frame(parameter = rownames(sm_tab), truth = vapply(rownames(sm_tab), truth_of, numeric(1)), mean = sm_tab[, "mean"], sd = sm_tab[, "sd"],
                  lo = sm_tab[, "5%"], hi = sm_tab[, "95%"], Rhat = sm_tab[, "Rhat"], n_eff = sm_tab[, "n_eff"], stringsAsFactors = FALSE)
rec$covered <- rec$truth >= rec$lo & rec$truth <= rec$hi
rec$bias <- rec$mean - rec$truth
rec$z <- rec$bias / rec$sd
rec$recovery_target <- !sub("\\[.*$", "", rec$parameter) %in% c("tau", "alpha0", "alpha", "sd_re", "sigma_b")   # stage 1 and the RE scales are fitted to the REAL y1 / a simulated a_i, not targets
rownames(rec) <- NULL
## mu_state (bilateral ideal points) of the states with events: correlation of the posterior mean with the truth, coverage per state
ms <- rstan::extract(fit, pars = "mu_state")$mu_state
with_ev <- sort(unique(dat_sim$sy_state[dat_sim$ev_sy]))
ms_mean <- apply(ms, c(2, 3), mean); ms_lo <- apply(ms, c(2, 3), quantile, 0.05); ms_hi <- apply(ms, c(2, 3), quantile, 0.95)
cov_ms <- mean(TRUTH$mu_state[with_ev, ] >= ms_lo[with_ev, ] & TRUTH$mu_state[with_ev, ] <= ms_hi[with_ev, ])
cor_ms <- c(cor(ms_mean[with_ev, 1], TRUTH$mu_state[with_ev, 1]), cor(ms_mean[with_ev, 2], TRUTH$mu_state[with_ev, 2]))
rec <- bind_rows(rec, data.frame(parameter = c("mu_state[material] cor(post mean, truth) over states with events", "mu_state[relational] cor(post mean, truth) over states with events", "mu_state 90% coverage over states with events x axes"),
                                 truth = c(1, 1, 0.9), mean = c(cor_ms, cov_ms), sd = NA, lo = NA, hi = NA, Rhat = NA, n_eff = NA, covered = NA, bias = NA, z = NA, recovery_target = TRUE))
rec$n_init <- N_INIT; rec$n_noninit <- N_NONINIT; rec$chains <- CHAINS; rec$iter <- ITER; rec$warmup <- WARMUP; rec$seed <- SEED; rec$elapsed_min <- elapsed
write.csv(rec, file.path(res_dir, "context_b_recovery.csv"), row.names = FALSE)
write.csv(design |> mutate(N1 = dat_sim$N1, N2 = dat_sim$N2, E = dat_sim$E, S = dat_sim$S, sum_y1 = sum(dat_sim$y1), n_states_with_events = length(with_ev)),
          file.path(res_dir, "context_b_recovery_design.csv"), row.names = FALSE)
saveRDS(list(fit = fit, dat = dat_sim, truth = TRUTH, design = design, keep = keep, args = list(n_init = N_INIT, n_noninit = N_NONINIT, chains = CHAINS, iter = ITER, warmup = WARMUP, seed = SEED, sample = SAMPLE, lag = LAG)),
        file.path(sim_dir, "fits", "context_b_sim.rds"))
tg <- rec |> filter(recovery_target, !is.na(covered))
logf("recovery over the %d stage-b / stage-2 parameters: 90%% coverage %.2f, mean |bias| %.3f, max |z| %.2f, max Rhat %.3f, min n_eff %.0f", nrow(tg), mean(tg$covered), mean(abs(tg$bias)), max(abs(tg$z)), max(tg$Rhat), min(tg$n_eff))
key <- rec |> filter(grepl("^(mu0|delta_lead|delta_join|sigma_g|rho_x)", parameter))
for (r in seq_len(nrow(key))) logf("  %-14s truth %6.2f  mean %6.2f  [%6.2f, %6.2f]  covered %-5s  bias %6.2f  Rhat %.3f", key$parameter[r], key$truth[r], key$mean[r], key$lo[r], key$hi[r], key$covered[r], key$bias[r], key$Rhat[r])
logf("mu_state recovery (states with events, n = %d): cor material %.3f, relational %.3f; 90%% coverage %.2f", length(with_ev), cor_ms[1], cor_ms[2], cov_ms)
logf("written: %s, %s", file.path(res_dir, "context_b_recovery.csv"), file.path(res_dir, "context_b_recovery_design.csv"))
