## ============================================================================
## run_sim_1d.R -- small recovery simulation for the one-dimensional model M1m (D_mode = 1)
##
## Usage: Rscript v2/sim/run_sim_1d.R [chains] [iter] [warmup]
## Truth: simulate_targeting("M1m") (stage-2 utility on dM only). Fits M1m and, for contrast, M1r
## (the wrong axis) and M1b (both axes) on the same data; reports recovery / coverage of beta[,1],
## sigma_x[1], Omega (c, a1), tau, gamma, the 3-level axis classification, and LOO.
## Writes v2/sim/results/sim1d_*.csv. Seed 123.
## ============================================================================
suppressPackageStartupMessages({ library(rstan); library(loo); library(dplyr) })
args <- commandArgs(trailingOnly = TRUE)
CHAINS <- if (length(args) >= 1) as.integer(args[1]) else 2L
ITER   <- if (length(args) >= 2) as.integer(args[2]) else 500L
WARMUP <- if (length(args) >= 3) as.integer(args[3]) else 250L
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
sim_dir <- normalizePath(dirname(this_file)); Sys.setenv(V2_ROOT = normalizePath(file.path(sim_dir, "..")))
source(file.path(sim_dir, "..", "R", "00_helpers.R")); source(file.path(sim_dir, "..", "R", "typology.R")); source(file.path(sim_dir, "simulate.R"))
log_open("run_sim_1d"); set.seed(V2_SEED); rstan_options(auto_write = FALSE); options(mc.cores = min(CHAINS, parallel::detectCores()))
res_dir <- file.path(sim_dir, "results"); dir.create(res_dir, showWarnings = FALSE)
sm <- readRDS(file.path(v2_root(), "stan", "build", "m1_latent_additive.rds"))
sim <- simulate_targeting("M1m", S = 30L, Tt = 12L, seed = V2_SEED)
logf("1-D world: S=%d N1=%d N2=%d E=%d", sim$meta$S, sim$meta$N1, sim$meta$N2, sim$meta$E)
tr <- sim$truth
W_cols <- c("W1", "W2"); X_cols <- paste0("X", 1:4); Z_cols <- c("Z1", "Z2")
excl <- c("z_re", "z_b", "V_raw", "U_free", "L_Om_joint", "L_Om", "re")
specs <- list(M1m = list(D_mode = 1L, joint = 0L), M1r = list(D_mode = 2L, joint = 0L), M1b = list(D_mode = 3L, joint = 1L))
rec <- list(); conv <- list(); cls_rows <- list(); loos <- list()
q <- function(x) c(mean = mean(x), lo = unname(quantile(x, 0.05)), hi = unname(quantile(x, 0.95)))
for (nm in names(specs)) {
  sp <- specs[[nm]]
  dat <- build_stan_data(sim$sy, sim$alts, sim$ev, W_cols, X_cols, Z_cols, joint = sp$joint, x_ref = rep(0, 4), D_mode = sp$D_mode)
  t0 <- Sys.time()
  fit <- sampling(sm, data = dat, chains = CHAINS, iter = ITER, warmup = WARMUP, seed = V2_SEED, refresh = 0,
                  control = list(adapt_delta = 0.9), pars = excl, include = FALSE)
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins")); logf("%s fitted in %.1f min", nm, el)
  sdiag <- get_sampler_params(fit, inc_warmup = FALSE)
  smry <- summary(fit, pars = c("alpha0", "alpha", "tau", "mu0", "beta", "sigma_x", "sd_re", "sigma_b", "gamma", "Omega"))$summary
  conv[[nm]] <- data.frame(model = nm, elapsed_min = round(el, 1), divergences = sum(sapply(sdiag, function(m) sum(m[, "divergent__"]))),
                           worst_rhat = max(smry[, "Rhat"], na.rm = TRUE), min_ess = min(smry[, "n_eff"], na.rm = TRUE))
  post <- rstan::extract(fit, pars = c("tau", "mu0", "beta", "Sigma", "Omega", "gamma", "mu_state"))
  rows <- list(); add <- function(name, draws, truth) rows[[length(rows) + 1]] <<- data.frame(param = name, truth = truth, t(q(draws)))
  add("tau", post$tau, tr$tau)
  D <- if (sp$D_mode == 3L) 2L else 1L
  if (D == 1L) {
    add("mu0[1]", as.numeric(post$mu0), if (sp$D_mode == 1L) tr$mu0[1] else NA)   # [draws, 1] for row_vector[1]
    b1 <- matrix(post$beta, nrow = dim(post$beta)[1])          # rstan keeps a [draws, P, 1] array for matrix[P, 1]
    for (p in 1:4) add(sprintf("beta[%d,1]", p), b1[, p], if (sp$D_mode == 1L) tr$beta[p, 1] else NA)
    add("Sigma[1,1]", as.numeric(post$Sigma), if (sp$D_mode == 1L) tr$sigma_x[1]^2 else NA)
    add("Omega[1,1] (c)", post$Omega[, 1, 1], tr$Omega[1, 1]); add("Omega[2,2] (a)", post$Omega[, 2, 2], if (sp$D_mode == 1L) tr$Omega[2, 2] else NA)
    add("Omega[1,2] (c,a)", post$Omega[, 1, 2], if (sp$D_mode == 1L) tr$Omega[1, 2] else NA)
  } else {
    add("mu0[1]", post$mu0[, 1], tr$mu0[1])
    for (p in 1:4) add(sprintf("beta[%d,1]", p), post$beta[, p, 1], tr$beta[p, 1])
    add("Sigma[1,1]", post$Sigma[, 1, 1], tr$sigma_x[1]^2); add("Sigma[1,2]", post$Sigma[, 1, 2], 0)
    add("Omega[1,1] (c)", post$Omega[, 1, 1], tr$Omega[1, 1]); add("Omega[2,2] (a1)", post$Omega[, 2, 2], tr$Omega[2, 2]); add("Omega[1,2] (c,a1)", post$Omega[, 1, 2], tr$Omega[1, 2])
  }
  for (p in 1:2) add(sprintf("gamma[%d]", p), post$gamma[, p], tr$gamma[p])
  rec[[nm]] <- bind_rows(rows) |> mutate(model = nm, covered = truth >= lo & truth <= hi, abs_err = abs(mean - truth))
  mu <- mu_state_draws(fit)
  if (D == 1L) {
    cls <- classify_types(mu, r0 = 0.8)
    acc <- if (sp$D_mode == 1L) mean(cls$map_type == sim$types_1d) else NA
    cls_rows[[nm]] <- data.frame(model = nm, axis = c("material", "relational")[sp$D_mode], accuracy_3level = acc,
                                 mean_abs_err_mu = if (sp$D_mode == 1L) mean(abs(cls$mu_mean - tr$mu_state_ref0[, 1])) else NA)
  } else {
    cls_rows[[nm]] <- data.frame(model = nm, axis = "both", accuracy_3level = mean(ifelse(cls_m <- apply(mu[, , 1], 2, mean) > 0.8, "weaker", ifelse(cls_m < -0.8, "stronger", "peer")) == sim$types_1d),
                                 mean_abs_err_mu = mean(abs(apply(mu[, , 1], 2, mean) - tr$mu_state_ref0[, 1])))
  }
  ll2 <- loo::extract_log_lik(fit, "log_lik_stage2", merge_chains = FALSE); l2 <- loo::loo(ll2, r_eff = loo::relative_eff(exp(ll2)))
  ll1 <- loo::extract_log_lik(fit, "log_lik_stage1", merge_chains = FALSE); l1 <- loo::loo(ll1, r_eff = loo::relative_eff(exp(ll1)))
  loos[[nm]] <- data.frame(model = nm, elpd_stage2 = l2$estimates["elpd_loo", 1], se_stage2 = l2$estimates["elpd_loo", 2],
                           elpd_stage1 = l1$estimates["elpd_loo", 1], se_stage1 = l1$estimates["elpd_loo", 2])
  rm(fit); gc()
}
write.csv(bind_rows(rec), file.path(res_dir, "sim1d_recovery.csv"), row.names = FALSE)
write.csv(bind_rows(conv), file.path(res_dir, "sim1d_convergence.csv"), row.names = FALSE)
write.csv(bind_rows(cls_rows), file.path(res_dir, "sim1d_classification.csv"), row.names = FALSE)
write.csv(bind_rows(loos), file.path(res_dir, "sim1d_loo.csv"), row.names = FALSE)
print(bind_rows(rec) |> filter(model == "M1m")); print(bind_rows(cls_rows)); print(bind_rows(loos))
logf("run_sim_1d finished")
