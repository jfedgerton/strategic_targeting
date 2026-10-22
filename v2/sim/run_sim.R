## ============================================================================
## run_sim.R -- fit M1a, M1b, M2a, M2b (K = 2) to data simulated from M1b and M2b
##
## Usage: Rscript v2/sim/run_sim.R [chains] [iter] [warmup]
## Defaults are sized for a 2-core / 8 GB container: chains = 2, iter = 700, warmup = 350.
## On Roar use chains = 4, iter = 2000, warmup = 1000.
## Writes: v2/sim/fits/<truth>_<model>.rds (gitignored), v2/sim/results/*.csv
## Seed 123 everywhere (audit F11).
## ============================================================================
suppressPackageStartupMessages({ library(rstan); library(loo); library(posterior); library(dplyr) })
args <- commandArgs(trailingOnly = TRUE)
CHAINS <- if (length(args) >= 1) as.integer(args[1]) else 2L
ITER   <- if (length(args) >= 2) as.integer(args[2]) else 700L
WARMUP <- if (length(args) >= 3) as.integer(args[3]) else 350L
ONLY   <- if (length(args) >= 4) strsplit(args[4], ",")[[1]] else NULL   # e.g. "M1b:M1a,M2b:M2b"

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
sim_dir <- normalizePath(dirname(this_file))
Sys.setenv(V2_ROOT = normalizePath(file.path(sim_dir, "..")))
source(file.path(sim_dir, "..", "R", "00_helpers.R"))
source(file.path(sim_dir, "..", "R", "typology.R"))
source(file.path(sim_dir, "simulate.R"))
log_open("run_sim")
set.seed(V2_SEED)
rstan_options(auto_write = FALSE)
options(mc.cores = min(CHAINS, parallel::detectCores()))
dir.create(file.path(sim_dir, "fits"), showWarnings = FALSE)
dir.create(file.path(sim_dir, "results"), showWarnings = FALSE)

models <- list(m1 = readRDS(file.path(v2_root(), "stan", "build", "m1_latent_additive.rds")),
               m2 = readRDS(file.path(v2_root(), "stan", "build", "m2_ame.rds")))
W_cols <- c("W1", "W2"); X_cols <- paste0("X", 1:4); Z_cols <- c("Z1", "Z2")
specs <- data.frame(model = c("M1a", "M1b", "M2a", "M2b"),
                    family = c("m1", "m1", "m2", "m2"),
                    joint = c(0L, 1L, 0L, 1L), stringsAsFactors = FALSE)

recovery_rows <- list(); conv_rows <- list(); class_rows <- list(); loo_rows <- list(); ame_rows <- list()
fits_by_truth <- list()

for (truth_family in c("M1b", "M2b")) {
  sim <- simulate_targeting(truth_family, seed = V2_SEED)
  logf("truth = %s : S=%d N1=%d N2=%d E=%d", truth_family, sim$meta$S, sim$meta$N1, sim$meta$N2, sim$meta$E)
  write.csv(data.frame(state_idx = seq_len(sim$meta$S), type = sim$types,
                       a1 = sim$truth$a_i[, 1], a2 = sim$truth$a_i[, 2], c = sim$truth$c_i),
            file.path(sim_dir, "results", paste0("truth_", truth_family, "_states.csv")), row.names = FALSE)
  for (k in seq_len(nrow(specs))) {
    sp <- specs[k, ]
    tag <- paste0(truth_family, ":", sp$model)
    if (!is.null(ONLY) && !(tag %in% ONLY)) next
    dat <- build_stan_data(sim$sy, sim$alts, sim$ev, W_cols, X_cols, Z_cols, joint = sp$joint,
                           x_ref = rep(0, 4), K = if (sp$family == "m2") 2L else NULL)
    logf("fitting %s on %s data (chains=%d iter=%d warmup=%d)", sp$model, truth_family, CHAINS, ITER, WARMUP)
    t0 <- Sys.time()
    fit <- sampling(models[[sp$family]], data = dat, chains = CHAINS, iter = ITER, warmup = WARMUP,
                    seed = V2_SEED, refresh = max(1L, ITER %/% 10L),
                    control = list(adapt_delta = 0.9, max_treedepth = 10),
                    pars = c("z_re", "z_b", "V_raw", "U_free", "L_Om_joint", "L_Om", "re"), include = FALSE)
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    logf("  done in %.1f min", el)
    saveRDS(list(fit = fit, dat = dat, spec = sp, truth_family = truth_family, elapsed_min = el),
            file.path(sim_dir, "fits", paste0(truth_family, "_", sp$model, ".rds")))

    ## ---- convergence ------------------------------------------------------------------
    sdiag <- get_sampler_params(fit, inc_warmup = FALSE)
    n_div <- sum(sapply(sdiag, function(m) sum(m[, "divergent__"])))
    n_maxtd <- sum(sapply(sdiag, function(m) sum(m[, "treedepth__"] >= 10)))
    sm <- summary(fit, pars = c("alpha0", "alpha", "tau", "beta", "sigma_x", "sd_re", "sigma_b", "gamma",
                                "Sigma", "Omega", "R_Om", "mu_state", if (sp$family == "m2") c("U", "V")))$summary
    sm <- sm[is.finite(sm[, "Rhat"]), , drop = FALSE]     # structural zeros in U anchor block -> NA
    conv_rows[[tag]] <- data.frame(truth = truth_family, model = sp$model, elapsed_min = round(el, 1),
                                   divergences = n_div, max_treedepth_hits = n_maxtd,
                                   worst_rhat = max(sm[, "Rhat"]), min_ess_bulk = min(sm[, "n_eff"]),
                                   frac_rhat_gt_1.05 = mean(sm[, "Rhat"] > 1.05))
    ## ---- recovery and 90% coverage ---------------------------------------------------------
    tr <- sim$truth
    post <- rstan::extract(fit, pars = c("alpha0", "alpha", "tau", "beta", "Sigma", "Omega", "gamma", "mu_state"))
    q <- function(x) c(mean = mean(x), lo = unname(quantile(x, 0.05)), hi = unname(quantile(x, 0.95)))
    rows <- list()
    add <- function(name, draws, truth) rows[[length(rows) + 1]] <<- data.frame(param = name, truth = truth, t(q(draws)))
    add("alpha0", post$alpha0, tr$alpha0)
    for (p in 1:2) add(sprintf("alpha[%d]", p), post$alpha[, p], tr$alpha[p])
    add("tau", post$tau, tr$tau)
    for (p in 1:4) for (d in 1:2) add(sprintf("beta[%d,%d]", p, d), post$beta[, p, d], tr$beta[p, d])
    for (p in 1:2) add(sprintf("gamma[%d]", p), post$gamma[, p], tr$gamma[p])
    add("Sigma[1,1]", post$Sigma[, 1, 1], tr$Sigma[1, 1]); add("Sigma[2,2]", post$Sigma[, 2, 2], tr$Sigma[2, 2])
    add("Sigma[1,2]", post$Sigma[, 1, 2], tr$Sigma[1, 2])
    for (i in 1:3) for (j in i:3) add(sprintf("Omega[%d,%d]", i, j), post$Omega[, i, j], tr$Omega[i, j])
    rec <- bind_rows(rows) |> mutate(truth_family = truth_family, model = sp$model,
                                      covered = truth >= lo & truth <= hi, abs_err = abs(mean - truth),
                                      fixed_zero = (model %in% c("M1a", "M2a")) & param %in% c("Sigma[1,2]", "Omega[2,3]"))
    recovery_rows[[tag]] <- rec
    ## ---- type classification from posterior mu_state (x_ref = 0 so mu_i = a_i) --------------
    cls <- classify_types(post$mu_state, r0 = 0.8)           # draws x S x 2 -> per-state table
    cls$true_type <- sim$types
    cls$state_idx <- seq_len(sim$meta$S)
    acc <- mean(cls$map_type == cls$true_type)
    class_rows[[tag]] <- data.frame(truth = truth_family, model = sp$model, accuracy = acc,
                                    n_states = sim$meta$S,
                                    mean_abs_err_mu = mean(abs(apply(post$mu_state, c(2, 3), mean) - tr$a_i)))
    write.csv(cls, file.path(sim_dir, "results", paste0("classification_", truth_family, "_", sp$model, ".csv")), row.names = FALSE)
    ## ---- LOO, stage 1 and stage 2 kept separate ----------------------------------------------
    ll2 <- loo::extract_log_lik(fit, "log_lik_stage2", merge_chains = FALSE)
    ll1 <- loo::extract_log_lik(fit, "log_lik_stage1", merge_chains = FALSE)
    loo2 <- loo::loo(ll2, r_eff = loo::relative_eff(exp(ll2)))
    loo1 <- loo::loo(ll1, r_eff = loo::relative_eff(exp(ll1)))
    loo_rows[[tag]] <- data.frame(truth = truth_family, model = sp$model,
                                  elpd_stage2 = loo2$estimates["elpd_loo", 1], se_stage2 = loo2$estimates["elpd_loo", 2],
                                  bad_k_stage2 = sum(loo2$diagnostics$pareto_k > 0.7),
                                  elpd_stage1 = loo1$estimates["elpd_loo", 1], se_stage1 = loo1$estimates["elpd_loo", 2],
                                  bad_k_stage1 = sum(loo1$diagnostics$pareto_k > 0.7))
    fits_by_truth[[truth_family]][[sp$model]] <- list(loo1 = loo1, loo2 = loo2)
    ## ---- AME recovery ---------------------------------------------------------------------------
    if (sp$family == "m2") {
      Ud <- rstan::extract(fit, pars = "U", permuted = FALSE)   # iter x chain x (S*K)
      Vd <- rstan::extract(fit, pars = "V", permuted = FALSE)
      S <- dat$S; K <- dat$K
      U4 <- array(Ud, c(dim(Ud)[1], dim(Ud)[2], S, K)); V4 <- array(Vd, c(dim(Vd)[1], dim(Vd)[2], S, K))
      al <- procrustes_align(U4, V4)
      rot_dev <- max(sapply(al$rotations, function(Rm) max(abs(Rm - diag(K)))))
      UV_draws <- array(NA_real_, c(dim(U4)[1] * dim(U4)[2], S, S)); ii <- 0
      for (ch in seq_len(dim(U4)[2])) for (it in seq_len(dim(U4)[1])) {
        ii <- ii + 1; UV_draws[ii, , ] <- matrix(al$U[it, ch, , ], S, K) %*% t(matrix(al$V[it, ch, , ], S, K))
      }
      UV_mean <- apply(UV_draws, c(2, 3), mean)
      off <- row(UV_mean) != col(UV_mean)
      if (!is.null(tr$UV)) {
        cor_uv <- cor(UV_mean[off], tr$UV[off])
        lo <- apply(UV_draws, c(2, 3), quantile, 0.05); hi <- apply(UV_draws, c(2, 3), quantile, 0.95)
        cov_uv <- mean((tr$UV >= lo & tr$UV <= hi)[off])
        ## alternative: u'v recovered only up to the b_j / a_i decomposition; also report sender-centred version
      } else { cor_uv <- NA; cov_uv <- NA }
      ame_rows[[tag]] <- data.frame(truth = truth_family, model = sp$model, cor_UV_offdiag = cor_uv,
                                    coverage90_UV_offdiag = cov_uv, sd_UV_posterior_mean = sd(UV_mean[off]),
                                    max_procrustes_rotation_dev = rot_dev)
    }
    rm(fit, post); gc()
  }
}

## ---- pairwise LOO comparisons a vs b within truth ---------------------------------------------
cmp_rows <- list()
for (tf in names(fits_by_truth)) {
  f <- fits_by_truth[[tf]]
  for (pair in list(c("M1a", "M1b"), c("M2a", "M2b"), c("M1b", "M2b"))) {
    if (!all(pair %in% names(f))) next
    c2 <- loo::loo_compare(f[[pair[1]]]$loo2, f[[pair[2]]]$loo2)
    c1 <- loo::loo_compare(f[[pair[1]]]$loo1, f[[pair[2]]]$loo1)
    cmp_rows[[length(cmp_rows) + 1]] <- data.frame(
      truth = tf, comparison = paste(pair, collapse = " vs "),
      stage2_elpd_diff_b_minus_a = f[[pair[2]]]$loo2$estimates["elpd_loo", 1] - f[[pair[1]]]$loo2$estimates["elpd_loo", 1],
      stage2_se_diff = c2[2, "se_diff"],
      stage1_elpd_diff_b_minus_a = f[[pair[2]]]$loo1$estimates["elpd_loo", 1] - f[[pair[1]]]$loo1$estimates["elpd_loo", 1],
      stage1_se_diff = c1[2, "se_diff"])
  }
}
## misclassification between a and b
mis_rows <- list()
for (tf in c("M1b", "M2b")) for (fam in c("M1", "M2")) {
  fa <- file.path(sim_dir, "results", sprintf("classification_%s_%sa.csv", tf, fam))
  fb <- file.path(sim_dir, "results", sprintf("classification_%s_%sb.csv", tf, fam))
  if (file.exists(fa) && file.exists(fb)) {
    a <- read.csv(fa); b <- read.csv(fb)
    mis_rows[[length(mis_rows) + 1]] <- data.frame(truth = tf, family = fam,
      disagreement_a_vs_b = mean(a$map_type != b$map_type),
      acc_a = mean(a$map_type == a$true_type), acc_b = mean(b$map_type == b$true_type))
  }
}
res_dir <- file.path(sim_dir, "results")
write.csv(bind_rows(recovery_rows), file.path(res_dir, "recovery.csv"), row.names = FALSE)
write.csv(bind_rows(conv_rows), file.path(res_dir, "convergence.csv"), row.names = FALSE)
write.csv(bind_rows(class_rows), file.path(res_dir, "classification_accuracy.csv"), row.names = FALSE)
write.csv(bind_rows(loo_rows), file.path(res_dir, "loo.csv"), row.names = FALSE)
if (length(cmp_rows)) write.csv(bind_rows(cmp_rows), file.path(res_dir, "loo_compare.csv"), row.names = FALSE)
if (length(mis_rows)) write.csv(bind_rows(mis_rows), file.path(res_dir, "misclassification_a_vs_b.csv"), row.names = FALSE)
if (length(ame_rows)) write.csv(bind_rows(ame_rows), file.path(res_dir, "ame_recovery.csv"), row.names = FALSE)
logf("run_sim finished; results in %s", res_dir)
