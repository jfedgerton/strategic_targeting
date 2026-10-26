## ============================================================================
## run_sim.R -- THE simulation driver: every simulation task in one file, selected by --task
##
## Usage: Rscript v2/sim/run_sim.R --task=main|persistent|1d|measurement|benchmarks [--chains=2] [--iter=700] [--warmup=350] [--only=]
##   main        fit M1a, M1b, M2a, M2b (K = 2) to data simulated from M1b and M2b under the CURRENT model (simulate_current(): non-zero
##               mu0, gap scaling, in_iv / in_s2 masks, minimal W with the power levels and log system size, X = power changes, event
##               weights, r0 = 0.5, x_ref = initiator mean); writes v2/sim/fits/<truth>_<model>.rds (gitignored) and v2/sim/results/*.csv
##               --only = comma-separated "truth:model" cells to (re)fit, e.g. "M1b:M1a,M2b:M2b"; cells not listed reuse their saved
##               fit when present and are skipped otherwise
##   persistent  the same with persistent positions (random-walk step SD 0.06 instead of 0.3); files get the suffix _persistent
##   1d          recovery of the one-dimensional model M1m (truth simulate_current("M1m"); fits M1m, M1r, M1b); sim1d_*.csv
##   measurement the measurement tests of measurement_tests.R on a 30-state world (LOO a vs b, misclassification, temporal holdout
##               with frozen constants, permutation in rescore mode, relational-level error inflation); measurement_*.csv
##   benchmarks  the archived benchmarks (archive/R/benchmarks.R: directed-dyad logit, old gap regression) on the simulated frames,
##               in its own task because the brms compile needs several GB; benchmarks_sim.csv
## Defaults are sized for a 2-core / 8 GB container (chains = 2, iter = 700, warmup = 350); on Roar slurm/05_simulations.sbatch
## passes 4 / 2000 / 1000 (1500 / 750 for measurement). write_sim_report.R assembles SIM_REPORT.md from the result CSVs.
## Seed 123 everywhere (audit F11). The old simulate_targeting() world remains in simulate.R for the unit tests only.
## ============================================================================
suppressPackageStartupMessages({ library(rstan); library(loo); library(posterior); library(dplyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
sim_dir <- normalizePath(dirname(this_file)); R_DIR <- file.path(sim_dir, "..", "R")
Sys.setenv(V2_ROOT = normalizePath(file.path(sim_dir, "..")))
source(file.path(R_DIR, "spec.R")); source(file.path(R_DIR, "helpers_data.R")); source(file.path(R_DIR, "helpers_model.R"))
source(file.path(R_DIR, "typology.R")); source(file.path(R_DIR, "measurement_tests.R")); source(file.path(sim_dir, "simulate.R"))
args <- commandArgs(trailingOnly = TRUE)
TASK <- get_arg(args, "task", "")
if (!TASK %in% c("main", "persistent", "1d", "measurement", "benchmarks")) stop("--task must be main | persistent | 1d | measurement | benchmarks; got '", TASK, "'")
CHAINS <- as.integer(get_arg(args, "chains", 2)); ITER <- as.integer(get_arg(args, "iter", if (TASK %in% c("main", "persistent")) 700 else 500))
WARMUP <- as.integer(get_arg(args, "warmup", if (TASK %in% c("main", "persistent")) 350 else 250))
ONLY <- { o <- get_arg(args, "only", ""); if (nzchar(o)) strsplit(o, ",")[[1]] else NULL }
stopifnot(CHAINS >= 1, ITER > WARMUP, WARMUP >= 0)
log_open(paste0("run_sim_", TASK)); set.seed(V2_SEED)
rstan_options(auto_write = FALSE); options(mc.cores = min(CHAINS, n_cores()))      # n_cores(): detectCores() with a fallback when it returns NA (helpers_model.R)
dir.create(file.path(sim_dir, "fits"), showWarnings = FALSE); res_dir <- file.path(sim_dir, "results"); dir.create(res_dir, showWarnings = FALSE)
logf("task %s: chains %d, iter %d, warmup %d", TASK, CHAINS, ITER, WARMUP)
R0 <- 0.5
excl <- c("z_re", "z_b", "V_raw", "U_free", "L_Om_joint", "L_Om", "re")

## ============================================================================
## main | persistent: M1a, M1b, M2a, M2b on M1b and M2b truths
## ============================================================================
if (TASK %in% c("main", "persistent")) {
  WORLD <- if (TASK == "persistent") "persistent" else ""
  STEP_SD <- if (WORLD == "persistent") 0.06 else 0.3
  SUF <- if (WORLD == "persistent") "_persistent" else ""
  models <- list(m1 = readRDS(file.path(v2_root(), "stan", "build", "m1_latent_additive.rds")),
                 m2 = readRDS(file.path(v2_root(), "stan", "build", "m2_ame.rds")))
  specs <- data.frame(model = c("M1a", "M1b", "M2a", "M2b"),
                      family = c("m1", "m1", "m2", "m2"),
                      joint = c(0L, 1L, 0L, 1L), stringsAsFactors = FALSE)

  recovery_rows <- list(); conv_rows <- list(); class_rows <- list(); loo_rows <- list(); ame_rows <- list()
  fits_by_truth <- list()

  for (truth_family in c("M1b", "M2b")) {
    sim <- simulate_current(truth_family, seed = V2_SEED, step_sd = STEP_SD)
    W_cols <- sim$cols$W; X_cols <- sim$cols$X; Z_cols <- sim$cols$Z
    logf("truth = %s (%s world, step SD %.2f): S=%d N1=%d N2=%d E=%d (events with weight 1/2: %d); masked state-years %d (IV-only rows %d, stage-2-only rows %d); mu0 = (%.2f, %.2f); x_ref = (%.3f, %.3f)",
         truth_family, if (WORLD == "") "default" else WORLD, STEP_SD, sim$meta$S, sim$meta$N1, sim$meta$N2, sim$meta$E, 2 * sim$meta$n_coalition,
         sim$meta$n_masked_sy, sim$meta$n_iv_only, sim$meta$n_s2_only, sim$truth$mu0[1], sim$truth$mu0[2], sim$x_ref[1], sim$x_ref[2])
    write.csv(data.frame(state_idx = seq_len(sim$meta$S), type = sim$types, planted = sim$truth$planted,
                         mu_ref1 = sim$truth$mu_state_ref[, 1], mu_ref2 = sim$truth$mu_state_ref[, 2],
                         a1 = sim$truth$a_i[, 1], a2 = sim$truth$a_i[, 2], c = sim$truth$c_i),
              file.path(sim_dir, "results", paste0("truth_", truth_family, "_states", SUF, ".csv")), row.names = FALSE)
    for (k in seq_len(nrow(specs))) {
      sp <- specs[k, ]
      tag <- paste0(truth_family, ":", sp$model)
      fit_path <- file.path(sim_dir, "fits", paste0(truth_family, "_", sp$model, SUF, ".rds"))
      dat <- build_stan_data(sim$sy, sim$alts, sim$ev, W_cols, X_cols, Z_cols, joint = sp$joint,
                             x_ref = sim$x_ref, K = if (sp$family == "m2") 2L else NULL)
      stopifnot(any(dat$in_iv == 0L), any(dat$in_s2 == 0L), any(dat$ev_w < 1))
      reused <- FALSE
      if (!is.null(ONLY) && !(tag %in% ONLY)) {
        ## not selected for (re)fitting: reuse the saved fit if one exists (its results are recomputed below)
        if (!file.exists(fit_path)) next
        saved <- readRDS(fit_path); fit <- saved$fit; el <- saved$elapsed_min; reused <- TRUE
        logf("reusing saved fit %s (%s)", tag, format(file.mtime(fit_path), "%Y-%m-%d %H:%M"))
      } else {
        logf("fitting %s on %s data (chains=%d iter=%d warmup=%d)", sp$model, truth_family, CHAINS, ITER, WARMUP)
        t0 <- Sys.time()
        fit <- sampling(models[[sp$family]], data = dat, chains = CHAINS, iter = ITER, warmup = WARMUP,
                        seed = V2_SEED, refresh = max(1L, ITER %/% 10L),
                        control = list(adapt_delta = 0.9, max_treedepth = 10),
                        pars = c("z_re", "z_b", "V_raw", "U_free", "L_Om_joint", "L_Om", "re"), include = FALSE)
        el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
        logf("  done in %.1f min", el)
        saveRDS(list(fit = fit, dat = dat, spec = sp, truth_family = truth_family, elapsed_min = el, fitted_at = Sys.time()),
                fit_path)
      }

      ## ---- convergence ------------------------------------------------------------------
      sdiag <- get_sampler_params(fit, inc_warmup = FALSE)
      n_div <- sum(sapply(sdiag, function(m) sum(m[, "divergent__"])))
      n_maxtd <- sum(sapply(sdiag, function(m) sum(m[, "treedepth__"] >= 10)))
      sm <- summary(fit, pars = c("alpha0", "alpha", "tau", "mu0", "beta", "sigma_x", "sd_re", "sigma_b", "gamma",
                                  "Sigma", "Omega", "R_Om", "mu_state", if (sp$family == "m2") c("U", "V")))$summary
      sm <- sm[is.finite(sm[, "Rhat"]), , drop = FALSE]     # structural zeros in U anchor block -> NA
      conv_rows[[tag]] <- data.frame(truth = truth_family, model = sp$model, reused_saved_fit = reused, elapsed_min = round(el, 1),
                                     divergences = n_div, max_treedepth_hits = n_maxtd,
                                     worst_rhat = max(sm[, "Rhat"]), min_ess_bulk = min(sm[, "n_eff"]),
                                     frac_rhat_gt_1.05 = mean(sm[, "Rhat"] > 1.05))
      ## ---- recovery and 90% coverage ---------------------------------------------------------
      tr <- sim$truth
      post <- rstan::extract(fit, pars = c("alpha0", "alpha", "tau", "mu0", "beta", "Sigma", "Omega", "gamma", "mu_state"))
      q <- function(x) c(mean = mean(x), lo = unname(quantile(x, 0.05)), hi = unname(quantile(x, 0.95)))
      rows <- list()
      add <- function(name, draws, truth) rows[[length(rows) + 1]] <<- data.frame(param = name, truth = truth, t(q(draws)))
      add("alpha0", post$alpha0, tr$alpha0)
      for (p in seq_along(tr$alpha)) add(sprintf("alpha[%d] %s", p, W_cols[p]), post$alpha[, p], tr$alpha[p])
      add("tau", post$tau, tr$tau)
      for (d in 1:2) add(sprintf("mu0[%d]", d), post$mu0[, d], tr$mu0[d])
      for (p in 1:2) for (d in 1:2) add(sprintf("beta[%d,%d] %s", p, d, X_cols[p]), post$beta[, p, d], tr$beta[p, d])
      for (p in seq_along(tr$gamma)) add(sprintf("gamma[%d] %s", p, Z_cols[p]), post$gamma[, p], tr$gamma[p])
      add("Sigma[1,1]", post$Sigma[, 1, 1], tr$Sigma[1, 1]); add("Sigma[2,2]", post$Sigma[, 2, 2], tr$Sigma[2, 2])
      add("Sigma[1,2]", post$Sigma[, 1, 2], tr$Sigma[1, 2])
      for (i in 1:3) for (j in i:3) add(sprintf("Omega[%d,%d]", i, j), post$Omega[, i, j], tr$Omega[i, j])
      rec <- bind_rows(rows) |> mutate(truth_family = truth_family, model = sp$model,
                                        covered = truth >= lo & truth <= hi, abs_err = abs(mean - truth),
                                        fixed_zero = (model %in% c("M1a", "M2a")) & param %in% c("Sigma[1,2]", "Omega[2,3]"))
      recovery_rows[[tag]] <- rec
      ## ---- type classification from posterior mu_state at x_ref (the initiator mean), r0 = 0.5; truth = the true mu_i classified alike --
      cls <- classify_types(post$mu_state, r0 = R0)            # draws x S x 2 -> per-state table
      cls$true_type <- sim$types
      cls$state_idx <- seq_len(sim$meta$S)
      acc <- mean(cls$map_type == cls$true_type)
      class_rows[[tag]] <- data.frame(truth = truth_family, model = sp$model, r0 = R0, accuracy = acc,
                                      n_states = sim$meta$S,
                                      mean_abs_err_mu = mean(abs(apply(post$mu_state, c(2, 3), mean) - tr$mu_state_ref)))
      write.csv(cls, file.path(sim_dir, "results", paste0("classification_", truth_family, "_", sp$model, SUF, ".csv")), row.names = FALSE)
      ## ---- LOO: stage 2, stage 1 (conditional) and the joint state-year LOO (B1) ------------------------------------------
      ll2 <- loo::extract_log_lik(fit, "log_lik_stage2", merge_chains = FALSE)
      ll1 <- loo::extract_log_lik(fit, "log_lik_stage1", merge_chains = FALSE)
      loo2 <- loo::loo(ll2, r_eff = loo::relative_eff(exp(ll2)))
      loo1 <- loo::loo(ll1, r_eff = loo::relative_eff(exp(ll1)))
      looj <- loo_joint_sy(fit, dat)
      loo_rows[[tag]] <- data.frame(truth = truth_family, model = sp$model,
                                    elpd_stage2 = loo2$estimates["elpd_loo", 1], se_stage2 = loo2$estimates["elpd_loo", 2],
                                    bad_k_stage2 = sum(loo2$diagnostics$pareto_k > 0.7),
                                    elpd_stage1_conditional = loo1$estimates["elpd_loo", 1], se_stage1_conditional = loo1$estimates["elpd_loo", 2],
                                    bad_k_stage1_conditional = sum(loo1$diagnostics$pareto_k > 0.7),
                                    elpd_joint_stateyear = looj$estimates["elpd_loo", 1], se_joint_stateyear = looj$estimates["elpd_loo", 2],
                                    bad_k_joint_stateyear = sum(looj$diagnostics$pareto_k > 0.7))
      fits_by_truth[[truth_family]][[sp$model]] <- list(loo1 = loo1, loo2 = loo2, looj = looj)
      ## ---- AME recovery ---------------------------------------------------------------------------
      if (sp$family == "m2") {
        Ud <- rstan::extract(fit, pars = "U", permuted = FALSE)   # iter x chain x (S*K)
        Vd <- rstan::extract(fit, pars = "V", permuted = FALSE)
        S <- dat$S; K <- dat$K
        U4 <- array(Ud, c(dim(Ud)[1], dim(Ud)[2], S, K)); V4 <- array(Vd, c(dim(Vd)[1], dim(Vd)[2], S, K))
        stopifnot(max(abs(apply(U4, c(1, 2, 4), mean))) < 1e-8, max(abs(apply(V4, c(1, 2, 4), mean))) < 1e-8)  # centred U, V
        al <- procrustes_align(U4, V4)
        rot_dev <- max(sapply(al$rotations, function(Rm) max(abs(Rm - diag(K)))))
        UV_draws <- array(NA_real_, c(dim(U4)[1] * dim(U4)[2], S, S)); ii <- 0
        for (ch in seq_len(dim(U4)[2])) for (it in seq_len(dim(U4)[1])) {
          ii <- ii + 1; UV_draws[ii, , ] <- matrix(al$U[it, ch, , ], S, K) %*% t(matrix(al$V[it, ch, , ], S, K))
        }
        UV_mean <- apply(UV_draws, c(2, 3), mean)
        off <- row(UV_mean) != col(UV_mean)
        if (!is.null(tr$UV)) {
          ## the identified quantity is the doubly-centred u'v (column means of U and V removed; the
          ## U-mean part is absorbed by b_j in the likelihood), so compare with that version of the truth
          UV_true <- sweep(tr$U, 2, colMeans(tr$U)) %*% t(sweep(tr$V, 2, colMeans(tr$V)))
          cor_uv <- cor(UV_mean[off], UV_true[off])
          lo <- apply(UV_draws, c(2, 3), quantile, 0.05); hi <- apply(UV_draws, c(2, 3), quantile, 0.95)
          cov_uv <- mean((UV_true >= lo & UV_true <= hi)[off])
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
      cj <- loo::loo_compare(f[[pair[1]]]$looj, f[[pair[2]]]$looj)
      cmp_rows[[length(cmp_rows) + 1]] <- data.frame(
        truth = tf, comparison = paste(pair, collapse = " vs "),
        stage2_elpd_diff_b_minus_a = f[[pair[2]]]$loo2$estimates["elpd_loo", 1] - f[[pair[1]]]$loo2$estimates["elpd_loo", 1],
        stage2_se_diff = c2[2, "se_diff"],
        stage1_conditional_elpd_diff_b_minus_a = f[[pair[2]]]$loo1$estimates["elpd_loo", 1] - f[[pair[1]]]$loo1$estimates["elpd_loo", 1],
        stage1_conditional_se_diff = c1[2, "se_diff"],
        joint_stateyear_elpd_diff_b_minus_a = f[[pair[2]]]$looj$estimates["elpd_loo", 1] - f[[pair[1]]]$looj$estimates["elpd_loo", 1],
        joint_stateyear_se_diff = cj[2, "se_diff"])
    }
  }
  ## misclassification between a and b
  mis_rows <- list()
  for (tf in c("M1b", "M2b")) for (fam in c("M1", "M2")) {
    fa <- file.path(sim_dir, "results", sprintf("classification_%s_%sa%s.csv", tf, fam, SUF))
    fb <- file.path(sim_dir, "results", sprintf("classification_%s_%sb%s.csv", tf, fam, SUF))
    if (file.exists(fa) && file.exists(fb)) {
      a <- read.csv(fa); b <- read.csv(fb)
      mis_rows[[length(mis_rows) + 1]] <- data.frame(truth = tf, family = fam,
        disagreement_a_vs_b = mean(a$map_type != b$map_type),
        acc_a = mean(a$map_type == a$true_type), acc_b = mean(b$map_type == b$true_type))
    }
  }
  write.csv(bind_rows(recovery_rows), file.path(res_dir, paste0("recovery", SUF, ".csv")), row.names = FALSE)
  write.csv(bind_rows(conv_rows), file.path(res_dir, paste0("convergence", SUF, ".csv")), row.names = FALSE)
  write.csv(bind_rows(class_rows), file.path(res_dir, paste0("classification_accuracy", SUF, ".csv")), row.names = FALSE)
  write.csv(bind_rows(loo_rows), file.path(res_dir, paste0("loo", SUF, ".csv")), row.names = FALSE)
  if (length(cmp_rows)) write.csv(bind_rows(cmp_rows), file.path(res_dir, paste0("loo_compare", SUF, ".csv")), row.names = FALSE)
  if (length(mis_rows)) write.csv(bind_rows(mis_rows), file.path(res_dir, paste0("misclassification_a_vs_b", SUF, ".csv")), row.names = FALSE)
  if (length(ame_rows)) write.csv(bind_rows(ame_rows), file.path(res_dir, paste0("ame_recovery", SUF, ".csv")), row.names = FALSE)
  logf("run_sim task %s finished; results in %s", TASK, res_dir)
}

## ============================================================================
## 1d: the one-dimensional model M1m (D_mode = 1) -- M1m (right axis), M1r (wrong axis), M1b (both axes) on the same data
## ============================================================================
if (TASK == "1d") {
  sm <- readRDS(file.path(v2_root(), "stan", "build", "m1_latent_additive.rds"))
  sim <- simulate_current("M1m", S = 30L, Tt = 12L, seed = V2_SEED)
  logf("1-D world: S=%d N1=%d N2=%d E=%d (weight-1/2 events %d; IV-only rows %d, stage-2-only rows %d)", sim$meta$S, sim$meta$N1, sim$meta$N2, sim$meta$E, 2 * sim$meta$n_coalition, sim$meta$n_iv_only, sim$meta$n_s2_only)
  tr <- sim$truth
  W_cols <- sim$cols$W; X_cols <- sim$cols$X; Z_cols <- sim$cols$Z
  specs <- list(M1m = list(D_mode = 1L, joint = 0L), M1r = list(D_mode = 2L, joint = 0L), M1b = list(D_mode = 3L, joint = 1L))
  rec <- list(); conv <- list(); cls_rows <- list(); loos <- list()
  q <- function(x) c(mean = mean(x), lo = unname(quantile(x, 0.05)), hi = unname(quantile(x, 0.95)))
  for (nm in names(specs)) {
    sp <- specs[[nm]]
    dat <- build_stan_data(sim$sy, sim$alts, sim$ev, W_cols, X_cols, Z_cols, joint = sp$joint, x_ref = sim$x_ref, D_mode = sp$D_mode)
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
      for (p in 1:2) add(sprintf("beta[%d,1]", p), b1[, p], if (sp$D_mode == 1L) tr$beta[p, 1] else NA)
      add("Sigma[1,1]", as.numeric(post$Sigma), if (sp$D_mode == 1L) tr$sigma_x[1]^2 else NA)
      add("Omega[1,1] (c)", post$Omega[, 1, 1], tr$Omega[1, 1]); add("Omega[2,2] (a)", post$Omega[, 2, 2], if (sp$D_mode == 1L) tr$Omega[2, 2] else NA)
      add("Omega[1,2] (c,a)", post$Omega[, 1, 2], if (sp$D_mode == 1L) tr$Omega[1, 2] else NA)
    } else {
      add("mu0[1]", post$mu0[, 1], tr$mu0[1])
      for (p in 1:2) add(sprintf("beta[%d,1]", p), post$beta[, p, 1], tr$beta[p, 1])
      add("Sigma[1,1]", post$Sigma[, 1, 1], tr$sigma_x[1]^2); add("Sigma[1,2]", post$Sigma[, 1, 2], 0)
      add("Omega[1,1] (c)", post$Omega[, 1, 1], tr$Omega[1, 1]); add("Omega[2,2] (a1)", post$Omega[, 2, 2], tr$Omega[2, 2]); add("Omega[1,2] (c,a1)", post$Omega[, 1, 2], tr$Omega[1, 2])
    }
    for (p in seq_along(tr$gamma)) add(sprintf("gamma[%d]", p), post$gamma[, p], tr$gamma[p])
    rec[[nm]] <- bind_rows(rows) |> mutate(model = nm, covered = truth >= lo & truth <= hi, abs_err = abs(mean - truth))
    mu <- mu_state_draws(fit)
    if (D == 1L) {
      cls <- classify_types(mu, r0 = R0)
      acc <- if (sp$D_mode == 1L) mean(cls$map_type == sim$types_1d) else NA
      cls_rows[[nm]] <- data.frame(model = nm, axis = c("material", "relational")[sp$D_mode], r0 = R0, accuracy_3level = acc,
                                   mean_abs_err_mu = if (sp$D_mode == 1L) mean(abs(cls$mu_mean - tr$mu_state_ref[, 1])) else NA)
    } else {
      ## the material axis of the 2-D fit, classified by the 3-level rule on the posterior MEAN (review fix 2026-10-25, B5: the old
      ## code compared a LOGICAL vector with -0.8, so "stronger" could never occur)
      mu_mean1 <- apply(mu[, , 1], 2, mean)
      cls3 <- ifelse(mu_mean1 > R0, "weaker", ifelse(mu_mean1 < -R0, "stronger", "peer"))
      stopifnot(is.character(cls3), all(cls3 %in% c("weaker", "peer", "stronger")))
      cls_rows[[nm]] <- data.frame(model = nm, axis = "both", r0 = R0, accuracy_3level = mean(cls3 == sim$types_1d),
                                   mean_abs_err_mu = mean(abs(mu_mean1 - tr$mu_state_ref[, 1])))
    }
    ll2 <- loo::extract_log_lik(fit, "log_lik_stage2", merge_chains = FALSE); l2 <- loo::loo(ll2, r_eff = loo::relative_eff(exp(ll2)))
    ll1 <- loo::extract_log_lik(fit, "log_lik_stage1", merge_chains = FALSE); l1 <- loo::loo(ll1, r_eff = loo::relative_eff(exp(ll1)))
    lj <- loo_joint_sy(fit, dat)
    loos[[nm]] <- data.frame(model = nm, elpd_stage2 = l2$estimates["elpd_loo", 1], se_stage2 = l2$estimates["elpd_loo", 2],
                             elpd_stage1_conditional = l1$estimates["elpd_loo", 1], se_stage1_conditional = l1$estimates["elpd_loo", 2],
                             elpd_joint_stateyear = lj$estimates["elpd_loo", 1], se_joint_stateyear = lj$estimates["elpd_loo", 2])
    rm(fit); gc()
  }
  write.csv(bind_rows(rec), file.path(res_dir, "sim1d_recovery.csv"), row.names = FALSE)
  write.csv(bind_rows(conv), file.path(res_dir, "sim1d_convergence.csv"), row.names = FALSE)
  write.csv(bind_rows(cls_rows), file.path(res_dir, "sim1d_classification.csv"), row.names = FALSE)
  write.csv(bind_rows(loos), file.path(res_dir, "sim1d_loo.csv"), row.names = FALSE)
  print(bind_rows(rec) |> filter(model == "M1m")); print(bind_rows(cls_rows)); print(bind_rows(loos))
  logf("run_sim task 1d finished")
}

## ============================================================================
## measurement: the measurement tests on a 30-state world (truth M1b)
## ============================================================================
if (TASK == "measurement") {
  models <- list(m1 = readRDS(file.path(v2_root(), "stan", "build", "m1_latent_additive.rds")),
                 m2 = readRDS(file.path(v2_root(), "stan", "build", "m2_ame.rds")))
  sim <- simulate_current("M1b", S = 30L, Tt = 12L, seed = V2_SEED)
  W_cols <- sim$cols$W; X_cols <- sim$cols$X; Z_cols <- sim$cols$Z; x_ref <- sim$x_ref
  logf("small world: S=%d N1=%d N2=%d E=%d (weight-1/2 events %d; IV-only rows %d, stage-2-only rows %d)", sim$meta$S, sim$meta$N1, sim$meta$N2, sim$meta$E, 2 * sim$meta$n_coalition, sim$meta$n_iv_only, sim$meta$n_s2_only)
  sy <- sim$sy; alts <- sim$alts; ev <- sim$ev

  ## ---- baseline fits a and b (M1) ------------------------------------------------------------------------------
  fit_ab <- function(joint) {
    dat <- build_stan_data(sy, alts, ev, W_cols, X_cols, Z_cols, joint = joint, x_ref = x_ref)
    list(fit = sampling(models$m1, data = dat, chains = CHAINS, iter = ITER, warmup = WARMUP, seed = V2_SEED, refresh = 0,
                        control = list(adapt_delta = 0.9), pars = excl, include = FALSE), dat = dat)
  }
  t0 <- Sys.time(); fa <- fit_ab(0L); fb <- fit_ab(1L); fit_a <- fa$fit; fit_b <- fb$fit
  logf("baseline a and b fitted in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

  ## ---- 1. LOO a vs b (stage 1 and 2 separately) -----------------------------------------------------------------
  lc <- loo_compare_ab(fit_a, fit_b, fa$dat, fb$dat); print(lc)                   # stage 2, stage 1 conditional, joint state-year (B1)
  write.csv(lc, file.path(res_dir, "measurement_loo_ab.csv"), row.names = FALSE)

  ## ---- 2. misclassification a vs b ------------------------------------------------------------------------------
  mc <- misclassification_ab(fit_a, fit_b, r0 = R0)
  cls_a <- classify_types(mu_state_draws(fit_a), r0 = R0); cls_b <- classify_types(mu_state_draws(fit_b), r0 = R0)
  mis <- data.frame(rate_a_vs_b = mc$rate, mean_abs_prob_shift = mc$mean_abs_prob_shift,
                    accuracy_a = mean(cls_a$map_type == sim$types), accuracy_b = mean(cls_b$map_type == sim$types))
  print(mis); write.csv(mis, file.path(res_dir, "measurement_misclassification.csv"), row.names = FALSE)
  write.csv(as.data.frame(mc$table), file.path(res_dir, "measurement_misclassification_table.csv"), row.names = FALSE)

  ## ---- 3. temporal holdout --------------------------------------------------------------------------------------
  cutoff <- sort(unique(sy$year))[9]              # train on the first 9 usable years, test on the last 3
  t0 <- Sys.time()
  freeze <- list(positions = sim$levels, gap_raw = c(dM = "dM_logratio", dR = "dR_pr_logratio"), levels = list(dM = "M_lvl", dR = "R_lvl"),
                 change_raw = c(dM_i = "dM_raw_h", dR_i = "dRpr_raw_h"), target_change = c(dM_i = "dM_j_h", dR_i = "dRpr_j_h"),
                 ctrl_cent = c("ongoing_tm1", "log_prior_mids", "dem_tm1", "log_n_system_tm1"), z_cent = c("log_dy_prior_mids", "dM_j_h", "dRpr_j_h"))
  th <- temporal_holdout(sy, alts, ev, W_cols, X_cols, Z_cols, cutoff = cutoff, family = "m1", stan_models = models,
                         chains = CHAINS, iter = ITER, warmup = WARMUP, n_draws = 100, x_ref = x_ref, freeze = freeze)
  logf("temporal holdout fitted in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
  print(th$summary); write.csv(th$summary, file.path(res_dir, "measurement_temporal_holdout.csv"), row.names = FALSE)

  ## ---- 4. permutation test (rescore mode; refit mode is available but too slow here) -----------------------------
  pt <- permutation_test(sy, alts, ev, W_cols, X_cols, Z_cols, fit_a, fit_b, n_perm = 20, mode = "rescore", n_draws = 100, x_ref = x_ref)
  print(pt$summary); write.csv(pt$summary, file.path(res_dir, "measurement_permutation.csv"), row.names = FALSE)
  write.csv(data.frame(perm = seq_along(pt$permuted_advantage), adv = pt$permuted_advantage, a = pt$permuted_a, b = pt$permuted_b),
            file.path(res_dir, "measurement_permutation_draws.csv"), row.names = FALSE)

  ## ---- 5. systematic measurement error ---------------------------------------------------------------------------
  ## "recent conflict" at the POSITION-year: the state initiated in that year, so its measured relational level that year is inflated
  ## (deltas in log points of R, like the real-data run); the gaps, changes and level column are rebuilt from the perturbed levels (B8)
  recent <- sy |> transmute(state_idx, year, recent = y1)
  t0 <- Sys.time()
  me <- measurement_error_check(sy, alts, ev, W_cols, X_cols, Z_cols, recent = recent, deltas = c(0.25, 0.5, 1.0),
                                baseline_cls = cls_b, levels = sim$levels, h = 1L, gap_scale_R = sim$truth$gap_R,
                                family = "m1", joint = 1L, stan_models = models,
                                chains = CHAINS, iter = ITER, warmup = WARMUP, r0 = R0, x_ref = x_ref)
  logf("measurement-error refits in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
  me$baseline_accuracy <- mean(cls_b$map_type == sim$types)
  print(me); write.csv(me, file.path(res_dir, "measurement_error_inflation.csv"), row.names = FALSE)

  ## Benchmarks on the simulated frames are --task=benchmarks (separate R session: the brms compile needs ~4 GB).
  logf("run_sim task measurement finished")
}

## ============================================================================
## benchmarks: the archived benchmarks (directed-dyad logit, old gap regression) on the simulated frames
## ============================================================================
if (TASK == "benchmarks") {
  suppressPackageStartupMessages(library(brms))
  source(file.path(sim_dir, "..", "archive", "R", "benchmarks.R"))   # archived option; the assembly functions are still exercised here
  sim <- simulate_current("M1b", S = 30L, Tt = 12L, seed = V2_SEED)          # the current-model world (2026-10-25, B5)
  sy <- sim$sy; alts <- sim$alts |> filter(in_s2 == 1L); ev <- sim$ev; X_cols <- sim$cols$X; Z_cols <- sim$cols$Z
  ev <- ev |> mutate(alt_row = match(paste(sim$alts$sy_row[alt_row], sim$alts$target_idx[alt_row]), paste(alts$sy_row, alts$target_idx)))
  stopifnot(!anyNA(ev$alt_row))
  ## ---- 6. benchmarks on the simulated frames (assembly + short brms fits) ---------------------------------------------
  ev_keys <- data.frame(state_idx = alts$state_idx[ev$alt_row], target_idx = alts$target_idx[ev$alt_row], year = alts$year[ev$alt_row])
  ## cinc_ratio = plogis(raw log material gap) (B10); simulate_current() carries dM_logratio like the real alternatives files
  d1 <- assemble_dyad_logit(alts, ev_keys, id_cols = c("state_idx", "target_idx", "year"), z_cols = Z_cols, logratio_col = "dM_logratio")
  stopifnot(sum(d1$y) == nrow(ev))
  d1$off <- 0
  t0 <- Sys.time()
  b1 <- brm(bf(as.formula(paste("y ~ cinc_ratio + dR +", paste(Z_cols, collapse = " + "), "+ offset(off) + (1 | state_idx) + (1 | target_idx)"))), data = d1, family = bernoulli(),
            prior = c(set_prior("normal(0, 1)", class = "b"), set_prior("student_t(3, 0, 1.5)", class = "sd")),
            chains = CHAINS, cores = CHAINS, iter = ITER, warmup = WARMUP, seed = V2_SEED, refresh = 0)
  d2 <- assemble_gap_regression(sy, alts, ev, X_cols)
  b2 <- brm(gap_formula(X_cols), data = d2, prior = set_prior("normal(0, 1)", class = "b"),
            chains = CHAINS, cores = CHAINS, iter = ITER, warmup = WARMUP, seed = V2_SEED, refresh = 0)
  logf("benchmarks fitted in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
  fe1 <- as.data.frame(fixef(b1)); fe1$model <- "dyad_logit"; fe1$term <- rownames(fe1)
  fe2 <- as.data.frame(fixef(b2)); fe2$model <- "gap_mvbind"; fe2$term <- rownames(fe2)
  write.csv(bind_rows(fe1, fe2), file.path(res_dir, "benchmarks_sim.csv"), row.names = FALSE)
  ## truth for the gap regression's dM ~ X1 slope is NOT beta[1,1]: the realised target position is a
  ## softmax draw around mu, so the slope recovers beta only up to shrinkage; recorded for the report.
  logf("run_sim task benchmarks finished")
}
