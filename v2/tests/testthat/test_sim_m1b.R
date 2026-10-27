## M1b simulation recovery mode of sim/sim_context_b.R (PI request 2026-10-27, after the M1c recovery bias): the --model switch, the
## summariser's --model=M1b path, the 07 sbatch, and the M1c default left unchanged (wiring checks; the real-data pilots run when the
## data and the compiled programs are present).
suppressPackageStartupMessages({ library(testthat); library(dplyr) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)

## synthetic M1b replicate files written the way sim_context_b.R --model=M1b writes them
m1b_pars <- c("mu0[1]", "mu0[2]", "sigma_x[1]", "sigma_x[2]", "rho_x[1]", "Sigma[1,1]", "Sigma[1,2]", "Sigma[2,2]", "beta[1,1]", "gamma[1]", "alpha0", "alpha[1]", "tau", "sd_re[1]", "sigma_b")
m1b_blk <- c(rep("stage_2_core", 8), "stage_2", "stage_2", rep("stage_1", 3), "re_scales", "re_scales")
m1b_truth <- c(-1, 0.3, 1.2, 0.9, 0.4, 1.44, 0.432, 0.81, 0.3, 1, -2, 0.2, 0.5, 0.9, 0.6)
write_m1b_rep <- function(r, dir, iter = 500L, covered_all = FALSE, prefix = "m1b", model_col = TRUE) {
  set.seed(2000 + r)
  m <- m1b_truth + rnorm(length(m1b_truth), 0, 0.2); sdv <- rep(0.3, length(m1b_truth)); lo <- m - 1.645 * sdv; hi <- m + 1.645 * sdv
  if (covered_all) { lo <- m1b_truth - 1; hi <- m1b_truth + 1 }
  rec <- data.frame(parameter = m1b_pars, truth = m1b_truth, mean = m, sd = sdv, lo = lo, hi = hi, Rhat = 1.001, n_eff = 400, covered = m1b_truth >= lo & m1b_truth <= hi, bias = m - m1b_truth, z = (m - m1b_truth) / sdv,
                    block = m1b_blk, recovery_target = TRUE, stringsAsFactors = FALSE)
  rec <- rbind(rec, data.frame(parameter = "mu_state 90% coverage over states with events x axes", truth = 0.9, mean = 0.91, sd = NA, lo = NA, hi = NA, Rhat = NA, n_eff = NA, covered = NA, bias = NA, z = NA, block = "mu_state", recovery_target = TRUE))
  rec$rep <- r; rec$seed <- 123L + r; rec$n_init <- 250L; rec$n_noninit <- 250L; rec$chains <- 4L; rec$iter <- iter; rec$warmup <- iter %/% 2L; rec$elapsed_min <- 20; rec$divergences <- 0L; rec$treedepth_hits <- 0L; rec$sample <- "fmid"; rec$lag <- 3L
  if (model_col) rec$model <- "M1b"
  write.csv(rec, file.path(dir, sprintf("%s_recovery_r%d.csv", prefix, r)), row.names = FALSE)
  write.csv(data.frame(group = "all", n_events = 400L, w_sum = 280, n_states = 120L, n_state_years = 260L, rep = r, seed = 123L + r,
                       N1 = 500L, N2 = 60000L, E = 400L, S = 217L, sum_y1 = 260L, alpha0_true = -2, n_states_with_events = 120L, n_pattern_redraws = 0L, n_dup_events = 9L),
            file.path(dir, sprintf("%s_design_r%d.csv", prefix, r)), row.names = FALSE)
}

test_that("sim_context_b.R --model: M1c is the default and its literals are untouched; the M1b branch has its own directories, files, truth, program and block", {
  src <- readLines(file.path(root, "sim", "sim_context_b.R"))
  expect_true(any(grepl('MODEL <- get_arg(args, "model", "M1c")', src, fixed = TRUE)))
  expect_true(any(grepl('stop("--model must be M1c (default) or M1b")', src, fixed = TRUE)))
  ## the M1c path: the same directories, file names and program as before the switch (byte-identical outputs)
  expect_true(any(grepl('fit_dir <- file.path(sim_dir, "fits", if (SMOKE) "context_b_smoke" else "context_b")', src, fixed = TRUE)))
  expect_true(any(grepl('if (M1C) stopifnot(rec_file == file.path(res_dir, sprintf("context_b_recovery_%s.csv", tag)), dsn_file == file.path(res_dir, sprintf("context_b_design_%s.csv", tag)))', src, fixed = TRUE)))
  expect_true(any(grepl('file.path(fit_dir, sprintf("context_b_sim_%s.rds", tag))', src, fixed = TRUE)))
  expect_true(any(grepl('u <- util_rows(rows_all, n, if (M1C) pat$grp[k] else 0L)', src, fixed = TRUE)))        # M1c keeps the group softmax, M1b the single one
  expect_true(any(grepl('shift <- rbind(c(0, 0), TRUTH$delta_lead, TRUTH$delta_join)', src, fixed = TRUE)))
  expect_true(any(grepl("if (!M1C) rec$model <- MODEL", src, fixed = TRUE)))                                    # the model column is M1b-only: the M1c CSV is unchanged
  expect_false(any(grepl("fit_stem", src, fixed = TRUE)))
  ## the production directories are created only after the PILOT exit
  expect_equal(sum(grepl("dir.create(", src, fixed = TRUE)), 1L)
  expect_true(which(grepl("dir.create(fit_dir, showWarnings = FALSE, recursive = TRUE); dir.create(res_dir", src, fixed = TRUE)) > which(grepl('if (PILOT) { logf("PILOT done; no recovery table")', src, fixed = TRUE)))
  ## the M1b path
  expect_true(any(grepl('fit_dir <- file.path(sim_dir, "fits", if (SMOKE) "m1b_recovery_smoke" else "m1b_recovery")', src, fixed = TRUE)))
  expect_true(any(grepl('res_dir <- file.path(sim_dir, "results", if (SMOKE) "m1b_recovery_smoke" else "m1b_recovery")', src, fixed = TRUE)))
  expect_true(any(grepl('prefix <- "m1b"', src, fixed = TRUE))); expect_true(any(grepl('sprintf("m1b_sim_%s.rds", tag)', src, fixed = TRUE)))
  expect_true(any(grepl('A <- assemble_cell_spec(SAMPLE, LAG, parse_spec(character(0)), model = MODEL)', src, fixed = TRUE)))
  expect_true(any(grepl('stopifnot(!"multi_mid_jan" %in% A$W_cols, identical(A$W_cols, c("dM_i", "dR_i", w_controls())))', src, fixed = TRUE)))
  ## paired design (audit of 39232d3): ONE leader-first selection for both models, no model branch around it
  expect_equal(sum(grepl("keep_init <- ", src, fixed = TRUE)), 1L); expect_equal(sum(grepl("keep_lead <- ", src, fixed = TRUE)), 1L)
  expect_true(any(grepl('keep_init <- sort(c(keep_lead, sample(setdiff(init_sy, keep_lead), N_INIT - length(keep_lead))))', src, fixed = TRUE)))
  expect_false(any(grepl("sample(init_sy", src, fixed = TRUE)))
  expect_true(any(grepl('TRUTH$sigma_x <- unname(TRUTH$sigma_g["bilateral", ])', src, fixed = TRUE)))             # the M1c bilateral truth
  expect_true(any(grepl('shift <- rbind(c(0, 0)); Sinv <- list(solve(TRUTH$Sigma))', src, fixed = TRUE)))
  expect_true(any(grepl('stage2_loglik_R(dat_sim, list(mu0 = TRUTH$mu0, beta = TRUTH$beta, a_i = TRUTH$a_i, b_j = TRUTH$b_j, gamma = TRUTH$gamma, sigma_x = TRUTH$sigma_x, rho = TRUTH$rho))', src, fixed = TRUE)))
  expect_true(any(grepl('family <- stan_program(MODEL)', src, fixed = TRUE)))
  expect_true(any(grepl('c("mu0", "sigma_x", "rho_x", "Sigma", "beta", "gamma", "tau", "alpha0", "alpha", "sd_re", "sigma_b")', src, fixed = TRUE)))
  expect_true(any(grepl('sm_tab <- sm_tab[rownames(sm_tab) != "Sigma[2,1]", , drop = FALSE]', src, fixed = TRUE)))
  expect_true(any(grepl('grepl("^(mu0|sigma_x|rho_x|Sigma)", rec$parameter) ~ "stage_2_core"', src, fixed = TRUE)))
  expect_true(any(grepl('data.frame(group = "all", n_events = dat_sim$E', src, fixed = TRUE)))
  expect_true(any(grepl('control = list(adapt_delta = SAMPLER$M1$adapt_delta, max_treedepth = SAMPLER$M1$max_treedepth)', src, fixed = TRUE)))
  ## the sbatch: open account, array 0-19, --model=M1b, the M1 production chains, the env overrides, the summariser's M1b call
  sb <- readLines(file.path(root, "slurm", "07_sim_m1b.sbatch"))
  expect_true(any(grepl("^#SBATCH --account=open", sb))); expect_true(any(grepl("^#SBATCH --array=0-19", sb))); expect_true(any(grepl("^#SBATCH --job-name=v2_sim_m1b", sb)))
  expect_true(any(grepl('Rscript v2/sim/sim_context_b.R --model=M1b --rep="\\$r"', sb))); expect_true(any(grepl("--chains=4 --iter=2000 --warmup=1000", sb, fixed = TRUE)))
  expect_true(any(grepl('--n_init="${N_INIT:-250}" --n_noninit="${N_NONINIT:-250}"', sb, fixed = TRUE)))
  expect_true(any(grepl("summarise_context_b.R --model=M1b", sb, fixed = TRUE))); expect_true(any(grepl("r=\\$\\{SLURM_ARRAY_TASK_ID", sb)))
  ## the 06 sbatch still runs M1c (no --model argument)
  sb6 <- readLines(file.path(root, "slurm", "06_sim_context_b.sbatch")); expect_false(any(grepl("--model", sb6, fixed = TRUE)))
  expect_true(any(grepl('Rscript v2/sim/sim_context_b.R --rep="\\$r"', sb6)))
})

test_that("summarise_context_b.R --model=M1b: the m1b_recovery files -> m1b_coverage*.csv with the stage_2_core block; M1c files are not read; M1c default unchanged; refusals", {
  tmp <- tempfile("simm1b_v2_"); v2 <- file.path(tmp, "v2"); dir.create(v2, recursive = TRUE)
  file.copy(file.path(root, "R"), v2, recursive = TRUE); dir.create(file.path(v2, "sim")); file.copy(list.files(file.path(root, "sim"), pattern = "\\.R$", full.names = TRUE), file.path(v2, "sim")); dir.create(file.path(v2, "logs"))
  res_b <- file.path(v2, "sim", "results", "m1b_recovery"); dir.create(res_b, recursive = TRUE)
  res_c <- file.path(v2, "sim", "results", "context_b"); dir.create(res_c, recursive = TRUE)
  run <- function(...) { out <- suppressWarnings(system2("Rscript", c(file.path(v2, "sim", "summarise_context_b.R"), ...), stdout = TRUE, stderr = TRUE)); list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  for (r in 0:4) write_m1b_rep(r, res_b)
  expect_true(run("--model=M1b")$status != 0L)                                                     # 5 < 20 replicates by default
  expect_true(run("--model=M1x", "--min_reps=5")$status != 0L)                                        # an unknown model is refused
  expect_true(run("--min_reps=5")$status != 0L)                                                     # the M1c default finds no context_b files: the M1b files are not read
  r5 <- run("--model=M1b", "--min_reps=5"); expect_equal(r5$status, 0L, info = paste(tail(r5$out, 20), collapse = "\n"))
  for (f in c("m1b_coverage.csv", "m1b_coverage_blocks.csv", "m1b_replicates.csv")) expect_true(file.exists(file.path(v2, "sim", "results", f)), info = f)
  for (f in c("context_b_coverage.csv", "context_b_coverage_blocks.csv", "context_b_replicates.csv")) expect_false(file.exists(file.path(v2, "sim", "results", f)), info = f)
  cov <- read.csv(file.path(v2, "sim", "results", "m1b_coverage.csv"), stringsAsFactors = FALSE)
  expect_equal(nrow(cov), length(m1b_pars) + 1L); expect_true(all(cov$n_reps == 5L))
  expect_equal(cov$parameter[seq_along(m1b_pars)], c("mu0[1]", "mu0[2]", "sigma_x[1]", "sigma_x[2]", "rho_x[1]", "Sigma[1,1]", "Sigma[1,2]", "Sigma[2,2]", "beta[1,1]", "gamma[1]", "alpha0", "alpha[1]", "tau", "sd_re[1]", "sigma_b"))
  expect_equal(cov$block[cov$parameter %in% c("mu0[1]", "sigma_x[2]", "rho_x[1]", "Sigma[1,2]")], rep("stage_2_core", 4)); expect_false("stage_b" %in% cov$block)
  all_rec <- do.call(rbind, lapply(0:4, function(r) read.csv(file.path(res_b, sprintf("m1b_recovery_r%d.csv", r)), stringsAsFactors = FALSE)))
  for (pp in m1b_pars) { d <- all_rec[all_rec$parameter == pp, ]; expect_equal(cov$coverage90[cov$parameter == pp], mean(d$covered)); expect_equal(cov$mean_bias[cov$parameter == pp], mean(d$bias)); expect_equal(cov$rmse[cov$parameter == pp], sqrt(mean(d$bias^2))) }
  expect_equal(cov$mean_over_reps[cov$block == "mu_state"], 0.91)
  bl <- read.csv(file.path(v2, "sim", "results", "m1b_coverage_blocks.csv"), stringsAsFactors = FALSE)
  expect_equal(bl$block, c("re_scales", "stage_1", "stage_2", "stage_2_core", "all")); expect_equal(bl$n_params[bl$block == "stage_2_core"], 8L); expect_equal(bl$n_params[bl$block == "all"], length(m1b_pars))
  rt <- read.csv(file.path(v2, "sim", "results", "m1b_replicates.csv")); expect_equal(rt$rep, 0:4); expect_equal(rt$seed, 123:127); expect_equal(rt$n_events, rep(400L, 5)); expect_false("n_leader" %in% names(rt))
  for (r in 0:4) write_m1b_rep(r, res_b, covered_all = TRUE); expect_equal(run("--model=M1b", "--min_reps=5")$status, 0L)
  ## a replicate file without the model column is refused in the M1b path
  write_m1b_rep(0, res_b, model_col = FALSE); expect_true(run("--model=M1b", "--min_reps=5")$status != 0L); write_m1b_rep(0, res_b)
  expect_true(all(read.csv(file.path(v2, "sim", "results", "m1b_coverage.csv"))$coverage90[seq_along(m1b_pars)] == 1))
  ## mixed settings and a duplicated index are refused in the M1b path too
  write_m1b_rep(5, res_b, iter = 2000L); expect_true(run("--model=M1b", "--min_reps=5")$status != 0L); file.remove(file.path(res_b, c("m1b_recovery_r5.csv", "m1b_design_r5.csv")))
  write_m1b_rep(0, res_b); file.copy(file.path(res_b, "m1b_recovery_r0.csv"), file.path(res_b, "m1b_recovery_r00.csv")); expect_true(run("--model=M1b", "--min_reps=5")$status != 0L); file.remove(file.path(res_b, "m1b_recovery_r00.csv"))
  ## M1b-format files dropped into the M1c directory are refused (block / group assertions), so the two designs cannot be pooled
  for (r in 0:4) write_m1b_rep(r, res_c, prefix = "context_b"); expect_true(run("--min_reps=5")$status != 0L)
  unlink(tmp, recursive = TRUE)
})

test_that("sim_context_b.R pilots on the real fmid lag-3 design: M1c (default) and M1b both run through sampling; the M1b design has no groups and no multi_mid_jan", {
  skip_if_not(file.exists(file.path(root, "data", "manifest.csv")) && file.exists(file.path(root, "data", "alts_lag3_all.rds")), "real data files not built")
  for (m in c("m1_context", "m1_latent_additive")) if (!file.exists(file.path(root, "stan", "build", paste0(m, ".rds")))) skip(paste("compiled model missing:", m))
  run <- function(...) { out <- suppressWarnings(system2("Rscript", c(file.path(root, "sim", "sim_context_b.R"), ...), stdout = TRUE, stderr = TRUE)); list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  ## a pilot writes no CSV / rds (quits after sampling), so the production directories stay untouched
  before <- list.files(file.path(root, "sim"), recursive = TRUE); before_dirs <- list.dirs(file.path(root, "sim"), recursive = TRUE)
  c1 <- run("--rep=0", "--pilot=TRUE", "--n_init=60", "--n_noninit=40"); expect_equal(c1$status, 0L, info = paste(tail(c1$out, 30), collapse = "\n"))
  expect_true(any(grepl("M1c simulation recovery on the observed fmid lag-3 design, replicate 0 (seed 123)", c1$out, fixed = TRUE)))
  expect_true(any(grepl("P1=9 P=2", c1$out, fixed = TRUE))); expect_true(any(grepl("groups: bilateral", c1$out, fixed = TRUE))); expect_true(any(grepl("slots (bilateral", c1$out, fixed = TRUE)))
  expect_true(any(grepl("PILOT done; no recovery table", c1$out, fixed = TRUE)))
  b1 <- run("--model=M1b", "--rep=0", "--pilot=TRUE", "--n_init=60", "--n_noninit=40"); expect_equal(b1$status, 0L, info = paste(tail(b1$out, 30), collapse = "\n"))
  expect_true(any(grepl("M1b simulation recovery on the observed fmid lag-3 design, replicate 0 (seed 123)", b1$out, fixed = TRUE)))
  expect_true(any(grepl("P1=8 P=2", b1$out, fixed = TRUE))); expect_true(any(grepl("groups: all", b1$out, fixed = TRUE))); expect_true(any(grepl("(M1b: weight multisets", b1$out, fixed = TRUE)))
  ## the paired design: the M1b pilot keeps exactly the state-years the M1c pilot kept (same design line, same pool size)
  dl <- function(o) sub("^\\[[0-9:]+\\] ", "", grep("^\\[[0-9:]+\\] design: ", o, value = TRUE))
  expect_length(dl(c1$out), 1L); expect_equal(dl(b1$out), dl(c1$out)); expect_true(grepl("with a leader event, 60 kept), 60 initiating kept", dl(b1$out), fixed = TRUE))
  expect_true(any(grepl("PILOT done; no recovery table", b1$out, fixed = TRUE)))
  expect_true(run("--model=M1x", "--rep=0", "--pilot=TRUE")$status != 0L)
  expect_equal(list.files(file.path(root, "sim"), recursive = TRUE), before); expect_equal(list.dirs(file.path(root, "sim"), recursive = TRUE), before_dirs)
})
