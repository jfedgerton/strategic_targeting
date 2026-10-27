## ============================================================================
## summarise_context_b.R -- the coverage table of the M1c simulation replicates (sim/sim_context_b.R, one CSV per replicate)
##
## Usage: Rscript v2/sim/summarise_context_b.R [--min_reps=20]
## Reads v2/sim/results/context_b/context_b_recovery_r*.csv (and the matching _design_ files), refuses fewer than --min_reps replicates or
## replicates with differing run settings, and writes
##   v2/sim/results/context_b_coverage.csv       per parameter: truth (mean over replicates; alpha0 is calibrated per replicate), n_reps,
##                                               coverage90 (share of replicates whose 90% interval covers the truth), mean bias, RMSE,
##                                               mean posterior SD, mean |z|, share of |z| > 2, max Rhat, min n_eff, block
##   v2/sim/results/context_b_coverage_blocks.csv the same pooled over the parameter blocks (stage_b, stage_2, stage_1, re_scales) and overall
##   v2/sim/results/context_b_replicates.csv     one row per replicate: seed, E, sum(y1), events per group, elapsed, divergences, max Rhat
## Numbers are read from the replicate files, never typed. A log goes to v2/logs.
## ============================================================================
suppressPackageStartupMessages({ library(dplyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
sim_dir <- normalizePath(dirname(this_file)); R_DIR <- file.path(sim_dir, "..", "R")
Sys.setenv(V2_ROOT = normalizePath(file.path(sim_dir, "..")))
source(file.path(R_DIR, "spec.R")); source(file.path(R_DIR, "helpers_data.R"))
args <- commandArgs(trailingOnly = TRUE)
MIN_REPS <- as.integer(get_arg(args, "min_reps", 20))
stopifnot(MIN_REPS >= 1)
log_open("summarise_context_b")
res_dir <- file.path(sim_dir, "results", "context_b")
files <- sort(list.files(res_dir, pattern = "^context_b_recovery_r[0-9]+\\.csv$", full.names = TRUE))
if (length(files) < MIN_REPS) stop(sprintf("%d replicate file(s) in %s; at least %d required (--min_reps)", length(files), res_dir, MIN_REPS))
rec <- bind_rows(lapply(files, read.csv, stringsAsFactors = FALSE))
stopifnot(all(c("parameter", "truth", "mean", "sd", "lo", "hi", "covered", "bias", "z", "Rhat", "n_eff", "block", "rep", "seed", "n_init", "n_noninit", "chains", "iter", "warmup",
                "elapsed_min", "divergences", "treedepth_hits", "sample", "lag") %in% names(rec)))
reps <- rec |> distinct(rep, seed, n_init, n_noninit, chains, iter, warmup, sample, lag, elapsed_min, divergences, treedepth_hits)
if (anyDuplicated(reps$rep)) stop("a replicate index appears in more than one file")
if (nrow(distinct(reps, n_init, n_noninit, chains, iter, warmup, sample, lag)) != 1L) stop("the replicates were run with different settings (n_init / n_noninit / chains / iter / warmup / sample / lag); summarise one design at a time")
stopifnot(all(reps$seed == V2_SEED + reps$rep))
## every replicate reports the same parameter set
par_sets <- rec |> group_by(rep) |> summarise(k = paste(sort(parameter), collapse = "|"), .groups = "drop")
if (n_distinct(par_sets$k) != 1L) stop("the replicates report different parameter sets")
dsn <- bind_rows(lapply(sub("recovery", "design", files), function(f) { if (!file.exists(f)) stop("missing design file ", f); read.csv(f, stringsAsFactors = FALSE) }))
rep_tab <- dsn |> group_by(rep, seed, N1, N2, E, S, sum_y1, alpha0_true, n_states_with_events, n_pattern_redraws) |>
  summarise(n_bilateral = n_events[group == "bilateral"], n_leader = n_events[group == "leader"], n_joiner = n_events[group == "joiner"], .groups = "drop") |>
  inner_join(rec |> filter(!is.na(Rhat)) |> group_by(rep) |> summarise(max_rhat = max(Rhat), min_n_eff = min(n_eff), .groups = "drop"), by = "rep") |>
  inner_join(reps |> select(rep, elapsed_min, divergences, treedepth_hits), by = "rep") |> arrange(rep)
stopifnot(nrow(rep_tab) == length(files))
## per parameter (the interval-based rows; the mu_state summary rows are reported as means)
tg <- rec |> filter(!is.na(covered))
cov <- tg |> group_by(parameter, block) |>
  summarise(n_reps = n(), truth = mean(truth), coverage90 = mean(covered), mean_bias = mean(bias), rmse = sqrt(mean(bias^2)), mean_post_sd = mean(sd), mean_abs_z = mean(abs(z)),
            share_abs_z_gt2 = mean(abs(z) > 2), max_rhat = max(Rhat), min_n_eff = min(n_eff), .groups = "drop") |>
  mutate(.o = match(sub("\\[.*$", "", parameter), c("mu0", "delta_lead", "delta_join", "mu0_lead", "mu0_join", "sigma_g", "rho_x", "beta", "gamma", "alpha0", "alpha", "tau", "sd_re", "sigma_b"))) |>
  arrange(.o, parameter) |> select(-.o)
stopifnot(all(cov$n_reps == length(files)))
ms <- rec |> filter(is.na(covered)) |> group_by(parameter, block) |> summarise(n_reps = n(), truth = mean(truth), mean_over_reps = mean(mean), sd_over_reps = sd(mean), .groups = "drop")
blocks <- bind_rows(tg |> group_by(block) |> summarise(n_params = n_distinct(parameter), n_reps = n_distinct(rep), coverage90 = mean(covered), mean_abs_bias = mean(abs(bias)), rmse = sqrt(mean(bias^2)),
                                                       mean_abs_z = mean(abs(z)), share_abs_z_gt2 = mean(abs(z) > 2), max_rhat = max(Rhat), min_n_eff = min(n_eff), .groups = "drop"),
                    tg |> summarise(block = "all", n_params = n_distinct(parameter), n_reps = n_distinct(rep), coverage90 = mean(covered), mean_abs_bias = mean(abs(bias)), rmse = sqrt(mean(bias^2)),
                                    mean_abs_z = mean(abs(z)), share_abs_z_gt2 = mean(abs(z) > 2), max_rhat = max(Rhat), min_n_eff = min(n_eff)))
write.csv(bind_rows(cov, ms), file.path(sim_dir, "results", "context_b_coverage.csv"), row.names = FALSE)
write.csv(blocks, file.path(sim_dir, "results", "context_b_coverage_blocks.csv"), row.names = FALSE)
write.csv(rep_tab, file.path(sim_dir, "results", "context_b_replicates.csv"), row.names = FALSE)
logf("%d replicates (seeds %d-%d; %s lag %d; n_init %d, n_noninit %d; %d chains x %d / %d): E %d-%d, sum(y1) %d-%d; divergences total %d; max Rhat %.3f; elapsed %.0f-%.0f min",
     nrow(rep_tab), min(rep_tab$seed), max(rep_tab$seed), reps$sample[1], reps$lag[1], reps$n_init[1], reps$n_noninit[1], reps$chains[1], reps$iter[1], reps$warmup[1],
     min(rep_tab$E), max(rep_tab$E), min(rep_tab$sum_y1), max(rep_tab$sum_y1), sum(rep_tab$divergences), max(rep_tab$max_rhat), min(rep_tab$elapsed_min), max(rep_tab$elapsed_min))
for (r in seq_len(nrow(blocks))) logf("  block %-9s %2d parameters: 90%% coverage %.3f, mean |bias| %.3f, RMSE %.3f, mean |z| %.2f, |z| > 2 in %.3f", blocks$block[r], blocks$n_params[r], blocks$coverage90[r], blocks$mean_abs_bias[r], blocks$rmse[r], blocks$mean_abs_z[r], blocks$share_abs_z_gt2[r])
low <- cov |> filter(coverage90 < 0.8)
if (nrow(low)) logf("parameters with 90%% coverage below 0.80 over %d replicates: %s", length(files), paste(sprintf("%s %.2f (bias %.2f)", low$parameter, low$coverage90, low$mean_bias), collapse = "; ")) else logf("every parameter's 90%% coverage is at least 0.80")
for (r in seq_len(nrow(ms))) logf("  %s: mean over replicates %.3f (SD %.3f)", ms$parameter[r], ms$mean_over_reps[r], ms$sd_over_reps[r])
logf("written: %s", paste(file.path(sim_dir, "results", c("context_b_coverage.csv", "context_b_coverage_blocks.csv", "context_b_replicates.csv")), collapse = ", "))
