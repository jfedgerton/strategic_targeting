## The permutation test on the contiguity + distance cell and its parallel refit mode (2026-10-29, PI request: the calibrated
## permutation test on the fmid lag-3 M1 _dist cell before the sweep). Driver argument handling and refusals, cell naming, the
## per-permutation part files and the combiner on stub data, and sequential == combined on a toy refit.
suppressPackageStartupMessages({ library(testthat); library(dplyr) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "spec.R")); source(file.path(root, "R", "helpers_data.R")); source(file.path(root, "R", "helpers_model.R")); source(file.path(root, "R", "typology.R"))
source(file.path(root, "R", "measurement_tests.R")); source(file.path(root, "sim", "simulate.R"))

test_that("cell naming: the distance cell's measurement outputs and fits carry _dist; contiguity = none is never a perm_error cell", {
  sp <- parse_spec("--distance=capital")
  expect_equal(measurement_cell_name("fmid", 3L, "m1", sp), "fmid_lag3_all_pr_M1_dist")
  expect_equal(measurement_cell_name("mid", 3L, "m1", sp), "mid_lag3_all_pr_M1_dist")
  expect_equal(measurement_cell_name("fmid", 3L, "m2", sp), "fmid_lag3_all_pr_M2_K2_dist")
  expect_equal(cell_name("fmid", 3L, "M1a", cell_spec("M1a", sp)), "fmid_lag3_all_pr_M1a_dist"); expect_equal(cell_name("fmid", 3L, "M1b", cell_spec("M1b", sp)), "fmid_lag3_all_pr_M1b_dist")
  expect_equal(cell_name("mid", 3L, "M1a", cell_spec("M1a", sp)), "mid_lag3_all_pr_M1a_dist"); expect_equal(cell_name("mid", 3L, "M1b", cell_spec("M1b", sp)), "mid_lag3_all_pr_M1b_dist")
  ## the main cell's names are unchanged
  expect_equal(measurement_cell_name("mid", 5L, "m1", parse_spec(character(0))), "mid_lag5_all_pr_M1")
  expect_equal(permutation_part_path("/x/measurement_fmid_lag3_all_pr_M1_dist_permutation_refit_parts", 7), "/x/measurement_fmid_lag3_all_pr_M1_dist_permutation_refit_parts/perm_0007.csv")
  ## the distance + major cell (2026-10-30): outputs and fits carry _dist_major, never the _dist names
  spm <- parse_spec(c("--distance=capital", "--major=joint"))
  expect_equal(measurement_cell_name("fmid", 3L, "m1", spm), "fmid_lag3_all_pr_M1_dist_major")
  expect_equal(measurement_cell_name("mid", 3L, "m1", spm), "mid_lag3_all_pr_M1_dist_major")
  expect_equal(cell_name("fmid", 3L, "M1a", cell_spec("M1a", spm)), "fmid_lag3_all_pr_M1a_dist_major"); expect_equal(cell_name("fmid", 3L, "M1b", cell_spec("M1b", spm)), "fmid_lag3_all_pr_M1b_dist_major")
  expect_equal(cell_name("mid", 3L, "M1a", cell_spec("M1a", spm)), "mid_lag3_all_pr_M1a_dist_major"); expect_equal(cell_name("mid", 3L, "M1b", cell_spec("M1b", spm)), "mid_lag3_all_pr_M1b_dist_major")
  expect_false(measurement_cell_name("fmid", 3L, "m1", spm) == measurement_cell_name("fmid", 3L, "m1", sp))
  sp_preset <- parse_spec(character(0)); for (nm in names(DIST_PRESET$variants$distance_major)) sp_preset[[nm]] <- DIST_PRESET$variants$distance_major[[nm]]
  expect_equal(cell_name("fmid", 3L, "M1b", cell_spec("M1b", spm)), cell_name("fmid", 3L, "M1b", sp_preset))   # = the dist preset's distance_major cell (the existing M1b_dist_major fits)
})

test_that("robustness driver: perm_error takes --distance=capital for the permutation test on an explicit --sample / --lag, refuses it for the other tests, keeps the other refusals, and checks the subset / combine arguments", {
  tmp <- tempfile("drv_v2_"); dir.create(file.path(tmp, "v2"), recursive = TRUE)
  file.copy(file.path(root, "R"), file.path(tmp, "v2"), recursive = TRUE); dir.create(file.path(tmp, "v2", "logs")); dir.create(file.path(tmp, "v2", "output"))
  drv <- file.path(tmp, "v2", "R", "05_robustness.R")
  run <- function(...) { out <- suppressWarnings(system2("Rscript", c(drv, "--preset=perm_error", "--family=m1", ...), stdout = TRUE, stderr = TRUE)); list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  refused <- function(r, msg) r$status != 0L && any(grepl(msg, r$out, fixed = TRUE))
  ## the throw-away root has no data directory, so a run that gets PAST the argument checks stops there: that line is the sentinel of acceptance
  passed_args <- function(r) r$status != 0L && !any(grepl("perm_error refuses", r$out, fixed = TRUE)) && any(grepl("does not exist", r$out, fixed = TRUE))
  log_names <- function() list.files(file.path(tmp, "v2", "logs"))
  DIST <- c("--sample=fmid", "--lag=3", "--distance=capital")
  ## accepted: the permutation test alone (both modes), the subset of the refit mode, the combiner; the log is named after the _dist cell
  r <- run(DIST, "--tests=permutation", "--perm_mode=rescore"); expect_true(passed_args(r), info = paste(tail(r$out, 5), collapse = "\n"))
  expect_true(any(grepl("^05_robustness_perm_error_fmid_lag3_all_pr_M1_dist_[0-9_]+\\.log$", log_names())))
  r <- run(DIST, "--tests=permutation", "--perm_mode=refit"); expect_true(passed_args(r))
  r <- run(DIST, "--tests=permutation", "--perm_mode=refit", "--perm_index=7", "--n_perm=20"); expect_true(passed_args(r), info = paste(tail(r$out, 5), collapse = "\n"))
  expect_true(any(grepl("^05_robustness_perm_error_fmid_lag3_all_pr_M1_dist_perm7_7_", log_names())))
  r <- run(DIST, "--tests=permutation", "--perm_mode=refit", "--perm_from=3", "--perm_to=5", "--n_perm=20"); expect_true(passed_args(r))
  expect_true(any(grepl("^05_robustness_perm_error_fmid_lag3_all_pr_M1_dist_perm3_5_", log_names())))
  r <- run(DIST, "--tests=permutation_combine", "--perm_mode=refit", "--n_perm=20"); expect_true(passed_args(r))
  expect_true(any(grepl("^05_robustness_perm_error_fmid_lag3_all_pr_M1_dist_combine_", log_names())))
  r <- run("--sample=mid", "--lag=3", "--distance=capital", "--tests=permutation", "--perm_mode=refit", "--perm_index=1"); expect_true(passed_args(r))
  expect_true(any(grepl("^05_robustness_perm_error_mid_lag3_all_pr_M1_dist_perm1_1_", log_names())))
  ## the main cell still runs with the default test list (the sentinel again, no refusal) and without --sample / --lag (ROBUSTNESS_MAIN_CELL)
  r <- run(); expect_true(passed_args(r)); expect_true(any(grepl("^05_robustness_perm_error_mid_lag5_all_pr_M1_[0-9_]+\\.log$", log_names())))
  ## refused: the distance with the holdout / inflation tests (the default list included)
  expect_true(refused(run(DIST), "perm_error refuses --distance=capital with --tests=holdout,permutation,inflation"))
  expect_true(refused(run(DIST, "--tests=holdout"), "perm_error refuses --distance=capital with --tests=holdout"))
  expect_true(refused(run(DIST, "--tests=permutation,inflation"), "perm_error refuses --distance=capital"))
  expect_true(passed_args(run("--distance=capital", "--tests=permutation")))   # the main cell's sample / lag with the distance: accepted (mid lag 5 _dist)
  ## refused whatever the test: contiguity = none, interact, winsor, pairs != all
  expect_true(refused(run(DIST, "--contiguity=none", "--tests=permutation"), "perm_error refuses --contiguity=none"))
  expect_true(refused(run(DIST, "--interact=prefixed", "--tests=permutation"), "perm_error refuses --interact=prefixed"))
  ## --major=joint (2026-10-30): accepted with the distance for the permutation test (both modes, the subset, the combiner); the log is named
  ## after the _dist_major cell; refused for the other tests and without the distance
  MAJ <- c(DIST, "--major=joint")
  r <- run(MAJ, "--tests=permutation", "--perm_mode=rescore"); expect_true(passed_args(r), info = paste(tail(r$out, 5), collapse = "\n"))
  expect_true(any(grepl("^05_robustness_perm_error_fmid_lag3_all_pr_M1_dist_major_[0-9_]+\\.log$", log_names())))
  r <- run(MAJ, "--tests=permutation", "--perm_mode=refit"); expect_true(passed_args(r))
  r <- run(MAJ, "--tests=permutation", "--perm_mode=refit", "--perm_index=4", "--n_perm=20", "--threads_per_chain=8"); expect_true(passed_args(r), info = paste(tail(r$out, 5), collapse = "\n"))
  expect_true(any(grepl("^05_robustness_perm_error_fmid_lag3_all_pr_M1_dist_major_perm4_4_", log_names())))
  r <- run(MAJ, "--tests=permutation_combine", "--perm_mode=refit", "--n_perm=20", "--threads_per_chain=8"); expect_true(passed_args(r))
  expect_true(any(grepl("^05_robustness_perm_error_fmid_lag3_all_pr_M1_dist_major_combine_", log_names())))
  r <- run("--sample=mid", "--lag=3", "--distance=capital", "--major=joint", "--tests=permutation", "--perm_mode=refit", "--perm_index=2"); expect_true(passed_args(r))
  expect_true(any(grepl("^05_robustness_perm_error_mid_lag3_all_pr_M1_dist_major_perm2_2_", log_names())))
  expect_true(refused(run(MAJ), "perm_error refuses --major=joint with --tests=holdout,permutation,inflation"))
  expect_true(refused(run(MAJ, "--tests=holdout"), "perm_error refuses --major=joint with --tests=holdout"))
  expect_true(refused(run(MAJ, "--tests=inflation"), "perm_error refuses --major=joint with --tests=inflation"))
  expect_true(refused(run(MAJ, "--tests=permutation,inflation"), "perm_error refuses --major=joint with --tests=permutation,inflation"))
  expect_true(refused(run("--sample=fmid", "--lag=3", "--major=joint", "--tests=permutation"), "perm_error refuses --major=joint without --distance=capital"))
  expect_true(refused(run("--sample=fmid", "--lag=3", "--major=joint", "--tests=permutation_combine", "--perm_mode=refit"), "perm_error refuses --major=joint without --distance=capital"))
  expect_true(refused(run(MAJ, "--contiguity=none", "--tests=permutation"), "perm_error refuses --contiguity=none"))
  expect_true(refused(run(MAJ, "--interact=prefixed", "--tests=permutation"), "perm_error refuses --interact=prefixed"))
  ## --threads_per_chain (2026-10-30): a positive integer; > 1 only with --perm_mode=refit and the permutation test or its combiner; 1 everywhere
  expect_true(passed_args(run(DIST, "--tests=permutation", "--perm_mode=refit", "--perm_index=1", "--threads_per_chain=1")))
  expect_true(passed_args(run(DIST, "--tests=permutation", "--perm_mode=refit", "--threads_per_chain=4")))                       # the sequential refit loop
  expect_true(passed_args(run("--tests=holdout", "--threads_per_chain=1")))                                                        # 1 is the default everywhere
  expect_true(refused(run(DIST, "--tests=permutation", "--perm_mode=rescore", "--threads_per_chain=8"), "--threads_per_chain=8 applies to the permutation refits only"))
  expect_true(refused(run("--tests=holdout", "--threads_per_chain=8"), "--threads_per_chain=8 applies to the permutation refits only"))
  expect_true(refused(run("--tests=inflation", "--perm_mode=refit", "--threads_per_chain=2"), "--threads_per_chain=2 applies to the permutation refits only"))
  expect_true(refused(run(DIST, "--tests=permutation", "--perm_mode=refit", "--threads_per_chain=0"), "must be a positive integer"))
  expect_true(refused(run(DIST, "--tests=permutation", "--perm_mode=refit", "--threads_per_chain=2.5"), "must be a positive integer"))
  expect_true(refused(run(DIST, "--tests=permutation", "--perm_mode=refit", "--threads_per_chain=-1"), "must be a positive integer"))
  expect_true(refused(run(DIST, "--tests=permutation", "--perm_mode=refit", "--threads_per_chain"), "--threads_per_chain needs a value"))
  expect_true(run(DIST, "--tests=permutation", "--perm_mode=refit", "--threads_per_chain=eight")$status != 0L)
  expect_true(refused(run(DIST, "--changes=winsor", "--tests=permutation"), "perm_error refuses --changes=winsor"))
  expect_true(refused(run(DIST, "--pairs=originators", "--tests=permutation"), "perm_error refuses --pairs=originators"))
  ## the subset arguments: refit mode and the permutation test alone, integers, from <= to, to <= n_perm, one form, no bare flag
  expect_true(refused(run(DIST, "--tests=permutation", "--perm_mode=rescore", "--perm_index=1"), "need --tests=permutation --perm_mode=refit alone"))
  expect_true(refused(run(DIST, "--tests=permutation,inflation", "--perm_mode=refit", "--perm_index=1"), "need --tests=permutation --perm_mode=refit alone"))
  expect_true(refused(run(DIST, "--tests=permutation", "--perm_mode=refit", "--perm_from=5", "--perm_to=3"), "must be integers with 1 <= from <= to"))
  expect_true(refused(run(DIST, "--tests=permutation", "--perm_mode=refit", "--perm_index=0"), "must be integers with 1 <= from <= to"))
  expect_true(refused(run(DIST, "--tests=permutation", "--perm_mode=refit", "--perm_index=2.5"), "must be integers with 1 <= from <= to"))
  expect_true(refused(run(DIST, "--tests=permutation", "--perm_mode=refit", "--perm_index=21", "--n_perm=20"), "exceeds --n_perm=20"))
  expect_true(refused(run(DIST, "--tests=permutation", "--perm_mode=refit", "--perm_from=1"), "--perm_from and --perm_to go together"))
  expect_true(refused(run(DIST, "--tests=permutation", "--perm_mode=refit", "--perm_index=1", "--perm_from=1", "--perm_to=1"), "pass one form, not both"))
  expect_true(refused(run(DIST, "--tests=permutation", "--perm_mode=refit", "--perm_index"), "--perm_index needs a value"))
  ## the combiner: alone and in refit mode
  expect_true(refused(run(DIST, "--tests=permutation_combine,permutation", "--perm_mode=refit"), "--tests=permutation_combine goes alone"))
  expect_true(refused(run(DIST, "--tests=permutation_combine", "--perm_mode=rescore"), "--tests=permutation_combine goes alone, with --perm_mode=refit"))
  expect_true(run(DIST, "--tests=permutation_combine", "--perm_mode=refit", "--perm_index=1")$status != 0L)
  expect_true(run(DIST, "--tests=nonsense")$status != 0L)
  unlink(tmp, recursive = TRUE)
})

test_that("part files: written atomically with the settings, read back, and the combiner refuses a missing, duplicated or extra index, another setting, a wrong seed, a file / row index mismatch", {
  settings <- data.frame(cell = "fmid_lag3_all_pr_M1_dist", family = "M1", ame_K = NA_integer_, sample = "fmid", lag = 3L, spec = "distance=capital", chains = 4L, iter = 1500L, warmup = 750L,
                         v2_seed = 123L, threads_per_chain = 8L, grainsize = choose_grainsize(14597L, 8L), E = 1055L, N1 = 14597L, sum_y1 = 695L,
                         fit_a = "fit_fmid_lag3_all_pr_M1a_dist.rds", fit_b = "fit_fmid_lag3_all_pr_M1b_dist.rds", statistic = "stage2_loo_elpd_of_refits", stringsAsFactors = FALSE)
  expect_equal(settings$grainsize, 457L)
  one <- function(p, a, b) list(perm = p, seed = V2_SEED + p, a = a, b = b, adv = b - a)
  pd <- file.path(tempfile("parts_"), "measurement_fmid_lag3_all_pr_M1_dist_permutation_refit_parts")
  expect_error(read_permutation_parts(pd, settings), "no parts directory")
  paths <- c(write_permutation_part(pd, one(1, -100, -98), settings), write_permutation_part(pd, one(3, -101, -99.5), settings), write_permutation_part(pd, one(2, -100.5, -100.2), settings))
  expect_equal(basename(paths), c("perm_0001.csv", "perm_0003.csv", "perm_0002.csv")); expect_false(any(grepl("\\.tmp$", list.files(pd))))
  parts <- read_permutation_parts(pd, settings)
  expect_equal(parts$perm, 1:3); expect_equal(parts$seed, 123L + 1:3); expect_equal(parts$a, c(-100, -100.5, -101)); expect_equal(parts$b, c(-98, -100.2, -99.5))
  expect_true(all(names(settings) %in% names(parts)))
  ## the combiner on the three parts: in index order, the sequential formula
  pt <- permutation_combine(parts, obs_a = -100, obs_b = -97, statistic = "stage2_loo_elpd_of_refits", n_perm = 3)
  expect_equal(pt$permuted_a, c(-100, -100.5, -101)); expect_equal(pt$permuted_b, c(-98, -100.2, -99.5)); expect_equal(pt$permuted_advantage, c(2, 0.3, 1.5))
  expect_equal(pt$observed_advantage_b_minus_a, 3); expect_equal(pt$summary$n_exceedances, 0L); expect_equal(pt$p_value, 1 / 4); expect_equal(pt$summary$n_perm, 3)
  expect_equal(pt$summary$statistic, "stage2_loo_elpd_of_refits"); expect_match(pt$summary$interpretation, "calibrated permutation p-value")
  pt2 <- permutation_combine(parts, obs_a = -100, obs_b = -99, statistic = "s", n_perm = 3)       # observed advantage 1: 2 and 1.5 exceed
  expect_equal(pt2$summary$n_exceedances, 2L); expect_equal(pt2$p_value, 3 / 4)
  ## missing index (n_perm larger than the parts), extra index (n_perm smaller), duplicate
  expect_error(permutation_combine(parts, -100, -97, "s", n_perm = 4), "missing permutation index 4")
  expect_error(permutation_combine(parts, -100, -97, "s", n_perm = 2), "index outside 1..2: 3")
  expect_error(permutation_combine(rbind(parts, parts[1, ]), -100, -97, "s", n_perm = 3), "duplicated permutation index 1")
  expect_error(permutation_combine(parts[-2, ], -100, -97, "s", n_perm = 3), "missing permutation index 2")
  ## another setting in one part: refused by file name and column
  s2 <- settings; s2$iter <- 1000L
  write_permutation_part(pd, one(4, -100, -99), s2)
  expect_error(read_permutation_parts(pd, settings), "perm_0004.csv were made with iter = 1000, this run has 1500")
  file.remove(permutation_part_path(pd, 4))
  s3 <- settings; s3$cell <- "mid_lag3_all_pr_M1_dist"; s3$fit_b <- "fit_mid_lag3_all_pr_M1b_dist.rds"
  write_permutation_part(pd, one(4, -100, -99), s3)
  expect_error(read_permutation_parts(pd, settings), "perm_0004.csv were made with cell = mid_lag3_all_pr_M1_dist")
  file.remove(permutation_part_path(pd, 4))
  s4 <- settings; s4$ame_K <- 2L
  write_permutation_part(pd, one(4, -100, -99), s4)
  expect_error(read_permutation_parts(pd, settings), "perm_0004.csv were made with ame_K = 2")
  file.remove(permutation_part_path(pd, 4))
  ## another threads / grainsize setting (2026-10-30): a serial part (threads 1, grainsize N1) is refused by an 8-thread combine run, and the
  ## reverse; a part whose grainsize alone differs (another N1 -> grainsize pairing) is refused as well
  s5 <- settings; s5$threads_per_chain <- 1L; s5$grainsize <- 14597L
  write_permutation_part(pd, one(4, -100, -99), s5)
  expect_error(read_permutation_parts(pd, settings), "perm_0004.csv were made with threads_per_chain = 1, this run has 8")
  expect_error(read_permutation_parts(pd, s5), "perm_0001.csv, perm_0002.csv, perm_0003.csv were made with threads_per_chain = 8, this run has 1")   # and the reverse
  file.remove(permutation_part_path(pd, 4))
  s6 <- settings; s6$grainsize <- 456L
  write_permutation_part(pd, one(4, -100, -99), s6)
  expect_error(read_permutation_parts(pd, settings), "perm_0004.csv were made with grainsize = 456, this run has 457")
  file.remove(permutation_part_path(pd, 4))
  s7 <- settings; s7$cell <- "fmid_lag3_all_pr_M1_dist_major"; s7$spec <- "distance=capital major=joint"; s7$fit_a <- "fit_fmid_lag3_all_pr_M1a_dist_major.rds"; s7$fit_b <- "fit_fmid_lag3_all_pr_M1b_dist_major.rds"
  write_permutation_part(pd, one(4, -100, -99), s7)
  expect_error(read_permutation_parts(pd, settings), "perm_0004.csv were made with cell = fmid_lag3_all_pr_M1_dist_major")   # a _dist_major part never combines into the _dist test
  file.remove(permutation_part_path(pd, 4))
  ## a wrong seed is refused at write time and at read time; a file whose index differs from its row is refused
  expect_error(write_permutation_part(pd, list(perm = 5, seed = 999, a = -1, b = -1, adv = 0), settings))
  d <- read.csv(permutation_part_path(pd, 3), stringsAsFactors = FALSE); d$seed <- 999L; write.csv(d, permutation_part_path(pd, 3), row.names = FALSE)
  expect_error(read_permutation_parts(pd, settings), "perm_0003.csv were not run with seed V2_SEED \\+ perm")
  d$seed <- 123L + 3L; d$perm <- 9L; write.csv(d, permutation_part_path(pd, 3), row.names = FALSE)
  expect_error(read_permutation_parts(pd, settings), "perm_0003.csv hold another permutation index")
  d$perm <- 3L; d$a <- NA; write.csv(d, permutation_part_path(pd, 3), row.names = FALSE)
  expect_error(read_permutation_parts(pd, settings), "non-finite statistic in part\\(s\\) perm_0003.csv")
  write_permutation_part(pd, one(3, -101, -99.5), settings); expect_equal(read_permutation_parts(pd, settings)$perm, 1:3)
  ## a non-finite statistic is refused at write time; a stray file that is not perm_<pppp>.csv is ignored
  expect_error(write_permutation_part(pd, list(perm = 6, seed = 129, a = Inf, b = -1, adv = -Inf), settings))
  writeLines("junk", file.path(pd, "perm_7.csv")); expect_equal(read_permutation_parts(pd, settings)$perm, 1:3)
  unlink(dirname(pd), recursive = TRUE)
})

test_that("sequential refit mode == per-permutation refits + combiner on a toy cell (same seeds, same sampler settings, same p-value); rescore mode unchanged", {
  p1 <- file.path(root, "stan", "build", "m1_latent_additive.rds")
  if (!file.exists(p1)) skip("compiled model missing")
  suppressPackageStartupMessages(library(rstan))
  sm <- list(m1 = readRDS(p1))
  sim <- simulate_targeting("M1b", S = 7L, Tt = 4L, seed = 123L, target_rate = 0.6)
  Wc <- c("W1", "W2"); Xc <- paste0("X", 1:4); Zc <- c("Z1", "Z2"); xr <- rep(0, 4)
  da <- build_stan_data(sim$sy, sim$alts, sim$ev, Wc, Xc, Zc, joint = 0L, x_ref = xr, allow_orphan_initiations = TRUE)
  db <- build_stan_data(sim$sy, sim$alts, sim$ev, Wc, Xc, Zc, joint = 1L, x_ref = xr, anchor = da$anchor, allow_orphan_initiations = TRUE)
  fa <- suppressWarnings(sampling(sm$m1, data = da, chains = 1, iter = 40, warmup = 20, seed = 123, refresh = 0, pars = EXCL_PARS, include = FALSE))
  fb <- suppressWarnings(sampling(sm$m1, data = db, chains = 1, iter = 40, warmup = 20, seed = 123, refresh = 0, pars = EXCL_PARS, include = FALSE))
  n_perm <- 3L; ch <- 1; it <- 40; wu <- 20
  seq_pt <- suppressWarnings(permutation_test(sim$sy, sim$alts, sim$ev, Wc, Xc, Zc, fa, fb, n_perm = n_perm, mode = "refit", family = "m1", stan_models = sm, chains = ch, iter = it, warmup = wu, x_ref = xr))
  expect_equal(length(seq_pt$permuted_a), n_perm); expect_true(all(is.finite(seq_pt$permuted_advantage))); expect_equal(seq_pt$summary$statistic, "stage2_loo_elpd_of_refits")
  expect_equal(seq_pt$p_value, (1 + seq_pt$summary$n_exceedances) / (1 + n_perm))
  ## the parts, one permutation at a time, in a scrambled order, through the part files
  obs <- suppressWarnings(permutation_observed(sim$sy, sim$alts, sim$ev, Wc, Xc, Zc, fa, fb, mode = "refit", family = "m1", x_ref = xr))
  expect_equal(obs$obs_a, seq_pt$observed_a); expect_equal(obs$obs_b, seq_pt$observed_b)
  settings <- data.frame(cell = "toy", family = "M1", ame_K = NA_integer_, chains = ch, iter = it, warmup = wu, E = da$E, N1 = da$N1, stringsAsFactors = FALSE)
  pd <- tempfile("toy_parts_")
  for (p in c(2L, 3L, 1L)) {
    one <- suppressWarnings(permutation_refit_one(p, sim$sy, sim$alts, sim$ev, Wc, Xc, Zc, obs, family = "m1", stan_models = sm, chains = ch, iter = it, warmup = wu))
    expect_equal(one$seed, V2_SEED + p); write_permutation_part(pd, one, settings)
  }
  parts <- read_permutation_parts(pd, settings)
  comb <- permutation_combine(parts, obs$obs_a, obs$obs_b, obs$statistic, n_perm = n_perm)
  ## identical statistics permutation by permutation (same shuffle seed, same sampling seed and settings) and the same p-value and summary
  expect_equal(comb$permuted_a, seq_pt$permuted_a, tolerance = 1e-8); expect_equal(comb$permuted_b, seq_pt$permuted_b, tolerance = 1e-8)
  expect_equal(comb$permuted_advantage, seq_pt$permuted_advantage, tolerance = 1e-8)
  expect_equal(comb$p_value, seq_pt$p_value); expect_equal(comb$summary$n_exceedances, seq_pt$summary$n_exceedances)
  expect_equal(comb$summary[, setdiff(names(comb$summary), c("mean_permuted_advantage", "mean_permuted_a", "mean_permuted_b"))],
               seq_pt$summary[, setdiff(names(seq_pt$summary), c("mean_permuted_advantage", "mean_permuted_a", "mean_permuted_b"))])
  expect_equal(comb$summary$mean_permuted_advantage, seq_pt$summary$mean_permuted_advantage, tolerance = 1e-8)
  ## the serial parts record threads 1 and grainsize N1 (the previous behaviour)
  expect_equal(one$threads_per_chain, 1L); expect_equal(one$grainsize, da$N1); expect_equal(rstan_options("threads_per_chain"), 1L)
  ## within-chain threads (2026-10-30): a refit with threads = 2 sets rstan's threads_per_chain and slices the state-years
  ## (grainsize = choose_grainsize(N1, 2) < N1), records both, is reproducible like with like (threads 2 against threads 2, the same seeds),
  ## and its parts combine only with parts of the same threads / grainsize setting
  gs2 <- choose_grainsize(da$N1, 2L); expect_true(gs2 < da$N1)
  one2 <- suppressWarnings(permutation_refit_one(1L, sim$sy, sim$alts, sim$ev, Wc, Xc, Zc, obs, family = "m1", stan_models = sm, chains = ch, iter = it, warmup = wu, threads = 2L))
  expect_equal(one2$threads_per_chain, 2L); expect_equal(one2$grainsize, gs2); expect_equal(one2$perm, 1L); expect_equal(one2$seed, V2_SEED + 1L)
  expect_true(is.finite(one2$a)); expect_true(is.finite(one2$b)); expect_equal(one2$adv, one2$b - one2$a); expect_equal(rstan_options("threads_per_chain"), 2L)
  one2b <- suppressWarnings(permutation_refit_one(1L, sim$sy, sim$alts, sim$ev, Wc, Xc, Zc, obs, family = "m1", stan_models = sm, chains = ch, iter = it, warmup = wu, threads = 2L))
  expect_equal(one2b$a, one2$a, tolerance = 1e-8); expect_equal(one2b$b, one2$b, tolerance = 1e-8)
  settings1 <- cbind(settings, threads_per_chain = 1L, grainsize = da$N1); settings2 <- cbind(settings, threads_per_chain = 2L, grainsize = gs2)
  pd2 <- tempfile("toy_parts2_"); write_permutation_part(pd2, one2, settings2)
  expect_error(read_permutation_parts(pd2, settings1), "threads_per_chain = 2, this run has 1")
  expect_equal(read_permutation_parts(pd2, settings2)$perm, 1L)
  one_s <- one; one_s$perm <- 2L; one_s$seed <- V2_SEED + 2L                  # a serial part beside it, under the next index
  write_permutation_part(pd2, one_s, settings1)
  expect_error(read_permutation_parts(pd2, settings2), "perm_0002.csv were made with threads_per_chain = 1, this run has 2")
  ## threads outside the refit mode, a non-integer or zero threads: refused
  expect_error(permutation_refit_one(1L, sim$sy, sim$alts, sim$ev, Wc, Xc, Zc, obs, family = "m1", stan_models = sm, chains = ch, iter = it, warmup = wu, threads = 0L))
  expect_error(permutation_refit_one(1L, sim$sy, sim$alts, sim$ev, Wc, Xc, Zc, obs, family = "m1", stan_models = sm, chains = ch, iter = it, warmup = wu, threads = 1.5))
  rstan_options(threads_per_chain = 1)
  ## a permutation outside the refit mode is refused; the rescore mode still runs and is labelled as the diagnostic with p NA; threads there are refused
  obs_r <- suppressWarnings(permutation_observed(sim$sy, sim$alts, sim$ev, Wc, Xc, Zc, fa, fb, mode = "rescore", family = "m1", n_draws = 5, x_ref = xr))
  expect_error(permutation_refit_one(1L, sim$sy, sim$alts, sim$ev, Wc, Xc, Zc, obs_r, family = "m1", stan_models = sm))
  rs <- suppressWarnings(permutation_test(sim$sy, sim$alts, sim$ev, Wc, Xc, Zc, fa, fb, n_perm = 2, mode = "rescore", family = "m1", n_draws = 5, x_ref = xr))
  expect_true(is.na(rs$p_value)); expect_match(rs$summary$interpretation, "RELIANCE DIAGNOSTIC"); expect_equal(rs$observed_a, obs_r$obs_a); expect_equal(length(rs$permuted_a), 2)
  expect_error(permutation_test(sim$sy, sim$alts, sim$ev, Wc, Xc, Zc, fa, fb, n_perm = 2, mode = "rescore", family = "m1", n_draws = 5, x_ref = xr, threads = 2L), "threads applies to the refit mode only")
  unlink(pd, recursive = TRUE); unlink(pd2, recursive = TRUE)
})

test_that("08_perm_refit.sbatch dispatches the three modes with the 1..N_PERM guard and the threads flag; submit_perm_refit.sh submits the rescore job (open), the refit array 1-N_PERM on basic with chains x threads cores plus the combine job afterok, or the combine job alone", {
  sb <- readLines(file.path(root, "slurm", "08_perm_refit.sbatch"))
  expect_true(any(grepl("^  rescore\\)", sb))); expect_true(any(grepl("^  refit\\)", sb))); expect_true(any(grepl("^  combine\\)", sb)))
  expect_true(any(grepl("--tests=permutation --perm_mode=rescore ;;$", sb)))                                              # the rescore run takes no threads flag (it samples nothing)
  expect_true(any(grepl("--tests=permutation --perm_mode=refit --perm_index=\\$p --threads_per_chain=\\$THR", sb)))   # the refit tasks and the combiner carry the SAME threads
  expect_true(any(grepl("--tests=permutation_combine --perm_mode=refit --threads_per_chain=\\$THR", sb))); expect_true(any(grepl("outside 1\\.\\.\\$NP", sb)))
  expect_true(any(grepl("^THR=\\$\\{THREADS_PER_CHAIN:-8\\}", sb)))
  expect_true(any(grepl("^#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32", sb))); expect_true(any(grepl("^#SBATCH --mem=128G", sb))); expect_true(any(grepl("^#SBATCH --time=2-00:00:00", sb)))
  expect_true(any(grepl("^#SBATCH --account=jfe4_cr_default", sb))); expect_true(any(grepl("^#SBATCH --partition=basic", sb)))
  ## the task body through a stub Rscript: the flags the driver receives in every mode, with the default and an explicit THREADS_PER_CHAIN, and
  ## a two-switch SWITCHES (distance + major) reaching the driver as two arguments
  tmp <- tempfile("sub_"); dir.create(file.path(tmp, "v2", "R"), recursive = TRUE); dir.create(file.path(tmp, "bin"))
  file.copy(file.path(root, "slurm"), file.path(tmp, "v2"), recursive = TRUE)
  writeLines("", file.path(tmp, "v2", "R", "01_build_data.R"))                   # 00_env.sh checks for it
  writeLines(c("#!/bin/bash", "echo \"MODULE $*\" >/dev/null"), file.path(tmp, "bin", "module")); Sys.chmod(file.path(tmp, "bin", "module"), "0755")
  env_stub <- file.path(tmp, "v2", "slurm", "00_env.sh")                        # the module / gcc lines of the real 00_env.sh cannot run here; keep V2_ROOT only
  writeLines(c("export V2_ROOT=\"${V2_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}/v2}\"", "mkdir -p \"$V2_ROOT/logs\" \"$V2_ROOT/output\""), env_stub)
  rec_r <- file.path(tmp, "rscript_calls.txt")
  writeLines(c("#!/bin/bash", sprintf("echo \"$*\" >> %s", rec_r)), file.path(tmp, "bin", "Rscript")); Sys.chmod(file.path(tmp, "bin", "Rscript"), "0755")
  run_body <- function(env) { out <- suppressWarnings(system2("bash", "v2/slurm/08_perm_refit.sbatch", stdout = TRUE, stderr = TRUE,
                                                                env = c(sprintf("PATH=%s:%s", file.path(tmp, "bin"), Sys.getenv("PATH")), sprintf("SLURM_SUBMIT_DIR=%s", tmp), env)))
    list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  old <- setwd(tmp); on.exit(setwd(old), add = TRUE)
  r <- run_body(c("MODE=refit", "SLURM_ARRAY_TASK_ID=3")); expect_equal(r$status, 0L, info = paste(r$out, collapse = "\n"))
  calls <- readLines(rec_r); expect_length(calls, 1L)
  expect_equal(calls[1], "v2/R/05_robustness.R --preset=perm_error --family=m1 --sample=fmid --lag=3 --r0=0.5 --chains=4 --iter=1500 --warmup=750 --n_perm=20 --distance=capital --tests=permutation --perm_mode=refit --perm_index=3 --threads_per_chain=8")
  r <- run_body(c("MODE=refit", "SLURM_ARRAY_TASK_ID=1", "THREADS_PER_CHAIN=4", "SWITCHES='--distance=capital --major=joint'", "SAMPLE=mid")); expect_equal(r$status, 0L, info = paste(r$out, collapse = "\n"))
  calls <- readLines(rec_r); expect_length(calls, 2L)
  expect_equal(calls[2], "v2/R/05_robustness.R --preset=perm_error --family=m1 --sample=mid --lag=3 --r0=0.5 --chains=4 --iter=1500 --warmup=750 --n_perm=20 --distance=capital --major=joint --tests=permutation --perm_mode=refit --perm_index=1 --threads_per_chain=4")
  r <- run_body(c("MODE=combine", "THREADS_PER_CHAIN=4", "SWITCHES='--distance=capital --major=joint'")); expect_equal(r$status, 0L)
  calls <- readLines(rec_r); expect_length(calls, 3L)
  expect_equal(calls[3], "v2/R/05_robustness.R --preset=perm_error --family=m1 --sample=fmid --lag=3 --r0=0.5 --chains=4 --iter=1500 --warmup=750 --n_perm=20 --distance=capital --major=joint --tests=permutation_combine --perm_mode=refit --threads_per_chain=4")
  r <- run_body(c("MODE=rescore", "THREADS_PER_CHAIN=8")); expect_equal(r$status, 0L)
  calls <- readLines(rec_r); expect_length(calls, 4L)
  expect_equal(calls[4], "v2/R/05_robustness.R --preset=perm_error --family=m1 --sample=fmid --lag=3 --r0=0.5 --chains=4 --iter=1500 --warmup=750 --n_perm=20 --distance=capital --tests=permutation --perm_mode=rescore")
  expect_true(run_body(c("MODE=refit", "SLURM_ARRAY_TASK_ID=21"))$status != 0L); expect_true(run_body(c("MODE=refit", "SLURM_ARRAY_TASK_ID=0"))$status != 0L)
  expect_true(run_body(c("MODE=refit", "SLURM_ARRAY_TASK_ID=1", "THREADS_PER_CHAIN=0"))$status != 0L); expect_true(run_body(c("MODE=refit", "SLURM_ARRAY_TASK_ID=1", "THREADS_PER_CHAIN=x"))$status != 0L)
  expect_length(readLines(rec_r), 4L)
  ## the submit script against a stub sbatch
  rec <- file.path(tmp, "sbatch_calls.txt")
  writeLines(c("#!/bin/bash", sprintf("echo \"MODE='${MODE-<unset>}' SAMPLE='${SAMPLE-<unset>}' LAG='${LAG-<unset>}' FAMILY='${FAMILY-<unset>}' SWITCHES='${SWITCHES-<unset>}' N_PERM='${N_PERM-<unset>}' CHAINS='${CHAINS-<unset>}' THREADS_PER_CHAIN='${THREADS_PER_CHAIN-<unset>}' args: $*\" >> %s", rec),
               sprintf("echo \"$(( 4000 + $(wc -l < %s) ));cluster\"", rec)), file.path(tmp, "bin", "sbatch"))
  Sys.chmod(file.path(tmp, "bin", "sbatch"), "0755")
  run <- function(env, ...) { out <- suppressWarnings(system2("bash", c("v2/slurm/submit_perm_refit.sh", ...), stdout = TRUE, stderr = TRUE,
                                                              env = c(sprintf("PATH=%s:%s", file.path(tmp, "bin"), Sys.getenv("PATH")), env)))
    list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  ## rescore: one open-account job on 4 cores, the defaults (fmid, 3, m1, --distance=capital, 20; threads 8 exported but not used by the rescore run)
  r <- run("MODE=rescore"); expect_equal(r$status, 0L, info = paste(r$out, collapse = "\n"))
  calls <- readLines(rec); expect_length(calls, 1L)
  expect_match(calls[1], "^MODE='rescore' SAMPLE='fmid' LAG='3' FAMILY='m1' SWITCHES='--distance=capital' N_PERM='20' CHAINS='4' THREADS_PER_CHAIN='8' args: --parsable --export=ALL --account=open --time=24:00:00 --cpus-per-task=4 --mem=32G v2/slurm/08_perm_refit.sbatch$")
  expect_true(any(grepl("perm rescore submitted: job 4001", r$out)))
  ## refit: the array 1-N_PERM with chains x threads cores (no account override: the sbatch's jfe4_cr_default / basic), then the combine job afterok on open
  r <- run(c("MODE=refit", "N_PERM=20"), "--time=1-00:00:00"); expect_equal(r$status, 0L, info = paste(r$out, collapse = "\n"))
  calls <- readLines(rec); expect_length(calls, 3L)
  expect_match(calls[2], "^MODE='refit' SAMPLE='fmid' LAG='3' FAMILY='m1' SWITCHES='--distance=capital' N_PERM='20' CHAINS='4' THREADS_PER_CHAIN='8' args: --parsable --export=ALL --array=1-20 --cpus-per-task=32 --mem=128G --time=1-00:00:00 v2/slurm/08_perm_refit.sbatch$")
  expect_match(calls[3], "^MODE='combine' SAMPLE='fmid' LAG='3' FAMILY='m1' SWITCHES='--distance=capital' N_PERM='20' CHAINS='4' THREADS_PER_CHAIN='8' args: --parsable --export=ALL --account=open --time=04:00:00 --cpus-per-task=4 --mem=32G --dependency=afterok:4002 --kill-on-invalid-dep=yes --time=1-00:00:00 v2/slurm/08_perm_refit.sbatch$")
  expect_true(any(grepl("perm refit array submitted: job 4002 \\(indices 1-20, one refit of a and b per task, 4 chains x 8 threads = 32 cores, 128G", r$out))); expect_true(any(grepl("perm combine submitted: job 4003 \\(afterok:4002", r$out)))
  ## mid, a throttle, per-job options, N_PERM = 5, 4 threads per chain (16 cores, 64 GB), the distance + major switches (2026-10-30)
  r <- run(c("MODE=refit", "SAMPLE=mid", "N_PERM=5", "THROTTLE=2", "THREADS_PER_CHAIN=4", "SWITCHES='--distance=capital --major=joint'", "ARRAY_OPTS=--time=1-00:00:00", "COMBINE_OPTS=--mem=16G")); expect_equal(r$status, 0L, info = paste(r$out, collapse = "\n"))
  calls <- readLines(rec); expect_length(calls, 5L)
  expect_match(calls[4], "^MODE='refit' SAMPLE='mid' LAG='3' FAMILY='m1' SWITCHES='--distance=capital --major=joint' N_PERM='5' CHAINS='4' THREADS_PER_CHAIN='4' args: --parsable --export=ALL --array=1-5%2 --cpus-per-task=16 --mem=64G --time=1-00:00:00 v2/slurm/08_perm_refit.sbatch$")
  expect_match(calls[5], "^MODE='combine' SAMPLE='mid' .* SWITCHES='--distance=capital --major=joint' .* THREADS_PER_CHAIN='4' args: .* --dependency=afterok:4004 --kill-on-invalid-dep=yes --mem=16G v2/slurm/08_perm_refit.sbatch$")
  expect_true(any(grepl("4 chains x 4 threads = 16 cores, 64G", r$out)))
  ## combine alone; the main-spec cell through SWITCHES=""; 2 chains x 8 threads = 16 cores for a later refit array
  r <- run(c("MODE=combine", "SWITCHES=", "SAMPLE=mid", "LAG=5")); expect_equal(r$status, 0L)
  calls <- readLines(rec); expect_length(calls, 6L)
  expect_match(calls[6], "^MODE='combine' SAMPLE='mid' LAG='5' FAMILY='m1' SWITCHES='' N_PERM='20' CHAINS='4' THREADS_PER_CHAIN='8' args: --parsable --export=ALL --account=open --time=04:00:00 --cpus-per-task=4 --mem=32G v2/slurm/08_perm_refit.sbatch$")
  r <- run(c("MODE=refit", "CHAINS=2", "N_PERM=3")); expect_equal(r$status, 0L)
  calls <- readLines(rec); expect_length(calls, 8L)
  expect_match(calls[7], "^MODE='refit' .* N_PERM='3' CHAINS='2' THREADS_PER_CHAIN='8' args: --parsable --export=ALL --array=1-3 --cpus-per-task=16 --mem=64G v2/slurm/08_perm_refit.sbatch$")
  ## refusals: --array, a bad N_PERM, a bad THREADS_PER_CHAIN / CHAINS, no MODE, an unknown MODE
  expect_true(run("MODE=refit", "--array=1-3")$status != 0L); expect_true(run(c("MODE=refit", "N_PERM=0"))$status != 0L)
  expect_true(run(c("MODE=refit", "N_PERM=x"))$status != 0L); expect_true(run("MODE=nonsense")$status != 0L); expect_true(run(character(0))$status != 0L)
  expect_true(run(c("MODE=refit", "THREADS_PER_CHAIN=0"))$status != 0L); expect_true(run(c("MODE=refit", "THREADS_PER_CHAIN=2.5"))$status != 0L); expect_true(run(c("MODE=refit", "CHAINS=0"))$status != 0L)
  expect_length(readLines(rec), 8L)
  unlink(tmp, recursive = TRUE)
})
