## End-to-end smoke test of 02_fit_model.R on the toy cell (audit fix 1, 2026-10-26): M1b, then M2b starting from the M1b fit.
## The audit found that 02_fit_model.R never loaded assert_fit_matches_cell(), so every M2 cell fell back to random inits with exit 0.
## This test runs the SCRIPT (not the helpers) as a subprocess in a throw-away v2 root: R/ copied, stan/ symlinked to the compiled
## programs, data/ = the toy cell of helper-toy-cell.R, output/ empty. Short chains (2 x 80), so it takes about a minute.
## Needs the compiled models (Rscript v2/R/compile_models.R); skipped otherwise.
suppressPackageStartupMessages({ library(testthat) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)

## a throw-away v2 root the fit script can treat as v2_root(): 02_fit_model.R sets V2_ROOT to the parent of its own directory
make_smoke_root <- function() {
  tmp <- tempfile("smoke_v2_"); dir.create(file.path(tmp, "v2"), recursive = TRUE)
  file.copy(file.path(root, "R"), file.path(tmp, "v2"), recursive = TRUE)
  ok <- file.symlink(file.path(root, "stan"), file.path(tmp, "v2", "stan")); if (!ok) file.copy(file.path(root, "stan"), file.path(tmp, "v2"), recursive = TRUE)
  write_toy_cell(file.path(tmp, "v2", "data"))
  dir.create(file.path(tmp, "v2", "output")); dir.create(file.path(tmp, "v2", "logs"))
  file.path(tmp, "v2")
}
run_fit <- function(v2, ...) {
  out <- suppressWarnings(system2("Rscript", c(file.path(v2, "R", "02_fit_model.R"), ...), stdout = TRUE, stderr = TRUE))
  list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out)
}
short <- c("--sample=mid", "--lag=3", "--chains=2", "--iter=80", "--warmup=40", "--max_treedepth=6")

test_that("02_fit_model.R end to end on the toy cell: M1b, then M2b from the M1b posterior means; an unusable M1 file stops the run; no M1 file falls back loudly", {
  for (m in c("m1_latent_additive", "m2_ame")) if (!file.exists(file.path(root, "stan", "build", paste0(m, ".rds")))) skip(paste("compiled model missing:", m))
  v2 <- make_smoke_root(); out_dir <- file.path(v2, "output")
  m1_file <- file.path(out_dir, "fit_mid_lag3_all_pr_M1b.rds"); m2_file <- file.path(out_dir, "fit_mid_lag3_all_pr_M2b_K2.rds")
  ## (1) M1b: Stan's random inits, recorded as such
  r1 <- run_fit(v2, short, "--model=M1b")
  expect_equal(r1$status, 0L, info = paste(tail(r1$out, 30), collapse = "\n"))
  expect_true(file.exists(m1_file))
  o1 <- readRDS(m1_file)
  expect_equal(o1$args$init_source, "stan_random"); expect_equal(o1$init$source, "stan_random"); expect_equal(o1$args$model, "M1b"); expect_false(o1$args$pilot)
  s1 <- read.csv(file.path(out_dir, "summary_mid_lag3_all_pr_M1b.csv")); expect_true(all(s1$init_source == "stan_random")); expect_true(all(is.na(s1$init_from)))
  expect_match(o1$args$stan_source_sha256, "^[0-9a-f]{64}$")
  ## (2) M2b: every chain starts at the M1b posterior means (the file is auto-located in output/)
  r2 <- run_fit(v2, short, "--model=M2b")
  expect_equal(r2$status, 0L, info = paste(tail(r2$out, 30), collapse = "\n"))
  expect_true(any(grepl("start at the posterior means of fit_mid_lag3_all_pr_M1b.rds", r2$out, fixed = TRUE)))
  expect_false(any(grepl("RANDOM inits", r2$out, fixed = TRUE)))
  expect_true(file.exists(m2_file))
  o2 <- readRDS(m2_file)
  expect_equal(o2$args$init_source, "m1_posterior_means"); expect_equal(o2$init$source, "m1_posterior_means")
  expect_equal(o2$args$init_from, "fit_mid_lag3_all_pr_M1b.rds"); expect_equal(o2$init$m1_cell, "mid_lag3_all_pr_M1b")
  expect_equal(o2$init$n_copied, 13L)                                                     # 11 shared blocks + rho_x + L_Om_joint for the b models
  expect_equal(o2$init$jitter_sd, INIT_JITTER_SD); expect_equal(o2$init$uv_range, INIT_UV_RANGE)
  s2 <- read.csv(file.path(out_dir, "summary_mid_lag3_all_pr_M2b_K2.csv"))
  expect_true(all(s2$init_source == "m1_posterior_means")); expect_true(all(s2$init_from == "fit_mid_lag3_all_pr_M1b.rds"))
  expect_equal(o2$dat$E, o1$dat$E); expect_equal(o2$dat$K, 2L); expect_equal(o2$fit@mode, 0L)
  ## (3) an M1 file that EXISTS but is unusable must stop the M2 run (non-zero exit, no M2 fit written), never fall back
  file.remove(m2_file); file.remove(file.path(out_dir, "summary_mid_lag3_all_pr_M2b_K2.csv"))
  keep <- file.path(out_dir, "fit_M1b_keep.rds"); file.copy(m1_file, keep, overwrite = TRUE)
  bad <- o1; bad$args$pilot <- TRUE; saveRDS(bad, m1_file)                               # a pilot fit under the production name
  r3 <- run_fit(v2, short, "--model=M2b")
  expect_true(r3$status != 0L, info = paste(tail(r3$out, 30), collapse = "\n"))
  expect_true(any(grepl("PILOT", r3$out, fixed = TRUE))); expect_false(file.exists(m2_file))
  bad <- o1; bad$dat$E <- bad$dat$E + 1L; saveRDS(bad, m1_file)                           # a fit of another cell (E differs)
  r3b <- run_fit(v2, short, "--model=M2b")
  expect_true(r3b$status != 0L); expect_true(any(grepl("has E = ", r3b$out, fixed = TRUE))); expect_false(file.exists(m2_file))
  writeLines("not an rds", m1_file)                                                        # an unreadable file under the production name
  r3c <- run_fit(v2, short, "--model=M2b")
  expect_true(r3c$status != 0L); expect_true(any(grepl("unreadable", r3c$out, fixed = TRUE))); expect_false(file.exists(m2_file))
  ## --init_from: a usable file under another name works, a missing one is an error
  r3d <- run_fit(v2, short, "--model=M2b", paste0("--init_from=", keep))
  expect_equal(r3d$status, 0L, info = paste(tail(r3d$out, 30), collapse = "\n"))
  expect_equal(readRDS(m2_file)$args$init_from, "fit_M1b_keep.rds"); expect_equal(readRDS(m2_file)$args$init_source, "m1_posterior_means"); file.remove(m2_file)
  expect_true(run_fit(v2, short, "--model=M2b", "--init_from=/nonexistent/fit.rds")$status != 0L)
  ## (4) NO M1 file: the fall back to random inits is allowed, loud and recorded
  file.remove(m1_file)
  r4 <- run_fit(v2, short, "--model=M2b")
  expect_equal(r4$status, 0L, info = paste(tail(r4$out, 30), collapse = "\n"))
  expect_true(any(grepl("WARNING: no matching M1 fit (fit_mid_lag3_all_pr_M1b.rds)", r4$out, fixed = TRUE)))
  expect_true(any(grepl("starts from Stan's RANDOM inits", r4$out, fixed = TRUE)))
  o4 <- readRDS(m2_file); expect_equal(o4$args$init_source, "stan_random"); expect_true(is.na(o4$args$init_from))
  expect_true(all(read.csv(file.path(out_dir, "summary_mid_lag3_all_pr_M2b_K2.csv"))$init_source == "stan_random"))
  unlink(dirname(v2), recursive = TRUE)
})
