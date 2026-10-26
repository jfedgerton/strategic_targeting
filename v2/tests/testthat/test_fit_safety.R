## Fit safety and provenance (review 2026-10-25, C1 / C2 / C3 / C4) and the summary labels (B2)
suppressPackageStartupMessages({ library(testthat) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "spec.R")); source(file.path(root, "R", "helpers_data.R")); source(file.path(root, "R", "helpers_model.R"))

test_that("C1: pilot fits are refused, by the saved flag or by the pilot settings; C2: stale fits need allow_stale", {
  prod <- list(cell = "x", args = list(pilot = FALSE, chains = 4L, iter = 2000L, warmup = 1000L))
  expect_true(assert_production_fit(prod))
  expect_error(assert_production_fit(list(cell = "p", args = list(pilot = TRUE, chains = 1L, iter = 60L))), "PILOT")
  expect_error(assert_production_fit(list(cell = "old_pilot", args = list(chains = 1L, iter = 60L))), "PILOT")          # a pre-flag fit with the pilot settings
  expect_true(assert_production_fit(list(cell = "old_prod", args = list(chains = 4L, iter = 2000L))))                  # a pre-flag production fit passes
  expect_error(assert_production_fit(list(cell = "noargs")), "no saved sampling configuration")
  st <- list(cell = "s", args = prod$args, stale = list(reason = "changed", changed_files = "events.rds"))
  expect_error(assert_production_fit(st), "STALE")
  expect_true(assert_production_fit(st, allow_stale = TRUE))
})

test_that("C2: an unreadable fit file is treated as missing with a log line", {
  td <- tempfile(); dir.create(td)
  writeLines("not an rds", file.path(td, "fit_bad.rds"))
  expect_null(read_fit_or_null(file.path(td, "fit_bad.rds")))
  saveRDS(list(cell = "no_fit"), file.path(td, "fit_partial.rds")); expect_null(read_fit_or_null(file.path(td, "fit_partial.rds")))
  saveRDS(list(cell = "ok", fit = 1, dat = list()), file.path(td, "fit_ok.rds")); expect_equal(read_fit_or_null(file.path(td, "fit_ok.rds"))$cell, "ok")
})

test_that("C3: the compiled program must be the checked-out source; the sidecar hash must match", {
  p <- file.path(root, "stan", "build", "m1_latent_additive.rds")
  if (!file.exists(p)) skip("compiled model missing")
  suppressPackageStartupMessages(library(rstan))
  sm <- readRDS(p); src <- file.path(root, "stan", "m1_latent_additive.stan"); side <- file.path(root, "stan", "build", "m1_latent_additive.sha256")
  prov <- assert_compiled_current(sm, src, side)
  expect_match(prov$stan_source_sha256, "^[0-9a-f]{64}$"); expect_equal(prov$stan_source_sha256, stan_source_sha256(src))
  expect_equal(trimws(readLines(side)[1]), prov$stan_source_sha256)
  ## an edited source (one changed prior) is detected
  td <- tempfile(); dir.create(td); alt <- file.path(td, "m1_latent_additive.stan")
  writeLines(sub("alpha0 ~ normal(0, 2.5);", "alpha0 ~ normal(0, 5);", readLines(src), fixed = TRUE), alt)
  expect_error(assert_compiled_current(sm, alt), "stale binary")
  ## trailing whitespace / line endings do not matter
  writeLines(paste0(readLines(src), "  "), alt); expect_equal(assert_compiled_current(sm, alt)$compiled_model_code_sha256, prov$compiled_model_code_sha256)
  ## a wrong sidecar is detected
  writeLines("deadbeef", file.path(td, "x.sha256")); expect_error(assert_compiled_current(sm, src, file.path(td, "x.sha256")), "sidecar")
  expect_error(assert_compiled_current(sm, src, file.path(td, "missing.sha256")), "missing")
})

test_that("B2: summary labels map parameters to axes and covariates in every D_mode", {
  W <- c("dM_i", "dR_i", "M_lvl", "R_lvl", "ongoing_tm1"); X <- c("dM_i", "dR_i"); Z <- c("contig_land", "joint_dem", "dM_j_h")
  l3 <- label_parameters(c("alpha0", "alpha[3]", "tau", "mu0[2]", "beta[2,1]", "gamma[3]", "sd_re[1]", "sd_re[3]", "Sigma[1,2]", "Omega[2,3]", "R_Om[1,3]", "sigma_x[1]", "u_diag[1]"), W, X, Z, 3L)
  expect_equal(l3$covariate[l3$parameter == "alpha[3]"], "M_lvl"); expect_equal(l3$axis[l3$parameter == "mu0[2]"], "relational")
  expect_equal(l3$covariate[l3$parameter == "beta[2,1]"], "dR_i"); expect_equal(l3$axis[l3$parameter == "beta[2,1]"], "material")
  expect_equal(l3$covariate[l3$parameter == "gamma[3]"], "dM_j_h"); expect_equal(l3$covariate[l3$parameter == "sd_re[3]"], "a_relational")
  expect_equal(l3$covariate[l3$parameter == "Omega[2,3]"], "a_material x a_relational"); expect_equal(l3$axis[l3$parameter == "Sigma[1,2]"], "material x relational")
  expect_equal(l3$covariate[l3$parameter == "u_diag[1]"], "latent_dim_1"); expect_equal(l3$covariate[l3$parameter == "tau"], "IV")
  l2 <- label_parameters(c("mu0[1]", "beta[1,1]", "sd_re[2]", "Sigma[1,1]"), W, X, Z, 2L)
  expect_equal(l2$axis[1], "relational"); expect_equal(l2$covariate[2], "dM_i"); expect_equal(l2$covariate[3], "a_relational"); expect_equal(l2$axis[4], "relational x relational")
  expect_equal(axis_labels(1L), "material"); expect_equal(axis_labels(3L), c("material", "relational"))
  expect_error(label_parameters("alpha[9]", W, X, Z, 3L))                                  # an index beyond W is an error, not a silent NA
})

test_that("C4: a fit must score exactly the assembled cell", {
  source(file.path(root, "R", "typology.R")); source(file.path(root, "R", "measurement_tests.R"))
  dat <- list(E = 3L, N1 = 5L, S = 4L, ev_sy = c(1L, 2L, 2L), ev_row = c(1L, 4L, 5L), alt_target = c(2L, 3L, 4L, 1L, 3L), ev_w = c(1, 0.5, 0.5))
  obj <- list(cell = "c", dat = dat, W_cols = "w", X_cols = "x", Z_cols = "z")
  expect_true(assert_fit_matches_cell(obj, dat, "w", "x", "z"))
  d2 <- dat; d2$E <- 2L; expect_error(assert_fit_matches_cell(obj, d2, "w", "x", "z"), "E = 3")
  d3 <- dat; d3$alt_target[5] <- 2L; expect_error(assert_fit_matches_cell(obj, d3, "w", "x", "z"), "events")
  expect_error(assert_fit_matches_cell(obj, dat, "w", "x", c("z", "z2")), "Z columns")
  d4 <- dat; d4$ev_w <- c(1, 1, 1); expect_error(assert_fit_matches_cell(obj, d4, "w", "x", "z"), "weights")
})
