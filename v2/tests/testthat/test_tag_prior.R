## Run tag and prior-sensitivity switches (2026-10-30; PI-approved checks on the jointness result Sigma[1,2] / rho of M1b).
## (1) --tag=<alnum>: _t<tag> is appended LAST to every name built from the spec; the default leaves every name byte-identical; tagged
##     and untagged names never collide; the M2 init lookup takes the tagged M1 fit only; 03_extract.R refuses a renamed fit;
##     05_robustness.R refuses the tag. (2) --rho_prior / --sigma_prior: the suffixes, the refusals (M2, M1c, rho with a non-joint M1
##     model), the program routing (m1_latent_prior only off the defaults), the hash guard on both programs, the two Stan sources differ
##     in exactly the hyperparameter lines, and (with the compiled programs) the default hyperparameters reproduce the main program's log
##     density and draws while flat / tight move the log density by exactly the intended term. (3) 02_fit_model.R end to end on the toy
##     cell with a tag, a tagged M2, and a prior switch (short chains; skipped without the compiled programs).
suppressPackageStartupMessages({ library(testthat) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "spec.R")); source(file.path(root, "R", "helpers_data.R")); source(file.path(root, "R", "helpers_model.R"))

main <- parse_spec(character(0))

test_that("tag: parsing, the _t<tag> suffix last, byte-identical names without it, no collision, the tag is not a switch", {
  expect_equal(main$tag, ""); expect_equal(spec_tag(main), ""); expect_equal(SPEC_DEFAULTS$tag, ""); expect_equal(TAG_PREFIX, "_t")
  expect_equal(spec_tag(list(ww2 = "drop")), "")                                           # a spec saved before the tag existed
  t1 <- parse_spec("--tag=long"); expect_equal(t1$tag, "long"); expect_equal(spec_tag(t1), "long")
  expect_equal(parse_spec("--tag=A1b2")$tag, "A1b2")
  for (bad in c("--tag=a-b", "--tag=a_b", "--tag=a.b", "--tag=a b", "--tag=lo/ng")) expect_error(parse_spec(bad), "letters and digits only", info = bad)
  expect_error(parse_spec("--tag="), "empty value is refused"); expect_error(parse_spec("--tag"), "--tag needs a value", fixed = TRUE)
  ## the suffix goes last: after every data switch, after the prior switches, before _pilot / _STALE (which 02_fit_model.R appends to the cell)
  expect_equal(cell_name("fmid", 3L, "M1b", t1), "fmid_lag3_all_pr_M1b_tlong")
  expect_equal(cell_name("fmid", 3L, "M2b", t1), "fmid_lag3_all_pr_M2b_K2_tlong")
  expect_equal(cell_name("fmid", 3L, "M1c", t1), "fmid_lag3_all_pr_M1c_tlong")
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec(c("--distance=capital", "--tag=long"))), "fmid_lag3_all_pr_M1b_dist_tlong")
  expect_equal(cell_name("mid", 3L, "M1b", parse_spec(c("--distance=capital", "--major=joint", "--tag=long"))), "mid_lag3_all_pr_M1b_dist_major_tlong")
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec(c("--rho_prior=flat", "--sigma_prior=ln05", "--tag=pri"))), "fmid_lag3_all_pr_M1b_rhoflat_sigln05_tpri")
  all_on <- parse_spec(c("--pairs=originators", "--context=coalition", "--coal_ties=strict", "--coal_drop=top3", "--ww2=drop", "--controls=full",
                         "--distance=capital", "--contiguity=none", "--major=joint", "--interact=prefixed", "--weights=none", "--changes=winsor", "--direction=sideA",
                         "--rho_prior=flat", "--sigma_prior=ln05", "--tag=x9"))
  expect_equal(spec_suffix(all_on), "_orig_co_strict_drop3_noww2_fullctrl_dist_nocontig_major_ix_w1_winsor_sideA_rhoflat_sigln05_tx9")
  expect_equal(compare_group_name("fmid", 3L, t1), "fmid_lag3_all_pr_tlong"); expect_equal(measurement_cell_name("fmid", 3L, "m1", t1), "fmid_lag3_all_pr_M1_tlong")
  ## without a tag every name is the documented one, byte for byte, for the main spec and every preset variant and model
  specs <- c(list(main = character(0)), lapply(ROBUSTNESS_PRESETS[c("noww2", "noweights", "fullcontrols", "winsor")], function(p) sprintf("--%s=%s", names(p), unlist(p))),
             list(dist = "--distance=capital", dist_major = c("--distance=capital", "--major=joint"), nocontig = c("--distance=capital", "--contiguity=none"), co = "--context=coalition"))
  untagged <- character(0); tagged <- character(0)
  for (nm in names(specs)) for (m in MODELS) {
    a <- cell_name("fmid", 3L, m, parse_spec(specs[[nm]])); b <- cell_name("fmid", 3L, m, parse_spec(c(specs[[nm]], "--tag=long"))); c2 <- cell_name("fmid", 3L, m, parse_spec(c(specs[[nm]], "--tag=long2")))
    expect_equal(b, paste0(a, "_tlong"), info = paste(nm, m)); expect_equal(c2, paste0(a, "_tlong2"), info = paste(nm, m))
    untagged <- c(untagged, a); tagged <- c(tagged, b, c2)
  }
  expect_equal(cell_name("fmid", 3L, "M1b", main), "fmid_lag3_all_pr_M1b"); expect_equal(cell_name("fmid", 3L, "M1b", parse_spec("--distance=capital")), "fmid_lag3_all_pr_M1b_dist")
  expect_length(intersect(untagged, tagged), 0L); expect_equal(anyDuplicated(c(untagged, tagged)), 0L)
  ## no switch suffix can be read as a tag: none starts with _t followed by letters / digits only
  expect_false(any(grepl("^_t[A-Za-z0-9]+$", unlist(SPEC_SUFFIX))))
  ## the tag is not a switch: a tagged main-spec run is still the main specification (02_fit_model.R logs it as such and the tag separately)
  expect_true(spec_is_main(t1)); expect_false(spec_is_main(parse_spec(c("--tag=long", "--ww2=drop"))))
  expect_true(grepl("tag=long", spec_string(t1)))
  ## the M2 init lookup takes the TAGGED M1 fit and never the untagged one (and the reverse)
  td <- tempfile("tag_"); dir.create(td)
  saveRDS(list(cell = "x"), file.path(td, "fit_mid_lag3_all_pr_M1b.rds"))
  expect_equal(m1_cell_for_m2("mid", 3L, "M2b", parse_spec("--tag=long")), "mid_lag3_all_pr_M1b_tlong")
  expect_true(is.na(locate_m1_fit("mid", 3L, "M2b", parse_spec("--tag=long"), td)))                     # only the untagged file exists: NA, not that file
  expect_equal(basename(locate_m1_fit("mid", 3L, "M2b", main, td)), "fit_mid_lag3_all_pr_M1b.rds")
  saveRDS(list(cell = "y"), file.path(td, "fit_mid_lag3_all_pr_M1b_tlong.rds")); file.remove(file.path(td, "fit_mid_lag3_all_pr_M1b.rds"))
  expect_equal(basename(locate_m1_fit("mid", 3L, "M2b", parse_spec("--tag=long"), td)), "fit_mid_lag3_all_pr_M1b_tlong.rds")
  expect_true(is.na(locate_m1_fit("mid", 3L, "M2b", main, td)))                                          # only the tagged file exists: the untagged lookup finds nothing
  unlink(td, recursive = TRUE)
})

test_that("prior switches: suffixes, refusals, program routing, hyperparameters, and old specs count as the main priors", {
  expect_equal(SPEC_CHOICES$rho_prior, c("beta25", "flat")); expect_equal(SPEC_CHOICES$sigma_prior, c("ln1", "ln05"))
  expect_equal(main$rho_prior, "beta25"); expect_equal(main$sigma_prior, "ln1"); expect_equal(PRIOR_SWITCHES, c("rho_prior", "sigma_prior"))
  expect_equal(tail(names(SPEC_SUFFIX), 2), c("rho_prior", "sigma_prior"))                                # after every data switch
  expect_equal(unname(SPEC_SUFFIX$rho_prior), "_rhoflat"); expect_equal(unname(SPEC_SUFFIX$sigma_prior), "_sigln05")
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec("--rho_prior=flat")), "fmid_lag3_all_pr_M1b_rhoflat")
  expect_equal(cell_name("mid", 3L, "M1b", parse_spec("--sigma_prior=ln05")), "mid_lag3_all_pr_M1b_sigln05")
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec(c("--rho_prior=flat", "--sigma_prior=ln05"))), "fmid_lag3_all_pr_M1b_rhoflat_sigln05")
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec(c("--distance=capital", "--rho_prior=flat"))), "fmid_lag3_all_pr_M1b_dist_rhoflat")
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec(c("--distance=capital", "--sigma_prior=ln05"))), "fmid_lag3_all_pr_M1b_dist_sigln05")
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec(c("--distance=capital", "--major=joint", "--sigma_prior=ln05"))), "fmid_lag3_all_pr_M1b_dist_major_sigln05")
  expect_equal(cell_name("fmid", 3L, "M1a", parse_spec("--sigma_prior=ln05")), "fmid_lag3_all_pr_M1a_sigln05")              # sigma_x exists in every M1 model
  expect_error(parse_spec("--rho_prior=uniform"), "--rho_prior must be one of beta25 | flat", fixed = TRUE)
  expect_error(parse_spec("--sigma_prior=0.5"), "--sigma_prior must be one of ln1 | ln05", fixed = TRUE)
  expect_error(parse_spec("--rho_prior"), "--rho_prior needs a value", fixed = TRUE); expect_error(parse_spec("--sigma_prior"), "--sigma_prior needs a value", fixed = TRUE)
  ## refused with the other programs (M2, M1c) and the rho prior with a non-joint M1 model
  for (m in c("M2m", "M2r", "M2a", "M2b", "M1c")) { expect_error(cell_name("fmid", 3L, m, parse_spec("--rho_prior=flat")), "M1 m / r / a / b programs only", info = m)
                                                   expect_error(cell_name("fmid", 3L, m, parse_spec("--sigma_prior=ln05")), "M1 m / r / a / b programs only", info = m) }
  for (m in c("M1m", "M1r", "M1a")) expect_error(cell_name("fmid", 3L, m, parse_spec("--rho_prior=flat")), "no Sigma correlation rho", info = m)
  expect_error(cell_name("fmid", 3L, "M1a", parse_spec(c("--rho_prior=flat", "--sigma_prior=ln05"))), "no Sigma correlation rho")
  ## program routing: the main program for the defaults, the prior copy off them; M2 / M1c never route there
  expect_equal(stan_program("M1b"), "m1_latent_additive"); expect_equal(stan_program("M1b", main), "m1_latent_additive"); expect_equal(stan_program("M1b", NULL), "m1_latent_additive")
  expect_equal(stan_program("M1b", parse_spec("--rho_prior=flat")), "m1_latent_prior"); expect_equal(stan_program("M1a", parse_spec("--sigma_prior=ln05")), "m1_latent_prior")
  expect_equal(stan_program("M1b", parse_spec("--tag=long")), "m1_latent_additive")                     # a tag changes no program
  expect_equal(stan_program("M2b", parse_spec("--rho_prior=flat")), "m2_ame"); expect_equal(stan_program("M1c", parse_spec("--rho_prior=flat")), "m1_context")
  expect_true("m1_latent_prior" %in% STAN_PROGRAMS); expect_true(file.exists(file.path(root, "stan", "m1_latent_prior.stan")))
  ## the hyperparameters: the defaults are the main program's literals
  expect_equal(prior_hyper(main), list(prior_rho_shape = 2.5, prior_sigma_x_sd = 1))
  expect_equal(prior_hyper(parse_spec("--rho_prior=flat")), list(prior_rho_shape = 1, prior_sigma_x_sd = 1))
  expect_equal(prior_hyper(parse_spec("--sigma_prior=ln05")), list(prior_rho_shape = 2.5, prior_sigma_x_sd = 0.5))
  expect_equal(prior_hyper(parse_spec(c("--rho_prior=flat", "--sigma_prior=ln05"))), list(prior_rho_shape = 1, prior_sigma_x_sd = 0.5))
  expect_equal(prior_hyper(NULL), list(prior_rho_shape = 2.5, prior_sigma_x_sd = 1)); expect_equal(prior_hyper(list(ww2 = "drop")), list(prior_rho_shape = 2.5, prior_sigma_x_sd = 1))
  ## a spec saved before 2026-10-30 (no rho_prior / sigma_prior / tag fields) is the main priors and gets its old suffix
  old <- main[setdiff(names(main), c("rho_prior", "sigma_prior", "tag"))]; old$ww2 <- "drop"
  expect_true(spec_prior_is_main(old)); expect_true(spec_prior_is_main(NULL)); expect_true(spec_prior_is_main(main)); expect_false(spec_prior_is_main(parse_spec("--sigma_prior=ln05")))
  expect_equal(spec_suffix(old), "_noww2"); expect_equal(stan_program("M1b", old), "m1_latent_additive")
})

test_that("prior switches: the hash guard checks a prior-switch fit against the prior program's sidecar, never the main one's, and a fit that claims another program is refused", {
  td <- tempfile("guard_"); dir.create(file.path(td, "stan", "build"), recursive = TRUE)
  for (m in c("m1_latent_additive", "m1_latent_prior")) {
    file.copy(file.path(root, "stan", paste0(m, ".stan")), file.path(td, "stan"))
    writeLines(stan_source_sha256(file.path(td, "stan", paste0(m, ".stan"))), file.path(td, "stan", "build", paste0(m, ".sha256")))
  }
  sha_main <- stan_source_sha256(file.path(td, "stan", "m1_latent_additive.stan")); sha_prior <- stan_source_sha256(file.path(td, "stan", "m1_latent_prior.stan"))
  expect_false(sha_main == sha_prior)
  flat <- parse_spec("--rho_prior=flat")
  ok_main <- list(cell = "m", args = list(model = "M1b", stan_source_sha256 = sha_main, stan_program = "m1_latent_additive"), spec = main)
  ok_old  <- list(cell = "o", args = list(model = "M1b", stan_source_sha256 = sha_main))                            # a pre-2026-10-30 fit: no spec fields, no stan_program
  ok_prior <- list(cell = "p", args = list(model = "M1b", stan_source_sha256 = sha_prior, stan_program = "m1_latent_prior"), spec = flat)
  expect_true(assert_fit_stan_current(ok_main, root = td)); expect_true(assert_fit_stan_current(ok_old, root = td)); expect_true(assert_fit_stan_current(ok_prior, root = td))
  ## a prior-switch fit stamped with the MAIN program's hash (or the reverse) is a fit of another program
  expect_error(assert_fit_stan_current(list(cell = "x", args = list(model = "M1b", stan_source_sha256 = sha_main, stan_program = "m1_latent_prior"), spec = flat), root = td), "different programs")
  expect_error(assert_fit_stan_current(list(cell = "y", args = list(model = "M1b", stan_source_sha256 = sha_prior, stan_program = "m1_latent_additive"), spec = main), root = td), "different programs")
  ## the recorded program must agree with what the model code and the switches imply
  expect_error(assert_fit_stan_current(list(cell = "z", args = list(model = "M1b", stan_source_sha256 = sha_prior, stan_program = "m1_latent_additive"), spec = flat), root = td), "not what its name says")
  expect_error(assert_fit_stan_current(list(cell = "w", args = list(model = "M1b", stan_source_sha256 = sha_main, stan_program = "m1_latent_prior"), spec = main), root = td), "not what its name says")
  ## editing the prior source after compilation is caught for prior fits and leaves main fits untouched, and the reverse
  src_p <- file.path(td, "stan", "m1_latent_prior.stan"); keep <- readLines(src_p)
  writeLines(c(keep, "// edited"), src_p); expect_error(assert_fit_stan_current(ok_prior, root = td), "changed after compilation"); expect_true(assert_fit_stan_current(ok_main, root = td))
  writeLines(keep, src_p); expect_true(assert_fit_stan_current(ok_prior, root = td))
  unlink(td, recursive = TRUE)
})

test_that("the prior program is the main program with the two hyperparameters as data and nothing else", {
  main_src <- readLines(file.path(root, "stan", "m1_latent_additive.stan")); prior_src <- readLines(file.path(root, "stan", "m1_latent_prior.stan"))
  ## (a) the main program is byte-identical to the one the sidecar was compiled from, when the sidecar is there (the existing fits stay valid)
  side <- file.path(root, "stan", "build", "m1_latent_additive.sha256")
  if (file.exists(side)) expect_equal(trimws(readLines(side)[1]), stan_source_sha256(file.path(root, "stan", "m1_latent_additive.stan")))
  ## (b) substitute the literals back and drop the added lines: the result is the main source line for line
  back <- prior_src
  i_sig <- which(back == "    sigma_x ~ lognormal(0, prior_sigma_x_sd);"); expect_length(i_sig, 1L); back[i_sig] <- "    sigma_x ~ lognormal(0, 1);"
  i_rho <- which(back == "      target += beta_lpdf((rho_x[1] + 1) / 2 | prior_rho_shape, prior_rho_shape);"); expect_length(i_rho, 1L); back[i_rho] <- "      target += beta_lpdf((rho_x[1] + 1) / 2 | 2.5, 2.5);"
  i_dat <- which(grepl("^  real<lower=0> prior_rho_shape;", back) | grepl("^  real<lower=0> prior_sigma_x_sd;", back)); expect_length(i_dat, 2L); expect_equal(diff(i_dat), 1L)
  back <- back[-i_dat]
  ## the header: the prior program's lines 2..10 replace the main program's line 2 (the file name line)
  expect_equal(main_src[1], prior_src[1]); expect_true(grepl("m1_latent_additive.stan", main_src[2], fixed = TRUE)); expect_true(grepl("m1_latent_prior.stan", prior_src[2], fixed = TRUE))
  hdr_end <- which(prior_src == "//")[2]       # the blank comment line that closes the inserted header block (= the main program's line 3)
  expect_true(all(grepl("^//", prior_src[2:hdr_end]))); expect_equal(main_src[3], "//")
  back <- c(back[1], main_src[2], back[hdr_end:length(back)])
  expect_identical(back, main_src)
  ## (c) the hyperparameters are used exactly where the literals were, and nowhere else
  expect_equal(sum(grepl("prior_rho_shape", prior_src)), 3L); expect_equal(sum(grepl("prior_sigma_x_sd", prior_src)), 3L)   # header, data, use
  expect_false(any(grepl("prior_rho_shape|prior_sigma_x_sd", main_src)))
  expect_true(any(grepl("beta_lpdf((r_c_a1 + 1) / 2 | 2.5, 2.5)", prior_src, fixed = TRUE)))        # the D = 1 corr(c, a) prior is NOT switched
})

test_that("compiled programs: the default hyperparameters reproduce the main program's log density and draws; flat / tight move the log density by exactly the intended term", {
  p_main <- file.path(root, "stan", "build", "m1_latent_additive.rds"); p_prior <- file.path(root, "stan", "build", "m1_latent_prior.rds")
  if (!file.exists(p_main) || !file.exists(p_prior)) skip("compiled model missing (m1_latent_additive / m1_latent_prior)")
  suppressPackageStartupMessages(library(rstan)); source(file.path(root, "sim", "simulate.R"))
  sm_main <- readRDS(p_main); sm_prior <- readRDS(p_prior)
  assert_compiled_current(sm_prior, file.path(root, "stan", "m1_latent_prior.stan"), file.path(root, "stan", "build", "m1_latent_prior.sha256"))
  sim <- simulate_targeting("M1b", S = 6L, Tt = 3L, seed = 123L, target_rate = 0.6)
  dat <- build_stan_data(sim$sy, sim$alts, sim$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"), joint = 1L, x_ref = rep(0, 4))
  expect_null(dat$prior_rho_shape); expect_null(dat$prior_sigma_x_sd)                          # build_stan_data() adds nothing: the main data are unchanged
  with_prior <- function(d, spec) { h <- prior_hyper(spec); d$prior_rho_shape <- h$prior_rho_shape; d$prior_sigma_x_sd <- h$prior_sigma_x_sd; d }
  d_def <- with_prior(dat, main); d_flat <- with_prior(dat, parse_spec("--rho_prior=flat")); d_tight <- with_prior(dat, parse_spec("--sigma_prior=ln05"))
  ## log densities at two random parameter points, through the empty stanfit objects (chains = 0)
  f_main <- sampling(sm_main, data = dat, chains = 0); f_def <- sampling(sm_prior, data = d_def, chains = 0)
  f_flat <- sampling(sm_prior, data = d_flat, chains = 0); f_tight <- sampling(sm_prior, data = d_tight, chains = 0)
  n_up <- get_num_upars(f_main); expect_equal(get_num_upars(f_def), n_up)
  set.seed(1); u1 <- rnorm(n_up, 0, 0.5); u2 <- rnorm(n_up, 0, 0.5)
  lp <- function(f, u) log_prob(f, u, adjust_transform = TRUE, gradient = FALSE)
  c1 <- constrain_pars(f_def, u1); c2 <- constrain_pars(f_def, u2)
  ## (i) defaults == the main program, term for term (same propto conventions: the literal 1 and 2.5 are data in both)
  expect_equal(lp(f_def, u1), lp(f_main, u1), tolerance = 1e-12); expect_equal(lp(f_def, u2), lp(f_main, u2), tolerance = 1e-12)
  g_main <- log_prob(f_main, u1, adjust_transform = TRUE, gradient = TRUE); g_def <- log_prob(f_def, u1, adjust_transform = TRUE, gradient = TRUE)
  expect_equal(attr(g_def, "gradient"), attr(g_main, "gradient"), tolerance = 1e-10)
  ## (ii) flat rho: the log density loses exactly the Beta(2.5, 2.5) term on (rho + 1) / 2 (Beta(1, 1) contributes 0); nothing else moves
  for (k in 1:2) { u <- list(u1, u2)[[k]]; cp <- list(c1, c2)[[k]]
    expect_equal(lp(f_flat, u) - lp(f_def, u), -dbeta((cp$rho_x[1] + 1) / 2, 2.5, 2.5, log = TRUE), tolerance = 1e-10) }
  expect_gt(abs(lp(f_flat, u1) - lp(f_def, u1)), 1e-3)
  ## (iii) tight sigma: lognormal(0, 0.5) - lognormal(0, 1) = -1.5 (log sigma_x)^2 per element plus a constant that cancels between two points
  tight_term <- function(cp) -1.5 * sum(log(cp$sigma_x)^2)
  expect_equal((lp(f_tight, u1) - lp(f_def, u1)) - (lp(f_tight, u2) - lp(f_def, u2)), tight_term(c1) - tight_term(c2), tolerance = 1e-10)
  expect_gt(abs(tight_term(c1) - tight_term(c2)), 1e-3)
  ## (iv) same seed, same program up to the data values: the draws of the prior program at the defaults are the main program's draws
  s_main <- suppressWarnings(sampling(sm_main, data = dat, chains = 1, iter = 120, warmup = 60, seed = 123, refresh = 0, pars = c("mu0", "sigma_x", "rho", "gamma", "tau")))
  s_def  <- suppressWarnings(sampling(sm_prior, data = d_def, chains = 1, iter = 120, warmup = 60, seed = 123, refresh = 0, pars = c("mu0", "sigma_x", "rho", "gamma", "tau")))
  expect_equal(as.array(s_def), as.array(s_main), tolerance = 1e-12)
  s_flat <- suppressWarnings(sampling(sm_prior, data = d_flat, chains = 1, iter = 120, warmup = 60, seed = 123, refresh = 0, pars = c("mu0", "sigma_x", "rho", "gamma", "tau")))
  expect_false(isTRUE(all.equal(as.array(s_flat), as.array(s_main), tolerance = 1e-8)))
})

## ---- 02_fit_model.R / 03_extract.R / 05_robustness.R end to end on the toy cell (short chains) -----------------------------------------------------
make_root <- function() {
  tmp <- tempfile("tagprior_v2_"); dir.create(file.path(tmp, "v2"), recursive = TRUE)
  file.copy(file.path(root, "R"), file.path(tmp, "v2"), recursive = TRUE)
  ok <- file.symlink(file.path(root, "stan"), file.path(tmp, "v2", "stan")); if (!ok) file.copy(file.path(root, "stan"), file.path(tmp, "v2"), recursive = TRUE)
  write_toy_cell(file.path(tmp, "v2", "data"))
  dir.create(file.path(tmp, "v2", "output")); dir.create(file.path(tmp, "v2", "logs"))
  file.path(tmp, "v2")
}
run_script <- function(v2, script, ...) {
  out <- suppressWarnings(system2("Rscript", c(file.path(v2, "R", script), ...), stdout = TRUE, stderr = TRUE))
  list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out)
}
short <- c("--sample=mid", "--lag=3", "--chains=2", "--iter=80", "--warmup=40", "--max_treedepth=6")

test_that("02_fit_model.R: a tagged run writes only tagged files, a tagged M2 starts from the tagged M1 only, 03_extract.R takes the tagged fit and refuses a renamed one", {
  for (m in c("m1_latent_additive", "m2_ame")) if (!file.exists(file.path(root, "stan", "build", paste0(m, ".rds")))) skip(paste("compiled model missing:", m))
  v2 <- make_root(); out_dir <- file.path(v2, "output"); log_dir <- file.path(v2, "logs")
  files <- function(pat) list.files(out_dir, pattern = pat)
  ## dry run with a tag: standata carries the tag
  d <- run_script(v2, "02_fit_model.R", short, "--model=M1b", "--tag=long", "--dry_run=TRUE"); expect_equal(d$status, 0L, info = paste(tail(d$out, 20), collapse = "\n"))
  expect_true(file.exists(file.path(out_dir, "standata_mid_lag3_all_pr_M1b_tlong.rds"))); expect_false(file.exists(file.path(out_dir, "standata_mid_lag3_all_pr_M1b.rds")))
  sd <- readRDS(file.path(out_dir, "standata_mid_lag3_all_pr_M1b_tlong.rds")); expect_equal(sd$cell, "mid_lag3_all_pr_M1b_tlong"); expect_equal(sd$tag, "long"); expect_equal(sd$spec$tag, "long"); expect_equal(sd$stan_program, "m1_latent_additive")
  ## the untagged fit first, then the tagged one with other sampler settings: both exist afterwards, the untagged file untouched
  r0 <- run_script(v2, "02_fit_model.R", short, "--model=M1b"); expect_equal(r0$status, 0L, info = paste(tail(r0$out, 20), collapse = "\n"))
  f0 <- file.path(out_dir, "fit_mid_lag3_all_pr_M1b.rds"); expect_true(file.exists(f0)); md5_0 <- tools::md5sum(f0)
  r1 <- run_script(v2, "02_fit_model.R", "--sample=mid", "--lag=3", "--chains=2", "--iter=100", "--warmup=50", "--max_treedepth=6", "--model=M1b", "--tag=long")
  expect_equal(r1$status, 0L, info = paste(tail(r1$out, 20), collapse = "\n"))
  f1 <- file.path(out_dir, "fit_mid_lag3_all_pr_M1b_tlong.rds"); expect_true(file.exists(f1)); expect_equal(tools::md5sum(f0), md5_0)
  expect_true(any(grepl("run tag long", r1$out, fixed = TRUE))); expect_true(any(grepl("cell = mid_lag3_all_pr_M1b_tlong; spec: ", r1$out, fixed = TRUE)))
  expect_true(any(grepl("(MAIN specification)", r1$out, fixed = TRUE)))                                   # the tag is not a switch
  expect_true(file.exists(file.path(out_dir, "summary_mid_lag3_all_pr_M1b_tlong.csv"))); expect_true(file.exists(file.path(out_dir, "summary_mid_lag3_all_pr_M1b.csv")))
  expect_length(list.files(log_dir, pattern = "^02_fit_mid_lag3_all_pr_M1b_tlong_"), 2L); expect_length(list.files(log_dir, pattern = "^02_fit_mid_lag3_all_pr_M1b_[0-9]"), 1L)   # tagged: the dry run and the fit
  o1 <- readRDS(f1); expect_equal(o1$cell, "mid_lag3_all_pr_M1b_tlong"); expect_equal(o1$args$tag, "long"); expect_equal(o1$spec$tag, "long"); expect_equal(o1$args$iter, 100L)
  expect_equal(o1$args$stan_program, "m1_latent_additive"); expect_equal(o1$prior$rho_prior, "beta25"); expect_equal(o1$prior$prior_rho_shape, 2.5); expect_equal(o1$prior$prior_sigma_x_sd, 1)
  expect_null(o1$dat$prior_rho_shape)                                                                      # the main program's data: no hyperparameter field
  o0 <- readRDS(f0); expect_equal(o0$args$tag, ""); expect_equal(o0$spec$tag, ""); expect_null(o0$dat$prior_rho_shape); expect_identical(o0$dat, o1$dat)
  s1 <- read.csv(file.path(out_dir, "summary_mid_lag3_all_pr_M1b_tlong.csv")); expect_true(all(s1$tag == "long")); expect_true(all(s1$rho_prior == "beta25")); expect_true(all(s1$sigma_prior == "ln1")); expect_true(all(s1$stan_program == "m1_latent_additive"))
  s0 <- read.csv(file.path(out_dir, "summary_mid_lag3_all_pr_M1b.csv")); expect_true(all(is.na(s0$tag) | s0$tag == ""))
  ## 03_extract.R on the tagged fit: every output carries the tag; the untagged outputs are not written
  e1 <- run_script(v2, "03_extract.R", paste0("--fit=", f1)); expect_equal(e1$status, 0L, info = paste(tail(e1$out, 20), collapse = "\n"))
  expect_true(file.exists(file.path(out_dir, "typology_mid_lag3_all_pr_M1b_tlong.csv"))); expect_false(file.exists(file.path(out_dir, "typology_mid_lag3_all_pr_M1b.csv")))
  expect_true(any(grepl("run tag long: this is the TAGGED fit", e1$out, fixed = TRUE)))
  ## a renamed fit is refused: the tagged fit under the untagged name, and the untagged fit under a tagged name
  ren <- file.path(out_dir, "fit_fmid_lag3_all_pr_M1b.rds"); file.copy(f1, ren)
  e2 <- run_script(v2, "03_extract.R", paste0("--fit=", ren)); expect_true(e2$status != 0L); expect_true(any(grepl("a renamed fit is refused", e2$out, fixed = TRUE))); file.remove(ren)
  ren2 <- file.path(out_dir, "fit_mid_lag3_all_pr_M1b_tother.rds"); file.copy(f0, ren2)
  e3 <- run_script(v2, "03_extract.R", paste0("--fit=", ren2)); expect_true(e3$status != 0L); expect_true(any(grepl("must be named fit_mid_lag3_all_pr_M1b.rds", e3$out, fixed = TRUE))); file.remove(ren2)
  expect_false(file.exists(file.path(out_dir, "typology_mid_lag3_all_pr_M1b.csv")))
  ## a tagged M2 starts from the TAGGED M1 fit; with a tag that has no M1 fit it falls back loudly and never takes the untagged or another tag's file
  r2 <- run_script(v2, "02_fit_model.R", short, "--model=M2b", "--tag=long"); expect_equal(r2$status, 0L, info = paste(tail(r2$out, 20), collapse = "\n"))
  expect_true(any(grepl("start at the posterior means of fit_mid_lag3_all_pr_M1b_tlong.rds", r2$out, fixed = TRUE)))
  f2 <- file.path(out_dir, "fit_mid_lag3_all_pr_M2b_K2_tlong.rds"); expect_true(file.exists(f2)); expect_false(file.exists(file.path(out_dir, "fit_mid_lag3_all_pr_M2b_K2.rds")))
  o2 <- readRDS(f2); expect_equal(o2$args$init_from, "fit_mid_lag3_all_pr_M1b_tlong.rds"); expect_equal(o2$init$m1_cell, "mid_lag3_all_pr_M1b_tlong"); expect_equal(o2$args$tag, "long")
  r3 <- run_script(v2, "02_fit_model.R", short, "--model=M2b", "--tag=other"); expect_equal(r3$status, 0L, info = paste(tail(r3$out, 20), collapse = "\n"))
  expect_true(any(grepl("WARNING: no matching M1 fit (fit_mid_lag3_all_pr_M1b_tother.rds)", r3$out, fixed = TRUE))); expect_true(any(grepl("RANDOM inits", r3$out, fixed = TRUE)))
  expect_equal(readRDS(file.path(out_dir, "fit_mid_lag3_all_pr_M2b_K2_tother.rds"))$args$init_source, "stan_random")
  expect_setequal(files("^fit_"), c("fit_mid_lag3_all_pr_M1b.rds", "fit_mid_lag3_all_pr_M1b_tlong.rds", "fit_mid_lag3_all_pr_M2b_K2_tlong.rds", "fit_mid_lag3_all_pr_M2b_K2_tother.rds"))
  ## bad tags are refused before anything is written
  for (bad in c("--tag=a-b", "--tag=", "--tag")) { rb <- run_script(v2, "02_fit_model.R", short, "--model=M1b", bad, "--dry_run=TRUE"); expect_true(rb$status != 0L, info = bad) }
  expect_length(files("^standata_"), 1L)
  unlink(dirname(v2), recursive = TRUE)
})

test_that("02_fit_model.R: a prior switch fits the prior program with the hyperparameters as data, names and records it, and 03_extract.R accepts the fit", {
  for (m in c("m1_latent_additive", "m1_latent_prior")) if (!file.exists(file.path(root, "stan", "build", paste0(m, ".rds")))) skip(paste("compiled model missing:", m))
  v2 <- make_root(); out_dir <- file.path(v2, "output")
  d <- run_script(v2, "02_fit_model.R", short, "--model=M1b", "--rho_prior=flat", "--sigma_prior=ln05", "--dry_run=TRUE"); expect_equal(d$status, 0L, info = paste(tail(d$out, 20), collapse = "\n"))
  sd <- readRDS(file.path(out_dir, "standata_mid_lag3_all_pr_M1b_rhoflat_sigln05.rds"))
  expect_equal(sd$stan_program, "m1_latent_prior"); expect_equal(sd$dat$prior_rho_shape, 1); expect_equal(sd$dat$prior_sigma_x_sd, 0.5); expect_equal(sd$prior$rho_prior, "flat"); expect_equal(sd$prior$sigma_prior, "ln05")
  d0 <- run_script(v2, "02_fit_model.R", short, "--model=M1b", "--dry_run=TRUE"); expect_equal(d0$status, 0L)
  sd0 <- readRDS(file.path(out_dir, "standata_mid_lag3_all_pr_M1b.rds")); expect_null(sd0$dat$prior_rho_shape); expect_equal(sd0$stan_program, "m1_latent_additive")
  sdx <- sd$dat; sdx$prior_rho_shape <- NULL; sdx$prior_sigma_x_sd <- NULL; expect_identical(sdx, sd0$dat)      # the same data apart from the two hyperparameters
  r <- run_script(v2, "02_fit_model.R", short, "--model=M1b", "--rho_prior=flat"); expect_equal(r$status, 0L, info = paste(tail(r$out, 20), collapse = "\n"))
  expect_true(any(grepl("Stan program: m1_latent_prior (prior sensitivity", r$out, fixed = TRUE))); expect_true(any(grepl("compiled m1_latent_prior matches m1_latent_prior.stan", r$out, fixed = TRUE)))
  f <- file.path(out_dir, "fit_mid_lag3_all_pr_M1b_rhoflat.rds"); expect_true(file.exists(f)); expect_false(file.exists(file.path(out_dir, "fit_mid_lag3_all_pr_M1b.rds")))
  o <- readRDS(f); expect_equal(o$args$stan_program, "m1_latent_prior"); expect_equal(o$args$rho_prior, "flat"); expect_equal(o$args$sigma_prior, "ln1"); expect_equal(o$args$prior_rho_shape, 1); expect_equal(o$args$prior_sigma_x_sd, 1)
  expect_equal(o$prior, list(rho_prior = "flat", sigma_prior = "ln1", prior_rho_shape = 1, prior_sigma_x_sd = 1, stan_program = "m1_latent_prior"))
  expect_equal(o$dat$prior_rho_shape, 1); expect_equal(o$dat$prior_sigma_x_sd, 1); expect_equal(o$spec$rho_prior, "flat")
  expect_equal(o$args$stan_source_sha256, trimws(readLines(file.path(root, "stan", "build", "m1_latent_prior.sha256"))[1]))
  expect_false(o$args$stan_source_sha256 == trimws(readLines(file.path(root, "stan", "build", "m1_latent_additive.sha256"))[1]))
  s <- read.csv(file.path(out_dir, "summary_mid_lag3_all_pr_M1b_rhoflat.csv")); expect_true(all(s$rho_prior == "flat")); expect_true(all(s$sigma_prior == "ln1")); expect_true(all(s$stan_program == "m1_latent_prior"))
  expect_true(assert_fit_stan_current(o, root = root))
  e <- run_script(v2, "03_extract.R", paste0("--fit=", f)); expect_equal(e$status, 0L, info = paste(tail(e$out, 20), collapse = "\n"))
  expect_true(file.exists(file.path(out_dir, "typology_mid_lag3_all_pr_M1b_rhoflat.csv"))); expect_true(any(grepl("prior sensitivity fit (m1_latent_prior): rho_prior = flat (Beta shape 1)", e$out, fixed = TRUE)))
  ## the fit stamped as the main program under the prior switch, or a prior fit claiming the main program, is refused by the guard
  doc <- o; doc$args$stan_program <- "m1_latent_additive"; df <- file.path(out_dir, "fit_mid_lag3_all_pr_M1b_rhoflat.rds"); saveRDS(doc, df)
  e2 <- run_script(v2, "03_extract.R", paste0("--fit=", df)); expect_true(e2$status != 0L); expect_true(any(grepl("not what its name says", e2$out, fixed = TRUE)))
  ## refusals through the script: M2 / M1c with a prior switch, rho with M1a
  expect_true(run_script(v2, "02_fit_model.R", short, "--model=M2b", "--rho_prior=flat", "--dry_run=TRUE")$status != 0L)
  expect_true(run_script(v2, "02_fit_model.R", short, "--model=M1c", "--sigma_prior=ln05", "--dry_run=TRUE")$status != 0L)
  expect_true(run_script(v2, "02_fit_model.R", short, "--model=M1a", "--rho_prior=flat", "--dry_run=TRUE")$status != 0L)
  unlink(dirname(v2), recursive = TRUE)
})

test_that("05_robustness.R refuses a run tag and the prior switches for every preset, before anything runs", {
  tmp <- tempfile("rob_"); dir.create(file.path(tmp, "v2"), recursive = TRUE)
  file.copy(file.path(root, "R"), file.path(tmp, "v2"), recursive = TRUE); write_toy_cell(file.path(tmp, "v2", "data")); dir.create(file.path(tmp, "v2", "output")); dir.create(file.path(tmp, "v2", "logs"))
  rec <- file.path(tmp, "calls.txt")
  writeLines(c("args <- commandArgs(trailingOnly = TRUE)", sprintf("cat(paste(args, collapse = ' '), '\\n', sep = '', file = '%s', append = TRUE)", rec)), file.path(tmp, "v2", "R", "02_fit_model.R"))
  v2 <- file.path(tmp, "v2")
  for (p in names(ROBUSTNESS_PRESETS)) {
    r <- run_script(v2, "05_robustness.R", paste0("--preset=", p), "--tag=long", "--dry_run=TRUE")
    expect_true(r$status != 0L, info = p); expect_true(any(grepl("--tag=long is not an option of 05_robustness.R", r$out, fixed = TRUE)), info = p)
    r <- run_script(v2, "05_robustness.R", paste0("--preset=", p), "--rho_prior=flat", "--dry_run=TRUE")
    expect_true(r$status != 0L, info = p); expect_true(any(grepl("--rho_prior / --sigma_prior are not options of 05_robustness.R", r$out, fixed = TRUE)), info = p)
    r <- run_script(v2, "05_robustness.R", paste0("--preset=", p), "--sigma_prior=ln05", "--dry_run=TRUE"); expect_true(r$status != 0L, info = p)
  }
  expect_false(file.exists(rec))                                                                           # no 02_fit_model.R call was made
  ## the dist preset still runs its four cells, with the prior switches at their defaults on every call
  r <- run_script(v2, "05_robustness.R", "--preset=dist", "--dry_run=TRUE"); expect_equal(r$status, 0L, info = paste(tail(r$out, 20), collapse = "\n"))
  calls <- readLines(rec); expect_length(calls, 4L)
  expect_true(all(grepl("--rho_prior=beta25", calls, fixed = TRUE))); expect_true(all(grepl("--sigma_prior=ln1", calls, fixed = TRUE))); expect_false(any(grepl("--tag", calls, fixed = TRUE)))
  unlink(tmp, recursive = TRUE)
})
