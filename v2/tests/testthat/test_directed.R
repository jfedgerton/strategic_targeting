## The directed comparison specification (PI-approved 2026-11-01; switch comparison = directed, suffix _cdir): the conventional comparison
## with the directed stage-2 Z (Z_DIRECTED: tgt_dem, tgt_dem_x_init_dem, tgt_major, tgt_major_x_init_major in place of joint_dem, one_major,
## both_major; Pz = 9) for M1a / M1b / M1c / M1cr and the AME model M2b (P = 0 legal in m2_ame.stan since 2026-11-01).
## (a) parsing, suffix, cell names, refusals; (b) the exact Z on the toy cell for every model, the identities, within-set variation, and the
## default / conventional cells byte-identical to the 4748bbf code (sourced from git); (c) m2_ame.stan: P = 0 samples for K = 1 and K = 2,
## and for P >= 1 the log density, gradient and generated quantities equal those of the 4748bbf program (compiled from git);
## (d) 02_fit_model.R end to end: dry runs, the M2b init from the M1b fit of the SAME directed cell (never _conv / default), refusals,
## 05_robustness.R refuses the switch; (e) the real lag-3 cells when built.
suppressPackageStartupMessages({ library(testthat); library(dplyr); library(rstan) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "assemble_cell.R"))
## the toy cell with the build's dyad-year status columns (dem1 / dem2, major1 / major2), consistent with the toy's joint_dem (democracies 10, 20),
## major_count (majors 10, 40) and the panel's dem_tm1 / cowmaj_tm1
toy_dirstat <- function(a) a |> mutate(dem1 = as.integer(ccode1 %in% c(10L, 20L)), dem2 = as.integer(ccode2 %in% c(10L, 20L)),
                                       major1 = as.integer(ccode1 %in% c(10L, 40L)), major2 = as.integer(ccode2 %in% c(10L, 40L)))
toy_joiners_dir <- function(e) e |> mutate(orig1 = ifelse(initiator == 30L & target %in% c(20L, 40L), 0L, 1L), is_joiner_pair = orig1 == 0L)   # every stage-2 group non-empty (M1c)
Z_DIR_EXPECTED <- c("contig_land", "log_capdist", "tgt_dem", "tgt_dem_x_init_dem", "alliance_def_off", "tgt_major", "tgt_major_x_init_major", "dy_mid_tm1", "log_dy_prior_mids")
## the 4748bbf spec.R + assemble_cell.R in their own environment (helpers_data.R / helpers_model.R are unchanged and come from the global env)
old_code_env <- function() {
  e <- new.env(parent = globalenv())
  for (f in c("spec.R", "assemble_cell.R")) {
    src <- tryCatch(suppressWarnings(system2("git", c("-C", root, "show", paste0("4748bbf:v2/R/", f)), stdout = TRUE, stderr = FALSE)), error = function(err) NULL)
    if (is.null(src) || !length(src) || !is.null(attr(src, "status"))) return(NULL)
    tf <- tempfile(fileext = ".R"); writeLines(src, tf); sys.source(tf, envir = e); unlink(tf)
  }
  e
}

test_that("spec: comparison = directed adds _cdir in the _conv position; the fixed switches are refused in any value; models M1a / M1b / M1c / M1cr / M2b only", {
  d <- parse_spec("--comparison=directed")
  expect_true(spec_is_directed(d)); expect_false(spec_is_conventional(d)); expect_false(spec_is_main(d)); expect_equal(spec_suffix(d), "_cdir")
  expect_false(spec_is_directed(parse_spec(character(0)))); expect_false(spec_is_directed(parse_spec("--comparison=conventional"))); expect_false(spec_is_directed(NULL)); expect_false(spec_is_directed(list()))
  expect_equal(DIRECTED_MODELS, c("M1a", "M1b", "M1c", "M1cr", "M2b")); expect_equal(Z_DIRECTED, Z_DIR_EXPECTED); expect_length(Z_DIRECTED, 9L)
  for (nm in CONVENTIONAL_FIXED) expect_equal(d[[nm]], SPEC_DEFAULTS[[nm]], info = nm)
  ## the 12 cells of the request
  d_ap <- parse_spec(c("--comparison=directed", "--roles=aprime"))
  for (s in c("fmid", "mid")) {
    expect_equal(cell_name(s, 3L, "M1a", d), sprintf("%s_lag3_all_pr_M1a_cdir", s)); expect_equal(cell_name(s, 3L, "M1b", d), sprintf("%s_lag3_all_pr_M1b_cdir", s))
    expect_equal(cell_name(s, 3L, "M1c", d_ap), sprintf("%s_lag3_all_pr_M1c_cdir_aprime", s)); expect_equal(cell_name(s, 3L, "M1cr", d_ap), sprintf("%s_lag3_all_pr_M1cr_cdir_aprime", s))
    expect_equal(cell_name(s, 3L, "M2b", parse_spec(c("--comparison=directed", "--ame_K=1"))), sprintf("%s_lag3_all_pr_M2b_K1_cdir", s))
    expect_equal(cell_name(s, 3L, "M2b", d), sprintf("%s_lag3_all_pr_M2b_K2_cdir", s))
  }
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec(c("--comparison=directed", "--ww2=drop", "--weights=none"))), "fmid_lag3_all_pr_M1b_noww2_cdir_w1")
  expect_false(any(duplicated(c(unlist(SPEC_SUFFIX), "_K2", "_pilot", "_STALE", "_samp10", "_noCH", "_polrel"))))
  ## the M2b init cell: the M1b cell of the SAME directed spec, for every K
  for (k in 1:3) expect_equal(m1_cell_for_m2("fmid", 3L, "M2b", parse_spec(c("--comparison=directed", sprintf("--ame_K=%d", k)))), "fmid_lag3_all_pr_M1b_cdir")
  ## the switches the spec fixes conflict in ANY value, the default included
  for (bad in c("--controls=minimal", "--controls=full", "--distance=none", "--distance=capital", "--contiguity=land", "--contiguity=none", "--major=none", "--major=joint", "--interact=none", "--interact=prefixed", "--changes=raw", "--changes=winsor"))
    expect_error(parse_spec(c("--comparison=directed", bad)), "conflicts with --comparison=directed", fixed = TRUE, info = bad)
  ## the other models, the role coding with a non-context model, the prior switches with M2b, a bad value
  for (m in c("M1m", "M1r", "M2a", "M2m", "M2r")) expect_error(cell_name("fmid", 3L, m, d), "the directed comparison is specified for M1a, M1b, M1c, M1cr, M2b only", fixed = TRUE, info = m)
  expect_error(cell_name("fmid", 3L, "M2b", d_ap), "the role coding applies to M1c, M1cr, M2c only", fixed = TRUE); expect_error(cell_name("fmid", 3L, "M1b", d_ap), "the role coding applies")
  expect_error(cell_name("fmid", 3L, "M2b", parse_spec(c("--comparison=directed", "--rho_prior=flat"))), "the prior switches apply to the M1 m / r / a / b programs only")
  expect_error(parse_spec("--comparison=cdir"), "--comparison must be one of none | conventional | directed", fixed = TRUE)
  ## the M1-family K rule still holds
  expect_error(cell_name("fmid", 3L, "M1b", parse_spec(c("--comparison=directed", "--ame_K=1"))), "K applies to the M2 family only")
  ## the permutation part string keeps the switch when on; the default and conventional strings are those of 4748bbf
  expect_true(grepl(" comparison=directed ", perm_part_spec_string(d), fixed = TRUE))
  ## (c) existing names and spec strings are those of 4748bbf
  old <- old_code_env(); if (is.null(old)) skip("git history with commit 4748bbf not available")
  combos <- list(character(0), "--comparison=conventional", c("--comparison=conventional", "--roles=aprime"), c("--distance=capital", "--major=joint"), "--ww2=drop", c("--comparison=conventional", "--ww2=drop", "--weights=none", "--tag=x1"),
                 c("--distance=capital", "--contiguity=none"), "--ame_K=1", c("--context=coalition", "--coal_drop=top3"))
  for (cb in combos) {
    sn <- parse_spec(cb); so <- old$parse_spec(cb)
    expect_identical(spec_string(sn), old$spec_string(so), info = paste(cb, collapse = " ")); expect_identical(perm_part_spec_string(sn), old$perm_part_spec_string(so), info = paste(cb, collapse = " "))
    for (m in c("M1a", "M1b", "M1c", "M1cr", "M2b", "M1m")) {
      nn <- tryCatch(cell_name("fmid", 3L, m, cell_spec(m, sn)), error = function(e) "ERR"); no <- tryCatch(old$cell_name("fmid", 3L, m, old$cell_spec(m, so)), error = function(e) "ERR")
      expect_identical(nn, no, info = paste(m, paste(cb, collapse = " ")))
    }
  }
})

test_that("assembly on the toy cell: exact Z (Pz = 9), W and no X under M1a, M1b, M1c, M1cr and M2b; the identities on every row; within-set variation; the refusals", {
  td <- write_toy_cell(tempfile("toycell"), alts_edit = toy_dirstat, events_edit = toy_joiners_dir)
  d <- parse_spec("--comparison=directed"); cv <- parse_spec("--comparison=conventional")
  out <- lapply(c(M1a = "M1a", M1b = "M1b", M1c = "M1c", M1cr = "M1cr", M2b = "M2b"), function(m) quiet(assemble_cell_spec("mid", 3L, cell_spec(m, d), data_dir = td, model = m)))
  C <- quiet(assemble_cell_spec("mid", 3L, cv, data_dir = td, model = "M1b"))
  for (m in names(out)) {
    A <- out[[m]]
    expect_identical(A$Z_cols, Z_DIR_EXPECTED, info = m); expect_identical(A$W_cols, W_CONVENTIONAL, info = m); expect_identical(A$X_cols, character(0), info = m); expect_length(A$x_ref, 0L)
    expect_true(A$directed); expect_false(A$conventional); expect_equal(A$distance, "capital"); expect_equal(A$major, "joint"); expect_equal(A$contiguity, "land")
    al <- A$alts
    ## the definitions and the identities, on every row
    expect_identical(al$tgt_dem, as.integer(al$dem2)); expect_identical(al$tgt_major, as.integer(al$major2))
    expect_identical(al$tgt_dem_x_init_dem, as.integer(al$dem2 * al$dem1)); expect_identical(al$tgt_major_x_init_major, as.integer(al$major2 * al$major1))
    expect_true(all(al$tgt_dem_x_init_dem == al$joint_dem))
    expect_true(all(al$one_major == al$tgt_major + al$major1 - 2L * al$tgt_major * al$major1)); expect_true(all(al$both_major == al$tgt_major * al$major1))
    ## the initiator status is the panel's (uncentred): democracies 10, 20; majors 10, 40
    expect_equal(al$dem1, as.integer(al$ccode1 %in% c(10L, 20L))); expect_equal(al$major1, as.integer(al$ccode1 %in% c(10L, 40L)))
    ## 0/1 and uncentred; the continuous Z columns centred over the IV rows as under the conventional spec
    for (cc in c("contig_land", "tgt_dem", "tgt_dem_x_init_dem", "alliance_def_off", "tgt_major", "tgt_major_x_init_major", "dy_mid_tm1")) expect_true(all(al[[cc]] %in% c(0L, 1L)), info = cc)
    expect_equal(al$log_dy_prior_mids, log1p(al$dy_prior_mids_h) - mean(log1p(al$dy_prior_mids_h)[al$in_iv == 1L]))
    expect_lt(abs(mean(al$log_capdist[al$in_iv == 1L])), 1e-12)
    ## the tgt_* columns vary within choice sets (IV sets and event-bearing stage-2 sets)
    v <- al |> group_by(sy_row) |> summarise(across(c(tgt_dem, tgt_dem_x_init_dem, tgt_major, tgt_major_x_init_major), ~ n_distinct(.x) > 1L))
    for (cc in c("tgt_dem", "tgt_dem_x_init_dem", "tgt_major", "tgt_major_x_init_major")) expect_true(any(v[[cc]]), info = cc)
    expect_true(all(A$directed_stats$var_iv_sets[Z_DIRECTED_TGT] > 0)); expect_true(all(A$directed_stats$var_s2_sets[Z_DIRECTED_TGT] > 0))
    ## the Stan data: P = 0, X N1 x 0, Pz = 9, P1 = 8
    dat <- build_stan_data(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, x_ref = A$x_ref, K = if (m == "M2b") 2L else NULL, allow_orphan_initiations = TRUE)
    expect_equal(dat$P, 0L); expect_equal(dim(dat$X), c(dat$N1, 0L)); expect_equal(dat$Pz, 9L); expect_equal(dat$P1, 8L); expect_length(dat$x_ref, 0L)
    ## everything but the Z swap is the conventional cell: rows, events, weights, W, and the shared Z columns
    expect_identical(A$sy, C$sy); expect_identical(A$ev[, c("sy_row", "alt_row", "w")], C$ev[, c("sy_row", "alt_row", "w")])
    for (cc in intersect(A$Z_cols, C$Z_cols)) expect_identical(al[[cc]], C$alts[[cc]], info = cc)
  }
  expect_identical(out$M1c$ev, out$M1cr$ev); expect_identical(out$M1a$alts, out$M2b$alts)
  ## the refusals: missing / inconsistent status columns, constant target columns, the switch requirements, lag != 3, exclusive flags
  td2 <- write_toy_cell(tempfile("toycell"))                                              # no dem1 / dem2 / major1 / major2
  expect_error(quiet(assemble_cell_spec("mid", 3L, d, data_dir = td2, model = "M1b")), "directed comparison: the alternatives lack dem1, dem2, major1, major2")
  td3 <- write_toy_cell(tempfile("toycell"), alts_edit = function(a) toy_dirstat(a) |> mutate(dem1 = ifelse(ccode1 == 30L, 1L, dem1)))
  expect_error(quiet(assemble_cell_spec("mid", 3L, d, data_dir = td3, model = "M1b")), "dem1 differs from the panel's dem_tm1")
  td4 <- write_toy_cell(tempfile("toycell"), alts_edit = function(a) toy_dirstat(a) |> mutate(major1 = ifelse(ccode1 == 20L, 1L, major1)))
  expect_error(quiet(assemble_cell_spec("mid", 3L, d, data_dir = td4, model = "M1b")), "major1 differs from the panel's cowmaj_tm1")
  td5 <- write_toy_cell(tempfile("toycell"), alts_edit = function(a) toy_dirstat(a) |> mutate(joint_dem = ifelse(ccode1 == 10L & ccode2 == 30L, 1L, joint_dem)))
  expect_error(quiet(assemble_cell_spec("mid", 3L, d, data_dir = td5, model = "M1b")), "tgt_dem_x_init_dem != joint_dem")
  td6 <- write_toy_cell(tempfile("toycell"), alts_edit = function(a) toy_dirstat(a) |> mutate(major2 = as.integer(ccode2 == 10L)))   # major2 inconsistent with major_count
  expect_error(quiet(assemble_cell_spec("mid", 3L, d, data_dir = td6, model = "M1b")), "one_major != tgt_major")
  td7 <- write_toy_cell(tempfile("toycell"), alts_edit = function(a) toy_dirstat(a) |> mutate(dem2 = 0L, joint_dem = 0L))           # no democratic target anywhere
  expect_error(quiet(assemble_cell_spec("mid", 3L, d, data_dir = td7, model = "M1b")), "a target column does not vary within any choice set")
  expect_error(assemble_cell("mid", 3L, "all", DIRECTED = TRUE, data_dir = td), "DIRECTED = TRUE needs DISTANCE = capital, MAJOR = joint")
  expect_error(assemble_cell("mid", 3L, "all", DIRECTED = TRUE, CONVENTIONAL = TRUE, DISTANCE = "capital", MAJOR = "joint", data_dir = td), "exclusive")
  expect_error(assemble_cell("mid", 3L, "all", DIRECTED = NA, data_dir = td))
  expect_error(assemble_cell_spec("mid", 1L, d, data_dir = td, model = "M1b"), "--comparison=directed is defined on the lag-3 data", fixed = TRUE)
  expect_error(assemble_cell_spec("mid", 3L, d, data_dir = td, model = "M2a"))
  unlink(c(td, td2, td3, td4, td5, td6, td7), recursive = TRUE)
})

test_that("existing cells are unchanged: default, M1c, _dist_major, conventional (M1b, M1c) and roles cells assemble byte-identically to the 4748bbf code", {
  old <- old_code_env(); if (is.null(old)) skip("git history with commit 4748bbf not available")
  for (edit in list(identity, toy_dirstat)) {                                             # with and without the status columns in the alternatives file
    td <- write_toy_cell(tempfile("toycell"), alts_edit = edit, events_edit = toy_joiners_dir)
    cases <- list(list(character(0), "M1b"), list(character(0), "M1c"), list(character(0), "M2b"), list(c("--distance=capital", "--major=joint"), "M1b"),
                  list("--comparison=conventional", "M1b"), list("--comparison=conventional", "M1c"), list("--comparison=conventional", "M1cr"), list("--ww2=drop", "M1a"), list("--controls=full", "M1b"))
    for (cs in cases) {
      sn <- parse_spec(cs[[1]]); so <- old$parse_spec(cs[[1]])
      An <- quiet(assemble_cell_spec("mid", 3L, sn, data_dir = td, model = cs[[2]])); Ao <- quiet(old$assemble_cell_spec("mid", 3L, so, data_dir = td, model = cs[[2]]))
      expect_identical(setdiff(names(An), names(Ao)), c("directed", "directed_stats", "joiner_set", "ev_sets", "jset_stats", "ev_joiner", "variant", "mid3_stats", "reach", "reach_stats", "subsample", "subsample_stats", "init_sub_stats"))   # joiner sets 2026-11-04; run 0 / mid3 2026-11-10 (all NULL / none / the variant name here)
      expect_equal(An$joiner_set, "all"); expect_null(An$ev_sets); expect_null(An$jset_stats)
      for (nm in names(Ao)) expect_identical(An[[nm]], Ao[[nm]], info = paste(nm, cs[[2]], paste(cs[[1]], collapse = " ")))
      expect_false(An$directed); expect_null(An$directed_stats)
    }
    unlink(td, recursive = TRUE)
  }
})

test_that("m2_ame.stan: P >= 1 gives the log density, gradient and generated quantities of the 4748bbf program; P = 0 samples for K = 1 and K = 2 and equals an all-zero X", {
  sm2 <- skip_if_no_build("m2_ame")
  old_src <- tryCatch(suppressWarnings(system2("git", c("-C", root, "show", "4748bbf:v2/stan/m2_ame.stan"), stdout = TRUE, stderr = FALSE)), error = function(e) NULL)
  if (is.null(old_src) || !length(old_src) || !is.null(attr(old_src, "status"))) skip("git history with commit 4748bbf not available")
  expect_true(any(grepl("^  int<lower=1> P;$", old_src))); expect_true(any(grepl("^  int<lower=0> P;", readLines(file.path(root, "stan", "m2_ame.stan")))))
  ## the source differs from 4748bbf in the P declaration line only
  new_src <- readLines(file.path(root, "stan", "m2_ame.stan")); expect_equal(length(new_src), length(old_src))
  expect_equal(which(new_src != old_src), which(grepl("^  int<lower=1> P;$", old_src)))
  sm_old <- stan_model(model_code = paste(old_src, collapse = "\n"), model_name = "m2_ame_4748bbf", auto_write = FALSE)
  rstan_options(threads_per_chain = 1)
  for (cfg in list(list(joint = 1L, K = 2L, D_mode = 3L), list(joint = 0L, K = 1L, D_mode = 3L), list(joint = 1L, K = 1L, D_mode = 1L))) {
    dat <- tiny_data("M2b", joint = cfg$joint, K = cfg$K, D_mode = cfg$D_mode); expect_gte(dat$P, 1L)
    init <- stan_inits(dat, make_pars(dat, with_uv = TRUE))$init
    f_old <- sampling(sm_old, data = dat, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(init), seed = 123, refresh = 0)
    f_new <- sampling(sm2, data = dat, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(init), seed = 123, refresh = 0)
    expect_equal(f_old@mode, 0L); expect_equal(f_new@mode, 0L)
    expect_equal(as.numeric(rstan::extract(f_new, "sigma_b")$sigma_b), init$sigma_b); expect_equal(as.numeric(rstan::extract(f_old, "sigma_b")$sigma_b), init$sigma_b)   # the fixed point is the init
    np <- get_num_upars(f_new); expect_equal(np, get_num_upars(f_old))
    set.seed(11)
    for (s in 1:5) {
      u <- rnorm(np, 0, 0.5)
      for (adj in c(TRUE, FALSE)) expect_identical(log_prob(f_old, u, adjust_transform = adj), log_prob(f_new, u, adjust_transform = adj))
      expect_identical(grad_log_prob(f_old, u), grad_log_prob(f_new, u))
    }
    for (par in c("log_lik_stage1", "log_lik_stage2", "mu_state", "mu_event", "Sigma", "Omega", "R_Om", "U", "V")) expect_identical(as.numeric(rstan::extract(f_old, par)[[par]]), as.numeric(rstan::extract(f_new, par)[[par]]), info = paste(par, cfg$K, cfg$joint))
  }
  ## P = 0 equals P = 1 with an all-zero X column: identical log-liks and ideal points (beta multiplies zeros)
  dat1 <- tiny_data("M2b", joint = 1L, K = 2L)
  dat0 <- dat1; dat0$P <- 0L; dat0$X <- matrix(numeric(0), dat1$N1, 0L); dat0$x_ref <- numeric(0)
  datz <- dat1; datz$P <- 1L; datz$X <- matrix(0, dat1$N1, 1L); datz$x_ref <- array(0, 1)
  i1 <- stan_inits(dat1, make_pars(dat1, with_uv = TRUE))$init; i0 <- i1; i0$beta <- matrix(numeric(0), 0L, 2L); iz <- i1; iz$beta <- matrix(c(0.3, -0.2), 1L, 2L)
  f0 <- sampling(sm2, data = dat0, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(i0), seed = 123, refresh = 0)
  fz <- sampling(sm2, data = datz, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(iz), seed = 123, refresh = 0)
  expect_equal(f0@mode, 0L); expect_equal(f0@par_dims$beta, c(0, 2))
  expect_equal(as.numeric(rstan::extract(f0, "sigma_b")$sigma_b), i1$sigma_b); expect_equal(as.numeric(rstan::extract(fz, "sigma_b")$sigma_b), i1$sigma_b)   # the inits were used
  for (par in c("log_lik_stage1", "log_lik_stage2", "mu_state", "mu_event")) expect_identical(as.numeric(rstan::extract(f0, par)[[par]]), as.numeric(rstan::extract(fz, par)[[par]]), info = par)
  ## P = 0 samples on the directed toy cell for K = 1 and K = 2 (M2b), mu_state = mu0 + a_i = mu_event of the event's state
  td <- write_toy_cell(tempfile("toycell"), alts_edit = toy_dirstat)
  for (K in 1:2) {
    A <- quiet(assemble_cell_spec("mid", 3L, parse_spec(c("--comparison=directed", sprintf("--ame_K=%d", K))), data_dir = td, model = "M2b"))
    dat <- build_stan_data(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, x_ref = A$x_ref, K = K, allow_orphan_initiations = TRUE)
    expect_equal(dat$P, 0L); expect_equal(dat$Pz, 9L); expect_equal(dat$K, K)
    fit <- sampling(sm2, data = dat, chains = 1, iter = 40, warmup = 20, seed = 123, refresh = 0, pars = c("z_re", "z_b", "V_raw", "U_free", "L_Om_joint", "L_Om", "re"), include = FALSE)
    expect_equal(fit@mode, 0L, info = paste("K", K)); expect_equal(fit@par_dims$beta, c(0, 2)); expect_equal(fit@par_dims$gamma, 9); expect_equal(fit@par_dims$U, c(dat$S, K))
    e <- rstan::extract(fit, pars = c("mu0", "a_i", "mu_state", "mu_event", "log_lik_stage2"), permuted = TRUE)
    expect_true(all(is.finite(e$mu_state))); expect_true(all(is.finite(e$log_lik_stage2)))
    for (s in seq_len(dat$S)) expect_equal(e$mu_state[, s, ], e$mu0 + e$a_i[, s, ], tolerance = 1e-12, check.attributes = FALSE)
    for (ev in seq_len(dat$E)) expect_equal(e$mu_event[, ev, ], e$mu_state[, dat$sy_state[dat$ev_sy[ev]], ])
  }
  unlink(td, recursive = TRUE)
})

test_that("02_fit_model.R --comparison=directed end to end on the toy cell: dry runs, M2b inits from the matching directed M1b only, refusals; 05_robustness.R refuses the switch", {
  for (m in c("m1_context", "m1_context_rho", "m1_latent_additive", "m2_ame")) if (!file.exists(file.path(root, "stan", "build", paste0(m, ".rds")))) skip(paste("compiled model missing:", m))
  tmp <- tempfile("smoke_dir_"); dir.create(file.path(tmp, "v2"), recursive = TRUE)
  file.copy(file.path(root, "R"), file.path(tmp, "v2"), recursive = TRUE)
  ok <- file.symlink(file.path(root, "stan"), file.path(tmp, "v2", "stan")); if (!ok) file.copy(file.path(root, "stan"), file.path(tmp, "v2"), recursive = TRUE)
  write_toy_cell(file.path(tmp, "v2", "data"), alts_edit = toy_dirstat, events_edit = toy_joiners_dir)
  dir.create(file.path(tmp, "v2", "output")); dir.create(file.path(tmp, "v2", "logs"))
  v2 <- file.path(tmp, "v2"); out_dir <- file.path(v2, "output")
  run <- function(script, ...) { out <- suppressWarnings(system2("Rscript", c(file.path(v2, "R", script), ...), stdout = TRUE, stderr = TRUE)); list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  prog_of <- c(M1a = "m1_latent_additive", M1b = "m1_latent_additive", M1c = "m1_context", M1cr = "m1_context_rho", M2b = "m2_ame")
  ## (a) dry runs: exact columns, P = 0, Pz = 9, cell names (M1c / M1cr with the default stage-b roles and with aprime)
  dry <- list(list("M1a", character(0), "M1a_cdir"), list("M1b", character(0), "M1b_cdir"), list("M1c", character(0), "M1c_cdir"), list("M1cr", character(0), "M1cr_cdir"),
              list("M2b", "--ame_K=1", "M2b_K1_cdir"), list("M2b", "--ame_K=2", "M2b_K2_cdir"))
  for (dd in dry) {
    r0 <- run("02_fit_model.R", "--sample=mid", "--lag=3", paste0("--model=", dd[[1]]), "--comparison=directed", dd[[2]], "--dry_run=TRUE")
    expect_equal(r0$status, 0L, info = paste(dd[[1]], paste(tail(r0$out, 30), collapse = "\n")))
    sd <- readRDS(file.path(out_dir, sprintf("standata_mid_lag3_all_pr_%s.rds", dd[[3]])))
    expect_identical(sd$Z_cols, Z_DIR_EXPECTED); expect_identical(sd$W_cols, W_CONVENTIONAL); expect_identical(sd$X_cols, character(0))
    expect_equal(sd$dat$P, 0L); expect_equal(sd$dat$Pz, 9L); expect_equal(sd$dat$P1, 8L); expect_true(sd$directed); expect_false(sd$conventional); expect_equal(sd$spec$comparison, "directed")
    expect_equal(sd$cell, sprintf("mid_lag3_all_pr_%s", dd[[3]])); expect_equal(sd$stan_program, prog_of[[dd[[1]]]])
    if (dd[[1]] == "M2b") expect_equal(sd$dat$K, as.integer(sub("--ame_K=", "", dd[[2]])))
    expect_true(any(grepl("directed comparison (2026-11-01, suffix _cdir): W = ", r0$out, fixed = TRUE)))
  }
  r0 <- run("02_fit_model.R", "--sample=mid", "--lag=3", "--model=M1cr", "--comparison=directed", "--roles=aprime", "--dry_run=TRUE")   # the toy has no roles file: refused loudly
  expect_true(r0$status != 0L)
  ## refusals: the other M2 models, the fixed switches, the roles switch with M2b, a bad value
  for (bad in list("--model=M2a", "--model=M2m", "--model=M2r", "--model=M1m", c("--model=M2b", "--major=joint"), c("--model=M1b", "--controls=minimal"), c("--model=M2b", "--distance=capital"), c("--model=M2b", "--roles=aprime"))) {
    r <- run("02_fit_model.R", "--sample=mid", "--lag=3", "--comparison=directed", "--dry_run=TRUE", bad); expect_true(r$status != 0L, info = paste(bad, collapse = " "))
  }
  ## (b) the conventional and default M1b fits exist; a directed M2b run must NOT start from either of them (random-init fallback, said loudly)
  for (a in list("--comparison=conventional", character(0))) {
    r <- run("02_fit_model.R", "--sample=mid", "--lag=3", "--model=M1b", a, "--chains=2", "--iter=60", "--warmup=30", "--max_treedepth=6"); expect_equal(r$status, 0L, info = paste(tail(r$out, 30), collapse = "\n"))
  }
  expect_true(file.exists(file.path(out_dir, "fit_mid_lag3_all_pr_M1b_conv.rds"))); expect_true(file.exists(file.path(out_dir, "fit_mid_lag3_all_pr_M1b.rds")))
  expect_true(is.na(locate_m1_fit("mid", 3L, "M2b", parse_spec(c("--comparison=directed", "--ame_K=1")), out_dir)))
  r <- run("02_fit_model.R", "--sample=mid", "--lag=3", "--model=M2b", "--comparison=directed", "--ame_K=1", "--chains=2", "--iter=60", "--warmup=30", "--max_treedepth=6")
  expect_equal(r$status, 0L, info = paste(tail(r$out, 30), collapse = "\n"))
  o <- readRDS(file.path(out_dir, "fit_mid_lag3_all_pr_M2b_K1_cdir.rds")); expect_equal(o$init$source, "stan_random"); expect_true(any(grepl("fit_mid_lag3_all_pr_M1b_cdir.rds", r$out, fixed = TRUE)))
  ## naming the conventional fit explicitly is an error, never an init
  r <- run("02_fit_model.R", "--sample=mid", "--lag=3", "--model=M2b", "--comparison=directed", "--ame_K=1", "--chains=2", "--iter=60", "--warmup=30", paste0("--init_from=", file.path(out_dir, "fit_mid_lag3_all_pr_M1b_conv.rds")))
  expect_true(r$status != 0L)
  ## (c) the directed M1b fit, then M2b K = 1 and K = 2 start from it
  r <- run("02_fit_model.R", "--sample=mid", "--lag=3", "--model=M1b", "--comparison=directed", "--chains=2", "--iter=60", "--warmup=30", "--max_treedepth=6"); expect_equal(r$status, 0L, info = paste(tail(r$out, 30), collapse = "\n"))
  o1 <- readRDS(file.path(out_dir, "fit_mid_lag3_all_pr_M1b_cdir.rds")); expect_true(o1$args$directed); expect_identical(o1$Z_cols, Z_DIR_EXPECTED); expect_equal(o1$dat$P, 0L)
  expect_equal(basename(locate_m1_fit("mid", 3L, "M2b", parse_spec(c("--comparison=directed", "--ame_K=2")), out_dir)), "fit_mid_lag3_all_pr_M1b_cdir.rds")
  for (K in 1:2) {
    r <- run("02_fit_model.R", "--sample=mid", "--lag=3", "--model=M2b", "--comparison=directed", sprintf("--ame_K=%d", K), "--chains=2", "--iter=60", "--warmup=30", "--max_treedepth=6")
    expect_equal(r$status, 0L, info = paste(K, paste(tail(r$out, 30), collapse = "\n")))
    o2 <- readRDS(file.path(out_dir, sprintf("fit_mid_lag3_all_pr_M2b_K%d_cdir.rds", K)))
    expect_equal(o2$init$source, "m1_posterior_means"); expect_equal(o2$init$m1_fit, "fit_mid_lag3_all_pr_M1b_cdir.rds"); expect_equal(o2$init$m1_cell, "mid_lag3_all_pr_M1b_cdir")
    expect_equal(o2$dat$P, 0L); expect_equal(o2$dat$K, K); expect_identical(o2$Z_cols, Z_DIR_EXPECTED); expect_true(o2$args$directed); expect_equal(o2$args$stan_program, "m2_ame")
    expect_equal(o2$fit@par_dims$beta, c(0, 2)); expect_equal(o2$args$iter, 60L)
    expect_true(any(grepl("the init fit is the M1 fit of the same directed cell (mid_lag3_all_pr_M1b_cdir)", r$out, fixed = TRUE)))
    s2 <- read.csv(file.path(out_dir, sprintf("summary_mid_lag3_all_pr_M2b_K%d_cdir.csv", K)))
    expect_false(any(grepl("^beta\\[", s2$parameter))); expect_equal(s2$covariate[s2$parameter == "gamma[3]"], "tgt_dem"); expect_equal(s2$covariate[s2$parameter == "gamma[7]"], "tgt_major_x_init_major")
  }
  ## the M2 sampler defaults are untouched (SAMPLER$M2)
  expect_equal(sampler_settings("M2b"), list(chains = 4L, iter = 4000L, warmup = 2000L, adapt_delta = 0.95, max_treedepth = 10L))
  ## (d) the robustness driver refuses the switch for every preset; no call is made
  rec <- file.path(tmp, "calls.txt")
  writeLines(c("args <- commandArgs(trailingOnly = TRUE)", sprintf("cat(paste(args, collapse = ' '), '\\n', sep = '', file = '%s', append = TRUE)", rec)), file.path(v2, "R", "02_fit_model.R"))
  for (ps in c("noww2", "dist", "context_b", "major_ix", "perm_error")) { dr <- run("05_robustness.R", paste0("--preset=", ps), "--comparison=directed", "--dry_run=TRUE"); expect_true(dr$status != 0L, info = ps); expect_true(any(grepl("comparison=directed is not an option of 05_robustness.R", dr$out, fixed = TRUE)), info = ps) }
  expect_false(file.exists(rec))
  unlink(tmp, recursive = TRUE)
})

test_that("real data (lag 3, when built): the directed fmid / mid cells assemble with Z_DIRECTED and the conventional rows; identities hold on every row", {
  man_path <- file.path(v2_root(), "data", "manifest.csv")
  skip_if_not(file.exists(man_path) && file.exists(file.path(v2_root(), "data", "alts_lag3_all.rds")), "real data files not built")
  d <- parse_spec("--comparison=directed")
  for (s in c("fmid", "mid")) {
    A <- quiet(assemble_cell_spec(s, 3L, d, model = "M2b"))
    expect_identical(A$Z_cols, Z_DIR_EXPECTED); expect_identical(A$W_cols, W_CONVENTIONAL); expect_identical(A$X_cols, character(0))
    expect_equal(nrow(A$sy), 14597L); expect_equal(nrow(A$alts), 1720110L); expect_equal(nrow(A$ev), if (s == "fmid") 1055L else 3258L)
    expect_true(all(A$alts$tgt_dem_x_init_dem == A$alts$joint_dem)); expect_true(all(A$alts$both_major == A$alts$tgt_major_x_init_major))
    expect_true(all(A$alts$one_major == A$alts$tgt_major + A$alts$major1 - 2L * A$alts$tgt_major * A$alts$major1))
    expect_true(all(A$directed_stats$var_s2_sets[Z_DIRECTED_TGT] > 0)); expect_false(anyNA(as.matrix(A$alts[, A$Z_cols])))
    rm(A); invisible(gc())
  }
})
