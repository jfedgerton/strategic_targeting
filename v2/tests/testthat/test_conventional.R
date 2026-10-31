## The conventional comparison specification (PI-approved 2026-10-31; switch comparison = conventional, suffix _conv): parsing and the refusals
## of the switches it fixes, the exact W / Z / X column sets under M1a / M1b / M1c / M1cr on the toy cell (helper-toy-cell.R), no peace years,
## P = 0 in the three M1 Stan programs (a tiny fit of m1_latent_additive and m1_context), the end-to-end 02_fit_model.R run on the toy cell,
## the robustness driver's refusal, and the invariance of every existing cell (spec string form, column sets, Stan data) at the default.
suppressPackageStartupMessages({ library(testthat); library(dplyr); library(rstan) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "assemble_cell.R"))
toy_joiners_conv <- function(e) e |> mutate(orig1 = ifelse(initiator == 30L & target %in% c(20L, 40L), 0L, 1L), is_joiner_pair = orig1 == 0L)   # 30 joins: every stage-2 group non-empty (M1c)

test_that("spec: comparison defaults to none; conventional adds _conv after _ix and before _w1; the switches it fixes are refused in any value; models M1a / M1b / M1c only", {
  main <- parse_spec(character(0))
  expect_equal(main$comparison, "none"); expect_true(spec_is_main(main)); expect_false(spec_is_conventional(main)); expect_false(spec_is_conventional(NULL)); expect_false(spec_is_conventional(list()))
  expect_equal(SPEC_CHOICES$comparison, c("none", "conventional")); expect_equal(SPEC_SUFFIX$comparison, c(conventional = "_conv"))
  expect_equal(CONVENTIONAL_MODELS, c("M1a", "M1b", "M1c", "M1cr")); expect_equal(CONVENTIONAL_FIXED, c("controls", "distance", "contiguity", "major", "interact", "changes"))
  conv <- parse_spec("--comparison=conventional")
  expect_true(spec_is_conventional(conv)); expect_false(spec_is_main(conv)); expect_equal(spec_suffix(conv), "_conv")
  for (nm in CONVENTIONAL_FIXED) expect_equal(conv[[nm]], SPEC_DEFAULTS[[nm]], info = nm)       # the fixed switches stay at their defaults in the spec
  expect_equal(cell_name("fmid", 3L, "M1a", conv), "fmid_lag3_all_pr_M1a_conv"); expect_equal(cell_name("mid", 3L, "M1b", conv), "mid_lag3_all_pr_M1b_conv")
  expect_equal(cell_name("fmid", 3L, "M1c", conv), "fmid_lag3_all_pr_M1c_conv"); expect_equal(cell_name("mid", 3L, "M1cr", conv), "mid_lag3_all_pr_M1cr_conv")
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec(c("--comparison=conventional", "--ww2=drop", "--weights=none"))), "fmid_lag3_all_pr_M1b_noww2_conv_w1")
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec(c("--comparison=conventional", "--rho_prior=flat", "--tag=x1"))), "fmid_lag3_all_pr_M1b_conv_rhoflat_tx1")
  expect_equal(which(names(SPEC_SUFFIX) == "comparison"), which(names(SPEC_SUFFIX) == "interact") + 1L); expect_equal(names(SPEC_SUFFIX)[which(names(SPEC_SUFFIX) == "comparison") + 2L], "weights")   # roles (2026-10-31) sits between
  expect_false(any(duplicated(c(unlist(SPEC_SUFFIX), "_K2", "_pilot", "_STALE", "_samp10", "_noCH", "_polrel"))))
  ## the switches the spec fixes conflict in ANY value, the default included
  for (bad in c("--controls=minimal", "--controls=full", "--distance=none", "--distance=capital", "--contiguity=land", "--contiguity=none", "--major=none", "--major=joint", "--interact=none", "--interact=prefixed", "--changes=raw", "--changes=winsor"))
    expect_error(parse_spec(c("--comparison=conventional", bad)), "conflicts with --comparison=conventional", fixed = TRUE, info = bad)
  expect_error(parse_spec(c("--distance=capital", "--major=joint", "--comparison=conventional")), "--distance=capital conflicts", fixed = TRUE)
  ## the data switches that do not touch W / Z / X combine with it
  ok <- parse_spec(c("--comparison=conventional", "--pairs=originators", "--context=coalition", "--coal_drop=top3", "--ww2=drop", "--weights=none", "--direction=sideA"))
  expect_true(spec_is_conventional(ok)); expect_equal(spec_suffix(ok), "_orig_co_drop3_noww2_conv_w1_sideA")
  ## bad values, the bare flag, and the other models
  expect_error(parse_spec("--comparison=TRUE"), "--comparison must be one of none | conventional", fixed = TRUE); expect_error(parse_spec("--comparison"), "--comparison needs a value", fixed = TRUE)
  for (m in c("M1m", "M1r", "M2a", "M2b", "M2m", "M2r")) expect_error(cell_name("fmid", 3L, m, conv), "the conventional comparison is specified for M1a, M1b, M1c, M1cr only", fixed = TRUE, info = m)
  expect_error(cell_name("fmid", 3L, "M1a", parse_spec(c("--comparison=conventional", "--rho_prior=flat"))), "rho prior applies to the joint model M1b only")   # the existing rule still holds
  ## the permutation part string: the 7c4645e form at the default, the switch kept when on
  sp <- parse_spec("--distance=capital")
  expect_identical(perm_part_spec_string(sp), "pairs=all context=all coal_ties=inclusive coal_drop=none direction=cowwar ww2=include weights=dispute controls=minimal major=none interact=none distance=capital contiguity=land changes=raw ame_K=2")
  expect_identical(perm_part_spec_string(parse_spec(c("--distance=capital", "--rho_prior=flat", "--tag=a"))), perm_part_spec_string(sp))
  expect_true(grepl(" comparison=conventional ", perm_part_spec_string(conv), fixed = TRUE)); expect_false(grepl("comparison", perm_part_spec_string(sp), fixed = TRUE))
  ## the column constants
  expect_equal(W_CONVENTIONAL, c("M_lvl", "R_lvl", "cowmaj_tm1", "dem_tm1", "ongoing_tm1", "log_prior_mids", "log_n_system_tm1", "multi_mid_jan"))
  expect_equal(Z_CONVENTIONAL, c("contig_land", "log_capdist", "joint_dem", "alliance_def_off", "one_major", "both_major", "dy_mid_tm1", "log_dy_prior_mids"))
  expect_false(any(CHANGE_COLUMNS %in% c(W_CONVENTIONAL, Z_CONVENTIONAL))); expect_false(any(grepl("peace", c(W_CONVENTIONAL, Z_CONVENTIONAL)))); expect_equal(PEACE_COLUMNS, c("log_peace_years", "log_dy_peace_years"))
})

test_that("assembly on the toy cell: exact W / Z / X under M1a, M1b, M1c and M1cr; no peace years; the main cells are unchanged", {
  td <- write_toy_cell(tempfile("toycell"), events_edit = toy_joiners_conv)
  conv <- parse_spec("--comparison=conventional"); main <- parse_spec(character(0))
  panel <- readRDS(file.path(td, "stage1_panel_lag3.rds")); alts0 <- readRDS(file.path(td, "alts_lag3_all.rds"))
  out <- lapply(c(M1a = "M1a", M1b = "M1b", M1c = "M1c", M1cr = "M1cr"), function(m) quiet(assemble_cell_spec("mid", 3L, conv, data_dir = td, model = m)))
  for (m in names(out)) {
    A <- out[[m]]
    expect_identical(A$W_cols, W_CONVENTIONAL, info = m); expect_identical(A$Z_cols, Z_CONVENTIONAL, info = m); expect_identical(A$X_cols, character(0), info = m)
    expect_true(A$conventional); expect_equal(A$distance, "capital"); expect_equal(A$major, "joint"); expect_equal(A$contiguity, "land"); expect_false(A$full_controls)
    expect_false(any(CHANGE_COLUMNS %in% c(A$W_cols, A$X_cols, A$Z_cols)))
    expect_length(A$x_ref, 0L)
    expect_false(anyNA(as.matrix(A$sy[, A$W_cols]))); expect_false(anyNA(as.matrix(A$alts[, A$Z_cols])))
    ## (d) no peace years in either stage (PI decision 2026-10-31); the prior-MID counts and lagged conflict stay
    expect_false(any(grepl("peace", c(A$W_cols, A$Z_cols)))); expect_true(all(c("ongoing_tm1", "log_prior_mids") %in% A$W_cols)); expect_true(all(c("dy_mid_tm1", "log_dy_prior_mids") %in% A$Z_cols))
    expect_equal(A$alts$log_dy_prior_mids, log1p(A$alts$dy_prior_mids_h) - mean(log1p(A$alts$dy_prior_mids_h)[A$alts$in_iv == 1L]))   # centred over the IV rows
    ## the 0/1 columns of Z stay 0/1; the W indicators and multi_mid_jan are grand-centred
    for (cc in c("contig_land", "joint_dem", "alliance_def_off", "one_major", "both_major", "dy_mid_tm1")) expect_true(all(A$alts[[cc]] %in% c(0L, 1L)), info = cc)
    for (cc in c("cowmaj_tm1", "dem_tm1", "ongoing_tm1", "multi_mid_jan", "log_prior_mids", "log_n_system_tm1")) expect_lt(abs(mean(A$sy[[cc]])), 1e-12)
    expect_equal(A$sy$multi_mid_jan, panel$multi_mid_jan[match(A$sy$ccode, panel$ccode)] - mean(panel$multi_mid_jan))
    ## the Stan data: P = 0, X is N1 x 0, x_ref empty, P1 / Pz the list lengths
    dat <- build_stan_data(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, x_ref = A$x_ref, allow_orphan_initiations = TRUE)
    expect_equal(dat$P, 0L); expect_equal(dim(dat$X), c(dat$N1, 0L)); expect_length(dat$x_ref, 0L); expect_equal(dat$P1, 8L); expect_equal(dat$Pz, 8L)
  }
  ## the three models assemble the same frames (the model changes nothing at assembly; M1c passes ev$grp to Stan)
  expect_identical(out$M1a$sy, out$M1b$sy); expect_identical(out$M1a$alts, out$M1c$alts); expect_identical(out$M1a$ev, out$M1c$ev); expect_identical(out$M1c$ev, out$M1cr$ev); expect_identical(out$M1c$sy, out$M1cr$sy); expect_true("grp" %in% names(out$M1c$ev))
  ## the refusals at the assembly level
  expect_error(assemble_cell("mid", 3L, "all", CONVENTIONAL = TRUE, data_dir = td), "CONVENTIONAL = TRUE needs DISTANCE = capital, MAJOR = joint")
  expect_error(assemble_cell("mid", 3L, "all", CONVENTIONAL = TRUE, DISTANCE = "capital", MAJOR = "joint", CHANGES = "winsor", data_dir = td), "CONVENTIONAL = TRUE needs")
  expect_error(assemble_cell("mid", 3L, "all", CONVENTIONAL = TRUE, DISTANCE = "capital", MAJOR = "joint", FULL_CONTROLS = TRUE, data_dir = td), "CONVENTIONAL = TRUE needs")
  expect_error(assemble_cell("mid", 3L, "all", CONVENTIONAL = "yes", data_dir = td))
  ## a missing or NA alliance column stops the run (the conventional Z needs it; no other cell does)
  td2 <- write_toy_cell(tempfile("toycell"), alts_edit = function(a) a |> select(-alliance_def_off))
  expect_error(quiet(assemble_cell_spec("mid", 3L, conv, data_dir = td2, model = "M1b")), "alliance_def_off")
  td3 <- write_toy_cell(tempfile("toycell"), alts_edit = function(a) a |> mutate(alliance_def_off = ifelse(ccode1 == 10L & ccode2 == 20L, NA_integer_, alliance_def_off)))
  expect_error(quiet(assemble_cell_spec("mid", 3L, conv, data_dir = td3, model = "M1b")), "alliance_def_off")
  ## (c) the existing cells are unchanged: the main spec, the M1c main spec, and the _dist_major cell assemble exactly as before the switch existed
  M <- quiet(assemble_cell_spec("mid", 3L, main, data_dir = td, model = "M1b")); M0 <- quiet(assemble_cell("mid", 3L, "all", data_dir = td))
  expect_identical(M$W_cols, c("dM_i", "dR_i", "M_lvl", "R_lvl", "ongoing_tm1", "log_prior_mids", "dem_tm1", "log_n_system_tm1")); expect_identical(M$X_cols, c("dM_i", "dR_i"))
  expect_identical(M$Z_cols, c("contig_land", "joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "dM_j_h", "dRpr_j_h", "dch_j_missing"))
  for (nm in c("sy", "alts", "ev", "W_cols", "X_cols", "Z_cols", "x_ref")) expect_identical(M[[nm]], M0[[nm]], info = nm)
  expect_false(M$conventional); expect_false(any(grepl("peace", c(M$W_cols, M$Z_cols))))
  C <- quiet(assemble_cell_spec("mid", 3L, main, data_dir = td, model = "M1c")); expect_identical(C$W_cols, c(M$W_cols, "multi_mid_jan")); expect_identical(C$Z_cols, M$Z_cols); expect_false(C$conventional)
  DM <- quiet(assemble_cell_spec("mid", 3L, parse_spec(c("--distance=capital", "--major=joint")), data_dir = td, model = "M1b"))
  expect_identical(DM$W_cols, c(M$W_cols, "cowmaj_tm1")); expect_identical(DM$Z_cols, c("contig_land", "joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "log_capdist", "dM_j_h", "dRpr_j_h", "dch_j_missing", "one_major", "both_major"))
  expect_false(DM$conventional)
  ## the conventional cell shares the rows, events, weights and groups of the _dist_major cell (same NA-distance drop: none on the toy) and differs in the columns only
  expect_identical(out$M1b$ev, DM$ev); expect_identical(nrow(out$M1b$alts), nrow(DM$alts)); expect_identical(out$M1b$alts$log_capdist, DM$alts$log_capdist); expect_identical(out$M1b$sy$cowmaj_tm1, DM$sy$cowmaj_tm1)
  unlink(c(td, td2, td3), recursive = TRUE)
})

test_that("the three M1 Stan programs take P = 0: a tiny fit of m1_latent_additive (a and b) and m1_context on the conventional toy cell; the existing programs are unchanged elsewhere", {
  sm_add <- skip_if_no_build("m1_latent_additive"); sm_ctx <- skip_if_no_build("m1_context"); sm_rho <- skip_if_no_build("m1_context_rho")
  ## the sources: P >= 0 in the three M1 programs, P >= 1 still in m2_ame.stan (not part of this build)
  for (m in c("m1_latent_additive", "m1_latent_prior", "m1_context", "m1_context_rho")) expect_true(any(grepl("^  int<lower=0> P;", readLines(file.path(root, "stan", paste0(m, ".stan"))))), info = m)
  expect_true(any(grepl("^  int<lower=1> P;", readLines(file.path(root, "stan", "m2_ame.stan")))))
  td <- write_toy_cell(tempfile("toycell"), events_edit = toy_joiners_conv)
  conv <- parse_spec("--comparison=conventional")
  A <- quiet(assemble_cell_spec("mid", 3L, conv, data_dir = td, model = "M1c"))
  for (jt in c(0L, 1L)) {
    dat <- build_stan_data(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = jt, x_ref = A$x_ref, allow_orphan_initiations = TRUE)
    expect_equal(dat$P, 0L)
    fit <- sampling(sm_add, data = dat, chains = 1, iter = 40, warmup = 20, seed = 123, refresh = 0, pars = c("z_re", "z_b", "L_Om_joint", "L_Om", "re"), include = FALSE)
    expect_equal(fit@mode, 0L); expect_equal(fit@par_dims$beta, c(0, 2)); expect_equal(fit@par_dims$alpha, 8); expect_equal(fit@par_dims$gamma, 8)
    ## mu_state = mu0 + a_i exactly (no X): check on one draw
    e <- rstan::extract(fit, pars = c("mu0", "mu_state", "Omega"), permuted = TRUE)
    expect_true(all(is.finite(e$mu_state))); expect_equal(dim(e$mu_state)[2:3], c(dat$S, 2L))
    ## with no X, mu_event equals mu_state of the event's state (both are mu0 + a_i)
    me <- rstan::extract(fit, pars = "mu_event", permuted = TRUE)$mu_event
    for (ev in seq_len(dat$E)) expect_equal(me[, ev, ], e$mu_state[, dat$sy_state[dat$ev_sy[ev]], ], info = paste("event", ev))
  }
  datc <- build_stan_data(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, x_ref = A$x_ref, allow_orphan_initiations = TRUE); datc$G <- 3L; datc$ev_grp <- as.integer(A$ev$grp)
  expect_true(all(0:2 %in% datc$ev_grp))
  for (prog in c("m1_context", "m1_context_rho")) {
    fitc <- sampling(if (prog == "m1_context") sm_ctx else sm_rho, data = datc, chains = 1, iter = 40, warmup = 20, seed = 123, refresh = 0, pars = c("z_re", "z_b", "L_Om", "re"), include = FALSE)
    expect_equal(fitc@mode, 0L, info = prog); expect_equal(fitc@par_dims$beta, c(0, 2)); expect_equal(fitc@par_dims$delta, c(2, 2)); expect_equal(fitc@par_dims$sigma_g, c(3, 2))
    if (prog == "m1_context_rho") { expect_equal(fitc@par_dims$rho_g, 3); expect_false("rho_x" %in% fitc@model_pars) } else expect_null(fitc@par_dims$rho_g)
    ec <- rstan::extract(fitc, pars = c("mu0", "delta", "mu_state", "mu_event"), permuted = TRUE)
    expect_true(all(is.finite(ec$mu_state))); expect_true(all(is.finite(ec$mu_event)))
    ## mu_event carries the group shift on top of mu0 + a_i (= mu_state, the reference ideal point) and nothing from X
    for (ev in seq_len(datc$E)) {
      g <- datc$ev_grp[ev]; shift <- if (g >= 1L) ec$delta[, g, ] else 0 * ec$delta[, 1, ]
      expect_equal(as.numeric(ec$mu_event[, ev, ]), as.numeric(ec$mu_state[, datc$sy_state[datc$ev_sy[ev]], ] + shift), info = paste(prog, "event", ev))
    }
  }
  unlink(td, recursive = TRUE)
})

test_that("02_fit_model.R --comparison=conventional end to end on the toy cell: dry runs and short fits of M1a, M1b, M1c; the refusals; the robustness driver refuses the switch", {
  for (m in c("m1_context", "m1_context_rho", "m1_latent_additive")) if (!file.exists(file.path(root, "stan", "build", paste0(m, ".rds")))) skip(paste("compiled model missing:", m))
  tmp <- tempfile("smoke_conv_"); dir.create(file.path(tmp, "v2"), recursive = TRUE)
  file.copy(file.path(root, "R"), file.path(tmp, "v2"), recursive = TRUE)
  ok <- file.symlink(file.path(root, "stan"), file.path(tmp, "v2", "stan")); if (!ok) file.copy(file.path(root, "stan"), file.path(tmp, "v2"), recursive = TRUE)
  write_toy_cell(file.path(tmp, "v2", "data"), events_edit = toy_joiners_conv)
  dir.create(file.path(tmp, "v2", "output")); dir.create(file.path(tmp, "v2", "logs"))
  v2 <- file.path(tmp, "v2"); out_dir <- file.path(v2, "output")
  run <- function(script, ...) { out <- suppressWarnings(system2("Rscript", c(file.path(v2, "R", script), ...), stdout = TRUE, stderr = TRUE)); list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  ## (a) the dry runs of the three models: the exact column sets, P = 0, the cell names
  prog_of <- function(m) c(M1a = "m1_latent_additive", M1b = "m1_latent_additive", M1c = "m1_context", M1cr = "m1_context_rho")[[m]]
  for (m in c("M1a", "M1b", "M1c", "M1cr")) {
    r0 <- run("02_fit_model.R", "--sample=mid", "--lag=3", paste0("--model=", m), "--comparison=conventional", "--dry_run=TRUE")
    expect_equal(r0$status, 0L, info = paste(m, paste(tail(r0$out, 30), collapse = "\n")))
    sd <- readRDS(file.path(out_dir, sprintf("standata_mid_lag3_all_pr_%s_conv.rds", m)))
    expect_identical(sd$W_cols, W_CONVENTIONAL); expect_identical(sd$Z_cols, Z_CONVENTIONAL); expect_identical(sd$X_cols, character(0))
    expect_equal(sd$dat$P, 0L); expect_equal(sd$dat$P1, 8L); expect_equal(sd$dat$Pz, 8L); expect_equal(dim(sd$dat$X), c(sd$dat$N1, 0L)); expect_true(sd$conventional); expect_equal(sd$spec$comparison, "conventional")
    expect_equal(sd$cell, sprintf("mid_lag3_all_pr_%s_conv", m)); expect_equal(sd$model, m); expect_equal(prog_of(m), sd$stan_program)
    expect_true(any(grepl("conventional comparison (2026-10-31, suffix _conv): W = ", r0$out, fixed = TRUE)))
    expect_true(any(grepl("X = none (P = 0", r0$out, fixed = TRUE)))
    if (m %in% CONTEXT_MODELS) { expect_equal(sd$dat$ev_grp, c(0L, 1L, 2L, 1L, 1L)); expect_equal(sd$dat$G, 3L); expect_true(any(grepl(sprintf("%s stage-2 groups (ev_grp, G = 3", m), r0$out, fixed = TRUE))) } else { expect_null(sd$dat$ev_grp); expect_null(sd$dat$G) }
  }
  ## the refusals: a fixed switch in any value, the other models, a bad value
  for (bad in list(c("--model=M1b", "--controls=minimal"), c("--model=M1b", "--distance=capital"), c("--model=M1c", "--major=joint"), c("--model=M1b", "--changes=winsor"), c("--model=M1b", "--contiguity=land"), c("--model=M1b", "--interact=none"),
                   c("--model=M1m"), c("--model=M2b"))) {
    r <- run("02_fit_model.R", "--sample=mid", "--lag=3", "--comparison=conventional", "--dry_run=TRUE", bad); expect_true(r$status != 0L, info = paste(bad, collapse = " "))
  }
  r <- run("02_fit_model.R", "--sample=mid", "--lag=3", "--model=M1b", "--comparison=yes", "--dry_run=TRUE"); expect_true(r$status != 0L)
  ## (b) short fits of the three models (2 x 60): P = 0 reaches Stan, the fit records the spec, the summary labels alpha / gamma by the lists
  for (m in c("M1a", "M1b", "M1c", "M1cr")) {
    r1 <- run("02_fit_model.R", "--sample=mid", "--lag=3", paste0("--model=", m), "--comparison=conventional", "--chains=2", "--iter=60", "--warmup=30", "--max_treedepth=6")
    expect_equal(r1$status, 0L, info = paste(m, paste(tail(r1$out, 30), collapse = "\n")))
    f <- file.path(out_dir, sprintf("fit_mid_lag3_all_pr_%s_conv.rds", m)); expect_true(file.exists(f))
    o <- readRDS(f)
    expect_equal(o$dat$P, 0L); expect_identical(o$W_cols, W_CONVENTIONAL); expect_identical(o$Z_cols, Z_CONVENTIONAL); expect_identical(o$X_cols, character(0)); expect_length(o$x_ref, 0L)
    expect_equal(o$args$comparison, "conventional"); expect_true(o$args$conventional); expect_equal(o$args$distance, "capital"); expect_equal(o$args$major, "joint"); expect_equal(o$spec$comparison, "conventional")
    expect_equal(o$spec$distance, "none"); expect_equal(o$spec$major, "none")           # the spec's own switch fields stay at their defaults; the spec carries them
    expect_equal(o$fit@par_dims$beta, c(0, 2)); expect_equal(o$args$model, m); expect_equal(o$args$stan_program, prog_of(m))
    expect_equal(o$args$stan_source_sha256, trimws(readLines(file.path(root, "stan", "build", paste0(o$args$stan_program, ".sha256")))[1]))
    s1 <- read.csv(file.path(out_dir, sprintf("summary_mid_lag3_all_pr_%s_conv.csv", m)))
    expect_false(any(grepl("^beta\\[", s1$parameter))); expect_equal(s1$covariate[s1$parameter == "alpha[8]"], "multi_mid_jan"); expect_equal(s1$covariate[s1$parameter == "alpha[6]"], "log_prior_mids")
    expect_equal(s1$covariate[s1$parameter == "gamma[2]"], "log_capdist"); expect_equal(s1$covariate[s1$parameter == "gamma[8]"], "log_dy_prior_mids"); expect_equal(s1$covariate[s1$parameter == "gamma[7]"], "dy_mid_tm1")
    expect_true(all(c("mu0[1]", "mu0[2]") %in% s1$parameter)); if (m == "M1c") expect_true(all(c("delta[1,1]", "delta[2,2]", "sigma_g[3,2]", "rho_x", "mu0_g[3,2]") %in% s1$parameter))
    if (m == "M1cr") { expect_true(all(c("delta[1,1]", "delta[2,2]", "sigma_g[3,2]", "rho_g[1]", "rho_g[3]", "mu0_g[3,2]") %in% s1$parameter)); expect_false("rho_x" %in% s1$parameter); expect_equal(s1$covariate[s1$parameter == "rho_g[2]"], "leader") }
    ## 03_extract.R takes the fit (the typology at x_ref = none is mu0 + a_i)
    e1 <- run("03_extract.R", paste0("--fit=", f)); expect_equal(e1$status, 0L, info = paste(tail(e1$out, 30), collapse = "\n"))
    expect_true(file.exists(file.path(out_dir, sprintf("typology_mid_lag3_all_pr_%s_conv.csv", m))))
    if (m %in% CONTEXT_MODELS) { cm <- read.csv(file.path(out_dir, sprintf("context_model_mid_lag3_all_pr_%s_conv.csv", m))); expect_equal(sum(cm$quantity == "rho_g (group tolerance correlation)"), if (m == "M1cr") 3L else 0L)
      for (g in STAGE2_GROUPS) expect_true(file.exists(file.path(out_dir, sprintf("typology_mid_lag3_all_pr_%s_conv_%s.csv", m, g)))) }
  }
  ## the main-spec M1b fit of the same toy data (P = 2) still runs through the recompiled program and the names never collide
  r2 <- run("02_fit_model.R", "--sample=mid", "--lag=3", "--model=M1b", "--dry_run=TRUE"); expect_equal(r2$status, 0L, info = paste(tail(r2$out, 30), collapse = "\n"))
  sd2 <- readRDS(file.path(out_dir, "standata_mid_lag3_all_pr_M1b.rds")); expect_equal(sd2$dat$P, 2L); expect_false(sd2$conventional); expect_equal(sd2$spec$comparison, "none")
  ## the robustness driver refuses the switch for every preset (a stub 02_fit_model.R would record the calls; none must be made)
  rec <- file.path(tmp, "calls.txt")
  writeLines(c("args <- commandArgs(trailingOnly = TRUE)", sprintf("cat(paste(args, collapse = ' '), '\\n', sep = '', file = '%s', append = TRUE)", rec)), file.path(v2, "R", "02_fit_model.R"))
  for (ps in c("noww2", "dist", "context_b", "major_ix", "perm_error")) { d <- run("05_robustness.R", paste0("--preset=", ps), "--comparison=conventional", "--dry_run=TRUE"); expect_true(d$status != 0L, info = ps); expect_true(any(grepl("comparison=conventional is not an option of 05_robustness.R", d$out, fixed = TRUE)), info = ps) }
  expect_false(file.exists(rec))
  d <- run("05_robustness.R", "--preset=dist", "--dry_run=TRUE"); expect_equal(d$status, 0L, info = paste(tail(d$out, 20), collapse = "\n"))
  calls <- readLines(rec); expect_length(calls, 4L); expect_true(all(grepl("--comparison=none", calls, fixed = TRUE)))   # every preset call passes the switch at its default
  unlink(tmp, recursive = TRUE)
})

test_that("real data (lag 3, when built): the conventional fmid / mid cells assemble for M1a, M1b and M1c with the exact column sets and no NA; the main cells are unchanged", {
  man_path <- file.path(v2_root(), "data", "manifest.csv")
  skip_if_not(file.exists(man_path) && file.exists(file.path(v2_root(), "data", "alts_lag3_all.rds")), "real data files not built")
  conv <- parse_spec("--comparison=conventional")
  for (s in c("fmid", "mid")) {
    A <- quiet(assemble_cell_spec(s, 3L, conv, model = "M1c"))
    expect_identical(A$W_cols, W_CONVENTIONAL); expect_identical(A$Z_cols, Z_CONVENTIONAL); expect_identical(A$X_cols, character(0)); expect_length(A$x_ref, 0L)
    expect_false(anyNA(as.matrix(A$sy[, A$W_cols]))); expect_false(anyNA(as.matrix(A$alts[, A$Z_cols])))
    expect_equal(nrow(A$sy), 14597L); expect_equal(nrow(A$alts), 1720110L); expect_equal(nrow(A$ev), if (s == "fmid") 1055L else 3258L); expect_equal(sum(A$sy$y1), if (s == "fmid") 695L else 2086L)
    expect_true(all(A$group_table$n_usable > 0L))
    B <- quiet(assemble_cell_spec(s, 3L, conv, model = "M1b")); expect_identical(B$W_cols, A$W_cols); expect_identical(B$Z_cols, A$Z_cols); expect_identical(B$ev[, c("sy_row", "alt_row", "w")], A$ev[, c("sy_row", "alt_row", "w")])
  }
  M <- quiet(assemble_cell_spec("fmid", 3L, parse_spec(character(0)), model = "M1b"))
  expect_identical(M$W_cols, c("dM_i", "dR_i", "M_lvl", "R_lvl", "ongoing_tm1", "log_prior_mids", "dem_tm1", "log_n_system_tm1")); expect_identical(M$X_cols, c("dM_i", "dR_i"))
  expect_identical(M$Z_cols, c("contig_land", "joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "dM_j_h", "dRpr_j_h", "dch_j_missing")); expect_false(M$conventional)
})
