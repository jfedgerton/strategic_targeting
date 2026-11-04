## Ideal-point shifters (switch mu_shift = traits, suffix _mst) and the rho-regression model M1bz (stan/m1_latent_rhoreg.stan), both under the
## directed comparison (PI-approved 2026-11-02).
## (a) parsing, suffix, cell names, refusals, the spec strings / names of every existing cell equal to the b43d87f code (sourced from git);
## (b) the exact X (P = 4) of the mu_shift cells for M1a / M1b / M1c / M1cr / M2b / M1bz and the exact Zr of M1bz on the toy cell, x_ref, the
##     assembly otherwise identical to the directed cell; (c) default / _conv / _cdir cells assembled byte-identically to the b43d87f code;
## (d) m1_latent_rhoreg: with g = 0 the log density equals the b program's at rho_x = tanh(g0) up to the two prior terms (accounted for
##     explicitly), the generated log-liks are identical, and with g != 0 every state-year's stage-1 and stage-2 log-liks equal the b program's at
##     rho_x = rho_n of that state-year; the toy cell compiles and samples; (e) 02_fit_model.R end to end (dry runs, the M2b init from the matching
##     M1b with the switch only, an M1bz fit and its labels); 05_robustness.R refuses the switch and the model; (f) the real lag-3 cells when built.
suppressPackageStartupMessages({ library(testthat); library(dplyr); library(rstan) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "assemble_cell.R"))
toy_dirstat <- function(a) a |> mutate(dem1 = as.integer(ccode1 %in% c(10L, 20L)), dem2 = as.integer(ccode2 %in% c(10L, 20L)),
                                       major1 = as.integer(ccode1 %in% c(10L, 40L)), major2 = as.integer(ccode2 %in% c(10L, 40L)))
toy_joiners_dir <- function(e) e |> mutate(orig1 = ifelse(initiator == 30L & target %in% c(20L, 40L), 0L, 1L), is_joiner_pair = orig1 == 0L)
X_EXP <- c("M_lvl", "R_lvl", "cowmaj_tm1", "dem_tm1")
## the b43d87f spec.R + assemble_cell.R in their own environment (helpers_data.R / helpers_model.R come from the global env; assemble_cell()
## does not call label_parameters(), the only helper changed in this commit)
old_code_env <- function() {
  e <- new.env(parent = globalenv())
  for (f in c("spec.R", "assemble_cell.R")) {
    src <- tryCatch(suppressWarnings(system2("git", c("-C", root, "show", paste0("b43d87f:v2/R/", f)), stdout = TRUE, stderr = FALSE)), error = function(err) NULL)
    if (is.null(src) || !length(src) || !is.null(attr(src, "status"))) return(NULL)
    tf <- tempfile(fileext = ".R"); writeLines(src, tf); sys.source(tf, envir = e); unlink(tf)
  }
  e
}

test_that("spec: mu_shift = traits adds _mst right after _cdir, directed only; M1bz is directed only; refusals; existing names and spec strings unchanged", {
  d <- parse_spec("--comparison=directed"); dm <- parse_spec(c("--comparison=directed", "--mu_shift=traits"))
  expect_equal(SPEC_CHOICES$mu_shift, c("none", "traits")); expect_equal(SPEC_DEFAULTS$mu_shift, "none"); expect_equal(SPEC_SUFFIX$mu_shift, c(traits = "_mst"))
  expect_equal(which(names(SPEC_SUFFIX) == "mu_shift"), which(names(SPEC_SUFFIX) == "comparison") + 1L)
  expect_equal(spec_suffix(dm), "_cdir_mst"); expect_false(spec_is_main(dm)); expect_equal(dm$mu_shift, "traits"); expect_equal(d$mu_shift, "none")
  expect_equal(RHOREG_MODELS, "M1bz"); expect_equal(MU_SHIFT_MODELS, c("M1a", "M1b", "M1c", "M1cr", "M2b", "M1bz")); expect_equal(X_MU_TRAITS, X_EXP)
  expect_true("M1bz" %in% ALL_MODELS); expect_false("M1bz" %in% MODELS); expect_false("M1bz" %in% DIRECTED_MODELS)
  expect_equal(stan_program("M1bz"), "m1_latent_rhoreg"); expect_equal(stan_program("M1bz", d), "m1_latent_rhoreg"); expect_true("m1_latent_rhoreg" %in% STAN_PROGRAMS)
  expect_equal(model_D_mode("M1bz"), 3L); expect_equal(model_joint("M1bz"), 1L); expect_equal(model_family("M1bz"), "M1")
  ## the cells of the request (lag 3, directed) and the M2b cell of the switch
  dap <- parse_spec(c("--comparison=directed", "--mu_shift=traits", "--roles=aprime"))
  for (s in c("fmid", "mid")) {
    expect_equal(cell_name(s, 3L, "M1b", dm), sprintf("%s_lag3_all_pr_M1b_cdir_mst", s))
    expect_equal(cell_name(s, 3L, "M1c", dap), sprintf("%s_lag3_all_pr_M1c_cdir_mst_aprime", s))
    expect_equal(cell_name(s, 3L, "M1bz", d), sprintf("%s_lag3_all_pr_M1bz_cdir", s)); expect_equal(cell_name(s, 3L, "M1bz", dm), sprintf("%s_lag3_all_pr_M1bz_cdir_mst", s))
    expect_equal(cell_name(s, 3L, "M2b", dm), sprintf("%s_lag3_all_pr_M2b_K2_cdir_mst", s))
    for (m in c("M1a", "M1cr")) expect_equal(cell_name(s, 3L, m, dm), sprintf("%s_lag3_all_pr_%s_cdir_mst", s, m))
  }
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec(c("--comparison=directed", "--mu_shift=traits", "--ww2=drop", "--weights=none", "--tag=x1"))), "fmid_lag3_all_pr_M1b_noww2_cdir_mst_w1_tx1")
  expect_false(any(duplicated(c(unlist(SPEC_SUFFIX), "_K2", "_pilot", "_STALE", "_samp10", "_noCH", "_polrel"))))
  ## the M2b init cell: the directed M1b WITH the switch (m1_spec() keeps it), for every K; without the switch the plain directed M1b
  for (k in 1:3) expect_equal(m1_cell_for_m2("fmid", 3L, "M2b", parse_spec(c("--comparison=directed", "--mu_shift=traits", sprintf("--ame_K=%d", k)))), "fmid_lag3_all_pr_M1b_cdir_mst")
  expect_equal(m1_cell_for_m2("fmid", 3L, "M2b", d), "fmid_lag3_all_pr_M1b_cdir")
  ## refusals: mu_shift without the directed comparison (default, conventional), a bad value, a bare flag
  expect_error(parse_spec("--mu_shift=traits"), "--mu_shift=traits needs --comparison=directed", fixed = TRUE)
  expect_error(parse_spec(c("--comparison=conventional", "--mu_shift=traits")), "--mu_shift=traits needs --comparison=directed", fixed = TRUE)
  expect_error(parse_spec(c("--comparison=directed", "--mu_shift=levels")), "--mu_shift must be one of none | traits", fixed = TRUE)
  expect_error(parse_spec(c("--comparison=directed", "--mu_shift")), "--mu_shift needs a value", fixed = TRUE)
  expect_true(identical(parse_spec(c("--comparison=directed", "--mu_shift=none"))$mu_shift, "none"))
  ## M1bz: directed only; no roles, no prior switches, no K; the other models with mu_shift are refused by the directed rule already
  for (sp in list(parse_spec(character(0)), parse_spec("--distance=capital")))
    expect_error(cell_name("fmid", 3L, "M1bz", sp), "the model M1bz needs --comparison=directed", fixed = TRUE)
  expect_error(cell_name("fmid", 3L, "M1bz", parse_spec("--comparison=conventional")), "the conventional comparison is specified for M1a, M1b, M1c, M1cr only", fixed = TRUE)
  expect_error(cell_name("fmid", 3L, "M1bz", parse_spec(c("--comparison=directed", "--roles=aprime"))), "the role coding applies to M1c, M1cr, M2c only", fixed = TRUE)
  expect_error(cell_name("fmid", 3L, "M1bz", parse_spec(c("--comparison=directed", "--rho_prior=flat"))), "the prior switches apply to the M1 m / r / a / b programs only")
  expect_error(cell_name("fmid", 3L, "M1bz", parse_spec(c("--comparison=directed", "--sigma_prior=ln05"))), "the prior switches apply to the M1 m / r / a / b programs only")
  expect_error(cell_name("fmid", 3L, "M1bz", parse_spec(c("--comparison=directed", "--ame_K=1"))), "K applies to the M2 family only")
  for (m in c("M1m", "M1r", "M2a", "M2m", "M2r")) expect_error(cell_name("fmid", 3L, m, dm), "the directed comparison is specified for M1a, M1b, M1c, M1cr, M2b only", fixed = TRUE, info = m)
  ## a mu_shift field set by hand on a non-directed spec is refused by cell_name() too
  bad <- parse_spec(character(0)); bad$mu_shift <- "traits"; expect_error(cell_name("fmid", 3L, "M1b", bad), "--mu_shift=traits needs --comparison=directed", fixed = TRUE)
  ## spec strings: mu_shift appears only when on; the permutation part string keeps it when on
  expect_false(grepl("mu_shift", spec_string(d), fixed = TRUE)); expect_true(grepl(" mu_shift=traits ", spec_string(dm), fixed = TRUE))
  expect_false(grepl("mu_shift", perm_part_spec_string(d), fixed = TRUE)); expect_true(grepl(" mu_shift=traits ", perm_part_spec_string(dm), fixed = TRUE))
  ## existing names and spec strings are those of b43d87f
  old <- old_code_env(); if (is.null(old)) skip("git history with commit b43d87f not available")
  combos <- list(character(0), "--comparison=conventional", "--comparison=directed", c("--comparison=directed", "--roles=aprime"), c("--comparison=directed", "--ame_K=1"),
                 c("--comparison=conventional", "--roles=aprime"), c("--distance=capital", "--major=joint"), "--ww2=drop", c("--comparison=directed", "--ww2=drop", "--weights=none", "--tag=x1"),
                 c("--distance=capital", "--contiguity=none"), "--ame_K=1", c("--context=coalition", "--coal_drop=top3"), c("--rho_prior=flat", "--sigma_prior=ln05"))
  for (cb in combos) {
    sn <- parse_spec(cb); so <- old$parse_spec(cb)
    expect_identical(spec_string(sn), old$spec_string(so), info = paste(cb, collapse = " ")); expect_identical(perm_part_spec_string(sn), old$perm_part_spec_string(so), info = paste(cb, collapse = " "))
    expect_identical(spec_suffix(sn), old$spec_suffix(so)); expect_identical(spec_is_main(sn), old$spec_is_main(so))
    for (m in c("M1a", "M1b", "M1c", "M1cr", "M2b", "M1m", "M2a")) {
      nn <- tryCatch(cell_name("fmid", 3L, m, cell_spec(m, sn)), error = function(e) "ERR"); no <- tryCatch(old$cell_name("fmid", 3L, m, old$cell_spec(m, so)), error = function(e) "ERR")
      expect_identical(nn, no, info = paste(m, paste(cb, collapse = " ")))
      expect_identical(tryCatch(stan_program(m, sn), error = function(e) "ERR"), tryCatch(old$stan_program(m, so), error = function(e) "ERR"))
    }
  }
})

test_that("assembly on the toy cell: mu_shift cells for M1a / M1b / M1c / M1cr / M2b / M1bz have X = M_lvl, R_lvl, cowmaj_tm1, dem_tm1 (P = 4) centred as in W; M1bz's Zr is the same four columns; the rest is the directed cell", {
  td <- write_toy_cell(tempfile("toycell"), alts_edit = toy_dirstat, events_edit = toy_joiners_dir)
  d <- parse_spec("--comparison=directed"); dm <- parse_spec(c("--comparison=directed", "--mu_shift=traits"))
  panel <- readRDS(file.path(td, "stage1_panel_lag3.rds"))
  for (m in c("M1a", "M1b", "M1c", "M1cr", "M2b", "M1bz")) {
    A0 <- quiet(assemble_cell_spec("mid", 3L, cell_spec(m, d), data_dir = td, model = m)); A <- quiet(assemble_cell_spec("mid", 3L, cell_spec(m, dm), data_dir = td, model = m))
    expect_identical(A$X_cols, X_EXP, info = m); expect_identical(A$W_cols, W_CONVENTIONAL); expect_identical(A$Z_cols, Z_DIRECTED); expect_true(A$directed)
    expect_identical(A0$X_cols, character(0)); expect_length(A0$x_ref, 0L)
    ## the X values are the W columns as they enter W: the levels within-year centred by the build (untouched here), the indicators grand-centred
    sy <- A$sy; pr <- panel[match(paste(sy$ccode, sy$year), paste(panel$ccode, panel$year)), ]
    expect_equal(sy$M_lvl, pr$M_lvl); expect_equal(sy$R_lvl, pr$R_lvl)
    expect_equal(sy$cowmaj_tm1, pr$cowmaj_tm1 - mean(pr$cowmaj_tm1)); expect_equal(sy$dem_tm1, pr$dem_tm1 - mean(pr$dem_tm1))
    ## x_ref: the mean X over the initiating state-years (the rule of every cell with an X)
    expect_identical(names(A$x_ref), X_EXP); expect_equal(unname(A$x_ref), unname(colMeans(as.matrix(sy[sy$y1 == 1L, X_EXP]))))
    ## everything but X / x_ref is the directed cell
    for (nm in setdiff(names(A0), c("X_cols", "x_ref"))) expect_identical(A[[nm]], A0[[nm]], info = paste(m, nm))
    expect_identical(names(A), names(A0))
    ## the Stan data: P = 4 and X equal to those W columns; W, Z, events identical to the directed cell's
    dat <- build_stan_data(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, x_ref = A$x_ref, K = if (m == "M2b") 2L else NULL, allow_orphan_initiations = TRUE)
    dat0 <- build_stan_data(A0$sy, A0$alts, A0$ev, A0$W_cols, A0$X_cols, A0$Z_cols, joint = 1L, x_ref = A0$x_ref, K = if (m == "M2b") 2L else NULL, allow_orphan_initiations = TRUE)
    expect_equal(dat$P, 4L); expect_equal(dim(dat$X), c(dat$N1, 4L)); expect_identical(dat$X, dat$W[, match(X_EXP, W_CONVENTIONAL)]); expect_length(dat$x_ref, 4L)
    for (nm in setdiff(names(dat0), c("P", "X", "x_ref"))) expect_identical(dat[[nm]], dat0[[nm]], info = paste(m, nm))
    ## M1bz: Zr = the same four columns, with and without the switch
    if (m == "M1bz") for (AA in list(A, A0)) {
      Zr <- rhoreg_Zr(AA)
      expect_equal(dim(Zr), c(nrow(AA$sy), 4L)); expect_null(dimnames(Zr))
      expect_identical(Zr, unname(as.matrix(AA$sy[, X_EXP]))); expect_identical(Zr, unname(dat0$W[, match(X_EXP, W_CONVENTIONAL)]))
    }
  }
  ## refusals: mu_shift without the directed comparison, M1bz outside it, Zr of a non-directed cell or of a doctored one
  expect_error(assemble_cell("mid", 3L, "all", MU_SHIFT = "traits", data_dir = td), "MU_SHIFT = traits needs DIRECTED = TRUE")
  expect_error(assemble_cell("mid", 3L, "all", MU_SHIFT = "traits", CONVENTIONAL = TRUE, DISTANCE = "capital", MAJOR = "joint", data_dir = td), "MU_SHIFT = traits needs DIRECTED = TRUE")
  expect_error(assemble_cell("mid", 3L, "all", MU_SHIFT = "yes", data_dir = td))
  expect_error(assemble_cell_spec("mid", 3L, parse_spec(character(0)), data_dir = td, model = "M1bz"), "defined under --comparison=directed only")
  bad <- parse_spec(character(0)); bad$mu_shift <- "traits"; expect_error(quiet(assemble_cell_spec("mid", 3L, bad, data_dir = td, model = "M1b")))
  Adef <- quiet(assemble_cell_spec("mid", 3L, parse_spec(character(0)), data_dir = td, model = "M1b")); expect_error(rhoreg_Zr(Adef))
  Ad <- quiet(assemble_cell_spec("mid", 3L, d, data_dir = td, model = "M1bz")); Ad2 <- Ad; Ad2$sy$dem_tm1 <- Ad2$sy$dem_tm1 + 1; expect_error(rhoreg_Zr(Ad2))
  unlink(td, recursive = TRUE)
})

test_that("existing cells are unchanged: default, M1c, M2b, _dist_major, conventional and directed cells assemble (and build their Stan data) byte-identically to the b43d87f code", {
  old <- old_code_env(); if (is.null(old)) skip("git history with commit b43d87f not available")
  td <- write_toy_cell(tempfile("toycell"), alts_edit = toy_dirstat, events_edit = toy_joiners_dir)
  cases <- list(list(character(0), "M1b"), list(character(0), "M1c"), list(character(0), "M2b"), list(c("--distance=capital", "--major=joint"), "M1b"), list("--ww2=drop", "M1a"), list("--controls=full", "M1b"),
                list("--comparison=conventional", "M1b"), list("--comparison=conventional", "M1c"), list("--comparison=conventional", "M1cr"),
                list("--comparison=directed", "M1a"), list("--comparison=directed", "M1b"), list("--comparison=directed", "M1c"), list("--comparison=directed", "M1cr"), list(c("--comparison=directed", "--ame_K=1"), "M2b"))
  for (cs in cases) {
    sn <- parse_spec(cs[[1]]); so <- old$parse_spec(cs[[1]])
    An <- quiet(assemble_cell_spec("mid", 3L, cell_spec(cs[[2]], sn), data_dir = td, model = cs[[2]])); Ao <- quiet(old$assemble_cell_spec("mid", 3L, old$cell_spec(cs[[2]], so), data_dir = td, model = cs[[2]]))
    expect_identical(names(An), c(names(Ao), "joiner_set", "ev_sets", "jset_stats", "ev_joiner"), info = paste(cs[[2]], paste(cs[[1]], collapse = " ")))   # joiner sets 2026-11-04
    expect_equal(An$joiner_set, "all"); expect_null(An$ev_sets)
    for (nm in names(Ao)) expect_identical(An[[nm]], Ao[[nm]], info = paste(nm, cs[[2]], paste(cs[[1]], collapse = " ")))
    K <- if (cs[[2]] == "M2b") sn$ame_K else NULL
    expect_identical(build_stan_data(An$sy, An$alts, An$ev, An$W_cols, An$X_cols, An$Z_cols, joint = model_joint(cs[[2]]), x_ref = An$x_ref, K = K, allow_orphan_initiations = TRUE),
                     build_stan_data(Ao$sy, Ao$alts, Ao$ev, Ao$W_cols, Ao$X_cols, Ao$Z_cols, joint = model_joint(cs[[2]]), x_ref = Ao$x_ref, K = K, allow_orphan_initiations = TRUE))
  }
  unlink(td, recursive = TRUE)
})

## ---- the compiled programs ----------------------------------------------------------------------------------------------------------------------
## constrained parameter list of the b program -> the rhoreg program's (rho_x replaced by g0 = atanh(rho_x) and g)
to_rhoreg <- function(p, g) { q <- p[setdiff(names(p), "rho_x")]; q$g0 <- atanh(as.numeric(p$rho_x)); q$g <- array(g, dim = length(g)); q }
fixed_fit <- function(sm, dat, init) sampling(sm, data = dat, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(init), seed = 123, refresh = 0)
pars_of <- function(fit, p) { out <- p[intersect(names(p), c(PARAM_BLOCK, "g0", "g"))]; complete_init(fit, out) }
gq <- function(fit, nm) { a <- rstan::extract(fit, nm)[[nm]]; stopifnot(dim(a)[1] == 1L); as.numeric(a) }   # one Fixed_param draw
## the directed toy cell (mixed choice sets: an IV-only and a stage-2-only row) as M1bz Stan data
toy_rhoreg_data <- function(mu_shift = FALSE) {
  td <- write_toy_cell(tempfile("toycell"), alts_edit = toy_dirstat, events_edit = toy_joiners_dir)
  A <- quiet(assemble_cell_spec("mid", 3L, parse_spec(c("--comparison=directed", if (mu_shift) "--mu_shift=traits")), data_dir = td, model = "M1bz"))
  dat <- build_stan_data(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, x_ref = A$x_ref, allow_orphan_initiations = TRUE)
  Zr <- rhoreg_Zr(A); dat$Q <- ncol(Zr); dat$Zr <- Zr
  unlink(td, recursive = TRUE); dat
}
## gen_init(dat, seed): a b-program init of the data's dimensions and the natural parameters behind it (for stage2_loglik_R())
gen_init <- function(dat, seed = 7) {
  set.seed(seed); S <- dat$S
  L <- t(chol(matrix(c(1, 0.2, 0.1, 0.2, 1, 0.3, 0.1, 0.3, 1), 3))); sd_re <- c(0.5, 1, 0.8)
  z_re <- matrix(rnorm(3 * S), 3, S); z_b <- rnorm(S); sigma_b <- 0.6; sigma_x <- c(0.9, 1.2); mu0 <- c(0.7, -0.4)
  beta <- matrix(rnorm(dat$P * 2, 0, 0.3), dat$P, 2); gamma <- rnorm(dat$Pz, 0, 0.3)
  init <- list(alpha0 = -0.5, alpha = array(rnorm(dat$P1, 0, 0.3), dim = dat$P1), tau = 0.4, mu0 = array(mu0, 2), beta = beta, sigma_x = array(sigma_x, 2),
               rho_x = array(0.4, 1), sd_re = array(sd_re, 3), L_Om_joint = array(L, c(1, 3, 3)), z_re = z_re, sigma_b = sigma_b, z_b = array(z_b, S), gamma = array(gamma, dim = dat$Pz))
  re <- t(diag(sd_re) %*% L %*% z_re)
  nat <- list(mu0 = mu0, beta = beta, a_i = re[, 2:3, drop = FALSE], b_j = sigma_b * (z_b - mean(z_b)), gamma = gamma, sigma_x = sigma_x, rho = 0.4)
  list(init = init, nat = nat)
}
## the simulated tiny cell (P = 4) with a random N1 x 4 Zr
tiny_rhoreg_data <- function() { dat <- tiny_data("M1b", joint = 1L); set.seed(31); dat$Q <- 4L; dat$Zr <- matrix(rnorm(dat$N1 * 4L), dat$N1, 4L); dat }

test_that("m1_latent_rhoreg: the source is the b variant with rho_n = tanh(g0 + Zr[n] g) in the IV and the stage-2 softmax, g0 / g priors in place of rho_x and its Beta prior", {
  src <- readLines(file.path(root, "stan", "m1_latent_rhoreg.stan")); b <- readLines(file.path(root, "stan", "m1_latent_additive.stan"))
  expect_false(any(grepl("rho_x", src[!grepl("^\\s*//", src)])))                     # no rho_x outside comments
  expect_true(any(grepl("^  real g0;", src))); expect_true(any(grepl("^  vector\\[Q\\] g;", src))); expect_true(any(grepl("^  int<lower=0> Q;", src))); expect_true(any(grepl("^  matrix\\[N1, Q\\] Zr;", src)))
  expect_true(any(grepl("^      g0 ~ normal\\(0, 1\\);", src))); expect_true(any(grepl("^      g ~ normal\\(0, 0.5\\);", src))); expect_false(any(grepl("beta_lpdf\\(\\(rho", src)))
  expect_equal(sum(grepl("real rho = tanh(g0 + Zr[n] * g);", src, fixed = TRUE)), 1L)  # one rho per state-year in the partial sum: IV and stage 2 of n
  ## undoing the documented edits gives the b variant's body line for line (header block aside)
  body <- src[(which(grepl("^// ---- the source program's header", src)) + 1L):length(src)]
  body <- body[!grepl("^  // ---- M1bz \\(2026-11-02\\)|^  int<lower=0> Q;|^  matrix\\[N1, Q\\] Zr;|real rho = tanh\\(g0 \\+ Zr\\[n\\] \\* g\\);|^  if \\(joint != 1 \\|\\| D_mode != 3\\) reject|^    real rho0 = tanh\\(g0\\);", body)]
  edits <- c("matrix X, matrix Zr," = "matrix X,", "real g0, vector g, vector b_j" = "real rho, vector b_j", "W, X, Zr, cs_start" = "W, X, cs_start", "sigma_x, g0, g, b_j" = "sigma_x, rho, b_j",
             "tanh(g0 + Zr[n] * g)" = "rho", "tanh(g0 + Zr[ev_sy[ev]] * g)" = "rho", "rho0 * sigma_x" = "rho * sigma_x", "  matrix[D, D] Sigma;         // M1bz: at Zr = 0, rho = tanh(g0)" = "  matrix[D, D] Sigma;")
  for (k in names(edits)) body <- gsub(k, edits[[k]], body, fixed = TRUE)
  ## the parameter / prior / transformed-parameter lines of rho_x come back
  i_g0 <- which(body == "  real g0;                                                       // M1bz: atanh(rho) at Zr = 0")
  expect_length(i_g0, 1L)
  body <- c(body[seq_len(i_g0 - 1L)], "  array[joint * (D - 1)] real<lower=-1, upper=1> rho_x;          // D = 2, joint only", body[(i_g0 + 2L):length(body)])
  i_pr <- which(body == "      // M1bz: the correlation regression on the atanh scale replaces rho_x and its Beta(2.5, 2.5) prior"); expect_length(i_pr, 1L)
  body <- c(body[seq_len(i_pr - 1L)], "      // Sigma's correlation (b only): explicit (rho + 1) / 2 ~ Beta(2.5, 2.5), the same marginal as the shared",
            "      // correlations (2026-10-25 P3; before this the prior was the implicit uniform(-1, 1) of the bounds)", "      target += beta_lpdf((rho_x[1] + 1) / 2 | 2.5, 2.5);", body[(i_pr + 3L):length(body)])
  i_lom <- which(body == "  matrix[Dre, Dre] L_Om;"); body <- append(body, "  real rho;", after = i_lom)
  i_j <- which(body == "      L_Om = L_Om_joint[1];"); body <- append(body, "      rho = rho_x[1];", after = i_j)
  for (pat in c("      L_Om = cholesky_decompose(R_a);", "    L_Om = cholesky_decompose(R_1);")) { i <- which(body == pat); body <- append(body, sub("L_Om = .*$", "rho = 0;", pat), after = i) }
  b_body <- b[(which(grepl("^// Two-stage spatial choice model", b)) - 1L):length(b)]   # the source after its title line
  expect_equal(body, b_body)
})

test_that("m1_latent_rhoreg: with g = 0 the log density is the b program's at rho_x = tanh(g0) minus its Beta prior plus the g0 / g priors; the log-liks are identical", {
  sm_z <- skip_if_no_build("m1_latent_rhoreg"); sm_b <- skip_if_no_build("m1_latent_additive")
  rstan_options(threads_per_chain = 1)
  for (dat in list(tiny_rhoreg_data(), toy_rhoreg_data(), toy_rhoreg_data(mu_shift = TRUE))) {
    expect_equal(dat$joint, 1L); expect_equal(dat$D_mode, 3L)
    p0 <- gen_init(dat)$init
    f_b <- fixed_fit(sm_b, dat, p0); f_z <- fixed_fit(sm_z, dat, to_rhoreg(p0, rep(0, dat$Q)))
    expect_equal(f_b@mode, 0L); expect_equal(f_z@mode, 0L); expect_equal(get_num_upars(f_z), get_num_upars(f_b) + dat$Q)
    set.seed(5)
    for (k in 1:5) {
      u_b <- rnorm(get_num_upars(f_b), 0, 0.6)
      pb <- pars_of(f_b, constrain_pars(f_b, u_b)); pz <- pars_of(f_z, to_rhoreg(pb, rep(0, dat$Q)))
      u_z <- unconstrain_pars(f_z, pz)
      g0 <- pz$g0; rho <- as.numeric(pb$rho_x); expect_equal(tanh(g0), rho, tolerance = 1e-12)
      ## propto log densities without Jacobians: the ~ statements drop their constants (g0 ~ N(0, 1) -> -g0^2 / 2; g ~ N(0, 0.5) -> -2 sum g^2 = 0 here), the
      ## b program's target += beta_lpdf keeps its constant: lp_z - lp_b = -g0^2 / 2 - dbeta((rho + 1) / 2, 2.5, 2.5, log = TRUE)
      d_exp <- -0.5 * g0^2 - dbeta((rho + 1) / 2, 2.5, 2.5, log = TRUE)
      lp_b <- log_prob(f_b, u_b, adjust_transform = FALSE); lp_z <- log_prob(f_z, u_z, adjust_transform = FALSE)
      expect_equal(lp_z - lp_b, d_exp, tolerance = 1e-9)
      ## the log-liks and every shared generated quantity are identical at the point
      gb <- fixed_fit(sm_b, dat, pb); gz <- fixed_fit(sm_z, dat, pz)
      for (nm in c("log_lik_stage1", "log_lik_stage2", "mu_state", "mu_event", "Sigma", "Omega", "R_Om")) expect_equal(gq(gz, nm), gq(gb, nm), tolerance = 1e-12, info = nm)
    }
    ## the g prior term: moving only g moves the log density by -2 sum(g^2) plus the likelihood change; at a point where Zr g only moves rho the prior term
    ## is checked through prior_only = 1 (both likelihoods skipped)
    dp <- dat; dp$prior_only <- 1L
    f_zp <- fixed_fit(sm_z, dp, to_rhoreg(p0, rep(0, dat$Q))); pz0 <- pars_of(f_zp, to_rhoreg(p0, rep(0, dat$Q))); g1 <- c(0.3, -0.7, 0.2, 1.1)[seq_len(dat$Q)]
    pz1 <- pz0; pz1$g <- array(g1, dim = dat$Q)
    expect_equal(log_prob(f_zp, unconstrain_pars(f_zp, pz1), adjust_transform = FALSE) - log_prob(f_zp, unconstrain_pars(f_zp, pz0), adjust_transform = FALSE), -0.5 * sum((g1 / 0.5)^2), tolerance = 1e-10)
  }
})

test_that("m1_latent_rhoreg: with g != 0 each state-year's stage-1 and stage-2 log-liks are the b program's at rho_x = tanh(g0 + Zr[n] g) (rho_n enters the IV and the softmax of n)", {
  sm_z <- skip_if_no_build("m1_latent_rhoreg"); sm_b <- skip_if_no_build("m1_latent_additive")
  rstan_options(threads_per_chain = 1)
  for (dat in list(toy_rhoreg_data(), tiny_rhoreg_data())) {
    gi <- gen_init(dat); p0 <- gi$init
    g <- c(0.4, -0.3, 0.5, -0.6); pz <- to_rhoreg(p0, g); pz$g0 <- 0.25
    gz <- fixed_fit(sm_z, dat, pz); expect_equal(gz@mode, 0L)
    rho_n <- tanh(0.25 + as.numeric(dat$Zr %*% g)); expect_gt(sd(rho_n), 0.05)
    ll1 <- gq(gz, "log_lik_stage1"); ll2 <- gq(gz, "log_lik_stage2")
    ## Sigma is reported at Zr = 0
    sg <- matrix(gq(gz, "Sigma"), 2, 2); expect_equal(sg[1, 2], tanh(0.25) * prod(p0$sigma_x), tolerance = 1e-12)
    ## the toy cell has 4 state-years; the tiny cell 24: check the first 8
    for (n in seq_len(min(dat$N1, 8L))) {
      pb <- p0; pb$rho_x <- array(rho_n[n], 1)
      gb <- fixed_fit(sm_b, dat, pb)
      expect_equal(ll1[n], gq(gb, "log_lik_stage1")[n], tolerance = 1e-12, info = paste("stage 1, state-year", n))
      evn <- which(dat$ev_sy == n)
      if (length(evn)) expect_equal(ll2[evn], gq(gb, "log_lik_stage2")[evn], tolerance = 1e-12, info = paste("stage 2, state-year", n))
    }
    ## and the brute-force stage-2 likelihood per state-year at its own rho
    for (n in unique(dat$ev_sy)[1:min(3, length(unique(dat$ev_sy)))]) {
      pr <- gi$nat; pr$rho <- rho_n[n]
      ref <- stage2_loglik_R(dat, pr)$loglik
      expect_equal(ll2[dat$ev_sy == n], ref[dat$ev_sy == n], tolerance = 1e-8, info = paste("brute force, state-year", n))
    }
  }
  ## the program refuses a non-joint or one-axis data set
  dat <- tiny_rhoreg_data(); p0 <- to_rhoreg(gen_init(dat)$init, rep(0, 4))
  for (bad in list(list(joint = 0L), list(D_mode = 1L))) { db <- modifyList(dat, bad); f <- suppressWarnings(tryCatch(fixed_fit(sm_z, db, p0), error = function(e) NULL)); expect_true(is.null(f) || f@mode != 0L) }
})

test_that("m1_latent_rhoreg: the directed toy cell compiles and samples (with and without mu_shift), Q = 4", {
  sm_z <- skip_if_no_build("m1_latent_rhoreg")
  rstan_options(threads_per_chain = 1)
  for (ms in c(FALSE, TRUE)) {
    dat <- toy_rhoreg_data(mu_shift = ms)
    expect_equal(dat$Q, 4L); expect_equal(dat$P, if (ms) 4L else 0L)
    fit <- sampling(sm_z, data = dat, chains = 1, iter = 200, warmup = 100, seed = 123, refresh = 0, pars = c("z_re", "z_b", "L_Om_joint", "L_Om", "re"), include = FALSE)
    expect_equal(fit@mode, 0L); expect_equal(fit@par_dims$g, 4L); expect_length(fit@par_dims$g0, 0L); expect_false("rho_x" %in% fit@model_pars)
    e <- rstan::extract(fit, pars = c("g0", "g", "log_lik_stage1", "log_lik_stage2", "mu_state", "Sigma"))
    expect_true(all(is.finite(e$g0))); expect_true(all(is.finite(e$log_lik_stage1))); expect_true(all(is.finite(e$log_lik_stage2)))
    expect_equal(as.numeric(e$Sigma[, 1, 2]), as.numeric(tanh(e$g0) * sqrt(e$Sigma[, 1, 1] * e$Sigma[, 2, 2])), tolerance = 1e-10)
    expect_equal(fit@par_dims$beta, c(dat$P, 2L))
  }
})

test_that("02_fit_model.R end to end on the toy cell: mu_shift dry runs for every model, M1bz with and without it, the M2b init from the matching M1b with the switch only, an M1bz fit; 05_robustness.R refuses both", {
  for (m in c("m1_context", "m1_context_rho", "m1_latent_additive", "m2_ame", "m1_latent_rhoreg")) if (!file.exists(file.path(root, "stan", "build", paste0(m, ".rds")))) skip(paste("compiled model missing:", m))
  tmp <- tempfile("smoke_mst_"); dir.create(file.path(tmp, "v2"), recursive = TRUE)
  file.copy(file.path(root, "R"), file.path(tmp, "v2"), recursive = TRUE)
  ok <- file.symlink(file.path(root, "stan"), file.path(tmp, "v2", "stan")); if (!ok) file.copy(file.path(root, "stan"), file.path(tmp, "v2"), recursive = TRUE)
  write_toy_cell(file.path(tmp, "v2", "data"), alts_edit = toy_dirstat, events_edit = toy_joiners_dir)
  dir.create(file.path(tmp, "v2", "output")); dir.create(file.path(tmp, "v2", "logs"))
  v2 <- file.path(tmp, "v2"); out_dir <- file.path(v2, "output")
  run <- function(script, ...) { out <- suppressWarnings(system2("Rscript", c(file.path(v2, "R", script), ...), stdout = TRUE, stderr = TRUE)); list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  prog_of <- c(M1a = "m1_latent_additive", M1b = "m1_latent_additive", M1c = "m1_context", M1cr = "m1_context_rho", M2b = "m2_ame", M1bz = "m1_latent_rhoreg")
  ## (a) dry runs: exact X (P = 4) for every mu_shift model; M1bz with and without the switch carries Zr (Q = 4)
  dry <- list(list("M1a", "--mu_shift=traits", "M1a_cdir_mst"), list("M1b", "--mu_shift=traits", "M1b_cdir_mst"), list("M1c", "--mu_shift=traits", "M1c_cdir_mst"),
              list("M1cr", "--mu_shift=traits", "M1cr_cdir_mst"), list("M2b", c("--mu_shift=traits", "--ame_K=1"), "M2b_K1_cdir_mst"),
              list("M1bz", "--mu_shift=traits", "M1bz_cdir_mst"), list("M1bz", character(0), "M1bz_cdir"))
  for (dd in dry) {
    r0 <- run("02_fit_model.R", "--sample=mid", "--lag=3", paste0("--model=", dd[[1]]), "--comparison=directed", dd[[2]], "--dry_run=TRUE")
    expect_equal(r0$status, 0L, info = paste(dd[[1]], paste(tail(r0$out, 30), collapse = "\n")))
    sd <- readRDS(file.path(out_dir, sprintf("standata_mid_lag3_all_pr_%s.rds", dd[[3]])))
    ms <- "--mu_shift=traits" %in% dd[[2]]
    expect_identical(sd$X_cols, if (ms) X_EXP else character(0)); expect_equal(sd$dat$P, if (ms) 4L else 0L); expect_identical(sd$W_cols, W_CONVENTIONAL); expect_identical(sd$Z_cols, Z_DIRECTED)
    expect_equal(sd$mu_shift, if (ms) "traits" else "none"); expect_equal(sd$cell, sprintf("mid_lag3_all_pr_%s", dd[[3]])); expect_equal(sd$stan_program, prog_of[[dd[[1]]]])
    if (ms) { expect_identical(sd$dat$X, sd$dat$W[, match(X_EXP, W_CONVENTIONAL)]); expect_true(any(grepl("mu_shift = traits (2026-11-02, suffix _mst): X = M_lvl,R_lvl,cowmaj_tm1,dem_tm1 (P = 4", r0$out, fixed = TRUE))) }
    if (dd[[1]] == "M1bz") { expect_equal(sd$dat$Q, 4L); expect_identical(sd$dat$Zr, sd$dat$W[, match(X_EXP, W_CONVENTIONAL)]); expect_identical(sd$Zr_cols, X_EXP) } else { expect_null(sd$dat$Q); expect_null(sd$dat$Zr) }
  }
  ## the default directed dry run carries neither X nor Zr (and its Stan data equal the directed toy's in test_directed.R)
  r0 <- run("02_fit_model.R", "--sample=mid", "--lag=3", "--model=M1b", "--comparison=directed", "--dry_run=TRUE"); expect_equal(r0$status, 0L)
  sd <- readRDS(file.path(out_dir, "standata_mid_lag3_all_pr_M1b_cdir.rds")); expect_equal(sd$dat$P, 0L); expect_null(sd$dat$Zr); expect_equal(sd$mu_shift, "none")
  ## refusals: the switch without the directed comparison, M1bz without it, M1bz with roles / a prior switch / a K, the switch with an M2 model other than M2b
  for (bad in list(c("--model=M1b", "--mu_shift=traits"), c("--model=M1b", "--comparison=conventional", "--mu_shift=traits"), "--model=M1bz", c("--model=M1bz", "--comparison=conventional"),
                   c("--model=M1bz", "--comparison=directed", "--roles=aprime"), c("--model=M1bz", "--comparison=directed", "--rho_prior=flat"), c("--model=M1bz", "--comparison=directed", "--ame_K=1"),
                   c("--model=M2a", "--comparison=directed", "--mu_shift=traits"))) {
    r <- run("02_fit_model.R", "--sample=mid", "--lag=3", "--dry_run=TRUE", bad); expect_true(r$status != 0L, info = paste(bad, collapse = " "))
  }
  ## (b) the directed M1b WITHOUT the switch exists; a mu_shift M2b must not start from it (auto: random fallback, loudly; named: an error)
  r <- run("02_fit_model.R", "--sample=mid", "--lag=3", "--model=M1b", "--comparison=directed", "--chains=2", "--iter=60", "--warmup=30", "--max_treedepth=6"); expect_equal(r$status, 0L, info = paste(tail(r$out, 30), collapse = "\n"))
  expect_true(is.na(locate_m1_fit("mid", 3L, "M2b", parse_spec(c("--comparison=directed", "--mu_shift=traits", "--ame_K=1")), out_dir)))
  r <- run("02_fit_model.R", "--sample=mid", "--lag=3", "--model=M2b", "--comparison=directed", "--mu_shift=traits", "--ame_K=1", "--chains=2", "--iter=60", "--warmup=30", "--max_treedepth=6")
  expect_equal(r$status, 0L, info = paste(tail(r$out, 30), collapse = "\n"))
  o <- readRDS(file.path(out_dir, "fit_mid_lag3_all_pr_M2b_K1_cdir_mst.rds")); expect_equal(o$init$source, "stan_random"); expect_true(any(grepl("fit_mid_lag3_all_pr_M1b_cdir_mst.rds", r$out, fixed = TRUE)))
  expect_equal(o$args$mu_shift, "traits"); expect_identical(o$X_cols, X_EXP); expect_equal(o$fit@par_dims$beta, c(4L, 2L))
  r <- run("02_fit_model.R", "--sample=mid", "--lag=3", "--model=M2b", "--comparison=directed", "--mu_shift=traits", "--ame_K=1", "--chains=2", "--iter=60", "--warmup=30", paste0("--init_from=", file.path(out_dir, "fit_mid_lag3_all_pr_M1b_cdir.rds")))
  expect_true(r$status != 0L)
  ## (c) the directed M1b WITH the switch; the mu_shift M2b starts from it; the plain directed M2b must not take it (named: an error)
  r <- run("02_fit_model.R", "--sample=mid", "--lag=3", "--model=M1b", "--comparison=directed", "--mu_shift=traits", "--chains=2", "--iter=60", "--warmup=30", "--max_treedepth=6"); expect_equal(r$status, 0L, info = paste(tail(r$out, 30), collapse = "\n"))
  o1 <- readRDS(file.path(out_dir, "fit_mid_lag3_all_pr_M1b_cdir_mst.rds")); expect_equal(o1$args$mu_shift, "traits"); expect_identical(o1$X_cols, X_EXP); expect_equal(o1$dat$P, 4L)
  s1 <- read.csv(file.path(out_dir, "summary_mid_lag3_all_pr_M1b_cdir_mst.csv")); expect_equal(s1$covariate[s1$parameter == "beta[3,2]"], "cowmaj_tm1"); expect_equal(s1$axis[s1$parameter == "beta[3,2]"], "relational")
  r <- run("02_fit_model.R", "--sample=mid", "--lag=3", "--model=M2b", "--comparison=directed", "--mu_shift=traits", "--ame_K=2", "--chains=2", "--iter=60", "--warmup=30", "--max_treedepth=6")
  expect_equal(r$status, 0L, info = paste(tail(r$out, 30), collapse = "\n"))
  o2 <- readRDS(file.path(out_dir, "fit_mid_lag3_all_pr_M2b_K2_cdir_mst.rds"))
  expect_equal(o2$init$source, "m1_posterior_means"); expect_equal(o2$init$m1_cell, "mid_lag3_all_pr_M1b_cdir_mst"); expect_equal(o2$dat$P, 4L)
  expect_true(any(grepl("the init fit is the M1 fit of the same directed cell (mid_lag3_all_pr_M1b_cdir_mst)", r$out, fixed = TRUE)))
  r <- run("02_fit_model.R", "--sample=mid", "--lag=3", "--model=M2b", "--comparison=directed", "--ame_K=1", "--chains=2", "--iter=60", "--warmup=30", paste0("--init_from=", file.path(out_dir, "fit_mid_lag3_all_pr_M1b_cdir_mst.rds")))
  expect_true(r$status != 0L)
  ## (d) an M1bz fit (with the switch): the program, g0 / g in the summary with their labels, Zr saved, the hash guard of 03_extract.R's helper passes
  r <- run("02_fit_model.R", "--sample=mid", "--lag=3", "--model=M1bz", "--comparison=directed", "--mu_shift=traits", "--chains=2", "--iter=60", "--warmup=30", "--max_treedepth=6"); expect_equal(r$status, 0L, info = paste(tail(r$out, 30), collapse = "\n"))
  oz <- readRDS(file.path(out_dir, "fit_mid_lag3_all_pr_M1bz_cdir_mst.rds"))
  expect_equal(oz$args$stan_program, "m1_latent_rhoreg"); expect_true(oz$args$rhoreg); expect_identical(oz$Zr_cols, X_EXP); expect_equal(oz$dat$Q, 4L); expect_false("rho_x" %in% oz$fit@model_pars)
  expect_silent(assert_fit_stan_current(oz, root = v2))
  sz <- read.csv(file.path(out_dir, "summary_mid_lag3_all_pr_M1bz_cdir_mst.csv"))
  expect_equal(sz$covariate[match(sprintf("g[%d]", 1:4), sz$parameter)], X_EXP); expect_equal(sz$covariate[sz$parameter == "g0"], "intercept (atanh rho at Zr = 0)"); expect_true(all(c("beta[1,1]", "Sigma[1,2]") %in% sz$parameter))
  ## (e) the robustness driver refuses the switch and the model for every preset; no call is made
  rec <- file.path(tmp, "calls.txt")
  writeLines(c("args <- commandArgs(trailingOnly = TRUE)", sprintf("cat(paste(args, collapse = ' '), '\\n', sep = '', file = '%s', append = TRUE)", rec)), file.path(v2, "R", "02_fit_model.R"))
  for (ps in c("noww2", "dist", "context_b", "major_ix", "perm_error")) {
    dr <- run("05_robustness.R", paste0("--preset=", ps), "--mu_shift=traits", "--dry_run=TRUE"); expect_true(dr$status != 0L, info = ps)
    dr <- run("05_robustness.R", paste0("--preset=", ps), "--comparison=directed", "--mu_shift=traits", "--dry_run=TRUE"); expect_true(dr$status != 0L, info = ps)
  }
  dr <- run("05_robustness.R", "--preset=noww2", "--cells=mid:3:M1bz", "--dry_run=TRUE"); expect_true(dr$status != 0L); expect_true(any(grepl("the model M1bz is not an option of 05_robustness.R", dr$out, fixed = TRUE)))
  dr <- run("05_robustness.R", "--preset=noww2", "--cells=mid:3:M1b,fmid:3:M1bz", "--dry_run=TRUE"); expect_true(dr$status != 0L)
  expect_false(file.exists(rec))
  ## a preset without the switch still passes it at its default
  dr <- run("05_robustness.R", "--preset=noww2", "--cells=mid:3:M1b", "--dry_run=TRUE"); expect_equal(dr$status, 0L)
  expect_true(all(grepl("--mu_shift=none", readLines(rec), fixed = TRUE)))
  unlink(tmp, recursive = TRUE)
})

test_that("real data (lag 3, when built): fmid / mid M1b_cdir_mst and M1bz_cdir assemble with X / Zr = the four W columns and the directed rows", {
  man_path <- file.path(v2_root(), "data", "manifest.csv")
  skip_if_not(file.exists(man_path) && file.exists(file.path(v2_root(), "data", "alts_lag3_all.rds")), "real data files not built")
  for (s in c("fmid", "mid")) {
    A <- quiet(assemble_cell_spec(s, 3L, parse_spec(c("--comparison=directed", "--mu_shift=traits")), model = "M1b"))
    expect_identical(A$X_cols, X_EXP); expect_identical(A$Z_cols, Z_DIRECTED); expect_identical(A$W_cols, W_CONVENTIONAL)
    expect_equal(nrow(A$sy), 14597L); expect_equal(nrow(A$alts), 1720110L); expect_equal(nrow(A$ev), if (s == "fmid") 1055L else 3258L)
    expect_equal(unname(A$x_ref), unname(colMeans(as.matrix(A$sy[A$sy$y1 == 1L, X_EXP]))))
    Zr <- rhoreg_Zr(A); expect_equal(dim(Zr), c(14597L, 4L))
    rm(A); invisible(gc())
  }
})
