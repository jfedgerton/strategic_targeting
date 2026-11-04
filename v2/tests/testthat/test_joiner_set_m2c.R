## The joiner choice-set restriction (switch joiner_set = all | opposing, suffix _jset) and the AME model with initiation types M2c
## (stan/m2_context.stan), both PI-approved 2026-11-04.
## (A) joiner_opposing_sets() on synthetic disputes; the switch, the suffix and the refusals; the assembly on the toy cell; the Stan programs:
##     (a) with every event's set = its stage-2 set (joiner_set = all) m1_context / m1_context_rho give the log density, gradient and log-liks of
##     the acec29f programs (compiled from git) bit for bit; (b) with restricted joiner events each restricted event's log_lik_stage2 is the
##     log-softmax over its allowed rows recomputed in R from the same parameters, every other event is unchanged; the programs' own checks;
##     (d) the real lag-3 cells (skipped without the data): every joiner set non-empty, holding the chosen target, opposite-side states only.
## (B) M2c: (a) with V = 0 the log-liks equal M1c's and the log density differs by exactly the AME priors; restricted sets as in (b) with u'v;
##     (b) refusals; (c) a toy fit recovers the planted sign of a delta; 02_fit_model.R / 03_extract.R end to end on the toy cell.
suppressPackageStartupMessages({ library(testthat); library(dplyr); library(rstan) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "assemble_cell.R")); source(file.path(root, "R", "context_model.R"))

mk_dd <- function(...) bind_rows(...) |> transmute(dispnum = as.integer(dispnum), ccode1 = as.integer(ccode), year = as.integer(year), sidea1 = as.integer(side))
P <- function(d, cc, side, years = 2000L) tibble(dispnum = d, ccode = cc, year = years, side = side)
J <- function(d, i, j, s1, grp = 2L, year = 2000L) tibble(dispnum = as.integer(d), initiator = as.integer(i), target = as.integer(j), year = as.integer(year), sidea1 = as.integer(s1), grp_aprime = as.integer(grp))

test_that("joiner_opposing_sets(): opposite-side participants of the onset year; flips, switchers, absent states; the target must be in the set", {
  dd <- mk_dd(P(1, 10, 1), P(1, 11, 1), P(1, 20, 0), P(1, 21, 0), P(1, 22, 0, 2001),       # 22 enters in 2001: not in the 2000 set
              P(2, 30, 1), P(2, 40, 0), P(2, 41, 0), P(2, 41, 1), P(2, 42, 1))             # 41 on both sides in 2000 (switcher)
  pr <- bind_rows(J(1, 11, 20, 1), J(1, 10, 20, 1, grp = 0L),                              # joiner 11 on side A; 10 is not a joiner (no set)
                  J(2, 40, 30, 0),                                                         # joiner 40 on COW side B (flipped: own side 0)
                  J(1, 22, 10, 0, year = 2001L))                                           # joiner 22 (side B) in 2001: opposite = side A members of 2001 (none but...)
  dd <- bind_rows(dd, mk_dd(P(1, 10, 1, 2001)))
  r <- joiner_opposing_sets(pr, dd)
  s <- split(r$opp_ccode, paste(r$initiator, r$target, r$year))
  expect_equal(sort(s[["11 20 2000"]]), c(20L, 21L))                                       # side B in 2000, never the joiner, never 22 (absent in 2000)
  expect_equal(sort(s[["40 30 2000"]]), c(30L, 41L, 42L))                                  # side A in 2000 incl. the switcher 41
  expect_equal(s[["22 10 2001"]], 10L)                                                     # a singleton set
  expect_false("10 20 2000" %in% names(s))                                                 # non-joiners get no set
  expect_true(r$opp_two_sides[r$initiator == 40L & r$opp_ccode == 41L]); expect_equal(attr(r, "diagnostics")$n_singleton, 1L); expect_equal(attr(r, "diagnostics")$n_events, 3L)
  ## the chosen target not on the opposite side (or absent in the onset year) fails loudly, listing the events
  expect_error(joiner_opposing_sets(J(1, 11, 10, 1), dd), "1 joiner events have a chosen target that is not on the opposite side.*11->10 2000", perl = TRUE)
  expect_error(joiner_opposing_sets(J(1, 11, 22, 1), dd), "not on the opposite side")
  expect_error(joiner_opposing_sets(J(1, 10, 20, 1, grp = 0L), dd), "no joiner event")
})

test_that("spec: joiner_set = all | opposing, suffix _jset right after _aprime, roles = aprime and the role models only; M2c directed + aprime only; strings unchanged at the default", {
  expect_equal(SPEC_CHOICES$joiner_set, c("all", "opposing")); expect_equal(SPEC_SUFFIX$joiner_set, c(opposing = "_jset"))
  expect_equal(which(names(SPEC_SUFFIX) == "joiner_set"), which(names(SPEC_SUFFIX) == "roles") + 1L)
  dj <- parse_spec(c("--comparison=directed", "--roles=aprime", "--joiner_set=opposing"))
  expect_equal(cell_name("fmid", 3L, "M1c", dj), "fmid_lag3_all_pr_M1c_cdir_aprime_jset"); expect_equal(cell_name("mid", 3L, "M1cr", dj), "mid_lag3_all_pr_M1cr_cdir_aprime_jset")
  expect_equal(cell_name("fmid", 3L, "M2c", parse_spec(c("--comparison=directed", "--roles=aprime", "--joiner_set=opposing", "--ame_K=1", "--tag=long"))), "fmid_lag3_all_pr_M2c_K1_cdir_aprime_jset_tlong")
  expect_equal(cell_name("mid", 3L, "M2c", parse_spec(c("--comparison=directed", "--roles=aprime"))), "mid_lag3_all_pr_M2c_K2_cdir_aprime")
  expect_equal(m1_cell_for_m2("fmid", 3L, "M2c", parse_spec(c("--comparison=directed", "--roles=aprime", "--joiner_set=opposing", "--ame_K=1", "--tag=long"))), "fmid_lag3_all_pr_M1c_cdir_aprime_jset_tlong")
  expect_equal(cell_name("fmid", 3L, "M1c", parse_spec(c("--roles=aprime", "--joiner_set=opposing"))), "fmid_lag3_all_pr_M1c_aprime_jset")
  expect_equal(stan_program("M2c"), "m2_context"); expect_true("m2_context" %in% STAN_PROGRAMS); expect_equal(model_family("M2c"), "M2"); expect_equal(model_D_mode("M2c"), 3L); expect_equal(model_joint("M2c"), 1L)
  expect_false("M2c" %in% MODELS); expect_equal(ROLE_MODELS, c("M1c", "M1cr", "M2c"))
  ## refusals
  expect_error(parse_spec("--joiner_set=opposing"), "--joiner_set=opposing needs --roles=aprime", fixed = TRUE)
  expect_error(parse_spec(c("--comparison=directed", "--joiner_set=opposing")), "needs --roles=aprime")
  expect_error(parse_spec("--joiner_set=yes"), "--joiner_set must be one of all | opposing", fixed = TRUE)
  expect_error(parse_spec("--joiner_set"), "needs a value")
  for (m in c("M1b", "M1a", "M2b", "M1bz")) expect_error(cell_name("fmid", 3L, m, dj), info = m)
  expect_error(cell_name("fmid", 3L, "M2c", parse_spec("--roles=aprime")), "needs --comparison=directed and --roles=aprime")
  expect_error(cell_name("fmid", 3L, "M2c", parse_spec("--comparison=directed")), "needs --comparison=directed and --roles=aprime")
  expect_error(cell_name("fmid", 3L, "M2c", parse_spec(c("--comparison=conventional", "--roles=aprime"))))
  expect_error(cell_name("fmid", 3L, "M2c", parse_spec(c("--comparison=directed", "--roles=aprime", "--mu_shift=traits"))), "specified for")
  expect_error(cell_name("fmid", 3L, "M2c", parse_spec(c("--comparison=directed", "--roles=aprime", "--rho_prior=flat"))), "prior switches")
  expect_error(cell_name("fmid", 3L, "M2c", parse_spec(character(0))))
  ## strings: unchanged at the default, the switch shown when on
  d0 <- parse_spec(c("--comparison=directed", "--roles=aprime"))
  expect_false(grepl("joiner_set", spec_string(d0), fixed = TRUE)); expect_false(grepl("joiner_set", perm_part_spec_string(d0), fixed = TRUE))
  expect_true(grepl("joiner_set=opposing", spec_string(dj), fixed = TRUE)); expect_true(grepl("joiner_set=opposing", perm_part_spec_string(dj), fixed = TRUE))
  expect_identical(spec_suffix(d0), "_cdir_aprime"); expect_false(any(duplicated(c(unlist(SPEC_SUFFIX), "_K2", "_pilot", "_STALE"))))
  ## the git version before this change gives the same names and strings for every existing combination
  e <- new.env(parent = globalenv())
  src <- tryCatch(suppressWarnings(system2("git", c("-C", root, "show", "acec29f:v2/R/spec.R"), stdout = TRUE, stderr = FALSE)), error = function(err) NULL)
  if (is.null(src) || !length(src) || !is.null(attr(src, "status"))) skip("git history with commit acec29f not available")
  tf <- tempfile(fileext = ".R"); writeLines(src, tf); sys.source(tf, envir = e); unlink(tf)
  for (cb in list(character(0), c("--comparison=directed", "--roles=aprime"), c("--comparison=directed", "--roles=aprime", "--tag=long"), c("--comparison=conventional", "--roles=aprime"), "--roles=aprime",
                  c("--comparison=directed", "--mu_shift=traits"), c("--distance=capital", "--major=joint"), "--ame_K=1")) {
    sn <- parse_spec(cb); so <- e$parse_spec(cb)
    expect_identical(spec_string(sn), e$spec_string(so), info = paste(cb, collapse = " ")); expect_identical(perm_part_spec_string(sn), e$perm_part_spec_string(so), info = paste(cb, collapse = " "))
    for (m in c("M1b", "M1c", "M1cr", "M2b", "M1bz")) {
      nn <- tryCatch(cell_name("fmid", 3L, m, cell_spec(m, sn)), error = function(err) "ERR"); no <- tryCatch(e$cell_name("fmid", 3L, m, e$cell_spec(m, so)), error = function(err) "ERR")
      expect_identical(nn, no, info = paste(m, paste(cb, collapse = " ")))
    }
  }
})

## toy roles / sets: 30 -> 20 (dispute 2) and 40 -> 20 (dispute 7) are joiners; 30's stage-2 set is {10, 20} -> restricted to {20} (size 1);
## 40's is {10, 20, 30} -> restricted to {20, 30}
toy_roles_j <- tibble(initiator = c(10L, 10L, 30L, 30L, 30L, 40L), target = c(20L, 30L, 20L, 40L, 10L, 20L), year = 1990L, dispnum = c(1L, 4L, 2L, 3L, 6L, 7L), grp_aprime = c(0L, 1L, 2L, 3L, 3L, 2L))
toy_jsets <- tibble(initiator = c(30L, 40L, 40L), target = 20L, year = 1990L, dispnum = c(2L, 7L, 7L), opp_ccode = c(20L, 20L, 30L))
toy_dirstat <- function(a) a |> mutate(dem1 = as.integer(ccode1 %in% c(10L, 20L)), dem2 = as.integer(ccode2 %in% c(10L, 20L)),
                                       major1 = as.integer(ccode1 %in% c(10L, 40L)), major2 = as.integer(ccode2 %in% c(10L, 40L)))

test_that("assembly: joiner_set = opposing reads joiner_sets_aprime.rds and returns the allowed rows of the joiner events; everything else unchanged; refusals", {
  td <- write_toy_cell(tempfile("toycell"), roles = toy_roles_j, jsets = toy_jsets)
  ap <- parse_spec("--roles=aprime"); aj <- parse_spec(c("--roles=aprime", "--joiner_set=opposing"))
  A0 <- quiet(assemble_cell_spec("mid", 3L, ap, data_dir = td, model = "M1c")); A1 <- quiet(assemble_cell_spec("mid", 3L, aj, data_dir = td, model = "M1c"))
  expect_null(A0$ev_sets); expect_equal(A0$joiner_set, "all"); expect_equal(A1$joiner_set, "opposing")
  for (nm in c("sy", "alts", "ev", "W_cols", "Z_cols", "X_cols", "group_table", "x_ref")) expect_identical(A1[[nm]], A0[[nm]], info = nm)
  key <- paste(A1$sy$ccode[A1$ev$sy_row], A1$alts$ccode2[A1$ev$alt_row])
  i30 <- which(key == "30 20"); i40 <- which(key == "40 20")
  expect_equal(A1$ev_sets$alt_row[A1$ev_sets$ev_idx == i30], A1$ev$alt_row[i30])
  expect_equal(sort(A1$alts$ccode2[A1$ev_sets$alt_row[A1$ev_sets$ev_idx == i40]]), c(20L, 30L))
  expect_setequal(unique(A1$ev_sets$ev_idx), c(i30, i40))
  expect_equal(A1$jset_stats$n_joiner_events, 2L); expect_equal(A1$jset_stats$n_singleton, 1L); expect_equal(A1$jset_stats$n_removed, 2L)   # 30: 1 removed, 40: 1 removed
  out <- capture.output(invisible(assemble_cell_spec("mid", 3L, aj, data_dir = td, model = "M1c")))
  expect_true(any(grepl("joiner_set = opposing (2026-11-04): 2 usable joiner events", out, fixed = TRUE))); expect_true(any(grepl("1 events with set size 1", out, fixed = TRUE)))
  ## Stan data: event_alt_sets() gives the full sets without ev_sets and the restricted ones with it
  d0 <- event_alt_sets(build_stan_data(A0$sy, A0$alts, A0$ev, A0$W_cols, A0$X_cols, A0$Z_cols, joint = 1L, x_ref = A0$x_ref, allow_orphan_initiations = TRUE))
  d1 <- event_alt_sets(build_stan_data(A1$sy, A1$alts, A1$ev, A1$W_cols, A1$X_cols, A1$Z_cols, joint = 1L, x_ref = A1$x_ref, allow_orphan_initiations = TRUE), A1$ev_sets)
  s2n <- tabulate(findInterval(which(d0$in_s2 == 1L), d0$cs_start), nbins = d0$N1)
  expect_equal(as.integer(d0$ev_alt_n), s2n[d0$ev_sy]); expect_equal(as.integer(d1$ev_alt_n)[-c(i30, i40)], s2n[d1$ev_sy][-c(i30, i40)]); expect_equal(as.integer(d1$ev_alt_n)[c(i30, i40)], c(1L, 2L))
  ## refusals: no joiner-set file, an event without a set, a chosen target outside its set, the switch without the aprime roles
  td0 <- write_toy_cell(tempfile("toycell"), roles = toy_roles_j)
  expect_error(quiet(assemble_cell_spec("mid", 3L, aj, data_dir = td0, model = "M1c")))
  td2 <- write_toy_cell(tempfile("toycell"), roles = toy_roles_j, jsets = toy_jsets[toy_jsets$initiator == 40L, ])
  expect_error(quiet(assemble_cell_spec("mid", 3L, aj, data_dir = td2, model = "M1c")), "1 joiner events have no row in joiner_sets_aprime.rds (rebuild the data): 30->20 1990", fixed = TRUE)
  td3 <- write_toy_cell(tempfile("toycell"), roles = toy_roles_j, jsets = toy_jsets |> mutate(opp_ccode = ifelse(initiator == 40L & opp_ccode == 20L, 10L, opp_ccode)))
  expect_error(quiet(assemble_cell_spec("mid", 3L, aj, data_dir = td3, model = "M1c")), "the chosen target of 1 joiner events is not in its restricted set: 40->20 1990")
  expect_error(quiet(assemble_cell("mid", 3L, "all", JOINER_SET = "opposing", data_dir = td, CONTEXT_MODEL = TRUE)), "needs ROLES = aprime")
  expect_error(quiet(assemble_cell_spec("mid", 3L, aj, data_dir = td, model = "M1b")), "the restriction applies to")
  expect_error(event_alt_sets(d0, data.frame(ev_idx = i40, alt_row = d0$ev_row[i30])), "outside its state-year's stage-2 set|not in its allowed set")
  unlink(c(td, td0, td2, td3), recursive = TRUE)
})

## ---- the Stan programs ----------------------------------------------------------------------------------------------------------------
tiny_role_data <- function(seed = 123L, K = NULL) {
  set.seed(seed)
  sim <- simulate_targeting("M1b", S = 6L, Tt = 6L, K = 2L, seed = seed, target_rate = 0.6)
  dat <- build_stan_data(sim$sy, sim$alts, sim$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"), joint = 1L, x_ref = rep(0, 4), D_mode = 3L, K = K)
  dat$G <- 4L; dat$ev_grp <- as.integer((seq_len(dat$E) - 1L) %% 4L)
  dat
}
## restricted sets for the group-2 (joiner) events: a random subset of the stage-2 rows holding the chosen row; the first joiner a singleton
restrict_joiners <- function(dat, seed = 11L) {
  set.seed(seed); rows <- list(); first <- TRUE
  for (e in which(dat$ev_grp == 2L)) {
    n <- dat$ev_sy[e]; r <- dat$cs_start[n]:dat$cs_end[n]; r <- r[dat$in_s2[r] == 1L]; oth <- setdiff(r, dat$ev_row[e])
    k <- if (first || length(oth) <= 1L) 0L else sample.int(length(oth) - 1L, 1L)
    keep <- oth[sample.int(length(oth), k)]
    rows[[length(rows) + 1L]] <- data.frame(ev_idx = e, alt_row = c(dat$ev_row[e], keep)); first <- FALSE
  }
  bind_rows(rows)
}
role_pars <- function(dat, rho_g = FALSE) {
  set.seed(7); S <- dat$S
  p <- list(mu0 = c(0.7, -0.4), delta = rbind(c(-0.6, 0.5), c(0.9, -0.8), c(0.2, 0.3)), beta = matrix(rnorm(dat$P * 2, 0, 0.3), dat$P, 2),
            a_i = matrix(rnorm(S * 2), S, 2), c_i = rnorm(S, 0, 0.5), b_j = { z <- rnorm(S); 0.6 * (z - mean(z)) }, gamma = c(0.5, -0.3),
            sigma_g = rbind(c(0.9, 1.2), c(1.6, 0.8), c(1.1, 1.9), c(0.7, 1.4)), rho = if (rho_g) c(0.4, -0.2, 0.6, 0.1) else 0.4, alpha0 = -0.5, alpha = c(0.3, -0.2), tau = 0.4)
  p
}
role_inits <- function(dat, p, rho_g = FALSE) {
  sd_re <- c(0.5, 1, 1); L <- t(chol(matrix(c(1, 0.2, 0.1, 0.2, 1, 0.3, 0.1, 0.3, 1), 3)))
  z_re <- solve(diag(sd_re) %*% L) %*% t(cbind(p$c_i, p$a_i))
  init <- list(alpha0 = p$alpha0, alpha = p$alpha, tau = p$tau, mu0 = p$mu0, delta = p$delta, beta = p$beta, sigma_g = p$sigma_g, sd_re = sd_re, L_Om = L, z_re = z_re,
               sigma_b = 0.6, z_b = p$b_j / 0.6, gamma = p$gamma)
  if (rho_g) init$rho_g <- p$rho else init$rho_x <- p$rho
  init
}
ame_inits <- function(dat, init, zero_v = FALSE) {
  set.seed(5); S <- dat$S; K <- dat$K
  init$u_diag <- as.array(runif(K, 0.3, 0.8)); init$u_lower <- array(rnorm(K * (K - 1) / 2, 0, 0.4), dim = K * (K - 1) / 2)
  init$U_free <- matrix(rnorm((S - K - 1) * K, 0, 0.5), S - K - 1, K); init$V_raw <- if (zero_v) matrix(0, S, K) else matrix(rnorm(S * K, 0, 0.5), S, K)
  init
}
## the U / V of the program from its parameters (anchor block, hard centring), for the R recomputation
ame_UV <- function(dat, init) {
  S <- dat$S; K <- dat$K; U <- matrix(0, S, K); idx <- 1L
  for (k in 1:K) for (l in 1:K) { if (l == k) U[dat$anchor[k], l] <- init$u_diag[k] else if (l < k) { U[dat$anchor[k], l] <- init$u_lower[idx]; idx <- idx + 1L } }
  nona <- setdiff(seq_len(S), dat$anchor); if (S - K > 1) U[nona[seq_len(S - K - 1)], ] <- init$U_free
  U[nona[S - K], ] <- -colSums(U[-nona[S - K], , drop = FALSE])
  list(U = U, V = sweep(init$V_raw, 2, colMeans(init$V_raw)), nonanchor = nona)
}
## R recomputation of every event's log_lik_stage2 over its ALLOWED rows (ev_alt_loc), from the same parameters (+ u_i'v_j when UV is given)
loglik_sets_R <- function(dat, p, UV = NULL) {
  G <- dat$G; rho <- rep(as.numeric(p$rho), length.out = G); shift <- rbind(c(0, 0), matrix(p$delta, ncol = 2))
  Mu <- rep(1, dat$N1) %o% p$mu0 + dat$X %*% p$beta + p$a_i[dat$sy_state, , drop = FALSE]; zg <- as.numeric(dat$Z %*% p$gamma)
  vapply(seq_len(dat$E), function(e) {
    n <- dat$ev_sy[e]; g <- dat$ev_grp[e]; rows <- dat$cs_start[n] - 1L + dat$ev_alt_loc[dat$ev_alt_first[e] + seq_len(dat$ev_alt_n[e]) - 1L]
    Sg <- matrix(c(p$sigma_g[g + 1, 1]^2, rho[g + 1] * prod(p$sigma_g[g + 1, ]), rho[g + 1] * prod(p$sigma_g[g + 1, ]), p$sigma_g[g + 1, 2]^2), 2); Si <- solve(Sg)
    u <- vapply(rows, function(r) { d <- dat$x[r, ] - Mu[n, ] - shift[g + 1, ]; -0.5 * sum(d * (Si %*% d)) + p$b_j[dat$alt_target[r]] + zg[r] +
      if (is.null(UV)) 0 else sum(UV$U[dat$sy_state[n], ] * UV$V[dat$alt_target[r], ]) }, 0)
    u[match(dat$ev_row[e], rows)] - log(sum(exp(u)))
  }, 0)
}
fixed_run <- function(sm, dat, init) suppressWarnings(sampling(sm, data = dat, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(init), seed = 123, refresh = 0))
gq <- function(fit, nm) as.numeric(rstan::extract(fit, nm)[[nm]])          # one draw: the whole array flattened
lp_grad <- function(fit, init, adjust = TRUE) { l <- log_prob(fit, unconstrain_pars(fit, init), adjust_transform = adjust, gradient = TRUE); list(lp = as.numeric(l), grad = as.numeric(attr(l, "gradient"))) }
old_program <- function(prog) {
  src <- tryCatch(suppressWarnings(system2("git", c("-C", root, "show", paste0("acec29f:v2/stan/", prog, ".stan")), stdout = TRUE, stderr = FALSE)), error = function(err) NULL)
  if (is.null(src) || !length(src) || !is.null(attr(src, "status"))) return(NULL)
  cache <- file.path(tempdir(), paste0("acec29f_", prog, ".rds"))
  if (file.exists(cache)) return(readRDS(cache))
  sm <- stan_model(model_code = paste(src, collapse = "\n"), model_name = paste0(prog, "_acec29f"), auto_write = FALSE); saveRDS(sm, cache); sm
}

test_that("(A a) joiner_set = all: m1_context and m1_context_rho give the log density, gradient and log-liks of the acec29f programs bit for bit", {
  for (prog in c("m1_context", "m1_context_rho")) {
    sm <- skip_if_no_build(prog); old <- old_program(prog); if (is.null(old)) skip("git history with commit acec29f not available")
    rg <- prog == "m1_context_rho"
    dat <- event_alt_sets(tiny_role_data()); p <- role_pars(dat, rg); init <- role_inits(dat, p, rg)
    expect_true(all(tabulate(dat$ev_grp + 1L, nbins = 4L) >= 2L))
    dat_old <- dat[setdiff(names(dat), c("n_ev_alt", "ev_alt_first", "ev_alt_n", "ev_alt_loc"))]
    rstan_options(threads_per_chain = 1)
    fn <- fixed_run(sm, dat, init); fo <- fixed_run(old, dat_old, init)
    for (nm in c("log_lik_stage1", "log_lik_stage2", "mu_event", "mu_state", "Sigma_g", "mu0_g")) expect_identical(gq(fn, nm), gq(fo, nm), info = paste(prog, nm))
    for (adj in c(TRUE, FALSE)) { a <- lp_grad(fn, init, adj); b <- lp_grad(fo, init, adj); expect_identical(a$lp, b$lp, info = prog); expect_identical(a$grad, b$grad, info = prog) }
    ## a second random point, and 2 threads with a small grainsize (the reduce_sum slicing)
    init2 <- init; init2$delta[2, ] <- c(-0.3, 1.1); init2$mu0 <- c(0.1, 0.2)
    expect_identical(lp_grad(fn, init2)$lp, lp_grad(fo, init2)$lp, info = prog)
    d2 <- dat; d2$grainsize <- 2L; d2o <- dat_old; d2o$grainsize <- 2L
    rstan_options(threads_per_chain = 2)
    f2 <- fixed_run(sm, d2, init); f2o <- fixed_run(old, d2o, init)
    expect_identical(lp_grad(f2, init)$lp, lp_grad(f2o, init)$lp, info = prog); expect_equal(lp_grad(f2, init)$lp, lp_grad(fn, init)$lp, tolerance = 1e-10)
    rstan_options(threads_per_chain = 1)
  }
})

test_that("(A b) joiner_set = opposing: each restricted joiner event's log_lik_stage2 is the log-softmax over its allowed rows (R recomputation); other events unchanged; the program's checks", {
  for (prog in c("m1_context", "m1_context_rho")) {
    sm <- skip_if_no_build(prog); rg <- prog == "m1_context_rho"
    dat0 <- event_alt_sets(tiny_role_data()); es <- restrict_joiners(dat0); dat <- event_alt_sets(tiny_role_data(), es)
    p <- role_pars(dat, rg); init <- role_inits(dat, p, rg)
    jn <- which(dat$ev_grp == 2L); expect_true(length(jn) >= 2L); expect_true(any(dat$ev_alt_n[jn] == 1L)); expect_true(all(dat$ev_alt_n[jn] < dat0$ev_alt_n[jn]))
    f1 <- fixed_run(sm, dat, init); f0 <- fixed_run(sm, dat0, init)
    ll1 <- gq(f1, "log_lik_stage2"); ll0 <- gq(f0, "log_lik_stage2")
    R1 <- loglik_sets_R(dat, p); R0 <- loglik_sets_R(dat0, p)
    expect_equal(ll1, R1, tolerance = 1e-10); expect_equal(ll0, R0, tolerance = 1e-10)
    expect_equal(R0, stage2_loglik_R_context(dat0, p)$loglik, tolerance = 1e-10)         # the full-set recomputation is the existing brute force
    expect_identical(ll1[-jn], ll0[-jn]); expect_true(all(ll1[jn] >= ll0[jn] - 1e-12)); expect_equal(ll1[jn][dat$ev_alt_n[jn] == 1L], rep(0, sum(dat$ev_alt_n[jn] == 1L)))
    expect_identical(gq(f1, "log_lik_stage1"), gq(f0, "log_lik_stage1"))                  # stage 1 / IV unchanged
    ## the log density moves by exactly the weighted stage-2 difference
    expect_equal(lp_grad(f1, init, FALSE)$lp - lp_grad(f0, init, FALSE)$lp, sum(dat$ev_w * (R1 - R0)), tolerance = 1e-8)
    ## thread invariance with restricted sets
    d2 <- dat; d2$grainsize <- 2L; rstan_options(threads_per_chain = 2)
    expect_equal(lp_grad(fixed_run(sm, d2, init), init)$lp, lp_grad(f1, init)$lp, tolerance = 1e-10); rstan_options(threads_per_chain = 1)
    ## the program rejects malformed sets
    rejected <- function(bad) { f <- tryCatch(suppressWarnings(fixed_run(sm, bad, init)), error = function(e) NULL); is.null(f) || f@mode != 0L }
    e1 <- jn[which(dat$ev_alt_n[jn] >= 2L)[1]]; pos <- dat$ev_alt_first[e1] + 0:1
    bad <- dat; bad$ev_alt_loc[pos] <- rev(bad$ev_alt_loc[pos]); expect_true(rejected(bad), info = "not increasing")
    bad <- dat; nr <- dat$ev_sy[e1]; bad$ev_alt_loc[dat$ev_alt_first[e1]] <- dat$cs_end[nr] - dat$cs_start[nr] + 2L; expect_true(rejected(bad), info = "outside the block")
    e2 <- jn[dat$ev_alt_n[jn] == 1L][1]; bad <- dat; bad$ev_row[e2] <- dat$cs_start[dat$ev_sy[e2]] - 1L + setdiff(dat0$ev_alt_loc[dat0$ev_alt_first[e2] + seq_len(dat0$ev_alt_n[e2]) - 1L], dat$ev_alt_loc[dat$ev_alt_first[e2]])[1]
    expect_true(rejected(bad), info = "chosen row not allowed")
    bad <- dat; bad$ev_alt_first[e1] <- bad$n_ev_alt; expect_true(rejected(bad), info = "index past the end")
  }
})

test_that("(B a) M2c: with V = 0 the log-liks equal M1c's and the log density differs by the AME priors only; with restricted sets and u'v != 0 the R recomputation holds", {
  smc <- skip_if_no_build("m1_context"); sm2 <- skip_if_no_build("m2_context")
  for (K in 1:2) {
    for (restr in c(FALSE, TRUE)) {
      base <- tiny_role_data(K = K); es <- if (restr) restrict_joiners(event_alt_sets(base)) else NULL
      dat <- event_alt_sets(base, es); datc <- dat[setdiff(names(dat), c("K", "anchor"))]
      p <- role_pars(dat); init <- role_inits(dat, p); i2 <- ame_inits(dat, init, zero_v = TRUE)
      fc <- fixed_run(smc, datc, init); f2 <- fixed_run(sm2, dat, i2)
      for (nm in c("log_lik_stage1", "log_lik_stage2", "mu_event", "mu_state", "Sigma_g", "mu0_g", "Sigma", "R_Om")) expect_equal(gq(f2, nm), gq(fc, nm), tolerance = 1e-12, info = paste(K, restr, nm))
      uv <- ame_UV(dat, i2)
      prior_ame <- -0.5 * (sum(i2$u_diag^2) + sum(i2$u_lower^2) + sum(uv$U[uv$nonanchor, ]^2) + sum(i2$V_raw^2))
      expect_equal(lp_grad(f2, i2, FALSE)$lp - lp_grad(fc, init, FALSE)$lp, prior_ame, tolerance = 1e-8, info = paste(K, restr))
      expect_equal(max(abs(rstan::extract(f2, "V")$V)), 0)
      ## u'v != 0: the generated log-liks are the R recomputation over each event's allowed rows with the latent term, stage 1 with the reference IV
      i3 <- ame_inits(dat, init); f3 <- fixed_run(sm2, dat, i3); uv3 <- ame_UV(dat, i3)
      expect_equal(gq(f3, "log_lik_stage2"), loglik_sets_R(dat, p, uv3), tolerance = 1e-10, info = paste(K, restr))
      U3 <- matrix(rstan::extract(f3, "U")$U[1, , ], dat$S, K); expect_equal(U3, uv3$U, tolerance = 1e-12); expect_equal(colSums(U3), rep(0, K), tolerance = 1e-12)
      Mu <- rep(1, dat$N1) %o% p$mu0 + dat$X %*% p$beta + p$a_i[dat$sy_state, ]; zg <- as.numeric(dat$Z %*% p$gamma); Si <- solve(matrix(c(p$sigma_g[1, 1]^2, p$rho * prod(p$sigma_g[1, ]), p$rho * prod(p$sigma_g[1, ]), p$sigma_g[1, 2]^2), 2))
      lse <- vapply(seq_len(dat$N1), function(n) { r <- dat$cs_start[n]:dat$cs_end[n]; r <- r[dat$in_iv[r] == 1L]
        log(sum(exp(vapply(r, function(rr) { d <- dat$x[rr, ] - Mu[n, ]; -0.5 * sum(d * (Si %*% d)) + p$b_j[dat$alt_target[rr]] + zg[rr] + sum(uv3$U[dat$sy_state[n], ] * uv3$V[dat$alt_target[rr], ]) }, 0)))) }, 0)
      eta <- p$alpha0 + as.numeric(dat$W %*% p$alpha) + p$c_i[dat$sy_state] + p$tau * lse
      expect_equal(gq(f3, "log_lik_stage1"), dbinom(dat$y1, 1, plogis(eta), log = TRUE), tolerance = 1e-10)
    }
  }
  ## the program's own checks
  dat <- event_alt_sets(tiny_role_data(K = 2L)); i2 <- ame_inits(dat, role_inits(dat, role_pars(dat)))
  rejected <- function(bad) { f <- tryCatch(suppressWarnings(fixed_run(sm2, bad, i2)), error = function(e) NULL); is.null(f) || f@mode != 0L }
  bad <- dat; bad$D_mode <- 1L; expect_true(rejected(bad)); bad <- dat; bad$joint <- 0L; expect_true(rejected(bad)); bad <- dat; bad$sampled <- 1L; expect_true(rejected(bad))
  bad <- dat; bad$anchor <- array(rep(dat$anchor[1], 2), dim = 2); expect_true(rejected(bad)); bad <- dat; bad$ev_grp[1] <- 4L; expect_true(rejected(bad))
  ## thread invariance
  d2 <- dat; d2$grainsize <- 2L; rstan_options(threads_per_chain = 2)
  expect_equal(lp_grad(fixed_run(sm2, d2, i2), i2)$lp, lp_grad(fixed_run(sm2, dat, i2), i2)$lp, tolerance = 1e-10); rstan_options(threads_per_chain = 1)
})

test_that("(B c) M2c toy fit (2 chains x 400): recovers the planted sign of delta (group 1, material axis)", {
  sm2 <- skip_if_no_build("m2_context")
  set.seed(123); S <- 12L; Tt <- 15L; G <- 2L; K <- 1L
  sy <- expand.grid(state_idx = seq_len(S), t = seq_len(Tt)) |> arrange(t, state_idx) |> mutate(y1 = rbinom(n(), 1, 0.6), W1 = rnorm(n()), W2 = rnorm(n()))
  pos <- array(rnorm(S * Tt * 2), c(S, Tt, 2))
  alts <- bind_rows(lapply(seq_len(nrow(sy)), function(n) { i <- sy$state_idx[n]; tt <- sy$t[n]; j <- setdiff(seq_len(S), i)
    data.frame(sy_row = n, target_idx = j, dM = pos[i, tt, 1] - pos[j, tt, 1], dR = pos[i, tt, 2] - pos[j, tt, 2], Z1 = rnorm(length(j)), Z2 = rnorm(length(j))) }))
  a_i <- matrix(rnorm(S * 2, 0, 0.3), S, 2); b_j <- rnorm(S, 0, 0.3); U <- rnorm(S, 0, 0.3); V <- rnorm(S, 0, 0.3); delta1 <- c(2, 0); gam <- c(0.3, -0.2)
  ev <- list(); grp <- integer(0)
  for (n in which(sy$y1 == 1L)) {
    r <- which(alts$sy_row == n); i <- sy$state_idx[n]; taken <- integer(0)
    for (k in 1:2) {
      g <- rbinom(1, 1, 0.5); mu <- a_i[i, ] + if (g == 1L) delta1 else 0
      u <- -0.5 * ((alts$dM[r] - mu[1])^2 + (alts$dR[r] - mu[2])^2) + b_j[alts$target_idx[r]] + gam[1] * alts$Z1[r] + gam[2] * alts$Z2[r] + U[i] * V[alts$target_idx[r]]
      pr <- exp(u - max(u)); pr[r %in% taken] <- 0
      ch <- r[sample.int(length(r), 1, prob = pr)]; taken <- c(taken, ch); ev[[length(ev) + 1L]] <- data.frame(sy_row = n, alt_row = ch); grp <- c(grp, g)
    }
  }
  ev <- bind_rows(ev)
  dat <- build_stan_data(sy, alts, ev, c("W1", "W2"), character(0), c("Z1", "Z2"), joint = 1L, x_ref = numeric(0), K = K, D_mode = 3L)
  dat$G <- G; dat$ev_grp <- as.integer(grp); dat <- event_alt_sets(dat)
  expect_true(dat$E > 150L)
  rstan_options(threads_per_chain = 1)
  fit <- suppressWarnings(sampling(sm2, data = dat, chains = 2, iter = 400, warmup = 200, seed = 123, init_r = 0.5, refresh = 0, cores = 2, control = list(adapt_delta = 0.9)))
  expect_equal(fit@mode, 0L); d <- rstan::extract(fit, "delta")$delta[, 1, 1]
  expect_gt(mean(d), 0); expect_gt(quantile(d, 0.05), 0)
  expect_lt(abs(mean(rstan::extract(fit, "delta")$delta[, 1, 2])), 1.5)                   # the relational shift (truth 0) is not pushed away
})

test_that("02_fit_model.R / 03_extract.R end to end on the toy cell: M2c refused without its M1c fit, inits from the matching M1c (with and without --joiner_set); refusals", {
  for (m in c("m1_context", "m2_context")) if (!file.exists(file.path(root, "stan", "build", paste0(m, ".rds")))) skip(paste("compiled model missing:", m))
  tmp <- tempfile("smoke_m2c_"); dir.create(file.path(tmp, "v2"), recursive = TRUE)
  file.copy(file.path(root, "R"), file.path(tmp, "v2"), recursive = TRUE)
  ok <- file.symlink(file.path(root, "stan"), file.path(tmp, "v2", "stan")); if (!ok) file.copy(file.path(root, "stan"), file.path(tmp, "v2"), recursive = TRUE)
  write_toy_cell(file.path(tmp, "v2", "data"), alts_edit = toy_dirstat, roles = toy_roles_j, jsets = toy_jsets)
  dir.create(file.path(tmp, "v2", "output")); dir.create(file.path(tmp, "v2", "logs"))
  v2 <- file.path(tmp, "v2"); out_dir <- file.path(v2, "output")
  run <- function(script, ...) { out <- suppressWarnings(system2("Rscript", c(file.path(v2, "R", script), ...), stdout = TRUE, stderr = TRUE)); list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  base <- c("--sample=mid", "--lag=3", "--comparison=directed", "--roles=aprime")
  short <- c("--chains=2", "--iter=60", "--warmup=30", "--max_treedepth=6")
  ## refusals (dry runs): M2c outside directed + aprime, with mu_shift; the joiner switch without aprime or with a non-role model
  for (bad in list(c("--sample=mid", "--lag=3", "--model=M2c", "--roles=aprime"), c("--sample=mid", "--lag=3", "--model=M2c", "--comparison=directed"),
                   c(base, "--model=M2c", "--mu_shift=traits"), c("--sample=mid", "--lag=3", "--model=M1c", "--comparison=directed", "--joiner_set=opposing"),
                   c("--sample=mid", "--lag=3", "--model=M1b", "--comparison=directed", "--roles=aprime", "--joiner_set=opposing"))) {
    r <- run("02_fit_model.R", bad, "--dry_run=TRUE"); expect_true(r$status != 0L, info = paste(bad, collapse = " "))
  }
  ## a dry run of M2c works without the M1c fit (no inits needed); sampling without it is refused
  r <- run("02_fit_model.R", base, "--model=M2c", "--ame_K=1", "--joiner_set=opposing", "--dry_run=TRUE"); expect_equal(r$status, 0L, info = paste(tail(r$out, 20), collapse = "\n"))
  sd <- readRDS(file.path(out_dir, "standata_mid_lag3_all_pr_M2c_K1_cdir_aprime_jset.rds")); expect_equal(sd$stan_program, "m2_context"); expect_equal(sd$dat$K, 1L); expect_equal(sd$dat$G, 4L)
  expect_equal(sd$joiner_set, "opposing"); expect_equal(sd$jset_stats$n_singleton, 1L); expect_true(sd$dat$n_ev_alt < sum(tabulate(findInterval(which(sd$dat$in_s2 == 1L), sd$dat$cs_start), nbins = sd$dat$N1)[sd$dat$ev_sy]))
  r <- run("02_fit_model.R", base, "--model=M2c", "--ame_K=1", short); expect_true(r$status != 0L); expect_true(any(grepl("M2c needs the matching M1c fit fit_mid_lag3_all_pr_M1c_cdir_aprime.rds", r$out, fixed = TRUE)))
  expect_false(file.exists(file.path(out_dir, "fit_mid_lag3_all_pr_M2c_K1_cdir_aprime.rds")))
  for (js in list(character(0), "--joiner_set=opposing")) {
    suf <- if (length(js)) "_jset" else ""
    r1 <- run("02_fit_model.R", base, "--model=M1c", js, short); expect_equal(r1$status, 0L, info = paste(tail(r1$out, 30), collapse = "\n"))
    o1 <- readRDS(file.path(out_dir, sprintf("fit_mid_lag3_all_pr_M1c_cdir_aprime%s.rds", suf)))
    expect_equal(o1$args$joiner_set, if (length(js)) "opposing" else "all"); expect_equal(!is.null(o1$jset_stats), length(js) > 0L)
    if (length(js)) { expect_true("joiner_sets_aprime.rds" %in% o1$manifest$file); expect_true(any(grepl("joiner_set = opposing (2026-11-04): 2 joiner events", r1$out, fixed = TRUE))) }
    else expect_true(any(grepl("joiner_set = all: every event is scored against its state-year's whole stage-2 set", r1$out, fixed = TRUE)))
    for (K in 1:2) {
      r2 <- run("02_fit_model.R", base, "--model=M2c", sprintf("--ame_K=%d", K), js, short); expect_equal(r2$status, 0L, info = paste(K, suf, paste(tail(r2$out, 30), collapse = "\n")))
      cell <- sprintf("mid_lag3_all_pr_M2c_K%d_cdir_aprime%s", K, suf)
      o2 <- readRDS(file.path(out_dir, paste0("fit_", cell, ".rds")))
      expect_equal(o2$init$source, "m1_posterior_means"); expect_equal(o2$init$m1_cell, o1$cell); expect_equal(o2$args$stan_program, "m2_context"); expect_equal(o2$dat$K, K)
      expect_equal(o2$fit@par_dims$delta, c(3L, 2L)); expect_equal(o2$fit@par_dims$sigma_g, c(4L, 2L)); expect_equal(o2$fit@par_dims$U, c(o2$dat$S, K))
      expect_identical(o2$dat$ev_alt_loc, o1$dat$ev_alt_loc); expect_identical(o2$dat$ev_grp, o1$dat$ev_grp); expect_equal(o2$group_labels, APRIME_GROUPS)
      s2 <- read.csv(file.path(out_dir, paste0("summary_", cell, ".csv"))); expect_equal(s2$covariate[s2$parameter == "delta[2,1]"], "joiner"); expect_true("u_diag[1]" %in% s2$parameter)
      e2 <- run("03_extract.R", paste0("--fit=", file.path(out_dir, paste0("fit_", cell, ".rds")))); expect_equal(e2$status, 0L, info = paste(tail(e2$out, 30), collapse = "\n"))
      for (f in c(sprintf("context_model_%s.csv", cell), sprintf("ame_U_%s.csv", cell), sprintf("ame_pairs_%s.csv", cell), sprintf("typology_%s_joiner.csv", cell))) expect_true(file.exists(file.path(out_dir, f)), info = f)
    }
  }
  ## the M1c fit of ANOTHER joiner_set is never an init: naming it explicitly is an error
  r <- run("02_fit_model.R", base, "--model=M2c", "--ame_K=1", "--joiner_set=opposing", short, paste0("--init_from=", file.path(out_dir, "fit_mid_lag3_all_pr_M1c_cdir_aprime.rds")))
  expect_true(r$status != 0L)
  ## 05_robustness.R refuses the switch and the model
  for (a in list(c("--preset=noww2", "--roles=aprime", "--joiner_set=opposing"), c("--preset=noww2", "--cells=mid:3:M2c"))) {
    dr <- run("05_robustness.R", a, "--dry_run=TRUE"); expect_true(dr$status != 0L, info = paste(a, collapse = " "))
  }
  unlink(tmp, recursive = TRUE)
})

test_that("(A d) real data (lag 3, when built): every joiner's restricted set is non-empty, holds the chosen target and only opposite-side states", {
  man_path <- file.path(v2_root(), "data", "manifest.csv")
  skip_if_not(file.exists(man_path) && file.exists(file.path(v2_root(), "data", "joiner_sets_aprime.rds")) && file.exists(file.path(v2_root(), "data", "alts_lag3_all.rds")), "real data files not built")
  skip_if_not(requireNamespace("peacesciencer", quietly = TRUE), "peacesciencer not installed")
  dd <- peacesciencer::cow_mid_dirdisps |> transmute(dispnum = as.integer(dispnum), ccode = as.integer(ccode1), year = as.integer(year), side = as.integer(sidea1)) |> distinct()
  dj <- parse_spec(c("--comparison=directed", "--roles=aprime", "--joiner_set=opposing"))
  for (s in c("fmid", "mid")) {
    A <- quiet(assemble_cell_spec(s, 3L, dj, model = "M1c"))
    jn <- which(A$ev$grp == 2L); expect_true(length(jn) > 0L); expect_setequal(unique(A$ev_sets$ev_idx), jn)
    n_set <- tabulate(A$ev_sets$ev_idx, nbins = nrow(A$ev))[jn]; expect_true(all(n_set >= 1L))
    expect_true(all(paste(jn, A$ev$alt_row[jn]) %in% paste(A$ev_sets$ev_idx, A$ev_sets$alt_row)))               # the chosen target is in every set
    ## opposite side, checked against COW directly: the initiator's side(s) in the onset year vs each allowed state's side(s) in that year, in the
    ## weight-group dispute of the event (events.rds dispnum)
    evf <- readRDS(file.path(v2_root(), "data", "events.rds")) |> filter(!escalation) |> select(initiator, target, year, dispnum)
    k <- data.frame(ev_idx = A$ev_sets$ev_idx, ini = A$sy$ccode[A$ev$sy_row[A$ev_sets$ev_idx]], tgt = A$alts$ccode2[A$ev$alt_row[A$ev_sets$ev_idx]],
                    year = A$sy$year[A$ev$sy_row[A$ev_sets$ev_idx]], opp = A$alts$ccode2[A$ev_sets$alt_row]) |>
      left_join(evf, by = c("ini" = "initiator", "tgt" = "target", "year"))
    expect_false(anyNA(k$dispnum)); expect_true(all(A$alts$in_s2[A$ev_sets$alt_row] == 1L)); expect_true(all(A$alts$sy_row[A$ev_sets$alt_row] == A$ev$sy_row[A$ev_sets$ev_idx]))
    si <- dd |> rename(ini = ccode, side_i = side); so <- dd |> rename(opp = ccode, side_o = side)
    chk <- k |> left_join(si, by = c("dispnum", "year", "ini"), relationship = "many-to-many") |> left_join(so, by = c("dispnum", "year", "opp"), relationship = "many-to-many") |>
      group_by(ev_idx, opp) |> summarise(present = !all(is.na(side_o)), differs = any(side_o != side_i, na.rm = TRUE), .groups = "drop")
    expect_true(all(chk$present), info = s); expect_true(all(chk$differs), info = s); expect_false(any(k$opp == k$ini))
    message(sprintf("real data %s: %d joiner events; restricted set size min %d / median %.1f / max %d; size 1: %d; removed %d", s, length(jn), min(n_set), median(n_set), max(n_set), sum(n_set == 1L), A$jset_stats$n_removed))
    rm(A); invisible(gc())
  }
})
