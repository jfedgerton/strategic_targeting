## Run 1 (built 2026-11-10; README section 14): the MID-propensity opportunity weight q and the offset lambda x ln q in the M1c model.
## (A) the switches opp_weight / opp_outcome / opp_draw: choices, suffixes, positions, names, refusals, the program; every existing name and string unchanged.
## (B) the auxiliary frame: undirected outcomes, hostility >= 3, the history window [t - 3, t - 1] (never t), sea contiguity, A4 by hand (and the January
##     presence on a toy dispute); NO target-strength variable can enter A (scrambling every strength column leaves the frame identical).
## (C) the expanding window: year t is predicted from the years < t ONLY (truncating or corrupting the data from t on changes nothing up to t; a window
##     table that trains on the prediction year is refused), the burn-in is fitted in sample, the decade rule, the fit equals stats::glm.
## (D) the draws: full coefficient vectors with common random numbers, coherent across rows and windows (never per-row draws), reproducible under the seed.
## (E) Stan: m1_context_offset with a zero offset is m1_context to the last bit (log density, gradient, generated quantities); a non-zero offset enters the
##     stage-2 utility and the stage-1 IV exactly as a Z column with coefficient 1 does.
## (F) 02_fit_model.R end to end on the toy cell: lambda = 0 gives the Stan data of the cell without the switch plus a zero offset; lambda and the draw
##     index pick the right column; composite weights unchanged; stale or incomplete weights stop the run; short fits and the mixture script.
## (G) N_eff, the per-cell diagnostics by hand, the mixture table; the choice probabilities under the weights still sum to one.
## (H) slurm: 13_opp_weights.sbatch and submit_opp_weights.sh against stubs. (I) the real builds (skipped without v2/data_dymid).
suppressPackageStartupMessages({ library(testthat); library(dplyr) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "assemble_cell.R")); source(file.path(root, "R", "opp_weights.R"))
spec_o <- function(...) parse_spec(c("--comparison=directed", "--roles=bcj", ...))

## ---- (A) the switches ----------------------------------------------------------------------------------------------------------------------------------
test_that("spec: opp_weight, opp_outcome, opp_draw: choices, suffixes, names, the offset program; existing names and strings unchanged", {
  expect_equal(SPEC_CHOICES$opp_weight, c("none", "0", "0.5", "1")); expect_equal(SPEC_CHOICES$opp_outcome, c("mid", "host3")); expect_equal(SPEC_CHOICES$opp_draw, as.character(0:5))
  expect_equal(unname(OPP_LAMBDA[c("0", "0.5", "1")]), c(0, 0.5, 1)); expect_equal(OPP_MODELS, "M1c"); expect_equal(OPPW_N_DRAWS, 5L); expect_equal(OPPW_BURN_IN, 1865L)
  expect_equal(SPEC_SUFFIX$opp_weight, c("0" = "_opw0", "0.5" = "_opw05", "1" = "_opw1")); expect_equal(SPEC_SUFFIX$opp_outcome, c(host3 = "_oh3"))
  expect_equal(SPEC_SUFFIX$opp_draw, c("1" = "_d1", "2" = "_d2", "3" = "_d3", "4" = "_d4", "5" = "_d5"))
  nm <- names(SPEC_SUFFIX); expect_equal(nm[which(nm == "events") + 1:3], c("opp_weight", "opp_outcome", "opp_draw")); expect_true(which(nm == "opp_draw") < which(nm == "rho_prior"))
  expect_false(any(duplicated(c(unlist(SPEC_SUFFIX), "_K2", "_pilot", "_STALE")))); expect_false(any(grepl("^_t[A-Za-z0-9]", unlist(SPEC_SUFFIX))))
  ## the run-1 cells
  expect_equal(cell_name("fmid", 3L, "M1c", spec_o("--entry=later", "--opp_weight=0", "--tag=long")), "fmid_lag3_all_pr_M1c_cdir_bcj_ent_opw0_tlong")
  expect_equal(cell_name("mid", 3L, "M1c", spec_o("--entry=later", "--opp_weight=0.5", "--tag=long")), "mid_lag3_all_pr_M1c_cdir_bcj_ent_opw05_tlong")
  expect_equal(cell_name("fmid", 3L, "M1c", spec_o("--events=dymid", "--opp_weight=1", "--tag=long")), "fmid_lag3_all_pr_M1c_cdir_bcj_dy_opw1_tlong")
  expect_equal(cell_name("fmid", 3L, "M1c", spec_o("--events=dymid", "--opp_weight=1", "--opp_draw=3", "--tag=long")), "fmid_lag3_all_pr_M1c_cdir_bcj_dy_opw1_d3_tlong")
  expect_equal(cell_name("mid", 3L, "M1c", spec_o("--events=dymid", "--opp_weight=1", "--opp_outcome=host3", "--opp_draw=5")), "mid_lag3_all_pr_M1c_cdir_bcj_dy_opw1_oh3_d5")
  ## the program: the offset copy only when the switch is on (lambda = 0 included); everything else as before
  expect_equal(stan_program("M1c", spec_o("--entry=later", "--opp_weight=0")), "m1_context_offset"); expect_equal(stan_program("M1c", spec_o("--entry=later")), "m1_context")
  expect_equal(stan_program("M1c"), "m1_context"); expect_equal(stan_program("M1c", parse_spec()), "m1_context"); expect_equal(stan_program("M1cr", spec_o("--entry=later")), "m1_context_rho")
  expect_true("m1_context_offset" %in% STAN_PROGRAMS); expect_true(file.exists(file.path(root, "stan", "m1_context_offset.stan")))
  expect_equal(spec_opp_weight(list()), "none"); expect_equal(spec_opp_outcome(list()), "mid"); expect_equal(spec_opp_draw(list()), "0"); expect_false(spec_opp_on(list())); expect_false(spec_opp_on(parse_spec()))
  expect_true(spec_opp_on(spec_o("--entry=later", "--opp_weight=0")))
  expect_true(grepl("opp_weight=1 opp_draw=3", spec_string(spec_o("--events=dymid", "--opp_weight=1", "--opp_draw=3")), fixed = TRUE))
  expect_true(grepl("opp_weight=0.5", perm_part_spec_string(spec_o("--entry=later", "--opp_weight=0.5")), fixed = TRUE))
  expect_equal(oppw_foundation(spec_o("--entry=later")), "data_entry"); expect_equal(oppw_foundation(spec_o("--events=dymid")), "data_dymid")
  expect_equal(basename(oppw_path(spec_o("--events=dymid", "--opp_weight=1", "--opp_outcome=host3"), "x")), "oppw_data_dymid_host3_lag3.rds")
  ## every existing name and string, hard-coded from the code before the change
  expect_identical(spec_string(parse_spec()), "pairs=all context=all coal_ties=inclusive coal_drop=none direction=cowwar ww2=include weights=dispute controls=minimal major=none interact=none distance=none contiguity=land changes=raw rho_prior=beta25 sigma_prior=ln1 comparison=none roles=stageb ame_K=2 tag=")
  expect_identical(spec_string(spec_o("--entry=later", "--tag=long")),
                   "pairs=all context=all coal_ties=inclusive coal_drop=none direction=cowwar ww2=include weights=dispute controls=minimal major=none interact=none distance=none contiguity=land changes=raw rho_prior=beta25 sigma_prior=ln1 comparison=directed roles=bcj entry=later ame_K=2 tag=long")
  expect_identical(perm_part_spec_string(parse_spec()), "pairs=all context=all coal_ties=inclusive coal_drop=none direction=cowwar ww2=include weights=dispute controls=minimal major=none interact=none distance=none contiguity=land changes=raw ame_K=2")
  expect_identical(spec_suffix(parse_spec()), ""); expect_true(spec_is_main(parse_spec()))
  expect_equal(cell_name("fmid", 3L, "M1c", spec_o("--entry=later", "--tag=long")), "fmid_lag3_all_pr_M1c_cdir_bcj_ent_tlong")
  expect_equal(cell_name("fmid", 3L, "M1c", spec_o("--events=dymid", "--tag=long")), "fmid_lag3_all_pr_M1c_cdir_bcj_dy_tlong")
  expect_equal(cell_name("fmid", 3L, "M1c", spec_o("--entry=later", "--reach=capdist", "--tag=long")), "fmid_lag3_all_pr_M1c_cdir_reach_bcj_ent_tlong")
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec()), "fmid_lag3_all_pr_M1b"); expect_equal(cell_name("war", 3L, "M2b", parse_spec()), "war_lag3_all_pr_M2b_K2")
})
test_that("spec: run-1 refusals (other specifications, other models, WWII dropped, a draw without lambda = 1, an outcome without a weight)", {
  for (bad in list(c("--opp_weight=1"), c("--comparison=directed", "--entry=later", "--roles=bc", "--opp_weight=1"), c("--comparison=directed", "--roles=aprime", "--opp_weight=1"),
                   c("--comparison=directed", "--roles=bcj", "--entry=later", "--opp_weight=1", "--ww2=drop"), c("--comparison=directed", "--roles=bcj", "--entry=later", "--opp_weight=1", "--variant=polrel"),
                   c("--comparison=directed", "--roles=bcj", "--entry=later", "--opp_weight=1", "--variant=mid3"), c("--comparison=directed", "--roles=bcj", "--entry=later", "--opp_weight=1", "--joiner_set=opposing"),
                   c("--comparison=directed", "--roles=bcj", "--entry=later", "--opp_weight=1", "--reach=capdist"), c("--comparison=directed", "--roles=bcj", "--entry=later", "--opp_weight=1", "--subsample=prtarget"),
                   c("--comparison=directed", "--roles=bcj", "--entry=later", "--opp_weight=1", "--mu_shift=traits"), c("--comparison=directed", "--roles=bcj", "--entry=later", "--opp_weight=1", "--pairs=joiners"),
                   c("--comparison=directed", "--roles=bcj", "--entry=later", "--opp_weight=0.5", "--opp_draw=2"), c("--comparison=directed", "--roles=bcj", "--entry=later", "--opp_weight=0", "--opp_draw=1"),
                   c("--comparison=directed", "--roles=bcj", "--entry=later", "--opp_draw=1"), c("--comparison=directed", "--roles=bcj", "--entry=later", "--opp_outcome=host3"),
                   c("--comparison=directed", "--roles=bcj", "--entry=later", "--opp_weight=2"), c("--comparison=directed", "--roles=bcj", "--entry=later", "--opp_weight=05"),
                   c("--comparison=directed", "--roles=bcj", "--entry=later", "--opp_weight=1", "--opp_draw=6"), c("--comparison=directed", "--roles=bcj", "--entry=later", "--opp_weight")))
    expect_error(parse_spec(bad), info = paste(bad, collapse = " "))
  expect_error(parse_spec(c("--comparison=directed", "--roles=bcj", "--entry=later", "--opp_weight=1", "--ww2=drop")), "Second World War stays in")
  ## no model but M1c takes the offset (no AME, no pooled model, no role-specific rho)
  s1 <- spec_o("--entry=later", "--opp_weight=1")
  for (m in c("M2c", "M1cr", "M1b", "M2b", "M1a")) expect_error(cell_name("fmid", 3L, m, s1), info = m)
  expect_error(cell_name("fmid", 3L, "M2c", s1), "no AME model takes the offset")
})

## ---- toy fixtures of the auxiliary model ------------------------------------------------------------------------------------------------------------
## a synthetic frame: 7 states, every directed dyad, 1850-1875, a logit outcome; a few stage-2-only rows (in_iv = 0)
toy_frame <- function(seed = 1L, years = 1850:1875) {
  set.seed(seed); st <- seq(10L, 70L, by = 10L)
  fr <- expand.grid(ccode2 = st, ccode1 = st, year = years) |> filter(ccode1 != ccode2) |> arrange(year, ccode1, ccode2)
  n <- nrow(fr); dy <- paste(pmin(fr$ccode1, fr$ccode2), pmax(fr$ccode1, fr$ccode2))
  ld <- stats::setNames(rnorm(length(unique(dy)), 6, 1), unique(dy)); cl <- stats::setNames(rbinom(length(unique(dy)), 1, 0.3), unique(dy))
  lvl <- stats::setNames(rnorm(length(st) * length(years)), paste(rep(st, length(years)), rep(years, each = length(st))))
  fr$in_iv <- 1L; fr$in_s2 <- 1L; fr$in_iv[sample(n, 15L)] <- 0L
  fr$contig_land <- as.integer(cl[dy]); fr$contig_sea <- as.integer((1L - fr$contig_land) * rbinom(n, 1, 0.3)); fr$log_capdist <- as.numeric(ld[dy])
  fr$M_lvl_i <- as.numeric(lvl[paste(fr$ccode1, fr$year)]); fr$dymid_hist3 <- rbinom(n, 1, 0.15); fr$a4_joinable <- rbinom(n, 1, 0.2); fr$decade <- 10L * (fr$year %/% 10L)
  eta <- -2.2 + 0.9 * fr$contig_land + 0.4 * fr$contig_sea - 0.5 * (fr$log_capdist - 6) + 0.2 * (fr$log_capdist - 6) * fr$M_lvl_i + 0.8 * fr$dymid_hist3 + 0.5 * fr$a4_joinable
  fr$y_mid <- rbinom(n, 1, stats::plogis(eta)); fr$y_host3 <- fr$y_mid * rbinom(n, 1, 0.8)
  fr
}
TOY_BURN <- 1860L

## ---- (B) the frame and A4 ------------------------------------------------------------------------------------------------------------------------------
test_that("oppw_frame: undirected outcomes, hostility >= 3, history in [t - 3, t - 1] and never t, sea contiguity, A4; no target-strength variable can enter A", {
  st <- c(10L, 20L, 30L, 40L); yrs <- 1990:1996
  alts <- expand.grid(ccode2 = st, ccode1 = st, year = yrs) |> filter(ccode1 != ccode2)
  set.seed(3); n <- nrow(alts)
  alts <- alts |> mutate(in_iv = 1L, in_s2 = 1L, contig_land = as.integer(ccode1 + ccode2 == 30L), contig_any = as.integer((ccode1 + ccode2) %in% c(30L, 50L)), log_capdist = log1p(runif(n, 100, 5000)),
                         M_j = rnorm(n), Rpr_j = rnorm(n), major2 = rbinom(n, 1, 0.5), major_count = rbinom(n, 2, 0.5), polrel = rbinom(n, 1, 0.5), dM = rnorm(n), dR_pr = rnorm(n), dM_j_h = rnorm(n), dRpr_j_h = rnorm(n))
  alts$in_iv[alts$ccode1 == 10L & alts$ccode2 == 40L & alts$year == 1994L] <- 0L
  alts <- alts[sample(n), ]                                          # the frame sorts by year itself
  panel <- expand.grid(ccode = st, year = yrs) |> mutate(M_lvl = ccode / 100 + (year - 1990) / 10, cowmaj_tm1 = 1L)
  ## onsets: 10-20 in 1992 (hostility 2 / 1: below 3) and 20 -> 10 in 1995 (4 / 1); 30 -> 40 in 1990 (5 / 5)
  ev <- data.frame(initiator = c(10L, 20L, 30L), target = c(20L, 10L, 40L), year = c(1992L, 1995L, 1990L), hostlev1 = c(2L, 4L, 5L), hostlev2 = c(1L, 1L, 5L))
  ## dispute 1 at 1 January 1992: 20, 30, 40 present; 10 not a party. In 1993 state 10 is a party as well.
  jp <- list(present = data.frame(dispnum = 1L, ccode = c(20L, 30L, 40L, 20L, 30L, 40L), year = rep(c(1992L, 1993L), each = 3L)),
             already_in = data.frame(dispnum = 1L, ccode = c(20L, 30L, 40L, 10L, 20L, 30L, 40L), year = c(rep(1992L, 3L), rep(1993L, 4L))),
             flags = panel |> transmute(ccode, year, multi_mid_jan = as.integer(ccode == 10L & year == 1992L)))
  fr <- oppw_frame(alts, panel, ev, jp)
  expect_equal(nrow(fr), n); expect_false(is.unsorted(fr$year)); expect_false(any(OPPW_TARGET_STRENGTH %in% names(fr)))
  k <- function(i, j, t) which(fr$ccode1 == i & fr$ccode2 == j & fr$year == t)
  ## outcome: both directions of the pair in the onset year, no other year
  expect_equal(sum(fr$y_mid), 6L); for (r in list(c(10, 20, 1992), c(20, 10, 1992), c(10, 20, 1995), c(20, 10, 1995), c(30, 40, 1990), c(40, 30, 1990))) expect_equal(fr$y_mid[k(r[1], r[2], r[3])], 1L)
  expect_equal(sum(fr$y_host3), 4L); expect_equal(fr$y_host3[k(10, 20, 1992)], 0L); expect_equal(fr$y_host3[k(10, 20, 1995)], 1L)
  ## history: the three years after an onset, never the onset year, never t - 4
  expect_equal(fr$dymid_hist3[c(k(10, 20, 1992), k(10, 20, 1993), k(20, 10, 1994), k(10, 20, 1995), k(10, 20, 1996))], c(0L, 1L, 1L, 1L, 1L))     # 1996: the 1995 onset
  expect_equal(fr$dymid_hist3[c(k(30, 40, 1990), k(30, 40, 1991), k(40, 30, 1993), k(30, 40, 1994))], c(0L, 1L, 1L, 0L))
  expect_equal(sum(fr$dymid_hist3), 2L * (3L + 3L + 1L)); expect_false(any(fr$dymid_hist3[fr$y_mid == 1L & fr$year %in% c(1990L, 1992L)] == 1L))
  ## A4: only the rows of the non-party 10 towards the present 20, 30, 40 in 1992
  expect_equal(sum(fr$a4_joinable), 3L); expect_equal(fr$a4_joinable[c(k(10, 20, 1992), k(10, 30, 1992), k(10, 40, 1992))], c(1L, 1L, 1L)); expect_equal(fr$a4_joinable[k(20, 30, 1992)], 0L)
  expect_equal(fr$a4_joinable[k(10, 20, 1993)], 0L)
  expect_equal(fr$contig_sea, as.integer((fr$ccode1 + fr$ccode2 == 50L))); expect_equal(fr$M_lvl_i, fr$ccode1 / 100 + (fr$year - 1990) / 10); expect_equal(fr$in_iv[k(10, 40, 1994)], 0L)
  ## no target-strength variable reaches the frame: scrambling (or deleting) every such column of the inputs leaves it identical
  a2 <- alts; for (cc in intersect(OPPW_TARGET_STRENGTH, names(a2))) a2[[cc]] <- rev(a2[[cc]]) * 7 + 1
  expect_identical(oppw_frame(a2, panel, ev, jp), fr)
  expect_identical(oppw_frame(alts[, OPPW_ALTS_COLS], panel[, OPPW_PANEL_COLS], ev, jp), fr)
  p2 <- panel; p2$cowmaj_tm1 <- 0L; expect_identical(oppw_frame(alts, p2, ev, jp), fr)
  expect_false(any(OPPW_TARGET_STRENGTH %in% c(OPPW_A_COLS, OPPW_A_BASELINE, OPPW_ALTS_COLS, OPPW_PANEL_COLS)))
  expect_true(all(c("M_j", "Rpr_j", "major2", "tgt_major", "tgt_major_x_init_major", "dM", "dR_pr", "polrel") %in% OPPW_TARGET_STRENGTH))
  expect_equal(OPPW_A_COLS, c("contig_land", "contig_sea", "log_capdist", "log_capdist_x_init_M", "dymid_hist3", "a4_joinable"))
  dl <- oppw_decade_levels(fr$year, fr$y_mid)
  for (bad in c("M_j", "major2", "tgt_major", "dM")) expect_error(oppw_design(cbind(fr, M_j = 1, major2 = 1, tgt_major = 1, dM = 1), c("contig_land", bad), dl), info = bad)
  ## a row without a capital distance is dropped and counted; an initiator state-year outside the panel stops
  a3 <- alts; a3$log_capdist[1:2] <- NA; f3 <- oppw_frame(a3, panel, ev, jp); expect_equal(nrow(f3), n - 2L); expect_equal(attr(f3, "n_nodist"), 2L)
  expect_error(oppw_frame(alts, panel[panel$ccode != 10L, ], ev, jp), "not in the stage-1 panel")
  ## A4 consistency: a joinable row whose initiator has multi_mid_jan = 0 stops
  jp_bad <- jp; jp_bad$flags$multi_mid_jan <- 0L; expect_error(oppw_frame(alts, panel, ev, jp_bad), "multi_mid_jan = 0")
})
test_that("oppw_jan_presence / oppw_a4 on a toy dispute: >= 3 participants PRESENT at 1 January, a coarse entry in t does not count, a party is never a joiner", {
  ## dispute 1: 10 (side A) v 20 (side B) from June 1990; 30 joins B in March 1991; 40 joins A day-dated 1 January 1992; 50 joins A in February 1992 (month
  ## resolution); all to 1993 (year-resolution exits). Dispute 2: 60 v 70, 1990-1992, never multilateral.
  mk <- function(dispnum, ccode, res, y, m, d, last) { e <- date_keys(res, y, m, d); x <- date_keys(3L, last, NA, NA)
    data.frame(dispnum = dispnum, ccode = ccode, spell = 1L, spell_first = y, spell_last = last, entry_res = e$res, entry_y = e$y, entry_mkey = e$mkey, entry_dkey = e$dkey,
               exit_res = x$res, exit_y = x$y, exit_mkey = x$mkey, exit_dkey = x$dkey) }
  psd <- rbind(mk(1L, 10L, 2L, 1990L, 6L, NA, 1993L), mk(1L, 20L, 2L, 1990L, 6L, NA, 1993L), mk(1L, 30L, 2L, 1991L, 3L, NA, 1993L), mk(1L, 40L, 1L, 1992L, 1L, 1L, 1993L),
               mk(1L, 50L, 2L, 1992L, 2L, NA, 1993L), mk(2L, 60L, 3L, 1990L, NA, NA, 1992L), mk(2L, 70L, 3L, 1990L, NA, NA, 1992L))
  side <- c("10" = 1L, "20" = 0L, "30" = 0L, "40" = 1L, "50" = 1L, "60" = 1L, "70" = 0L)
  dd <- do.call(rbind, lapply(seq_len(nrow(psd)), function(r) data.frame(dispnum = psd$dispnum[r], ccode1 = psd$ccode[r], year = psd$spell_first[r]:psd$spell_last[r], sidea1 = side[[as.character(psd$ccode[r])]])))
  sy <- expand.grid(ccode = seq(10L, 70L, by = 10L), year = 1990:1993)
  jp <- oppw_jan_presence(psd, dd, sy)
  expect_equal(jp$multi$year, c(1992L, 1993L)); expect_equal(jp$multi$n_side_a + jp$multi$n_side_b, c(4L, 5L))       # January 1991: two present, not multilateral
  expect_setequal(jp$present$ccode[jp$present$year == 1992L], c(10L, 20L, 30L, 40L)); expect_setequal(jp$present$ccode[jp$present$year == 1993L], c(10L, 20L, 30L, 40L, 50L))
  keys <- expand.grid(ccode2 = seq(10L, 70L, by = 10L), ccode1 = seq(10L, 70L, by = 10L), year = 1990:1993) |> filter(ccode1 != ccode2)
  a4 <- oppw_a4(keys, jp); on <- keys[a4 == 1L, ]
  expect_true(all(on$year %in% 1992:1993)); expect_equal(sum(a4), 3L * 4L + 2L * 5L)                                 # 1992: 50, 60, 70 x four present; 1993: 60, 70 x five
  expect_setequal(unique(on$ccode1[on$year == 1992L]), c(50L, 60L, 70L)); expect_setequal(unique(on$ccode2[on$year == 1992L]), c(10L, 20L, 30L, 40L))   # 50 entered in February 1992: not yet a party, and not present
  expect_setequal(unique(on$ccode1[on$year == 1993L]), c(60L, 70L)); expect_true(50L %in% on$ccode2[on$year == 1993L])
  expect_false(any(a4[keys$year %in% 1990:1991] == 1L)); expect_false(any(on$ccode1 %in% c(10L, 20L, 30L, 40L)))
})

## ---- (C) the expanding window ------------------------------------------------------------------------------------------------------------------------
test_that("oppw_windows / oppw_decade_levels: the burn-in window is in sample, every later year is predicted from the years before it; the decade rule", {
  w <- oppw_windows(1850:1875, TOY_BURN)
  expect_equal(w$train_max, c(1860L, 1861:1874)); expect_equal(w$pred_from, c(1850L, 1862:1875)); expect_equal(w$pred_to, c(1861L, 1862:1875))
  expect_true(all(w$train_max[-1] < w$pred_from[-1])); expect_true(w$pred_to[1] > w$train_max[1]); expect_true(w$pred_from[1] <= w$train_max[1])   # window 1: the burn-in in sample + the first held-out year
  expect_setequal(unlist(Map(seq, w$pred_from, w$pred_to)), 1850:1875); expect_false(any(duplicated(unlist(Map(seq, w$pred_from, w$pred_to)))))
  expect_error(oppw_windows(1870:1880, 1865L)); expect_error(oppw_windows(1850:1860, 1865L))
  ## decades: 1850s (10 years, onsets), 1860s (2 years: merged down), then 3 years of the 1860s; a decade without an onset is merged down; an early one up
  yr <- c(rep(1850:1859, each = 2), 1860, 1860, 1861, 1861); y <- rep(c(0, 1), length(yr) / 2)
  d <- oppw_decade_levels(yr, y); expect_equal(d$levels, 1850L); expect_equal(d$map(c(1850L, 1860L, 1870L)), c(1850L, 1850L, 1850L))
  d <- oppw_decade_levels(c(yr, 1862, 1862), c(y, 0, 1)); expect_equal(d$levels, c(1850L, 1860L)); expect_equal(d$map(c(1850L, 1860L, 1870L, 1990L)), c(1850L, 1860L, 1860L, 1860L))
  d <- oppw_decade_levels(c(yr, 1862, 1862), c(y[1:20], 0, 0, 0, 0, 0, 0)); expect_equal(d$levels, 1850L)                                         # three years, no onset: merged
  d <- oppw_decade_levels(c(1848, 1849, yr), c(1, 0, y)); expect_equal(d$levels, 1850L); expect_equal(d$map(1840L), 1850L)                       # two years of the 1840s: into the first level
  expect_error(oppw_decade_levels(c(1850, 1851), c(0, 1)))
})
test_that("the expanding window never uses year t or later: truncating or corrupting the data from t on changes no prediction up to t; a leaky window table is refused", {
  fr <- toy_frame()
  r <- oppw_expanding(fr, "y_mid", burn_in = TOY_BURN)
  expect_equal(dim(r$lq), c(nrow(fr), 6L)); expect_equal(colnames(r$lq), paste0("d", 0:5)); expect_true(all(r$lq < 0)); expect_false(anyNA(r$lq))
  expect_equal(r$windows$train_year_max, r$windows$train_max); expect_true(all(r$windows$train_year_max[-1] < r$windows$pred_from[-1]))
  expect_equal(r$held_out, fr$year > TOY_BURN); expect_true(all(r$windows$glm_max_diff[c(1, nrow(r$windows))] < 1e-6))
  ## (1) the data cut after year T give exactly the predictions of the full run up to T, for the MLE and every draw
  for (T in c(1861L, 1866L, 1870L, 1874L)) { rT <- oppw_expanding(fr[fr$year <= T, ], "y_mid", burn_in = TOY_BURN); expect_identical(rT$lq, r$lq[fr$year <= T, ], info = paste("truncated at", T)) }
  ## (2) the outcome of year t and everything after it (outcome and predictors) can be anything: the predictions up to t do not move; the year after does
  t0 <- 1868L; f2 <- fr; f2$y_mid[f2$year == t0] <- 1L - f2$y_mid[f2$year == t0]
  late <- f2$year > t0; set.seed(9); f2$y_mid[late] <- rbinom(sum(late), 1, 0.5); f2$log_capdist[late] <- f2$log_capdist[late] + rnorm(sum(late)); f2$a4_joinable[late] <- 1L - f2$a4_joinable[late]
  r2 <- oppw_expanding(f2, "y_mid", burn_in = TOY_BURN)
  expect_identical(r2$lq[fr$year <= t0, ], r$lq[fr$year <= t0, ])
  expect_false(isTRUE(all.equal(r2$lq[fr$year == t0 + 1L, 1], r$lq[fr$year == t0 + 1L, 1])))             # the test has teeth: year t0's outcome does move year t0 + 1
  ## (3) inside the burn-in the fit IS in sample: changing a burn-in outcome moves the burn-in predictions
  f3 <- fr; i <- which(f3$year == 1855L & f3$in_iv == 1L)[1:10]; f3$y_mid[i] <- 1L - f3$y_mid[i]
  expect_false(isTRUE(all.equal(oppw_expanding(f3, "y_mid", burn_in = TOY_BURN, n_draws = 0L)$lq[fr$year == 1852L, 1], r$lq[fr$year == 1852L, 1])))
  ## (4) only the IV rows are fitted: the outcome of a stage-2-only row is never used, but the row is predicted
  f4 <- fr; f4$y_mid[f4$in_iv == 0L] <- 1L - f4$y_mid[f4$in_iv == 0L]; expect_identical(oppw_expanding(f4, "y_mid", burn_in = TOY_BURN)$lq, r$lq); expect_true(all(is.finite(r$lq[fr$in_iv == 0L, ])))
  ## (5) a window table that trains on its own prediction year (or later) is refused by the loop's own check
  good <- oppw_windows; on.exit(assign("oppw_windows", good, envir = environment(oppw_expanding)), add = TRUE)
  assign("oppw_windows", function(years, burn_in) { w <- good(years, burn_in); w$train_max[5] <- w$pred_from[5]; w }, envir = environment(oppw_expanding))
  expect_error(oppw_expanding(fr, "y_mid", burn_in = TOY_BURN, n_draws = 0L))
  assign("oppw_windows", function(years, burn_in) { w <- good(years, burn_in); w$train_max[8] <- w$pred_from[8] + 2L; w }, envir = environment(oppw_expanding))
  expect_error(oppw_expanding(fr, "y_mid", burn_in = TOY_BURN, n_draws = 0L))
  assign("oppw_windows", good, envir = environment(oppw_expanding))
  expect_error(oppw_expanding(fr[order(-fr$year), ], "y_mid", burn_in = TOY_BURN))                        # an unsorted frame is refused
  expect_error(oppw_expanding(fr, "y_mid", a_cols = c("contig_land", "M_j"), burn_in = TOY_BURN))
})
test_that("the window fit is the logit MLE of stats::glm on the standardised design (coefficients, covariance, predictions); training constants only", {
  fr <- toy_frame(2L); tr <- fr[fr$year <= 1866L & fr$in_iv == 1L, ]
  dl <- oppw_decade_levels(tr$year, tr$y_mid); dz <- oppw_design(tr, OPPW_A_COLS, dl); f <- oppw_fit(dz$X, tr$y_mid)
  g <- stats::glm(tr$y_mid ~ dz$X - 1, family = stats::binomial())
  expect_equal(f$coef, unname(stats::coef(g)), tolerance = 1e-7); expect_equal(unname(f$V), unname(stats::vcov(g)), tolerance = 1e-5)
  fg <- oppw_fit_glm(dz$X, tr$y_mid); expect_equal(fg$coef, f$coef, tolerance = 1e-8); expect_equal(fg$V, f$V, tolerance = 1e-6)
  ## the design: training means / SDs, the product centred on the TRAINING mean of log distance, decade columns after the first level
  expect_equal(colnames(dz$X), c("(Intercept)", OPPW_A_COLS, "decade_1860")); expect_equal(unname(colMeans(dz$X[, OPPW_A_COLS])), rep(0, 6), tolerance = 1e-12)
  expect_equal(unname(apply(dz$X[, OPPW_A_COLS], 2, sd)), rep(1, 6), tolerance = 1e-12); expect_equal(dz$centre$ld_mean, mean(tr$log_capdist))
  prod_raw <- (tr$log_capdist - mean(tr$log_capdist)) * tr$M_lvl_i; expect_equal(unname(dz$X[, "log_capdist_x_init_M"]), (prod_raw - mean(prod_raw)) / sd(prod_raw))
  pr <- fr[fr$year == 1867L, ]; Xp <- oppw_design(pr, OPPW_A_COLS, dl, centre = dz$centre)$X
  expect_equal(unname(Xp[, "log_capdist"]), (pr$log_capdist - mean(tr$log_capdist)) / sd(tr$log_capdist))        # the prediction rows are scaled with the TRAINING constants
  r <- oppw_expanding(fr, "y_mid", burn_in = TOY_BURN)
  expect_equal(r$lq[fr$year == 1867L, "d0"], as.numeric(stats::plogis(Xp %*% f$coef, log.p = TRUE)), tolerance = 1e-8)
  expect_equal(r$coefs$coef[r$coefs$train_max == 1866L], f$coef, tolerance = 1e-8)
  ## the fitted probabilities do not depend on the centring and scaling (the logit MLE is invariant): raw columns give the same q
  raw <- cbind(1, as.matrix(tr[, c("contig_land", "contig_sea", "log_capdist")]), prod_raw, as.matrix(tr[, c("dymid_hist3", "a4_joinable")]), dz$X[, "decade_1860"])
  expect_equal(as.numeric(raw %*% stats::coef(stats::glm(tr$y_mid ~ raw - 1, family = stats::binomial()))), as.numeric(dz$X %*% f$coef), tolerance = 1e-6)
  expect_error(oppw_fit(dz$X, rep(0L, nrow(tr)))); tr0 <- tr; tr0$a4_joinable <- 0L; expect_error(oppw_design(tr0, OPPW_A_COLS, dl), "constant predictor")
  ## separation stops loudly
  Xs <- cbind(1, c(rep(0, 50), rep(1, 50))); expect_error(oppw_fit(Xs, c(rep(0L, 50), rep(1L, 50))))
})

## ---- (D) the draws ---------------------------------------------------------------------------------------------------------------------------------------
test_that("the draws are draws of the FULL coefficient vector with common random numbers: coherent across rows and windows, reproducible, never per-row", {
  fr <- toy_frame(4L); r <- oppw_expanding(fr, "y_mid", burn_in = TOY_BURN)
  expect_identical(oppw_expanding(fr, "y_mid", burn_in = TOY_BURN)$lq, r$lq)                                 # seed 123: reproducible
  expect_false(isTRUE(all.equal(oppw_expanding(fr, "y_mid", burn_in = TOY_BURN, seed = 7L)$lq[, 2], r$lq[, 2])))
  expect_identical(oppw_expanding(fr, "y_mid", burn_in = TOY_BURN, seed = 7L)$lq[, 1], r$lq[, 1])            # the MLE column does not depend on the seed
  expect_identical(oppw_expanding(fr, "y_mid", burn_in = TOY_BURN, n_draws = 0L)$lq[, 1], r$lq[, 1])
  set.seed(123); p_max <- 1L + length(OPPW_A_COLS) + length(unique(fr$decade)); expect_identical(r$z, matrix(rnorm(5L * p_max), nrow = 5L)); expect_equal(dim(r$z), c(5L, p_max))
  expect_equal(length(unique(apply(r$lq, 2, function(v) digest::digest(v)))), 6L)                             # five draws, all different from each other and from the MLE
  lgt <- function(lq) lq - log1p(-exp(lq))                                                                    # logit of q
  for (wi in c(1L, 4L, nrow(r$windows))) {
    w <- r$windows[wi, ]; tr <- fr[fr$year <= w$train_max & fr$in_iv == 1L, ]; pr <- fr$year >= w$pred_from & fr$year <= w$pred_to
    dl <- oppw_decade_levels(tr$year, tr$y_mid); dz <- oppw_design(tr, OPPW_A_COLS, dl); f <- oppw_fit(dz$X, tr$y_mid); Xp <- oppw_design(fr[pr, ], OPPW_A_COLS, dl, centre = dz$centre)$X
    L <- t(chol(f$V)); p <- length(f$coef)
    for (d in 1:5) {
      b_d <- f$coef + drop(L %*% r$z[d, seq_len(p)])                                                          # ONE coefficient vector per draw and window, from the draw's own z
      expect_equal(r$lq[pr, d + 1L], as.numeric(stats::plogis(Xp %*% b_d, log.p = TRUE)), tolerance = 1e-8, info = paste("window", wi, "draw", d))
      ## coherence: the shift of the linear predictor against the MLE is exactly linear in the design (a per-row draw would leave a residual)
      shift <- lgt(r$lq[pr, d + 1L]) - lgt(r$lq[pr, 1L]); expect_lt(max(abs(stats::lm.fit(Xp, shift)$residuals)), 1e-6)
    }
    ## the draws scatter around the MLE with the asymptotic covariance: the implied z is the stored one
    expect_equal(drop(solve(L, (f$coef + drop(L %*% r$z[2, seq_len(p)])) - f$coef)), r$z[2, seq_len(p)], tolerance = 1e-8)
  }
  ## a per-row marginal draw (what the procedure must never be) fails the coherence check, so the check has teeth
  pr <- fr$year == 1870L; set.seed(1); bad <- lgt(r$lq[pr, 1L]) + rnorm(sum(pr), 0, 0.3)
  tr <- fr[fr$year <= 1869L & fr$in_iv == 1L, ]; dl <- oppw_decade_levels(tr$year, tr$y_mid); dz <- oppw_design(tr, OPPW_A_COLS, dl)
  expect_gt(max(abs(stats::lm.fit(oppw_design(fr[pr, ], OPPW_A_COLS, dl, centre = dz$centre)$X, bad - lgt(r$lq[pr, 1L]))$residuals)), 0.05)
})
test_that("auxiliary diagnostics: log loss, precision-recall, calibration by hand", {
  y <- c(1L, 0L, 0L, 1L, 0L); q <- c(0.8, 0.1, 0.4, 0.3, 0.2)
  expect_equal(oppw_logloss(y, log(q)), -mean(c(log(0.8), log(0.9), log(0.6), log(0.3), log(0.8))))
  p <- oppw_pr(y, log(q)); expect_equal(p$prevalence, 0.4); expect_equal(p$avg_precision, (1 / 1 + 2 / 3) / 2); expect_equal(p$prec_at_rec50, 1); expect_equal(p$prec_at_rec75, 2 / 3)
  cal <- oppw_calibration(y, log(q), c("a", "a", "b", "b", "b"), "g")
  expect_equal(cal$n, c(2L, 3L)); expect_equal(cal$n_y, c(1L, 1L)); expect_equal(cal$obs_rate, c(0.5, 1 / 3)); expect_equal(cal$mean_q, c(0.45, 0.3)); expect_equal(cal$obs_over_pred, c(0.5 / 0.45, (1 / 3) / 0.3))
  expect_false(any(grepl("auc", names(c(p, cal)), ignore.case = TRUE)))
})

## ---- (E) Stan: the offset program -------------------------------------------------------------------------------------------------------------------------
toy_ctx_o <- function() {
  source(file.path(root, "sim", "simulate.R"), local = TRUE)
  sim <- simulate_targeting("M1b", S = 6L, Tt = 4L, K = 2L, seed = 123L, target_rate = 0.6)
  dat <- build_stan_data(sim$sy, sim$alts, sim$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"), joint = 1L, x_ref = rep(0.1, 4), D_mode = 3L)
  dat$G <- 3L; dat$ev_grp <- as.integer((seq_len(dat$E) - 1L) %% 3L)
  ## mixed choice sets (IV and stage-2 sets differ) so that the offset is exercised in both log-sum-exps separately
  nb <- dat$cs_end - dat$cs_start + 1L; big <- which(nb >= 4L & !(seq_len(dat$N1) %in% dat$ev_sy))[1:2]
  dat$in_iv[dat$cs_start[big[1]]] <- 0L; dat$in_s2[dat$cs_end[big[2]]] <- 0L
  event_alt_sets(dat)
}
rand_point_o <- function(dat, seed, Pz = dat$Pz) {
  set.seed(seed); S <- dat$S
  list(alpha0 = rnorm(1), alpha = rnorm(dat$P1, 0, 0.5), tau = rnorm(1, 0, 0.5), mu0 = rnorm(2), beta = matrix(rnorm(dat$P * 2, 0, 0.3), dat$P, 2),
       sigma_g = matrix(exp(rnorm(6, 0, 0.3)), 3, 2), rho_x = runif(1, -0.8, 0.8), sd_re = exp(rnorm(3, 0, 0.3)), L_Om = t(chol(matrix(c(1, 0.2, 0.1, 0.2, 1, 0.3, 0.1, 0.3, 1), 3))),
       z_re = matrix(rnorm(3 * S), 3, S), sigma_b = exp(rnorm(1, 0, 0.3)), z_b = rnorm(S), gamma = rnorm(Pz, 0, 0.5), delta = rbind(rnorm(2, 0, 0.7), rnorm(2, 0, 0.7)))
}
gq_o <- function(fit, par) as.numeric(rstan::extract(fit, par)[[par]])
fixed_o <- function(sm, dat, init) rstan::sampling(sm, data = dat, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(init), seed = 123, refresh = 0)
test_that("m1_context_offset with a zero offset IS m1_context: identical log density, gradient and generated quantities; the source differs by the offset only", {
  suppressPackageStartupMessages(library(rstan))
  p_c <- file.path(root, "stan", "build", "m1_context.rds"); p_o <- file.path(root, "stan", "build", "m1_context_offset.rds")
  ## the source: m1_context.stan plus the offset lines (removing them gives the original program text)
  src_c <- readLines(file.path(root, "stan", "m1_context.stan")); src_o <- readLines(file.path(root, "stan", "m1_context_offset.stan"))
  back <- src_o[!grepl("^  vector\\[N2\\] opp_off;", src_o)]
  back <- gsub(" + segment(opp_off, s, e - s + 1)", "", back, fixed = TRUE); back <- gsub(", vector opp_off,", ",", back, fixed = TRUE); back <- gsub("Z, Pz, opp_off,", "Z, Pz,", back, fixed = TRUE)
  back <- gsub("m1_context_offset:", "m1_context:", back, fixed = TRUE); back <- sub("[[:space:]]*// (run 1|zg_local = Z gamma).*$", "", back)
  code <- function(x) { x <- x[!grepl("^//", x)]; x[nzchar(trimws(x))] }
  expect_identical(code(back), code(src_c)); expect_false(any(grepl("opp_off", src_c)))
  if (!file.exists(p_c) || !file.exists(p_o)) skip("compiled m1_context / m1_context_offset missing")
  smc <- readRDS(p_c); smo <- readRDS(p_o); rstan_options(threads_per_chain = 1)
  dat <- toy_ctx_o(); dat0 <- c(dat, list(opp_off = numeric(dat$N2)))
  expect_identical(dat0[names(dat)], dat)                                                                   # the offset program's data = the current data + a zero vector
  p1 <- rand_point_o(dat, 1L); f_c <- fixed_o(smc, dat, p1); f_o <- fixed_o(smo, dat0, p1); expect_equal(f_c@mode, 0L); expect_equal(f_o@mode, 0L)
  for (s in 1:6) {
    p <- rand_point_o(dat, 100L + s); u_c <- unconstrain_pars(f_c, p); u_o <- unconstrain_pars(f_o, p); expect_identical(u_c, u_o)
    for (adj in c(TRUE, FALSE)) expect_identical(log_prob(f_o, u_o, adjust_transform = adj), log_prob(f_c, u_c, adjust_transform = adj))    # bit-identical
    expect_identical(as.numeric(grad_log_prob(f_o, u_o)), as.numeric(grad_log_prob(f_c, u_c)))
  }
  for (par in c("log_lik_stage1", "log_lik_stage2", "mu_state", "mu_event", "Sigma", "Sigma_g", "mu0_g", "Omega")) expect_identical(gq_o(f_o, par), gq_o(f_c, par), info = par)
  bad <- tryCatch(suppressWarnings(suppressMessages(quiet(fixed_o(smo, dat, p1)))), error = function(e) NULL)   # the offset program refuses data without opp_off
  expect_true(is.null(bad) || bad@mode != 0L)
})
test_that("a non-zero offset enters the stage-2 utility and the stage-1 IV identically: it is exactly a Z column with coefficient 1 (and lambda scales it)", {
  suppressPackageStartupMessages(library(rstan))
  p_c <- file.path(root, "stan", "build", "m1_context.rds"); p_o <- file.path(root, "stan", "build", "m1_context_offset.rds")
  if (!file.exists(p_c) || !file.exists(p_o)) skip("compiled m1_context / m1_context_offset missing")
  smc <- readRDS(p_c); smo <- readRDS(p_o); rstan_options(threads_per_chain = 1)
  dat <- toy_ctx_o(); set.seed(11); lq <- -rexp(dat$N2, 0.5)                                                  # ln q of every alternative row
  for (lambda in c(0.5, 1)) {
    dat_o <- c(dat, list(opp_off = lambda * lq)); dat_z <- dat; dat_z$Z <- cbind(dat$Z, lq); dat_z$Pz <- dat$Pz + 1L
    p1 <- rand_point_o(dat, 1L); pz1 <- p1; pz1$gamma <- c(p1$gamma, lambda)
    f_o <- fixed_o(smo, dat_o, p1); f_z <- fixed_o(smc, dat_z, pz1); f_c <- fixed_o(smc, dat, p1)
    ## both stages' log-likelihood terms: the IV (stage 1) and every event's choice probability (stage 2) are those of the Z-column model
    for (par in c("log_lik_stage1", "log_lik_stage2")) expect_equal(gq_o(f_o, par), gq_o(f_z, par), tolerance = 1e-12, info = paste(par, lambda))
    expect_gt(max(abs(gq_o(f_o, "log_lik_stage1") - gq_o(f_c, "log_lik_stage1"))), 1e-3)                      # and the offset does move the IV ...
    expect_gt(max(abs(gq_o(f_o, "log_lik_stage2") - gq_o(f_c, "log_lik_stage2"))), 1e-3)                      # ... and the stage-2 probabilities
    d <- vapply(1:5, function(s) { p <- rand_point_o(dat, 200L + s); pz <- p; pz$gamma <- c(p$gamma, lambda)
      log_prob(f_o, unconstrain_pars(f_o, p)) - log_prob(f_z, unconstrain_pars(f_z, pz)) }, numeric(1))
    expect_lt(diff(range(d)), 1e-9)                                                                           # the log densities differ by the prior term of the extra gamma only (a constant)
    expect_equal(d[1], 0.5 * lambda^2, tolerance = 1e-9)                                                      # -log N(lambda | 0, 1) up to the dropped constant
  }
  ## the IV of a state-year by hand: log sum over its in_iv rows of q^lambda exp(V) = the stage-1 linear predictor's lse (through log_lik_stage1 above);
  ## the stage-2 probabilities under the weights sum to one over a state-year's stage-2 set (the softmax of V + lambda ln q)
  set.seed(5); V <- rnorm(7); q <- runif(7, 0.01, 0.9)
  for (lambda in c(0, 0.5, 1)) { P <- q^lambda * exp(V) / sum(q^lambda * exp(V)); expect_equal(sum(P), 1); expect_equal(log(P), (V + lambda * log(q)) - log(sum(exp(V + lambda * log(q)))), tolerance = 1e-12) }
})

## ---- (F) 02_fit_model.R end to end on the toy cell ---------------------------------------------------------------------------------------------------------
toy_dir_o <- function(a) { a <- a |> mutate(dem1 = as.integer(ccode1 %in% c(10L, 20L)), dem2 = as.integer(ccode2 %in% c(10L, 20L)), major1 = as.integer(ccode1 %in% c(10L, 40L)), major2 = as.integer(ccode2 %in% c(10L, 40L)))
  a$dy_mid_tm1[a$ccode1 == 10L & a$ccode2 == 40L] <- 1L; a }
toy_roles_o <- tibble(initiator = c(10L, 10L, 30L, 30L, 30L, 40L), target = c(20L, 30L, 20L, 40L, 10L, 20L), year = 1990L, dispnum = c(1L, 4L, 2L, 3L, 6L, 7L), grp_aprime = c(0L, 1L, 2L, 0L, 0L, 2L))
toy_lq <- function(c1, c2, d) -(c1 / 10 + c2 / 100) * (1 + d / 10)                 # ln q of a toy row, by draw: a known function of the key
toy_weights <- function(v2, foundation = "data_entry", outcome = "mid") {
  a <- readRDS(file.path(v2, foundation, "alts_lag3_all.rds")); man <- read.csv(file.path(v2, foundation, "manifest.csv"), stringsAsFactors = FALSE)
  keys <- as.data.frame(a[, c("ccode1", "ccode2", "year")]); lq <- vapply(0:5, function(d) toy_lq(keys$ccode1, keys$ccode2, d), numeric(nrow(keys))); colnames(lq) <- paste0("d", 0:5)
  f <- c("alts_lag3_all.rds", "stage1_panel_lag3.rds")
  list(keys = keys, in_iv = a$in_iv, lq = lq, input = data.frame(file = f, sha256 = man$sha256[match(f, man$file)], role = "foundation"), a_cols = OPPW_A_COLS, outcome = outcome, foundation = foundation,
       lag = 3L, burn_in = OPPW_BURN_IN, n_draws = 5L, seed = 123L, git_commit = "toy", built = "toy")
}
test_that("02_fit_model.R on the toy cell: lambda = 0 = the current Stan data + a zero offset; lambda / draw pick the column; weights unchanged; stale weights stop; fits and the mixture", {
  tmp <- tempfile("oppw_e2e"); dir.create(file.path(tmp, "v2"), recursive = TRUE); file.copy(file.path(root, "R"), file.path(tmp, "v2"), recursive = TRUE)
  ok <- file.symlink(file.path(root, "stan"), file.path(tmp, "v2", "stan")); if (!ok) file.copy(file.path(root, "stan"), file.path(tmp, "v2"), recursive = TRUE)
  v2 <- file.path(tmp, "v2"); out_dir <- file.path(v2, "output"); dir.create(out_dir); dir.create(file.path(v2, "logs"))
  write_toy_cell(file.path(v2, "data_entry"), alts_edit = toy_dir_o, roles = toy_roles_o)
  base_args <- c("--sample=mid", "--lag=3", "--model=M1c", "--comparison=directed", "--entry=later", "--roles=bcj")
  run <- function(script, ...) { out <- suppressWarnings(system2("Rscript", c(file.path(v2, "R", script), ...), stdout = TRUE, stderr = TRUE, env = "V2_ROOT=")); list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  dry <- function(...) run("02_fit_model.R", base_args, "--dry_run=TRUE", ...)
  rd <- function(cell) readRDS(file.path(out_dir, paste0("standata_", cell, ".rds")))
  r <- dry(); expect_equal(r$status, 0L, info = paste(tail(r$out, 15), collapse = "\n")); base <- rd("mid_lag3_all_pr_M1c_cdir_bcj_ent")
  expect_false("opp" %in% names(base)); expect_false("opp_off" %in% names(base$dat)); expect_equal(base$stan_program, "m1_context")
  ## without the weights file the run stops and says how to build it
  r <- dry("--opp_weight=1"); expect_true(r$status != 0L); expect_true(any(grepl("08_opportunity_weights.R --entry=later --opp_outcome=mid", r$out, fixed = TRUE)))
  W <- toy_weights(v2); saveRDS(W, file.path(out_dir, "oppw_data_entry_mid_lag3.rds"))
  A <- suppressWarnings(quiet(assemble_cell_spec("mid", 3L, spec_o("--entry=later"), data_dir = file.path(v2, "data_entry"), model = "M1c"))); c1 <- A$alts$ccode1; c2 <- A$alts$ccode2
  ## lambda = 0: the Stan data of the current cell, field by field, plus opp_off = 0
  r <- dry("--opp_weight=0"); expect_equal(r$status, 0L, info = paste(tail(r$out, 15), collapse = "\n")); o0 <- rd("mid_lag3_all_pr_M1c_cdir_bcj_ent_opw0")
  expect_identical(o0$dat[names(base$dat)], base$dat); expect_identical(setdiff(names(o0$dat), names(base$dat)), "opp_off"); expect_identical(o0$dat$opp_off, numeric(base$dat$N2))
  expect_equal(o0$stan_program, "m1_context_offset"); expect_equal(o0$opp$lambda, 0); expect_equal(o0$opp$draw, 0L); expect_identical(o0$Z_cols, base$Z_cols); expect_identical(o0$W_cols, base$W_cols)
  expect_identical(o0$group_table, base$group_table)
  ## lambda = 0.5 and 1: lambda x ln q of the MLE column, row by row; the draw index picks its own column; everything else unchanged
  r <- dry("--opp_weight=0.5"); expect_equal(r$status, 0L, info = paste(tail(r$out, 15), collapse = "\n")); o05 <- rd("mid_lag3_all_pr_M1c_cdir_bcj_ent_opw05")
  expect_equal(o05$dat$opp_off, 0.5 * toy_lq(c1, c2, 0)); expect_identical(o05$dat[names(base$dat)], base$dat)
  r <- dry("--opp_weight=1"); expect_equal(r$status, 0L); o1 <- rd("mid_lag3_all_pr_M1c_cdir_bcj_ent_opw1"); expect_equal(o1$dat$opp_off, toy_lq(c1, c2, 0))
  r <- dry("--opp_weight=1", "--opp_draw=3"); expect_equal(r$status, 0L); o13 <- rd("mid_lag3_all_pr_M1c_cdir_bcj_ent_opw1_d3"); expect_equal(o13$dat$opp_off, toy_lq(c1, c2, 3)); expect_equal(o13$opp$draw, 3L)
  expect_identical(o13$dat[names(base$dat)], base$dat)
  ## the composite (D8) weights still sum as before, per role and in total; the per-role diagnostics carry them
  expect_identical(o1$dat$ev_w, base$dat$ev_w); br <- o1$opp$by_role
  expect_equal(br$w_events[br$role == "all_events"], sum(base$dat$ev_w)); expect_equal(br$w_events[1:3], base$group_table$w_usable); expect_equal(br$n_events[1:3], base$group_table$n_usable)
  expect_equal(sum(br$w_events[1:3]), sum(base$dat$ev_w)); expect_equal(br$role, c(BCJ_GROUPS, "all_events", "non_initiators", "all_state_years"))
  expect_true(file.exists(file.path(out_dir, "oppw_cell_mid_lag3_all_pr_M1c_cdir_bcj_ent_opw1.csv")))
  expect_true(any(grepl("SENSITIVITY specification, q is an opportunity weight", r$out, fixed = TRUE)))
  ## refusals at the command line and stale / incomplete weights
  for (bad in list(c("--opp_weight=0.5", "--opp_draw=2"), c("--opp_weight=1", "--ww2=drop"), c("--opp_weight=1", "--variant=polrel"), c("--opp_weight=1", "--reach=capdist"), c("--opp_draw=1"), c("--opp_weight"),
                   c("--opp_weight=1", "--opp_outcome=host3")))                                             # the last: no host3 weights file
    expect_true(dry(bad)$status != 0L, info = paste(bad, collapse = " "))
  expect_true(run("02_fit_model.R", "--sample=mid", "--lag=3", "--model=M1b", "--comparison=directed", "--entry=later", "--opp_weight=1", "--dry_run=TRUE")$status != 0L)
  W2 <- W; W2$input$sha256[1] <- paste(rep("0", 64), collapse = ""); saveRDS(W2, file.path(out_dir, "oppw_data_entry_mid_lag3.rds"))
  r <- dry("--opp_weight=1"); expect_true(r$status != 0L); expect_true(any(grepl("built from other data files", r$out, fixed = TRUE)))
  W3 <- W; W3$keys <- W$keys[-1, ]; W3$lq <- W$lq[-1, ]; saveRDS(W3, file.path(out_dir, "oppw_data_entry_mid_lag3.rds"))
  r <- dry("--opp_weight=1"); expect_true(r$status != 0L); expect_true(any(grepl("have no opportunity weight", r$out, fixed = TRUE)))
  W4 <- W; W4$outcome <- "host3"; saveRDS(W4, file.path(out_dir, "oppw_data_entry_mid_lag3.rds")); expect_true(dry("--opp_weight=1")$status != 0L)
  saveRDS(W, file.path(out_dir, "oppw_data_entry_mid_lag3.rds"))
  ## short fits: the point fit and the five draws at lambda = 1, then the mixture script; 03_extract.R accepts an offset fit
  if (!file.exists(file.path(root, "stan", "build", "m1_context_offset.rds"))) { unlink(tmp, recursive = TRUE); skip("compiled m1_context_offset missing") }
  fit <- function(...) run("02_fit_model.R", base_args, "--chains=2", "--iter=80", "--warmup=40", "--max_treedepth=6", ...)
  r <- fit("--opp_weight=1"); expect_equal(r$status, 0L, info = paste(tail(r$out, 25), collapse = "\n"))
  o <- readRDS(file.path(out_dir, "fit_mid_lag3_all_pr_M1c_cdir_bcj_ent_opw1.rds"))
  expect_equal(o$args$stan_program, "m1_context_offset"); expect_equal(o$opp$lambda, 1); expect_equal(o$dat$opp_off, toy_lq(c1, c2, 0))
  expect_equal(o$args$stan_source_sha256, trimws(readLines(file.path(root, "stan", "build", "m1_context_offset.sha256"))[1]))
  e <- run("03_extract.R", paste0("--fit=", file.path(out_dir, "fit_mid_lag3_all_pr_M1c_cdir_bcj_ent_opw1.rds"))); expect_equal(e$status, 0L, info = paste(tail(e$out, 25), collapse = "\n"))
  comb <- function(...) run("09_opp_combine.R", "--sample=mid", "--comparison=directed", "--entry=later", "--roles=bcj", ...)
  r <- comb("--opp_weight=1"); expect_true(r$status != 0L); expect_true(any(grepl("missing fit", r$out, fixed = TRUE)))                   # the draws are not there yet
  for (d in 1:5) { r <- fit("--opp_weight=1", paste0("--opp_draw=", d)); expect_equal(r$status, 0L, info = paste(tail(r$out, 25), collapse = "\n")) }
  r <- comb("--opp_weight=1"); expect_equal(r$status, 0L, info = paste(tail(r$out, 25), collapse = "\n"))
  tab <- read.csv(file.path(out_dir, "oppw_mixture_mid_lag3_all_pr_M1c_cdir_bcj_ent_opw1.csv"), stringsAsFactors = FALSE)
  expect_equal(nrow(tab), 6L + 6L + 1L); expect_setequal(tab$family, c("mu", "sigma", "rho")); expect_equal(sum(tab$family == "rho"), 1L); expect_true(all(tab$n_aux_draws == 5L)); expect_true(all(tab$n_post_draws_each == 80L))
  expect_equal(tab$sd_ratio, tab$mix_sd / tab$point_sd); expect_identical(tab$exceeds_10pct, tab$sd_ratio > 1.10); expect_setequal(tab$role[tab$family == "mu"], BCJ_GROUPS)
  mx <- readRDS(file.path(out_dir, "oppw_mixture_mid_lag3_all_pr_M1c_cdir_bcj_ent_opw1.rds")); expect_equal(nrow(mx$mix), 5L * 80L)
  expect_true(any(grepl("VERDICT: the mixture SDs", r$out, fixed = TRUE)))
  expect_true(comb("--opp_weight=0.5")$status != 0L); expect_true(comb("--opp_weight=1", "--opp_draw=2")$status != 0L)
  ## a draw fitted with another weights file is refused
  o5 <- readRDS(file.path(out_dir, "fit_mid_lag3_all_pr_M1c_cdir_bcj_ent_opw1_d5.rds")); o5$opp$weights_sha256 <- "x"; saveRDS(o5, file.path(out_dir, "fit_mid_lag3_all_pr_M1c_cdir_bcj_ent_opw1_d5.rds"))
  r <- comb("--opp_weight=1"); expect_true(r$status != 0L); expect_true(any(grepl("another weights file", r$out, fixed = TRUE)))
  unlink(tmp, recursive = TRUE)
})

## ---- (G) N_eff, the per-cell diagnostics, the mixture --------------------------------------------------------------------------------------------------
test_that("oppw_neff and oppw_cell_stats by hand; oppw_offset; oppw_mixture", {
  expect_equal(unname(oppw_neff(rep(0.2, 4))), c(0.8, 4)); expect_equal(unname(oppw_neff(c(1, 1e-12, 1e-12)))[2], 1, tolerance = 1e-9)
  expect_equal(unname(oppw_neff(c(0.5, 0.25, 0.25))), c(1, 1 / (0.25 + 0.0625 + 0.0625))); expect_error(oppw_neff(c(0.5, 0)))
  ## two state-years: n = 1 (rows 1-4; row 4 stage-2 only) initiates twice (roles a and b), n = 2 (rows 5-7) does not
  dat <- list(N1 = 2L, N2 = 7L, E = 2L, G = 2L, cs_start = c(1L, 5L), cs_end = c(4L, 7L), in_iv = c(1L, 1L, 1L, 0L, 1L, 1L, 1L), in_s2 = rep(1L, 7L), y1 = c(1L, 0L),
              ev_sy = c(1L, 1L), ev_row = c(2L, 4L), ev_grp = c(0L, 1L), ev_w = c(1, 0.5))
  q <- c(0.4, 0.004, 0.2, 0.1, 0.3, 0.3, 0.3); s <- oppw_cell_stats(dat, log(q), c("a", "b"), low_prob = 0.1)
  expect_equal(s$threshold, unname(quantile(q, 0.1))); expect_equal(s$threshold, 0.004 + 0.6 * 0.096); b <- s$by_role; expect_equal(b$role, c("a", "b", "all_events", "non_initiators", "all_state_years"))
  expect_equal(b$n_events, c(1L, 1L, 2L, NA, NA)); expect_equal(b$w_events, c(1, 0.5, 1.5, NA, NA)); expect_equal(b$n_low, c(1L, 0L, 1L, NA, NA))        # the chosen 0.004 is low; 0.1 is not
  expect_equal(b$share_low_n, c(1, 0, 0.5, NA, NA)); expect_equal(b$w_low, c(1, 0, 1, NA, NA)); expect_equal(b$share_low_w, c(1, 0, 1 / 1.5, NA, NA)); expect_equal(b$q_chosen_median[1:2], c(0.004, 0.1))
  sq1 <- 0.4 + 0.004 + 0.2; ne1 <- sq1^2 / (0.4^2 + 0.004^2 + 0.2^2)                                                 # IV rows only (row 4 is stage-2 only)
  expect_equal(b$sum_q_median, c(sq1, sq1, sq1, 0.9, mean(c(sq1, 0.9)))); expect_equal(b$n_eff_median, c(ne1, ne1, ne1, 3, mean(c(ne1, 3)))); expect_equal(b$n_sy, c(1L, 1L, 1L, 1L, 2L))
  expect_equal(b$n_iv_median, c(3, 3, 3, 3, 3))
  ## oppw_offset: key match in the cell's row order, lambda and draw, refusals
  W <- list(keys = data.frame(ccode1 = c(10L, 10L, 20L), ccode2 = c(20L, 30L, 10L), year = 1990L), lq = cbind(d0 = c(-1, -2, -3), d1 = c(-1.5, -2.5, -3.5), d2 = 0, d3 = 0, d4 = 0, d5 = 0),
            input = data.frame(file = "alts_lag3_all.rds", sha256 = "abc", role = "foundation"))
  alts <- data.frame(ccode1 = c(20L, 10L), ccode2 = c(10L, 30L), year = 1990L)
  expect_equal(oppw_offset(alts, W, 1, 0L)$off, c(-3, -2)); expect_equal(oppw_offset(alts, W, 0.5, 0L)$off, c(-1.5, -1)); expect_identical(oppw_offset(alts, W, 0, 0L)$off, c(0, 0))
  expect_equal(oppw_offset(alts, W, 0, 0L)$lq, c(-3, -2)); expect_equal(oppw_offset(alts, W, 1, 1L)$off, c(-3.5, -2.5))
  expect_error(oppw_offset(alts, W, 0.5, 1L)); expect_error(oppw_offset(alts, W, 0.3, 0L)); expect_error(oppw_offset(rbind(alts, data.frame(ccode1 = 30L, ccode2 = 10L, year = 1990L)), W, 1, 0L), "no opportunity weight")
  expect_silent(oppw_offset(alts, W, 1, 0L, manifest = data.frame(file = "alts_lag3_all.rds", sha256 = "abc")))
  expect_error(oppw_offset(alts, W, 1, 0L, manifest = data.frame(file = "alts_lag3_all.rds", sha256 = "abd")), "other data files")
  ## the mixture: equal weights, the SD ratio and the 10% flag
  set.seed(2); pt <- cbind(a = rnorm(4000), b = rnorm(4000)); dl <- lapply(1:5, function(d) cbind(a = rnorm(4000, mean = (d - 3) * 0.6), b = rnorm(4000)))
  m <- oppw_mixture(dl, pt); t <- m$table
  expect_equal(nrow(m$mix), 20000L); expect_equal(t$mix_mean, unname(colMeans(do.call(rbind, dl)))); expect_equal(t$mix_sd, unname(apply(do.call(rbind, dl), 2, sd)))
  expect_equal(t$exceeds_10pct, c(TRUE, FALSE)); expect_gt(t$sd_ratio[1], 1.2); expect_lt(abs(t$sd_ratio[2] - 1), 0.05); expect_equal(t$between_draw_sd[1], sd(vapply(dl, function(x) mean(x[, "a"]), numeric(1))))
  expect_error(oppw_mixture(c(dl[1:4], list(dl[[5]][1:100, ])), pt), "equally weighted"); expect_error(oppw_mixture(dl[1], pt))
})

## ---- (H) slurm ---------------------------------------------------------------------------------------------------------------------------------------------
test_that("slurm: 13_opp_weights.sbatch issues the auxiliary runs, the 32 valid M1c cells and the four mixtures (stub Rscript); submit_opp_weights.sh chains them afterok (stub sbatch)", {
  tmp <- tempfile("oppw_slurm"); dir.create(file.path(tmp, "v2"), recursive = TRUE); dir.create(file.path(tmp, "bin"))
  for (d in c("slurm", "env")) file.copy(file.path(root, d), file.path(tmp, "v2"), recursive = TRUE)
  dir.create(file.path(tmp, "v2", "R")); file.copy(file.path(root, "R", "01_build_data.R"), file.path(tmp, "v2", "R"))
  rec <- file.path(tmp, "calls.txt"); ctr <- file.path(tmp, "ctr.txt")
  writeLines(c("#!/bin/bash", sprintf("echo \"$*\" >> %s", rec)), file.path(tmp, "bin", "Rscript")); writeLines(c("#!/bin/bash", "exit 0"), file.path(tmp, "bin", "module"))
  writeLines(c("#!/bin/bash", sprintf("n=$(cat %s 2>/dev/null || echo 1000); n=$((n + 1)); echo $n > %s", ctr, ctr), sprintf("echo \"sbatch job=$n args: $*\" >> %s", rec), "echo \"$n;cluster\""), file.path(tmp, "bin", "sbatch"))
  for (b in c("Rscript", "module", "sbatch")) Sys.chmod(file.path(tmp, "bin", b), "0755")
  env0 <- c(sprintf("PATH=%s:%s", file.path(tmp, "bin"), Sys.getenv("PATH")), sprintf("SLURM_SUBMIT_DIR=%s", tmp), "V2_ROOT=", "BASH_ENV=")
  task <- function(i, env = character(0)) { unlink(rec); out <- suppressWarnings(system2("bash", file.path(tmp, "v2", "slurm", "13_opp_weights.sbatch"), stdout = TRUE, stderr = TRUE, env = c(env0, if (!is.null(i)) sprintf("SLURM_ARRAY_TASK_ID=%d", i), env)))
    list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out, call = if (file.exists(rec)) readLines(rec) else character(0)) }
  cell_of <- function(call) { a <- strsplit(call, " ")[[1]]; a <- a[grepl("^--", a)]; cell_name(get_arg(a, "sample", ""), 3L, get_arg(a, "model", ""), parse_spec(a)) }
  want <- unlist(lapply(c("ent", "dy"), function(ds) unlist(lapply(c("fmid", "mid"), function(s) sprintf("%s_lag3_all_pr_M1c_cdir_bcj_%s_%s_tlong", s, ds, c("opw0", "opw05", "opw1", paste0("opw1_d", 1:5)))))))
  expect_length(unique(want), 32L)
  for (i in 0:31) { r <- task(i, "DRY_RUN=TRUE"); expect_equal(r$status, 0L, info = paste(i, paste(r$out, collapse = "\n"))); expect_length(r$call, 1L)
    expect_equal(cell_of(r$call), want[i + 1L], info = i); expect_match(r$call, "02_fit_model.R", fixed = TRUE); expect_match(r$call, "--dry_run=TRUE", fixed = TRUE)
    expect_match(r$call, "--iter=6000 --warmup=3000", fixed = TRUE); expect_false(grepl("ww2", r$call)) }
  expect_equal(want[c(1, 14, 19)], c("fmid_lag3_all_pr_M1c_cdir_bcj_ent_opw0_tlong", "mid_lag3_all_pr_M1c_cdir_bcj_ent_opw1_d3_tlong", "fmid_lag3_all_pr_M1c_cdir_bcj_dy_opw1_tlong"))   # the header's examples
  expect_true(task(32L)$status != 0L); expect_true(task(NULL)$status != 0L); expect_match(task(2L)$call, "--dry_run=FALSE", fixed = TRUE)
  r <- task(5L, "OPP_OUTCOME=host3"); expect_equal(cell_of(r$call), "fmid_lag3_all_pr_M1c_cdir_bcj_ent_opw1_oh3_d3_tlong"); expect_true(task(5L, "OPP_OUTCOME=war")$status != 0L)
  ## the auxiliary step: both foundations x both outcomes; the mixture step: the four lambda = 1 cells
  r <- task(0L, "STEP=aux"); expect_equal(r$status, 0L); expect_length(r$call, 4L); expect_true(all(grepl("08_opportunity_weights.R", r$call, fixed = TRUE)))
  expect_equal(sub(".*08_opportunity_weights.R ", "", r$call), c("--entry=later --opp_outcome=mid", "--events=dymid --opp_outcome=mid", "--entry=later --opp_outcome=host3", "--events=dymid --opp_outcome=host3"))
  r <- task(0L, c("STEP=aux", "OPP_OUTCOMES=mid")); expect_length(r$call, 2L)
  r <- task(0L, "STEP=combine"); expect_equal(r$status, 0L); expect_length(r$call, 4L); expect_true(all(grepl("09_opp_combine.R", r$call, fixed = TRUE)))
  expect_equal(vapply(r$call, function(cl) { a <- strsplit(cl, " ")[[1]]; a <- a[grepl("^--", a)]; cell_name(get_arg(a, "sample", ""), 3L, "M1c", parse_spec(a)) }, character(1), USE.NAMES = FALSE),
               c("fmid_lag3_all_pr_M1c_cdir_bcj_ent_opw1_tlong", "mid_lag3_all_pr_M1c_cdir_bcj_ent_opw1_tlong", "fmid_lag3_all_pr_M1c_cdir_bcj_dy_opw1_tlong", "mid_lag3_all_pr_M1c_cdir_bcj_dy_opw1_tlong"))
  expect_true(task(0L, "STEP=other")$status != 0L)
  hdr <- readLines(file.path(root, "slurm", "13_opp_weights.sbatch")); hdr0 <- readLines(file.path(root, "slurm", "12_run0.sbatch"))
  expect_true(all(c("#SBATCH --account=open", "#SBATCH --partition=open", "#SBATCH --requeue", "#SBATCH --time=48:00:00", "#SBATCH --array=0-31", "#SBATCH --mem=128G") %in% hdr))
  for (k in c("--account", "--partition", "--requeue", "--time", "--mem", "--nodes")) expect_identical(grep(paste0("^#SBATCH ", k), hdr, value = TRUE), grep(paste0("^#SBATCH ", k), hdr0, value = TRUE), info = k)   # the M1 resources of run 0
  ## the submit helper: aux, then the array afterok on it, then the mixture afterok on the array; nothing else
  owd <- setwd(tmp); on.exit(setwd(owd), add = TRUE); unlink(rec); unlink(ctr)
  out <- suppressWarnings(system2("bash", "v2/slurm/submit_opp_weights.sh", stdout = TRUE, stderr = TRUE, env = env0)); expect_null(attr(out, "status"), info = paste(out, collapse = "\n"))
  calls <- readLines(rec); expect_length(calls, 3L)
  expect_match(calls[1], "^sbatch job=1001 args: --parsable --job-name=v2_oppw_aux --array=0 --export=ALL,STEP=aux --time=06:00:00 --mem=32G --cpus-per-task=1 v2/slurm/13_opp_weights.sbatch$")
  expect_match(calls[2], "^sbatch job=1002 args: --parsable --array=0-31 --export=ALL,STEP=fit --dependency=afterok:1001 --kill-on-invalid-dep=yes v2/slurm/13_opp_weights.sbatch$")
  expect_match(calls[3], "^sbatch job=1003 args: --parsable --job-name=v2_oppw_mix --array=0 --export=ALL,STEP=combine --time=06:00:00 --mem=64G --cpus-per-task=1 --dependency=afterok:1002 --kill-on-invalid-dep=yes v2/slurm/13_opp_weights.sbatch$")
  expect_false(any(grepl("jfe4_cr_default", calls)))                                                         # everything on the open account unless FIT_OPTS says otherwise
  unlink(rec); out <- suppressWarnings(system2("bash", "v2/slurm/submit_opp_weights.sh", stdout = TRUE, stderr = TRUE, env = c(env0, "A_JOB=777", "FIT_IDX=16-23", "COMBINE=FALSE")))
  calls <- readLines(rec); expect_length(calls, 1L); expect_match(calls[1], "--array=16-23 --export=ALL,STEP=fit --dependency=afterok:777", fixed = TRUE)
  out <- suppressWarnings(system2("bash", c("v2/slurm/submit_opp_weights.sh", "--array=0-3"), stdout = TRUE, stderr = TRUE, env = env0)); expect_false(is.null(attr(out, "status")))
  setwd(owd); unlink(tmp, recursive = TRUE)
})

## ---- (I) the real builds -----------------------------------------------------------------------------------------------------------------------------------
test_that("real lag-3 dymid build: the January presence reproduces multi_mid_jan; outcome, history and A4 recomputed independently; no year-t onset in the history", {
  d <- file.path(root, "data_dymid"); if (!file.exists(file.path(d, "alts_lag3_all.rds"))) skip("v2/data_dymid not built")
  suppressPackageStartupMessages(library(peacesciencer))
  panel <- readRDS(file.path(d, "stage1_panel_lag3.rds")); ev <- readRDS(file.path(d, "events.rds"))
  alts <- readRDS(file.path(d, "alts_lag3_all.rds")); alts <- alts[alts$year %in% c(1914:1920, 1939:1946, 1990:1993), ]
  dd <- cow_mid_dirdisps |> mutate(across(c(dispnum, ccode1, ccode2, year), as.integer))
  psd <- participant_spell_dates(dd, dispute_date_bounds(cow_mid_disps |> mutate(dispnum = as.integer(dispnum)), dd), part = gml_part)
  jp <- oppw_jan_presence(psd, dd, panel |> select(ccode, year))
  chk <- panel |> select(ccode, year, multi_mid_jan) |> inner_join(jp$flags, by = c("ccode", "year")); expect_equal(nrow(chk), nrow(panel)); expect_identical(chk$multi_mid_jan.x, chk$multi_mid_jan.y)
  fr <- oppw_frame(alts, panel, ev, jp); expect_equal(nrow(fr), sum(!is.na(alts$log_capdist))); expect_false(any(OPPW_TARGET_STRENGTH %in% names(fr)))
  ## the outcome, independently: an events row of the unordered pair in the row's own year
  pk <- paste(pmin(ev$initiator, ev$target), pmax(ev$initiator, ev$target), ev$year); fk <- paste(pmin(fr$ccode1, fr$ccode2), pmax(fr$ccode1, fr$ccode2), fr$year)
  expect_identical(fr$y_mid, as.integer(fk %in% pk)); expect_identical(fr$y_host3, as.integer(fk %in% pk[pmax(ev$hostlev1, ev$hostlev2) >= 3L])); expect_gt(sum(fr$y_mid), sum(fr$y_host3))
  ## the history, independently: an onset in t - 1, t - 2 or t - 3; a pair whose ONLY onset in [t - 3, t] is in year t has no history
  on <- unique(data.frame(a = pmin(ev$initiator, ev$target), b = pmax(ev$initiator, ev$target), y = ev$year)); fa <- pmin(fr$ccode1, fr$ccode2); fb <- pmax(fr$ccode1, fr$ccode2)
  h <- Reduce(`|`, lapply(1:3, function(k) paste(fa, fb, fr$year - k) %in% paste(on$a, on$b, on$y))); expect_identical(fr$dymid_hist3, as.integer(h))
  expect_true(any(fr$y_mid == 1L & fr$dymid_hist3 == 0L)); expect_true(any(fr$y_mid == 0L & fr$dymid_hist3 == 1L))
  ## A4: the United States in January 1941 is not a party to the Second World War (dispute 258) and Germany is present in it; Germany is a party
  k <- function(i, j, t) which(fr$ccode1 == i & fr$ccode2 == j & fr$year == t)
  expect_true(258L %in% jp$multi$dispnum[jp$multi$year == 1941L]); expect_true(255L %in% jp$present$ccode[jp$present$dispnum == 258L & jp$present$year == 1941L])
  expect_false(2L %in% jp$already_in$ccode[jp$already_in$dispnum == 258L & jp$already_in$year == 1941L]); expect_equal(fr$a4_joinable[k(2L, 255L, 1941L)], 1L)
  ## independent recomputation of A4 on the 1941 rows: j present in a multilateral dispute i is not already in
  r41 <- fr[fr$year == 1941L, ]; pres <- jp$present[jp$present$year == 1941L, ]; inn <- jp$already_in[jp$already_in$year == 1941L, ]
  a4 <- vapply(seq_len(nrow(r41)), function(r) { dj <- pres$dispnum[pres$ccode == r41$ccode2[r]]; as.integer(length(setdiff(dj, inn$dispnum[inn$ccode == r41$ccode1[r]])) > 0L) }, integer(1))
  expect_identical(r41$a4_joinable, a4); expect_gt(sum(a4), 0L)
  expect_true(all(fr$a4_joinable[fr$year == 1990L] %in% 0:1)); expect_equal(fr$contig_sea + fr$contig_land, alts$contig_any[match(paste(fr$ccode1, fr$ccode2, fr$year), paste(alts$ccode1, alts$ccode2, alts$year))])
})
