## Run 0 and the lag-3 MID-history choice set (built 2026-11-10; README section 13).
## Revised the same day after the audit of d60da78 (PI decisions): --variant=mid3 is REFUSED by parse_spec() and cell_name() (tested here); its construction
## stays in the code and the blocks (B), (C) and the mid3 half of (F) are kept as UNIT TESTS of mid3_history() and of the assembly called directly with a
## hand-set spec (spec_mid3() below), not as tests of a runnable cell. Run 0b carries the polrel cell's events, roles and weights (no D8 rescale).
## (A) the switches reach / subsample and the variant value mid3: choices, suffixes, positions, names, refusals; every existing name and spec string unchanged
##     (hard-coded, and against the git version before the change); defaults byte-identical: the toy assembly of the existing cells is identical() to the
##     assembly by the code of commit aeb876b.
## (B) mid3_history(): the window [t - 3, t - 1], never year t (a test that fails when t is included), both directions, every row of the events file.
## (C) the mid3 assembly on a toy cell: IV / stage-2 sets, the chosen target always present, event-level sets = history rows + the event's OWN chosen target,
##     chosen-only events dropped and counted by role, state-years without a history IV alternative leave stage 1; Stan data (ragged sets).
## (D) run 0: reach (the two columns recomputed by hand, everything else unchanged), prtarget (the events polrel drops, full set, D8 weights renormalised to
##     sum 1 per split group), bilateral (M1b, stage 1 unchanged), majinit / mininit (both stages, W / Z columns), refusals.
## (E) 02_fit_model.R --dry_run end to end on the toy cell for every new cell type; the two sbatch files and the submit helper against stubs.
## (F) the real lag-3 builds (skipped without v2/data_entry / v2/data_dymid): mid3 history recomputed independently, no year-t onset, chosen targets present;
##     prtarget drops exactly the events --variant=polrel drops.
suppressPackageStartupMessages({ library(testthat); library(dplyr) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "assemble_cell.R"))
BASE_COMMIT <- "aeb876b"                                   # review-dymid tip the run-0 build started from

## ---- toy fixtures ------------------------------------------------------------------------------------------------------------------------------------
## the directed columns of the toy alternatives (as in test_entry_roles.R): democracies 10, 20; major powers 10, 40
toy_dir <- function(a) a |> mutate(dem1 = as.integer(ccode1 %in% c(10L, 20L)), dem2 = as.integer(ccode2 %in% c(10L, 20L)),
                                   major1 = as.integer(ccode1 %in% c(10L, 40L)), major2 = as.integer(ccode2 %in% c(10L, 40L)))
## aprime roles of the six toy events: 10->20 alone, 10->30 co-originator, 30->20 joiner, 30->40 alone, 30->10 alone, 40->20 joiner
toy_roles <- tibble(initiator = c(10L, 10L, 30L, 30L, 30L, 40L), target = c(20L, 30L, 20L, 40L, 10L, 20L), year = 1990L, dispnum = c(1L, 4L, 2L, 3L, 6L, 7L), grp_aprime = c(0L, 1L, 2L, 0L, 0L, 2L))
## MID-history rows appended to the toy events: earlier onsets in NO sample (flagged as escalations, as the real file's non-sample rows are), all columns
## copied from the first toy row. hist = data.frame(initiator, target, year).
add_history <- function(hist) function(ev) {
  h <- ev[rep(1L, nrow(hist)), ]
  h$initiator <- as.integer(hist$initiator); h$target <- as.integer(hist$target); h$year <- as.integer(hist$year)
  h$dispnum <- 100L + seq_len(nrow(hist)); h$dispnums <- as.character(h$dispnum); h$n_disputes <- 1L; h$escalation <- TRUE; h$is_ww2 <- FALSE
  for (cc in grep("^in_", names(h), value = TRUE)) h[[cc]] <- FALSE
  bind_rows(ev, h)
}
## toy A history: 10-20 (onset 1988), 10-40 (1989), 20-40 (1987, the first year of the window); 10-30 has an onset in 1986 (t - 4: outside the window)
## and its sample events 10->30 / 30->10 in 1990 (year t: never history). 30-40 has history through dy_mid_tm1 only (the toy's (30, 40) row).
HIST_A <- data.frame(initiator = c(20L, 40L, 20L, 30L), target = c(10L, 10L, 40L, 10L), year = c(1988L, 1989L, 1987L, 1986L))
key_ev  <- function(A) paste(A$sy$ccode[A$ev$sy_row], A$alts$ccode2[A$ev$alt_row])
key_evy <- function(A) paste(A$sy$ccode[A$ev$sy_row], A$alts$ccode2[A$ev$alt_row], A$sy$year[A$ev$sy_row])           # with the year (the real builds)
## the toy's only lagged dyadic MID is (30, 40); a second one on (10, 40) keeps dy_mid_tm1 non-constant when state-year 30 is cut (run 0d, majors)
toy_dir2 <- function(a) { a <- toy_dir(a); a$dy_mid_tm1[a$ccode1 == 10L & a$ccode2 == 40L] <- 1L; a }
key_row <- function(A, rows) paste(A$sy$ccode[A$alts$sy_row[rows]], A$alts$ccode2[rows])
spec_dir <- function(...) parse_spec(c("--comparison=directed", "--entry=later", ...))
## the mid3 spec set by hand: parse_spec() refuses --variant=mid3, so the unit tests of the construction build the spec without the command line
spec_mid3 <- function(data_arg = "--entry=later") { s <- parse_spec(c("--comparison=directed", data_arg, "--roles=bcj")); s$variant <- "mid3"; s }

## ---- (A) the switches ----------------------------------------------------------------------------------------------------------------------------------
test_that("spec: reach, subsample and variant = mid3: choices, suffixes, positions, names; existing names and strings unchanged (hard-coded)", {
  expect_equal(SPEC_CHOICES$reach, c("none", "capdist")); expect_equal(SPEC_CHOICES$subsample, c("none", "prtarget", "bilateral", "majinit", "mininit"))
  expect_equal(SPEC_CHOICES$variant, c("all", "polrel", "mid3")); expect_equal(MID3_WINDOW, 3L)
  expect_equal(SPEC_SUFFIX$reach, c(capdist = "_reach")); expect_equal(SPEC_SUFFIX$subsample, c(prtarget = "_sprt", bilateral = "_sbil", majinit = "_smaj", mininit = "_smin"))
  expect_equal(which(names(SPEC_SUFFIX) == "reach"), which(names(SPEC_SUFFIX) == "mu_shift") + 1L)            # right after _cdir / _mst
  expect_equal(which(names(SPEC_SUFFIX) == "subsample"), which(names(SPEC_SUFFIX) == "joiner_set") + 1L)      # right after the role / joiner-set suffixes
  expect_false(any(duplicated(c(unlist(SPEC_SUFFIX), "_K2", "_pilot", "_STALE")))); expect_false(any(grepl("^_t[A-Za-z0-9]", unlist(SPEC_SUFFIX))))
  ## the new cells
  expect_equal(cell_name("fmid", 3L, "M1c", spec_dir("--roles=bcj", "--reach=capdist", "--tag=long")), "fmid_lag3_all_pr_M1c_cdir_reach_bcj_ent_tlong")
  expect_equal(cell_name("fmid", 3L, "M1c", spec_dir("--roles=bcj", "--subsample=prtarget")), "fmid_lag3_all_pr_M1c_cdir_bcj_sprt_ent")
  expect_equal(cell_name("mid", 3L, "M1b", spec_dir("--subsample=bilateral")), "mid_lag3_all_pr_M1b_cdir_sbil_ent")
  expect_equal(cell_name("fmid", 3L, "M1c", parse_spec(c("--comparison=directed", "--events=dymid", "--roles=bcj", "--subsample=majinit"))), "fmid_lag3_all_pr_M1c_cdir_bcj_smaj_dy")
  expect_equal(cell_name("mid", 3L, "M1c", spec_dir("--roles=bcj", "--subsample=mininit")), "mid_lag3_all_pr_M1c_cdir_bcj_smin_ent")
  ## mid3 (audit P1-1, PI decision 2026-11-10): refused on the command line and as a cell name; the value stays in SPEC_CHOICES so the refusal can say why
  expect_error(spec_dir("--roles=bcj", "--variant=mid3", "--tag=long"), "--variant=mid3 is refused \\(audit P1-1, PI decision 2026-11-10\\).*ONLY non-history member.*year t would define the set by the outcome")
  s3 <- spec_mid3()
  for (m in c("M1c", "M2c", "M1b")) expect_error(cell_name("fmid", 3L, m, s3), "cell_name: --variant=mid3 is refused", info = m)
  expect_error(m1_cell_for_m2("fmid", 3L, "M2c", s3), "is refused")
  expect_equal(spec_alts_variant(s3), "all"); expect_equal(spec_alts_variant(spec_dir("--roles=bcj", "--variant=polrel")), "polrel"); expect_equal(spec_alts_variant(parse_spec()), "all")
  expect_true(grepl("reach=capdist", spec_string(spec_dir("--reach=capdist")), fixed = TRUE))
  expect_true(grepl("subsample=prtarget", perm_part_spec_string(spec_dir("--subsample=prtarget")), fixed = TRUE))
  ## every existing name and string, hard-coded from the code before the change
  expect_identical(spec_string(parse_spec()), "pairs=all context=all coal_ties=inclusive coal_drop=none direction=cowwar ww2=include weights=dispute controls=minimal major=none interact=none distance=none contiguity=land changes=raw rho_prior=beta25 sigma_prior=ln1 comparison=none roles=stageb ame_K=2 tag=")
  expect_identical(spec_string(spec_dir("--roles=bcj", "--variant=polrel", "--tag=long")),
                   "pairs=all context=all coal_ties=inclusive coal_drop=none direction=cowwar ww2=include weights=dispute controls=minimal major=none interact=none distance=none contiguity=land changes=raw rho_prior=beta25 sigma_prior=ln1 comparison=directed roles=bcj entry=later variant=polrel ame_K=2 tag=long")
  expect_identical(spec_suffix(parse_spec()), ""); expect_true(spec_is_main(parse_spec())); expect_equal(spec_reach(parse_spec()), "none"); expect_equal(spec_subsample(parse_spec()), "none")
  expect_equal(spec_reach(list()), "none"); expect_equal(spec_subsample(list()), "none")            # a spec saved before 2026-11-10
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec()), "fmid_lag3_all_pr_M1b"); expect_equal(cell_name("war", 3L, "M2b", parse_spec()), "war_lag3_all_pr_M2b_K2")
  expect_equal(cell_name("fmid", 3L, "M1c", spec_dir("--roles=bcj", "--tag=long")), "fmid_lag3_all_pr_M1c_cdir_bcj_ent_tlong")
  expect_equal(cell_name("fmid", 3L, "M2c", spec_dir("--roles=bcj", "--variant=polrel", "--tag=long")), "fmid_lag3_polrel_pr_M2c_K2_cdir_bcj_ent_tlong")
  expect_equal(cell_name("mid", 3L, "M1c", parse_spec(c("--comparison=directed", "--events=dymid", "--roles=bcj", "--joiner_set=opposing"))), "mid_lag3_all_pr_M1c_cdir_bcj_jset_dy")
  expect_equal(directed_w_cols(), W_CONVENTIONAL); expect_equal(directed_z_cols(), Z_DIRECTED); expect_equal(length(Z_DIRECTED), 9L)
  expect_equal(directed_w_cols("majinit"), setdiff(W_CONVENTIONAL, "cowmaj_tm1")); expect_equal(directed_z_cols("none", "mininit"), setdiff(Z_DIRECTED, "tgt_major_x_init_major"))
  expect_equal(directed_z_cols("capdist"), c(Z_DIRECTED, "log_capdist_sq", "log_capdist_x_init_M")); expect_false(any(Z_REACH %in% Z_CENTRED))
  expect_equal(directed_w_cols("prtarget"), W_CONVENTIONAL); expect_equal(directed_z_cols("none", "bilateral"), Z_DIRECTED)
})

test_that("spec: refusals of the invalid combinations (parse_spec and cell_name), bare flags and bad values", {
  expect_error(parse_spec("--reach=capdist"), "needs --comparison=directed"); expect_error(parse_spec("--subsample=prtarget"), "needs --comparison=directed")
  expect_error(parse_spec(c("--variant=mid3", "--roles=aprime")), "--variant=mid3 is refused"); expect_error(parse_spec("--variant=mid3"), "--variant=mid3 is refused")
  expect_error(parse_spec("--reach"), "needs a value"); expect_error(parse_spec("--subsample"), "needs a value"); expect_error(parse_spec("--reach=TRUE"), "must be one of")
  expect_error(parse_spec(c("--comparison=directed", "--subsample=majors")), "must be one of"); expect_error(parse_spec(c("--comparison=directed", "--variant=mid5")), "must be one of")
  d <- "--comparison=directed"
  expect_error(parse_spec(c(d, "--reach=capdist", "--mu_shift=traits")), "not specified"); expect_error(parse_spec(c(d, "--entry=later", "--roles=bcj", "--reach=capdist", "--variant=polrel")), "full choice set")
  expect_error(parse_spec(c(d, "--reach=capdist", "--subsample=prtarget")), "one switch each")
  expect_error(parse_spec(c(d, "--entry=later", "--roles=bcj", "--subsample=prtarget", "--variant=polrel")), "full choice set")
  expect_error(parse_spec(c(d, "--entry=later", "--roles=bcj", "--subsample=prtarget", "--variant=mid3")), "--variant=mid3 is refused")
  expect_error(parse_spec(c(d, "--entry=later", "--roles=bcj", "--reach=capdist", "--variant=mid3")), "--variant=mid3 is refused")
  for (sw in c("--pairs=joiners", "--context=bilateral", "--ww2=drop", "--mu_shift=traits")) {
    expect_error(parse_spec(c(d, "--subsample=majinit", sw)), "one event / sample filter at a time", info = sw)
    expect_error(parse_spec(c(d, "--entry=later", "--roles=bcj", "--variant=mid3", sw)), "--variant=mid3 is refused", info = sw)
  }
  expect_error(parse_spec(c(d, "--entry=later", "--roles=bcj", "--joiner_set=opposing", "--subsample=prtarget")), "one event / sample filter at a time")
  expect_error(parse_spec(c(d, "--entry=later", "--roles=bcj", "--joiner_set=opposing", "--variant=mid3")), "--variant=mid3 is refused")
  expect_error(parse_spec(c(d, "--entry=later", "--roles=bcj", "--subsample=bilateral")), "pooled model M1b and no --roles")
  expect_error(parse_spec(c(d, "--subsample=bilateral")), "needs --entry=later or --events=dymid")
  expect_silent(parse_spec(c(d, "--events=dymid", "--subsample=bilateral")))
  ## models
  expect_error(cell_name("fmid", 3L, "M1b", spec_dir("--reach=capdist")), "run 0a is specified for M1c only"); expect_error(cell_name("fmid", 3L, "M2c", spec_dir("--roles=bcj", "--reach=capdist")), "specified for M1c only")
  expect_error(cell_name("fmid", 3L, "M1c", spec_dir("--subsample=bilateral")), "specified for M1b only.*M1c stops on an empty group")
  for (sb in c("prtarget", "majinit", "mininit")) { expect_error(cell_name("fmid", 3L, "M1b", spec_dir(paste0("--subsample=", sb))), "specified for M1c only", info = sb)
    expect_error(cell_name("fmid", 3L, "M2c", spec_dir("--roles=bcj", paste0("--subsample=", sb))), "specified for M1c only", info = sb) }
  for (m in c("M1b", "M1cr", "M2b", "M1a", "M1c", "M2c")) expect_error(cell_name("fmid", 3L, m, spec_mid3()), "--variant=mid3 is refused", info = m)
  for (m in c("M1b", "M2b", "M1a")) expect_error(cell_name("fmid", 3L, m, parse_spec("--variant=polrel")), "the politically relevant choice set is specified for", info = m)   # unchanged
  ## the assembly refuses the same combinations when called directly
  td <- file.path(tempfile("r0ref"), "data_entry"); write_toy_cell(td, alts_edit = toy_dir, roles = toy_roles)
  expect_error(quiet(assemble_cell("mid", 3L, "all", REACH = "capdist", data_dir = td)), "REACH = capdist needs DIRECTED")
  expect_error(quiet(assemble_cell("mid", 3L, "all", SUBSAMPLE = "prtarget", data_dir = td)), "SUBSAMPLE = prtarget needs DIRECTED")
  expect_error(quiet(assemble_cell("mid", 3L, "mid3", data_dir = td)), "VARIANT = mid3 needs DIRECTED")
  expect_error(quiet(assemble_cell("mid", 3L, "polrel", DISTANCE = "capital", MAJOR = "joint", DIRECTED = TRUE, SUBSAMPLE = "majinit", data_dir = td)), "SUBSAMPLE = majinit needs")
  expect_error(quiet(assemble_cell_spec("mid", 3L, spec_dir("--reach=capdist"), data_dir = td, model = "M1b")))
  unlink(dirname(td), recursive = TRUE)
})

test_that("defaults: names, strings and the toy assembly of the existing cells are identical to the code before the change (git aeb876b)", {
  e <- new.env(parent = globalenv())
  for (f in c("spec.R", "helpers_data.R", "helpers_model.R", "assemble_cell.R")) {
    src <- tryCatch(suppressWarnings(system2("git", c("-C", root, "show", sprintf("%s:v2/R/%s", BASE_COMMIT, f)), stdout = TRUE, stderr = FALSE)), error = function(err) NULL)
    if (is.null(src) || !length(src) || !is.null(attr(src, "status"))) skip(sprintf("git history with commit %s not available", BASE_COMMIT))
    tf <- tempfile(fileext = ".R"); writeLines(src, tf); sys.source(tf, envir = e); unlink(tf)
  }
  combos <- list(character(0), "--ame_K=1", c("--comparison=directed", "--entry=later", "--roles=bcj"), c("--comparison=directed", "--entry=later", "--roles=bcj", "--tag=long"),
                 c("--comparison=directed", "--entry=later", "--roles=bcj", "--variant=polrel"), c("--comparison=directed", "--entry=later", "--roles=bc", "--joiner_set=opposing"),
                 c("--comparison=directed", "--events=dymid", "--roles=bcj"), c("--comparison=directed", "--roles=aprime"), c("--comparison=conventional"), c("--comparison=directed", "--mu_shift=traits"),
                 c("--distance=capital", "--major=joint"), "--pairs=joiners", c("--context=coalition", "--coal_drop=top3"), "--ww2=drop")
  for (cb in combos) {
    sn <- parse_spec(cb); so <- e$parse_spec(cb)
    expect_identical(sn[names(so)], so, info = paste(cb, collapse = " "))                                   # the spec itself: the old fields unchanged, two new ones at none
    expect_identical(setdiff(names(sn), names(so)), c("reach", "subsample")); expect_identical(unlist(sn[c("reach", "subsample")], use.names = FALSE), c("none", "none"))
    expect_identical(spec_string(sn), e$spec_string(so), info = paste(cb, collapse = " ")); expect_identical(perm_part_spec_string(sn), e$perm_part_spec_string(so), info = paste(cb, collapse = " "))
    expect_identical(spec_suffix(sn), e$spec_suffix(so)); expect_identical(basename(spec_data_dir(sn, "x")), basename(e$spec_data_dir(so, "x")))
    for (m in ALL_MODELS) for (s in c("fmid", "mid")) {
      nn <- tryCatch(cell_name(s, 3L, m, cell_spec(m, sn)), error = function(err) "ERR"); no <- tryCatch(e$cell_name(s, 3L, m, e$cell_spec(m, so)), error = function(err) "ERR")
      expect_identical(nn, no, info = paste(m, s, paste(cb, collapse = " ")))
    }
  }
  ## the assembled toy cells: every field the old code returned is identical(); the new code only ADDS fields (all NULL / "none" / the variant name)
  base <- tempfile("r0base"); pr_cut <- function(a) a |> filter(!(ccode1 == 10L & ccode2 == 30L), !(ccode1 == 40L & ccode2 == 30L))
  jsets <- tibble(initiator = c(30L, 40L, 40L), target = 20L, year = 1990L, dispnum = c(2L, 7L, 7L), opp_ccode = c(20L, 20L, 30L))
  for (dd in c("data", "data_entry", "data_dymid")) write_toy_cell(file.path(base, dd), alts_edit = toy_dir, roles = toy_roles, jsets = jsets, polrel = pr_cut)
  cells <- list(list(character(0), "M1b"), list(character(0), "M1c"), list("--pairs=originators", "M1b"), list("--ww2=drop", "M1b"), list("--context=coalition", "M1b"),
                list(c("--distance=capital", "--major=joint"), "M1b"), list("--comparison=conventional", "M1c"), list("--comparison=directed", "M1b"),
                list(c("--comparison=directed", "--mu_shift=traits"), "M1b"), list(c("--comparison=directed", "--roles=aprime", "--joiner_set=opposing"), "M1c"),
                list(c("--comparison=directed", "--entry=later", "--roles=bcj"), "M1c"), list(c("--comparison=directed", "--entry=later", "--roles=bcj"), "M2c"),
                list(c("--comparison=directed", "--entry=later", "--roles=bc", "--joiner_set=opposing"), "M1c"), list(c("--comparison=directed", "--entry=later", "--roles=bcj", "--variant=polrel"), "M1c"),
                list(c("--comparison=directed", "--events=dymid", "--roles=bcj"), "M1c"), list(c("--comparison=directed", "--events=dymid", "--roles=bcj", "--variant=polrel"), "M2c"))
  for (cl in cells) for (s in c("mid", "fmid")) {
    sn <- parse_spec(cl[[1]]); so <- e$parse_spec(cl[[1]]); td <- spec_data_dir(sn, base); info <- paste(s, cl[[2]], paste(cl[[1]], collapse = " "))
    ln <- capture.output(An <- assemble_cell_spec(s, 3L, sn, data_dir = td, model = cl[[2]])); lo <- capture.output(Ao <- e$assemble_cell_spec(s, 3L, so, data_dir = td, model = cl[[2]]))
    expect_identical(An[names(Ao)], Ao, info = info)
    expect_identical(setdiff(names(An), names(Ao)), c("variant", "mid3_stats", "reach", "reach_stats", "subsample", "subsample_stats", "init_sub_stats"), info = info)
    expect_true(is.null(An$mid3_stats) && is.null(An$reach_stats) && is.null(An$subsample_stats) && is.null(An$init_sub_stats) && An$reach == "none" && An$subsample == "none", info = info)
    expect_identical(sub("^\\[[0-9:]+\\] ", "", ln), sub("^\\[[0-9:]+\\] ", "", lo), info = info)           # the log lines too (timestamps stripped)
  }
  unlink(base, recursive = TRUE)
})

## ---- (B) the history window ----------------------------------------------------------------------------------------------------------------------------
test_that("mid3_history(): onsets in [t - 3, t - 1] only, NEVER year t, both directions, every row of the events file", {
  ev <- tibble(initiator = c(1L, 3L, 5L), target = c(2L, 1L, 6L), year = c(2000L, 2000L, 1995L), in_mid = c(TRUE, FALSE, TRUE), escalation = c(FALSE, TRUE, FALSE))
  h <- mid3_history(ev, 3L); k <- paste(h$ccode1, h$ccode2, h$year)
  expect_setequal(h$year[h$ccode1 == 1L & h$ccode2 == 2L], 2001:2003)                              # t = onset + 1, + 2, + 3
  expect_false("1 2 2000" %in% k); expect_false("2 1 2000" %in% k)                                  # the onset year itself is NOT history: a window that included t fails here
  expect_false(any(h$year == 2000L)); expect_false(any(h$year == 1995L))                            # no key at any onset year of its own pair
  expect_false("1 2 2004" %in% k); expect_false("1 2 1999" %in% k)                                  # t - 4 is outside; the future is not history
  expect_setequal(h$year[h$ccode1 == 2L & h$ccode2 == 1L], 2001:2003)                              # either direction
  expect_true(all(c("3 1 2001", "1 3 2003") %in% k))                                                # a row in no sample (an escalation row) counts: all hostility levels
  expect_equal(nrow(h), 2L * 3L * 3L); expect_false(any(duplicated(k)))
  expect_setequal(mid3_history(ev, 1L)$year[1:2], 2001L); expect_equal(nrow(mid3_history(ev, 1L)), 6L)
  ## two onsets of one pair: the union of their windows, still never an onset year unless ANOTHER onset precedes it in the window
  ev2 <- tibble(initiator = c(1L, 2L), target = c(2L, 1L), year = c(2000L, 2002L)); h2 <- mid3_history(ev2, 3L)
  expect_setequal(h2$year[h2$ccode1 == 1L], 2001:2005)                                              # 2002 is history because of the 2000 onset, not its own
  expect_false(2000L %in% h2$year)
  expect_error(mid3_history(tibble(initiator = 1L, target = 1L, year = 2000L)), "initiator != events\\$target"); expect_error(mid3_history(ev, 0L))
})

## ---- (C) the mid3 assembly -----------------------------------------------------------------------------------------------------------------------------
test_that("mid3 assembly (toy): IV = history rows, each event scored against history rows + its OWN chosen target, chosen-only events dropped by role, year t excluded", {
  base <- tempfile("r0mid3"); td <- file.path(base, "data_entry")
  write_toy_cell(td, alts_edit = toy_dir, events_edit = add_history(HIST_A), roles = toy_roles)
  s3 <- spec_mid3(); sa <- spec_dir("--roles=bcj")
  out <- capture.output(A <- assemble_cell_spec("mid", 3L, s3, data_dir = td, model = "M1c")); Aa <- quiet(assemble_cell_spec("mid", 3L, sa, data_dir = td, model = "M1c"))
  expect_equal(A$variant, "mid3"); expect_equal(nrow(A$sy), 4L)                                       # every state-year has a history IV alternative here
  ## history per (i, j): 10-20, 10-40, 20-40 by an onset in 1987-1989 (both directions); 30 -> 40 by dy_mid_tm1; 10-30 NOT (its onsets are 1986 and 1990 = t)
  hist <- setNames(A$alts$mid3_hist, key_row(A, seq_len(nrow(A$alts))))
  expect_setequal(names(hist)[hist == 1L], c("10 20", "10 40", "20 10", "20 40", "30 40", "40 10", "40 20"))
  expect_false("30 10" %in% names(hist)[hist == 1L]); expect_equal(unname(hist[["10 30"]]), 0L)       # the year-t onset 10 -> 30 / 30 -> 10 (1990) is not history
  ## the IV set = the January rows with history; the stage-2 rows = history rows + the sample's chosen targets
  expect_setequal(key_row(A, which(A$alts$in_iv == 1L)), c("10 20", "20 10", "20 40", "30 40", "40 10", "40 20"))     # 10 -> 40 is stage-2 only in the toy
  expect_setequal(key_row(A, which(A$alts$in_s2 == 1L)), c("10 20", "10 40", "10 30", "20 10", "20 40", "30 20", "30 10", "40 10", "40 20"))
  expect_true(all(A$alts$in_iv == 0L | A$alts$mid3_hist == 1L)); expect_true(all(A$alts$in_s2 == 0L | A$alts$mid3_hist == 1L | A$alts$mid3_chosen == 1L))
  ## events: 30 -> 20 and 30 -> 10 have no other history alternative in stage 2 (30's only history row, 40, is IV only) and are dropped; 30 stays in stage 1
  expect_setequal(key_ev(A), c("10 20", "10 30", "40 20")); expect_equal(A$mid3_stats$n_set1, 2L); expect_equal(A$mid3_stats$n_events, 5L); expect_equal(A$mid3_stats$n_usable, 3L)
  expect_equal(A$sy$y1[A$sy$ccode == 30L], 1L); expect_equal(A$n_orphan, 1L); expect_identical(A$sy$y1, Aa$sy$y1)
  ## the chosen target is always present in its own set; set = history stage-2 rows of the state-year + the own chosen target ONLY
  sets <- split(A$alts$ccode2[A$ev_sets$alt_row], A$ev_sets$ev_idx); names(sets) <- key_ev(A)[as.integer(names(sets))]
  expect_setequal(sets[["10 20"]], c(20L, 40L))                                                       # NOT 30: the other event's chosen target without history
  expect_setequal(sets[["10 30"]], c(20L, 40L, 30L))                                                  # history rows + its own chosen target (no history)
  expect_setequal(sets[["40 20"]], c(10L, 20L))
  for (i in seq_len(nrow(A$ev))) expect_true(A$ev$alt_row[i] %in% A$ev_sets$alt_row[A$ev_sets$ev_idx == i])
  expect_equal(A$mid3_stats$ev$chosen_hist[match(c("10 20", "10 30", "40 20"), paste(A$mid3_stats$ev$ccode1, A$mid3_stats$ev$ccode2))], c(1L, 0L, 1L))
  ## by role (bcj: 10->20 bilateral, 10->30 coalition originator, 30->20 / 40->20 joiners, 30->10 bilateral): events, chosen-has-history share, size-1 share
  br <- A$mid3_stats$by_role
  expect_equal(br$group, c("bilateral", "coalition_originator", "joiner")); expect_equal(br$n_events, c(2, 1, 2)); expect_equal(br$share_chosen_hist, c(0.5, 0, 0.5))
  expect_equal(br$n_set1, c(1, 0, 1)); expect_equal(br$share_set1, c(0.5, 0, 0.5)); expect_equal(br$n_lost_with_sy, c(0L, 0L, 0L))
  expect_equal(A$group_table$n_usable, c(1L, 1L, 1L)); expect_true(any(grepl("2 events whose set holds only the chosen target are dropped from stage 2", out, fixed = TRUE)))
  ## the Stan data: ragged sets, every event with >= 2 allowed rows; the event restricted in Stan is the one whose set is smaller than its state-year's stage-2 rows
  dat <- build_stan_data(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, x_ref = A$x_ref, allow_orphan_initiations = TRUE, D_mode = 3L)
  dat <- event_alt_sets(dat, A$ev_sets)
  expect_equal(as.integer(dat$ev_alt_n), c(2L, 3L, 2L)[match(key_ev(A), c("10 20", "10 30", "40 20"))]); expect_equal(dat$n_ev_alt, 7L); expect_true(all(dat$ev_alt_n >= 2L))
  ## W, X, Z columns and the kept rows' values are those of the full-set directed cell (Z centring constants differ: other IV rows)
  expect_identical(A$W_cols, Aa$W_cols); expect_identical(A$Z_cols, Z_DIRECTED); expect_identical(A$sy[, c("ccode", "year", "y1")], Aa$sy[, c("ccode", "year", "y1")])
  unlink(base, recursive = TRUE)
})

test_that("mid3 assembly (toy): a state-year without a history alternative in its IV set leaves stage 1 with its events; left-censored history is refused", {
  base <- tempfile("r0mid3b"); td <- file.path(base, "data_entry")
  ## 10-30 history added (onset 1988); state-year 20's two history rows (10, 40) are made stage-2 only, so its IV set has no history row
  histB <- rbind(HIST_A, data.frame(initiator = 10L, target = 30L, year = 1988L))
  cutB <- function(a) { a <- toy_dir(a); a$in_iv[a$ccode1 == 20L & a$ccode2 %in% c(10L, 40L)] <- 0L; a }
  write_toy_cell(td, alts_edit = cutB, events_edit = add_history(histB), roles = toy_roles)
  out <- capture.output(A <- assemble_cell_spec("mid", 3L, spec_mid3(), data_dir = td, model = "M1c"))
  expect_equal(A$sy$ccode, c(10L, 30L, 40L)); expect_equal(A$mid3_stats$n_sy_dropped, 1L); expect_equal(A$mid3_stats$y1_dropped, 0L); expect_false(20L %in% A$sy$ccode[A$alts$sy_row])
  expect_true(any(grepl("1 of 4 panel state-years have NO history alternative in their IV set and leave stage 1", out, fixed = TRUE)))
  expect_true(any(grepl("dropped from stage 1: 1 panel state-years without a MID-history alternative", out, fixed = TRUE)))
  ## with the 10-30 history both events of 10 see {20, 30, 40}; 30 -> 20 is scored against {10, 20}; 30 -> 10 (chosen WITH history, no other history
  ## stage-2 row: 40 is IV only) is a chosen-only set and is dropped
  expect_setequal(key_ev(A), c("10 20", "10 30", "30 20", "40 20")); expect_equal(A$mid3_stats$n_set1, 1L); expect_equal(A$mid3_stats$by_role$n_set1, c(1, 0, 0))
  sets <- split(A$alts$ccode2[A$ev_sets$alt_row], A$ev_sets$ev_idx); names(sets) <- key_ev(A)[as.integer(names(sets))]
  expect_setequal(sets[["30 20"]], c(10L, 20L)); expect_setequal(sets[["10 20"]], c(20L, 30L, 40L)); expect_setequal(sets[["10 30"]], c(20L, 30L, 40L))
  unlink(base, recursive = TRUE)
  ## an INITIATING state-year dropped: 30 loses its only history (dy_mid_tm1 moved from (30, 40) to (20, 40)), so it leaves stage 1 with y1 = 1 and its
  ## two expected events are counted by role (30 -> 10 bilateral, 30 -> 20 joiner; 30 -> 40 was stage 1 only)
  td2 <- file.path(base, "data_entry")
  cutC <- function(a) { a <- toy_dir(a); a$dy_mid_tm1 <- as.integer(a$ccode1 == 20L & a$ccode2 == 40L); a }
  write_toy_cell(td2, alts_edit = cutC, events_edit = add_history(HIST_A), roles = toy_roles)
  out <- capture.output(A2 <- assemble_cell_spec("mid", 3L, spec_mid3(), data_dir = td2, model = "M1c"))
  expect_equal(A2$sy$ccode, c(10L, 20L, 40L)); expect_equal(A2$mid3_stats$n_sy_dropped, 1L); expect_equal(A2$mid3_stats$y1_dropped, 1L); expect_equal(sum(A2$sy$y1), 2L)
  expect_equal(A2$mid3_stats$by_role$n_lost_with_sy, c(1L, 0L, 1L)); expect_equal(A2$mid3_stats$n_lost_with_sy, 2L); expect_setequal(key_ev(A2), c("10 20", "10 30", "40 20"))
  expect_true(any(grepl("2 expected stage-2 events of sample mid leave with the 1 state-years", out, fixed = TRUE)))
  ## history that would be left-censored (alternatives of 1990, events file starting in 1990) is refused, not treated as "no history"
  td3 <- file.path(tempfile("r0mid3c"), "data_entry"); write_toy_cell(td3, alts_edit = toy_dir, roles = toy_roles)
  expect_error(quiet(assemble_cell_spec("mid", 3L, spec_mid3(), data_dir = td3, model = "M1c")), "min\\(alts\\$year\\) - MID3_WINDOW >= min\\(events\\$year\\)")
  unlink(c(base, dirname(td3)), recursive = TRUE)
})

## ---- (D) run 0 -----------------------------------------------------------------------------------------------------------------------------------------
test_that("run 0a (toy): reach = capdist appends the squared centred distance and its product with the initiator's M_lvl, centred over the IV rows; nothing else changes", {
  base <- tempfile("r0reach"); td <- file.path(base, "data_entry"); write_toy_cell(td, alts_edit = toy_dir, roles = toy_roles)
  A0 <- quiet(assemble_cell_spec("mid", 3L, spec_dir("--roles=bcj"), data_dir = td, model = "M1c")); out <- capture.output(A <- assemble_cell_spec("mid", 3L, spec_dir("--roles=bcj", "--reach=capdist"), data_dir = td, model = "M1c"))
  expect_identical(A$Z_cols, c(Z_DIRECTED, "log_capdist_sq", "log_capdist_x_init_M")); expect_equal(length(A$Z_cols), 11L); expect_identical(A$W_cols, W_CONVENTIONAL); expect_equal(A$reach, "capdist")
  expect_identical(A$sy, A0$sy); expect_identical(A$ev, A0$ev); expect_identical(A$alts[, names(A0$alts)], A0$alts)       # same rows, same events, every old column identical
  iv <- A0$alts$in_iv == 1L; cd <- A0$alts$log_capdist                                                # already centred over the IV rows in the base cell
  expect_lt(abs(mean(cd[iv])), 1e-12)
  sq <- cd^2; xm <- cd * A0$sy$M_lvl[A0$alts$sy_row]
  expect_equal(A$alts$log_capdist_sq, sq - mean(sq[iv]), tolerance = 1e-12); expect_equal(A$alts$log_capdist_x_init_M, xm - mean(xm[iv]), tolerance = 1e-12)
  expect_lt(abs(mean(A$alts$log_capdist_sq[iv])), 1e-10); expect_lt(abs(mean(A$alts$log_capdist_x_init_M[iv])), 1e-10)
  expect_equal(A$reach_stats$cor_lin_sq, cor(cd[iv], sq[iv]), tolerance = 1e-10); expect_gt(A$reach_stats$sd_sq, 0); expect_gt(A$reach_stats$sd_x, 0)
  ## the product varies WITHIN a choice set (an initiator-only term would cancel in the conditional logit)
  expect_true(all(tapply(A$alts$log_capdist_x_init_M, A$alts$sy_row, function(v) length(unique(v)) > 1L)))
  expect_true(any(grepl("reach = capdist (run 0a", out, fixed = TRUE)))
  dat <- build_stan_data(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, x_ref = A$x_ref, allow_orphan_initiations = TRUE, D_mode = 3L)
  expect_equal(dat$Pz, 11L); expect_false(anyNA(dat$Z))
  unlink(base, recursive = TRUE)
})

test_that("run 0b (toy): prtarget = the polrel cell's events, roles and weights (no D8 rescale) on the FULL choice set, stage 1 unchanged", {
  base <- tempfile("r0prt"); td <- file.path(base, "data_entry")
  ## 10 -> 30 and 40 -> 30 are not politically relevant. 10's two events are made ONE dispute (dispnum 4, weights 0.5 each, both co-originator events), so
  ## dropping 10 -> 30 splits a D8 group (which is NOT rescaled, as in the polrel cell); 30 -> 20 and 30 -> 10 are made one dispute (2) with weights 0.25 / 0.75 and stay whole
  ev_ed <- function(ev) { ev$dispnum[ev$initiator == 10L & ev$target == 20L] <- 4L; ev$dispnums[ev$initiator == 10L & ev$target == 20L] <- "4"
    ev$w_mid[ev$initiator == 10L] <- 0.5; ev$dispnum[ev$initiator == 30L & ev$target == 10L] <- 2L; ev$dispnums[ev$initiator == 30L & ev$target == 10L] <- "2"
    ev$w_mid[ev$initiator == 30L & ev$target == 20L] <- 0.25; ev$w_mid[ev$initiator == 30L & ev$target == 10L] <- 0.75; ev }
  rl <- toy_roles |> mutate(dispnum = ifelse(initiator == 10L & target == 20L, 4L, ifelse(initiator == 30L & target == 10L, 2L, dispnum)),
                            grp_aprime = ifelse(initiator == 10L, 1L, ifelse(initiator == 30L & target == 10L, 2L, grp_aprime)))
  al_ed <- function(a) { a <- toy_dir(a); a$polrel[(a$ccode1 == 10L & a$ccode2 == 30L) | (a$ccode1 == 40L & a$ccode2 == 30L)] <- 0L; a }
  write_toy_cell(td, alts_edit = al_ed, events_edit = ev_ed, roles = rl, polrel = function(a) a |> filter(polrel == 1L))
  A0 <- quiet(assemble_cell_spec("mid", 3L, spec_dir("--roles=bcj"), data_dir = td, model = "M1c"))
  out <- capture.output(A <- assemble_cell_spec("mid", 3L, spec_dir("--roles=bcj", "--subsample=prtarget"), data_dir = td, model = "M1c"))
  Ap <- quiet(assemble_cell_spec("mid", 3L, spec_dir("--roles=bcj", "--variant=polrel"), data_dir = td, model = "M1c"))
  expect_setequal(key_ev(A), setdiff(key_ev(A0), "10 30")); expect_setequal(key_ev(A), key_ev(Ap))    # the events of the polrel cell ...
  expect_identical(A$alts, A0$alts); expect_equal(nrow(Ap$alts), nrow(A0$alts) - 2L)                  # ... on the FULL choice set (polrel loses two rows)
  expect_identical(A$sy, A0$sy); expect_equal(sum(A$sy$y1), 3L)                                       # stage 1 unchanged: y1 not recomputed
  expect_identical(A$W_cols, A0$W_cols); expect_identical(A$Z_cols, A0$Z_cols)
  ## 0b = the polrel cell on the full choice set (audit of d60da78, PI decision 2026-11-10): identical events, roles AND weights. No D8 rescale: the
  ## split group (10, dispute 4, 1990) keeps the build's 0.5 on 10 -> 20, exactly as in the polrel cell
  o <- order(key_ev(A)); op <- order(key_ev(Ap))
  expect_identical(key_ev(A)[o], key_ev(Ap)[op]); expect_identical(A$ev$grp[o], Ap$ev$grp[op]); expect_identical(A$ev$w[o], Ap$ev$w[op])
  expect_identical(A$group_table[, c("group", "n_usable", "w_usable")], Ap$group_table[, c("group", "n_usable", "w_usable")])
  w <- setNames(A$ev$w, key_ev(A)); w0 <- setNames(A0$ev$w, key_ev(A0))
  expect_equal(unname(w[["10 20"]]), 0.5); expect_identical(w[names(w)], w0[names(w)])                 # every kept event keeps its weight of the unfiltered cell
  expect_equal(unname(w[c("30 20", "30 10")]), c(0.25, 0.75)); expect_lt(sum(w[c("10 20")]), 1)       # the split group no longer sums to 1 (logged, as in polrel)
  st <- A$subsample_stats
  expect_equal(st$labels, c("bilateral", "coalition_originator", "joiner")); expect_equal(st$n_dropped, c(0L, 1L, 0L)); expect_equal(st$n_events_rescaled, 0L); expect_equal(st$n_groups_split, 1L)
  expect_equal(st$w_kept_before, 3.5, tolerance = 1e-12); expect_identical(st$w_kept_after, st$w_kept_before)       # the five kept sample events (30 -> 40 is stage 1 only): 0.5 + 0.25 + 1 + 0.75 + 1
  expect_true(any(grepl("subsample = prtarget", out, fixed = TRUE)))
  ## weights = none: no rescale, weight 1 everywhere
  Aw <- quiet(assemble_cell_spec("mid", 3L, spec_dir("--roles=bcj", "--subsample=prtarget", "--weights=none"), data_dir = td, model = "M1c")); expect_true(all(Aw$ev$w == 1))
  unlink(base, recursive = TRUE)
})

test_that("run 0c / 0d (toy): bilateral keeps the bcj bilateral events under M1b with stage 1 unchanged; majinit / mininit cut BOTH stages and drop the constant columns", {
  base <- tempfile("r0sub"); td <- file.path(base, "data_entry"); write_toy_cell(td, alts_edit = toy_dir2, roles = toy_roles)
  A0 <- quiet(assemble_cell_spec("mid", 3L, spec_dir(), data_dir = td, model = "M1b"))
  out <- capture.output(Ab <- assemble_cell_spec("mid", 3L, spec_dir("--subsample=bilateral"), data_dir = td, model = "M1b"))
  expect_setequal(key_ev(Ab), c("10 20", "30 10"))                                                    # aprime 0 (alone originator) events with a stage-2 row; 30 -> 40 is stage 1 only
  expect_identical(Ab$sy, A0$sy); expect_identical(Ab$alts, A0$alts); expect_identical(Ab$W_cols, W_CONVENTIONAL); expect_identical(Ab$Z_cols, Z_DIRECTED)
  expect_equal(sum(Ab$sy$y1), 3L); expect_equal(Ab$n_orphan, 1L)                                      # 40 keeps y1 = 1 without a stage-2 event
  expect_equal(Ab$subsample_stats$n_dropped, c(0L, 1L, 2L)); expect_equal(Ab$subsample_stats$n_groups_split, 0L); expect_equal(Ab$subsample_stats$n_events_rescaled, 0L)
  w0 <- setNames(A0$ev$w, key_ev(A0)); expect_equal(setNames(Ab$ev$w, key_ev(Ab))[c("10 20", "30 10")], w0[c("10 20", "30 10")])
  ## 0d majors (10, 40): the state-years, the events, W without cowmaj_tm1, Z without the major product
  out <- capture.output(Am <- assemble_cell_spec("mid", 3L, spec_dir("--roles=bcj", "--subsample=majinit"), data_dir = td, model = "M1c"))
  A1 <- quiet(assemble_cell_spec("mid", 3L, spec_dir("--roles=bcj"), data_dir = td, model = "M1c"))
  expect_equal(Am$sy$ccode, c(10L, 40L)); expect_equal(sum(Am$sy$y1), 2L); expect_setequal(key_ev(Am), c("10 20", "10 30", "40 20")); expect_setequal(unique(Am$sy$ccode[Am$alts$sy_row]), c(10L, 40L))
  expect_identical(Am$W_cols, setdiff(W_CONVENTIONAL, "cowmaj_tm1")); expect_identical(Am$Z_cols, setdiff(Z_DIRECTED, "tgt_major_x_init_major")); expect_equal(length(Am$Z_cols), 8L)
  expect_setequal(Am$alts$ccode2[Am$alts$sy_row == 1L], c(20L, 30L, 40L))                              # targets of either status stay in a kept choice set
  expect_equal(Am$init_sub_stats$n_kept, 2L); expect_equal(Am$init_sub_stats$n_panel, 4L); expect_equal(Am$init_sub_stats$y1_kept, 2L); expect_equal(Am$init_sub_stats$y1_panel, 3L)
  expect_equal(Am$subsample_stats$n_dropped, c(2L, 0L, 1L)); expect_equal(Am$subsample_stats$n_events_rescaled, 0L)
  wm <- setNames(Am$ev$w, key_ev(Am)); w1 <- setNames(A1$ev$w, key_ev(A1)); expect_equal(wm[names(wm)], w1[names(wm)])       # D8 groups kept whole: weights as built
  expect_true(any(grepl("cowmaj_tm1 is constant (1)", out, fixed = TRUE))); expect_true(any(grepl("tgt_major_x_init_major equals tgt_major", out, fixed = TRUE)))
  dat <- build_stan_data(Am$sy, Am$alts, Am$ev, Am$W_cols, Am$X_cols, Am$Z_cols, joint = 1L, x_ref = Am$x_ref, allow_orphan_initiations = TRUE, D_mode = 3L)
  expect_equal(dat$N1, 2L); expect_equal(dat$P1, 7L); expect_equal(dat$Pz, 8L); expect_equal(sum(dat$y1), 2L)
  ## 0d minors (20, 30): the toy's only minor initiator (30, a non-democracy with two democratic stage-2 targets) leaves no within-set variation in the
  ## democracy columns, so the directed identification check stops the run, as it must; the real cells are covered by the dry runs / block (F)
  expect_error(quiet(assemble_cell_spec("mid", 3L, spec_dir("--roles=bcj", "--subsample=mininit"), data_dir = td, model = "M1c")), "a target column does not vary within any choice set")
  unlink(base, recursive = TRUE)
})

## ---- (E) end to end ------------------------------------------------------------------------------------------------------------------------------------
test_that("02_fit_model.R --dry_run end to end on a toy cell for every new cell type; refusals at the command line", {
  tmp <- tempfile("r0e2e"); dir.create(file.path(tmp, "v2"), recursive = TRUE); file.copy(file.path(root, "R"), file.path(tmp, "v2"), recursive = TRUE)
  ok <- file.symlink(file.path(root, "stan"), file.path(tmp, "v2", "stan")); if (!ok) file.copy(file.path(root, "stan"), file.path(tmp, "v2"), recursive = TRUE)
  write_toy_cell(file.path(tmp, "v2", "data_entry"), alts_edit = toy_dir2, events_edit = add_history(HIST_A), roles = toy_roles)   # (10, 40) has history through its 1989 onset already
  dir.create(file.path(tmp, "v2", "output")); dir.create(file.path(tmp, "v2", "logs")); v2 <- file.path(tmp, "v2"); out_dir <- file.path(v2, "output")
  run <- function(...) { out <- suppressWarnings(system2("Rscript", c(file.path(v2, "R", "02_fit_model.R"), "--sample=mid", "--lag=3", "--comparison=directed", "--entry=later", "--dry_run=TRUE", ...), stdout = TRUE, stderr = TRUE,
                                                 env = "V2_ROOT=")); list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  rd <- function(cell) readRDS(file.path(out_dir, paste0("standata_", cell, ".rds")))
  ## mid3 is refused at the command line for both models, with the reason, and writes nothing
  for (m in list(c("--model=M1c"), c("--model=M2c", "--ame_K=2"))) { r <- run(m, "--roles=bcj", "--variant=mid3"); expect_true(r$status != 0L); expect_true(any(grepl("--variant=mid3 is refused (audit P1-1", r$out, fixed = TRUE))) }
  expect_length(list.files(out_dir, pattern = "mid3"), 0L)
  r <- run("--model=M1c", "--roles=bcj", "--reach=capdist"); expect_equal(r$status, 0L, info = paste(tail(r$out, 15), collapse = "\n"))
  o <- rd("mid_lag3_all_pr_M1c_cdir_reach_bcj_ent"); expect_equal(o$dat$Pz, 11L); expect_identical(o$Z_cols, c(Z_DIRECTED, Z_REACH)); expect_false(is.null(o$reach_stats)); expect_null(o$mid3_stats)
  r <- run("--model=M1c", "--roles=bcj", "--subsample=prtarget"); expect_equal(r$status, 0L, info = paste(tail(r$out, 15), collapse = "\n")); expect_equal(rd("mid_lag3_all_pr_M1c_cdir_bcj_sprt_ent")$subsample_stats$n_events_rescaled, 0L)
  r <- run("--model=M1b", "--subsample=bilateral"); expect_equal(r$status, 0L, info = paste(tail(r$out, 15), collapse = "\n"))
  o <- rd("mid_lag3_all_pr_M1b_cdir_sbil_ent"); expect_equal(o$dat$E, 2L); expect_equal(sum(o$dat$y1), 3L); expect_true("roles_aprime.rds" %in% o$manifest$file)
  r <- run("--model=M1c", "--roles=bcj", "--subsample=majinit"); expect_equal(r$status, 0L, info = paste(tail(r$out, 15), collapse = "\n"))
  o <- rd("mid_lag3_all_pr_M1c_cdir_bcj_smaj_ent"); expect_equal(o$dat$N1, 2L); expect_equal(o$dat$P1, 7L); expect_equal(o$dat$Pz, 8L); expect_identical(o$W_cols, setdiff(W_CONVENTIONAL, "cowmaj_tm1"))
  ## an existing cell saves exactly the fields it saved before (no run-0 field)
  r <- run("--model=M1c", "--roles=bcj"); expect_equal(r$status, 0L); o <- rd("mid_lag3_all_pr_M1c_cdir_bcj_ent")
  expect_false(any(c("mid3_stats", "reach_stats", "subsample_stats", "init_sub_stats") %in% names(o))); expect_equal(o$dat$Pz, 9L); expect_identical(o$Z_cols, Z_DIRECTED)
  ## refusals
  for (bad in list(c("--model=M1b", "--reach=capdist"), c("--model=M1c", "--roles=bcj", "--subsample=bilateral"), c("--model=M1c", "--subsample=bilateral"), c("--model=M1b", "--variant=mid3"),
                   c("--model=M1c", "--roles=bcj", "--subsample=prtarget", "--variant=polrel"), c("--model=M1c", "--roles=bcj", "--reach")))
    expect_true(run(bad)$status != 0L, info = paste(bad, collapse = " "))
  unlink(tmp, recursive = TRUE)
})

test_that("slurm: 11_cset_grid.sbatch and 12_run0.sbatch issue valid cells (stub Rscript); submit_cset_grid.sh chains the M2c array afterok on the M1c array (stub sbatch)", {
  tmp <- tempfile("r0slurm"); dir.create(file.path(tmp, "v2"), recursive = TRUE); dir.create(file.path(tmp, "bin"))
  for (d in c("slurm", "env")) file.copy(file.path(root, d), file.path(tmp, "v2"), recursive = TRUE)
  dir.create(file.path(tmp, "v2", "R")); file.copy(file.path(root, "R", "01_build_data.R"), file.path(tmp, "v2", "R"))
  rec <- file.path(tmp, "calls.txt"); ctr <- file.path(tmp, "counter")
  writeLines(c("#!/bin/bash", sprintf("echo \"$*\" >> %s", rec)), file.path(tmp, "bin", "Rscript")); writeLines(c("#!/bin/bash", "exit 0"), file.path(tmp, "bin", "module"))
  writeLines(c("#!/bin/bash", sprintf("n=$(cat %s 2>/dev/null || echo 1000); n=$((n + 1)); echo $n > %s", ctr, ctr), sprintf("echo \"sbatch job=$n args: $*\" >> %s", rec), "echo \"$n;cluster\""), file.path(tmp, "bin", "sbatch"))
  for (b in c("Rscript", "module", "sbatch")) Sys.chmod(file.path(tmp, "bin", b), "0755")
  env0 <- c(sprintf("PATH=%s:%s", file.path(tmp, "bin"), Sys.getenv("PATH")), sprintf("SLURM_SUBMIT_DIR=%s", tmp), "V2_ROOT=")
  task <- function(file, i, env = character(0)) { unlink(rec); out <- suppressWarnings(system2("bash", file.path(tmp, "v2", "slurm", file), stdout = TRUE, stderr = TRUE, env = c(env0, sprintf("SLURM_ARRAY_TASK_ID=%d", i), env)))
    list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out, call = if (file.exists(rec)) readLines(rec) else character(0)) }
  cell_of <- function(call) { a <- strsplit(sub("^v2/R/02_fit_model.R +", "", call), " +")[[1]]; sp <- parse_spec(a); m <- get_arg(a, "model", ""); cell_name(get_arg(a, "sample", ""), as.integer(get_arg(a, "lag", 0)), m, cell_spec(m, sp)) }
  ## the choice-set grid: 8 cells ({all, polrel} x {pairs, dymid} x {M1c, M2c K2}; mid3 dropped 2026-11-10), dry run (no init fit needed)
  want11 <- c("fmid_lag3_all_pr_M1c_cdir_bcj_ent_tlong", "fmid_lag3_polrel_pr_M1c_cdir_bcj_ent_tlong", "fmid_lag3_all_pr_M1c_cdir_bcj_dy_tlong", "fmid_lag3_polrel_pr_M1c_cdir_bcj_dy_tlong",
              "fmid_lag3_all_pr_M2c_K2_cdir_bcj_ent_tlong", "fmid_lag3_polrel_pr_M2c_K2_cdir_bcj_ent_tlong", "fmid_lag3_all_pr_M2c_K2_cdir_bcj_dy_tlong", "fmid_lag3_polrel_pr_M2c_K2_cdir_bcj_dy_tlong")
  for (i in 0:7) { r <- task("11_cset_grid.sbatch", i, "DRY_RUN=TRUE"); expect_equal(r$status, 0L, info = paste(i, paste(r$out, collapse = "\n"))); expect_length(r$call, 1L)
    expect_equal(cell_of(r$call), want11[i + 1L], info = i); expect_match(r$call, "--dry_run=TRUE", fixed = TRUE); expect_false(grepl("mid3", r$call))
    expect_equal(grepl("--iter=6000 --warmup=3000 --adapt_delta=0.9", r$call, fixed = TRUE), i <= 3L, info = i); expect_equal(grepl("--ame_K=2", r$call, fixed = TRUE), i >= 4L, info = i) }
  ## the pairs cells are the run-7 cells (their names are the run-7 names): not in the default array. Open-account convention of the repository: --partition=basic
  hdr <- readLines(file.path(root, "slurm", "11_cset_grid.sbatch")); expect_true("#SBATCH --array=2-3" %in% hdr); expect_true(all(c("#SBATCH --account=open", "#SBATCH --requeue", "#SBATCH --time=48:00:00") %in% hdr))
  expect_true(any(grepl("^#SBATCH --partition=basic( |$)", hdr))); expect_false(any(grepl("partition=open", hdr))); expect_false(any(grepl("vr=mid3", hdr)))
  expect_identical(grep("^#SBATCH --(account|partition|requeue)", hdr, value = TRUE)[1:3], grep("^#SBATCH --(account|partition|requeue)", readLines(file.path(root, "slurm", "01_build.sbatch")), value = TRUE)[1:3])
  ## a real (non-dry) M2c task: refused on the open account, and refused without its M1c fit on the paid account; an M1c task runs
  r <- task("11_cset_grid.sbatch", 6L); expect_true(r$status != 0L); expect_true(any(grepl("M2c cell on the open account", r$out))); expect_length(r$call, 0L)
  r <- task("11_cset_grid.sbatch", 6L, "SLURM_JOB_ACCOUNT=jfe4_cr_default"); expect_true(r$status != 0L); expect_true(any(grepl("needs v2/output/fit_fmid_lag3_all_pr_M1c_cdir_bcj_dy_tlong.rds", r$out, fixed = TRUE)))
  dir.create(file.path(tmp, "v2", "output"), showWarnings = FALSE); file.create(file.path(tmp, "v2", "output", "fit_fmid_lag3_all_pr_M1c_cdir_bcj_dy_tlong.rds"))
  r <- task("11_cset_grid.sbatch", 6L, "SLURM_JOB_ACCOUNT=jfe4_cr_default"); expect_equal(r$status, 0L, info = paste(r$out, collapse = "\n")); expect_equal(cell_of(r$call), want11[7])
  r <- task("11_cset_grid.sbatch", 2L); expect_equal(r$status, 0L); expect_match(r$call, "--dry_run=FALSE", fixed = TRUE)
  expect_true(task("11_cset_grid.sbatch", 8L)$status != 0L)
  ## run 0: 20 cells
  sw <- c("reach", "sprt", "sbil", "smaj", "smin")
  for (i in 0:19) { r <- task("12_run0.sbatch", i, "DRY_RUN=TRUE"); expect_equal(r$status, 0L, info = paste(i, paste(r$out, collapse = "\n"))); expect_length(r$call, 1L)
    s <- c("fmid", "mid")[i %% 2L + 1L]; ds <- c("_ent", "_dy")[(i %% 4L) %/% 2L + 1L]; k <- i %/% 4L + 1L
    want <- switch(sw[k], reach = sprintf("%s_lag3_all_pr_M1c_cdir_reach_bcj%s_tlong", s, ds), sbil = sprintf("%s_lag3_all_pr_M1b_cdir_sbil%s_tlong", s, ds), sprintf("%s_lag3_all_pr_M1c_cdir_bcj_%s%s_tlong", s, sw[k], ds))
    expect_equal(cell_of(r$call), want, info = i); expect_match(r$call, "--iter=6000 --warmup=3000 --adapt_delta=0.9", fixed = TRUE) }
  expect_true(task("12_run0.sbatch", 20L)$status != 0L)
  hdr <- readLines(file.path(root, "slurm", "12_run0.sbatch")); expect_true(all(c("#SBATCH --account=open", "#SBATCH --requeue", "#SBATCH --array=0-19") %in% hdr))
  expect_true(any(grepl("^#SBATCH --partition=basic( |$)", hdr))); expect_false(any(grepl("partition=open", hdr)))
  ## the submit helper: the M1c array on the header's open account, the M2c array on the paid allocation, afterok on it
  unlink(rec); old <- setwd(tmp); on.exit(setwd(old), add = TRUE)
  out <- suppressWarnings(system2("bash", "v2/slurm/submit_cset_grid.sh", stdout = TRUE, stderr = TRUE, env = env0)); expect_null(attr(out, "status"), info = paste(out, collapse = "\n"))
  calls <- readLines(rec); expect_length(calls, 2L)
  expect_match(calls[1], "^sbatch job=1001 args: --parsable --array=2-3 --export=ALL v2/slurm/11_cset_grid.sbatch$")
  expect_match(calls[2], "^sbatch job=1002 args: --parsable --array=6-7 --export=ALL --account=jfe4_cr_default --partition=basic --time=7-00:00:00 --dependency=afterok:1001 --kill-on-invalid-dep=yes v2/slurm/11_cset_grid.sbatch$")
  out <- suppressWarnings(system2("bash", c("v2/slurm/submit_cset_grid.sh", "--array=0-3"), stdout = TRUE, stderr = TRUE, env = env0)); expect_false(is.null(attr(out, "status")))
  out <- suppressWarnings(system2("bash", "v2/slurm/submit_cset_grid.sh", stdout = TRUE, stderr = TRUE, env = c(env0, "A_JOB=777", "M2_IDX=4-7"))); calls <- readLines(rec); expect_length(calls, 3L)
  expect_match(calls[3], "--array=4-7 .*--dependency=afterok:777 --kill-on-invalid-dep=yes")
  unlink(tmp, recursive = TRUE)
})

## ---- (F) the real builds ---------------------------------------------------------------------------------------------------------------------------------
test_that("real lag-3 builds: mid3 history recomputed independently, never a year-t onset, chosen targets present; prtarget = the events polrel drops", {
  for (dd in c("data_entry", "data_dymid")) {
    td <- file.path(root, dd); if (!file.exists(file.path(td, "manifest.csv"))) skip(paste("no build in", td))
    da <- if (dd == "data_entry") "--entry=later" else "--events=dymid"
    A <- quiet(assemble_cell_spec("fmid", 3L, spec_mid3(da), data_dir = td, model = "M1c"))             # the refused construction, called directly (unit test)
    ev <- readRDS(file.path(td, "events.rds"))
    ## independent recomputation: for every IV row, an onset of the unordered pair with 1 <= t - onset <= 3, or dy_mid_tm1 = 1
    on <- unique(data.frame(a = pmin(ev$initiator, ev$target), b = pmax(ev$initiator, ev$target), onset = ev$year))
    rows <- data.frame(a = pmin(A$sy$ccode[A$alts$sy_row], A$alts$ccode2), b = pmax(A$sy$ccode[A$alts$sy_row], A$alts$ccode2), t = A$sy$year[A$alts$sy_row], r = seq_len(nrow(A$alts)))
    j <- merge(rows, on, by = c("a", "b")); gap <- j$t - j$onset
    has_onset <- rows$r %in% j$r[gap >= 1L & gap <= 3L]
    expect_identical(as.integer(has_onset), as.integer(A$alts$mid3_onset), info = dd); expect_identical(A$alts$mid3_hist, as.integer(has_onset | A$alts$dy_mid_tm1 == 1L), info = dd)
    expect_true(all(A$alts$mid3_hist[A$alts$in_iv == 1L] == 1L), info = dd)
    ## a pair whose ONLY onsets are in year t (or later, or before t - 3) and that was not in a MID at t - 1 never has history
    only_t <- !rows$r %in% j$r[gap >= 1L & gap <= 3L] & rows$r %in% j$r[gap == 0L] & A$alts$dy_mid_tm1 == 0L
    expect_gt(sum(only_t), 0L); expect_true(all(A$alts$mid3_hist[only_t] == 0L), info = dd); expect_true(all(A$alts$in_iv[only_t] == 0L), info = dd)
    ## every event: its chosen target in its set, at least one other alternative, every other allowed row with history
    n_set <- tabulate(A$ev_sets$ev_idx, nbins = nrow(A$ev)); expect_true(all(n_set >= 2L), info = dd)
    expect_equal(nrow(anti_join(data.frame(ev_idx = seq_len(nrow(A$ev)), alt_row = A$ev$alt_row), A$ev_sets, by = c("ev_idx", "alt_row"))), 0L, info = dd)
    expect_true(all(A$alts$mid3_hist[A$ev_sets$alt_row] == 1L | A$ev_sets$alt_row == A$ev$alt_row[A$ev_sets$ev_idx]), info = dd)
    expect_equal(A$mid3_stats$n_usable + A$mid3_stats$n_set1, A$mid3_stats$n_events, info = dd); expect_true(all(A$sy$y1[A$ev$sy_row] == 1L))
    ## prtarget: the events of the polrel cell, on the full choice set
    Ap <- quiet(assemble_cell_spec("fmid", 3L, parse_spec(c("--comparison=directed", da, "--roles=bcj", "--variant=polrel")), data_dir = td, model = "M1c"))
    As <- quiet(assemble_cell_spec("fmid", 3L, parse_spec(c("--comparison=directed", da, "--roles=bcj", "--subsample=prtarget")), data_dir = td, model = "M1c"))
    expect_false(any(duplicated(key_evy(As)))); expect_setequal(key_evy(As), key_evy(Ap)); expect_identical(As$ev$grp[order(key_evy(As))], Ap$ev$grp[order(key_evy(Ap))], info = dd); expect_gt(nrow(As$alts), nrow(Ap$alts))
    expect_identical(As$ev$w[order(key_evy(As))], Ap$ev$w[order(key_evy(Ap))], info = dd)             # no D8 rescale: the polrel cell's weights exactly
    expect_identical(As$group_table[, c("group", "n_usable", "w_usable")], Ap$group_table[, c("group", "n_usable", "w_usable")], info = dd)
    if (dd == "data_entry") expect_equal(round(As$group_table$w_usable[3], 1), 238.9)                    # pairs fmid joiner weight = the polrel cell's             # the rescale only raises weights
    rm(A, Ap, As, ev, rows, j); gc()
    if (dd == "data_entry") {
      ## 0d on the real cell: minor-power initiators only, in both stages
      An <- quiet(assemble_cell_spec("fmid", 3L, parse_spec(c("--comparison=directed", da, "--roles=bcj", "--subsample=mininit")), data_dir = td, model = "M1c"))
      expect_true(all(An$sy$cowmaj_tm1 == 0L)); expect_true(all(An$alts$major1 == 0L)); expect_true(any(An$alts$major2 == 1L)); expect_true(all(An$sy$y1[An$ev$sy_row] == 1L))
      expect_identical(An$W_cols, setdiff(W_CONVENTIONAL, "cowmaj_tm1")); expect_identical(An$Z_cols, setdiff(Z_DIRECTED, "tgt_major_x_init_major")); expect_lt(nrow(An$sy), 14597L)
      rm(An); gc()
    }
  }
})
