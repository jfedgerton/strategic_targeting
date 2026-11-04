## The joiner entry direction (switch entry = keep | later, suffix _ent, build variant v2/data_entry) and the role codings bc / bcj (suffixes _bc /
## _bcj; need entry = later), PI decisions of 2026-11-04.
## (a) entry_later_direction() on a hand-built dispute: originator vs originator, originator vs joiner, joiner vs joiner by entry date, ties (also at
## mixed resolution); the collapse with and without the entry columns; (b) the switches, suffixes, data directory and refusals; strings and names
## unchanged at the defaults; (c) the bc / bcj mapping at assembly, joiner_set under bc / bcj restricts exactly the joiner events; (d) 02_fit_model.R
## end to end on a toy entry cell (M1c bc / bcj, M1cr bcj, M2c from the matching M1c); (e) the real entry build (skipped when absent).
suppressPackageStartupMessages({ library(testthat); library(dplyr) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "assemble_cell.R"))

## events after the direction rule (Side A initiators here), dispute 1 in 2000. A1, B1 originators; joiners with entry months: B2 May, A2 August,
## A3 May (ties B2), A4 March (earlier than B2), A5 day-dated 2000-05-20 (ties B2's month-resolution May at the coarser resolution)
ev_e <- function(i, j, o1, o2) tibble(dispnum = 1L, initiator = as.integer(i), target = as.integer(j), year = 2000L, sidea1 = 1L, sidea2 = 0L, orig1 = as.integer(o1), orig2 = as.integer(o2),
                                      fatality1 = 1L, fatality2 = 2L, hostlev1 = 3L, hostlev2 = 4L)
psd_e <- function(cc, res, m, d = NA_integer_) { k <- date_keys(res, 2000L, m, d)
  tibble(dispnum = 1L, ccode = as.integer(cc), spell = 1L, spell_first = 2000L, spell_last = 2000L, entry_res = as.integer(res), entry_y = 2000L, entry_mkey = k$mkey, entry_dkey = k$dkey) }

test_that("entry_later_direction(): the PI's rule on a hand-built dispute", {
  ev <- bind_rows(ev_e(11, 21, 1, 1), ev_e(11, 22, 1, 0), ev_e(12, 21, 0, 1), ev_e(12, 22, 0, 0), ev_e(13, 22, 0, 0), ev_e(14, 22, 0, 0), ev_e(15, 22, 0, 0))
  psd <- bind_rows(psd_e(11, 2L, 1L), psd_e(21, 2L, 1L), psd_e(22, 2L, 5L), psd_e(12, 2L, 8L), psd_e(13, 2L, 5L), psd_e(14, 2L, 3L), psd_e(15, 1L, 5L, 20L))
  r <- entry_later_direction(ev, psd)
  key <- paste(r$initiator, r$target)
  expect_equal(key, c("11 21", "22 11", "12 21", "12 22", "13 22", "22 14", "15 22"))
  expect_equal(r$entry_rule, c("orig_orig", "joiner_orig", "joiner_orig", "joiner_joiner_later", "joiner_joiner_tie", "joiner_joiner_later", "joiner_joiner_tie"))
  expect_equal(r$entry_reversed, c(FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE))
  ## a reversed event swaps every pair column; no late-target event remains
  expect_equal(unlist(r[2, c("sidea1", "sidea2", "orig1", "orig2", "fatality1", "fatality2", "hostlev1", "hostlev2")]), c(sidea1 = 0, sidea2 = 1, orig1 = 0, orig2 = 1, fatality1 = 2, fatality2 = 1, hostlev1 = 4, hostlev2 = 3))
  expect_false(any(r$orig1 == 1L & r$orig2 == 0L))
  d <- attr(r, "diagnostics"); expect_equal(d$n_reversed, 2L); expect_equal(d$n_reversed_joiner_orig, 1L); expect_equal(d$n_reversed_jj, 1L); expect_equal(d$n_spell_fallback, 0L)
  ## a participant without a spell in its dispute fails loudly
  expect_error(entry_later_direction(ev, psd |> filter(ccode != 22L)), "without a spell")
  ## the collapse: unchanged without the entry columns; with them entry_reversed is ORed and entry_rule pasted
  expect_true(all(c("entry_reversed", "entry_rule") %in% names(r)))
})

test_that("spec: entry = keep | later (suffix _ent, v2/data_entry), roles bc / bcj (suffixes _bc / _bcj) need entry = later; refusals; defaults unchanged", {
  expect_equal(SPEC_CHOICES$entry, c("keep", "later")); expect_equal(SPEC_SUFFIX$entry, c(later = "_ent")); expect_equal(SPEC_SUFFIX$roles, c(aprime = "_aprime", bc = "_bc", bcj = "_bcj"))
  expect_equal(which(names(SPEC_SUFFIX) == "entry"), which(names(SPEC_SUFFIX) == "direction") + 1L)
  el <- parse_spec("--entry=later")
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec(c("--comparison=directed", "--entry=later"))), "fmid_lag3_all_pr_M1b_cdir_ent")
  b <- parse_spec(c("--comparison=directed", "--roles=bc", "--entry=later")); bj <- parse_spec(c("--comparison=directed", "--roles=bcj", "--entry=later", "--joiner_set=opposing"))
  expect_equal(cell_name("mid", 3L, "M1c", b), "mid_lag3_all_pr_M1c_cdir_bc_ent"); expect_equal(cell_name("mid", 3L, "M1cr", bj), "mid_lag3_all_pr_M1cr_cdir_bcj_jset_ent")
  expect_equal(cell_name("fmid", 3L, "M2c", parse_spec(c("--comparison=directed", "--roles=bcj", "--entry=later", "--ame_K=1", "--tag=long"))), "fmid_lag3_all_pr_M2c_K1_cdir_bcj_ent_tlong")
  expect_equal(m1_cell_for_m2("fmid", 3L, "M2c", parse_spec(c("--comparison=directed", "--roles=bc", "--entry=later", "--joiner_set=opposing", "--ame_K=2"))), "fmid_lag3_all_pr_M1c_cdir_bc_jset_ent")
  expect_equal(basename(spec_data_dir(el, root)), "data_entry"); expect_equal(basename(spec_data_dir(parse_spec(character(0)), root)), "data")
  expect_equal(stage2_group_labels("bc"), c("bilateral", "coalition")); expect_equal(stage2_group_labels("bcj"), c("bilateral", "coalition_originator", "joiner"))
  expect_error(parse_spec("--roles=bc"), "--roles=bc needs --entry=later", fixed = TRUE); expect_error(parse_spec(c("--roles=bcj", "--entry=keep")), "needs --entry=later")
  expect_error(parse_spec(c("--entry=later", "--direction=sideA")), "cowwar direction only")
  expect_error(parse_spec("--entry=yes"), "--entry must be one of keep | later", fixed = TRUE)
  for (m in c("M1b", "M2b")) expect_error(cell_name("fmid", 3L, m, b), "the role coding applies", info = m)
  expect_false(grepl("entry", spec_string(parse_spec(character(0))), fixed = TRUE)); expect_true(grepl("entry=later", spec_string(el), fixed = TRUE))
  expect_false(grepl("entry", perm_part_spec_string(parse_spec("--roles=aprime")), fixed = TRUE)); expect_true(grepl("entry=later", perm_part_spec_string(el), fixed = TRUE))
})

## toy entry cell: no late-target event; groups (aprime) 10 -> 20 alone, 10 -> 30 co-originator, 30 -> 20 joiner, 30 -> 40 alone (fighting at ref),
## 30 -> 10 alone, 40 -> 20 joiner
toy_roles_e <- tibble(initiator = c(10L, 10L, 30L, 30L, 30L, 40L), target = c(20L, 30L, 20L, 40L, 10L, 20L), year = 1990L, dispnum = c(1L, 4L, 2L, 3L, 6L, 7L), grp_aprime = c(0L, 1L, 2L, 0L, 0L, 2L))
toy_jsets_e <- tibble(initiator = c(30L, 40L, 40L), target = 20L, year = 1990L, dispnum = c(2L, 7L, 7L), opp_ccode = c(20L, 20L, 30L))
toy_dirstat <- function(a) a |> mutate(dem1 = as.integer(ccode1 %in% c(10L, 20L)), dem2 = as.integer(ccode2 %in% c(10L, 20L)),
                                       major1 = as.integer(ccode1 %in% c(10L, 40L)), major2 = as.integer(ccode2 %in% c(10L, 40L)))

test_that("assembly: bc / bcj map the aprime groups (cross-tab logged); joiner_set restricts exactly the joiner events (bc: the joiner subset of group 1); late-target events refused", {
  base <- tempfile("entrytoy"); td <- file.path(base, "data_entry"); write_toy_cell(td, roles = toy_roles_e, jsets = toy_jsets_e)
  sp <- function(...) parse_spec(c("--entry=later", ...))
  Aa <- quiet(assemble_cell_spec("mid", 3L, sp("--roles=aprime"), data_dir = td, model = "M1c"))
  Ab <- quiet(assemble_cell_spec("mid", 3L, sp("--roles=bc"), data_dir = td, model = "M1c")); Aj <- quiet(assemble_cell_spec("mid", 3L, sp("--roles=bcj"), data_dir = td, model = "M1c"))
  expect_equal(Ab$group_table$group, BC_GROUPS); expect_equal(Aj$group_table$group, BCJ_GROUPS)
  expect_equal(Ab$ev$grp, c(0L, 1L, 1L, NA)[Aa$ev$grp + 1L]); expect_equal(Aj$ev$grp, Aa$ev$grp)
  for (nm in c("sy", "alts", "W_cols", "Z_cols")) { expect_identical(Ab[[nm]], Aa[[nm]], info = nm); expect_identical(Aj[[nm]], Aa[[nm]], info = nm) }
  expect_equal(Ab$ev_joiner, Aa$ev$grp == 2L); expect_equal(Aj$ev_joiner, Aa$ev$grp == 2L)
  out <- capture.output(invisible(assemble_cell_spec("mid", 3L, sp("--roles=bc"), data_dir = td, model = "M1c")))
  expect_true(any(grepl("roles = bc: cross-tab of the 6 sample events, aprime 4 groups (rows) x bc 2 groups (columns)", out, fixed = TRUE)))
  for (r in c("bc", "bcj")) {
    A <- quiet(assemble_cell_spec("mid", 3L, sp(paste0("--roles=", r), "--joiner_set=opposing"), data_dir = td, model = "M1c"))
    expect_setequal(unique(A$ev_sets$ev_idx), which(A$ev_joiner))                                  # exactly the joiner events
    if (r == "bc") expect_true(any(A$ev$grp == 1L & !A$ev_joiner))                                 # a co-originator event in group 1 stays unrestricted
    expect_equal(A$jset_stats$n_joiner_events, 2L)
  }
  ## late-target events (an entry = keep coding) are refused under bc / bcj; the data directory must be data_entry
  td2 <- file.path(tempfile("entrytoy"), "data_entry"); write_toy_cell(td2, roles = toy_roles_e |> mutate(grp_aprime = ifelse(initiator == 30L & target == 10L, 3L, grp_aprime)))
  expect_error(quiet(assemble_cell_spec("mid", 3L, sp("--roles=bcj"), data_dir = td2, model = "M1c")), "aprime group 3 (target joined late)", fixed = TRUE)
  td3 <- tempfile("toycell"); write_toy_cell(td3, roles = toy_roles_e)
  expect_error(quiet(assemble_cell_spec("mid", 3L, sp("--roles=bc"), data_dir = td3, model = "M1c")))
  unlink(c(base, dirname(td2), td3), recursive = TRUE)
})

test_that("02_fit_model.R end to end on a toy entry cell: M1c bc / bcj (with joiner_set), M1cr bcj, M2c bcj from the matching M1c, M1b; refusals", {
  for (m in c("m1_context", "m1_context_rho", "m2_context", "m1_latent_additive")) if (!file.exists(file.path(root, "stan", "build", paste0(m, ".rds")))) skip(paste("compiled model missing:", m))
  tmp <- tempfile("smoke_ent_"); dir.create(file.path(tmp, "v2"), recursive = TRUE)
  file.copy(file.path(root, "R"), file.path(tmp, "v2"), recursive = TRUE)
  ok <- file.symlink(file.path(root, "stan"), file.path(tmp, "v2", "stan")); if (!ok) file.copy(file.path(root, "stan"), file.path(tmp, "v2"), recursive = TRUE)
  write_toy_cell(file.path(tmp, "v2", "data_entry"), alts_edit = toy_dirstat, roles = toy_roles_e, jsets = toy_jsets_e)
  dir.create(file.path(tmp, "v2", "output")); dir.create(file.path(tmp, "v2", "logs"))
  v2 <- file.path(tmp, "v2"); out_dir <- file.path(v2, "output")
  run <- function(script, ...) { out <- suppressWarnings(system2("Rscript", c(file.path(v2, "R", script), ...), stdout = TRUE, stderr = TRUE)); list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  base <- c("--sample=mid", "--lag=3", "--comparison=directed", "--entry=later"); short <- c("--chains=2", "--iter=60", "--warmup=30", "--max_treedepth=6")
  ## refusals: bc without entry = later (and the keep data directory does not exist here), M2c without its M1c fit
  expect_true(run("02_fit_model.R", "--sample=mid", "--lag=3", "--comparison=directed", "--roles=bc", "--model=M1c", "--dry_run=TRUE")$status != 0L)
  r <- run("02_fit_model.R", base, "--roles=bcj", "--model=M2c", "--ame_K=1", short); expect_true(r$status != 0L)
  r <- run("02_fit_model.R", base, "--model=M1b", short); expect_equal(r$status, 0L, info = paste(tail(r$out, 20), collapse = "\n"))
  expect_true(file.exists(file.path(out_dir, "fit_mid_lag3_all_pr_M1b_cdir_ent.rds"))); expect_equal(readRDS(file.path(out_dir, "fit_mid_lag3_all_pr_M1b_cdir_ent.rds"))$data_dir, "data_entry")
  for (cs in list(list("M1c", "bc", character(0)), list("M1c", "bc", "--joiner_set=opposing"), list("M1c", "bcj", character(0)), list("M1c", "bcj", "--joiner_set=opposing"), list("M1cr", "bcj", character(0)))) {
    r <- run("02_fit_model.R", base, paste0("--model=", cs[[1]]), paste0("--roles=", cs[[2]]), cs[[3]], short); expect_equal(r$status, 0L, info = paste(unlist(cs), paste(tail(r$out, 30), collapse = "\n")))
    cell <- sprintf("mid_lag3_all_pr_%s_cdir_%s%s_ent", cs[[1]], cs[[2]], if (length(cs[[3]])) "_jset" else "")
    o <- readRDS(file.path(out_dir, paste0("fit_", cell, ".rds")))
    expect_equal(o$dat$G, if (cs[[2]] == "bc") 2L else 3L); expect_equal(o$group_labels, stage2_group_labels(cs[[2]])); expect_equal(o$spec$entry, "later")
    if (length(cs[[3]])) { s2n <- tabulate(findInterval(which(o$dat$in_s2 == 1L), o$dat$cs_start), nbins = o$dat$N1); jn <- if (cs[[2]] == "bcj") which(o$dat$ev_grp == 2L) else NULL
      expect_equal(sum(o$dat$ev_alt_n < s2n[o$dat$ev_sy]), 2L); if (!is.null(jn)) expect_identical(which(o$dat$ev_alt_n < s2n[o$dat$ev_sy]), jn) }
  }
  r <- run("02_fit_model.R", base, "--model=M2c", "--roles=bcj", "--ame_K=1", "--joiner_set=opposing", short); expect_equal(r$status, 0L, info = paste(tail(r$out, 30), collapse = "\n"))
  o2 <- readRDS(file.path(out_dir, "fit_mid_lag3_all_pr_M2c_K1_cdir_bcj_jset_ent.rds")); expect_equal(o2$init$m1_cell, "mid_lag3_all_pr_M1c_cdir_bcj_jset_ent"); expect_equal(o2$dat$G, 3L)
  e2 <- run("03_extract.R", paste0("--fit=", file.path(out_dir, "fit_mid_lag3_all_pr_M2c_K1_cdir_bcj_jset_ent.rds"))); expect_equal(e2$status, 0L, info = paste(tail(e2$out, 30), collapse = "\n"))
  for (g in BCJ_GROUPS) expect_true(file.exists(file.path(out_dir, sprintf("typology_mid_lag3_all_pr_M2c_K1_cdir_bcj_jset_ent_%s.csv", g))), info = g)
  ## the bc M1c fit is never the init of a bcj M2c (another cell): explicitly named -> refused
  r <- run("02_fit_model.R", base, "--model=M2c", "--roles=bcj", "--ame_K=2", short, paste0("--init_from=", file.path(out_dir, "fit_mid_lag3_all_pr_M1c_cdir_bc_ent.rds"))); expect_true(r$status != 0L)
  dr <- run("05_robustness.R", "--preset=noww2", "--entry=later", "--dry_run=TRUE"); expect_true(dr$status != 0L)
  unlink(tmp, recursive = TRUE)
})

test_that("real entry build (when built): UK -> Germany 1939 and not the reverse, no late-target event, bc / bcj assemble, joiner sets restrict exactly the joiner events", {
  dir <- file.path(v2_root(), "data_entry")
  skip_if_not(file.exists(file.path(dir, "manifest.csv")) && file.exists(file.path(dir, "alts_lag3_all.rds")), "entry build (v2/data_entry) not present")
  ev <- readRDS(file.path(dir, "events.rds")) |> filter(!escalation)
  expect_true(any(ev$initiator == 200L & ev$target == 255L & ev$year == 1939L)); expect_false(any(ev$initiator == 255L & ev$target == 200L & ev$year == 1939L))
  expect_true(any(ev$initiator == 255L & ev$target == 290L & ev$year == 1939L))
  ro <- readRDS(file.path(dir, "roles_aprime.rds")); expect_false(any(ro$grp_aprime == 3L)); expect_false(any(ev$orig1 == 1L & ev$orig2 == 0L))
  for (s in c("fmid", "mid")) {
    for (r in c("bc", "bcj")) {
      A <- quiet(assemble_cell_spec(s, 3L, parse_spec(c("--comparison=directed", "--entry=later", paste0("--roles=", r), "--joiner_set=opposing")), model = "M1c"))
      expect_setequal(unique(A$ev_sets$ev_idx), which(A$ev_joiner)); expect_true(all(paste(which(A$ev_joiner), A$ev$alt_row[A$ev_joiner]) %in% paste(A$ev_sets$ev_idx, A$ev_sets$alt_row)))
      expect_equal(nrow(A$group_table), if (r == "bc") 2L else 3L); expect_true(all(A$group_table$n_usable > 0L))
      message(sprintf("entry build %s %s: %s; joiner events %d, size 1 %d, removed %d", s, r, paste(sprintf("%s %d (%.1f)", A$group_table$group, A$group_table$n_usable, A$group_table$w_usable), collapse = ", "),
                      A$jset_stats$n_joiner_events, A$jset_stats$n_singleton, A$jset_stats$n_removed))
      rm(A); invisible(gc())
    }
  }
})
