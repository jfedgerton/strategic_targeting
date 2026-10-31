## The `distance` switch of 2026-10-28 (PI request: capital distance as a stage-2 control, and distance + major powers; contiguity stays in):
## parsing and the _dist suffix, the Z columns with and without major = joint, the IV-row centring, the NA-distance row drop with its event
## reconciliation on the toy cell (helper-toy-cell.R), the dist preset through the robustness driver, the sbatch and the submit script, and the
## real build when present. The default (distance = none) must leave the main specification byte-identical.
suppressPackageStartupMessages({ library(testthat); library(dplyr) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "assemble_cell.R"))

test_that("spec: distance defaults to none, adds _dist after _fullctrl and before _major / _ix, is a switch (not archived), and the preset table is fixed", {
  main <- parse_spec(character(0))
  expect_equal(main$distance, "none"); expect_true(spec_is_main(main)); expect_equal(spec_suffix(main), "")
  expect_equal(SPEC_CHOICES$distance, c("none", "capital")); expect_equal(SPEC_SUFFIX$distance, c(capital = "_dist"))
  expect_equal(parse_spec("--distance=capital")$distance, "capital"); expect_false(spec_is_main(parse_spec("--distance=capital")))
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec("--distance=capital")), "fmid_lag3_all_pr_M1b_dist")
  expect_equal(cell_name("mid", 3L, "M1b", parse_spec(c("--distance=capital", "--major=joint"))), "mid_lag3_all_pr_M1b_dist_major")
  expect_equal(cell_name("mid", 3L, "M1b", parse_spec(c("--major=joint", "--distance=capital"))), "mid_lag3_all_pr_M1b_dist_major")   # order fixed by SPEC_SUFFIX, not by the command line
  expect_equal(cell_name("mid", 3L, "M2b", parse_spec(c("--distance=capital", "--interact=prefixed"))), "mid_lag3_all_pr_M2b_K2_dist_ix")
  expect_equal(cell_name("mid", 3L, "M1b", parse_spec(c("--controls=full", "--distance=capital", "--major=joint", "--interact=prefixed"))), "mid_lag3_all_pr_M1b_fullctrl_dist_major_ix")
  all_on <- parse_spec(c("--pairs=originators", "--context=coalition", "--coal_ties=strict", "--coal_drop=top3", "--ww2=drop", "--controls=full",
                         "--distance=capital", "--major=joint", "--interact=prefixed", "--weights=none", "--changes=winsor", "--direction=sideA"))
  expect_equal(spec_suffix(all_on), "_orig_co_strict_drop3_noww2_fullctrl_dist_major_ix_w1_winsor_sideA")
  expect_equal(names(SPEC_SUFFIX), c("pairs", "context", "coal_ties", "coal_drop", "ww2", "controls", "distance", "contiguity", "major", "interact", "comparison", "weights", "changes", "direction", "rho_prior", "sigma_prior"))   # rho_prior, sigma_prior appended 2026-10-30 (after every data switch)   # contiguity added the same day (test_contiguity.R)
  expect_false(any(duplicated(c(unlist(SPEC_SUFFIX), "_K2", "_pilot", "_STALE", "_samp10", "_noCH", "_polrel"))))
  ## the archived spelling --distance=TRUE is refused as a bad VALUE (parse_spec), not as an archived option; the other archived flags still are
  expect_error(parse_spec("--distance=TRUE"), "--distance must be one of none | capital", fixed = TRUE)
  expect_error(parse_spec("--distance=FALSE"), "--distance must be one of"); expect_error(parse_spec("--distance=1"), "--distance must be one of")
  ## a bare --<switch> (no =value) is refused for every switch, so `02_fit_model.R --distance` can never run the main cell silently (audit of ec9d8ae)
  for (nm in names(SPEC_CHOICES)) expect_error(parse_spec(paste0("--", nm)), sprintf("--%s needs a value", nm), fixed = TRUE)
  expect_error(parse_spec("--ame_K"), "--ame_K needs a value", fixed = TRUE)
  expect_error(parse_spec(c("--major=joint", "--distance")), "--distance needs a value", fixed = TRUE); expect_error(parse_spec("--major"), "--major needs a value", fixed = TRUE)
  expect_equal(parse_spec("--distance=capital")$distance, "capital")                                   # the valued form is unaffected
  expect_false("distance" %in% ARCHIVED_OPTIONS)
  expect_true(refuse_archived_options("--distance=capital", "02_fit_model.R"))
  for (old in c("variant", "sampled_m", "centrality", "conflict_history", "dy_peace_years", "drop_ww2", "full_controls", "event_weights", "K")) {   # contiguity left the list the same day
    expect_true(old %in% ARCHIVED_OPTIONS); expect_error(refuse_archived_options(paste0("--", old, "=1"), "02_fit_model.R"), "is not an option")
  }
  ## the build still takes --direction only, so --distance is refused there through the switch list
  expect_error(refuse_archived_options("--distance=capital", "01_build_data.R", also = setdiff(names(SPEC_CHOICES), "direction")), "--distance is not an option of 01_build_data.R")
  ## the old global constant is gone: nothing but the switch can set the distance
  spec_src <- readLines(file.path(root, "R", "spec.R")); expect_false(any(grepl("^DISTANCE\\s*<-", spec_src)))
  for (f in c("02_fit_model.R", "05_robustness.R", "assemble_cell.R")) expect_false(any(grepl("DISTANCE = DISTANCE|distance = DISTANCE, contiguity|!DISTANCE|if \\(DISTANCE\\)", readLines(file.path(root, "R", f)))), info = f)
  ## the preset table: two variants x two cells, M1b, lag 3, no interactions, no war, no M2
  expect_true("dist" %in% names(ROBUSTNESS_PRESETS)); expect_equal(ROBUSTNESS_PRESETS$dist, list())
  expect_equal(DIST_PRESET$cells, c("fmid:3:M1b", "mid:3:M1b")); expect_equal(names(DIST_PRESET$variants), c("distance", "distance_major"))
  expect_equal(DIST_PRESET$variants$distance, list(distance = "capital")); expect_equal(DIST_PRESET$variants$distance_major, list(distance = "capital", major = "joint"))
  for (v in DIST_PRESET$variants) { for (nm in names(v)) expect_true(v[[nm]] %in% SPEC_CHOICES[[nm]]); expect_false("interact" %in% names(v)) }
  cells <- unlist(lapply(DIST_PRESET$variants, function(v) vapply(DIST_PRESET$cells, function(c) { p <- strsplit(c, ":")[[1]]
    s <- parse_spec(character(0)); for (nm in names(v)) s[[nm]] <- v[[nm]]; cell_name(p[1], as.integer(p[2]), p[3], s) }, "")))
  expect_equal(unname(cells), c("fmid_lag3_all_pr_M1b_dist", "mid_lag3_all_pr_M1b_dist", "fmid_lag3_all_pr_M1b_dist_major", "mid_lag3_all_pr_M1b_dist_major"))
  expect_equal(anyDuplicated(cells), 0L)
  expect_false(any(cells %in% c(vapply(MAIN_SAMPLES, function(s) cell_name(s, 3L, "M1b", main), ""), "fmid_lag3_all_pr_M1b_major", "mid_lag3_all_pr_M1b_major", "mid_lag3_all_pr_M1b_major_ix")))
})

test_that("column lists: log_capdist after the dyadic controls and before the target changes, centred; the major indicators follow it; the logical flag is refused", {
  expect_equal(z_columns(), c("contig_land", "joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "dM_j_h", "dRpr_j_h"))
  expect_equal(z_columns(DISTANCE = "none"), z_columns())
  expect_equal(z_columns(DISTANCE = "capital"), c("contig_land", "joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "log_capdist", "dM_j_h", "dRpr_j_h"))
  expect_equal(z_columns(DISTANCE = "capital", MAJOR = "joint"), c("contig_land", "joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "log_capdist", "dM_j_h", "dRpr_j_h", "one_major", "both_major"))
  expect_equal(z_columns(DISTANCE = "capital", TARGET_MISSING = TRUE, MAJOR = "joint"), c("contig_land", "joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "log_capdist", "dM_j_h", "dRpr_j_h", "dch_j_missing", "one_major", "both_major"))
  expect_equal(z_columns(DISTANCE = "capital", FULL_CONTROLS = TRUE), c("contig_land", "joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "alliance_def_off", "log_dy_peace_years", "log_capdist", "dM_j_h", "dRpr_j_h"))
  expect_true("contig_land" %in% z_columns(DISTANCE = "capital"))                 # contiguity stays in
  expect_true("log_capdist" %in% Z_CENTRED)                                      # continuous, centred over the IV rows
  for (bad in list(TRUE, FALSE, "yes", c("none", "capital"), 1L)) expect_error(z_columns(DISTANCE = bad), info = deparse(bad))
  ## W is untouched by the switch (the distance is dyadic)
  expect_equal(w_controls(), c("M_lvl", "R_lvl", "ongoing_tm1", "log_prior_mids", "dem_tm1", "log_n_system_tm1"))
})

test_that("assembly on the toy cell: log_capdist centred over the IV rows, everything else identical; with major the indicators follow; defaults identical", {
  td <- write_toy_cell(tempfile("toycell"))
  A  <- quiet(assemble_cell("mid", 3L, "all", data_dir = td))
  D  <- quiet(assemble_cell("mid", 3L, "all", DISTANCE = "capital", data_dir = td))
  DM <- quiet(assemble_cell("mid", 3L, "all", DISTANCE = "capital", MAJOR = "joint", data_dir = td))
  d <- function(X) build_stan_data(X$sy, X$alts, X$ev, X$W_cols, X$X_cols, X$Z_cols, joint = 1L, x_ref = X$x_ref, allow_orphan_initiations = TRUE)
  ## (a) defaults: the plain call, the explicit none and assemble_cell_spec() of the main spec are identical; no distance statistics
  A0 <- quiet(assemble_cell("mid", 3L, "all", DISTANCE = "none", data_dir = td)); S0 <- quiet(assemble_cell_spec("mid", 3L, parse_spec(character(0)), data_dir = td))
  expect_identical(d(A), d(A0)); expect_identical(d(A), d(S0)); expect_identical(A$Z_cols, A0$Z_cols)
  expect_equal(A$distance, "none"); expect_null(A$distance_stats); expect_false("log_capdist" %in% A$Z_cols)
  expect_error(quiet(assemble_cell("mid", 3L, "all", DISTANCE = TRUE, data_dir = td)))     # the old logical flag is refused
  ## (b) distance = capital on a toy without NA distances: same rows, same events, same W / X; Z gains exactly log_capdist before the target changes
  expect_equal(D$Z_cols, c("contig_land", "joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "log_capdist", "dM_j_h", "dRpr_j_h", "dch_j_missing"))
  expect_equal(D$W_cols, A$W_cols); expect_equal(D$X_cols, A$X_cols); expect_identical(D$ev, A$ev); expect_identical(D$sy, A$sy)
  expect_equal(nrow(D$alts), nrow(A$alts)); expect_identical(D$alts[, A$Z_cols], A$alts[, A$Z_cols])
  expect_equal(d(D)$Pz, d(A)$Pz + 1L); expect_equal(d(D)$E, d(A)$E); expect_equal(d(D)$N2, d(A)$N2); expect_identical(d(D)$Z[, match(A$Z_cols, D$Z_cols)], d(A)$Z)
  ## the centring: raw log1p(km) of the file minus its mean over the IV rows; the stage-2-only row (10 -> 40) gets the same constant
  alts0 <- readRDS(file.path(td, "alts_lag3_all.rds")); key <- function(a) paste(a$ccode1, a$ccode2)
  raw <- alts0$log_capdist[match(key(D$alts), key(alts0))]; iv <- D$alts$in_iv == 1L
  expect_equal(D$alts$log_capdist, raw - mean(raw[iv])); expect_lt(abs(mean(D$alts$log_capdist[iv])), 1e-12); expect_false(abs(mean(D$alts$log_capdist)) < 1e-12)   # not centred over all rows
  expect_equal(unname(d(D)$Z[, match("log_capdist", D$Z_cols)]), D$alts$log_capdist)
  ## the statistics: nothing dropped, cor with contiguity over the IV rows, chosen vs stage-2 means on the UNCENTRED scale
  st <- D$distance_stats
  expect_equal(st$n_rows_dropped, 0L); expect_equal(st$n_rows_dropped_iv, 0L); expect_equal(st$n_rows_dropped_s2, 0L); expect_equal(st$n_sy_affected, 0L); expect_equal(st$n_events_dropped, 0L)
  expect_equal(st$mean_iv, mean(raw[iv])); expect_equal(st$sd_iv, sd(raw[iv])); expect_equal(st$cor_contig_iv, cor(raw[iv], D$alts$contig_land[iv]))
  expect_equal(st$mean_chosen, mean(raw[D$ev$alt_row])); expect_equal(st$mean_s2, mean(raw[D$alts$in_s2 == 1L])); expect_equal(length(D$ev$alt_row), 5L)
  ## (c) distance + major: log_capdist, then the two indicators; the indicators and W as in the major-only cell
  M <- quiet(assemble_cell("mid", 3L, "all", MAJOR = "joint", data_dir = td))
  expect_equal(DM$Z_cols, c(D$Z_cols, "one_major", "both_major")); expect_equal(DM$W_cols, M$W_cols); expect_equal(DM$W_cols, c(A$W_cols, "cowmaj_tm1"))
  expect_identical(DM$alts[, Z_MAJOR], M$alts[, Z_MAJOR]); expect_identical(DM$alts$log_capdist, D$alts$log_capdist); expect_identical(DM$ev, A$ev)
  expect_equal(d(DM)$Pz, d(A)$Pz + 3L); expect_equal(d(DM)$P1, d(A)$P1 + 1L)
  expect_identical(DM$distance_stats, D$distance_stats); expect_identical(DM$major_stats, M$major_stats)
  ## assemble_cell_spec() passes the switch through, alone and with major
  S <- quiet(assemble_cell_spec("mid", 3L, parse_spec("--distance=capital"), data_dir = td)); expect_identical(d(S), d(D))
  SM <- quiet(assemble_cell_spec("mid", 3L, parse_spec(c("--distance=capital", "--major=joint")), data_dir = td)); expect_identical(d(SM), d(DM))
  ## M1c (CONTEXT_MODEL) with the distance is not specified and refused
  expect_error(quiet(assemble_cell("mid", 3L, "all", DISTANCE = "capital", CONTEXT_MODEL = TRUE, data_dir = td)), "not specified")
})

test_that("NA distances: rows leave both sets at once, state-years and events are counted, a chosen target without a distance drops its event and is reconciled, an emptied IV set is an error", {
  d <- function(X) build_stan_data(X$sy, X$alts, X$ev, X$W_cols, X$X_cols, X$Z_cols, joint = 1L, x_ref = X$x_ref, allow_orphan_initiations = TRUE)
  na_at <- function(pairs) function(a) { for (p in pairs) a$log_capdist[a$ccode1 == p[1] & a$ccode2 == p[2]] <- NA_real_; a }
  td0 <- write_toy_cell(tempfile("toycell")); A <- quiet(assemble_cell("mid", 3L, "all", data_dir = td0))
  ## the default never touches a NA distance: the main specification on a file WITH NA distances is identical to the clean file's
  td_na <- write_toy_cell(tempfile("toycell"), alts_edit = na_at(list(c(20L, 30L), c(40L, 20L))))
  A_na <- quiet(assemble_cell("mid", 3L, "all", data_dir = td_na)); expect_identical(d(A_na), d(A))
  ## (a) one non-chosen alternative (20 -> 30, in both sets; 20 never initiates): one row dropped from both sets, one state-year affected, no event lost
  td1 <- write_toy_cell(tempfile("toycell"), alts_edit = na_at(list(c(20L, 30L))))
  D1 <- quiet(assemble_cell("mid", 3L, "all", DISTANCE = "capital", data_dir = td1)); st <- D1$distance_stats
  expect_equal(nrow(D1$alts), nrow(A$alts) - 1L); expect_false(any(D1$alts$ccode1 == 20L & D1$alts$ccode2 == 30L)); expect_false(anyNA(D1$alts$log_capdist))
  expect_equal(st$n_rows_dropped, 1L); expect_equal(st$n_rows_dropped_iv, 1L); expect_equal(st$n_rows_dropped_s2, 1L); expect_equal(st$n_sy_affected, 1L); expect_equal(st$n_events_dropped, 0L)
  expect_equal(d(D1)$E, d(A)$E); expect_equal(d(D1)$N2, d(A)$N2 - 1L); expect_equal(d(D1)$N1, d(A)$N1); expect_identical(d(D1)$y1, d(A)$y1)
  expect_equal(sum(D1$alts$in_iv[D1$alts$ccode1 == 20L]), 2L)                                       # 20's IV set shrank from 3 to 2
  ## the centring constant is the mean over the REMAINING IV rows
  alts1 <- readRDS(file.path(td1, "alts_lag3_all.rds")); key <- function(a) paste(a$ccode1, a$ccode2)
  raw1 <- alts1$log_capdist[match(key(D1$alts), key(alts1))]; expect_false(anyNA(raw1)); expect_equal(D1$alts$log_capdist, raw1 - mean(raw1[D1$alts$in_iv == 1L]))
  ## the events are the same five, their alt_row re-indexed into the smaller block
  expect_equal(nrow(D1$ev), 5L); expect_true(all(D1$alts$in_s2[D1$ev$alt_row] == 1L)); expect_true(all(D1$alts$sy_row[D1$ev$alt_row] == D1$ev$sy_row))
  ## (b) a chosen target without a distance: 40 -> 20 is 40's only event (dispute 7), so the event leaves stage 2, is counted, and 40 keeps y1 = 1 as an orphan
  td2 <- write_toy_cell(tempfile("toycell"), alts_edit = na_at(list(c(40L, 20L))))
  D2 <- quiet(assemble_cell("mid", 3L, "all", DISTANCE = "capital", data_dir = td2)); st <- D2$distance_stats
  expect_equal(st$n_rows_dropped, 1L); expect_equal(st$n_events_dropped, 1L); expect_equal(st$n_sy_affected, 1L)
  expect_equal(nrow(D2$ev), 4L); expect_equal(d(D2)$E, d(A)$E - 1L); expect_equal(d(D2)$N1, d(A)$N1); expect_identical(d(D2)$y1, d(A)$y1)   # stage 1 unchanged
  expect_equal(D2$n_orphan, A$n_orphan + 1L)
  sy40 <- which(D2$sy$ccode == 40L); expect_equal(D2$sy$y1[sy40], 1L); expect_false(sy40 %in% D2$ev$sy_row)
  expect_false(any(D2$alts$ccode1 == 40L & D2$alts$ccode2 == 20L))
  expect_equal(st$mean_chosen, mean(readRDS(file.path(td2, "alts_lag3_all.rds")) |> semi_join(D2$alts[D2$ev$alt_row, c("ccode1", "ccode2")], by = c("ccode1", "ccode2")) |> pull(log_capdist)))
  ## with the chosen target AND a non-chosen one missing the counts add up row by row
  td3 <- write_toy_cell(tempfile("toycell"), alts_edit = na_at(list(c(40L, 20L), c(20L, 30L), c(10L, 40L))))      # 10 -> 40 is stage-2 only (in_iv = 0)
  D3 <- quiet(assemble_cell("mid", 3L, "all", DISTANCE = "capital", data_dir = td3)); st <- D3$distance_stats
  expect_equal(st$n_rows_dropped, 3L); expect_equal(st$n_rows_dropped_iv, 2L); expect_equal(st$n_rows_dropped_s2, 3L); expect_equal(st$n_sy_affected, 3L); expect_equal(st$n_events_dropped, 1L)
  expect_equal(nrow(D3$ev), 4L); expect_equal(d(D3)$N2, d(A)$N2 - 3L)
  ## (c) every IV-set alternative of a state-year without a distance: an error, not a silent drop (20's IV set is 10, 30, 40)
  td4 <- write_toy_cell(tempfile("toycell"), alts_edit = na_at(list(c(20L, 10L), c(20L, 30L), c(20L, 40L))))
  expect_error(quiet(assemble_cell("mid", 3L, "all", DISTANCE = "capital", data_dir = td4)), "would lose every IV-set alternative")
  ## (d) a stage-2 set reduced to the chosen target alone: 30's stage-2 set is {10, 20} with events to both; dropping 30 -> 10 leaves 30 -> 20 as the
  ## only alternative, which the existing chosen-only check refuses
  td5 <- write_toy_cell(tempfile("toycell"), alts_edit = na_at(list(c(30L, 10L))))
  expect_error(quiet(assemble_cell("mid", 3L, "all", DISTANCE = "capital", data_dir = td5)), "holds only the chosen target")
  ## the same NA file under distance = none is untouched by all of this
  expect_identical(d(quiet(assemble_cell("mid", 3L, "all", data_dir = td4))), d(A)); expect_identical(d(quiet(assemble_cell("mid", 3L, "all", data_dir = td5))), d(A))
})

test_that("robustness driver: the dist preset expands to the four (variant, cell) calls of 02_fit_model.R; the fixed switches are refused; perm_error refuses --distance", {
  tmp <- tempfile("drv_v2_"); dir.create(file.path(tmp, "v2"), recursive = TRUE)
  file.copy(file.path(root, "R"), file.path(tmp, "v2"), recursive = TRUE); dir.create(file.path(tmp, "v2", "logs")); dir.create(file.path(tmp, "v2", "output"))
  rec <- file.path(tmp, "calls.txt")
  writeLines(c("args <- commandArgs(trailingOnly = TRUE)", sprintf("cat(paste(args, collapse = ' '), '\\n', sep = '', file = '%s', append = TRUE)", rec)), file.path(tmp, "v2", "R", "02_fit_model.R"))
  drv <- file.path(tmp, "v2", "R", "05_robustness.R")
  run <- function(...) { out <- suppressWarnings(system2("Rscript", c(drv, ...), stdout = TRUE, stderr = TRUE)); list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  r <- run("--preset=dist", "--dry_run=TRUE", "--threads_per_chain=2")
  expect_equal(r$status, 0L, info = paste(tail(r$out, 20), collapse = "\n"))
  calls <- readLines(rec); expect_length(calls, 4L)
  arg <- function(call, nm) unname(sub(paste0("^--", nm, "="), "", grep(paste0("^--", nm, "="), strsplit(call, " ")[[1]], value = TRUE)))
  args_of <- function(nm) unname(vapply(calls, arg, "", nm))
  ## variant slowest (distance, distance_major), cell fastest (fmid, mid) = the --cell_index order 0-3 of DIST_PRESET
  expect_equal(args_of("sample"), rep(c("fmid", "mid"), 2))
  expect_equal(args_of("distance"), rep("capital", 4))
  expect_equal(args_of("major"), rep(c("none", "joint"), each = 2))
  expect_true(all(args_of("interact") == "none"))
  expect_true(all(args_of("lag") == "3")); expect_true(all(args_of("model") == "M1b")); expect_true(all(args_of("dry_run") == "TRUE")); expect_true(all(args_of("threads_per_chain") == "2"))
  for (nm in setdiff(names(SPEC_CHOICES), c("distance", "major"))) expect_true(all(args_of(nm) == SPEC_DEFAULTS[[nm]]), info = nm)
  expect_true(any(grepl("cell fmid_lag3_all_pr_M1b_dist:", r$out, fixed = TRUE))); expect_true(any(grepl("cell mid_lag3_all_pr_M1b_dist:", r$out, fixed = TRUE)))
  expect_true(any(grepl("cell fmid_lag3_all_pr_M1b_dist_major:", r$out, fixed = TRUE))); expect_true(any(grepl("cell mid_lag3_all_pr_M1b_dist_major:", r$out, fixed = TRUE)))
  ## one task of the Slurm array: --cell_index=3 is mid, distance_major; a pass-through switch adds its suffix in SPEC_SUFFIX order
  file.remove(rec); r3 <- run("--preset=dist", "--cell_index=3", "--ww2=drop")
  expect_equal(r3$status, 0L); calls <- readLines(rec); expect_length(calls, 1L)
  expect_equal(arg(calls, "sample"), "mid"); expect_equal(arg(calls, "distance"), "capital"); expect_equal(arg(calls, "major"), "joint"); expect_equal(arg(calls, "ww2"), "drop")
  expect_true(any(grepl("cell mid_lag3_all_pr_M1b_noww2_dist_major:", r3$out, fixed = TRUE)))
  expect_true(run("--preset=dist", "--cell_index=4")$status != 0L)
  for (bad in c("--cells=mid:3:M1b", "--distance=capital", "--distance=none", "--major=joint")) {
    rb <- run("--preset=dist", bad); expect_true(rb$status != 0L, info = bad); expect_true(any(grepl("the dist preset fixes --cells, --distance and --major", rb$out, fixed = TRUE)), info = bad)
  }
  ## the other fixed-table presets refuse --distance only where they fix every switch (context_b); major_ix and context pass it through
  expect_true(run("--preset=context_b", "--distance=capital")$status != 0L)
  ## --distance=capital through a refit preset gives the combined suffix (the switch is a pass-through there)
  file.remove(rec); rn <- run("--preset=noww2", "--cells=fmid:3:M1b", "--distance=capital", "--dry_run=TRUE")
  expect_equal(rn$status, 0L); expect_true(any(grepl("cell fmid_lag3_all_pr_M1b_noww2_dist:", rn$out, fixed = TRUE)))
  ## perm_error refuses --distance (the refit-based tests have not been wired for the NA-row drop), before any fit or data file is touched
  rb <- run("--preset=perm_error", "--family=m1", "--distance=capital")
  expect_true(rb$status != 0L); expect_true(any(grepl("perm_error refuses --distance=capital", rb$out, fixed = TRUE)), info = paste(tail(rb$out, 5), collapse = "\n"))
  ## the archived spelling is refused by the driver as a bad value; a bare --distance (no value) is refused too, by parse_spec(), before any call
  ra <- run("--preset=dist", "--distance=TRUE"); expect_true(ra$status != 0L)
  file.remove(rec); rbare <- run("--preset=noww2", "--distance"); expect_true(rbare$status != 0L); expect_true(any(grepl("--distance needs a value", rbare$out, fixed = TRUE))); expect_false(file.exists(rec))
  ## the sbatch dispatches the preset with the 0-3 range guard
  sb <- readLines(file.path(root, "slurm", "04_robustness.sbatch"))
  expect_true(any(grepl("^  dist)", sb))); expect_true(any(grepl("--preset=dist --cell_index=\\$i", sb))); expect_true(any(grepl('\\[ "\\$i" -le 3 \\] \\|\\| \\{ echo "dist: array index', sb)))
  expect_true(any(grepl("major_ix\\|dist\\|dist_nocontig\\|perm_error", sb)))
  unlink(tmp, recursive = TRUE)
})

test_that("submit_robustness.sh: the dist preset goes as ONE array 0-3 with no dependency; CELLS is refused", {
  tmp <- tempfile("sub_"); dir.create(file.path(tmp, "v2"), recursive = TRUE); dir.create(file.path(tmp, "bin"))
  file.copy(file.path(root, "slurm"), file.path(tmp, "v2"), recursive = TRUE)
  rec <- file.path(tmp, "sbatch_calls.txt")
  writeLines(c("#!/bin/bash", sprintf("echo \"PRESET='${PRESET-<unset>}' CELLS='${CELLS-<unset>}' args: $*\" >> %s", rec), "echo \"3001;cluster\""), file.path(tmp, "bin", "sbatch"))
  Sys.chmod(file.path(tmp, "bin", "sbatch"), "0755")
  run <- function(env, ...) { out <- suppressWarnings(system2("bash", c("v2/slurm/submit_robustness.sh", ...), stdout = TRUE, stderr = TRUE,
                                                              env = c(sprintf("PATH=%s:%s", file.path(tmp, "bin"), Sys.getenv("PATH")), env)))
    list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  old <- setwd(tmp); on.exit(setwd(old), add = TRUE)
  r <- run("PRESET=dist", "--account=open"); expect_equal(r$status, 0L, info = paste(r$out, collapse = "\n"))
  calls <- readLines(rec); expect_length(calls, 1L)
  expect_match(calls[1], "^PRESET='dist' CELLS='<unset>' args: --parsable --export=ALL --array=0-3 --account=open v2/slurm/04_robustness.sbatch$")
  expect_true(any(grepl("dist array submitted: job 3001 \\(indices 0-3: 0-1 distance, 2-3 distance_major; fmid, mid; M1b only, no dependency\\)", r$out)))
  expect_true(run(c("PRESET=dist", "CELLS=mid:3:M1b"))$status != 0L); expect_length(readLines(rec), 1L)
  expect_true(run("PRESET=dist", "--array=0-1")$status != 0L); expect_length(readLines(rec), 1L)
})

test_that("real build (when present): the four preset cells assemble at lag 3 with the expected dimensions, no row is dropped (no NA distance in the build), the defaults reproduce the saved dry runs byte for byte", {
  dd <- file.path(v2_root(), "data")
  skip_if_not(file.exists(file.path(dd, "events.rds")) && file.exists(file.path(dd, "alts_lag3_all.rds")) && file.exists(file.path(dd, "manifest.csv")), "real data files not built")
  expect_equal(sum(is.na(readRDS(file.path(dd, "alts_lag3_all.rds"))$log_capdist)), 0L)
  d <- function(X) build_stan_data(X$sy, X$alts, X$ev, X$W_cols, X$X_cols, X$Z_cols, joint = 1L, x_ref = X$x_ref, allow_orphan_initiations = TRUE, D_mode = 3L)
  for (s in c("fmid", "mid")) {
    A <- quiet(assemble_cell_spec(s, 3L, parse_spec(character(0)))); dA <- d(A); A_Z <- A$Z_cols; A_W <- A$W_cols; iv <- A$alts$in_iv == 1L
    expect_equal(c(dA$P1, dA$P, dA$Pz), c(8L, 2L, 7L)); expect_null(A$distance_stats)
    ## main-spec invariance: the saved dry run of the main cell, when made from these inputs, is reproduced identically by the defaults
    sd_path <- file.path(v2_root(), "output", sprintf("standata_%s_lag3_all_pr_M1b.rds", s))
    if (file.exists(sd_path)) {
      sd <- readRDS(sd_path); man <- read.csv(file.path(dd, "manifest.csv"), stringsAsFactors = FALSE)
      if (all(sd$manifest$sha256 == man$sha256[match(sd$manifest$file, man$file)])) {
        dat <- build_stan_data(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, x_ref = A$x_ref, allow_orphan_initiations = TRUE, D_mode = 3L, grainsize = sd$dat$grainsize)
        expect_identical(dat, sd$dat); rm(dat, sd)
      }
    }
    rm(A); invisible(gc())
    dims <- list(distance = c(8L, 2L, 8L), distance_major = c(9L, 2L, 10L))
    for (v in names(DIST_PRESET$variants)) {
      sp <- parse_spec(character(0)); for (nm in names(DIST_PRESET$variants[[v]])) sp[[nm]] <- DIST_PRESET$variants[[v]][[nm]]
      X <- quiet(assemble_cell_spec(s, 3L, sp)); dX <- d(X)
      expect_equal(c(dX$P1, dX$P, dX$Pz), dims[[v]], info = paste(s, v))
      expect_equal(dX$E, dA$E); expect_equal(dX$N1, dA$N1); expect_equal(dX$N2, dA$N2); expect_identical(dX$y1, dA$y1); expect_identical(dX$ev_w, dA$ev_w); expect_identical(dX$ev_row, dA$ev_row)
      expect_identical(dX$Z[, match(A_Z, X$Z_cols)], dA$Z); expect_identical(dX$W[, match(A_W, X$W_cols)], dA$W); expect_identical(dX$X, dA$X)
      st <- X$distance_stats
      expect_equal(st$n_rows_dropped, 0L); expect_equal(st$n_events_dropped, 0L); expect_equal(st$n_sy_affected, 0L)
      expect_lt(abs(mean(dX$Z[iv, match("log_capdist", X$Z_cols)])), 1e-8)
      expect_true(is.finite(st$cor_contig_iv) && st$cor_contig_iv < 0)            # nearer capitals are the contiguous ones
      expect_true(is.finite(st$mean_chosen) && is.finite(st$mean_s2) && st$mean_chosen < st$mean_s2)   # chosen targets are nearer than the average alternative
      if (sp$major == "joint") { expect_equal(sum(X$W_cols == "cowmaj_tm1"), 1L); expect_equal(tail(X$Z_cols, 2), Z_MAJOR); expect_true(all(X$major_stats$n_events > 0)) }
      rm(X, dX); invisible(gc())
    }
    rm(dA); invisible(gc())
  }
})
