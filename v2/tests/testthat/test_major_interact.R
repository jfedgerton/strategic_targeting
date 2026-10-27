## The PI-approved stage-2 switches of 2026-10-27: major (none | joint) and interact (none | prefixed), their suffixes, the assembly on the
## toy cell (helper-toy-cell.R: states 10 and 40 are majors), the major_ix preset through the robustness driver and the submit script, and
## the real build when present. The defaults must leave the main specification byte-identical.
suppressPackageStartupMessages({ library(testthat); library(dplyr) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "assemble_cell.R"))

test_that("spec: major and interact default to none, add _major / _ix after _fullctrl, combine to _major_ix, and every suffix stays unique", {
  main <- parse_spec(character(0))
  expect_equal(main$major, "none"); expect_equal(main$interact, "none"); expect_true(spec_is_main(main)); expect_equal(spec_suffix(main), "")
  expect_equal(SPEC_CHOICES$major, c("none", "joint")); expect_equal(SPEC_CHOICES$interact, c("none", "prefixed"))
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec("--major=joint")), "fmid_lag3_all_pr_M1b_major")
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec("--interact=prefixed")), "fmid_lag3_all_pr_M1b_ix")
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec(c("--major=joint", "--interact=prefixed"))), "fmid_lag3_all_pr_M1b_major_ix")
  expect_equal(cell_name("mid", 3L, "M2b", parse_spec(c("--interact=prefixed", "--major=joint"))), "mid_lag3_all_pr_M2b_K2_major_ix")   # order fixed by SPEC_SUFFIX, not by the command line
  ## position in the suffix order: after controls, before weights / changes / direction
  all_on <- parse_spec(c("--pairs=originators", "--context=coalition", "--coal_ties=strict", "--coal_drop=top3", "--ww2=drop", "--controls=full",
                         "--major=joint", "--interact=prefixed", "--weights=none", "--changes=winsor", "--direction=sideA"))
  expect_equal(spec_suffix(all_on), "_orig_co_strict_drop3_noww2_fullctrl_major_ix_w1_winsor_sideA")
  expect_equal(names(SPEC_SUFFIX), c("pairs", "context", "coal_ties", "coal_drop", "ww2", "controls", "major", "interact", "weights", "changes", "direction"))
  ## uniqueness against every other suffix, the fixed ones and the archived ones; no suffix is a prefix of another followed by a letter run
  all_suffixes <- c(unlist(SPEC_SUFFIX), "_K2", "_pilot", "_STALE", "_samp10", "_dist", "_noCH", "_polrel")
  expect_false(any(duplicated(all_suffixes)))
  expect_false(spec_is_main(parse_spec("--major=joint"))); expect_false(spec_is_main(parse_spec("--interact=prefixed")))
  expect_error(parse_spec("--major=TRUE"), "--major must be one of"); expect_error(parse_spec("--interact=all"), "--interact must be one of")
  ## the six preset cells have six distinct names, none of them a main-grid or context name
  cells <- unlist(lapply(MAJOR_IX_PRESET$variants, function(v) vapply(MAJOR_IX_PRESET$cells, function(c) { p <- strsplit(c, ":")[[1]]
    s <- parse_spec(character(0)); for (nm in names(v)) s[[nm]] <- v[[nm]]; cell_name(p[1], as.integer(p[2]), p[3], s) }, "")))
  expect_equal(unname(cells), c("fmid_lag3_all_pr_M1b_major", "mid_lag3_all_pr_M1b_major", "fmid_lag3_all_pr_M1b_ix", "mid_lag3_all_pr_M1b_ix",
                                "fmid_lag3_all_pr_M1b_major_ix", "mid_lag3_all_pr_M1b_major_ix"))
  expect_equal(anyDuplicated(cells), 0L)
  expect_false(any(cells %in% c(vapply(MAIN_SAMPLES, function(s) cell_name(s, 3L, "M1b", main), ""), "mid_lag3_all_pr_M1b_bi", "mid_lag3_all_pr_M1b_co", "mid_lag3_all_pr_M1b_co_drop3")))
  ## the preset table
  expect_true("major_ix" %in% names(ROBUSTNESS_PRESETS)); expect_equal(ROBUSTNESS_PRESETS$major_ix, list())
  expect_equal(MAJOR_IX_PRESET$cells, c("fmid:3:M1b", "mid:3:M1b")); expect_equal(names(MAJOR_IX_PRESET$variants), c("major", "interact", "both"))
  for (v in MAJOR_IX_PRESET$variants) for (nm in names(v)) expect_true(v[[nm]] %in% SPEC_CHOICES[[nm]])
  expect_equal(MAJOR_IX_PRESET$variants$both, list(major = "joint", interact = "prefixed"))
  expect_false(any(grepl("war|M2", MAJOR_IX_PRESET$cells)))
})

test_that("column lists: W gains exactly cowmaj_tm1 only under major = joint; Z appends the indicators and the products in a fixed order", {
  expect_equal(w_controls(MAJOR = "joint"), c(w_controls(), "cowmaj_tm1"))
  expect_equal(w_controls(), c("M_lvl", "R_lvl", "ongoing_tm1", "log_prior_mids", "dem_tm1", "log_n_system_tm1"))
  expect_equal(w_controls(FULL_CONTROLS = TRUE, MAJOR = "joint"), w_controls(FULL_CONTROLS = TRUE))      # the full set already holds it, not twice
  expect_equal(sum(w_controls(FULL_CONTROLS = TRUE, MAJOR = "joint") == "cowmaj_tm1"), 1L)
  expect_equal(x_controls(), character(0))
  expect_equal(z_columns(TARGET_MISSING = TRUE), c("contig_land", "joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "dM_j_h", "dRpr_j_h", "dch_j_missing"))
  expect_equal(z_columns(TARGET_MISSING = TRUE, MAJOR = "joint"), c(z_columns(TARGET_MISSING = TRUE), "one_major", "both_major"))
  expect_equal(z_columns(TARGET_MISSING = TRUE, INTERACT = "prefixed"), c(z_columns(TARGET_MISSING = TRUE), "ix_dMj_dRj", "ix_dMi_dMj", "ix_dRi_dRj", "ix_dMi_dRj"))
  expect_equal(z_columns(MAJOR = "joint", INTERACT = "prefixed"), c(z_columns(), Z_MAJOR, Z_INTERACT))
  expect_equal(z_columns(FULL_CONTROLS = TRUE, MAJOR = "joint"), c(z_columns(FULL_CONTROLS = TRUE), Z_MAJOR))   # major_count (the count) stays, the indicators follow
  expect_true(all(Z_INTERACT %in% Z_CENTRED)); expect_false(any(Z_MAJOR %in% Z_CENTRED))                      # continuous products centred, 0/1 indicators not
  expect_error(w_controls(MAJOR = "both")); expect_error(z_columns(INTERACT = "all"))
})

test_that("assembly on the toy cell: indicators from the 0 / 1 / 2 major counts, uncentred; cowmaj_tm1 in W, centred; products equal the uncentred change products before their IV-row centring; defaults identical", {
  td <- write_toy_cell(tempfile("toycell"))
  A  <- quiet(assemble_cell("mid", 3L, "all", data_dir = td))
  M  <- quiet(assemble_cell("mid", 3L, "all", MAJOR = "joint", data_dir = td))
  I  <- quiet(assemble_cell("mid", 3L, "all", INTERACT = "prefixed", data_dir = td))
  MI <- quiet(assemble_cell("mid", 3L, "all", MAJOR = "joint", INTERACT = "prefixed", data_dir = td))
  d <- function(X) build_stan_data(X$sy, X$alts, X$ev, X$W_cols, X$X_cols, X$Z_cols, joint = 1L, x_ref = X$x_ref, allow_orphan_initiations = TRUE)
  ## (a) defaults: the plain call and the explicit none / none are identical, and assemble_cell_spec() of the main spec too
  A0 <- quiet(assemble_cell("mid", 3L, "all", MAJOR = "none", INTERACT = "none", data_dir = td)); S0 <- quiet(assemble_cell_spec("mid", 3L, parse_spec(character(0)), data_dir = td))
  expect_identical(d(A), d(A0)); expect_identical(d(A), d(S0)); expect_identical(A$W_cols, A0$W_cols); expect_identical(A$Z_cols, A0$Z_cols)
  expect_equal(A$major, "none"); expect_equal(A$interact, "none"); expect_null(A$major_stats)
  expect_false("cowmaj_tm1" %in% A$W_cols); expect_false(any(c(Z_MAJOR, Z_INTERACT) %in% A$Z_cols))
  ## (b) major = joint: W gains exactly cowmaj_tm1 (last), X unchanged, Z gains exactly the two indicators (last)
  expect_equal(M$W_cols, c(A$W_cols, "cowmaj_tm1")); expect_equal(M$X_cols, A$X_cols); expect_equal(M$Z_cols, c(A$Z_cols, "one_major", "both_major"))
  expect_equal(setdiff(M$W_cols, A$W_cols), "cowmaj_tm1"); expect_equal(setdiff(I$W_cols, A$W_cols), character(0))
  ## cowmaj_tm1 is grand-centred like dem_tm1 (the panel has majors 10 and 40 -> raw values 1, 0, 0, 1 in ccode order)
  panel <- readRDS(file.path(td, "stage1_panel_lag3.rds"))
  expect_equal(M$sy$cowmaj_tm1, panel$cowmaj_tm1[match(M$sy$ccode, panel$ccode)] - mean(panel$cowmaj_tm1)); expect_lt(abs(mean(M$sy$cowmaj_tm1)), 1e-12)
  expect_equal(d(M)$P1, d(A)$P1 + 1L); expect_equal(unname(d(M)$W[, match("cowmaj_tm1", M$W_cols)]), M$sy$cowmaj_tm1)
  ## the indicators: one_major = (major_count == 1), both_major = (major_count == 2), zero the reference, NOT centred
  alts0 <- readRDS(file.path(td, "alts_lag3_all.rds")); key <- function(a) paste(a$ccode1, a$ccode2)
  mc <- alts0$major_count[match(key(M$alts), key(alts0))]
  expect_equal(M$alts$one_major, as.integer(mc == 1L)); expect_equal(M$alts$both_major, as.integer(mc == 2L)); expect_true(all(M$alts$one_major %in% 0:1))
  expect_equal(sort(unique(mc)), 0:2)                                                                            # the toy has all three categories
  expect_equal(d(M)$Pz, d(A)$Pz + 2L); expect_equal(unname(d(M)$Z[, match("one_major", M$Z_cols)]), M$alts$one_major)
  ## shares and within-choice-set variation: 12 rows, 2 with zero majors, 8 with one, 2 with both; events 10->20 (1), 10->30 (1), 30->20 (0), 30->10 (1), 40->20 (1)
  expect_equal(unname(M$major_stats$share_alts), c(2, 8, 2) / 12); expect_equal(unname(M$major_stats$n_events), c(1L, 4L, 0L)); expect_equal(unname(M$major_stats$share_events), c(1, 4, 0) / 5)
  ## IV sets: 10 {20, 30} no variation; 20 {10, 30, 40} one_major varies; 30 {10, 20, 40} varies; 40 {10, 20, 30} one and both vary
  expect_equal(unname(M$major_stats$var_iv_sets), c(4L, 3L, 1L))
  ## stage-2 sets of the initiating state-years with a usable event (10, 30, 40): 10 {20, 30, 40} both vary; 30 {10, 20} one_major; 40 both
  expect_equal(unname(M$major_stats$var_s2_sets), c(3L, 3L, 2L))
  ## everything else is untouched by the indicators: same rows, events, weights, other Z columns
  expect_identical(M$ev, A$ev); expect_identical(M$sy$y1, A$sy$y1); expect_identical(M$alts[, A$Z_cols], A$alts[, A$Z_cols]); expect_identical(d(M)$Z[, match(A$Z_cols, M$Z_cols)], d(A)$Z)
  ## no within-set variation -> error (every dyad gets the same count)
  td_c <- write_toy_cell(tempfile("toycell"), alts_edit = function(a) { a$major_count <- 1L; a })
  expect_error(quiet(assemble_cell("mid", 3L, "all", MAJOR = "joint", data_dir = td_c)), "does not vary within any choice set")
  ## a major_count outside 0..2 is refused
  td_b <- write_toy_cell(tempfile("toycell"), alts_edit = function(a) { a$major_count[1] <- 3L; a })
  expect_error(quiet(assemble_cell("mid", 3L, "all", MAJOR = "joint", data_dir = td_b)))
  ## (c) interact = prefixed: exactly the four products, after the base Z; W unchanged
  expect_equal(I$Z_cols, c(A$Z_cols, Z_INTERACT)); expect_equal(I$W_cols, A$W_cols); expect_equal(d(I)$Pz, d(A)$Pz + 4L)
  ## before their centring the products are the uncentred dM_i / dR_i of the row's state-year (as in X) times the uncentred dM_j_h / dRpr_j_h
  ## of the target (the data columns as they enter Z before the IV-row centring)
  raw_Mj <- alts0$dM_j_h[match(key(I$alts), key(alts0))]; raw_Rj <- alts0$dRpr_j_h[match(key(I$alts), key(alts0))]
  dMi <- I$sy$dM_i[I$alts$sy_row]; dRi <- I$sy$dR_i[I$alts$sy_row]
  expect_equal(dMi, panel$dM_h[match(I$alts$ccode1, panel$ccode)])                                          # dM_i is the panel change, uncentred
  iv <- I$alts$in_iv == 1L
  prod <- list(ix_dMj_dRj = raw_Mj * raw_Rj, ix_dMi_dMj = dMi * raw_Mj, ix_dRi_dRj = dRi * raw_Rj, ix_dMi_dRj = dMi * raw_Rj)
  for (nm in Z_INTERACT) {
    expect_equal(I$alts[[nm]], prod[[nm]] - mean(prod[[nm]][iv]), info = nm)                              # product, then centred over the IV rows
    expect_lt(abs(mean(I$alts[[nm]][iv])), 1e-12)
    expect_equal(unname(d(I)$Z[, match(nm, I$Z_cols)]), I$alts[[nm]])
  }
  ## the products are NOT built from the centred Z columns (the toy IV mean of dM_j_h is not 0, so the two constructions differ)
  expect_false(isTRUE(all.equal(I$alts$ix_dMj_dRj, I$alts$dM_j_h * I$alts$dRpr_j_h - mean((I$alts$dM_j_h * I$alts$dRpr_j_h)[iv]))))
  ## a target without a change value contributes 0 to every product before the centring
  miss <- I$alts$dch_j_missing == 1L
  for (nm in Z_INTERACT) expect_true(all(abs(I$alts[[nm]][miss] + mean(prod[[nm]][iv])) < 1e-12), info = nm)
  expect_identical(I$ev, A$ev); expect_identical(I$alts[, A$Z_cols], A$alts[, A$Z_cols])
  ## with winsor the products use the CLIPPED changes on both sides
  Iw <- quiet(assemble_cell("mid", 3L, "all", INTERACT = "prefixed", CHANGES = "winsor", data_dir = td))
  capM <- Iw$winsor$M; capR <- Iw$winsor$R
  cMi <- Iw$sy$dM_i[Iw$alts$sy_row]; cMj <- pmin(pmax(raw_Mj, capM[1]), capM[2]); cRj <- pmin(pmax(raw_Rj, capR[1]), capR[2])
  expect_true(any(cMj != raw_Mj)); expect_true(any(cMi != dMi))
  p <- cMi * cMj; expect_equal(Iw$alts$ix_dMi_dMj, p - mean(p[iv]))
  p <- cMj * cRj; expect_equal(Iw$alts$ix_dMj_dRj, p - mean(p[iv]))
  ## (d) both switches: the union, in the documented order, and the base columns untouched
  expect_equal(MI$Z_cols, c(A$Z_cols, Z_MAJOR, Z_INTERACT)); expect_equal(MI$W_cols, M$W_cols)
  expect_identical(MI$alts[, Z_MAJOR], M$alts[, Z_MAJOR]); expect_identical(MI$alts[, Z_INTERACT], I$alts[, Z_INTERACT]); expect_identical(MI$ev, A$ev)
  expect_equal(d(MI)$Pz, d(A)$Pz + 6L); expect_equal(d(MI)$P1, d(A)$P1 + 1L)
  ## assemble_cell_spec() passes both switches through
  S <- quiet(assemble_cell_spec("mid", 3L, parse_spec(c("--major=joint", "--interact=prefixed")), data_dir = td)); expect_identical(d(S), d(MI))
  ## under full controls cowmaj_tm1 is not duplicated and the count stays beside the indicators
  F <- quiet(assemble_cell("mid", 3L, "all", MAJOR = "joint", FULL_CONTROLS = TRUE, data_dir = td))
  expect_equal(sum(F$W_cols == "cowmaj_tm1"), 1L); expect_true(all(c("major_count", Z_MAJOR) %in% F$Z_cols))
  expect_equal(F$W_cols, quiet(assemble_cell("mid", 3L, "all", FULL_CONTROLS = TRUE, data_dir = td))$W_cols)
})

test_that("robustness driver: the major_ix preset expands to the six (variant, cell) calls of 02_fit_model.R; the fixed switches are refused", {
  tmp <- tempfile("drv_v2_"); dir.create(file.path(tmp, "v2"), recursive = TRUE)
  file.copy(file.path(root, "R"), file.path(tmp, "v2"), recursive = TRUE); dir.create(file.path(tmp, "v2", "logs")); dir.create(file.path(tmp, "v2", "output"))
  rec <- file.path(tmp, "calls.txt")
  writeLines(c("args <- commandArgs(trailingOnly = TRUE)", sprintf("cat(paste(args, collapse = ' '), '\\n', sep = '', file = '%s', append = TRUE)", rec)), file.path(tmp, "v2", "R", "02_fit_model.R"))
  drv <- file.path(tmp, "v2", "R", "05_robustness.R")
  run <- function(...) { out <- suppressWarnings(system2("Rscript", c(drv, ...), stdout = TRUE, stderr = TRUE)); list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  r <- run("--preset=major_ix", "--dry_run=TRUE", "--threads_per_chain=2")
  expect_equal(r$status, 0L, info = paste(tail(r$out, 20), collapse = "\n"))
  calls <- readLines(rec); expect_length(calls, 6L)
  arg <- function(call, nm) unname(sub(paste0("^--", nm, "="), "", grep(paste0("^--", nm, "="), strsplit(call, " ")[[1]], value = TRUE)))
  args_of <- function(nm) unname(vapply(calls, arg, "", nm))
  ## variant slowest (major, interact, both), cell fastest (fmid, mid) = the --cell_index order 0-5 of MAJOR_IX_PRESET
  expect_equal(args_of("sample"), rep(c("fmid", "mid"), 3))
  expect_equal(args_of("major"), rep(c("joint", "none", "joint"), each = 2))
  expect_equal(args_of("interact"), rep(c("none", "prefixed", "prefixed"), each = 2))
  expect_true(all(args_of("lag") == "3")); expect_true(all(args_of("model") == "M1b")); expect_true(all(args_of("dry_run") == "TRUE")); expect_true(all(args_of("threads_per_chain") == "2"))
  for (nm in setdiff(names(SPEC_CHOICES), c("major", "interact"))) expect_true(all(args_of(nm) == SPEC_DEFAULTS[[nm]]), info = nm)
  expect_true(any(grepl("cell fmid_lag3_all_pr_M1b_major:", r$out, fixed = TRUE))); expect_true(any(grepl("cell mid_lag3_all_pr_M1b_ix:", r$out, fixed = TRUE)))
  expect_true(any(grepl("cell mid_lag3_all_pr_M1b_major_ix:", r$out, fixed = TRUE)))
  ## one task of the Slurm array: --cell_index=5 is mid, both; a pass-through switch adds its suffix
  file.remove(rec); r5 <- run("--preset=major_ix", "--cell_index=5", "--ww2=drop")
  expect_equal(r5$status, 0L); calls <- readLines(rec); expect_length(calls, 1L)
  expect_equal(arg(calls, "sample"), "mid"); expect_equal(arg(calls, "major"), "joint"); expect_equal(arg(calls, "interact"), "prefixed"); expect_equal(arg(calls, "ww2"), "drop")
  expect_true(any(grepl("cell mid_lag3_all_pr_M1b_noww2_major_ix:", r5$out, fixed = TRUE)))      # _noww2 precedes _major_ix in SPEC_SUFFIX order
  expect_true(run("--preset=major_ix", "--cell_index=6")$status != 0L)
  for (bad in c("--cells=mid:3:M1b", "--major=joint", "--interact=prefixed")) expect_true(run("--preset=major_ix", bad)$status != 0L, info = bad)
  ## the sbatch dispatches the preset
  sb <- readLines(file.path(root, "slurm", "04_robustness.sbatch")); expect_true(any(grepl("^  major_ix)", sb))); expect_true(any(grepl("--preset=major_ix --cell_index=\\$i", sb)))
  unlink(tmp, recursive = TRUE)
})

test_that("submit_robustness.sh: the major_ix preset goes as ONE array 0-5 with no dependency, like context; CELLS is refused", {
  tmp <- tempfile("sub_"); dir.create(file.path(tmp, "v2"), recursive = TRUE); dir.create(file.path(tmp, "bin"))
  file.copy(file.path(root, "slurm"), file.path(tmp, "v2"), recursive = TRUE)
  rec <- file.path(tmp, "sbatch_calls.txt")
  writeLines(c("#!/bin/bash", sprintf("echo \"PRESET='${PRESET-<unset>}' CELLS='${CELLS-<unset>}' args: $*\" >> %s", rec), "echo \"2001;cluster\""), file.path(tmp, "bin", "sbatch"))
  Sys.chmod(file.path(tmp, "bin", "sbatch"), "0755")
  run <- function(env, ...) { out <- suppressWarnings(system2("bash", c("v2/slurm/submit_robustness.sh", ...), stdout = TRUE, stderr = TRUE,
                                                              env = c(sprintf("PATH=%s:%s", file.path(tmp, "bin"), Sys.getenv("PATH")), env)))
    list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  old <- setwd(tmp); on.exit(setwd(old), add = TRUE)
  r <- run("PRESET=major_ix", "--account=open"); expect_equal(r$status, 0L, info = paste(r$out, collapse = "\n"))
  calls <- readLines(rec); expect_length(calls, 1L)
  expect_match(calls[1], "^PRESET='major_ix' CELLS='<unset>' args: --parsable --export=ALL --array=0-5 --account=open v2/slurm/04_robustness.sbatch$")
  expect_true(any(grepl("major_ix array submitted: job 2001 \\(indices 0-5: 0-1 major, 2-3 interact, 4-5 both; fmid, mid; M1b only, no dependency\\)", r$out)))
  expect_true(run(c("PRESET=major_ix", "CELLS=mid:3:M1b"))$status != 0L); expect_length(readLines(rec), 1L)
})

test_that("real build (when present): the six preset cells assemble at lag 3 with the expected dimensions, both indicators vary within choice sets, the defaults reproduce the saved dry runs", {
  dd <- file.path(v2_root(), "data")
  skip_if_not(file.exists(file.path(dd, "events.rds")) && file.exists(file.path(dd, "alts_lag3_all.rds")) && file.exists(file.path(dd, "manifest.csv")), "real data files not built")
  ## one assembly in memory at a time (1.7 million alternative rows per cell): the base cell's Stan data is kept, each variant is built,
  ## compared and dropped
  d <- function(X) build_stan_data(X$sy, X$alts, X$ev, X$W_cols, X$X_cols, X$Z_cols, joint = 1L, x_ref = X$x_ref, allow_orphan_initiations = TRUE, D_mode = 3L)
  for (s in c("fmid", "mid")) {
    A <- quiet(assemble_cell_spec(s, 3L, parse_spec(character(0)))); dA <- d(A); A_Z <- A$Z_cols; A_W <- A$W_cols; iv <- A$alts$in_iv == 1L
    expect_equal(c(dA$P1, dA$P, dA$Pz), c(8L, 2L, 7L))
    ## the saved dry run of the main cell, when made from these inputs, is reproduced by the defaults
    sd_path <- file.path(v2_root(), "output", sprintf("standata_%s_lag3_all_pr_M1b.rds", s))
    if (file.exists(sd_path)) {
      sd <- readRDS(sd_path); man <- read.csv(file.path(dd, "manifest.csv"), stringsAsFactors = FALSE)
      if (all(sd$manifest$sha256 == man$sha256[match(sd$manifest$file, man$file)])) {
        dat <- build_stan_data(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, x_ref = A$x_ref, allow_orphan_initiations = TRUE, D_mode = 3L, grainsize = sd$dat$grainsize)
        expect_identical(dat, sd$dat); rm(dat, sd)
      }
    }
    rm(A); invisible(gc())
    dims <- list(major = c(9L, 2L, 9L), interact = c(8L, 2L, 11L), both = c(9L, 2L, 13L))
    for (v in names(MAJOR_IX_PRESET$variants)) {
      sp <- parse_spec(character(0)); for (nm in names(MAJOR_IX_PRESET$variants[[v]])) sp[[nm]] <- MAJOR_IX_PRESET$variants[[v]][[nm]]
      X <- quiet(assemble_cell_spec(s, 3L, sp)); dX <- d(X)
      expect_equal(c(dX$P1, dX$P, dX$Pz), dims[[v]], info = paste(s, v))
      expect_equal(dX$E, dA$E); expect_equal(dX$N1, dA$N1); expect_equal(dX$N2, dA$N2); expect_identical(dX$y1, dA$y1); expect_identical(dX$ev_w, dA$ev_w); expect_identical(dX$ev_row, dA$ev_row)
      expect_identical(dX$Z[, match(A_Z, X$Z_cols)], dA$Z); expect_identical(dX$W[, match(A_W, X$W_cols)], dA$W); expect_identical(dX$X, dA$X)
      if (sp$major == "joint") {
        ## both indicators vary within some IV set and some event-bearing stage-2 set; every category holds events
        st <- X$major_stats
        expect_true(st$var_iv_sets[["one_major"]] > 0); expect_true(st$var_iv_sets[["both_major"]] > 0); expect_true(st$var_s2_sets[["one_major"]] > 0); expect_true(st$var_s2_sets[["both_major"]] > 0)
        expect_true(all(st$n_events > 0)); expect_equal(sum(st$n_events), dA$E); expect_equal(sum(st$share_alts), 1, tolerance = 1e-12)
        expect_true(all(dX$Z[, match("one_major", X$Z_cols)] + dX$Z[, match("both_major", X$Z_cols)] <= 1))
        expect_equal(sum(X$W_cols == "cowmaj_tm1"), 1L)
      }
      if (sp$interact == "prefixed") for (nm in Z_INTERACT) expect_lt(abs(mean(dX$Z[iv, match(nm, X$Z_cols)])), 1e-8)
      rm(X, dX); invisible(gc())
    }
    rm(dA); invisible(gc())
  }
})
