## The `contiguity` switch of 2026-10-28 (PI request, same day as `distance`: the four dist cells WITHOUT contiguity, so the two can be compared):
## parsing and the _nocontig suffix (adjacent to _dist), the refusal without distance = capital, the archived spellings, the Z columns (no
## contig_land; log_capdist centred over the IV rows), the toy-cell assembly against the distance cell, the dist_nocontig preset through the
## robustness driver, the sbatch and the submit script, and the real build when present. The default (contiguity = land) must leave the main
## specification AND every existing preset byte-identical. Follows test_distance.R.
suppressPackageStartupMessages({ library(testthat); library(dplyr) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "assemble_cell.R"))

test_that("spec: contiguity defaults to land, none needs distance = capital, _nocontig follows _dist, old spellings are refused, archived list without it", {
  main <- parse_spec(character(0))
  expect_equal(main$contiguity, "land"); expect_true(spec_is_main(main)); expect_equal(spec_suffix(main), "")
  expect_equal(SPEC_CHOICES$contiguity, c("land", "none")); expect_equal(SPEC_SUFFIX$contiguity, c(none = "_nocontig"))
  expect_equal(SPEC_DEFAULTS$contiguity, "land")
  ## none is only meaningful with the distance: refused loudly whatever else is set
  expect_error(parse_spec("--contiguity=none"), "--contiguity=none needs --distance=capital", fixed = TRUE)
  expect_error(parse_spec(c("--contiguity=none", "--distance=none")), "needs --distance=capital", fixed = TRUE)
  expect_error(parse_spec(c("--contiguity=none", "--major=joint")), "needs --distance=capital", fixed = TRUE)
  expect_error(parse_spec(c("--distance=none", "--contiguity=none", "--controls=full")), "needs --distance=capital", fixed = TRUE)
  s <- parse_spec(c("--distance=capital", "--contiguity=none")); expect_equal(s$contiguity, "none"); expect_equal(s$distance, "capital"); expect_false(spec_is_main(s))
  expect_equal(parse_spec(c("--distance=capital", "--contiguity=land"))$contiguity, "land")
  expect_equal(parse_spec("--contiguity=land"), main)                                              # the explicit default IS the main specification
  ## cell names: _nocontig right after _dist, before _major / _ix; the command-line order does not matter
  expect_equal(cell_name("fmid", 3L, "M1b", s), "fmid_lag3_all_pr_M1b_dist_nocontig")
  expect_equal(cell_name("mid", 3L, "M1b", parse_spec(c("--distance=capital", "--contiguity=none", "--major=joint"))), "mid_lag3_all_pr_M1b_dist_nocontig_major")
  expect_equal(cell_name("mid", 3L, "M1b", parse_spec(c("--major=joint", "--contiguity=none", "--distance=capital"))), "mid_lag3_all_pr_M1b_dist_nocontig_major")
  expect_equal(cell_name("mid", 3L, "M2b", parse_spec(c("--contiguity=none", "--distance=capital", "--interact=prefixed"))), "mid_lag3_all_pr_M2b_K2_dist_nocontig_ix")
  expect_equal(cell_name("mid", 3L, "M1b", parse_spec(c("--controls=full", "--distance=capital", "--contiguity=none", "--major=joint", "--interact=prefixed"))), "mid_lag3_all_pr_M1b_fullctrl_dist_nocontig_major_ix")
  all_on <- parse_spec(c("--pairs=originators", "--context=coalition", "--coal_ties=strict", "--coal_drop=top3", "--ww2=drop", "--controls=full",
                         "--distance=capital", "--contiguity=none", "--major=joint", "--interact=prefixed", "--weights=none", "--changes=winsor", "--direction=sideA"))
  expect_equal(spec_suffix(all_on), "_orig_co_strict_drop3_noww2_fullctrl_dist_nocontig_major_ix_w1_winsor_sideA")
  expect_equal(names(SPEC_SUFFIX), c("pairs", "context", "coal_ties", "coal_drop", "ww2", "controls", "distance", "contiguity", "major", "interact", "weights", "changes", "direction"))
  expect_equal(which(names(SPEC_SUFFIX) == "contiguity"), which(names(SPEC_SUFFIX) == "distance") + 1L)      # adjacent to _dist
  expect_false(any(duplicated(c(unlist(SPEC_SUFFIX), "_K2", "_pilot", "_STALE", "_samp10", "_noCH", "_polrel"))))
  ## every existing cell name is unchanged: the switch adds nothing at its default
  for (nm in setdiff(names(SPEC_CHOICES), "contiguity")) for (v in SPEC_CHOICES[[nm]][-1]) {
    a <- c(sprintf("--%s=%s", nm, v), if (nm %in% c("coal_drop", "coal_ties")) "--context=coalition")     # coal_drop needs the coalition context; coal_ties is inert without one
    sp <- parse_spec(a); sp2 <- parse_spec(c(a, "--contiguity=land"))
    expect_equal(cell_name("mid", 3L, "M1b", sp), cell_name("mid", 3L, "M1b", sp2)); expect_false(grepl("nocontig", cell_name("mid", 3L, "M1b", sp)))
  }
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec("--distance=capital")), "fmid_lag3_all_pr_M1b_dist")
  expect_equal(cell_name("mid", 3L, "M1b", parse_spec(c("--distance=capital", "--major=joint"))), "mid_lag3_all_pr_M1b_dist_major")
  ## the archived spellings --contiguity=TRUE / FALSE are refused as bad VALUES (parse_spec), the bare flag as a flag without a value; the option
  ## is no longer in ARCHIVED_OPTIONS (so --contiguity=land passes the archived-flag check), the other archived flags still are
  expect_error(parse_spec("--contiguity=TRUE"), "--contiguity must be one of land | none", fixed = TRUE)
  expect_error(parse_spec("--contiguity=FALSE"), "--contiguity must be one of land | none", fixed = TRUE)
  expect_error(parse_spec("--contiguity=1"), "--contiguity must be one of"); expect_error(parse_spec("--contiguity=off"), "--contiguity must be one of")
  expect_error(parse_spec("--contiguity"), "--contiguity needs a value (--contiguity=land | none); a bare flag is refused", fixed = TRUE)
  expect_error(parse_spec(c("--distance=capital", "--contiguity")), "--contiguity needs a value", fixed = TRUE)
  for (nm in names(SPEC_CHOICES)) expect_error(parse_spec(paste0("--", nm)), sprintf("--%s needs a value", nm), fixed = TRUE)   # the f4a8a98 refusal holds for every switch
  expect_false("contiguity" %in% ARCHIVED_OPTIONS); expect_false("distance" %in% ARCHIVED_OPTIONS)
  expect_equal(ARCHIVED_OPTIONS, c("variant", "sampled_m", "centrality", "conflict_history", "dy_peace_years", "drop_ww2", "full_controls", "event_weights", "K"))
  expect_true(refuse_archived_options(c("--contiguity=none", "--distance=capital"), "02_fit_model.R")); expect_true(refuse_archived_options("--contiguity=land", "02_fit_model.R"))
  for (old in ARCHIVED_OPTIONS) expect_error(refuse_archived_options(paste0("--", old, "=1"), "02_fit_model.R"), "is not an option")
  ## the build still takes --direction only, so --contiguity is refused there through the switch list
  expect_error(refuse_archived_options("--contiguity=none", "01_build_data.R", also = setdiff(names(SPEC_CHOICES), "direction")), "--contiguity is not an option of 01_build_data.R")
  ## the old global constant is gone: nothing but the switch can set the contiguity
  spec_src <- readLines(file.path(root, "R", "spec.R")); expect_false(any(grepl("^CONTIGUITY\\s*<-", spec_src)))
  ## (the scripts pass spec$contiguity; the logical idioms of the old flag are gone from assemble_cell.R, whose return list names its own CONTIGUITY argument)
  for (f in c("02_fit_model.R", "05_robustness.R")) expect_false(any(grepl("CONTIGUITY", readLines(file.path(root, "R", f)))), info = f)
  expect_false(any(grepl("if \\(CONTIGUITY\\)|!CONTIGUITY|CONTIGUITY &&|CONTIGUITY = TRUE|CONTIGUITY = FALSE", grep("^\\s*#", readLines(file.path(root, "R", "assemble_cell.R")), value = TRUE, invert = TRUE))))
  ## the preset table: the dist table with contiguity = none on top; same cells, same variant names, same order; no interactions, no war, no M2
  expect_true("dist_nocontig" %in% names(ROBUSTNESS_PRESETS)); expect_equal(ROBUSTNESS_PRESETS$dist_nocontig, list())
  expect_equal(DIST_NOCONTIG_PRESET$cells, DIST_PRESET$cells); expect_equal(DIST_NOCONTIG_PRESET$cells, c("fmid:3:M1b", "mid:3:M1b"))
  expect_equal(names(DIST_NOCONTIG_PRESET$variants), c("distance", "distance_major")); expect_equal(names(DIST_NOCONTIG_PRESET$variants), names(DIST_PRESET$variants))
  expect_equal(DIST_NOCONTIG_PRESET$variants$distance, list(distance = "capital", contiguity = "none"))
  expect_equal(DIST_NOCONTIG_PRESET$variants$distance_major, list(distance = "capital", contiguity = "none", major = "joint"))
  for (v in names(DIST_PRESET$variants)) expect_equal(DIST_NOCONTIG_PRESET$variants[[v]][names(DIST_PRESET$variants[[v]])], DIST_PRESET$variants[[v]])   # dist + contiguity = none, nothing else
  for (v in DIST_NOCONTIG_PRESET$variants) { for (nm in names(v)) expect_true(v[[nm]] %in% SPEC_CHOICES[[nm]]); expect_false("interact" %in% names(v)); expect_equal(v$contiguity, "none") }
  cells <- unlist(lapply(DIST_NOCONTIG_PRESET$variants, function(v) vapply(DIST_NOCONTIG_PRESET$cells, function(c) { p <- strsplit(c, ":")[[1]]
    s <- parse_spec(character(0)); for (nm in names(v)) s[[nm]] <- v[[nm]]; cell_name(p[1], as.integer(p[2]), p[3], s) }, "")))
  expect_equal(unname(cells), c("fmid_lag3_all_pr_M1b_dist_nocontig", "mid_lag3_all_pr_M1b_dist_nocontig", "fmid_lag3_all_pr_M1b_dist_nocontig_major", "mid_lag3_all_pr_M1b_dist_nocontig_major"))
  expect_equal(anyDuplicated(cells), 0L)
  dist_cells <- unlist(lapply(DIST_PRESET$variants, function(v) vapply(DIST_PRESET$cells, function(c) { p <- strsplit(c, ":")[[1]]
    s <- parse_spec(character(0)); for (nm in names(v)) s[[nm]] <- v[[nm]]; cell_name(p[1], as.integer(p[2]), p[3], s) }, "")))
  expect_false(any(cells %in% dist_cells)); expect_equal(unname(sub("_nocontig", "", cells, fixed = TRUE)), unname(dist_cells))    # cell-by-cell counterparts
  expect_false(any(cells %in% c(vapply(MAIN_SAMPLES, function(s) cell_name(s, 3L, "M1b", main), ""), "fmid_lag3_all_pr_M1b_major", "mid_lag3_all_pr_M1b_major_ix")))
})

test_that("column lists: contiguity = none drops contig_land and nothing else; log_capdist keeps its place and its centring; the logical flag is refused", {
  expect_equal(z_columns(), c("contig_land", "joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "dM_j_h", "dRpr_j_h"))
  expect_equal(z_columns(CONTIGUITY = "land"), z_columns())
  expect_equal(z_columns(DISTANCE = "capital", CONTIGUITY = "none"), c("joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "log_capdist", "dM_j_h", "dRpr_j_h"))
  expect_equal(z_columns(DISTANCE = "capital", CONTIGUITY = "none"), setdiff(z_columns(DISTANCE = "capital"), "contig_land"))
  expect_equal(z_columns(DISTANCE = "capital", CONTIGUITY = "none", MAJOR = "joint"), c("joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "log_capdist", "dM_j_h", "dRpr_j_h", "one_major", "both_major"))
  expect_equal(z_columns(DISTANCE = "capital", CONTIGUITY = "none", TARGET_MISSING = TRUE, MAJOR = "joint"), c("joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "log_capdist", "dM_j_h", "dRpr_j_h", "dch_j_missing", "one_major", "both_major"))
  expect_equal(z_columns(DISTANCE = "capital", CONTIGUITY = "none", FULL_CONTROLS = TRUE), c("joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "alliance_def_off", "log_dy_peace_years", "log_capdist", "dM_j_h", "dRpr_j_h"))
  expect_equal(z_columns(DISTANCE = "capital", CONTIGUITY = "none", INTERACT = "prefixed"), c("joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "log_capdist", "dM_j_h", "dRpr_j_h", Z_INTERACT))
  for (di in c("none", "capital")) for (fc in c(TRUE, FALSE)) for (mj in c("none", "joint")) for (ix in c("none", "prefixed")) for (tm in c(TRUE, FALSE)) {
    land <- z_columns(DISTANCE = di, CONTIGUITY = "land", FULL_CONTROLS = fc, MAJOR = mj, INTERACT = ix, TARGET_MISSING = tm)
    none <- z_columns(DISTANCE = di, CONTIGUITY = "none", FULL_CONTROLS = fc, MAJOR = mj, INTERACT = ix, TARGET_MISSING = tm)
    expect_equal(land[1], "contig_land"); expect_equal(none, land[-1])                      # exactly the first column leaves (the list itself does not enforce the distance)
  }
  expect_false("contig_land" %in% Z_CENTRED); expect_true("log_capdist" %in% Z_CENTRED)      # the 0/1 column was never centred; the distance is
  for (bad in list(TRUE, FALSE, "yes", c("land", "none"), 1L, NA)) expect_error(z_columns(CONTIGUITY = bad), info = deparse(bad))
  ## W is untouched by the switch (the contiguity is dyadic)
  expect_equal(w_controls(), c("M_lvl", "R_lvl", "ongoing_tm1", "log_prior_mids", "dem_tm1", "log_n_system_tm1"))
})

test_that("assembly on the toy cell: the nocontig Z is the distance Z without contig_land, everything else identical; refused without the distance; defaults identical", {
  td <- write_toy_cell(tempfile("toycell"))
  A  <- quiet(assemble_cell("mid", 3L, "all", data_dir = td))
  D  <- quiet(assemble_cell("mid", 3L, "all", DISTANCE = "capital", data_dir = td))
  N  <- quiet(assemble_cell("mid", 3L, "all", DISTANCE = "capital", CONTIGUITY = "none", data_dir = td))
  NM <- quiet(assemble_cell("mid", 3L, "all", DISTANCE = "capital", CONTIGUITY = "none", MAJOR = "joint", data_dir = td))
  DM <- quiet(assemble_cell("mid", 3L, "all", DISTANCE = "capital", MAJOR = "joint", data_dir = td))
  d <- function(X) build_stan_data(X$sy, X$alts, X$ev, X$W_cols, X$X_cols, X$Z_cols, joint = 1L, x_ref = X$x_ref, allow_orphan_initiations = TRUE)
  ## (a) defaults: the plain call, the explicit land and assemble_cell_spec() of the main spec are identical
  A0 <- quiet(assemble_cell("mid", 3L, "all", CONTIGUITY = "land", data_dir = td)); S0 <- quiet(assemble_cell_spec("mid", 3L, parse_spec(character(0)), data_dir = td))
  expect_identical(d(A), d(A0)); expect_identical(d(A), d(S0)); expect_identical(A$Z_cols, A0$Z_cols); expect_equal(A$contiguity, "land"); expect_equal(A$distance, "none")
  for (bad in list(TRUE, FALSE)) expect_error(quiet(assemble_cell("mid", 3L, "all", DISTANCE = "capital", CONTIGUITY = bad, data_dir = td)))   # the old logical flag is refused
  ## (b) none without the distance is refused by the assembly itself (parse_spec() refuses it earlier on the command line)
  expect_error(quiet(assemble_cell("mid", 3L, "all", CONTIGUITY = "none", data_dir = td)), "CONTIGUITY = none needs DISTANCE = capital", fixed = TRUE)
  expect_error(quiet(assemble_cell("mid", 3L, "all", DISTANCE = "none", CONTIGUITY = "none", MAJOR = "joint", data_dir = td)), "needs DISTANCE = capital", fixed = TRUE)
  ## (c) distance + none: Z = the distance Z minus contig_land (its first column), same rows, same events, same W / X, same distance statistics
  expect_equal(N$Z_cols, c("joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "log_capdist", "dM_j_h", "dRpr_j_h", "dch_j_missing")); expect_equal(N$Z_cols, D$Z_cols[-1]); expect_equal(D$Z_cols[1], "contig_land")
  expect_equal(N$contiguity, "none"); expect_equal(N$distance, "capital")
  expect_equal(N$W_cols, A$W_cols); expect_equal(N$X_cols, A$X_cols); expect_identical(N$ev, D$ev); expect_identical(N$sy, D$sy); expect_identical(N$ev, A$ev)
  expect_equal(nrow(N$alts), nrow(D$alts)); expect_identical(N$alts, D$alts)                 # the alternatives frame is the same (contig_land stays in it, out of Z)
  expect_identical(N$distance_stats, D$distance_stats)
  dN <- d(N); dD <- d(D); dA <- d(A)
  expect_equal(dN$Pz, dD$Pz - 1L); expect_equal(dN$Pz, dA$Pz); expect_equal(dN$E, dA$E); expect_equal(dN$N1, dA$N1); expect_equal(dN$N2, dA$N2)
  expect_identical(dN$Z, dD$Z[, -1, drop = FALSE]); expect_identical(dN$W, dD$W); expect_identical(dN$X, dD$X); expect_identical(dN$y1, dA$y1); expect_identical(dN$ev_w, dA$ev_w)
  expect_false("contig_land" %in% N$Z_cols); expect_true("contig_land" %in% names(N$alts))
  ## log_capdist: the same column as in the distance cell, centred over the IV rows
  iv <- N$alts$in_iv == 1L; alts0 <- readRDS(file.path(td, "alts_lag3_all.rds")); key <- function(a) paste(a$ccode1, a$ccode2)
  raw <- alts0$log_capdist[match(key(N$alts), key(alts0))]
  expect_equal(N$alts$log_capdist, raw - mean(raw[iv])); expect_lt(abs(mean(N$alts$log_capdist[iv])), 1e-12); expect_identical(N$alts$log_capdist, D$alts$log_capdist)
  expect_equal(unname(dN$Z[, match("log_capdist", N$Z_cols)]), N$alts$log_capdist); expect_lt(abs(mean(dN$Z[iv, match("log_capdist", N$Z_cols)])), 1e-12)
  expect_equal(N$distance_stats$cor_contig_iv, cor(raw[iv], N$alts$contig_land[iv]))        # still reported: what the distance replaces
  ## (d) distance + none + major: log_capdist, then the two indicators; W as in the major cell
  expect_equal(NM$Z_cols, c(N$Z_cols, "one_major", "both_major")); expect_equal(NM$Z_cols, DM$Z_cols[-1]); expect_equal(NM$W_cols, DM$W_cols); expect_equal(NM$W_cols, c(A$W_cols, "cowmaj_tm1"))
  dNM <- d(NM); dDM <- d(DM)
  expect_identical(dNM$Z, dDM$Z[, -1, drop = FALSE]); expect_identical(dNM$W, dDM$W); expect_equal(dNM$Pz, dDM$Pz - 1L); expect_equal(dNM$P1, dA$P1 + 1L)
  expect_identical(NM$distance_stats, D$distance_stats); expect_identical(NM$major_stats, DM$major_stats)
  ## assemble_cell_spec() passes the switch through, alone and with major
  S <- quiet(assemble_cell_spec("mid", 3L, parse_spec(c("--distance=capital", "--contiguity=none")), data_dir = td)); expect_identical(d(S), dN)
  SM <- quiet(assemble_cell_spec("mid", 3L, parse_spec(c("--distance=capital", "--contiguity=none", "--major=joint")), data_dir = td)); expect_identical(d(SM), dNM)
  ## (e) the NA-distance drop is the distance cell's under either contiguity value (same rows leave, same event counted)
  na_at <- function(pairs) function(a) { for (p in pairs) a$log_capdist[a$ccode1 == p[1] & a$ccode2 == p[2]] <- NA_real_; a }
  td2 <- write_toy_cell(tempfile("toycell"), alts_edit = na_at(list(c(40L, 20L), c(20L, 30L))))
  D2 <- quiet(assemble_cell("mid", 3L, "all", DISTANCE = "capital", data_dir = td2)); N2 <- quiet(assemble_cell("mid", 3L, "all", DISTANCE = "capital", CONTIGUITY = "none", data_dir = td2))
  expect_identical(N2$distance_stats, D2$distance_stats); expect_equal(N2$distance_stats$n_rows_dropped, 2L); expect_equal(N2$distance_stats$n_events_dropped, 1L)
  expect_identical(N2$alts, D2$alts); expect_identical(N2$ev, D2$ev); expect_identical(d(N2)$Z, d(D2)$Z[, -1, drop = FALSE]); expect_equal(N2$n_orphan, A$n_orphan + 1L)
  ## M1c (CONTEXT_MODEL) with the switch is not specified and refused
  expect_error(quiet(assemble_cell("mid", 3L, "all", DISTANCE = "capital", CONTIGUITY = "none", CONTEXT_MODEL = TRUE, data_dir = td)), "not specified")
})

test_that("robustness driver: the dist_nocontig preset expands to the four (variant, cell) calls of 02_fit_model.R; fixed switches refused; dist refuses --contiguity; perm_error and context_b refuse", {
  tmp <- tempfile("drv_v2_"); dir.create(file.path(tmp, "v2"), recursive = TRUE)
  file.copy(file.path(root, "R"), file.path(tmp, "v2"), recursive = TRUE); dir.create(file.path(tmp, "v2", "logs")); dir.create(file.path(tmp, "v2", "output"))
  rec <- file.path(tmp, "calls.txt")
  writeLines(c("args <- commandArgs(trailingOnly = TRUE)", sprintf("cat(paste(args, collapse = ' '), '\\n', sep = '', file = '%s', append = TRUE)", rec)), file.path(tmp, "v2", "R", "02_fit_model.R"))
  drv <- file.path(tmp, "v2", "R", "05_robustness.R")
  run <- function(...) { out <- suppressWarnings(system2("Rscript", c(drv, ...), stdout = TRUE, stderr = TRUE)); list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  rmrec <- function() if (file.exists(rec)) file.remove(rec)
  r <- run("--preset=dist_nocontig", "--dry_run=TRUE", "--threads_per_chain=2")
  expect_equal(r$status, 0L, info = paste(tail(r$out, 20), collapse = "\n"))
  calls <- readLines(rec); expect_length(calls, 4L)
  arg <- function(call, nm) unname(sub(paste0("^--", nm, "="), "", grep(paste0("^--", nm, "="), strsplit(call, " ")[[1]], value = TRUE)))
  args_of <- function(nm) unname(vapply(calls, arg, "", nm))
  ## variant slowest (distance, distance_major), cell fastest (fmid, mid) = the --cell_index order 0-3, the same as the dist preset's
  expect_equal(args_of("sample"), rep(c("fmid", "mid"), 2))
  expect_equal(args_of("distance"), rep("capital", 4)); expect_equal(args_of("contiguity"), rep("none", 4))
  expect_equal(args_of("major"), rep(c("none", "joint"), each = 2))
  expect_true(all(args_of("interact") == "none"))
  expect_true(all(args_of("lag") == "3")); expect_true(all(args_of("model") == "M1b")); expect_true(all(args_of("dry_run") == "TRUE")); expect_true(all(args_of("threads_per_chain") == "2"))
  for (nm in setdiff(names(SPEC_CHOICES), c("distance", "contiguity", "major"))) expect_true(all(args_of(nm) == SPEC_DEFAULTS[[nm]]), info = nm)
  for (cell in c("fmid_lag3_all_pr_M1b_dist_nocontig", "mid_lag3_all_pr_M1b_dist_nocontig", "fmid_lag3_all_pr_M1b_dist_nocontig_major", "mid_lag3_all_pr_M1b_dist_nocontig_major"))
    expect_true(any(grepl(paste0("cell ", cell, ":"), r$out, fixed = TRUE)), info = cell)
  ## one task of the Slurm array: --cell_index=3 is mid, distance_major; a pass-through switch adds its suffix in SPEC_SUFFIX order
  rmrec(); r3 <- run("--preset=dist_nocontig", "--cell_index=3", "--ww2=drop")
  expect_equal(r3$status, 0L); calls <- readLines(rec); expect_length(calls, 1L)
  expect_equal(arg(calls, "sample"), "mid"); expect_equal(arg(calls, "distance"), "capital"); expect_equal(arg(calls, "contiguity"), "none"); expect_equal(arg(calls, "major"), "joint"); expect_equal(arg(calls, "ww2"), "drop")
  expect_true(any(grepl("cell mid_lag3_all_pr_M1b_noww2_dist_nocontig_major:", r3$out, fixed = TRUE)))
  rmrec(); r0 <- run("--preset=dist_nocontig", "--cell_index=0"); expect_equal(r0$status, 0L); expect_true(any(grepl("cell fmid_lag3_all_pr_M1b_dist_nocontig:", r0$out, fixed = TRUE)))
  expect_true(run("--preset=dist_nocontig", "--cell_index=4")$status != 0L)
  for (bad in list("--cells=mid:3:M1b", "--distance=capital", "--distance=none", c("--distance=capital", "--contiguity=none"), "--contiguity=land", "--major=joint")) {
    rmrec(); rb <- run("--preset=dist_nocontig", bad); expect_true(rb$status != 0L, info = paste(bad, collapse = " "))
    expect_true(any(grepl("the dist_nocontig preset fixes --cells, --distance, --contiguity and --major", rb$out, fixed = TRUE)), info = paste(bad, collapse = " ")); expect_false(file.exists(rec))
  }
  ## --contiguity=none alone is refused earlier, by parse_spec() (it needs --distance=capital), whatever the preset
  rmrec(); rb <- run("--preset=dist_nocontig", "--contiguity=none"); expect_true(rb$status != 0L); expect_true(any(grepl("needs --distance=capital", rb$out, fixed = TRUE))); expect_false(file.exists(rec))
  ## the dist preset refuses --contiguity (its cells with contiguity = none ARE this preset's); its own four calls are unchanged (contiguity = land)
  rmrec(); rb <- run("--preset=dist", "--contiguity=land"); expect_true(rb$status != 0L); expect_true(any(grepl("the dist preset fixes --contiguity=land", rb$out, fixed = TRUE))); expect_false(file.exists(rec))
  rmrec(); rb <- run("--preset=dist", "--contiguity=none"); expect_true(rb$status != 0L); expect_true(any(grepl("needs --distance=capital", rb$out, fixed = TRUE))); expect_false(file.exists(rec))
  rmrec(); rd <- run("--preset=dist", "--dry_run=TRUE"); expect_equal(rd$status, 0L); calls <- readLines(rec); expect_length(calls, 4L)
  expect_true(all(args_of("contiguity") == "land")); expect_true(all(args_of("distance") == "capital"))
  for (cell in c("fmid_lag3_all_pr_M1b_dist:", "mid_lag3_all_pr_M1b_dist:", "fmid_lag3_all_pr_M1b_dist_major:", "mid_lag3_all_pr_M1b_dist_major:")) expect_true(any(grepl(paste0("cell ", cell), rd$out, fixed = TRUE)), info = cell)
  ## every other preset passes the default through as --contiguity=land (the switch list is the whole of SPEC_CHOICES); major_ix with the pair passes it through
  rmrec(); rn <- run("--preset=noww2", "--cells=fmid:3:M1b", "--dry_run=TRUE"); expect_equal(rn$status, 0L); expect_equal(arg(readLines(rec), "contiguity"), "land")
  rmrec(); rp <- run("--preset=noww2", "--cells=fmid:3:M1b", "--distance=capital", "--contiguity=none", "--dry_run=TRUE")
  expect_equal(rp$status, 0L); expect_true(any(grepl("cell fmid_lag3_all_pr_M1b_noww2_dist_nocontig:", rp$out, fixed = TRUE)))
  rmrec(); rm_ <- run("--preset=major_ix", "--cell_index=0", "--distance=capital", "--contiguity=none"); expect_equal(rm_$status, 0L)
  expect_true(any(grepl("cell fmid_lag3_all_pr_M1b_dist_nocontig_major:", rm_$out, fixed = TRUE)))
  ## --contiguity=none without --distance=capital is refused by parse_spec() before any call, through any preset
  rmrec(); rx <- run("--preset=noww2", "--contiguity=none"); expect_true(rx$status != 0L); expect_true(any(grepl("needs --distance=capital", rx$out, fixed = TRUE))); expect_false(file.exists(rec))
  ## context_b fixes every switch; perm_error refuses the switch before any fit or data file is touched
  expect_true(run("--preset=context_b", "--distance=capital", "--contiguity=none")$status != 0L)
  rb <- run("--preset=perm_error", "--family=m1", "--distance=capital", "--contiguity=none")
  expect_true(rb$status != 0L); expect_true(any(grepl("perm_error refuses --contiguity=none", rb$out, fixed = TRUE)), info = paste(tail(rb$out, 5), collapse = "\n"))
  ## the archived spellings and the bare flag are refused by the driver before any call
  for (bad in c("--contiguity=TRUE", "--contiguity=FALSE", "--contiguity")) { rmrec(); ra <- run("--preset=noww2", bad); expect_true(ra$status != 0L, info = bad); expect_false(file.exists(rec)) }
  expect_true(any(grepl("--contiguity needs a value", run("--preset=noww2", "--contiguity")$out, fixed = TRUE)))
  ## the sbatch dispatches the preset with the 0-3 range guard, and the dist case is unchanged
  sb <- readLines(file.path(root, "slurm", "04_robustness.sbatch"))
  expect_true(any(grepl("^  dist_nocontig)", sb))); expect_true(any(grepl("--preset=dist_nocontig --cell_index=\\$i", sb))); expect_true(any(grepl('\\[ "\\$i" -le 3 \\] \\|\\| \\{ echo "dist_nocontig: array index', sb)))
  expect_true(any(grepl("^  dist)", sb))); expect_true(any(grepl("--preset=dist --cell_index=\\$i", sb)))
  expect_true(any(grepl("major_ix\\|dist\\|dist_nocontig\\|perm_error", sb)))
  unlink(tmp, recursive = TRUE)
})

test_that("submit_robustness.sh: the dist_nocontig preset goes as ONE array 0-3 with no dependency; CELLS is refused; dist unchanged", {
  tmp <- tempfile("sub_"); dir.create(file.path(tmp, "v2"), recursive = TRUE); dir.create(file.path(tmp, "bin"))
  file.copy(file.path(root, "slurm"), file.path(tmp, "v2"), recursive = TRUE)
  rec <- file.path(tmp, "sbatch_calls.txt")
  writeLines(c("#!/bin/bash", sprintf("echo \"PRESET='${PRESET-<unset>}' CELLS='${CELLS-<unset>}' args: $*\" >> %s", rec), "echo \"3002;cluster\""), file.path(tmp, "bin", "sbatch"))
  Sys.chmod(file.path(tmp, "bin", "sbatch"), "0755")
  run <- function(env, ...) { out <- suppressWarnings(system2("bash", c("v2/slurm/submit_robustness.sh", ...), stdout = TRUE, stderr = TRUE,
                                                              env = c(sprintf("PATH=%s:%s", file.path(tmp, "bin"), Sys.getenv("PATH")), env)))
    list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  old <- setwd(tmp); on.exit(setwd(old), add = TRUE)
  r <- run("PRESET=dist_nocontig", "--account=open"); expect_equal(r$status, 0L, info = paste(r$out, collapse = "\n"))
  calls <- readLines(rec); expect_length(calls, 1L)
  expect_match(calls[1], "^PRESET='dist_nocontig' CELLS='<unset>' args: --parsable --export=ALL --array=0-3 --account=open v2/slurm/04_robustness.sbatch$")
  expect_true(any(grepl("dist_nocontig array submitted: job 3002 \\(indices 0-3: 0-1 distance, 2-3 distance_major; fmid, mid; contiguity=none; M1b only, no dependency\\)", r$out)))
  expect_true(run(c("PRESET=dist_nocontig", "CELLS=mid:3:M1b"))$status != 0L); expect_length(readLines(rec), 1L)
  expect_true(run("PRESET=dist_nocontig", "--array=0-1")$status != 0L); expect_length(readLines(rec), 1L)
  r2 <- run("PRESET=dist"); expect_equal(r2$status, 0L); calls <- readLines(rec); expect_length(calls, 2L)
  expect_match(calls[2], "^PRESET='dist' CELLS='<unset>' args: --parsable --export=ALL --array=0-3 v2/slurm/04_robustness.sbatch$")
})

test_that("real build (when present): the four nocontig cells are the dist cells without contig_land; the main spec and every existing preset are byte-identical to the saved dry runs", {
  dd <- file.path(v2_root(), "data")
  skip_if_not(file.exists(file.path(dd, "events.rds")) && file.exists(file.path(dd, "alts_lag3_all.rds")) && file.exists(file.path(dd, "manifest.csv")), "real data files not built")
  d <- function(X) build_stan_data(X$sy, X$alts, X$ev, X$W_cols, X$X_cols, X$Z_cols, joint = 1L, x_ref = X$x_ref, allow_orphan_initiations = TRUE, D_mode = 3L)
  man <- read.csv(file.path(dd, "manifest.csv"), stringsAsFactors = FALSE)
  ## a saved dry run (standata_<cell>.rds of 02_fit_model.R --dry_run=TRUE) is AVAILABLE when it was made from the current inputs (same manifest
  ## digests); every available one must be reproduced identically by the current code. Files from an earlier build are not comparable and are skipped.
  saved_available <- function(cell) {
    sd_path <- file.path(v2_root(), "output", sprintf("standata_%s.rds", cell))
    if (!file.exists(sd_path)) return(FALSE)
    sd <- readRDS(sd_path); on.exit(rm(sd))
    all(sd$manifest$sha256 == man$sha256[match(sd$manifest$file, man$file)])
  }
  check_saved <- function(X, cell) {
    sd <- readRDS(file.path(v2_root(), "output", sprintf("standata_%s.rds", cell)))
    dat <- build_stan_data(X$sy, X$alts, X$ev, X$W_cols, X$X_cols, X$Z_cols, joint = 1L, x_ref = X$x_ref, allow_orphan_initiations = TRUE, D_mode = 3L, grainsize = sd$dat$grainsize)
    expect_identical(dat, sd$dat, info = cell); expect_identical(X$Z_cols, sd$Z_cols, info = cell); expect_identical(X$W_cols, sd$W_cols, info = cell)
    rm(dat, sd); invisible(TRUE)
  }
  n_avail <- 0L; n_checked <- 0L
  for (s in c("fmid", "mid")) {
    main <- parse_spec(character(0))
    A <- quiet(assemble_cell_spec(s, 3L, main)); dA <- d(A)
    expect_equal(c(dA$P1, dA$P, dA$Pz), c(8L, 2L, 7L)); expect_equal(A$contiguity, "land"); expect_equal(A$Z_cols[1], "contig_land")
    if (saved_available(cell_name(s, 3L, "M1b", main))) { n_avail <- n_avail + 1L; n_checked <- n_checked + isTRUE(check_saved(A, cell_name(s, 3L, "M1b", main))) }
    ## every existing preset variant of this sample (M1b, lag 3) carries contiguity = land, no _nocontig in its name, and is byte-identical to its
    ## saved dry run where one is available from these inputs (the dist cells are checked below, with their nocontig counterparts)
    existing <- c(ROBUSTNESS_PRESETS[c("noww2", "noweights", "fullcontrols", "winsor")], CONTEXT_PRESET$variants, MAJOR_IX_PRESET$variants)
    for (v in names(existing)) {
      sp <- parse_spec(character(0)); for (nm in names(existing[[v]])) sp[[nm]] <- existing[[v]][[nm]]
      expect_equal(sp$contiguity, "land", info = v)
      cell <- cell_name(s, 3L, "M1b", sp); expect_false(grepl("nocontig", cell), info = v)
      if (!saved_available(cell)) next
      n_avail <- n_avail + 1L
      X <- quiet(assemble_cell_spec(s, 3L, sp)); n_checked <- n_checked + isTRUE(check_saved(X, cell)); rm(X); invisible(gc())
    }
    ## the four nocontig cells against their dist counterparts, assembled now
    for (v in names(DIST_NOCONTIG_PRESET$variants)) {
      sp <- parse_spec(character(0)); for (nm in names(DIST_PRESET$variants[[v]])) sp[[nm]] <- DIST_PRESET$variants[[v]][[nm]]
      sn <- parse_spec(character(0)); for (nm in names(DIST_NOCONTIG_PRESET$variants[[v]])) sn[[nm]] <- DIST_NOCONTIG_PRESET$variants[[v]][[nm]]
      expect_equal(cell_name(s, 3L, "M1b", sn), sub("_dist", "_dist_nocontig", cell_name(s, 3L, "M1b", sp), fixed = TRUE))
      D <- quiet(assemble_cell_spec(s, 3L, sp)); dD <- d(D)
      if (saved_available(cell_name(s, 3L, "M1b", sp))) { n_avail <- n_avail + 1L; n_checked <- n_checked + isTRUE(check_saved(D, cell_name(s, 3L, "M1b", sp))) }
      N <- quiet(assemble_cell_spec(s, 3L, sn)); dN <- d(N)
      expect_equal(N$Z_cols, D$Z_cols[-1], info = paste(s, v)); expect_equal(D$Z_cols[1], "contig_land"); expect_false("contig_land" %in% N$Z_cols)
      expect_equal(c(dN$P1, dN$P, dN$Pz), c(dD$P1, dD$P, dD$Pz - 1L), info = paste(s, v))
      expect_equal(c(dN$P1, dN$P, dN$Pz), list(distance = c(8L, 2L, 7L), distance_major = c(9L, 2L, 9L))[[v]], info = paste(s, v))
      expect_equal(dN$E, dA$E); expect_equal(dN$N1, dA$N1); expect_equal(dN$N2, dA$N2); expect_identical(dN$y1, dA$y1); expect_identical(dN$ev_w, dD$ev_w); expect_identical(dN$ev_row, dD$ev_row)
      expect_identical(dN$Z, dD$Z[, -1, drop = FALSE]); expect_identical(dN$W, dD$W); expect_identical(dN$X, dD$X)
      expect_identical(N$distance_stats, D$distance_stats); expect_equal(N$distance_stats$n_rows_dropped, 0L); expect_equal(N$distance_stats$n_events_dropped, 0L)
      iv <- N$alts$in_iv == 1L; expect_lt(abs(mean(dN$Z[iv, match("log_capdist", N$Z_cols)])), 1e-8)
      expect_true(is.finite(N$distance_stats$cor_contig_iv) && N$distance_stats$cor_contig_iv < 0)
      if (sn$major == "joint") { expect_equal(sum(N$W_cols == "cowmaj_tm1"), 1L); expect_equal(tail(N$Z_cols, 2), Z_MAJOR); expect_identical(N$major_stats, D$major_stats) }
      if (saved_available(cell_name(s, 3L, "M1b", sn))) { n_avail <- n_avail + 1L; n_checked <- n_checked + isTRUE(check_saved(N, cell_name(s, 3L, "M1b", sn))) }   # the saved nocontig dry runs, once written
      rm(D, N, dD, dN); invisible(gc())
    }
    rm(A, dA); invisible(gc())
  }
  expect_equal(n_checked, n_avail)     # every available saved dry run was compared (how many depends on which dry runs were saved from the current build)
})
