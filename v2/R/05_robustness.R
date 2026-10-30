## ============================================================================
## 05_robustness.R -- ONE driver for the robustness checks, as named presets (spec.R, ROBUSTNESS_PRESETS)
##
## Refit presets (a switch flipped on a list of cells; each cell is a call of 02_fit_model.R, so the fit carries the switch
## suffix and never overwrites a main-spec output):
##   noww2        ww2 = drop          the Second World War events and their initiations removed            suffix _noww2
##   sideA        direction = sideA   MID Side A as the initiator everywhere (needs v2/data_sideA, see below)   suffix _sideA
##   noweights    weights = none      every stage-2 event with weight 1                                       suffix _w1
##   fullcontrols controls = full     main + log_peace_years (W), alliance_def_off + log_dy_peace_years (Z)    suffix _fullctrl
##   winsor       changes = winsor    power changes winsorised at the 1st / 99th panel percentiles             suffix _winsor
## Context preset (dual typology stage a, 2026-10-25; CONTEXT_PRESET in spec.R): M1b, lag 3, samples mid and fmid, default tie convention,
##   three switch combinations: context = bilateral (_bi), context = coalition (_co), context = coalition + coal_drop = top3 (_co_drop3).
##   --cell_index = 0..5 maps variant (slowest) x cell: 0-1 bilateral, 2-3 coalition, 4-5 coalition_drop3; without it all six run in sequence.
##   Post-processing: 06_context_compare.R --sample=mid --lag=3 (needs the _bi and _co fits; _co_drop3 optional).
## Context_b preset (dual typology stage b, 2026-10-27; CONTEXT_B_PRESET in spec.R): the M1c model (one Stan program, m1_context.stan, with
##   group-specific ideal-point shifts and tolerance scales over bilateral / coalition leader / joiner events) at lag 3 on mid and fmid with
##   every switch at its default (the model itself is the variant; cells <sample>_lag3_all_pr_M1c). --cell_index = 0 mid, 1 fmid.
##   Post-processing: 03_extract.R --fit=v2/output/fit_<sample>_lag3_all_pr_M1c.rds (context_model_*.csv, typology_*_<group>.csv).
## major_ix preset (PI-approved stage-2 switches, 2026-10-27; MAJOR_IX_PRESET in spec.R): M1b, lag 3, samples fmid and mid, three switch
##   combinations: major = joint (_major), interact = prefixed (_ix), both (_major_ix). --cell_index = 0..5 maps variant (slowest) x cell:
##   0-1 major (fmid, mid), 2-3 interact, 4-5 both; without it all six run in sequence. M1b only: no M1 -> M2 dependency.
## dist preset (PI request 2026-10-28: distance and distance + major powers; contiguity stays in; DIST_PRESET in spec.R): M1b, lag 3, samples
##   fmid and mid, two switch combinations: distance = capital (_dist) and distance = capital + major = joint (_dist_major); no interactions.
##   --cell_index = 0..3 maps variant (slowest) x cell: 0-1 distance (fmid, mid), 2-3 distance_major; without it all four run in sequence. M1b only.
## dist_nocontig preset (PI request 2026-10-28, same day: the dist cells WITHOUT contiguity, distance replaces it; DIST_NOCONTIG_PRESET in spec.R):
##   the same four cells with contiguity = none on top, distance = capital + contiguity = none (_dist_nocontig) and + major = joint
##   (_dist_nocontig_major); the same --cell_index mapping 0-1 distance (fmid, mid), 2-3 distance_major, so the two presets compare cell by cell.
##   --cells, --distance, --contiguity and --major are refused (the dist preset refuses --contiguity as well: its cells with contiguity = none ARE
##   this preset's cells). M1b only.
## Test preset (no switch; the existing permutation and systematic-measurement-error tests on the main cell's a / b fits):
##   perm_error   the temporal holdout (B9, constants frozen on the training years; --cutoff), the permutation test (rescore =
##                reliance diagnostic, refit = calibrated p-value) and the relational-level error inflation, per family, on
##                ROBUSTNESS_MAIN_CELL (mid, lag 5); default --tests=holdout,permutation,inflation (04_robustness.sbatch runs
##                permutation + inflation in one task per family and the holdout in one task per cutoff (1990, 2000) and family)
##                Refused with --changes=winsor or --pairs != all (audit fix 15): the inflation and holdout refits rebuild the
##                changes and the frozen constants from the raw levels, which does not reproduce a winsorised cell, and the pairs
##                cells recompute y1 from a subset of the events; neither case has been checked for the refit-based tests.
##                Refused with --major != none or --interact != none as well (audit fix 4, 2026-10-27): the refits would leave the
##                ix_* products stale and cowmaj_tm1 un-re-centred. --distance=capital (2026-10-29, PI request: the permutation test on
##                the contiguity + distance main cell, fmid lag 3) is accepted for the PERMUTATION test only (--tests=permutation or
##                permutation_combine): the shuffle runs on the assembled frames as they are (the NA-distance row drop and the event
##                reconciliation already happened at assembly, and the _dist a / b fits are asserted to score exactly those frames), so
##                nothing is rebuilt from raw levels. The holdout and inflation tests still refuse it (they re-centre z_cent from A$Z_cols
##                but the NA-distance row drop and the event reconciliation have not been wired into their rebuilds). --contiguity != land
##                stays refused for every test (contiguity must stay in). The outputs carry the _dist cell name (measurement_cell_name()).
##                Parallel refit mode (2026-10-29): --perm_index=p (or --perm_from=p1 --perm_to=p2) with --tests=permutation
##                --perm_mode=refit refits ONLY those permutations (seed V2_SEED + p each, the same chains / iter / warmup) and writes one
##                row per permutation to output/measurement_<cell>_permutation_refit_parts/perm_<pppp>.csv (a, b, adv, seed, settings);
##                --tests=permutation_combine --perm_mode=refit --n_perm=N reads the parts 1..N (refuses a missing, duplicated or extra
##                index and parts whose settings differ), computes the observed statistic from the a / b fits as the sequential loop does,
##                and writes the SAME measurement_<cell>_permutation_refit.csv and _draws.csv with p = (1 + exceedances) / (1 + N).
##                The sequential loop (no --perm_index / --perm_from) is unchanged and gives the same numbers (tested on a toy cell).
##                --threads_per_chain=k (2026-10-30, PI request): k reduce_sum threads per refit chain (rstan_options(threads_per_chain = k),
##                grainsize = choose_grainsize(N1, k), as 02_fit_model.R; Slurm cpus = chains x k). Default 1 = the previous serial chains
##                (grainsize N1) and the previous output. k > 1 is taken with --perm_mode=refit and --tests=permutation or
##                permutation_combine only (the other tests and the rescore mode keep grainsize N1 and would ignore it). The part files
##                carry threads_per_chain and grainsize; the combiner refuses a part made with other values (reduce_sum sums the slices in
##                another order, so the statistics are not bit-identical across grainsizes), so pass the SAME --threads_per_chain to the
##                refit tasks and the combiner.
##                --major=joint (2026-10-30, PI request): accepted with --distance=capital for --tests=permutation / permutation_combine
##                only, like the distance: the shuffle acts on the assembled frames (cowmaj_tm1 centred, one_major / both_major in Z as
##                assembled) and the refits take them as they are, nothing is re-centred, and the _dist_major a / b fits
##                (fit_<s>_lag<h>_all_pr_M1a_dist_major / _M1b_dist_major) are asserted to score exactly those frames. The outputs carry
##                _dist_major (measurement_<s>_lag<h>_all_pr_M1_dist_major_*). --major without --distance=capital stays refused, as do the
##                holdout and inflation tests with it (their rebuilds re-centre W from w_controls() without the switch).
##
## Usage:
##   Rscript v2/R/05_robustness.R --preset=noww2 [--cells=mid:3:M1b,fmid:3:M1b,...] [--cell_index=i] [--dry_run=TRUE]
##   Rscript v2/R/05_robustness.R --preset=context [--cell_index=0..5] [--dry_run=TRUE] [--threads_per_chain=k]
##   Rscript v2/R/05_robustness.R --preset=context_b [--cell_index=0..1] [--dry_run=TRUE] [--threads_per_chain=k]
##   Rscript v2/R/05_robustness.R --preset=major_ix [--cell_index=0..5] [--dry_run=TRUE] [--threads_per_chain=k]
##   Rscript v2/R/05_robustness.R --preset=dist [--cell_index=0..3] [--dry_run=TRUE] [--threads_per_chain=k]
##   Rscript v2/R/05_robustness.R --preset=dist_nocontig [--cell_index=0..3] [--dry_run=TRUE] [--threads_per_chain=k]
##          [--threads_per_chain=k] [--chains= --iter= --warmup=] [other switches]
##       cells      : "sample:lag:model" list; default ROBUSTNESS_CELLS (spec.R). --cell_index=i (0-based) runs ONE cell of the
##                    list (the Slurm array of 04_robustness.sbatch); without it every cell runs in sequence.
##       other switches are passed through, so --preset=noww2 --controls=full fits the _noww2_fullctrl cells
##   Rscript v2/R/05_robustness.R --preset=perm_error --family=m1|m2 [--tests=holdout,permutation,inflation] [--n_perm=50]
##          [--perm_mode=rescore|refit] [--deltas=0.25,0.5,1.0] [--cutoff=1990] [--r0=0.5] [--chains=4] [--iter=1500] [--warmup=750]
##          [--sample=mid --lag=5] [--allow_stale=FALSE] [switches]
##   Rscript v2/R/05_robustness.R --preset=perm_error --family=m1 --sample=fmid --lag=3 --distance=capital --tests=permutation --perm_mode=rescore --n_perm=50
##   Rscript v2/R/05_robustness.R --preset=perm_error --family=m1 --sample=fmid --lag=3 --distance=capital --tests=permutation --perm_mode=refit --perm_index=7 --threads_per_chain=8
##   Rscript v2/R/05_robustness.R --preset=perm_error --family=m1 --sample=fmid --lag=3 --distance=capital --tests=permutation_combine --perm_mode=refit --n_perm=20 --threads_per_chain=8
##   Rscript v2/R/05_robustness.R --preset=perm_error --family=m1 --sample=fmid --lag=3 --distance=capital --major=joint --tests=permutation --perm_mode=refit --perm_index=7 --threads_per_chain=8
##       needs the production a / b fits of the cell (fit_<cell>_M1a / M1b or M2a_K<K> / M2b_K<K> with the switch suffix) from
##       02_fit_model.R; refuses pilots, stale fits and fits that do not score the assembled cell (C1 / C4)
##       writes v2/output/measurement_<cell>_*.csv with <cell> = <sample>_lag<h>_all_pr_<M1|M2>[_K<K>]<suffix> (measurement_cell_name(),
##       spec.R: the M2 name always carries its K, so --ame_K=3 never overwrites the K=2 files; review 2026-10-26 P2); every table
##       carries family and ame_K columns
## --direction=sideA (and the sideA preset) read v2/data_sideA, built by `Rscript v2/R/01_build_data.R --direction=sideA`.
## ============================================================================
suppressPackageStartupMessages({ library(dplyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
R_DIR <- dirname(this_file)
source(file.path(R_DIR, "spec.R")); source(file.path(R_DIR, "helpers_data.R")); source(file.path(R_DIR, "helpers_model.R"))
Sys.setenv(V2_ROOT = normalizePath(file.path(R_DIR, "..")))
args <- commandArgs(trailingOnly = TRUE)
refuse_archived_options(args, "05_robustness.R")         # the same archived flags 02_fit_model.R refuses (audit fix 14)
PRESET <- get_arg(args, "preset", "")
if (!PRESET %in% names(ROBUSTNESS_PRESETS)) stop("--preset must be one of ", paste(names(ROBUSTNESS_PRESETS), collapse = " | "), "; got '", PRESET, "'")
## the switches: whatever the command line says, with the preset's switch forced on top
spec <- parse_spec(args)
for (nm in names(ROBUSTNESS_PRESETS[[PRESET]])) spec[[nm]] <- ROBUSTNESS_PRESETS[[PRESET]][[nm]]
switch_args <- c(sprintf("--%s=%s", names(SPEC_CHOICES), unlist(spec[names(SPEC_CHOICES)])), sprintf("--ame_K=%d", spec$ame_K))

if (PRESET != "perm_error") {
  ## ---- refit presets: one 02_fit_model.R call per cell --------------------------------------------------------------------------
  if (PRESET == "context") {
    ## the context preset is a fixed table of (variant, cell): the three switch combinations of CONTEXT_PRESET x its two cells
    if (length(grep("^--(cells|context|coal_ties|coal_drop)=", args))) stop("the context preset fixes --cells, --context, --coal_ties and --coal_drop (CONTEXT_PRESET in spec.R); pass other switches only")
    parts <- strsplit(CONTEXT_PRESET$cells, ":", fixed = TRUE)
    cells_tab <- expand.grid(k = seq_along(CONTEXT_PRESET$cells), v = seq_along(CONTEXT_PRESET$variants))      # cell fastest, variant slowest
    cells_tab <- data.frame(sample = vapply(parts, `[`, "", 1)[cells_tab$k], lag = as.integer(vapply(parts, `[`, "", 2))[cells_tab$k], model = vapply(parts, `[`, "", 3)[cells_tab$k],
                            variant = names(CONTEXT_PRESET$variants)[cells_tab$v], stringsAsFactors = FALSE)
  } else if (PRESET == "context_b") {
    ## the M1c model on its fixed cells, every switch at its default: no --cells and no switch may be passed (the cell names would still be
    ## unique through the model code, but the stage-b run is defined on the main specification)
    if (length(grep("^--(cells|pairs|context|coal_ties|coal_drop|direction|ww2|weights|controls|major|interact|distance|contiguity|changes|ame_K)=", args)))
      stop("the context_b preset fixes its cells (CONTEXT_B_PRESET in spec.R) and runs M1c on the main specification; no --cells and no switch may be passed")
    parts <- strsplit(CONTEXT_B_PRESET$cells, ":", fixed = TRUE)
    cells_tab <- data.frame(sample = vapply(parts, `[`, "", 1), lag = as.integer(vapply(parts, `[`, "", 2)), model = vapply(parts, `[`, "", 3), variant = NA_character_, stringsAsFactors = FALSE)
    stopifnot(all(cells_tab$model %in% CONTEXT_MODELS))
  } else if (PRESET == "major_ix") {
    ## the major_ix preset (2026-10-27) is a fixed table as well: the three switch combinations of MAJOR_IX_PRESET x its two cells
    if (length(grep("^--(cells|major|interact)=", args))) stop("the major_ix preset fixes --cells, --major and --interact (MAJOR_IX_PRESET in spec.R); pass other switches only")
    parts <- strsplit(MAJOR_IX_PRESET$cells, ":", fixed = TRUE)
    cells_tab <- expand.grid(k = seq_along(MAJOR_IX_PRESET$cells), v = seq_along(MAJOR_IX_PRESET$variants))      # cell fastest, variant slowest
    cells_tab <- data.frame(sample = vapply(parts, `[`, "", 1)[cells_tab$k], lag = as.integer(vapply(parts, `[`, "", 2))[cells_tab$k], model = vapply(parts, `[`, "", 3)[cells_tab$k],
                            variant = names(MAJOR_IX_PRESET$variants)[cells_tab$v], stringsAsFactors = FALSE)
  } else if (PRESET == "dist") {
    ## the dist preset (2026-10-28) is a fixed table as well: the two switch combinations of DIST_PRESET x its two cells
    if (length(grep("^--(cells|distance|major)=", args))) stop("the dist preset fixes --cells, --distance and --major (DIST_PRESET in spec.R); pass other switches only")
    ## --contiguity is fixed too (2026-10-28, same day): dist with contiguity = none would be the dist_nocontig preset's cells under another route
    if (length(grep("^--contiguity=", args))) stop("the dist preset fixes --contiguity=land (its cells with contiguity = none are the dist_nocontig preset); pass other switches only")
    parts <- strsplit(DIST_PRESET$cells, ":", fixed = TRUE)
    cells_tab <- expand.grid(k = seq_along(DIST_PRESET$cells), v = seq_along(DIST_PRESET$variants))      # cell fastest, variant slowest
    cells_tab <- data.frame(sample = vapply(parts, `[`, "", 1)[cells_tab$k], lag = as.integer(vapply(parts, `[`, "", 2))[cells_tab$k], model = vapply(parts, `[`, "", 3)[cells_tab$k],
                            variant = names(DIST_PRESET$variants)[cells_tab$v], stringsAsFactors = FALSE)
  } else if (PRESET == "dist_nocontig") {
    ## the dist_nocontig preset (2026-10-28, same day) is the dist table with contiguity = none: the two switch combinations of DIST_NOCONTIG_PRESET x its two cells
    if (length(grep("^--(cells|distance|contiguity|major)=", args))) stop("the dist_nocontig preset fixes --cells, --distance, --contiguity and --major (DIST_NOCONTIG_PRESET in spec.R); pass other switches only")
    parts <- strsplit(DIST_NOCONTIG_PRESET$cells, ":", fixed = TRUE)
    cells_tab <- expand.grid(k = seq_along(DIST_NOCONTIG_PRESET$cells), v = seq_along(DIST_NOCONTIG_PRESET$variants))      # cell fastest, variant slowest
    cells_tab <- data.frame(sample = vapply(parts, `[`, "", 1)[cells_tab$k], lag = as.integer(vapply(parts, `[`, "", 2))[cells_tab$k], model = vapply(parts, `[`, "", 3)[cells_tab$k],
                            variant = names(DIST_NOCONTIG_PRESET$variants)[cells_tab$v], stringsAsFactors = FALSE)
  } else {
    cells <- strsplit(get_arg(args, "cells", paste(ROBUSTNESS_CELLS, collapse = ",")), ",")[[1]]
    parts <- strsplit(cells, ":", fixed = TRUE)
    if (any(lengths(parts) != 3L)) stop("--cells entries must be sample:lag:model, e.g. mid:3:M1b")
    cells_tab <- data.frame(sample = vapply(parts, `[`, "", 1), lag = as.integer(vapply(parts, `[`, "", 2)), model = vapply(parts, `[`, "", 3), variant = NA_character_, stringsAsFactors = FALSE)
  }
  stopifnot(all(cells_tab$sample %in% ALL_SAMPLES), all(cells_tab$lag %in% LAGS), all(cells_tab$model %in% if (PRESET == "context_b") CONTEXT_MODELS else MODELS))
  idx <- get_arg(args, "cell_index", -1)
  if (idx >= 0) { if (idx >= nrow(cells_tab)) stop(sprintf("--cell_index=%d is outside the %d-cell list", idx, nrow(cells_tab))); cells_tab <- cells_tab[idx + 1, , drop = FALSE] }
  log_open(paste0("05_robustness_", PRESET, if (idx >= 0) sprintf("_cell%d", idx) else ""))
  logf("preset %s -> switches %s (suffix %s); %d cell(s): %s", PRESET, spec_string(spec), spec_suffix(spec), nrow(cells_tab),
       paste(sprintf("%s:%d:%s%s", cells_tab$sample, cells_tab$lag, cells_tab$model, ifelse(is.na(cells_tab$variant), "", paste0(" [", cells_tab$variant, "]"))), collapse = " "))
  if (!PRESET %in% c("context", "context_b", "major_ix", "dist", "dist_nocontig") && spec_is_main(spec)) stop("the preset left every switch at its default; refusing to refit main-spec cells from the robustness driver")
  ## pass-through of the fit script's own options
  passthrough <- grep("^--(dry_run|threads_per_chain|chains|iter|warmup|adapt_delta|max_treedepth|pilot)=", args, value = TRUE)
  fit_script <- file.path(R_DIR, "02_fit_model.R")
  for (k in seq_len(nrow(cells_tab))) {
    cl <- cells_tab[k, ]
    spec_k <- spec
    if (PRESET == "context") for (nm in names(CONTEXT_PRESET$variants[[cl$variant]])) spec_k[[nm]] <- CONTEXT_PRESET$variants[[cl$variant]][[nm]]
    if (PRESET == "major_ix") for (nm in names(MAJOR_IX_PRESET$variants[[cl$variant]])) spec_k[[nm]] <- MAJOR_IX_PRESET$variants[[cl$variant]][[nm]]
    if (PRESET == "dist") for (nm in names(DIST_PRESET$variants[[cl$variant]])) spec_k[[nm]] <- DIST_PRESET$variants[[cl$variant]][[nm]]
    if (PRESET == "dist_nocontig") for (nm in names(DIST_NOCONTIG_PRESET$variants[[cl$variant]])) spec_k[[nm]] <- DIST_NOCONTIG_PRESET$variants[[cl$variant]][[nm]]
    if (spec_is_main(spec_k) && PRESET != "context_b") stop("refusing to refit a main-spec cell from the robustness driver")   # context_b: the M1c model IS the variant
    switch_args <- c(sprintf("--%s=%s", names(SPEC_CHOICES), unlist(spec_k[names(SPEC_CHOICES)])), sprintf("--ame_K=%d", spec_k$ame_K))
    cell <- cell_name(cl$sample, cl$lag, cl$model, spec_k)
    cmd_args <- c(fit_script, sprintf("--sample=%s", cl$sample), sprintf("--lag=%d", cl$lag), sprintf("--model=%s", cl$model), switch_args, passthrough)
    logf("[%d/%d] cell %s: Rscript %s", k, nrow(cells_tab), cell, paste(cmd_args, collapse = " "))
    status <- system2("Rscript", cmd_args)
    if (status != 0L) stop(sprintf("02_fit_model.R failed for cell %s (exit status %d)", cell, status))
    logf("[%d/%d] cell %s done", k, nrow(cells_tab), cell)
  }
  logf("05_robustness preset %s finished", PRESET)
  quit(save = "no", status = 0)
}

## ============================================================================
## perm_error: the permutation test and the systematic measurement error (plus the optional temporal holdout) on one cell's
## a / b fits. This is the former run_measurement_real.R (archive/R/), unchanged in substance.
## ============================================================================
suppressPackageStartupMessages({ library(rstan); library(loo) })
source(file.path(R_DIR, "typology.R")); source(file.path(R_DIR, "measurement_tests.R")); source(file.path(R_DIR, "assemble_cell.R"))
SAMPLE <- get_arg(args, "sample", ROBUSTNESS_MAIN_CELL$sample); LAG <- as.integer(get_arg(args, "lag", ROBUSTNESS_MAIN_CELL$lag))
FAMILY <- get_arg(args, "family", "m1"); K <- spec$ame_K; CUTOFF <- as.integer(get_arg(args, "cutoff", 1990))
N_PERM <- as.integer(get_arg(args, "n_perm", 50)); PERM_MODE <- get_arg(args, "perm_mode", "rescore")
DELTAS <- as.numeric(strsplit(get_arg(args, "deltas", "0.25,0.5,1.0"), ",")[[1]]); R0 <- get_arg(args, "r0", R0)
CHAINS <- as.integer(get_arg(args, "chains", 4)); ITER <- as.integer(get_arg(args, "iter", 1500)); WARMUP <- as.integer(get_arg(args, "warmup", 750))
TESTS <- strsplit(get_arg(args, "tests", "holdout,permutation,inflation"), ",")[[1]]   # the pre-restructure default (run_measurement_real.R); audit fix 6
ALLOW_STALE <- get_arg(args, "allow_stale", FALSE)
stopifnot(SAMPLE %in% ALL_SAMPLES, LAG %in% LAGS, all(TESTS %in% c("holdout", "permutation", "inflation", "permutation_combine")), FAMILY %in% c("m1", "m2"),
          PERM_MODE %in% c("rescore", "refit"), all(is.finite(DELTAS)), N_PERM >= 1, length(TESTS) >= 1, !anyDuplicated(TESTS))
## the parallel refit mode (2026-10-29): a subset of permutation indices per task, and the combiner; both refit mode + the permutation test alone
for (nm in c("perm_index", "perm_from", "perm_to")) if (length(grep(paste0("^--", nm, "$"), args))) stop(sprintf("--%s needs a value; a bare flag is refused", nm))
PERM_INDEX <- get_arg(args, "perm_index", NA_real_); PERM_FROM <- get_arg(args, "perm_from", NA_real_); PERM_TO <- get_arg(args, "perm_to", NA_real_)
if (!is.na(PERM_INDEX) && (!is.na(PERM_FROM) || !is.na(PERM_TO))) stop("--perm_index=p is the shorthand of --perm_from=p --perm_to=p; pass one form, not both")
if (!is.na(PERM_INDEX)) { PERM_FROM <- PERM_INDEX; PERM_TO <- PERM_INDEX }
if (xor(is.na(PERM_FROM), is.na(PERM_TO))) stop("--perm_from and --perm_to go together (or --perm_index=p for one permutation)")
PERM_SUBSET <- !is.na(PERM_FROM)
if (PERM_SUBSET) {
  if (PERM_FROM != as.integer(PERM_FROM) || PERM_TO != as.integer(PERM_TO) || PERM_FROM < 1 || PERM_TO < PERM_FROM) stop(sprintf("--perm_from=%s --perm_to=%s must be integers with 1 <= from <= to", PERM_FROM, PERM_TO))
  if (PERM_TO > N_PERM) stop(sprintf("--perm_to=%d exceeds --n_perm=%d (the combiner reads the parts 1..n_perm)", as.integer(PERM_TO), N_PERM))
  if (!identical(TESTS, "permutation") || PERM_MODE != "refit") stop("--perm_index / --perm_from / --perm_to need --tests=permutation --perm_mode=refit alone (one Slurm task per permutation; the other tests have their own tasks)")
  PERM_IDX <- seq.int(as.integer(PERM_FROM), as.integer(PERM_TO))
}
if ("permutation_combine" %in% TESTS) {
  if (!identical(TESTS, "permutation_combine") || PERM_MODE != "refit") stop("--tests=permutation_combine goes alone, with --perm_mode=refit (it combines the per-permutation refit parts)")
}
## within-chain threads for the permutation refits (2026-10-30): k >= 1; k > 1 only where a refit of the permutation test samples (refit mode,
## the permutation test or its combiner, which must carry the same setting to accept the parts); the other tests keep grainsize = N1
if (length(grep("^--threads_per_chain$", args))) stop("--threads_per_chain needs a value; a bare flag is refused")
THREADS <- get_arg(args, "threads_per_chain", 1)
if (length(THREADS) != 1 || !is.finite(THREADS) || THREADS < 1 || THREADS != as.integer(THREADS)) stop(sprintf("--threads_per_chain=%s must be a positive integer", paste(THREADS, collapse = ",")))
THREADS <- as.integer(THREADS)
if (THREADS > 1 && (PERM_MODE != "refit" || !all(TESTS %in% c("permutation", "permutation_combine"))))
  stop(sprintf("--threads_per_chain=%d applies to the permutation refits only: --perm_mode=refit with --tests=permutation or permutation_combine (got --tests=%s --perm_mode=%s); the holdout, inflation and rescore paths keep grainsize = N1", THREADS, paste(TESTS, collapse = ","), PERM_MODE))
## audit fix 15 (2026-10-26): the refit-based tests are not defined for winsorised changes or a pairs subset (see the header)
if (spec$changes == "winsor") stop("perm_error refuses --changes=winsor: the inflation and holdout refits rebuild the changes from the raw levels and would not reproduce the winsorised cell")
if (spec$pairs != "all") stop(sprintf("perm_error refuses --pairs=%s: the refit-based tests are defined on the full event sample (y1 of every initiating state-year)", spec$pairs))
## audit fix 4 (2026-10-27): the refits (holdout, inflation, permutation refit mode) rebuild the changes and re-centre the controls from the
## raw levels (freeze$ctrl_cent / z_cent below), which would leave the ix_* products of --interact stale and cowmaj_tm1 of --major un-re-centred
## --major=joint (2026-10-30, PI request): the permutation test only, and only on the distance cell (the _dist_major cells of the dist preset):
## like the distance, the shuffle and the refits take the assembled frames as they are (cowmaj_tm1 centred at assembly, one_major / both_major
## in Z), nothing is re-centred, and the _dist_major fits are asserted below to score exactly these frames. The holdout and inflation refits
## re-centre the W controls from w_controls() without the switch (cowmaj_tm1 would not be re-centred) and keep refusing it.
if (spec$major != "none" && !all(TESTS %in% c("permutation", "permutation_combine")))
  stop(sprintf("perm_error refuses --major=%s with --tests=%s: only the permutation test (--tests=permutation or permutation_combine) takes the major-power switch; the holdout and inflation refits re-centre the W controls from w_controls() without the switch, so cowmaj_tm1 would not be re-centred", spec$major, paste(TESTS, collapse = ",")))
if (spec$major != "none" && spec$distance != "capital")
  stop(sprintf("perm_error refuses --major=%s without --distance=capital: the permutation test takes the major-power switch on the distance cell only (the _dist_major a / b fits of the dist preset; 2026-10-30)", spec$major))
if (spec$contiguity != "land") stop(sprintf("perm_error refuses --contiguity=%s: contiguity stays in Z for every test of this preset (2026-10-29; the permutation test takes --distance=capital with contiguity = land only)", spec$contiguity))   # before the distance check, so the message names the switch
## --distance=capital (2026-10-29): the permutation test only. It shuffles dR within the assembled stage-2 sets and never rebuilds a column from
## the raw levels, so the NA-distance row drop and the event reconciliation of assembly are the data; the _dist fits are asserted below to score
## exactly these frames. The holdout and inflation refits still re-centre Z from the assembled columns without that wiring and keep refusing.
if (spec$distance != "none" && !all(TESTS %in% c("permutation", "permutation_combine")))
  stop(sprintf("perm_error refuses --distance=%s with --tests=%s: only the permutation test (--tests=permutation or permutation_combine) takes the distance cell; the holdout and inflation refits re-centre Z from the assembled columns but the NA-distance row drop and the event reconciliation have not been wired into them (2026-10-28)", spec$distance, paste(TESTS, collapse = ",")))
if (spec$interact != "none") stop(sprintf("perm_error refuses --interact=%s: the refits rebuild the changes but not the ix_* products in Z, which would stay at their original values (the switch has not been wired into the refit-based tests)", spec$interact))
cell <- measurement_cell_name(SAMPLE, LAG, FAMILY, spec)     # the M2 name carries its K (review 2026-10-26, P2); spec.R: the _dist cell gives measurement_<..>_M1_dist_*
log_open(paste0("05_robustness_perm_error_", cell, if (PERM_SUBSET) sprintf("_perm%d_%d", PERM_IDX[1], PERM_IDX[length(PERM_IDX)]) else if ("permutation_combine" %in% TESTS) "_combine" else "")); set.seed(V2_SEED)
rstan_options(auto_write = FALSE)
rstan_options(threads_per_chain = THREADS)                 # as 02_fit_model.R; 1 = the previous serial chains (the other tests keep grainsize = N1 anyway)
options(mc.cores = min(CHAINS, max(1L, n_cores() %/% THREADS)))
out_dir <- file.path(v2_root(), "output"); data_dir <- spec_data_dir(spec, v2_root())
if (!dir.exists(data_dir)) stop("data directory ", data_dir, " does not exist")
logf("perm_error on cell %s; tests %s; spec %s; data %s; threads_per_chain = %d; mc.cores = %d", cell, paste(TESTS, collapse = ","), spec_string(spec), basename(data_dir), THREADS, getOption("mc.cores"))
## compiled models, checked against the checked-out sources (C3)
models <- list(m1 = readRDS(file.path(v2_root(), "stan", "build", "m1_latent_additive.rds")),
               m2 = readRDS(file.path(v2_root(), "stan", "build", "m2_ame.rds")))
for (nm in names(models)) assert_compiled_current(models[[nm]], file.path(v2_root(), "stan", paste0(c(m1 = "m1_latent_additive", m2 = "m2_ame")[nm], ".stan")),
                                                  file.path(v2_root(), "stan", "build", paste0(c(m1 = "m1_latent_additive", m2 = "m2_ame")[nm], ".sha256")))
A <- assemble_cell_spec(SAMPLE, LAG, spec, data_dir = data_dir)
## the a / b fits of the cell, named by 02_fit_model.R (cell_name())
ma <- paste0(toupper(FAMILY), "a"); mb <- paste0(toupper(FAMILY), "b")    # cell_spec(): the M1 fits carry no K whatever --ame_K says (P2)
fa_path <- file.path(out_dir, paste0("fit_", cell_name(SAMPLE, LAG, ma, cell_spec(ma, spec)), ".rds"))
fb_path <- file.path(out_dir, paste0("fit_", cell_name(SAMPLE, LAG, mb, cell_spec(mb, spec)), ".rds"))
if (ALLOW_STALE) { if (!file.exists(fa_path)) fa_path <- sub("\\.rds$", "_STALE.rds", fa_path); if (!file.exists(fb_path)) fb_path <- sub("\\.rds$", "_STALE.rds", fb_path) }
if (!file.exists(fa_path) || !file.exists(fb_path)) stop("the a / b fits of the cell are missing: ", basename(fa_path), ", ", basename(fb_path), " (run 02_fit_model.R first)")
obj_a <- readRDS(fa_path); obj_b <- readRDS(fb_path)
## C1 / C2 / C4: production fits of the CURRENT data that score exactly this assembled cell
dat_chk <- build_stan_data(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 0L, x_ref = A$x_ref, K = if (FAMILY == "m2") K else NULL, allow_orphan_initiations = TRUE)
for (o in list(obj_a, obj_b)) {
  assert_production_fit(o, allow_stale = ALLOW_STALE)
  assert_fit_current(o, file.path(data_dir, "manifest.csv"))
  assert_fit_stan_current(o)                                   # sampled with the compiled program's source (audit fix 5)
  assert_fit_matches_cell(o, dat_chk, A$W_cols, A$X_cols, A$Z_cols)
  stopifnot(o$dat$D_mode == 3L)
}
stopifnot(obj_a$dat$joint == 0L, obj_b$dat$joint == 1L)
logf("fits %s and %s verified: production, current data, same E = %d, same events, same W / X / Z", basename(fa_path), basename(fb_path), dat_chk$E)
## every written table names the family and, for M2, the AME K it was computed with (metadata of the K-specific file names; P2)
with_meta <- function(tab) cbind(family = toupper(FAMILY), ame_K = if (FAMILY == "m2") K else NA_integer_, tab)
fit_a <- obj_a$fit; fit_b <- obj_b$fit
## positions (relational LEVELS) and the gap constant for the holdout freezing (B9) and the level perturbation (B8)
pos <- assert_matches_manifest(file.path(data_dir, "sy_positions.rds"), file.path(data_dir, "manifest.csv"))
levels <- pos |> transmute(state_idx = match(ccode, A$state_codes), year = as.integer(year), M, R = R_pr) |> filter(!is.na(state_idx))
scale_tbl <- read.csv(file.path(data_dir, "change_scale.csv"), stringsAsFactors = FALSE)
GAP_R <- scale_tbl$pooled_sd[is.na(scale_tbl$h) & scale_tbl$variable == "gap_R_pr"]; GAP_M <- scale_tbl$pooled_sd[is.na(scale_tbl$h) & scale_tbl$variable == "gap_M"]
stopifnot(length(GAP_R) == 1, length(GAP_M) == 1, is.finite(GAP_R), is.finite(GAP_M))
FULL <- spec$controls == "full"
## ---- temporal holdout (constants frozen on the training years, B9) ---------------------------------------------------------------
if ("holdout" %in% TESTS) {
  freeze <- list(positions = levels, gap_raw = c(dM = "dM_logratio", dR = "dR_pr_logratio"), levels = list(dM = "M_lvl", dR = "R_lvl"),
                 change_raw = c(dM_i = "dM_raw_h", dR_i = "dRpr_raw_h"), target_change = c(dM_i = "dM_j_h", dR_i = "dRpr_j_h"),
                 ctrl_cent = setdiff(union(w_controls(FULL_CONTROLS = FULL), x_controls(FULL_CONTROLS = FULL)), W_LEVELS), z_cent = intersect(A$Z_cols, Z_CENTRED))
  th <- temporal_holdout(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, cutoff = CUTOFF, family = FAMILY, K = K,
                         stan_models = models, chains = CHAINS, iter = ITER, warmup = WARMUP, x_ref = A$x_ref, freeze = freeze)
  write.csv(with_meta(th$summary), file.path(out_dir, sprintf("measurement_%s_temporal_holdout_%d.csv", cell, CUTOFF)), row.names = FALSE)
  logf("temporal holdout written (cutoff %d)", CUTOFF)
}
## ---- permutation (B7: refit mode = calibrated p-value; rescore mode = reliance diagnostic) ----------------------------------------
## the per-permutation parts of the parallel refit mode (2026-10-29): one row per permutation with its settings, in a per-cell directory
parts_dir <- file.path(out_dir, sprintf("measurement_%s_permutation_refit_parts", cell))
## the settings every part must share, and that the combiner checks against its own run: the cell and its spec, the family / K, the sampler
## settings of the refits (chains / iter / warmup and, since 2026-10-30, the reduce_sum threads per chain and the grainsize they imply), the
## data dimensions the fits were verified against, and the fit files
PERM_GRAINSIZE <- choose_grainsize(dat_chk$N1, THREADS)      # = the grainsize of every refit of this run (permutation_refit_one() recomputes and asserts it)
part_settings <- data.frame(cell = cell, family = toupper(FAMILY), ame_K = if (FAMILY == "m2") K else NA_integer_, sample = SAMPLE, lag = LAG,
                            spec = spec_string(spec), chains = CHAINS, iter = ITER, warmup = WARMUP, v2_seed = V2_SEED,
                            threads_per_chain = THREADS, grainsize = PERM_GRAINSIZE,
                            E = dat_chk$E, N1 = dat_chk$N1, sum_y1 = sum(dat_chk$y1), fit_a = basename(fa_path), fit_b = basename(fb_path),
                            statistic = "stage2_loo_elpd_of_refits", stringsAsFactors = FALSE)
if ("permutation" %in% TESTS && PERM_SUBSET) {
  ## ---- parallel refit mode, one task: the permutations PERM_IDX only; no summary, no p-value (the combiner does that) ----------------------
  obs <- permutation_observed(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, fit_a, fit_b, mode = "refit", family = FAMILY, K = K, x_ref = A$x_ref)
  stopifnot(obs$dat_a$E == dat_chk$E, obs$dat_a$N1 == dat_chk$N1)
  logf("permutation refit, permutations %d..%d of n_perm %d (seeds %d..%d); parts in %s; threads_per_chain %d, grainsize %d (N1 = %d); observed a %.3f b %.3f (for the log only, the combiner recomputes them)",
       PERM_IDX[1], PERM_IDX[length(PERM_IDX)], N_PERM, V2_SEED + PERM_IDX[1], V2_SEED + PERM_IDX[length(PERM_IDX)], basename(parts_dir), THREADS, PERM_GRAINSIZE, dat_chk$N1, obs$obs_a, obs$obs_b)
  for (p in PERM_IDX) {
    t0 <- Sys.time()
    one <- permutation_refit_one(p, A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, obs, family = FAMILY, stan_models = models, chains = CHAINS, iter = ITER, warmup = WARMUP, threads = THREADS)
    stopifnot(one$perm == p, one$seed == V2_SEED + p, one$threads_per_chain == THREADS, one$grainsize == PERM_GRAINSIZE)
    path <- write_permutation_part(parts_dir, one, part_settings)
    logf("permutation %d written: a %.3f b %.3f adv %.3f (%.1f min) -> %s", p, one$a, one$b, one$adv, as.numeric(difftime(Sys.time(), t0, units = "mins")), basename(path))
  }
  logf("permutation refit parts %d..%d done; combine with --tests=permutation_combine --perm_mode=refit --n_perm=%d once every part 1..%d exists", PERM_IDX[1], PERM_IDX[length(PERM_IDX)], N_PERM, N_PERM)
} else if ("permutation" %in% TESTS) {
  ## ---- the sequential test (unchanged): both modes, all n_perm permutations in this process -----------------------------------------------
  pt <- permutation_test(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, fit_a, fit_b, n_perm = N_PERM, mode = PERM_MODE,
                         family = FAMILY, K = K, stan_models = models, chains = CHAINS, iter = ITER, warmup = WARMUP, x_ref = A$x_ref, threads = THREADS)
  write.csv(with_meta(pt$summary), file.path(out_dir, sprintf("measurement_%s_permutation_%s.csv", cell, PERM_MODE)), row.names = FALSE)
  write.csv(data.frame(perm = seq_along(pt$permuted_advantage), adv = pt$permuted_advantage, a = pt$permuted_a, b = pt$permuted_b),
            file.path(out_dir, sprintf("measurement_%s_permutation_%s_draws.csv", cell, PERM_MODE)), row.names = FALSE)
  logf("permutation (%s) written: %s", PERM_MODE, pt$summary$interpretation)
}
if ("permutation_combine" %in% TESTS) {
  ## ---- the combiner: the parts 1..N_PERM plus the observed statistic from the a / b fits, the sequential output files --------------------
  parts <- read_permutation_parts(parts_dir, part_settings)       # refuses parts with other settings, index / seed mismatches, non-finite statistics
  logf("%d part(s) read from %s (threads_per_chain %d, grainsize %d, as this run): indices %s", nrow(parts), basename(parts_dir), THREADS, PERM_GRAINSIZE, paste(parts$perm, collapse = " "))
  obs <- permutation_observed(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, fit_a, fit_b, mode = "refit", family = FAMILY, K = K, x_ref = A$x_ref)
  pt <- permutation_combine(parts[, c("perm", "a", "b")], obs$obs_a, obs$obs_b, obs$statistic, n_perm = N_PERM, mode = "refit")
  stopifnot(pt$summary$n_perm == N_PERM, pt$summary$p_value == (1 + pt$summary$n_exceedances) / (1 + N_PERM))
  write.csv(with_meta(pt$summary), file.path(out_dir, sprintf("measurement_%s_permutation_%s.csv", cell, "refit")), row.names = FALSE)
  write.csv(data.frame(perm = seq_along(pt$permuted_advantage), adv = pt$permuted_advantage, a = pt$permuted_a, b = pt$permuted_b),
            file.path(out_dir, sprintf("measurement_%s_permutation_%s_draws.csv", cell, "refit")), row.names = FALSE)
  logf("permutation (refit, combined from %d parts in %s): observed advantage %.3f, exceedances %d, p = %.4f; %s", N_PERM, basename(parts_dir),
       pt$observed_advantage_b_minus_a, pt$summary$n_exceedances, pt$p_value, pt$summary$interpretation)
}
## ---- systematic measurement error (B8): inflate the relational LEVEL of state-years with recent conflict ---------------------------
## recent conflict at the POSITION-year: the state initiated in that year (its measured centrality that year is inflated)
if ("inflation" %in% TESTS) {
  recent <- A$sy |> transmute(state_idx, year, recent = y1)
  cls_b <- classify_types(mu_state_draws(fit_b), r0 = R0)
  me <- measurement_error_check(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, recent = recent, deltas = DELTAS,
                                baseline_cls = cls_b, levels = levels, h = LAG, gap_scale_R = GAP_R, family = FAMILY, joint = 1L, K = K, stan_models = models,
                                chains = CHAINS, iter = ITER, warmup = WARMUP, r0 = R0, x_ref = A$x_ref)
  write.csv(with_meta(me), file.path(out_dir, sprintf("measurement_%s_error_inflation.csv", cell)), row.names = FALSE)
  logf("measurement-error inflation written (deltas %s log points of R)", paste(DELTAS, collapse = ","))
}
logf("05_robustness preset perm_error finished")
