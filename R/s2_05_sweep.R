## R/s2_05_sweep.R -- the stage-2 parallel multiverse (PI 2026-11-22): the joint two-stage FE model (R/lib/s2_cell.R) for every relational
## measure in the given positions files, on one window x sample (one stage-2 build + its stage-1 init dump). Each measure uses its OWN
## positions (stage-1 dR, stage-2 surface gaps) and its OWN edges (stage-2 ties; option A, measure-specific ties).
## Point fits only (two-way clustered SEs); the bootstrapped sub-multiverse is R/s2_07_boot_cells.R (cells picked by R/s2_06_pick_boot_cells.R).
## Resumable: one CSV per (positions file x chunk); a chunk whose CSV exists is skipped. A failing cell is written as status ERROR, never dropped.
## Usage:
##   Rscript R/s2_05_sweep.R --build=DIR --dump=long_init_...rds --edges_dir=DIR --out=DIR --positions=f1.rds,f2.rds,...
##          [--validation=measure_validation.csv --set=all|plausible] [--chunk=50] [--cores=4] [--max_cells=0]
suppressPackageStartupMessages({ library(data.table); library(parallel) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); RD <- dirname(normalizePath(this_file))
for (f in c("s2.R", "s2_fit.R", "sry.R", "s2_cell.R")) source(file.path(RD, "lib", f))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
known <- c("build", "dump", "edges_dir", "out", "positions", "validation", "set", "chunk", "cores", "max_cells")
fl <- sub("^--([a-z_]+)=.*$", "\\1", args); if (any(!fl %in% known) || anyDuplicated(fl)) stop("unknown or duplicated option: ", paste(fl, collapse = " "))
BUILD <- opt("build"); DUMPF <- opt("dump"); EDIR <- opt("edges_dir"); OUT <- opt("out"); POSF <- strsplit(opt("positions", ""), ",")[[1]]
VAL <- opt("validation"); SET <- opt("set", "all"); CHUNK <- as.integer(opt("chunk", "50")); CORES <- as.integer(opt("cores", "4")); MAXC <- as.integer(opt("max_cells", "0"))
stopifnot(!is.null(BUILD), !is.null(DUMPF), !is.null(EDIR), !is.null(OUT), length(POSF) > 0, all(file.exists(POSF)), SET %in% c("all", "plausible"), SET == "all" || !is.null(VAL))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT, "s2_sweep.log")
logf <- function(...) { s <- sprintf(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOG, append = TRUE) }

B <- s2_base(BUILD); W_REF <- B$window; dump <- readRDS(DUMPF)
## the stage-1 dump must be the init frame of the same window and sample as the stage-2 build
cm <- if (is.null(dump$coorig_months)) 0L else as.integer(dump$coorig_months)
smp <- c(stage1 = "fmid", gml = "fmid", mid = "mid")[[B$meta$fatal_def]]
if (dump$types != "init" || cm != W_REF || dump$sample != smp) stop(sprintf("dump (%s, window %d, %s) does not match the build (init, window %d, %s)", dump$types, cm, dump$sample, W_REF, smp))
logf("stage-2 sweep %s | build %s (window %d, fatal_def %s) | dump %s | set %s | cores %d", format(Sys.time()), BUILD, W_REF, B$meta$fatal_def, basename(DUMPF), SET, CORES)
arm_of <- function(f) { b <- basename(f); if (grepl("^positions_igo_arm", b)) "igo" else if (grepl("^positions_core_", b)) "core" else if (grepl("^positions_dipl_", b)) "dipl" else stop("cannot tell arm of ", b) }
edge_file_of <- function(f) { b <- sub("\\.rds$", "", basename(f)); if (startsWith(b, "positions_core_G")) file.path(EDIR, "grid_layers.rds") else if (startsWith(b, "positions_core_BOOKc12_")) file.path(EDIR, "edges_book_c12.rds") else if (startsWith(b, "positions_core_BOOKraw_")) file.path(EDIR, "edges_book_raw.rds") else if (arm_of(f) == "igo") file.path(EDIR, "edges_igo_arm.rds") else file.path(EDIR, paste0(sub("^positions_", "edges_", b), ".rds")) }
FPF <- c(s2_code_files(RD), script = normalizePath(this_file), s2_build_files(BUILD, "build"), dump = DUMPF, setNames(POSF, paste0("pos:", basename(POSF))),
         setNames(vapply(POSF, edge_file_of, ""), paste0("edges:", basename(vapply(POSF, edge_file_of, "")))), if (!is.null(VAL)) c(validation = VAL))
s2_fingerprint_check(OUT, FPF[!duplicated(names(FPF))])
PLAUS <- if (SET == "plausible") { v <- fread(VAL); v[plausible == TRUE, paste(basename(file), spec_id)] } else NULL

one_cell <- function(spec, R, key, arm) {
  t0 <- proc.time()[[3]]
  r <- tryCatch({
    E <- s2_spec_edges(arm, spec, EDIR); D <- s1_add_tie(s1_prepare(dump, key, R), E, dump$roles); X <- s2_prepare(B, key, R, E)
    o <- joint_fe_cell(D, dump, X, B$TS); o[, `:=`(status = "OK", s2_rows_dropped = attr(X, "dropped"), msg = "")]; o
  }, error = function(e) data.table(stage = NA_character_, alt = NA_character_, axis = NA_character_, gap = NA_real_, est = NA_real_, se = NA_real_,
                                    status = if (grepl("^position constant", conditionMessage(e))) "SKIPPED" else if (grepl("fit invalid|no convergence", conditionMessage(e))) "INVALID" else "ERROR",
                                    msg = conditionMessage(e)))   # SKIPPED and INVALID are final results (as in R/07); ERROR is retried once
  r[, `:=`(spec_id = spec, secs = round(proc.time()[[3]] - t0, 1))]; r
}

n_done <- 0L
for (pf in POSF) {
  arm <- arm_of(pf); pz <- readRDS(pf); key <- as.data.table(pz$key); specs <- colnames(pz$R)
  if (!is.null(PLAUS)) specs <- specs[paste(basename(pf), specs) %in% PLAUS]
  if (MAXC > 0) specs <- head(specs, MAXC)
  if (!length(specs)) { logf("%s: no measures selected", basename(pf)); next }
  ch <- split(specs, ceiling(seq_along(specs) / CHUNK)); tag <- sub("\\.rds$", "", basename(pf))
  logf("%s (%s arm): %d measures in %d chunks", basename(pf), arm, length(specs), length(ch))
  rm(list = ls(.edge_cache), envir = .edge_cache); invisible(s2_spec_edges(arm, specs[1], EDIR))   # one edge file per positions file: load once, before forking
  run_specs <- function(sp) { res <- mclapply(sp, function(s) one_cell(s, pz$R[, s], key, arm), mc.cores = CORES, mc.preschedule = FALSE)
    bad <- vapply(res, function(x) !is.data.table(x), TRUE)
    if (any(bad)) res[bad] <- lapply(which(bad), function(i) data.table(status = "ERROR", msg = paste("worker failed:", as.character(res[[i]])), spec_id = sp[i]))
    rbindlist(res, fill = TRUE) }
  for (ci in seq_along(ch)) {
    f <- file.path(OUT, sprintf("cells_%s_w%d_%s_c%03d.csv", tag, W_REF, smp, ci)); fs <- sub("\\.csv$", ".specs", f)
    if (file.exists(f)) { if (!file.exists(fs) || !identical(readLines(fs), ch[[ci]])) stop("chunk ", basename(f), " exists but holds a different spec list (changed set, --chunk or --max_cells); use a new --out")
      prev <- fread(f); er <- unique(prev[status == "ERROR", spec_id]); if (!length(er)) next
      ## one retry of ERROR cells (e.g. a worker killed for memory); SKIPPED cells (constant positions) are final
      logf("  chunk %d: retrying %d ERROR cells", ci, length(er)); rr <- run_specs(er)
      ## an error that repeats with the same message is deterministic: FAILED (final, reported), not retried forever (recheck P1)
      pm <- prev[status == "ERROR", .(pmsg = msg[1]), by = spec_id]; rr <- merge(rr, pm, by = "spec_id", all.x = TRUE)
      rr[status == "ERROR" & !is.na(pmsg) & msg == pmsg, status := "FAILED"]; rr[, pmsg := NULL]
      rr[, `:=`(positions_file = basename(pf), arm = arm, window = W_REF, sample = smp, fatal_def = B$meta$fatal_def)]
      res <- rbind(prev[!spec_id %in% er], rr, fill = TRUE); s2_write_atomic(res, f); next }
    t0 <- Sys.time()
    res <- run_specs(ch[[ci]])
    res[, `:=`(positions_file = basename(pf), arm = arm, window = W_REF, sample = smp, fatal_def = B$meta$fatal_def)]
    writeLines(ch[[ci]], paste0(fs, ".tmp")); if (!file.rename(paste0(fs, ".tmp"), fs)) stop("could not write ", fs); s2_write_atomic(res, f)
    n_err <- uniqueN(res[status == "ERROR", spec_id]); n_done <- n_done + length(ch[[ci]])
    logf("  chunk %d/%d: %d measures, %d errors, %d skipped (constant positions), %d invalid fits, %.1f min", ci, length(ch), length(ch[[ci]]), n_err, uniqueN(res[status == "SKIPPED", spec_id]), uniqueN(res[status == "INVALID", spec_id]), as.numeric(difftime(Sys.time(), t0, units = "mins")))
  }
}
## completion (GPT P2-4): DONE only if no ERROR cell remains in this call's chunks; otherwise exit 3 ("finished with gaps"; rerun retries them)
ALLC <- rbindlist(lapply(list.files(OUT, pattern = "^cells_.*\\.csv$", full.names = TRUE), fread), fill = TRUE)
ALLC <- ALLC[positions_file %in% basename(POSF) & window == W_REF & sample == smp]
st <- ALLC[, .(status = status[1]), by = .(positions_file, spec_id)]
logf("census: %s", paste(sprintf("%s %d", names(table(st$status)), table(st$status)), collapse = ", "))
if (any(st$status == "ERROR")) { logf("S2_SWEEP_GAPS: %d ERROR cells remain (rerun to retry)", sum(st$status == "ERROR")); quit(status = 3) }
logf("S2_SWEEP_DONE (%d measures fitted this call)", n_done)
