## R/s2_07_boot_cells.R -- two-stage dispute bootstrap of the cells listed by R/s2_06_pick_boot_cells.R or R/grid_pick_boot_cells.R (stage-2 appendix, layer 2).
## Each draw resamples disputes and refits BOTH steps (R/lib/s2_cell.R joint_fe_draw; identical to R/s2_03_joint_fe.R). Draws are common
## across cells that share a build (seed 123; draw b is the same dispute-weight vector for every cell of a window x sample, and the first 500
## draws of a 2000-draw cell are the 500 draws of the others). Resumable in blocks of draws; percentile 95% intervals.
## Usage: Rscript R/s2_07_boot_cells.R --cells=cells.csv --posdir=DIR --builds_root=DIR --dumps_dir=DIR --edges_dir=DIR --out=DIR [--cores=4] [--block=50] [--only=b001,b002]
suppressPackageStartupMessages({ library(data.table); library(parallel) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); RD <- dirname(normalizePath(this_file))
for (f in c("s2.R", "s2_fit.R", "sry.R", "s2_cell.R", "grid.R")) source(file.path(RD, "lib", f))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
CELLS <- fread(opt("cells")); PD <- opt("posdir"); BR <- opt("builds_root"); DD <- opt("dumps_dir"); EDIR <- opt("edges_dir"); OUT <- opt("out")
CORES <- as.integer(opt("cores", "4")); BLOCK <- as.integer(opt("block", "50")); ONLY <- opt("only")
if (!is.null(ONLY)) CELLS <- CELLS[cell_id %in% strsplit(ONLY, ",")[[1]]]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); LOG <- file.path(OUT, "s2_boot.log")
logf <- function(...) { s <- sprintf(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOG, append = TRUE) }
wtag <- function(W) if (W > 0) sprintf("_w%dd", W) else ""
arm_of <- function(b) if (grepl("^positions_igo_arm", b)) "igo" else if (grepl("^positions_core_", b)) "core" else if (grepl("^positions_dipl_", b)) "dipl" else stop("arm of ", b)
logf("two-stage bootstrap of %d cells %s | cores %d | block %d", nrow(CELLS), format(Sys.time()), CORES, BLOCK)
## provenance (GPT P1-2): code, cell list and every input the listed cells use
edge_file_of <- function(f, spec) if (startsWith(f, "positions_core_G")) file.path(EDIR, "grid_layers.rds") else if (startsWith(f, "positions_core_BOOKc12_")) file.path(EDIR, "edges_book_c12.rds") else if (startsWith(f, "positions_core_BOOKraw_")) file.path(EDIR, "edges_book_raw.rds") else if (arm_of(f) == "igo") file.path(EDIR, "edges_igo_arm.rds") else file.path(EDIR, paste0(sub("^positions_", "edges_", sub("\\.rds$", "", f)), ".rds"))
bk <- unique(CELLS[, .(window, sample)])
FPF <- c(s2_code_files(RD), script = normalizePath(this_file),   # cells.csv itself is not fingerprinted: cell_identity.txt protects each cell, so a regenerated list can add cells
         unlist(lapply(seq_len(nrow(bk)), function(i) s2_build_files(file.path(BR, sprintf("w%d_%s", bk$window[i], bk$sample[i])), sprintf("build_w%d_%s", bk$window[i], bk$sample[i])))),
         unlist(lapply(seq_len(nrow(bk)), function(i) { f <- file.path(DD, sprintf("long_init_bon%s_relyear_%s.rds", wtag(bk$window[i]), bk$sample[i])); setNames(f, paste0("dump:", basename(f))) })),
         setNames(file.path(PD, unique(CELLS$file)), paste0("pos:", unique(CELLS$file))),
         { e <- unique(vapply(unique(CELLS$file), edge_file_of, "", spec = "")); setNames(e, paste0("edges:", basename(e))) })
s2_fingerprint_check(OUT, FPF)
QMIN <- 0.8   # prespecified usable-draw threshold, per quantity

run_cell <- function(cl) {
  cd <- file.path(OUT, cl$cell_id); dir.create(cd, showWarnings = FALSE)
  ## a cell folder belongs to exactly one (file, spec, window, sample, n_boot, block); a regenerated cells.csv or changed --block stops loudly
  idf <- file.path(cd, "cell_identity.txt"); idv <- paste(cl$file, cl$spec_id, cl$window, cl$sample, cl$n_boot, BLOCK, sep = " | ")
  if (file.exists(idf)) { if (!identical(readLines(idf), idv)) stop(sprintf("%s holds a different cell (%s); use a new --out", cd, readLines(idf))) } else writeLines(idv, idf)
  sf <- file.path(cd, "summary.csv")   # GPT P2-3: a finished cell is final only if no draw was lost to a worker; otherwise its blocks are retried below
  if (file.exists(sf)) { S0 <- fread(sf); if (isTRUE(all(S0$n_draws_other_fail == 0))) return(invisible(NULL))
    logf("%s: %d draws were lost to workers; retrying them", cl$cell_id, S0$n_draws_other_fail[1]); unlink(sf) }
  B <- s2_base(file.path(BR, sprintf("w%d_%s", cl$window, cl$sample))); W_REF <<- B$window; stopifnot(W_REF == cl$window)
  dump <- readRDS(file.path(DD, sprintf("long_init_bon%s_relyear_%s.rds", wtag(cl$window), cl$sample)))
  cm <- if (is.null(dump$coorig_months)) 0L else as.integer(dump$coorig_months); smp <- c(stage1 = "fmid", gml = "fmid", mid = "mid")[[B$meta$fatal_def]]
  if (dump$types != "init" || cm != cl$window || dump$sample != cl$sample || smp != cl$sample) stop("dump/build/cell mismatch for ", cl$cell_id)
  pz <- readRDS(file.path(PD, cl$file)); key <- as.data.table(pz$key); R <- pz$R[, cl$spec_id]; arm <- arm_of(cl$file)
  if (startsWith(cl$file, "positions_core_G")) grid_check_provenance(pz, file.path(EDIR, "grid_layers.rds"))   # GPT P2-5
  E <- s2_spec_edges(arm, cl$spec_id, EDIR); D <- s1_add_tie(s1_prepare(dump, key, R), E, dump$roles); X <- s2_prepare(B, key, R, E)
  ref_order <- names(sort(table(D[chosen == 1L, c2]), decreasing = TRUE))
  pf <- file.path(cd, "point.csv")
  if (!file.exists(pf)) { P <- joint_fe_cell(D, dump, X, B$TS); s2_write_atomic(P, pf) } else P <- fread(pf)   # an invalid point fit stops here (cell INVALID)
  ## common draws: Bayesian-bootstrap dispute weights over the build's disputes (measure-independent), generated sequentially from seed 123 (boot_weights)
  disp_all <- sort(unique(c(as.data.table(dump$D)$disp, B$X$dispnum)))
  MM <- boot_weights(disp_all, cl$n_boot)   # Bayesian bootstrap (lib/s2_cell.R)
  nb <- ceiling(cl$n_boot / BLOCK); t0 <- Sys.time()
  for (k in seq_len(nb)) { bf <- file.path(cd, sprintf("draws_%03d.rds", k)); ix <- ((k - 1) * BLOCK + 1):min(k * BLOCK, cl$n_boot)
    one <- function(m) tryCatch(joint_fe_draw(D, dump, X, B$TS, m, ref_order), error = function(e) paste("numerical:", conditionMessage(e)))
    if (file.exists(bf)) { old <- readRDS(bf); if (!identical(old$ix, ix)) stop("draw block ", bf, " does not hold draws ", min(ix), "-", max(ix))
      ## GPT P2-4: draws lost to a worker (killed fork: try-error / NULL) are retried; numerical failures are final
      rt <- which(!vapply(old$draws, function(x) (is.numeric(x) && identical(names(x), boot_names(dump))) || (is.character(x) && !inherits(x, "try-error") && startsWith(x, "numerical:")), TRUE))
      if (!length(rt)) next
      old$draws[rt] <- mclapply(MM[ix[rt]], one, mc.cores = CORES, mc.preschedule = FALSE); DR <- old$draws
    } else DR <- mclapply(MM[ix], one, mc.cores = CORES, mc.preschedule = FALSE)
    saveRDS(list(ix = ix, draws = DR), paste0(bf, ".tmp")); if (!file.rename(paste0(bf, ".tmp"), bf)) stop("could not write ", bf) }
  DRA <- unlist(lapply(sort(list.files(cd, pattern = "^draws_.*\\.rds$", full.names = TRUE)), function(f) readRDS(f)$draws), recursive = FALSE)
  stopifnot(length(DRA) == cl$n_boot)
  okd <- vapply(DRA, function(x) is.numeric(x) && identical(names(x), boot_names(dump)), TRUE)
  numf <- vapply(DRA, function(x) is.character(x) && !inherits(x, "try-error") && startsWith(x, "numerical:"), TRUE)
  if (!any(okd)) stop("no usable draws for ", cl$cell_id, ": ", paste(head(unique(unlist(DRA[numf])), 3), collapse = " | "))
  M <- do.call(rbind, DRA[okd])
  pt <- c(setNames(P[stage == "S1", est], paste("S1", P[stage == "S1", alt], P[stage == "S1", axis], P[stage == "S1", gap], sep = ".")),
          setNames(P[stage == "S1NF", est], paste("S1NF", P[stage == "S1NF", alt], P[stage == "S1NF", axis], P[stage == "S1NF", gap], sep = ".")),
          setNames(P[stage == "S1NFT" & axis != "tie", est], paste("S1NFT", P[stage == "S1NFT" & axis != "tie", alt], P[stage == "S1NFT" & axis != "tie", axis], P[stage == "S1NFT" & axis != "tie", gap], sep = ".")),
          setNames(P[stage == "S1NFT" & axis == "tie", est], paste0("S1NFT.", P[stage == "S1NFT" & axis == "tie", alt], ".tie")),
          setNames(P[stage == "S2" & axis %in% c("M", "R"), est], paste("S2", P[stage == "S2" & axis %in% c("M", "R"), alt], P[stage == "S2" & axis %in% c("M", "R"), axis], P[stage == "S2" & axis %in% c("M", "R"), gap], sep = ".")),
          setNames(P[axis == "theta", est], paste0(P[axis == "theta", alt], ":bbar")),
          setNames(P[axis == "tie_sum", est], c("A:tie_sum_A", "B:tie_sum_B", "O:tie_sum_A")))
  q <- boot_names(dump); stopifnot(all(q %in% names(pt)))
  ## GPT P1-3: usable draws counted PER QUANTITY; an interval is reported only with a finite point estimate and >= QMIN x n_boot finite draws
  nf <- colSums(is.finite(M[, q, drop = FALSE]))
  qok <- is.finite(pt[q]) & nf >= QMIN * cl$n_boot
  qf <- function(k, pr) if (qok[k]) unname(quantile(M[is.finite(M[, k]), k], pr)) else NA_real_
  S <- data.table(cell_id = cl$cell_id, file = cl$file, spec_id = cl$spec_id, window = cl$window, sample = cl$sample, rule = cl$rule, q = q, est = unname(pt[q]),
                  lo = vapply(q, qf, 0, pr = .025), hi = vapply(q, qf, 0, pr = .975),
                  share_pos = vapply(q, function(k) if (qok[k]) mean(M[is.finite(M[, k]), k] > 0) else NA_real_, 0),
                  n_finite = unname(nf), n_draws_ok = sum(okd), n_draws_numerical_fail = sum(numf), n_draws_other_fail = sum(!okd & !numf), n_boot = cl$n_boot,
                  flag = ifelse(qok, "", ifelse(is.finite(pt[q]), sprintf("fewer than %d%% finite draws", round(100 * QMIN)), "no finite point estimate")),
                  fail_msgs = paste(head(unique(unlist(DRA[numf])), 5), collapse = " | "))
  s2_write_atomic(S, file.path(cd, "summary.csv")); unlink(file.path(cd, "error.csv"))
  logf("%s %s w%d %s: %d/%d draws ok, %.1f min | theta A %.2f [%.2f, %.2f]", cl$cell_id, cl$spec_id, cl$window, cl$sample, sum(okd), cl$n_boot,
       as.numeric(difftime(Sys.time(), t0, units = "mins")), S[q == "A:bbar", est], S[q == "A:bbar", lo], S[q == "A:bbar", hi])
}
## a failing cell is recorded (ERROR summary) and the run moves on; rerunning retries it
for (r in seq_len(nrow(CELLS))) { cl <- CELLS[r]
  tryCatch(run_cell(cl), error = function(e) { m <- conditionMessage(e)
    if (grepl("holds a different|does not hold draws|mismatch|could not write", m)) stop(m)   # identity, resume and I/O errors stop the run
    stt <- if (grepl("^position constant", m)) "SKIPPED" else if (grepl("fit invalid|no convergence|no usable draws", m)) "INVALID" else "ERROR"
    ef <- file.path(OUT, cl$cell_id, "error.csv")   # an ERROR repeating with the same message is deterministic: FAILED (final, reported)
    if (stt == "ERROR" && file.exists(ef) && identical(fread(ef)$msg[1], m)) stt <- "FAILED"
    logf("%s %s: %s %s", cl$cell_id, cl$spec_id, stt, m)
    s2_write_atomic(data.table(cell_id = cl$cell_id, file = cl$file, spec_id = cl$spec_id, window = cl$window, sample = cl$sample, rule = cl$rule, status = stt, msg = m),
                    file.path(OUT, cl$cell_id, "error.csv")) }) }
ALL <- rbindlist(lapply(file.path(OUT, CELLS$cell_id), function(d) if (file.exists(file.path(d, "summary.csv"))) fread(file.path(d, "summary.csv"))[, status := "OK"]
  else if (file.exists(file.path(d, "error.csv"))) fread(file.path(d, "error.csv")) else data.table(cell_id = basename(d), status = "ERROR", msg = "not run")), fill = TRUE)
s2_write_atomic(ALL, file.path(OUT, "boot_cells_summary.csv"))
st <- ALL[, .(status = status[1]), by = cell_id]; logf("census: %s", paste(sprintf("%s %d", names(table(st$status)), table(st$status)), collapse = ", "))
## GPT P2-8: every bootstrapped cell is prespecified, so any status other than OK is a gap (ERROR is retried on rerun; INVALID / FAILED / SKIPPED are reported)
if (any(st$status != "OK")) { logf("S2_BOOT_GAPS: %d cells not OK (%s)", sum(st$status != "OK"), paste(sprintf("%s %s", st[status != "OK", cell_id], st[status != "OK", status]), collapse = ", ")); quit(status = 3) }
logf("S2_BOOT_DONE")
