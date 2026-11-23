## R/s2_07_boot_cells.R -- two-stage dispute bootstrap of the cells listed by R/s2_06_pick_boot_cells.R (stage-2 appendix, layer 2).
## Each draw resamples disputes and refits BOTH steps (R/lib/s2_cell.R joint_fe_draw; identical to R/s2_03_joint_fe.R). Draws are common
## across cells that share a build (seed 123; draw b is the same multiplicity vector for every cell of a window x sample, and the first 500
## draws of a 2000-draw cell are the 500 draws of the others). Resumable in blocks of draws; percentile 95% intervals.
## Usage: Rscript R/s2_07_boot_cells.R --cells=cells.csv --posdir=DIR --builds_root=DIR --dumps_dir=DIR --edges_dir=DIR --out=DIR [--cores=4] [--block=50] [--only=b001,b002]
suppressPackageStartupMessages({ library(data.table); library(parallel) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); RD <- dirname(normalizePath(this_file))
for (f in c("s2.R", "s2_fit.R", "sry.R", "s2_cell.R")) source(file.path(RD, "lib", f))
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

for (r in seq_len(nrow(CELLS))) {
  cl <- CELLS[r]; cd <- file.path(OUT, cl$cell_id); dir.create(cd, showWarnings = FALSE)
  if (file.exists(file.path(cd, "summary.csv"))) next
  B <- s2_base(file.path(BR, sprintf("w%d_%s", cl$window, cl$sample))); W_REF <- B$window; stopifnot(W_REF == cl$window)
  dump <- readRDS(file.path(DD, sprintf("long_init_bon%s_relyear_%s.rds", wtag(cl$window), cl$sample)))
  pz <- readRDS(file.path(PD, cl$file)); key <- as.data.table(pz$key); R <- pz$R[, cl$spec_id]; arm <- arm_of(cl$file)
  D <- s1_prepare(dump, key, R); X <- s2_prepare(B, key, R, s2_spec_edges(arm, cl$spec_id, EDIR))
  ref_order <- names(sort(table(D[chosen == 1L, c2]), decreasing = TRUE))
  if (!file.exists(file.path(cd, "point.csv"))) { P <- joint_fe_cell(D, dump, X, B$TS); fwrite(P, file.path(cd, "point.csv")) } else P <- fread(file.path(cd, "point.csv"))
  ## common draws: multiplicities over the build's disputes (measure-independent), generated sequentially from seed 123
  disp_all <- sort(unique(c(as.data.table(dump$D)$disp, B$X$dispnum)))
  set.seed(123); MM <- lapply(seq_len(cl$n_boot), function(b) { m <- tabulate(sample.int(length(disp_all), length(disp_all), replace = TRUE), length(disp_all)); names(m) <- disp_all; m })
  nb <- ceiling(cl$n_boot / BLOCK); t0 <- Sys.time()
  for (k in seq_len(nb)) { bf <- file.path(cd, sprintf("draws_%03d.rds", k)); if (file.exists(bf)) next
    ix <- ((k - 1) * BLOCK + 1):min(k * BLOCK, cl$n_boot)
    DR <- mclapply(MM[ix], function(m) tryCatch(joint_fe_draw(D, dump, X, B$TS, m, ref_order), error = function(e) conditionMessage(e)), mc.cores = CORES, mc.preschedule = FALSE)
    saveRDS(list(ix = ix, draws = DR), paste0(bf, ".tmp")); file.rename(paste0(bf, ".tmp"), bf) }
  DRA <- unlist(lapply(sort(list.files(cd, pattern = "^draws_.*\\.rds$", full.names = TRUE)), function(f) readRDS(f)$draws), recursive = FALSE)
  okd <- vapply(DRA, is.numeric, TRUE); M <- do.call(rbind, DRA[okd])
  pt <- c(setNames(P[stage == "S1", est], paste("S1", P[stage == "S1", alt], P[stage == "S1", axis], P[stage == "S1", gap], sep = ".")),
          setNames(P[stage == "S2" & axis %in% c("M", "R"), est], paste("S2", P[stage == "S2" & axis %in% c("M", "R"), alt], P[stage == "S2" & axis %in% c("M", "R"), axis], P[stage == "S2" & axis %in% c("M", "R"), gap], sep = ".")),
          setNames(P[axis == "theta", est], paste0(P[axis == "theta", alt], ":bbar")),
          setNames(P[axis == "tie_sum", est], c("A:tie_sum_A", "B:tie_sum_B", "O:tie_sum_A")))
  q <- intersect(boot_names(dump), names(pt))
  S <- data.table(cell_id = cl$cell_id, file = cl$file, spec_id = cl$spec_id, window = cl$window, sample = cl$sample, rule = cl$rule, q = q, est = pt[q],
                  lo = apply(M[, q, drop = FALSE], 2, quantile, .025, na.rm = TRUE), hi = apply(M[, q, drop = FALSE], 2, quantile, .975, na.rm = TRUE),
                  share_pos = colMeans(M[, q, drop = FALSE] > 0, na.rm = TRUE), n_ok = sum(okd), n_boot = cl$n_boot,
                  fail_msgs = paste(unique(unlist(DRA[!okd])), collapse = " | "))
  if (sum(okd) < 0.8 * cl$n_boot) S[, flag := "fewer than 80% usable draws"] else S[, flag := ""]
  fwrite(S, file.path(cd, "summary.csv"))
  logf("%s %s w%d %s: %d/%d draws ok, %.1f min | theta A %.2f [%.2f, %.2f]", cl$cell_id, cl$spec_id, cl$window, cl$sample, sum(okd), cl$n_boot,
       as.numeric(difftime(Sys.time(), t0, units = "mins")), S[q == "A:bbar", est], S[q == "A:bbar", lo], S[q == "A:bbar", hi])
}
fwrite(rbindlist(lapply(file.path(OUT, CELLS$cell_id, "summary.csv"), function(f) if (file.exists(f)) fread(f)), fill = TRUE), file.path(OUT, "boot_cells_summary.csv"))
logf("S2_BOOT_DONE")
