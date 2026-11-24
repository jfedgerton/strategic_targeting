## R/main_table.R -- the main-text table (PI 2026-11-24): stage 1, stage 2 and the joint two-stage model for the MAIN measure, one window x sample.
##   (1) Stage 1  : between-target conditional logit (no target FE, S1NF), the main stage-1 comparison; contrasts vs a peer target (gap 0)
##   (2) Stage 2  : the joining model alone (categorical A / B / O vs out, S2_REF_DESIGN), no target trait
##   (3) Joint    : two-stage FE model: stage-1 target FE (within targets), bbar = mean target FE of the dispute's onset targets, stage 2 with bbar
## Inference: the two-stage dispute bootstrap with the common draws of R/s2_07_boot_cells.R (seed 123). Columns (1) and (3) are read from that
## cell's draws; column (2) is fitted here on the same draws (same multiplicities, same seed), so all three use identical resamples.
## 95% percentile intervals; a quantity is reported only with a finite point estimate and >= 80% finite draws (as s2_07).
## Usage: Rscript R/main_table.R --build=DIR --dump=FILE --positions=FILE --spec=ID --edges_dir=DIR --cell_dir=DIR --out=DIR [--cores=2]
suppressPackageStartupMessages({ library(data.table); library(parallel) })
for (f in c("s2.R", "s2_fit.R", "sry.R", "s2_cell.R")) source(file.path("R/lib", f))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
BUILD <- opt("build"); DUMPF <- opt("dump"); POSF <- opt("positions"); SPEC <- opt("spec"); EDIR <- opt("edges_dir"); CELL <- opt("cell_dir"); OUT <- opt("out")
CORES <- as.integer(opt("cores", "2")); stopifnot(!is.null(BUILD), !is.null(DUMPF), !is.null(POSF), !is.null(SPEC), !is.null(EDIR), !is.null(CELL), !is.null(OUT))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); set.seed(123)
B <- s2_base(BUILD); W_REF <- B$window; dump <- readRDS(DUMPF)
pz <- readRDS(POSF); key <- as.data.table(pz$key); stopifnot(SPEC %in% colnames(pz$R)); R <- pz$R[, SPEC]
E <- s2_spec_edges("core", SPEC, EDIR); X <- s2_prepare(B, key, R, E)
idv <- readLines(file.path(CELL, "cell_identity.txt")); stopifnot(grepl(SPEC, idv, fixed = TRUE))
NB <- as.integer(trimws(strsplit(idv, "|", fixed = TRUE)[[1]][5]))
P <- fread(file.path(CELL, "point.csv"))
DRA <- unlist(lapply(sort(list.files(CELL, pattern = "^draws_.*\\.rds$", full.names = TRUE)), function(f) readRDS(f)$draws), recursive = FALSE)
stopifnot(length(DRA) == NB); okd <- vapply(DRA, function(x) is.numeric(x) && identical(names(x), boot_names(dump)), TRUE); MJ <- do.call(rbind, DRA[okd])
## (2) stage 2 alone: point and the same draws
F0 <- fit_mnl(X, S2_REF_DESIGN, "stage2 alone"); S0 <- surf_contrasts(F0)
tie0 <- c("A:tie_sum_A", "B:tie_sum_B", "O:tie_sum_A"); stopifnot(all(tie0 %in% names(F0$b)))
nm0 <- c(paste("S2A", S0$alt, S0$axis, S0$gap, sep = "."), paste0("S2A.", tie0))
pt0 <- setNames(c(S0$est, F0$b[tie0]), nm0)
disp_all <- sort(unique(c(as.data.table(dump$D)$disp, B$X$dispnum)))
set.seed(123); MM <- lapply(seq_len(NB), function(b) { m <- tabulate(sample.int(length(disp_all), length(disp_all), replace = TRUE), length(disp_all)); names(m) <- disp_all; m })
D0 <- mclapply(MM, function(m) tryCatch({ wt <- m[as.character(X$dispnum)]; Xb <- X[wt > 0]; w <- as.numeric(wt[wt > 0])
  f <- fit_mnl(Xb, S2_REF_DESIGN, "boot", vc = FALSE, w = w); s <- surf_contrasts(f); setNames(c(s$est, f$b[tie0]), c(paste("S2A", s$alt, s$axis, s$gap, sep = "."), paste0("S2A.", tie0)))[nm0] },
  error = function(e) NULL), mc.cores = CORES, mc.preschedule = FALSE)
M0 <- do.call(rbind, D0[!vapply(D0, is.null, TRUE)])
## assemble
ci <- function(M, k, est) { x <- M[is.finite(M[, k]), k]; if (!is.finite(est) || length(x) < 0.8 * NB) c(NA, NA) else unname(quantile(x, c(.025, .975))) }
row <- function(col, label, q, est, M) { b <- ci(M, q, est); data.table(column = col, row = label, est = est, lo = b[1], hi = b[2], n_draws = sum(is.finite(M[, q]))) }
getP <- function(st, al, ax, g) P[stage == st & alt == al & axis == ax & (is.na(g) & is.na(gap) | gap == g), est]
T <- list()
lab_s1 <- c("M.1" = "weaker target (material gap +1)", "M.-1" = "stronger target (material gap -1)", "R.1" = "less central target (relational gap +1)", "R.-1" = "more central target (relational gap -1)")
for (r in dump$roles) for (k in names(lab_s1)) { ax <- sub("\\..*$", "", k); g <- as.numeric(sub("^[MR]\\.", "", k))
  T[[length(T) + 1]] <- row("(1) Stage 1, between targets", paste(r, lab_s1[[k]]), paste("S1NF", r, ax, g, sep = "."), getP("S1NF", r, ax, g), MJ)
  T[[length(T) + 1]] <- row("(3) Joint: stage 1 within targets", paste(r, lab_s1[[k]]), paste("S1", r, ax, g, sep = "."), getP("S1", r, ax, g), MJ) }
lab_s2 <- c("M.1" = "attacker weaker (material gap +1)", "M.-1" = "attacker stronger (material gap -1)", "R.1" = "attacker less central (gap +1)", "R.-1" = "attacker more central (gap -1)")
for (al in c("A", "B", "O")) for (k in names(lab_s2)) { ax <- sub("\\..*$", "", k); g <- as.numeric(sub("^[MR]\\.", "", k)); q0 <- paste("S2A", al, ax, g, sep = ".")
  T[[length(T) + 1]] <- row("(2) Stage 2 alone", paste("join", al, lab_s2[[k]]), q0, pt0[[q0]], M0)
  T[[length(T) + 1]] <- row("(3) Joint: stage 2", paste("join", al, lab_s2[[k]]), paste("S2", al, ax, g, sep = "."), getP("S2", al, ax, g), MJ) }
for (i in 1:3) { al <- c("A", "B", "O")[i]
  T[[length(T) + 1]] <- row("(2) Stage 2 alone", paste("join", al, "recruitment ties"), paste0("S2A.", tie0[i]), pt0[[paste0("S2A.", tie0[i])]], M0)
  T[[length(T) + 1]] <- row("(3) Joint: stage 2", paste("join", al, "recruitment ties"), tie0[i], getP("S2", al, "tie_sum", NA), MJ)
  T[[length(T) + 1]] <- row("(3) Joint: stage 2", paste("join", al, "theta (target trait)"), paste0(al, ":bbar"), getP("S2", al, "theta", NA), MJ) }
TAB <- rbindlist(T)
TAB[, cell := ifelse(is.na(lo), sprintf("%.2f", est), sprintf("%.2f [%.2f, %.2f]%s", est, lo, hi, ifelse(lo > 0 | hi < 0, "*", "")))]
fwrite(TAB, file.path(OUT, "main_table_long.csv"))
N <- data.table(quantity = c("stage-1 events", "stage-2 candidate rows", "joined A", "joined B", "joined O", "bootstrap draws (joint usable)", "bootstrap draws (stage 2 alone usable)"),
                value = c(P$n_s1_events[1], P$n_s2_rows[1], P$n_s2_A[1], P$n_s2_B[1], P$n_s2_O[1], nrow(MJ), nrow(M0)))
fwrite(N, file.path(OUT, "main_table_n.csv"))
W <- dcast(TAB, row ~ column, value.var = "cell", fill = "")
fwrite(W, file.path(OUT, "main_table_wide.csv"))
cat(sprintf("MAIN_TABLE_DONE %s | joint draws %d, stage-2-alone draws %d of %d\n", OUT, nrow(MJ), nrow(M0), NB))
