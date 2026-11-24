## R/main_table.R -- the main-text table (PI 2026-11-24): stage 1, stage 2 and the joint two-stage model for the MAIN measure, one window x sample.
##   (1) Stage 1  : between-target conditional logit (no target FE, S1NF), the main stage-1 comparison; contrasts vs a peer target (gap 0)
##   (2) Stage 2  : the joining model alone (categorical A / B / O vs out, S2_REF_DESIGN), no target trait, on the joint model's disputes
##                  (those with a finite bbar), so (2) vs (3) differs only by the target trait (PI 2026-11-24)
##   (3) Joint    : two-stage FE model: stage-1 target FE (within targets), bbar = mean target FE of the dispute's onset targets, stage 2 with bbar
## Inference: the two-stage dispute bootstrap with the common draws of R/s2_07_boot_cells.R (seed 123). Columns (1) and (3) are read from that
## cell's draws; column (2) is fitted here on the same draws (same Bayesian-bootstrap dispute weights, same seed; boot_weights), so all three use identical draws.
## 95% percentile intervals; a quantity is reported only with a finite point estimate and >= 80% finite draws (as s2_07).
## Missing-data robustness (PI 2026-11-24), written to main_table_robust_*.csv on the first N_ROBUST draws of the same weights:
##   (2b) stage 2 alone on all stage-2 disputes (incl. those without bbar: every onset target chosen < 2 times in stage 1)
##   (3b) joint model with b = 0 (the mean target) for targets without a repeat-target FE, so every stage-2 dispute has a bbar
##   plus a census of the rows the t-1 timing rule and the bbar rule remove (main_table_missingness.csv).
## Usage: Rscript R/main_table.R --build=DIR --dump=FILE --positions=FILE --spec=ID --edges_dir=DIR --cell_dir=DIR --out=DIR [--cores=2] [--n_robust=500]
suppressPackageStartupMessages({ library(data.table); library(parallel) })
for (f in c("s2.R", "s2_fit.R", "sry.R", "s2_cell.R")) source(file.path("R/lib", f))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
BUILD <- opt("build"); DUMPF <- opt("dump"); POSF <- opt("positions"); SPEC <- opt("spec"); EDIR <- opt("edges_dir"); CELL <- opt("cell_dir"); OUT <- opt("out")
CORES <- as.integer(opt("cores", "2")); N_ROBUST <- as.integer(opt("n_robust", "500")); stopifnot(!is.null(BUILD), !is.null(DUMPF), !is.null(POSF), !is.null(SPEC), !is.null(EDIR), !is.null(CELL), !is.null(OUT))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); set.seed(123)
B <- s2_base(BUILD); W_REF <- B$window; dump <- readRDS(DUMPF)
pz <- readRDS(POSF); key <- as.data.table(pz$key); stopifnot(SPEC %in% colnames(pz$R)); R <- pz$R[, SPEC]
E <- s2_spec_edges("core", SPEC, EDIR); X <- s2_prepare(B, key, R, E)
idv <- readLines(file.path(CELL, "cell_identity.txt")); stopifnot(grepl(SPEC, idv, fixed = TRUE))
idp <- trimws(strsplit(idv, " | ", fixed = TRUE)[[1]]); stopifnot(length(idp) == 6L, idp[2] == SPEC); NB <- as.integer(idp[5])   # file | spec | window | sample | n_boot | block (spec ids contain "|")
P <- fread(file.path(CELL, "point.csv"))
DRA <- unlist(lapply(sort(list.files(CELL, pattern = "^draws_.*\\.rds$", full.names = TRUE)), function(f) readRDS(f)$draws), recursive = FALSE)
stopifnot(length(DRA) == NB); okd <- vapply(DRA, function(x) is.numeric(x) && identical(names(x), boot_names(dump)), TRUE); MJ <- do.call(rbind, DRA[okd])
## stage-1 data and the joint model's disputes (finite bbar), from the same point fit as s2_07
D <- s1_add_tie(s1_prepare(dump, key, R), E, dump$roles); ref_order <- names(sort(table(D[chosen == 1L, c2]), decreasing = TRUE))
PJ <- joint_fe_cell(D, dump, X, B$TS); JS <- attr(PJ, "disp_joint"); stopifnot(length(JS) > 0, all(JS %in% X$dispnum))
stopifnot(isTRUE(all.equal(PJ[stage == "S2" & axis == "theta" & alt == "A", est], P[stage == "S2" & axis == "theta" & alt == "A", est], tolerance = 1e-6)))   # same point fit as the cell
XJ <- X[dispnum %in% JS]; stopifnot(nrow(XJ) == P$n_s2_rows[1])
tie0 <- c("A:tie_sum_A", "B:tie_sum_B", "O:tie_sum_A")
s2a <- function(Xs, w = NULL, tag = "S2A", vc = TRUE) { f <- fit_mnl(Xs, S2_REF_DESIGN, tag, vc = vc, w = w); s <- surf_contrasts(f)
  setNames(c(s$est, f$b[tie0]), c(paste(tag, s$alt, s$axis, s$gap, sep = "."), paste0(tag, ".", tie0))) }
nm0 <- names(s2a(XJ)); pt0 <- s2a(XJ)                       # (2) on the joint sample
nmb <- names(s2a(X, tag = "S2Ab")); ptb <- s2a(X, tag = "S2Ab")   # (2b) on all stage-2 disputes
disp_all <- sort(unique(c(as.data.table(dump$D)$disp, B$X$dispnum)))
MM <- boot_weights(disp_all, NB)   # Bayesian bootstrap, the same draws as s2_07 (lib/s2_cell.R)
D0 <- mclapply(MM, function(m) tryCatch(s2a(XJ, as.numeric(m[as.character(XJ$dispnum)]), vc = FALSE)[nm0], error = function(e) NULL), mc.cores = CORES, mc.preschedule = FALSE)
M0 <- do.call(rbind, D0[!vapply(D0, is.null, TRUE)])
## robustness draws: the first N_ROBUST draws (boot_weights is prefix-stable, so these are draws 1..N_ROBUST of the main run)
NR <- min(N_ROBUST, NB); MR <- MM[seq_len(NR)]
Db <- mclapply(MR, function(m) tryCatch(s2a(X, as.numeric(m[as.character(X$dispnum)]), tag = "S2Ab", vc = FALSE)[nmb], error = function(e) NULL), mc.cores = CORES, mc.preschedule = FALSE)
MB <- do.call(rbind, Db[!vapply(Db, is.null, TRUE)])
PF <- joint_fe_cell(D, dump, X, B$TS, fill_b = TRUE); stopifnot(length(attr(PF, "disp_joint")) == uniqueN(X$dispnum))
DF <- mclapply(MR, function(m) tryCatch(joint_fe_draw(D, dump, X, B$TS, m, ref_order, fill_b = TRUE), error = function(e) NULL), mc.cores = CORES, mc.preschedule = FALSE)
MF <- do.call(rbind, DF[vapply(DF, function(x) is.numeric(x) && identical(names(x), boot_names(dump)), TRUE)])
## assemble
ci <- function(M, k, est, nb) { x <- if (is.null(M)) numeric(0) else M[is.finite(M[, k]), k]; if (!is.finite(est) || length(x) < 0.8 * nb) c(NA, NA) else unname(quantile(x, c(.025, .975))) }
row <- function(col, label, q, est, M, nb = NB) { b <- ci(M, q, est, nb); data.table(column = col, row = label, est = est, lo = b[1], hi = b[2], n_draws = if (is.null(M)) 0L else sum(is.finite(M[, q]))) }
getP <- function(st, al, ax, g) P[stage == st & alt == al & axis == ax & (is.na(g) & is.na(gap) | gap == g), est]
T <- list()
lab_s1 <- c("M.1" = "weaker target (material gap +1)", "M.-1" = "stronger target (material gap -1)", "R.1" = "less central target (relational gap +1)", "R.-1" = "more central target (relational gap -1)")
for (r in dump$roles) for (k in names(lab_s1)) { ax <- sub("\\..*$", "", k); g <- as.numeric(sub("^[MR]\\.", "", k))
  T[[length(T) + 1]] <- row("(1) Stage 1, between targets", paste(r, lab_s1[[k]]), paste("S1NF", r, ax, g, sep = "."), getP("S1NF", r, ax, g), MJ)
  T[[length(T) + 1]] <- row("(3) Joint: stage 1 within targets", paste(r, lab_s1[[k]]), paste("S1", r, ax, g, sep = "."), getP("S1", r, ax, g), MJ) }
## stage-2 gaps are candidate minus the side the candidate would fight (A and O fight the target side B; B fights the attacking side A), so +1 = that side is weaker / less central (audit 2026-11-24)
lab_s2 <- c("M.1" = "%s weaker (material gap +1)", "M.-1" = "%s stronger (material gap -1)", "R.1" = "%s less central (relational gap +1)", "R.-1" = "%s more central (relational gap -1)")
opp <- c(A = "target side", B = "attacking side", O = "target side")
for (al in c("A", "B", "O")) for (k in names(lab_s2)) { ax <- sub("\\..*$", "", k); g <- as.numeric(sub("^[MR]\\.", "", k)); q0 <- paste("S2A", al, ax, g, sep = ".")
  lb <- paste("join", al, sprintf(lab_s2[[k]], opp[[al]]))
  T[[length(T) + 1]] <- row("(2) Stage 2 alone", lb, q0, pt0[[q0]], M0)
  T[[length(T) + 1]] <- row("(3) Joint: stage 2", lb, paste("S2", al, ax, g, sep = "."), getP("S2", al, ax, g), MJ) }
for (i in 1:3) { al <- c("A", "B", "O")[i]
  T[[length(T) + 1]] <- row("(2) Stage 2 alone", paste("join", al, "recruitment ties"), paste0("S2A.", tie0[i]), pt0[[paste0("S2A.", tie0[i])]], M0)
  T[[length(T) + 1]] <- row("(3) Joint: stage 2", paste("join", al, "recruitment ties"), tie0[i], getP("S2", al, "tie_sum", NA), MJ)
  T[[length(T) + 1]] <- row("(3) Joint: stage 2", paste("join", al, "theta (target trait)"), paste0(al, ":bbar"), getP("S2", al, "theta", NA), MJ) }
TAB <- rbindlist(T)
fmt <- function(TB) TB[, cell := ifelse(is.na(lo), sprintf("%.2f", est), sprintf("%.2f [%.2f, %.2f]%s", est, lo, hi, ifelse(lo > 0 | hi < 0, "*", "")))]
fmt(TAB); fwrite(TAB, file.path(OUT, "main_table_long.csv"))
N <- data.table(quantity = c("stage-1 events", "stage-2 disputes, (2) and (3)", "stage-2 candidate rows, (2) and (3)", "joined A", "joined B", "joined O",
                             "stage-2 disputes, all (2b, 3b)", "stage-2 candidate rows, all", "joined A, all", "joined B, all", "joined O, all",
                             "bootstrap draws, joint usable", "bootstrap draws, stage 2 alone usable", "robustness draws attempted", "robustness draws usable (2b / 3b)"),
                value = c(P$n_s1_events[1], length(JS), nrow(XJ), sum(XJ$y == "A"), sum(XJ$y == "B"), sum(XJ$y == "orig"),
                          uniqueN(X$dispnum), nrow(X), sum(X$y == "A"), sum(X$y == "B"), sum(X$y == "orig"), nrow(MJ), nrow(M0), NR,
                          paste(if (is.null(MB)) 0L else nrow(MB), if (is.null(MF)) 0L else nrow(MF), sep = " / ")))
fwrite(N, file.path(OUT, "main_table_n.csv"))
W <- dcast(TAB, row ~ column, value.var = "cell", fill = ""); W <- W[match(unique(TAB$row), row)]   # keep the intended row order
fwrite(W, file.path(OUT, "main_table_wide.csv"))
## ---- missing-data robustness: (2) vs (2b) and (3) vs (3b) side by side ----
gF <- function(st, al, ax, g) PF[stage == st & alt == al & axis == ax & (is.na(g) & is.na(gap) | gap == g), est]
RB <- list()
for (al in c("A", "B", "O")) for (k in names(lab_s2)) { ax <- sub("\\..*$", "", k); g <- as.numeric(sub("^[MR]\\.", "", k)); lb <- paste("join", al, sprintf(lab_s2[[k]], opp[[al]]))
  RB[[length(RB) + 1]] <- row("(2) Stage 2 alone, joint sample", lb, paste("S2A", al, ax, g, sep = "."), pt0[[paste("S2A", al, ax, g, sep = ".")]], M0)
  RB[[length(RB) + 1]] <- row("(2b) Stage 2 alone, all disputes", lb, paste("S2Ab", al, ax, g, sep = "."), ptb[[paste("S2Ab", al, ax, g, sep = ".")]], MB, NR)
  RB[[length(RB) + 1]] <- row("(3) Joint", lb, paste("S2", al, ax, g, sep = "."), getP("S2", al, ax, g), MJ)
  RB[[length(RB) + 1]] <- row("(3b) Joint, b = 0 for one-time targets", lb, paste("S2", al, ax, g, sep = "."), gF("S2", al, ax, g), MF, NR) }
for (al in c("A", "B", "O")) {
  RB[[length(RB) + 1]] <- row("(3) Joint", paste("join", al, "theta (target trait)"), paste0(al, ":bbar"), getP("S2", al, "theta", NA), MJ)
  RB[[length(RB) + 1]] <- row("(3b) Joint, b = 0 for one-time targets", paste("join", al, "theta (target trait)"), paste0(al, ":bbar"), gF("S2", al, "theta", NA), MF, NR) }
RT <- rbindlist(RB); fmt(RT); fwrite(RT, file.path(OUT, "main_table_robust_long.csv"))
RW <- dcast(RT, row ~ column, value.var = "cell", fill = ""); RW <- RW[match(unique(RT$row), row)]; fwrite(RW, file.path(OUT, "main_table_robust_wide.csv"))
## census of what the t-1 timing rule (build reference cohort) and the bbar rule remove
O1 <- as.data.table(readRDS(file.path(BUILD, sprintf("s2_oneshot_w%d.rds", W_REF))))
MC <- rbind(O1[excl_reason != "", .(rows = .N, disputes_lost_entirely = NA_integer_, joins = sum(y != "out")), by = .(rule = paste0("build: ", excl_reason))],
            data.table(rule = "build: disputes with no reference-cohort row", rows = 0L, disputes_lost_entirely = uniqueN(O1$dispnum) - uniqueN(O1[in_ref == TRUE]$dispnum), joins = NA_integer_),
            data.table(rule = "joint: no repeat target (bbar undefined)", rows = nrow(X) - nrow(XJ), disputes_lost_entirely = uniqueN(X$dispnum) - length(JS), joins = sum(X[!dispnum %in% JS]$y != "out")))
fwrite(MC, file.path(OUT, "main_table_missingness.csv"))
cat(sprintf("MAIN_TABLE_DONE %s | joint draws %d, stage-2-alone draws %d of %d | robustness %d draws (2b %d, 3b %d usable)\n", OUT, nrow(MJ), nrow(M0), NB, NR,
            if (is.null(MB)) 0L else nrow(MB), if (is.null(MF)) 0L else nrow(MF)))
