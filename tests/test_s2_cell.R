## tests/test_s2_cell.R -- acceptance for the stage-2 parallel-multiverse cell (R/lib/s2_cell.R) on real builds:
## (1) the main measure's edges rebuilt from the per-layer edge file equal the builder's edges; (2) the cell reproduces R/s2_03_joint_fe.R's
## point estimates for the main measure. Usage: Rscript tests/test_s2_cell.R BUILD_DIR EDGES_DIR DUMP POSITIONS JOINT_FE_DIR
suppressPackageStartupMessages(library(data.table))
source("R/lib/s2.R"); source("R/lib/s2_fit.R"); source("R/lib/sry.R"); source("R/lib/s2_cell.R")
a <- commandArgs(TRUE); ok_all <- TRUE
check <- function(name, cond) { cat(sprintf("%-92s %s\n", name, if (isTRUE(cond)) "PASS" else "FAIL")); if (!isTRUE(cond)) ok_all <<- FALSE }
B <- s2_base(a[1]); W_REF <- B$window
MAIN <- "W2contig_ptaall|igo+def+off+neu+nonagg+consul+dca+un|pr"
E <- s2_spec_edges("core", MAIN, a[2])
cmp <- merge(E, B$main_edges, by = c("a", "b", "year"), all = TRUE)
check("1 main-measure edges from the per-layer edge file equal the builder's edges (same dyad-years, max |diff| < 1e-12)", !anyNA(cmp) && max(abs(cmp$w.x - cmp$w.y)) < 1e-12)
pz <- readRDS(a[4]); col <- if (MAIN %in% colnames(pz$R)) MAIN else "W2contig|noPTA|pr"
key <- as.data.table(pz$key); R <- pz$R[, col]
dump <- readRDS(a[3]); D <- s1_prepare(dump, key, R)
X <- s2_prepare(B, key, R, E)
check("2 stage-2 inputs rebuilt from the main measure reproduce the build's ties and gaps", { o <- merge(X, B$X[, .(dispnum, ccode, t0 = tie_sum_A, g0 = dR_kB)], by = c("dispnum", "ccode"))
  max(abs(o$tie_sum_A - o$t0)) < 1e-10 && max(abs(o$dR_kB - o$g0)) < 1e-10 && nrow(X) == nrow(B$X) })
res <- joint_fe_cell(D, dump, X, B$TS)
J <- readRDS(file.path(a[5], "joint_fe.rds"))
th_ref <- J$F_t$b[c("A:bbar", "B:bbar", "O:bbar")]
check("3 cell theta equals R/s2_03_joint_fe.R (main measure)", isTRUE(all.equal(unname(res[axis == "theta", est]), unname(th_ref), tolerance = 1e-5)))
s1r <- J$S1C[model == "target FE"]
check("3 cell stage-1 FE contrasts equal R/s2_03_joint_fe.R", isTRUE(all.equal(res[stage == "S1", est], s1r$est, tolerance = 1e-5)))
bf <- file.path(a[5], "boot_draws_raw.rds")
if (file.exists(bf)) { DR0 <- readRDS(bf); disp_all <- sort(unique(c(as.data.table(dump$D)$disp, B$X$dispnum)))
  set.seed(123); MM <- lapply(1:2, function(b) { m <- tabulate(sample.int(length(disp_all), length(disp_all), replace = TRUE), length(disp_all)); names(m) <- disp_all; m })
  ro <- names(sort(table(D[chosen == 1L, c2]), decreasing = TRUE))
  dd <- lapply(MM, function(m) joint_fe_draw(D, dump, X, B$TS, m, ro))
  check("4 bootstrap draws 1-2 (weights, common seed) equal R/s2_03_joint_fe.R's draws (row replication)",
        all(vapply(1:2, function(i) isTRUE(all.equal(unname(dd[[i]]), unname(DR0[[i]][names(dd[[i]])]), tolerance = 1e-5)), TRUE))) }
print(res[, .(stage, alt, axis, gap, est = round(est, 3), se = round(se, 3))])
cat(if (ok_all) "ALL PASS\n" else "SOME TESTS FAILED\n"); quit(status = if (ok_all) 0 else 1)
