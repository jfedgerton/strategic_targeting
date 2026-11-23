## tests/test_s2_cell.R -- acceptance for the stage-2 parallel-multiverse cell (R/lib/s2_cell.R) on real builds:
## (1) the main measure's edges rebuilt from the per-layer edge file equal the builder's edges; (2) the cell reproduces R/s2_03_joint_fe.R's
## point estimates for the main measure. Usage: Rscript tests/test_s2_cell.R BUILD_DIR EDGES_DIR DUMP POSITIONS JOINT_FE_DIR [POSITIONS_WITH_EDGE_SPEC_IDS]
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
## audit P1: a non-main measure must carry its OWN centrality controls and ties, not the main measure's
alt <- setdiff(colnames(pz$R), col)[1]; Xa <- s2_prepare(B, key, pz$R[, alt], E)
oa <- merge(Xa[, .(dispnum, ccode, oc = own_cent, tc = tgt_cent_mean)], B$X[, .(dispnum, ccode, oc0 = own_cent, tc0 = tgt_cent_mean)], by = c("dispnum", "ccode"))
RA <- data.table(ccode = key$ccode, ty = key$year, Rv = pz$R[, alt]); chk <- merge(Xa[, .(dispnum, ccode, oc = own_cent)], unique(B$row0[, .(dispnum, ccode, ty = onset_year - 1L)]), by = c("dispnum", "ccode"))
chk <- merge(chk, RA, by = c("ccode", "ty"))
check(sprintf("5 measure '%s': own_cent equals its own position at onset - 1 and differs from the main measure's", alt), max(abs(chk$oc - chk$Rv)) < 1e-12 && max(abs(oa$oc - oa$oc0)) > 1e-6)
check("5 main measure: rebuilt own_cent and tgt_cent_mean equal the build's", { o <- merge(X[, .(dispnum, ccode, oc = own_cent, tc = tgt_cent_mean)], B$X[, .(dispnum, ccode, oc0 = own_cent, tc0 = tgt_cent_mean)], by = c("dispnum", "ccode")); max(abs(o$oc - o$oc0)) < 1e-12 && max(abs(o$tc - o$tc0)) < 1e-12 })
## GPT P2-8: alternate-measure TIES from the measure's own edges (positions file in a[6] with spec ids matching the edge files)
if (length(a) >= 6) { pa <- readRDS(a[6]); ka <- as.data.table(pa$key); sa <- grep("\\|igo\\|pr$", colnames(pa$R), value = TRUE)[1]
  Ea <- s2_spec_edges("core", sa, a[2]); Xt <- s2_prepare(B, ka, pa$R[, sa], Ea)
  T1 <- s2_ties(B$row0[, .(dispnum, ccode, m = 0L, year_tie = onset_year - 1L)], B$R0, Ea, B$cover)
  o <- merge(Xt[, .(dispnum, ccode, t = tie_sum_A)], T1[, .(dispnum, ccode, t1 = tie_sum_A)], by = c("dispnum", "ccode"))
  o0 <- merge(Xt[, .(dispnum, ccode, t = tie_sum_A)], B$X[, .(dispnum, ccode, t0 = tie_sum_A)], by = c("dispnum", "ccode"))
  check(sprintf("6 measure '%s': ties come from its own edges (equal to s2_ties on them; differ from the main measure's)", sa), max(abs(o$t - o$t1)) < 1e-12 && max(abs(o0$t - o0$t0)) > 1e-6)
  ## every weighting family with an edge file present: edges in [0,1]; W1raw / Euclidean families reach 1 in every covered year
  for (ef in list.files(a[2], pattern = "^edges_core_.*\\.rds$")) { wt <- sub("^edges_core_(.*)_pta[a-z]+\\.rds$", "\\1", ef); tg <- sub("^.*_(pta[a-z]+)\\.rds$", "\\1", ef)
    Ef <- s2_spec_edges("core", sprintf("%s_%s|igo+def|pr", wt, tg), a[2]); mx <- Ef[, max(w), by = year]
    check(sprintf("6 edges %s: weights in (0,1]%s", ef, if (!grepl("sry", wt)) ", year max 1" else ""), all(Ef$w > 0 & Ef$w <= 1) && (grepl("sry", wt) || all(abs(mx$V1 - 1) < 1e-12))) }
  ## constant and coverage guards
  Rc <- pa$R[, sa]; yy <- 1950L; Rc[ka$year == yy] <- 0.5
  check("6 a measure constant in a used year is refused with 'position constant' (stage 1 and stage 2)",
        grepl("^position constant", tryCatch({ s1_prepare(dump, ka, Rc); "" }, error = conditionMessage)) && grepl("^position constant", tryCatch({ s2_prepare(B, ka, Rc, Ea); "" }, error = conditionMessage)))
  Rn <- pa$R[, sa]; Rn[which(ka$year == 1990L)[1]] <- NA
  check("6 a measure covering different state-years than the main measure is refused (main-defined cohort)", grepl("different state-years", tryCatch({ s2_prepare(B, ka, Rn, Ea); "" }, error = conditionMessage)))
} else cat("6 alternate-measure ties, edge families, guards                                              SKIP (no positions file given)\n")
## GPT P1-1: an invalid stage-1 fit is refused, never returned as estimates (20 events cannot identify the target FE)
Dsm <- D[ev_id %in% head(unique(D$ev_id), 20)]
check("7 an invalid stage-1 FE fit stops with 'stage-1 fit invalid'", grepl("stage-1 fit invalid", tryCatch({ s1_fit_fe(Dsm, dump, names(sort(table(Dsm[chosen == 1L, c2]), decreasing = TRUE))); "" }, error = conditionMessage)))
bf <- file.path(a[5], "boot_draws_raw.rds")
if (file.exists(bf)) { DR0 <- readRDS(bf); disp_all <- sort(unique(c(as.data.table(dump$D)$disp, B$X$dispnum)))
  set.seed(123); MM <- lapply(1:2, function(b) { m <- tabulate(sample.int(length(disp_all), length(disp_all), replace = TRUE), length(disp_all)); names(m) <- disp_all; m })
  ro <- names(sort(table(D[chosen == 1L, c2]), decreasing = TRUE))
  dd <- lapply(MM, function(m) joint_fe_draw(D, dump, X, B$TS, m, ro))
  check("4 bootstrap draws 1-2 (weights, common seed) equal R/s2_03_joint_fe.R's draws (row replication)",
        all(vapply(1:2, function(i) { k <- intersect(names(dd[[i]]), names(DR0[[i]])); length(k) >= 26 && isTRUE(all.equal(unname(dd[[i]][k]), unname(DR0[[i]][k]), tolerance = 1e-5)) }, TRUE))) }   # S1NF (added 2026-11-23) is not in the reference draws else cat("4 bootstrap draws vs R/s2_03_joint_fe.R                                                    SKIP (no boot_draws_raw.rds)\n")
print(res[, .(stage, alt, axis, gap, est = round(est, 3), se = round(se, 3))])
cat(if (ok_all) "ALL PASS\n" else "SOME TESTS FAILED\n"); quit(status = if (ok_all) 0 else 1)
