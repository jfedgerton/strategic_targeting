## R/lib/s2_cell.R -- one cell of the stage-2 parallel multiverse (PI 2026-11-22): the joint two-stage FE model for a given relational measure,
## window and sample, with the measure's OWN positions (stage-1 dR, stage-2 surface) and OWN edges (stage-2 ties). Shared by R/s2_05_sweep.R
## and its tests; mirrors R/s2_03_joint_fe.R exactly for the main measure (checked by tests/test_s2_cell.R).
## Requires R/lib/s2.R, R/lib/s2_fit.R, R/lib/sry.R and W_REF set by the caller.
suppressPackageStartupMessages(library(survival))

## ---- edges of any measure from the per-layer edge files written by R/02_build_positions.R --edges_out ---------------------------------
.edge_cache <- new.env()
s2_spec_edges <- function(arm, spec_id, edges_dir) {
  p <- strsplit(spec_id, "|", fixed = TRUE)[[1]]
  if (startsWith(p[1], "G")) {   # grid measures (R/lib/grid.R, PI 2026-11-24): edges rebuilt from grid_layers.rds by the same function as the positions
    f <- file.path(edges_dir, "grid_layers.rds"); if (is.null(.edge_cache[[f]])) { if (!file.exists(f)) stop("missing edge file ", basename(f)); g <- readRDS(f); setDT(g); setorder(g, year, ccode1, ccode2); .edge_cache[[f]] <- g }
    if (!exists("grid_edges", mode = "function")) source("R/lib/grid.R")
    return(grid_edges(.edge_cache[[f]], spec_id))
  }
  if (startsWith(p[1], "BOOKc12_") || startsWith(p[1], "BOOKraw_")) {   # book measures (R/book_measure_c12.R weighted; R/book_measure_raw.R unweighted = main, PI 2026-11-24): edges stored per spec
    f <- file.path(edges_dir, if (startsWith(p[1], "BOOKraw_")) "edges_book_raw.rds" else "edges_book_c12.rds"); if (is.null(.edge_cache[[f]])) { if (!file.exists(f)) stop("missing edge file ", basename(f)); .edge_cache[[f]] <- readRDS(f) }
    E <- .edge_cache[[f]]$edges[[spec_id]]; if (is.null(E)) stop("no book-measure edges for ", spec_id)
    stopifnot(all(E$w > 0 & E$w <= 1 + 1e-12)); return(E[, .(a, b, year, w = pmin(w, 1))])
  }
  if (arm == "igo") { f <- file.path(edges_dir, "edges_igo_arm.rds"); key <- f
    if (is.null(.edge_cache[[key]])) .edge_cache[[key]] <- readRDS(f); e <- .edge_cache[[key]]
    E <- e$E[, p[2]]; dy <- e$dy
    if (p[2] %in% c("raw_count", "bip_newman")) { yr <- split(seq_len(nrow(dy)), dy$year); En <- numeric(length(E)); for (ix in yr) { m <- max(E[ix]); if (m > 0) En[ix] <- E[ix] / m }; E <- En }   # raw_count and bip_newman (sum of 1/(k_g - 1), can exceed 1) are scaled by the year maximum for ties; cosine is already in [0,1]
  } else {
    wt <- sub("_pta.*$", "", p[1]); tag <- sub("^.*_(pta[a-z]+)$", "\\1", p[1]); L <- strsplit(p[2], "+", fixed = TRUE)[[1]]
    f <- file.path(edges_dir, sprintf("edges_%s_%s_%s.rds", arm, wt, tag)); key <- f
    if (is.null(.edge_cache[[key]])) { if (!file.exists(f)) stop("missing edge file ", basename(f)); .edge_cache[[key]] <- readRDS(f) }
    e <- .edge_cache[[key]]; stopifnot(all(L %in% colnames(e$EL))); dy <- e$dy
    yr <- split(seq_len(nrow(dy)), dy$year)
    if (wt %in% c("W6srybook", "W7sryyear", "W6srybookmm", "W7sryyearmm")) E <- as.numeric(sry_combo(e$EL[, L, drop = FALSE], yr, edge = if (grepl("mm$", wt)) "minmax" else "clip"))
    else if (wt == "W1raw") { S <- rowSums(e$EL[, L, drop = FALSE]); E <- numeric(length(S)); for (ix in yr) { m <- max(S[ix]); if (m > 0) E[ix] <- S[ix] / m } }
    else { S <- sqrt(rowSums(e$EL[, L, drop = FALSE]^2)); E <- numeric(length(S)); for (ix in yr) { r <- range(S[ix]); if (r[2] - r[1] >= 1e-12) E[ix] <- (S[ix] - r[1]) / (r[2] - r[1]) } }
  }
  stopifnot(all(E >= -1e-12 & E <= 1 + 1e-12))
  data.table(a = dy$ccode1, b = dy$ccode2, year = dy$year, w = pmin(pmax(E, 0), 1))[w > 0]
}

## ---- stage 1 (as R/07 and R/s2_03): dR from the measure's positions at t-1, scaled by the gap scale of year t; quadratic gap terms by type
s1_prepare <- function(dump, key, R) {
  D <- as.data.table(dump$D); ROLES <- dump$roles; kk <- paste(key$ccode, key$year)
  sc <- tapply(R, key$year, s2_gap_scale)
  i <- match(paste(D$c1, D$yr - 1L), kk); j <- match(paste(D$c2, D$yr - 1L), kk)
  if (anyNA(i) || anyNA(j)) stop("stage 1: positions missing for event rows")
  D[, dR := (R[i] - R[j]) / sc[as.character(yr)]]
  ## constant positions (GPT P2-5): the R/07 rule, sd < 1e-8 in any year the fit uses (positions at t-1, scale of year t)
  uy <- sort(unique(c(D$yr, D$yr - 1L))); uy <- uy[uy >= min(key$year)]; sds <- vapply(split(R, key$year)[as.character(uy)], sd, 0)
  if (any(!(sds >= 1e-8))) stop(sprintf("position constant in stage-1 year(s) %s", paste(uy[!(sds >= 1e-8)], collapse = " ")))
  if (any(!is.finite(D$dR))) stop("stage 1: non-finite dR")
  D[, `:=`(dM2 = dM^2, dR2 = dR^2, dMR = dM * dR)]
  for (g in 0:1) for (v in c("dM", "dM2", "dR", "dR2", "dMR")) D[, paste0(v, "_", ROLES[g + 1]) := get(v) * (type == g)]
  D[]
}
## robustness term (PI 2026-11-23): the initiator-target composite tie from the measure's OWN edges at t-1 (the year of the positions),
## by initiator type (tie_bil, tie_coal); a pair with no positive edge has tie 0. Used only by the S1NFT variant (between targets + tie).
s1_add_tie <- function(D, EDGE, roles) {
  D <- copy(D); E <- EDGE[, .(k = paste(pmin(a, b), pmax(a, b), year), w)]; stopifnot(!anyDuplicated(E$k))
  D[, tie_it := E$w[match(paste(pmin(c1, c2), pmax(c1, c2), yr - 1L), E$k)]][is.na(tie_it), tie_it := 0]
  for (g in 0:1) D[, paste0("tie_", roles[g + 1]) := tie_it * (type == g)]
  D[]
}
s1_terms <- function(dump) as.vector(outer(c("dM", "dM2", "dR", "dR2", "dMR"), dump$roles, paste, sep = "_"))
s1_fit_fe <- function(D, dump, ref_order) {   # target FE on the rows of chosen targets only (exact limit; never-chosen FE -> -Inf)
  tgb <- unique(D[chosen == 1L, c2]); dd <- D[c2 %in% tgb]; dd <- dd[dd[, .(k = .N), by = ev_id][k > 1L], on = "ev_id"]
  dd[, tfe := droplevels(factor(as.character(c2), levels = intersect(ref_order, as.character(tgb))))]
  f <- as.formula(paste("Surv(rep(1, nrow(dd)), chosen) ~", paste(c(s1_terms(dump), dump$Z_cols, "tfe"), collapse = " + "), "+ strata(ev_id)"))
  ## validity (GPT P1-1), the R/07 policy: a convergence / infinite / singular warning, a non-finite coefficient or the iteration limit makes
  ## the fit INVALID; callers stop with "stage-1 fit invalid" (cell status INVALID; a bootstrap draw counts as a numerical failure)
  w_ <- character(0)
  m <- withCallingHandlers(coxph(f, data = dd, weights = w, method = "breslow", control = coxph.control(iter.max = 50)),
                           warning = function(e) { w_ <<- c(w_, conditionMessage(e)); invokeRestart("muffleWarning") })
  bad <- c(if (any(grepl("converge|infinite|singular", w_, ignore.case = TRUE))) paste(unique(sub(" variable.*$", "", w_)), collapse = " | "),
           if (!all(is.finite(coef(m)))) "non-finite coefficients", if (m$iter >= 50L) "iteration limit")
  if (length(bad)) stop("stage-1 fit invalid: ", paste(bad, collapse = "; "))
  list(m = m, dd = dd, warnings = w_)
}
## stage 1 WITHOUT target FE (PI 2026-11-23: the between-target comparison is the main stage-1 question; the target-FE fit is the
## within-target appendix test and supplies bbar). Same validity policy as s1_fit_fe. Inference: three-way clustered covariance over
## initiator, dispute and target (inclusion-exclusion; repeated targets make target clustering the honest between-target choice)
s1_fit_nofe <- function(D, dump, extra = character(0)) {
  f <- as.formula(paste("Surv(rep(1, nrow(D)), chosen) ~", paste(c(s1_terms(dump), extra, dump$Z_cols), collapse = " + "), "+ strata(ev_id)"))
  w_ <- character(0)
  m <- withCallingHandlers(coxph(f, data = D, weights = w, method = "breslow", control = coxph.control(iter.max = 50)),
                           warning = function(e) { w_ <<- c(w_, conditionMessage(e)); invokeRestart("muffleWarning") })
  bad <- c(if (any(grepl("converge|infinite|singular", w_, ignore.case = TRUE))) paste(unique(sub(" variable.*$", "", w_)), collapse = " | "),
           if (!all(is.finite(coef(m)))) "non-finite coefficients", if (m$iter >= 50L) "iteration limit")
  if (length(bad)) stop("stage-1 fit invalid (no target FE): ", paste(bad, collapse = "; "))
  m
}
s1_vc3w <- function(m, dd) { cl <- function(g) { Dd <- rowsum(residuals(m, type = "dfbeta", weighted = TRUE), g, reorder = FALSE); crossprod(Dd) }
  V <- cl(dd$init) + cl(dd$disp) + cl(dd$c2) - cl(paste(dd$init, dd$disp)) - cl(paste(dd$init, dd$c2)) - cl(paste(dd$disp, dd$c2)) + cl(paste(dd$init, dd$disp, dd$c2))
  e <- eigen(V, symmetric = TRUE); if (any(e$values < 0)) V <- e$vectors %*% diag(pmax(e$values, 0)) %*% t(e$vectors)
  dimnames(V) <- list(names(coef(m)), names(coef(m))); V }
s1_vc2w <- function(m, dd) { cl <- function(g) { Dd <- rowsum(residuals(m, type = "dfbeta", weighted = TRUE), g, reorder = FALSE); crossprod(Dd) }
  V <- cl(dd$init) + cl(dd$disp) - cl(dd$init_disp); e <- eigen(V, symmetric = TRUE)
  if (any(e$values < 0)) V <- e$vectors %*% diag(pmax(e$values, 0)) %*% t(e$vectors); dimnames(V) <- list(names(coef(m)), names(coef(m))); V }
s1_contrasts <- function(b, V, roles) { out <- list()
  for (r in roles) for (ax in c("M", "R")) for (sg in c(1, -1)) { a <- setNames(numeric(length(b)), names(b)); a[paste0("d", ax, "_", r)] <- sg; a[paste0("d", ax, "2_", r)] <- 1
    out[[length(out) + 1L]] <- data.table(stage = "S1", alt = r, axis = ax, gap = sg, est = sum(a * b), se = if (is.null(V)) NA_real_ else sqrt(max(drop(t(a) %*% V %*% a), 0))) }
  rbindlist(out) }

## ---- stage 2 inputs for a measure: gaps from its positions, ties from its edges (one-shot, month-0 rosters, onset year - 1) ----------
s2_prepare <- function(B, key, R, EDGE) {
  PV <- merge(as.data.table(B$SI$PV)[, .(ccode, year, M)], data.table(ccode = key$ccode, year = key$year, R = R), by = c("ccode", "year"), all = TRUE)
  ## common cohort (GPT P2-7): B$X is the main build's eligible cohort. Position-dependent eligibility is identical for this measure only if it
  ## has a finite position on exactly the same state-years as the main measure; stop otherwise (never a silently different cohort)
  pk <- as.data.table(B$SI$PV)[is.finite(R), paste(ccode, year)]; mk <- paste(key$ccode, key$year)[is.finite(R)]
  if (!setequal(pk, mk)) stop(sprintf("measure covers different state-years than the main measure (%d vs %d); the main-defined cohort would not apply", length(mk), length(pk)))
  ## constant positions (GPT P2-5): sd < 1e-8 in a year that eligible rows use (positions at ty, scale of year t); unused years are dropped
  ER <- merge(as.data.table(B$SI$rows), unique(B$X[, .(dispnum, ccode)]), by = c("dispnum", "ccode"))
  uy <- intersect(as.character(unique(c(ER$ty, ER$t))), unique(as.character(key$year)))
  sds <- vapply(split(R, key$year)[uy], sd, 0); if (any(!(sds >= 1e-8))) stop("position constant in stage-2 year(s) ", paste(uy[!(sds >= 1e-8)], collapse = " "))
  sR <- vapply(split(R, key$year), function(v) if (sd(v) >= 1e-8) s2_gap_scale(v) else NA_real_, 0); sR <- sR[is.finite(sR)]
  sv <- s2_surface(ER, B$SI$mem_roster, PV, B$SI$sM, sR, "mean")   # eligible rows only
  ## the measure's OWN centrality controls (audit P1, 2026-11-22): own_cent = R of the candidate at onset - 1; tgt_cent_mean / _min over the
  ## baseline target set at onset - 1 (NA if any target lacks a position), exactly as R/s2_01_build_riskset.R builds them from the main measure
  RK <- data.table(ccode = key$ccode, ty = key$year, Rv = R)
  own <- merge(unique(B$row0[, .(dispnum, ccode, ty = onset_year - 1L)]), RK, by = c("ccode", "ty"), all.x = TRUE)[, .(dispnum, ccode, own_cent = Rv)]
  tg <- merge(merge(B$TS, unique(B$row0[, .(dispnum, ty = onset_year - 1L)]), by = "dispnum"), RK, by = c("ccode", "ty"), all.x = TRUE)
  tc <- tg[, .(tgt_cent_mean = if (anyNA(Rv)) NA_real_ else mean(Rv), tgt_cent_min = if (anyNA(Rv)) NA_real_ else min(Rv)), by = dispnum]
  drop <- intersect(c("dM_kA", "dR_kA", "dM_kB", "dR_kB", "tie_sum_A", "tie_sum_B", "tie_max_A", "tie_max_B", "own_cent", "tgt_cent_mean", "tgt_cent_min"), names(B$X))
  X <- merge(B$X[, !drop, with = FALSE], sv[, .(dispnum, ccode, dM_kA, dR_kA, dM_kB, dR_kB)], by = c("dispnum", "ccode"))
  X <- merge(merge(X, own, by = c("dispnum", "ccode"), all.x = TRUE), tc, by = "dispnum", all.x = TRUE)
  T1 <- s2_ties(B$row0[, .(dispnum, ccode, m = 0L, year_tie = onset_year - 1L)], B$R0, EDGE, B$cover)
  X <- merge(X, T1[, .(dispnum, ccode, tie_sum_A, tie_sum_B, tie_max_A, tie_max_B)], by = c("dispnum", "ccode"))
  n0 <- nrow(X); X <- X[is.finite(dM_kA) & is.finite(dR_kA) & is.finite(dM_kB) & is.finite(dR_kB) & is.finite(tie_sum_A) & is.finite(tie_sum_B) & is.finite(own_cent) & is.finite(tgt_cent_mean)]
  add_surface(X); attr(X, "dropped") <- n0 - nrow(X); X
}

## ---- the joint two-stage FE model on one cell ------------------------------------------------------------------------------------------
## returns a long table: stage-1 FE contrasts (two-way SE), stage-2 surface contrasts and theta (two-way SE conditional on bbar), counts
joint_fe_cell <- function(D, dump, X, TS, design = S2_REF_DESIGN, fill_b = FALSE) {   # fill_b: robustness, targets without a repeat-target FE get b = 0 (the mean)
  ch <- D[chosen == 1L]; NCH0 <- sort(table(ch$c2), decreasing = TRUE); REP <- names(NCH0)[NCH0 >= 2]
  s1 <- s1_fit_fe(D, dump, names(NCH0)); b1 <- coef(s1$m); keep <- !grepl("^tfe", names(b1)); V1 <- s1_vc2w(s1$m, s1$dd)
  bfe <- b1[!keep]; names(bfe) <- sub("^tfe", "", names(bfe)); bfe <- c(setNames(0, levels(s1$dd$tfe)[1]), bfe)   # the reference level actually used by the fit
  b <- bfe[names(bfe) %in% REP & is.finite(bfe)]; b <- pmin(pmax(b - mean(b), -5), 5)
  if (fill_b) { nt <- setdiff(as.character(unique(TS$ccode)), names(b)); b <- c(b, setNames(rep(0, length(nt)), nt)) }   # robustness (PI 2026-11-24): one-time / unchosen targets at the mean
  bb <- TS[, .(dispnum, b = b[as.character(ccode)])][, .(bbar = mean(b, na.rm = TRUE)), by = dispnum][is.finite(bbar)]
  X1 <- merge(X, bb, by = "dispnum")
  F_t <- fit_mnl(X1, lapply(design, function(v) c(v, "bbar")), "cell")
  sc <- surf_contrasts(F_t)[, .(stage = "S2", alt, axis, gap, est, se)]
  th <- data.table(stage = "S2", alt = c("A", "B", "O"), axis = "theta", gap = NA_real_, est = F_t$b[c("A:bbar", "B:bbar", "O:bbar")], se = sqrt(diag(F_t$V))[c("A:bbar", "B:bbar", "O:bbar")])
  tie <- data.table(stage = "S2", alt = c("A", "B", "O"), axis = "tie_sum", gap = NA_real_, est = F_t$b[c("A:tie_sum_A", "B:tie_sum_B", "O:tie_sum_A")],
                    se = sqrt(diag(F_t$V))[c("A:tie_sum_A", "B:tie_sum_B", "O:tie_sum_A")])
  mN <- s1_fit_nofe(D, dump); kN <- names(coef(mN)) %in% s1_terms(dump)
  nf <- s1_contrasts(coef(mN)[kN], s1_vc3w(mN, D)[kN, kN], dump$roles)[, stage := "S1NF"]   # between targets (main stage-1 comparison)
  TT <- paste0("tie_", dump$roles); nft <- NULL
  if (all(TT %in% names(D))) {   # robustness: between targets + the initiator-target composite tie
    mT <- s1_fit_nofe(D, dump, TT); VT <- s1_vc3w(mT, D); kT <- names(coef(mT)) %in% s1_terms(dump)
    nft <- rbind(s1_contrasts(coef(mT)[kT], VT[kT, kT], dump$roles)[, stage := "S1NFT"],
                 data.table(stage = "S1NFT", alt = dump$roles, axis = "tie", gap = NA_real_, est = coef(mT)[TT], se = sqrt(diag(VT))[TT])) }
  out <- rbind(s1_contrasts(b1[keep], V1[keep, keep], dump$roles), nf, nft, sc, th, tie)
  out[, `:=`(n_s1_events = nrow(ch), n_s2_rows = F_t$n, n_s2_A = sum(X1$y == "A"), n_s2_B = sum(X1$y == "B"), n_s2_O = sum(X1$y == "orig"))]
  setattr(out, "disp_joint", sort(unique(X1$dispnum)))   # the joint model's disputes (finite bbar); main_table fits stage 2 alone on these
  out
}

## ---- one two-stage bootstrap draw of a cell (as R/s2_03_joint_fe.R): dispute weights m (named by dispnum; Bayesian bootstrap, boot_weights) reweight BOTH steps;
## step 1 refits the target-FE model (reference = most-chosen target present in the draw, b = 0), bbar is rebuilt, step 2 refits with row
## weights (equal to row replication, tests/test_s2_fit.R 4). Returns a named vector of the cell's quantities (names as boot_names()).
## ---- dispute weights of the two-stage bootstrap (PI 2026-11-24: Bayesian bootstrap, Rubin 1981): each draw gives every dispute an
## Exp(1) weight scaled to mean 1 (Dirichlet weights x n), seed 123. Every dispute stays in every draw (no zero weights), so no target, joiner
## or sparse covariate cell drops out of a draw. Shared by R/s2_07_boot_cells.R, R/main_table.R and R/marginal_effects.R so all three use
## identical draws. Replaces multinomial resampling (tabulate(sample.int(...))).
boot_weights <- function(disp_all, NB) {
  stopifnot(length(disp_all) > 1L, !anyDuplicated(disp_all), NB >= 1L)
  set.seed(123)
  lapply(seq_len(NB), function(b) { g <- rexp(length(disp_all)); m <- g / mean(g); names(m) <- disp_all; m })
}
boot_names <- function(dump) { r <- dump$roles; ax <- rep(c("M", "R"), each = 2); sg <- c(1, -1)
  c(paste("S1", rep(r, each = 4), ax, sg, sep = "."), paste("S1NF", rep(r, each = 4), ax, sg, sep = "."), paste("S1NFT", rep(r, each = 4), ax, sg, sep = "."),
    paste0("S1NFT.", r, ".tie"), paste("S2", rep(c("A", "B", "O"), each = 4), ax, sg, sep = "."),
    "A:bbar", "B:bbar", "O:bbar", "A:tie_sum_A", "B:tie_sum_B", "O:tie_sum_A") }
joint_fe_draw <- function(D, dump, X, TS, m, ref_order, design = S2_REF_DESIGN, fill_b = FALSE) {
  ch <- D[chosen == 1L]; NCH <- table(ch$c2); REP <- names(NCH)[NCH >= 2]   # repeat targets of the full sample (as R/s2_03)
  Db <- copy(D); Db[, w := w * m[as.character(disp)]]; Db <- Db[w > 0]
  s1 <- s1_fit_fe(Db, dump, ref_order); bb <- coef(s1$m); be <- bb[grepl("^tfe", names(bb))]; names(be) <- sub("^tfe", "", names(be))
  refb <- setdiff(intersect(ref_order, as.character(unique(Db[chosen == 1L, c2]))), names(be))[1]; be <- c(setNames(0, refb), be)
  c1 <- s1_contrasts(bb[!grepl("^tfe", names(bb))], NULL, dump$roles)
  mN <- s1_fit_nofe(Db, dump); kN <- names(coef(mN)) %in% s1_terms(dump); cN <- s1_contrasts(coef(mN)[kN], NULL, dump$roles)
  TT <- paste0("tie_", dump$roles); vT <- NULL
  if (all(TT %in% names(Db))) { mT <- s1_fit_nofe(Db, dump, TT); kT <- names(coef(mT)) %in% s1_terms(dump); cT <- s1_contrasts(coef(mT)[kT], NULL, dump$roles)
    vT <- c(setNames(cT$est, paste("S1NFT", cT$alt, cT$axis, cT$gap, sep = ".")), setNames(coef(mT)[TT], paste0("S1NFT.", dump$roles, ".tie"))) }
  b <- be[names(be) %in% REP & is.finite(be)]; b <- pmin(pmax(b - mean(b), -5), 5)
  if (fill_b) { nt <- setdiff(as.character(unique(TS$ccode)), names(b)); b <- c(b, setNames(rep(0, length(nt)), nt)) }
  bbt <- TS[, .(dispnum, b = b[as.character(ccode)])][, .(bbar = mean(b, na.rm = TRUE)), by = dispnum][is.finite(bbar)]
  Xb <- merge(X, bbt, by = "dispnum"); wt <- m[as.character(Xb$dispnum)]; Xb <- Xb[wt > 0]; wt <- wt[wt > 0]
  ft <- fit_mnl(Xb, lapply(design, function(v) c(v, "bbar")), "boot", vc = FALSE, w = as.numeric(wt))
  if (!"A:bbar" %in% names(ft$b)) stop("bbar dropped by the separation rule in this draw")
  sc <- surf_contrasts(ft)
  v <- c(setNames(c1$est, paste("S1", c1$alt, c1$axis, c1$gap, sep = ".")), setNames(cN$est, paste("S1NF", cN$alt, cN$axis, cN$gap, sep = ".")), vT, setNames(sc$est, paste("S2", sc$alt, sc$axis, sc$gap, sep = ".")),
         ft$b[c("A:bbar", "B:bbar", "O:bbar", "A:tie_sum_A", "B:tie_sum_B", "O:tie_sum_A")])
  bn <- boot_names(dump); setNames(v[bn], bn)   # a dropped term gives NA for that quantity only, never a misnamed draw (recheck P2)
}

## base build of one window x sample (measure-independent parts), from an R/s2_01_build_riskset.R output folder
s2_base <- function(build_dir) {
  META <- readRDS(file.path(build_dir, "s2_meta.rds")); W <- as.integer(META$window)
  ONE <- as.data.table(readRDS(file.path(build_dir, sprintf("s2_oneshot_w%d.rds", W))))
  w0 <- W; CL <- as.data.table(readRDS(file.path(build_dir, "s2_classification_w0_6.rds")))[W == w0]   # (audit P2: `W == W` filtered nothing)
  X <- ONE[in_ref == TRUE]; for (v in names(X)) if (is.logical(X[[v]])) X[, (v) := as.numeric(get(v))]
  ROS <- as.data.table(readRDS(file.path(build_dir, "s2_rosters_carry.rds")))
  E <- as.data.table(readRDS(file.path(build_dir, "s2_main_edges.rds")))
  list(meta = META, window = W, X = X, SI = readRDS(file.path(build_dir, "s2_surface_inputs.rds")), row0 = ONE[, .(dispnum, ccode, onset_year)],
       R0 = ROS[m == 0L], cover = if (!is.null(META$cover)) sort(META$cover) else { warning("build predates META$cover; using main-edge years"); sort(unique(E$year)) }, TS = CL[baseline == TRUE & side_rel == "B", .(dispnum, ccode)], main_edges = E)
}

## ---- provenance (GPT P1-2): an output folder is bound to the contents of its code and inputs. fingerprint.tsv lists (role, file, md5);
## on a rerun every listed role must have the same md5 (new roles may be added, e.g. more positions files), otherwise the run stops.
s2_fingerprint_check <- function(out, files) {   # files: named character vector, names = role, values = paths
  stopifnot(length(files) > 0, !anyDuplicated(names(files)), all(file.exists(files)))
  cur <- data.table(role = names(files), file = unname(files), md5 = unname(tools::md5sum(files)))
  f <- file.path(out, "fingerprint.tsv")
  if (!file.exists(f) && length(setdiff(list.files(out), c("s2_sweep.log", "s2_boot.log")))) stop(out, " holds outputs but no fingerprint (made before provenance checks); use a new output folder")
  if (file.exists(f)) { old <- fread(f, sep = "\t"); m <- merge(old, cur, by = "role", suffixes = c("_old", "_new"))
    bad <- m[md5_old != md5_new]
    if (nrow(bad)) stop(sprintf("%s was produced from different code or inputs (%s); use a new output folder", out, paste(head(bad$role, 5), collapse = ", ")))
    cur <- unique(rbind(old, cur[!role %in% old$role]), by = "role") }
  fwrite(cur, paste0(f, ".tmp"), sep = "\t"); if (!file.rename(paste0(f, ".tmp"), f)) stop("could not write ", f)
  invisible(cur)
}
s2_code_files <- function(rd) c(lib_s2 = file.path(rd, "lib", "s2.R"), lib_s2_fit = file.path(rd, "lib", "s2_fit.R"), lib_sry = file.path(rd, "lib", "sry.R"), lib_s2_cell = file.path(rd, "lib", "s2_cell.R"), lib_grid = file.path(rd, "lib", "grid.R"))   # grid.R rebuilds grid ties at run time (audit 2026-11-24)
s2_build_files <- function(bd, tag) { f <- sort(list.files(bd, pattern = "\\.rds$", full.names = TRUE)); setNames(f, paste0(tag, ":", basename(f))) }
s2_write_atomic <- function(x, f) { fwrite(x, paste0(f, ".tmp")); if (!file.rename(paste0(f, ".tmp"), f)) stop("could not write ", f) }
