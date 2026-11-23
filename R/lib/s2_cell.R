## R/lib/s2_cell.R -- one cell of the stage-2 parallel multiverse (PI 2026-11-22): the joint two-stage FE model for a given relational measure,
## window and sample, with the measure's OWN positions (stage-1 dR, stage-2 surface) and OWN edges (stage-2 ties). Shared by R/s2_05_sweep.R
## and its tests; mirrors R/s2_03_joint_fe.R exactly for the main measure (checked by tests/test_s2_cell.R).
## Requires R/lib/s2.R, R/lib/s2_fit.R, R/lib/sry.R and W_REF set by the caller.
suppressPackageStartupMessages(library(survival))

## ---- edges of any measure from the per-layer edge files written by R/02_build_positions.R --edges_out ---------------------------------
.edge_cache <- new.env()
s2_spec_edges <- function(arm, spec_id, edges_dir) {
  p <- strsplit(spec_id, "|", fixed = TRUE)[[1]]
  if (arm == "igo") { f <- file.path(edges_dir, "edges_igo_arm.rds"); key <- f
    if (is.null(.edge_cache[[key]])) .edge_cache[[key]] <- readRDS(f); e <- .edge_cache[[key]]
    E <- e$E[, p[2]]; dy <- e$dy
    if (p[2] == "raw_count") { yr <- split(seq_len(nrow(dy)), dy$year); En <- numeric(length(E)); for (ix in yr) { m <- max(E[ix]); if (m > 0) En[ix] <- E[ix] / m }; E <- En }   # the IGO arm's centrality takes raw weights; scale to [0,1] for ties
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
  D[, `:=`(dM2 = dM^2, dR2 = dR^2, dMR = dM * dR)]
  for (g in 0:1) for (v in c("dM", "dM2", "dR", "dR2", "dMR")) D[, paste0(v, "_", ROLES[g + 1]) := get(v) * (type == g)]
  D[]
}
s1_terms <- function(dump) as.vector(outer(c("dM", "dM2", "dR", "dR2", "dMR"), dump$roles, paste, sep = "_"))
s1_fit_fe <- function(D, dump, ref_order) {   # target FE on the rows of chosen targets only (exact limit; never-chosen FE -> -Inf)
  tgb <- unique(D[chosen == 1L, c2]); dd <- D[c2 %in% tgb]; dd <- dd[dd[, .(k = .N), by = ev_id][k > 1L], on = "ev_id"]
  dd[, tfe := droplevels(factor(as.character(c2), levels = intersect(ref_order, as.character(tgb))))]
  f <- as.formula(paste("Surv(rep(1, nrow(dd)), chosen) ~", paste(c(s1_terms(dump), dump$Z_cols, "tfe"), collapse = " + "), "+ strata(ev_id)"))
  list(m = suppressWarnings(coxph(f, data = dd, weights = w, method = "breslow", control = coxph.control(iter.max = 50))), dd = dd)
}
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
  sR <- vapply(split(R, key$year), s2_gap_scale, 0)
  sv <- s2_surface(B$SI$rows, B$SI$mem_roster, PV, B$SI$sM, sR, "mean")
  X <- merge(B$X[, !c("dM_kA", "dR_kA", "dM_kB", "dR_kB", "tie_sum_A", "tie_sum_B", "tie_max_A", "tie_max_B"), with = FALSE], sv[, .(dispnum, ccode, dM_kA, dR_kA, dM_kB, dR_kB)], by = c("dispnum", "ccode"))
  cover <- sort(unique(B$cover))
  T1 <- s2_ties(B$row0[, .(dispnum, ccode, m = 0L, year_tie = onset_year - 1L)], B$R0, EDGE, cover)
  X <- merge(X, T1[, .(dispnum, ccode, tie_sum_A, tie_sum_B, tie_max_A, tie_max_B)], by = c("dispnum", "ccode"))
  n0 <- nrow(X); X <- X[is.finite(dM_kA) & is.finite(dR_kA) & is.finite(dM_kB) & is.finite(dR_kB) & is.finite(tie_sum_A) & is.finite(tie_sum_B)]
  add_surface(X); attr(X, "dropped") <- n0 - nrow(X); X
}

## ---- the joint two-stage FE model on one cell ------------------------------------------------------------------------------------------
## returns a long table: stage-1 FE contrasts (two-way SE), stage-2 surface contrasts and theta (two-way SE conditional on bbar), counts
joint_fe_cell <- function(D, dump, X, TS, design = S2_REF_DESIGN) {
  ch <- D[chosen == 1L]; NCH0 <- sort(table(ch$c2), decreasing = TRUE); REP <- names(NCH0)[NCH0 >= 2]
  s1 <- s1_fit_fe(D, dump, names(NCH0)); b1 <- coef(s1$m); keep <- !grepl("^tfe", names(b1)); V1 <- s1_vc2w(s1$m, s1$dd)
  bfe <- b1[!keep]; names(bfe) <- sub("^tfe", "", names(bfe)); bfe <- c(setNames(0, names(NCH0)[1]), bfe)
  b <- bfe[names(bfe) %in% REP & is.finite(bfe)]; b <- pmin(pmax(b - mean(b), -5), 5)
  bb <- TS[, .(dispnum, b = b[as.character(ccode)])][, .(bbar = mean(b, na.rm = TRUE)), by = dispnum][is.finite(bbar)]
  X1 <- merge(X, bb, by = "dispnum")
  F_t <- fit_mnl(X1, lapply(design, function(v) c(v, "bbar")), "cell")
  sc <- surf_contrasts(F_t)[, .(stage = "S2", alt, axis, gap, est, se)]
  th <- data.table(stage = "S2", alt = c("A", "B", "O"), axis = "theta", gap = NA_real_, est = F_t$b[c("A:bbar", "B:bbar", "O:bbar")], se = sqrt(diag(F_t$V))[c("A:bbar", "B:bbar", "O:bbar")])
  tie <- data.table(stage = "S2", alt = c("A", "B", "O"), axis = "tie_sum", gap = NA_real_, est = F_t$b[c("A:tie_sum_A", "B:tie_sum_B", "O:tie_sum_A")],
                    se = sqrt(diag(F_t$V))[c("A:tie_sum_A", "B:tie_sum_B", "O:tie_sum_A")])
  out <- rbind(s1_contrasts(b1[keep], V1[keep, keep], dump$roles), sc, th, tie)
  out[, `:=`(n_s1_events = nrow(ch), n_s2_rows = F_t$n, n_s2_A = sum(X1$y == "A"), n_s2_B = sum(X1$y == "B"), n_s2_O = sum(X1$y == "orig"))]
  out[]
}

## base build of one window x sample (measure-independent parts), from an R/s2_01_build_riskset.R output folder
s2_base <- function(build_dir) {
  META <- readRDS(file.path(build_dir, "s2_meta.rds")); W <- as.integer(META$window)
  ONE <- as.data.table(readRDS(file.path(build_dir, sprintf("s2_oneshot_w%d.rds", W))))
  CL <- as.data.table(readRDS(file.path(build_dir, "s2_classification_w0_6.rds")))[W == W]
  X <- ONE[in_ref == TRUE]; for (v in names(X)) if (is.logical(X[[v]])) X[, (v) := as.numeric(get(v))]
  ROS <- as.data.table(readRDS(file.path(build_dir, "s2_rosters_carry.rds")))
  E <- as.data.table(readRDS(file.path(build_dir, "s2_main_edges.rds")))
  list(meta = META, window = W, X = X, SI = readRDS(file.path(build_dir, "s2_surface_inputs.rds")), row0 = ONE[, .(dispnum, ccode, onset_year)],
       R0 = ROS[m == 0L], cover = sort(unique(E$year)), TS = CL[baseline == TRUE & side_rel == "B", .(dispnum, ccode)], main_edges = E)
}
