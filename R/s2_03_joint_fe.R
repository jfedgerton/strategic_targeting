## R/s2_03_joint_fe.R -- two-step fixed-effects joint model (PI 2026-11-22): the frequentist stand-in for the AME's shared target effect.
##   Step 1  stage-1 conditional logit (init types, W = 1, fatal, main measure; type-specific quadratic gap surface + Z_DIRECTED), without and
##           with target fixed effects b_j (targets chosen at least once; never-chosen alternatives pooled as the reference).
##   Step 2  stage-2 reference model (R/lib/s2_fit.R fit_mnl, frozen S2_REF_DESIGN, availability mask) plus bbar = mean b_j over the dispute's
##           original target set, entering A, B and O (theta_A, theta_B, theta_O). b_j is used only for targets chosen >= 2 times and clipped
##           to [-5, 5] (rule fixed in explore/x_s2_fe_twostep.R before these fits); disputes whose target set has no usable b_j are dropped.
##   Inference: stage 1 two-way CGM (initiator x dispute, R/07 vc_2w); stage 2 two-way CGM conditional on bbar; and a dispute-cluster
##   bootstrap refitting BOTH steps on the same dispute draw (generated-regressor uncertainty), percentile 95% intervals, seed 123.
## Usage: Rscript R/s2_03_joint_fe.R --dump=long_init_*.rds --positions=positions.rds --poscol=COL --build=DIR --out=DIR [--B=200]
suppressPackageStartupMessages({ library(data.table); library(survival) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file)), "lib", "s2.R")); source(file.path(dirname(normalizePath(this_file)), "lib", "s2_fit.R"))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
DUMP <- opt("dump"); POSF <- opt("positions"); POSCOL <- opt("poscol"); BUILD <- opt("build"); OUT <- opt("out"); B <- as.integer(opt("B", "200"))
stopifnot(!is.null(DUMP), !is.null(POSF), !is.null(POSCOL), !is.null(BUILD), !is.null(OUT), B >= 0)
W_REF <- 1L; dir.create(OUT, showWarnings = FALSE, recursive = TRUE); set.seed(123)
LOG <- file.path(OUT, "s2_joint_fe.log"); cat("", file = LOG)
logf <- function(...) { s <- sprintf(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOG, append = TRUE) }
logf("two-step FE joint model %s; dump %s; build %s; B = %d", format(Sys.time()), basename(DUMP), BUILD, B)

## ---- step 1 data (as R/07: positions at t-1, dR scaled by the gap scale of year t) ---------------------------------------------------
d <- readRDS(DUMP); stopifnot(d$types == "init", as.integer(d$coorig_months) == 1L)
D <- as.data.table(d$D); Z <- d$Z_cols; ROLES <- d$roles
pz <- readRDS(POSF); key <- as.data.table(pz$key); R <- pz$R[, POSCOL]; kk <- paste(key$ccode, key$year)
sc <- tapply(R, key$year, s2_gap_scale)
i <- match(paste(D$c1, D$yr - 1L), kk); j <- match(paste(D$c2, D$yr - 1L), kk); stopifnot(!anyNA(i), !anyNA(j))
D[, dR := (R[i] - R[j]) / sc[as.character(yr)]]
D[, `:=`(dM2 = dM^2, dR2 = dR^2, dMR = dM * dR)]
GAP <- c("dM", "dM2", "dR", "dR2", "dMR")
for (g in 0:1) for (v in GAP) D[, paste0(v, "_", ROLES[g + 1]) := get(v) * (type == g)]
GT <- as.vector(outer(GAP, ROLES, paste, sep = "_"))
ch <- D[chosen == 1L]; tg <- sort(unique(ch$c2))
## target FE (fix 2026-11-22): a target never chosen in the sample has fitted FE -> -Inf (quasi-separation), which made the FE level drift and the
## CGM variance explode. Its rows contribute nothing in the limit, so the FE model is fitted on the rows of chosen targets only (the exact
## limiting likelihood), with the most-chosen target as the reference level; the no-FE model keeps every alternative
NCH0 <- sort(table(ch$c2), decreasing = TRUE)
D[, tfe := factor(as.character(c2), levels = names(NCH0))]
NCH <- table(ch$c2); REP <- names(NCH)[NCH >= 2]
logf("step 1: events %d (bilateral %d, coalition %d), rows %d; targets chosen %d (>= 2: %d)", nrow(ch), sum(ch$type == 0), sum(ch$type == 1), nrow(D), length(tg), length(REP))
fit1 <- function(dd, fe) { if (fe) { tgb <- unique(dd[chosen == 1L, c2]); dd <- dd[c2 %in% tgb]; dd[, tfe := droplevels(factor(as.character(c2), levels = intersect(names(NCH0), as.character(tgb))))]
    dd <- dd[dd[, .(k = .N), by = ev_id][k > 1L], on = "ev_id"] }
  f <- as.formula(paste("Surv(rep(1, nrow(dd)), chosen) ~", paste(c(GT, Z, if (fe) "tfe"), collapse = " + "), "+ strata(ev_id)"))
  suppressWarnings(coxph(f, data = dd, weights = w, method = "breslow", control = coxph.control(iter.max = 50))) }
vc_cl <- function(m, g) { Dd <- rowsum(residuals(m, type = "dfbeta", weighted = TRUE), g, reorder = FALSE); crossprod(Dd) }
vc_2w <- function(m, dd) { V <- vc_cl(m, dd$init) + vc_cl(m, dd$disp) - vc_cl(m, dd$init_disp); e <- eigen(V, symmetric = TRUE)
  if (any(e$values < 0)) V <- e$vectors %*% diag(pmax(e$values, 0)) %*% t(e$vectors); dimnames(V) <- list(names(coef(m)), names(coef(m))); attr(V, "n_neg") <- sum(e$values < 0); V }
s1_contr <- function(b, V = NULL) { out <- list()
  for (r in ROLES) for (ax in c("M", "R")) for (sg in c(1, -1)) { a <- setNames(numeric(length(b)), names(b))
    a[paste0("d", ax, "_", r)] <- sg; a[paste0("d", ax, "2_", r)] <- 1
    out[[length(out) + 1L]] <- data.table(type = r, axis = ax, gap = sg, est = sum(a * b), se = if (is.null(V)) NA_real_ else sqrt(max(drop(t(a) %*% V %*% a), 0))) }
  rbindlist(out) }
m0 <- fit1(D, FALSE); m1 <- fit1(D, TRUE)
D1 <- D[c2 %in% tg]; D1 <- D1[D1[, .(k = .N), by = ev_id][k > 1L], on = "ev_id"]
logf("step 1 target-FE model: rows of chosen targets only %d of %d (events kept %d)", nrow(D1), nrow(D), uniqueN(D1$ev_id))
b0 <- coef(m0); V0 <- vc_2w(m0, D); b1 <- coef(m1); keep1 <- !grepl("^tfe", names(b1))
V1 <- vc_2w(m1, D1)
S1 <- rbind(data.table(model = "no FE", term = names(b0), est = b0, se = sqrt(diag(V0))),
            data.table(model = "target FE", term = names(b1)[keep1], est = b1[keep1], se = sqrt(diag(V1))[keep1]))
S1C <- rbind(s1_contr(b0, V0)[, model := "no FE"], s1_contr(b1[keep1], V1[keep1, keep1])[, model := "target FE"])
logf("step 1 two-way CGM negative eigenvalues: no FE %d, target FE %d", attr(V0, "n_neg"), attr(V1, "n_neg"))
bfe <- b1[!keep1]; names(bfe) <- sub("^tfe", "", names(bfe)); bfe <- c(setNames(0, names(NCH0)[1]), bfe)   # the reference target has b = 0
logf("step 1 target FE (relative to the most-chosen target): range %.2f to %.2f; centred on repeat targets, |b_j| > 5 (clipped) %d", min(bfe, na.rm = TRUE), max(bfe, na.rm = TRUE), { bb <- bfe[names(bfe) %in% REP]; sum(abs(bb - mean(bb, na.rm = TRUE)) > 5, na.rm = TRUE) })

## ---- step 2 -----------------------------------------------------------------------------------------------------------------------
ONE <- as.data.table(readRDS(file.path(BUILD, "s2_oneshot_w1.rds"))); CL <- as.data.table(readRDS(file.path(BUILD, "s2_classification_w0_6.rds")))[W == 1L]
TS <- CL[baseline == TRUE & side_rel == "B", .(dispnum, ccode)]
X <- ONE[in_ref == TRUE]
for (v in names(X)) if (is.logical(X[[v]]) && v %in% unlist(S2_REF_DESIGN)) X[, (v) := as.numeric(get(v))]
## only differences among chosen targets are identified: b_j is centred on the repeat targets' mean before clipping (location is absorbed by
## the stage-2 intercepts)
bbar_of <- function(b) { b <- b[names(b) %in% REP & is.finite(b)]; b <- pmin(pmax(b - mean(b), -5), 5); tb <- TS[, .(dispnum, b = b[as.character(ccode)])]; tb[, .(bbar = mean(b, na.rm = TRUE)), by = dispnum][is.finite(bbar)] }
DES_T <- lapply(S2_REF_DESIGN, function(v) c(v, "bbar"))
X1 <- merge(X, bbar_of(bfe), by = "dispnum")
logf("step 2: reference rows %d; with a usable bbar %d (%d disputes dropped); outcomes %s", nrow(X), nrow(X1), uniqueN(X$dispnum) - uniqueN(X1$dispnum),
     paste(sprintf("%s %d", names(table(X1$y)), table(X1$y)), collapse = ", "))
F_nt <- fit_mnl(X1, S2_REF_DESIGN, "no theta (same rows)"); F_t <- fit_mnl(X1, DES_T, "theta")
for (al in names(F_t$drops)) if (length(F_t$drops[[al]])) logf("separation rule: dropped from %s: %s", al, paste(F_t$drops[[al]], collapse = ", "))
logf("step 2: LR theta block (3 df) %.2f; theta (two-way SE conditional on bbar): A %.3f (%.3f), B %.3f (%.3f), O %.3f (%.3f)", 2 * (F_t$ll - F_nt$ll),
     F_t$b["A:bbar"], sqrt(F_t$V["A:bbar", "A:bbar"]), F_t$b["B:bbar"], sqrt(F_t$V["B:bbar", "B:bbar"]), F_t$b["O:bbar"], sqrt(F_t$V["O:bbar", "O:bbar"]))

## ---- bootstrap over disputes (both steps) ----------------------------------------------------------------------------------------
disp_all <- sort(unique(c(D$disp, X$dispnum)))
MM <- lapply(seq_len(B), function(b) { m <- tabulate(sample.int(length(disp_all), length(disp_all), replace = TRUE), length(disp_all)); names(m) <- disp_all; m })
fitw <- function(Xw, des, wt) {   # weighted refit: replicate rows by their dispute multiplicity
  Xr <- Xw[rep(seq_len(.N), wt)]; fit_mnl(Xr, des, "boot", vc = FALSE) }
one_draw <- function(m) tryCatch({
  Db <- copy(D); Db[, w := w * m[as.character(disp)]]; Db <- Db[w > 0]
  mb <- fit1(Db, TRUE); bb <- coef(mb); be <- bb[grepl("^tfe", names(bb))]; names(be) <- sub("^tfe", "", names(be))
  refb <- setdiff(intersect(names(NCH0), as.character(unique(Db[chosen == 1L, c2]))), names(be))[1]; be <- c(setNames(0, refb), be)
  s1 <- s1_contr(bb[!grepl("^tfe", names(bb))])
  Xb <- merge(X, bbar_of(be), by = "dispnum"); wt <- m[as.character(Xb$dispnum)]; Xb <- Xb[wt > 0]; wt <- wt[wt > 0]
  ft <- fitw(Xb, DES_T, wt)
  if (!"A:bbar" %in% names(ft$b)) stop("bbar dropped by the separation rule in this draw")
  c(setNames(s1$est, paste("S1", s1$type, s1$axis, s1$gap, sep = ".")), setNames(surf_contrasts(ft)$est, paste("S2", surf_contrasts(ft)$alt, surf_contrasts(ft)$axis, surf_contrasts(ft)$gap, sep = ".")),
    ft$b[c("A:bbar", "B:bbar", "O:bbar", "A:tie_sum_A", "B:tie_sum_B", "O:tie_sum_A", "B:share_B")])
}, error = function(e) { message("draw failed: ", conditionMessage(e)); NULL })
t0 <- Sys.time(); DR <- parallel::mclapply(MM, one_draw, mc.cores = 2L)
saveRDS(DR, file.path(OUT, "boot_draws_raw.rds"))
fails <- sum(vapply(DR, is.null, TRUE)); DR <- do.call(rbind, DR[!vapply(DR, is.null, TRUE)])
logf("bootstrap: %d of %d draws succeeded (%.1f min)", nrow(DR), B, as.numeric(difftime(Sys.time(), t0, units = "mins")))
s1p <- s1_contr(b1[keep1]); s2p <- surf_contrasts(F_t)
PT <- c(setNames(s1p$est, paste("S1", s1p$type, s1p$axis, s1p$gap, sep = ".")), setNames(s2p$est, paste("S2", s2p$alt, s2p$axis, s2p$gap, sep = ".")),
        F_t$b[c("A:bbar", "B:bbar", "O:bbar", "A:tie_sum_A", "B:tie_sum_B", "O:tie_sum_A", "B:share_B")])
miss <- setdiff(names(PT), colnames(DR)); if (length(miss)) stop("bootstrap draws lack: ", paste(miss, collapse = ", "), " | draw names: ", paste(head(colnames(DR), 40), collapse = ", "))
BOOT <- data.table(q = names(PT), est = PT, lo = apply(DR[, names(PT), drop = FALSE], 2, quantile, .025), hi = apply(DR[, names(PT), drop = FALSE], 2, quantile, .975), share_pos = colMeans(DR[, names(PT), drop = FALSE] > 0))
for (r in seq_len(nrow(BOOT))) logf("  %-22s %8.3f  [%7.3f, %7.3f]  share>0 %.2f", BOOT$q[r], BOOT$est[r], BOOT$lo[r], BOOT$hi[r], BOOT$share_pos[r])
saveRDS(list(S1 = S1, S1C = S1C, F_t = F_t, F_nt = F_nt, BOOT = BOOT, n1 = nrow(ch), n1_types = table(ch$type), n_rows1 = nrow(D), B = B, fails = fails), file.path(OUT, "joint_fe.rds"))
logf("S2_JOINT_FE_DONE")
