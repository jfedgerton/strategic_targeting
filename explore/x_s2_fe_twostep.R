## explore/x_s2_fe_twostep.R -- EXPLORATORY, PRE-FREEZE (PI 2026-11-22: "just to see"). Not part of the pre-specified pipeline.
## Two-step fixed-effects stand-in for the joint AME's shared target effect:
##   step 1  stage-1 conditional logit (init types, W = 1, fatal, main measure; M1 gap terms + Z) with target fixed effects -> b_j
##   step 2  stage-2 four-category logit (A, B, orig, out) on one-shot complete-horizon rows (stage-1 fatal set, W = 1), with
##           bbar = mean b_j over the dispute's original target set entering A, B and orig (theta_A, theta_B, theta_O)
##   inference: dispute-cluster bootstrap refitting BOTH steps on the same dispute draw (generated-regressor uncertainty), seed 123
## Usage: Rscript explore/x_s2_fe_twostep.R DUMP_INIT POSITIONS POSCOL S2_DIR OUT_DIR [B]
suppressPackageStartupMessages({ library(data.table); library(survival) })
a <- commandArgs(TRUE); DUMP <- a[1]; POSF <- a[2]; POSCOL <- a[3]; S2 <- a[4]; OUT <- a[5]; B <- if (length(a) >= 6) as.integer(a[6]) else 200L
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); set.seed(123)
LOG <- file.path(OUT, "x_s2_fe_twostep.log"); cat("", file = LOG)
logf <- function(...) { s <- sprintf(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOG, append = TRUE) }
logf("EXPLORATORY, PRE-FREEZE two-step FE run %s; B = %d", format(Sys.time()), B)

## ---- step-1 data: main-measure dR (positions at t-1, scaled by the gap scale of year t; as R/07) ---------------------------------
d <- readRDS(DUMP); stopifnot(d$types == "init", as.integer(d$coorig_months) == 1L)
D <- as.data.table(d$D); Z <- d$Z_cols
pz <- readRDS(POSF); key <- as.data.table(pz$key); R <- pz$R[, POSCOL]
gs <- function(v) { n <- length(v); sqrt(2 * n * sum((v - mean(v))^2) / (n * (n - 1))) }
sc <- tapply(R, key$year, gs)
kk <- paste(key$ccode, key$year)
i <- match(paste(D$c1, D$yr - 1L), kk); j <- match(paste(D$c2, D$yr - 1L), kk)
stopifnot(!anyNA(i), !anyNA(j))
D[, dR := (R[i] - R[j]) / sc[as.character(yr)]]
D[, `:=`(dM2 = dM^2, dR2 = dR^2, dMR = dM * dR)]
GAP <- c("dM", "dM2", "dR", "dR2", "dMR"); ROLES <- d$roles
for (g in 0:1) for (v in GAP) D[, paste0(v, "_", ROLES[g + 1]) := get(v) * (type == g)]
GT <- as.vector(outer(GAP, ROLES, paste, sep = "_"))
ch <- D[chosen == 1L]
## target FE for every target chosen at least once; never-chosen alternatives are the reference (their b_j is not estimable)
tg_ch <- sort(unique(ch$c2))
D[, tfe := ifelse(c2 %in% tg_ch, as.character(c2), "ref")]; D[, tfe := relevel(factor(tfe), "ref")]
logf("step 1: events %d, rows %d, targets chosen %d (chosen >= 2: %d), never-chosen alternatives pooled as reference %d", nrow(ch), nrow(D), length(tg_ch), sum(table(ch$c2) >= 2), uniqueN(D[tfe == "ref", c2]))

fit1 <- function(dd, fe = TRUE) {
  f <- as.formula(paste("Surv(rep(1, nrow(dd)), chosen) ~", paste(c(GT, Z, if (fe) "tfe"), collapse = " + "), "+ strata(ev_id)"))
  suppressWarnings(coxph(f, data = dd, weights = w, method = "breslow", control = coxph.control(iter.max = 50)))
}
m0 <- fit1(D, fe = FALSE); m1 <- fit1(D, fe = TRUE)
bfe <- coef(m1)[grepl("^tfe", names(coef(m1)))]; names(bfe) <- sub("^tfe", "", names(bfe))
logf("step 1 without FE: dR_bil %.3f, dR_coal %.3f | with target FE: dR_bil %.3f, dR_coal %.3f; FE NA (not identified) %d; FE range %.2f to %.2f",
     coef(m0)["dR_bil"], coef(m0)["dR_coal"], coef(m1)["dR_bil"], coef(m1)["dR_coal"], sum(is.na(bfe)), min(bfe, na.rm = TRUE), max(bfe, na.rm = TRUE))
logf("step 1 with FE: |b_j| > 10 (quasi-separated targets) %d", sum(abs(bfe) > 10, na.rm = TRUE))

## ---- step-2 data -------------------------------------------------------------------------------------------------------------------
ONE <- as.data.table(readRDS(file.path(S2, "s2_oneshot_w1.rds"))); CL <- as.data.table(readRDS(file.path(S2, "s2_classification_w0_6.rds")))[W == 1L]
TS <- CL[baseline == TRUE & side_rel == "B", .(dispnum, ccode)]
X <- ONE[complete == TRUE]
X[, `:=`(lcap_A = log1p(min_capdist_A), lcap_B = log1p(min_capdist_B), lcinc = log(cinc_k))]
need <- c("tie_sum_A", "tie_sum_B", "tgt_cent_mean", "lcap_A", "lcap_B", "ally_side_A", "ally_side_B", "likely", "lcinc")
n0 <- nrow(X); X <- X[complete.cases(X[, ..need]) & is.finite(lcinc)]
logf("step 2: complete-horizon rows %d; dropped for missing covariates %d; outcomes %s", n0, n0 - nrow(X), paste(sprintf("%s %d", names(table(X$y)), table(X$y)), collapse = ", "))

NCH <- table(ch$c2); REP <- names(NCH)[NCH >= 2]
## b_j used only for targets chosen >= 2 times (one-off targets are not usefully identified) and clipped to [-5, 5]
bbar_of <- function(b) { b <- pmin(pmax(b[names(b) %in% REP], -5), 5); tb <- TS[, .(dispnum, b = b[as.character(ccode)])]; tb[, .(bbar = mean(b, na.rm = TRUE), n_b = sum(!is.na(b))), by = dispnum] }
## four-category logit: utilities relative to out
VA <- c("tie_sum_A", "tgt_cent_mean", "lcap_A", "ally_side_A", "likely", "lcinc")
VB <- c("tie_sum_B", "tgt_cent_mean", "lcap_B", "ally_side_B", "likely", "lcinc")
VO <- c("tie_sum_A", "lcap_A", "ally_side_A")   # likely dropped: every originator entrant is in the likely tier (separation)
fit2 <- function(X, wt, with_b) {
  va <- c(VA, if (with_b) "bbar"); vb <- c(VB, if (with_b) "bbar"); vo <- c(VO, if (with_b) "bbar")
  ## covariates centred and scaled by fixed constants (CS, from the full point-fit data) so BFGS is well conditioned; slopes are reported per
  ## original unit (divided back), intercepts are at the covariate means
  cs <- function(M) sweep(sweep(M, 2, CS$m[colnames(M)]), 2, CS$s[colnames(M)], "/")
  XA <- cbind(1, cs(as.matrix(X[, ..va]))); XB <- cbind(1, cs(as.matrix(X[, ..vb]))); XO <- cbind(1, cs(as.matrix(X[, ..vo])))
  yA <- X$y == "A"; yB <- X$y == "B"; yO <- X$y == "orig"
  ka <- ncol(XA); kb <- ncol(XB); ko <- ncol(XO)
  nll <- function(p) { uA <- XA %*% p[1:ka]; uB <- XB %*% p[ka + 1:kb]; uO <- XO %*% p[ka + kb + 1:ko]
    mx <- pmax(0, uA, uB, uO); lse <- mx + log(exp(-mx) + exp(uA - mx) + exp(uB - mx) + exp(uO - mx))
    -sum(wt * (yA * uA + yB * uB + yO * uO - lse)) }
  gr <- function(p) { uA <- XA %*% p[1:ka]; uB <- XB %*% p[ka + 1:kb]; uO <- XO %*% p[ka + kb + 1:ko]
    mx <- pmax(0, uA, uB, uO); den <- exp(-mx) + exp(uA - mx) + exp(uB - mx) + exp(uO - mx)
    pA <- exp(uA - mx) / den; pB <- exp(uB - mx) / den; pO <- exp(uO - mx) / den
    -c(crossprod(XA, wt * (yA - pA)), crossprod(XB, wt * (yB - pB)), crossprod(XO, wt * (yO - pO))) }
  st <- c(log(mean(yA)), rep(0, ka - 1), log(mean(yB)), rep(0, kb - 1), log(mean(yO)), rep(0, ko - 1))
  o <- optim(st, nll, gr, method = "BFGS", control = list(maxit = 5000, reltol = 1e-14))
  gmax <- max(abs(gr(o$par)))
  sc <- c(1, CS$s[va], 1, CS$s[vb], 1, CS$s[vo]); par <- o$par / sc
  nm <- c(paste0("A:", c("(int@mean)", va)), paste0("B:", c("(int@mean)", vb)), paste0("O:", c("(int@mean)", vo)))
  list(coef = setNames(par, nm), conv = if (o$convergence == 0 && gmax < 1e-3 * max(1, sum(wt) / 1000)) 0L else 1L, ll = -o$value, gmax = gmax)
}
bb <- bbar_of(bfe); X0 <- merge(X, bb, by = "dispnum", all.x = TRUE)
allv <- unique(c(VA, VB, VO, "bbar")); CS <- list(m = sapply(allv, function(v) mean(X0[[v]], na.rm = TRUE)), s = sapply(allv, function(v) sd(X0[[v]], na.rm = TRUE)))
logf("step 2: repeat targets used for b_j %d; clipped at +-5 %d; disputes with no repeat-target b_j in the target set (dropped) %d", length(REP), sum(abs(bfe[names(bfe) %in% REP]) > 5), 0L); logf("step 2: disputes dropped (no usable b_j) %d", uniqueN(X0[is.na(bbar), dispnum]))
X0 <- X0[!is.na(bbar)]
f0 <- fit2(X0, rep(1, nrow(X0)), FALSE); f1 <- fit2(X0, rep(1, nrow(X0)), TRUE)
logf("point fits converged: no-theta %d (max |grad| %.1e), theta %d (%.1e); LR (theta block, 3 df) %.2f", f0$conv, f0$gmax, f1$conv, f1$gmax, 2 * (f1$ll - f0$ll))

## ---- bootstrap over disputes (both steps) ----------------------------------------------------------------------------------------
disp_all <- sort(unique(c(D$disp, X0$dispnum)))
MM <- lapply(seq_len(B), function(b) { m <- tabulate(sample.int(length(disp_all), length(disp_all), replace = TRUE), length(disp_all)); names(m) <- disp_all; m })
t0 <- Sys.time()
one_draw <- function(m) {
  Db <- copy(D); Db[, w := w * m[as.character(disp)]]; Db <- Db[w > 0]
  r <- tryCatch({
    mb <- fit1(Db, TRUE); be <- coef(mb)[grepl("^tfe", names(coef(mb)))]; names(be) <- sub("^tfe", "", names(be))
    Xb <- merge(X, bbar_of(be), by = "dispnum", all.x = TRUE); Xb <- Xb[!is.na(bbar)]
    wt <- m[as.character(Xb$dispnum)]; keep <- wt > 0; Xb <- Xb[keep]; wt <- wt[keep]
    g0 <- fit2(Xb, wt, FALSE); g1 <- fit2(Xb, wt, TRUE)
    if (g1$conv != 0 || g0$conv != 0) stop("no convergence")
    c(g1$coef, setNames(g0$coef, paste0("noTheta.", names(g0$coef))), dR_bil_fe = unname(coef(mb)["dR_bil"]), dR_coal_fe = unname(coef(mb)["dR_coal"]))
  }, error = function(e) NULL)
  r
}
draws <- parallel::mclapply(MM, one_draw, mc.cores = 2L)
fails <- sum(vapply(draws, is.null, TRUE)); draws <- draws[!vapply(draws, is.null, TRUE)]
logf("bootstrap time %.1f min, failures %d", as.numeric(difftime(Sys.time(), t0, units = "mins")), fails)
DR <- do.call(rbind, draws)
pt <- c(f1$coef, setNames(f0$coef, paste0("noTheta.", names(f0$coef))), dR_bil_fe = unname(coef(m1)["dR_bil"]), dR_coal_fe = unname(coef(m1)["dR_coal"]))
res <- data.table(term = names(pt), est = pt, lo = apply(DR[, names(pt)], 2, quantile, .025), hi = apply(DR[, names(pt)], 2, quantile, .975),
                  share_pos = colMeans(DR[, names(pt)] > 0))
fwrite(res, file.path(OUT, "x_s2_fe_twostep_coefs.csv"))
logf("bootstrap: %d successful of %d", nrow(DR), B)
for (k in seq_len(nrow(res))) logf("%-28s %8.3f  [%7.3f, %7.3f]  share>0 %.2f", res$term[k], res$est[k], res$lo[k], res$hi[k], res$share_pos[k])
logf("X_DONE")
