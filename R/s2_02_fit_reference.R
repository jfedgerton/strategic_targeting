## R/s2_02_fit_reference.R -- stage-2 reference one-shot model and its robustness checks (memo section 13, surface amendment; rev 5).
## Four-category alternative-specific logit (A = join initiating side, B = join target side, O = originator entry; stay out = reference) on the
## reference estimation cohort (in_ref) with the frozen design S2_REF_DESIGN: a RESTRICTED opponent-surface model (the stage-1 quadratic
## surface in the joiner's gaps to the side it fights, which constrains own and opponent position effects to be equal and opposite), the
## stage-1 directed controls toward the opponent, and recruitment terms toward the joined side. Inference: two-way cluster-robust (candidate
## state x dispute; covariance-level PSD repair, as stage 1) and a two-way pigeonhole bootstrap.
## Reported: coefficients; surface utility contrasts (log-odds relative to staying out, +1 / -1 gap on one axis); predicted-probability
## contrasts (25th -> 75th percentile of each gap, parametric simulation intervals); the 15-df Wald test that all surface coefficients are 0
## (not a test of equal surfaces across alternatives).
## Robustness (GPT surface review, rev 5):
##   A  opposing-side variants: 40 nominal surface definitions (aggregation x members x lag x dR scale). Membership and timing variants use
##      opposing-side controls rebuilt for the same members and timing; aggregation variants change the surface only (min/max describe a
##      composite opponent, lead and mean actual quantities). Each variant is compared with the reference refitted on the same rows.
##   L  the linear design (a different restricted model, NOT nested in the surface model; no LR test is reported between them).
##   B1 / B2 artificial-assignment DIAGNOSTICS (not formal tests): candidates' gap vectors shuffled within dispute (surface-block LR); A/B
##      labels shuffled among third-party entrants within dispute (side-choice statistic = LR of the A and B surface terms). Reported only
##      when >= 80% of draws are usable; otherwise "unavailable".
##   B4 model-based null (formal test): outcomes simulated from the fitted no-surface model with the same availability, both models refitted;
##      p = share of simulated surface-block LRs >= observed.
##   B3 position-level contamination (positions perturbed with alliance- / contiguity-correlated noise, both signs, k = 0.25 / 0.5 / 1, plus
##      20 random-noise draws at each k), surface rebuilt, under two scalings: fixed baseline gap scales and gap scales recomputed from the
##      perturbed positions (as stage 1). Recruitment edges are unchanged.
## Usage: Rscript R/s2_02_fit_reference.R --build=DIR --out=DIR [--n_perm=100] [--n_null=200] [--n_boot=200] [--n_noise=20] [--cores=1]
##          [--checks=all|none]
suppressPackageStartupMessages(library(data.table))
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file)), "lib", "s2.R")); source(file.path(dirname(normalizePath(this_file)), "lib", "s2_fit.R"))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
fl <- sub("^--([a-z_]+)=.*$", "\\1", args); if (any(!fl %in% c("build", "out", "n_perm", "n_null", "n_boot", "n_noise", "cores", "checks")) || anyDuplicated(fl)) stop("unknown or duplicated option")
BUILD <- opt("build"); OUT <- opt("out"); N_PERM <- as.integer(opt("n_perm", "100")); N_NULL <- as.integer(opt("n_null", "200"))
N_BOOT <- as.integer(opt("n_boot", "200")); N_NOISE <- as.integer(opt("n_noise", "20")); CORES <- as.integer(opt("cores", "1")); CHECKS <- opt("checks", "all")
stopifnot(!is.null(BUILD), !is.null(OUT), CHECKS %in% c("all", "none"), N_PERM >= 0, N_NULL >= 0, N_BOOT >= 0, CORES >= 1)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); set.seed(123); RNGkind("L'Ecuyer-CMRG")
LOG <- file.path(OUT, "s2_fit.log"); cat("", file = LOG)
logf <- function(...) { s <- sprintf(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOG, append = TRUE) }
GIT <- tryCatch(system2("git", c("-C", dirname(normalizePath(this_file)), "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
logf("stage-2 reference fit %s; build %s; git %s; R %s; cores %d", format(Sys.time()), BUILD, GIT, getRversion(), CORES)
pmap <- function(X, FUN) if (CORES > 1L) parallel::mclapply(X, FUN, mc.cores = CORES, mc.set.seed = TRUE) else lapply(X, FUN)
usable <- function(v, n) { k <- sum(!is.na(v)); list(ok = k >= 0.8 * n, k = k) }

## ---- data ----------------------------------------------------------------------------------------------------------------------------
META <- readRDS(file.path(BUILD, "s2_meta.rds")); stopifnot(META$mode == "release")
W_REF <- as.integer(META$window)
ONE <- as.data.table(readRDS(file.path(BUILD, sprintf("s2_oneshot_w%d.rds", W_REF))))
VARS <- as.data.table(readRDS(file.path(BUILD, "s2_surface_variants.rds"))); SI <- readRDS(file.path(BUILD, "s2_surface_inputs.rds"))
OCV <- readRDS(file.path(BUILD, "s2_opp_control_variants.rds"))
X0 <- ONE[in_ref == TRUE]
for (v in names(X0)) if (is.logical(X0[[v]]) && v %in% unlist(c(S2_REF_DESIGN, S2_REF_DESIGN_LINEAR))) X0[, (v) := as.numeric(get(v))]
logf("build: fatal_def %s, window %d; reference cohort %d rows, %d disputes, %d candidate states; outcomes %s", META$fatal_def, W_REF, nrow(X0), uniqueN(X0$dispnum), uniqueN(X0$ccode),
     paste(sprintf("%s %d", names(table(X0$y)), table(X0$y)), collapse = ", "))
GCOLS <- c("dM_kA", "dR_kA", "dM_kB", "dR_kB"); SURFC <- c(paste0("sA_", c("dM", "dM2", "dR", "dR2", "dMR")), paste0("sB_", c("dM", "dM2", "dR", "dR2", "dMR")))
X0 <- add_surface(X0)

## ---- 1. reference fit --------------------------------------------------------------------------------------------------------------
F0 <- fit_mnl(X0, S2_REF_DESIGN, "reference")
for (al in names(F0$drops)) if (length(F0$drops[[al]])) logf("separation rule: dropped from %s: %s (restricted specification)", al, paste(F0$drops[[al]], collapse = ", "))
se0 <- sqrt(diag(F0$V))
COEF <- data.table(term = names(F0$b), est = F0$b, se = se0, z = F0$b / se0, p = 2 * pnorm(-abs(F0$b / se0)))
fwrite(COEF, file.path(OUT, "reference_coefs.csv"))
SC0 <- surf_contrasts(F0); SC0[, `:=`(lo = est - 1.96 * se, hi = est + 1.96 * se)]
fwrite(SC0, file.path(OUT, "reference_surface_contrasts.csv"))
wb <- wald(F0, surf_block)
logf("reference fit: n %d, logLik %.2f, max |grad| %.1e, Hessian condition %.1e, covariance negative eigenvalues %d (share of absolute eigenvalue mass %.2e); Wald (all 15 surface coefficients = 0) chi2 %.2f, p = %.4f",
     F0$n, F0$ll, F0$gmax, F0$cond, attr(F0$V, "n_neg"), attr(F0$V, "neg_share"), wb["chi2"], wb["p"])
logf("surface utility contrasts (log-odds vs staying out; gap = joiner minus the side it fights; 95%% two-way intervals):")
for (r in seq_len(nrow(SC0))) logf("  %s  %s %+d: %7.3f  [%7.3f, %7.3f]  %s", SC0$alt[r], SC0$axis[r], SC0$gap[r], SC0$est[r], SC0$lo[r], SC0$hi[r], SC0$status[r])
PC0 <- prob_contrasts(F0, X0, 500L); fwrite(PC0, file.path(OUT, "reference_prob_contrasts.csv"))
logf("predicted-probability contrasts (25th -> 75th percentile of the gap; base = mean predicted probability; 95%% parametric-simulation intervals):")
for (r in seq_len(nrow(PC0))) logf("  %s  %s [%.2f -> %.2f]: base %.5f, change %+.5f [%+.5f, %+.5f]", PC0$alt[r], PC0$gap[r], PC0$q25[r], PC0$q75[r], PC0$base_p[r], PC0$contrast[r], PC0$lo[r], PC0$hi[r])
supp <- X0[, .(p01_kB_R = quantile(dR_kB, .01), p99_kB_R = quantile(dR_kB, .99), share_R_in_pm1 = mean(abs(dR_kB) <= 1), share_M_in_pm1 = mean(abs(dM_kB) <= 1))]
logf("support of the +-1 endpoints (gap to B): dR 1st-99th pct [%.2f, %.2f]; share of rows with |dR| <= 1 %.2f, |dM| <= 1 %.2f", supp$p01_kB_R, supp$p99_kB_R, supp$share_R_in_pm1, supp$share_M_in_pm1)
logf("coefficients (two-way SE):"); for (r in seq_len(nrow(COEF))) logf("  %-22s %8.3f  (%.3f)  p %.3f", COEF$term[r], COEF$est[r], COEF$se[r], COEF$p[r])
saveRDS(list(F0 = F0, design = S2_REF_DESIGN, window = W_REF, meta = META), file.path(OUT, "reference_fit.rds"))
if (CHECKS == "none") { logf("S2_FIT_DONE (checks skipped)"); quit(status = 0) }

## ---- 2. two-way pigeonhole bootstrap -------------------------------------------------------------------------------------------------
BT <- pmap(seq_len(N_BOOT), function(b) tryCatch({ f <- fit_mnl(X0, S2_REF_DESIGN, "boot", vc = FALSE, w = pigeonhole_w(X0), start = F0$par_scaled)
  sc <- surf_contrasts(f); c(setNames(sc$est, paste(sc$alt, sc$axis, sc$gap, sep = ".")), f$b) }, error = function(e) NULL))
nb <- sum(!vapply(BT, is.null, TRUE))
if (nb >= 0.8 * N_BOOT && nb > 0) { BTm <- do.call(rbind, lapply(BT[!vapply(BT, is.null, TRUE)], function(v) v[intersect(names(v), c(paste(SC0$alt, SC0$axis, SC0$gap, sep = "."), names(F0$b)))]))
  qq <- data.table(q = colnames(BTm), lo = apply(BTm, 2, quantile, .025, na.rm = TRUE), hi = apply(BTm, 2, quantile, .975, na.rm = TRUE), share_pos = colMeans(BTm > 0, na.rm = TRUE))
  fwrite(qq, file.path(OUT, "pigeonhole_bootstrap.csv"))
  logf("pigeonhole bootstrap: %d of %d draws usable; surface contrasts [2.5, 97.5] pct:", nb, N_BOOT)
  for (r in which(grepl("^[ABO]\\.[MR]\\.", qq$q))) logf("  %-10s [%7.3f, %7.3f]  share>0 %.2f", qq$q[r], qq$lo[r], qq$hi[r], qq$share_pos[r])
} else logf("pigeonhole bootstrap: UNAVAILABLE (%d of %d draws usable; rule >= 80%%)", nb, N_BOOT)

## ---- 3. A: opposing-side variants (controls rebuilt for membership/timing), common-cohort comparison; L: linear design --------------------
keys <- sub("^dR_kB__", "", grep("^dR_kB__", names(VARS), value = TRUE))
OPPC <- unique(c(sub("_B$", "", S2_Z_TOWARD("B")), sub("_A$", "", S2_Z_TOWARD("A"))))
oc_cols <- function(W) intersect(names(W), unlist(lapply(c("_A", "_B"), function(sf) paste0(c("land_side", "lcap", "dem_side", "dem_x", "defoff_side", "major_side", "major_x", "mid_tm1_side", "prior_side"), sf))))
VR <- pmap(keys, function(key) tryCatch({
  p <- strsplit(key, "_")[[1]]; ock <- paste(p[2], p[3], sep = "_")                       # agg_members_lag_scale -> members_lag
  Xv <- merge(X0[, !c(GCOLS, SURFC), with = FALSE], VARS[, c("dispnum", "ccode", paste0(GCOLS, "__", key)), with = FALSE], by = c("dispnum", "ccode"))
  setnames(Xv, paste0(GCOLS, "__", key), GCOLS)
  if (ock != "roster_l1") { W <- OCV[[ock]]; cc <- oc_cols(W); Xv[, (cc) := NULL]; Xv <- merge(Xv, W[, c("dispnum", "ccode", cc), with = FALSE], by = c("dispnum", "ccode")) }
  need <- unique(unlist(S2_REF_DESIGN)); need <- intersect(c(GCOLS, need), names(Xv))
  n_all <- nrow(Xv); Xv <- Xv[Reduce(`&`, lapply(c(GCOLS, intersect(need, oc_cols(Xv))), function(v) is.finite(Xv[[v]])))]; add_surface(Xv)
  f <- fit_mnl(Xv, S2_REF_DESIGN, key)
  fr <- fit_mnl(X0[paste(dispnum, ccode) %in% paste(Xv$dispnum, Xv$ccode)], S2_REF_DESIGN, paste(key, "ref-rows"))      # reference on the same rows
  merge(surf_contrasts(f)[, .(alt, axis, gap, est, se, status)], surf_contrasts(fr)[, .(alt, axis, gap, est_ref_same_rows = est)], by = c("alt", "axis", "gap"))[,
    `:=`(variant = key, controls = if (ock == "roster_l1") "reference" else paste("rebuilt", ock), n = f$n, dropped_rows = n_all - f$n)]
}, error = function(e) data.table(variant = key, status = paste("failed:", conditionMessage(e)))))
VRt <- rbindlist(VR, fill = TRUE); fwrite(VRt, file.path(OUT, "variants_surface_contrasts.csv"))
nfail <- sum(grepl("^failed", VRt$status)); ok <- VRt[status == "ok"]
ok[, `:=`(lo = est - 1.96 * se, hi = est + 1.96 * se)]
sumv <- ok[, .(variants = uniqueN(variant), same_sign = mean(sign(est) == sign(SC0$est[SC0$alt == alt[1] & SC0$axis == axis[1] & SC0$gap == gap[1]])), excl0 = mean(lo > 0 | hi < 0),
               same_sign_as_ref_same_rows = mean(sign(est) == sign(est_ref_same_rows)), min = min(est), max = max(est)), by = .(alt, axis, gap)]
logf("A: %d variants requested, %d failed (listed in variants_surface_contrasts.csv); per contrast: share with the reference sign / share whose interval excludes 0 / share agreeing with the reference refitted on the same rows / range:",
     length(keys), uniqueN(VRt[grepl("^failed", status), variant]))
for (r in seq_len(nrow(sumv))) logf("  %s %s %+d: same sign %.2f, excludes 0 %.2f, agrees on same rows %.2f, range [%.2f, %.2f]", sumv$alt[r], sumv$axis[r], sumv$gap[r], sumv$same_sign[r], sumv$excl0[r],
                                    sumv$same_sign_as_ref_same_rows[r], sumv$min[r], sumv$max[r])
FL <- tryCatch(fit_mnl(X0, S2_REF_DESIGN_LINEAR, "linear"), error = function(e) { logf("linear design failed: %s", conditionMessage(e)); NULL })
if (!is.null(FL)) { cl <- data.table(term = names(FL$b), est = FL$b, se = sqrt(diag(FL$V))); fwrite(cl, file.path(OUT, "linear_design_coefs.csv"))
  logf("L (linear design; a different restricted model, not nested): logLik %.2f vs surface %.2f (no LR test); position terms %s", FL$ll, F0$ll,
       paste(sprintf("%s %.2f (%.2f)", cl$term, cl$est, cl$se)[grepl("cent|ratio", cl$term)], collapse = "; ")) }

## ---- 4. B1 / B2 diagnostics and B4 model-based null ------------------------------------------------------------------------------------
NULL_DESIGN <- lapply(S2_REF_DESIGN, function(v) setdiff(v, c(S2_SURF_A, S2_SURF_B)))
F_null <- fit_mnl(X0, NULL_DESIGN, "null", vc = FALSE)
LR_OBS <- 2 * (F0$ll - F_null$ll)
logf("surface-block LR (15 df) %.2f", LR_OBS)
P1 <- unlist(pmap(seq_len(N_PERM), function(b) tryCatch({ Xp <- copy(X0); ord <- Xp[, .(i = .I, j = .I[sample.int(.N)]), by = dispnum]
  Xp[ord$i, (GCOLS) := X0[ord$j, ..GCOLS]]; add_surface(Xp); 2 * (fit_mnl(Xp, S2_REF_DESIGN, "B1", vc = FALSE, start = F0$par_scaled)$ll - F_null$ll) }, error = function(e) NA_real_)))
## B2 side-choice statistic (GPT P2-6): among third-party entrants (A or B) whose two sides were both available, a binary logit of A vs B on
## both surface blocks (sA_*, sB_*); statistic = LR of those 10 terms. Labels are shuffled among those entrants within dispute.
ENT <- X0[y %in% c("A", "B") & H > W_REF]
side_lr <- function(E) { E <- copy(E); E[, yA := as.integer(y == "A")]
  g1 <- suppressWarnings(glm(reformulate(SURFC, "yA"), data = E, family = binomial())); g0 <- glm(yA ~ 1, data = E, family = binomial())
  if (!g1$converged) return(NA_real_); 2 * as.numeric(logLik(g1) - logLik(g0)) }
LR_SIDE_OBS <- side_lr(ENT)
logf("B2 side-choice statistic: %d entrants with both sides available (A %d, B %d); observed LR (10 df) %.2f", nrow(ENT), sum(ENT$y == "A"), sum(ENT$y == "B"), LR_SIDE_OBS)
P2 <- unlist(pmap(seq_len(N_PERM), function(b) tryCatch({ E <- copy(ENT); E[, y := sample(y), by = dispnum]; side_lr(E) }, error = function(e) NA_real_)))
fwrite(data.table(draw = seq_len(N_PERM), B1_gaps_within_dispute = P1, B2_side_labels = P2), file.path(OUT, "diagnostic_permutations.csv"))
for (pp in list(list("B1 (gap vectors shuffled within dispute; surface-block LR)", P1, LR_OBS), list("B2 (A/B labels shuffled among entrants within dispute; side-choice LR)", P2, LR_SIDE_OBS))) {
  u <- usable(pp[[2]], N_PERM)
  if (u$ok) logf("%s DIAGNOSTIC: observed %.2f; share of %d usable draws >= observed %.3f; draw median %.2f", pp[[1]], pp[[3]], u$k, mean(pp[[2]] >= pp[[3]], na.rm = TRUE), median(pp[[2]], na.rm = TRUE))
  else logf("%s DIAGNOSTIC: UNAVAILABLE (%d of %d usable)", pp[[1]], u$k, N_PERM) }
PN <- mnl_predict(F_null, X0)
P4 <- unlist(pmap(seq_len(N_NULL), function(b) tryCatch({ Xs <- copy(X0); u <- runif(nrow(Xs)); cp <- t(apply(PN, 1, cumsum))
  Xs[, y := ifelse(u < cp[, "A"], "A", ifelse(u < cp[, "B"], "B", ifelse(u < cp[, "O"], "orig", "out")))]
  2 * (fit_mnl(Xs, S2_REF_DESIGN, "B4", vc = FALSE)$ll - fit_mnl(Xs, NULL_DESIGN, "B4-null", vc = FALSE)$ll) }, error = function(e) NA_real_)))
fwrite(data.table(draw = seq_len(N_NULL), lr = P4), file.path(OUT, "model_based_null.csv"))
u4 <- usable(P4, N_NULL)
if (u4$ok) logf("B4 model-based null (outcomes simulated from the fitted no-surface model): observed LR %.2f; p = %.4f (%d usable draws; null LR median %.2f, 95th pct %.2f)", LR_OBS,
                (1 + sum(P4 >= LR_OBS, na.rm = TRUE)) / (1 + u4$k), u4$k, median(P4, na.rm = TRUE), quantile(P4, .95, na.rm = TRUE)) else logf("B4 model-based null: UNAVAILABLE (%d of %d usable)", u4$k, N_NULL)

## ---- 5. B3 position-level contamination (fixed vs renormalized gap scales) ------------------------------------------------------------
PV <- as.data.table(SI$PV); CONF <- as.data.table(SI$conf)
PV <- merge(PV, CONF, by = c("ccode", "year"), all.x = TRUE); PV[is.na(n_def), n_def := 0L]; PV[is.na(n_contig), n_contig := 0L]
zyr <- function(x, yr) ave(x, yr, FUN = function(v) if (length(v) > 1 && sd(v) > 0) (v - mean(v)) / sd(v) else 0 * v)
PV <- PV[!is.na(R)]; PV[, `:=`(z_def = zyr(log1p(n_def), year), z_contig = zyr(log1p(n_contig), year), sdR = ave(R, year, FUN = sd))]
grid <- rbind(CJ(conf = c("def", "contig"), sgn = c(1, -1), k = c(0.25, 0.5, 1), draw = 1L), CJ(conf = "rand", sgn = 1, k = c(0.25, 0.5, 1), draw = seq_len(N_NOISE)))
grid <- rbind(grid[, scale := "fixed"], copy(grid)[, scale := "renormalized"])
MER <- pmap(seq_len(nrow(grid)), function(i) tryCatch({ g <- grid[i]; Pp <- copy(PV)
  z <- if (g$conf == "rand") { set.seed(123 + 1000 * g$draw + round(100 * g$k)); rnorm(nrow(Pp)) } else Pp[[paste0("z_", g$conf)]]
  Pp[, R := R + g$sgn * g$k * sdR * z]
  sR2 <- if (g$scale == "fixed") SI$sR else vapply(split(Pp$R, Pp$year), s2_gap_scale, 0)
  sv <- s2_surface(SI$rows, SI$mem_roster, Pp[, .(ccode, year, M, R)], SI$sM, sR2, "mean")
  Xm <- merge(X0[, !c(GCOLS, SURFC), with = FALSE], sv[, c("dispnum", "ccode", GCOLS), with = FALSE], by = c("dispnum", "ccode"))
  Xm <- Xm[is.finite(dR_kA) & is.finite(dR_kB) & is.finite(dM_kA) & is.finite(dM_kB)]; add_surface(Xm)
  rk <- suppressWarnings(cor(PV$R, Pp$R, method = "spearman"))
  f <- fit_mnl(Xm, S2_REF_DESIGN, "B3")
  surf_contrasts(f)[, `:=`(confounder = g$conf, sign = g$sgn, k = g$k, draw = g$draw, scale = g$scale, rank_corr = rk, n = f$n)]
}, error = function(e) data.table(confounder = grid$conf[i], sign = grid$sgn[i], k = grid$k[i], draw = grid$draw[i], scale = grid$scale[i], status = paste("failed:", conditionMessage(e)))))
MERt <- rbindlist(MER, fill = TRUE); fwrite(MERt, file.path(OUT, "contamination_surface_contrasts.csv"))
logf("B3: %d perturbed fits requested, %d failed; R +1 contrasts by condition (mean over random draws), A / B / O:", nrow(grid), sum(grepl("^failed", MERt$status)))
sm <- MERt[status == "ok" & axis == "R" & gap == 1, .(est = mean(est), rk = mean(rank_corr)), by = .(scale, confounder, sign, k, alt)]
for (g in split(sm, by = c("scale", "confounder", "sign", "k"))) logf("  %-12s %-6s %+d k %.2f (rank corr %.2f): A %.2f, B %.2f, O %.2f", g$scale[1], g$confounder[1], g$sign[1], g$k[1], g$rk[1],
  g[alt == "A", est], g[alt == "B", est], g[alt == "O", est])
logf("S2_FIT_DONE")
