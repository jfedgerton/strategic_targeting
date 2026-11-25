## R/explore_specs.R -- exploration of alternative specifications for the main cell (PI 2026-11-24; test code, run on the laptop).
## Every spec is refit on the SAME Bayesian-bootstrap draws as the main table (boot_weights, seed 123), so intervals are comparable.
##   Stage 1   S1NF  between targets (main)        | S1FE  target FE (within)       | S1MUND Mundlak (gaps + target-mean gaps)
##             S1TE  target x era FE               | every S1 spec also reports coalition - bilateral differences (GPT 1)
##   Stage 2   S2FULL multinomial, all disputes (main for GPT 3) | S2JS same on the joint-model disputes (matched cohort)
##             S2SEP  own-side and opposing-side ties together (GPT 2)
##             S2BIN  per-option binary logit, no FE (reference for the FE versions) | S2CFE candidate FE | S2DFE dispute FE
##             S2MON  monthly competing-risk multinomial (GPT 5; clustered SEs, bootstrap only with --monthly_boot=1)
##   Joint     JMAIN  stage-1 target FE -> bbar -> stage 2 (main) | JSEP joint with own + opposing ties (GPT 4: theta as dependence)
## Output (--out): explore_long.csv (every quantity), explore_summary.csv (each variant vs its reference: robust / attenuated / flips / new)
## Usage: Rscript R/explore_specs.R --build=DIR --dump=FILE --positions=FILE --spec=ID --edges_dir=DIR --out=DIR [--n_boot=200] [--cores=8] [--monthly_boot=0]
suppressPackageStartupMessages({ library(data.table); library(parallel); library(survival) })
for (f in c("s2.R", "s2_fit.R", "sry.R", "s2_cell.R")) source(file.path("R/lib", f))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
BUILD <- opt("build"); DUMPF <- opt("dump"); POSF <- opt("positions"); SPEC <- opt("spec"); EDIR <- opt("edges_dir"); OUT <- opt("out")
NB <- as.integer(opt("n_boot", "200")); CORES <- as.integer(opt("cores", "8")); MONB <- opt("monthly_boot", "0") == "1"
stopifnot(!is.null(BUILD), !is.null(DUMPF), !is.null(POSF), !is.null(SPEC), !is.null(EDIR), !is.null(OUT))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); set.seed(123); t0 <- Sys.time()

## ---- 1. inputs (as main_table.R) ----
B <- s2_base(BUILD); W_REF <- B$window; dump <- readRDS(DUMPF); ROLES <- dump$roles
pz <- readRDS(POSF); key <- as.data.table(pz$key); stopifnot(SPEC %in% colnames(pz$R)); R <- pz$R[, SPEC]
E <- s2_spec_edges("core", SPEC, EDIR); X <- s2_prepare(B, key, R, E)
D <- s1_prepare(dump, key, R); ref_order <- names(sort(table(D[chosen == 1L, c2]), decreasing = TRUE))

## ---- 2. extra columns for the variants ----
## Mundlak: each candidate target's mean material / relational gap over all its stage-1 rows, by initiator type
D[, `:=`(mdM = mean(dM), mdR = mean(dR)), by = c2]
for (g in 0:1) { D[, paste0("mdM_", ROLES[g + 1]) := mdM * (type == g)]; D[, paste0("mdR_", ROLES[g + 1]) := mdR * (type == g)] }
MUND <- as.vector(outer(c("mdM", "mdR"), ROLES, paste, sep = "_"))
## eras for the target x era FE
D[, era := cut(yr, c(-Inf, 1913, 1945, 1989, Inf), labels = c("e1816", "e1914", "e1946", "e1990"))]
## stage-2 designs: separate own-side and opposing-side ties
DSEP <- list(A = c(S2_REF_DESIGN$A, "tie_sum_B"), B = c(S2_REF_DESIGN$B, "tie_sum_A"), O = c(S2_REF_DESIGN$O, "tie_sum_B"))
TIE_OWN <- c(A = "tie_sum_A", B = "tie_sum_B", O = "tie_sum_A"); TIE_OPP <- c(A = "tie_sum_B", B = "tie_sum_A", O = "tie_sum_B")
## joint-model disputes from the point fit (fixed across draws, as main_table column 2)
REP0 <- { n <- table(D[chosen == 1L, c2]); names(n)[n >= 2] }

## ---- 3. helpers ----
s1c <- function(b, tag) {   # stage-1 contrasts at gap +-1 per role, plus coalition - bilateral
  v <- c(); for (r in ROLES) for (ax in c("M", "R")) for (sg in c(1, -1)) v[paste(tag, r, ax, sg, sep = ".")] <- sg * b[[paste0("d", ax, "_", r)]] + b[[paste0("d", ax, "2_", r)]]
  for (ax in c("M", "R")) for (sg in c(1, -1)) v[paste(tag, "diff", ax, sg, sep = ".")] <- v[[paste(tag, ROLES[2], ax, sg, sep = ".")]] - v[[paste(tag, ROLES[1], ax, sg, sep = ".")]]
  v }
s2c <- function(b, tag, ties = TRUE) {   # stage-2 surface contrasts (A, O fight the target side: sA_; B fights the attackers: sB_), ties, theta
  v <- c(); for (al in c("A", "B", "O")) { pre <- SURF_TERMS[[al]]
    for (ax in c("M", "R")) for (sg in c(1, -1)) { t1 <- paste0(al, ":", pre, "d", ax); t2 <- paste0(al, ":", pre, "d", ax, "2")
      v[paste(tag, al, ax, sg, sep = ".")] <- if (all(c(t1, t2) %in% names(b))) sg * b[[t1]] + b[[t2]] else NA_real_ }
    if (ties) { v[paste(tag, al, "tie_own", sep = ".")] <- b[paste0(al, ":", TIE_OWN[[al]])]; v[paste(tag, al, "tie_opp", sep = ".")] <- b[paste0(al, ":", TIE_OPP[[al]])] }
    v[paste(tag, al, "theta", sep = ".")] <- b[paste0(al, ":bbar")] }
  v[!is.na(names(v))] }
blk <- function(expr) tryCatch(expr, error = function(e) NULL)   # a failing spec gives NA for its quantities only
## per-option binary conditional logit (Efron approximation) with strata = candidate or dispute; terms constant within strata are dropped
bin_fit <- function(Xs, w, al, strat = NULL) {
  av <- mnl_avail(Xs, al)[, 1]; d <- Xs[av]; d[, `:=`(yy = as.integer(y == S2_CATS_FIT[[al]]), ww = w[av])]
  tm <- S2_REF_DESIGN[[al]]
  if (!is.null(strat)) { d[, st := get(strat)]; tm <- tm[vapply(tm, function(v) d[, sd(get(v)), by = st][, any(V1 > 1e-10, na.rm = TRUE)], TRUE)]
    f <- as.formula(paste("Surv(rep(1, nrow(d)), yy) ~", paste(tm, collapse = " + "), "+ strata(st)"))
    m <- coxph(f, data = d, weights = ww, method = "efron") } else {
    tm <- tm[vapply(tm, function(v) sd(d[[v]]) > 1e-10, TRUE)]
    m <- suppressWarnings(glm(as.formula(paste("yy ~", paste(tm, collapse = " + "))), data = d, weights = ww, family = quasibinomial())) }
  b <- coef(m); b <- b[is.finite(b)]; setNames(b, paste0(al, ":", names(b))) }
bin_all <- function(Xs, w, strat, tag) { b <- unlist(lapply(c("A", "B", "O"), function(al) blk(bin_fit(Xs, w, al, strat)))); if (is.null(b)) return(NULL); s2c(b, tag) }
## stage-1 target FE -> bbar per dispute (as joint_fe_draw); returns list(contrasts, bbar table)
s1_bbar <- function(Db) {
  s1 <- s1_fit_fe(Db, dump, ref_order); bb <- coef(s1$m); be <- bb[grepl("^tfe", names(bb))]; names(be) <- sub("^tfe", "", names(be))
  refb <- setdiff(intersect(ref_order, as.character(unique(Db[chosen == 1L, c2]))), names(be))[1]; be <- c(setNames(0, refb), be)
  b <- be[names(be) %in% REP0 & is.finite(be)]; b <- pmin(pmax(b - mean(b), -5), 5)
  bbt <- B$TS[, .(dispnum, b = b[as.character(ccode)])][, .(bbar = mean(b, na.rm = TRUE)), by = dispnum][is.finite(bbar)]
  list(c = s1c(bb[!grepl("^tfe", names(bb))], "S1FE"), bbt = bbt) }

## ---- 4. one draw: every spec on dispute weights m (m = 1 gives the point estimates) ----
one_draw <- function(m, monthly = FALSE) {
  Db <- copy(D)[, w := w * m[as.character(disp)]]; wX <- as.numeric(m[as.character(X$dispnum)])
  v <- list()
  ## stage 1
  v$s1nf <- blk({ mm <- s1_fit_nofe(Db, dump); s1c(coef(mm), "S1NF") })
  v$s1mu <- blk({ mm <- s1_fit_nofe(Db, dump, MUND); b <- coef(mm); c(s1c(b, "S1MUND"), setNames(b[MUND], paste0("S1MUND.between.", MUND))) })
  v$s1te <- blk({ De <- copy(Db)[, c2 := paste(c2, era)]; ro <- names(sort(table(De[chosen == 1L, c2]), decreasing = TRUE))
                  s1 <- s1_fit_fe(De, dump, ro); bb <- coef(s1$m); s1c(bb[!grepl("^tfe", names(bb))], "S1TE") })
  J <- blk(s1_bbar(Db)); v$s1fe <- J$c
  ## stage 2 alone
  JS <- if (is.null(J)) integer(0) else J$bbt$dispnum
  v$full <- blk(s2c(fit_mnl(X, S2_REF_DESIGN, "full", vc = FALSE, w = wX)$b, "S2FULL"))
  v$js   <- blk({ k <- X$dispnum %in% JS0; s2c(fit_mnl(X[k], S2_REF_DESIGN, "js", vc = FALSE, w = wX[k])$b, "S2JS") })
  v$sep  <- blk(s2c(fit_mnl(X, DSEP, "sep", vc = FALSE, w = wX)$b, "S2SEP"))
  v$bin  <- blk(bin_all(X, wX, NULL, "S2BIN")); v$cfe <- blk(bin_all(X, wX, "ccode", "S2CFE")); v$dfe <- blk(bin_all(X, wX, "dispnum", "S2DFE"))
  ## joint
  if (!is.null(J)) { X1 <- merge(X, J$bbt, by = "dispnum"); w1 <- as.numeric(m[as.character(X1$dispnum)])
    v$jm <- blk(s2c(fit_mnl(X1, lapply(S2_REF_DESIGN, c, "bbar"), "jm", vc = FALSE, w = w1)$b, "JMAIN"))
    v$jsep <- blk(s2c(fit_mnl(X1, lapply(DSEP, c, "bbar"), "jsep", vc = FALSE, w = w1)$b, "JSEP")) }
  if (monthly) v$mon <- blk(s2c(fit_mon(wmon = as.numeric(m[as.character(MX$dispnum)]), vc = FALSE)$b, "S2MON"))
  unlist(unname(v)) }

## ---- 5. monthly competing-risk data: monthly rows of the reference cohort, onset covariates from X, ties by month, piecewise duration ----
MX <- as.data.table(readRDS(file.path(BUILD, sprintf("s2_monthly_w%d.rds", W_REF))))
onset_cols <- setdiff(unique(c(unlist(S2_REF_DESIGN), "likely")), c("tie_sum_A", "tie_sum_B"))
MX <- merge(MX[, .(dispnum, ccode, m, y_m, avail_A, avail_B, avail_orig, tie_sum_A, tie_sum_B)], X[, c("dispnum", "ccode", onset_cols), with = FALSE], by = c("dispnum", "ccode"))
MX[, `:=`(y = y_m, H = Inf, dur_1_2 = as.numeric(m %in% 1:2), dur_3_5 = as.numeric(m %in% 3:5), dur_6_11 = as.numeric(m %in% 6:11), dur_12p = as.numeric(m >= 12))]
DMON <- lapply(S2_REF_DESIGN, c, "dur_1_2", "dur_3_5", "dur_6_11", "dur_12p")
fit_mon <- function(wmon, vc) {   # fit_mnl with the monthly availability masks (restored on exit)
  old <- mnl_avail; on.exit(assign("mnl_avail", old, envir = globalenv()))
  assign("mnl_avail", function(X, alts) cbind(A = X$avail_A, B = X$avail_B, O = X$avail_orig)[, alts, drop = FALSE], envir = globalenv())
  fit_mnl(MX, DMON, "monthly", vc = vc, w = wmon) }

## ---- 6. point estimates ----
JS0 <- s1_bbar(D)$bbt$dispnum
cat("points: FE specs\n"); PT <- one_draw(setNames(rep(1, length(unique(c(D$disp, X$dispnum)))), unique(c(D$disp, X$dispnum))))
## monthly: clustered SEs when the Hessian inverts; otherwise point only (singular Hessian under investigation; flagged in the log, never silent)
cat("points: monthly model\n")
FM <- tryCatch(fit_mon(rep(1, nrow(MX)), vc = TRUE), error = function(e) { cat("WARNING monthly SEs unavailable (", conditionMessage(e), "); point estimates only\n"); NULL })
if (is.null(FM)) FM <- tryCatch(fit_mon(rep(1, nrow(MX)), vc = FALSE), error = function(e) { cat("WARNING monthly model failed:", conditionMessage(e), "\n"); NULL })
PM <- if (is.null(FM)) c() else s2c(FM$b, "S2MON"); SEM <- c()
if (!is.null(FM) && !is.null(FM$V)) for (q in grep("\\.(M|R)\\.", names(PM), value = TRUE)) { p <- strsplit(q, ".", fixed = TRUE)[[1]]; al <- p[2]; ax <- p[3]; sg <- as.numeric(p[4])
  a <- setNames(numeric(length(FM$b)), names(FM$b)); t1 <- paste0(al, ":", SURF_TERMS[[al]], "d", ax); a[t1] <- sg; a[paste0(t1, "2")] <- 1; SEM[q] <- sqrt(drop(t(a) %*% FM$V %*% a)) }
if (!is.null(FM) && !is.null(FM$V)) for (al in c("A", "B", "O")) for (k in c("tie_own", "tie_opp")) { q <- paste("S2MON", al, k, sep = "."); t <- paste0(al, ":", if (k == "tie_own") TIE_OWN[[al]] else TIE_OPP[[al]]); if (q %in% names(PM) && t %in% rownames(FM$V)) SEM[q] <- sqrt(FM$V[t, t]) }
cat(sprintf("points done (%.1f min): %d quantities; monthly rows %d\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), length(PT), nrow(MX)))

## ---- 7. bootstrap draws (the main table's draws) ----
disp_all <- sort(unique(c(as.data.table(dump$D)$disp, B$X$dispnum)))
MM <- boot_weights(disp_all, NB)
DR <- mclapply(MM, function(m) tryCatch(one_draw(m, MONB), error = function(e) NULL), mc.cores = CORES, mc.preschedule = FALSE)
nm <- unique(c(names(PT), if (MONB) names(PM)))
DM <- do.call(rbind, lapply(DR, function(x) if (is.numeric(x)) x[nm] else rep(NA_real_, length(nm)))); colnames(DM) <- nm
cat(sprintf("draws done (%.1f min): %d of %d returned\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), sum(vapply(DR, is.numeric, TRUE)), NB))

## ---- 8. table: percentile intervals (>= 80% finite draws); monthly without bootstrap uses its two-way clustered SE ----
L <- rbindlist(lapply(c(names(PT), names(PM)), function(q) {
  est <- if (q %in% names(PT)) PT[[q]] else PM[[q]]; x <- if (q %in% colnames(DM)) DM[, q] else numeric(0); x <- x[is.finite(x)]
  if (length(x) >= 0.8 * NB) { lo <- quantile(x, .025); hi <- quantile(x, .975); src <- "bootstrap" }
  else if (q %in% names(SEM) && is.finite(SEM[[q]])) { lo <- est - 1.96 * SEM[[q]]; hi <- est + 1.96 * SEM[[q]]; src <- "clustered SE" }
  else { lo <- hi <- NA_real_; src <- "none" }
  p <- strsplit(q, ".", fixed = TRUE)[[1]]
  data.table(spec = p[1], q = paste(p[-1], collapse = "."), est = est, lo = unname(lo), hi = unname(hi), n_draws = length(x), ci = src) }))
L <- L[is.finite(est)][, sig := !is.na(lo) & (lo > 0 | hi < 0)]
fwrite(L, file.path(OUT, "explore_long.csv"))

## ---- 9. summary: each variant vs its reference on the shared quantities ----
REF <- c(S1FE = "S1NF", S1MUND = "S1NF", S1TE = "S1NF", S2JS = "S2FULL", S2SEP = "S2FULL", S2MON = "S2FULL", S2CFE = "S2BIN", S2DFE = "S2BIN", JSEP = "JMAIN")
S <- rbindlist(lapply(names(REF), function(s) { a <- L[spec == s]; r <- L[spec == REF[[s]], .(q, est_ref = est, sig_ref = sig)]
  m <- merge(a[, .(spec, q, est, lo, hi, sig)], r, by = "q"); if (!nrow(m)) return(NULL)
  m[, verdict := fifelse(sig_ref & sig & sign(est) == sign(est_ref), "robust", fifelse(sig_ref & !sig & sign(est) == sign(est_ref), "attenuated",
                 fifelse(sig & sign(est) != sign(est_ref), "flips", fifelse(!sig_ref & sig, "new", fifelse(sig_ref, "lost (sign change, n.s.)", "null both")))))]
  m[, ref := REF[[s]]][] }))
fwrite(S, file.path(OUT, "explore_summary.csv"))
print(dcast(S, spec ~ verdict, fun.aggregate = length))
cat(sprintf("EXPLORE_DONE %s (%.1f min)\n", OUT, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
