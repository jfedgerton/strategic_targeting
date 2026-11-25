## R/explore_specs.R -- exploration of alternative specifications for the main cell (PI 2026-11-24; test code, run on the laptop).
## Every spec is refit on the SAME Bayesian-bootstrap draws as the main table (boot_weights, seed 123), so intervals are comparable.
##   Stage 1   S1NF  between targets (main)        | S1FE  target FE (within)       | S1MUND Mundlak (gaps + target-mean gaps)
##             S1TE  target x era FE               | every S1 spec also reports coalition - bilateral differences (GPT 1)
##   Stage 2   S2FULL multinomial, all disputes (main for GPT 3) | S2JS same on the joint-model disputes (matched cohort)
##             S2SEP  own-side ties only (the main design has ties to both sides since 2026-11-24)
##             S2BIN  per-option binary logit, no FE (reference for the FE versions) | S2CFE candidate FE | S2DFE dispute FE
##             S2MON  monthly competing-risk multinomial (GPT 5; clustered SEs, bootstrap only with --monthly_boot=1)
##   Joint     JMAIN  stage-1 target FE -> bbar -> stage 2 (main) | JSEP joint with own-side ties only
## Added (PI 2026-11-24, robustness for this paper):
##             S1RIV  initiator x target (rivalry) FE   (within-era fits dropped, PI 2026-11-24; era enters only as target x era FE)
##             S2TWFE per-option binary logit with candidate + dispute FE together (fixest)
##             JTE    joint with bbar from target x era FE (bbar of a dispute = its targets' FE in the dispute's onset era)
##             S1ME_* / S2ME_* measurement error in the relational gap, lambda 0.5 and 1: noise (classical) and conf (relational gap
##                    contaminated by the material gap, mimicking power-driven ties); compared with S1NF / S2FULL
##             explore_permutation.csv: relational gaps permuted within each stage-1 event / stage-2 dispute (point fits, --n_perm)
## Output (--out): explore_long.csv (every quantity), explore_summary.csv (each variant vs its reference: robust / attenuated / flips / new)
## Usage: Rscript R/explore_specs.R --build=DIR --dump=FILE --positions=FILE --spec=ID --edges_dir=DIR --out=DIR [--n_boot=200] [--cores=8] [--monthly_boot=0] [--n_perm=200]
suppressPackageStartupMessages({ library(data.table); library(parallel); library(survival); library(fixest) })
for (f in c("s2.R", "s2_fit.R", "sry.R", "s2_cell.R")) source(file.path("R/lib", f))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
BUILD <- opt("build"); DUMPF <- opt("dump"); POSF <- opt("positions"); SPEC <- opt("spec"); EDIR <- opt("edges_dir"); OUT <- opt("out")
NB <- as.integer(opt("n_boot", "200")); CORES <- as.integer(opt("cores", "8")); MONB <- opt("monthly_boot", "0") == "1"; NPERM <- as.integer(opt("n_perm", "200"))
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
ERA_BR <- c(-Inf, 1913, 1945, 1989, Inf); ERA_LB <- c("e1816", "e1914", "e1946", "e1990")
D[, era := as.character(cut(yr, ERA_BR, labels = ERA_LB))]
stopifnot("onset_year" %in% names(X))
DERA <- unique(X[, .(dispnum, era = as.character(cut(onset_year, ERA_BR, labels = ERA_LB)))]); stopifnot(!anyDuplicated(DERA$dispnum), !anyNA(DERA$era))
TSE <- merge(B$TS, DERA, by = "dispnum")[, key := paste(ccode, era)]   # onset targets keyed by target x onset era
REP_TE <- { n <- table(D[chosen == 1L, paste(c2, era)]); names(n)[n >= 2] }
## rivalry FE: dyadic controls that never change within an initiator-target pair are absorbed by the pair FE
DUMP_RIV <- dump; DUMP_RIV$Z_cols <- setdiff(dump$Z_cols, c("contig_land", "log_capdist"))
## measurement error: fixed noise vectors (seed 123, identical in every draw); conf = relational gap contaminated by the material gap
set.seed(123); NZ1 <- rnorm(nrow(D)); NZ2 <- rnorm(nrow(X)); NZ3 <- rnorm(nrow(X)); ME_L <- c(0.5, 1)
stopifnot(all(c("dR_kA", "dR_kB", "dM_kA", "dM_kB") %in% names(X)))
SD1 <- c(R = sd(D$dR), M = sd(D$dM)); SD2 <- c(RA = sd(X$dR_kA), RB = sd(X$dR_kB), MA = sd(X$dM_kA), MB = sd(X$dM_kB))
## stage-2 designs: separate own-side and opposing-side ties
DSEP <- S2_REF_DESIGN_OWN   # main design now has ties to both sides (PI 2026-11-24); S2SEP / JSEP are now the own-side-only versions
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
  if (identical(strat, "twfe")) {   # candidate + dispute FE together; fixest drops collinear terms and all-0 / all-1 FE groups
    m <- suppressMessages(suppressWarnings(feglm(as.formula(paste("yy ~", paste(tm, collapse = " + "), "| ccode + dispnum")), data = d, weights = ~ww, family = binomial(), notes = FALSE)))
    b <- coef(m); b <- b[is.finite(b)]; return(setNames(b, paste0(al, ":", names(b)))) }
  if (!is.null(strat)) { d[, st := get(strat)]; tm <- tm[vapply(tm, function(v) d[, sd(get(v)), by = st][, any(V1 > 1e-10, na.rm = TRUE)], TRUE)]
    f <- as.formula(paste("Surv(rep(1, nrow(d)), yy) ~", paste(tm, collapse = " + "), "+ strata(st)"))
    m <- coxph(f, data = d, weights = ww, method = "efron") } else {
    tm <- tm[vapply(tm, function(v) sd(d[[v]]) > 1e-10, TRUE)]
    m <- suppressWarnings(glm(as.formula(paste("yy ~", paste(tm, collapse = " + "))), data = d, weights = ww, family = quasibinomial())) }
  b <- coef(m); b <- b[is.finite(b)]; setNames(b, paste0(al, ":", names(b))) }
bin_all <- function(Xs, w, strat, tag) { b <- unlist(lapply(c("A", "B", "O"), function(al) blk(bin_fit(Xs, w, al, strat)))); if (is.null(b)) return(NULL); s2c(b, tag) }
## stage-1 target FE -> bbar per dispute (as joint_fe_draw); returns list(contrasts, bbar table)
## finer FE (target x era, initiator x target): a level that appears only when chosen has an infinite FE, so its events are dropped
## (iterated until none is left); the number of events dropped is returned as attribute "dropped"
fe_prune <- function(Dx) { n0 <- uniqueN(Dx$ev_id)
  repeat { tgb <- unique(Dx[chosen == 1L, c2]); dd <- Dx[c2 %in% tgb]; dd <- dd[dd[, .(k = .N), by = ev_id][k > 1L], on = "ev_id"]
    bad <- setdiff(unique(dd[chosen == 1L, c2]), dd[chosen == 0L, unique(c2)])
    if (!length(bad)) break; Dx <- Dx[!ev_id %in% Dx[chosen == 1L & c2 %in% bad, ev_id]] }
  setattr(Dx, "dropped", n0 - uniqueN(Dx$ev_id)); Dx }
## finer-FE fit: like s1_fit_fe, but FE levels that are aliased (NA) or not converged are tolerated; the gap and control coefficients must
## be finite and below 15 in absolute value (else the fit is invalid, as in s1_fit_fe). Aliased / diverging levels get no FE (no bbar).
s1_fit_fe_tol <- function(Dx, dmp, ro) {
  tgb <- unique(Dx[chosen == 1L, c2]); dd <- Dx[c2 %in% tgb]; dd <- dd[dd[, .(k = .N), by = ev_id][k > 1L], on = "ev_id"]
  dd[, tfe := droplevels(factor(as.character(c2), levels = intersect(ro, as.character(tgb))))]
  f <- as.formula(paste("Surv(rep(1, nrow(dd)), chosen) ~", paste(c(s1_terms(dmp), dmp$Z_cols, "tfe"), collapse = " + "), "+ strata(ev_id)"))
  m <- suppressWarnings(coxph(f, data = dd, weights = w, method = "breslow", control = coxph.control(iter.max = 50)))
  b <- coef(m); g <- b[!grepl("^tfe", names(b))]
  if (!all(is.finite(g)) || any(abs(g) > 15)) stop("finer-FE fit invalid: gap / control coefficients not finite")
  b[grepl("^tfe", names(b)) & (!is.finite(b) | abs(b) > 15)] <- NA; list(b = b, dd = dd) }
## te = TRUE: target x era FE; a dispute's bbar averages its onset targets' FE in the dispute's onset era (repeat target-eras only)
s1_bbar <- function(Db, te = FALSE) {
  if (te) Db <- fe_prune(copy(Db)[, c2 := paste(c2, era)])
  ro <- if (te) names(sort(table(Db[chosen == 1L, c2]), decreasing = TRUE)) else ref_order
  bb <- if (te) s1_fit_fe_tol(Db, dump, ro)$b else coef(s1_fit_fe(Db, dump, ro)$m); be <- bb[grepl("^tfe", names(bb))]; names(be) <- sub("^tfe", "", names(be))
  refb <- setdiff(intersect(ro, as.character(unique(Db[chosen == 1L, c2]))), names(be))[1]; be <- c(setNames(0, refb), be)
  b <- be[names(be) %in% (if (te) REP_TE else REP0) & is.finite(be)]; b <- pmin(pmax(b - mean(b), -5), 5)
  bbt <- if (te) TSE[, .(dispnum, b = b[key])][, .(bbar = mean(b, na.rm = TRUE)), by = dispnum][is.finite(bbar)]
         else B$TS[, .(dispnum, b = b[as.character(ccode)])][, .(bbar = mean(b, na.rm = TRUE)), by = dispnum][is.finite(bbar)]
  list(c = s1c(bb[!grepl("^tfe", names(bb))], if (te) "S1TE" else "S1FE"), bbt = bbt) }
## stage-1 relational terms rebuilt after dR changes (measurement error, permutation)
s1_redo <- function(Dx) { Dx[, `:=`(dR2 = dR^2, dMR = dM * dR)]
  for (g in 0:1) for (v in c("dR", "dR2", "dMR")) Dx[, paste0(v, "_", ROLES[g + 1]) := get(v) * (type == g)]; Dx }
me_D <- function(Db, kind, lam) { Dx <- copy(Db)
  if (kind == "noise") Dx[, dR := dR + lam * SD1[["R"]] * NZ1] else Dx[, dR := dR + lam * SD1[["R"]] / SD1[["M"]] * dM]; s1_redo(Dx) }
me_X <- function(kind, lam) { Xx <- copy(X)
  if (kind == "noise") Xx[, `:=`(dR_kA = dR_kA + lam * SD2[["RA"]] * NZ2, dR_kB = dR_kB + lam * SD2[["RB"]] * NZ3)]
  else Xx[, `:=`(dR_kA = dR_kA + lam * SD2[["RA"]] / SD2[["MA"]] * dM_kA, dR_kB = dR_kB + lam * SD2[["RB"]] / SD2[["MB"]] * dM_kB)]
  add_surface(Xx) }
lamtag <- function(l) sprintf("%02d", round(10 * l))   # 0.5 -> "05", 1 -> "10" (spec names carry no dots)

## ---- 4. one draw: every spec on dispute weights m (m = 1 gives the point estimates) ----
one_draw <- function(m, monthly = FALSE) {
  Db <- copy(D)[, w := w * m[as.character(disp)]]; wX <- as.numeric(m[as.character(X$dispnum)])
  v <- list()
  ## stage 1
  v$s1nf <- blk({ mm <- s1_fit_nofe(Db, dump); s1c(coef(mm), "S1NF") })
  v$s1mu <- blk({ mm <- s1_fit_nofe(Db, dump, MUND); b <- coef(mm); c(s1c(b, "S1MUND"), setNames(b[MUND], paste0("S1MUND.between.", MUND))) })
  JT <- blk(s1_bbar(Db, te = TRUE)); v$s1te <- JT$c
  J <- blk(s1_bbar(Db)); v$s1fe <- J$c
  v$s1riv <- blk({ Dr <- fe_prune(copy(Db)[, c2 := paste(c1, c2, sep = "_")]); ro <- names(sort(table(Dr[chosen == 1L, c2]), decreasing = TRUE))
                   bb <- s1_fit_fe_tol(Dr, DUMP_RIV, ro)$b; s1c(bb[!grepl("^tfe", names(bb))], "S1RIV") })
  for (k in c("noise", "conf")) for (l in ME_L) {
    v[[paste0("me1", k, l)]] <- blk({ mm <- s1_fit_nofe(me_D(Db, k, l), dump); s1c(coef(mm), paste0("S1ME_", k, lamtag(l))) })
    v[[paste0("me2", k, l)]] <- blk(s2c(fit_mnl(XME[[paste(k, l)]], S2_REF_DESIGN, "me", vc = FALSE, w = wX)$b, paste0("S2ME_", k, lamtag(l)))) }
  ## stage 2 alone
  JS <- if (is.null(J)) integer(0) else J$bbt$dispnum
  v$full <- blk(s2c(fit_mnl(X, S2_REF_DESIGN, "full", vc = FALSE, w = wX)$b, "S2FULL"))
  v$js   <- blk({ k <- X$dispnum %in% JS0; s2c(fit_mnl(X[k], S2_REF_DESIGN, "js", vc = FALSE, w = wX[k])$b, "S2JS") })
  v$sep  <- blk(s2c(fit_mnl(X, DSEP, "sep", vc = FALSE, w = wX)$b, "S2SEP"))
  v$bin  <- blk(bin_all(X, wX, NULL, "S2BIN")); v$cfe <- blk(bin_all(X, wX, "ccode", "S2CFE")); v$dfe <- blk(bin_all(X, wX, "dispnum", "S2DFE"))
  v$twfe <- blk(bin_all(X, wX, "twfe", "S2TWFE"))
  ## joint
  if (!is.null(J)) { X1 <- merge(X, J$bbt, by = "dispnum"); w1 <- as.numeric(m[as.character(X1$dispnum)])
    v$jm <- blk(s2c(fit_mnl(X1, lapply(S2_REF_DESIGN, c, "bbar"), "jm", vc = FALSE, w = w1)$b, "JMAIN"))
    v$jsep <- blk(s2c(fit_mnl(X1, lapply(DSEP, c, "bbar"), "jsep", vc = FALSE, w = w1)$b, "JSEP")) }
  if (!is.null(JT)) { X2 <- merge(X, JT$bbt, by = "dispnum"); w2 <- as.numeric(m[as.character(X2$dispnum)])
    v$jte <- blk(s2c(fit_mnl(X2, lapply(S2_REF_DESIGN, c, "bbar"), "jte", vc = FALSE, w = w2)$b, "JTE")) }
  if (monthly) v$mon <- blk(s2c(fit_mon(wmon = as.numeric(m[as.character(MX$dispnum)]), vc = FALSE)$b, "S2MON"))
  unlist(unname(v)) }

XME <- list(); for (k in c("noise", "conf")) for (l in ME_L) XME[[paste(k, l)]] <- me_X(k, l)   # contaminated stage-2 data (fixed across draws)

## ---- 5. monthly competing-risk data: monthly rows of the reference cohort, onset covariates from X, ties by month, piecewise duration ----
MX <- as.data.table(readRDS(file.path(BUILD, sprintf("s2_monthly_w%d.rds", W_REF))))
onset_cols <- setdiff(unique(c(unlist(S2_REF_DESIGN), "likely")), c("tie_sum_A", "tie_sum_B"))
MX <- merge(MX[, .(dispnum, ccode, m, y_m, avail_A, avail_B, avail_orig, tie_sum_A, tie_sum_B)], X[, c("dispnum", "ccode", onset_cols), with = FALSE], by = c("dispnum", "ccode"))
MX[, `:=`(y = y_m, H = Inf, dur_1_2 = as.numeric(m %in% 1:2), dur_3_5 = as.numeric(m %in% 3:5), dur_6_11 = as.numeric(m %in% 6:11), dur_12p = as.numeric(m >= 12))]
DMON <- lapply(S2_REF_DESIGN, c, "dur_1_2", "dur_3_5", "dur_6_11", "dur_12p")
DMON$A <- setdiff(DMON$A, "dur_1_2")   # joining the attackers is available only after the originator window, so its four duration dummies sum to 1 (collinear with the intercept); months 1-2 are its reference
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

## ---- 7b. permutation (point fits): relational gaps shuffled across candidates within each stage-1 event and within each stage-2
## dispute (both gaps of a candidate move together). Seed 123 + i per permutation, so results do not depend on --cores.
## p_perm is two-sided against the permutation distribution itself (not against 0): stage-2 permuted contrasts are not centred on 0.
perm_one <- function(i) { set.seed(123 + i)
  Dp <- copy(D); Dp[, dR := dR[sample.int(.N)], by = ev_id]; Dp <- s1_redo(Dp)
  Xp <- copy(X); Xp[, c("dR_kA", "dR_kB") := { o <- sample.int(.N); list(dR_kA[o], dR_kB[o]) }, by = dispnum]; Xp <- add_surface(Xp)
  c(blk(s1c(coef(s1_fit_nofe(Dp, dump)), "S1NF")), blk(s2c(fit_mnl(Xp, S2_REF_DESIGN, "perm", vc = FALSE)$b, "S2FULL"))) }
if (NPERM > 0) {
  PR <- mclapply(seq_len(NPERM), function(i) tryCatch(perm_one(i), error = function(e) NULL), mc.cores = CORES, mc.preschedule = FALSE)
  qp <- grep("^(S1NF|S2FULL)\\.", names(PT), value = TRUE)
  PMX <- do.call(rbind, lapply(PR, function(x) if (is.numeric(x)) x[qp] else rep(NA_real_, length(qp)))); colnames(PMX) <- qp
  PRM <- rbindlist(lapply(qp, function(q) { x <- PMX[, q]; x <- x[is.finite(x)]; p <- strsplit(q, ".", fixed = TRUE)[[1]]
    data.table(spec = p[1], q = paste(p[-1], collapse = "."), est = PT[[q]], perm_mean = mean(x), perm_sd = sd(x), p_perm = min(1, 2 * min((1 + sum(x >= PT[[q]])) / (1 + length(x)), (1 + sum(x <= PT[[q]])) / (1 + length(x)))), n_perm = length(x)) }))
  fwrite(PRM, file.path(OUT, "explore_permutation.csv"))
  cat(sprintf("permutations done (%.1f min): %d of %d returned\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), sum(vapply(PR, is.numeric, TRUE)), NPERM)) }

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
REF <- c(S1FE = "S1NF", S1MUND = "S1NF", S1TE = "S1NF", S1RIV = "S1NF", S2JS = "S2FULL", S2SEP = "S2FULL", S2MON = "S2FULL",
         S2CFE = "S2BIN", S2DFE = "S2BIN", S2TWFE = "S2BIN", JSEP = "JMAIN", JTE = "JMAIN",
         setNames(rep("S1NF", 4), paste0("S1ME_", rep(c("noise", "conf"), each = 2), lamtag(ME_L))),
         setNames(rep("S2FULL", 4), paste0("S2ME_", rep(c("noise", "conf"), each = 2), lamtag(ME_L))))
S <- rbindlist(lapply(names(REF), function(s) { a <- L[spec == s]; r <- L[spec == REF[[s]], .(q, est_ref = est, sig_ref = sig)]
  m <- merge(a[, .(spec, q, est, lo, hi, sig)], r, by = "q"); if (!nrow(m)) return(NULL)
  m[, verdict := fifelse(sig_ref & sig & sign(est) == sign(est_ref), "robust", fifelse(sig_ref & !sig & sign(est) == sign(est_ref), "attenuated",
                 fifelse(sig & sign(est) != sign(est_ref), "flips", fifelse(!sig_ref & sig, "new", fifelse(sig_ref, "lost (sign change, n.s.)", "null both")))))]
  m[, ref := REF[[s]]][] }))
fwrite(S, file.path(OUT, "explore_summary.csv"))
print(dcast(S, spec ~ verdict, fun.aggregate = length))
cat(sprintf("EXPLORE_DONE %s (%.1f min)\n", OUT, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
