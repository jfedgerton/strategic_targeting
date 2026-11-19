## ============================================================================
## 22_rel_sweep.R  --  Part 2d of the relational-measure sweep (PI design 2026-11-18)
## For every candidate relational position (21_rel_measures.R) and every type scheme (long frames dumped by
## 10_clogit_roles.R --dump_long=TRUE), rebuild the relational gap exactly as 10_clogit_roles.R does and refit
##   M0  pooled conditional logit   (gap terms dM, dM2, dR, dR2, dMR + Z)        -> is there ANY relational effect?
##   M1  type-interacted            (gap terms by type + Z)                       -> does it differ by type?
## with the same strata, D8 weights, Z and two-way CGM (initiator x dispute) variance as 10_clogit_roles.R section 4.
##
## Gap: dR = (R_initiator,t-1 - R_target,t-1) / s, as in 01_build_data.R (positive = target LESS central).
##   --rel_scale=year in the dump: s = gap_scale(R over all states in the EVENT year t) (10_clogit_roles.R section 3 convention);
##   pooled: s = gap_scale(R, year) over all state-years.
## Validation (always run first): with R = log PageRank from the dump's sy_positions.rds, the rebuilt dR must equal the dumped dR
## (max abs diff < 1e-8; the script stops otherwise). Because the gap is identical, the BASELINE fit is the fit on the dumped gap
## (checked once by hand against 10_clogit_roles.R output: coefficients and two-way SEs identical).
## A spec is refused (row model = SKIPPED, with the reason) when its position is constant in any year this fit uses (t-1 for positions;
## also t under the year-t scale): the year scale would be ~0 and the gap ~1e13 (independent review 2026-11-18, P1-1; GPT P2-5).
## A fit with separation / non-convergence warnings, non-finite coefficients or the iteration limit is model = INVALID: it stays in the
## census with its iterations and max |dR| but has no estimate rows (GPT review P1-3). Each sweep log line is the census for that file.
## Validity is per model (full audit 2026-11-19 P2-11): an invalid M1 leaves M0 rows in place. M0 pools type-specific controls (P2-6); M0t
## (pooled gaps, type-specific controls) is added when such controls exist. The scheme tag carries the dump's bilateral rule (_bon / _ever).
## Per spec the output also records coxph iterations, any coxph warning, and max |dR| (review P2-7).
##
## Output (one row per spec x scheme x quantity): <out_dir>/sweep_<scheme tag>_<positions file stem>.csv with
##   spec_id, scheme, model, quantity, est, se, lo, hi, p, df   (95% intervals, two-way CGM)
##   M0: b_dR, b_dR2, b_dMR, u_R_less (dR=+1: b_dR + b_dR2), u_R_more (dR=-1: -b_dR + b_dR2), the same for M, axis_R (3 df), axis_M (3 df)
##   M1: per type the same; per pair of types the differences (diff_*), diff_axis_R, diff_axis_M (3 df), diff_surface (5 df)
## Usage: Rscript v2/R/22_rel_sweep.R <long_dump.rds> <positions_*.rds> <out_dir> [--yscale=t|tm1] [--cores=N] [--max_specs=K]
## ============================================================================
suppressMessages({ library(survival); library(parallel) })
set.seed(123)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
DUMP <- args[1]; POSF <- args[2]; OUT_DIR <- args[3]
opt <- function(name, default) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
CORES <- as.integer(opt("cores", "1")); MAX_SPECS <- as.integer(opt("max_specs", "0"))
YSCALE <- opt("yscale", "t"); stopifnot(YSCALE %in% c("t", "tm1"))
POS_LAG <- as.integer(opt("pos_lag", "1")); stopifnot(POS_LAG >= -5L, POS_LAG <= 20L)     # PI 2026-11-18: positions at t-L (punishment-vs-selection check); L = 0 is the event year, L < 0 a lead (positions AFTER the event year: t+|L|), 2026-11-18 test 2
DECOMP <- as.integer(opt("decomp", "0")); stopifnot(DECOMP %in% c(0L, 5L, 10L))
if (DECOMP > 0L && POS_LAG < 1L) stop("--decomp needs positions at t-1 or earlier (pos_lag >= 1)")   # >0: dR = long-run baseline gap, plus type-specific recent change
SAMPLE_LAG <- as.integer(opt("sample_lag", "0")); stopifnot(SAMPLE_LAG >= -5L, SAMPLE_LAG <= 20L)   # review P2-5: restrict to the strata / alternatives available at t-SAMPLE_LAG (comparison sample for pos_lag)
## 2026-11-19 additions (GPT full audit P2-7, P2-8; multiverse loop):
##  --scale=year|pooled : year (default) divides the gap by the year scale (event year t at pos_lag 1 under --yscale=t, else the positions' year);
##                        pooled divides by ONE scale over all years, so lag and sample-lag runs differ only in information (P2-7)
##  --decomp_terms=linear|full : full adds type-specific dRc^2, dM x dRc and dR x dRc (the terms the quadratic implies; P2-8), tested jointly
##  --merror=none|def|contig|random --merror_k=K --layers=<layers_dy.rds> : systematic measurement error that mimics a confounder, added to
##                        every position before the gap: R' = R + K * sd_year(R) * z_year(m), m = log(1 + defense partners) | log(1 + land
##                        neighbours) (both from the layers panel, so --layers is required) | N(0, 1) (seed 123)
##  --stepB=TRUE : failed-bargain step B on a mid dump -- weighted logit of ev_fatal on the CHOSEN row's type-specific gaps, type intercepts and Z
SCALE <- opt("scale", "year"); stopifnot(SCALE %in% c("year", "pooled"))
DECOMP_TERMS <- opt("decomp_terms", "linear"); stopifnot(DECOMP_TERMS %in% c("linear", "full"))
MERROR <- opt("merror", "none"); stopifnot(MERROR %in% c("none", "def", "contig", "random")); MERROR_K <- as.numeric(opt("merror_k", "0"))
if (MERROR != "none") stopifnot(is.finite(MERROR_K), MERROR_K > 0)
if (MERROR %in% c("def", "contig")) stopifnot(file.exists(opt("layers", "")))
STEPB <- as.logical(opt("stepB", "FALSE")); stopifnot(!is.na(STEPB))
KNOWN <- c("cores", "max_specs", "yscale", "pos_lag", "decomp", "sample_lag", "bilateral", "scale", "decomp_terms", "merror", "merror_k", "layers", "stepB")
fl <- sub("^--([A-Za-z_]+)=.*$", "\\1", grep("^--", args, value = TRUE)); if (any(!fl %in% KNOWN)) stop("unknown option(s): ", paste(setdiff(fl, KNOWN), collapse = ", "))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
logf <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
tlab <- function(L) if (L < 0L) paste0("t+", -L) else paste0("t-", L)   # readable position year for log lines
logf("R %s | survival %s | cores %d", getRversion(), packageVersion("survival"), CORES)

## gap_scale(): copied verbatim from helpers_data.R (pooled within-year SD of pairwise differences)
gap_scale <- function(v, yr) {
  stopifnot(length(v) == length(yr), !anyNA(v), all(is.finite(v)))
  ss <- tapply(v, yr, function(z) sum((z - mean(z))^2)); n <- tapply(v, yr, length)
  out <- sqrt(sum(2 * n * ss) / sum(n * (n - 1)))
  stopifnot(is.finite(out), out > 0)
  out
}

## ---- 0. inputs ----------------------------------------------------------------------------------------------------
dump <- readRDS(DUMP); D <- dump$D; ROLES <- dump$roles; Z_cols <- dump$Z_cols; G <- length(ROLES)
stopifnot(all(c("ev_id", "chosen", "type", "w", "init", "disp", "init_disp", "dM", "dR", "c1", "c2", "yr", Z_cols) %in% names(D)),
          all(D$type %in% 0:(G - 1L)), all(tapply(D$chosen, D$ev_id, sum) == 1L))
## GPT full audit 2026-11-19 P1-1: the bilateral rule must be recorded in the dump, must match --bilateral when given, and is part of the scheme tag
BILAT_ARG <- opt("bilateral", "")
if (is.null(dump$bilateral_rule) || !(dump$bilateral_rule %in% c("onset", "ever"))) stop("dump has no bilateral_rule (rebuild it with 10_clogit_roles.R)")
if (nzchar(BILAT_ARG) && BILAT_ARG != dump$bilateral_rule) stop(sprintf("--bilateral=%s but the dump was built with %s", BILAT_ARG, dump$bilateral_rule))
SCHEME <- paste0(dump$types, if (dump$bilateral_rule == "onset") "_bon" else "_ever", if (SAMPLE_LAG != 0L) paste0("_samp", if (SAMPLE_LAG < 0L) paste0("lead", -SAMPLE_LAG) else SAMPLE_LAG) else "", if (POS_LAG != 1L) (if (POS_LAG < 0L) paste0("_lead", -POS_LAG) else paste0("_lag", POS_LAG)) else "", if (DECOMP > 0L) paste0("_decomp", DECOMP, if (DECOMP_TERMS == "full") "full" else "") else "", if (SCALE == "pooled") "_pooled" else "", if (MERROR != "none") paste0("_me", MERROR, MERROR_K) else "", if (STEPB) "_stepB" else "",
                 if (!is.null(dump$coorig_months) && dump$coorig_months > 0L) paste0("_w", dump$coorig_months) else "", if (dump$rel_scale == "year") paste0("_relyear", if (YSCALE == "tm1") "tm1" else "") else "", if (isTRUE(dump$unit_weights)) "_uw" else "", "_", dump$sample)
if (YSCALE == "tm1" && dump$rel_scale != "year") stop("--yscale=tm1 needs a dump with rel_scale = year")
logf("dump %s: scheme %s, %d rows, %d strata, types %s, rel_scale %s, commit %s", basename(DUMP), SCHEME, nrow(D), length(unique(D$ev_id)),
     paste(ROLES, collapse = "/"), dump$rel_scale, dump$git_commit)
if (STEPB) { stopifnot(dump$sample == "mid", "ev_fatal" %in% names(D), DECOMP == 0L); logf("step B: %d of %d events fatal", sum(D$ev_fatal[D$chosen == 1L]), sum(D$chosen == 1L)) }
if (DECOMP > 0L && SCALE == "pooled") stop("--decomp with --scale=pooled is not supported")
P <- readRDS(POSF); key <- P$key; Rm <- P$R; specs <- P$specs
stopifnot(identical(colnames(Rm), specs$spec_id), nrow(key) == nrow(Rm), !anyDuplicated(key[, c("ccode", "year")]))
kidx <- setNames(seq_len(nrow(key)), paste(key$ccode, key$year))
if (SAMPLE_LAG != 0L) {   # SAMPLE_LAG < 0: restrict to the strata / alternatives that have positions at t+|SAMPLE_LAG| (comparison sample for a lead)
  mi <- is.na(kidx[paste(D$c1, D$yr - SAMPLE_LAG)]); mj <- is.na(kidx[paste(D$c2, D$yr - SAMPLE_LAG)])
  bad <- unique(D$ev_id[(mi | mj) & (D$chosen == 1L | mi)]); keep <- !(mi | mj) & !(D$ev_id %in% bad)
  logf("sample_lag = %d: restricted to the %s sample (dropped %d strata, %d other alternatives); positions at %s", SAMPLE_LAG, tlab(SAMPLE_LAG), length(bad), sum((mi | mj) & !(D$ev_id %in% bad)), tlab(POS_LAG))
  D <- D[keep, ]
}
i_row <- kidx[paste(D$c1, D$yr - POS_LAG)]; j_row <- kidx[paste(D$c2, D$yr - POS_LAG)]
if (POS_LAG != 1L) {   # states / years without a position at t-L: drop those alternatives, and the whole stratum if the chosen target or initiator lacks one
  miss <- is.na(i_row) | is.na(j_row)
  bad_ev <- unique(D$ev_id[miss & (D$chosen == 1L | is.na(i_row))])
  keep <- !miss & !(D$ev_id %in% bad_ev)
  logf("pos_lag = %d: dropped %d strata (chosen target or initiator without a position at %s) and %d other alternatives", POS_LAG, length(bad_ev), tlab(POS_LAG),
       sum(miss & !(D$ev_id %in% bad_ev)))
  w_before <- tapply(D$w[D$chosen == 1L], D$type[D$chosen == 1L], sum)
  D <- D[keep, ]; i_row <- i_row[keep]; j_row <- j_row[keep]
  stopifnot(all(tapply(D$chosen, D$ev_id, sum) == 1L))
  w_after <- tapply(D$w[D$chosen == 1L], D$type[D$chosen == 1L], sum)
  logf("pos_lag = %d: event weight kept by type %s (old weights retained, not renormalised)", POS_LAG, paste(sprintf("%s %.1f of %.1f", ROLES[as.integer(names(w_before)) + 1L], w_after[names(w_before)], w_before), collapse = "; "))
}
if (anyNA(i_row) || anyNA(j_row)) stop(sprintf("positions missing at t-%d for %d initiator / %d target rows", POS_LAG, sum(is.na(i_row)), sum(is.na(j_row))))
yr_ev <- as.character(D$yr)
key_years <- split(seq_len(nrow(key)), key$year)
stopifnot(all(unique(c(yr_ev, as.character(D$yr - POS_LAG))) %in% names(key_years)))

## GPT review P2-5: the years a fit actually needs are the positions' year (t-1) plus, under the year-t scale, year t
## audit P2-7: year t is needed only when it is the scale year (pos_lag 1, --yscale=t, --scale=year)
used_years <- sort(unique(if (dump$rel_scale == "year" && YSCALE == "t" && POS_LAG == 1L && SCALE == "year") c(D$yr, D$yr - POS_LAG) else D$yr - POS_LAG))
if (DECOMP > 0L) used_years <- sort(unique(c(used_years, unlist(lapply(used_years, function(y) (y - DECOMP):(y - 1L))))))
used_years <- used_years[used_years >= min(key$year)]
## DECOMP (audit P2-8): for key row (c, y), the same state's rows in the EARLIER years y-DECOMP..y-1 (a baseline that does not contain the
## positions' year itself); change = centred position at y minus that baseline. A state with no earlier year has baseline = its own year (change 0)
## and is flagged decomp_short.
if (DECOMP > 0L) {
  base_idx <- lapply(seq_len(nrow(key)), function(k) { v <- kidx[paste(key$ccode[k], (key$year[k] - DECOMP):(key$year[k] - 1L))]; v <- v[!is.na(v)]; if (length(v)) v else k })
  ## review P1-2: positions are centred within year before averaging (raw log PageRank drifts over time, so an uncentred baseline of a young
  ## state is mechanically different); a target or initiator with fewer than DECOMP years gets the indicator decomp_short (in Z)
  yr_of_key <- split(seq_len(nrow(key)), key$year)
  centre <- function(R) { out <- R; for (ix in yr_of_key) out[ix] <- R[ix] - mean(R[ix]); out }
  rbase <- function(Rc) vapply(base_idx, function(v) mean(Rc[v]), 0)
  short_key <- vapply(seq_along(base_idx), function(k) length(base_idx[[k]]) < DECOMP || identical(base_idx[[k]], k), TRUE)
  D$decomp_short <- as.numeric(short_key[j_row]) - as.numeric(short_key[i_row])   # target short minus initiator short (initiator constant within stratum)
  Z_cols <- c(Z_cols, "decomp_short")
  logf("decomp = %d: rows whose target has a short baseline window %.3f", DECOMP, mean(short_key[j_row]))
}
## measurement-error confounder per key row (z within year), built once
if (MERROR != "none") {
  if (MERROR == "random") { set.seed(123); mvar <- rnorm(nrow(key)) } else {
    LY <- readRDS(opt("layers", "")); stopifnot(all(c("ccode1", "ccode2", "year") %in% names(LY)))
    e <- if (MERROR == "def") LY[LY$def == 1L, c("ccode1", "ccode2", "year")] else LY[LY$contig_land == 1L, c("ccode1", "ccode2", "year")]   # both from the layers panel (20_build_layers.R)
    e2 <- rbind(data.frame(a = e$ccode1, y = e$year), data.frame(a = e$ccode2, y = e$year))
    cnt <- table(paste(e2$a, e2$y)); mvar <- log1p(as.numeric(cnt[paste(key$ccode, key$year)])); mvar[is.na(mvar)] <- 0
  }
  mz <- mvar; for (ix in key_years) { s_ <- sd(mvar[ix]); mz[ix] <- if (s_ > 0) (mvar[ix] - mean(mvar[ix])) / s_ else 0 }
  logf("merror %s, k = %.2f: confounder z built for %d state-years", MERROR, MERROR_K, length(mz))
}
add_merror <- function(R) { if (MERROR == "none") return(R); out <- R; for (ix in key_years) out[ix] <- R[ix] + MERROR_K * sd(R[ix]) * mz[ix]; out }
wsd <- function(x, w) { m <- sum(w * x) / sum(w); sqrt(sum(w * (x - m)^2) / sum(w)) }
make_dR <- function(R) {                                          # R: position per key row
  R <- add_merror(R)
  sds <- vapply(key_years[as.character(used_years)], function(ix) sd(R[ix]), 0)
  if (any(sds < 1e-8)) stop(sprintf("position constant in event-used year(s) %s", paste(used_years[sds < 1e-8], collapse = " ")))
  if (DECOMP > 0L) { Rc_ <- centre(R); Rb_ <- rbase(Rc_); num <- Rb_[i_row] - Rb_[j_row]; numc <- (Rc_[i_row] - Rb_[i_row]) - (Rc_[j_row] - Rb_[j_row]) } else num <- R[i_row] - R[j_row]
  if (dump$rel_scale == "year" && SCALE == "year") {
    sc_yr <- if (YSCALE == "t" && POS_LAG == 1L) yr_ev else as.character(D$yr - POS_LAG)   # review P2-4: with a lag, scale by the positions' own year
    ev_years <- unique(sc_yr)
    ys <- vapply(key_years[ev_years], function(ix) gap_scale(R[ix], rep(1L, length(ix))), 0)
    stopifnot(all(ys > 1e-8))
    sc <- ys[sc_yr]
  } else sc <- gap_scale(R, key$year)
  if (DECOMP > 0L) { dRc <- numc / sc; list(dR = num / sc, dRc = dRc / wsd(dRc, D$w)) } else num / sc   # dRc per event-weighted SD of the change gap (audit P2-8)
}

## ---- 1. model machinery (as 10_clogit_roles.R section 4) ----------------------------------------------------------
GAP <- c("dM", "dM2", "dR", "dR2", "dMR")
GAP_G <- as.vector(outer(GAP, ROLES, paste, sep = "_"))
fit_cl <- function(d, terms, zc = Z_cols) coxph(as.formula(paste("Surv(rep(1, nrow(d)), chosen) ~", paste(c(terms, zc), collapse = " + "), "+ strata(ev_id)")),
                                   data = d, weights = w, method = "breslow", cluster = init)
vc_cl <- function(m, g) { Dd <- rowsum(residuals(m, type = "dfbeta", weighted = TRUE), g, reorder = FALSE); V <- crossprod(Dd); dimnames(V) <- list(names(coef(m)), names(coef(m))); V }
vc_2w <- function(m, d) {
  V <- vc_cl(m, d$init) + vc_cl(m, d$disp) - vc_cl(m, d$init_disp)
  e <- eigen(V, symmetric = TRUE); nn <- sum(e$values <= 0)
  if (nn > 0) { V <- e$vectors %*% diag(pmax(e$values, 0)) %*% t(e$vectors); dimnames(V) <- list(names(coef(m)), names(coef(m))) }
  attr(V, "n_neg") <- nn; V
}
lin <- function(b, V, a) { a <- a[names(b)]; a[is.na(a)] <- 0; est <- sum(a * b); se <- sqrt(max(drop(t(a) %*% V %*% a), 0)); c(est = est, se = se) }
wald <- function(b, V, R) { d <- R %*% b; S <- R %*% V %*% t(R); s <- tryCatch(as.numeric(t(d) %*% solve(S, d)), error = function(e) NA_real_)
  c(chi2 = s, df = nrow(R), p = if (is.na(s)) NA_real_ else pchisq(s, nrow(R), lower.tail = FALSE)) }
Rsel <- function(nm, terms, minus = NULL) { R <- matrix(0, length(terms), length(nm), dimnames = list(NULL, nm))
  for (i in seq_along(terms)) { R[i, terms[i]] <- 1; if (!is.null(minus)) R[i, minus[i]] <- -1 }; R }
row_out <- function(spec, model, q, est = NA, se = NA, p = NA, df = 1L)
  data.frame(spec_id = spec, scheme = SCHEME, model = model, quantity = q, est = est, se = se, lo = est - 1.96 * se, hi = est + 1.96 * se,
             p = if (is.na(p) && !is.na(se) && se > 0) 2 * pnorm(-abs(est / se)) else p, df = df, stringsAsFactors = FALSE)

## the quantities for one set of gap coefficients with suffix sfx ("" for M0, "_<role>" for M1)
axis_rows <- function(spec, model, b, V, sfx, lab, minus_sfx = NULL) {
  nm <- names(b); out <- list()
  t_ <- function(v) paste0(v, sfx); m_ <- function(v) if (is.null(minus_sfx)) NULL else paste0(v, minus_sfx)
  for (v in c("dR", "dR2", "dMR", "dM", "dM2")) {
    a <- setNames(rep(0, length(nm)), nm); a[t_(v)] <- 1; if (!is.null(minus_sfx)) a[m_(v)] <- -1
    l <- lin(b, V, a); out[[length(out) + 1]] <- row_out(spec, model, paste0(lab, "b_", v), l[["est"]], l[["se"]])
  }
  for (q in list(c("R", "dR", 1, "u_R_less"), c("R", "dR", -1, "u_R_more"), c("M", "dM", 1, "u_M_weaker"), c("M", "dM", -1, "u_M_stronger"))) {
    r <- as.numeric(q[3]); a <- setNames(rep(0, length(nm)), nm); a[t_(q[2])] <- r; a[t_(paste0(q[2], "2"))] <- r^2
    if (!is.null(minus_sfx)) { a[m_(q[2])] <- -r; a[m_(paste0(q[2], "2"))] <- -r^2 }
    l <- lin(b, V, a); out[[length(out) + 1]] <- row_out(spec, model, paste0(lab, q[4]), l[["est"]], l[["se"]])
  }
  for (ax in list(c("axis_R", "dR", "dR2", "dMR"), c("axis_M", "dM", "dM2", "dMR"))) {
    w <- wald(b, V, Rsel(nm, t_(ax[2:4]), m_(ax[2:4]))); out[[length(out) + 1]] <- row_out(spec, model, paste0(lab, ax[1]), w[["chi2"]], NA, w[["p"]], 3L)
  }
  if (!is.null(minus_sfx)) { w <- wald(b, V, Rsel(nm, t_(GAP), m_(GAP))); out[[length(out) + 1]] <- row_out(spec, model, paste0(lab, "surface"), w[["chi2"]], NA, w[["p"]], 5L) }
  do.call(rbind, out)
}

## ---- step B (2026-11-19): conditional on a recorded MID onset, is the dispute fatal? One row per event (the chosen target). Weighted logit
## (quasibinomial, event weights w) of ev_fatal on type intercepts, type-specific gap terms of the CHOSEN target and its Z; two-way CR0
## (initiator x dispute, sandwich multiway). B0 = pooled gap terms, B1 = type-specific. Quantities as axis_rows (log-odds scale).
run_stepB <- function(spec, d, zc) {
  dd <- d[d$chosen == 1L, ]; stopifnot(!anyDuplicated(dd$ev_id))
  for (g in 1:(G - 1L)) dd[[paste0("ty_", ROLES[g + 1])]] <- as.numeric(dd$type == g)
  TY <- paste0("ty_", ROLES[-1])
  fitB <- function(terms) { f <- as.formula(paste("ev_fatal ~", paste(c(TY, terms, zc), collapse = " + ")))
    w_ <- character(0); m <- withCallingHandlers(glm(f, data = dd, weights = w, family = quasibinomial()), warning = function(e) { w_ <<- c(w_, conditionMessage(e)); invokeRestart("muffleWarning") })
    bad <- !m$converged || any(!is.finite(coef(m))) || any(grepl("fitted probabilities numerically 0 or 1", w_))
    V <- if (bad) NULL else { v <- sandwich::vcovCL(m, cluster = ~ init + disp, type = "HC0", cadjust = FALSE, multi0 = FALSE, fix = TRUE); v }
    list(m = m, V = V, bad = bad, wrn = w_) }
  out <- list()
  for (mm in list(list("B0", GAP), list("B1", GAP_G))) {
    f <- fitB(mm[[2]])
    if (f$bad) { out[[length(out) + 1]] <- row_out(spec, "INVALID", paste0(mm[[1]], " fit invalid: ", paste(unique(c(f$wrn, if (!f$m$converged) "not converged")), collapse = " | "))); next }
    b <- coef(f$m); b <- b[names(b) != "(Intercept)"]; V <- f$V[names(b), names(b)]
    if (mm[[1]] == "B0") out[[length(out) + 1]] <- axis_rows(spec, "B0", b, V, "", "") else {
      for (k in ROLES) out[[length(out) + 1]] <- axis_rows(spec, "B1", b, V, paste0("_", k), paste0(k, ":"))
      for (pr in lapply(combn(ROLES, 2, simplify = FALSE), rev))
        out[[length(out) + 1]] <- axis_rows(spec, "B1", b, V, paste0("_", pr[1]), paste0(pr[1], "-", pr[2], ":"), minus_sfx = paste0("_", pr[2]))
      w <- wald(b, V, Rsel(names(b), unlist(lapply(ROLES[-1], function(o) paste0(GAP, "_", o))), rep(paste0(GAP, "_", ROLES[1]), G - 1L)))
      out[[length(out) + 1]] <- row_out(spec, "B1", "surfaces_differ", w[["chi2"]], NA, w[["p"]], as.integer(w[["df"]]))
    }
    out[[length(out) + 1]] <- row_out(spec, mm[[1]], "n_events", nrow(dd)); out[[length(out) + 1]] <- row_out(spec, mm[[1]], "n_fatal", sum(dd$ev_fatal))
  }
  out[[length(out) + 1]] <- row_out(spec, "B1", "max_abs_dR", max(abs(dd$dR)))
  do.call(rbind, out)
}

run_spec <- function(spec, dR) {
  dRc <- NULL; if (is.list(dR)) { dRc <- dR$dRc; dR <- dR$dR }
  d <- D; d$dR <- dR
  zc <- Z_cols
  extra <- character(0)
  if (!is.null(dRc)) {
    for (g in 0:(G - 1L)) d[[paste0("dRc_", ROLES[g + 1])]] <- dRc * (d$type == g); zc <- c(Z_cols, paste0("dRc_", ROLES))
    if (DECOMP_TERMS == "full") {   # audit P2-8: the terms the quadratic utility implies once dR = baseline + change; tested jointly (decomp_extra)
      for (g in 0:(G - 1L)) { k <- ROLES[g + 1]; tg <- d$type == g
        d[[paste0("dRc2_", k)]] <- dRc^2 * tg; d[[paste0("dMdRc_", k)]] <- d$dM * dRc * tg; d[[paste0("dRdRc_", k)]] <- dR * dRc * tg
        extra <- c(extra, paste0(c("dRc2_", "dMdRc_", "dRdRc_"), k)) }
      zc <- c(zc, extra)
    }
  }
  if (!all(is.finite(d$dR))) stop("non-finite dR")
  d$dM2 <- d$dM^2; d$dR2 <- d$dR^2; d$dMR <- d$dM * d$dR
  for (g in 0:(G - 1L)) for (v in GAP) d[[paste0(v, "_", ROLES[g + 1])]] <- d[[v]] * (d$type == g)
  if (STEPB) return(run_stepB(spec, d, zc))
  ## GPT full audit 2026-11-19 P2-6: M0 is the FULLY pooled benchmark -- type-specific controls (x_<role> = x * 1[type == role], from
  ## 24_rel_variants.R or the decomposition) are summed back into one pooled column, so M0 is the same model under every type scheme.
  ## When such controls exist, M0t keeps them type-specific with pooled gap terms (the model nested in M1, for the M0t-vs-M1 comparison).
  zc0 <- zc; ts_groups <- list()
  for (z in zc) for (k in ROLES) if (endsWith(z, paste0("_", k))) { pre <- substr(z, 1, nchar(z) - nchar(k) - 1L); ts_groups[[pre]] <- c(ts_groups[[pre]], z) }
  ts_groups <- ts_groups[vapply(ts_groups, length, 0L) == G]
  for (pre in names(ts_groups)) { pc <- paste0(pre, "_pooled"); d[[pc]] <- Reduce(`+`, lapply(ts_groups[[pre]], function(z) d[[z]])); zc0 <- c(setdiff(zc0, ts_groups[[pre]]), pc) }
  ## GPT review P1-3 + full audit P2-11: each model has its own validity; an invalid M1 no longer removes a usable M0 (and vice versa)
  fit_w <- function(gap, z) { w_ <- character(0); m <- withCallingHandlers(fit_cl(d, gap, z), warning = function(w) { w_ <<- c(w_, conditionMessage(w)); invokeRestart("muffleWarning") })
    bad <- any(grepl("converge|infinite|singular", w_, ignore.case = TRUE)) || !all(is.finite(coef(m))) || m$iter >= 20L
    list(m = m, wrn = w_, bad = bad) }
  fits <- list(M0 = fit_w(GAP, zc0), M1 = fit_w(GAP_G, zc))
  if (length(ts_groups)) fits$M0t <- fit_w(GAP, zc)
  out <- list()
  for (nm_ in names(fits)) {
    f <- fits[[nm_]]
    if (f$bad) { out[[length(out) + 1]] <- row_out(spec, "INVALID", paste0(nm_, " fit invalid: ", if (length(f$wrn)) paste(unique(f$wrn), collapse = " | ") else "non-finite coefficients or iteration limit"))
                 out[[length(out) + 1]] <- row_out(spec, nm_, "iter", f$m$iter); next }
    b <- coef(f$m); V <- vc_2w(f$m, d); zz <- if (nm_ == "M0") zc0 else zc
    if (nm_ == "M1") {
      for (k in ROLES) out[[length(out) + 1]] <- axis_rows(spec, "M1", b, V, paste0("_", k), paste0(k, ":"))
      if (G >= 2L) {
        for (pr in lapply(combn(ROLES, 2, simplify = FALSE), rev))
          out[[length(out) + 1]] <- axis_rows(spec, "M1", b, V, paste0("_", pr[1]), paste0(pr[1], "-", pr[2], ":"), minus_sfx = paste0("_", pr[2]))
        w <- wald(b, V, Rsel(names(b), unlist(lapply(ROLES[-1], function(o) paste0(GAP, "_", o))), rep(paste0(GAP, "_", ROLES[1]), G - 1L)))
        out[[length(out) + 1]] <- row_out(spec, "M1", "surfaces_differ", w[["chi2"]], NA, w[["p"]], as.integer(w[["df"]]))
      }
      if (length(extra)) { w <- wald(b, V, Rsel(names(b), extra)); out[[length(out) + 1]] <- row_out(spec, "M1", "decomp_extra", w[["chi2"]], NA, w[["p"]], as.integer(w[["df"]])) }
    } else out[[length(out) + 1]] <- axis_rows(spec, nm_, b, V, "", "")
    for (z in zz) { a <- setNames(rep(0, length(b)), names(b)); a[z] <- 1; l <- lin(b, V, a)   # control coefficients (variants: aggr_*, def_*, dRc_*)
      out[[length(out) + 1]] <- row_out(spec, nm_, paste0("z:", z), l[["est"]], l[["se"]]) }
    out[[length(out) + 1]] <- row_out(spec, nm_, "loglik", f$m$loglik[2]); out[[length(out) + 1]] <- row_out(spec, nm_, "cgm_neg_eigen", attr(V, "n_neg"))
    out[[length(out) + 1]] <- row_out(spec, nm_, "iter", f$m$iter)
    out[[length(out) + 1]] <- row_out(spec, nm_, "n_warnings", length(f$wrn))   # only benign warnings can remain here
    if (length(f$wrn)) out[[length(out)]]$quantity <- paste0("n_warnings: ", paste(unique(f$wrn), collapse = " | "))
  }
  out[[length(out) + 1]] <- row_out(spec, "M1", "max_abs_dR", max(abs(d$dR)))
  out <- do.call(rbind, out)
  n_na <- sum(out$df > 1 & is.na(out$p))                          # GPT review P1-3: Wald tests that could not be computed (singular)
  rbind(out, row_out(spec, "ALL", "n_wald_na", n_na))
}

## ---- 2. validation against the dumped gap (baseline log IGO PageRank from the data build) --------------------------
posy <- readRDS(file.path(dump$data_dir, "sy_positions.rds"))
stopifnot(all(c("ccode", "year", "R_pr") %in% names(posy)))
Rb <- posy$R_pr[match(paste(key$ccode, key$year), paste(posy$ccode, posy$year))]
## the build's positions cover the state-years it uses; rows absent there are never looked up (asserted below)
need <- unique(c(i_row, j_row))
BASE_OK <- !anyNA(Rb[need])
if (!BASE_OK && POS_LAG == 1L) stop("build positions missing for state-years the fit uses")
if (!BASE_OK) logf("pos_lag %d: build positions (sy_positions) missing for %d used state-years; BASELINE row skipped", POS_LAG, sum(is.na(Rb[need])))
## gap_scale for the year convention uses ALL states of the year in sy_positions (as 10_clogit_roles.R); recompute on posy directly
if (dump$rel_scale == "year") {
  ysb <- vapply(split(posy$R_pr, posy$year), function(v) gap_scale(v, rep(1L, length(v))), 0)
  dRb <- (Rb[i_row] - Rb[j_row]) / ysb[yr_ev]
} else dRb <- (Rb[i_row] - Rb[j_row]) / gap_scale(posy$R_pr, posy$year)
mx <- if (BASE_OK) max(abs(dRb - D$dR)) else NA_real_; logf("validation: rebuilt baseline dR vs dumped dR, max abs diff %.2e", mx)
if (isTRUE(mx > 1e-8) && POS_LAG == 1L && DECOMP == 0L) stop("validation failed: the sweep does not reproduce the dumped relational gap")
if (POS_LAG != 1L || DECOMP > 0L) logf("validation skipped (pos_lag %d, decomp %d): the BASELINE row uses the build positions at %s without decomposition", POS_LAG, DECOMP, tlab(POS_LAG))
## the key table (COW state-years) and sy_positions should hold the same state-years; differences change the year scale and are reported
n_key_only <- sum(is.na(Rb)); n_pos_only <- sum(!paste(posy$ccode, posy$year) %in% paste(key$ccode, key$year))
logf("state-years in the measure key but not in sy_positions %d; in sy_positions but not in the key %d (year scales use each measure's own key)", n_key_only, n_pos_only)
if (YSCALE == "tm1" || POS_LAG != 1L) {   # validation above is for the dumped (year t) convention; the BASELINE row uses this run's convention (scale year = positions' year)
  ysb1 <- vapply(split(posy$R_pr, posy$year), function(v) gap_scale(v, rep(1L, length(v))), 0)
  dRb <- (Rb[i_row] - Rb[j_row]) / ysb1[as.character(D$yr - POS_LAG)]
  logf("BASELINE rebuilt with the positions-year scale (correlation with the dumped gap %.4f)", cor(dRb, D$dR))
}
## audit P2-8: the BASELINE row is not transformed like the candidates under --decomp / --merror / --scale=pooled, so it is skipped there
BASE_RUN <- BASE_OK && DECOMP == 0L && MERROR == "none" && SCALE == "year"
if (BASE_OK && !BASE_RUN) logf("BASELINE row skipped (decomp / merror / pooled scale apply to candidate measures only)")
base_rows <- if (BASE_RUN) run_spec("BASELINE|igo_pagerank_build|pr", dRb) else NULL

## ---- 3. sweep -------------------------------------------------------------------------------------------------------
todo <- seq_len(ncol(Rm)); if (MAX_SPECS > 0) todo <- todo[seq_len(min(MAX_SPECS, length(todo)))]
## GPT review P2-5: eligibility is decided per fit by make_dR() on the years this fit needs; the all-years flag from 21 is metadata only
ok_id <- rep(TRUE, length(todo))
logf("positions %s: %d specs (%d constant in >= 5%% of all years, metadata only; eligibility checked on the years each fit uses)", basename(POSF), length(todo), sum(!specs$identified[todo]))
t0 <- Sys.time()
one <- function(k) tryCatch(run_spec(colnames(Rm)[k], make_dR(Rm[, k])),
                            error = function(e) row_out(colnames(Rm)[k], if (grepl("^position constant", conditionMessage(e))) "SKIPPED" else "ERROR", conditionMessage(e)))
res <- if (CORES > 1) mclapply(todo[ok_id], one, mc.cores = CORES, mc.preschedule = FALSE) else lapply(todo[ok_id], one)
res <- c(if (!is.null(base_rows)) list(base_rows), res)
bad_res <- vapply(res, function(z) !is.data.frame(z), TRUE)
if (any(bad_res)) stop(sprintf("%d worker results are not data frames (killed fork?)", sum(bad_res)))
out <- do.call(rbind, res)
stopifnot(setequal(unique(out$spec_id), c(if (BASE_RUN) "BASELINE|igo_pagerank_build|pr", colnames(Rm)[todo])))   # review P2-8: no spec may vanish
n_err <- sum(out$model == "ERROR")
f_out <- file.path(OUT_DIR, sprintf("sweep_%s_%s.csv", SCHEME, sub("\\.rds$", "", basename(POSF))))
write.csv(out, f_out, row.names = FALSE)
## census per model (full audit P2-11): valid = the model's estimate rows exist; candidate specs only (BASELINE excluded)
cand <- out[out$spec_id != "BASELINE|igo_pagerank_build|pr", ]
nv <- function(m) { m2 <- if (STEPB) c(M0 = "B0", M1 = "B1")[[m]] else m; length(unique(cand$spec_id[cand$model == m2 & cand$quantity %in% c("u_R_less", paste0(ROLES[1], ":u_R_less"))])) }
inv <- function(m) { m2 <- if (STEPB) c(M0 = "B0", M1 = "B1")[[m]] else m; length(unique(cand$spec_id[cand$model == "INVALID" & startsWith(cand$quantity, paste0(m2, " fit invalid"))])) }
logf("CENSUS %s: specs %d | valid M0 %d, M1 %d | skipped (constant in a used year) %d | invalid fits M0 %d, M1 %d | errors %d | valid fits with a Wald test not computable %d | in %.1f min",
     f_out, length(todo), nv("M0"), nv("M1"),
     length(unique(cand$spec_id[cand$model == "SKIPPED"])), inv("M0"), inv("M1"), n_err,
     length(unique(out$spec_id[out$quantity == "n_wald_na" & out$est > 0])), as.numeric(difftime(Sys.time(), t0, units = "mins")))
if (n_err > 0) logf("ERRORS: %s", paste(unique(out$quantity[out$model == "ERROR"]), collapse = " | "))
