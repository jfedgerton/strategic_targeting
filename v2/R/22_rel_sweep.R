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
## A spec is refused (row model = SKIPPED, with the reason) when its position is constant in any year an event uses (t for the year
## scale, t-1 for positions): the year scale would be ~0 and the gap ~1e13 (independent review 2026-11-18, P1-1).
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
YSCALE <- opt("yscale", "t"); stopifnot(YSCALE %in% c("t", "tm1"))   # PI 2026-11-18: both conventions (scale from year t, as 10_clogit_roles.R, or t-1, the positions' year)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
logf <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
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
SCHEME <- paste0(dump$types, if (dump$rel_scale == "year") paste0("_relyear", if (YSCALE == "tm1") "tm1" else "") else "", if (isTRUE(dump$unit_weights)) "_uw" else "", "_", dump$sample)
if (YSCALE == "tm1" && dump$rel_scale != "year") stop("--yscale=tm1 needs a dump with rel_scale = year")
logf("dump %s: scheme %s, %d rows, %d strata, types %s, rel_scale %s, commit %s", basename(DUMP), SCHEME, nrow(D), length(unique(D$ev_id)),
     paste(ROLES, collapse = "/"), dump$rel_scale, dump$git_commit)
P <- readRDS(POSF); key <- P$key; Rm <- P$R; specs <- P$specs
stopifnot(identical(colnames(Rm), specs$spec_id), nrow(key) == nrow(Rm), !anyDuplicated(key[, c("ccode", "year")]))
kidx <- setNames(seq_len(nrow(key)), paste(key$ccode, key$year))
i_row <- kidx[paste(D$c1, D$yr - 1L)]; j_row <- kidx[paste(D$c2, D$yr - 1L)]
if (anyNA(i_row) || anyNA(j_row)) stop(sprintf("positions missing at t-1 for %d initiator / %d target rows", sum(is.na(i_row)), sum(is.na(j_row))))
yr_ev <- as.character(D$yr)
key_years <- split(seq_len(nrow(key)), key$year)
stopifnot(all(unique(c(yr_ev, as.character(D$yr - 1L))) %in% names(key_years)))

used_years <- sort(unique(c(D$yr, D$yr - 1L)))
make_dR <- function(R) {                                          # R: position per key row
  sds <- vapply(key_years[as.character(used_years)], function(ix) sd(R[ix]), 0)
  if (any(sds < 1e-8)) stop(sprintf("position constant in event-used year(s) %s", paste(used_years[sds < 1e-8], collapse = " ")))
  num <- R[i_row] - R[j_row]
  if (dump$rel_scale == "year") {
    sc_yr <- if (YSCALE == "t") yr_ev else as.character(D$yr - 1L)
    ev_years <- unique(sc_yr)
    ys <- vapply(key_years[ev_years], function(ix) gap_scale(R[ix], rep(1L, length(ix))), 0)
    stopifnot(all(ys > 1e-8))
    num / ys[sc_yr]
  } else num / gap_scale(R, key$year)
}

## ---- 1. model machinery (as 10_clogit_roles.R section 4) ----------------------------------------------------------
GAP <- c("dM", "dM2", "dR", "dR2", "dMR")
GAP_G <- as.vector(outer(GAP, ROLES, paste, sep = "_"))
fit_cl <- function(d, terms) coxph(as.formula(paste("Surv(rep(1, nrow(d)), chosen) ~", paste(c(terms, Z_cols), collapse = " + "), "+ strata(ev_id)")),
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

run_spec <- function(spec, dR) {
  d <- D; d$dR <- dR
  if (!all(is.finite(d$dR))) stop("non-finite dR")
  d$dM2 <- d$dM^2; d$dR2 <- d$dR^2; d$dMR <- d$dM * d$dR
  for (g in 0:(G - 1L)) for (v in GAP) d[[paste0(v, "_", ROLES[g + 1])]] <- d[[v]] * (d$type == g)
  wrn <- character(0)
  M0 <- withCallingHandlers(fit_cl(d, GAP), warning = function(w) { wrn <<- c(wrn, conditionMessage(w)); invokeRestart("muffleWarning") })
  M1 <- withCallingHandlers(fit_cl(d, GAP_G), warning = function(w) { wrn <<- c(wrn, conditionMessage(w)); invokeRestart("muffleWarning") })
  for (m in list(M0, M1)) if (!all(is.finite(coef(m)))) stop("non-finite coefficients")
  V0 <- vc_2w(M0, d); V1 <- vc_2w(M1, d); b0 <- coef(M0); b1 <- coef(M1)
  out <- list(axis_rows(spec, "M0", b0, V0, "", ""))
  for (k in ROLES) out[[length(out) + 1]] <- axis_rows(spec, "M1", b1, V1, paste0("_", k), paste0(k, ":"))
  for (pr in lapply(combn(ROLES, 2, simplify = FALSE), rev))
    out[[length(out) + 1]] <- axis_rows(spec, "M1", b1, V1, paste0("_", pr[1]), paste0(pr[1], "-", pr[2], ":"), minus_sfx = paste0("_", pr[2]))
  w <- wald(b1, V1, Rsel(names(b1), unlist(lapply(ROLES[-1], function(o) paste0(GAP, "_", o))), rep(paste0(GAP, "_", ROLES[1]), G - 1L)))
  out[[length(out) + 1]] <- row_out(spec, "M1", "surfaces_differ", w[["chi2"]], NA, w[["p"]], as.integer(w[["df"]]))
  out[[length(out) + 1]] <- row_out(spec, "M0", "loglik", M0$loglik[2]); out[[length(out) + 1]] <- row_out(spec, "M1", "loglik", M1$loglik[2])
  out[[length(out) + 1]] <- row_out(spec, "M0", "cgm_neg_eigen", attr(V0, "n_neg")); out[[length(out) + 1]] <- row_out(spec, "M1", "cgm_neg_eigen", attr(V1, "n_neg"))
  out[[length(out) + 1]] <- row_out(spec, "M0", "iter", M0$iter); out[[length(out) + 1]] <- row_out(spec, "M1", "iter", M1$iter)
  out[[length(out) + 1]] <- row_out(spec, "M1", "max_abs_dR", max(abs(d$dR)))
  out[[length(out) + 1]] <- row_out(spec, "M1", "n_warnings", length(wrn))
  if (length(wrn)) out[[length(out)]]$quantity <- paste0("n_warnings: ", paste(unique(wrn), collapse = " | "))
  do.call(rbind, out)
}

## ---- 2. validation against the dumped gap (baseline log IGO PageRank from the data build) --------------------------
posy <- readRDS(file.path(dump$data_dir, "sy_positions.rds"))
stopifnot(all(c("ccode", "year", "R_pr") %in% names(posy)))
Rb <- posy$R_pr[match(paste(key$ccode, key$year), paste(posy$ccode, posy$year))]
## the build's positions cover the state-years it uses; rows absent there are never looked up (asserted below)
need <- unique(c(i_row, j_row)); stopifnot(!anyNA(Rb[need]))
## gap_scale for the year convention uses ALL states of the year in sy_positions (as 10_clogit_roles.R); recompute on posy directly
if (dump$rel_scale == "year") {
  ysb <- vapply(split(posy$R_pr, posy$year), function(v) gap_scale(v, rep(1L, length(v))), 0)
  dRb <- (Rb[i_row] - Rb[j_row]) / ysb[yr_ev]
} else dRb <- (Rb[i_row] - Rb[j_row]) / gap_scale(posy$R_pr, posy$year)
mx <- max(abs(dRb - D$dR)); logf("validation: rebuilt baseline dR vs dumped dR, max abs diff %.2e", mx)
if (mx > 1e-8) stop("validation failed: the sweep does not reproduce the dumped relational gap")
## the key table (COW state-years) and sy_positions should hold the same state-years; differences change the year scale and are reported
n_key_only <- sum(is.na(Rb)); n_pos_only <- sum(!paste(posy$ccode, posy$year) %in% paste(key$ccode, key$year))
logf("state-years in the measure key but not in sy_positions %d; in sy_positions but not in the key %d (year scales use each measure's own key)", n_key_only, n_pos_only)
if (YSCALE == "tm1") {   # validation above is for the dumped (year t) convention; the BASELINE row uses this run's convention
  ysb1 <- vapply(split(posy$R_pr, posy$year), function(v) gap_scale(v, rep(1L, length(v))), 0)
  dRb <- (Rb[i_row] - Rb[j_row]) / ysb1[as.character(D$yr - 1L)]
  logf("yscale = tm1: BASELINE rebuilt with the year t-1 scale (correlation with the dumped gap %.4f)", cor(dRb, D$dR))
}
base_rows <- run_spec("BASELINE|igo_pagerank_build|pr", dRb)

## ---- 3. sweep -------------------------------------------------------------------------------------------------------
todo <- seq_len(ncol(Rm)); if (MAX_SPECS > 0) todo <- todo[seq_len(min(MAX_SPECS, length(todo)))]
ok_id <- specs$identified[todo]
logf("positions %s: %d specs (%d marked not identified, skipped)", basename(POSF), length(todo), sum(!ok_id))
t0 <- Sys.time()
one <- function(k) tryCatch(run_spec(colnames(Rm)[k], make_dR(Rm[, k])),
                            error = function(e) row_out(colnames(Rm)[k], if (grepl("^position constant", conditionMessage(e))) "SKIPPED" else "ERROR", conditionMessage(e)))
res <- if (CORES > 1) mclapply(todo[ok_id], one, mc.cores = CORES, mc.preschedule = FALSE) else lapply(todo[ok_id], one)
res <- c(list(base_rows), res, lapply(todo[!ok_id], function(k) row_out(colnames(Rm)[k], "SKIPPED", "not identified (constant position within year)")))
bad_res <- vapply(res, function(z) !is.data.frame(z), TRUE)
if (any(bad_res)) stop(sprintf("%d worker results are not data frames (killed fork?)", sum(bad_res)))
out <- do.call(rbind, res)
stopifnot(setequal(unique(out$spec_id), c("BASELINE|igo_pagerank_build|pr", colnames(Rm)[todo])))   # review P2-8: no spec may vanish
n_err <- sum(out$model == "ERROR")
f_out <- file.path(OUT_DIR, sprintf("sweep_%s_%s.csv", SCHEME, sub("\\.rds$", "", basename(POSF))))
write.csv(out, f_out, row.names = FALSE)
logf("wrote %s: %d specs in %.1f min; errors %d; skipped %d (constant in a used year %d, flagged in 21 %d); specs with coxph warnings %d",
     f_out, length(todo), as.numeric(difftime(Sys.time(), t0, units = "mins")), n_err, length(unique(out$spec_id[out$model == "SKIPPED"])),
     sum(out$model == "SKIPPED" & grepl("^position constant", out$quantity)), sum(!ok_id), sum(grepl("^n_warnings: ", out$quantity)))
if (n_err > 0) logf("ERRORS: %s", paste(unique(out$quantity[out$model == "ERROR"]), collapse = " | "))
