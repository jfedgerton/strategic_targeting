## 10_clogit_roles.R (2026-11-16, PI request): the simplest test of whether target selection differs by initiation type.
## Stage 2 only, as a weighted conditional logit (survival::coxph, Breslow = exact conditional logit with one chosen row per stratum; the
## event weight is constant within the stratum, so it weights the event's log-likelihood), on EXACTLY the events, stage-2 choice sets, gaps
## (dM, dR), directed Z and D8 event weights that M1c sees: the cell is assembled with assemble_cell_spec() and build_stan_data() under the
## grid flags (directed, roles bcj, joiner_set all, full choice set) and the three bcj groups are then recoded here (switch --types).
## Omitted relative to M1c: stage 1 (initiation, tau), the initiator ideal-point intercept a_i and the target effects b_j. The magnitudes
## (mu, sigma, rho) are therefore NOT directly comparable with M1c's (audit A3 P1-4: omitted a_i inflates sigma and pulls rho toward 0).
## coxph's log-likelihood differs from Stan's weighted stage-2 log-likelihood by the constant -sum(w log w); coefficients and tests are unaffected.
##
## Types (--types, PI decision 2026-11-16: bilateral vs coalition is the main coding; the 14-initiator coalition-originator group is too thin):
##   bc          (main)       bil = bcj 0 (alone originator); coal = bcj 1 (coalition originator) + bcj 2 (joiner)
##   bc_nocoorig (check)      bil = bcj 0; coal = bcj 2 only; the bcj-1 events are DROPPED from stage 2 (does merging them matter?)
##   bcj         (appendix)   bil, coal (= coalition originator), join
##
## The spatial utility -1/2 (x - mu)' P (x - mu) = b'x + 1/2 x'Hx + const with H = -P and b = P mu, so for each type g
##   u = b1 dM + b2 dM^2 + b3 dR + b4 dR^2 + b5 dM*dR    ->    H = [[2 b2, b5], [b5, 2 b4]],  Sigma = -H^{-1},  mu = Sigma b.
## An ideal point exists only where H is negative definite (a peak). Reporting rule (audit F1): mu, sigma, rho and their intervals are printed
## for a type only when (i) the point estimate has a peak, (ii) at least 95% of the bootstrap draws have a peak, and (iii) the point mu lies
## inside the range of that type's observed alternatives on both axes; otherwise the type is reported as "no peak" with its bootstrap share.
## A contrast between two types is printed only when both pass. Percentile intervals use the peaked draws only, which is why rule (ii) exists.
##
## Models (all with the directed Z and strata = event):
##   M0  pooled: the five gap terms shared by all types
##   M1  type-specific: the five gap terms separately by type                   <- the main model
##   M1s type-specific linear and squared terms, ONE shared cross coefficient (a common utility interaction P12; NOT a shared rho)
## Tests (Wald, CR0 cluster-robust). PRIMARY clustering (PI decision 2026-11-16, joint audit M1): two-way initiator x dispute (Cameron, Gelbach
##   and Miller 2011, negative eigenvalues set to zero); initiator-only and target-year are reported as sensitivities. Tests:
##   M1 vs M0 (5 (G - 1) df) = "the utility surfaces differ by type"; equal cross coefficients (G - 1 df) = "the utility interaction differs by
##   type" (NOT a test of equal rho, audit A3 P1-5); per pair of types, a 3-df test per axis restricting linear, squared AND cross terms
##   (GPT review P1: the axis effect is b_k + 2 b_kk x_k + b_MR x_other) and the 5-df surface test. Term-by-term contrasts use the two-way SE.
## Intervals: cluster bootstrap, seed 123, percentile, resampling DISPUTES by default (--boot_unit=dispute; PI decision M2: initiators and disputes
##   form one giant connected component, so no unit nests both, and resampling initiator x dispute cells would ignore both dependencies; the
##   dispute is the dependence that changed step B). --boot_unit=initiator resamples initiating states. Intervals are reported only when B >= 200
##   and >= 90% of the draws succeed (joint audit B3); a draw whose fit errors or warns is dropped for that model and logged with its reason.
## --loo=TRUE adds leave-one-dispute-out and leave-one-initiator-out refits of M1 (joint audit M3), written to <cell>_loo.csv.
##
## Usage (same flags as the grid cell; --roles=bcj is required because the types are recoded from it):
##   Rscript v2/R/10_clogit_roles.R --sample=fmid --lag=3 --model=M1c --comparison=directed --entry=later --variant=all --roles=bcj --joiner_set=all --types=bc
##   --events=dymid (instead of --entry=later) for jointly fatal, --reach=capdist for the reach terms, --boot=1000, --boot_unit=dispute|initiator,
##   --loo=TRUE. Unknown, retired or duplicated options are refused (joint audit B2).
suppressPackageStartupMessages({ library(dplyr); library(survival) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
R_DIR <- dirname(normalizePath(this_file))
source(file.path(R_DIR, "spec.R")); source(file.path(R_DIR, "helpers_data.R")); source(file.path(R_DIR, "helpers_model.R"))
source(file.path(R_DIR, "assemble_cell.R"))
Sys.setenv(V2_ROOT = normalizePath(file.path(R_DIR, "..")))   # audit F5 (A1 P2-3): read and write this checkout only, as 02_fit_model.R does

args <- commandArgs(trailingOnly = TRUE)
## joint audit B2: refuse retired options, unknown flags, bare flags and duplicated flags (get_arg would silently take the first)
refuse_archived_options(args, "10_clogit_roles.R")
SCRIPT_FLAGS <- c("sample", "lag", "model", "boot", "types", "boot_unit", "loo")
KNOWN_FLAGS <- c(names(SPEC_CHOICES), "ame_K", "tag", SCRIPT_FLAGS)
bad <- args[!grepl("^--[A-Za-z_0-9]+=.+$", args)]
if (length(bad)) stop(sprintf("arguments must be --name=value; refused: %s", paste(bad, collapse = " ")))
flag_names <- sub("^--([A-Za-z_0-9]+)=.*$", "\\1", args)
if (any(!flag_names %in% KNOWN_FLAGS)) stop(sprintf("unknown option(s) for 10_clogit_roles.R: %s", paste(unique(flag_names[!flag_names %in% KNOWN_FLAGS]), collapse = ", ")))
if (anyDuplicated(flag_names)) stop(sprintf("option(s) given more than once: %s", paste(unique(flag_names[duplicated(flag_names)]), collapse = ", ")))
SAMPLE <- get_arg(args, "sample", "fmid"); LAG <- as.integer(get_arg(args, "lag", 3)); MODEL <- get_arg(args, "model", "M1c")
B <- as.integer(get_arg(args, "boot", 1000))
TYPES <- get_arg(args, "types", "bc"); BOOT_UNIT <- get_arg(args, "boot_unit", "dispute"); LOO <- get_arg(args, "loo", FALSE)
stopifnot(SAMPLE %in% ALL_SAMPLES, LAG == 3L, MODEL == "M1c", B >= 0, TYPES %in% c("bc", "bc_nocoorig", "bcj"), BOOT_UNIT %in% c("dispute", "initiator"), is.logical(LOO))
spec <- parse_spec(args)
## audit F5 (A1 P2-1): the cells M1c is defined on, and nothing else (02_fit_model.R refuses the same)
stopifnot(spec_is_directed(spec), spec$roles == "bcj", spec_variant(spec) == "all", (if (is.null(spec$joiner_set)) "all" else spec$joiner_set) == "all",
          spec$controls == "minimal", spec$context == "all", (if (is.null(spec$mu_shift)) "none" else spec$mu_shift) == "none")
cell <- paste0("clogit_", cell_name(SAMPLE, LAG, MODEL, spec), "_", TYPES)
log_open(paste0("10_", cell))
logf("cell %s; spec %s; types %s; bootstrap B = %d over %ss, seed 123; leave-one-out %s", cell, spec_string(spec), TYPES, B, BOOT_UNIT, LOO)

## ---- 1. assemble exactly as 02_fit_model.R does for M1c ----------------------------------------------------
data_dir <- spec_data_dir(spec, v2_root())
manifest <- read.csv(file.path(data_dir, "manifest.csv"), stringsAsFactors = FALSE)
A <- assemble_cell_spec(SAMPLE, LAG, spec, data_dir = data_dir, manifest = manifest, model = MODEL)
sy <- A$sy; alts <- A$alts; ev <- A$ev; Z_cols <- A$Z_cols
dat <- build_stan_data(sy, alts, ev, A$W_cols, A$X_cols, Z_cols, joint = 1L, sampled = A$sampled, n_full = A$n_full,
                       x_ref = A$x_ref, allow_orphan_initiations = TRUE)
gt <- A$group_table
## audit F5 (A1 P2-2): the bcj labels and order are asserted, as in 02_fit_model.R
stopifnot(A$sampled == 0L, identical(gt$grp, 0:2), identical(gt$group, stage2_group_labels("bcj")), identical(gt$group, c("bilateral", "coalition_originator", "joiner")),
          "grp" %in% names(ev), length(ev$grp) == dat$E, all(ev$grp %in% 0:2), identical(as.integer(dat$ev_sy), as.integer(ev$sy_row)))
logf("bcj events %d (bilateral %d, coalition originator %d, joiner %d); weight sums %.1f / %.1f / %.1f; Z = %s",
     dat$E, sum(ev$grp == 0), sum(ev$grp == 1), sum(ev$grp == 2),
     sum(dat$ev_w[ev$grp == 0]), sum(dat$ev_w[ev$grp == 1]), sum(dat$ev_w[ev$grp == 2]), paste(Z_cols, collapse = ","))

## ---- 2. types: recode the bcj groups (0-based type index t per event; keep = events in stage 2) ------------------
if (TYPES == "bcj") {
  ROLES <- c("bil", "coal", "join"); t_ev <- ev$grp; keep <- rep(TRUE, dat$E)
} else if (TYPES == "bc") {
  ROLES <- c("bil", "coal"); t_ev <- as.integer(ev$grp > 0L); keep <- rep(TRUE, dat$E)
} else {
  ROLES <- c("bil", "coal"); t_ev <- as.integer(ev$grp > 0L); keep <- ev$grp != 1L
}
G <- length(ROLES)
stopifnot(length(t_ev) == dat$E, all(t_ev %in% 0:(G - 1L)), all(tabulate(t_ev[keep] + 1L, nbins = G) > 0L))   # every type has events
logf("types %s: %s; events kept %d of %d (dropped %d)", TYPES,
     paste(sprintf("%s %d (weights %.1f)", ROLES, tabulate(t_ev[keep] + 1L, nbins = G), vapply(0:(G - 1L), function(g) sum(dat$ev_w[keep & t_ev == g]), 0)), collapse = "; "),
     sum(keep), dat$E, sum(!keep))

## ---- 3. long frame: one stratum per kept event = its state-year's stage-2 set (in_s2 = 1), the chosen row flagged ---------
init_cc <- A$state_codes[sy$state_idx]                  # initiating state's ccode per state-year (cluster id)
rows <- lapply(which(keep), function(e) {
  r <- dat$ev_start[e]:dat$ev_end[e]; r <- r[dat$in_s2[r] == 1L]
  data.frame(ev_id = e, row = r, chosen = as.integer(r == dat$ev_row[e]), type = t_ev[e], w = dat$ev_w[e], init = init_cc[dat$ev_sy[e]])
})
L <- bind_rows(rows)
stopifnot(all(tapply(L$chosen, L$ev_id, sum) == 1L), all(tapply(L$row, L$ev_id, length) >= 2L), length(unique(L$ev_id)) == sum(keep))
L$dM <- dat$x[L$row, 1]; L$dR <- dat$x[L$row, 2]
Zm <- dat$Z[L$row, , drop = FALSE]; colnames(Zm) <- Z_cols; L <- cbind(L, Zm)
stopifnot(!anyNA(L[, c("dM", "dR", Z_cols)]))
L$dM2 <- L$dM^2; L$dR2 <- L$dR^2; L$dMR <- L$dM * L$dR
GAP <- c("dM", "dM2", "dR", "dR2", "dMR")
for (g in 0:(G - 1L)) for (v in GAP) L[[paste0(v, "_", ROLES[g + 1])]] <- L[[v]] * (L$type == g)
GAP_G <- as.vector(outer(GAP, ROLES, paste, sep = "_"))           # dM_bil, dM2_bil, ..., dMR_<last type>
## data support per type (reporting rule iii): the range of dM and dR over the type's alternatives
SUPPORT <- lapply(0:(G - 1L), function(g) { s <- L$type == g; c(dM_lo = min(L$dM[s]), dM_hi = max(L$dM[s]), dR_lo = min(L$dR[s]), dR_hi = max(L$dR[s])) })
names(SUPPORT) <- ROLES
logf("long frame: %d rows, %d strata, median set size %.0f", nrow(L), length(unique(L$ev_id)), median(table(L$ev_id)))
## cluster ids for the dependence sensitivity (GPT review 2026-11-16 P2-5): the event's dispute (events.rds dispnum of its non-escalation row,
## matched on (initiator, chosen target, year)) and its chosen target-year. Two-way = initiator x dispute; target-year is a separate sensitivity.
ek <- data.frame(ev_id = which(keep), initiator = as.integer(sy$ccode[ev$sy_row[keep]]), target = as.integer(alts$ccode2[ev$alt_row[keep]]),
                 year = as.integer(sy$year[ev$sy_row[keep]]))
stopifnot(!anyNA(ek), all(alts$ccode1[ev$alt_row[keep]] == ek$initiator), all(alts$year[ev$alt_row[keep]] == ek$year), !anyDuplicated(ek[, c("initiator", "target", "year")]),
          all(ek$initiator == init_cc[dat$ev_sy[ek$ev_id]]))
ev_rds <- readRDS(file.path(data_dir, "events.rds"))
stopifnot(all(c("initiator", "target", "year", "escalation", "dispnum") %in% names(ev_rds)))
ne <- ev_rds[!ev_rds$escalation, c("initiator", "target", "year", "dispnum")]
stopifnot(!anyDuplicated(ne[, c("initiator", "target", "year")]))
ek <- merge(ek, ne, by = c("initiator", "target", "year"), all.x = TRUE, sort = FALSE)
if (anyNA(ek$dispnum)) stop(sprintf("%d events have no non-escalation row with a dispnum in events.rds", sum(is.na(ek$dispnum))))
ek <- ek[match(unique(L$ev_id), ek$ev_id), ]
stopifnot(nrow(ek) == length(unique(L$ev_id)))
ek$disp <- ek$dispnum; ek$tgtyr <- paste(ek$target, ek$year, sep = "_"); ek$init_disp <- paste(ek$initiator, ek$dispnum, sep = "_")
L$disp <- ek$disp[match(L$ev_id, ek$ev_id)]; L$tgtyr <- ek$tgtyr[match(L$ev_id, ek$ev_id)]; L$init_disp <- ek$init_disp[match(L$ev_id, ek$ev_id)]
stopifnot(!anyNA(L[, c("disp", "tgtyr", "init_disp")]))
logf("clusters: %d initiators, %d disputes, %d initiator-disputes, %d target-years", length(unique(ek$initiator)), length(unique(ek$disp)),
     length(unique(ek$init_disp)), length(unique(ek$tgtyr)))

## ---- 4. fits -----------------------------------------------------------------------------------------------
fit_cl <- function(d, terms, robust = TRUE) {
  f <- as.formula(paste("Surv(rep(1, nrow(d)), chosen) ~", paste(c(terms, Z_cols), collapse = " + "), "+ strata(ev_id)"))
  if (robust) coxph(f, data = d, weights = w, method = "breslow", cluster = init)
  else coxph(f, data = d, weights = w, method = "breslow")
}
TERMS_M1S <- c(setdiff(GAP_G, paste0("dMR_", ROLES)), "dMR")
M0 <- fit_cl(L, GAP)
M1 <- fit_cl(L, GAP_G)
M1s <- fit_cl(L, TERMS_M1S)
for (m in list(M0, M1, M1s)) stopifnot(all(is.finite(coef(m))), !anyNA(coef(m)))

## cluster-robust variances from the per-row dfbeta residuals summed within clusters (CR0; coxph's own robust variance is the initiator case)
## two-way (initiator x dispute) = V(initiator) + V(dispute) - V(initiator-dispute) (Cameron, Gelbach and Miller 2011)
vc_cl <- function(m, g) { D <- rowsum(residuals(m, type = "dfbeta", weighted = TRUE), g, reorder = FALSE); V <- crossprod(D); dimnames(V) <- list(names(coef(m)), names(coef(m))); V }
## Cameron, Gelbach and Miller (2011): the two-way variance need not be positive definite; negative eigenvalues are set to zero (as sandwich fix = TRUE)
vc_set <- function(m, d = L) {
  V_init <- vc_cl(m, d$init); V_disp <- vc_cl(m, d$disp); V_id <- vc_cl(m, d$init_disp)
  V_2w <- V_init + V_disp - V_id
  e <- eigen(V_2w, symmetric = TRUE); n_neg <- sum(e$values <= 0)
  if (n_neg > 0) V_2w <- e$vectors %*% diag(pmax(e$values, 0)) %*% t(e$vectors)
  dimnames(V_2w) <- dimnames(V_init)
  list(init = V_init, twoway = V_2w, tgtyr = vc_cl(m, d$tgtyr), twoway_psd = TRUE, twoway_n_neg_fixed = n_neg)
}
wald <- function(m, R, r = rep(0, nrow(R)), V = vcov(m)) {            # Wald test of R beta = r with variance V (default: coxph robust, initiator)
  b <- coef(m); d <- R %*% b - r
  s <- as.numeric(t(d) %*% solve(R %*% V %*% t(R)) %*% d); c(chi2 = s, df = nrow(R), p = pchisq(s, nrow(R), lower.tail = FALSE))
}
nm1 <- names(coef(M1))
Rm <- function(pairs) { R <- matrix(0, length(pairs), length(nm1), dimnames = list(NULL, nm1)); for (i in seq_along(pairs)) { R[i, pairs[[i]][1]] <- 1; R[i, pairs[[i]][2]] <- -1 }; R }
REF <- ROLES[1]; OTHERS <- ROLES[-1]
PAIRS_REF <- unlist(lapply(GAP, function(v) lapply(OTHERS, function(o) c(paste0(v, "_", o), paste0(v, "_", REF)))), recursive = FALSE)
stopifnot(length(PAIRS_REF) == 5L * (G - 1L))
VS <- vc_set(M1)
stopifnot(isTRUE(all.equal(unname(VS$init), unname(vcov(M1)), tolerance = 1e-6)))   # the dfbeta construction reproduces coxph's initiator-robust variance
logf("two-way (initiator x dispute) variance: %d non-positive eigenvalues set to zero (CGM fix)", VS$twoway_n_neg_fixed)
R_roles <- Rm(PAIRS_REF); R_cross <- Rm(lapply(OTHERS, function(o) c(paste0("dMR_", o), paste0("dMR_", REF))))
## PRIMARY = two-way initiator x dispute (chi2 and p); initiator-only and target-year p-values as sensitivities
wald3 <- function(R, m = M1, V3 = VS) { b <- wald(m, R, V = V3$twoway); a <- wald(m, R, V = V3$init); t <- wald(m, R, V = V3$tgtyr)
  c(chi2 = b[["chi2"]], df = b[["df"]], p = b[["p"]], p_initiator = a[["p"]], p_tgtyr = t[["p"]]) }
W_roles <- wald3(R_roles); W_cross <- wald3(R_cross)
logf("Wald (CR0; PRIMARY initiator x dispute; sensitivities initiator / target-year): utility surfaces differ by type (M1 vs M0) chi2 %.2f df %d p %.4g (%.4g / %.4g); cross coefficient differs by type chi2 %.2f df %d p %.4g (%.4g / %.4g) (NOT a test of equal rho)",
     W_roles[["chi2"]], W_roles[["df"]], W_roles[["p"]], W_roles[["p_initiator"]], W_roles[["p_tgtyr"]], W_cross[["chi2"]], W_cross[["df"]], W_cross[["p"]], W_cross[["p_initiator"]], W_cross[["p_tgtyr"]])
## pairwise contrasts, term by term: difference of the type coefficients, two-way SE (initiator x dispute), 95% interval, p; plus per-axis 3-df and
## whole-surface 5-df tests. GPT review 2026-11-16 P1: an axis test must include the cross term (the axis effect is b_k + 2 b_kk x_k + b_MR x_other);
## the per-axis test below restricts linear, squared AND cross, i.e. "the two types' utility does not differ along this axis anywhere".
V1 <- VS$twoway; b1 <- coef(M1)
PAIRS_R <- combn(ROLES, 2, simplify = FALSE); PAIRS_R <- lapply(PAIRS_R, rev)          # (later type, earlier type), e.g. coal - bil
pw <- NULL
for (pr in PAIRS_R) {
  for (v in GAP) {
    a <- paste0(v, "_", pr[1]); z <- paste0(v, "_", pr[2]); d <- b1[[a]] - b1[[z]]; se <- sqrt(V1[a, a] + V1[z, z] - 2 * V1[a, z])
    pw <- rbind(pw, data.frame(pair = paste(pr, collapse = "-"), term = v, diff = d, se = se, lo = d - 1.96 * se, hi = d + 1.96 * se, p = 2 * pnorm(-abs(d / se)), df = 1,
                               p_initiator = NA_real_, p_tgtyr = NA_real_))
  }
  for (blk in list(c("dM", "dM2", "dMR"), c("dR", "dR2", "dMR"), GAP)) {
    w2 <- wald3(Rm(lapply(blk, function(v) c(paste0(v, "_", pr[1]), paste0(v, "_", pr[2])))))
    pw <- rbind(pw, data.frame(pair = paste(pr, collapse = "-"), term = paste0("joint:", paste(blk, collapse = "+")), diff = NA, se = NA, lo = NA, hi = NA,
                               p = w2[["p"]], df = w2[["df"]], p_initiator = w2[["p_initiator"]], p_tgtyr = w2[["p_tgtyr"]]))
  }
}
pw$note <- ifelse(pw$term == "dMR", "utility interaction P12, not rho",
           ifelse(pw$term == "dM", "difference at dM = dR = 0 only (squared and cross terms can change its sign elsewhere)",
           ifelse(grepl("^joint:dM\\+dM2\\+dMR$|^joint:dR\\+dR2\\+dMR$", pw$term), "axis effect differs anywhere (linear, squared, cross); not mu or sigma", "")))

## ---- 5. derived quantities ------------------------------------------------------------------------------------
QN <- c("muM", "muR", "sdM", "sdR", "rho")
derive <- function(b, pooled = FALSE) {                   # b: named coefficient vector; per type mu, sigma, rho (NA if no peak) and peak flag
  out <- c()
  for (g in if (pooled) "pooled" else ROLES) {
    gb <- function(v) { nm <- if (pooled) v else paste0(v, "_", g); if (!pooled && !nm %in% names(b)) nm <- v; stopifnot(nm %in% names(b)); b[[nm]] }   # M1s: shared dMR
    H <- matrix(c(2 * gb("dM2"), gb("dMR"), gb("dMR"), 2 * gb("dR2")), 2)
    ev_H <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    if (all(ev_H < 0)) {
      S <- -solve(H); mu <- as.vector(S %*% c(gb("dM"), gb("dR")))
      out <- c(out, setNames(c(mu, sqrt(diag(S)), S[1, 2] / sqrt(S[1, 1] * S[2, 2]), 1), paste0(c(QN, "peak"), "_", g)))
    } else out <- c(out, setNames(c(rep(NA, 5), 0), paste0(c(QN, "peak"), "_", g)))
  }
  if (!pooled) for (q in QN) for (pr in PAIRS_R) out[paste0(q, "_", pr[1], "-", pr[2])] <- out[paste0(q, "_", pr[1])] - out[paste0(q, "_", pr[2])]
  out
}
pt <- list(M0 = derive(coef(M0), pooled = TRUE), M1 = derive(coef(M1)), M1s = derive(coef(M1s)))

## ---- 6. cluster bootstrap (unit: dispute by default, or initiating state) --------------------------------------
## audit F5 (A2/A3 P2): a fit that errors OR warns (separation, non-convergence) is dropped for that model only; failures counted per model
fit_or_null <- function(d, terms) {
  warned <- FALSE
  res <- tryCatch(withCallingHandlers(coef(fit_cl(d, terms, robust = FALSE)), warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") }),
                  error = function(e) NULL)
  if (warned || is.null(res) || any(!is.finite(res))) NULL else res
}
set.seed(123)
ev_unit <- unique(L[, c("ev_id", if (BOOT_UNIT == "dispute") "disp" else "init")]); names(ev_unit) <- c("ev_id", "unit")
stopifnot(!anyDuplicated(ev_unit$ev_id))
units <- unique(ev_unit$unit); ev_by_unit <- split(ev_unit$ev_id, ev_unit$unit)
logf("bootstrap unit: %s (%d units)", BOOT_UNIT, length(units))
Lsplit <- split(L, L$ev_id)
bs <- list(M1 = NULL, M1s = NULL); bs_coef <- NULL; n_fail <- c(M1 = 0L, M1s = 0L, no_type = 0L)
fail_log <- data.frame(draw = integer(0), model = character(0), reason = character(0))   # GPT review P2-6: replicate IDs and failure reasons
for (bb in seq_len(B)) {
  draw <- sample(units, length(units), replace = TRUE)
  evs <- unlist(ev_by_unit[as.character(draw)], use.names = FALSE)
  Lb <- bind_rows(lapply(seq_along(evs), function(i) { x <- Lsplit[[as.character(evs[i])]]; x$ev_id <- i; x }))
  if (any(tabulate(Lb$type[Lb$chosen == 1L] + 1L, nbins = G) == 0L)) {     # a type with no events in this draw
    n_fail["no_type"] <- n_fail["no_type"] + 1L; fail_log <- rbind(fail_log, data.frame(draw = bb, model = "both", reason = "a type has no events")); next }
  c1 <- fit_or_null(Lb, GAP_G); c2 <- fit_or_null(Lb, TERMS_M1S)
  if (is.null(c1)) { n_fail["M1"] <- n_fail["M1"] + 1L; fail_log <- rbind(fail_log, data.frame(draw = bb, model = "M1", reason = "error, warning or non-finite fit"))
  } else { bs$M1 <- rbind(bs$M1, derive(c1)); bs_coef <- rbind(bs_coef, c(draw = bb, c1)) }
  if (is.null(c2)) { n_fail["M1s"] <- n_fail["M1s"] + 1L; fail_log <- rbind(fail_log, data.frame(draw = bb, model = "M1s", reason = "error, warning or non-finite fit"))
  } else bs$M1s <- rbind(bs$M1s, derive(c2))
  if (bb %% 100 == 0) logf("bootstrap %d / %d", bb, B)
}
logf("bootstrap: %d draws; failed or warned M1 %d, M1s %d; draws with a type missing %d (dropped)", B, n_fail["M1"], n_fail["M1s"], n_fail["no_type"])
## GPT review P2-6: intervals only when at least MIN_OK of the draws succeeded for that model
MIN_OK <- 0.90; MIN_B <- 200L                                    # joint audit B3: intervals need B >= 200 and >= 90% successful draws
n_ok <- c(M1 = if (is.null(bs$M1)) 0L else nrow(bs$M1), M1s = if (is.null(bs$M1s)) 0L else nrow(bs$M1s))
boot_ok <- if (B >= MIN_B) n_ok / B >= MIN_OK else c(M1 = FALSE, M1s = FALSE)
logf("bootstrap successful draws: M1 %d, M1s %d of %d; intervals reported for %s", n_ok[["M1"]], n_ok[["M1s"]], B,
     paste(names(boot_ok)[boot_ok], collapse = ", "))

## ---- 6b. leave-one-dispute-out and leave-one-initiator-out (joint audit M3; --loo=TRUE) ------------------------------
## Refit M1 without each dispute (and each initiating state) and record the primary (two-way) tests and the coefficient differences of each
## non-reference type from the reference type on dM and dR. Drops that leave a type without events are skipped and logged.
loo <- NULL
if (LOO) {
  ev_d <- unique(L[, c("ev_id", "disp", "init")])
  for (kind in c("dispute", "initiator")) {
    ids <- unique(if (kind == "dispute") ev_d$disp else ev_d$init)
    for (id in ids) {
      drop_ev <- if (kind == "dispute") ev_d$ev_id[ev_d$disp == id] else ev_d$ev_id[ev_d$init == id]
      Ld <- L[!L$ev_id %in% drop_ev, ]
      if (any(tabulate(Ld$type[Ld$chosen == 1L] + 1L, nbins = G) == 0L)) { logf("loo %s %s skipped: a type has no events", kind, id); next }
      md <- tryCatch(fit_cl(Ld, GAP_G, robust = FALSE), error = function(e) NULL)
      if (is.null(md) || any(!is.finite(coef(md)))) { logf("loo %s %s: fit failed", kind, id); next }
      Vd <- vc_set(md, Ld); nmd <- names(coef(md))
      Rd <- function(pairs) { R <- matrix(0, length(pairs), length(nmd), dimnames = list(NULL, nmd)); for (i in seq_along(pairs)) { R[i, pairs[[i]][1]] <- 1; R[i, pairs[[i]][2]] <- -1 }; R }
      w_s <- wald3(Rd(PAIRS_REF), md, Vd)
      o <- OTHERS[1]
      w_m <- wald3(Rd(lapply(c("dM", "dM2", "dMR"), function(v) c(paste0(v, "_", o), paste0(v, "_", REF)))), md, Vd)
      w_r <- wald3(Rd(lapply(c("dR", "dR2", "dMR"), function(v) c(paste0(v, "_", o), paste0(v, "_", REF)))), md, Vd)
      bd <- coef(md)
      loo <- rbind(loo, data.frame(drop = kind, id = id, events_dropped = length(drop_ev), surfaces_chi2 = w_s[["chi2"]], surfaces_p = w_s[["p"]],
                                   material_p = w_m[["p"]], relational_p = w_r[["p"]],
                                   dM_diff = bd[[paste0("dM_", o)]] - bd[[paste0("dM_", REF)]], dR_diff = bd[[paste0("dR_", o)]] - bd[[paste0("dR_", REF)]]))
    }
  }
  base <- c(surfaces_p = W_roles[["p"]], dM_diff = b1[[paste0("dM_", OTHERS[1])]] - b1[[paste0("dM_", REF)]])
  logf("leave-one-out: %d refits; surfaces p (primary) ranges %.3g to %.3g (full sample %.3g); %s - %s dM difference ranges %.3f to %.3f (full %.3f)",
       nrow(loo), min(loo$surfaces_p), max(loo$surfaces_p), base[["surfaces_p"]], OTHERS[1], REF, min(loo$dM_diff), max(loo$dM_diff), base[["dM_diff"]])
}

## ---- 7. tables (reporting rule F1) ------------------------------------------------------------------------
in_support <- function(p, g) { s <- SUPPORT[[g]]; !is.na(p[[paste0("muM_", g)]]) && p[[paste0("muM_", g)]] >= s[["dM_lo"]] && p[[paste0("muM_", g)]] <= s[["dM_hi"]] &&
                               p[[paste0("muR_", g)]] >= s[["dR_lo"]] && p[[paste0("muR_", g)]] <= s[["dR_hi"]] }
fmt_iv <- function(est, d, star = TRUE) {
  q <- quantile(d, c(0.025, 0.975), na.rm = TRUE)
  sprintf("%.2f [%.2f, %.2f]%s", est, q[1], q[2], ifelse(star & (q[1] > 0 | q[2] < 0), "*", ""))
}
tab <- NULL
for (mm in c("M1", "M1s")) {
  p <- pt[[mm]]; D <- if (boot_ok[[mm]]) bs[[mm]] else NULL        # too few successful draws: no intervals, no peak shares
  tab <- rbind(tab, data.frame(model = mm, quantity = "bootstrap_successful_draws", estimate = sprintf("%d of %d", n_ok[[mm]], B),
                               status = if (boot_ok[[mm]]) "intervals reported" else if (B < MIN_B) sprintf("B < %d: no intervals", MIN_B) else sprintf("fewer than %.0f%% succeeded: no intervals", 100 * MIN_OK), stringsAsFactors = FALSE))
  share <- setNames(vapply(ROLES, function(g) if (is.null(D)) NA_real_ else mean(D[, paste0("peak_", g)]), 0), ROLES)
  ok <- setNames(vapply(ROLES, function(g) p[[paste0("peak_", g)]] == 1 && !is.na(share[[g]]) && share[[g]] >= 0.95 && in_support(p, g), TRUE), ROLES)
  for (g in ROLES) {
    reason <- if (p[[paste0("peak_", g)]] != 1) "no peak at the point estimate" else if (is.na(share[[g]])) "no bootstrap" else if (share[[g]] < 0.95) "peak in < 95% of draws" else if (!in_support(p, g)) "ideal point outside the data" else ""
    tab <- rbind(tab, data.frame(model = mm, quantity = paste0("peak_share_", g), estimate = sprintf("%.3f", share[[g]]), status = if (ok[[g]]) "reported" else reason, stringsAsFactors = FALSE))
    for (q in QN) tab <- rbind(tab, data.frame(model = mm, quantity = paste0(q, "_", g),
                                               estimate = if (ok[[g]]) fmt_iv(p[[paste0(q, "_", g)]], D[, paste0(q, "_", g)], star = !q %in% c("sdM", "sdR")) else "no peak",
                                               status = if (ok[[g]]) "reported" else reason, stringsAsFactors = FALSE))
  }
  for (pr in PAIRS_R) for (q in QN) {
    both <- ok[[pr[1]]] && ok[[pr[2]]]; nm <- paste0(q, "_", pr[1], "-", pr[2])
    tab <- rbind(tab, data.frame(model = mm, quantity = nm, estimate = if (both) fmt_iv(p[[nm]], D[, nm]) else "not reported",
                                 status = if (both) "reported" else "a type without a clear peak", stringsAsFactors = FALSE))
  }
}
tab <- rbind(tab, data.frame(model = "M0", quantity = c(paste0(QN, "_pooled"), "peak_pooled"), estimate = sprintf("%.2f", pt$M0[c(paste0(QN, "_pooled"), "peak_pooled")]),
                             status = "point estimate only", stringsAsFactors = FALSE))
coefs <- bind_rows(lapply(c("M0", "M1", "M1s"), function(mm) { m <- get(mm); s <- summary(m)$coefficients
  data.frame(model = mm, term = rownames(s), coef = s[, "coef"], robust_se = s[, "robust se"], p = s[, ncol(s)], row.names = NULL) }))
tests <- data.frame(test = c(sprintf("utility surfaces differ by type (M1 vs M0, %d df)", W_roles[["df"]]),
                             sprintf("utility interaction (cross coefficient) differs by type (%d df; NOT equal rho)", W_cross[["df"]])),
                    chi2 = c(W_roles[["chi2"]], W_cross[["chi2"]]), df = c(W_roles[["df"]], W_cross[["df"]]),
                    p_init_x_dispute = c(W_roles[["p"]], W_cross[["p"]]), p_initiator = c(W_roles[["p_initiator"]], W_cross[["p_initiator"]]),
                    p_target_year = c(W_roles[["p_tgtyr"]], W_cross[["p_tgtyr"]]))
od <- file.path(v2_root(), "output"); dir.create(od, showWarnings = FALSE)
write.csv(tab, file.path(od, paste0(cell, "_derived.csv")), row.names = FALSE)
write.csv(coefs, file.path(od, paste0(cell, "_coefs.csv")), row.names = FALSE)
write.csv(tests, file.path(od, paste0(cell, "_tests.csv")), row.names = FALSE)
write.csv(pw, file.path(od, paste0(cell, "_pairwise.csv")), row.names = FALSE)
if (!is.null(loo)) write.csv(loo, file.path(od, paste0(cell, "_loo.csv")), row.names = FALSE)
saveRDS(list(pt = pt, bs = bs, bs_coef = bs_coef, fail_log = fail_log, n_ok = n_ok, n_fail = n_fail, B = B, spec = spec, types = TYPES, roles = ROLES,
             support = SUPPORT, n_rows = nrow(L), E = sum(keep), V = VS, boot_unit = BOOT_UNIT, loo = loo),
        file.path(od, paste0(cell, "_boot.rds")))
print(tab, right = FALSE, row.names = FALSE); print(tests, row.names = FALSE)
logf("done: %s_{derived,coefs,tests,pairwise}.csv and _boot.rds", cell)
