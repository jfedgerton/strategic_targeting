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
## Tests (robust Wald, cluster = initiator, CR0): M1 vs M0 (5 (G - 1) df) = "the utility surfaces differ by type" (valid);
##   equal cross coefficients (G - 1 df) = "the utility interaction of the two gaps differs by type". This is NOT a test of equal rho
##   (audit A3 P1-5: equal b5 does not imply equal rho); the per-axis blocks test (P mu)_k and P_kk jointly, not ideal points or tolerances.
## Intervals: cluster bootstrap over initiating states, B draws, seed 123, percentile. A draw whose fit errors or warns (e.g. separation)
## is dropped for that model only; failures are counted per model.
##
## Usage (same flags as the grid cell; --roles=bcj is required because the types are recoded from it):
##   Rscript v2/R/10_clogit_roles.R --sample=fmid --lag=3 --model=M1c --comparison=directed --entry=later --variant=all --roles=bcj --joiner_set=all --types=bc
##   --events=dymid (instead of --entry=later) for jointly fatal, --reach=capdist for the reach terms, --boot=200 for a shorter run.
suppressPackageStartupMessages({ library(dplyr); library(survival) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
R_DIR <- dirname(normalizePath(this_file))
source(file.path(R_DIR, "spec.R")); source(file.path(R_DIR, "helpers_data.R")); source(file.path(R_DIR, "helpers_model.R"))
source(file.path(R_DIR, "assemble_cell.R"))
Sys.setenv(V2_ROOT = normalizePath(file.path(R_DIR, "..")))   # audit F5 (A1 P2-3): read and write this checkout only, as 02_fit_model.R does

args <- commandArgs(trailingOnly = TRUE)
SAMPLE <- get_arg(args, "sample", "fmid"); LAG <- as.integer(get_arg(args, "lag", 3)); MODEL <- get_arg(args, "model", "M1c")
B <- as.integer(get_arg(args, "boot", 1000))
TYPES <- get_arg(args, "types", "bc")
stopifnot(SAMPLE %in% ALL_SAMPLES, LAG == 3L, MODEL == "M1c", B >= 0, TYPES %in% c("bc", "bc_nocoorig", "bcj"))
spec <- parse_spec(args)
## audit F5 (A1 P2-1): the cells M1c is defined on, and nothing else (02_fit_model.R refuses the same)
stopifnot(spec_is_directed(spec), spec$roles == "bcj", spec_variant(spec) == "all", (if (is.null(spec$joiner_set)) "all" else spec$joiner_set) == "all",
          spec$controls == "minimal", spec$context == "all", (if (is.null(spec$mu_shift)) "none" else spec$mu_shift) == "none")
cell <- paste0("clogit_", cell_name(SAMPLE, LAG, MODEL, spec), "_", TYPES)
log_open(paste0("10_", cell))
logf("cell %s; spec %s; types %s; bootstrap B = %d, seed 123", cell, spec_string(spec), TYPES, B)

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

wald <- function(m, R, r = rep(0, nrow(R))) {            # robust Wald test of R beta = r
  b <- coef(m); V <- vcov(m); d <- R %*% b - r
  s <- as.numeric(t(d) %*% solve(R %*% V %*% t(R)) %*% d); c(chi2 = s, df = nrow(R), p = pchisq(s, nrow(R), lower.tail = FALSE))
}
nm1 <- names(coef(M1))
Rm <- function(pairs) { R <- matrix(0, length(pairs), length(nm1), dimnames = list(NULL, nm1)); for (i in seq_along(pairs)) { R[i, pairs[[i]][1]] <- 1; R[i, pairs[[i]][2]] <- -1 }; R }
REF <- ROLES[1]; OTHERS <- ROLES[-1]
PAIRS_REF <- unlist(lapply(GAP, function(v) lapply(OTHERS, function(o) c(paste0(v, "_", o), paste0(v, "_", REF)))), recursive = FALSE)
stopifnot(length(PAIRS_REF) == 5L * (G - 1L))
W_roles <- wald(M1, Rm(PAIRS_REF))
W_cross <- wald(M1, Rm(lapply(OTHERS, function(o) c(paste0("dMR_", o), paste0("dMR_", REF)))))
logf("Wald (robust, cluster = initiator): utility surfaces differ by type (M1 vs M0) chi2 %.2f df %d p %.4g; utility interaction (cross coefficient) differs by type chi2 %.2f df %d p %.4g (NOT a test of equal rho)",
     W_roles[1], W_roles[2], W_roles[3], W_cross[1], W_cross[2], W_cross[3])
## pairwise contrasts, term by term: difference of the type coefficients, robust SE, 95% interval, p; plus per-axis 2-df and whole-surface 5-df tests
V1 <- vcov(M1); b1 <- coef(M1)
PAIRS_R <- combn(ROLES, 2, simplify = FALSE); PAIRS_R <- lapply(PAIRS_R, rev)          # (later type, earlier type), e.g. coal - bil
pw <- NULL
for (pr in PAIRS_R) {
  for (v in GAP) {
    a <- paste0(v, "_", pr[1]); z <- paste0(v, "_", pr[2]); d <- b1[[a]] - b1[[z]]; se <- sqrt(V1[a, a] + V1[z, z] - 2 * V1[a, z])
    pw <- rbind(pw, data.frame(pair = paste(pr, collapse = "-"), term = v, diff = d, se = se, lo = d - 1.96 * se, hi = d + 1.96 * se, p = 2 * pnorm(-abs(d / se)), df = 1))
  }
  for (blk in list(c("dM", "dM2"), c("dR", "dR2"), GAP)) {
    w2 <- wald(M1, Rm(lapply(blk, function(v) c(paste0(v, "_", pr[1]), paste0(v, "_", pr[2])))))
    pw <- rbind(pw, data.frame(pair = paste(pr, collapse = "-"), term = paste0("joint:", paste(blk, collapse = "+")), diff = NA, se = NA, lo = NA, hi = NA, p = w2[["p"]], df = w2[["df"]]))
  }
}
pw$note <- ifelse(pw$term == "dMR", "utility interaction P12, not rho", ifelse(grepl("^joint:dM\\+dM2$|^joint:dR\\+dR2$", pw$term), "tests (P mu)_k and P_kk jointly, not mu or sigma", ""))

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

## ---- 6. cluster bootstrap over initiating states ---------------------------------------------------------
## audit F5 (A2/A3 P2): a fit that errors OR warns (separation, non-convergence) is dropped for that model only; failures counted per model
fit_or_null <- function(d, terms) {
  warned <- FALSE
  res <- tryCatch(withCallingHandlers(coef(fit_cl(d, terms, robust = FALSE)), warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") }),
                  error = function(e) NULL)
  if (warned || is.null(res) || any(!is.finite(res))) NULL else res
}
set.seed(123)
ev_init <- unique(L[, c("ev_id", "init")]); inits <- unique(ev_init$init); ev_by_init <- split(ev_init$ev_id, ev_init$init)
Lsplit <- split(L, L$ev_id)
bs <- list(M1 = NULL, M1s = NULL); n_fail <- c(M1 = 0L, M1s = 0L, no_type = 0L)
for (bb in seq_len(B)) {
  draw <- sample(inits, length(inits), replace = TRUE)
  evs <- unlist(ev_by_init[as.character(draw)], use.names = FALSE)
  Lb <- bind_rows(lapply(seq_along(evs), function(i) { x <- Lsplit[[as.character(evs[i])]]; x$ev_id <- i; x }))
  if (any(tabulate(Lb$type[Lb$chosen == 1L] + 1L, nbins = G) == 0L)) { n_fail["no_type"] <- n_fail["no_type"] + 1L; next }   # a type with no events in this draw
  c1 <- fit_or_null(Lb, GAP_G); c2 <- fit_or_null(Lb, TERMS_M1S)
  if (is.null(c1)) n_fail["M1"] <- n_fail["M1"] + 1L else bs$M1 <- rbind(bs$M1, derive(c1))
  if (is.null(c2)) n_fail["M1s"] <- n_fail["M1s"] + 1L else bs$M1s <- rbind(bs$M1s, derive(c2))
  if (bb %% 100 == 0) logf("bootstrap %d / %d", bb, B)
}
logf("bootstrap: %d draws; failed or warned M1 %d, M1s %d; draws with a type missing %d (dropped)", B, n_fail["M1"], n_fail["M1s"], n_fail["no_type"])

## ---- 7. tables (reporting rule F1) ------------------------------------------------------------------------
in_support <- function(p, g) { s <- SUPPORT[[g]]; !is.na(p[[paste0("muM_", g)]]) && p[[paste0("muM_", g)]] >= s[["dM_lo"]] && p[[paste0("muM_", g)]] <= s[["dM_hi"]] &&
                               p[[paste0("muR_", g)]] >= s[["dR_lo"]] && p[[paste0("muR_", g)]] <= s[["dR_hi"]] }
fmt_iv <- function(est, d, star = TRUE) {
  q <- quantile(d, c(0.025, 0.975), na.rm = TRUE)
  sprintf("%.2f [%.2f, %.2f]%s", est, q[1], q[2], ifelse(star & (q[1] > 0 | q[2] < 0), "*", ""))
}
tab <- NULL
for (mm in c("M1", "M1s")) {
  p <- pt[[mm]]; D <- bs[[mm]]
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
                    chi2 = c(W_roles[1], W_cross[1]), df = c(W_roles[2], W_cross[2]), p = c(W_roles[3], W_cross[3]))
od <- file.path(v2_root(), "output"); dir.create(od, showWarnings = FALSE)
write.csv(tab, file.path(od, paste0(cell, "_derived.csv")), row.names = FALSE)
write.csv(coefs, file.path(od, paste0(cell, "_coefs.csv")), row.names = FALSE)
write.csv(tests, file.path(od, paste0(cell, "_tests.csv")), row.names = FALSE)
write.csv(pw, file.path(od, paste0(cell, "_pairwise.csv")), row.names = FALSE)
saveRDS(list(pt = pt, bs = bs, n_fail = n_fail, B = B, spec = spec, types = TYPES, roles = ROLES, support = SUPPORT, n_rows = nrow(L), E = sum(keep)),
        file.path(od, paste0(cell, "_boot.rds")))
print(tab, right = FALSE, row.names = FALSE); print(tests, row.names = FALSE)
logf("done: %s_{derived,coefs,tests,pairwise}.csv and _boot.rds", cell)
