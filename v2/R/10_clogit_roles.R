## 10_clogit_roles.R (2026-11-16, PI request): the simplest test of "target selection changes by initiation type as a joint consideration".
## Stage 2 only, as a weighted conditional logit (survival::coxph, Breslow = exact conditional logit with one chosen row per stratum),
## on EXACTLY the events, choice sets, gaps (dM, dR), directed Z and 1/n event weights that M1c sees: the cell is assembled with
## assemble_cell_spec() and build_stan_data() under the same flags as the grid (directed, roles bcj, joiner_set all, full choice set).
## Omitted relative to M1c: stage 1 (initiation, tau), the initiator ideal-point intercept a_i and the target effects b_j.
##
## The spatial utility -1/2 (x - mu)' P (x - mu) = b'x + 1/2 x'Hx + const with H = -P and b = P mu, so for each role g
##   u = b1 dM + b2 dM^2 + b3 dR + b4 dR^2 + b5 dM*dR    ->    H = [[2 b2, b5], [b5, 2 b4]],  Sigma = -H^{-1},  mu = Sigma b.
## Sigma is in the same units as M1c's Sigma_g (logit error scale 1), so mu, sigma and rho are directly comparable. A draw whose H is not
## negative definite has no ideal point (the surface has no peak); those draws are counted and excluded from the derived quantities.
##
## Models (all with the directed Z and strata = event):
##   M0  pooled: the five gap terms shared by the three roles
##   M1  role-specific: the five gap terms separately for bilateral (0), coalition originator (1), joiner (2)   <- the main test
##   M1s role-specific linear and squared terms, ONE cross term (the clogit analogue of M1c's shared rho)
## Tests: robust Wald (cluster = initiator) for M1 vs M0 (10 df) and for equal cross terms in M1 (2 df). Intervals for mu, sigma, rho and
## the role contrasts: cluster bootstrap over initiating states, B = 1000, seed 123, 95% percentile intervals.
##
## Usage (same flags as the grid cell):
##   Rscript v2/R/10_clogit_roles.R --sample=fmid --lag=3 --model=M1c --comparison=directed --entry=later --variant=all --roles=bcj --joiner_set=all
##   add --events=dymid (instead of --entry=later) for jointly fatal, --reach=capdist for the reach terms, --boot=200 for a shorter run.
suppressPackageStartupMessages({ library(dplyr); library(survival) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
R_DIR <- dirname(normalizePath(this_file))
source(file.path(R_DIR, "spec.R")); source(file.path(R_DIR, "helpers_data.R")); source(file.path(R_DIR, "helpers_model.R"))
source(file.path(R_DIR, "assemble_cell.R"))

args <- commandArgs(trailingOnly = TRUE)
SAMPLE <- get_arg(args, "sample", "fmid"); LAG <- as.integer(get_arg(args, "lag", 3)); MODEL <- get_arg(args, "model", "M1c")
B <- as.integer(get_arg(args, "boot", 1000))
stopifnot(SAMPLE %in% ALL_SAMPLES, LAG == 3L, MODEL == "M1c", B >= 0)
spec <- parse_spec(args)
stopifnot(spec_is_directed(spec), spec$roles == "bcj", spec_variant(spec) == "all", (if (is.null(spec$joiner_set)) "all" else spec$joiner_set) == "all")
cell <- paste0("clogit_", cell_name(SAMPLE, LAG, MODEL, spec))
log_open(paste0("10_", cell))
logf("cell %s; spec %s; bootstrap B = %d, seed 123", cell, spec_string(spec), B)

## ---- 1. assemble exactly as 02_fit_model.R does for M1c ----------------------------------------------------
data_dir <- spec_data_dir(spec, v2_root())
manifest <- read.csv(file.path(data_dir, "manifest.csv"), stringsAsFactors = FALSE)
A <- assemble_cell_spec(SAMPLE, LAG, spec, data_dir = data_dir, manifest = manifest, model = MODEL)
sy <- A$sy; alts <- A$alts; ev <- A$ev; Z_cols <- A$Z_cols
dat <- build_stan_data(sy, alts, ev, A$W_cols, A$X_cols, Z_cols, joint = 1L, sampled = A$sampled, n_full = A$n_full,
                       x_ref = A$x_ref, allow_orphan_initiations = TRUE)
stopifnot(A$sampled == 0L, "grp" %in% names(ev), length(ev$grp) == dat$E, all(ev$grp %in% 0:2), identical(as.integer(dat$ev_sy), as.integer(ev$sy_row)))
logf("events %d (bilateral %d, coalition originator %d, joiner %d); weight sums %.1f / %.1f / %.1f; Z = %s",
     dat$E, sum(ev$grp == 0), sum(ev$grp == 1), sum(ev$grp == 2), sum(ev$w[ev$grp == 0]), sum(ev$w[ev$grp == 1]), sum(ev$w[ev$grp == 2]), paste(Z_cols, collapse = ","))

## ---- 2. long frame: one stratum per event = its state-year's stage-2 set (in_s2 = 1), the chosen row flagged ---------
init_cc <- A$state_codes[sy$state_idx]                  # initiating state's ccode per state-year (cluster id)
rows <- lapply(seq_len(dat$E), function(e) {
  r <- dat$ev_start[e]:dat$ev_end[e]; r <- r[dat$in_s2[r] == 1L]
  data.frame(ev_id = e, row = r, chosen = as.integer(r == dat$ev_row[e]), grp = ev$grp[e], w = dat$ev_w[e], init = init_cc[dat$ev_sy[e]])
})
L <- bind_rows(rows)
stopifnot(all(tapply(L$chosen, L$ev_id, sum) == 1L), all(tapply(L$row, L$ev_id, length) >= 2L))
L$dM <- dat$x[L$row, 1]; L$dR <- dat$x[L$row, 2]
Zm <- dat$Z[L$row, , drop = FALSE]; colnames(Zm) <- Z_cols; L <- cbind(L, Zm)
stopifnot(!anyNA(L[, c("dM", "dR", Z_cols)]))
## gap terms, pooled and by role
L$dM2 <- L$dM^2; L$dR2 <- L$dR^2; L$dMR <- L$dM * L$dR
GAP <- c("dM", "dM2", "dR", "dR2", "dMR"); ROLES <- c("bil", "coal", "join")
for (g in 0:2) for (v in GAP) L[[paste0(v, "_", ROLES[g + 1])]] <- L[[v]] * (L$grp == g)
GAP_G <- as.vector(outer(GAP, ROLES, paste, sep = "_"))           # dM_bil, dM2_bil, ..., dMR_join
logf("long frame: %d rows, %d strata, median set size %.0f", nrow(L), dat$E, median(table(L$ev_id)))

## ---- 3. fits -----------------------------------------------------------------------------------------------
fit_cl <- function(d, terms, robust = TRUE) {
  f <- as.formula(paste("Surv(rep(1, nrow(d)), chosen) ~", paste(c(terms, Z_cols), collapse = " + "), "+ strata(ev_id)"))
  if (robust) coxph(f, data = d, weights = w, method = "breslow", cluster = init)
  else coxph(f, data = d, weights = w, method = "breslow")
}
M0 <- fit_cl(L, GAP)
M1 <- fit_cl(L, GAP_G)
M1s <- fit_cl(L, c(setdiff(GAP_G, paste0("dMR_", ROLES)), "dMR"))
for (m in list(M0, M1, M1s)) stopifnot(all(is.finite(coef(m))), !anyNA(coef(m)))

wald <- function(m, R, r = rep(0, nrow(R))) {            # robust Wald test of R beta = r
  b <- coef(m); V <- vcov(m); d <- R %*% b - r
  s <- as.numeric(t(d) %*% solve(R %*% V %*% t(R)) %*% d); c(chi2 = s, df = nrow(R), p = pchisq(s, nrow(R), lower.tail = FALSE))
}
nm1 <- names(coef(M1))
Rm <- function(pairs) { R <- matrix(0, length(pairs), length(nm1), dimnames = list(NULL, nm1)); for (i in seq_along(pairs)) { R[i, pairs[[i]][1]] <- 1; R[i, pairs[[i]][2]] <- -1 }; R }
W_roles <- wald(M1, Rm(unlist(lapply(GAP, function(v) list(c(paste0(v, "_coal"), paste0(v, "_bil")), c(paste0(v, "_join"), paste0(v, "_bil")))), recursive = FALSE)))
W_cross <- wald(M1, Rm(list(c("dMR_coal", "dMR_bil"), c("dMR_join", "dMR_bil"))))
logf("Wald (robust, cluster = initiator): role-specific gap terms vs pooled chi2 %.2f df %d p %.4f; equal cross terms chi2 %.2f df %d p %.4f",
     W_roles[1], W_roles[2], W_roles[3], W_cross[1], W_cross[2], W_cross[3])

## ---- 4. derived quantities ------------------------------------------------------------------------------------
derive <- function(b) {                                   # b: named coefficient vector; returns mu, sigma, rho per role (NA if no peak)
  out <- c()
  for (g in ROLES) {
    gb <- function(v) if (paste0(v, "_", g) %in% names(b)) b[[paste0(v, "_", g)]] else b[[v]]   # M1s: shared dMR; M0: pooled names
    H <- matrix(c(2 * gb("dM2"), gb("dMR"), gb("dMR"), 2 * gb("dR2")), 2)
    ev_H <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    if (all(ev_H < 0)) {
      S <- -solve(H); mu <- as.vector(S %*% c(gb("dM"), gb("dR")))
      out <- c(out, setNames(c(mu, sqrt(diag(S)), S[1, 2] / sqrt(S[1, 1] * S[2, 2]), 1), paste0(c("muM", "muR", "sdM", "sdR", "rho", "peak"), "_", g)))
    } else out <- c(out, setNames(c(rep(NA, 5), 0), paste0(c("muM", "muR", "sdM", "sdR", "rho", "peak"), "_", g)))
  }
  for (q in c("muM", "muR", "sdM", "sdR", "rho")) {
    out[paste0(q, "_coal-bil")] <- out[paste0(q, "_coal")] - out[paste0(q, "_bil")]
    out[paste0(q, "_join-bil")] <- out[paste0(q, "_join")] - out[paste0(q, "_bil")]
    out[paste0(q, "_join-coal")] <- out[paste0(q, "_join")] - out[paste0(q, "_coal")]
  }
  out
}
pt <- list(M0 = derive(c(setNames(coef(M0)[GAP], paste0(GAP, "")))), M1 = derive(coef(M1)), M1s = derive(coef(M1s)))

## ---- 5. cluster bootstrap over initiating states ---------------------------------------------------------
set.seed(123)
inits <- unique(L$init); ev_by_init <- split(unique(L[, c("ev_id", "init")])$ev_id, unique(L[, c("ev_id", "init")])$init)
Lsplit <- split(L, L$ev_id)
bs <- list(M1 = NULL, M1s = NULL); n_fail <- 0L
for (bb in seq_len(B)) {
  draw <- sample(inits, length(inits), replace = TRUE)
  evs <- unlist(ev_by_init[as.character(draw)], use.names = FALSE)
  Lb <- bind_rows(lapply(seq_along(evs), function(i) { x <- Lsplit[[as.character(evs[i])]]; x$ev_id <- i; x }))
  r1 <- tryCatch(derive(coef(fit_cl(Lb, GAP_G, robust = FALSE))), error = function(e) NULL)
  r2 <- tryCatch(derive(coef(fit_cl(Lb, c(setdiff(GAP_G, paste0("dMR_", ROLES)), "dMR"), robust = FALSE))), error = function(e) NULL)
  if (is.null(r1) || is.null(r2)) { n_fail <- n_fail + 1L; next }
  bs$M1 <- rbind(bs$M1, r1); bs$M1s <- rbind(bs$M1s, r2)
  if (bb %% 100 == 0) logf("bootstrap %d / %d", bb, B)
}
logf("bootstrap: %d of %d replicates failed to fit (dropped)", n_fail, B)

## ---- 6. tables ------------------------------------------------------------------------------------------
fmt <- function(est, d) {
  if (is.null(d) || all(is.na(d))) return(sprintf("%.2f", est))
  q <- quantile(d, c(0.025, 0.975), na.rm = TRUE)
  sprintf("%.2f [%.2f, %.2f]%s", est, q[1], q[2], ifelse(q[1] > 0 | q[2] < 0, "*", ""))
}
tab <- NULL
for (mm in c("M1", "M1s")) {
  for (q in setdiff(names(pt[[mm]]), grep("^peak", names(pt[[mm]]), value = TRUE))) {
    d <- if (!is.null(bs[[mm]])) bs[[mm]][, q] else NULL
    tab <- rbind(tab, data.frame(model = mm, quantity = q, estimate = fmt(pt[[mm]][[q]], d), stringsAsFactors = FALSE))
  }
  for (g in ROLES) tab <- rbind(tab, data.frame(model = mm, quantity = paste0("share_draws_with_peak_", g),
                                                estimate = if (!is.null(bs[[mm]])) sprintf("%.3f", mean(bs[[mm]][, paste0("peak_", g)])) else "", stringsAsFactors = FALSE))
}
tab <- rbind(tab, data.frame(model = "M0", quantity = names(pt$M0)[1:6], estimate = sprintf("%.2f", pt$M0[1:6]), stringsAsFactors = FALSE))
coefs <- bind_rows(lapply(c("M0", "M1", "M1s"), function(mm) { m <- get(mm); s <- summary(m)$coefficients
  data.frame(model = mm, term = rownames(s), coef = s[, "coef"], robust_se = s[, "robust se"], p = s[, ncol(s)], row.names = NULL) }))
tests <- data.frame(test = c("role-specific gap terms vs pooled (M1 vs M0)", "equal cross terms across roles (M1)"),
                    chi2 = c(W_roles[1], W_cross[1]), df = c(W_roles[2], W_cross[2]), p = c(W_roles[3], W_cross[3]))
od <- file.path(v2_root(), "output"); dir.create(od, showWarnings = FALSE)
write.csv(tab, file.path(od, paste0(cell, "_derived.csv")), row.names = FALSE)
write.csv(coefs, file.path(od, paste0(cell, "_coefs.csv")), row.names = FALSE)
write.csv(tests, file.path(od, paste0(cell, "_tests.csv")), row.names = FALSE)
saveRDS(list(pt = pt, bs = bs, n_fail = n_fail, B = B, spec = spec, n_rows = nrow(L), E = dat$E), file.path(od, paste0(cell, "_boot.rds")))
print(tab, right = FALSE, row.names = FALSE); print(tests, row.names = FALSE)
logf("done: %s_{derived,coefs,tests}.csv and _boot.rds", cell)
