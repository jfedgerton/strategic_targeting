## 11_two_step_simple.R (2026-11-16, PI request): the simple two-step check of "selection versus failed bargaining", step B.
## Step A (selection, whom states confront) = 10_clogit_roles.R on the MID onsets of the interacting-dyads build:
##   Rscript v2/R/10_clogit_roles.R --sample=mid --events=dymid ... --types=bc
## Step B (failure, which confrontations turn fatal) = this script, on the SAME events (same assembly, same flags, same types): a weighted logit
## of the event's fatal flag on the chosen target's gaps, type-specific dM, dM^2, dR, dR^2, dM*dR, a type intercept and the chosen row's directed Z.
## Outcome: the event's in_<outcome> flag of events.rds (the dymid build's dyadic fatal coding; fmid = jointly fatal, the default), matched on
## (initiator, target, year) exactly as fail_outcomes() of the two-step branch does (every event found once, a mid event, not an escalation).
## Weights: the D8 event weights (dat$ev_w), as in step A. Family quasibinomial (same coefficients as binomial; non-integer weights allowed).
## Inference: cluster-robust variance by initiating state (sandwich::vcovCL, HC0, as step A's CR0); readings from a cluster bootstrap over
## initiating states (B draws, seed 123, percentile 95% intervals). Omitted relative to the Stan step B (two-step branch): the initiator random
## intercept and the joint posterior. This is the transparent baseline; the Stan two-step is the main estimate once its reading rules are fixed.
##
## Readings per type and axis (no labels, no thresholds; PI asked for the evidence, the decision table waits for the GPT review):
##   dP_<axis>_<type> = the average, over the type's events (each at its own covariates), of the predicted fatality probability with that
##   axis's gap set to the type's 90th percentile minus the same with it set to the 10th percentile (the other gap and Z as observed).
##   A positive value: confrontations against WEAKER (dM) or less central (dR) targets are more likely to turn fatal (gap = initiator minus target).
## Tests: robust Wald, gap terms differ by type (5 (G - 1) df), and per type whether the axis's terms (linear + squared, 2 df) are zero.
##
## Usage:
##   Rscript v2/R/11_two_step_simple.R --sample=mid --lag=3 --model=M1c --comparison=directed --events=dymid --variant=all --roles=bcj --joiner_set=all --types=bc
##   --outcome=fmid|fmidany|fmidub (default fmid), --reach=capdist, --boot=1000
suppressPackageStartupMessages({ library(dplyr); library(sandwich) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
R_DIR <- dirname(normalizePath(this_file))
source(file.path(R_DIR, "spec.R")); source(file.path(R_DIR, "helpers_data.R")); source(file.path(R_DIR, "helpers_model.R"))
source(file.path(R_DIR, "assemble_cell.R"))
Sys.setenv(V2_ROOT = normalizePath(file.path(R_DIR, "..")))

args <- commandArgs(trailingOnly = TRUE)
SAMPLE <- get_arg(args, "sample", "mid"); LAG <- as.integer(get_arg(args, "lag", 3)); MODEL <- get_arg(args, "model", "M1c")
B <- as.integer(get_arg(args, "boot", 1000)); TYPES <- get_arg(args, "types", "bc"); OUTCOME <- get_arg(args, "outcome", "fmid")
stopifnot(SAMPLE == "mid", LAG == 3L, MODEL == "M1c", B >= 0, TYPES %in% c("bc", "bc_nocoorig", "bcj"), OUTCOME %in% c("fmid", "fmidany", "fmidub"))
spec <- parse_spec(args)
stopifnot(spec_is_directed(spec), spec$roles == "bcj", spec_variant(spec) == "all", (if (is.null(spec$joiner_set)) "all" else spec$joiner_set) == "all",
          spec$controls == "minimal", spec$context == "all", (if (is.null(spec$mu_shift)) "none" else spec$mu_shift) == "none",
          spec_events(spec) == "dymid", spec$ww2 == "include")
cell <- paste0("stepB_", cell_name(SAMPLE, LAG, MODEL, spec), "_", TYPES, "_", OUTCOME)
log_open(paste0("11_", cell))
logf("cell %s; spec %s; types %s; outcome in_%s; bootstrap B = %d, seed 123", cell, spec_string(spec), TYPES, OUTCOME, B)

## ---- 1. assemble exactly as step A (10_clogit_roles.R) and 02_fit_model.R do ----------------------------------
data_dir <- spec_data_dir(spec, v2_root())
manifest <- read.csv(file.path(data_dir, "manifest.csv"), stringsAsFactors = FALSE)
A <- assemble_cell_spec(SAMPLE, LAG, spec, data_dir = data_dir, manifest = manifest, model = MODEL)
sy <- A$sy; alts <- A$alts; ev <- A$ev; Z_cols <- A$Z_cols
dat <- build_stan_data(sy, alts, ev, A$W_cols, A$X_cols, Z_cols, joint = 1L, sampled = A$sampled, n_full = A$n_full,
                       x_ref = A$x_ref, allow_orphan_initiations = TRUE)
gt <- A$group_table
stopifnot(A$sampled == 0L, identical(gt$grp, 0:2), identical(gt$group, c("bilateral", "coalition_originator", "joiner")),
          length(ev$grp) == dat$E, all(ev$grp %in% 0:2), identical(as.integer(dat$ev_sy), as.integer(ev$sy_row)))

## ---- 2. types (identical rule to 10_clogit_roles.R) ---------------------------------------------------------
if (TYPES == "bcj") { ROLES <- c("bil", "coal", "join"); t_ev <- ev$grp; keep <- rep(TRUE, dat$E)
} else if (TYPES == "bc") { ROLES <- c("bil", "coal"); t_ev <- as.integer(ev$grp > 0L); keep <- rep(TRUE, dat$E)
} else { ROLES <- c("bil", "coal"); t_ev <- as.integer(ev$grp > 0L); keep <- ev$grp != 1L }
G <- length(ROLES)

## ---- 3. outcome: the event's fatal flag in events.rds, matched on (initiator, target, year) ------------------------
keys <- data.frame(e = seq_len(dat$E), initiator = as.integer(sy$ccode[ev$sy_row]), target = as.integer(alts$ccode2[ev$alt_row]), year = as.integer(sy$year[ev$sy_row]))
stopifnot(!anyNA(keys), all(alts$ccode1[ev$alt_row] == keys$initiator), all(alts$year[ev$alt_row] == keys$year), all(keys$initiator != keys$target),
          !anyDuplicated(keys[, c("initiator", "target", "year")]))
events <- readRDS(file.path(data_dir, "events.rds"))
need <- c("initiator", "target", "year", "escalation", "in_mid", "in_fmid", "in_fmidany", "in_fmidub")
stopifnot(all(need %in% names(events)))
stopifnot("dispnum" %in% names(events))
ne <- events |> filter(!escalation) |> select(all_of(c(need[-4], "dispnum")))
stopifnot(!anyDuplicated(ne[, c("initiator", "target", "year")]))
j <- keys |> left_join(ne, by = c("initiator", "target", "year"))
if (anyNA(j$in_mid)) stop(sprintf("%d events have no non-escalation row in events.rds", sum(is.na(j$in_mid))))
stopifnot(nrow(j) == dat$E, identical(j$e, keys$e), all(j$in_mid), all(j$in_fmidany[j$in_fmid]), all(j$in_fmidub[j$in_fmid]))
y_all <- as.integer(j[[paste0("in_", OUTCOME)]])
stopifnot(all(y_all %in% 0:1))

## ---- 4. event frame: the chosen row's gaps and Z, type, weight, cluster -----------------------------------------
stopifnot(!anyNA(j$dispnum))
D <- data.frame(e = which(keep), y = y_all[keep], type = t_ev[keep], w = dat$ev_w[keep], init = keys$initiator[keep],
                disp = j$dispnum[keep], tgtyr = paste(keys$target[keep], keys$year[keep], sep = "_"))
D$dM <- dat$x[dat$ev_row[D$e], 1]; D$dR <- dat$x[dat$ev_row[D$e], 2]
Zm <- dat$Z[dat$ev_row[D$e], , drop = FALSE]; colnames(Zm) <- Z_cols; D <- cbind(D, Zm)
stopifnot(!anyNA(D))
GAP <- c("dM", "dM2", "dR", "dR2", "dMR")
D$dM2 <- D$dM^2; D$dR2 <- D$dR^2; D$dMR <- D$dM * D$dR
for (g in 0:(G - 1L)) { D[[paste0("int_", ROLES[g + 1])]] <- as.integer(D$type == g); for (v in GAP) D[[paste0(v, "_", ROLES[g + 1])]] <- D[[v]] * (D$type == g) }
GAP_G <- as.vector(outer(GAP, ROLES, paste, sep = "_")); INT <- paste0("int_", ROLES[-1])     # bilateral = reference intercept
cnt <- vapply(0:(G - 1L), function(g) c(n = sum(D$type == g), fatal = sum(D$y[D$type == g]), w = sum(D$w[D$type == g]), w_fatal = sum(D$w[D$type == g & D$y == 1L])), numeric(4))
logf("step-B events %d (of %d assembled; dropped %d); by type (events, fatal, weight, fatal weight): %s", nrow(D), dat$E, sum(!keep),
     paste(sprintf("%s %d / %d / %.1f / %.1f", ROLES, cnt[1, ], cnt[2, ], cnt[3, ], cnt[4, ]), collapse = "; "))
stopifnot(all(cnt[2, ] > 0), all(cnt[2, ] < cnt[1, ]))         # every type has fatal and non-fatal events

## ---- 5. fit --------------------------------------------------------------------------------------------------
FORM <- as.formula(paste("y ~", paste(c(INT, GAP_G, Z_cols), collapse = " + ")))
fitB <- function(d) glm(FORM, data = d, weights = w, family = quasibinomial())
m <- fitB(D)
stopifnot(m$converged, all(is.finite(coef(m))))
## cluster-robust variances, CR0 (no small-sample factor, matching step A's coxph robust variance; GPT review P3):
## by initiator (main), two-way initiator x dispute (sandwich multiway), and by target-year (sensitivity)
V <- vcovCL(m, cluster = D$init, type = "HC0", cadjust = FALSE)
V2w_raw <- vcovCL(m, cluster = D[, c("init", "disp")], type = "HC0", cadjust = FALSE, multi0 = FALSE)
neg_eig <- sum(eigen(V2w_raw, symmetric = TRUE, only.values = TRUE)$values <= 0)
## Cameron, Gelbach and Miller (2011): a two-way variance need not be positive definite; replace its negative eigenvalues with zero (sandwich fix = TRUE)
V2w <- vcovCL(m, cluster = D[, c("init", "disp")], type = "HC0", cadjust = FALSE, multi0 = FALSE, fix = TRUE)
Vty <- vcovCL(m, cluster = D$tgtyr, type = "HC0", cadjust = FALSE)
psd2w <- TRUE
logf("two-way (initiator x dispute) variance: %d non-positive eigenvalues of %d before the CGM fix%s", neg_eig, ncol(V2w_raw), if (neg_eig > 0) " (fixed: set to zero)" else "")
wald <- function(R, VV = V) { b <- coef(m); d <- R %*% b; s <- as.numeric(t(d) %*% solve(R %*% VV %*% t(R)) %*% d); c(chi2 = s, df = nrow(R), p = pchisq(s, nrow(R), lower.tail = FALSE)) }
wald3 <- function(R) { a <- wald(R); data.frame(chi2 = a[["chi2"]], df = a[["df"]], p = a[["p"]], p_twoway_init_x_dispute = if (psd2w) wald(R, V2w)[["p"]] else NA_real_,
                                                p_target_year = wald(R, Vty)[["p"]]) }
nmb <- names(coef(m))
Rsel <- function(pos, neg = NULL) { R <- matrix(0, length(pos), length(nmb), dimnames = list(NULL, nmb)); for (i in seq_along(pos)) { R[i, pos[i]] <- 1; if (!is.null(neg)) R[i, neg[i]] <- -1 }; R }
## GPT review 2026-11-16 P1: an axis effect is b_k + 2 b_kk x_k + b_MR x_other, so "no effect along an axis" restricts linear, squared AND cross (3 df).
## Failing to reject is not evidence of no effect (no equivalence test here).
AX <- list(material = c("dM", "dM2", "dMR"), relational = c("dR", "dR2", "dMR"))
tests <- NULL
if (G > 1) {
  for (blk in c(list(all = GAP), AX)) for (o in ROLES[-1]) {
    lab <- if (identical(blk, GAP)) sprintf("%s - %s: gap terms differ (%d df)", o, ROLES[1], length(blk)) else
           sprintf("%s - %s: %s axis effect differs anywhere (linear, squared, cross; %d df)", o, ROLES[1], if (blk[1] == "dM") "material" else "relational", length(blk))
    tests <- rbind(tests, cbind(test = lab, wald3(Rsel(paste0(blk, "_", o), paste0(blk, "_", ROLES[1]))))) }
}
for (g in ROLES) for (ax in names(AX)) {
  tests <- rbind(tests, cbind(test = sprintf("%s: %s gap has no effect on fatality anywhere (linear, squared, cross = 0; 3 df)", g, ax),
                              wald3(Rsel(paste0(AX[[ax]], "_", g)))))
}
logf("tests:\n%s", paste(capture.output(print(tests, row.names = FALSE)), collapse = "\n"))

## ---- 6. readings: average predicted change across each type's 10th-90th percentile of each gap ----------------------
q1090 <- lapply(0:(G - 1L), function(g) list(dM = quantile(D$dM[D$type == g], c(0.1, 0.9), names = FALSE), dR = quantile(D$dR[D$type == g], c(0.1, 0.9), names = FALSE)))
names(q1090) <- ROLES
dP <- function(mod, d) {
  out <- c()
  for (gi in seq_along(ROLES)) { g <- ROLES[gi]; s <- d$type == (gi - 1L)
    for (ax in c("dM", "dR")) {
      p_at <- function(val) { nd <- d[s, ]; nd[[ax]] <- val; nd$dM2 <- nd$dM^2; nd$dR2 <- nd$dR^2; nd$dMR <- nd$dM * nd$dR
        for (v in GAP) for (h in ROLES) nd[[paste0(v, "_", h)]] <- nd[[v]] * (nd$type == (match(h, ROLES) - 1L))
        weighted.mean(predict(mod, newdata = nd, type = "response"), nd$w) }
      out[paste0("dP_", ax, "_", g)] <- p_at(q1090[[g]][[ax]][2]) - p_at(q1090[[g]][[ax]][1])
    }
    out[paste0("p_fatal_", g)] <- weighted.mean(fitted(mod)[s], d$w[s])
  }
  if (G > 1) for (ax in c("dM", "dR")) out[paste0("dP_", ax, "_coal-bil")] <- out[paste0("dP_", ax, "_coal")] - out[paste0("dP_", ax, "_bil")]
  out
}
pt <- dP(m, D)
set.seed(123)
inits <- unique(D$init); by_init <- split(seq_len(nrow(D)), D$init)
bs <- NULL; bs_coef <- NULL; n_fail <- 0L
fail_log <- data.frame(draw = integer(0), reason = character(0))   # GPT review P2-6: replicate IDs and failure reasons
for (bb in seq_len(B)) {
  rows <- unlist(by_init[as.character(sample(inits, length(inits), replace = TRUE))], use.names = FALSE)
  Db <- D[rows, ]
  tt <- table(factor(Db$type, levels = 0:(G - 1L))); ft <- tapply(Db$y, factor(Db$type, levels = 0:(G - 1L)), sum)
  if (any(tt == 0L) || any(ft == 0L, na.rm = TRUE) || any(ft == tt, na.rm = TRUE)) {
    n_fail <- n_fail + 1L; fail_log <- rbind(fail_log, data.frame(draw = bb, reason = "a type lacks events or lacks fatal / non-fatal events")); next }
  warned <- FALSE
  mb <- tryCatch(withCallingHandlers(fitB(Db), warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") }), error = function(e) NULL)
  if (is.null(mb) || warned || !mb$converged || any(!is.finite(coef(mb)))) {
    n_fail <- n_fail + 1L; fail_log <- rbind(fail_log, data.frame(draw = bb, reason = if (is.null(mb)) "error" else if (warned) "warning (e.g. separation)" else "non-convergence or non-finite")); next }
  bs <- rbind(bs, dP(mb, Db)); bs_coef <- rbind(bs_coef, c(draw = bb, coef(mb)))
  if (bb %% 100 == 0) logf("bootstrap %d / %d", bb, B)
}
MIN_OK <- 0.90
n_ok <- if (is.null(bs)) 0L else nrow(bs)
boot_ok <- B > 0 && n_ok / B >= MIN_OK
logf("bootstrap: %d draws, %d successful, %d dropped; intervals %s", B, n_ok, n_fail, if (boot_ok) "reported" else sprintf("NOT reported (fewer than %.0f%% succeeded)", 100 * MIN_OK))
if (!boot_ok) bs <- NULL

## ---- 7. tables -------------------------------------------------------------------------------------------------
fmt <- function(est, d) { if (is.null(d)) return(sprintf("%.3f", est)); q <- quantile(d, c(0.025, 0.975), na.rm = TRUE)
  sprintf("%.3f [%.3f, %.3f]%s", est, q[1], q[2], ifelse(q[1] > 0 | q[2] < 0, "*", "")) }
rd <- data.frame(quantity = names(pt), estimate = vapply(names(pt), function(q) fmt(pt[[q]], if (is.null(bs)) NULL else bs[, q]), ""), row.names = NULL)
rd$range <- vapply(rd$quantity, function(q) { if (!grepl("^dP_d[MR]_(bil|coal|join)$", q)) return(""); p <- strsplit(q, "_")[[1]]; r <- q1090[[p[3]]][[p[2]]]; sprintf("%.2f to %.2f", r[1], r[2]) }, "")
## GPT review P1-2 / P2-4: these are NET ENDPOINT CHANGES for each type's own population and range; they cannot show a flat curve (an interior
## peak or offsetting interactions give zero), and the coal-bil difference mixes ranges, baselines and covariate distributions. The common-scale
## comparison is in 12_two_step_joint.R.
rd$note <- ifelse(grepl("^dP_d[MR]_(bil|coal|join)$", rd$quantity), "net endpoint change (10th to 90th percentile of the type's own gap); not a test of flatness",
           ifelse(grepl("-", rd$quantity), "difference of own-population endpoint changes (different ranges and populations); not a type comparison on a common scale", ""))
rd <- rbind(data.frame(quantity = "bootstrap_successful_draws", estimate = sprintf("%d of %d", n_ok, B), range = "",
                       note = if (boot_ok) "intervals reported" else sprintf("fewer than %.0f%% succeeded: no intervals", 100 * MIN_OK)), rd)
se <- sqrt(diag(V)); coefs <- data.frame(term = nmb, coef = coef(m), robust_se = se, p = 2 * pnorm(-abs(coef(m) / se)), row.names = NULL)
od <- file.path(v2_root(), "output"); dir.create(od, showWarnings = FALSE)
write.csv(rd, file.path(od, paste0(cell, "_readings.csv")), row.names = FALSE)
write.csv(coefs, file.path(od, paste0(cell, "_coefs.csv")), row.names = FALSE)
write.csv(tests, file.path(od, paste0(cell, "_tests.csv")), row.names = FALSE)
saveRDS(list(pt = pt, bs = bs, bs_coef = bs_coef, fail_log = fail_log, n_ok = n_ok, n_fail = n_fail, B = B, spec = spec, types = TYPES, outcome = OUTCOME,
             q1090 = q1090, counts = cnt, E = nrow(D), V = list(init = V, twoway = V2w, tgtyr = Vty, twoway_psd = psd2w)),
        file.path(od, paste0(cell, "_boot.rds")))
print(rd, row.names = FALSE); print(tests, row.names = FALSE)
logf("done: %s_{readings,coefs,tests}.csv and _boot.rds", cell)
