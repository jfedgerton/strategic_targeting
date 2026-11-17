## 11_two_step_simple.R (2026-11-16, PI request): the simple two-step check of "selection versus failed bargaining", step B.
## Step A (selection, whom states confront) = 10_clogit_roles.R on the MID onsets of the interacting-dyads build:
##   Rscript v2/R/10_clogit_roles.R --sample=mid --events=dymid ... --types=bc
## Step B (the conditional probability that a recorded dyadic dispute is jointly fatal) = this script, on the SAME events (same assembly, same flags, same types): a weighted logit
## of the event's fatal flag on the chosen target's gaps, type-specific dM, dM^2, dR, dR^2, dM*dR, a type intercept and the chosen row's directed Z.
## Outcome: the event's in_<outcome> flag of events.rds (the dymid build's dyadic fatal coding; fmid = jointly fatal, the default), matched on
## (initiator, target, year) exactly as fail_outcomes() of the two-step branch does (every event found once, a mid event, not an escalation).
## Weights: the D8 event weights (dat$ev_w), as in step A. Family quasibinomial (same coefficients as binomial; non-integer weights allowed).
## Inference (CR0, sandwich::vcovCL with cadjust = FALSE). PRIMARY clustering (PI decision 2026-11-16, joint audit M1): two-way initiator x dispute
## (sandwich multiway, fix = TRUE sets negative eigenvalues to zero); initiator-only and target-year as sensitivities. Readings from a pigeonhole
## two-way bootstrap by default (GPT recheck 2026-11-16, PI-approved: initiators and disputes resampled independently, each event's D8 weight
## times the product of its two multiplicities; --boot_unit=dispute or initiator are one-way sensitivities; the unit is part of the output
## name), seed 123, percentile 95% intervals, reported only when B >= 200
## and >= 90% of draws succeed. --loo=TRUE: leave-one-dispute-out and leave-one-initiator-out refits (M3). --firth=TRUE: a bias-reduced
## (Firth-type, brglm2 AS_mean) refit as a small-sample sensitivity (M4; its sandwich variance uses the ML score, an approximation).
## Omitted relative to the Stan step B (two-step branch): the initiator random intercept and the joint posterior.
##
## Readings per type and axis (no labels, no thresholds; PI asked for the evidence, the decision table waits for the GPT review):
##   dP_<axis>_<type> = the average, over the type's events (each at its own covariates), of the predicted fatality probability with that
##   axis's gap set to the type's 90th percentile minus the same with it set to the 10th percentile (the other gap and Z as observed).
##   A positive value: confrontations against WEAKER (dM) or less central (dR) targets are more likely to be jointly fatal (gap = initiator minus target).
## Tests: Wald, between types (5-df surface, 3-df per axis: linear, squared, cross) and per type whether an axis has no effect anywhere (3 df).
##
## Usage:
##   Rscript v2/R/11_two_step_simple.R --sample=mid --lag=3 --model=M1c --comparison=directed --events=dymid --variant=all --roles=bcj --joiner_set=all --types=bc
##   --outcome=fmid|fmidany|fmidub (default fmid), --reach=capdist, --boot=1000, --boot_unit=pigeonhole|dispute|initiator, --loo=TRUE, --firth=TRUE.
##   Unknown, retired or duplicated options are refused (joint audit B2).
suppressPackageStartupMessages({ library(dplyr); library(sandwich) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
R_DIR <- dirname(normalizePath(this_file))
source(file.path(R_DIR, "spec.R")); source(file.path(R_DIR, "helpers_data.R")); source(file.path(R_DIR, "helpers_model.R"))
source(file.path(R_DIR, "assemble_cell.R"))
Sys.setenv(V2_ROOT = normalizePath(file.path(R_DIR, "..")))

args <- commandArgs(trailingOnly = TRUE)
refuse_archived_options(args, "11_two_step_simple.R")
SCRIPT_FLAGS <- c("sample", "lag", "model", "boot", "types", "outcome", "boot_unit", "loo", "firth")
KNOWN_FLAGS <- c(names(SPEC_CHOICES), "ame_K", "tag", SCRIPT_FLAGS)
bad <- args[!grepl("^--[A-Za-z_0-9]+=.+$", args)]
if (length(bad)) stop(sprintf("arguments must be --name=value; refused: %s", paste(bad, collapse = " ")))
flag_names <- sub("^--([A-Za-z_0-9]+)=.*$", "\\1", args)
if (any(!flag_names %in% KNOWN_FLAGS)) stop(sprintf("unknown option(s) for 11_two_step_simple.R: %s", paste(unique(flag_names[!flag_names %in% KNOWN_FLAGS]), collapse = ", ")))
if (anyDuplicated(flag_names)) stop(sprintf("option(s) given more than once: %s", paste(unique(flag_names[duplicated(flag_names)]), collapse = ", ")))
SAMPLE <- get_arg(args, "sample", "mid"); LAG <- as.integer(get_arg(args, "lag", 3)); MODEL <- get_arg(args, "model", "M1c")
B <- as.integer(get_arg(args, "boot", 1000)); TYPES <- get_arg(args, "types", "bc"); OUTCOME <- get_arg(args, "outcome", "fmid")
BOOT_UNIT <- get_arg(args, "boot_unit", "pigeonhole"); LOO <- get_arg(args, "loo", FALSE); FIRTH <- get_arg(args, "firth", FALSE)
stopifnot(SAMPLE == "mid", LAG == 3L, MODEL == "M1c", B >= 0, TYPES %in% c("bc", "bc_nocoorig", "bcj"), OUTCOME %in% c("fmid", "fmidany", "fmidub"),
          BOOT_UNIT %in% c("pigeonhole", "dispute", "initiator"), is.logical(LOO), is.logical(FIRTH))
spec <- parse_spec(args)
stopifnot(spec_is_directed(spec), spec$roles == "bcj", spec_variant(spec) == "all", (if (is.null(spec$joiner_set)) "all" else spec$joiner_set) == "all",
          spec$controls == "minimal", spec$context == "all", (if (is.null(spec$mu_shift)) "none" else spec$mu_shift) == "none",
          spec_events(spec) == "dymid", spec$ww2 == "include")
cell <- paste0("stepB_", cell_name(SAMPLE, LAG, MODEL, spec), "_", TYPES, "_", OUTCOME, "_b", c(pigeonhole = "pig", dispute = "disp", initiator = "init")[[BOOT_UNIT]])   # GPT recheck P2
log_open(paste0("11_", cell))
logf("cell %s; spec %s; types %s; outcome in_%s; bootstrap B = %d over %ss, seed 123; leave-one-out %s; Firth %s", cell, spec_string(spec), TYPES, OUTCOME, B, BOOT_UNIT, LOO, FIRTH)

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
## cluster-robust variances, CR0 (no small-sample factor, matching step A): PRIMARY two-way initiator x dispute (Cameron, Gelbach and Miller 2011,
## negative eigenvalues set to zero, sandwich fix = TRUE); initiator-only and target-year as sensitivities
vset <- function(mod, d) {
  raw <- vcovCL(mod, cluster = d[, c("init", "disp")], type = "HC0", cadjust = FALSE, multi0 = FALSE)
  list(twoway = vcovCL(mod, cluster = d[, c("init", "disp")], type = "HC0", cadjust = FALSE, multi0 = FALSE, fix = TRUE),
       init = vcovCL(mod, cluster = d$init, type = "HC0", cadjust = FALSE), tgtyr = vcovCL(mod, cluster = d$tgtyr, type = "HC0", cadjust = FALSE),
       n_neg = sum(eigen(raw, symmetric = TRUE, only.values = TRUE)$values <= 0))
}
VS <- vset(m, D)
logf("two-way (initiator x dispute) variance: %d non-positive eigenvalues of %d before the CGM fix%s", VS$n_neg, ncol(VS$twoway), if (VS$n_neg > 0) " (fixed: set to zero)" else "")
wald <- function(R, VV, mod = m) { b <- coef(mod); d <- R %*% b; s <- as.numeric(t(d) %*% solve(R %*% VV %*% t(R)) %*% d); c(chi2 = s, df = nrow(R), p = pchisq(s, nrow(R), lower.tail = FALSE)) }
wald3 <- function(R, mod = m, V3 = VS) { b <- wald(R, V3$twoway, mod); data.frame(chi2 = b[["chi2"]], df = b[["df"]], p_init_x_dispute = b[["p"]],
                                                p_initiator = wald(R, V3$init, mod)[["p"]], p_target_year = wald(R, V3$tgtyr, mod)[["p"]]) }
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
u_init <- unique(D$init); u_disp <- unique(D$disp); i_idx <- match(D$init, u_init); d_idx <- match(D$disp, u_disp)
draw_mult <- function() {                                              # GPT recheck P2: pigeonhole two-way (default) or one-way multiplicities
  mi <- if (BOOT_UNIT == "dispute") rep(1L, length(u_init)) else tabulate(sample.int(length(u_init), length(u_init), replace = TRUE), length(u_init))
  md <- if (BOOT_UNIT == "initiator") rep(1L, length(u_disp)) else tabulate(sample.int(length(u_disp), length(u_disp), replace = TRUE), length(u_disp))
  mi[i_idx] * md[d_idx] }
logf("bootstrap: %s (%d initiators, %d disputes); largest initiator / dispute share of event weight %.3f / %.3f", BOOT_UNIT, length(u_init), length(u_disp),
     max(tapply(D$w, D$init, sum)) / sum(D$w), max(tapply(D$w, D$disp, sum)) / sum(D$w))
bs <- NULL; bs_coef <- NULL; n_fail <- 0L
fail_log <- data.frame(draw = integer(0), reason = character(0))   # GPT review P2-6: replicate IDs and failure reasons
for (bb in seq_len(B)) {
  mult <- draw_mult(); Db <- D[mult > 0L, ]; Db$w <- Db$w * mult[mult > 0L]
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
MIN_OK <- 0.90; MIN_B <- 200L                                     # joint audit B3
n_ok <- if (is.null(bs)) 0L else nrow(bs)
boot_ok <- B >= MIN_B && n_ok / B >= MIN_OK
logf("bootstrap: %d draws, %d successful, %d dropped; intervals %s", B, n_ok, n_fail,
     if (boot_ok) "reported" else if (B < MIN_B) sprintf("NOT reported (B < %d)", MIN_B) else sprintf("NOT reported (fewer than %.0f%% succeeded)", 100 * MIN_OK))
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
                       note = if (boot_ok) "intervals reported" else if (B < MIN_B) sprintf("B < %d: no intervals", MIN_B) else sprintf("fewer than %.0f%% succeeded: no intervals", 100 * MIN_OK)), rd)
coef_table <- function(mod, V3) { se <- sqrt(diag(V3$twoway)); se1 <- sqrt(diag(V3$init))
  data.frame(term = names(coef(mod)), coef = coef(mod), se_init_x_dispute = se, p_init_x_dispute = 2 * pnorm(-abs(coef(mod) / se)),
             se_initiator = se1, p_initiator = 2 * pnorm(-abs(coef(mod) / se1)), row.names = NULL) }
coefs <- coef_table(m, VS)

## ---- 7b. leave-one-dispute-out and leave-one-initiator-out (joint audit M3) and the Firth-type sensitivity (M4) ----------------------
loo <- NULL
if (LOO && G > 1) {
  o <- ROLES[2]; r0 <- ROLES[1]
  for (kind in c("dispute", "initiator")) for (id in unique(if (kind == "dispute") D$disp else D$init)) {
    drop <- if (kind == "dispute") D$disp == id else D$init == id
    Dd <- D[!drop, ]
    tt <- table(factor(Dd$type, levels = 0:(G - 1L))); ft <- tapply(Dd$y, factor(Dd$type, levels = 0:(G - 1L)), sum)
    if (any(tt == 0L) || any(ft == 0L, na.rm = TRUE) || any(ft == tt, na.rm = TRUE)) { logf("loo %s %s skipped: a type lacks fatal or non-fatal events", kind, id); next }
    md <- tryCatch(suppressWarnings(fitB(Dd)), error = function(e) NULL)
    if (is.null(md) || !md$converged || any(!is.finite(coef(md)))) { logf("loo %s %s: fit failed", kind, id); next }
    Vd <- vset(md, Dd); nmd <- names(coef(md))
    Rd <- function(pos, neg = NULL) { R <- matrix(0, length(pos), length(nmd), dimnames = list(NULL, nmd)); for (i in seq_along(pos)) { R[i, pos[i]] <- 1; if (!is.null(neg)) R[i, neg[i]] <- -1 }; R }
    loo <- rbind(loo, data.frame(drop = kind, id = id, events_dropped = sum(drop),
      surfaces_p = wald3(Rd(paste0(GAP, "_", o), paste0(GAP, "_", r0)), md, Vd)$p_init_x_dispute,
      material_diff_p = wald3(Rd(paste0(AX$material, "_", o), paste0(AX$material, "_", r0)), md, Vd)$p_init_x_dispute,
      coal_material_p = wald3(Rd(paste0(AX$material, "_", o)), md, Vd)$p_init_x_dispute,
      bil_material_p = wald3(Rd(paste0(AX$material, "_", r0)), md, Vd)$p_init_x_dispute,
      dM_diff = coef(md)[[paste0("dM_", o)]] - coef(md)[[paste0("dM_", r0)]]))
  }
  logf("leave-one-out: %d refits; coal-bil surfaces p (primary) ranges %.3g to %.3g; coal material p %.3g to %.3g", nrow(loo),
       min(loo$surfaces_p), max(loo$surfaces_p), min(loo$coal_material_p), max(loo$coal_material_p))
}
firth <- NULL
if (FIRTH) {
  if (!requireNamespace("brglm2", quietly = TRUE)) stop("--firth=TRUE needs the brglm2 package")
  ## brglmFit does not take quasi families; binomial gives the same mean model (the non-integer-weight warning is expected and muffled)
  mf <- suppressWarnings(glm(FORM, data = D, weights = w, family = binomial(), method = brglm2::brglmFit, type = "AS_mean"))
  stopifnot(mf$converged, all(is.finite(coef(mf))))
  VF <- vset(mf, D)
  ft <- NULL
  if (G > 1) for (blk in c(list(all = GAP), AX)) for (o in ROLES[-1])
    ft <- rbind(ft, cbind(test = sprintf("%s - %s: %s (Firth)", o, ROLES[1], paste(blk, collapse = "+")),
                          wald3(Rsel(paste0(blk, "_", o), paste0(blk, "_", ROLES[1])), mf, VF)))
  for (g in ROLES) for (ax in names(AX)) ft <- rbind(ft, cbind(test = sprintf("%s: %s axis no effect (Firth)", g, ax), wald3(Rsel(paste0(AX[[ax]], "_", g)), mf, VF)))
  firth <- list(coefs = coef_table(mf, VF), tests = ft)
  logf("Firth-type refit: max |coef change| %.3f", max(abs(coef(mf) - coef(m))))
}
od <- file.path(v2_root(), "output"); dir.create(od, showWarnings = FALSE)
write.csv(rd, file.path(od, paste0(cell, "_readings.csv")), row.names = FALSE)
write.csv(coefs, file.path(od, paste0(cell, "_coefs.csv")), row.names = FALSE)
write.csv(tests, file.path(od, paste0(cell, "_tests.csv")), row.names = FALSE)
if (!is.null(loo)) write.csv(loo, file.path(od, paste0(cell, "_loo.csv")), row.names = FALSE)
if (!is.null(firth)) { write.csv(firth$coefs, file.path(od, paste0(cell, "_firth_coefs.csv")), row.names = FALSE)
  write.csv(firth$tests, file.path(od, paste0(cell, "_firth_tests.csv")), row.names = FALSE) }
saveRDS(list(pt = pt, bs = bs, bs_coef = bs_coef, fail_log = fail_log, n_ok = n_ok, n_fail = n_fail, B = B, spec = spec, types = TYPES, outcome = OUTCOME,
             q1090 = q1090, counts = cnt, E = nrow(D), V = VS, boot_unit = BOOT_UNIT, loo = loo, firth = firth),
        file.path(od, paste0(cell, "_boot.rds")))
written <- c(paste0(cell, c("_readings.csv", "_coefs.csv", "_tests.csv", "_boot.rds")), if (!is.null(loo)) paste0(cell, "_loo.csv"),
             if (!is.null(firth)) paste0(cell, c("_firth_coefs.csv", "_firth_tests.csv")))   # GPT recheck P2: manifest of THIS run's artifacts
gitc <- tryCatch(system2("git", c("-C", R_DIR, "rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE), error = function(e) NA_character_, warning = function(w) NA_character_)
write.csv(data.frame(file = written, run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), git = if (length(gitc)) gitc[1] else NA_character_,
                     args = paste(args, collapse = " ")), file.path(od, paste0(cell, "_manifest.csv")), row.names = FALSE)
print(rd, row.names = FALSE); print(tests, row.names = FALSE)
logf("done: %s_{readings,coefs,tests}.csv and _boot.rds", cell)
