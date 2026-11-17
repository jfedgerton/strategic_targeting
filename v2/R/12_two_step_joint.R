## 12_two_step_joint.R (2026-11-16, PI-approved after the GPT review of the simple two-step): selection and fatality on ONE scale.
## Step A (selection) and step B (fatality given confrontation) are fitted on the SAME MID onsets of the interacting-dyads build, with the same
## assembly, types, gaps, directed Z and D8 weights as 10_clogit_roles.R (step A) and 11_two_step_simple.R (step B), and resampled TOGETHER
## (one cluster bootstrap drives both fits), so the uncertainty of the decomposition is joint.
##
## Decomposition (GPT review). For two candidate targets j and k that face the same initiator, year, type and choice set, the model-implied
## odds of a FATAL confrontation with j rather than k split additively:
##     log P_F(j)/P_F(k) = log P_A(j)/P_A(k)  +  log r_B(j)/r_B(k)
##                         [selection]           [fatality: a log RISK ratio, not the step-B logit coefficient]
## with P_A the step-A choice probability and r_B the step-B probability that the confrontation turns (jointly) fatal.
## Here j is a candidate with gaps x = (dM, dR) and k a reference candidate x0; the two candidates share everything else (Z), so in step A the
## Z terms cancel and selection_g(x) = u_g(x) - u_g(x0) exactly. Step B's risk depends on the target's Z and the type intercept, so the fatality
## component is averaged over a COMMON reference distribution (all events of both types, each at its own chosen-target Z, D8-weighted):
##     fatality_g(x) = sum_e w_e [log r_g(x, Z_e) - log r_g(x0, Z_e)] / sum_e w_e.
##     total_g(x)    = selection_g(x) + fatality_g(x) = the weighted AVERAGE over reference events of the log P_F ratio (not the log of an
##                     averaged P_F); and the coal - bil difference of each component.
## ESTIMAND (joint audit B4): each type's fitted models applied to the SAME pooled reference events (a common-scale comparison: what each type's
## selection and fatality rules imply for one set of confrontations). Each curve is a one-dimensional slice through the surface with the other
## gap held at x0 (the dM x dR term makes slices at other values differ). The reference events, x0, the grids and the overlap flags are fixed
## from the original sample, so the bootstrap intervals are conditional on them.
## Step-B fatality curves on the probability scale use the same reference distribution: p_g(x) = sum_e w_e r_g(x, Z_e) / sum_e w_e.
## Own-type reference (joint audit M5): fatalityOwn_g, pfatalOwn_g and totalOwn_g average over type g's own events instead; their coal - bil
## differences mix populations by construction and are reported for comparison only.
##
## Common supported profiles: x0 = the pooled median (over all events' chosen targets) of dM and dR. Material curve: dM over 21 points from the
## pooled 5th to 95th percentile of chosen dM, dR = x0's. Relational curve: likewise for dR, dM = x0's. Overlap (joint audit B6): a box of
## +-0.25 SD on BOTH axes around each profile point; per type, thin_B = fewer than 20 of the type's CHOSEN targets in the box (step B is fit on
## chosen dyads, so a thin_B point is extrapolation for fatality and total) and thin_A = fewer than 20 of the type's ALTERNATIVES in the box.
## Intervals: pointwise 95% percentile, matched cluster bootstrap resampling DISPUTES by default (--boot_unit=dispute; or initiator), seed 123;
## reported only when B >= 200 and >= 90% of draws succeed for BOTH steps; replicate IDs, failure reasons and both steps' coefficient draws
## are saved. No labels (selection / failure / both) are produced: the curves and intervals are the evidence (GPT review).
## Identification (GPT review): step A is target choice conditional on a militarized confrontation and its type; step B is the probability that
## a confrontation turns jointly fatal given it occurred (escalation severity, not bargaining failure in general).
##
## Usage:
##   Rscript v2/R/12_two_step_joint.R --sample=mid --lag=3 --model=M1c --comparison=directed --events=dymid --variant=all --roles=bcj --joiner_set=all --types=bc
##   --outcome=fmid|fmidany|fmidub (default fmid), --reach=capdist, --boot=1000, --boot_unit=dispute|initiator,
##   --sim_check=TRUE --sim_seed=N (simulation check; outputs carry _SIMCHECK_seedN and never overwrite real outputs).
##   Unknown, retired or duplicated options are refused (joint audit B2).
suppressPackageStartupMessages({ library(dplyr); library(survival) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
R_DIR <- dirname(normalizePath(this_file))
source(file.path(R_DIR, "spec.R")); source(file.path(R_DIR, "helpers_data.R")); source(file.path(R_DIR, "helpers_model.R"))
source(file.path(R_DIR, "assemble_cell.R"))
Sys.setenv(V2_ROOT = normalizePath(file.path(R_DIR, "..")))

args <- commandArgs(trailingOnly = TRUE)
refuse_archived_options(args, "12_two_step_joint.R")
SCRIPT_FLAGS <- c("sample", "lag", "model", "boot", "types", "outcome", "boot_unit", "sim_check", "sim_seed")
KNOWN_FLAGS <- c(names(SPEC_CHOICES), "ame_K", "tag", SCRIPT_FLAGS)
bad <- args[!grepl("^--[A-Za-z_0-9]+=.+$", args)]
if (length(bad)) stop(sprintf("arguments must be --name=value; refused: %s", paste(bad, collapse = " ")))
flag_names <- sub("^--([A-Za-z_0-9]+)=.*$", "\\1", args)
if (any(!flag_names %in% KNOWN_FLAGS)) stop(sprintf("unknown option(s) for 12_two_step_joint.R: %s", paste(unique(flag_names[!flag_names %in% KNOWN_FLAGS]), collapse = ", ")))
if (anyDuplicated(flag_names)) stop(sprintf("option(s) given more than once: %s", paste(unique(flag_names[duplicated(flag_names)]), collapse = ", ")))
SAMPLE <- get_arg(args, "sample", "mid"); LAG <- as.integer(get_arg(args, "lag", 3)); MODEL <- get_arg(args, "model", "M1c")
B <- as.integer(get_arg(args, "boot", 1000)); TYPES <- get_arg(args, "types", "bc"); OUTCOME <- get_arg(args, "outcome", "fmid")
BOOT_UNIT <- get_arg(args, "boot_unit", "dispute"); SIM <- get_arg(args, "sim_check", FALSE); SIM_SEED <- as.integer(get_arg(args, "sim_seed", 123))
stopifnot(SAMPLE == "mid", LAG == 3L, MODEL == "M1c", B >= 0, TYPES %in% c("bc", "bc_nocoorig", "bcj"), OUTCOME %in% c("fmid", "fmidany", "fmidub"),
          BOOT_UNIT %in% c("dispute", "initiator"), is.logical(SIM), !is.na(SIM_SEED))
spec <- parse_spec(args)
stopifnot(spec_is_directed(spec), spec$roles == "bcj", spec_variant(spec) == "all", (if (is.null(spec$joiner_set)) "all" else spec$joiner_set) == "all",
          spec$controls == "minimal", spec$context == "all", (if (is.null(spec$mu_shift)) "none" else spec$mu_shift) == "none",
          spec_events(spec) == "dymid", spec$ww2 == "include")
cell <- paste0("joint_", cell_name(SAMPLE, LAG, MODEL, spec), "_", TYPES, "_", OUTCOME, if (SIM) paste0("_SIMCHECK_seed", SIM_SEED) else "")   # joint audit B1
log_open(paste0("12_", cell))
logf("cell %s; spec %s; types %s; outcome in_%s; matched bootstrap B = %d over %ss, seed 123%s", cell, spec_string(spec), TYPES, OUTCOME, B, BOOT_UNIT,
     if (SIM) sprintf("; SIMULATION CHECK seed %d", SIM_SEED) else "")

## ---- 1. assemble (identical to 10 and 11) ----------------------------------------------------------------------------
data_dir <- spec_data_dir(spec, v2_root())
manifest <- read.csv(file.path(data_dir, "manifest.csv"), stringsAsFactors = FALSE)
A <- assemble_cell_spec(SAMPLE, LAG, spec, data_dir = data_dir, manifest = manifest, model = MODEL)
sy <- A$sy; alts <- A$alts; ev <- A$ev; Z_cols <- A$Z_cols
dat <- build_stan_data(sy, alts, ev, A$W_cols, A$X_cols, Z_cols, joint = 1L, sampled = A$sampled, n_full = A$n_full,
                       x_ref = A$x_ref, allow_orphan_initiations = TRUE)
gt <- A$group_table
stopifnot(A$sampled == 0L, identical(gt$grp, 0:2), identical(gt$group, c("bilateral", "coalition_originator", "joiner")),
          length(ev$grp) == dat$E, all(ev$grp %in% 0:2), identical(as.integer(dat$ev_sy), as.integer(ev$sy_row)))

## ---- 2. types (identical rule to 10 and 11) ----------------------------------------------------------------------------
if (TYPES == "bcj") { ROLES <- c("bil", "coal", "join"); t_ev <- ev$grp; keep <- rep(TRUE, dat$E)
} else if (TYPES == "bc") { ROLES <- c("bil", "coal"); t_ev <- as.integer(ev$grp > 0L); keep <- rep(TRUE, dat$E)
} else { ROLES <- c("bil", "coal"); t_ev <- as.integer(ev$grp > 0L); keep <- ev$grp != 1L }
G <- length(ROLES)

## ---- 3. outcome (identical to 11) ------------------------------------------------------------------------------------------
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

## ---- 4. frames: L (step A, one row per alternative) and D (step B, one row per event), the same kept events --------------------
GAP <- c("dM", "dM2", "dR", "dR2", "dMR")
add_terms <- function(d) { d$dM2 <- d$dM^2; d$dR2 <- d$dR^2; d$dMR <- d$dM * d$dR
  for (g in 0:(G - 1L)) for (v in GAP) d[[paste0(v, "_", ROLES[g + 1])]] <- d[[v]] * (d$type == g); d }
GAP_G <- as.vector(outer(GAP, ROLES, paste, sep = "_"))
init_cc <- A$state_codes[sy$state_idx]
L <- bind_rows(lapply(which(keep), function(e) {
  r <- dat$ev_start[e]:dat$ev_end[e]; r <- r[dat$in_s2[r] == 1L]
  data.frame(ev_id = e, row = r, chosen = as.integer(r == dat$ev_row[e]), type = t_ev[e], w = dat$ev_w[e], init = init_cc[dat$ev_sy[e]])
}))
stopifnot(all(tapply(L$chosen, L$ev_id, sum) == 1L), all(tapply(L$row, L$ev_id, length) >= 2L), length(unique(L$ev_id)) == sum(keep))
L$one <- 1; L$dM <- dat$x[L$row, 1]; L$dR <- dat$x[L$row, 2]
Zm <- dat$Z[L$row, , drop = FALSE]; colnames(Zm) <- Z_cols; L <- cbind(L, Zm); L <- add_terms(L)
stopifnot(!anyNA(j$dispnum))
D <- data.frame(ev_id = which(keep), y = y_all[keep], type = t_ev[keep], w = dat$ev_w[keep], init = keys$initiator[keep], disp = j$dispnum[keep])
D$dM <- dat$x[dat$ev_row[D$ev_id], 1]; D$dR <- dat$x[dat$ev_row[D$ev_id], 2]
ZD <- dat$Z[dat$ev_row[D$ev_id], , drop = FALSE]; colnames(ZD) <- Z_cols; D <- cbind(D, ZD); D <- add_terms(D)
for (g in 0:(G - 1L)) D[[paste0("int_", ROLES[g + 1])]] <- as.integer(D$type == g)
INT <- paste0("int_", ROLES[-1])
stopifnot(!anyNA(L), !anyNA(D), identical(sort(unique(L$ev_id)), sort(D$ev_id)))
chA <- L[L$chosen == 1L, ]; chA <- chA[match(D$ev_id, chA$ev_id), ]
for (v in c("dM", "dR", "init", "type", "w", Z_cols)) stopifnot(isTRUE(all.equal(chA[[v]], D[[v]], check.attributes = FALSE)))   # B5: step B = step A's chosen rows
cnt <- vapply(0:(G - 1L), function(g) c(n = sum(D$type == g), fatal = sum(D$y[D$type == g])), numeric(2))
stopifnot(all(cnt[2, ] > 0), all(cnt[2, ] < cnt[1, ]))
logf("events %d (dropped %d); by type (events, fatal): %s; long frame %d rows", nrow(D), sum(!keep),
     paste(sprintf("%s %d / %d", ROLES, cnt[1, ], cnt[2, ]), collapse = "; "), nrow(L))

## the REAL chosen targets fix the profiles, the grids and the overlap counts (joint audit B7: identical across simulation seeds and in real runs)
D_real <- D[, c("ev_id", "type", "dM", "dR")]
## ---- 4b. optional simulation check (--sim_check=TRUE): replace the observed choices and fatal flags with draws from KNOWN coefficients on the
## real design (real choice sets, gaps, Z, weights, types); the decomposition must then recover decomp(truth). Types bc and bc_nocoorig only.
## Never used for results; its outputs carry _SIMCHECK_seedN. Z effects are nonzero in both steps (B7), so the averaging over Z is tested.
if (SIM) {
  stopifnot(TYPES %in% c("bc", "bc_nocoorig"), length(Z_cols) >= 3)
  set.seed(SIM_SEED)
  TRUE_A <- c(dM_bil = -0.5, dM2_bil = -0.2, dR_bil = 0.3, dR2_bil = 0, dMR_bil = 0.1, dM_coal = -1.2, dM2_coal = -0.1, dR_coal = 0.2, dR2_coal = -0.1, dMR_coal = -0.2,
              setNames(c(0.3, -0.2), Z_cols[1:2]))
  TRUE_B <- c(`(Intercept)` = -2, int_coal = 0.5, dM_bil = 0, dM2_bil = 0, dR_bil = 0.4, dR2_bil = 0, dMR_bil = 0, dM_coal = -0.6, dM2_coal = 0.1, dR_coal = 0, dR2_coal = 0, dMR_coal = 0.2,
              setNames(c(0.25, -0.3), Z_cols[c(1, 3)]))
  uA <- as.vector(as.matrix(L[, names(TRUE_A)]) %*% TRUE_A)
  L$chosen <- 0L
  for (idx in split(seq_len(nrow(L)), L$ev_id)) { p <- exp(uA[idx] - max(uA[idx])); L$chosen[idx[sample.int(length(idx), 1L, prob = p)]] <- 1L }
  ch <- L[L$chosen == 1L, ]; ch <- ch[match(D$ev_id, ch$ev_id), ]
  for (v in c("dM", "dR", Z_cols)) D[[v]] <- ch[[v]]
  D <- add_terms(D)
  etaB <- TRUE_B[["(Intercept)"]] + TRUE_B[["int_coal"]] * D$int_coal + as.vector(as.matrix(D[, names(TRUE_B)[-(1:2)]]) %*% TRUE_B[-(1:2)])
  D$y <- rbinom(nrow(D), 1L, plogis(etaB))
  logf("SIMULATION CHECK: choices and fatal flags replaced by draws from known coefficients (%d fatal of %d)", sum(D$y), nrow(D))
}

## ---- 5. fits -------------------------------------------------------------------------------------------------------------
FORM_A <- as.formula(paste("Surv(one, chosen) ~", paste(c(GAP_G, Z_cols), collapse = " + "), "+ strata(ev_id)"))
FORM_B <- as.formula(paste("y ~", paste(c(INT, GAP_G, Z_cols), collapse = " + ")))
fitA <- function(d) coxph(FORM_A, data = d, weights = w, method = "breslow")
fitB <- function(d) glm(FORM_B, data = d, weights = w, family = quasibinomial())
mA <- fitA(L); mB <- fitB(D)
bA <- coef(mA); bB <- coef(mB)
stopifnot(all(is.finite(bA)), all(is.finite(bB)), mB$converged, all(c(GAP_G, Z_cols) %in% names(bA)), all(c("(Intercept)", INT, GAP_G, Z_cols) %in% names(bB)))

## ---- 6. profiles, overlap and the decomposition ---------------------------------------------------------------------------
x0 <- c(dM = median(D_real$dM), dR = median(D_real$dR))
NG <- 21L; H_SD <- 0.25; MIN_OVL <- 20L
grid <- list(dM = seq(quantile(D_real$dM, 0.05), quantile(D_real$dM, 0.95), length.out = NG), dR = seq(quantile(D_real$dR, 0.05), quantile(D_real$dR, 0.95), length.out = NG))
profiles <- rbind(data.frame(axis = "material", dM = grid$dM, dR = x0[["dR"]]), data.frame(axis = "relational", dM = x0[["dM"]], dR = grid$dR))
profiles$x <- ifelse(profiles$axis == "material", profiles$dM, profiles$dR)
hM <- H_SD * sd(D_real$dM); hR <- H_SD * sd(D_real$dR)                     # B6: two-axis box, +-0.25 SD on each axis
alt_rows <- L[, c("ev_id", "type", "dM", "dR")]                            # every alternative (the real design, also under --sim_check), for thin_A
for (g in 0:(G - 1L)) {
  sB <- D_real$type == g; sA <- alt_rows$type == g
  profiles[[paste0("nB_", ROLES[g + 1])]] <- vapply(seq_len(nrow(profiles)), function(i)
    sum(abs(D_real$dM[sB] - profiles$dM[i]) <= hM & abs(D_real$dR[sB] - profiles$dR[i]) <= hR), 0L)
  profiles[[paste0("nA_", ROLES[g + 1])]] <- vapply(seq_len(nrow(profiles)), function(i)
    sum(abs(alt_rows$dM[sA] - profiles$dM[i]) <= hM & abs(alt_rows$dR[sA] - profiles$dR[i]) <= hR), 0L)
}
profiles$thin_B <- apply(profiles[, paste0("nB_", ROLES), drop = FALSE] < MIN_OVL, 1, any)
profiles$thin_A <- apply(profiles[, paste0("nA_", ROLES), drop = FALSE] < MIN_OVL, 1, any)
logf("reference x0 = (dM %.3f, dR %.3f); material grid %.2f to %.2f, relational grid %.2f to %.2f; overlap box +-%.3f (dM) x +-%.3f (dR); thin for step B (any type < %d chosen): %d of %d; thin for step A (alternatives): %d",
     x0[["dM"]], x0[["dR"]], min(grid$dM), max(grid$dM), min(grid$dR), max(grid$dR), hM, hR, MIN_OVL, sum(profiles$thin_B), nrow(profiles), sum(profiles$thin_A))
gterms <- function(dM, dR) cbind(dM, dM^2, dR, dR^2, dM * dR)            # in GAP order
ZD_mat <- as.matrix(D[, Z_cols, drop = FALSE]); wD <- D$w / sum(D$w)
decomp <- function(bA, bB) {
  out <- list()
  T0 <- gterms(x0[["dM"]], x0[["dR"]]); TP <- gterms(profiles$dM, profiles$dR)
  zlin <- as.vector(ZD_mat %*% bB[Z_cols])                                # each reference event's Z contribution in step B
  for (gi in seq_len(G)) { g <- ROLES[gi]
    ba <- bA[paste0(GAP, "_", g)]; bb <- bB[paste0(GAP, "_", g)]
    sel <- as.vector(TP %*% ba) - as.vector(T0 %*% ba)
    c0 <- bB[["(Intercept)"]] + (if (gi > 1) bB[[paste0("int_", g)]] else 0)
    eta0 <- c0 + zlin + as.vector(T0 %*% bb)                                # per reference event, at x0
    lr0 <- plogis(eta0, log.p = TRUE)
    wO <- ifelse(D$type == gi - 1L, D$w, 0); wO <- wO / sum(wO)             # M5: the type's own events
    fat <- pB <- fatO <- pBO <- numeric(nrow(TP))
    for (i in seq_len(nrow(TP))) { eta <- c0 + zlin + sum(TP[i, ] * bb); dl <- plogis(eta, log.p = TRUE) - lr0; pe <- plogis(eta)
      fat[i] <- sum(wD * dl); pB[i] <- sum(wD * pe); fatO[i] <- sum(wO * dl); pBO[i] <- sum(wO * pe) }
    out[[paste0("selection_", g)]] <- sel; out[[paste0("fatality_", g)]] <- fat; out[[paste0("total_", g)]] <- sel + fat; out[[paste0("pfatal_", g)]] <- pB
    out[[paste0("fatalityOwn_", g)]] <- fatO; out[[paste0("totalOwn_", g)]] <- sel + fatO; out[[paste0("pfatalOwn_", g)]] <- pBO
  }
  if (G > 1) for (o in ROLES[-1]) for (cmp in c("selection", "fatality", "total", "pfatal", "fatalityOwn", "totalOwn", "pfatalOwn"))
    out[[paste0(cmp, "_", o, "-", ROLES[1])]] <- out[[paste0(cmp, "_", o)]] - out[[paste0(cmp, "_", ROLES[1])]]
  out
}
pt <- decomp(bA, bB)
if (SIM) {
  tA <- setNames(rep(0, length(bA)), names(bA)); tA[names(TRUE_A)] <- TRUE_A
  tB <- setNames(rep(0, length(bB)), names(bB)); tB[names(TRUE_B)] <- TRUE_B
  truth <- decomp(tA, tB)
  err <- vapply(names(pt), function(k) max(abs(pt[[k]] - truth[[k]])), 0)
  logf("SIMULATION CHECK: max |estimate - truth| over the grid, by curve:\n%s", paste(sprintf("  %-24s %.3f  (truth range %.2f to %.2f)", names(err), err,
       vapply(truth, min, 0), vapply(truth, max, 0)), collapse = "\n"))
  saveRDS(list(pt = pt, truth = truth, err = err, profiles = profiles), file.path(v2_root(), "output", paste0(cell, ".rds")))
  if (B == 0) quit(save = "no")
}

## ---- 7. matched cluster bootstrap (unit: dispute by default, or initiating state) ------------------------------------------
quiet_fit <- function(f, d) { warned <- FALSE
  res <- tryCatch(withCallingHandlers(f(d), warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") }), error = function(e) NULL)
  list(fit = res, warned = warned) }
set.seed(123)
unit_of <- if (BOOT_UNIT == "dispute") D$disp else D$init
units <- unique(unit_of); ev_by_unit <- split(D$ev_id, unit_of)
logf("bootstrap unit: %s (%d units)", BOOT_UNIT, length(units))
Lsplit <- split(L, L$ev_id); Drow <- setNames(seq_len(nrow(D)), D$ev_id)
bs <- vector("list", length(pt)); names(bs) <- names(pt)
coefA <- NULL; coefB <- NULL; fail_log <- data.frame(draw = integer(0), reason = character(0))
for (bb in seq_len(B)) {
  evs <- unlist(ev_by_unit[as.character(sample(units, length(units), replace = TRUE))], use.names = FALSE)
  Lb <- bind_rows(lapply(seq_along(evs), function(i) { x <- Lsplit[[as.character(evs[i])]]; x$ev_id <- i; x }))
  Db <- D[Drow[as.character(evs)], ]; Db$ev_id <- seq_along(evs)
  tt <- table(factor(Db$type, levels = 0:(G - 1L))); ft <- tapply(Db$y, factor(Db$type, levels = 0:(G - 1L)), sum)
  if (any(tt == 0L) || any(ft == 0L, na.rm = TRUE) || any(ft == tt, na.rm = TRUE)) {
    fail_log <- rbind(fail_log, data.frame(draw = bb, reason = "a type lacks events or lacks fatal / non-fatal events")); next }
  fa <- quiet_fit(fitA, Lb); fb <- quiet_fit(fitB, Db)
  bad <- c(if (is.null(fa$fit)) "A error", if (fa$warned) "A warning", if (!is.null(fa$fit) && any(!is.finite(coef(fa$fit)))) "A non-finite",
           if (is.null(fb$fit)) "B error", if (fb$warned) "B warning (e.g. separation)", if (!is.null(fb$fit) && (!fb$fit$converged || any(!is.finite(coef(fb$fit))))) "B non-convergence")
  if (length(bad)) { fail_log <- rbind(fail_log, data.frame(draw = bb, reason = paste(bad, collapse = "; "))); next }
  r <- decomp(coef(fa$fit), coef(fb$fit))
  for (k in names(r)) bs[[k]] <- rbind(bs[[k]], r[[k]])
  coefA <- rbind(coefA, c(draw = bb, coef(fa$fit))); coefB <- rbind(coefB, c(draw = bb, coef(fb$fit)))
  if (bb %% 100 == 0) logf("bootstrap %d / %d", bb, B)
}
MIN_OK <- 0.90; MIN_B <- 200L                                      # joint audit B3
n_ok <- if (is.null(coefA)) 0L else nrow(coefA)
boot_ok <- B >= MIN_B && n_ok / B >= MIN_OK
logf("matched bootstrap: %d draws, %d successful (both steps), %d dropped; intervals %s", B, n_ok, nrow(fail_log),
     if (boot_ok) "reported" else if (B < MIN_B) sprintf("NOT reported (B < %d)", MIN_B) else sprintf("NOT reported (fewer than %.0f%% succeeded)", 100 * MIN_OK))

## ---- 8. tables and figure ---------------------------------------------------------------------------------------------------
curves <- NULL
for (k in names(pt)) {
  sp <- strsplit(k, "_", fixed = TRUE)[[1]]; comp <- sp[1]; who <- paste(sp[-1], collapse = "_")
  lo <- hi <- rep(NA_real_, nrow(profiles))
  if (boot_ok) { lo <- apply(bs[[k]], 2, quantile, 0.025); hi <- apply(bs[[k]], 2, quantile, 0.975) }
  curves <- rbind(curves, data.frame(component = comp, type = who, axis = profiles$axis, x = profiles$x, dM = profiles$dM, dR = profiles$dR,
                                     estimate = pt[[k]], lo = lo, hi = hi, thin_A = profiles$thin_A, thin_B = profiles$thin_B,
                                     profiles[, c(paste0("nA_", ROLES), paste0("nB_", ROLES)), drop = FALSE]))
}
od <- file.path(v2_root(), "output"); dir.create(od, showWarnings = FALSE)
write.csv(curves, file.path(od, paste0(cell, "_curves.csv")), row.names = FALSE)
write.csv(data.frame(term = names(bA), coef = bA), file.path(od, paste0(cell, "_coefA.csv")), row.names = FALSE)
write.csv(data.frame(term = names(bB), coef = bB), file.path(od, paste0(cell, "_coefB.csv")), row.names = FALSE)
saveRDS(list(pt = pt, bs = if (boot_ok) bs else NULL, coefA = coefA, coefB = coefB, fail_log = fail_log, n_ok = n_ok, B = B, x0 = x0, profiles = profiles,
             spec = spec, types = TYPES, outcome = OUTCOME, roles = ROLES, counts = cnt, E = nrow(D), boot_unit = BOOT_UNIT), file.path(od, paste0(cell, "_boot.rds")))
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  cd <- curves[curves$component %in% c("selection", "fatality", "total") & curves$type %in% ROLES, ]
  cd$thin <- ifelse(cd$component == "selection", cd$thin_A, cd$thin_B | cd$thin_A)          # selection needs alternatives; fatality and total need chosen
  cd$component <- factor(cd$component, levels = c("selection", "fatality", "total"),
                         labels = c("Selection\nlog P_A ratio", "Fatality\nlog risk ratio", "Total\nlog P_F ratio"))
  cd$axis <- factor(cd$axis, levels = c("material", "relational"), labels = c("Material gap (initiator - target)", "Relational gap (initiator - target)"))
  p <- ggplot(cd, aes(x = x, y = estimate, colour = type, fill = type)) + geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    { if (boot_ok) geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) } + geom_line(linewidth = 0.7) +
    geom_point(data = cd[cd$thin, ], shape = 4, size = 1.5, show.legend = FALSE) +
    facet_grid(component ~ axis, scales = "free") + labs(x = NULL, y = "Relative to the reference profile", colour = "Type", fill = "Type",
         caption = "x = thin support: a type has < 20 targets (alternatives for selection, chosen for fatality and total) within 0.25 SD on both gaps") +
    theme_bw(base_size = 10) + theme(legend.position = "bottom", strip.text.y = element_text(angle = 0))
  ggsave(file.path(od, paste0(cell, "_decomposition.png")), p, width = 8, height = 7, dpi = 200)
}
logf("done: %s_{curves,coefA,coefB}.csv, _boot.rds, _decomposition.png", cell)
