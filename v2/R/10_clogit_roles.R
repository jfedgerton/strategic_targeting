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
## Intervals: percentile, seed 123. PRIMARY (GPT recheck 2026-11-16, PI-approved): pigeonhole two-way bootstrap, resampling initiators and
##   disputes independently and multiplying each event's D8 weight by the product of the two multiplicities (--boot_unit=pigeonhole, default).
##   --boot_unit=dispute or initiator are one-way sensitivities. The unit is part of the output name (_bpig, _bdisp, _binit). Intervals are reported only when B >= 200
##   and >= 90% of the draws succeed (joint audit B3); a draw whose fit errors or warns is dropped for that model and logged with its reason.
## --loo=TRUE adds leave-one-dispute-out and leave-one-initiator-out refits of M1 (joint audit M3), written to <cell>_loo.csv.
## Relational-measure diagnostics (PI request 2026-11-16, measurement check of log IGO PageRank):
##   --rel_scale=year rescales the relational gap by the WITHIN-YEAR spread of log PageRank instead of the pooled spread (the measure compresses
##   over time: within-year SD 0.85 in 1880, 0.25 in 2000), i.e. dR_year = dR_pooled * GAP_SCALE_R / gap_scale(R_pr in the event's year);
##   --era=pre1946|from1946 keeps only events before 1946 or from 1946 on. Both are diagnostics, not the main specification.
## Per-type axis tests (3 df: linear, squared, cross = 0) are also reported: "does this axis matter for this type at all".
##
## Usage (same flags as the grid cell; --roles=bcj is required because the types are recoded from it):
##   Rscript v2/R/10_clogit_roles.R --sample=fmid --lag=3 --model=M1c --comparison=directed --entry=later --variant=all --roles=bcj --joiner_set=all --types=bc
##   --events=dymid (instead of --entry=later) for jointly fatal, --reach=capdist for the reach terms, --boot=1000, --boot_unit=pigeonhole|dispute|initiator,
##   --loo=TRUE, --coorig_months=1|2|3 (with --coorig_bound=definite|possible) or --entry_mod=TRUE (entry-timing robustness, section 1b). Unknown, retired or duplicated options are refused (joint audit B2).
suppressPackageStartupMessages({ library(dplyr); library(survival) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
R_DIR <- dirname(normalizePath(this_file))
source(file.path(R_DIR, "spec.R")); source(file.path(R_DIR, "helpers_data.R")); source(file.path(R_DIR, "helpers_model.R"))
source(file.path(R_DIR, "assemble_cell.R"))
Sys.setenv(V2_ROOT = normalizePath(file.path(R_DIR, "..")))   # audit F5 (A1 P2-3): read and write this checkout only, as 02_fit_model.R does

args <- commandArgs(trailingOnly = TRUE)
## joint audit B2: refuse retired options, unknown flags, bare flags and duplicated flags (get_arg would silently take the first)
refuse_archived_options(args, "10_clogit_roles.R")
SCRIPT_FLAGS <- c("sample", "lag", "model", "boot", "types", "boot_unit", "loo", "rel_scale", "era", "coorig_months", "coorig_bound", "entry_mod")
KNOWN_FLAGS <- c(names(SPEC_CHOICES), "ame_K", "tag", SCRIPT_FLAGS)
bad <- args[!grepl("^--[A-Za-z_0-9]+=.+$", args)]
if (length(bad)) stop(sprintf("arguments must be --name=value; refused: %s", paste(bad, collapse = " ")))
flag_names <- sub("^--([A-Za-z_0-9]+)=.*$", "\\1", args)
if (any(!flag_names %in% KNOWN_FLAGS)) stop(sprintf("unknown option(s) for 10_clogit_roles.R: %s", paste(unique(flag_names[!flag_names %in% KNOWN_FLAGS]), collapse = ", ")))
if (anyDuplicated(flag_names)) stop(sprintf("option(s) given more than once: %s", paste(unique(flag_names[duplicated(flag_names)]), collapse = ", ")))
SAMPLE <- get_arg(args, "sample", "fmid"); LAG <- as.integer(get_arg(args, "lag", 3)); MODEL <- get_arg(args, "model", "M1c")
B <- as.integer(get_arg(args, "boot", 1000))
TYPES <- get_arg(args, "types", "bc"); BOOT_UNIT <- get_arg(args, "boot_unit", "pigeonhole"); LOO <- get_arg(args, "loo", FALSE)
REL_SCALE <- get_arg(args, "rel_scale", "pooled"); ERA <- get_arg(args, "era", "all")
COORIG_W <- as.integer(get_arg(args, "coorig_months", 0)); ENTRY_MOD <- get_arg(args, "entry_mod", FALSE)   # PI request 2026-11-17 (robustness; see section 1b)
stopifnot(COORIG_W %in% c(0L, 1L, 2L, 3L), is.logical(ENTRY_MOD), !(ENTRY_MOD && COORIG_W > 0L))
stopifnot(REL_SCALE %in% c("pooled", "year"), ERA %in% c("all", "pre1946", "from1946"))
stopifnot(SAMPLE %in% ALL_SAMPLES, LAG == 3L, MODEL == "M1c", B >= 0, TYPES %in% c("bc", "bc_nocoorig", "bcj"), BOOT_UNIT %in% c("pigeonhole", "dispute", "initiator"), is.logical(LOO))
spec <- parse_spec(args)
## audit F5 (A1 P2-1): the cells M1c is defined on, and nothing else (02_fit_model.R refuses the same)
stopifnot(spec_is_directed(spec), spec$roles == "bcj", spec_variant(spec) == "all", (if (is.null(spec$joiner_set)) "all" else spec$joiner_set) == "all",
          spec$controls == "minimal", spec$context == "all", (if (is.null(spec$mu_shift)) "none" else spec$mu_shift) == "none")
cell <- paste0("clogit_", cell_name(SAMPLE, LAG, MODEL, spec), "_", TYPES, if (REL_SCALE == "year") "_relyear" else "", if (ERA != "all") paste0("_", ERA) else "",
              if (COORIG_W > 0L) paste0("_com", COORIG_W, if (get_arg(args, "coorig_bound", "definite") == "possible") "p" else "d") else "", if (ENTRY_MOD) "_entrymod" else "",
              "_b", c(pigeonhole = "pig", dispute = "disp", initiator = "init")[[BOOT_UNIT]])   # GPT recheck P2: bootstrap unit is part of the run identity
log_open(paste0("10_", cell))
logf("cell %s; spec %s; types %s; bootstrap B = %d, unit %s, seed 123; leave-one-out %s; relational scale %s; era %s", cell, spec_string(spec), TYPES, B, BOOT_UNIT, LOO, REL_SCALE, ERA)

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

## ---- 1b. entry timing (PI request 2026-11-17; robustness for early vs late coordination; rebuilt after the GPT and independent reviews) -----
## The default coding makes every state that enters a dispute after its onset a joiner, however soon after. Two robustness options:
##   --coorig_months=1|2|3 : an "early same-side participant" sensitivity (not a claim that late entrants were originators). A joiner whose entry
##       is within W month-index steps of the dispute onset is recoded coalition originator, and an alone originator becomes a coalition originator
##       when ANOTHER participant on its own side at entry enters within W. W = 0 is the default coding. Month-index windows: W = 1 can span almost
##       two months of elapsed time.
##   --coorig_bound=definite|possible : dates known only to the year make the classification uncertain. definite (default) counts an entry as
##       early only when its LATEST possible month is within W of the EARLIEST possible onset month; possible counts it when its earliest possible
##       month is within W of the latest possible onset. Every window is run under both; the events whose role differs are written out.
##   --entry_mod=TRUE (types bcj only): the joiners' five gap terms are interacted with m = log(1 + entry lag in months, midpoint of its bounds),
##       centred at the joiners' weighted median, so the joiner coefficients and every joiner-based quantity are read at the median joiner. The lag is
##       descriptive heterogeneity, not a pre-treatment moderator (entry timing and the target are often chosen together).
## Source: participants_timing.rds of the data build, one row per (dispute, participant, spell) of the event participation table BEFORE any key
## collapse: COW originator flag, side in the spell's first year (side_conflict, side_switch), participant entry month bounds (GML participant dates
## at day or month resolution, else the whole first year) and the dispute onset month bounds (COW start). A fitted event is matched to its weight-
## group dispute (dispnum of its non-escalation row; the role comes from the same dispute) and to the initiator's spell covering the event year.
ev_rds <- readRDS(file.path(data_dir, "events.rds"))
lag_ev <- rep(NA_real_, dat$E); lag_lo_ev <- lag_hi_ev <- rep(NA_real_, dat$E); role_alt <- NULL
COORIG_BOUND <- get_arg(args, "coorig_bound", "definite"); stopifnot(COORIG_BOUND %in% c("definite", "possible"))
if (COORIG_W > 0L || ENTRY_MOD) {
  stopifnot(TYPES %in% c("bc", "bcj"))
  pt_path <- file.path(data_dir, "participants_timing.rds")
  if (!file.exists(pt_path)) stop("entry timing: ", pt_path, " is missing; rebuild the data with the 2026-11-17 build (01_build_data.R)")
  P <- readRDS(pt_path)
  stopifnot(all(c("dispnum", "ccode", "spell", "spell_first", "spell_last", "orig", "side_entry", "side_conflict", "entry_lo", "entry_hi", "onset_lo", "onset_hi") %in% names(P)))
  ne_t <- ev_rds[!ev_rds$escalation, c("initiator", "target", "year", "dispnum", "dispnums")]
  stopifnot(!anyDuplicated(ne_t[, c("initiator", "target", "year")]))
  k_all <- data.frame(e = seq_len(dat$E), initiator = as.integer(sy$ccode[ev$sy_row]), target = as.integer(alts$ccode2[ev$alt_row]), year = as.integer(sy$year[ev$sy_row]))
  k_all <- merge(k_all, ne_t, by = c("initiator", "target", "year"), all.x = TRUE, sort = FALSE); k_all <- k_all[order(k_all$e), ]
  if (anyNA(k_all$dispnum)) stop(sprintf("entry timing: %d stage-2 events have no non-escalation row in events.rds", sum(is.na(k_all$dispnum))))
  stopifnot(identical(k_all$e, seq_len(dat$E)))
  ## the initiator's spell in its weight-group dispute that covers the event year (the latest spell starting by then)
  sp_of <- function(d, cc, y) { r <- which(P$dispnum == d & P$ccode == cc & P$spell_first <= y); if (!length(r)) NA_integer_ else r[which.max(P$spell_first[r])] }
  prow <- mapply(sp_of, k_all$dispnum, k_all$initiator, k_all$year)
  if (anyNA(prow)) stop(sprintf("entry timing: %d fitted events have no participant spell for their initiator in their weight-group dispute", sum(is.na(prow))))
  ## the build's originator flag must agree with the role: joiners (bcj 2) are exactly the non-originators of the weight-group dispute
  if (any((P$orig[prow] == 0L) != (ev$grp == 2L))) stop(sprintf("entry timing: originator flag and role disagree for %d events", sum((P$orig[prow] == 0L) != (ev$grp == 2L))))
  lag_lo_ev <- pmax(P$entry_lo[prow] - P$onset_hi[prow], 0)                      # earliest possible lag
  lag_hi_ev <- pmax(P$entry_hi[prow] - P$onset_lo[prow], 0)                      # latest possible lag
  lag_ev <- (lag_lo_ev + lag_hi_ev) / 2                                           # midpoint (moderator)
  lag_ev[ev$grp != 2L] <- 0; lag_lo_ev[ev$grp != 2L] <- 0; lag_hi_ev[ev$grp != 2L] <- 0   # originators enter at onset by definition
  jl <- lag_ev[ev$grp == 2L]
  logf("entry timing: joiner entry lag (months, midpoint) quantiles 10/25/50/75/90%% = %s; joiners whose lag bounds differ %d of %d; definitely / possibly within 1, 2, 3 months: %s",
       paste(round(quantile(jl, c(.1, .25, .5, .75, .9)), 1), collapse = " / "), sum(lag_lo_ev[ev$grp == 2L] != lag_hi_ev[ev$grp == 2L]), length(jl),
       paste(sprintf("%d/%d", sapply(1:3, function(w) sum(lag_hi_ev[ev$grp == 2L] <= w)), sapply(1:3, function(w) sum(lag_lo_ev[ev$grp == 2L] <= w))), collapse = ", "))
  if (COORIG_W > 0L) {
    W <- COORIG_W
    ## the recoding rule under one bound convention; partners: other participants of the dispute on the same side at entry, side not conflicted
    recode <- function(bound) {
      early_lag <- function(lo, hi) if (bound == "definite") hi <= W else lo <= W
      Pl <- P; Pl$lag_lo <- pmax(Pl$entry_lo - Pl$onset_hi, 0); Pl$lag_hi <- pmax(Pl$entry_hi - Pl$onset_lo, 0)
      Pl$early_partner <- Pl$orig == 0L & !Pl$side_conflict & early_lag(Pl$lag_lo, Pl$lag_hi)
      g <- ev$grp
      early_j <- ev$grp == 2L & early_lag(lag_lo_ev, lag_hi_ev)
      my_side <- P$side_entry[prow]; my_conf <- P$side_conflict[prow]
      alone_up <- ev$grp == 0L & !my_conf & vapply(seq_len(dat$E), function(e) any(Pl$early_partner & Pl$dispnum == k_all$dispnum[e] & Pl$side_entry == my_side[e] & Pl$ccode != k_all$initiator[e]), TRUE)
      g[early_j | alone_up] <- 1L
      list(grp = g, early_j = early_j, alone_up = alone_up)
    }
    rc <- recode(COORIG_BOUND); rc_alt <- recode(setdiff(c("definite", "possible"), COORIG_BOUND))
    grp_old <- ev$grp; ev$grp <- rc$grp
    diff_e <- which(rc$grp != rc_alt$grp)
    role_alt <- data.frame(initiator = k_all$initiator[diff_e], target = k_all$target[diff_e], year = k_all$year[diff_e], dispnum = k_all$dispnum[diff_e],
                           grp_this = rc$grp[diff_e], grp_other_bound = rc_alt$grp[diff_e], w = dat$ev_w[diff_e])
    logf("coorig window %d months, %s bound: joiners recoded coalition originator %d; alone originators recoded %d; bcj now %d / %d / %d (was %d / %d / %d); events whose role differs under the %s bound %d (weight %.1f)",
         W, COORIG_BOUND, sum(rc$early_j), sum(rc$alone_up), sum(ev$grp == 0L), sum(ev$grp == 1L), sum(ev$grp == 2L), sum(grp_old == 0L), sum(grp_old == 1L), sum(grp_old == 2L),
         setdiff(c("definite", "possible"), COORIG_BOUND), length(diff_e), sum(dat$ev_w[diff_e]))
    g8 <- tapply(ev$grp, paste(k_all$initiator, k_all$dispnum, k_all$year), function(x) length(unique(x)))
    if (any(g8 > 1L)) stop(sprintf("coorig window: roles not constant within %d (initiator, dispute, year) groups", sum(g8 > 1L)))
    ## regression checks (independent review 2026-11-17)
    gp <- which(k_all$initiator == 255L & k_all$target == 290L & k_all$year == 1939L)
    if (length(gp)) { stopifnot(ev$grp[gp] == 0L); logf("check: Germany -> Poland 1939 stays bilateral") }
    for (chk in list(c(2L, 645L, 1997L), c(2L, 731L, 1950L))) { r <- which(k_all$initiator == chk[1] & k_all$target == chk[2] & k_all$year == chk[3])
      if (length(r)) logf("check: %d -> %d %d coded %d (was %d)", chk[1], chk[2], chk[3], ev$grp[r[1]], grp_old[r[1]]) }
  }
}

## ---- 2. types: recode the bcj groups (0-based type index t per event; keep = events in stage 2) ------------------
if (TYPES == "bcj") {
  ROLES <- c("bil", "coal", "join"); t_ev <- ev$grp; keep <- rep(TRUE, dat$E)
} else if (TYPES == "bc") {
  ROLES <- c("bil", "coal"); t_ev <- as.integer(ev$grp > 0L); keep <- rep(TRUE, dat$E)
} else {
  ROLES <- c("bil", "coal"); t_ev <- as.integer(ev$grp > 0L); keep <- ev$grp != 1L
}
G <- length(ROLES)
ev_year <- as.integer(sy$year[ev$sy_row])
if (ERA == "pre1946") keep <- keep & ev_year < 1946L
if (ERA == "from1946") keep <- keep & ev_year >= 1946L
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
if (REL_SCALE == "year") {
  ## the build divides every relational log-ratio by ONE pooled scale (gap_scale over all years); here it is replaced by that year's own scale
  posy <- readRDS(file.path(data_dir, "sy_positions.rds"))
  stopifnot(all(c("year", "R_pr") %in% names(posy)))
  GSR <- gap_scale(posy$R_pr, posy$year)
  ys <- vapply(split(posy$R_pr, posy$year), function(v) gap_scale(v, rep(1L, length(v))), 0)
  yr_row <- as.character(sy$year[dat$ev_sy[L$ev_id]])
  stopifnot(all(yr_row %in% names(ys)), all(is.finite(ys)), all(ys > 0))
  L$dR <- L$dR * GSR / ys[yr_row]
  logf("relational gap rescaled within year: pooled scale %.3f; year scales %.3f to %.3f (median %.3f)", GSR, min(ys), max(ys), median(ys))
}
Zm <- dat$Z[L$row, , drop = FALSE]; colnames(Zm) <- Z_cols; L <- cbind(L, Zm)
stopifnot(!anyNA(L[, c("dM", "dR", Z_cols)]))
L$dM2 <- L$dM^2; L$dR2 <- L$dR^2; L$dMR <- L$dM * L$dR
GAP <- c("dM", "dM2", "dR", "dR2", "dMR")
for (g in 0:(G - 1L)) for (v in GAP) L[[paste0(v, "_", ROLES[g + 1])]] <- L[[v]] * (L$type == g)
GAP_G <- as.vector(outer(GAP, ROLES, paste, sep = "_"))           # dM_bil, dM2_bil, ..., dMR_<last type>
MOD <- character(0)
if (ENTRY_MOD) {                                                  # section 1b: joiner gap terms x centred log(1 + entry lag)
  stopifnot(TYPES == "bcj")
  jw <- dat$ev_w[keep & t_ev == 2L]; jm <- log1p(lag_ev[keep & t_ev == 2L]); o <- order(jm)
  M_CTR <- jm[o][which(cumsum(jw[o]) >= sum(jw) / 2)[1]]          # weighted median over joiner events
  L$m_lag <- ifelse(L$type == 2L, log1p(lag_ev[L$ev_id]) - M_CTR, 0)
  stopifnot(!anyNA(L$m_lag))
  for (v in GAP) L[[paste0(v, "_jlag")]] <- L[[paste0(v, "_join")]] * L$m_lag
  MOD <- paste0(GAP, "_jlag"); GAP_G <- c(GAP_G, MOD)
  logf("entry moderator: centred at the weighted median joiner lag %.0f months (log1p %.3f); joiner m range %.2f to %.2f", expm1(M_CTR), M_CTR, min(L$m_lag[L$type == 2L]), max(L$m_lag[L$type == 2L]))
}
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
W_roles_med <- NULL
if (ENTRY_MOD) {   # GPT review 2026-11-17 P1: with the moderator the full nested null (M1 vs M0) also sets the five *_jlag terms to zero
  R_mod <- diag(length(nm1))[match(MOD, nm1), , drop = FALSE]; colnames(R_mod) <- nm1; stopifnot(!anyNA(match(MOD, nm1)))
  W_roles_med <- W_roles; W_roles <- wald3(rbind(R_roles, R_mod))
  logf("entry moderator: surfaces differ at the median joiner lag only (%d df) chi2 %.2f p %.4g; the full M1 vs M0 test below includes the moderator terms", W_roles_med[["df"]], W_roles_med[["chi2"]], W_roles_med[["p"]])
}
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
           ifelse(pw$term %in% c("dM", "dR"), "difference at dM = dR = 0 only (squared and cross terms can change its sign elsewhere)",
           ifelse(grepl("^joint:dM\\+dM2\\+dMR$|^joint:dR\\+dR2\\+dMR$", pw$term), "axis effect differs anywhere (linear, squared, cross); not mu or sigma", "")))
if (ENTRY_MOD) pw$note <- ifelse(grepl("join", pw$pair), paste0(pw$note, ifelse(nchar(pw$note), "; ", ""), "joiner side evaluated at the median joiner entry lag"), pw$note)

if (ENTRY_MOD) {                                                  # section 1b: do joiners' gap effects change with entry lag? (5 df, primary two-way)
  W_mod <- wald3(R_mod)
  logf("entry moderator: joiner gap terms x log(1 + lag) jointly zero chi2 %.2f df %d p %.4g (initiator %.4g / target-year %.4g)",
       W_mod[["chi2"]], W_mod[["df"]], W_mod[["p"]], W_mod[["p_initiator"]], W_mod[["p_tgtyr"]])
}
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

## ---- 6. cluster bootstrap (pigeonhole two-way by default; dispute or initiator one-way sensitivities) --------------------------------------
## audit F5 (A2/A3 P2): a fit that errors OR warns (separation, non-convergence) is dropped for that model only; failures counted per model
fit_or_null <- function(d, terms) {
  warned <- FALSE
  res <- tryCatch(withCallingHandlers(coef(fit_cl(d, terms, robust = FALSE)), warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") }),
                  error = function(e) NULL)
  if (warned || is.null(res) || any(!is.finite(res))) NULL else res
}
## GPT recheck 2026-11-16 P2 (PI-approved): PRIMARY = pigeonhole two-way bootstrap (Owen and Eckles 2012; Davezies, D'Haultfoeuille and Guyonvarch
## 2021): initiators and disputes are resampled INDEPENDENTLY with replacement and each event's fixed D8 weight is multiplied by the product of its
## initiator's and its dispute's multiplicities (applied to every alternative in its choice set; zero-multiplier events drop out). Dispute-only
## and initiator-only resampling (multiplier = that unit's multiplicity) are labelled one-way sensitivities. A stratum's weight is constant
## across its rows, so a multiplier k is identical to duplicating the event k times.
set.seed(123)
ev_ui <- unique(L[, c("ev_id", "init", "disp")]); stopifnot(!anyDuplicated(ev_ui$ev_id))
u_init <- unique(ev_ui$init); u_disp <- unique(ev_ui$disp)
i_idx <- match(ev_ui$init, u_init); d_idx <- match(ev_ui$disp, u_disp)
draw_mult <- function() {
  mi <- if (BOOT_UNIT == "dispute") rep(1L, length(u_init)) else tabulate(sample.int(length(u_init), length(u_init), replace = TRUE), length(u_init))
  md <- if (BOOT_UNIT == "initiator") rep(1L, length(u_disp)) else tabulate(sample.int(length(u_disp), length(u_disp), replace = TRUE), length(u_disp))
  setNames(mi[i_idx] * md[d_idx], ev_ui$ev_id)
}
logf("bootstrap: %s (%d initiators, %d disputes)", c(pigeonhole = "pigeonhole two-way (initiators x disputes, independent multiplicities)",
     dispute = "one-way, disputes (sensitivity)", initiator = "one-way, initiators (sensitivity)")[[BOOT_UNIT]], length(u_init), length(u_disp))
bs <- list(M1 = NULL, M1s = NULL); bs_coef <- NULL; n_fail <- c(M1 = 0L, M1s = 0L, no_type = 0L)
fail_log <- data.frame(draw = integer(0), model = character(0), reason = character(0))   # GPT review P2-6: replicate IDs and failure reasons
for (bb in seq_len(B)) {
  mult <- draw_mult()
  Lb <- L[mult[as.character(L$ev_id)] > 0L, ]; Lb$w <- Lb$w * mult[as.character(Lb$ev_id)]
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
## GPT recheck P2: the export uses the PRIMARY two-way variance, with initiator-only and target-year columns labelled as sensitivities
coefs <- bind_rows(lapply(c("M0", "M1", "M1s"), function(mm) { m <- get(mm); b <- coef(m); V3 <- if (mm == "M1") VS else vc_set(m)
  se2 <- sqrt(diag(V3$twoway)); se1 <- sqrt(diag(V3$init)); se3 <- sqrt(diag(V3$tgtyr))
  data.frame(model = mm, term = names(b), coef = b, se_init_x_dispute = se2, p_init_x_dispute = 2 * pnorm(-abs(b / se2)),
             se_initiator = se1, p_initiator = 2 * pnorm(-abs(b / se1)), se_target_year = se3, p_target_year = 2 * pnorm(-abs(b / se3)), row.names = NULL) }))
## per-type axis tests: does this axis matter for this type at all (3 df: linear, squared, cross = 0)
Rz <- function(terms) { R <- matrix(0, length(terms), length(nm1), dimnames = list(NULL, nm1)); for (i in seq_along(terms)) R[i, terms[i]] <- 1; R }
axis_rows <- NULL
for (g in ROLES) for (ax in list(material = c("dM", "dM2", "dMR"), relational = c("dR", "dR2", "dMR"))) {
  tm <- paste0(ax, "_", g); lab_med <- ""
  if (ENTRY_MOD && g == "join") {   # GPT review P1: the joiner axis has no effect at ANY lag only if its three *_jlag terms are zero too
    w_med <- wald3(Rz(tm))
    axis_rows <- rbind(axis_rows, data.frame(test = sprintf("%s: %s axis has no effect at the median joiner lag (3 df)", g, if (ax[1] == "dM") "material" else "relational"),
                                             chi2 = w_med[["chi2"]], df = w_med[["df"]], p_init_x_dispute = w_med[["p"]], p_initiator = w_med[["p_initiator"]], p_target_year = w_med[["p_tgtyr"]]))
    tm <- c(tm, paste0(ax, "_jlag")); lab_med <- ", any entry lag"
  }
  w0 <- wald3(Rz(tm))
  axis_rows <- rbind(axis_rows, data.frame(test = sprintf("%s: %s axis has no effect on target choice%s (%d df)", g, if (ax[1] == "dM") "material" else "relational", lab_med, length(tm)),
                                           chi2 = w0[["chi2"]], df = w0[["df"]], p_init_x_dispute = w0[["p"]], p_initiator = w0[["p_initiator"]], p_target_year = w0[["p_tgtyr"]]))
}
tests <- data.frame(test = c(sprintf("utility surfaces differ by type (M1 vs M0, %d df%s)", W_roles[["df"]], if (ENTRY_MOD) ", incl. the joiner entry-lag terms" else ""),
                             sprintf("utility interaction (cross coefficient) differs by type (%d df; NOT equal rho)", W_cross[["df"]])),
                    chi2 = c(W_roles[["chi2"]], W_cross[["chi2"]]), df = c(W_roles[["df"]], W_cross[["df"]]),
                    p_init_x_dispute = c(W_roles[["p"]], W_cross[["p"]]), p_initiator = c(W_roles[["p_initiator"]], W_cross[["p_initiator"]]),
                    p_target_year = c(W_roles[["p_tgtyr"]], W_cross[["p_tgtyr"]]))
tests <- rbind(tests, axis_rows)
if (!is.null(W_roles_med)) tests <- rbind(tests, data.frame(test = sprintf("utility surfaces differ by type at the median joiner lag (%d df)", W_roles_med[["df"]]), chi2 = W_roles_med[["chi2"]],
                                                            df = W_roles_med[["df"]], p_init_x_dispute = W_roles_med[["p"]], p_initiator = W_roles_med[["p_initiator"]], p_target_year = W_roles_med[["p_tgtyr"]]))
emod <- NULL
if (ENTRY_MOD) {
  tests <- rbind(tests, data.frame(test = "joiner gap terms do not change with log(1 + entry lag in months) (5 df)", chi2 = W_mod[["chi2"]], df = W_mod[["df"]],
                                   p_init_x_dispute = W_mod[["p"]], p_initiator = W_mod[["p_initiator"]], p_target_year = W_mod[["p_tgtyr"]]))
  ## joiner endpoint utilities (other gap at 0) at entry lags of 0, 1, 3, 12 and 36 months, and joiner - bilateral; percentile intervals when boot_ok
  ep <- function(b, lagd) { mm <- log1p(lagd) - M_CTR; jb <- function(v) b[[paste0(v, "_join")]] + b[[paste0(v, "_jlag")]] * mm
    c(join_M_weaker = jb("dM") + jb("dM2"), join_M_stronger = -jb("dM") + jb("dM2"), join_R_less_central = jb("dR") + jb("dR2"), join_R_more_central = -jb("dR") + jb("dR2"),
      diff_M_weaker = jb("dM") + jb("dM2") - b[["dM_bil"]] - b[["dM2_bil"]], diff_M_stronger = -jb("dM") + jb("dM2") + b[["dM_bil"]] - b[["dM2_bil"]],
      diff_R_less_central = jb("dR") + jb("dR2") - b[["dR_bil"]] - b[["dR2_bil"]], diff_R_more_central = -jb("dR") + jb("dR2") + b[["dR_bil"]] - b[["dR2_bil"]]) }
  for (lagd in c(0, 1, 3, 12, 36)) {
    e0 <- ep(coef(M1), lagd)
    q <- if (isTRUE(boot_ok[["M1"]])) apply(t(apply(bs_coef[, -1, drop = FALSE], 1, function(r) ep(r, lagd))), 2, quantile, c(.025, .975)) else matrix(NA, 2, length(e0))
    jlag <- lag_ev[keep & t_ev == 2L]; nsup <- sum(abs(log1p(jlag) - log1p(lagd)) <= 0.5)   # GPT review P2: joiner events near this lag (|log(1+lag) diff| <= 0.5)
    emod <- rbind(emod, data.frame(lag_months = lagd, quantity = names(e0), estimate = e0, lo = q[1, ], hi = q[2, ], n_joiners_near_lag = nsup, row.names = NULL))
  }
}
od <- file.path(v2_root(), "output"); dir.create(od, showWarnings = FALSE)
write.csv(tab, file.path(od, paste0(cell, "_derived.csv")), row.names = FALSE)
write.csv(coefs, file.path(od, paste0(cell, "_coefs.csv")), row.names = FALSE)
write.csv(tests, file.path(od, paste0(cell, "_tests.csv")), row.names = FALSE)
write.csv(pw, file.path(od, paste0(cell, "_pairwise.csv")), row.names = FALSE)
if (!is.null(loo)) write.csv(loo, file.path(od, paste0(cell, "_loo.csv")), row.names = FALSE)
if (!is.null(emod)) write.csv(emod, file.path(od, paste0(cell, "_entrymod.csv")), row.names = FALSE)
if (!is.null(role_alt)) write.csv(role_alt, file.path(od, paste0(cell, "_roleswitch.csv")), row.names = FALSE)   # events whose role differs under the other date bound
saveRDS(list(pt = pt, bs = bs, bs_coef = bs_coef, fail_log = fail_log, n_ok = n_ok, n_fail = n_fail, B = B, spec = spec, types = TYPES, roles = ROLES,
             support = SUPPORT, n_rows = nrow(L), E = sum(keep), V = VS, boot_unit = BOOT_UNIT, loo = loo, coorig_months = COORIG_W, entry_mod = ENTRY_MOD,
             lag_ev = lag_ev, lag_lo_ev = lag_lo_ev, lag_hi_ev = lag_hi_ev, role_alt = role_alt, coorig_bound = COORIG_BOUND, m_ctr = if (ENTRY_MOD) M_CTR else NA_real_),
        file.path(od, paste0(cell, "_boot.rds")))
## GPT recheck P2: a manifest names exactly the artifacts written by THIS run (an older optional file, e.g. _loo.csv, is not part of it)
written <- c(paste0(cell, c("_derived.csv", "_coefs.csv", "_tests.csv", "_pairwise.csv", "_boot.rds")), if (!is.null(loo)) paste0(cell, "_loo.csv"),
             if (!is.null(emod)) paste0(cell, "_entrymod.csv"), if (!is.null(role_alt)) paste0(cell, "_roleswitch.csv"))
gitc <- tryCatch(system2("git", c("-C", R_DIR, "rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE), error = function(e) NA_character_, warning = function(w) NA_character_)
dirty <- tryCatch(length(system2("git", c("-C", R_DIR, "status", "--porcelain", "--untracked-files=no", "--", "."), stdout = TRUE, stderr = FALSE)) > 0, error = function(e) NA, warning = function(w) NA)
write.csv(data.frame(file = written, run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), git = if (length(gitc)) gitc[1] else NA_character_, git_dirty = dirty,
                     args = paste(args, collapse = " ")), file.path(od, paste0(cell, "_manifest.csv")), row.names = FALSE)
print(tab, right = FALSE, row.names = FALSE); print(tests, row.names = FALSE)
logf("done: %s_{derived,coefs,tests,pairwise}.csv and _boot.rds", cell)
