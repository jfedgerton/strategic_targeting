## R/s2_01_build_riskset.R -- stage-2 (joining) risk set on real data. Contract: docs/stage2_design_memo.md section 13 (rev 6).
## Builds the dispute table, participant classification (windows 0-6), prospective cohort with horizons, one-shot rows and monthly rows at
## the reference window, rosters (carry-forward; end-policy census), main-measure ties, tiers, folds, side-level jointly-fatal status and
## covariates; runs the R1/N2 acceptance tests against a stage-1 bcjs dump. NO MODEL IS FIT.
## Usage:
##   Rscript R/s2_01_build_riskset.R --data_entry=DIR --layers=layers_dy.rds --out=DIR [--dump_bcjs=long_bcjs_bon_w1d_relyear_fmid.rds]
##          [--verify_positions=positions.rds --verify_col='W2contig|noPTA|pr'] [--window=1]
suppressPackageStartupMessages({ library(data.table); library(peacesciencer) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file)), "lib", "s2.R"))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
known <- c("data_entry", "layers", "out", "dump_bcjs", "verify_positions", "verify_col", "window", "fatal_def")
fl <- sub("^--([a-z_]+)=.*$", "\\1", args); if (any(!fl %in% known) || anyDuplicated(fl)) stop("unknown or duplicated option: ", paste(fl, collapse = " "))
DE <- opt("data_entry"); LAYERS_RDS <- opt("layers"); OUT <- opt("out"); W_REF <- as.integer(opt("window", "1")); FATAL_DEF <- opt("fatal_def", "gml")
stopifnot(!is.null(DE), !is.null(LAYERS_RDS), !is.null(OUT), W_REF %in% 0:6, FATAL_DEF %in% c("gml", "stage1"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT, "s2_build.log"); cat("", file = LOG)
logf <- function(...) { s <- sprintf(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOG, append = TRUE) }
logf("stage-2 risk-set build %s; data_entry %s; layers %s; reference window %d", format(Sys.time()), DE, LAYERS_RDS, W_REF)
logf("git commit %s", tryCatch(system2("git", c("-C", dirname(normalizePath(this_file)), "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown"))
set.seed(123)

## ---- 0. inputs ---------------------------------------------------------------------------------------------------------------------
P <- as.data.table(readRDS(file.path(DE, "participants_timing.rds")))
EV <- readRDS(file.path(DE, "events.rds"))
DY <- as.data.table(readRDS(LAYERS_RDS))
GP <- as.data.table(gml_part); GD <- as.data.table(gml_mid_disps); GDY <- as.data.table(gml_mid_ddydisps)
logf("participants_timing %d spells / %d disputes; gml_part %d rows; gml_mid_disps %d disputes", nrow(P), uniqueN(P$dispnum), nrow(GP), nrow(GD))

## ---- 1. dispute sample (B5) ----------------------------------------------------------------------------------------------------
## --fatal_def=gml (memo B5 text): GML dispute fatality > 0. --fatal_def=stage1: the primary disputes of the stage-1 fatal events (events.rds
## in_fmid, the COW-based source chain stage 1 samples on). The memo assumed these coincide; they do not (overlap logged), so both are built.
fat_g <- GD[fatality > 0, dispnum]; fat_s1 <- unique(as.integer(as.data.table(EV)[in_fmid == TRUE, dispnum]))
logf("fatal dispute sets: GML %d; stage-1 fatal-event primary disputes %d; both %d; GML only %d; stage-1 only %d; using --fatal_def=%s",
     length(fat_g), length(fat_s1), length(intersect(fat_g, fat_s1)), length(setdiff(fat_g, fat_s1)), length(setdiff(fat_s1, fat_g)), FATAL_DEF)
fat_d <- if (FATAL_DEF == "gml") fat_g else fat_s1
logf("GML disputes: fatal %d, nonfatal %d, unknown (-9) %d (unknown dropped, counted)", length(fat_g), sum(GD$fatality == 0), sum(GD$fatality < 0))
miss <- setdiff(fat_d, P$dispnum); if (length(miss)) logf("fatal disputes absent from participants_timing (dropped): %s", paste(miss, collapse = " "))
fat_d <- intersect(fat_d, P$dispnum)
## consistency with the stage-1 event file's dispute fatality (events.rds fatality_disp)
ev_f <- unique(as.data.table(EV)[!EV$escalation & n_disputes == 1L, .(dispnum, fatality_disp)])
cmp <- merge(ev_f, GD[, .(dispnum, fatality)], by = "dispnum")
logf("check: events.rds fatality_disp vs gml_mid_disps fatality on %d single-source disputes: %d disagree", nrow(cmp), sum(cmp$fatality_disp != cmp$fatality))

## ---- 2. dispute table: reference onset, end month, initiating side, fold ---------------------------------------------------------------
ref <- P[orig == 1L, .(onset_lo = min(entry_lo), onset_hi = min(entry_hi)), by = dispnum]           # stage-1 TIMING_REF = originator
GP[, end_mk := s2_mk(endyear, endmon)]; GP[, st_mk := s2_mk(styear, stmon)]
gend <- GP[, .(end_gml = max(end_mk)), by = dispnum]
P[, exit_fallback := s2_mk(spell_last, 12L)]
pend <- P[, .(end_p = max(exit_fallback)), by = dispnum]
DT <- merge(merge(data.table(dispnum = fat_d), ref, by = "dispnum", all.x = TRUE), gend, by = "dispnum", all.x = TRUE)
DT <- merge(DT, pend, by = "dispnum", all.x = TRUE)
stopifnot(!anyNA(DT$onset_lo))
DT[, end_mk := ifelse(is.na(end_gml), end_p, end_gml)]
logf("dispute end month: GML participant end dates %d disputes; fallback (Dec of last spell year) %d", sum(!is.na(DT$end_gml)), sum(is.na(DT$end_gml)))
bad_end <- DT[end_mk < onset_lo]; if (nrow(bad_end)) { logf("disputes ending before their reference onset (dropped): %s", paste(bad_end$dispnum, collapse = " ")); DT <- DT[end_mk >= onset_lo] }
IS <- as.data.table(s2_init_side(P, EV, cow_war_inter, DT$dispnum))
DT <- merge(DT, IS, by = "dispnum")
logf("initiating side by rule: 1 (originator events) %d, 2 (Side A, not war-linked) %d, 3 (war rule) %d, 4 (war ambiguous -> Side A) %d, 5 (mixed-source war link, unresolved) %d; reviewed %d (hard %d)",
     sum(DT$side_rule == 1), sum(DT$side_rule == 2), sum(DT$side_rule == 3), sum(DT$side_rule == 4), sum(DT$side_rule == 5), sum(DT$side_case != ""), sum(DT$side_case == "hard"))
if (any(DT$side_rule == 5L)) { logf("dropped (unresolved side): %s", paste(DT[side_rule == 5L, dispnum], collapse = " ")); DT <- DT[side_rule != 5L] }
DT[, onset_year := s2_mk_year(onset_lo)]; DT[, fold := s2_fold(onset_year)]
logf("disputes in the stage-2 sample: %d; by fold: %s", nrow(DT), paste(sprintf("%d:%d", 1:6, tabulate(DT$fold, 6)), collapse = " "))

## ---- 3. participant spells with exit months -------------------------------------------------------------------------------------
PT <- P[dispnum %in% DT$dispnum, .(dispnum, ccode, spell, spell_first, spell_last, side = side_entry, orig, entry_lo, entry_hi, side_switch, exit_fallback)]
gx <- GP[, .(dispnum, ccode, g_sty = styear, g_endy = endyear, end_mk)]
mx <- gx[PT, on = .(dispnum, ccode, g_sty <= spell_last, g_endy >= spell_first), nomatch = NA, .(dispnum, ccode, spell = i.spell, end_mk = x.end_mk)]
mx <- mx[, .(exit_gml = if (all(is.na(end_mk))) NA_integer_ else max(end_mk, na.rm = TRUE)), by = .(dispnum, ccode, spell)]
PT <- merge(PT, mx, by = c("dispnum", "ccode", "spell"))
PT[, exit_mk := ifelse(is.na(exit_gml), exit_fallback, exit_gml)]
n_inv <- sum(PT$exit_mk < PT$entry_hi); PT[exit_mk < entry_hi, exit_mk := entry_hi]
logf("participant spells %d (%d participants); exit month from GML %d, fallback Dec of last year %d; exit before definite entry (set to entry) %d; side switches %d (first side kept); second spells %d (ignored as entries, kept in rosters)",
     nrow(PT), uniqueN(PT[, .(dispnum, ccode)]), sum(!is.na(PT$exit_gml)), sum(is.na(PT$exit_gml)), n_inv, sum(PT$side_switch), sum(PT$spell > 1L))
## target set and initiating side nonempty (B3)
CLr <- s2_classify(PT, DT, W_REF)
ts <- CLr[baseline == TRUE, .(nA = sum(side_rel == "A"), nB = sum(side_rel == "B")), by = dispnum]
no_ts <- setdiff(DT$dispnum, ts[nA > 0 & nB > 0, dispnum])
if (length(no_ts)) { logf("disputes without a baseline participant on one side (dropped, B3): %s", paste(no_ts, collapse = " ")); DT <- DT[!dispnum %in% no_ts]; PT <- PT[!dispnum %in% no_ts] }
logf("baseline sets: disputes with >1 initiating-side originator %d, >1 target-set member %d", sum(ts$nA > 1), sum(ts$nB > 1))

## ---- 4. classification at windows 0-6 and acceptance tests (R1, N2) ----------------------------------------------------------------
CLW <- rbindlist(lapply(0:6, function(w) s2_classify(PT, DT, w)[, W := w]))
tab <- CLW[baseline == FALSE, .N, by = .(W, cat)]; tab <- dcast(tab, W ~ cat, value.var = "N", fill = 0)
logf("non-baseline participants by window (A / B / orig): %s", paste(sprintf("W%d %d/%d/%d", tab$W, tab$A, tab$B, if ("orig" %in% names(tab)) tab$orig else 0L), collapse = "; "))
CL <- CLW[W == W_REF][, W := NULL]
sm <- CL[baseline == FALSE & m_entry == 0L, .N, by = side_rel]
logf("same-month (onset-month) non-original entrants by side (R1 count): %s", paste(sprintf("%s %d", sm$side_rel, sm$N), collapse = ", "))
DUMP <- opt("dump_bcjs")
if (!is.null(DUMP)) {
  d <- readRDS(DUMP); stopifnot(d$types == "bcjs", as.integer(d$coorig_months) == W_REF)
  ch <- unique(as.data.table(d$D)[chosen == 1L, .(c1, disp, yr, type, orig_joiner)])
  ## stage 1 classifies each event by the participant spell covering the event year; stage 2 classifies the FIRST spell (re-entries ignored,
  ## memo section 4). Events from later spells are a stated projection difference, reported separately and excluded from the tests.
  sp <- PT[ch, on = .(dispnum = disp, ccode = c1, spell_first <= yr, spell_last >= yr), nomatch = NA, .(disp = i.disp, c1 = i.c1, yr = i.yr, spell = x.spell)]
  sp <- sp[, .(spell = min(spell)), by = .(disp, c1, yr)]
  ch <- merge(ch, sp, by = c("disp", "c1", "yr"), all.x = TRUE)
  logf("acceptance: stage-1 events from a later participant spell (excluded; re-entries are not stage-2 entries): %d [%s]", ch[!is.na(spell) & spell > 1L, .N],
       paste(ch[!is.na(spell) & spell > 1L, sprintf("%d in %d (%d)", c1, disp, yr)], collapse = "; "))
  ch <- unique(ch[is.na(spell) | spell == 1L, .(c1, disp, type, orig_joiner)])
  X <- merge(ch, CL, by.x = c("disp", "c1"), by.y = c("dispnum", "ccode"))
  logf("acceptance: stage-1 bcjs W%d events %d (unique initiator-dispute); in the stage-2 sample %d (others: dispute not in the GML-fatal sample)", W_REF, nrow(ch), nrow(X))
  a_base <- X[, sum((!orig_joiner) != baseline)]
  a_prom <- X[orig_joiner == TRUE, sum((type == 1L) != (cat == "orig"))]
  a_side <- X[orig_joiner == TRUE & type %in% 2:3, sum(ifelse(type == 2L, "A", "B") != side_rel)]
  logf("acceptance (a) baseline roles: mismatches %d of %d", a_base, nrow(X))
  logf("acceptance (b) promotions at W%d: mismatches %d of %d original joiners (stage-1 promoted %d, stage-2 orig %d)", W_REF, a_prom, sum(X$orig_joiner), X[orig_joiner == TRUE, sum(type == 1L)], X[orig_joiner == TRUE, sum(cat == "orig")])
  logf("acceptance (c) joiner sides (stage-1 joinA/joinB vs stage-2 A/B): mismatches %d of %d", a_side, X[orig_joiner == TRUE & type %in% 2:3, .N])
  nofit <- CL[baseline == FALSE & !paste(dispnum, ccode) %in% paste(ch$disp, ch$c1)]
  logf("N2: non-baseline participants without a fitted stage-1 event, classified by the participant-level rule: %d (A %d, B %d, orig %d)", nrow(nofit), sum(nofit$cat == "A"), sum(nofit$cat == "B"), sum(nofit$cat == "orig"))
  if (a_base + a_prom + a_side > 0) { mm <- X[((!orig_joiner) != baseline) | (orig_joiner & (type == 1L) != (cat == "orig")) | (orig_joiner & type %in% 2:3 & ifelse(type == 2L, "A", "B") != side_rel)]
    fwrite(mm, file.path(OUT, "acceptance_mismatches.csv")); logf("ACCEPTANCE FAILED: mismatches written to acceptance_mismatches.csv") }
}

## ---- 5. cohort, horizons, one-shot and monthly rows (B1, B2, N1) ----------------------------------------------------------------------
CS <- as.data.table(cow_states)[, .(ccode = as.integer(ccode), st_mk = s2_mk(styear, stmonth), end_mk = s2_mk(endyear, endmonth))]
COH <- s2_cohort(DT, CL, CS)
DYN <- s2_dynamic_entrants(COH, CL)
OS <- s2_oneshot(COH, CL)
logf("cohort rows (state x dispute) %d; complete horizon %d; truncated: administrative %d, system exit %d; dynamic entrants (excluded) %d; entries after the horizon (treated as out) %d",
     nrow(COH), sum(COH$complete), sum(COH$trunc_admin), sum(COH$trunc_exit), nrow(DYN), sum(OS$late))
logf("one-shot outcomes W%d, all rows: %s; complete-horizon rows: %s", W_REF, paste(sprintf("%s %d", names(table(OS$y)), table(OS$y)), collapse = ", "),
     paste(sprintf("%s %d", names(table(OS[complete == TRUE]$y)), table(OS[complete == TRUE]$y)), collapse = ", "))
for (w in 0:6) { o <- s2_oneshot(COH, CLW[W == w][, W := NULL]); logf("  window %d: one-shot outcomes %s; retrospective cohort drops %d", w, paste(sprintf("%s %d", names(table(o$y)), table(o$y)), collapse = ", "), sum(o$y == "orig")) }
MON <- s2_monthly(OS, W_REF)
logf("monthly rows W%d: %d (%d state-disputes); outcomes %s; censored last rows %d", W_REF, nrow(MON), nrow(OS), paste(sprintf("%s %d", names(table(MON$y_m)), table(MON$y_m)), collapse = ", "),
     nrow(OS[y == "out"]))
logf("dispute duration (months, reference to end): median %.0f, 90%% %.0f, max %d", median(DT$end_mk - DT$onset_lo), quantile(DT$end_mk - DT$onset_lo, .9), max(DT$end_mk - DT$onset_lo))

## ---- 6. rosters (B3) ---------------------------------------------------------------------------------------------------------------
PS <- merge(PT[, .(dispnum, ccode, spell, exit_mk, entry_hi)], CL[, .(dispnum, ccode, side_rel, baseline)], by = c("dispnum", "ccode"))
PS <- merge(PS, DT[, .(dispnum, onset_lo)], by = "dispnum")
PS[, m_in := ifelse(baseline & spell == 1L, 0L, pmax(entry_hi - onset_lo, 0L))]; PS[, m_out := exit_mk - onset_lo]
MO <- unique(MON[, .(dispnum, m)])
ROS <- s2_rosters(PS, MO, "carry"); ROE <- s2_rosters(PS, MO, "end")
emp <- attr(ROE, "empty")
logf("rosters: dispute-months %d; side-months carried forward (empty side) %d in %d disputes (A %d, B %d)", nrow(MO), nrow(emp), uniqueN(emp$dispnum), sum(emp$side_rel == "A"), sum(emp$side_rel == "B"))

## ---- 7. main-measure ties (section 4, N4) ------------------------------------------------------------------------------------------
## main measure (config.sh MAIN_SPEC): W2contig, layers igo+def+off+neu+nonagg+consul+dca+un, clipped residuals / layer-year max, Euclidean combo
## rescaled within year. Same code path as R/02_build_positions.R edge_w() for W2contig (copied; verified below by recomputing PageRank).
MAIN_L <- c("igo", "def", "off", "neu", "nonagg", "consul", "dca", "un")
setorder(DY, year, ccode1, ccode2)
V <- sapply(MAIN_L, function(l) { v <- as.numeric(DY[[l]]); if (l == "igo") v <- log1p(v); v })
first_year <- sapply(MAIN_L, function(l) { v <- V[, l]; min(DY$year[!is.na(v) & v != 0]) })
yr_split <- split(seq_len(nrow(DY)), DY$year)
yr_has <- sapply(MAIN_L, function(l) { v <- V[, l]; a <- tapply(!is.na(v) & v != 0, DY$year, any); a[as.character(DY$year)] })
EL <- sapply(MAIN_L, function(l) {
  v <- V[, l]; ex <- DY$year >= first_year[[l]] & yr_has[, l]; ok <- if (l == "un") ex & !is.na(v) else ex
  d <- data.frame(v = v, year = DY$year, contig_land = DY$contig_land)[ok, ]
  r <- resid(lm(v ~ contig_land + year, data = d)); x <- numeric(nrow(DY)); x[ok] <- r
  tol <- 1e-8 * max(1, max(abs(v[ok]))); x[abs(x) < tol] <- 0; x <- pmax(x, 0)
  e <- numeric(nrow(DY)); for (ix in yr_split) { j <- ix[ok[ix]]; if (length(j)) { m <- max(x[j]); if (m > 0) e[j] <- x[j] / m } }
  e[!ok] <- 0; e })
S <- sqrt(rowSums(EL^2)); EW <- numeric(length(S))
for (ix in yr_split) { r <- range(S[ix]); if (r[2] - r[1] >= 1e-12) EW[ix] <- (S[ix] - r[1]) / (r[2] - r[1]) }
EDGE <- data.table(a = DY$ccode1, b = DY$ccode2, year = DY$year, w = EW)[w > 0]
logf("main-measure edges: %d positive dyad-years of %d", nrow(EDGE), nrow(DY))
VP <- opt("verify_positions")
if (!is.null(VP)) {
  suppressPackageStartupMessages(library(igraph))
  pz <- readRDS(VP); vc <- opt("verify_col"); stopifnot(!is.null(vc), vc %in% colnames(pz$R))
  ys <- c(1850L, 1914L, 1938L, 1960L, 1990L, 2010L); dev <- c()
  for (y in ys) { k <- pz$key$year == y; vs <- pz$key$ccode[k]; e <- EDGE[year == y]
    g <- graph_from_data_frame(e[, .(as.character(a), as.character(b), weight = w)], directed = FALSE, vertices = data.frame(name = as.character(vs)))
    p <- page_rank(g, directed = FALSE)$vector; dev <- c(dev, max(abs(log(p[as.character(vs)]) - pz$R[k, vc]))) }
  logf("verify: log PageRank from the stage-2 edges vs %s (%s), max abs deviation by year %s: %s", basename(VP), vc, paste(ys, collapse = "/"), paste(sprintf("%.1e", dev), collapse = " "))
  if (any(dev > 1e-6)) stop("stage-2 edges do not reproduce the stage-1 main-measure positions")
}
## one-shot ties: month-0 rosters, edges of onset year - 1 (lag 1) and onset year - 5 (lag 5)
R0 <- ROS[m == 0L]
row0 <- merge(OS[, .(dispnum, ccode)], DT[, .(dispnum, onset_year)], by = "dispnum")[, m := 0L]
T1 <- s2_ties(row0[, .(dispnum, ccode, m, year_tie = onset_year - 1L)], R0, EDGE)
T5 <- s2_ties(row0[, .(dispnum, ccode, m, year_tie = onset_year - 5L)], R0, EDGE)
setnames(T5, c("tie_sum_A", "tie_max_A", "tie_sum_B", "tie_max_B"), paste0(c("tie_sum_A", "tie_max_A", "tie_sum_B", "tie_max_B"), "_l5"))
n_pre <- sum(row0$onset_year - 1L < min(DY$year)); logf("one-shot rows whose lag-1 tie year precedes the layers (ties 0, flagged): %d", n_pre)
## monthly ties: rosters at m - 1 (lagged), edges of (calendar year of month m) - 1
MON <- merge(MON, DT[, .(dispnum, onset_lo)], by = "dispnum")
MON[, year_tie := s2_mk_year(onset_lo + m) - 1L]
TM <- s2_ties(MON[, .(dispnum, ccode, m, year_tie)], ROS, EDGE)
MON <- merge(MON, TM[, .(dispnum, ccode, m, tie_sum_A, tie_max_A, tie_sum_B, tie_max_B)], by = c("dispnum", "ccode", "m"))
logf("ties: one-shot rows with any positive tie to A %d / to B %d of %d; monthly rows %d / %d of %d", sum(T1$tie_max_A > 0), sum(T1$tie_max_B > 0), nrow(T1),
     sum(MON$tie_max_A > 0), sum(MON$tie_max_B > 0), nrow(MON))

## ---- 8. tiers (frozen at onset), covariates, jointly fatal ---------------------------------------------------------------------------
BASE <- CL[baseline == TRUE, .(dispnum, ccode)]
DYk <- rbind(DY[, .(k = ccode1, j = ccode2, year, conttype, capdist, ally)], DY[, .(k = ccode2, j = ccode1, year, conttype, capdist, ally)])
pr <- merge(merge(COH[, .(dispnum, ccode)], BASE[, .(dispnum, j = ccode)], by = "dispnum", allow.cartesian = TRUE), DT[, .(dispnum, onset_year)], by = "dispnum")
pr <- merge(pr, DYk[, .(k, j, year, conttype)], by.x = c("ccode", "j", "onset_year"), by.y = c("k", "j", "year"), all.x = TRUE)
pr[, ty := onset_year - 1L]
pr <- merge(pr, DYk[, .(k, j, year, ally)], by.x = c("ccode", "j", "ty"), by.y = c("k", "j", "year"), all.x = TRUE)
MJ <- as.data.table(cow_majors)[, .(ccode, st_mk = s2_mk(styear, stmonth), end_mk = s2_mk(endyear, endmonth))]
maj <- MJ[merge(COH[, .(dispnum, ccode)], DT[, .(dispnum, onset_lo)], by = "dispnum"), on = .(ccode, st_mk <= onset_lo, end_mk >= onset_lo), nomatch = 0L, .(dispnum = i.dispnum, ccode)]
TIER <- s2_tiers(COH, BASE, CONTIG = pr[!is.na(conttype) & conttype %in% 1:5, .(dispnum, ccode, j)], MAJ = unique(maj), ALLY = pr[!is.na(ally) & ally > 0, .(dispnum, ccode, j)])
logf("relevance tiers (frozen at onset): likely %d of %d cohort rows (contiguous %d, same region %d, major %d, allied at t-1 %d)", sum(TIER$likely), nrow(TIER), sum(TIER$any_contig), sum(TIER$any_reg), sum(TIER$major), sum(TIER$any_ally))
## side covariates at onset (one-shot): CINC ratio and side share (t - 1), minimum capital distance, alliance with any side member (t - 1)
NMC <- as.data.table(cow_nmc)[, .(ccode, year, cinc)]
sx <- merge(merge(row0[, .(dispnum, ccode, onset_year)], R0[, .(dispnum, j = ccode, side_rel)], by = "dispnum", allow.cartesian = TRUE)[ccode != j], NMC, by.x = c("j", "onset_year"), by.y = c("ccode", "year"), all.x = TRUE)
sx[, ty := onset_year - 1L]
sx <- merge(sx, DYk[, .(k, j, year, capdist)], by.x = c("ccode", "j", "onset_year"), by.y = c("k", "j", "year"), all.x = TRUE)
sx <- merge(sx, DYk[, .(k, j, year, ally_t1 = ally)], by.x = c("ccode", "j", "ty"), by.y = c("k", "j", "year"), all.x = TRUE)
## side CINC uses t - 1: re-merge cinc at ty
sx[, cinc := NULL]; sx <- merge(sx, NMC, by.x = c("j", "ty"), by.y = c("ccode", "year"), all.x = TRUE)
SC <- sx[, .(side_cinc = if (all(is.na(cinc))) NA_real_ else sum(cinc, na.rm = TRUE), side_cinc_na = sum(is.na(cinc)), min_capdist = suppressWarnings(min(capdist, na.rm = TRUE)), ally_side = as.integer(any(ally_t1 > 0, na.rm = TRUE))), by = .(dispnum, ccode, side_rel)]
SC[!is.finite(min_capdist), min_capdist := NA_real_]
SC <- dcast(SC, dispnum + ccode ~ side_rel, value.var = c("side_cinc", "side_cinc_na", "min_capdist", "ally_side"))
own <- merge(row0[, .(dispnum, ccode, ty = onset_year - 1L)], NMC, by.x = c("ccode", "ty"), by.y = c("ccode", "year"), all.x = TRUE)
SC <- merge(SC, own[, .(dispnum, ccode, cinc_k = cinc)], by = c("dispnum", "ccode"))
SC[, `:=`(cinc_ratio_A = ifelse(side_cinc_A > 0, cinc_k / side_cinc_A, NA_real_), cinc_ratio_B = ifelse(side_cinc_B > 0, cinc_k / side_cinc_B, NA_real_),
          share_A = ifelse(side_cinc_A + side_cinc_B > 0, side_cinc_A / (side_cinc_A + side_cinc_B), NA_real_))]
logf("covariates: side CINC zero or wholly missing (ratio NA, never 0/0) A %d / B %d", sum(!(SC$side_cinc_A > 0), na.rm = TRUE) + sum(is.na(SC$side_cinc_A)), sum(!(SC$side_cinc_B > 0), na.rm = TRUE) + sum(is.na(SC$side_cinc_B)))
SC[, share_B := 1 - share_A]
logf("covariates: rows with missing own CINC %d; side CINC with a missing member %d (A) / %d (B); missing min capital distance A %d / B %d",
     sum(is.na(SC$cinc_k)), sum(SC$side_cinc_na_A > 0), sum(SC$side_cinc_na_B > 0), sum(is.na(SC$min_capdist_A)), sum(is.na(SC$min_capdist_B)))
## target centrality (main measure, t - 1): mean and minimum over the original target set
POSV <- if (!is.null(VP)) data.table(pz$key, R = pz$R[, opt("verify_col")]) else NULL
if (!is.null(POSV)) {
  tgt <- merge(merge(BASE, CL[baseline == TRUE, .(dispnum, ccode, side_rel)], by = c("dispnum", "ccode"))[side_rel == "B"], DT[, .(dispnum, ty = onset_year - 1L)], by = "dispnum")
  tgt <- merge(tgt, POSV, by.x = c("ccode", "ty"), by.y = c("ccode", "year"), all.x = TRUE)
  TC <- tgt[, .(tgt_cent_mean = mean(R), tgt_cent_min = min(R)), by = dispnum]
  logf("target centrality: disputes with a missing t-1 position for a target-set member %d", sum(is.na(TC$tgt_cent_mean)))
} else TC <- data.table(dispnum = DT$dispnum, tgt_cent_mean = NA_real_, tgt_cent_min = NA_real_)
## jointly fatal (side level) from participant fatality codes (gml_mid_ddydisps, participant = ccode1; max over its rows; -9 unknown)
PF <- GDY[, .(fat = as.integer(if (any(fatality1 > 0)) max(fatality1) else if (all(fatality1 == 0)) 0 else -9)), by = .(dispnum, ccode = ccode1)]
CLall <- CL[, .(dispnum, ccode, side_rel)]
JF <- s2_jfatal(CLall, PF)
logf("jointly fatal (side level): yes %d, no %d, unresolved %d; participants without a fatality code (unknown) %d",
     sum(JF$jfatal == "yes"), sum(JF$jfatal == "no"), sum(JF$jfatal == "unresolved"), nrow(CLall[!PF, on = .(dispnum, ccode)]))
DT <- merge(merge(DT, JF, by = "dispnum"), TC, by = "dispnum", all.x = TRUE)
DT <- merge(DT, GD[, .(dispnum, fat_level = fatality)], by = "dispnum", all.x = TRUE)

## ---- 9. assemble and write ------------------------------------------------------------------------------------------------------
ONE <- Reduce(function(a, b) merge(a, b, by = c("dispnum", "ccode"), all.x = TRUE),
              list(OS, T1[, .(dispnum, ccode, tie_sum_A, tie_max_A, tie_sum_B, tie_max_B)], T5[, .(dispnum, ccode, tie_sum_A_l5, tie_max_A_l5, tie_sum_B_l5, tie_max_B_l5)],
                   TIER[, .(dispnum, ccode, likely, any_contig, any_reg, major, any_ally)], SC))
ONE <- merge(ONE, DT[, .(dispnum, fold, onset_year, jfatal, fat_level, tgt_cent_mean, tgt_cent_min, side_case)], by = "dispnum")
stopifnot(nrow(ONE) == nrow(OS), !anyDuplicated(ONE[, .(dispnum, ccode)]))
MON <- merge(MON, TIER[, .(dispnum, ccode, likely)], by = c("dispnum", "ccode"))
MON <- merge(MON, DT[, .(dispnum, fold)], by = "dispnum")
saveRDS(DT, file.path(OUT, "s2_disputes.rds")); saveRDS(PT, file.path(OUT, "s2_participant_spells.rds")); saveRDS(CLW, file.path(OUT, "s2_classification_w0_6.rds"))
saveRDS(ONE, file.path(OUT, sprintf("s2_oneshot_w%d.rds", W_REF))); saveRDS(MON, file.path(OUT, sprintf("s2_monthly_w%d.rds", W_REF)))
saveRDS(ROS, file.path(OUT, "s2_rosters_carry.rds")); saveRDS(EDGE, file.path(OUT, "s2_main_edges.rds"))
logf("written: s2_disputes (%d), s2_oneshot_w%d (%d), s2_monthly_w%d (%d), rosters (%d), edges (%d)", nrow(DT), W_REF, nrow(ONE), W_REF, nrow(MON), nrow(ROS), nrow(EDGE))
logf("S2_BUILD_DONE")
