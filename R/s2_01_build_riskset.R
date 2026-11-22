## R/s2_01_build_riskset.R -- stage-2 (joining) risk set on real data. Contract: docs/stage2_design_memo.md section 13 (rev 6).
## Builds the dispute table, participant classification (windows 0-6), prospective cohort with horizons, one-shot rows and monthly rows at
## the reference window, rosters (carry-forward; end-policy census), main-measure ties, tiers, folds, side-level jointly-fatal status and the
## one-shot reference covariates; reconciles every stage-1 projection record. NO MODEL IS FIT.
## Revision 2026-11-22 (GPT build review, docs/REVIEW_stage2_build_gpt.md): P1-1 spell exits capped at each spell's last year; P1-2 required
## positions input and the full one-shot reference covariate set with a schema assertion; P2-3 fail-closed acceptance with explicit
## reconciliation reasons (--mode=release); P2-4 missing vs zero (tie coverage, strict side aggregates, t-1 geography); P2-5 PI decision:
## the stage-1 fatal set is the reference (--fatal_def=stage1 default), ID lists, fatality-code decomposition, imputed-end flags; P2-6 broader
## acceptance (full-universe baselines, promotion boundaries, uncertain dates, all-year PageRank check).
## Usage:
##   Rscript R/s2_01_build_riskset.R --data_entry=DIR --layers=layers_dy.rds --positions=positions.rds --poscol='W2contig|noPTA|pr' --out=DIR
##          --dump_bcjs=long_bcjs_bon_w1d_relyear_fmid.rds [--mode=release|exploratory] [--window=1] [--fatal_def=stage1|gml]
suppressPackageStartupMessages({ library(data.table); library(peacesciencer); library(igraph) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file)), "lib", "s2.R"))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
known <- c("data_entry", "layers", "positions", "poscol", "out", "dump_bcjs", "mode", "window", "fatal_def")
fl <- sub("^--([a-z_]+)=.*$", "\\1", args); if (any(!fl %in% known) || anyDuplicated(fl)) stop("unknown or duplicated option: ", paste(fl, collapse = " "))
DE <- opt("data_entry"); LAYERS_RDS <- opt("layers"); POSF <- opt("positions"); POSCOL <- opt("poscol"); OUT <- opt("out")
W_REF <- as.integer(opt("window", "1")); FATAL_DEF <- opt("fatal_def", "stage1"); MODE <- opt("mode", "release"); DUMP <- opt("dump_bcjs")
stopifnot(!is.null(DE), !is.null(LAYERS_RDS), !is.null(POSF), !is.null(POSCOL), !is.null(OUT), W_REF %in% 0:6, FATAL_DEF %in% c("gml", "stage1"),
          MODE %in% c("release", "exploratory"))
if (MODE == "release" && is.null(DUMP)) stop("--mode=release requires --dump_bcjs (acceptance is mandatory for a release build)")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT, "s2_build.log"); cat("", file = LOG)
logf <- function(...) { s <- sprintf(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOG, append = TRUE) }
fail <- function(...) { s <- sprintf(...); logf("FAILED: %s", s); if (MODE == "release") stop(s, call. = FALSE) else logf("(exploratory mode: continuing)") }
GIT <- tryCatch(system2("git", c("-C", dirname(normalizePath(this_file)), "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
logf("stage-2 risk-set build %s; mode %s; fatal_def %s; reference window %d; git %s", format(Sys.time()), MODE, FATAL_DEF, W_REF, GIT)
logf("inputs: data_entry %s; layers %s; positions %s [%s]; dump %s", DE, LAYERS_RDS, POSF, POSCOL, if (is.null(DUMP)) "none" else DUMP)
set.seed(123)

## ---- 0. inputs ---------------------------------------------------------------------------------------------------------------------
P <- as.data.table(readRDS(file.path(DE, "participants_timing.rds")))
EV <- as.data.table(readRDS(file.path(DE, "events.rds")))
DY <- as.data.table(readRDS(LAYERS_RDS))
pz <- readRDS(POSF); stopifnot(POSCOL %in% colnames(pz$R)); POS <- data.table(pz$key, R = pz$R[, POSCOL]); stopifnot(!anyDuplicated(POS[, .(ccode, year)]))
GP <- as.data.table(gml_part); GD <- as.data.table(gml_mid_disps); GDY <- as.data.table(gml_mid_ddydisps)
logf("participants_timing %d spells / %d disputes; gml_part %d rows; gml_mid_disps %d disputes; positions %d state-years", nrow(P), uniqueN(P$dispnum), nrow(GP), nrow(GD), nrow(POS))
GP[, end_mk := s2_mk(endyear, endmon)]                           # s2_mk stops on invalid month codes (P2-5)

## ---- 1. dispute population (B5; PI 2026-11-22: stage-1 fatal set is the reference, GML a labelled alternative) ---------------------
fat_g <- GD[fatality > 0, dispnum]
fat_s1 <- unique(as.integer(EV[in_fmid == TRUE, dispnum]))      # "primary disputes represented by stage-1 fatal-event flags"
SETS <- data.table(dispnum = sort(unique(c(fat_g, fat_s1))))[, `:=`(in_gml = dispnum %in% fat_g, in_stage1 = dispnum %in% fat_s1)]
fwrite(SETS, file.path(OUT, "fatal_sets.csv"))
logf("population: GML-fatal %d; stage-1 fatal-event primary disputes %d; intersection %d; GML only %d; stage-1 only %d (ID lists in fatal_sets.csv); using %s",
     length(fat_g), length(fat_s1), sum(SETS$in_gml & SETS$in_stage1), sum(SETS$in_gml & !SETS$in_stage1), sum(!SETS$in_gml & SETS$in_stage1), FATAL_DEF)
fat_d <- if (FATAL_DEF == "gml") fat_g else fat_s1
miss <- setdiff(fat_d, P$dispnum); if (length(miss)) logf("availability exclusion: disputes absent from participants_timing %d [%s]", length(miss), paste(miss, collapse = " "))
fat_d <- intersect(fat_d, P$dispnum)
## fatality-code comparison: events.rds dispute fatality (stage-1 source chain) vs GML dispute fatality, single-source disputes
cmp <- merge(unique(EV[escalation == FALSE & n_disputes == 1L, .(dispnum, f_ev = fatality_disp)]), GD[, .(dispnum, f_g = fatality)], by = "dispnum")
cmp <- cmp[f_ev != f_g]
cmp[, kind := ifelse(f_ev < 0 | f_g < 0, "unknown code on one side", ifelse((f_ev > 0) != (f_g > 0), "fatal/nonfatal status differs", "severity level differs"))]
logf("fatality codes, events.rds vs GML (single-source disputes, %d differ): %s", nrow(cmp), paste(sprintf("%s %d", names(table(cmp$kind)), table(cmp$kind)), collapse = "; "))
fwrite(cmp, file.path(OUT, "fatality_code_differences.csv"))

## ---- 2. participant spells and exits (P1-1) --------------------------------------------------------------------------------------
PT <- P[dispnum %in% fat_d, .(dispnum, ccode, spell, spell_first, spell_last, side = side_entry, orig, entry_lo, entry_hi, side_switch)]
SE <- s2_spell_exits(PT, GP[, .(dispnum, ccode, styear, endyear, end_mk)])
PT <- merge(PT, SE, by = c("dispnum", "ccode", "spell"))
n_inv <- sum(PT$exit_mk < PT$entry_hi)
if (n_inv) { logf("spells whose exit precedes the definite entry month (exit set to entry month, flagged) %d", n_inv); PT[, exit_before_entry := exit_mk < entry_hi]; PT[exit_mk < entry_hi, exit_mk := entry_hi] } else PT[, exit_before_entry := FALSE]
logf("participant spells %d (%d participants): exit from GML %d, imputed (Dec of last spell year) %d, GML end capped at the spell's last year %d, >1 matching GML episode %d; side switches %d (first side kept)",
     nrow(PT), uniqueN(PT[, .(dispnum, ccode)]), sum(!PT$exit_imputed), sum(PT$exit_imputed), sum(PT$exit_capped), sum(PT$n_gml > 1L), sum(PT$side_switch))
for (r in which(PT$spell > 1L)) { k <- PT[dispnum == PT$dispnum[r] & ccode == PT$ccode[r]][order(spell)]
  logf("  second-spell audit: dispute %d, state %d: %s", PT$dispnum[r], PT$ccode[r], paste(sprintf("spell %d %d-%d side %d entry %d exit %d%s%s", k$spell, k$spell_first, k$spell_last, k$side, k$entry_hi, k$exit_mk,
       ifelse(k$exit_capped, " capped", ""), ifelse(k$exit_imputed, " imputed", "")), collapse = " | ")) }
gap_viol <- PT[order(dispnum, ccode, spell)][, .(viol = any(shift(exit_mk) >= s2_mk(spell_first, 1L), na.rm = TRUE)), by = .(dispnum, ccode)][viol == TRUE]
if (nrow(gap_viol)) fail("%d participants have a spell exit reaching into the next spell's first year", nrow(gap_viol))

## ---- 3. dispute table: reference onset, end (max spell exit), initiating side, fold ----------------------------------------------
ref <- P[orig == 1L, .(onset_lo = min(entry_lo), onset_hi = min(entry_hi)), by = dispnum]
dend <- s2_dispute_end(PT)                                        # recheck finding 1: endpoint provenance (exact / capped / imputed / adjusted)
DT <- merge(merge(data.table(dispnum = fat_d), ref, by = "dispnum", all.x = TRUE), dend, by = "dispnum", all.x = TRUE)
stopifnot(!anyNA(DT$onset_lo), !anyNA(DT$end_mk))
logf("dispute end = latest participant exit, by source: %s (exact = an unmodified GML exit attains it)", paste(sprintf("%s %d", names(table(DT$end_source)), table(DT$end_source)), collapse = ", "))
bad_end <- DT[end_mk < onset_lo]; if (nrow(bad_end)) { logf("availability exclusion: disputes ending before their reference onset %s", paste(bad_end$dispnum, collapse = " ")); DT <- DT[end_mk >= onset_lo] }
IS <- as.data.table(s2_init_side(P, EV, cow_war_inter, DT$dispnum))
DT <- merge(DT, IS, by = "dispnum")
logf("initiating side by rule: 1 %d, 2 %d, 3 %d, 4 %d, 5 (unresolved) %d; reviewed %d (hard %d)", sum(DT$side_rule == 1), sum(DT$side_rule == 2), sum(DT$side_rule == 3),
     sum(DT$side_rule == 4), sum(DT$side_rule == 5), sum(DT$side_case != ""), sum(DT$side_case == "hard"))
if (any(DT$side_rule == 5L)) { logf("availability exclusion: unresolved initiating side %s", paste(DT[side_rule == 5L, dispnum], collapse = " ")); DT <- DT[side_rule != 5L] }
PT <- PT[dispnum %in% DT$dispnum]
CLr <- s2_classify(PT, DT, W_REF)
ts <- CLr[baseline == TRUE, .(nA = sum(side_rel == "A"), nB = sum(side_rel == "B")), by = dispnum]
no_ts <- setdiff(DT$dispnum, ts[nA > 0 & nB > 0, dispnum])
if (length(no_ts)) { logf("availability exclusion: no baseline participant on one side (B3) %s", paste(no_ts, collapse = " ")); DT <- DT[!dispnum %in% no_ts]; PT <- PT[!dispnum %in% no_ts] }
DT[, onset_year := s2_mk_year(onset_lo)]; DT[, fold := s2_fold(onset_year)]
logf("disputes in the stage-2 sample: %d; by fold: %s; baseline sets with >1 initiating-side originator %d, >1 target-set member %d", nrow(DT),
     paste(sprintf("%d:%d", 1:6, tabulate(DT$fold, 6)), collapse = " "), sum(ts$nA > 1), sum(ts$nB > 1))

## ---- 4. classification (windows 0-6) and acceptance (R1, N2; P2-3, P2-6) ------------------------------------------------------------
CLW <- rbindlist(lapply(0:6, function(w) s2_classify(PT, DT, w)[, W := w]))
tab <- dcast(CLW[baseline == FALSE, .N, by = .(W, cat)], W ~ cat, value.var = "N", fill = 0); if (!"orig" %in% names(tab)) tab[, orig := 0L]
logf("non-baseline participants by window (A / B / orig): %s", paste(sprintf("W%d %d/%d/%d", tab$W, tab$A, tab$B, tab$orig), collapse = "; "))
CL <- CLW[W == W_REF][, W := NULL]
sm <- CL[baseline == FALSE & m_entry == 0L, .N, by = side_rel]
logf("same-month non-original entrants by side (R1 count): %s", paste(sprintf("%s %d", sm$side_rel, sm$N), collapse = ", "))
## (i) full-universe baseline equality: per dispute, the builder's baseline set = every participant with orig == 1 in participants_timing
bu <- merge(CL[baseline == TRUE, .(dispnum, ccode, b = TRUE)], unique(P[dispnum %in% DT$dispnum & orig == 1L, .(dispnum, ccode, o = TRUE)]), by = c("dispnum", "ccode"), all = TRUE)
if (anyNA(bu$b) || anyNA(bu$o)) fail("baseline sets differ from orig == 1 on the full participant universe (%d rows)", sum(is.na(bu$b) | is.na(bu$o))) else logf("acceptance (i): baseline sets equal orig == 1 on all %d baseline participants", nrow(bu))
## (ii) promotion boundaries: orig iff A side, W >= 1, m_entry <= W (definite bound); uncertain dates counted where the possible bound differs
PB <- merge(CLW[baseline == FALSE], PT[spell == 1L, .(dispnum, ccode, entry_lo)], by = c("dispnum", "ccode"))
PB <- merge(PB, DT[, .(dispnum, onset_hi)], by = "dispnum")
bnd <- PB[, sum((cat == "orig") != (side_rel == "A" & W >= 1L & m_entry <= W))]
if (bnd) fail("promotion boundary violated in %d participant-windows", bnd) else logf("acceptance (ii): promotion boundary holds in all %d non-baseline participant-windows", nrow(PB))
PB[, poss := side_rel == "A" & W >= 1L & pmax(entry_lo - onset_hi, 0L) <= W]
logf("uncertain dates: initiating-side entrants promoted under the possible bound but not the definite bound, by window: %s",
     paste(PB[, .(n = sum(poss & cat != "orig")), by = W][order(W), sprintf("W%d %d", W, n)], collapse = " "))
## (iii) reconciliation of every stage-1 projection record (initiator x dispute x event year x type)
if (!is.null(DUMP)) {
  d <- readRDS(DUMP); stopifnot(d$types == "bcjs", as.integer(d$coorig_months) == W_REF)
  pr1 <- unique(as.data.table(d$D)[chosen == 1L, .(c1, disp, yr, type, orig_joiner)])
  cov <- PT[pr1, on = .(dispnum = disp, ccode = c1, spell_first <= yr, spell_last >= yr), nomatch = NA, .(disp = i.disp, c1 = i.c1, yr = i.yr, type = i.type, spell = x.spell)]
  cov <- cov[, .(n_spell = sum(!is.na(spell)), spell = if (all(is.na(spell))) NA_integer_ else min(spell, na.rm = TRUE)), by = .(disp, c1, yr, type)]
  pr1 <- merge(pr1, cov, by = c("disp", "c1", "yr", "type"))
  pr1[, reason := ifelse(!disp %in% DT$dispnum, "outside selected disputes", ifelse(n_spell == 0L, "missing participant spell", ifelse(n_spell > 1L, "ambiguous spell", ifelse(spell > 1L, "later spell (re-entry)", "matched"))))]
  logf("acceptance (iii): stage-1 projection records %d: %s", nrow(pr1), paste(sprintf("%s %d", names(table(pr1$reason)), table(pr1$reason)), collapse = "; "))
  if (any(pr1$reason %in% c("missing participant spell", "ambiguous spell"))) fail("projection records with a missing or ambiguous covering spell: %d", sum(pr1$reason %in% c("missing participant spell", "ambiguous spell")))
  lr <- pr1[reason == "later spell (re-entry)"]; if (nrow(lr)) logf("  later-spell records (kept in rosters, excluded from the first-entry projection): %s", paste(sprintf("%d in %d (%d)", lr$c1, lr$disp, lr$yr), collapse = "; "))
  X <- merge(pr1[reason == "matched"], CL, by.x = c("disp", "c1"), by.y = c("dispnum", "ccode"), all.x = TRUE)
  if (anyNA(X$cat)) fail("matched projection records without a stage-2 classification: %d", sum(is.na(X$cat)))
  a_base <- X[, sum((!orig_joiner) != baseline)]; a_prom <- X[orig_joiner == TRUE, sum((type == 1L) != (cat == "orig"))]
  a_side <- X[orig_joiner == TRUE & type %in% 2:3, sum(ifelse(type == 2L, "A", "B") != side_rel)]
  logf("acceptance (iii): matched records %d; baseline-role mismatches %d; promotion mismatches %d of %d original-joiner records (stage-1 promoted %d); joiner-side mismatches %d of %d",
       nrow(X), a_base, a_prom, X[orig_joiner == TRUE, .N], X[orig_joiner == TRUE, sum(type == 1L)], a_side, X[orig_joiner == TRUE & type %in% 2:3, .N])
  nofit <- CL[baseline == FALSE & !paste(dispnum, ccode) %in% paste(pr1$disp, pr1$c1)]
  logf("N2: non-baseline participants without a stage-1 projection record, classified by the participant-level rule: %d (A %d, B %d, orig %d)", nrow(nofit), sum(nofit$cat == "A"), sum(nofit$cat == "B"), sum(nofit$cat == "orig"))
  if (a_base + a_prom + a_side > 0) { fwrite(X[((!orig_joiner) != baseline) | (orig_joiner & (type == 1L) != (cat == "orig")) | (orig_joiner & type %in% 2:3 & ifelse(type == 2L, "A", "B") != side_rel)],
                                             file.path(OUT, "acceptance_mismatches.csv")); fail("acceptance mismatches %d (acceptance_mismatches.csv)", a_base + a_prom + a_side) }
}

## ---- 5. cohort, horizons, one-shot and monthly rows (B1, B2, N1) ----------------------------------------------------------------------
CS <- as.data.table(cow_states)[, .(ccode = as.integer(ccode), st_mk = s2_mk(styear, stmonth), end_mk = s2_mk(endyear, endmonth))]
COH <- s2_cohort(DT, CL, CS)                                     # asserts one system spell per (dispute, state)
DYN <- s2_dynamic_entrants(COH, CL); OS <- s2_oneshot(COH, CL)
logf("cohort rows %d; complete horizon %d; truncated: administrative %d, system exit %d; dynamic entrants (excluded) %d; entries after the horizon (out) %d",
     nrow(COH), sum(COH$complete), sum(COH$trunc_admin), sum(COH$trunc_exit), nrow(DYN), sum(OS$late))
tabs <- function(y) paste(sprintf("%s %d", names(table(y)), table(y)), collapse = ", ")
logf("one-shot outcomes W%d: all rows %s; complete-horizon rows %s", W_REF, tabs(OS$y), tabs(OS[complete == TRUE]$y))
for (w in 0:6) { o <- s2_oneshot(COH, CLW[W == w][, W := NULL]); logf("  window %d: %s; retrospective cohort drops %d", w, tabs(o$y), sum(o$y == "orig")) }
MON <- s2_monthly(OS, W_REF)
logf("monthly rows W%d: %d (%d state-disputes); outcomes %s", W_REF, nrow(MON), nrow(OS), tabs(MON$y_m))
logf("dispute duration (months): median %.0f, 90%% %.0f, max %d", median(DT$end_mk - DT$onset_lo), quantile(DT$end_mk - DT$onset_lo, .9), max(DT$end_mk - DT$onset_lo))

## ---- 6. rosters (B3) ---------------------------------------------------------------------------------------------------------------
PS <- merge(PT[, .(dispnum, ccode, spell, exit_mk, entry_hi)], CL[, .(dispnum, ccode, side_rel, baseline)], by = c("dispnum", "ccode"))
PS <- merge(PS, DT[, .(dispnum, onset_lo)], by = "dispnum")
PS[, m_in := ifelse(baseline & spell == 1L, 0L, pmax(entry_hi - onset_lo, 0L))]; PS[, m_out := exit_mk - onset_lo]
MO <- unique(MON[, .(dispnum, m)])
ROS <- s2_rosters(PS, MO, "carry"); emp <- attr(s2_rosters(PS, MO, "end"), "empty")
logf("rosters: dispute-months %d; carried side-months %d in %d disputes (A %d, B %d); end policy would drop %d monthly rows", nrow(MO), nrow(emp), uniqueN(emp$dispnum),
     sum(emp$side_rel == "A"), sum(emp$side_rel == "B"), nrow(MON) - nrow(s2_truncate_end(MON, emp)))

## ---- 7. main-measure ties (section 4, N4; P2-4 coverage) ---------------------------------------------------------------------------
## main measure (config.sh MAIN_SPEC) via R/02_build_positions.R's W2contig code path (copied), checked against the positions in EVERY year
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
COVER <- sort(unique(DY$year))
dev <- vapply(COVER, function(y) { k <- POS$year == y; vs <- POS$ccode[k]; e <- EDGE[year == y]
  g <- graph_from_data_frame(e[, .(as.character(a), as.character(b), weight = w)], directed = FALSE, vertices = data.frame(name = as.character(vs)))
  max(abs(log(page_rank(g, directed = FALSE)$vector[as.character(vs)]) - POS$R[k])) }, 0)
logf("main-measure edges: %d positive dyad-years; covered years %d-%d (%d); log PageRank vs %s: max abs deviation over all years %.1e (worst year %d)",
     nrow(EDGE), min(COVER), max(COVER), length(COVER), POSCOL, max(dev), COVER[which.max(dev)])
if (max(dev) > 1e-6) stop("stage-2 edges do not reproduce the stage-1 main-measure positions")
R0 <- ROS[m == 0L]
row0 <- merge(OS[, .(dispnum, ccode)], DT[, .(dispnum, onset_year)], by = "dispnum")[, m := 0L]
T1 <- s2_ties(row0[, .(dispnum, ccode, m, year_tie = onset_year - 1L)], R0, EDGE, COVER)
T5 <- s2_ties(row0[, .(dispnum, ccode, m, year_tie = onset_year - 5L)], R0, EDGE, COVER)
setnames(T5, c("tie_sum_A", "tie_max_A", "tie_sum_B", "tie_max_B", "tie_nocov"), paste0(c("tie_sum_A", "tie_max_A", "tie_sum_B", "tie_max_B", "tie_nocov"), "_l5"))
logf("tie coverage: one-shot rows without a covered tie year, lag 1 %d, lag 5 %d (NA, flagged)", sum(T1$tie_nocov), sum(T5$tie_nocov_l5))
MON <- merge(MON, DT[, .(dispnum, onset_lo)], by = "dispnum")
MON[, year_tie := s2_mk_year(onset_lo + m) - 1L]
TM <- s2_ties(MON[, .(dispnum, ccode, m, year_tie)], ROS, EDGE, COVER)
MON <- merge(MON, TM[, .(dispnum, ccode, m, tie_sum_A, tie_max_A, tie_sum_B, tie_max_B, tie_nocov)], by = c("dispnum", "ccode", "m"))
logf("ties: one-shot rows with a positive tie to A %d / to B %d of %d; monthly rows without tie coverage %d", sum(T1$tie_max_A > 0, na.rm = TRUE), sum(T1$tie_max_B > 0, na.rm = TRUE), nrow(T1), sum(MON$tie_nocov))

## ---- 8. tiers (frozen at onset) ------------------------------------------------------------------------------------------------------
## tiers use contiguity in the onset year ("frozen at onset", memo section 4) and alliances at t - 1; side covariates below use t - 1 throughout
BASE <- CL[baseline == TRUE, .(dispnum, ccode)]
DYk <- rbind(DY[, .(k = ccode1, j = ccode2, year, conttype, capdist, ally, defoff = as.integer(def > 0 | off > 0))], DY[, .(k = ccode2, j = ccode1, year, conttype, capdist, ally, defoff = as.integer(def > 0 | off > 0))])
pr <- merge(merge(COH[, .(dispnum, ccode)], BASE[, .(dispnum, j = ccode)], by = "dispnum", allow.cartesian = TRUE), DT[, .(dispnum, onset_year)], by = "dispnum")
pr <- merge(pr, DYk[, .(k, j, year, conttype)], by.x = c("ccode", "j", "onset_year"), by.y = c("k", "j", "year"), all.x = TRUE)
pr[, ty := onset_year - 1L]
pr <- merge(pr, DYk[, .(k, j, year, ally)], by.x = c("ccode", "j", "ty"), by.y = c("k", "j", "year"), all.x = TRUE)
MJ <- as.data.table(cow_majors)[, .(ccode, st_mk = s2_mk(styear, stmonth), end_mk = s2_mk(endyear, endmonth))]
maj_at <- function(X) unique(MJ[X, on = .(ccode, st_mk <= onset_lo, end_mk >= onset_lo), nomatch = 0L, .(dispnum = i.dispnum, ccode)])
maj <- maj_at(merge(COH[, .(dispnum, ccode)], DT[, .(dispnum, onset_lo)], by = "dispnum"))
TIER <- s2_tiers(COH, BASE, CONTIG = pr[!is.na(conttype) & conttype %in% 1:5, .(dispnum, ccode, j)], MAJ = maj, ALLY = pr[!is.na(ally) & ally > 0, .(dispnum, ccode, j)])
logf("relevance tiers: likely %d of %d (contiguous %d, same region %d, major %d, allied at t-1 %d); candidate-baseline pairs without a dyad row in the onset year %d",
     sum(TIER$likely), nrow(TIER), sum(TIER$any_contig), sum(TIER$any_reg), sum(TIER$major), sum(TIER$any_ally), sum(is.na(pr$conttype)))

## ---- 9. one-shot reference covariates (memo section 7; t - 1; P1-2, P2-4) ----------------------------------------------------------
## Strict aggregation: a side aggregate is NA when any member's input is missing (reason flags kept); never a partial sum or a false "no".
NMC <- as.data.table(cow_nmc)[, .(ccode, year, cinc)]
sx <- merge(row0[, .(dispnum, ccode, ty = onset_year - 1L)], R0[, .(dispnum, j = ccode, side_rel)], by = "dispnum", allow.cartesian = TRUE)[ccode != j]
sx <- merge(sx, NMC[, .(j = ccode, ty = year, cinc_j = cinc)], by = c("j", "ty"), all.x = TRUE)
sx <- merge(sx, DYk[, .(ccode = k, j, ty = year, capdist, conttype, ally, defoff)], by = c("ccode", "j", "ty"), all.x = TRUE)
## stage-1 directed controls mirrored toward each side (assemble_cell.R Z_DIRECTED, t - 1): democracy = polity2 >= 6 with missing -> 0 (the
## build's POLITY_DEM_CUT is assumed to be 6: verify against v2/R/01_build_data.R section 6); MID ongoing between k and a member in t - 1
DEMY <- as.data.table(ccode_democracy)[, .(ccode, year, dem = as.integer(!is.na(polity2) & polity2 >= 6))]
sx <- merge(sx, DEMY[, .(j = ccode, ty = year, dem_j = dem)], by = c("j", "ty"), all.x = TRUE); sx[is.na(dem_j), dem_j := 0L]
## conflict histories (GPT surface review P2-4): COW-MID histories mirror stage 1 (the main build keeps COW-based histories) and are the
## reference; GML histories are kept as a labelled alternative (*_gml)
CDY <- as.data.table(cow_mid_ddydisps); CD <- as.data.table(cow_mid_disps)
mid_tm1_of <- function(DDY) { g <- unique(DDY[, .(pdisp = dispnum, a = ccode1, b = ccode2, year)])
  merge(sx[, .(dispnum, ccode, j, ty)], g, by.x = c("ccode", "j", "ty"), by.y = c("a", "b", "year"))[pdisp != dispnum, .(v = 1L), by = .(dispnum, ccode, j)] }
mt <- mid_tm1_of(CDY); setnames(mt, "v", "mid_tm1"); mtg <- mid_tm1_of(GDY); setnames(mtg, "v", "mid_tm1_gml")
sx <- merge(merge(sx, mt, by = c("dispnum", "ccode", "j"), all.x = TRUE), mtg, by = c("dispnum", "ccode", "j"), all.x = TRUE)
sx[is.na(mid_tm1), mid_tm1 := 0L]; sx[is.na(mid_tm1_gml), mid_tm1_gml := 0L]
## prior MIDs (recheck finding 4): DISTINCT GML disputes between k and any member of side s starting no later than t - 1 (current excluded)
pm <- unique(merge(CDY[, .(pdisp = dispnum, a = ccode1, b = ccode2)], CD[, .(pdisp = dispnum, styear)], by = "pdisp"))
pmg <- unique(merge(GDY[, .(pdisp = dispnum, a = ccode1, b = ccode2)], GD[, .(pdisp = dispnum, styear)], by = "pdisp"))
PRI <- s2_prior_side(sx[, .(dispnum, ccode, j, side_rel, ty)], pm)
PRIg <- s2_prior_side(sx[, .(dispnum, ccode, j, side_rel, ty)], pmg)
logf("conflict histories: COW-MID (reference, mirrors stage 1) and GML (alternative); candidate-member pairs with a MID at t-1 COW %d / GML %d; distinct prior-dispute counts > 0 by side COW %d / GML %d",
     sum(sx$mid_tm1), sum(sx$mid_tm1_gml), sum(PRI$prior_n > 0), sum(PRIg$prior_n > 0))
## side major-power status at t - 1 (recheck finding 2): major in December of the year before onset; the candidate's tier stays onset-based
maj_s <- maj_at(merge(R0[, .(dispnum, ccode)], DT[, .(dispnum, onset_lo = s2_mk(onset_year - 1L, 12L))], by = "dispnum"))[, maj := TRUE]
sx <- merge(sx, maj_s[, .(dispnum, j = ccode, maj)], by = c("dispnum", "j"), all.x = TRUE); sx[is.na(maj), maj := FALSE]
## side size counts the whole roster (including members the candidate cannot tie to: none, since cohort members are never baselines)
SC <- sx[, .(n_side = .N, side_cinc = s2_strict_sum(cinc_j), cinc_miss = anyNA(cinc_j), min_capdist = s2_strict_min(capdist),
             contig_side = s2_strict_any(ifelse(is.na(conttype), NA, conttype %in% 1:5)), ally_side = s2_strict_any(ifelse(is.na(ally), NA, ally > 0)),
             dyad_miss = anyNA(capdist) | anyNA(conttype) | anyNA(ally), major_side = as.integer(any(maj)),
             land_side = s2_strict_any(ifelse(is.na(conttype), NA, conttype == 1L)), defoff_side = s2_strict_any(ifelse(is.na(defoff), NA, defoff > 0)),
             dem_side = mean(dem_j), mid_tm1_side = as.integer(any(mid_tm1 == 1L)), mid_tm1_side_gml = as.integer(any(mid_tm1_gml == 1L))), by = .(dispnum, ccode, side_rel)]
SC <- merge(SC, PRI[, .(dispnum, ccode, side_rel, prior_side = log1p(prior_n))], by = c("dispnum", "ccode", "side_rel"))
SC <- merge(SC, PRIg[, .(dispnum, ccode, side_rel, prior_side_gml = log1p(prior_n))], by = c("dispnum", "ccode", "side_rel"))
SC <- dcast(SC, dispnum + ccode ~ side_rel, value.var = c("n_side", "side_cinc", "cinc_miss", "min_capdist", "contig_side", "ally_side", "dyad_miss", "prior_side", "major_side", "land_side", "defoff_side", "dem_side", "mid_tm1_side", "mid_tm1_side_gml", "prior_side_gml"))
own <- merge(row0[, .(dispnum, ccode, ty = onset_year - 1L)], NMC[, .(ccode, ty = year, cinc_k = cinc)], by = c("ccode", "ty"), all.x = TRUE)
own <- merge(own, POS[, .(ccode, ty = year, own_cent = R)], by = c("ccode", "ty"), all.x = TRUE)
own <- merge(own, DEMY[, .(ccode, ty = year, dem_k = dem)], by = c("ccode", "ty"), all.x = TRUE); own[is.na(dem_k), dem_k := 0L]
majk <- maj_at(merge(row0[, .(dispnum, ccode)], DT[, .(dispnum, onset_lo = s2_mk(onset_year - 1L, 12L))], by = "dispnum"))[, major_k := 1L]
own <- merge(own, majk, by = c("dispnum", "ccode"), all.x = TRUE); own[is.na(major_k), major_k := 0L]
SC <- merge(SC, own[, .(dispnum, ccode, cinc_k, own_cent, dem_k, major_k)], by = c("dispnum", "ccode"))
SC[, `:=`(dem_x_A = dem_k * dem_side_A, dem_x_B = dem_k * dem_side_B, major_x_A = major_k * major_side_A, major_x_B = major_k * major_side_B)]
SC[, `:=`(log_cinc_k = ifelse(cinc_k > 0, log(cinc_k), NA_real_),
          log_cinc_ratio_A = ifelse(cinc_k > 0 & side_cinc_A > 0, log(cinc_k / side_cinc_A), NA_real_), log_cinc_ratio_B = ifelse(cinc_k > 0 & side_cinc_B > 0, log(cinc_k / side_cinc_B), NA_real_),
          share_A = ifelse(side_cinc_A + side_cinc_B > 0, side_cinc_A / (side_cinc_A + side_cinc_B), NA_real_))]
SC[, share_B := 1 - share_A]
tgt <- merge(CL[baseline == TRUE & side_rel == "B", .(dispnum, ccode)], DT[, .(dispnum, ty = onset_year - 1L)], by = "dispnum")
tgt <- merge(tgt, POS[, .(ccode, ty = year, R)], by = c("ccode", "ty"), all.x = TRUE)
TC <- tgt[, .(tgt_cent_mean = if (anyNA(R)) NA_real_ else mean(R), tgt_cent_min = if (anyNA(R)) NA_real_ else min(R), tgt_pos_miss = anyNA(R)), by = dispnum]
logf("covariate missingness (rows): own CINC %d, own centrality %d; side CINC incomplete A %d / B %d; side dyad rows incomplete A %d / B %d; target set without t-1 position %d disputes",
     sum(is.na(SC$log_cinc_k)), sum(is.na(SC$own_cent)), sum(SC$cinc_miss_A), sum(SC$cinc_miss_B), sum(SC$dyad_miss_A), sum(SC$dyad_miss_B), sum(TC$tgt_pos_miss))

## ---- 10. jointly fatal; assemble; reference schema and estimation cohort ------------------------------------------------------------
PF <- GDY[, .(fat = as.integer(if (any(fatality1 > 0)) max(fatality1) else if (all(fatality1 == 0)) 0 else -9)), by = .(dispnum, ccode = ccode1)]
JF <- s2_jfatal(CL[, .(dispnum, ccode, side_rel)], PF)
logf("jointly fatal (side level): yes %d, no %d, unresolved %d; participants without a fatality code %d", sum(JF$jfatal == "yes"), sum(JF$jfatal == "no"), sum(JF$jfatal == "unresolved"),
     nrow(CL[!PF, on = .(dispnum, ccode)]))
DT <- merge(merge(DT, JF, by = "dispnum"), TC, by = "dispnum", all.x = TRUE)
DT <- merge(DT, GD[, .(dispnum, fat_level = fatality)], by = "dispnum", all.x = TRUE)
ONE <- Reduce(function(a, b) merge(a, b, by = c("dispnum", "ccode"), all.x = TRUE),
              list(OS, T1[, .(dispnum, ccode, tie_sum_A, tie_max_A, tie_sum_B, tie_max_B, tie_nocov)], T5[, .(dispnum, ccode, tie_sum_A_l5, tie_max_A_l5, tie_sum_B_l5, tie_max_B_l5, tie_nocov_l5)],
                   TIER[, .(dispnum, ccode, likely, any_contig, any_reg, major, any_ally)], SC))
ONE <- merge(ONE, DT[, .(dispnum, fold, onset_year, jfatal, fat_level, tgt_cent_mean, tgt_cent_min, side_case, end_source, end_imputed, end_exact)], by = "dispnum")
ONE[, `:=`(lcap_A = log1p(min_capdist_A), lcap_B = log1p(min_capdist_B))]          # fixed transforms named in S2_REF_DESIGN
## stage-2 targeting surface (PI 2026-11-22): joiner's gaps to the opposing side, stage-1 scaling (dM pooled scale of log CINC share over all
## state-years of sy_positions; dR scale of the main-measure positions in the onset year); side value = mean of the month-0 roster ("min" kept
## as sweep columns)
SYP <- as.data.table(readRDS(file.path(DE, "sy_positions.rds")))[, .(ccode, year, M)]
PV <- merge(SYP, POS[, .(ccode, year, R)], by = c("ccode", "year"), all = TRUE)
sM <- s2_gap_scale_pooled(SYP$M, SYP$year); sR <- vapply(split(POS$R, POS$year), s2_gap_scale, 0)
rows_s <- row0[, .(dispnum, ccode, ty = onset_year - 1L, t = onset_year)]
mem_s <- R0[, .(dispnum, side_rel, j = ccode)]
SURF <- s2_surface(rows_s, mem_s, PV, sM, sR, "mean")
SURFm <- s2_surface(rows_s, mem_s, PV, sM, sR, "min")[, .(dispnum, ccode, dM_kA_min = dM_kA, dR_kA_min = dR_kA, dM_kB_min = dM_kB, dR_kB_min = dR_kB)]
ONE <- merge(merge(ONE, SURF, by = c("dispnum", "ccode")), SURFm, by = c("dispnum", "ccode"))
logf("surface: dM scale %.3f (pooled log CINC share), dR scale by onset year (median %.3f); rows with a missing gap to A %d / to B %d; gap to B (dM, dR) quartiles %s / %s",
     sM, median(sR), sum(is.na(ONE$dR_kA) | is.na(ONE$dM_kA)), sum(is.na(ONE$dR_kB) | is.na(ONE$dM_kB)),
     paste(round(quantile(ONE$dM_kB, c(.25, .5, .75), na.rm = TRUE), 2), collapse = " "), paste(round(quantile(ONE$dR_kB, c(.25, .5, .75), na.rm = TRUE), 2), collapse = " "))
stopifnot(nrow(ONE) == nrow(OS))
## opposing-side robustness variants (A1 aggregation x A2 members x A3 timing x A4 dR scale); gaps only (the fitter forms the surface terms)
PP <- s2_primary_pair(EV, CL)
logf("primary pairs (first originator-originator event, initiator on the initiating side): %d of %d disputes", nrow(PP), nrow(DT))
mem_p <- rbind(PP[, .(dispnum, side_rel = "A", j = prim_init)], PP[, .(dispnum, side_rel = "B", j = prim_tgt)])
sRp <- setNames(rep(s2_gap_scale_pooled(POS$R, POS$year), length(sR)), names(sR))
VAR <- CJ(agg = c("mean", "min", "max", "lead", "total"), members = c("roster", "primary"), lag = c(1L, 5L), scale = c("year", "pooled"), sorted = FALSE)
VARS <- rows_s[, .(dispnum, ccode)]
for (i in seq_len(nrow(VAR))) { v <- VAR[i]; key <- paste(v$agg, v$members, paste0("l", v$lag), v$scale, sep = "_")
  sv <- s2_surface(rows_s[, .(dispnum, ccode, ty = t - v$lag, t)], if (v$members == "roster") mem_s else mem_p, PV, sM, if (v$scale == "year") sR else sRp, v$agg)
  setnames(sv, c("dM_kA", "dR_kA", "dM_kB", "dR_kB"), paste0(c("dM_kA", "dR_kA", "dM_kB", "dR_kB"), "__", key))
  VARS <- merge(VARS, sv[, c("dispnum", "ccode", paste0(c("dM_kA", "dR_kA", "dM_kB", "dR_kB"), "__", key)), with = FALSE], by = c("dispnum", "ccode"), all.x = TRUE) }
ref_key <- "mean_roster_l1_year"
stopifnot(isTRUE(all.equal(VARS[order(dispnum, ccode)][[paste0("dR_kB__", ref_key)]], ONE[order(dispnum, ccode)]$dR_kB)))   # reference variant = the built surface
nv <- VARS[, lapply(.SD, function(x) sum(is.na(x))), .SDcols = grep("^dR_kB__", names(VARS), value = TRUE)]
logf("surface variants: %d (aggregation x members x lag x dR scale); missing gap to B by variant ranges %d-%d rows (primary-pair and lag-5 coverage)", nrow(VAR), min(unlist(nv)), max(unlist(nv)))
saveRDS(VARS, file.path(OUT, "s2_surface_variants.rds"))
## inputs the fitter needs to rebuild the surface under perturbed positions (robustness B3): positions, rosters, scales, confounder counts
CONF <- rbind(DY[, .(ccode = ccode1, year, def, contig_land)], DY[, .(ccode = ccode2, year, def, contig_land)])[, .(n_def = sum(def > 0), n_contig = sum(contig_land > 0)), by = .(ccode, year)]
saveRDS(list(PV = PV, rows = rows_s, mem_roster = mem_s, mem_primary = mem_p, sM = sM, sR = sR, sR_pooled = sRp, conf = CONF), file.path(OUT, "s2_surface_inputs.rds"))
stopifnot(nrow(ONE) == nrow(OS), !anyDuplicated(ONE[, .(dispnum, ccode)]))
## reference one-shot: availability columns (S2_AVAIL) define the estimation cohort; S2_REF_DESIGN is the frozen model (recheck finding 3)
TERMS <- unique(unlist(S2_REF_DESIGN)); stopifnot(all(S2_AVAIL %in% names(ONE)), all(TERMS %in% names(ONE)))
ONE <- s2_ref_cohort(ONE, S2_AVAIL); census <- attr(ONE, "census")
stopifnot(all(ONE$y %in% S2_CATS), !anyDuplicated(ONE[, .(dispnum, ccode)]), all(ONE$H >= 0), all(ONE[y != "out", m_event <= H]))
ex <- ONE[in_ref == FALSE, .N, by = .(excl_reason, y)]
logf("reference estimation cohort: %d of %d rows; outcomes %s", sum(ONE$in_ref), nrow(ONE), tabs(ONE[in_ref == TRUE]$y))
logf("  exclusions by FIRST failing reason (mutually exclusive partition; horizon first, then S2_AVAIL order):")
for (r in unique(ex$excl_reason)) logf("    %-28s %s", r, paste(ex[excl_reason == r, sprintf("%s %d", y, N)], collapse = ", "))
logf("  nonexclusive missingness (all rows / complete-horizon rows): %s", paste(census[n_missing_any_row > 0, sprintf("%s %d/%d", column, n_missing_any_row, n_missing_complete_rows)], collapse = "; "))
fwrite(census, file.path(OUT, "missingness_census.csv"))
logf("  excluded rows by fold: %s; reference rows by endpoint source: %s", paste(ONE[in_ref == FALSE, .N, by = fold][order(fold), sprintf("%d:%d", fold, N)], collapse = " "),
     paste(ONE[in_ref == TRUE, .N, by = end_source][, sprintf("%s %d", end_source, N)], collapse = ", "))
stopifnot(!anyNA(ONE[in_ref == TRUE, ..TERMS]), all(vapply(TERMS, function(cc) all(is.finite(as.numeric(ONE[in_ref == TRUE][[cc]]))), TRUE)))
MON <- merge(MON, TIER[, .(dispnum, ccode, likely)], by = c("dispnum", "ccode"))
MON <- merge(MON, DT[, .(dispnum, fold, end_source, end_imputed, onset_year)], by = "dispnum")
stopifnot(all(MON$fold == s2_fold(MON$onset_year)))                  # every monthly row carries its dispute's onset fold (recheck finding 6)
MON[, onset_year := NULL]
META <- list(built = format(Sys.time()), git = GIT, mode = MODE, fatal_def = FATAL_DEF, fatal_label = if (FATAL_DEF == "stage1") "primary disputes represented by stage-1 fatal-event flags" else "GML dispute fatality > 0",
             window = W_REF, positions = POSF, poscol = POSCOL, avail_cols = S2_AVAIL, ref_design = S2_REF_DESIGN, ref_design_linear = S2_REF_DESIGN_LINEAR, monthly_covariates = "ties and frozen tier only (monthly non-tie covariates not built)")
for (o in c("DT", "ONE", "MON")) { x <- get(o); setattr(x, "s2_meta", META); assign(o, x) }
saveRDS(META, file.path(OUT, "s2_meta.rds")); saveRDS(DT, file.path(OUT, "s2_disputes.rds")); saveRDS(PT, file.path(OUT, "s2_participant_spells.rds"))
saveRDS(CLW, file.path(OUT, "s2_classification_w0_6.rds")); saveRDS(ONE, file.path(OUT, sprintf("s2_oneshot_w%d.rds", W_REF))); saveRDS(MON, file.path(OUT, sprintf("s2_monthly_w%d.rds", W_REF)))
saveRDS(ROS, file.path(OUT, "s2_rosters_carry.rds")); saveRDS(EDGE, file.path(OUT, "s2_main_edges.rds"))
logf("written (%s build): s2_disputes %d, s2_oneshot_w%d %d (reference cohort %d), s2_monthly_w%d %d, rosters %d, edges %d", MODE, nrow(DT), W_REF, nrow(ONE), sum(ONE$in_ref), W_REF, nrow(MON), nrow(ROS), nrow(EDGE))
logf("S2_BUILD_DONE")
