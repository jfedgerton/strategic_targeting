## R/lib/s2.R -- stage-2 (joining) risk-set library. Contract: docs/stage2_design_memo.md section 13 (revision 6; N1-N4 > R1-R3 > B1-C2).
## Pure functions on data.frames/data.tables so every rule is testable on fixtures (tests/test_s2.R). No model is fit here.
##
## Conventions
##   month key  mk = 12 * year + month - 1 (the data build's key; May 1902 = 22828)
##   reference  onset_lo / onset_hi = the dispute's earliest ORIGINATOR entry bounds (stage-1 TIMING_REF = "originator", R05 section 1b)
##   month idx  m = mk - onset_lo; an entry is dated by its definite (latest) bound: m_entry = max(entry_hi - onset_lo, 0) (memo section 2)
##   sides      COW side codes 0/1 per spell (initial side, side switches keep the first side); "A" = the dispute's initiating side (stage-1
##              case-reviewed map, s2_init_side), "B" = the other side
##   baseline   participants with orig == 1 (R1); they enter at month 0 by definition (as in stage 1, R05: "originators enter at onset")
suppressPackageStartupMessages(library(data.table))

S2_MK_END <- 12L * 2014L + 11L                                  # Dec 2014, administrative end of follow-up
S2_FOLDS <- data.frame(fold = 1:6, lo = c(1816L, 1914L, 1919L, 1939L, 1946L, 1990L), hi = c(1913L, 1918L, 1938L, 1945L, 1989L, 2014L))
S2_CATS <- c("A", "B", "orig", "out")

s2_mk <- function(y, m) {                                       # GPT build review P2-5: invalid month codes never become valid keys
  y <- as.integer(y); m <- as.integer(m)
  if (anyNA(y) || anyNA(m) || any(m < 1L | m > 12L) || any(y < 1800L | y > 2100L)) stop(sprintf("s2_mk: invalid year/month codes (%d rows)", sum(is.na(y) | is.na(m) | !(m %in% 1:12) | !(y %in% 1800:2100))))
  12L * y + m - 1L
}
s2_mk_year <- function(mk) as.integer(mk %/% 12L)

## ---- folds (B4): six disjoint era folds by the dispute's onset year --------------------------------------------------------------------
s2_fold <- function(year) {
  f <- findInterval(year, S2_FOLDS$lo)
  stopifnot(!anyNA(year), all(year >= 1816L & year <= 2014L), all(year <= S2_FOLDS$hi[f]))
  as.integer(f)
}

## ---- availability masks (R1, N1). Prospective: a function of the month index and the window ONLY -------------------------------------
## W = 0: no window promotion in stage 1 -> no originator alternative; A available in every month.
## W >= 1: months 0..W an initiating-side entry is always an originator entry (A masked, orig available); from W + 1, A available.
s2_mask <- function(m, W) {
  stopifnot(length(W) == 1L, W %in% 0:6, all(m >= 0))
  win <- W >= 1L & m <= W
  data.frame(avail_A = !win, avail_B = rep(TRUE, length(m)), avail_orig = win)
}
## auxiliary three-category q model / sequential first equation (R2): third-party entry T, originator entry O, out; same mask
s2_mask_aux <- function(m, W) { k <- s2_mask(m, W); data.frame(avail_T = k$avail_A | k$avail_B, avail_O = k$avail_orig) }

## ---- participant-level classification (R1, N2) -------------------------------------------------------------------------------------
## PT: one row per (dispnum, ccode, spell) with side (0/1), orig, entry_lo, entry_hi; DT: dispnum, onset_lo, init_side.
## Returns one row per (dispnum, ccode): first spell, side_rel A/B, baseline flag, m_entry, category at window W
## (baseline / A / B / orig). Promotion = initiating side, not baseline, W >= 1, definite bound m_entry <= W (R05 recode(), definite),
## applied to every participant whether or not stage 1 fitted an event for it.
s2_classify <- function(PT, DT, W) {
  stopifnot(W %in% 0:6, !anyDuplicated(DT$dispnum), all(PT$dispnum %in% DT$dispnum))
  P <- as.data.table(PT)[order(dispnum, ccode, spell)]
  P1 <- P[!duplicated(P[, .(dispnum, ccode)])]
  P1 <- merge(P1, as.data.table(DT)[, .(dispnum, onset_lo, init_side)], by = "dispnum", sort = FALSE)
  stopifnot(!anyNA(P1$init_side), all(P1$init_side %in% 0:1), all(P1$side %in% 0:1))
  P1[, side_rel := ifelse(side == init_side, "A", "B")]
  P1[, baseline := orig == 1L]
  P1[, m_entry := ifelse(baseline, 0L, pmax(entry_hi - onset_lo, 0L))]
  P1[, cat := ifelse(baseline, "baseline", ifelse(side_rel == "B", "B", ifelse(W >= 1L & m_entry <= W, "orig", "A")))]
  P1[order(dispnum, ccode), .(dispnum, ccode, side, side_rel, baseline, m_entry, cat)]
}

## ---- initiating side map (stage-1 R05 rules 1-2 and the case review), for every dispute ------------------------------------------------
## P: participants_timing (all spells); EV: events.rds; WI: peacesciencer::cow_war_inter; disputes: dispnums to map.
## rule 1 = COW side of the initiator in originator-originator non-escalation events (unique side); rule 2 = not war-linked -> Side A (1);
## rule 3 = war rule (MID side of the war's initiating states); rule 4 = war mapping ambiguous -> Side A; rule 5 = war link only from
## mixed-source rows (R05 stops on these for fitted events; here they are returned as unresolved, listed and dropped by the builder).
S2_SIDE_REVIEW <- c("112" = 1L, "362" = 1L, "1035" = 0L, "1293" = 1L, "1363" = 1L, "3098" = 1L, "200" = 0L, "611" = 1L, "1046" = 1L,
                    "3826" = 0L, "1379" = 1L, "2150" = 1L, "3131" = 1L)
S2_SIDE_HARD <- c(200L, 1379L)
s2_init_side <- function(P, EV, WI, disputes) {
  P <- as.data.frame(P); EV <- as.data.frame(EV); WI <- as.data.frame(WI)
  po <- unique(P[, c("dispnum", "ccode", "orig")])
  stopifnot(!anyDuplicated(po[, c("dispnum", "ccode")]))
  oe <- EV[!EV$escalation, c("dispnum", "initiator", "target", "war_link", "n_disputes")]
  ps0 <- P[order(P$dispnum, P$ccode, P$spell_first), c("dispnum", "ccode", "side_entry")]; ps0 <- ps0[!duplicated(ps0[, c("dispnum", "ccode")]), ]
  oe$side_init <- ps0$side_entry[match(paste(oe$dispnum, oe$initiator), paste(ps0$dispnum, ps0$ccode))]
  oe$pi <- po$orig[match(paste(oe$dispnum, oe$initiator), paste(po$dispnum, po$ccode))]
  oe$pt <- po$orig[match(paste(oe$dispnum, oe$target), paste(po$dispnum, po$ccode))]
  o1 <- oe[!is.na(oe$pi) & !is.na(oe$pt) & oe$pi == 1L & oe$pt == 1L & !is.na(oe$side_init), ]
  isd1 <- tapply(o1$side_init, o1$dispnum, function(s) if (length(unique(s)) == 1L) as.integer(s[1]) else NA_integer_)
  war_init_states <- function(w, side_members = FALSE) { x <- WI[WI$warnum == w, ]; side <- unique(c(x$sidea1[x$initiator1 == 1], x$sidea2[x$initiator2 == 1])); stopifnot(length(side) == 1L)
    if (!side_members) unique(c(x$ccode1[x$initiator1 == 1], x$ccode2[x$initiator2 == 1])) else unique(c(x$ccode1[x$sidea1 == side], x$ccode2[x$sidea2 == side])) }
  side_dir <- function(d) {
    rows <- oe$dispnum == d; single <- any(rows & oe$n_disputes == 1L); if (single) rows <- rows & oe$n_disputes == 1L
    wl <- unique(na.omit(oe$war_link[rows]))
    if (!single && length(wl)) return(c(side = NA_integer_, rule = 5L))
    if (!length(wl)) return(c(side = 1L, rule = 2L))
    st <- unique(unlist(lapply(wl, war_init_states)))
    if (!any(P$ccode[P$dispnum == d] %in% st)) st <- unique(unlist(lapply(wl, war_init_states, side_members = TRUE)))
    s <- unique(P$side_entry[P$dispnum == d & P$ccode %in% st])
    if (length(s) == 1L) c(side = as.integer(s), rule = 3L) else c(side = 1L, rule = 4L)
  }
  out <- data.frame(dispnum = as.integer(disputes), init_side = NA_integer_, side_rule = 1L)
  r1 <- isd1[as.character(out$dispnum)]; out$init_side <- as.integer(r1)
  for (i in which(is.na(out$init_side))) { r <- side_dir(out$dispnum[i]); out$init_side[i] <- r[["side"]]; out$side_rule[i] <- r[["rule"]] }
  hit <- intersect(names(S2_SIDE_REVIEW), as.character(out$dispnum))
  got <- out$init_side[match(as.integer(hit), out$dispnum)]
  if (any(is.na(got) | got != S2_SIDE_REVIEW[hit])) stop("s2_init_side disagrees with the case review for disputes: ", paste(hit[is.na(got) | got != S2_SIDE_REVIEW[hit]], collapse = " "))
  out$side_case <- ifelse(out$dispnum %in% S2_SIDE_HARD, "hard", ifelse(as.character(out$dispnum) %in% names(S2_SIDE_REVIEW), "reviewed", ""))
  out
}

## ---- cohort and horizons (B1, B2) ----------------------------------------------------------------------------------------------------
## SYS: system spells (ccode, st_mk, end_mk). Cohort of dispute d = members at onset_lo (st_mk <= onset_lo <= end_mk) that are not baseline
## participants of d. Horizon end (month index) = min(dispute end, Dec 2014, system exit). complete = dispute ended by Dec 2014 and the
## state did not exit the system before the dispute ended.
s2_cohort <- function(DT, CL, SYS) {
  DT <- as.data.table(DT); SYS <- as.data.table(SYS); CL <- as.data.table(CL)
  stopifnot(all(c("dispnum", "onset_lo", "end_mk") %in% names(DT)), all(DT$end_mk >= DT$onset_lo))
  C <- SYS[DT, on = .(st_mk <= onset_lo, end_mk >= onset_lo), nomatch = 0L, allow.cartesian = TRUE,
           .(dispnum = i.dispnum, ccode = x.ccode, sys_end = x.end_mk, onset_lo = i.onset_lo, d_end = i.end_mk)]
  stopifnot(!anyDuplicated(C[, .(dispnum, ccode)]))
  base <- CL[baseline == TRUE, .(dispnum, ccode)]
  C <- C[!base, on = .(dispnum, ccode)]
  C[, trunc_admin := d_end > S2_MK_END]
  C[, trunc_exit := sys_end < pmin(d_end, S2_MK_END)]
  C[, H := pmin(d_end, S2_MK_END, sys_end) - onset_lo]
  C[, complete := !trunc_admin & !trunc_exit]
  stopifnot(all(C$H >= 0))
  C[order(dispnum, ccode), .(dispnum, ccode, H, complete, trunc_admin, trunc_exit)]
}

## ---- one-shot labels (B1, B2; prospective cohort) -------------------------------------------------------------------------------------
## First entry (first spell) of a cohort member within its horizon: A / B / orig; else "out". Non-baseline participants outside the cohort
## (not system members at onset = dynamic entrants) are returned separately by s2_dynamic_entrants().
s2_oneshot <- function(COH, CL, cohort = "prospective") {
  stopifnot(cohort %in% c("prospective", "retrospective"))
  X <- merge(as.data.table(COH), as.data.table(CL)[baseline == FALSE, .(dispnum, ccode, m_entry, cat)], by = c("dispnum", "ccode"), all.x = TRUE)
  X[, late := !is.na(m_entry) & m_entry > H]                       # entry after the horizon (system exit or Dec 2014) -> not observed
  X[, y := ifelse(is.na(cat) | late, "out", cat)]
  X[, m_event := ifelse(y == "out", H, m_entry)]
  if (cohort == "retrospective") X <- X[y != "orig"]                # conditional on not becoming an originator (B1 alternative)
  X[order(dispnum, ccode), .(dispnum, ccode, y, m_event, H, complete, trunc_admin, trunc_exit, late)]
}
s2_dynamic_entrants <- function(COH, CL) {
  X <- as.data.table(CL)[baseline == FALSE]
  X[!as.data.table(COH)[, .(dispnum, ccode)], on = .(dispnum, ccode)]
}

## ---- monthly rows (B1, N1) ---------------------------------------------------------------------------------------------------------
## One row per cohort member and month 0..m_event; outcome in the event month, "out" before; masks from s2_mask(); rows after a first
## entry are not at risk (absorbing). A last row with y = "out" at H is censoring.
s2_monthly <- function(OS, W) {
  OS <- as.data.table(OS)
  M <- OS[rep(seq_len(.N), m_event + 1L)]
  M[, m := seq_len(.N) - 1L, by = .(dispnum, ccode)]
  M[, y_m := ifelse(m == m_event, y, "out")]
  mk <- s2_mask(M$m, W); M[, `:=`(avail_A = mk$avail_A, avail_B = mk$avail_B, avail_orig = mk$avail_orig)]
  if (any(M$y_m == "A" & !M$avail_A) || any(M$y_m == "orig" & !M$avail_orig)) stop("s2_monthly: an observed outcome is masked (classification and mask disagree)")
  M[, .(dispnum, ccode, m, y_m, avail_A, avail_B, avail_orig, H, complete)]
}

## ---- rosters (B3) ------------------------------------------------------------------------------------------------------------------
## Month 0: A = baseline initiating-side participants, B = original target set (baseline target side). Month m >= 1: participants whose
## spell covers month m - 1 (entry by the definite bound; baseline spells start at month 0), on the side of their FIRST spell. Empty side:
## "carry" keeps the last nonempty roster (flagged), "end" marks the month so follow-up can end (sweep level).
## PS: one row per (dispnum, ccode, spell) with m_in, m_out (month indices) and side_rel (first-spell side); months: data.frame(dispnum, m).
s2_rosters <- function(PS, months, policy = "carry") {
  stopifnot(policy %in% c("carry", "end"))
  PS <- as.data.table(PS); MO <- unique(as.data.table(months)[, .(dispnum, m)])[order(dispnum, m)]
  b0 <- PS[baseline == TRUE, .(dispnum, ccode, side_rel)]
  R0 <- unique(merge(MO[m == 0L], b0, by = "dispnum", allow.cartesian = TRUE))
  M1 <- MO[m >= 1L][, mm1 := m - 1L]
  R1 <- PS[M1, on = .(dispnum, m_in <= mm1, m_out >= mm1), nomatch = 0L, allow.cartesian = TRUE, .(dispnum, m = i.m, ccode = x.ccode, side_rel = x.side_rel)]
  R <- unique(rbind(R0[, .(dispnum, m, ccode, side_rel)], R1))
  ## per (dispute, month, side): fill empty sides
  grid <- MO[, .(side_rel = c("A", "B")), by = .(dispnum, m)]
  cnt <- R[, .(n = .N), by = .(dispnum, m, side_rel)]
  grid <- merge(grid, cnt, by = c("dispnum", "m", "side_rel"), all.x = TRUE); grid[is.na(n), n := 0L]
  setorder(grid, dispnum, side_rel, m)
  grid[, src_m := { s <- ifelse(n > 0L, m, NA_integer_); nafill(s, type = "locf") }, by = .(dispnum, side_rel)]
  if (any(is.na(grid$src_m))) stop("s2_rosters: a side is empty at month 0 (target set or initiating side unresolvable)")
  grid[, carried := n == 0L]
  if (policy == "end") { out <- R[, carried := FALSE]; attr(out, "empty") <- grid[carried == TRUE, .(dispnum, m, side_rel)]; return(out[]) }
  fill <- grid[carried == TRUE]
  RC <- R[fill, on = .(dispnum, side_rel, m = src_m), nomatch = 0L, allow.cartesian = TRUE, .(dispnum, m = i.m, ccode = x.ccode, side_rel)]
  out <- rbind(R[, carried := FALSE], RC[, carried := TRUE])
  stopifnot(!anyDuplicated(out[, .(dispnum, m, ccode, side_rel)]))
  out[order(dispnum, m, side_rel, ccode)]
}

## ---- ties (section 4, N4): raw edge weights on the measure's [0,1] scale --------------------------------------------------------------
## ROW: (dispnum, ccode, m, year_tie) candidate rows; ROS: rosters (dispnum, m, ccode, side_rel); E: (a, b, year, w), a < b, undirected.
## Returns tie_sum_A, tie_max_A, tie_sum_B, tie_max_B: 0 where no link in a COVERED year; NA (flag tie_nocov) where year_tie is not in
## `covered` (GPT build review P2-4: unavailable coverage is never an observed zero). A candidate never ties to itself.
s2_ties <- function(ROW, ROS, E, covered) {
  E <- as.data.table(E); stopifnot(all(E$a < E$b), !anyDuplicated(E[, .(a, b, year)]), all(E$w >= 0 & E$w <= 1), all(E$year %in% covered))
  EE <- rbind(E[, .(k = a, j = b, year, w)], E[, .(k = b, j = a, year, w)])
  X <- as.data.table(ROW)[as.data.table(ROS)[, .(dispnum, m, j = ccode, side_rel)], on = .(dispnum, m), allow.cartesian = TRUE, nomatch = 0L]
  X <- X[ccode != j]
  X <- merge(X, EE, by.x = c("ccode", "j", "year_tie"), by.y = c("k", "j", "year"), all.x = TRUE); X[is.na(w), w := 0]
  A <- X[, .(tie_sum = sum(w), tie_max = max(w)), by = .(dispnum, ccode, m, side_rel)]
  A <- dcast(A, dispnum + ccode + m ~ side_rel, value.var = c("tie_sum", "tie_max"), fill = 0)
  out <- merge(as.data.table(ROW), A, by = c("dispnum", "ccode", "m"), all.x = TRUE)
  for (cc in c("tie_sum_A", "tie_max_A", "tie_sum_B", "tie_max_B")) { if (!cc %in% names(out)) out[, (cc) := 0]; set(out, which(is.na(out[[cc]])), cc, 0) }
  out[, tie_nocov := !(year_tie %in% covered)]
  for (cc in c("tie_sum_A", "tie_max_A", "tie_sum_B", "tie_max_B")) set(out, which(out$tie_nocov), cc, NA_real_)
  out[order(dispnum, ccode, m)]
}

## ---- relevance tiers (section 4), frozen at onset ----------------------------------------------------------------------------------
s2_region <- function(cc) as.integer(cut(cc, c(1, 199, 399, 626, 699, 999), labels = FALSE))   # 1 Americas, 2 Europe, 3 Africa, 4 Middle East, 5 Asia/Oceania
## COH: cohort (dispnum, ccode); BASE: baseline participants (dispnum, ccode); CONTIG: (dispnum, ccode, j) pairs contiguous (<= 400 mi) at onset;
## MAJ: (dispnum, ccode) majors at onset; ALLY: (dispnum, ccode, j) any ATOP pact at t - 1.
s2_tiers <- function(COH, BASE, CONTIG, MAJ, ALLY) {
  COH <- as.data.table(COH)[, .(dispnum, ccode)]; BASE <- as.data.table(BASE)[, .(dispnum, j = ccode)]
  pr <- merge(COH, BASE, by = "dispnum", allow.cartesian = TRUE)[ccode != j]
  pr[, same_reg := s2_region(ccode) == s2_region(j)]
  ct <- unique(as.data.table(CONTIG)[, .(dispnum, ccode, j, contig = TRUE)]); al <- unique(as.data.table(ALLY)[, .(dispnum, ccode, j, ally = TRUE)])
  pr <- merge(merge(pr, ct, by = c("dispnum", "ccode", "j"), all.x = TRUE), al, by = c("dispnum", "ccode", "j"), all.x = TRUE)
  pr[is.na(contig), contig := FALSE]; pr[is.na(ally), ally := FALSE]
  T <- pr[, .(any_contig = any(contig), any_reg = any(same_reg), any_ally = any(ally)), by = .(dispnum, ccode)]
  T <- merge(COH, T, by = c("dispnum", "ccode"), all.x = TRUE)
  T[, major := paste(dispnum, ccode) %in% paste(MAJ$dispnum, MAJ$ccode)]
  for (cc in c("any_contig", "any_reg", "any_ally")) set(T, which(is.na(T[[cc]])), cc, FALSE)
  T[, likely := any_contig | any_reg | major | any_ally]
  T[order(dispnum, ccode)]
}

## ---- side-level jointly fatal (B5, R3 logical bounds) --------------------------------------------------------------------------------
## CL: classified participants (dispnum, ccode, side_rel); PF: participant fatality codes (dispnum, ccode, fat; -9 or NA = unknown).
## yes = each side has a known fatal participant; no = some side's participants are all known nonfatal; unresolved otherwise.
s2_jfatal <- function(CL, PF) {
  X <- merge(as.data.table(CL)[, .(dispnum, ccode, side_rel)], as.data.table(PF)[, .(dispnum, ccode, fat)], by = c("dispnum", "ccode"), all.x = TRUE)
  X[is.na(fat), fat := -9L]
  S <- X[, .(any_fatal = any(fat > 0), all_nonfatal = all(fat == 0)), by = .(dispnum, side_rel)]
  S <- dcast(S, dispnum ~ side_rel, value.var = c("any_fatal", "all_nonfatal"))
  stopifnot(!anyNA(S))
  S[, jfatal := ifelse(any_fatal_A & any_fatal_B, "yes", ifelse(all_nonfatal_A | all_nonfatal_B, "no", "unresolved"))]
  S[order(dispnum), .(dispnum, jfatal)]
}

## ---- incidence (B2, R2) ----------------------------------------------------------------------------------------------------------
## Flat model: p is a months x 4 matrix of category probabilities (A, B, orig, out) for one row, months 0..H, masked entries already 0.
## F_s(H) = sum_m S(m - 1) p_s(m), S(m) = prod_{r <= m} p_out(r).
s2_cif <- function(p) {
  stopifnot(is.matrix(p), ncol(p) == 4L, all(abs(rowSums(p) - 1) < 1e-10), all(p >= 0))
  S <- c(1, cumprod(p[, 4]))[seq_len(nrow(p))]
  setNames(colSums(S * p[, 1:3, drop = FALSE]), c("A", "B", "orig"))
}
## Sequential model: hT, hO hazards and pA = P(A | third-party entry) by month. F_A = sum_m S(m - 1) hT(m) pA(m), S = prod (1 - hT - hO).
s2_seq_cif <- function(hT, hO, pA) {
  stopifnot(length(hT) == length(hO), length(hT) == length(pA), all(hT + hO <= 1 + 1e-12), all(hT >= 0 & hO >= 0))
  S <- c(1, cumprod(1 - hT - hO))[seq_along(hT)]
  c(A = sum(S * hT * pA), B = sum(S * hT * (1 - pA)), orig = sum(S * hO))
}

## ---- legality table (C1, R3): machine-readable rules the task generator applies ------------------------------------------------------
S2_DIMS <- list(boundary = c("onset", "late"), estimator = c("oneshot", "monthly"), form = c("flat", "sequential"),
                cohort = c("prospective", "retrospective"), riskset = c("full_tier", "likely", "full_lambda1", "full_lambda0"),
                tier_timing = c("onset", "monthly"), tie_agg = c("sum", "max"), tie_timing = c("lag1", "onset"),
                target = c("mean", "min"), switch = c("first", "censor"), reentry = c("ignore", "new"), dynamic = c("exclude", "include"),
                alliance = c("on", "off"), lag = c("1", "5"), window = as.character(0:6), sample = c("fmid", "mid", "jfatal", paste0("lev", 2:6)),
                empty = c("carry", "end"), horizon = c("complete", "truncated"), agg2 = c("no", "yes"))
S2_REF <- vapply(S2_DIMS, `[`, "", 1)
## returns NA if illegal, else the duplicate-cell key (ineffective levels mapped to the reference)
s2_cell_key <- function(cell) {
  stopifnot(setequal(names(cell), names(S2_DIMS)), all(vapply(names(cell), function(d) cell[[d]] %in% S2_DIMS[[d]], TRUE)))
  c1 <- cell
  if (c1$estimator == "oneshot" && (c1$reentry == "new" || c1$tier_timing == "monthly" || c1$dynamic == "include" || c1$agg2 == "yes")) return(NA_character_)
  if (c1$estimator == "monthly" && c1$horizon == "truncated") return(NA_character_)
  if (c1$tier_timing == "monthly" && c1$riskset %in% c("full_lambda1", "full_lambda0")) return(NA_character_)   # updated tier x full set without tier terms: duplicate
  if (c1$window == "0" && c1$cohort == "retrospective") c1$cohort <- "prospective"           # W = 0: no originator category either way
  ## not mapped here (data-dependent): side-switch censoring is de-duplicated by the builder when no switcher affects a roster
  paste(vapply(names(S2_DIMS), function(d) c1[[d]], ""), collapse = "|")
}

## ---- spell exits (GPT build review P1-1) ----------------------------------------------------------------------------------------
## PT spells (dispnum, ccode, spell, spell_first, spell_last); G: GML participant episodes (dispnum, ccode, styear, endyear, end_mk).
## A GML episode overlapping a spell by year gives the exit month, CAPPED at December of the spell's last year so a GML episode spanning
## two separated stage-1 spells never bridges the inactive gap. No match -> December of spell_last (imputed). Flags: exit_imputed,
## exit_capped (GML end beyond the spell's last year), n_gml (matching episodes; > 1 counted).
s2_spell_exits <- function(PT, G) {
  PT <- as.data.table(PT); G <- as.data.table(G)
  mx <- G[PT, on = .(dispnum, ccode, styear <= spell_last, endyear >= spell_first), nomatch = NA,
          .(dispnum, ccode, spell = i.spell, spell_last = i.spell_last, end_mk = x.end_mk)]
  mx <- mx[, .(exit_gml = if (all(is.na(end_mk))) NA_integer_ else max(end_mk, na.rm = TRUE), n_gml = sum(!is.na(end_mk)), spell_last = spell_last[1]), by = .(dispnum, ccode, spell)]
  mx[, cap := s2_mk(spell_last, 12L)]
  mx[, exit_imputed := is.na(exit_gml)]
  mx[, exit_capped := !exit_imputed & exit_gml > cap]
  mx[, exit_mk := ifelse(exit_imputed, cap, pmin(exit_gml, cap))]
  stopifnot(nrow(mx) == nrow(PT), !anyNA(mx$exit_mk))
  mx[order(dispnum, ccode, spell), .(dispnum, ccode, spell, exit_mk, exit_imputed, exit_capped, n_gml)]
}

## ---- end-of-follow-up policy (B3 sweep level "end") -------------------------------------------------------------------------------
## MON monthly rows; EMPTY (dispnum, m, side_rel) empty side-months. Follow-up ends before the first empty month of either side.
s2_truncate_end <- function(MON, EMPTY) {
  MON <- as.data.table(MON); f <- as.data.table(EMPTY)[, .(m_stop = min(m)), by = dispnum]
  X <- merge(MON, f, by = "dispnum", all.x = TRUE)
  X <- X[is.na(m_stop) | m < m_stop]; X[, m_stop := NULL]
  X[]
}

## ---- conditional side equation (sequential model) respects the A mask (GPT build review P3-7) ---------------------------------------
## P(A | third-party entry, m) = exp(uA) / (exp(uA) + exp(uB)) where A is available, exactly 0 where avail_A is FALSE.
s2_side_prob <- function(uA, uB, avail_A) {
  stopifnot(length(uA) == length(uB), length(uA) == length(avail_A))
  p <- plogis(uA - uB); p[!avail_A] <- 0; p
}

## ---- recheck fixes (GPT build recheck 2026-11-22, docs/REVIEW_stage2_build_recheck_gpt.md) ------------------------------------------
## Endpoint provenance (finding 1): the dispute end is the latest participant exit; it is "exact" only if an UNMODIFIED GML exit attains it
## (not imputed, not capped, not raised to the entry month). end_source = exact | capped | imputed | adjusted (precedence imputed > capped >
## adjusted when no exact spell attains the maximum).
s2_dispute_end <- function(PT) {
  PT <- as.data.table(PT); stopifnot(all(c("exit_imputed", "exit_capped", "exit_before_entry") %in% names(PT)))
  PT[, exact := !exit_imputed & !exit_capped & !exit_before_entry]
  E <- PT[, { mx <- max(exit_mk); at <- exit_mk == mx
    .(end_mk = mx, end_source = if (any(exact[at])) "exact" else if (any(exit_imputed[at])) "imputed" else if (any(exit_capped[at])) "capped" else "adjusted") }, by = dispnum]
  E[, end_imputed := end_source == "imputed"]; E[, end_exact := end_source == "exact"]
  E[order(dispnum)]
}
## Prior conflict per side (finding 4): DISTINCT historical disputes between candidate k and any member of side s, starting no later than
## year ty (the year before onset), the current dispute excluded. PAIRS: (dispnum, ccode, j, side_rel, ty); PM: (pdisp, a, b, styear) with
## both orientations. Returns (dispnum, ccode, side_rel, prior_n).
s2_prior_side <- function(PAIRS, PM) {
  PAIRS <- as.data.table(PAIRS); PM <- as.data.table(PM)
  X <- merge(PAIRS, PM, by.x = c("ccode", "j"), by.y = c("a", "b"), allow.cartesian = TRUE)[pdisp != dispnum & styear <= ty]
  out <- X[, .(prior_n = uniqueN(pdisp)), by = .(dispnum, ccode, side_rel)]
  out <- merge(unique(PAIRS[, .(dispnum, ccode, side_rel)]), out, by = c("dispnum", "ccode", "side_rel"), all.x = TRUE)
  out[is.na(prior_n), prior_n := 0L]; out[order(dispnum, ccode, side_rel)]
}
## Strict side aggregation (memo B5 amendment): NA when any member's input is missing; never a partial sum or a false "no".
s2_strict_sum <- function(x) if (anyNA(x)) NA_real_ else sum(x)
s2_strict_min <- function(x) if (anyNA(x)) NA_real_ else min(x)
s2_strict_any <- function(x) if (anyNA(x)) NA_integer_ else as.integer(any(x))
## Reference estimation cohort (finding 5): incomplete horizon first, then the first missing availability column in list order (a
## mutually exclusive partition); the nonexclusive missingness census is returned separately.
s2_ref_cohort <- function(ONE, avail) {
  ONE <- as.data.table(ONE); stopifnot(all(avail %in% names(ONE)), "complete" %in% names(ONE))
  ONE[, excl_reason := ifelse(!complete, "incomplete horizon", "")]
  for (cc in avail) ONE[excl_reason == "" & !is.finite(as.numeric(get(cc))), excl_reason := paste("missing", cc)]
  ONE[, in_ref := excl_reason == ""]
  census <- data.table(column = avail, n_missing_any_row = vapply(avail, function(cc) sum(!is.finite(as.numeric(ONE[[cc]]))), 0L),
                       n_missing_complete_rows = vapply(avail, function(cc) sum(ONE$complete & !is.finite(as.numeric(ONE[[cc]]))), 0L))
  attr(ONE, "census") <- census; ONE
}
## Frozen reference design (finding 3): availability columns vs model terms. The fitter builds exactly these alternative-specific terms.
##   A (join initiating side, third party), B (join target side), O (originator entry, own intercept, A-side covariates), out = reference.
##   Every term has its own coefficient per alternative. "likely" (relevance tier, frozen at onset) enters as an alternative-specific
##   intercept shift (the "full set with tier interactions" cell). Ties are raw sums (N4); tie maxima are a sweep level, not reference terms.
S2_REF_DESIGN <- list(
  A = c("likely", "tie_sum_A", "tgt_cent_mean", "own_cent", "log_cinc_ratio_A", "share_A", "lcap_A", "contig_side_A", "ally_side_A", "prior_side_A",
        "n_side_A", "n_side_B", "major_side_A", "major_side_B"),
  B = c("likely", "tie_sum_B", "tgt_cent_mean", "own_cent", "log_cinc_ratio_B", "share_B", "lcap_B", "contig_side_B", "ally_side_B", "prior_side_B",
        "n_side_A", "n_side_B", "major_side_A", "major_side_B"),
  O = c("likely", "tie_sum_A", "tgt_cent_mean", "own_cent", "log_cinc_ratio_A", "share_A", "lcap_A", "contig_side_A", "ally_side_A", "prior_side_A",
        "n_side_A", "n_side_B", "major_side_A", "major_side_B"))
S2_AVAIL <- c("tie_sum_A", "tie_sum_B", "tgt_cent_mean", "own_cent", "log_cinc_ratio_A", "log_cinc_ratio_B", "share_A", "min_capdist_A", "min_capdist_B",
              "contig_side_A", "contig_side_B", "ally_side_A", "ally_side_B", "prior_side_A", "prior_side_B", "n_side_A", "n_side_B", "major_side_A", "major_side_B", "likely")

## ---- stage-2 targeting surface (PI 2026-11-22: "stage 2 should have the same surface") ------------------------------------------------
## Joining a side means fighting the other side, so each join alternative gets the stage-1 quadratic surface in the joiner's gaps to the
## OPPOSING side: join A and originator entry -> gaps to B (the original target set); join B -> gaps to A (the initiating side at onset).
## Gaps follow stage 1 exactly (R/05, R/07): dM = (M_k - M_side) / s_M with M = log CINC share and s_M the pooled gap scale over all
## state-years; dR = (R_k - R_side) / s_R(t) with R the relational position and s_R(t) the gap scale of positions in year t; positions and
## shares at t - 1, t = onset year. Side value = mean over the side's month-0 roster (default; "min" = the most isolated / weakest member,
## a sweep level). Strict: NA if any member or the candidate lacks a value.
## ROWS: (dispnum, ccode, ty, t); MEM: (dispnum, side_rel, j); PV: (ccode, year, M, R); sM scalar; sR named by year (character).
s2_gap_scale <- function(v) { n <- length(v); stopifnot(n > 1, !anyNA(v)); sqrt(2 * n * sum((v - mean(v))^2) / (n * (n - 1))) }
s2_gap_scale_pooled <- function(v, yr) { ss <- tapply(v, yr, function(z) sum((z - mean(z))^2)); n <- tapply(v, yr, length); sqrt(sum(2 * n * ss) / sum(n * (n - 1))) }
## agg (opposing-side robustness A1): "mean" (reference); "min" = weakest M and most isolated R (taken separately); "max" = strongest M and
## most central R; "lead" = the member with the largest CINC share (its M and R); "total" = the side's combined share and combined PageRank
## (M and R are logs, so total = log(sum(exp(.)))). MEM picks the members (A2: onset roster, or the primary initiator / target only); ROWS$ty
## picks the timing (A3); sR the dR scale (A4: by year, or one pooled value repeated for every year).
s2_surface <- function(ROWS, MEM, PV, sM, sR, agg = "mean") {
  stopifnot(agg %in% c("mean", "min", "max", "lead", "total"), is.finite(sM), sM > 0, all(sR > 0))
  ROWS <- as.data.table(ROWS); MEM <- as.data.table(MEM); PV <- as.data.table(PV)
  ms <- merge(merge(MEM, unique(ROWS[, .(dispnum, ty)]), by = "dispnum", allow.cartesian = TRUE), PV[, .(j = ccode, ty = year, Mj = M, Rj = R)], by = c("j", "ty"), all.x = TRUE)
  agg_f <- function(Mj, Rj) {
    if (anyNA(Mj) || anyNA(Rj)) return(list(M_side = NA_real_, R_side = NA_real_))
    switch(agg, mean = list(M_side = mean(Mj), R_side = mean(Rj)), min = list(M_side = min(Mj), R_side = min(Rj)), max = list(M_side = max(Mj), R_side = max(Rj)),
           lead = { i <- which.max(Mj); list(M_side = Mj[i], R_side = Rj[i]) },
           total = list(M_side = log(sum(exp(Mj))), R_side = log(sum(exp(Rj))))) }
  SV <- ms[, agg_f(Mj, Rj), by = .(dispnum, ty, side_rel)]
  SV <- dcast(SV, dispnum + ty ~ side_rel, value.var = c("M_side", "R_side"))
  X <- merge(merge(ROWS, PV[, .(ccode, ty = year, Mk = M, Rk = R)], by = c("ccode", "ty"), all.x = TRUE), SV, by = c("dispnum", "ty"), all.x = TRUE)
  s_r <- unname(sR[as.character(X$t)]); stopifnot(!anyNA(s_r))
  X[, `:=`(dM_kA = (Mk - M_side_A) / sM, dR_kA = (Rk - R_side_A) / s_r, dM_kB = (Mk - M_side_B) / sM, dR_kB = (Rk - R_side_B) / s_r)]
  ## per alternative: sA_* (join A: vs B), sB_* (join B: vs A); the originator alternative uses sA_*
  X[, `:=`(sA_dM = dM_kB, sA_dM2 = dM_kB^2, sA_dR = dR_kB, sA_dR2 = dR_kB^2, sA_dMR = dM_kB * dR_kB,
           sB_dM = dM_kA, sB_dM2 = dM_kA^2, sB_dR = dR_kA, sB_dR2 = dR_kA^2, sB_dMR = dM_kA * dR_kA)]
  X[order(dispnum, ccode), .(dispnum, ccode, dM_kA, dR_kA, dM_kB, dR_kB, sA_dM, sA_dM2, sA_dR, sA_dR2, sA_dMR, sB_dM, sB_dM2, sB_dR, sB_dR2, sB_dMR)]
}
## primary pair (A2): the dispute's first originator-originator event with the initiator on the initiating side (earliest ev_start_lo, then
## year, initiator, target). EV: events.rds; CL: classification (baseline, side_rel). Returns (dispnum, prim_init, prim_tgt); disputes without
## such an event are absent (counted by the builder).
s2_primary_pair <- function(EV, CL) {
  EV <- as.data.table(EV)[escalation == FALSE, .(dispnum, initiator, target, year, ev_start_lo)]
  B <- as.data.table(CL)[baseline == TRUE, .(dispnum, ccode, side_rel)]
  X <- merge(merge(EV, B[, .(dispnum, initiator = ccode, si = side_rel)], by = c("dispnum", "initiator")), B[, .(dispnum, target = ccode, st = side_rel)], by = c("dispnum", "target"))
  X <- X[si == "A" & st == "B"][order(dispnum, ev_start_lo, year, initiator, target)]
  X[!duplicated(dispnum), .(dispnum, prim_init = initiator, prim_tgt = target)]
}
S2_SURF_A <- c("sA_dM", "sA_dM2", "sA_dR", "sA_dR2", "sA_dMR"); S2_SURF_B <- c("sB_dM", "sB_dM2", "sB_dR", "sB_dR2", "sB_dMR")
## Surface reference design (PI 2026-11-22: "stage 2 should mirror stage 1"; supersedes the linear S2_REF_DESIGN as the reference, which
## stays a labelled sweep level "linear"). Each join alternative mirrors the stage-1 directed utility toward the side it FIGHTS:
##   surface (5 terms) + the stage-1 directed controls Z_DIRECTED toward the opposing side o (land contiguity, log distance, o's democracy,
##   democracy of k x o, defense/offense pact with o, major power on o, major k x o, MID with o at t - 1, log prior MIDs with o)
## plus the stage-2 recruitment terms toward the JOINED side s (tie sum, any ATOP pact) and the dispute-level balance and sizes.
## Side aggregation: geography = closest member; democracy = share of democratic members; binary = any member. The originator alternative
## mirrors join A. The linear position and relative-power terms (own and target-set centrality, CINC ratio) are nested in the surface.
S2_REF_DESIGN_LINEAR <- S2_REF_DESIGN
S2_Z_TOWARD <- function(o) c(paste0("land_side_", o), paste0("lcap_", o), paste0("dem_side_", o), paste0("dem_x_", o), paste0("defoff_side_", o),
                              paste0("major_side_", o), paste0("major_x_", o), paste0("mid_tm1_side_", o), paste0("prior_side_", o))
S2_REF_DESIGN <- list(
  A = c("likely", S2_SURF_A, S2_Z_TOWARD("B"), "tie_sum_A", "ally_side_A", "share_A", "n_side_A", "n_side_B"),
  B = c("likely", S2_SURF_B, S2_Z_TOWARD("A"), "tie_sum_B", "ally_side_B", "share_B", "n_side_A", "n_side_B"),
  O = c(S2_SURF_A, setdiff(S2_Z_TOWARD("B"), "dem_x_B"), "tie_sum_A", "ally_side_A", "share_A", "n_side_A", "n_side_B"))
## O restrictions (GPT surface review P2-3, declared): originator entry is structurally unavailable to unlikely candidates (fit_mnl mask), so
## "likely" is not an O term; dem_x_B is excluded from O because no originator entrant is a democracy facing a democratic target side (its MLE
## does not exist) -- a declared restricted specification, reported as such
S2_AVAIL <- unique(c(S2_AVAIL, "dM_kA", "dR_kA", "dM_kB", "dR_kB", "land_side_A", "land_side_B", "defoff_side_A", "defoff_side_B",
                     "dem_side_A", "dem_side_B", "mid_tm1_side_A", "mid_tm1_side_B"))
## sweep dimensions added by the surface amendment: opposing-side position (mean of the side / its most isolated, weakest member) and design
S2_DIMS$target <- NULL
S2_DIMS$opp_side <- c("mean", "min"); S2_DIMS$design <- c("surface", "linear")
S2_REF <- vapply(S2_DIMS, `[`, "", 1)
