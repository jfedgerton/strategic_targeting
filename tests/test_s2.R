## tests/test_s2.R -- fixtures for R/lib/s2.R, the stage-2 risk-set library (docs/stage2_design_memo.md section 13 fixture list).
## Run: Rscript tests/test_s2.R (from the repo root). Synthetic data only; seed 123.
source("R/lib/s2.R")
set.seed(123)
ok_all <- TRUE
check <- function(name, cond) { cat(sprintf("%-92s %s\n", name, if (isTRUE(cond)) "PASS" else "FAIL")); if (!isTRUE(cond)) ok_all <<- FALSE }
expect_error <- function(expr) inherits(tryCatch({ force(expr); NULL }, error = function(e) e), "error")

## ---- fixture world ----------------------------------------------------------------------------------------------------------------
## dispute 1: onset month key 100, initiating side = COW side 1, ends at 110
##   10 baseline A; 20 baseline B (target), exits at 105; 30 non-original B entering in the onset month (exits 104);
##   40 non-original A entering in the onset month; 50 non-original A entering at month 3; 60 never participates; 70 leaves the system at 104;
##   80 enters the system at 102 and joins side B at 103 (dynamic entrant)
## dispute 2: onset 200, ends 201, only the two originators (no observed promotion); 60 is in its cohort
## dispute 3: onset Nov 2014, ends after Dec 2014 (administrative truncation)
DT <- data.frame(dispnum = 1:3, onset_lo = c(100L, 200L, S2_MK_END - 1L), end_mk = c(110L, 201L, S2_MK_END + 5L), init_side = c(1L, 1L, 0L))
PT <- data.frame(dispnum = c(1, 1, 1, 1, 1, 1, 2, 2, 3, 3), ccode = c(10, 20, 30, 40, 50, 80, 11, 21, 12, 22), spell = 1L,
                 side = c(1, 0, 0, 1, 1, 0, 1, 0, 0, 1), orig = c(1, 1, 0, 0, 0, 0, 1, 1, 1, 1),
                 entry_lo = c(100, 100, 100, 100, 103, 103, 200, 200, S2_MK_END - 1L, S2_MK_END - 1L),
                 entry_hi = c(100, 100, 100, 100, 103, 103, 200, 200, S2_MK_END - 1L, S2_MK_END - 1L),
                 exit_mk = c(110, 105, 104, 110, 110, 110, 201, 201, S2_MK_END + 5L, S2_MK_END + 5L))
PT[] <- lapply(PT, as.integer)
SYS <- data.frame(ccode = c(10L, 20L, 30L, 40L, 50L, 60L, 70L, 80L, 11L, 21L, 12L, 22L), st_mk = 0L, end_mk = 99999L)
SYS$end_mk[SYS$ccode == 70L] <- 104L; SYS$st_mk[SYS$ccode == 80L] <- 102L

CL0 <- s2_classify(PT, DT, 0L); CL1 <- s2_classify(PT, DT, 1L); CL3 <- s2_classify(PT, DT, 3L)
cat_of <- function(CL, d, k) CL$cat[CL$dispnum == d & CL$ccode == k]

## 1. baseline roles from orig, not month equality (R1)
check("1 baseline = orig == 1 (10, 20 in dispute 1); same-month non-originators are not baseline", all(CL1[dispnum == 1 & baseline == TRUE, ccode] == c(10, 20)) && !any(CL1[dispnum == 1 & ccode %in% c(30, 40), baseline]))
## 2. same-month non-original B entrant retained as B
COH <- s2_cohort(DT, CL1, SYS)
OS1 <- s2_oneshot(COH, CL1)
check("2 same-month non-original B entrant (30): category B at every window, entry month 0", cat_of(CL0, 1, 30) == "B" && cat_of(CL3, 1, 30) == "B" && OS1[dispnum == 1 & ccode == 30, y == "B" & m_event == 0L])
## 3. one eligible same-month A entrant: third-party at W = 0, originator at W = 1
check("3 same-month A entrant (40): A at W = 0, orig at W = 1", cat_of(CL0, 1, 40) == "A" && cat_of(CL1, 1, 40) == "orig")
check("3 month-3 A entrant (50): A at W = 1, orig at W = 3 (definite bound <= W)", cat_of(CL1, 1, 50) == "A" && cat_of(CL3, 1, 50) == "orig")
## 4. originator probability exactly 0 after its window; A masked inside it (N1)
M1 <- s2_monthly(OS1, 1L)
r60 <- M1[dispnum == 1 & ccode == 60]
check("4 W = 1: avail_orig TRUE in months 0-1, FALSE from month 2; avail_A the complement; avail_B always", identical(r60$avail_orig, r60$m <= 1L) && identical(r60$avail_A, r60$m > 1L) && all(r60$avail_B))
check("4 W = 0: no originator alternative in any month, A available from month 0", { k <- s2_mask(0:10, 0L); !any(k$avail_orig) && all(k$avail_A) })
## 5. identical pre-decision eligibility, different eventual outcomes -> identical masks
r40 <- M1[dispnum == 1 & ccode == 40]; r50 <- M1[dispnum == 1 & ccode == 50]
check("5 candidates 40 (orig at m0), 50 (A at m3), 60 (never) share masks month by month", identical(r40[, .(avail_A, avail_B, avail_orig)], r60[m == 0, .(avail_A, avail_B, avail_orig)]) &&
        identical(r50[, .(avail_A, avail_B, avail_orig)], r60[m <= 3, .(avail_A, avail_B, avail_orig)]))
check("5 the mask function takes only the month index and window (no candidate or outcome argument)", identical(names(formals(s2_mask)), c("m", "W")))
## 6. dispute with no observed promotion still offers the originator option at W >= 1
r2 <- M1[dispnum == 2 & ccode == 60]
check("6 dispute 2 (no promotion observed): originator available to 60 in months 0-1 at W = 1", nrow(r2) == 2L && all(r2$avail_orig))
## 7. participant-level rule (N2): classification uses no stage-1 event input
check("7 s2_classify takes participants, disputes and W only (no fitted-event argument); 50 classified without an event", identical(names(formals(s2_classify)), c("PT", "DT", "W")) && length(cat_of(CL3, 1, 50)) == 1L)
## 8. same masks in the auxiliary q model and the sequential first equation
ka <- s2_mask_aux(0:8, 2L); km <- s2_mask(0:8, 2L)
check("8 auxiliary/sequential mask: T always available, O exactly where the flat model's orig is", all(ka$avail_T) && identical(ka$avail_O, km$avail_orig))
## 9. outcome rows: the entry month carries the outcome; earlier months are "out"; absorbing after entry
check("9 monthly rows of 50 at W = 1: months 0-3, outcome A at month 3 only", identical(r50$m, 0:3) && identical(r50$y_m, c("out", "out", "out", "A")))
check("9 an observed outcome in a masked cell stops the build", expect_error(s2_monthly(data.frame(dispnum = 9L, ccode = 1L, y = "A", m_event = 0L, H = 3L, complete = TRUE), 1L)))
## 10. system exit and administrative truncation (B2)
check("10 system exit: 70 has H = 4, trunc_exit, not complete; outcome out", COH[dispnum == 1 & ccode == 70, H == 4L & trunc_exit & !complete] && OS1[dispnum == 1 & ccode == 70, y == "out"])
check("10 administrative truncation: dispute 3 rows trunc_admin, H to Dec 2014 (1 month), not complete", all(COH[dispnum == 3, trunc_admin & !complete & H == 1L]))
check("10 complete-horizon one-shot reference excludes truncated rows; monthly keeps them censored", nrow(OS1[complete == TRUE & dispnum == 3]) == 0L && nrow(M1[dispnum == 3]) > 0L)
## 11. dynamic entrants are not in the reference cohort
check("11 80 (entered the system after onset) is not in dispute 1's cohort; listed as a dynamic entrant", !any(COH$dispnum == 1 & COH$ccode == 80) && any(s2_dynamic_entrants(COH, CL1)$ccode == 80))
check("11 baselines are not in their own dispute's cohort", !any(COH$dispnum == 1 & COH$ccode %in% c(10, 20)))
## 12. retrospective cohort drops eventual window originators
OSr <- s2_oneshot(COH, CL1, cohort = "retrospective")
check("12 retrospective cohort: 40 dropped at W = 1, no orig outcome left", !any(OSr$dispnum == 1 & OSr$ccode == 40) && !any(OSr$y == "orig"))

## ---- rosters (B3) ---------------------------------------------------------------------------------------------------------------
PS <- merge(as.data.frame(PT), as.data.frame(CL1[, .(dispnum, ccode, side_rel, baseline, m_entry)]), by = c("dispnum", "ccode"))
PS <- merge(PS, DT[, c("dispnum", "onset_lo")], by = "dispnum"); PS$m_in <- PS$m_entry; PS$m_out <- PS$exit_mk - PS$onset_lo
MO <- data.frame(dispnum = 1L, m = 0:10)
RC <- s2_rosters(PS, MO, "carry"); RE <- s2_rosters(PS, MO, "end")
ros <- function(R, m, s) sort(R$ccode[R$dispnum == 1 & R$m == m & R$side_rel == s])
check("13 month 0: A = {10}, B = {20} (baselines only; same-month entrants 30, 40 not yet in a roster)", identical(ros(RC, 0, "A"), 10L) && identical(ros(RC, 0, "B"), 20L))
check("13 month 1: A = {10, 40}, B = {20, 30} (active in month 0)", identical(ros(RC, 1, "A"), c(10L, 40L)) && identical(ros(RC, 1, "B"), c(20L, 30L)))
check("13 50 (enters month 3) is not in A's roster before month 4", !(50L %in% ros(RC, 3, "A")) && 50L %in% ros(RC, 4, "A"))
check("13 80 (side B from month 3) joins B's roster at month 4", identical(ros(RC, 4, "B"), c(20L, 30L, 80L)))
## B: 30 out after m 4, 20 out after m 5, 80 to m 10 -> B never empties here; build an emptying case
PSe <- PS[!(PS$dispnum == 1 & PS$ccode == 80), ]
RCe <- s2_rosters(PSe, MO, "carry"); REe <- s2_rosters(PSe, MO, "end")
check("14 empty side, carry: B at month 7 carries month 6's roster {20}, flagged", identical(sort(RCe$ccode[RCe$m == 7 & RCe$side_rel == "B"]), 20L) && all(RCe$carried[RCe$m == 7 & RCe$side_rel == "B"]))
check("14 empty side, end: no B roster at month 7; empty months listed from 7", nrow(REe[m == 7 & side_rel == "B"]) == 0L && min(attr(REe, "empty")[side_rel == "B", m]) == 7L)
check("14 empty side at month 0 (no target set) stops", expect_error(s2_rosters(PSe[!(PSe$ccode == 20), ], MO, "carry")))

## ---- ties (section 4) ---------------------------------------------------------------------------------------------------------------
E <- data.frame(a = c(10L, 20L, 40L, 10L), b = c(60L, 60L, 60L, 40L), year = 8L, w = c(0.5, 0.2, 0.9, 0.7))
ROW <- data.frame(dispnum = 1L, ccode = c(60L, 40L), m = 1L, year_tie = 8L)
TI <- s2_ties(ROW, RC, E, covered = 8L)
check("15 ties of 60 at month 1: A = {10, 40} -> sum 1.4, max 0.9; B = {20, 30} -> sum 0.2, max 0.2 (no link = 0)",
      TI[ccode == 60, isTRUE(all.equal(c(tie_sum_A, tie_max_A, tie_sum_B, tie_max_B), c(1.4, 0.9, 0.2, 0.2)))])
check("15 roster member 40 in A: only the 10-40 link counts (edges are a < b, so self-ties cannot arise; the filter is defensive)", TI[ccode == 40, isTRUE(all.equal(c(tie_sum_A, tie_max_A), c(0.7, 0.7)))])
check("15 covered year without links gives observed zeros (not NA)", { t2 <- s2_ties(transform(ROW, year_tie = 7L), RC, E, covered = 7:8); all(t2$tie_sum_A == 0) && !any(t2$tie_nocov) })
check("15 uncovered year gives NA ties and tie_nocov (never an observed zero; GPT P2-4)", { t3 <- s2_ties(transform(ROW, year_tie = 6L), RC, E, covered = 7:8); all(is.na(t3$tie_sum_A)) && all(t3$tie_nocov) })

## ---- folds (B4) ---------------------------------------------------------------------------------------------------------------------
check("16 folds exhaustive and disjoint over 1816-2014; 1914 dispute spanning to 1919 takes fold 2 (onset)", { f <- s2_fold(1816:2014); !anyNA(f) && identical(as.integer(table(f)), S2_FOLDS$hi - S2_FOLDS$lo + 1L) && s2_fold(1914L) == 2L })
check("16 year outside 1816-2014 stops", expect_error(s2_fold(2015L)))

## ---- tiers (section 4) --------------------------------------------------------------------------------------------------------------
COHt <- data.frame(dispnum = 1L, ccode = c(60L, 260L, 710L, 900L, 365L))
BASEt <- data.frame(dispnum = 1L, ccode = c(10L, 20L))      # region 1 (Americas)
T <- s2_tiers(COHt, BASEt, CONTIG = data.frame(dispnum = 1L, ccode = 260L, j = 20L), MAJ = data.frame(dispnum = 1L, ccode = 710L), ALLY = data.frame(dispnum = 1L, ccode = 900L, j = 10L))
check("17 tiers: same region (60), contiguous (260), major (710), allied (900) likely; 365 unlikely", identical(T$ccode, c(60L, 260L, 365L, 710L, 900L)) && identical(T$likely, c(TRUE, TRUE, FALSE, TRUE, TRUE)))

## ---- jointly fatal (B5, R3) ---------------------------------------------------------------------------------------------------------
CLf <- data.frame(dispnum = c(1, 1, 1, 2, 2, 3, 3, 3), ccode = c(1, 2, 3, 1, 2, 1, 2, 3), side_rel = c("A", "A", "B", "A", "B", "A", "B", "B"))
PF <- data.frame(dispnum = c(1, 1, 1, 2, 2, 3, 3, 3), ccode = c(1, 2, 3, 1, 2, 1, 2, 3), fat = c(0, -9, 3, 2, 1, 0, 4, NA))
J <- s2_jfatal(CLf, PF)
check("18 mixed side-level fatalities: A {0, unknown} vs B fatal -> unresolved", J[dispnum == 1, jfatal] == "unresolved")
check("18 both sides with a known fatal participant -> yes", J[dispnum == 2, jfatal] == "yes")
check("18 one side all known nonfatal -> no (unknown on the other side does not matter)", J[dispnum == 3, jfatal] == "no")

## ---- incidence (B2, R2) -------------------------------------------------------------------------------------------------------------
p <- rbind(c(0, 0.1, 0.2, 0.7), c(0.3, 0.1, 0, 0.6))
check("19 flat CIF: F_A = 0.7 * 0.3, F_B = 0.1 + 0.7 * 0.1, F_orig = 0.2", isTRUE(all.equal(unname(s2_cif(p)), c(0.21, 0.17, 0.2))))
check("19 sequential: a positive originator hazard lowers third-party A incidence through survival", {
  a <- s2_seq_cif(c(0.1, 0.2), c(0, 0), c(0.5, 0.5)); b <- s2_seq_cif(c(0.1, 0.2), c(0.3, 0), c(0.5, 0.5))
  isTRUE(all.equal(unname(b["A"]), 0.05 + 0.6 * 0.2 * 0.5)) && b["A"] < a["A"] })

## ---- initiating side (stage-1 rules) ------------------------------------------------------------------------------------------------
Pi <- data.frame(dispnum = c(5, 5, 5, 6, 6), ccode = c(1, 2, 3, 1, 2), spell = 1, spell_first = 1900, orig = c(1, 1, 0, 1, 1), side_entry = c(0, 1, 0, 1, 0))
EVi <- data.frame(dispnum = c(5, 5), initiator = c(1, 3), target = c(2, 2), war_link = NA, n_disputes = 1L, escalation = FALSE)
WIi <- data.frame(warnum = integer(0), ccode1 = integer(0), ccode2 = integer(0), sidea1 = integer(0), sidea2 = integer(0), initiator1 = integer(0), initiator2 = integer(0))
IS <- s2_init_side(Pi, EVi, WIi, c(5L, 6L))
check("20 rule 1: originator-originator event gives the initiator's side (0); joiner-initiated events ignored", IS$init_side[1] == 0L && IS$side_rule[1] == 1L)
check("20 rule 2: no originator event, not war-linked -> Side A (1)", IS$init_side[2] == 1L && IS$side_rule[2] == 2L)
check("20 disagreement with the case review stops (dispute 1035 must be side 0)", expect_error(s2_init_side(transform(Pi, dispnum = ifelse(dispnum == 6, 1035, dispnum)), EVi, WIi, c(5L, 1035L))))

## ---- legality table (C1, R3) --------------------------------------------------------------------------------------------------------
ref <- as.list(S2_REF)
check("21 reference cell is legal", !is.na(s2_cell_key(ref)))
check("21 one-shot x re-entries as new decisions is illegal", is.na(s2_cell_key(modifyList(ref, list(reentry = "new")))))
check("21 one-shot x monthly tiers / dynamic entrants / two-month bins illegal", all(is.na(c(s2_cell_key(modifyList(ref, list(tier_timing = "monthly"))), s2_cell_key(modifyList(ref, list(dynamic = "include"))), s2_cell_key(modifyList(ref, list(agg2 = "yes")))))))
check("21 monthly x truncated one-shot horizon illegal", is.na(s2_cell_key(modifyList(ref, list(estimator = "monthly", horizon = "truncated")))))
check("21 W = 0: retrospective and prospective cohorts share one key (deduplicated)", identical(s2_cell_key(modifyList(ref, list(window = "0", cohort = "retrospective"))), s2_cell_key(modifyList(ref, list(window = "0")))))

## ---- GPT build review fixtures (2026-11-22) -----------------------------------------------------------------------------------------
## 22. P1-1: stage-1 spells 1940 and 1942 (inactive 1941); one GML episode 1940-1942 must not bridge the gap
PTg <- data.frame(dispnum = 7L, ccode = 99L, spell = 1:2, spell_first = c(1940L, 1942L), spell_last = c(1940L, 1942L))
Gg <- data.frame(dispnum = 7L, ccode = 99L, styear = 1940L, endyear = 1942L, end_mk = s2_mk(1942L, 6L))
SE <- s2_spell_exits(PTg, Gg)
check("22 GML episode spanning two spells: spell 1 exit capped at Dec 1940 (flagged), spell 2 exit Jun 1942", SE$exit_mk[1] == s2_mk(1940L, 12L) && SE$exit_capped[1] && SE$exit_mk[2] == s2_mk(1942L, 6L) && !SE$exit_capped[2])
PSg <- data.frame(dispnum = 7L, ccode = c(99L, 99L, 1L, 2L), spell = c(1L, 2L, 1L, 1L), side_rel = c("A", "A", "A", "B"), baseline = c(FALSE, FALSE, TRUE, TRUE),
                  m_in = c(1L, 24L, 0L, 0L), m_out = c(SE$exit_mk - s2_mk(1940L, 1L), 40L, 40L))   # onset Jan 1940
RG <- s2_rosters(PSg, data.frame(dispnum = 7L, m = 0:40), "carry")
check("22 the participant is absent from A's roster throughout 1941 (months 13-24 use m - 1 = 12..23) and back from month 25", !any(RG$ccode == 99L & RG$m %in% 13:24) && all(25:30 %in% RG$m[RG$ccode == 99L]))
check("22 no GML match -> exit imputed at Dec of the last spell year (flagged)", { s0 <- s2_spell_exits(PTg, Gg[0, ]); all(s0$exit_imputed) && s0$exit_mk[2] == s2_mk(1942L, 12L) })
## 23. month keys: January / December, invalid codes; inclusive horizon rows
check("23 Jan 1900 and Dec 1900 keys are 22800 and 22811; Jan 1901 follows Dec 1900", s2_mk(1900L, 1L) == 22800L && s2_mk(1900L, 12L) == 22811L && s2_mk(1901L, 1L) == 22812L && s2_mk_year(22811L) == 1900L)
check("23 invalid month code (-9, 0, 13) stops", expect_error(s2_mk(1900L, -9L)) && expect_error(s2_mk(1900L, 0L)) && expect_error(s2_mk(1900L, 13L)))
check("23 a never-entering row has months 0..H inclusive (H + 1 rows, last row censored)", nrow(r60) == OS1[dispnum == 1 & ccode == 60, H] + 1L && tail(r60$y_m, 1) == "out")
## 24. jointly fatal: known fatal on both sides plus unknown participants -> yes
J2 <- s2_jfatal(data.frame(dispnum = 4, ccode = 1:4, side_rel = c("A", "A", "B", "B")), data.frame(dispnum = 4, ccode = 1:4, fat = c(2, -9, NA, 1)))
check("24 both sides with a known fatal participant and further unknowns -> yes", J2$jfatal == "yes")
## 25. a dispute starting Aug 1914 and ending 1919 takes fold 2 for all its rows
DT25 <- data.frame(dispnum = 25L, onset_lo = s2_mk(1914L, 8L), end_mk = s2_mk(1919L, 6L))
check("25 1914-1919 dispute: fold from the onset key = 2", s2_fold(s2_mk_year(DT25$onset_lo)) == 2L)
## 26. end policy truncates follow-up before the first empty side-month
MON26 <- data.table(dispnum = 1L, ccode = 60L, m = 0:10)
check("26 end policy: follow-up stops before month 7 (first empty B month)", identical(s2_truncate_end(MON26, attr(REe, "empty"))$m, 0:6))
## 27. conditional side equation respects the A mask (P(A | T) exactly 0 inside the window)
check("27 side probability is 0 where avail_A is FALSE, logistic elsewhere", { pp <- s2_side_prob(c(1, 1), c(0, 0), c(FALSE, TRUE)); pp[1] == 0 && abs(pp[2] - plogis(1)) < 1e-12 })

## ---- GPT recheck fixtures (2026-11-22) -------------------------------------------------------------------------------------------
## 28. endpoint provenance: a capped exit attaining the maximum is not "exact"
PTe <- data.table(dispnum = c(1L, 1L, 2L, 2L, 3L), exit_mk = c(s2_mk(1940L, 12L), s2_mk(1940L, 11L), 500L, 500L, 600L),
                  exit_imputed = c(FALSE, FALSE, FALSE, TRUE, TRUE), exit_capped = c(TRUE, FALSE, FALSE, FALSE, FALSE), exit_before_entry = FALSE)
DE <- s2_dispute_end(PTe)
check("28 capped GML exit (Jun 1942 capped to Dec 1940) sets the end: end_source capped, not exact", DE[dispnum == 1, end_source == "capped" & !end_exact & end_mk == s2_mk(1940L, 12L)])
check("28 tie at the maximum between an exact and an imputed exit -> exact; imputed only -> imputed", DE[dispnum == 2, end_source] == "exact" && DE[dispnum == 3, end_source] == "imputed")
## 29. prior conflict counts DISTINCT disputes per side: one historical MID with two roster members counts once; current and future excluded
PAIRS <- data.table(dispnum = 100L, ccode = 5L, j = c(1L, 2L), side_rel = "A", ty = 1950L)
PM <- data.table(pdisp = c(7L, 7L, 8L, 100L, 9L), a = 5L, b = c(1L, 2L, 1L, 1L, 2L), styear = c(1940L, 1940L, 1945L, 1951L, 1952L))
PM <- rbind(PM, PM[, .(pdisp, a = b, b = a, styear)])
check("29 prior MIDs: dispute 7 (both members) counts once, 8 counts, current (100) and future (9) excluded -> 2", s2_prior_side(PAIRS, PM)$prior_n == 2L)
check("29 no history -> 0, not NA", s2_prior_side(transform(PAIRS, ccode = 6L), PM)$prior_n == 0L)
## 30. strict aggregation
check("30 strict sum/min/any: a missing member gives NA; complete inputs aggregate", is.na(s2_strict_sum(c(0.02, NA))) && is.na(s2_strict_any(c(FALSE, NA))) &&
      is.na(s2_strict_min(c(3, NA))) && s2_strict_sum(c(0.02, 0.03)) == 0.05 && s2_strict_any(c(FALSE, TRUE)) == 1L && s2_strict_any(c(FALSE, FALSE)) == 0L)
## 31. reference-cohort exclusion precedence and the nonexclusive census
O31 <- data.table(complete = c(FALSE, TRUE, TRUE, TRUE), a = c(NA, NA, 1, 1), b = c(1, NA, NA, 1))
R31 <- s2_ref_cohort(O31, c("a", "b"))
check("31 precedence: incomplete horizon first, then the first missing column in list order", identical(R31$excl_reason, c("incomplete horizon", "missing a", "missing b", "")) && identical(R31$in_ref, c(FALSE, FALSE, FALSE, TRUE)))
check("31 census counts every missing value (a: 2 rows, 1 complete; b: 2 rows, 2 complete)", identical(attr(R31, "census")$n_missing_any_row, c(2L, 2L)) && identical(attr(R31, "census")$n_missing_complete_rows, c(1L, 2L)))
check("32 frozen design: every model term is an availability column or a fixed transform of one (lcap_s, share_B, surface squares/products); no tie maxima, no standalone own CINC",
      all(setdiff(unique(unlist(S2_REF_DESIGN)), c("lcap_A", "lcap_B", "share_B", S2_SURF_A, S2_SURF_B, "dem_x_A", "dem_x_B", "major_x_A", "major_x_B")) %in% S2_AVAIL) && !any(grepl("tie_max|log_cinc_k", unlist(S2_REF_DESIGN))))

## 33. stage-2 surface: gaps to the opposing side, stage-1 scaling, strict, side mean vs min
ROWS33 <- data.table(dispnum = 1L, ccode = 9L, ty = 1950L, t = 1951L)
MEM33 <- data.table(dispnum = 1L, side_rel = c("A", "A", "B"), j = c(1L, 2L, 3L))
PV33 <- data.table(ccode = c(9L, 1L, 2L, 3L), year = 1950L, M = c(-3, -2, -4, -5), R = c(-4, -3, -5, -6))
S33 <- s2_surface(ROWS33, MEM33, PV33, sM = 2, sR = c("1951" = 0.5))
check("33 dM_kA = (-3 - mean(-2, -4)) / 2 = 0; dR_kA = (-4 - (-4)) / 0.5 = 0; dM_kB = (-3 + 5) / 2 = 1; dR_kB = (-4 + 6) / 0.5 = 4",
      isTRUE(all.equal(c(S33$dM_kA, S33$dR_kA, S33$dM_kB, S33$dR_kB), c(0, 0, 1, 4))))
check("33 join A uses gaps to B (sA_dR = 4, sA_dR2 = 16, sA_dMR = 4); join B uses gaps to A (all 0)", S33$sA_dR == 4 && S33$sA_dR2 == 16 && S33$sA_dMR == 4 && S33$sB_dM == 0 && S33$sB_dR == 0)
check("33 min aggregation: weakest / most isolated A member (-4, -5) -> dM_kA = 0.5, dR_kA = 2", { s <- s2_surface(ROWS33, MEM33, PV33, 2, c("1951" = 0.5), agg = "min"); isTRUE(all.equal(c(s$dM_kA, s$dR_kA), c(0.5, 2))) })
check("33 a side member without a position makes that side's gaps NA (strict)", { pv <- copy(PV33); pv[ccode == 2L, R := NA]; s <- s2_surface(ROWS33, MEM33, pv, 2, c("1951" = 0.5)); is.na(s$dR_kA) && !is.na(s$dR_kB) })
check("33 frozen surface design: each alternative carries its 5 surface terms, no linear position terms (own_cent, tgt_cent_mean, log_cinc_ratio)",
      all(S2_SURF_A %in% S2_REF_DESIGN$A) && all(S2_SURF_B %in% S2_REF_DESIGN$B) && all(S2_SURF_A %in% S2_REF_DESIGN$O) && !any(grepl("own_cent|tgt_cent|log_cinc_ratio", unlist(S2_REF_DESIGN))))
check("33 mirroring: join A / originator take the stage-1 controls toward B (the side they fight); join B toward A; recruitment terms toward the joined side",
      all(S2_Z_TOWARD("B") %in% S2_REF_DESIGN$A) && all(S2_Z_TOWARD("A") %in% S2_REF_DESIGN$B) && identical(S2_REF_DESIGN$O, S2_REF_DESIGN$A) && "tie_sum_A" %in% S2_REF_DESIGN$A && "tie_sum_B" %in% S2_REF_DESIGN$B)
check("33 sweep dimensions: opposing-side position mean/min and design surface/linear; reference = mean, surface", S2_REF[["opp_side"]] == "mean" && S2_REF[["design"]] == "surface" && is.null(S2_DIMS$target))

cat(if (ok_all) "ALL PASS\n" else "SOME TESTS FAILED\n"); quit(status = if (ok_all) 0 else 1)
