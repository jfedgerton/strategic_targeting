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
TI <- s2_ties(ROW, RC, E)
check("15 ties of 60 at month 1: A = {10, 40} -> sum 1.4, max 0.9; B = {20, 30} -> sum 0.2, max 0.2 (no link = 0)",
      TI[ccode == 60, isTRUE(all.equal(c(tie_sum_A, tie_max_A, tie_sum_B, tie_max_B), c(1.4, 0.9, 0.2, 0.2)))])
check("15 roster member 40 in A: only the 10-40 link counts (edges are a < b, so self-ties cannot arise; the filter is defensive)", TI[ccode == 40, isTRUE(all.equal(c(tie_sum_A, tie_max_A), c(0.7, 0.7)))])
check("15 wrong tie year gives zeros, not NA", { t2 <- s2_ties(transform(ROW, year_tie = 7L), RC, E); all(t2$tie_sum_A == 0) && !anyNA(t2) })

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

cat(if (ok_all) "ALL PASS\n" else "SOME TESTS FAILED\n"); quit(status = if (ok_all) 0 else 1)
