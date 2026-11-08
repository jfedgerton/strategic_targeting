## The interacting-dyads event foundation (switch events = pairs | dymid, suffix _dy, build variant v2/data_dymid; PI-approved 2026-11-04; README 6.9).
## (a) helpers_dymid.R on hand-built DYMID and war rows: the dyad key, max-over-years fatality, re-entries, linking (merged / war-only / DYMID-only
## war), recip validation (2 accepted explicitly as not reciprocated, other values refused), the samples strict / any / ub (definition of
## 2026-11-05), the direction rules in the order (a) / (b) / (c) / fallback with ties (PI decision 2026-11-05), the participation table and the spell override,
## the copy validation and the adjudication table (2026-11-08: one toy case per rule, an uncovered disagreement, stale / inconsistent rows), the re-entry
## spells (2026-11-08: a re-entering dyad gives two spells, the gap year is not already fighting);
## (b) the switch, suffix, data directory, samples and refusals; strings and names unchanged at the default; (c) 02_fit_model.R end to end on a toy
## dymid cell (M1b, M1c bcj with joiner_set, M2c from the matching M1c); (d) the real dymid build (skipped when absent).
suppressPackageStartupMessages({ library(testthat); library(dplyr) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "assemble_cell.R")); source(file.path(root, "R", "helpers_dymid.R"))

## ---- toy DYMID rows (A < B copies; read_dymid() checks the mirror on a CSV below) -----------------------------------------------------------
## dispute 1 (2000-2001): 10 (side A) and 20 (side B) originators; joiners 30 (side B) and 40 (side A)
##   1: 10-20 two years; fatleva -9 in 2000 (listed only in the last year), 2 in 2001; fatlevb 1 in 2001 -> max 2 / 1, strict fatal; dyad_rolea 3 (20 initiated)
##   1: 10-30 fatleva 1, fatlevb -9 -> strict FALSE, ub TRUE, any TRUE (fatlev 1)
##   1: 30-40 both joiners, non-fatal; 1: 20-40 DYMID war == 1 without a war-data row (DYMID-only war dyad)
## dispute 2 (1990): 50-60 both originators, war dyad also in the war data (merged); dyad_rolea 1 (50) but D1 says 60 (rule (b) decides: 50 -> 60)
## dispute 5: 13-14 with two spells (2000-03 and 2002-06): one event dated 2000-03
## dispute 4: 11-12 with dyad_rolea 2 (simultaneous): no rule decides
dyrow <- function(disno, a, b, year, sy, sm, ey, em, fl = 0L, fa = 0L, fb = 0L, war = 0L, ra = 1L, oa = 1L, ob = 1L, sa = 1L, recip = 0L)
  data.frame(disno = disno, dyindex = disno, statea = a, stateb = b, year = year, strtday = 1L, strtmnth = sm, strtyr = sy, endday = 28L, endmnth = em, endyear = ey,
             fatlev = fl, hihost = 4L, recip = recip, sideaa = sa, sideab = 1L - sa, fatleva = fa, fatlevb = fb, hihosta = 4L, hihostb = 4L, orignata = oa, orignatb = ob,
             rolea = 1L, roleb = 3L, dyad_rolea = ra, dyad_roleb = 4L - ra, war = war)
toy_dy <- function() {
  r <- rbind(dyrow(1L, 10L, 20L, 2000L, 2000L, 3L, 2001L, 6L, fl = -9L, fa = -9L, fb = 0L, ra = 3L),
             dyrow(1L, 10L, 20L, 2001L, 2000L, 3L, 2001L, 6L, fl = 3L, fa = 2L, fb = 1L, ra = 3L),
             dyrow(1L, 10L, 30L, 2000L, 2000L, 5L, 2000L, 9L, fl = 1L, fa = 1L, fb = -9L, ob = 0L),
             dyrow(1L, 30L, 40L, 2000L, 2000L, 8L, 2000L, 9L, oa = 0L, ob = 0L, sa = 0L),
             dyrow(1L, 20L, 40L, 2000L, 2000L, 8L, 2000L, 12L, fl = 6L, fa = 6L, fb = 6L, war = 1L, ob = 0L, sa = 0L),
             dyrow(2L, 50L, 60L, 1990L, 1990L, 2L, 1990L, 7L, fl = 6L, fa = 6L, fb = 6L, war = 1L, ra = 1L),
             dyrow(5L, 13L, 14L, 2000L, 2000L, 3L, 2000L, 4L), dyrow(5L, 13L, 14L, 2002L, 2002L, 6L, 2002L, 8L, fl = 1L, fa = 1L, fb = 1L),
             dyrow(4L, 11L, 12L, 1995L, 1995L, 1L, 1995L, 2L, ra = 2L))
  mirror_dy(r)
}
## the mirror copies (A > B) of hand-built A < B rows, as the file lists them
mirror_dy <- function(r) {
  m <- r; m$statea <- r$stateb; m$stateb <- r$statea; m$fatleva <- r$fatlevb; m$fatlevb <- r$fatleva; m$sideaa <- r$sideab; m$sideab <- r$sideaa
  m$orignata <- r$orignatb; m$orignatb <- r$orignata; m$dyad_rolea <- r$dyad_roleb; m$dyad_roleb <- r$dyad_rolea; m$hihosta <- r$hihostb; m$hihostb <- r$hihosta
  rbind(r, m)
}
## toy war rows: 2: 50-60 (war 7, linked to the DYMID dyad); 3: 90-95 (war 9, no DYMID dyad: war-only, wardyadrolea 1 for 90)
toy_wr <- function() {
  w <- data.frame(warnum = c(7L, 9L), disno = c(2L, 3L), dyindex = c(2.001, 3.001), statea = c(50L, 90L), stateb = c(60L, 95L), warstrtmnth = c(4L, 6L), warstrtday = 1L,
                  warstrtyr = c(1990L, 1980L), warendmnth = c(7L, 3L), warenday = 1L, warendyr = c(1990L, 1981L), year = c(1990L, 1980L), warolea = c(1L, 1L), waroleb = c(3L, 3L),
                  wardyadrolea = c(1L, 1L), wardyadroleb = c(3L, 3L), batdtha = 600L, batdthb = 600L, batdths = 1200L)
  w <- rbind(w, transform(w[2, ], year = 1981L))
  m <- w; m$statea <- w$stateb; m$stateb <- w$statea; m$warolea <- w$waroleb; m$waroleb <- w$warolea; m$wardyadrolea <- w$wardyadroleb; m$wardyadroleb <- w$wardyadrolea
  rbind(w, m)
}

test_that("read_dymid / read_ddiw: structure checks, the mirror check, the day / month repair", {
  td <- tempfile("dymidraw"); dir.create(td)
  write.csv(toy_dy(), file.path(td, "dy.csv"), row.names = FALSE); write.csv(toy_wr(), file.path(td, "wr.csv"), row.names = FALSE)
  dy <- read_dymid(file.path(td, "dy.csv")); expect_equal(nrow(dy), 18L); expect_equal(attr(dy, "n_start_mismatch"), 0L)
  wr <- read_ddiw(file.path(td, "wr.csv")); expect_equal(attr(wr, "n_swapped"), 0L)
  bad <- toy_dy(); write.csv(bad[-nrow(bad), ], file.path(td, "bad.csv"), row.names = FALSE)
  expect_error(read_dymid(file.path(td, "bad.csv")), "do not mirror")
  b2 <- toy_dy(); b2$dyindex[1] <- 1001L; write.csv(b2, file.path(td, "bad2.csv"), row.names = FALSE); expect_error(read_dymid(file.path(td, "bad2.csv")), "dyindex differs")
  sw <- toy_wr(); sw$warstrtmnth[c(1, 4)] <- 24L; sw$warstrtday[c(1, 4)] <- 4L; sw$warendmnth[c(1, 4)] <- 1L; sw$warenday[c(1, 4)] <- 7L
  write.csv(sw, file.path(td, "sw.csv"), row.names = FALSE); w2 <- read_ddiw(file.path(td, "sw.csv"))
  expect_equal(attr(w2, "n_swapped"), 2L); expect_equal(w2$warstrtmnth[1], 4L); expect_equal(w2$warstrtday[1], 24L); expect_equal(w2$warendmnth[1], 7L)
  ## recip: 0 / 1 per the codebook; 2 (undocumented, 4 rows of 4.03) accepted and treated as NOT reciprocated, explicitly; anything else refused
  r2 <- toy_dy(); r2$recip[r2$disno == 4L] <- 2L; write.csv(r2, file.path(td, "r2.csv"), row.names = FALSE)
  d2 <- read_dymid(file.path(td, "r2.csv")); expect_equal(nrow(attr(d2, "recip2_rows")), 2L); expect_equal(nrow(attr(dy, "recip2_rows")), 0L)
  dd2 <- dymid_dyads(d2); expect_false(dd2$recip[dd2$disno == 4L]); expect_true(dd2$recip_code2[dd2$disno == 4L]); expect_false(any(dd2$recip_code2[dd2$disno != 4L]))
  r1 <- toy_dy(); r1$recip[r1$disno == 4L] <- 1L; write.csv(r1, file.path(td, "r1.csv"), row.names = FALSE)
  dd1 <- dymid_dyads(read_dymid(file.path(td, "r1.csv"))); expect_true(dd1$recip[dd1$disno == 4L]); expect_false(dd1$recip_code2[dd1$disno == 4L])
  r3 <- toy_dy(); r3$recip[1] <- 3L; write.csv(r3, file.path(td, "r3.csv"), row.names = FALSE); expect_error(read_dymid(file.path(td, "r3.csv")), "recip takes values")
  rn <- toy_dy(); rn$recip[1] <- -9L; write.csv(rn, file.path(td, "rn.csv"), row.names = FALSE); expect_error(read_dymid(file.path(td, "rn.csv")), "recip takes values")
  ra <- toy_dy(); ra$recip[1] <- NA; write.csv(ra, file.path(td, "ra.csv"), row.names = FALSE); expect_error(read_dymid(file.path(td, "ra.csv")), "recip takes values")
  unlink(td, recursive = TRUE)
})

test_that("dymid_dyads / ddiw_dyads / dymid_link: dyad key, max-over-years fatality, re-entry, merged / war-only / DYMID-only war", {
  d <- dymid_dyads(toy_dy()); w <- ddiw_dyads(toy_wr())
  expect_equal(nrow(d), 7L); expect_equal(nrow(w), 2L)
  r <- d[d$disno == 1L & d$a == 10L & d$b == 20L, ]
  expect_equal(c(r$fatlev_max, r$fatleva_max, r$fatlevb_max), c(3L, 2L, 1L)); expect_equal(r$onset_year, 2000L); expect_equal(r$start_ym, ym_index(2000L, 3L))
  expect_equal(r$init_dyadic, 20L)                                                         # dyad_rolea 3 for a = 10: 20 initiated
  rr <- d[d$disno == 5L, ]; expect_equal(rr$n_spells, 2L); expect_equal(rr$start_ym, ym_index(2000L, 3L)); expect_equal(rr$end_ym, ym_index(2002L, 8L)); expect_equal(rr$fatleva_max, 1L)
  expect_true(is.na(d$init_dyadic[d$disno == 4L]))                                        # simultaneous
  expect_true(is.na(d$fatlevb_max[d$disno == 1L & d$b == 30L]))
  e <- dymid_link(d, w); lc <- attr(e, "counts")
  expect_equal(nrow(e), 8L); expect_equal(lc$n_merged, 1L); expect_equal(lc$n_ddiw_only, 1L); expect_equal(lc$n_dymid_only_war, 1L)
  m <- e[e$disno == 2L, ]; expect_equal(m$source, "merged"); expect_true(m$war); expect_equal(m$warnum_ddiw, 7L); expect_equal(m$onset_year, 1990L); expect_equal(m$start_ym, ym_index(1990L, 2L))  # DYMID dates, not the war start
  wo <- e[e$disno == 3L, ]; expect_equal(wo$source, "ddiw_only"); expect_equal(wo$onset_year, 1980L); expect_equal(wo$start_ym, ym_index(1980L, 6L)); expect_equal(wo$init_b, 90L)
  expect_true(e$war[e$disno == 1L & e$a == 20L & e$b == 40L]); expect_false(e$war[e$disno == 1L & e$a == 10L & e$b == 30L])
  y <- attr(e, "years"); expect_equal(sum(y$disno == 3L), 2L); expect_equal(sum(y$disno == 1L & y$a == 10L & y$b == 20L), 2L)
  ## two war dyads on one (disno, pair) stop the link
  w2 <- bind_rows(w, w[w$disno == 2L, ] |> mutate(warnum = 8L)); expect_error(dymid_link(d, w2), "two war dyads")
})

test_that("dymid_fatal_codes: strict / any / ub (definition of 2026-11-05), missing sides, wars", {
  ## 1 both sides >= 1; 2 one side >= 1, other missing, fatlev 1; 3 all missing; 4 non-fatal; 5 war; 6 only the dyadic fatlev;
  ## 7 one side >= 1, other missing, fatlev 0; 8 both sides missing, fatlev 0; 9 both missing, fatlev 2; 10 one side 0, other missing, fatlev 1;
  ## 11 both sides >= 1 while fatlev 0 (inconsistent dyadic codes: strict decides); 12 one side missing, other >= 1, fatlev missing
  fc <- dymid_fatal_codes(fatleva = c(2L, 1L, NA, 0L, NA, 0L, 1L, NA, NA, 0L, 1L, NA), fatlevb = c(1L, NA, NA, 0L, NA, 0L, NA, NA, NA, NA, 1L, 2L),
                          fatlev = c(3L, 1L, NA, 0L, NA, 2L, 0L, 0L, 2L, 1L, 0L, NA), war = c(FALSE, FALSE, FALSE, FALSE, TRUE, rep(FALSE, 7)))
  expect_equal(fc$strict, c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE))
  expect_equal(fc$any,    c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE))       # 6: only the dyadic dispute's fatlev
  expect_equal(fc$ub,     c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE))
  expect_false(fc$ub[7]); expect_false(fc$ub[8])     # dyadic fatlev = 0 (coded non-fatal): a missing per-state field never makes it fatal
  expect_true(fc$ub[2]); expect_true(fc$ub[12])      # one side >= 1 and the other missing IS fatal (fatlev >= 1 or missing)
  expect_equal(fc$ub_branch, c("strict", "one_side_missing", "both_missing", "none", "strict", "none", "none", "none", "both_missing", "none", "strict", "one_side_missing"))
  expect_true(all(fc$ub[fc$strict])); expect_true(all(fc$any[fc$strict]))
  ## a fatlev = 0 dyad with BOTH per-state fields missing is not fatal (the previous definition counted it)
  expect_false(dymid_fatal_codes(NA_integer_, NA_integer_, 0L, FALSE)$ub)
})

## psd rows for the entry comparison (as in test_entry_roles.R)
psd_d <- function(disp, cc, res, m, y = 2000L, d = NA_integer_) { k <- date_keys(res, y, m, d)
  tibble(dispnum = disp, ccode = as.integer(cc), spell = 1L, spell_first = y, spell_last = y, entry_res = as.integer(res), entry_y = y, entry_mkey = k$mkey, entry_dkey = k$dkey) }
ev_d <- function(disp, i, j, o1, o2, war = FALSE, d1 = NA_integer_, b = NA_integer_, y = 2000L)
  tibble(dispnum = disp, initiator = as.integer(i), target = as.integer(j), year = y, sidea1 = 1L, sidea2 = 0L, orig1 = as.integer(o1), orig2 = as.integer(o2),
         fatality1 = 1L, fatality2 = 2L, hostlev1 = 3L, hostlev2 = 4L, war = war, d1_init = as.integer(d1), init_b = as.integer(b))

test_that("dymid_direction: order (a) -> (b) -> (c): (a) one joiner / two joiners by entry date, ties -> (b); (b) dyadic initiator, also for war dyads where D1 disagrees; (c) D1 only without (b); fallback", {
  ev <- bind_rows(ev_d(1L, 10L, 20L, 1, 1, b = 20L),                      # both originators, non-war: (b) -> 20 -> 10
                  ev_d(1L, 10L, 30L, 1, 0, b = 10L),                      # one joiner: (a) -> 30 -> 10 whatever (b) says
                  ev_d(1L, 40L, 30L, 0, 0, b = 30L),                      # two joiners: 30 May, 40 August -> 40 -> 30 (later entrant) although (b) says 30
                  ev_d(1L, 40L, 20L, 0, 1, war = TRUE, d1 = 20L),         # one joiner, war: (a) beats (c) -> 40 -> 20
                  ev_d(2L, 50L, 60L, 1, 1, war = TRUE, d1 = 60L, b = 50L, y = 1990L),   # war, both originators, D1 says 60 but (b) says 50: (b) decides -> 50 -> 60
                  ev_d(2L, 70L, 80L, 0, 0, b = 80L, y = 1990L),           # two joiners, same entry month (tie): (b) -> 80 -> 70
                  ev_d(2L, 71L, 81L, 0, 0, war = TRUE, d1 = 71L, b = 81L, y = 1990L),   # tie in a war dyad, D1 and (b) disagree: (b) -> 81 -> 71
                  ev_d(4L, 11L, 12L, 1, 1, y = 1995L),                    # simultaneous, no D1: Side A kept
                  ev_d(6L, 15L, 16L, 1, 1, war = FALSE, d1 = 16L, b = 15L),   # non-war dyad in a war-level dispute: (b) 15 -> 16
                  ev_d(7L, 17L, 18L, 1, 1, war = TRUE, d1 = 18L))         # war dyad without a dyadic initiator: (c) D1 -> 18 -> 17
  psd <- bind_rows(psd_d(1L, 10, 2L, 1L), psd_d(1L, 20, 2L, 1L), psd_d(1L, 30, 2L, 5L), psd_d(1L, 40, 2L, 8L),
                   psd_d(2L, 50, 2L, 2L, 1990L), psd_d(2L, 60, 2L, 2L, 1990L), psd_d(2L, 70, 2L, 3L, 1990L), psd_d(2L, 80, 2L, 3L, 1990L), psd_d(2L, 71, 2L, 3L, 1990L), psd_d(2L, 81, 2L, 3L, 1990L),
                   psd_d(4L, 11, 2L, 1L, 1995L), psd_d(4L, 12, 2L, 1L, 1995L), psd_d(6L, 15, 2L, 1L), psd_d(6L, 16, 2L, 1L), psd_d(7L, 17, 2L, 1L), psd_d(7L, 18, 2L, 1L))
  r <- dymid_direction(ev, psd)
  expect_equal(paste(r$initiator, r$target), c("20 10", "30 10", "40 30", "40 20", "50 60", "80 70", "81 71", "11 12", "15 16", "18 17"))
  expect_equal(r$dir_source, c("b_dyadic", "a_joiner_orig", "a_joiner_later", "a_joiner_orig", "b_dyadic", "b_dyadic", "b_dyadic", "fallback_sideA", "b_dyadic", "c_cowwar"))
  expect_equal(r$dir_rule[6:7], c("b_dyadic (joiner tie)", "b_dyadic (joiner tie)"))
  expect_equal(r$dir_flipped, r$sidea1 == 0L); expect_equal(r$dir_flipped, c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE))
  expect_equal(r$b_c_disagree, c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE))   # D1 decides and differs from (b), any dyad
  expect_equal(r$dir_changed, c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE))   # war dyads the previous (c)-first order gave to D1
  expect_equal(names(DYMID_DIR_LABELS), c("a_joiner_orig", "a_joiner_later", "b_dyadic", "c_cowwar", "fallback_sideA"))
  ## reoriented rows swap every pair column; no late-target event remains
  expect_equal(unlist(r[1, c("orig1", "orig2", "fatality1", "fatality2", "hostlev1", "hostlev2")]), c(orig1 = 1, orig2 = 1, fatality1 = 2, fatality2 = 1, hostlev1 = 4, hostlev2 = 3))
  expect_false(any(r$orig1 == 1L & r$orig2 == 0L))
  expect_error(dymid_direction(ev |> mutate(d1_init = ifelse(dispnum == 2L, 99L, d1_init)), psd))   # a D1 initiator outside the dyad
})

test_that("dymid_participation / dymid_disputes / dymid_override_pdi: COW pairs kept or replaced, DYMID-only pairs and disputes added, sides, spells re-dated", {
  ## COW: dispute 1 in 2000 with 10, 40 (side A) and 20, 30, 99 (side B); 30, 40 joiners; 99 never interacts. The opposing COW pairs are written both ways
  sides <- c(`10` = 1L, `40` = 1L, `20` = 0L, `30` = 0L, `99` = 0L); orig <- c(`10` = 1L, `40` = 0L, `20` = 1L, `30` = 0L, `99` = 1L)
  pr <- expand.grid(c1 = c(10L, 40L), c2 = c(20L, 30L, 99L))
  dd <- bind_rows(tibble(dispnum = 1L, ccode1 = pr$c1, ccode2 = pr$c2), tibble(dispnum = 1L, ccode1 = pr$c2, ccode2 = pr$c1)) |>
    mutate(year = 2000L, sidea1 = sides[as.character(ccode1)], sidea2 = sides[as.character(ccode2)], orig1 = orig[as.character(ccode1)], orig2 = orig[as.character(ccode2)])
  w0 <- ddiw_dyads(toy_wr()); w0e <- w0[0, ]; attr(w0e, "years") <- attr(w0, "years")[0, ]
  ## dispute 1 plus a DYMID-only dispute 7 (15 vs 16, unknown to COW)
  d <- dymid_dyads(bind_rows(toy_dy() |> filter(disno == 1L), mirror_dy(dyrow(7L, 15L, 16L, 2000L, 2000L, 4L, 2000L, 5L, sa = 0L)))); e <- dymid_link(d, w0e)
  p <- dymid_participation(dd, e, attr(e, "years")); dg <- attr(p, "diagnostics")
  ## COW pairs 10-20, 10-30, 10-99, 20-40, 30-40, 40-99: the four event dyads replaced, 10-99 and 40-99 (no interaction) kept as COW has them
  expect_equal(dg$n_cow_pairs, 6L); expect_equal(dg$n_cow_pairs_replaced, 4L); expect_equal(dg$n_cow_pairs_not_event, 2L); expect_equal(dg$n_event_pairs_added, 1L); expect_equal(dg$n_disputes_added, 1L)
  expect_true(all(c("10 99", "40 99") %in% paste(p$ccode1, p$ccode2)))
  expect_equal(p$sidea1[p$dispnum == 7L & p$ccode1 == 15L], 0L)                          # not a COW participant: DYMID's side (sideaa 0)
  expect_true(all(p$sidea1 + p$sidea2 == 1L))
  expect_equal(unique(p$orig1[p$ccode1 == 40L]), 0L); expect_equal(unique(p$sidea1[p$ccode1 == 40L & p$year == 2000L]), 1L)
  expect_equal(sum(p$ccode1 == 10L & p$ccode2 == 20L), 2L)                                 # the event dyad's own years (2000, 2001), both directions
  ds <- dymid_disputes(tibble(dispnum = 1L, styear = 2000L, stmon = 3L, mindur = 400, maxdur = 460, hiact = 18L, hostlev = 4L, fatality = 3L),
                       bind_rows(toy_dy() |> filter(disno == 1L), mirror_dy(dyrow(7L, 15L, 16L, 2000L, 2000L, 4L, 2000L, 5L, sa = 0L))), e, p)
  expect_equal(ds$hostlev_ev[ds$dispnum == 1L], 5L)                                        # 20 - 40 is a war dyad
  db <- dispute_date_bounds(ds, p); pdi <- pair_dispute_intervals(p, db, part = NULL)
  po <- dymid_override_pdi(pdi, attr(e, "spells"))
  s <- po[po$ccode1 == 10L & po$ccode2 == 30L, ]; expect_equal(c(s$start_lo, s$end_hi), c(ym_index(2000L, 5L), ym_index(2000L, 9L))); expect_equal(s$date_source, "dymid")
  expect_equal(po$date_source[po$ccode1 == 10L & po$ccode2 == 99L], "dispute")            # a co-participant pair keeps the rule
  expect_equal(attr(po, "n_dymid"), 2L * 5L)                                               # 5 event dyads, both directions, one spell each
  expect_equal(attr(po, "n_split_rows"), 0L); expect_equal(attr(po, "n_split_pairs"), 0L)    # no re-entry here: no row added (2026-11-08)
  expect_equal(ds$hostlev_ev[ds$dispnum == 7L], 4L); expect_equal(nrow(attr(ds, "constructed")), 1L)
})

## ---- copy disagreements and the adjudication table (2026-11-08) --------------------------------------------------------------------------------
## toy_dy() with the A > B copy of chosen rows altered: the (B, A) row of a dyad-year says something else than the (A, B) row
disagree_dy <- function(what = c("direction", "war", "fatal", "fatal_missing", "start", "end_last", "end_both", "sides", "hostility", "orig")) {
  what <- match.arg(what); d <- toy_dy(); m <- d$statea > d$stateb
  if (what == "direction") { k <- m & d$disno == 1L & d$statea == 20L & d$stateb == 10L; d$dyad_rolea[k] <- 3L; d$dyad_roleb[k] <- 1L }   # 20-10 copy: 10 initiated; the 10-20 copy (dyad_rolea 3) says 20 did
  if (what == "war") d$war[m & d$disno == 1L & d$statea == 30L & d$stateb == 10L] <- 1L                          # 30-10 copy flags war, 10-30 copy does not; no war-data dyad
  if (what == "fatal") { d$fatlevb[m & d$disno == 5L & d$year == 2002L] <- 0L; d$fatlev[m & d$disno == 5L & d$year == 2002L] <- 0L }   # 14-13 copy 2002: 13's deaths (its fatlevb) 0 (copy1: 1), fatlev 0 (copy1: 1)
  if (what == "fatal_missing") d$fatlevb[m & d$disno == 5L] <- -9L                                             # 14-13 copy, both years: 13's deaths missing (copy1: max 1)
  if (what == "start") d$strtmnth[m & d$disno == 2L] <- 4L                                                       # 60-50 copy: start April (copy1: February)
  if (what == "end_last") d$endmnth[m & d$disno == 5L & d$year == 2002L] <- 9L                                   # 14-13 copy, second spell: end September (copy1: August) = the dyad's end
  if (what == "end_both") d$endmnth[m & d$disno == 5L] <- c(5L, 9L)[match(d$year[m & d$disno == 5L], c(2000L, 2002L))]   # ... and the FIRST spell's end May (copy1: April), not visible at the dyad level
  if (what == "sides") { d$sideaa[m & d$disno == 1L & d$statea == 40L & d$stateb == 30L] <- 0L; d$sideab[m & d$disno == 1L & d$statea == 40L & d$stateb == 30L] <- 1L }   # 40-30 copy: sides inverted
  if (what == "hostility") d$hihosta[m & d$disno == 2L] <- 5L                                                   # 60-50 copy: 60 at hostility 5 (copy1: 4)
  if (what == "orig") d$orignata[m & d$disno == 2L] <- 0L                                                       # 60-50 copy: 60 a joiner (copy1: originator): no rule for orig_a_dy / orig_b_dy
  d
}
adj_row <- function(disno, a, b, field, v1, v2, res, rule, applies_to = "all", note = "")
  data.frame(disno = disno, a = a, b = b, field = field, value_copy1 = v1, value_copy2 = v2, resolved_value = res, applies_to = applies_to, rule = rule, note = note, stringsAsFactors = FALSE)
war_pairs_toy <- ddiw_dyads(toy_wr()) |> distinct(disno, a, b)

test_that("copy validation: identical copies pass with an empty table; dymid_copies mirrors; the disagreement list names dyad, field and both values", {
  d <- dymid_dyads(toy_dy()); expect_false(any(d$copy_adjudicated)); expect_equal(nrow(attr(d, "adjudication")), 0L)
  expect_true(all(c("fatlev_ub", "fatleva_ub", "fatlevb_ub", "adj_fields") %in% names(d))); expect_equal(d$fatleva_ub, d$fatleva_max); expect_equal(d$fatlevb_ub, d$fatlevb_max)
  sp <- attr(d, "spells"); expect_equal(nrow(sp), 8L); expect_equal(sp$spell[sp$disno == 5L], 1:2); expect_equal(sp$y_first[sp$disno == 5L], c(2000L, 2002L))
  cp <- dymid_copies(toy_dy()); expect_equal(nrow(cp$copy1), 9L); expect_true(all(cp$copy2$statea < cp$copy2$stateb)); expect_equal(cp$copy2$fatleva, cp$copy1$fatleva)
  dis <- dymid_copy_disagreements(dymid_dyads_copy(cp$copy1), dymid_dyads_copy(cp$copy2)); expect_equal(nrow(dis), 0L)
  cp <- dymid_copies(disagree_dy("direction")); dis <- dymid_copy_disagreements(dymid_dyads_copy(cp$copy1), dymid_dyads_copy(cp$copy2))
  expect_equal(nrow(dis), 1L); expect_equal(dis$field, "init_dyadic"); expect_equal(c(dis$value_copy1, dis$value_copy2), c("20", "10"))
  cp <- dymid_copies(disagree_dy("fatal_missing")); dis <- dymid_copy_disagreements(dymid_dyads_copy(cp$copy1), dymid_dyads_copy(cp$copy2))
  expect_equal(dis$field, "fatleva_max"); expect_equal(c(dis$value_copy1, dis$value_copy2), c("1", "-9"))
})

test_that("copy validation: an uncovered disagreement STOPS; so do a stale row, wrong copy values, a wrong resolved value, a rule that is not the field's", {
  expect_error(dymid_dyads(disagree_dy("direction")), "not covered by the adjudication table.*dispute 1 10-20 init_dyadic \\(copy1 20, copy2 10\\)")
  expect_error(dymid_dyads(disagree_dy("war"), war_pairs = war_pairs_toy), "not covered.*war_dymid")
  expect_error(dymid_dyads(disagree_dy("fatal")), "not covered")
  expect_error(dymid_dyads(disagree_dy("orig")), "not covered.*orig_b_dy")                      # a field with no rule: no resolution, stop
  ok <- adj_row(1L, 10L, 20L, "init_dyadic", "20", "10", "NA", "direction_uninformative")
  expect_error(dymid_dyads(toy_dy(), adj = ok), "match no copy disagreement")                     # stale: the data agree
  expect_error(dymid_dyads(disagree_dy("direction"), adj = adj_row(1L, 10L, 20L, "init_dyadic", "20", "20", "NA", "direction_uninformative")), "copy values differ from the data")
  expect_error(dymid_dyads(disagree_dy("direction"), adj = adj_row(1L, 10L, 20L, "init_dyadic", "20", "10", "20", "direction_uninformative")), "does not follow from the rule")
  expect_error(dymid_dyads(disagree_dy("direction"), adj = adj_row(1L, 10L, 20L, "init_dyadic", "20", "10", "NA", "hostility_higher")), "not a rule for field init_dyadic")
  expect_error(dymid_dyads(disagree_dy("war"), adj = adj_row(1L, 10L, 30L, "war_dymid", "0", "1", "0", "war_defer_to_war_data")), "needs war_pairs")
  ## a fatality disagreement needs both rules; the applies_to must be the rule's
  f1 <- adj_row(5L, 13L, 14L, "fatleva_max", "1", "0", "0", "fatal_strict_lower", "fmid"); f2 <- adj_row(5L, 13L, 14L, "fatleva_max", "1", "0", "1", "fatal_ub_higher", "fmidany,fmidub")
  l1 <- adj_row(5L, 13L, 14L, "fatlev_max", "1", "0", "0", "fatal_strict_lower", "fmid"); l2 <- adj_row(5L, 13L, 14L, "fatlev_max", "1", "0", "1", "fatal_ub_higher", "fmidany,fmidub")
  expect_error(dymid_dyads(disagree_dy("fatal"), adj = bind_rows(f1, l1, l2)), "needs one fatal_strict_lower row \\(fmid\\) and one fatal_ub_higher row")
  expect_error(dymid_dyads(disagree_dy("fatal"), adj = bind_rows(f1 |> mutate(applies_to = "all"), f2, l1, l2)), "needs applies_to = fmid")
  expect_silent(invisible(dymid_dyads(disagree_dy("fatal"), adj = bind_rows(f1, f2, l1, l2))))
  ## the table reader: columns, unknown rule, unknown field, duplicated rows
  td <- tempfile("adj"); dir.create(td)
  write.csv(ok, file.path(td, "a.csv"), row.names = FALSE); a <- dymid_adjudication_read(file.path(td, "a.csv")); expect_equal(a$resolved_value, "NA"); expect_type(a$disno, "integer")
  write.csv(ok |> mutate(rule = "coin_flip"), file.path(td, "b.csv"), row.names = FALSE); expect_error(dymid_adjudication_read(file.path(td, "b.csv")), "unknown rule")
  write.csv(ok |> mutate(field = "outcome"), file.path(td, "c.csv"), row.names = FALSE); expect_error(dymid_adjudication_read(file.path(td, "c.csv")), "does not validate")
  write.csv(bind_rows(ok, ok), file.path(td, "d.csv"), row.names = FALSE); expect_error(dymid_adjudication_read(file.path(td, "d.csv")), "duplicated")
  write.csv(ok |> select(-note), file.path(td, "e.csv"), row.names = FALSE); expect_error(dymid_adjudication_read(file.path(td, "e.csv")), "columns absent")
  expect_error(dymid_adjudication_read(file.path(td, "none.csv")), "does not exist")
  unlink(td, recursive = TRUE)
})

test_that("adjudication rule direction_uninformative: init_dyadic becomes NA, so (b) is skipped and the direction falls to (c) D1 or the fallback", {
  d <- dymid_dyads(disagree_dy("direction"), adj = adj_row(1L, 10L, 20L, "init_dyadic", "20", "10", "NA", "direction_uninformative"))
  r <- d[d$disno == 1L & d$a == 10L & d$b == 20L, ]; expect_true(is.na(r$init_dyadic)); expect_true(r$copy_adjudicated); expect_equal(r$adj_fields, "init_dyadic")
  expect_equal(sum(d$copy_adjudicated), 1L); expect_false(is.na(d$init_dyadic[d$disno == 2L]))
  e <- dymid_link(d, ddiw_dyads(toy_wr())); expect_true(is.na(e$init_b[e$disno == 1L & e$a == 10L & e$b == 20L]))
  ## through dymid_direction(): with no D1 decision the Side A fallback; with one, (c) decides a war dyad
  ev <- bind_rows(ev_d(1L, 10L, 20L, 1, 1, b = NA_integer_), ev_d(2L, 50L, 60L, 1, 1, war = TRUE, d1 = 60L, b = NA_integer_, y = 1990L))
  psd <- bind_rows(psd_d(1L, 10, 2L, 1L), psd_d(1L, 20, 2L, 1L), psd_d(2L, 50, 2L, 2L, 1990L), psd_d(2L, 60, 2L, 2L, 1990L))
  rr <- dymid_direction(ev, psd); expect_equal(rr$dir_source, c("fallback_sideA", "c_cowwar")); expect_equal(paste(rr$initiator, rr$target), c("10 20", "60 50"))
})

test_that("adjudication rule war_defer_to_war_data: the war data decides the contradicted war flag, both ways", {
  dd <- disagree_dy("war")
  d <- dymid_dyads(dd, adj = adj_row(1L, 10L, 30L, "war_dymid", "0", "1", "0", "war_defer_to_war_data"), war_pairs = war_pairs_toy)
  expect_false(d$war_dymid[d$disno == 1L & d$b == 30L]); expect_true(d$copy_adjudicated[d$disno == 1L & d$b == 30L])
  e <- dymid_link(d, ddiw_dyads(toy_wr())); expect_false(e$war[e$disno == 1L & e$a == 10L & e$b == 30L])
  ## with a war-data dyad for (1, 10, 30) the resolved value is 1 (and the table must say so)
  wp <- bind_rows(war_pairs_toy, tibble(disno = 1L, a = 10L, b = 30L))
  expect_error(dymid_dyads(dd, adj = adj_row(1L, 10L, 30L, "war_dymid", "0", "1", "0", "war_defer_to_war_data"), war_pairs = wp), "does not follow from the rule")
  d2 <- dymid_dyads(dd, adj = adj_row(1L, 10L, 30L, "war_dymid", "0", "1", "1", "war_defer_to_war_data"), war_pairs = wp); expect_true(d2$war_dymid[d2$disno == 1L & d2$b == 30L])
})

test_that("adjudication rules fatal_strict_lower / fatal_ub_higher: fmid takes the lower code, fmidany / fmidub the higher; missing sits between 0 and 1", {
  ## 13-14 in dispute 5: copy1 fatleva 1 / fatlevb 1 / fatlev 1 (strict fatal); copy2 says 14 had no deaths and the dyad was non-fatal
  tbl <- bind_rows(adj_row(5L, 13L, 14L, "fatleva_max", "1", "0", "0", "fatal_strict_lower", "fmid"), adj_row(5L, 13L, 14L, "fatleva_max", "1", "0", "1", "fatal_ub_higher", "fmidany,fmidub"),
                   adj_row(5L, 13L, 14L, "fatlev_max", "1", "0", "0", "fatal_strict_lower", "fmid"), adj_row(5L, 13L, 14L, "fatlev_max", "1", "0", "1", "fatal_ub_higher", "fmidany,fmidub"))
  d <- dymid_dyads(disagree_dy("fatal"), adj = tbl); r <- d[d$disno == 5L, ]
  expect_equal(c(r$fatleva_max, r$fatleva_ub, r$fatlev_max, r$fatlev_ub, r$fatlevb_max, r$fatlevb_ub), c(0L, 1L, 0L, 1L, 1L, 1L)); expect_equal(r$adj_fields, "fatlev_max;fatleva_max")
  lo <- dymid_fatal_codes(r$fatleva_max, r$fatlevb_max, r$fatlev_max, FALSE); hi <- dymid_fatal_codes(r$fatleva_ub, r$fatlevb_ub, r$fatlev_ub, FALSE)
  expect_false(lo$strict); expect_true(hi$any); expect_true(hi$ub)                                   # in fmidany and fmidub, not in fmid
  expect_true(lo$any); expect_false(lo$ub)                                                           # with the lower codes 14's deaths (1) still make it 'any'; not ub (no missing side)
  ## missing (-9) against a positive code: the strict sample takes the missing code (lower), the upper bound the positive one
  tbm <- bind_rows(adj_row(5L, 13L, 14L, "fatleva_max", "1", "-9", "-9", "fatal_strict_lower", "fmid"), adj_row(5L, 13L, 14L, "fatleva_max", "1", "-9", "1", "fatal_ub_higher", "fmidany,fmidub"))
  dm <- dymid_dyads(disagree_dy("fatal_missing"), adj = tbm); rm_ <- dm[dm$disno == 5L, ]
  expect_true(is.na(rm_$fatleva_max)); expect_equal(rm_$fatleva_ub, 1L)
  expect_false(dymid_fatal_codes(rm_$fatleva_max, rm_$fatlevb_max, rm_$fatlev_max, FALSE)$strict); expect_true(dymid_fatal_codes(rm_$fatleva_ub, rm_$fatlevb_ub, rm_$fatlev_ub, FALSE)$strict)
  ## missing against 0: the upper bound takes the missing code (0 < missing < 1)
  expect_equal(fat_rank(c(0L, NA, 1L, 6L)), c(0, 0.5, 1, 6))
  expect_error(dymid_dyads(disagree_dy("fatal_missing"), adj = tbm |> mutate(resolved_value = c("1", "1"))), "does not follow from the rule")
})

test_that("adjudication rules earliest_start / latest_end / hostility_higher / sides_uninformative", {
  d <- dymid_dyads(disagree_dy("start"), adj = adj_row(2L, 50L, 60L, "start_ym", "1990-02", "1990-04", "1990-02", "earliest_start"))
  expect_equal(d$start_ym[d$disno == 2L], ym_index(1990L, 2L)); expect_equal(attr(d, "spells")$spell_start_ym[attr(d, "spells")$disno == 2L], ym_index(1990L, 2L))
  expect_error(dymid_dyads(disagree_dy("start"), adj = adj_row(2L, 50L, 60L, "start_ym", "1990-02", "1990-04", "1990-04", "earliest_start")), "does not follow")
  ## latest_end resolves the spell carrying the dyad's latest end only; an end disagreement on another spell of the dyad has no row and stops (P2-2)
  er <- adj_row(5L, 13L, 14L, "end_ym", "2002-08", "2002-09", "2002-09", "latest_end")
  d <- dymid_dyads(disagree_dy("end_last"), adj = er); sp <- attr(d, "spells")[attr(d, "spells")$disno == 5L, ]
  expect_equal(d$end_ym[d$disno == 5L], ym_index(2002L, 9L)); expect_equal(sp$spell_end_ym, c(ym_index(2000L, 4L), ym_index(2002L, 9L)))
  expect_error(dymid_dyads(disagree_dy("end_both"), adj = er), "dates of a spell that is not the dyad's start or end.*dispute 5 13-14 spell 1")
  d <- dymid_dyads(disagree_dy("hostility"), adj = adj_row(2L, 50L, 60L, "hihostb", "4", "5", "5", "hostility_higher")); expect_equal(d$hihostb[d$disno == 2L], 5L)
  d <- dymid_dyads(disagree_dy("sides"), adj = bind_rows(adj_row(1L, 30L, 40L, "side_a_dy", "0", "1", "NA", "sides_uninformative"), adj_row(1L, 30L, 40L, "side_b_dy", "1", "0", "NA", "sides_uninformative")))
  y <- attr(d, "years"); expect_true(all(is.na(y$side_a_dy[y$disno == 1L & y$a == 30L]))); expect_false(anyNA(y$side_a_dy[!(y$disno == 1L & y$a == 30L)]))
  expect_equal(d$adj_fields[d$disno == 1L & d$a == 30L], "side_a_dy;side_b_dy"); expect_equal(sum(d$copy_adjudicated), 1L)
  expect_error(dymid_dyads(disagree_dy("sides"), adj = adj_row(1L, 30L, 40L, "side_a_dy", "0", "1", "NA", "sides_uninformative")), "not covered.*side_b_dy")
})

test_that("re-entry spells (2026-11-08): a dyad leaving and re-entering its dispute gives two spells; the gap year is not already fighting; one-spell dyads unchanged", {
  ## dispute 8: 17 (side A, originator) and 18 (side B, originator), Dec 2000 - Jun 2001 and Nov 2002 - Apr 2003: four consecutive participation years, ONE run
  dy <- bind_rows(toy_dy(), {
    r <- rbind(dyrow(8L, 17L, 18L, 2000L, 2000L, 12L, 2001L, 6L), dyrow(8L, 17L, 18L, 2001L, 2000L, 12L, 2001L, 6L),
               dyrow(8L, 17L, 18L, 2002L, 2002L, 11L, 2003L, 4L), dyrow(8L, 17L, 18L, 2003L, 2002L, 11L, 2003L, 4L))
    mirror_dy(r) })
  d <- dymid_dyads(dy); expect_equal(d$n_spells[d$disno == 8L], 2L); expect_equal(d$start_ym[d$disno == 8L], ym_index(2000L, 12L)); expect_equal(d$end_ym[d$disno == 8L], ym_index(2003L, 4L))
  w0 <- ddiw_dyads(toy_wr()); w0e <- w0[0, ]; attr(w0e, "years") <- attr(w0, "years")[0, ]
  e <- dymid_link(d, w0e); sp <- attr(e, "spells")[attr(e, "spells")$disno == 8L, ]
  expect_equal(sp$spell, 1:2); expect_equal(sp$spell_start_ym, c(ym_index(2000L, 12L), ym_index(2002L, 11L))); expect_equal(sp$spell_end_ym, c(ym_index(2001L, 6L), ym_index(2003L, 4L)))
  dd <- tibble(dispnum = 8L, ccode1 = c(17L, 18L), ccode2 = c(18L, 17L), sidea1 = c(1L, 0L), sidea2 = c(0L, 1L), orig1 = 1L, orig2 = 1L) |> tidyr::crossing(year = 2000:2003)
  p <- dymid_participation(dd, e[e$disno == 8L, ], attr(e, "years")[attr(e, "years")$disno == 8L, ])
  ds <- dymid_disputes(tibble(dispnum = 8L, styear = 2000L, stmon = 12L, mindur = 850, maxdur = 870, hiact = 16L, hostlev = 4L, fatality = 0L), dy |> filter(disno == 8L), e[e$disno == 8L, ], p)
  pdi <- pair_dispute_intervals(p, dispute_date_bounds(ds, p), part = NULL)
  expect_equal(nrow(pdi), 2L); expect_equal(pdi$spell_first, c(2000L, 2000L)); expect_equal(pdi$spell_last, c(2003L, 2003L))   # one bridged run per direction before the override
  expect_error(dymid_override_pdi(pdi, attr(e, "spells")), "belongs to no pair-dispute")   # spells of dyads the participation table lacks: stop
  po <- dymid_override_pdi(pdi, attr(e, "spells") |> filter(disno == 8L))
  expect_equal(nrow(po), 4L); expect_equal(attr(po, "n_split_pairs"), 2L); expect_equal(attr(po, "n_split_rows"), 2L); expect_equal(attr(po, "n_dymid"), 4L)
  s <- po[po$ccode1 == 17L, ]; expect_equal(s$spell, 1:2); expect_equal(s$spell_first, c(2000L, 2002L)); expect_equal(s$spell_last, c(2001L, 2003L)); expect_equal(s$pair_onset, c(2000L, 2000L))
  expect_equal(s$start_lo, c(ym_index(2000L, 12L), ym_index(2002L, 11L))); expect_equal(s$end_hi, c(ym_index(2001L, 6L), ym_index(2003L, 4L))); expect_true(all(s$start_exact)); expect_equal(unique(s$date_source), "dymid")
  ## already fighting at January: 2001 yes (inside spell 1), 2002 NO (the gap), 2003 yes; the bridged spell said yes in 2002
  af <- already_fighting_keys_month(po, NULL) |> filter(ccode1 == 17L) |> arrange(year); expect_equal(af$year, 2000:2003); expect_equal(af$drop, c(FALSE, TRUE, FALSE, TRUE))
  af_old <- already_fighting_keys_month(pdi, NULL) |> filter(ccode1 == 17L) |> arrange(year); expect_equal(af_old$drop, c(FALSE, TRUE, TRUE, TRUE))
  ## the first spell is the event's date; the year rule (year > pair onset) is unchanged
  expect_equal(af$drop_year, c(FALSE, TRUE, TRUE, TRUE))
  ev <- tibble(initiator = 17L, target = 18L, year = 2000L, dispnum = 8L); fe <- flag_escalations_month(ev, po); expect_equal(fe$ev_start_lo, ym_index(2000L, 12L)); expect_false(fe$escalation)
})

test_that("spec: events = pairs | dymid (suffix _dy after _ent, v2/data_dymid), samples per foundation, refusals; defaults unchanged", {
  expect_equal(SPEC_CHOICES$events, c("pairs", "dymid")); expect_equal(SPEC_SUFFIX$events, c(dymid = "_dy"))
  expect_equal(which(names(SPEC_SUFFIX) == "events"), which(names(SPEC_SUFFIX) == "entry") + 1L)
  dy <- parse_spec("--events=dymid"); expect_equal(basename(spec_data_dir(dy, root)), "data_dymid"); expect_equal(basename(spec_data_dir(parse_spec(character(0)), root)), "data")
  expect_equal(cell_name("fmid", 3L, "M1b", dy), "fmid_lag3_all_pr_M1b_dy")
  b <- parse_spec(c("--events=dymid", "--comparison=directed", "--roles=bcj", "--joiner_set=opposing"))
  expect_equal(cell_name("fmid", 3L, "M1c", b), "fmid_lag3_all_pr_M1c_cdir_bcj_jset_dy"); expect_equal(cell_name("mid", 3L, "M1cr", b), "mid_lag3_all_pr_M1cr_cdir_bcj_jset_dy")
  expect_equal(cell_name("fmid", 3L, "M2c", parse_spec(c("--events=dymid", "--comparison=directed", "--roles=bc", "--ame_K=1"))), "fmid_lag3_all_pr_M2c_K1_cdir_bc_dy")
  expect_equal(m1_cell_for_m2("fmid", 3L, "M2c", parse_spec(c("--events=dymid", "--comparison=directed", "--roles=bcj", "--ame_K=2"))), "fmid_lag3_all_pr_M1c_cdir_bcj_dy")
  expect_equal(cell_name("fmid", 3L, "M1c", parse_spec(c("--events=dymid", "--comparison=directed", "--roles=bcj", "--variant=polrel"))), "fmid_lag3_polrel_pr_M1c_cdir_bcj_dy")
  expect_equal(cell_name("fmidany", 3L, "M1b", dy), "fmidany_lag3_all_pr_M1b_dy")
  expect_error(cell_name("fmidany", 3L, "M1b", parse_spec(character(0))), "does not exist under --events=pairs"); expect_error(cell_name("fmidcow", 3L, "M1b", dy), "does not exist under --events=dymid")
  expect_error(parse_spec(c("--events=dymid", "--entry=later")), "redundant"); expect_error(parse_spec(c("--events=dymid", "--direction=sideA")), "cowwar direction only")
  expect_error(parse_spec(c("--events=dymid", "--roles=aprime")), "--roles=aprime is refused under --events=dymid", fixed = TRUE)
  expect_error(parse_spec("--events=yes"), "--events must be one of pairs | dymid", fixed = TRUE); expect_error(parse_spec("--events"), "needs a value")
  expect_silent(parse_spec(c("--events=dymid", "--roles=bc")))                             # bc / bcj need no --entry under dymid
  expect_error(parse_spec("--roles=bc"), "needs --entry=later or --events=dymid", fixed = TRUE)
  expect_false(grepl("events", spec_string(parse_spec(character(0))), fixed = TRUE)); expect_true(grepl("events=dymid", spec_string(dy), fixed = TRUE))
  expect_false(grepl("events", perm_part_spec_string(parse_spec(character(0))), fixed = TRUE)); expect_true(grepl("events=dymid", perm_part_spec_string(dy), fixed = TRUE))
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec(character(0))), "fmid_lag3_all_pr_M1b"); expect_equal(spec_events(list()), "pairs")
  expect_equal(EVENTS_SAMPLES$pairs, c("mid", "fmid", "fmidcow", "fmidub", "war")); expect_equal(EVENTS_SAMPLES$dymid, c("mid", "fmid", "fmidany", "fmidub", "war"))
})

## a toy dymid cell: the toy cell of helper-toy-cell.R with the dymid sample columns (no fmidcow); 30 a joiner on its events
dymid_edit <- function(ev) ev |> select(-in_fmidcow, -in_fmidcow_noww2, -w_fmidcow) |> mutate(in_fmidany = in_fmid, in_fmidany_noww2 = in_fmid_noww2, w_fmidany = w_fmid)
toy_roles_d <- tibble(initiator = c(10L, 10L, 30L, 30L, 30L, 40L), target = c(20L, 30L, 20L, 40L, 10L, 20L), year = 1990L, dispnum = c(1L, 4L, 2L, 3L, 6L, 7L), grp_aprime = c(0L, 1L, 2L, 0L, 0L, 2L))
toy_jsets_d <- tibble(initiator = c(30L, 40L, 40L), target = 20L, year = 1990L, dispnum = c(2L, 7L, 7L), opp_ccode = c(20L, 20L, 30L))
toy_dirstat_d <- function(a) a |> mutate(dem1 = as.integer(ccode1 %in% c(10L, 20L)), dem2 = as.integer(ccode2 %in% c(10L, 20L)),
                                         major1 = as.integer(ccode1 %in% c(10L, 40L)), major2 = as.integer(ccode2 %in% c(10L, 40L)))

test_that("assembly on a toy dymid cell: the data directory must be data_dymid, bcj without --entry, fmidcow refused", {
  base <- tempfile("dymidtoy"); td <- file.path(base, "data_dymid"); write_toy_cell(td, events_edit = dymid_edit, roles = toy_roles_d, jsets = toy_jsets_d)
  sp <- function(...) parse_spec(c("--events=dymid", ...))
  A <- quiet(assemble_cell_spec("mid", 3L, sp("--roles=bcj", "--joiner_set=opposing"), data_dir = td, model = "M1c"))
  expect_equal(A$group_table$group, BCJ_GROUPS); expect_setequal(unique(A$ev_sets$ev_idx), which(A$ev_joiner))
  expect_error(quiet(assemble_cell_spec("fmidcow", 3L, sp("--roles=bcj"), data_dir = td, model = "M1c")), "does not exist under --events=dymid")
  td2 <- file.path(base, "data"); write_toy_cell(td2, events_edit = dymid_edit, roles = toy_roles_d)
  expect_error(quiet(assemble_cell_spec("mid", 3L, sp("--roles=bcj"), data_dir = td2, model = "M1c")))
  unlink(base, recursive = TRUE)
})

test_that("02_fit_model.R end to end on a toy dymid cell: M1b, M1c bcj with joiner_set, M2c bcj from the matching M1c; refusals", {
  for (m in c("m1_context", "m2_context", "m1_latent_additive")) if (!file.exists(file.path(root, "stan", "build", paste0(m, ".rds")))) skip(paste("compiled model missing:", m))
  tmp <- tempfile("smoke_dy_"); dir.create(file.path(tmp, "v2"), recursive = TRUE)
  file.copy(file.path(root, "R"), file.path(tmp, "v2"), recursive = TRUE)
  ok <- file.symlink(file.path(root, "stan"), file.path(tmp, "v2", "stan")); if (!ok) file.copy(file.path(root, "stan"), file.path(tmp, "v2"), recursive = TRUE)
  write_toy_cell(file.path(tmp, "v2", "data_dymid"), alts_edit = toy_dirstat_d, events_edit = dymid_edit, roles = toy_roles_d, jsets = toy_jsets_d)
  dir.create(file.path(tmp, "v2", "output")); dir.create(file.path(tmp, "v2", "logs"))
  v2 <- file.path(tmp, "v2"); out_dir <- file.path(v2, "output")
  run <- function(script, ...) { out <- suppressWarnings(system2("Rscript", c(file.path(v2, "R", script), ...), stdout = TRUE, stderr = TRUE)); list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  base <- c("--sample=mid", "--lag=3", "--comparison=directed", "--events=dymid"); short <- c("--chains=2", "--iter=60", "--warmup=30", "--max_treedepth=6")
  expect_true(run("02_fit_model.R", base, "--roles=bcj", "--model=M1c", "--entry=later", "--dry_run=TRUE")$status != 0L)
  expect_true(run("02_fit_model.R", base, "--roles=aprime", "--model=M1c", "--dry_run=TRUE")$status != 0L)
  expect_true(run("02_fit_model.R", "--sample=fmidcow", "--lag=3", "--model=M1b", "--events=dymid", "--dry_run=TRUE")$status != 0L)
  r <- run("02_fit_model.R", base, "--model=M1b", short); expect_equal(r$status, 0L, info = paste(tail(r$out, 20), collapse = "\n"))
  o <- readRDS(file.path(out_dir, "fit_mid_lag3_all_pr_M1b_cdir_dy.rds")); expect_equal(o$data_dir, "data_dymid"); expect_equal(o$spec$events, "dymid")
  r <- run("02_fit_model.R", base, "--model=M1c", "--roles=bcj", "--joiner_set=opposing", short); expect_equal(r$status, 0L, info = paste(tail(r$out, 30), collapse = "\n"))
  o <- readRDS(file.path(out_dir, "fit_mid_lag3_all_pr_M1c_cdir_bcj_jset_dy.rds")); expect_equal(o$dat$G, 3L); expect_equal(o$group_labels, BCJ_GROUPS)
  r <- run("02_fit_model.R", base, "--model=M2c", "--roles=bcj", "--ame_K=1", "--joiner_set=opposing", short); expect_equal(r$status, 0L, info = paste(tail(r$out, 30), collapse = "\n"))
  o2 <- readRDS(file.path(out_dir, "fit_mid_lag3_all_pr_M2c_K1_cdir_bcj_jset_dy.rds")); expect_equal(o2$init$m1_cell, "mid_lag3_all_pr_M1c_cdir_bcj_jset_dy")
  dr <- run("05_robustness.R", "--preset=noww2", "--events=dymid", "--dry_run=TRUE"); expect_true(dr$status != 0L)
  unlink(tmp, recursive = TRUE)
})

test_that("real dymid build (when built): the PI's cases, counts per sample and bcj group, the stage-2 sets hold every chosen target, joiner sets", {
  dir <- file.path(v2_root(), "data_dymid")
  skip_if_not(file.exists(file.path(dir, "manifest.csv")) && file.exists(file.path(dir, "alts_lag3_all.rds")), "dymid build (v2/data_dymid) not present")
  ev <- readRDS(file.path(dir, "events.rds")); ne <- ev |> filter(!escalation)
  expect_true(any(ne$initiator == 255L & ne$target == 290L & ne$year == 1939L)); expect_false(any(ev$initiator == 290L & ev$target == 255L & ev$year == 1939L))
  expect_false(any((ev$initiator == 740L & ev$target == 290L) | (ev$initiator == 290L & ev$target == 740L)))
  expect_true(all(c("in_fmidany", "w_fmidany", "fatal_any", "ev_source") %in% names(ev))); expect_false("in_fmidcow" %in% names(ev))
  expect_true(all(ne$in_fmidany[ne$in_fmid])); expect_true(all(ne$in_fmidub[ne$in_fmid])); expect_true(all(ne$in_fmid[ne$in_war]))
  expect_false(any(ne$orig1 == 1L & ne$orig2 == 0L))
  ## fmidub (2026-11-05): a dyad coded non-fatal by its dyadic fatlev (0) is in fmidub only through the strict rule or war
  expect_false(any(ev$fatal_ub & !ev$fatal_main & ev$fatality_disp %in% 0L))
  ## README 6.9 examples: Morocco -> Iraq 1990 (dispute 3957) is a mid event by rule (a), not fmid; the war-only dyad 816 -> 811 1970 (1617) is an escalation
  mi <- ev |> filter(initiator == 600L, target == 645L, year == 1990L); expect_equal(nrow(mi), 1L); expect_true(mi$in_mid); expect_false(mi$in_fmid); expect_equal(mi$dir_source, "a_joiner_orig")
  wo <- ev |> filter(ev_source == "ddiw_only"); expect_equal(paste(wo$dispnum, wo$initiator, wo$target, wo$year), "1617 816 811 1970"); expect_true(wo$escalation)
  ## copy disagreements (2026-11-08): the adjudication table is checked in and covers the real data; the PI's examples
  adj <- dymid_adjudication_read(file.path(v2_root(), "data_raw", "dymid", "adjudication.csv")); expect_equal(nrow(adj), 75L); expect_equal(nrow(distinct(adj, disno, a, b)), 29L)
  expect_true(any(adj$disno == 2970L & adj$field == "init_dyadic")); expect_true(any(adj$disno == 4339L & adj$field == "war_dymid")); expect_true(any(adj$disno == 1353L & adj$a == 365L & adj$field == "fatlev_max"))
  nd <- ev |> filter(dispnum == 2970L, year == 1981L); expect_equal(nrow(nd), 1L); expect_equal(nd$dir_source, "fallback_sideA"); expect_equal(c(nd$initiator, nd$target), c(385L, 390L))   # Norway - Denmark: (b) uninformative, Side A kept
  cr <- ev |> filter(dispnum == 4339L, (initiator == 483L & target == 517L) | (initiator == 517L & target == 483L)); expect_equal(nrow(cr), 1L); expect_false(cr$war); expect_false(cr$in_war)   # Chad - Rwanda: no war-data dyad
  rl <- ev |> filter(dispnum == 1353L, initiator == 365L, target == 812L, year == 1962L); expect_equal(nrow(rl), 1L); expect_true(rl$in_fmidany); expect_false(rl$in_fmid); expect_false(rl$in_fmidub)   # Russia - Laos: fatlev 0 (fmid) / 6 (fmidany)
  gs <- ev |> filter(dispnum == 3833L, year == 1939L); expect_equal(nrow(gs), 1L); expect_equal(c(gs$initiator, gs$target), c(165L, 200L)); expect_equal(gs$dir_source, "fallback_sideA")   # Uruguay -> UK (Graf Spee): the A < B copy said 200 -> 165
  si <- ev |> filter(dispnum == 4182L, initiator == 652L, target == 666L, year == 1993L); expect_equal(nrow(si), 1L); expect_false(si$in_fmid); expect_true(si$in_fmidany); expect_true(si$in_fmidub)   # Syria -> Israel 1993: lower fatlevb 0
  expect_false(any(ne$in_fmid & ne$dispnum == 4182L & ne$year == 1993L & ((ne$initiator == 652L & ne$target == 666L) | (ne$initiator == 666L & ne$target == 652L))))
  expect_equal(c(sum(ne$in_mid), sum(ne$in_fmid), sum(ne$in_fmidany), sum(ne$in_fmidub), sum(ne$in_war)), c(3208L, 532L, 771L, 549L, 278L))   # post-collapse counts of 2026-11-08
  ## re-entry spells (2026-11-08): USA - Bulgaria (258) is not already fighting in 1943 (the gap Jun 1942 - Nov 1943), so both directions are lag-3 alternatives in both sets
  al <- readRDS(file.path(dir, "alts_lag3_all.rds")); ub <- al |> filter((ccode1 == 2L & ccode2 == 355L) | (ccode1 == 355L & ccode2 == 2L), year == 1943L)
  expect_equal(nrow(ub), 2L); expect_true(all(ub$in_s2 == 1L)); expect_true(all(ub$in_iv == 1L))
  ub42 <- al |> filter(ccode1 == 2L, ccode2 == 355L, year == 1942L); expect_true(nrow(ub42) == 0L || ub42$in_iv == 0L)   # inside the first spell at January 1942: already fighting
  ro <- readRDS(file.path(dir, "roles_aprime.rds")); expect_false(any(ro$grp_aprime == 3L))
  pl <- readRDS(file.path(dir, "stage1_panel_lag3.rds")); expect_true(all(paste0("y1_", EVENTS_SAMPLES$dymid) %in% names(pl)))
  for (s in EVENTS_SAMPLES$dymid) {
    jr <- ne |> filter(.data[[paste0("in_", s)]]) |> inner_join(ro, by = c("initiator", "target", "year"))
    message(sprintf("dymid build, sample %-7s: %4d events; bcj bilateral %d / coalition originator %d / joiner %d", s, sum(ne[[paste0("in_", s)]]), sum(jr$grp_aprime == 0L), sum(jr$grp_aprime == 1L), sum(jr$grp_aprime == 2L)))
  }
  for (s in c("fmid", "mid")) {
    A <- quiet(assemble_cell_spec(s, 3L, parse_spec(c("--comparison=directed", "--events=dymid", "--roles=bcj", "--joiner_set=opposing")), model = "M1c"))
    expect_true(all(A$alts$in_s2[A$ev$alt_row] == 1L))                                       # every usable event's chosen target is a stage-2 row of its state-year
    expect_true(all(A$alts$sy_row[A$ev$alt_row] == A$ev$sy_row))                             # ... of its own state-year block
    expect_setequal(unique(A$ev_sets$ev_idx), which(A$ev_joiner)); expect_true(all(paste(which(A$ev_joiner), A$ev$alt_row[A$ev_joiner]) %in% paste(A$ev_sets$ev_idx, A$ev_sets$alt_row)))
    expect_equal(nrow(A$group_table), 3L); expect_true(all(A$group_table$n_usable > 0L))
    message(sprintf("dymid build %s bcj jset: %s; joiner events %d, size 1 %d, removed %d; E = %d, sum(y1) = %d", s, paste(sprintf("%s %d (%.1f)", A$group_table$group, A$group_table$n_usable, A$group_table$w_usable), collapse = ", "),
                    A$jset_stats$n_joiner_events, A$jset_stats$n_singleton, A$jset_stats$n_removed, nrow(A$ev), sum(A$sy$y1)))
    rm(A); invisible(gc())
  }
  Ap <- quiet(assemble_cell_spec("fmid", 3L, parse_spec(c("--comparison=directed", "--events=dymid", "--roles=bcj", "--variant=polrel")), model = "M1c"))
  expect_true(nrow(Ap$ev) > 0L); rm(Ap)
  Ab <- quiet(assemble_cell_spec("fmid", 3L, parse_spec("--events=dymid"), model = "M1b")); expect_true(nrow(Ab$ev) > 0L)
  message(sprintf("dymid build fmid M1b main spec: E = %d, sum(y1) = %d, weights %.1f", nrow(Ab$ev), sum(Ab$sy$y1), sum(Ab$ev$w)))
})
