## The interacting-dyads event foundation (switch events = pairs | dymid, suffix _dy, build variant v2/data_dymid; PI-approved 2026-11-04; README 6.9).
## (a) helpers_dymid.R on hand-built DYMID and war rows: the dyad key, max-over-years fatality, re-entries, linking (merged / war-only / DYMID-only
## war), recip validation (2 accepted explicitly as not reciprocated, other values refused), the samples strict / any / ub (definition of
## 2026-11-05), the direction rules in the order (a) / (b) / (c) / fallback with ties (PI decision 2026-11-05), the participation table and the spell override;
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
  ## mirror copies (A > B), as the file lists them
  m <- r; m$statea <- r$stateb; m$stateb <- r$statea; m$fatleva <- r$fatlevb; m$fatlevb <- r$fatleva; m$sideaa <- r$sideab; m$sideab <- r$sideaa
  m$orignata <- r$orignatb; m$orignatb <- r$orignata; m$dyad_rolea <- r$dyad_roleb; m$dyad_roleb <- r$dyad_rolea
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
  d <- dymid_dyads(bind_rows(toy_dy() |> filter(disno == 1L), dyrow(7L, 15L, 16L, 2000L, 2000L, 4L, 2000L, 5L, sa = 0L))); e <- dymid_link(d, w0e)
  p <- dymid_participation(dd, e, attr(e, "years")); dg <- attr(p, "diagnostics")
  ## COW pairs 10-20, 10-30, 10-99, 20-40, 30-40, 40-99: the four event dyads replaced, 10-99 and 40-99 (no interaction) kept as COW has them
  expect_equal(dg$n_cow_pairs, 6L); expect_equal(dg$n_cow_pairs_replaced, 4L); expect_equal(dg$n_cow_pairs_not_event, 2L); expect_equal(dg$n_event_pairs_added, 1L); expect_equal(dg$n_disputes_added, 1L)
  expect_true(all(c("10 99", "40 99") %in% paste(p$ccode1, p$ccode2)))
  expect_equal(p$sidea1[p$dispnum == 7L & p$ccode1 == 15L], 0L)                          # not a COW participant: DYMID's side (sideaa 0)
  expect_true(all(p$sidea1 + p$sidea2 == 1L))
  expect_equal(unique(p$orig1[p$ccode1 == 40L]), 0L); expect_equal(unique(p$sidea1[p$ccode1 == 40L & p$year == 2000L]), 1L)
  expect_equal(sum(p$ccode1 == 10L & p$ccode2 == 20L), 2L)                                 # the event dyad's own years (2000, 2001), both directions
  ds <- dymid_disputes(tibble(dispnum = 1L, styear = 2000L, stmon = 3L, mindur = 400, maxdur = 460, hiact = 18L, hostlev = 4L, fatality = 3L),
                       bind_rows(toy_dy() |> filter(disno == 1L), dyrow(7L, 15L, 16L, 2000L, 2000L, 4L, 2000L, 5L, sa = 0L)), e, p)
  expect_equal(ds$hostlev_ev[ds$dispnum == 1L], 5L)                                        # 20 - 40 is a war dyad
  db <- dispute_date_bounds(ds, p); pdi <- pair_dispute_intervals(p, db, part = NULL)
  po <- dymid_override_pdi(pdi, attr(e, "years"))
  s <- po[po$ccode1 == 10L & po$ccode2 == 30L, ]; expect_equal(c(s$start_lo, s$end_hi), c(ym_index(2000L, 5L), ym_index(2000L, 9L))); expect_equal(s$date_source, "dymid")
  expect_equal(po$date_source[po$ccode1 == 10L & po$ccode2 == 99L], "dispute")            # a co-participant pair keeps the rule
  expect_equal(attr(po, "n_dymid"), 2L * 5L)                                               # 5 event dyads, both directions, one spell each
  expect_equal(ds$hostlev_ev[ds$dispnum == 7L], 4L); expect_equal(nrow(attr(ds, "constructed")), 1L)
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
