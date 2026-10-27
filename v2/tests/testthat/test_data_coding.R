## Unit tests of the data-coding logic in 01_build_data.R using small hand-built toy frames
## that mimic the peacesciencer columns. The full build needs the bundled peacesciencer data
## and is run by Rscript v2/R/01_build_data.R (on Roar or locally); these tests exercise the
## rules in isolation so a regression in the coding logic fails loudly without the full build.
suppressPackageStartupMessages({ library(testthat); library(dplyr); library(tidyr) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "spec.R")); source(file.path(root, "R", "helpers_data.R")); source(file.path(root, "R", "helpers_model.R"))

## ---- the coding rules under test, lifted verbatim from 01_build_data.R ---------------------
## disp needs the dispute-level columns fatality_disp, hostlev_disp and (month rule, 2026-10-24) styear, stmon, mindur, maxdur
## war_inter: COW inter-state war dyad-years for the direction rule (D1, 2026-10-25); the default is an empty table, so every
## event keeps the Side A rule (dir_source "sideA"), as before. part: participant-level dates (D5), NULL = dispute-level bounds.
empty_war_inter <- tibble(warnum = integer(), ccode1 = integer(), ccode2 = integer(), year = integer(), cowinterongoing = integer(), cowinteronset = integer(),
                          sidea1 = integer(), sidea2 = integer(), initiator1 = integer(), initiator2 = integer())
code_events <- function(dd, disp, include_joiners = FALSE,
                        mic = tibble(dispnum = integer(), mic_fatalmin = integer()),
                        gml = tibble(dispnum = integer(), initiator = integer(), target = integer(), fatal_gml = logical()),
                        war_inter = empty_war_inter, part = NULL, direction = "cowwar") {
  pair <- build_pair_table(dd, disp)                     # the build's pair table (earliest year's sides; P3 fix 2026-10-25)
  events <- pair |> filter(sidea1 == 1L, sidea2 == 0L)
  if (!include_joiners) events <- events |> filter(orig1 == 1L, orig2 == 1L)
  events <- events |> rename(initiator = ccode1, target = ccode2)
  ## D1: direction follows COW's war-initiating side in linked war-level disputes (window t..t+1, min 2 linked events, as the build)
  events <- cow_war_direction(events, war_dyad_years(war_inter), window = 1L, min_link = 2L)
  if (direction == "sideA") events <- events |> mutate(dir_flipped = FALSE, dir_source = "sideA", dir_rule = "sideA_build", war_link = NA_integer_)
  events <- apply_direction(events) |>
    mutate(is_joiner_pair = orig1 == 0L | orig2 == 0L,
           fatal1 = ifelse(fatality1 < 0, NA, fatality1 >= 1),
           fatal2 = ifelse(fatality2 < 0, NA, fatality2 >= 1),
           fatal_disp = ifelse(fatality_disp < 0, NA, fatality_disp >= 1),
           fatal = case_when(fatal1 %in% TRUE | fatal2 %in% TRUE | fatal_disp %in% TRUE ~ TRUE,
                             is.na(fatal1) & is.na(fatal2) & is.na(fatal_disp) ~ NA,
                             TRUE ~ FALSE),
           war_disp = hostlev_disp == 5L,
           war = hostlev_disp == 5L & (hostlev1 == 5L | hostlev2 == 5L))      # pair-level war (2026-10-24, at least one side at 5)
  ## ongoing-dispute rule at month resolution (2026-10-24): escalation = already fighting in another dispute at the
  ## event's dispute start date (year rule where undecided); excluded from every sample flag
  pdi <- pair_dispute_intervals(dd, dispute_date_bounds(disp, dd), part = part)
  events <- flag_escalations_month(events, pdi)
  events |>
    left_join(mic, by = "dispnum") |> left_join(gml, by = c("dispnum", "initiator", "target")) |>
    mutate(fatal_mic = ifelse(is.na(mic_fatalmin), NA, mic_fatalmin >= 1),
           fatal_main = code_fatal_main(fatal, fatal_mic, war, fatal_gml),
           in_mid = !escalation,
           in_fmid = fatal_main & !escalation,
           in_fmidub = (fatal_main | (is.na(fatal) & is.na(fatal_mic) & is.na(fatal_gml))) & !escalation,
           in_fmidcow = (fatal %in% TRUE) & !escalation,
           in_war = (war %in% TRUE) & !escalation) |>
    event_weights()                                        # D8: w_<sample> on the pre-collapse rows
}
## war_dyad_years() needs a complete table; the empty one yields no match
war_dyad_years_empty <- function() war_dyad_years(empty_war_inter)
## dispute-level toy rows with dates; mindur = maxdur means every COW day is known (the end stays bounded because
## peacesciencer drops the start day)
mk_disp <- function(dispnum, fatality_disp, hostlev_disp, styear, stmon, dur)
  tibble(dispnum = dispnum, fatality_disp = fatality_disp, hostlev_disp = hostlev_disp, styear = styear, stmon = stmon, mindur = dur, maxdur = dur)

## toy directed dispute-years: dispute 1 is bilateral A=10 vs B=20 (both originators), 1950-1951;
## dispute 2 is multilateral: side A originators 10, 11; side B originator 20, side B joiner 21 (orig = 0)
toy_dd <- tribble(
  ~dispnum, ~ccode1, ~ccode2, ~year, ~disponset, ~dispongoing, ~sidea1, ~sidea2, ~orig1, ~orig2, ~fatality1, ~fatality2, ~hostlev1, ~hostlev2,
  1L, 10L, 20L, 1950L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, 0L, 4L, 4L,
  1L, 10L, 20L, 1951L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 0L, 4L, 4L,
  1L, 20L, 10L, 1950L, 1L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 4L, 4L,
  1L, 20L, 10L, 1951L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 4L, 4L,
  2L, 10L, 20L, 1960L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, -9L, 5L, 5L,
  2L, 10L, 21L, 1960L, 1L, 1L, 1L, 0L, 1L, 0L, 1L, 0L, 5L, 3L,
  2L, 11L, 20L, 1960L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, -9L, 3L, 5L,
  2L, 11L, 21L, 1960L, 1L, 1L, 1L, 0L, 1L, 0L, 0L, 0L, 3L, 3L,
  2L, 20L, 10L, 1960L, 1L, 1L, 0L, 1L, 1L, 1L, -9L, 1L, 5L, 5L,
  2L, 20L, 11L, 1960L, 1L, 1L, 0L, 1L, 1L, 1L, -9L, 0L, 5L, 3L,
  2L, 21L, 10L, 1960L, 1L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 3L, 5L,
  2L, 21L, 11L, 1960L, 1L, 1L, 0L, 1L, 0L, 1L, 0L, 0L, 3L, 3L)
toy_disp <- bind_rows(mk_disp(1L, 0L, 4L, 1950L, 3L, 400), mk_disp(2L, 1L, 5L, 1960L, 7L, 30))

test_that("F4: initiation is Side A originator -> Side B originator, one direction per dispute pair", {
  ev <- code_events(toy_dd, toy_disp)
  expect_equal(nrow(ev), 3)                                     # (10->20, d1), (10->20, d2), (11->20, d2)
  expect_true(all(ev$initiator %in% c(10L, 11L)))
  expect_true(all(ev$target == 20L))                             # 21 is a joiner -> excluded
  expect_false(any(paste(ev$dispnum, ev$target, ev$initiator) %in% paste(ev$dispnum, ev$initiator, ev$target)))
  expect_equal(ev$year[ev$dispnum == 1L], 1950L)                # onset year, not every dispute-year
  ## multilateral: each A-originator x B-originator pair is an event
  expect_equal(sort(ev$initiator[ev$dispnum == 2L]), c(10L, 11L))
})

test_that("joiners (main build since 2026-10-24): every Side A x Side B pair is an event, is_joiner_pair marks pairs with a joiner", {
  ev <- code_events(toy_dd, toy_disp, include_joiners = TRUE)
  expect_equal(nrow(ev), 5)                                     # d1: 10->20; d2: 10->20, 10->21, 11->20, 11->21
  expect_true(21L %in% ev$target)
  expect_equal(sum(ev$is_joiner_pair), 2L)                      # the two pairs with the Side B joiner 21
  expect_equal(sort(ev$initiator[ev$is_joiner_pair]), c(10L, 11L))
  expect_false(any(ev$is_joiner_pair[ev$target == 20L]))
  ## one direction per dispute pair still holds with joiners, and the originator-only set is exactly the old event base
  expect_false(any(paste(ev$dispnum, ev$target, ev$initiator) %in% paste(ev$dispnum, ev$initiator, ev$target)))
  old <- code_events(toy_dd, toy_disp, include_joiners = FALSE)
  expect_equal(nrow(old), sum(!ev$is_joiner_pair))
  ## the joiner pairs take the dispute-level fatality like every other pair; war is pair-level (2026-10-24, at least one
  ## participant at hostlev 5): in dispute 2 every pair involving 10 or 20 (hostlev 5) is a war pair; a pair of 11 and 21
  ## (both hostlev 3) would be war_disp but not war
  expect_true(all(ev$in_fmid[ev$dispnum == 2L])); expect_true(all(ev$war_disp[ev$dispnum == 2L]))
  d2 <- ev$dispnum == 2L
  expect_equal(ev$in_war[d2], ev$initiator[d2] %in% c(10L, 20L) | ev$target[d2] %in% c(10L, 20L))
})

test_that("fatal / war samples: negative fatality codes -> NA, fatal if either side >= 1; war = hostlev 5", {
  ev <- code_events(toy_dd, toy_disp)
  d1 <- ev[ev$dispnum == 1L, ]; d2 <- ev[ev$dispnum == 2L & ev$initiator == 10L, ]; d2b <- ev[ev$dispnum == 2L & ev$initiator == 11L, ]
  expect_false(d1$in_fmid); expect_false(d1$in_war)
  expect_true(is.na(d2$fatal2)); expect_true(d2$in_fmid); expect_true(d2$in_war)
  expect_true(d2b$in_fmid)       # side 1 fatality 0, side 2 NA, dispute-level fatality 1 -> fatal
  ## all-NA fatality -> NA, not FALSE
  dd3 <- toy_dd |> filter(dispnum == 1L) |> mutate(fatality1 = -9L, fatality2 = -9L)
  ev3 <- code_events(dd3, mk_disp(1L, -9L, 4L, 1950L, 3L, 400))
  expect_true(is.na(ev3$fatal)); expect_false(ev3$in_fmid); expect_false(ev3$in_fmidcow)
})

test_that("fatal coding (2026-10-23): wars are fatal; MIC fills COW-missing only; fmidcow = COW only", {
  ## truth table of the helper
  cow <- c(TRUE, FALSE, NA, NA, NA, FALSE, NA)
  mic <- c(FALSE, TRUE, TRUE, FALSE, NA, NA, NA)
  war <- c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)
  gml <- c(FALSE, TRUE, FALSE, TRUE, NA, NA, NA)
  expect_identical(code_fatal_main(cow, mic, war, gml), c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE))
  ## GML fills only when COW and MIC are both missing
  expect_identical(code_fatal_main(c(NA, NA, NA, FALSE), c(NA, NA, FALSE, NA), rep(FALSE, 4), c(TRUE, FALSE, TRUE, TRUE)),
                   c(TRUE, FALSE, FALSE, FALSE))
  expect_error(code_fatal_main(cow, mic[-1], war, gml))
  ## all-NA COW fatality, not a war: MIC decides
  dd3 <- toy_dd |> filter(dispnum == 1L) |> mutate(fatality1 = -9L, fatality2 = -9L)
  disp3 <- mk_disp(1L, -9L, 4L, 1950L, 3L, 400)
  expect_true(code_events(dd3, disp3, mic = tibble(dispnum = 1L, mic_fatalmin = 30L))$in_fmid)
  expect_false(code_events(dd3, disp3, mic = tibble(dispnum = 1L, mic_fatalmin = 0L))$in_fmid)
  expect_false(code_events(dd3, disp3, mic = tibble(dispnum = 1L, mic_fatalmin = 30L))$in_fmidcow)
  ## all-NA COW fatality, a war-level dispute, not in MIC: fatal by the war rule ONLY when at least one side of the pair is
  ## at hostlev 5 (pair-level war, 2026-10-24); COW-only coding drops it either way
  ev_w <- code_events(dd3 |> mutate(hostlev1 = 5L, hostlev2 = 5L), mk_disp(1L, -9L, 5L, 1950L, 3L, 400))
  expect_true(ev_w$in_war); expect_true(ev_w$in_fmid); expect_false(ev_w$in_fmidcow)
  ev_w2 <- code_events(dd3 |> mutate(hostlev1 = 5L, hostlev2 = 4L), mk_disp(1L, -9L, 5L, 1950L, 3L, 400))   # one side at 5
  expect_true(ev_w2$in_war); expect_true(ev_w2$war_disp); expect_true(ev_w2$in_fmid); expect_equal(ev_w2$fatal_main, TRUE)
  ev_w3 <- code_events(dd3 |> mutate(hostlev1 = 4L, hostlev2 = 4L), mk_disp(1L, -9L, 5L, 1950L, 3L, 400))   # neither side at 5
  expect_false(ev_w3$in_war); expect_true(ev_w3$war_disp); expect_false(ev_w3$in_fmid); expect_equal(ev_w3$fatal_main, FALSE)
  ## MIC never overrides a COW code: dispute 1 is COW non-fatal, MIC says fatal -> stays non-fatal
  ev_o <- code_events(toy_dd, toy_disp, mic = tibble(dispnum = c(1L, 2L), mic_fatalmin = c(5L, 0L)))
  expect_false(ev_o$in_fmid[ev_o$dispnum == 1L])
  expect_true(all(ev_o$in_fmid[ev_o$dispnum == 2L]))   # COW fatal, MIC 0 -> stays fatal
  expect_true(all(ev_o$in_fmid[ev_o$in_war])); expect_true(all(ev_o$in_fmid[ev_o$in_fmidcow]))
  ## GML third: all-missing COW, not in MIC -> GML decides; upper bound counts only events no source codes
  g <- function(f) tibble(dispnum = 1L, initiator = 10L, target = 20L, fatal_gml = f)
  ev_g <- code_events(dd3, disp3, gml = g(TRUE)); expect_true(ev_g$in_fmid); expect_true(ev_g$in_fmidub)
  ev_g0 <- code_events(dd3, disp3, gml = g(FALSE)); expect_false(ev_g0$in_fmid); expect_false(ev_g0$in_fmidub)
  ev_none <- code_events(dd3, disp3); expect_false(ev_none$in_fmid); expect_true(ev_none$in_fmidub)
  ## MIC beats GML when both present
  expect_false(code_events(dd3, disp3, mic = tibble(dispnum = 1L, mic_fatalmin = 0L), gml = g(TRUE))$in_fmid)
})

test_that("F5: calendar-year lags via joins, never row lag()", {
  ## state 10 has positions in 1990, 1991, 1993 (1992 missing). Row lag would pair 1993 with 1991.
  pos <- tribble(~ccode, ~year, ~M, 10L, 1990L, 0.0, 10L, 1991L, 0.5, 10L, 1993L, 2.0)
  sy <- tibble(ccode = 10L, year = c(1992L, 1994L))
  h <- 1L
  panel <- sy |>
    inner_join(pos |> transmute(ccode, year = year + 1L, M_tm1 = M), by = c("ccode", "year")) |>
    inner_join(pos |> transmute(ccode, year = year + 1L + h, M_w0 = M), by = c("ccode", "year")) |>
    mutate(dM_h = M_tm1 - M_w0)
  expect_equal(nrow(panel), 1)                 # 1994 needs 1993 and 1992; 1992 is missing -> dropped, not bridged
  expect_equal(panel$year, 1992L); expect_equal(panel$dM_h, 0.5)
  ## row-lag version would have produced a (wrong) 1994 row: confirm the behaviours differ
  rowlag <- pos |> arrange(year) |> mutate(dM_row = M - lag(M))
  expect_equal(rowlag$dM_row[3], 1.5)          # the artefact v2 avoids
})

test_that("conflict-history windows: counts in [t-1-h-10, t-1-h-1] and peace years at t-1-h", {
  ## monadic onsets for state 10 in years 1970, 1975, 1979, 1980, 1985, 1988; in system 1960-2000
  grid <- tibble(ccode = 10L, year = 1960:2000) |>
    mutate(onsets = as.integer(year %in% c(1970L, 1975L, 1979L, 1980L, 1985L, 1988L)), ongoing = onsets) |>
    arrange(year) |>
    mutate(cum_onsets = cumsum(onsets),
           last_mid_year = cummax(ifelse(ongoing == 1L, year, min(year) - 1L)),
           mon_peace_years = year - coalesce(lag(last_mid_year), min(year) - 1L))
  ## t = 1991, h = 10: window [1980, 1990]; pre-window = [1970, 1979] -> onsets 1970, 1975, 1979 = 3
  t <- 1991L; h <- 10L; win_start <- t - 1L - h
  cum_at <- function(y) grid$cum_onsets[grid$year == y]
  expect_equal(win_start, 1980L)
  expect_equal(cum_at(win_start - 1L) - cum_at(win_start - 1L - 10L), 3L)
  ## peace years AS OF 1980 = years since last MID strictly before 1980 (1979) = 1
  expect_equal(grid$mon_peace_years[grid$year == 1980L], 1L)
  ## as of 1985: last before is 1980 -> 5; as of 1961 (no MID yet): 1961 - (1960 - 1) = 2
  expect_equal(grid$mon_peace_years[grid$year == 1985L], 5L)
  expect_equal(grid$mon_peace_years[grid$year == 1961L], 2L)
  ## windows for h = 1 vs h = 10 differ (controls are lag-specific)
  h1 <- t - 1L - 1L                                             # win_start for h = 1 is 1989
  expect_equal(cum_at(h1 - 1L) - cum_at(h1 - 1L - 10L), 4L)   # [1979, 1988] -> 1979, 1980, 1985, 1988
})

test_that("choice set: all other system members at t with positions at t-1; chosen target inside", {
  states <- tibble(ccode = c(10L, 20L, 30L, 40L), year = 1991L)
  pos <- tibble(ccode = c(10L, 20L, 30L), year = 1990L, M = c(1, 0, -1), R = c(0.5, -0.5, 0))   # 40 entered in 1991: no t-1 position
  panel <- tibble(ccode = 10L, year = 1991L, M_i = 1, R_i = 0.5)
  alts <- panel |> inner_join(states |> rename(ccode2 = ccode), by = "year", relationship = "many-to-many") |>
    filter(ccode != ccode2) |>
    inner_join(pos |> transmute(ccode2 = ccode, year = year + 1L, M_j = M, R_j = R), by = c("ccode2", "year")) |>
    mutate(dM = M_i - M_j, dR = R_i - R_j)
  expect_equal(sort(alts$ccode2), c(20L, 30L))
  expect_equal(alts$dM[alts$ccode2 == 30L], 2)
  expect_equal(alts$dR[alts$ccode2 == 20L], 1)
  ## event 10 -> 20 is inside; event 10 -> 40 would not be and must be dropped with a log line, not silently kept
  ev <- tibble(ccode1 = 10L, ccode2 = c(20L, 40L), year = 1991L)
  inside <- ev |> semi_join(alts, by = c("ccode1" = "ccode", "ccode2", "year"))
  expect_equal(inside$ccode2, 20L)
})

test_that("politically relevant choice set keeps contiguous or major-power dyads", {
  a <- tibble(contig_any = c(1L, 0L, 0L, 0L), major1 = c(0L, 1L, 0L, 0L), major2 = c(0L, 0L, 1L, 0L)) |>
    mutate(polrel = as.integer(contig_any == 1L | major1 == 1L | major2 == 1L))
  expect_equal(a$polrel, c(1L, 1L, 1L, 0L))
})

test_that("within-year standardisation invariant is asserted", {
  df <- tibble(year = rep(1:2, each = 3), v = c(1, 2, 3, 10, 20, 60)) |> group_by(year) |> mutate(z = as.numeric(scale(v))) |> ungroup()
  expect_true(assert_standardized(df, "z"))
  expect_error(assert_standardized(df, "v"))
})

test_that("manifest assertion catches a changed file", {
  tmp <- tempfile(fileext = ".rds"); saveRDS(data.frame(a = 1:3), tmp)
  man <- tempfile(fileext = ".csv")
  write_manifest(list(data.frame(file = basename(tmp), nrow = 3, ncol = 1, sha256 = file_digest(tmp), note = "t")), man)
  expect_equal(nrow(assert_matches_manifest(tmp, man)), 3)
  saveRDS(data.frame(a = 1:4), tmp)
  expect_error(assert_matches_manifest(tmp, man))
})

## ---- conflict-history counting: one count per state-dispute (and dyad-dispute) at the first onset year
count_mon <- function(dd) {
  mon_first <- dd |> filter(orig1 == 1L) |> group_by(ccode = ccode1, dispnum) |> summarise(year = min(year), .groups = "drop") |>
    count(ccode, year, name = "onsets")
  mon_first
}
count_dy <- function(dd) {
  dd |> filter(orig1 == 1L & orig2 == 1L) |> group_by(ccode1, ccode2, dispnum) |> summarise(year = min(year), .groups = "drop") |>
    count(ccode1, ccode2, year, name = "onsets")
}
test_that("a joiner entering a dispute in a later year does not create a second onset for the originator", {
  ## dispute 50: 10 vs 20 starts 1953; 30 joins side B in 1954 (disponset = 1 for the new pair 10-30 in 1954)
  dd <- tribble(
    ~dispnum, ~ccode1, ~ccode2, ~year, ~disponset, ~dispongoing, ~sidea1, ~sidea2, ~orig1, ~orig2,
    50L, 10L, 20L, 1953L, 1L, 1L, 1L, 0L, 1L, 1L,
    50L, 10L, 20L, 1954L, 0L, 1L, 1L, 0L, 1L, 1L,
    50L, 10L, 30L, 1954L, 1L, 1L, 1L, 0L, 1L, 0L,
    50L, 20L, 10L, 1953L, 1L, 1L, 0L, 1L, 1L, 1L,
    50L, 20L, 10L, 1954L, 0L, 1L, 0L, 1L, 1L, 1L,
    50L, 30L, 10L, 1954L, 1L, 1L, 0L, 1L, 0L, 1L,
    51L, 10L, 20L, 1954L, 1L, 1L, 1L, 0L, 1L, 1L,
    51L, 20L, 10L, 1954L, 1L, 1L, 0L, 1L, 1L, 1L)
  m <- count_mon(dd)
  expect_equal(m$onsets[m$ccode == 10L & m$year == 1953L], 1L)   # dispute 50 once, in 1953
  expect_equal(m$onsets[m$ccode == 10L & m$year == 1954L], 1L)   # dispute 51 only; NOT dispute 50 again
  expect_equal(sum(m$onsets[m$ccode == 10L]), 2L)                 # two distinct disputes in total
  ## the naive per-year n_distinct(dispnum[disponset == 1]) rule would give 2 for 1954
  naive <- dd |> filter(orig1 == 1L, ccode1 == 10L, year == 1954L) |> summarise(n = n_distinct(dispnum[disponset == 1L]))
  expect_equal(naive$n, 2L)
  d <- count_dy(dd)
  expect_equal(d$onsets[d$ccode1 == 10L & d$ccode2 == 20L & d$year == 1953L], 1L)
  expect_equal(d$onsets[d$ccode1 == 10L & d$ccode2 == 20L & d$year == 1954L], 1L)   # dispute 51
  expect_false(any(d$ccode2 == 30L))                                                   # joiner pair excluded (originators only)
})

## ---- distance covariate and the Z-column flags -------------------------------------------------------------
source(file.path(root, "R", "assemble_cell.R"))
test_that("capital distance enters as log1p(km) measured at t-1, like the other dyadic controls", {
  dy <- tibble(ccode1 = 10L, ccode2 = 20L, year = c(1989L, 1990L), capdist = c(1000, 1200)) |> mutate(log_capdist = log1p(capdist))
  alts <- tibble(ccode1 = 10L, ccode2 = 20L, year = 1991L) |> left_join(dy |> mutate(year = year + 1L), by = c("ccode1", "ccode2", "year"))
  expect_equal(alts$log_capdist, log1p(1200))            # the 1990 (t-1) value, not 1989
})
test_that("Z columns: main spec has land contiguity, the target power changes and no distance; flags add distance / drop contiguity / add the missing indicator", {
  ## 2026-10-24: dyadic peace years left the main Z, the lagged DV dy_mid_tm1 entered; --dy_peace_years restores the former;
  ## later that day alliance and the major-power count left the main Z as well (minimal controls; --full_controls restores them)
  expect_equal(z_columns(), c("contig_land", "joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "dM_j_h", "dRpr_j_h"))
  expect_false("log_dy_peace_years" %in% z_columns())
  expect_equal(z_columns(DY_PEACE_YEARS = TRUE), c("contig_land", "joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "log_dy_peace_years", "dM_j_h", "dRpr_j_h"))
  for (ch in c(TRUE, FALSE)) expect_true("dy_mid_tm1" %in% z_columns(CONFLICT_HISTORY = ch))    # the lagged DV is not a --conflict_history control
  expect_false("dy_mid_tm1" %in% Z_CENTRED)                                                  # 0/1, not centred
  expect_equal(SY_INDICATORS, c("cowmaj_tm1", "dem_tm1", "ongoing_tm1"))
  expect_true("log_capdist" %in% z_columns(DISTANCE = TRUE)); expect_true("contig_land" %in% z_columns(DISTANCE = TRUE))
  expect_false("contig_land" %in% z_columns(CONTIGUITY = FALSE))
  expect_false("log_dy_prior_mids" %in% z_columns(CONFLICT_HISTORY = FALSE))
  ## the target changes are in every specification, including the one the 1-D models use (z_columns is model-agnostic)
  for (ch in c(TRUE, FALSE)) for (di in c(TRUE, FALSE)) for (co in c(TRUE, FALSE))
    expect_true(all(c("dM_j_h", "dRpr_j_h") %in% z_columns(ch, di, co)))
  expect_false("dch_j_missing" %in% z_columns())
  expect_equal(tail(z_columns(TARGET_MISSING = TRUE), 1), "dch_j_missing")
  expect_true(all(c("dM_j_h", "dRpr_j_h") %in% Z_CENTRED)); expect_false("dch_j_missing" %in% Z_CENTRED)
})

## ---- positions = raw log shares (PI decision 2026-10-23); gaps are log ratios --------------------------------------
test_that("M = log CINC share and R_pr = log PageRank without within-year standardisation; gaps are log ratios", {
  pos <- tibble(ccode = rep(1:4, 2), year = rep(c(1990L, 1991L), each = 4),
                cinc_share = c(0.5, 0.3, 0.15, 0.05, 0.4, 0.4, 0.1, 0.1),
                pagerank = c(0.4, 0.3, 0.2, 0.1, 0.25, 0.25, 0.25, 0.25)) |>
    mutate(M = log(cinc_share), R_pr = log(pagerank))      # the build's mutate
  expect_equal(pos$M, log(pos$cinc_share)); expect_equal(pos$R_pr, log(pos$pagerank))
  expect_error(assert_standardized(pos, "M"))              # deliberately NOT standardised within year
  expect_true(assert_finite(pos, c("M", "R_pr")))
  y <- pos |> filter(year == 1990L)
  expect_equal(y$M[y$ccode == 1] - y$M[y$ccode == 2], log(0.5 / 0.3))            # log ratio, no year scale
  expect_equal(y$M[y$ccode == 1] - y$M[y$ccode == 4], log(0.5 / 0.05))
  expect_equal(y$R_pr[y$ccode == 1] - y$R_pr[y$ccode == 4], log(0.4 / 0.1))
  ## an equal-share year gives identical R_pr for every state (a within-year z-score would be undefined / 0 here)
  expect_equal(sd(pos$R_pr[pos$year == 1991L]), 0)
  ## the raw change of state 3 is the log ratio of its own share across the window
  h1 <- pos |> filter(ccode == 3) |> arrange(year)
  expect_equal(h1$M[2] - h1$M[1], log(0.1 / 0.15))
})

test_that("a zero, NA or non-finite PageRank stops the build like a zero CINC share", {
  pos <- tibble(ccode = 1:4, year = 1990L, pagerank = c(0.5, 0, NA, Inf))
  bad <- pos |> filter(!is.finite(pagerank) | pagerank <= 0)
  expect_equal(sort(bad$ccode), 2:4)
  expect_error(if (nrow(bad) > 0) stop("state-years have pagerank <= 0 or non-finite"), "pagerank")
  src <- readLines(file.path(v2_root(), "R", "01_build_data.R"))
  expect_true(any(grepl("bad_pr <- pos", src, fixed = TRUE)))
  expect_false(any(grepl("assert_standardized(pos, \"M\")", src, fixed = TRUE)))     # no z-score assertion on M / R_pr any more
  expect_false(any(grepl("assert_standardized(pos, \"R_pr\")", src, fixed = TRUE)))
})

## ---- power changes: one table per lag, year-centred, one pooled SD (build_change_table in helpers_data.R) ----------
## toy levels (log units): states 10 and 40 have every year 1988-1994; state 20 has no level in 1992 (gap);
## state 30 enters in 1992. The state-year panel adds 20 in 1992 (in the system, no position that year).
toy_pos <- tribble(
  ~ccode, ~year, ~M,   ~R_pr,
  10L, 1988L, -0.3, -1.2,   10L, 1989L, -0.1, -1.1,   10L, 1990L, 0.0, -1.0,   10L, 1991L, 0.2, -0.9,
  10L, 1992L, 0.5, -0.7,    10L, 1993L, 0.4, -0.6,    10L, 1994L, 0.9, -0.2,
  40L, 1988L, 1.0, 0.4,     40L, 1989L, 1.1, 0.3,     40L, 1990L, 0.9, 0.2,    40L, 1991L, 1.2, 0.0,
  40L, 1992L, 1.0, 0.1,     40L, 1993L, 1.3, 0.2,     40L, 1994L, 1.1, 0.6,
  20L, 1990L, -1.0, 0.0,    20L, 1991L, -1.4, 0.1,    20L, 1993L, -1.1, 0.3,   20L, 1994L, -1.0, 0.5,
  30L, 1992L, 2.0, 1.0,     30L, 1993L, 2.5, 1.4,     30L, 1994L, 2.2, 1.1)
toy_states <- toy_pos |> select(ccode, year) |> bind_rows(tibble(ccode = 20L, year = 1992L))

test_that("change table: raw change is level(t-1) - level(t-1-h) by calendar join; rows without both levels are dropped", {
  tb <- build_change_table(toy_pos, toy_states, h = 1L)$tbl
  expect_equal(unique(tb$h), 1L)
  ## state 10: change at t needs levels at t-1 and t-2 -> t = 1990..1994
  expect_equal(tb$year[tb$ccode == 10L], 1990:1994)
  expect_equal(tb$dM_raw_h[tb$ccode == 10L & tb$year == 1993L], 0.5 - 0.2)             # M(1992) - M(1991)
  expect_equal(tb$dRpr_raw_h[tb$ccode == 10L & tb$year == 1994L], -0.6 - (-0.7))      # R(1993) - R(1992)
  ## state 20 (no 1992 level): t = 1992 uses 1991 and 1990 (present); t = 1993 and 1994 need 1992 -> dropped, never bridged
  expect_equal(tb$year[tb$ccode == 20L], 1992L)
  expect_equal(tb$dM_raw_h[tb$ccode == 20L], -1.4 - (-1.0))
  ## state 30 entered 1992: its first change is at t = 1994 (levels 1993 and 1992)
  expect_equal(tb$year[tb$ccode == 30L], 1994L)
  expect_equal(tb$dM_raw_h[tb$ccode == 30L], 2.5 - 2.0)
  ## a row-lag version would have bridged state 20's gap (1993 - 1991); the calendar join does not
  rowlag <- toy_pos |> filter(ccode == 20L) |> arrange(year) |> mutate(d = M - lag(M))
  expect_equal(rowlag$d[3], -1.1 - (-1.4))
  expect_false(any(tb$ccode == 20L & tb$year == 1994L))
  ## h = 3: t needs t-1 and t-4; state 20 qualifies at t = 1994 (levels 1993 and 1990), 30 never does
  tb3 <- build_change_table(toy_pos, toy_states, h = 3L)$tbl
  expect_equal(tb3$dM_raw_h[tb3$ccode == 10L & tb3$year == 1994L], 0.4 - 0.0)
  expect_equal(tb3$year[tb3$ccode == 20L], 1994L); expect_equal(tb3$dRpr_raw_h[tb3$ccode == 20L], 0.3 - 0.0)
  expect_false(30L %in% tb3$ccode)
  expect_equal(unique(tb3$h), 3L)
})

test_that("change table: centred within year (mean of all states with a change that year), then ONE pooled SD per lag and variable", {
  ch <- build_change_table(toy_pos, toy_states, h = 1L)
  tb <- ch$tbl
  ## within-year mean of the transformed change is 0 in every year
  yr <- tb |> group_by(year) |> summarise(mM = mean(dM_h), mR = mean(dRpr_h), n = n(), .groups = "drop")
  expect_true(all(abs(yr$mM) < 1e-12)); expect_true(all(abs(yr$mR) < 1e-12))
  ## the year mean is taken over the states that HAVE a change that year: 1992 has 10, 40 and 20 (not 30)
  expect_equal(yr$n[yr$year == 1992L], 3L)
  ## pooled SD after scaling is exactly 1 for each variable, and the constants are the SDs of the year-centred raw change
  expect_equal(sd(tb$dM_h), 1); expect_equal(sd(tb$dRpr_h), 1)
  cM <- tb |> group_by(year) |> mutate(c = dM_raw_h - mean(dM_raw_h)) |> ungroup()
  expect_equal(ch$scale$pooled_sd[ch$scale$variable == "dM_h"], sd(cM$c))
  cR <- tb |> group_by(year) |> mutate(c = dRpr_raw_h - mean(dRpr_raw_h)) |> ungroup()
  expect_equal(ch$scale$pooled_sd[ch$scale$variable == "dRpr_h"], sd(cR$c))
  ## ONE constant: the ratio (transformed / centred raw) is the same number on every row (no per-year scaling)
  ratio <- tb$dM_h / cM$c
  expect_true(all(abs(ratio[is.finite(ratio)] - 1 / sd(cM$c)) < 1e-12))
  ## hand check one row: state 10 at t = 1993: raw 0.3; the other 1993 changes are 40: 1.0 - 1.2 = -0.2 (20 and 30 have none)
  expect_equal(tb$dM_h[tb$ccode == 10L & tb$year == 1993L], (0.3 - mean(c(0.3, -0.2))) / sd(cM$c))
  ## scale frame shape and the raw columns are kept
  expect_equal(ch$scale$h, c(1L, 1L)); expect_equal(ch$scale$variable, c("dM_h", "dRpr_h")); expect_equal(unique(ch$scale$n_rows), nrow(tb))
  expect_true(all(c("dM_raw_h", "dRpr_raw_h", "dM_h", "dRpr_h") %in% names(tb)))
  ## the year-centring removes any variable that is constant within a year, e.g. system growth: cor is exactly 0
  g <- c(`1990` = 0.05, `1991` = 0.10, `1992` = 0.00, `1993` = 0.20, `1994` = -0.02)[as.character(tb$year)]
  expect_true(abs(cor(tb$dM_h, g)) < 1e-10); expect_true(abs(cor(tb$dRpr_h, g)) < 1e-10)
  ## fail loudly: a non-finite level, a duplicate key, or too few rows
  expect_error(build_change_table(toy_pos |> mutate(M = ifelse(ccode == 10L & year == 1990L, NA, M)), toy_states, 1L), "non-finite")
  expect_error(build_change_table(bind_rows(toy_pos, toy_pos[1, ]), toy_states, 1L), "duplicate")
  expect_error(build_change_table(toy_pos, tibble(ccode = 10L, year = 1990L), 1L), "only 1 state-years")
})

test_that("target changes join from the same table; a target with no change value gets 0 and dch_j_missing = 1", {
  h <- 1L
  ch <- build_change_table(toy_pos, toy_states, h)
  panel <- ch$tbl                                               # initiator rows = the change table rows
  ## choice set at t = 1993: initiators 10 and 40 (20 has no change at 1993, 30 neither); targets = every other state at t
  alts <- panel |> filter(year == 1993L) |> select(ccode1 = ccode, year) |>
    inner_join(toy_states |> rename(ccode2 = ccode), by = "year", relationship = "many-to-many") |> filter(ccode1 != ccode2) |>
    left_join(ch$tbl |> transmute(ccode2 = ccode, year, dM_j_h = dM_h, dRpr_j_h = dRpr_h), by = c("ccode2", "year"))
  expect_true(all(is.na(alts$dM_j_h) == is.na(alts$dRpr_j_h)))
  alts <- alts |> mutate(dch_j_missing = as.integer(is.na(dM_j_h)), dM_j_h = coalesce(dM_j_h, 0), dRpr_j_h = coalesce(dRpr_j_h, 0))
  expect_equal(sort(alts$ccode2[alts$ccode1 == 10L]), c(20L, 30L, 40L))
  ## targets 20 and 30 have no change at 1993 -> 0 and indicator 1; target 40 gets ITS OWN transformed change (same table, same key)
  expect_equal(alts$dch_j_missing[alts$ccode1 == 10L & alts$ccode2 %in% c(20L, 30L)], c(1L, 1L))
  expect_equal(alts$dM_j_h[alts$ccode1 == 10L & alts$ccode2 == 20L], 0); expect_equal(alts$dRpr_j_h[alts$ccode1 == 10L & alts$ccode2 == 30L], 0)
  expect_equal(alts$dch_j_missing[alts$ccode1 == 10L & alts$ccode2 == 40L], 0L)
  expect_equal(alts$dM_j_h[alts$ccode1 == 10L & alts$ccode2 == 40L], panel$dM_h[panel$ccode == 40L & panel$year == 1993L])
  expect_equal(alts$dM_j_h[alts$ccode1 == 40L & alts$ccode2 == 10L], panel$dM_h[panel$ccode == 10L & panel$year == 1993L])
  ## the indicator decision used by assemble_cell(): in Z only when non-constant in the assembled cell
  n_miss <- sum(alts$dch_j_missing); target_missing <- n_miss > 0 && n_miss < nrow(alts)
  expect_true(target_missing); expect_true("dch_j_missing" %in% z_columns(TARGET_MISSING = target_missing))
  alts_all_ok <- alts |> filter(dch_j_missing == 0L); n0 <- sum(alts_all_ok$dch_j_missing)
  expect_false(n0 > 0 && n0 < nrow(alts_all_ok)); expect_false("dch_j_missing" %in% z_columns(TARGET_MISSING = n0 > 0 && n0 < nrow(alts_all_ok)))
})
test_that("a zero CINC share stops the build rather than being floored", {
  pos <- tibble(ccode = 1:3, year = 1990L, cinc_share = c(0.6, 0.4, 0))
  bad <- pos |> filter(is.na(cinc_share) | cinc_share <= 0)
  expect_equal(nrow(bad), 1); expect_equal(bad$ccode, 3L)
  expect_error(if (nrow(bad) > 0) stop("state-years have cinc_share <= 0"), "cinc_share")
})

test_that("a missing CINC share stops the build check before any filtering (Codex review of 36c0f7b)", {
  ## mirrors 01_build_data.R section 3: the zero/NA check runs on the full panel, before filter()
  pos <- tibble::tibble(ccode = 1:4, year = 1990L, cinc_share = c(0.6, 0.3, 0.1, NA))
  bad <- pos |> dplyr::filter(is.na(cinc_share) | cinc_share <= 0)
  expect_equal(nrow(bad), 1L)
  expect_error(if (nrow(bad) > 0) stop("state-years have cinc_share <= 0 or NA"), "cinc_share")
  src <- readLines(file.path(v2_root(), "R", "01_build_data.R"))
  i_check <- grep("bad_cinc <- pos", src, fixed = TRUE)[1]
  i_filter <- grep("filter(!is.na(cinc_share))", src, fixed = TRUE)
  expect_true(!is.na(i_check))
  expect_true(length(i_filter) == 0 || all(i_filter > i_check))   # no NA filter ahead of the check
})

test_that("gap scale (2026-10-23): one constant per axis equals the SD of within-year pair gaps", {
  set.seed(123)
  d <- data.frame(year = rep(c(1900L, 1950L, 2000L), times = c(4L, 7L, 12L)))
  d$v <- rnorm(nrow(d), mean = c(0, 3, -2)[match(d$year, c(1900L, 1950L, 2000L))], sd = 1.5)
  ## brute force over ordered pairs within year
  pr <- do.call(rbind, lapply(split(d, d$year), function(z) { g <- outer(z$v, z$v, "-"); data.frame(g = g[row(g) != col(g)]) }))
  expect_equal(gap_scale(d$v, d$year), sqrt(mean(pr$g^2)), tolerance = 1e-12)
  expect_equal(mean(pr$g), 0, tolerance = 1e-12)
  ## year means (e.g. log PageRank falling as the system grows) do not enter the constant
  d2 <- d; d2$v <- d2$v + 10 * as.integer(factor(d2$year))
  expect_equal(gap_scale(d2$v, d2$year), gap_scale(d$v, d$year), tolerance = 1e-12)
  expect_error(gap_scale(c(1, NA), c(1L, 1L)))
})

test_that("assert_fit_current (2026-10-23): a fit made from older data files is rejected", {
  td <- tempfile(); dir.create(td); mp <- file.path(td, "manifest.csv")
  write.csv(data.frame(file = c("events.rds", "alts_lag3_all.rds"), nrow = 1, ncol = 1, sha256 = c("aaa", "bbb"), note = ""), mp, row.names = FALSE)
  ok <- list(cell = "x", manifest = data.frame(file = c("events.rds", "alts_lag3_all.rds"), sha256 = c("aaa", "bbb")))
  stale <- list(cell = "x", manifest = data.frame(file = c("events.rds", "alts_lag3_all.rds"), sha256 = c("aaa", "OLD")))
  gone <- list(cell = "x", manifest = data.frame(file = "stage1_panel_lag9.rds", sha256 = "ccc"))
  expect_true(assert_fit_current(ok, mp))
  expect_error(assert_fit_current(stale, mp), "alts_lag3_all.rds")
  expect_error(assert_fit_current(gone, mp), "stage1_panel_lag9.rds")
})

## ---- joiners, ongoing disputes, lagged involvement (PI decisions 2026-10-24) ------------------------------------------
## Toy directed dispute-years (both directions, as in cow_mid_dirdisps):
##   dispute 70: 10 (A) vs 20 (B), 1980-1982; 30 joins Side B in 1982 (orig 0) -> pair (10, 30) onset 1982 (a LATER JOINER)
##   dispute 71: 10 (A) vs 20 (B), onset 1981, while dispute 70 is still active -> an ESCALATION of the ongoing fight
##   dispute 72: 40 (A) vs 50 (B), 1979 only; dispute 73: 40 (A) vs 50 (B), onset 1980 -> the earlier dispute ENDED at t-1
##   dispute 74: 10, 11 (A) vs 20, 21 (B), 1990, all originators -> multilateral A x B expansion (4 events)
mk_dd <- function(dispnum, a, b, years, orig_a = 1L, orig_b = 1L) {
  f <- tibble(dispnum = dispnum, ccode1 = a, ccode2 = b, year = years, disponset = as.integer(years == min(years)), dispongoing = 1L,
              sidea1 = 1L, sidea2 = 0L, orig1 = orig_a, orig2 = orig_b, fatality1 = 0L, fatality2 = 0L, hostlev1 = 3L, hostlev2 = 3L)
  bind_rows(f, f |> mutate(ccode1 = b, ccode2 = a, sidea1 = 0L, sidea2 = 1L, orig1 = orig_b, orig2 = orig_a))
}
toy_dd2 <- bind_rows(mk_dd(70L, 10L, 20L, 1980:1982), mk_dd(70L, 10L, 30L, 1982L, orig_b = 0L),
                     mk_dd(71L, 10L, 20L, 1981L), mk_dd(72L, 40L, 50L, 1979L), mk_dd(73L, 40L, 50L, 1980L),
                     mk_dd(74L, 10L, 20L, 1990L), mk_dd(74L, 10L, 21L, 1990L), mk_dd(74L, 11L, 20L, 1990L), mk_dd(74L, 11L, 21L, 1990L))
## dispute dates: 70 starts 1980-05 and lasts 900 days (ends in 1982); 71 starts 1981-03 while 70 is active; 72 is 1979-06,
## 10 days; 73 starts 1980-02; 74 is 1990-04
toy_disp2 <- bind_rows(mk_disp(70L, 0L, 3L, 1980L, 5L, 900), mk_disp(71L, 0L, 3L, 1981L, 3L, 20), mk_disp(72L, 0L, 3L, 1979L, 6L, 10),
                       mk_disp(73L, 0L, 3L, 1980L, 2L, 15), mk_disp(74L, 0L, 3L, 1990L, 4L, 10))

test_that("already-fighting keys: pair in a dispute active at t with pair onset before t, both directions, joiners included", {
  af <- already_fighting_keys(toy_dd2)
  k <- paste(af$ccode1, af$ccode2, af$year, af$dispnum)
  expect_true(all(c("10 20 1981 70", "20 10 1981 70", "10 20 1982 70", "20 10 1982 70") %in% k))
  expect_false("10 20 1980 70" %in% k)                            # onset year itself is not "already fighting"
  expect_false(any(af$dispnum %in% c(71L, 72L, 73L, 74L)))        # one-year disputes contribute nothing
  expect_false("10 30 1982 70" %in% k)                            # the later joiner's pair onset IS 1982: a selection, not ongoing
  expect_false("40 50 1980 73" %in% k); expect_false("40 50 1980 72" %in% k)   # dispute 72 ended in 1979
  expect_equal(nrow(af), 4L)
})

test_that("events with joiners: later joiner pair has its own onset year; multilateral A x B expansion; one direction per pair", {
  ev <- code_events(toy_dd2, toy_disp2, include_joiners = TRUE)
  e70 <- ev |> filter(dispnum == 70L) |> arrange(target)
  expect_equal(e70$target, c(20L, 30L)); expect_equal(e70$year, c(1980L, 1982L))
  expect_equal(e70$is_joiner_pair, c(FALSE, TRUE))
  e74 <- ev |> filter(dispnum == 74L)
  expect_equal(nrow(e74), 4L); expect_setequal(paste(e74$initiator, e74$target), c("10 20", "10 21", "11 20", "11 21"))
  expect_false(any(e74$is_joiner_pair))
  expect_false(any(paste(ev$dispnum, ev$target, ev$initiator) %in% paste(ev$dispnum, ev$initiator, ev$target)))
  expect_true(all(ev$initiator %in% c(10L, 11L, 40L)))            # Side A -> Side B only
  ## without joiners the (10, 30) pair is absent and nothing else changes
  old <- code_events(toy_dd2, toy_disp2, include_joiners = FALSE)
  expect_equal(nrow(old), nrow(ev) - 1L); expect_false(30L %in% old$target)
})

test_that("escalation: a second dispute between an already-fighting pair leaves stage 2 AND stage 1; a pair whose dispute ended at t-1 does not", {
  ev <- code_events(toy_dd2, toy_disp2, include_joiners = TRUE)
  e71 <- ev |> filter(dispnum == 71L)
  expect_true(e71$escalation); expect_false(e71$in_mid); expect_false(e71$in_fmid); expect_false(e71$in_war)
  expect_true(nrow(ev |> filter(dispnum == 71L)) == 1L)            # the row is kept, flagged, not deleted
  e73 <- ev |> filter(dispnum == 73L)
  expect_false(e73$escalation); expect_true(e73$in_mid)            # dispute 72 ended in 1979: 1980 is a fresh selection
  expect_false(any(ev$escalation[ev$dispnum != 71L]))
  ## stage 1: y1 from the in_mid flag; state 10 in 1981 had only the escalation -> y1 = 0; 1980 and 1982 -> 1
  panel <- tibble(ccode = c(10L, 10L, 10L, 40L, 40L), year = c(1980L, 1981L, 1982L, 1979L, 1980L))
  ev_sy <- ev |> group_by(ccode = initiator, year) |> summarise(n_ev_mid = sum(in_mid), n_all = n(), .groups = "drop")
  p <- panel |> left_join(ev_sy, by = c("ccode", "year")) |> mutate(across(c(n_ev_mid, n_all), ~ coalesce(.x, 0L)), y1_mid = as.integer(n_ev_mid > 0))
  expect_equal(p$y1_mid, c(1L, 0L, 1L, 1L, 1L))
  expect_equal(p$n_all[p$ccode == 10L & p$year == 1981L], 1L)      # the state-year had an event row, but it was an escalation
  ## recount from the non-escalation events agrees with y1
  re <- ev |> filter(!escalation) |> distinct(ccode = initiator, year) |> mutate(y = 1L)
  expect_equal(coalesce((panel |> left_join(re, by = c("ccode", "year")))$y, 0L), p$y1_mid)
  ## stage 2: choice sets at the state-year's reference date drop already-fighting alternatives (month rule, year fallback);
  ## the escalation's own row goes with them, no other event's does
  states <- tibble(ccode = c(10L, 11L, 20L, 21L, 30L, 40L, 50L, 60L))
  alts <- tidyr::expand_grid(ccode1 = c(10L, 40L), year = 1979:1982) |>
    cross_join(states |> rename(ccode2 = ccode)) |> filter(ccode1 != ccode2)
  pdi <- pair_dispute_intervals(toy_dd2, dispute_date_bounds(toy_disp2, toy_dd2))
  ref <- reference_dates(ev)
  expect_equal(ref$year[ref$ccode == 10L], c(1980L, 1982L, 1990L)) # 1981 had only the escalation: January is the reference
  expect_equal(ref$ref_lo[ref$ccode == 10L & ref$year == 1980L], ym_index(1980L, 5L))
  afm <- already_fighting_keys_month(pdi, ref)
  af <- afm |> filter(drop) |> distinct(ccode1, ccode2, year) |> mutate(af = 1L)
  expect_equal(afm$drop, afm$drop_year)                            # on this toy the month rule and the year rule agree on every key
  expect_true(afm$fallback[afm$ccode1 == 10L & afm$ccode2 == 20L & afm$year == 1982L & afm$dispnum == 70L])   # reference interval = all of 1982 (joiner entry month unknown)
  expect_false(afm$fallback[afm$ccode1 == 10L & afm$ccode2 == 20L & afm$year == 1981L & afm$dispnum == 70L])
  alts2 <- alts |> left_join(af, by = c("ccode1", "ccode2", "year")) |> filter(is.na(af)) |> select(-af)
  expect_equal(nrow(alts) - nrow(alts2), 2L)                       # (10, 20, 1981) and (10, 20, 1982)
  expect_false(any(alts2$ccode1 == 10L & alts2$ccode2 == 20L & alts2$year %in% 1981:1982))
  expect_true(any(alts2$ccode1 == 10L & alts2$ccode2 == 20L & alts2$year == 1980L))
  expect_true(any(alts2$ccode1 == 10L & alts2$ccode2 == 30L & alts2$year == 1982L))   # the joiner event's alternative stays
  expect_true(any(alts2$ccode1 == 40L & alts2$ccode2 == 50L & alts2$year == 1980L))   # ended-at-t-1 pair stays
  kept <- ev |> filter(in_mid, year <= 1982L) |> transmute(ccode1 = initiator, ccode2 = target, year)
  expect_equal(nrow(kept |> anti_join(alts2, by = c("ccode1", "ccode2", "year"))), 0L)   # every kept event has its row
  expect_equal(nrow(ev |> filter(escalation) |> transmute(ccode1 = initiator, ccode2 = target, year) |> semi_join(alts2, by = c("ccode1", "ccode2", "year"))), 0L)
})

test_that("ongoing_tm1 (stage 1) and dy_mid_tm1 (stage 2): involvement in any MID active at t-1, any role, joiners included", {
  act_s <- active_state_years(toy_dd2); act_p <- active_pair_years(toy_dd2)
  panel <- tibble(ccode = c(10L, 10L, 10L, 30L, 40L, 40L, 60L), year = c(1980L, 1981L, 1983L, 1983L, 1980L, 1981L, 1980L)) |>
    left_join(act_s |> transmute(ccode, year = year + 1L, ongoing_tm1 = active), by = c("ccode", "year")) |> mutate(ongoing_tm1 = coalesce(ongoing_tm1, 0L))
  expect_equal(panel$ongoing_tm1, c(0L, 1L, 1L, 1L, 1L, 1L, 0L))
  ## 10 in 1980: nothing in 1979 -> 0; 10 in 1981 (dispute 70 active 1980) -> 1; 10 and the JOINER 30 in 1983 (active 1982) -> 1
  ## 40 in 1980 (dispute 72 in 1979) -> 1; 40 in 1981 (dispute 73 in 1980) -> 1; 60 never -> 0
  expect_true(all(panel$ongoing_tm1 %in% c(0L, 1L))); expect_false(anyNA(panel$ongoing_tm1))
  alts <- tibble(ccode1 = c(40L, 40L, 10L, 30L, 20L, 10L), ccode2 = c(50L, 60L, 30L, 10L, 10L, 20L), year = c(1980L, 1980L, 1983L, 1983L, 1983L, 1980L)) |>
    left_join(act_p |> transmute(ccode1, ccode2, year = year + 1L, dy_mid_tm1 = active), by = c("ccode1", "ccode2", "year")) |> mutate(dy_mid_tm1 = coalesce(dy_mid_tm1, 0L))
  expect_equal(alts$dy_mid_tm1, c(1L, 0L, 1L, 1L, 1L, 0L))
  ## (40, 50, 1980): dispute 72 ended at t-1 -> 1 and the row is NOT already fighting; (40, 60): unrelated -> 0;
  ## (10, 30) and (30, 10) in 1983: the joiner pair was active in 1982 -> 1 in both directions; (10, 20, 1980): nothing in 1979 -> 0
  expect_false(any(paste(alts$ccode1, alts$ccode2, alts$year) %in% with(already_fighting_keys(toy_dd2), paste(ccode1, ccode2, year))))
  ## a remaining alternative with dy_mid_tm1 = 1 implies every dispute of the pair active at t-1 had ended by t (otherwise it would be already fighting)
  expect_false("dy_mid_tm1" %in% Z_CENTRED); expect_true("dy_mid_tm1" %in% z_columns())
})

test_that("WWII dispute list is derived from the COW data and asserted against WW2_DISPNUMS", {
  skip_if_not_installed("peacesciencer")
  suppressPackageStartupMessages(library(peacesciencer))
  dd <- cow_mid_dirdisps |> mutate(across(c(dispnum, ccode1, ccode2, year), as.integer))
  tbl <- ww2_disputes(cow_mid_disps |> mutate(dispnum = as.integer(dispnum)), dd, cow_war_inter)
  expect_setequal(tbl$dispnum[tbl$is_ww2], WW2_DISPNUMS)
  expect_true(258L %in% tbl$dispnum[tbl$is_ww2])
  expect_true(all(c(157L, 179L, 183L, 613L) %in% tbl$dispnum[!tbl$is_ww2]))   # the Sino-Japanese, Winter, Nomonhan and Franco-Thai wars
  expect_true(all(tbl$hiact[tbl$is_ww2 & tbl$dispnum != 258L] == WW2_HIACT_JOIN_WAR))
  expect_true(all(tbl$n_dyadyears_in_war139[tbl$is_ww2] > 0))
  expect_error(ww2_disputes(cow_mid_disps |> mutate(dispnum = as.integer(dispnum)) |> filter(dispnum != 258L), dd, cow_war_inter), "WW2_DISPNUMS")
})

## ---- month-level escalation rule, pair-level war, duplicate collapse, control sets (PI decisions 2026-10-24) -------------
## Toy disputes between 10 (A) and 20 (B), all originators, hostlev 3 unless stated:
##   80: 1989-11, 120 days -> ends 1990-02/03 (a 1990 row exists)        81: 1990-06 between the same pair
##       -> the YEAR rule calls 81 an escalation (80 active in 1990 with onset 1989); the MONTH rule does not (80 ended by June)
##   82: 1991-02, 200 days -> ends 1991-08/09                            83: 1991-06 between the same pair
##       -> the year rule misses it (82 onset 1991 is not before 1991); the month rule calls 83 an escalation (82 still active)
##   84: 1993-01, 150 days -> ends between 1993-05 and 1993-06           85: 1993-06 between the same pair
##       -> undecided at month resolution (the end bound straddles June): year fallback -> not an escalation, esc_fallback TRUE
##   86: 1995-03 by 30 (A) vs 20 (B), 700 days; 10 joins side A in 1996 (entry month unknown) and stays to the end (1997)
##   87: 1996-09 between 10 and 20 -> 10's entry into 86 is somewhere in 1996: undecided, year fallback (onset 1996 not before
##       1996) -> not an escalation; 88: 1997-02 between 10 and 20 -> entry in 1996 is certainly before 1997-02 and 86 is
##       certainly active (ends 1997-01/02 .. no: 700 days from 1995-03 end in 1997-01/02, straddling) -> see the test
toy_dd3 <- bind_rows(mk_dd(80L, 10L, 20L, 1989:1990), mk_dd(81L, 10L, 20L, 1990L), mk_dd(82L, 10L, 20L, 1991L), mk_dd(83L, 10L, 20L, 1991L),
                     mk_dd(84L, 10L, 20L, 1993L), mk_dd(85L, 10L, 20L, 1993L),
                     mk_dd(86L, 30L, 20L, 1995:1997), mk_dd(86L, 10L, 20L, 1996:1997, orig_a = 0L), mk_dd(87L, 10L, 20L, 1996L), mk_dd(88L, 10L, 20L, 1997L))
toy_disp3 <- bind_rows(mk_disp(80L, 0L, 3L, 1989L, 11L, 120), mk_disp(81L, 0L, 3L, 1990L, 6L, 5), mk_disp(82L, 0L, 3L, 1991L, 2L, 200), mk_disp(83L, 0L, 3L, 1991L, 6L, 5),
                       mk_disp(84L, 0L, 3L, 1993L, 1L, 150), mk_disp(85L, 0L, 3L, 1993L, 6L, 5), mk_disp(86L, 0L, 3L, 1995L, 3L, 700),
                       mk_disp(87L, 0L, 3L, 1996L, 9L, 5), mk_disp(88L, 0L, 3L, 1997L, 2L, 5))

test_that("dispute date bounds: exact start month, end bounded by the unknown start day and the duration, clamped to the last dispute-year", {
  db <- dispute_date_bounds(toy_disp3, toy_dd3)
  d80 <- db[db$dispnum == 80L, ]
  expect_equal(d80$start_lo, ym_index(1989L, 11L)); expect_equal(d80$start_hi, d80$start_lo); expect_true(d80$start_month_known)
  expect_equal(d80$end_lo, ym_index(1990L, 2L)); expect_equal(d80$end_hi, ym_index(1990L, 3L))   # Nov 1 + 119 d = Feb 28; Nov 30 + 120 d = Mar 30
  expect_equal(d80$end_rule, "duration")
  d84 <- db[db$dispnum == 84L, ]
  expect_equal(c(d84$end_lo, d84$end_hi), c(ym_index(1993L, 5L), ym_index(1993L, 6L)))         # Jan 1 + 149 d = May 30; Jan 31 + 150 d = Jun 30
  ## a dispute whose rows run past the duration-implied end keeps the whole last year (end_rule year_only)
  dd_long <- bind_rows(toy_dd3, mk_dd(80L, 10L, 20L, 1991L) |> mutate(disponset = 0L))
  db2 <- dispute_date_bounds(toy_disp3, dd_long)
  expect_equal(db2$end_rule[db2$dispnum == 80L], "year_only"); expect_equal(c(db2$end_lo[db2$dispnum == 80L], db2$end_hi[db2$dispnum == 80L]), ym_index(1991L, c(1L, 12L)))
  ## a start year that disagrees with the first dispute-year makes the start month unknown (whole first year)
  db3 <- dispute_date_bounds(toy_disp3 |> mutate(styear = ifelse(dispnum == 80L, 1990L, styear)), toy_dd3)
  expect_false(db3$start_month_known[db3$dispnum == 80L]); expect_equal(c(db3$start_lo[db3$dispnum == 80L], db3$start_hi[db3$dispnum == 80L]), ym_index(1989L, c(1L, 12L)))
  expect_error(dispute_date_bounds(toy_disp3 |> filter(dispnum != 80L), toy_dd3), "no dispute-level row")
  expect_error(dispute_date_bounds(toy_disp3 |> mutate(stmon = ifelse(dispnum == 80L, 13L, stmon)), toy_dd3), "start month")
  ## interval comparisons are three-valued
  expect_equal(interval_before(1, 1, 2, 2), TRUE); expect_equal(interval_before(2, 2, 2, 2), FALSE); expect_true(is.na(interval_before(1, 3, 2, 2)))
  expect_equal(interval_not_before(5, 6, 5, 5), TRUE); expect_equal(interval_not_before(3, 4, 5, 5), FALSE); expect_true(is.na(interval_not_before(4, 5, 5, 5)))
  expect_equal(month_rule(1, 1, 5, 6, 3, 3), TRUE); expect_equal(month_rule(1, 1, 2, 2, 3, 3), FALSE); expect_equal(month_rule(3, 3, 5, 6, 3, 3), FALSE)
  expect_equal(month_rule(1, 1, 3, 4, 3, 3), TRUE)                 # an end in the reference month or later counts as still active
  expect_true(is.na(month_rule(1, 1, 2, 4, 3, 3)))                 # an end bound straddling the reference month is undecided
})

test_that("pair-dispute intervals: later entrants get the whole entry year, early leavers the whole exit year, re-entries are separate spells", {
  pdi <- pair_dispute_intervals(toy_dd3, dispute_date_bounds(toy_disp3, toy_dd3))
  p86_30 <- pdi |> filter(dispnum == 86L, ccode1 == 30L, ccode2 == 20L); p86_10 <- pdi |> filter(dispnum == 86L, ccode1 == 10L, ccode2 == 20L)
  expect_true(p86_30$start_exact); expect_equal(p86_30$start_lo, ym_index(1995L, 3L))
  expect_false(p86_10$start_exact); expect_equal(c(p86_10$start_lo, p86_10$start_hi), ym_index(1996L, c(1L, 12L)))   # joined in 1996, month unknown
  expect_equal(c(p86_10$pair_onset, p86_10$pair_last), c(1996L, 1997L))
  expect_equal(c(p86_10$end_lo, p86_10$end_hi), c(p86_30$end_lo, p86_30$end_hi))                                    # both stay to the dispute's end
  ## an early leaver: 30 leaves 86 after 1996 -> its exit interval is all of 1996
  dd_leave <- toy_dd3 |> filter(!(dispnum == 86L & (ccode1 == 30L | ccode2 == 30L) & year == 1997L))
  pl <- pair_dispute_intervals(dd_leave, dispute_date_bounds(toy_disp3, dd_leave)) |> filter(dispnum == 86L, ccode1 == 30L)
  expect_equal(c(pl$end_lo, pl$end_hi), ym_index(1996L, c(1L, 12L)))
  ## a re-entry (pair absent in a middle year) is a second spell with its own unknown entry month
  dd_re <- bind_rows(toy_dd3, mk_dd(86L, 30L, 20L, 1995:1999) |> filter(year == 1999L) |> mutate(disponset = 0L)) |> filter(!(dispnum == 86L & ccode1 %in% c(30L, 20L) & ccode2 %in% c(30L, 20L) & year == 1998L))
  pr <- pair_dispute_intervals(dd_re, dispute_date_bounds(toy_disp3 |> mutate(mindur = ifelse(dispnum == 86L, 1500, mindur), maxdur = mindur), dd_re)) |> filter(dispnum == 86L, ccode1 == 30L)
  expect_equal(pr$spell, c(1L, 2L)); expect_equal(pr$spell_first, c(1995L, 1999L)); expect_equal(unique(pr$pair_onset), 1995L)
  expect_equal(pr$start_lo[2], ym_index(1999L, 1L)); expect_equal(pr$start_hi[2], pr$end_hi[2]); expect_true(pr$start_clipped[2])   # entry bound clipped to the feasible exit (D5)
  expect_equal(c(pr$end_lo[1], pr$end_hi[1]), ym_index(1997L, c(1L, 12L)))                       # spell 1 exits in 1997, month unknown
})

test_that("month rule: ended-before-the-new-start is not an escalation, still-active is, an undecided month falls back to the year rule", {
  ev <- code_events(toy_dd3, toy_disp3, include_joiners = TRUE)
  e <- function(d) ev[ev$dispnum == d, ]
  ## (a) 80 (1989-11, ends by 1990-03) then 81 (1990-06): the year rule says escalation, the month rule says a fresh selection
  expect_true(e(81L)$escalation_year); expect_false(e(81L)$escalation); expect_false(e(81L)$esc_fallback); expect_true(e(81L)$in_mid)
  ## (b) 82 (1991-02, ends 1991-08/09) then 83 (1991-06): still active -> escalation; the year rule missed it (within-year leak)
  expect_false(e(83L)$escalation_year); expect_true(e(83L)$escalation); expect_false(e(83L)$esc_fallback); expect_false(e(83L)$in_mid); expect_false(e(83L)$in_fmid)
  ## (c) 84 (1993-01, ends between 1993-05 and 1993-06) then 85 (1993-06): undecided -> year fallback (84 onset 1993 is not before 1993)
  expect_false(e(85L)$escalation); expect_true(e(85L)$esc_fallback); expect_false(e(85L)$escalation_year)
  ## (d) 10 joined 86 some time in 1996; 87 (1996-09) is undecided -> year fallback -> not an escalation, flagged as a fallback decision
  expect_false(e(87L)$escalation); expect_true(e(87L)$esc_fallback)
  expect_equal(c(e(87L)$ev_start_lo, e(87L)$ev_start_hi), rep(ym_index(1996L, 9L), 2))
  ## the joiner pair's own event in 86 has an unknown month (whole of 1996) and is not an escalation of anything
  expect_equal(c(e(86L)$ev_start_lo[e(86L)$initiator == 10L], e(86L)$ev_start_hi[e(86L)$initiator == 10L]), ym_index(1996L, c(1L, 12L)))
  expect_false(any(e(86L)$escalation))
  ## (e) 88 (1997-02): 10's entry in 1996 is certainly before 1997-02; 86 ends 1997-01/02 (700 days from 1995-03), straddling -> fallback;
  ##     the year rule says escalation (onset 1996 < 1997, 86 active in 1997)
  expect_true(e(88L)$escalation); expect_true(e(88L)$esc_fallback); expect_true(e(88L)$escalation_year)
  ## the first dispute of each pair is never an escalation; the comparison table is exposed for logging
  expect_false(any(ev$escalation[ev$dispnum %in% c(80L, 82L, 84L)]))
  cmp <- attr(ev, "comparisons"); expect_true(all(c("month_dec", "year_dec", "fallback", "decision") %in% names(cmp)))
  expect_equal(sum(cmp$fallback), 4L)                            # 85 vs 84, 87 vs 86, 88 vs 86, and 86 (10->20, entry month unknown) vs 87; 86 vs 88 is a certain FALSE
  ## the year rule inside the month helper reproduces flag_escalations() exactly
  expect_equal(ev$escalation_year, flag_escalations(ev |> select(initiator, target, year, dispnum), already_fighting_keys(toy_dd3))$escalation)
  ## stage 2: choice sets at the reference date. 10 initiates in 1990 at 1990-06 (81; 80 is the previous year's event), so 20 is
  ## not already fighting at the reference date (80 ended by March) and stays; in 1991 the reference is 82 (February), 20 stays
  ## (82 itself starts at the reference, not before it); in 1993 the reference is 84 (January), 20 stays
  pdi <- pair_dispute_intervals(toy_dd3, dispute_date_bounds(toy_disp3, toy_dd3)); ref <- reference_dates(ev)
  afm <- already_fighting_keys_month(pdi, ref)
  dropped <- afm |> filter(drop) |> distinct(ccode1, ccode2, year)
  expect_false(any(dropped$ccode1 == 10L & dropped$ccode2 == 20L & dropped$year %in% c(1990L, 1991L, 1993L)))
  expect_true(any(afm$drop_year & afm$ccode1 == 10L & afm$ccode2 == 20L & afm$year == 1990L))   # the year rule would have dropped 20 in 1990
  ## every kept event's chosen alternative is in its choice set
  kept <- ev |> filter(!escalation) |> transmute(ccode1 = initiator, ccode2 = target, year)
  expect_equal(nrow(kept |> semi_join(dropped, by = c("ccode1", "ccode2", "year"))), 0L)
  ## a state-year without a non-escalation event uses January: 20 (never an initiator) in 1996 has 30 dropped for sure (86 began
  ## 1995-03) and 10 kept for sure (10 entered 86 some time in 1996, which cannot be before January); in 1997 (reference January)
  ## 10's entry in 1996 is certainly before and 86 is certainly active in January 1997 (ends 1997-01/02) -> dropped, no fallback
  k <- afm |> filter(ccode1 == 20L, year == 1996L)
  expect_true(k$drop[k$ccode2 == 30L]); expect_false(k$fallback[k$ccode2 == 30L])
  expect_false(k$fallback[k$ccode2 == 10L & k$dispnum == 86L]); expect_false(k$drop[k$ccode2 == 10L & k$dispnum == 86L])
  k7 <- afm |> filter(ccode1 == 20L, year == 1997L, ccode2 == 10L, dispnum == 86L)
  expect_true(k7$drop); expect_false(k7$fallback)
})

test_that("war is pair-level: dispute war-level AND at least one participant at hostility level 5; the fatal war rule follows it", {
  ## dispute 90, hostlev 5, four pairs: (10, 20) both at 5; (10, 21) target at 3; (11, 20) initiator at 3; (11, 21) neither at 5;
  ## fatality all missing
  dd <- bind_rows(mk_dd(90L, 10L, 20L, 1970L), mk_dd(90L, 10L, 21L, 1970L), mk_dd(90L, 11L, 20L, 1970L), mk_dd(90L, 11L, 21L, 1970L)) |>
    mutate(hostlev1 = ifelse(ccode1 %in% c(10L, 20L), 5L, 3L), hostlev2 = ifelse(ccode2 %in% c(10L, 20L), 5L, 3L), fatality1 = -9L, fatality2 = -9L)
  ev <- code_events(dd, mk_disp(90L, -9L, 5L, 1970L, 4L, 100), include_joiners = TRUE) |> arrange(initiator, target)
  expect_equal(ev$war_disp, c(TRUE, TRUE, TRUE, TRUE))
  expect_equal(ev$war, c(TRUE, TRUE, TRUE, FALSE))
  expect_equal(ev$in_war, c(TRUE, TRUE, TRUE, FALSE))
  expect_equal(ev$fatal_main, c(TRUE, TRUE, TRUE, FALSE))         # the war rule makes the war pairs fatal, not the pair with neither side at 5
  expect_equal(ev$in_fmidub, c(TRUE, TRUE, TRUE, TRUE))           # the uncoded pair is fatal only in the upper bound
  ## participant hostlev is the max over the pair's years
  dd2 <- bind_rows(mk_dd(90L, 10L, 20L, 1970:1971)) |> mutate(hostlev1 = ifelse(year == 1971L, 5L, 4L), hostlev2 = ifelse(year == 1971L, 5L, 4L))
  expect_true(code_events(dd2, mk_disp(90L, 0L, 5L, 1970L, 4L, 400))$war)
})

## the extra columns the build adds after code_events() and before the collapse (so that collapse_duplicate_events() sees
## every column it classifies): COW-only / upper-bound codings, the source string, the no-WWII sample flags (review P2-3)
add_collapse_cols <- function(ev, ww2_disputes = integer(0)) {
  ev |> mutate(is_ww2 = dispnum %in% ww2_disputes, fatal_cow = fatal %in% TRUE, fatal_ub = fatal_main | (is.na(fatal) & is.na(fatal_mic) & is.na(fatal_gml)),
               fatal_source = ifelse(fatal %in% TRUE, "cow", "missing_nonfatal"), target_has_pos_tm1 = TRUE, fighting_at_ref = FALSE,
               in_mid_noww2 = in_mid & !is_ww2, in_fmid_noww2 = in_fmid & !is_ww2, in_fmidcow_noww2 = in_fmidcow & !is_ww2,
               in_fmidub_noww2 = in_fmidub & !is_ww2, in_war_noww2 = in_war & !is_ww2) |>
    select(-hostlev_disp, -fatality_disp, -hostlev1, -hostlev2, -fatality1, -fatality2, hostlev1, hostlev2, fatality1, fatality2, hostlev_disp, fatality_disp)
}
test_that("same-year duplicates collapse into one event per (initiator, target, year) with ORed flags; escalations untouched", {
  ## 10 -> 20 twice in 1990 (disputes 80, January, 10 days, a joiner pair, and 81, May, 60 days, originators); 81 is not an
  ## escalation under the month rule (80 ended in January for sure: the joiner pair's entry is unknown within 1990 but its end is
  ## January); one is fatal + war. 82 (June) starts while 81 is active -> escalation, kept as its own row
  dd <- bind_rows(mk_dd(80L, 10L, 20L, 1990L, orig_b = 0L), mk_dd(81L, 10L, 20L, 1990L) |> mutate(hostlev1 = 5L, hostlev2 = 5L, fatality1 = 5L),
                  mk_dd(82L, 10L, 20L, 1990L), mk_dd(83L, 11L, 20L, 1990L))
  disp <- bind_rows(mk_disp(80L, 0L, 3L, 1990L, 1L, 10), mk_disp(81L, 5L, 5L, 1990L, 5L, 60), mk_disp(82L, 0L, 3L, 1990L, 6L, 100), mk_disp(83L, 0L, 3L, 1990L, 7L, 5))
  ev <- code_events(dd, disp, include_joiners = TRUE) |> add_collapse_cols(ww2_disputes = 82L)
  expect_equal(nrow(ev), 4L)
  expect_true(ev$escalation[ev$dispnum == 82L])                   # 82 (June) starts while 81 (May, 60 days) is active -> escalation
  expect_false(ev$escalation[ev$dispnum == 81L])                  # 80 ended in January
  out <- collapse_duplicate_events(ev)
  expect_equal(nrow(out), 3L)                                     # (10, 20, 1990) once, the escalation row, (11, 20, 1990)
  k <- out |> filter(initiator == 10L, target == 20L, !escalation)
  expect_equal(nrow(k), 1L); expect_equal(k$dispnum, 80L); expect_equal(k$dispnums, "80;81"); expect_equal(k$n_disputes, 2L)
  expect_true(k$in_mid); expect_true(k$in_fmid); expect_true(k$in_war); expect_true(k$war); expect_true(k$fatal_main); expect_true(k$is_ww2 == FALSE)
  expect_false(k$is_joiner_pair)                                  # 80 is a joiner pair, 81 an originator pair -> not a joiner pair
  expect_equal(k$fatal_source, "cow;missing_nonfatal"); expect_equal(k$hostlev1, 5); expect_equal(k$ev_start_lo, ym_index(1990L, 1L))
  expect_true(k$in_mid_noww2); expect_true(k$in_war_noww2)        # neither dispute is WWII
  esc <- out |> filter(escalation); expect_equal(esc$dispnum, 82L); expect_equal(esc$dispnums, "82"); expect_equal(esc$n_disputes, 1L); expect_true(esc$is_ww2)
  expect_equal(out$dispnums[out$initiator == 11L], "83")
  expect_error(collapse_duplicate_events(ev |> mutate(extra = 1)), "unclassified")
  expect_error(collapse_duplicate_events(ev |> select(-in_war)), "absent")
  expect_error(collapse_duplicate_events(ev |> select(-in_war_noww2)), "absent")
  expect_error(collapse_duplicate_events(ev |> mutate(target_has_pos_tm1 = dispnum == 80L)), "differs")
  ## one chosen alternative per (i, j, t) afterwards
  expect_equal(nrow(out |> filter(!escalation) |> count(initiator, target, year) |> filter(n > 1)), 0L)
})

test_that("review P2-3: a collapsed key mixing a WWII and a non-WWII dispute keeps the no-WWII flags of the non-WWII dispute only", {
  ## 10 -> 20 in 1990 through dispute 80 (January, 10 days, non-fatal, hostlev 3) and dispute 81 (May, war-level, fatal), 81 a WWII
  ## dispute in this toy; 11 -> 20 through the WWII dispute 81 alone; 12 -> 20 through the non-WWII dispute 83 alone
  dd <- bind_rows(mk_dd(80L, 10L, 20L, 1990L), mk_dd(81L, 10L, 20L, 1990L) |> mutate(hostlev1 = 5L, hostlev2 = 5L, fatality1 = 5L),
                  mk_dd(81L, 11L, 20L, 1990L) |> mutate(hostlev1 = 5L, hostlev2 = 5L, fatality1 = 5L), mk_dd(83L, 12L, 20L, 1990L))
  disp <- bind_rows(mk_disp(80L, 0L, 3L, 1990L, 1L, 10), mk_disp(81L, 5L, 5L, 1990L, 5L, 60), mk_disp(83L, 0L, 3L, 1990L, 7L, 5))
  ev <- code_events(dd, disp, include_joiners = TRUE) |> add_collapse_cols(ww2_disputes = 81L)
  expect_false(any(ev$escalation))
  ## before the collapse the no-WWII flag is the sample flag of the non-WWII rows
  expect_equal(ev$in_mid_noww2, !ev$is_ww2); expect_equal(ev$in_fmid_noww2[ev$dispnum == 81L], c(FALSE, FALSE))
  out <- collapse_duplicate_events(ev)
  k <- out |> filter(initiator == 10L)
  expect_equal(k$dispnums, "80;81"); expect_true(k$is_ww2)                     # the ORed WWII flag marks the mixed key
  expect_true(k$in_mid); expect_true(k$in_fmid); expect_true(k$in_war)        # the sample flags take the WWII dispute's fatal / war coding
  expect_true(k$in_mid_noww2)                                                 # the non-WWII dispute 80 qualifies for mid ...
  expect_false(k$in_fmid_noww2); expect_false(k$in_war_noww2)                  # ... but is neither fatal nor war: no WWII coding leaks in
  expect_false(k$in_mid & !k$is_ww2)                                          # the old rule in_mid & !is_ww2 would have dropped the event
  k11 <- out |> filter(initiator == 11L); expect_true(k11$in_mid); expect_false(k11$in_mid_noww2)      # pure WWII: dropped
  k12 <- out |> filter(initiator == 12L); expect_true(k12$in_mid_noww2); expect_equal(k12$in_mid_noww2, k12$in_mid & !k12$is_ww2)   # pure non-WWII: unchanged
  ## the flag is always a subset of the sample flag
  for (s in c("mid", "fmid", "fmidcow", "fmidub", "war")) expect_true(all(out[[paste0("in_", s, "_noww2")]] <= out[[paste0("in_", s)]]))
})

test_that("control sets: minimal main W / X / Z, and --full_controls restores the previous specification exactly", {
  ## main W since 2026-10-25: the power levels (D2) after dM_i, dR_i and the centred log system size (D3) last; X unchanged
  expect_equal(w_controls(), c("M_lvl", "R_lvl", "ongoing_tm1", "log_prior_mids", "dem_tm1", "log_n_system_tm1"))
  expect_equal(c("dM_i", "dR_i", w_controls()), c("dM_i", "dR_i", "M_lvl", "R_lvl", "ongoing_tm1", "log_prior_mids", "dem_tm1", "log_n_system_tm1"))
  expect_equal(x_controls(), character(0))
  expect_equal(w_controls(CONFLICT_HISTORY = FALSE), c("M_lvl", "R_lvl", "ongoing_tm1", "dem_tm1", "log_n_system_tm1"))
  expect_equal(W_LEVELS, c("M_lvl", "R_lvl")); expect_false(any(W_LEVELS %in% x_controls()))
  expect_false("cowmaj_tm1" %in% w_controls()); expect_false("log_peace_years" %in% w_controls())
  expect_equal(z_columns(), c("contig_land", "joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "dM_j_h", "dRpr_j_h"))
  expect_false(any(c("alliance_def_off", "major_count") %in% z_columns()))
  ## controls = full (PI decision 2026-10-27, item 7): the MAIN sets plus the controls the minimal sets dropped, WITHOUT major-power status
  ## (its own switch): W = main W + log_peace_years; X = main X (no controls); Z = main Z + alliance_def_off + log_dy_peace_years after
  ## dy_mid_tm1. Until then `full` reproduced the 8bf0ce1 base (cowmaj_tm1 / major_count, no D2 / D3 columns, no ongoing_tm1, no dy_mid_tm1).
  full_w <- c("M_lvl", "R_lvl", "ongoing_tm1", "log_prior_mids", "dem_tm1", "log_n_system_tm1", "log_peace_years")
  full_z <- function(CONFLICT_HISTORY = TRUE, DISTANCE = FALSE, CONTIGUITY = TRUE, TARGET_MISSING = FALSE)
    c(if (CONTIGUITY) "contig_land", "joint_dem", if (CONFLICT_HISTORY) "log_dy_prior_mids", "dy_mid_tm1", "alliance_def_off", "log_dy_peace_years",
      if (DISTANCE) "log_capdist", "dM_j_h", "dRpr_j_h", if (TARGET_MISSING) "dch_j_missing")
  expect_equal(w_controls(FULL_CONTROLS = TRUE), full_w); expect_equal(x_controls(FULL_CONTROLS = TRUE), character(0))
  expect_equal(w_controls(FULL_CONTROLS = TRUE), c(w_controls(), "log_peace_years"))                                   # main + one column, appended
  expect_equal(w_controls(CONFLICT_HISTORY = FALSE, FULL_CONTROLS = TRUE), w_controls(CONFLICT_HISTORY = FALSE))      # the peace years are conflict history
  expect_false(any(c("cowmaj_tm1", "major_count") %in% c(w_controls(FULL_CONTROLS = TRUE), z_columns(FULL_CONTROLS = TRUE))))   # no major-power status in `full`
  expect_true(all(c("ongoing_tm1", "M_lvl", "R_lvl", "log_n_system_tm1") %in% w_controls(FULL_CONTROLS = TRUE))); expect_true("dy_mid_tm1" %in% z_columns(FULL_CONTROLS = TRUE))
  for (ch in c(TRUE, FALSE)) for (di in c(TRUE, FALSE)) for (co in c(TRUE, FALSE)) for (tm in c(TRUE, FALSE)) {
    expect_equal(z_columns(ch, di, co, TARGET_MISSING = tm, FULL_CONTROLS = TRUE), full_z(ch, di, co, tm))
    expect_equal(setdiff(z_columns(ch, di, co, TARGET_MISSING = tm, FULL_CONTROLS = TRUE), z_columns(ch, di, co, TARGET_MISSING = tm)), Z_FULL_EXTRA)   # exactly the two extras
  }
  expect_equal(z_columns(FULL_CONTROLS = TRUE), c("contig_land", "joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "alliance_def_off", "log_dy_peace_years", "dM_j_h", "dRpr_j_h"))
  ## under the full set the dyadic peace years are already in Z; DY_PEACE_YEARS adds nothing then and never duplicates the column
  expect_equal(z_columns(FULL_CONTROLS = TRUE, DY_PEACE_YEARS = TRUE), z_columns(FULL_CONTROLS = TRUE))
  expect_equal(z_columns(DY_PEACE_YEARS = TRUE), c("contig_land", "joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "log_dy_peace_years", "dM_j_h", "dRpr_j_h"))
  expect_true("log_dy_peace_years" %in% Z_CENTRED); expect_false("alliance_def_off" %in% Z_CENTRED)                 # continuous centred, 0/1 not
  ## the option flags combine with either control set the same way
  for (fc in c(TRUE, FALSE)) {
    expect_false("contig_land" %in% z_columns(CONTIGUITY = FALSE, FULL_CONTROLS = fc))
    expect_true("log_capdist" %in% z_columns(DISTANCE = TRUE, FULL_CONTROLS = fc))
    expect_false("log_dy_prior_mids" %in% z_columns(CONFLICT_HISTORY = FALSE, FULL_CONTROLS = fc))
    expect_equal(tail(z_columns(TARGET_MISSING = TRUE, FULL_CONTROLS = fc), 1), "dch_j_missing")
    expect_true(all(x_controls(FULL_CONTROLS = fc) %in% w_controls(FULL_CONTROLS = fc)))
  }
  ## every W / X control is one the build writes and the centring list covers the 0/1 indicators that enter
  expect_true(all(w_controls(FULL_CONTROLS = TRUE) %in% c(SY_INDICATORS, W_LEVELS, "log_prior_mids", "log_n_system_tm1", "log_peace_years")))
  expect_true(all(w_controls() %in% c(SY_INDICATORS, W_LEVELS, "log_prior_mids", "log_n_system_tm1")))
  src <- readLines(file.path(v2_root(), "R", "02_fit_model.R")); spec_src <- readLines(file.path(v2_root(), "R", "spec.R"))
  ## the full-controls option is the `controls` switch of spec.R (suffix _fullctrl) and reaches assemble_cell() through assemble_cell_spec()
  expect_true(any(grepl('controls  = c("minimal", "full")', spec_src, fixed = TRUE))); expect_true(any(grepl('full = "_fullctrl"', spec_src, fixed = TRUE)))
  expect_true(any(grepl("assemble_cell_spec(SAMPLE, LAG, spec", src, fixed = TRUE)))
  expect_true(any(grepl("manifest_at_load", src, fixed = TRUE))); expect_true(any(grepl("SLURM_ARRAY_TASK_ID", src, fixed = TRUE)))
  expect_true(any(grepl("changed while sampling", src, fixed = TRUE))); expect_true(any(grepl("_STALE", src, fixed = TRUE)))   # C2: saved as STALE, not discarded
  expect_true(any(grepl("assert_compiled_current", src, fixed = TRUE))); expect_true(any(grepl('"_pilot"', src, fixed = TRUE)))   # C3, C1
  expect_true(any(grepl("file.rename(tmp_path, fit_path)", src, fixed = TRUE)))                                                  # C2 atomic save
})

## ---- Codex review fixes (2026-10-24): joiner entry month (P1-2), two choice sets (P1-1), no-WWII flags (P2-3), base control
## ---- set (P2-4), manifest snapshot (P2-5), empty / chosen-only sets --------------------------------------------------------------
test_that("P1-2 same-year joiner: a pair with a joiner entering in the dispute's first year has entry interval [dispute start, December]", {
  ## dispute 90: 10 (A, originator) vs 20 (B, originator) from 1990-03 for 400 days (into 1991); 30 joins Side B in 1990 (orig 0)
  ## dispute 91: 10 -> 30 in 1990-07; dispute 92: 10 -> 30 in 1991-02
  dd <- bind_rows(mk_dd(90L, 10L, 20L, 1990:1991), mk_dd(90L, 10L, 30L, 1990:1991, orig_b = 0L), mk_dd(91L, 10L, 30L, 1990L), mk_dd(92L, 10L, 30L, 1991L))
  disp <- bind_rows(mk_disp(90L, 0L, 3L, 1990L, 3L, 400), mk_disp(91L, 0L, 3L, 1990L, 7L, 5), mk_disp(92L, 0L, 3L, 1991L, 2L, 5))
  db <- dispute_date_bounds(disp, dd)
  pdi <- pair_dispute_intervals(dd, db)
  p_orig <- pdi |> filter(dispnum == 90L, ccode1 == 10L, ccode2 == 20L); p_join <- pdi |> filter(dispnum == 90L, ccode1 == 10L, ccode2 == 30L)
  expect_true(p_orig$both_orig); expect_true(p_orig$start_exact); expect_equal(c(p_orig$start_lo, p_orig$start_hi), rep(ym_index(1990L, 3L), 2))
  expect_false(p_join$both_orig); expect_false(p_join$start_exact)
  expect_equal(c(p_join$start_lo, p_join$start_hi), c(ym_index(1990L, 3L), ym_index(1990L, 12L)))   # [dispute start, December of y0]
  expect_equal(p_join$spell_first, 1990L); expect_equal(c(p_join$end_lo, p_join$end_hi), c(p_orig$end_lo, p_orig$end_hi))   # exit bounds unchanged
  ## both directions of the joiner pair carry the widened entry
  p_rev <- pdi |> filter(dispnum == 90L, ccode1 == 30L, ccode2 == 10L); expect_equal(c(p_rev$start_lo, p_rev$start_hi), c(p_join$start_lo, p_join$start_hi))
  ## the pre-review rule gave the joiner pair the exact start month
  pdi_old <- pair_dispute_intervals(dd, db, joiner_entry_unknown = FALSE)
  p_join_old <- pdi_old |> filter(dispnum == 90L, ccode1 == 10L, ccode2 == 30L)
  expect_true(p_join_old$start_exact); expect_equal(c(p_join_old$start_lo, p_join_old$start_hi), rep(ym_index(1990L, 3L), 2))
  expect_identical(pdi_old$start_lo, pdi$start_lo); expect_identical(pdi_old$end_hi, pdi$end_hi)
  ## escalation flags follow: 91 (July 1990) against the joiner pair's dispute 90 is undecided (entry somewhere in March..December)
  ## -> year fallback (onset 1990 is not before 1990) -> NOT an escalation, flagged as a fallback decision; the old rule called it an
  ## escalation (March < July, still active). 92 (February 1991): the 1990 entry is certainly before and 90 is certainly active -> escalation
  ev <- code_events(dd, disp, include_joiners = TRUE)
  e91 <- ev |> filter(dispnum == 91L); e92 <- ev |> filter(dispnum == 92L)
  expect_false(e91$escalation); expect_true(e91$esc_fallback); expect_false(e91$escalation_year)
  expect_true(e92$escalation); expect_false(e92$esc_fallback)
  ev_old <- flag_escalations_month(ev |> select(initiator, target, year, dispnum), pdi_old)
  expect_true(ev_old$escalation[ev_old$dispnum == 91L]); expect_true(ev_old$escalation[ev_old$dispnum == 92L])
  ## the joiner pair's own event in 90 has the widened date interval
  e90j <- ev |> filter(dispnum == 90L, target == 30L); expect_equal(c(e90j$ev_start_lo, e90j$ev_start_hi), c(ym_index(1990L, 3L), ym_index(1990L, 12L)))
  ## originator status must not differ across a pair's dispute-years
  expect_error(pair_dispute_intervals(dd |> mutate(orig2 = ifelse(dispnum == 90L & ccode2 == 30L & year == 1991L, 1L, orig2)), db), "originator status")
  expect_error(pair_dispute_intervals(dd |> select(-orig1), db))
})

## toy for the two choice sets: dispute X = 10 (A) vs 20 (B), 1989-11 for 120 days (ends 1990-02/03, so a 1990 row exists);
## dispute Z = 40 (A) vs 10 (B), 1990-04 for 200 days (ends 1990-10/11); dispute Y = 10 (A) vs 30 (B) in 1990, start month per scenario
two_set_toy <- function(y_month) {
  dd <- bind_rows(mk_dd(1L, 10L, 20L, 1989:1990), mk_dd(2L, 40L, 10L, 1990L), if (!is.null(y_month)) mk_dd(3L, 10L, 30L, 1990L))
  disp <- bind_rows(mk_disp(1L, 0L, 3L, 1989L, 11L, 120), mk_disp(2L, 0L, 3L, 1990L, 4L, 200), if (!is.null(y_month)) mk_disp(3L, 0L, 3L, 1990L, y_month, 5))
  ev <- code_events(dd, disp, include_joiners = TRUE)
  pdi <- pair_dispute_intervals(dd, dispute_date_bounds(disp, dd))
  ref <- reference_dates(ev)
  af_iv <- already_fighting_keys_month(pdi, NULL) |> filter(drop) |> distinct(ccode1, ccode2, year)
  af_s2 <- already_fighting_keys_month(pdi, ref) |> filter(drop) |> distinct(ccode1, ccode2, year)
  grid <- tidyr::expand_grid(ccode1 = c(10L, 20L, 40L), year = 1990L, ccode2 = c(10L, 20L, 30L, 40L)) |> filter(ccode1 != ccode2)
  list(ev = ev, ref = ref, af_iv = af_iv, af_s2 = af_s2, alts = choice_set_flags(grid, af_iv, af_s2))
}
test_that("P1-1 regression: changing or adding a state-year's events never changes its stage-1 IV row set; only the stage-2 set moves", {
  A <- two_set_toy(6L)        # 10's earliest 1990 event in June: X has ended (20 back in the stage-2 set), Z is active (40 out of it)
  B <- two_set_toy(NULL)      # 10 has no 1990 event: reference January, both sets identical
  C <- two_set_toy(2L)        # 10's event in February: X's end straddles February -> undecided -> year rule drops 20 from the stage-2 set too
  iv_rows <- function(S) S$alts |> filter(ccode1 == 10L, in_iv == 1L) |> pull(ccode2) |> sort()
  s2_rows <- function(S) S$alts |> filter(ccode1 == 10L, in_s2 == 1L) |> pull(ccode2) |> sort()
  ## the IV set of (10, 1990) is {30, 40} in every scenario (20 is fighting 10 at January in dispute X; Z has not started)
  expect_equal(iv_rows(A), c(30L, 40L)); expect_equal(iv_rows(B), c(30L, 40L)); expect_equal(iv_rows(C), c(30L, 40L))
  expect_identical(A$af_iv, B$af_iv); expect_identical(A$af_iv, C$af_iv)                    # the January keys never see the events
  ## the stage-2 set does depend on the reference date
  expect_equal(s2_rows(A), c(20L, 30L)); expect_equal(s2_rows(B), c(30L, 40L)); expect_equal(s2_rows(C), c(30L, 40L))
  expect_equal(A$ref$ref_lo[A$ref$ccode == 10L & A$ref$year == 1990L], ym_index(1990L, 6L)); expect_equal(nrow(B$ref |> filter(ccode == 10L, year == 1990L)), 0L)
  ## row flags in scenario A: 20 is stage-2 only, 40 is IV only, 30 in both; in B and C the two sets coincide and 20 has no row
  a10 <- A$alts |> filter(ccode1 == 10L) |> arrange(ccode2)
  expect_equal(a10$ccode2, c(20L, 30L, 40L)); expect_equal(a10$in_iv, c(0L, 1L, 1L)); expect_equal(a10$in_s2, c(1L, 1L, 0L))
  for (S in list(B, C)) { a <- S$alts |> filter(ccode1 == 10L); expect_equal(sort(a$ccode2), c(30L, 40L)); expect_true(all(a$in_iv == 1L & a$in_s2 == 1L)) }
  expect_equal(attr(A$alts, "n_neither"), 1L); expect_equal(attr(B$alts, "n_neither"), 2L)   # (20, 10) for 20 in every scenario; (10, 20) in B and C
  ## non-initiating state-years always have identical sets: 20 never initiates
  for (S in list(A, B, C)) { a <- S$alts |> filter(ccode1 == 20L); expect_true(all(a$in_iv == a$in_s2)); expect_false(10L %in% a$ccode2) }
  ## the chosen targets of the kept events are stage-2 rows
  for (S in list(A, C)) { k <- S$ev |> filter(in_mid, year == 1990L) |> transmute(ccode1 = initiator, ccode2 = target, year)
    expect_equal(nrow(k |> anti_join(S$alts |> filter(in_s2 == 1L), by = c("ccode1", "ccode2", "year"))), 0L) }
  ## choice_set_flags() refuses duplicate keys and pre-existing flag columns; the n_neither attribute counts the removed rows
  expect_error(choice_set_flags(bind_rows(A$alts |> select(-in_iv, -in_s2), A$alts[1, ] |> select(-in_iv, -in_s2)), A$af_iv, A$af_s2), "duplicate")
  expect_error(choice_set_flags(A$alts, A$af_iv, A$af_s2))
})

## the toy cell on disk (write_toy_cell(), quiet()) lives in helper-toy-cell.R since 2026-10-25 so that test_context.R shares it

test_that("assembly: fighting_at_ref events leave stage 2 only; chosen rows are stage-2 rows; the two sets reach build_stan_data() and the R likelihood", {
  td <- write_toy_cell(tempfile("toycell"))
  A <- quiet(assemble_cell("mid", 3L, "all", data_dir = td))
  expect_equal(A$sy$ccode, c(10L, 20L, 30L, 40L)); expect_equal(A$sy$y1, c(1L, 0L, 1L, 1L))
  expect_equal(nrow(A$alts), 12L); expect_equal(sum(A$alts$in_iv), 11L); expect_equal(sum(A$alts$in_s2), 11L)
  ## usable events: 10->20, 10->30, 30->20, 30->10, 40->20; 30->40 is fighting_at_ref (its row exists with in_s2 = 0) and is not an event
  expect_equal(nrow(A$ev), 5L); expect_equal(A$n_fighting_at_ref, 1L); expect_equal(A$n_orphan, 0L)
  chosen <- A$alts[A$ev$alt_row, ]
  expect_true(all(chosen$in_s2 == 1L)); expect_false(any(chosen$ccode1 == 30L & chosen$ccode2 == 40L))
  expect_true(any(A$alts$ccode1 == 30L & A$alts$ccode2 == 40L & A$alts$in_s2 == 0L & A$alts$in_iv == 1L))   # the IV-only row stays in the block
  expect_true(all(A$alts$sy_row[A$ev$alt_row] == A$ev$sy_row))
  dat <- build_stan_data(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, x_ref = A$x_ref, allow_orphan_initiations = TRUE)
  expect_equal(dat$N1, 4L); expect_equal(dat$N2, 12L); expect_equal(dat$E, 5L)
  expect_equal(dat$in_iv, A$alts$in_iv); expect_equal(dat$in_s2, A$alts$in_s2); expect_equal(dat$P1, 8L); expect_equal(dat$P, 2L); expect_equal(dat$Pz, 7L)   # P1 = 8 since D2 / D3
  ## D8: the event weights reach the Stan data from the events file (30 -> 10, dispute 6, carries 0.5); EVENT_WEIGHTS = "none" gives 1
  expect_equal(dat$ev_w, A$ev$w); expect_true(any(dat$ev_w == 0.5)); expect_equal(sum(dat$ev_w == 1), 4L)
  A1 <- quiet(assemble_cell("mid", 3L, "all", EVENT_WEIGHTS = "none", data_dir = td)); expect_true(all(A1$ev$w == 1))
  ## D2: the levels are NOT grand-centred in assemble_cell() (0 = the system-average level that year); the other controls are
  expect_equal(A$sy$M_lvl, c(0.8, -0.3, -0.1, -0.4)); expect_lt(abs(mean(A$sy$ongoing_tm1)), 1e-12); expect_lt(abs(mean(A$sy$log_n_system_tm1)), 1e-12)
  expect_equal(A$n_flipped, 1L)
  expect_equal(dat$cs_end - dat$cs_start + 1L, rep(3L, 4))                    # one union block of 3 rows per state-year
  expect_equal(dat$n_full, c(2L, 3L, 3L, 3L))                                 # n_full counts the IV rows (10 has 2)
  expect_true(all(dat$in_s2[dat$ev_row] == 1L))
  ## the brute-force likelihood uses the IV rows for lse and the stage-2 rows for the event softmax
  set.seed(3); S <- dat$S
  pars <- list(mu0 = c(0.3, -0.2), beta = matrix(rnorm(4, 0, 0.3), 2, 2), a_i = matrix(rnorm(2 * S), S, 2), c_i = rnorm(S), b_j = rnorm(S),
               gamma = rnorm(dat$Pz), sigma_x = c(1, 1.3), rho = 0.2)
  bf <- stage2_loglik_R(dat, pars)
  r10 <- dat$cs_start[1]:dat$cs_end[1]; r30 <- dat$cs_start[3]:dat$cs_end[3]
  expect_equal(bf$lse[1], log(sum(exp(bf$util[r10][dat$in_iv[r10] == 1L]))))                       # IV of 10 excludes its stage-2-only row (40)
  expect_equal(bf$lse[3], log(sum(exp(bf$util[r30]))))                                             # IV of 30 is its whole block
  e10 <- which(dat$ev_sy == 1L)
  for (e in e10) expect_equal(bf$loglik[e], bf$util[dat$ev_row[e]] - log(sum(exp(bf$util[r10]))))  # stage 2 of 10 is its whole block
  e30 <- which(dat$ev_sy == 3L)
  for (e in e30) expect_equal(bf$loglik[e], bf$util[dat$ev_row[e]] - log(sum(exp(bf$util[r30][dat$in_s2[r30] == 1L]))))   # stage 2 of 30 excludes its IV-only row
  for (n in c(2L, 4L)) { r <- dat$cs_start[n]:dat$cs_end[n]; expect_equal(bf$lse[n], log(sum(exp(bf$util[r])))) }
  ## with every flag 1 the data equal the flag-free build exactly
  A1 <- A; A1$alts$in_iv <- 1L; A1$alts$in_s2 <- 1L
  d1 <- build_stan_data(A1$sy, A1$alts, A1$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, x_ref = A$x_ref, allow_orphan_initiations = TRUE)
  d0 <- build_stan_data(A1$sy, A1$alts |> select(-in_iv, -in_s2), A1$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, x_ref = A$x_ref, allow_orphan_initiations = TRUE)
  expect_identical(d1, d0); expect_true(all(d0$in_iv == 1L)); expect_equal(d0$n_full, rep(3L, 4))
  ## stage2_loglik_R refuses data without the flags
  expect_error(stage2_loglik_R(dat[setdiff(names(dat), "in_s2")], pars), "in_s2 are missing")
})

test_that("assembly: --drop_ww2 uses the no-WWII flags (mixed key kept, pure WWII dropped, y1 recounted); fighting_at_ref unchanged", {
  td <- write_toy_cell(tempfile("toycell"))
  A <- quiet(assemble_cell("mid", 3L, "all", DROP_WW2 = TRUE, data_dir = td))
  ## kept: 10->20, 10->30 (mixed: non-WWII dispute 4 or 5 qualifies), 30->20; dropped: 30->10 and 40->20 (pure WWII); 30->40 fighting_at_ref
  expect_equal(nrow(A$ev), 3L)
  chosen <- A$alts[A$ev$alt_row, ] |> arrange(ccode1, ccode2)
  expect_equal(paste(chosen$ccode1, chosen$ccode2), c("10 20", "10 30", "30 20"))
  expect_equal(A$sy$y1, c(1L, 0L, 1L, 0L))                                    # 40's only event was WWII -> y1 recounted to 0
  expect_equal(A$n_fighting_at_ref, 1L); expect_equal(A$n_ww2, 3L)
  ## the old rule in_mid & !is_ww2 would have dropped the mixed event as well
  ev <- readRDS(file.path(td, "events.rds")); expect_equal(sum(ev$in_mid & !ev$is_ww2), 3L); expect_equal(sum(ev$in_mid_noww2), 4L)   # 30 -> 40 (fighting_at_ref) counts in both
  ## fmid under --drop_ww2: the mixed key's fatal coding came from the WWII dispute only -> not in the no-WWII fatal sample
  F <- quiet(assemble_cell("fmid", 3L, "all", DROP_WW2 = TRUE, data_dir = td))
  expect_equal(nrow(F$ev), 2L); expect_equal(F$sy$y1, c(1L, 0L, 1L, 0L))
  ## without the option everything is kept
  B <- quiet(assemble_cell("mid", 3L, "all", data_dir = td)); expect_equal(nrow(B$ev), 5L); expect_equal(B$sy$y1, c(1L, 0L, 1L, 1L))
  ## a build without the no-WWII flags is refused
  td2 <- write_toy_cell(tempfile("toycell")); ev2 <- readRDS(file.path(td2, "events.rds")) |> select(-in_mid_noww2); saveRDS(ev2, file.path(td2, "events.rds"))
  man <- read.csv(file.path(td2, "manifest.csv"), stringsAsFactors = FALSE); man$sha256[man$file == "events.rds"] <- file_digest(file.path(td2, "events.rds")); man$ncol[man$file == "events.rds"] <- ncol(ev2)
  write.csv(man, file.path(td2, "manifest.csv"), row.names = FALSE)
  expect_error(quiet(assemble_cell("mid", 3L, "all", data_dir = td2)))
})

test_that("assembly stops on an empty stage-1 IV set and on an empty or chosen-only stage-2 set; so does build_stan_data()", {
  ## (20, 1990) keeps its rows but none in the IV set
  td <- write_toy_cell(tempfile("toycell"), alts_edit = function(a) { a$in_iv[a$ccode1 == 20L] <- 0L; a })
  expect_error(quiet(assemble_cell("mid", 3L, "all", data_dir = td)), "EMPTY stage-1 IV set")
  ## (10, 1990): only its chosen rows are stage-2 rows ... make 30 and 40 IV only: 10 -> 30 is an event, so the set {20, 30} would hold
  ## two chosen rows; drop the 10 -> 30 chosen row too by making the set {20} alone
  td <- write_toy_cell(tempfile("toycell"), alts_edit = function(a) { a$in_s2[a$ccode1 == 10L & a$ccode2 %in% c(30L, 40L)] <- 0L; a$in_iv[a$ccode1 == 10L] <- 1L; a })
  ## the chosen row of 10 -> 30 now has in_s2 = 0: that event is not usable, so the "all events have their row" check fires first
  expect_error(quiet(assemble_cell("mid", 3L, "all", data_dir = td)))
  ## a legal variant: 10 -> 30 is made fighting_at_ref so that only 10 -> 20 remains, and its set is the chosen row alone
  td <- write_toy_cell(tempfile("toycell"), alts_edit = function(a) { a$in_s2[a$ccode1 == 10L & a$ccode2 %in% c(30L, 40L)] <- 0L; a$in_iv[a$ccode1 == 10L] <- 1L; a })
  ev <- readRDS(file.path(td, "events.rds")); ev$fighting_at_ref[ev$initiator == 10L & ev$target == 30L] <- TRUE; saveRDS(ev, file.path(td, "events.rds"))
  man <- read.csv(file.path(td, "manifest.csv"), stringsAsFactors = FALSE); man$sha256[man$file == "events.rds"] <- file_digest(file.path(td, "events.rds")); write.csv(man, file.path(td, "manifest.csv"), row.names = FALSE)
  expect_error(quiet(assemble_cell("mid", 3L, "all", data_dir = td)), "empty or holds only the chosen target")
  ## build_stan_data() has the same two guards and the chosen-row guard
  td <- write_toy_cell(tempfile("toycell")); A <- quiet(assemble_cell("mid", 3L, "all", data_dir = td))
  a_empty <- A$alts; a_empty$in_iv[a_empty$sy_row == 2L] <- 0L
  expect_error(build_stan_data(A$sy, a_empty, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, allow_orphan_initiations = TRUE), "EMPTY stage-1 IV set")
  a_only <- A$alts; a_only$in_s2[a_only$sy_row == 4L & a_only$ccode2 != 20L] <- 0L      # 40's set = its chosen row (20) alone
  expect_error(build_stan_data(A$sy, a_only, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, allow_orphan_initiations = TRUE), "CHOSEN-ONLY")
  a_out <- A$alts; a_out$in_s2[A$ev$alt_row[1]] <- 0L
  expect_error(build_stan_data(A$sy, a_out, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, allow_orphan_initiations = TRUE), "outside the stage-2 set")
  a_bad <- A$alts; a_bad$in_iv[1] <- 0L; a_bad$in_s2[1] <- 0L
  expect_error(build_stan_data(A$sy, a_bad, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, allow_orphan_initiations = TRUE))
})

test_that("P2-5: a rebuild between the manifest snapshot and the load is caught; without a snapshot the live manifest is read", {
  td <- write_toy_cell(tempfile("toycell"))
  snapshot <- read.csv(file.path(td, "manifest.csv"), stringsAsFactors = FALSE)
  A <- quiet(assemble_cell("mid", 3L, "all", data_dir = td, manifest = snapshot)); expect_equal(nrow(A$ev), 5L)
  ## "rebuild": events.rds changes and the manifest on disk is updated accordingly
  ev <- readRDS(file.path(td, "events.rds")); ev$rebuilt <- 1L; saveRDS(ev, file.path(td, "events.rds"))
  man <- read.csv(file.path(td, "manifest.csv"), stringsAsFactors = FALSE); man$sha256[man$file == "events.rds"] <- file_digest(file.path(td, "events.rds")); man$ncol[man$file == "events.rds"] <- ncol(ev)
  write.csv(man, file.path(td, "manifest.csv"), row.names = FALSE)
  expect_error(quiet(assemble_cell("mid", 3L, "all", data_dir = td, manifest = snapshot)), "events.rds: sha256 digest does not match the manifest snapshot")
  B <- quiet(assemble_cell("mid", 3L, "all", data_dir = td)); expect_equal(nrow(B$ev), 5L)              # live manifest: consistent again
  ## a rebuild that changes a row count is caught by the nrow check first
  pn <- readRDS(file.path(td, "stage1_panel_lag3.rds")); saveRDS(bind_rows(pn, pn[1, ] |> mutate(ccode = 50L)), file.path(td, "stage1_panel_lag3.rds"))
  man$sha256[man$file == "stage1_panel_lag3.rds"] <- file_digest(file.path(td, "stage1_panel_lag3.rds")); man$nrow[man$file == "stage1_panel_lag3.rds"] <- 5L
  write.csv(man, file.path(td, "manifest.csv"), row.names = FALSE)
  expect_error(quiet(assemble_cell("mid", 3L, "all", data_dir = td, manifest = snapshot)), "stage1_panel_lag3.rds: nrow 5 != manifest 4")
  ## the helper accepts a path or a data.frame
  expect_error(assert_matches_manifest(file.path(td, "events.rds"), snapshot), "snapshot")
  expect_equal(nrow(assert_matches_manifest(file.path(td, "events.rds"), file.path(td, "manifest.csv"))), 6L)
  ## 02_fit_model.R passes its load-time snapshot into the assembly (through assemble_cell_spec())
  src <- readLines(file.path(v2_root(), "R", "02_fit_model.R"))
  expect_true(any(grepl("manifest = manifest_at_load", src, fixed = TRUE)))
  expect_true(any(grepl("manifest = manifest)", readLines(file.path(v2_root(), "R", "assemble_cell.R")), fixed = TRUE)))
  expect_true(any(grepl("assert_compiled_current(sm, stan_src", src, fixed = TRUE)))        # C3: the marker checks became the source / model_code identity check
})

test_that("Codex re-check P2: Z centring uses the IV rows only, so removing a stage-2-only row changes no stage-1 IV value", {
  td  <- write_toy_cell(tempfile("toycell"))
  td2 <- write_toy_cell(tempfile("toycell"), alts_edit = function(a) a |> filter(!(ccode1 == 10L & ccode2 == 40L)))   # 10 -> 40 is stage-2 only, unchosen
  A <- quiet(assemble_cell("mid", 3L, "all", data_dir = td)); B <- quiet(assemble_cell("mid", 3L, "all", data_dir = td2))
  expect_equal(nrow(A$alts) - 1L, nrow(B$alts)); expect_equal(A$Z_cols, B$Z_cols); expect_equal(A$sy, B$sy)
  ## centred Z values on the IV rows are identical
  key <- c("ccode1", "ccode2", "year")
  ia <- A$alts |> filter(in_iv == 1L) |> arrange(ccode1, ccode2); ib <- B$alts |> filter(in_iv == 1L) |> arrange(ccode1, ccode2)
  expect_equal(ia[, c(key, A$Z_cols)], ib[, c(key, B$Z_cols)])
  zc <- intersect(A$Z_cols, c("log_dy_prior_mids", "log_dy_peace_years", "log_capdist", "dM_j_h", "dRpr_j_h"))
  expect_true(all(abs(colMeans(as.matrix(ia[, zc]))) < 1e-8))
  ## stage-1 IV (lse) at fixed parameters is identical for every state-year
  da <- build_stan_data(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, x_ref = A$x_ref, allow_orphan_initiations = TRUE)
  db <- build_stan_data(B$sy, B$alts, B$ev, B$W_cols, B$X_cols, B$Z_cols, joint = 1L, x_ref = A$x_ref, allow_orphan_initiations = TRUE)
  set.seed(11); S <- da$S
  pars <- list(mu0 = c(0.3, -0.2), beta = matrix(rnorm(4, 0, 0.3), 2, 2), a_i = matrix(rnorm(2 * S), S, 2), c_i = rnorm(S), b_j = rnorm(S),
               gamma = rnorm(da$Pz), sigma_x = c(1, 1.3), rho = 0.2)
  expect_equal(stage2_loglik_R(da, pars)$lse, stage2_loglik_R(db, pars)$lse, tolerance = 1e-12)
  ## the old rule (centring over the whole union block) would have moved them
  old_centre <- function(a) { for (cc in zc) a[[cc]] <- a[[cc]] - mean(a[[cc]]); a }   # re-centring over all rows = the old rule
  dao <- build_stan_data(A$sy, old_centre(A$alts), A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, x_ref = A$x_ref, allow_orphan_initiations = TRUE)
  dbo <- build_stan_data(B$sy, old_centre(B$alts), B$ev, B$W_cols, B$X_cols, B$Z_cols, joint = 1L, x_ref = A$x_ref, allow_orphan_initiations = TRUE)
  expect_false(isTRUE(all.equal(stage2_loglik_R(dao, pars)$lse, stage2_loglik_R(dbo, pars)$lse, tolerance = 1e-12)))
})

test_that("Codex re-check P2: the manifest check hashes and loads the same bytes (file replaced after the read is not returned)", {
  td <- write_toy_cell(tempfile("toycell"))
  snapshot <- read.csv(file.path(td, "manifest.csv"), stringsAsFactors = FALSE)
  p <- file.path(td, "events.rds"); ev_a <- readRDS(p)
  ev_b <- ev_a; ev_b$rebuilt <- 1L
  ## file A is read and hashed; the file is replaced by B right after the read: the helper returns A and passes
  got <- assert_matches_manifest(p, snapshot, .after_read = function() saveRDS(ev_b, p))
  expect_identical(got, ev_a); expect_false("rebuilt" %in% names(got))
  ## now B is on disk: the same snapshot rejects it
  expect_error(assert_matches_manifest(p, snapshot), "snapshot")
  ## and with B read, replacing it by A after the read still returns B and fails the A digest
  expect_error(assert_matches_manifest(p, snapshot, .after_read = function() saveRDS(ev_a, p)), "snapshot")
  ## the digest of the in-memory bytes equals the file digest used when the manifest is written
  b <- readBin(p, "raw", n = file.info(p)$size)
  expect_equal(digest::digest(b, algo = "sha256", serialize = FALSE), file_digest(p))
})

test_that("assembly with controls = full returns exactly the main lists plus log_peace_years (W) and alliance_def_off, log_dy_peace_years (Z), in order", {
  td <- write_toy_cell(tempfile("toycell"))
  A <- quiet(assemble_cell("mid", 3L, "all", FULL_CONTROLS = TRUE, data_dir = td))
  expect_equal(A$W_cols, c("dM_i", "dR_i", "M_lvl", "R_lvl", "ongoing_tm1", "log_prior_mids", "dem_tm1", "log_n_system_tm1", "log_peace_years"))
  expect_equal(A$X_cols, c("dM_i", "dR_i"))
  expect_equal(A$Z_cols, c("contig_land", "joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "alliance_def_off", "log_dy_peace_years", "dM_j_h", "dRpr_j_h", "dch_j_missing"))
  expect_true(all(abs(colMeans(as.matrix(A$sy[, c("ongoing_tm1", "dem_tm1", "log_prior_mids", "log_n_system_tm1", "log_peace_years")]))) < 1e-8))   # grand-centred like the main controls
  expect_equal(length(A$x_ref), 2L)
  iv <- A$alts$in_iv == 1L; expect_lt(abs(mean(A$alts$log_dy_peace_years[iv])), 1e-8); expect_true(all(A$alts$alliance_def_off %in% 0:1))   # continuous extra centred over the IV rows, 0/1 extra not
  M <- quiet(assemble_cell("mid", 3L, "all", data_dir = td))
  ## the full assembly differs from the main one in the three extra columns only
  expect_identical(A$ev, M$ev); expect_identical(A$sy$y1, M$sy$y1); expect_identical(A$sy[, M$W_cols], M$sy[, M$W_cols]); expect_identical(A$alts[, M$Z_cols], M$alts[, M$Z_cols])
  expect_equal(M$W_cols, c("dM_i", "dR_i", "M_lvl", "R_lvl", "ongoing_tm1", "log_prior_mids", "dem_tm1", "log_n_system_tm1")); expect_equal(M$X_cols, c("dM_i", "dR_i"))
  expect_equal(M$Z_cols, c("contig_land", "joint_dem", "log_dy_prior_mids", "dy_mid_tm1", "dM_j_h", "dRpr_j_h", "dch_j_missing"))
})

test_that("sample_alternatives(): IV rows come from the IV set, event blocks from the stage-2 set, the output flags are resolved", {
  td <- write_toy_cell(tempfile("toycell")); A <- quiet(assemble_cell("mid", 3L, "all", data_dir = td))
  sa <- sample_alternatives(A$alts, A$ev, m = 1L, seed = 11L)
  iv <- sa$alts[seq_len(sa$n_iv), ]
  expect_equal(sa$n_full, c(2L, 3L, 3L, 3L))                                   # the IV-set sizes, not the block sizes
  expect_false(any(iv$sy_row == 1L & iv$ccode2 == 40L))                        # 10's stage-2-only row never enters the IV sample
  blocks <- sa$alts[-seq_len(sa$n_iv), ]
  expect_false(any(blocks$sy_row == 3L & blocks$ccode2 == 40L))               # 30's IV-only row never enters an event block
  expect_true(all(sa$alts$in_iv == 1L), all(sa$alts$in_s2 == 1L))
  expect_true(all(sa$ev$ev_end - sa$ev$ev_start + 1L == 2L))
  dat <- build_stan_data(A$sy, sa$alts, sa$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, sampled = 1L, n_full = sa$n_full, x_ref = A$x_ref, allow_orphan_initiations = TRUE)
  expect_equal(dat$sampled, 1L); expect_equal(dat$n_full, c(2L, 3L, 3L, 3L))
  ## sampled data with an unresolved flag is refused; a chosen row outside the stage-2 set is refused by the sampler
  a_bad <- sa$alts; a_bad$in_s2[1] <- 0L
  expect_error(build_stan_data(A$sy, a_bad, sa$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, sampled = 1L, n_full = sa$n_full, allow_orphan_initiations = TRUE), "in_iv = in_s2 = 1")
  a_out <- A$alts; a_out$in_s2[A$ev$alt_row[1]] <- 0L
  expect_error(sample_alternatives(a_out, A$ev, m = 1L), "not in the stage-2 set")
  ## repeated draws over many seeds never pick an excluded row
  for (sd in 1:30) { s2 <- sample_alternatives(A$alts, A$ev, m = 2L, seed = sd)
    expect_false(any(s2$alts$sy_row[seq_len(s2$n_iv)] == 1L & s2$alts$ccode2[seq_len(s2$n_iv)] == 40L))
    expect_false(any(s2$alts$sy_row[-seq_len(s2$n_iv)] == 3L & s2$alts$ccode2[-seq_len(s2$n_iv)] == 40L)) }
})

test_that("real build (when present): the alternatives carry the flags, no empty IV set, no chosen-only stage-2 set, snapshot check works", {
  man_path <- file.path(v2_root(), "data", "manifest.csv")
  skip_if_not(file.exists(man_path) && file.exists(file.path(v2_root(), "data", "alts_lag3_all.rds")), "real data files not built")
  man <- read.csv(man_path, stringsAsFactors = FALSE)
  alts <- readRDS(file.path(v2_root(), "data", "alts_lag3_all.rds"))
  expect_true(all(c("in_iv", "in_s2") %in% names(alts))); expect_true(all(alts$in_iv == 1L | alts$in_s2 == 1L))
  n_iv <- alts |> group_by(ccode1, year) |> summarise(n = sum(in_iv), .groups = "drop"); expect_true(all(n_iv$n >= 1L))
  ev <- readRDS(file.path(v2_root(), "data", "events.rds"))
  expect_true(all(c("in_mid_noww2", "in_war_noww2") %in% names(ev))); expect_true(all(ev$in_mid_noww2 <= ev$in_mid))
  ## a tampered snapshot (the panel digest) is caught before the large files are read
  bad <- man; bad$sha256[bad$file == "stage1_panel_lag3.rds"] <- "0000"
  expect_error(quiet(assemble_cell("mid", 3L, "all", manifest = bad)), "stage1_panel_lag3.rds: sha256 digest does not match the manifest snapshot")
})

## ---- initiator direction for war-level disputes (PI decision 2026-10-25, D1) ------------------------------------------------------
## mk_war(): COW inter-state war dyad-years, both directions, for one war; init = the war side (1 or 2) with the initiator
mk_war <- function(warnum, a, b, years, init = 2L, onset_year = min(years)) {
  f <- tibble(warnum = warnum, ccode1 = a, ccode2 = b, year = years, cowinterongoing = 1L, cowinteronset = as.integer(years == onset_year),
              sidea1 = 1L, sidea2 = 2L, initiator1 = as.integer(init == 1L), initiator2 = as.integer(init == 2L))
  bind_rows(f, f |> mutate(ccode1 = b, ccode2 = a, sidea1 = 2L, sidea2 = 1L, initiator1 = as.integer(init == 2L), initiator2 = as.integer(init == 1L)))
}
## disputes: 100 (war level) 10, 11 (A) vs 20 (B, originator), 21 (B, joiner) in 1950 -- COW war 500 has 10, 11 vs 20 (20 initiated), 21 absent
##           101 (war level) 30 (A) vs 40 (B) in 1960 -- two COW wars with OPPOSITE initiators match (501 in 1960, 502 in 1961): ambiguous
##           102 (hostlev 3) 50 (A) vs 60 (B) in 1970 -- war 503 matches but the dispute is not war level: untouched
##           103 (war level) 70 (A) vs 80, 81 (B) in 1980 -- war 504 covers (70, 80) only, 80 initiated: 70 -> 80 flips, 70 -> 81 has a single
##               pair-decided event behind its dispute (below the minimum of 2): kept as Side A
##           104 (war level) 90 (A) vs 91 (B) in 1990 -- war 505, 90 initiated: confirmed, not flipped
##           105 (war level) 92 (A) vs 93 (B) in 1995 -- war 506 begins in 1996 (t + 1), 93 initiated: flips through the window
war5 <- function(x) x |> mutate(hostlev1 = 5L, hostlev2 = 5L)
toy_dd_dir <- bind_rows(war5(mk_dd(100L, 10L, 20L, 1950L)), war5(mk_dd(100L, 10L, 21L, 1950L, orig_b = 0L)), war5(mk_dd(100L, 11L, 20L, 1950L)), war5(mk_dd(100L, 11L, 21L, 1950L, orig_b = 0L)),
                        war5(mk_dd(101L, 30L, 40L, 1960L)), mk_dd(102L, 50L, 60L, 1970L), war5(mk_dd(103L, 70L, 80L, 1980L)), war5(mk_dd(103L, 70L, 81L, 1980L)),
                        war5(mk_dd(104L, 90L, 91L, 1990L)), war5(mk_dd(105L, 92L, 93L, 1995L)))
toy_disp_dir <- bind_rows(mk_disp(100L, 5L, 5L, 1950L, 3L, 200), mk_disp(101L, 5L, 5L, 1960L, 4L, 100), mk_disp(102L, 0L, 3L, 1970L, 5L, 10),
                          mk_disp(103L, 5L, 5L, 1980L, 6L, 100), mk_disp(104L, 5L, 5L, 1990L, 7L, 100), mk_disp(105L, 5L, 5L, 1995L, 8L, 300))
toy_war <- bind_rows(mk_war(500L, 10L, 20L, 1950L), mk_war(500L, 11L, 20L, 1950L),
                     mk_war(501L, 30L, 40L, 1960L, init = 1L), mk_war(502L, 30L, 40L, 1961L, init = 2L),
                     mk_war(503L, 50L, 60L, 1970L), mk_war(504L, 70L, 80L, 1980L), mk_war(505L, 90L, 91L, 1990L, init = 1L), mk_war(506L, 92L, 93L, 1996L))

test_that("D1: war dyad-years carry the initiating side; a war without exactly one initiating side is rejected", {
  wd <- war_dyad_years(toy_war)
  expect_equal(nrow(wd), nrow(toy_war)); expect_true(all(wd$side1_init != wd$side2_init))
  expect_true(wd$side2_init[wd$warnum == 500L & wd$ccode1 == 10L]); expect_false(wd$side2_init[wd$warnum == 505L & wd$ccode1 == 90L])
  expect_error(war_dyad_years(toy_war |> mutate(initiator1 = 0L, initiator2 = 0L)), "exactly one initiating side")
  expect_error(war_dyad_years(toy_war |> mutate(initiator1 = 1L, initiator2 = 1L)), "exactly one initiating side")
})

test_that("D1: pair-level flip, dispute-level flip, ambiguous and single-link pairs kept, non-war disputes untouched, the t+1 window", {
  ev <- code_events(toy_dd_dir, toy_disp_dir, include_joiners = TRUE, war_inter = toy_war)
  e <- function(d) ev |> filter(dispnum == d) |> arrange(initiator, target)
  ## dispute 100: the two matched pairs flip (pair match), the two joiner pairs flip through the dispute-level side mapping
  d100 <- e(100L)
  expect_equal(nrow(d100), 4L); expect_true(all(d100$dir_flipped)); expect_true(all(d100$initiator %in% c(20L, 21L))); expect_true(all(d100$target %in% c(10L, 11L)))
  expect_equal(sort(d100$dir_source), c("cow_war_dispute", "cow_war_dispute", "cow_war_pair", "cow_war_pair"))
  expect_equal(d100$dir_source[d100$initiator == 20L], c("cow_war_pair", "cow_war_pair")); expect_equal(d100$dir_source[d100$initiator == 21L], c("cow_war_dispute", "cow_war_dispute"))
  expect_true(all(d100$war_link == 500L)); expect_true(all(d100$sidea1 == 0L & d100$sidea2 == 1L))
  ## originator status moved with the state: 20 (originator) -> 10 is an originator pair, 21 (joiner) -> 10 a joiner pair
  expect_false(d100$is_joiner_pair[d100$initiator == 20L & d100$target == 10L]); expect_true(d100$is_joiner_pair[d100$initiator == 21L & d100$target == 10L])
  expect_equal(d100$orig1[d100$initiator == 21L], c(0L, 0L)); expect_equal(d100$orig2[d100$initiator == 21L], c(1L, 1L))
  ## dispute 101: two wars with opposite initiators -> ambiguous, Side A kept
  expect_false(e(101L)$dir_flipped); expect_equal(e(101L)$dir_source, "sideA"); expect_equal(e(101L)$dir_rule, "pair_ambiguous"); expect_equal(e(101L)$initiator, 30L)
  ## dispute 102 is not war level: untouched although a war dyad matches
  expect_false(e(102L)$dir_flipped); expect_equal(e(102L)$dir_rule, "non_war_dispute"); expect_true(is.na(e(102L)$war_link))
  ## dispute 103: 70 -> 80 flips (pair match); 70 -> 81 is kept because only one pair-decided event backs the dispute mapping
  d103 <- e(103L)
  expect_true(any(d103$initiator == 80L & d103$target == 70L & d103$dir_flipped & d103$dir_source == "cow_war_pair"))
  expect_true(any(d103$initiator == 70L & d103$target == 81L & !d103$dir_flipped & d103$dir_rule == "dispute_no_majority_link"))
  ## dispute 104: the Side A state initiated the war -> confirmed, not flipped, still "cow_war_pair"
  expect_false(e(104L)$dir_flipped); expect_equal(e(104L)$dir_source, "cow_war_pair"); expect_equal(e(104L)$initiator, 90L)
  ## dispute 105: the war dyad-year is t + 1 -> matched and flipped; with window 0 it is unmatched
  expect_true(e(105L)$dir_flipped); expect_equal(e(105L)$initiator, 93L)
  pair0 <- build_pair_table(toy_dd_dir, toy_disp_dir) |> filter(sidea1 == 1L, sidea2 == 0L) |> rename(initiator = ccode1, target = ccode2)
  d0 <- cow_war_direction(pair0, war_dyad_years(toy_war), window = 0L)
  expect_equal(d0$dir_rule[d0$dispnum == 105L], "dispute_no_pair_match"); expect_false(d0$dir_flipped[d0$dispnum == 105L])
  ## one direction per dispute pair still holds, every flipped row has sidea1 = 0, and the attributes list the matches and links
  expect_false(any(paste(ev$dispnum, ev$target, ev$initiator) %in% paste(ev$dispnum, ev$initiator, ev$target)))
  expect_equal(ev$sidea1 == 0L, ev$dir_flipped)
  dl <- attr(cow_war_direction(pair0, war_dyad_years(toy_war)), "dispute_links")
  expect_equal(dl$n_pair_decided[dl$dispnum == 100L], 2L); expect_true(dl$map_ok[dl$dispnum == 100L]); expect_false(dl$map_ok[dl$dispnum == 103L])
  ## the sideA build option leaves every direction alone and marks the source
  evA <- code_events(toy_dd_dir, toy_disp_dir, include_joiners = TRUE, war_inter = toy_war, direction = "sideA")
  expect_false(any(evA$dir_flipped)); expect_true(all(evA$dir_source == "sideA")); expect_true(all(evA$initiator %in% c(10L, 11L, 30L, 50L, 70L, 90L, 92L)))
  ## without war data nothing is flipped (the pre-2026-10-25 coding)
  ev0 <- code_events(toy_dd_dir, toy_disp_dir, include_joiners = TRUE)
  expect_false(any(ev0$dir_flipped)); expect_true(all(ev0$dir_rule %in% c("non_war_dispute", "dispute_no_pair_match")))
  ## downstream: y1 follows the applied direction (20 initiates in 1950, 10 does not)
  y1 <- ev |> filter(in_mid) |> distinct(initiator, year)
  expect_true(any(y1$initiator == 20L & y1$year == 1950L)); expect_false(any(y1$initiator == 10L & y1$year == 1950L))
  ## apply_direction() swaps every per-participant column and refuses a table that is not Side A -> Side B
  expect_error(apply_direction(ev |> mutate(sidea1 = 1L)), "sidea")
})

test_that("D1 on the real data (when present): Germany -> Poland 1939 exists, Poland -> Germany 1939 does not", {
  f <- file.path(v2_root(), "data", "events.rds")
  if (!file.exists(f)) skip("real build not present")
  ev <- readRDS(f)
  expect_true(any(ev$initiator == 255L & ev$target == 290L & ev$year == 1939L))
  expect_false(any(ev$initiator == 290L & ev$target == 255L & ev$year == 1939L))
  g <- ev |> filter(initiator == 255L, target == 290L, year == 1939L)
  expect_true(g$dir_flipped); expect_equal(g$dir_source, "cow_war_pair"); expect_equal(g$war_link, 139L); expect_true(g$in_war)
  expect_true(all(c("dir_flipped", "dir_source", "dir_rule", "war_link", "w_mid", "w_war") %in% names(ev)))
  expect_equal(ev$dir_flipped, ev$sidea2 == 1L)
  ## USSR -> Japan 1939 (Nomonhan) and China -> USSR 1929 are reversed as the review found
  expect_true(any(ev$initiator == 740L & ev$target == 365L & ev$year == 1939L & ev$dir_flipped))
  expect_true(any(ev$initiator == 365L & ev$target == 710L & ev$year == 1929L & ev$dir_flipped))
  expect_false(any(ev$dir_flipped & ev$dir_source == "sideA")); expect_false(any(ev$dir_flipped & ev$hostlev_disp != 5L))
})

## ---- participant-level dates (PI decision 2026-10-25, D5) --------------------------------------------------------------------------
test_that("D5: participant start / end months give exact pair entry and exit; inconsistent participants fall back; entry bounds are clipped", {
  ## dispute 90: 10 (A, orig) vs 20 (B, orig) 1990-1991, 400 days from 1990-03; 30 joins Side B in 1990 (orig 0). Participant dates:
  ## 10 and 20 1990-03 .. 1991-04; 30 1990-09 .. 1991-04 -> pair (10, 30) entry 1990-09 exactly (was [1990-03, 1990-12])
  dd <- bind_rows(mk_dd(90L, 10L, 20L, 1990:1991), mk_dd(90L, 10L, 30L, 1990:1991, orig_b = 0L), mk_dd(91L, 10L, 30L, 1990L), mk_dd(92L, 10L, 30L, 1991L))
  disp <- bind_rows(mk_disp(90L, 0L, 3L, 1990L, 3L, 400), mk_disp(91L, 0L, 3L, 1990L, 7L, 5), mk_disp(92L, 0L, 3L, 1991L, 2L, 5))
  part <- tibble(dispnum = c(90L, 90L, 90L, 91L, 91L), ccode = c(10L, 20L, 30L, 10L, 30L),
                 styear = c(1990L, 1990L, 1990L, 1990L, 1990L), stmon = c(3L, 3L, 9L, 7L, 7L), stday = -9L,
                 endyear = c(1991L, 1991L, 1991L, 1990L, 1990L), endmon = c(4L, 4L, 4L, 7L, 7L), endday = -9L)
  db <- dispute_date_bounds(disp, dd)
  p0 <- pair_dispute_intervals(dd, db); p1 <- pair_dispute_intervals(dd, db, part = part)
  j0 <- p0 |> filter(dispnum == 90L, ccode1 == 10L, ccode2 == 30L); j1 <- p1 |> filter(dispnum == 90L, ccode1 == 10L, ccode2 == 30L)
  expect_equal(c(j0$start_lo, j0$start_hi), c(ym_index(1990L, 3L), ym_index(1990L, 12L))); expect_equal(j0$date_source, "dispute")
  expect_equal(c(j1$start_lo, j1$start_hi), rep(ym_index(1990L, 9L), 2)); expect_equal(c(j1$end_lo, j1$end_hi), rep(ym_index(1991L, 4L), 2))
  expect_true(j1$start_exact); expect_equal(j1$date_source, "participant"); expect_equal(j1$part_reason, "participant dates used")
  o1 <- p1 |> filter(dispnum == 90L, ccode1 == 10L, ccode2 == 20L); expect_equal(c(o1$start_lo, o1$end_hi), c(ym_index(1990L, 3L), ym_index(1991L, 4L)))
  ## dispute 92 is not in part: dispute-level bounds stay; both directions of a pair get the same dates
  expect_equal(p1$date_source[p1$dispnum == 92L], c("dispute", "dispute")); expect_equal(p1$part_reason[p1$dispnum == 92L][1], "dispute not in gml_part")
  expect_equal(p1 |> filter(dispnum == 90L, ccode1 == 30L, ccode2 == 10L) |> pull(start_lo), j1$start_lo)
  ## the escalation decision follows: 91 (July 1990) vs the joiner pair's dispute 90 was undecided (year fallback -> not an escalation);
  ## with the entry known to be September, July is certainly BEFORE the entry -> certainly not an escalation, no fallback
  ev0 <- code_events(dd, disp, include_joiners = TRUE); ev1 <- code_events(dd, disp, include_joiners = TRUE, part = part)
  expect_true(ev0$esc_fallback[ev0$dispnum == 91L]); expect_false(ev1$esc_fallback[ev1$dispnum == 91L]); expect_false(ev1$escalation[ev1$dispnum == 91L])
  expect_true(ev1$escalation[ev1$dispnum == 92L])                                             # February 1991: entered 1990-09, still active (ends 1991-04)
  ## a participant whose GML entry is AFTER the pair's first COW year contradicts the COW years: dispute bounds, reason recorded
  part_bad <- part |> mutate(styear = ifelse(ccode == 30L & dispnum == 90L, 1991L, styear), stmon = ifelse(ccode == 30L & dispnum == 90L, 2L, stmon))
  pb <- pair_dispute_intervals(dd, db, part = part_bad) |> filter(dispnum == 90L, ccode1 == 10L, ccode2 == 30L)
  expect_equal(pb$date_source, "dispute"); expect_match(pb$part_reason, "entry after the spell's first COW year")
  expect_equal(c(pb$start_lo, pb$start_hi), c(ym_index(1990L, 3L), ym_index(1990L, 12L)))
  ## a GML entry BEFORE the first COW year is clamped to January of that year (the pair was in the dispute by then)
  part_early <- part |> mutate(styear = ifelse(ccode == 30L & dispnum == 90L, 1989L, styear))
  pe <- pair_dispute_intervals(dd, db, part = part_early) |> filter(dispnum == 90L, ccode1 == 10L, ccode2 == 30L)
  expect_equal(pe$date_source, "participant"); expect_equal(c(pe$start_lo, pe$start_hi), rep(ym_index(1990L, 3L), 2))   # the pair entry is the LATER start, 10's 1990-03, so nothing to clamp
  part_both <- part |> mutate(styear = ifelse(dispnum == 90L, 1989L, styear))
  pc2 <- pair_dispute_intervals(dd, db, part = part_both) |> filter(dispnum == 90L, ccode1 == 10L, ccode2 == 20L)
  expect_match(pc2$part_reason, "entry clamped"); expect_equal(c(pc2$start_lo, pc2$start_hi), rep(ym_index(1990L, 1L), 2))   # both GML starts in 1989: clamped to January 1990
  ## a participant with no overlapping gml_part spell, and one with two spells, fall back too
  pn <- pair_dispute_intervals(dd, db, part = part |> filter(!(dispnum == 90L & ccode == 30L))) |> filter(dispnum == 90L, ccode1 == 10L, ccode2 == 30L)
  expect_match(pn$part_reason, "no gml_part spell")
  p2 <- pair_dispute_intervals(dd, db, part = bind_rows(part, part |> filter(dispnum == 90L, ccode == 30L) |> mutate(stmon = 11L))) |> filter(dispnum == 90L, ccode1 == 10L, ccode2 == 30L)
  expect_match(p2$part_reason, "several gml_part spells")
  ## clipping: a joiner pair entering some time in 1990 when the dispute ends in 1990-05 had start_hi (December) > end_hi (May)
  dd_c <- bind_rows(mk_dd(95L, 10L, 20L, 1990L), mk_dd(95L, 10L, 30L, 1990L, orig_b = 0L)); disp_c <- mk_disp(95L, 0L, 3L, 1990L, 2L, 60)
  pc <- pair_dispute_intervals(dd_c, dispute_date_bounds(disp_c, dd_c)) |> filter(ccode1 == 10L, ccode2 == 30L)
  expect_true(pc$start_clipped); expect_equal(pc$start_hi, pc$end_hi); expect_true(pc$start_lo <= pc$start_hi); expect_lt(pc$start_hi, ym_index(1990L, 12L))
  expect_false(any(pair_dispute_intervals(dd_c, dispute_date_bounds(disp_c, dd_c))$start_hi > pair_dispute_intervals(dd_c, dispute_date_bounds(disp_c, dd_c))$end_hi))
})

test_that("D5 on the real data: Egypt <-> Israel 1967 (disputes 1035 and 1480) get participant-level dates and the right escalation call", {
  skip_if_not_installed("peacesciencer")
  suppressPackageStartupMessages(library(peacesciencer))
  dd <- cow_mid_dirdisps |> mutate(across(c(dispnum, ccode1, ccode2, year), as.integer)) |> filter(dispnum %in% c(1035L, 1480L))
  disp <- cow_mid_disps |> mutate(dispnum = as.integer(dispnum)) |> filter(dispnum %in% c(1035L, 1480L))
  pdi <- pair_dispute_intervals(dd, dispute_date_bounds(disp, dd), part = gml_part)
  ## 1035 (the Six-Day War): Egypt 651 and Israel 666 are 1966-12 .. 1967-06 in gml_part but their COW pair-years begin in 1967, so the
  ## entry is clamped to January 1967 (exact; under the dispute-level bounds the pair, a later entrant, had the whole of 1967) and the
  ## exit is the exact 1967-06 (was the duration bound)
  p <- pdi |> filter(dispnum == 1035L, ccode1 == 651L, ccode2 == 666L)
  expect_equal(p$date_source, "participant"); expect_match(p$part_reason, "entry clamped")
  expect_equal(c(p$start_lo, p$start_hi), rep(ym_index(1967L, 1L), 2)); expect_equal(c(p$end_lo, p$end_hi), rep(ym_index(1967L, 6L), 2))
  p0 <- pair_dispute_intervals(dd, dispute_date_bounds(disp, dd)) |> filter(dispnum == 1035L, ccode1 == 651L, ccode2 == 666L)
  expect_equal(p0$start_lo, ym_index(1967L, 1L)); expect_equal(p0$start_hi, p0$end_hi); expect_true(p0$start_clipped); expect_lt(p0$end_lo, p0$end_hi)
  ## Jordan 663 joined 1967-05: the pair (663, 666) enters in May 1967, exits June
  pj <- pdi |> filter(dispnum == 1035L, ccode1 == 663L, ccode2 == 666L)
  expect_equal(c(pj$start_lo, pj$start_hi), rep(ym_index(1967L, 5L), 2)); expect_equal(pj$end_hi, ym_index(1967L, 6L)); expect_equal(pj$date_source, "participant")
  ## 1480 (the War of Attrition): Israel and Egypt 1967-06-29 .. 1970-08 -> pair entry 1967-06; the USSR's 1969 spell gives the (365, 666) pair a 1969 entry
  pw <- pdi |> filter(dispnum == 1480L, ccode1 == 666L, ccode2 == 651L)
  expect_equal(pw$date_source, "participant"); expect_equal(c(pw$start_lo, pw$start_hi), rep(ym_index(1967L, 6L), 2)); expect_equal(pw$end_hi, ym_index(1970L, 8L))
  pu <- pdi |> filter(dispnum == 1480L, ccode1 == 365L, ccode2 == 666L); expect_equal(pu$spell_first, 1970L); expect_equal(pu$date_source, "dispute")   # the USSR's gml_part spell (1969-04) does not overlap its COW year 1970
  ## Israel -> Egypt 1967 in 1480 (GML start June 29) against 1035 (GML end June 10 / 17): at MONTH resolution both are June 1967, a tie,
  ## and the same-month convention counts an end in the reference month as still active, so the event is an escalation WITHOUT a fallback.
  ## Under the dispute-level bounds the comparison was undecided and the year rule (onset 1967 not before 1967) said not an escalation.
  ## The day fields of gml_part would separate the two; the pipeline works at month resolution (documented limitation, README).
  ev_pairs <- build_pair_table(dd, disp |> transmute(dispnum, fatality_disp = fatality, hostlev_disp = hostlev)) |> filter(sidea1 == 1L, sidea2 == 0L) |> rename(initiator = ccode1, target = ccode2)
  e1 <- flag_escalations_month(ev_pairs, pdi); e0 <- flag_escalations_month(ev_pairs, pair_dispute_intervals(dd, dispute_date_bounds(disp, dd)))
  k <- e1$dispnum == 1480L & e1$initiator == 666L & e1$target == 651L
  expect_true(any(k)); expect_false(e1$esc_fallback[k]); expect_true(e1$escalation[k]); expect_true(e0$esc_fallback[k]); expect_false(e0$escalation[k])
  expect_equal(c(e1$ev_start_lo[k], e1$ev_start_hi[k]), rep(ym_index(1967L, 6L), 2))
  expect_false(any(e1$escalation[e1$dispnum == 1035L]))
})

## ---- stage-2 event weights (PI decision 2026-10-25, D8) --------------------------------------------------------------------------------
test_that("D8: w_<sample> = 1 / the initiator's targets in the same dispute and year, per sample; the collapse sums and caps at 1", {
  ## dispute 200: 10 (A) vs 20, 21, 22 (B) in 1990, hostlev 5 but only (10, 20) at participant level 5 -> war sample has one target (w_war 1),
  ## mid has three (w_mid 1/3); dispute 201: 10 -> 20 again in 1990 (collapses with the 200 pair: w_mid 1/3 + 1 capped at 1)
  dd <- bind_rows(mk_dd(200L, 10L, 20L, 1990L) |> mutate(hostlev1 = 5L, hostlev2 = 5L), mk_dd(200L, 10L, 21L, 1990L), mk_dd(200L, 10L, 22L, 1990L),
                  mk_dd(201L, 10L, 20L, 1990L), mk_dd(202L, 30L, 40L, 1991L), mk_dd(202L, 31L, 40L, 1991L))
  disp <- bind_rows(mk_disp(200L, 1L, 5L, 1990L, 3L, 100), mk_disp(201L, 0L, 3L, 1990L, 8L, 5), mk_disp(202L, 0L, 3L, 1991L, 1L, 5))
  ev <- code_events(dd, disp, include_joiners = TRUE)
  e200 <- ev |> filter(dispnum == 200L)
  expect_equal(e200$w_mid, rep(1 / 3, 3)); expect_equal(sum(e200$w_mid), 1)
  expect_equal(e200$w_war[e200$target == 20L], 1); expect_equal(e200$w_war[e200$target != 20L], c(0, 0))     # outside the war sample: 0
  expect_equal(ev$w_mid[ev$dispnum == 201L], 1)
  expect_equal(ev$w_mid[ev$dispnum == 202L], c(1, 1))                      # two initiators against one target: one decision each
  expect_true(all(ev$w_fmid[!ev$in_fmid] == 0)); expect_true(all(ev$w_mid[ev$in_mid] > 0))
  ## escalations carry 0 everywhere
  dd_e <- bind_rows(dd, mk_dd(203L, 10L, 21L, 1990L)); disp_e <- bind_rows(disp, mk_disp(203L, 0L, 3L, 1990L, 4L, 5))   # April, while 200 (March, 100 days) is active
  ev_e <- code_events(dd_e, disp_e, include_joiners = TRUE)
  expect_true(ev_e$escalation[ev_e$dispnum == 203L]); expect_equal(ev_e$w_mid[ev_e$dispnum == 203L], 0)
  ## the collapse: (10, 20, 1990) from 200 (1/3) and 201 (1) -> min(1, 4/3) = 1; (10, 21) and (10, 22) keep 1/3
  out <- collapse_duplicate_events(ev |> add_collapse_cols())
  k <- out |> filter(initiator == 10L, target == 20L); expect_equal(k$n_disputes, 2L); expect_equal(k$w_mid, 1); expect_equal(k$w_war, 1)
  expect_equal(out$w_mid[out$initiator == 10L & out$target == 21L], 1 / 3)
  expect_equal(out$dir_flipped, rep(FALSE, nrow(out))); expect_true(all(out$dir_source == "sideA"))
  expect_error(collapse_duplicate_events(ev |> add_collapse_cols() |> mutate(w_mid = 1.5)), "outside")
  expect_error(collapse_duplicate_events(ev |> add_collapse_cols() |> select(-w_war)), "absent")
  expect_error(event_weights(ev |> select(-in_war)))
})

test_that("P3: the pair table takes the sides of the pair's EARLIEST year and reports side changes", {
  ## pair (10, 20): 10 on Side A in 1939 and on Side B in 1941 (opposed again after switching sides, like France and the USSR in 258)
  dd <- bind_rows(mk_dd(300L, 10L, 20L, 1939L), mk_dd(300L, 20L, 10L, 1941L) |> mutate(disponset = 0L))
  pair <- build_pair_table(dd, mk_disp(300L, 5L, 5L, 1939L, 9L, 2000))
  p <- pair |> filter(ccode1 == 10L, ccode2 == 20L)
  expect_equal(p$year, 1939L); expect_equal(c(p$sidea1, p$sidea2), c(1L, 0L))                   # the 1939 sides, whatever the row order
  pair_rev <- build_pair_table(dd |> arrange(desc(year)), mk_disp(300L, 5L, 5L, 1939L, 9L, 2000)) |> filter(ccode1 == 10L, ccode2 == 20L)
  expect_equal(c(pair_rev$sidea1, pair_rev$sidea2), c(1L, 0L))
  expect_equal(nrow(attr(pair, "side_changes")), 2L)                                               # both directions of the one pair
  expect_equal(nrow(attr(build_pair_table(toy_dd, toy_disp), "side_changes")), 0L)
  expect_error(build_pair_table(dd |> mutate(orig1 = ifelse(year == 1941L, 0L, orig1)), mk_disp(300L, 5L, 5L, 1939L, 9L, 2000)), "originator status")
})

## ---- switches of spec.R at assembly (2026-10-25): pairs = originators | joiners, changes = winsor -----------------------------------
test_that("pairs switch: originators / joiners filter the events and recompute y1 from the kept events; all = unchanged", {
  ## 10 -> 30 (dispute 4;5), 30 -> 10 (dispute 6) and 40 -> 20 (dispute 7) become joiner pairs
  joiners <- function(e) e |> mutate(is_joiner_pair = dispnum %in% c(4L, 6L, 7L))
  td <- write_toy_cell(tempfile("toycell"), events_edit = joiners)
  A  <- quiet(assemble_cell("mid", 3L, "all", data_dir = td))
  Ao <- quiet(assemble_cell("mid", 3L, "all", PAIRS = "originators", data_dir = td))
  Aj <- quiet(assemble_cell("mid", 3L, "all", PAIRS = "joiners", data_dir = td))
  expect_equal(A$sy$y1, c(1L, 0L, 1L, 1L)); expect_equal(nrow(A$ev), 5L); expect_equal(A$pairs, "all")
  ## originators: 10 -> 20, 30 -> 20 usable, 30 -> 40 stage 1 only; state 40 (only a joiner event) gets y1 = 0
  expect_equal(Ao$sy$y1, c(1L, 0L, 1L, 0L)); expect_equal(nrow(Ao$ev), 2L); expect_equal(Ao$n_fighting_at_ref, 1L)
  chosen_o <- Ao$alts[Ao$ev$alt_row, ]; expect_setequal(paste(chosen_o$ccode1, chosen_o$ccode2), c("10 20", "30 20"))
  ## joiners: 10 -> 30, 30 -> 10, 40 -> 20; state 20 never initiates
  expect_equal(Aj$sy$y1, c(1L, 0L, 1L, 1L)); expect_equal(nrow(Aj$ev), 3L); expect_equal(Aj$n_fighting_at_ref, 0L)
  chosen_j <- Aj$alts[Aj$ev$alt_row, ]; expect_setequal(paste(chosen_j$ccode1, chosen_j$ccode2), c("10 30", "30 10", "40 20"))
  expect_true(any(Aj$ev$w == 0.5))                                             # the weights travel with the kept events
  ## x_ref follows the recomputed y1 (mean X over the initiating state-years of the kept group)
  expect_equal(unname(Ao$x_ref["dM_i"]), mean(c(0.5, 0.1))); expect_equal(unname(Aj$x_ref["dM_i"]), mean(c(0.5, 0.1, -0.4)))
  ## alternatives, W / X / Z and the Z centring do not depend on the pairs switch
  expect_equal(Ao$alts$dM_j_h, A$alts$dM_j_h); expect_equal(Ao$W_cols, A$W_cols); expect_equal(Ao$Z_cols, A$Z_cols)
  ## the stage-1 data of both variants pass build_stan_data(), and a group that leaves no event stops
  dat_o <- build_stan_data(Ao$sy, Ao$alts, Ao$ev, Ao$W_cols, Ao$X_cols, Ao$Z_cols, joint = 1L, x_ref = Ao$x_ref, allow_orphan_initiations = TRUE)
  expect_equal(sum(dat_o$y1), 2L); expect_equal(dat_o$E, 2L)
  td0 <- write_toy_cell(tempfile("toycell"))                                   # no joiner pairs at all
  expect_error(quiet(assemble_cell("mid", 3L, "all", PAIRS = "joiners", data_dir = td0)), "leaves no event")
  expect_error(quiet(assemble_cell("mid", 3L, "all", PAIRS = "both", data_dir = td0)))
})

## ---- D8 weights under pairs != all (PI decision 2026-10-26): split groups rescaled to sum 1 --------------------------------------------
## toy events for the pairs x D8 interaction (all 1990, mid): dispute 1: 10 -> 20 (originator, collapsed with dispute 9, w 1), 10 -> 30 and
## 10 -> 40 (joiners, 1/3 each; the group sums to 5/3 under pairs = all); dispute 2: 30 -> 20 (originator, 1/2), 30 -> 10 (joiner, 1/2);
## dispute 7: 40 -> 20 (joiner, 1); dispute 8: 40 -> 10 (originator, 1). Groups (10, 1), (30, 2) are split by both filters; (40, 7) and
## (40, 8) are kept whole by one filter and dropped whole by the other.
pairs_d8_events <- function(e) {
  stopifnot(e$initiator[1] == 10L, e$target[1] == 20L, e$dispnum[1] == 1L)           # the template row: a plain dispute-1 event
  e[rep(1L, 7L), ] |>
    mutate(initiator = c(10L, 10L, 10L, 30L, 30L, 40L, 40L), target = c(20L, 30L, 40L, 20L, 10L, 20L, 10L),
           dispnum = c(1L, 1L, 1L, 2L, 2L, 7L, 8L), dispnums = c("1;9", "1", "1", "2", "2", "7", "8"), n_disputes = c(2L, 1L, 1L, 1L, 1L, 1L, 1L),
           is_joiner_pair = c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE),
           w_mid = c(1, 1 / 3, 1 / 3, 0.5, 0.5, 1, 1), w_fmid = w_mid, w_fmidub = w_mid)
}
test_that("D8 x pairs (helpers): split groups are rescaled to sum 1, whole groups keep the build's weights, a collapsed key follows its dispnum", {
  td <- write_toy_cell(tempfile("toycell"), events_edit = pairs_d8_events)
  before <- readRDS(file.path(td, "events.rds"))
  key <- function(x) paste(x$initiator, x$target)
  ## originators: (10, 1) keeps the collapsed key alone (kept sum 1, factor 1), (30, 2) keeps 30 -> 20 (0.5 -> 1), (40, 8) whole, (40, 7) dropped
  ko <- before |> filter(!is_joiner_pair); ro <- d8_rescale_split_groups(before, ko, "w_mid")
  expect_equal(ro$w_mid[key(ro) == "10 20"], 1); expect_equal(ro$w_mid[key(ro) == "30 20"], 1); expect_equal(ro$w_mid[key(ro) == "40 10"], 1)
  expect_equal(attr(ro, "n_groups_split"), 2L); expect_equal(attr(ro, "n_groups_rescaled"), 1L); expect_equal(attr(ro, "n_events_rescaled"), 1L)
  expect_equal(attr(ro, "w_before"), 2.5); expect_equal(attr(ro, "w_after"), 3); expect_equal(attr(ro, "n_split_with_collapsed"), 1L)
  ## joiners: (10, 1) keeps 10 -> 30 and 10 -> 40 (1/3 + 1/3 -> 1/2 each), (30, 2) keeps 30 -> 10 (0.5 -> 1), (40, 7) whole at 1
  kj <- before |> filter(is_joiner_pair); rj <- d8_rescale_split_groups(before, kj, "w_mid")
  expect_equal(rj$w_mid[key(rj) == "10 30"], 0.5); expect_equal(rj$w_mid[key(rj) == "10 40"], 0.5); expect_equal(rj$w_mid[key(rj) == "30 10"], 1); expect_equal(rj$w_mid[key(rj) == "40 20"], 1)
  expect_equal(attr(rj, "n_groups_split"), 2L); expect_equal(attr(rj, "n_groups_rescaled"), 2L); expect_equal(attr(rj, "n_events_rescaled"), 3L)
  g <- rj |> group_by(initiator, dispnum, year) |> summarise(s = sum(w_mid), .groups = "drop"); expect_true(all(abs(g$s - 1) < 1e-12))
  ## the other weight columns are untouched; the row set and the column set are those of the kept frame
  expect_equal(rj$w_fmid, kj$w_fmid); expect_equal(names(rj), names(kj)); expect_equal(nrow(rj), nrow(kj))
  ## nothing filtered: no split group, weights identical (the (10, 1) group keeps its sum of 5/3)
  r0 <- d8_rescale_split_groups(before, before, "w_mid"); expect_equal(r0$w_mid, before$w_mid); expect_equal(attr(r0, "n_groups_split"), 0L)
  ## the group table and the whole-group assertion: dropping dispute 2 whole passes; dropping one single-dispute member of dispute 1 stops;
  ## dropping the collapsed key alone (as the no-WWII flag can) is a split through a collapsed key, counted and not an error
  g0 <- d8_group_table(before, before, "w_mid"); expect_equal(nrow(g0), 4L); expect_equal(sort(g0$w_before), c(1, 1, 1, 5 / 3))
  ok <- d8_assert_groups_whole(before, before |> filter(dispnum != 2L), "w_mid", "toy"); expect_equal(c(ok$n_whole_kept, ok$n_whole_dropped, ok$n_split_collapsed), c(3L, 1L, 0L))
  expect_error(d8_assert_groups_whole(before, before |> filter(!(initiator == 10L & target == 30L)), "w_mid", "toy"), "split among their single-dispute members")
  cl <- d8_assert_groups_whole(before, before |> filter(!(initiator == 10L & target == 20L)), "w_mid", "toy"); expect_equal(cl$n_split_collapsed, 1L)
  expect_error(d8_group_table(before |> filter(dispnum != 7L), before, "w_mid"), "not a subset")
})

test_that("D8 x pairs (assembly): originators / joiners rescale the split groups to sum 1; pairs = all and weights = none are unchanged", {
  td <- write_toy_cell(tempfile("toycell"), events_edit = pairs_d8_events)
  A  <- quiet(assemble_cell("mid", 3L, "all", data_dir = td))
  Ao <- quiet(assemble_cell("mid", 3L, "all", PAIRS = "originators", data_dir = td))
  Aj <- quiet(assemble_cell("mid", 3L, "all", PAIRS = "joiners", data_dir = td))
  A1 <- quiet(assemble_cell("mid", 3L, "all", PAIRS = "joiners", EVENT_WEIGHTS = "none", data_dir = td))
  chosen <- function(X) paste(X$alts$ccode1[X$ev$alt_row], X$alts$ccode2[X$ev$alt_row])
  ## pairs = all: the build's weights, the (10, 1) group summing to 5/3 included
  expect_equal(nrow(A$ev), 7L); expect_equal(sum(A$ev$w), 14 / 3); expect_equal(A$pairs_weights$n_groups_split, 0L)
  expect_equal(A$ev$w[chosen(A) == "10 30"], 1 / 3); expect_equal(A$ev$w[chosen(A) == "30 10"], 0.5)
  ## originators: 10 -> 20 (1), 30 -> 20 (0.5 -> 1), 40 -> 10 (1)
  expect_setequal(chosen(Ao), c("10 20", "30 20", "40 10")); expect_equal(Ao$ev$w, rep(1, 3)); expect_equal(Ao$sy$y1, c(1L, 0L, 1L, 1L))
  expect_equal(Ao$pairs_weights$n_groups_split, 2L); expect_equal(Ao$pairs_weights$n_groups_rescaled, 1L); expect_equal(Ao$pairs_weights$w_before, 2.5); expect_equal(Ao$pairs_weights$w_after, 3)
  ## joiners: 10 -> 30 and 10 -> 40 (1/3 -> 1/2), 30 -> 10 (0.5 -> 1), 40 -> 20 (1)
  expect_setequal(chosen(Aj), c("10 30", "10 40", "30 10", "40 20"))
  expect_equal(Aj$ev$w[chosen(Aj) == "10 30"], 0.5); expect_equal(Aj$ev$w[chosen(Aj) == "10 40"], 0.5); expect_equal(Aj$ev$w[chosen(Aj) == "30 10"], 1); expect_equal(Aj$ev$w[chosen(Aj) == "40 20"], 1)
  expect_equal(Aj$pairs_weights$n_groups_rescaled, 2L); expect_equal(Aj$pairs_weights$n_events_rescaled, 3L); expect_equal(Aj$pairs_weights$w_after, 3)
  ## weights = none: weight 1 everywhere, no rescale reported
  expect_true(all(A1$ev$w == 1)); expect_equal(A1$pairs_weights$n_groups_split, 0L); expect_true(is.na(A1$pairs_weights$w_before))
  ## the rescale does not touch the alternatives, the stage-1 frame beyond y1, or the chosen rows
  expect_equal(Aj$alts, A$alts); expect_equal(Aj$sy |> select(-y1), A$sy |> select(-y1)); expect_equal(sort(Aj$ev$alt_row), sort(A$ev$alt_row[chosen(A) %in% chosen(Aj)]))
  ## ww2 = drop on the toy drops disputes 6 and 7 whole (the default fixture); the group assertion passes and reports no split
  td0 <- write_toy_cell(tempfile("toycell"))
  expect_error(quiet(assemble_cell("mid", 3L, "all", DROP_WW2 = TRUE, data_dir = td0)), NA)
})

test_that("D8 x pairs (real data, when present): fmid lag 3 split groups sum to 1 under originators and joiners; pairs = all carries the build's weights", {
  man_path <- file.path(v2_root(), "data", "manifest.csv")
  skip_if_not(file.exists(man_path) && file.exists(file.path(v2_root(), "data", "alts_lag3_all.rds")), "real data files not built")
  events <- readRDS(file.path(v2_root(), "data", "events.rds")); e <- events |> filter(in_fmid)
  chosen_keys <- function(X) X$alts[X$ev$alt_row, c("ccode1", "ccode2", "year")] |> mutate(w = X$ev$w)
  check <- function(PAIRS) {
    X <- quiet(assemble_cell_spec("fmid", 3L, parse_spec(paste0("--pairs=", PAIRS))))
    kept <- if (PAIRS == "originators") e |> filter(!is_joiner_pair) else e |> filter(is_joiner_pair)
    g <- d8_group_table(e, kept, "w_fmid"); split <- g |> filter(n_kept > 0L, n_kept < n)
    ## the Stan-data weights of the kept events, matched back to the events file, sum to 1 in every split group whose members are all usable
    ck <- chosen_keys(X) |> inner_join(kept |> select(ccode1 = initiator, ccode2 = target, year, dispnum, w_build = w_fmid), by = c("ccode1", "ccode2", "year"))
    expect_equal(nrow(ck), nrow(X$ev))
    s <- ck |> group_by(ccode1, dispnum, year) |> summarise(n_usable = n(), s = sum(w), .groups = "drop") |>
      inner_join(split |> transmute(ccode1 = initiator, dispnum, year, n_kept), by = c("ccode1", "dispnum", "year"))
    expect_true(nrow(s) > 0); expect_true(all(abs(s$s[s$n_usable == s$n_kept] - 1) < 1e-9))
    ## untouched groups carry the build's weights
    whole <- ck |> anti_join(split |> transmute(ccode1 = initiator, dispnum, year), by = c("ccode1", "dispnum", "year"))
    expect_equal(whole$w, whole$w_build)
    expect_equal(X$pairs_weights$n_groups_split, nrow(split)); expect_equal(X$pairs_weights$n_groups, nrow(g))
    X
  }
  Xo <- check("originators"); Xj <- check("joiners")
  expect_equal(Xo$pairs_weights$n_groups_split, 49L); expect_equal(Xj$pairs_weights$n_groups_split, 49L)
  expect_equal(nrow(Xo$ev), 545L); expect_equal(round(sum(Xo$ev$w), 1), 529.1); expect_equal(round(Xo$pairs_weights$w_after, 1), 580.1)
  expect_equal(nrow(Xj$ev), 510L); expect_equal(round(sum(Xj$ev$w), 1), 253.2); expect_equal(round(Xj$pairs_weights$w_after, 1), 266.5)
  ## pairs = all: every usable event carries its build weight (the main-specification Stan data are unchanged)
  Xa <- quiet(assemble_cell_spec("fmid", 3L, parse_spec(character(0))))
  ca <- chosen_keys(Xa) |> inner_join(e |> select(ccode1 = initiator, ccode2 = target, year, w_build = w_fmid), by = c("ccode1", "ccode2", "year"))
  expect_equal(nrow(ca), 1055L); expect_equal(ca$w, ca$w_build); expect_equal(round(sum(ca$w), 1), 741.0)
})

test_that("changes switch: winsor clips dM_i, dR_i, dM_j_h, dRpr_j_h at the panel's 1st / 99th percentiles before the Z centring; raw = unchanged", {
  td <- write_toy_cell(tempfile("toycell"))
  A  <- quiet(assemble_cell("mid", 3L, "all", data_dir = td))
  Aw <- quiet(assemble_cell("mid", 3L, "all", CHANGES = "winsor", data_dir = td))
  panel <- readRDS(file.path(td, "stage1_panel_lag3.rds")); alts0 <- readRDS(file.path(td, "alts_lag3_all.rds"))
  capM <- unname(quantile(panel$dM_h, c(0.01, 0.99))); capR <- unname(quantile(panel$dRpr_h, c(0.01, 0.99)))
  expect_null(A$winsor); expect_equal(Aw$winsor$M, capM); expect_equal(Aw$winsor$R, capR); expect_equal(Aw$changes, "winsor")
  ## initiator changes: the extremes are clipped, the interior values untouched
  expect_equal(Aw$sy$dM_i, pmin(pmax(panel$dM_h, capM[1]), capM[2])); expect_equal(Aw$sy$dR_i, pmin(pmax(panel$dRpr_h, capR[1]), capR[2]))
  expect_equal(sum(Aw$sy$dM_i != A$sy$dM_i), 2L)                              # 0.5 and -0.4 are the 4-point extremes
  expect_true(all(Aw$sy$dM_i >= capM[1] & Aw$sy$dM_i <= capM[2]))
  ## target changes: clipped at the SAME caps, then centred over the IV rows (so compare the de-centred values)
  key <- function(a) paste(a$ccode1, a$ccode2)
  raw_j <- alts0$dM_j_h[match(key(Aw$alts), key(alts0))]
  clipped <- pmin(pmax(raw_j, capM[1]), capM[2])
  expect_equal(Aw$alts$dM_j_h, clipped - mean(clipped[Aw$alts$in_iv == 1L]))
  expect_true(any(raw_j < capM[1] | raw_j > capM[2]))                          # the toy alternatives do have values outside the caps
  expect_lt(abs(mean(Aw$alts$dM_j_h[Aw$alts$in_iv == 1L])), 1e-12)
  ## a target without a change value keeps 0 before the centring: its centred value is minus the IV-row mean of the clipped column
  miss <- Aw$alts$dch_j_missing == 1L
  expect_true(all(abs(Aw$alts$dM_j_h[miss] - (0 - mean(clipped[Aw$alts$in_iv == 1L]))) < 1e-12))
  ## everything else is identical to the raw assembly
  expect_equal(Aw$sy$y1, A$sy$y1); expect_equal(Aw$ev, A$ev); expect_equal(Aw$Z_cols, A$Z_cols); expect_equal(Aw$alts$contig_land, A$alts$contig_land)
  expect_error(quiet(assemble_cell("mid", 3L, "all", CHANGES = "trim", data_dir = td)))
})

test_that("assemble_cell_spec() maps the switches of spec.R onto assemble_cell() and the main spec is the plain call", {
  td <- write_toy_cell(tempfile("toycell"))
  main <- parse_spec(character(0))
  A <- quiet(assemble_cell("mid", 3L, "all", data_dir = td)); B <- quiet(assemble_cell_spec("mid", 3L, main, data_dir = td))
  expect_equal(A$sy, B$sy); expect_equal(A$alts, B$alts); expect_equal(A$ev, B$ev); expect_equal(A$W_cols, B$W_cols); expect_equal(A$Z_cols, B$Z_cols)
  C <- quiet(assemble_cell_spec("mid", 3L, parse_spec(c("--ww2=drop", "--weights=none", "--controls=full")), data_dir = td))
  D <- quiet(assemble_cell("mid", 3L, "all", DROP_WW2 = TRUE, EVENT_WEIGHTS = "none", FULL_CONTROLS = TRUE, data_dir = td))
  expect_equal(C$sy, D$sy); expect_equal(C$ev, D$ev); expect_equal(C$W_cols, D$W_cols); expect_equal(C$Z_cols, D$Z_cols)
  expect_true(all(C$ev$w == 1)); expect_equal(C$full_controls, TRUE); expect_equal(C$drop_ww2, TRUE)
})
