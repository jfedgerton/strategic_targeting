## Unit tests of the data-coding logic in 01_build_data.R using small hand-built toy frames
## that mimic the peacesciencer columns. The full build needs the bundled peacesciencer data
## and is run by Rscript v2/R/01_build_data.R (on Roar or locally); these tests exercise the
## rules in isolation so a regression in the coding logic fails loudly without the full build.
suppressPackageStartupMessages({ library(testthat); library(dplyr); library(tidyr) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "00_helpers.R"))

## ---- the coding rules under test, lifted verbatim from 01_build_data.R ---------------------
code_events <- function(dd, disp, include_joiners = FALSE,
                        mic = tibble(dispnum = integer(), mic_fatalmin = integer()),
                        gml = tibble(dispnum = integer(), initiator = integer(), target = integer(), fatal_gml = logical())) {
  pair <- dd |> group_by(dispnum, ccode1, ccode2) |>
    summarise(year = min(year), sidea1 = first(sidea1), sidea2 = first(sidea2),
              orig1 = first(orig1), orig2 = first(orig2),
              fatality1 = max(fatality1), fatality2 = max(fatality2),
              hostlev1 = max(hostlev1), hostlev2 = max(hostlev2), .groups = "drop") |>
    left_join(disp, by = "dispnum")
  events <- pair |> filter(sidea1 == 1L, sidea2 == 0L)
  if (!include_joiners) events <- events |> filter(orig1 == 1L, orig2 == 1L)
  events <- events |> rename(initiator = ccode1, target = ccode2) |>
    mutate(is_joiner_pair = orig1 == 0L | orig2 == 0L,
           fatal1 = ifelse(fatality1 < 0, NA, fatality1 >= 1),
           fatal2 = ifelse(fatality2 < 0, NA, fatality2 >= 1),
           fatal_disp = ifelse(fatality_disp < 0, NA, fatality_disp >= 1),
           fatal = case_when(fatal1 %in% TRUE | fatal2 %in% TRUE | fatal_disp %in% TRUE ~ TRUE,
                             is.na(fatal1) & is.na(fatal2) & is.na(fatal_disp) ~ NA,
                             TRUE ~ FALSE),
           war = hostlev_disp == 5L)
  ## ongoing-dispute rule (2026-10-24): escalation = already fighting in another dispute at t; excluded from every sample flag
  events <- flag_escalations(events, already_fighting_keys(dd))
  events |>
    left_join(mic, by = "dispnum") |> left_join(gml, by = c("dispnum", "initiator", "target")) |>
    mutate(fatal_mic = ifelse(is.na(mic_fatalmin), NA, mic_fatalmin >= 1),
           fatal_main = code_fatal_main(fatal, fatal_mic, war, fatal_gml),
           in_mid = !escalation,
           in_fmid = fatal_main & !escalation,
           in_fmidub = (fatal_main | (is.na(fatal) & is.na(fatal_mic) & is.na(fatal_gml))) & !escalation,
           in_fmidcow = (fatal %in% TRUE) & !escalation,
           in_war = (war %in% TRUE) & !escalation)
}

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
toy_disp <- tribble(~dispnum, ~fatality_disp, ~hostlev_disp,
                    1L, 0L, 4L,
                    2L, 1L, 5L)

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
  ## the joiner pairs are coded with the dispute-level fatality / hostility like every other pair (war rule, COW codes)
  expect_true(all(ev$in_war[ev$dispnum == 2L])); expect_true(all(ev$in_fmid[ev$dispnum == 2L]))
})

test_that("fatal / war samples: negative fatality codes -> NA, fatal if either side >= 1; war = hostlev 5", {
  ev <- code_events(toy_dd, toy_disp)
  d1 <- ev[ev$dispnum == 1L, ]; d2 <- ev[ev$dispnum == 2L & ev$initiator == 10L, ]; d2b <- ev[ev$dispnum == 2L & ev$initiator == 11L, ]
  expect_false(d1$in_fmid); expect_false(d1$in_war)
  expect_true(is.na(d2$fatal2)); expect_true(d2$in_fmid); expect_true(d2$in_war)
  expect_true(d2b$in_fmid)       # side 1 fatality 0, side 2 NA, dispute-level fatality 1 -> fatal
  ## all-NA fatality -> NA, not FALSE
  dd3 <- toy_dd |> filter(dispnum == 1L) |> mutate(fatality1 = -9L, fatality2 = -9L)
  ev3 <- code_events(dd3, tibble(dispnum = 1L, fatality_disp = -9L, hostlev_disp = 4L))
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
  disp3 <- tibble(dispnum = 1L, fatality_disp = -9L, hostlev_disp = 4L)
  expect_true(code_events(dd3, disp3, mic = tibble(dispnum = 1L, mic_fatalmin = 30L))$in_fmid)
  expect_false(code_events(dd3, disp3, mic = tibble(dispnum = 1L, mic_fatalmin = 0L))$in_fmid)
  expect_false(code_events(dd3, disp3, mic = tibble(dispnum = 1L, mic_fatalmin = 30L))$in_fmidcow)
  ## all-NA COW fatality, a war, not in MIC: fatal by the war rule; COW-only coding drops it
  ev_w <- code_events(dd3, tibble(dispnum = 1L, fatality_disp = -9L, hostlev_disp = 5L))
  expect_true(ev_w$in_war); expect_true(ev_w$in_fmid); expect_false(ev_w$in_fmidcow)
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
  ## 2026-10-24: dyadic peace years left the main Z, the lagged DV dy_mid_tm1 entered; --dy_peace_years restores the former
  expect_equal(z_columns(), c("contig_land", "alliance_def_off", "joint_dem", "major_count", "log_dy_prior_mids", "dy_mid_tm1",
                              "dM_j_h", "dRpr_j_h"))
  expect_false("log_dy_peace_years" %in% z_columns())
  expect_equal(z_columns(DY_PEACE_YEARS = TRUE), c("contig_land", "alliance_def_off", "joint_dem", "major_count", "log_dy_prior_mids", "dy_mid_tm1",
                                                   "log_dy_peace_years", "dM_j_h", "dRpr_j_h"))
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

## ---- power changes: one table per lag, year-centred, one pooled SD (build_change_table in 00_helpers.R) ----------
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
toy_disp2 <- tibble(dispnum = 70:74, fatality_disp = 0L, hostlev_disp = 3L)

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
  ## stage 2: choice sets at t drop already-fighting alternatives; the escalation's own row goes with them, no other event's does
  states <- tibble(ccode = c(10L, 11L, 20L, 21L, 30L, 40L, 50L, 60L))
  alts <- tidyr::expand_grid(ccode1 = c(10L, 40L), year = 1979:1982) |>
    cross_join(states |> rename(ccode2 = ccode)) |> filter(ccode1 != ccode2)
  af <- already_fighting_keys(toy_dd2) |> distinct(ccode1, ccode2, year) |> mutate(af = 1L)
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
