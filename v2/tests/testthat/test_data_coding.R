## Unit tests of the data-coding logic in 01_build_data.R using small hand-built toy frames
## that mimic the peacesciencer columns. The full build needs the bundled peacesciencer data
## and is run by Rscript v2/R/01_build_data.R (on Roar or locally); these tests exercise the
## rules in isolation so a regression in the coding logic fails loudly without the full build.
suppressPackageStartupMessages({ library(testthat); library(dplyr); library(tidyr) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "00_helpers.R"))

## ---- the coding rules under test, lifted verbatim from 01_build_data.R ---------------------
code_events <- function(dd, disp, include_joiners = FALSE) {
  pair <- dd |> group_by(dispnum, ccode1, ccode2) |>
    summarise(year = min(year), sidea1 = first(sidea1), sidea2 = first(sidea2),
              orig1 = first(orig1), orig2 = first(orig2),
              fatality1 = max(fatality1), fatality2 = max(fatality2),
              hostlev1 = max(hostlev1), hostlev2 = max(hostlev2), .groups = "drop") |>
    left_join(disp, by = "dispnum")
  events <- pair |> filter(sidea1 == 1L, sidea2 == 0L)
  if (!include_joiners) events <- events |> filter(orig1 == 1L, orig2 == 1L)
  events |> rename(initiator = ccode1, target = ccode2) |>
    mutate(fatal1 = ifelse(fatality1 < 0, NA, fatality1 >= 1),
           fatal2 = ifelse(fatality2 < 0, NA, fatality2 >= 1),
           fatal_disp = ifelse(fatality_disp < 0, NA, fatality_disp >= 1),
           fatal = case_when(fatal1 %in% TRUE | fatal2 %in% TRUE | fatal_disp %in% TRUE ~ TRUE,
                             is.na(fatal1) & is.na(fatal2) & is.na(fatal_disp) ~ NA,
                             TRUE ~ FALSE),
           war = hostlev_disp == 5L,
           in_mid = TRUE, in_fmid = fatal %in% TRUE, in_war = war %in% TRUE)
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

test_that("joiner flag adds joiner pairs for robustness", {
  ev <- code_events(toy_dd, toy_disp, include_joiners = TRUE)
  expect_equal(nrow(ev), 5)
  expect_true(21L %in% ev$target)
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
  expect_true(is.na(ev3$fatal)); expect_false(ev3$in_fmid)
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
test_that("Z columns: main spec has land contiguity and no distance; flags add distance / drop contiguity", {
  expect_equal(z_columns(), c("contig_land", "alliance_def_off", "joint_dem", "major_count", "log_dy_prior_mids", "log_dy_peace_years"))
  expect_true("log_capdist" %in% z_columns(DISTANCE = TRUE)); expect_true("contig_land" %in% z_columns(DISTANCE = TRUE))
  expect_false("contig_land" %in% z_columns(CONTIGUITY = FALSE))
  expect_false("log_dy_prior_mids" %in% z_columns(CONFLICT_HISTORY = FALSE))
})
