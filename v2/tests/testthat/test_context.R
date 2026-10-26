## Context split (dual typology stage a, 2026-10-25): the coalition flag of the build, the context / coal_ties / coal_drop switches, the
## assembly filter, and the comparison functions. Toy dispute-years and the toy cell of helper-toy-cell.R; the real build when present.
suppressPackageStartupMessages({ library(testthat); library(dplyr); library(tidyr) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "spec.R")); source(file.path(root, "R", "helpers_data.R")); source(file.path(root, "R", "helpers_model.R"))
source(file.path(root, "R", "assemble_cell.R")); source(file.path(root, "R", "typology.R")); source(file.path(root, "R", "context_compare.R"))

## ---- toy dispute-years: directed rows both ways, one dispute ---------------------------------------------------------------------------
dd_rows <- function(dispnum, a, b, years, orig_a = 1L, orig_b = 1L) bind_rows(
  tibble(dispnum = dispnum, ccode1 = a, ccode2 = b, year = years, sidea1 = 1L, sidea2 = 0L, orig1 = orig_a, orig2 = orig_b),
  tibble(dispnum = dispnum, ccode1 = b, ccode2 = a, year = years, sidea1 = 0L, sidea2 = 1L, orig1 = orig_b, orig2 = orig_a))
disp_row <- function(dispnum, styear, stmon, dur) tibble(dispnum = dispnum, styear = styear, stmon = stmon, mindur = dur, maxdur = dur)
part_row <- function(dispnum, ccode, sy, sm, sd, ey, em, ed) tibble(dispnum = dispnum, ccode = ccode, styear = sy, stmon = sm, stday = sd, endyear = ey, endmon = em, endday = ed)

test_that("date keys and comparisons: coarser resolution wins, ties at that resolution, the earliest-by rule degrades to the coarsest tied date", {
  k <- date_keys(c(1L, 2L, 3L), c(1939L, 1939L, 1939L), c(3L, 3L, NA), c(25L, NA, NA))
  expect_equal(k$mkey, c(ym_index(1939L, 3L), ym_index(1939L, 3L), NA)); expect_equal(k$dkey, c(as.integer(as.Date("1939-03-25")), NA, NA))
  d25 <- date_keys(1L, 1939L, 3L, 25L); d31 <- date_keys(1L, 1939L, 3L, 31L); mar <- date_keys(2L, 1939L, 3L, NA); yr <- date_keys(3L, 1939L, NA, NA); apr <- date_keys(2L, 1939L, 4L, NA)
  expect_equal(date_cmp(d25, d31)$cmp, -1L); expect_equal(date_cmp(d25, d31)$res, 1L)
  expect_equal(date_cmp(d25, mar)$cmp, 0L); expect_equal(date_cmp(d25, mar)$res, 2L)        # a day against a month in the same month: tie
  expect_equal(date_cmp(d25, apr)$cmp, -1L); expect_equal(date_cmp(yr, d25)$cmp, 0L); expect_equal(date_cmp(yr, d25)$res, 3L)
  expect_error(date_keys(1L, 1939L, 2L, 30L), "invalid calendar day")
  l <- date_later(d25, d31); expect_equal(l$res, 1L); expect_equal(l$dkey, d31$dkey)
  l <- date_later(d25, mar); expect_equal(l$res, 2L); expect_equal(l$mkey, mar$mkey); expect_true(is.na(l$dkey))     # tie -> the coarser date
  e <- date_earliest_by(bind_rows(d31, d25, apr), c(1L, 1L, 1L)); expect_equal(e$res, 1L); expect_equal(e$dkey, d25$dkey)
  e <- date_earliest_by(bind_rows(d31, mar, apr), c(1L, 1L, 1L)); expect_equal(e$res, 2L); expect_equal(e$mkey, mar$mkey)
  e <- date_earliest_by(bind_rows(d25, yr), c(1L, 1L)); expect_equal(e$res, 3L)
  e <- date_earliest_by(bind_rows(d25, apr, d31), c(2L, 1L, 2L)); expect_equal(e$g, c(1L, 2L)); expect_equal(e$dkey, c(apr$dkey * NA, d25$dkey))
})

test_that("participant spell dates: GML day / month where the years agree, clamping, fallbacks, re-entries as separate spells, membership disagreements", {
  dd <- bind_rows(dd_rows(1L, 10L, 20L, 1939:1940), dd_rows(1L, 10L, 30L, 1939:1940, orig_b = 0L), dd_rows(1L, 10L, 40L, c(1939L, 1941L), orig_b = 0L),
                  dd_rows(2L, 50L, 60L, 1950L), dd_rows(3L, 70L, 80L, 1960L))
  disp <- bind_rows(disp_row(1L, 1939L, 3L, 800), disp_row(2L, 1950L, 5L, 0), disp_row(3L, 1960L, 7L, 10))   # dispute 2 ends within its start month (exact end month)
  db <- dispute_date_bounds(disp, dd)
  part <- bind_rows(part_row(1L, 10L, 1939L, 3L, 25L, 1941L, 5L, 7L), part_row(1L, 20L, 1939L, 3L, 25L, 1940L, 6L, -9L),      # 20: exit month only
                    part_row(1L, 30L, 1938L, 12L, 1L, 1940L, 12L, 31L),                                                     # 30: GML start before the COW year -> clamped January
                    part_row(1L, 40L, 1939L, 9L, 3L, 1939L, 9L, 20L), part_row(1L, 40L, 1941L, 6L, -9L, 1941L, 8L, 1L),     # 40: two COW spells, one GML spell each
                    part_row(1L, 99L, 1939L, 4L, 1L, 1939L, 4L, 2L),                                                        # 99 in GML, not in COW
                    part_row(2L, 50L, 1950L, 5L, 2L, 1950L, 5L, 12L))                                                        # 60 in COW, not in GML; dispute 3 absent from GML
  psd <- participant_spell_dates(dd, db, part)
  g <- function(c, sp = 1L) psd[psd$dispnum == 1L & psd$ccode == c & psd$spell == sp, ]
  expect_equal(g(10L)$entry_res, 1L); expect_equal(g(10L)$entry_dkey, as.integer(as.Date("1939-03-25"))); expect_equal(g(10L)$exit_res, 1L)
  expect_equal(g(20L)$exit_res, 2L); expect_equal(g(20L)$exit_mkey, ym_index(1940L, 6L)); expect_equal(g(20L)$exit_src, "participant")
  expect_equal(g(30L)$entry_res, 2L); expect_equal(g(30L)$entry_mkey, ym_index(1939L, 1L)); expect_match(g(30L)$entry_src, "clamped")
  expect_equal(nrow(psd[psd$dispnum == 1L & psd$ccode == 40L, ]), 2L); expect_equal(g(40L, 1L)$entry_dkey, as.integer(as.Date("1939-09-03"))); expect_equal(g(40L, 2L)$entry_res, 2L)
  ## fallbacks: dispute 2 participant 60 (no GML row, originator in the first year) -> dispute start month; dispute 3 (absent from GML) the same
  p60 <- psd[psd$dispnum == 2L & psd$ccode == 60L, ]; expect_equal(p60$entry_res, 2L); expect_equal(p60$entry_mkey, ym_index(1950L, 5L)); expect_match(p60$entry_src, "no overlapping")
  expect_equal(p60$exit_res, 2L); expect_equal(p60$exit_mkey, ym_index(1950L, 5L))   # the dispute's end bounds are one exact month
  p70 <- psd[psd$dispnum == 3L & psd$ccode == 70L, ]; expect_equal(p70$entry_res, 2L); expect_match(p70$entry_src, "not in gml_part")
  expect_equal(p70$exit_res, 3L)                                      # 10 days from an unknown July day: the end spans July / August -> year
  mem <- attr(psd, "membership")
  expect_equal(nrow(mem), 2L); expect_true(any(mem$ccode == 99L & mem$issue == "gml_part participant absent from COW")); expect_true(any(mem$ccode == 60L & mem$issue == "COW participant absent from gml_part"))
  expect_true(all(psd$entry_y == psd$spell_first)); expect_true(all(psd$exit_y == psd$spell_last))
  ## without GML: every entry is a dispute-level date
  p0 <- participant_spell_dates(dd, db, part = NULL); expect_true(all(!startsWith(p0$entry_src, "participant"))); expect_equal(nrow(attr(p0, "membership")), 0L)
})

## a Germany-1939-like dispute: 1 = initiator 10 (orig) vs 20 (orig, same day), 30 and 40 join on day 31 of the same month, 50 joins in
## September; the initiator's 1939 group has pairs with 20 (March 25), 30 and 40 (March 31) and 50 (September); 60 attacks 20 in 1940
toy_context <- function(part_days = TRUE) {
  dd <- bind_rows(dd_rows(1L, 10L, 20L, 1939:1940), dd_rows(1L, 10L, 30L, 1939:1940, orig_b = 0L), dd_rows(1L, 10L, 40L, 1939:1940, orig_b = 0L),
                  dd_rows(1L, 10L, 50L, 1939:1940, orig_b = 0L), dd_rows(1L, 60L, 20L, 1940L, orig_a = 0L), dd_rows(1L, 60L, 30L, 1940L, orig_a = 0L, orig_b = 0L),
                  dd_rows(2L, 70L, 80L, 1950L), dd_rows(3L, 70L, 80L, 1950L))
  disp <- bind_rows(disp_row(1L, 1939L, 3L, 600), disp_row(2L, 1950L, 2L, 10), disp_row(3L, 1950L, 8L, 10))
  sd <- function(d) if (part_days) d else -9L
  part <- bind_rows(part_row(1L, 10L, 1939L, 3L, sd(25L), 1940L, 12L, 1L), part_row(1L, 20L, 1939L, 3L, sd(25L), 1940L, 12L, 1L),
                    part_row(1L, 30L, 1939L, 3L, sd(31L), 1940L, 12L, 1L), part_row(1L, 40L, 1939L, 3L, sd(31L), 1940L, 12L, 1L),
                    part_row(1L, 50L, 1939L, 9L, sd(3L), 1940L, 12L, 1L), part_row(1L, 60L, 1940L, 5L, sd(10L), 1940L, 12L, 1L))
  db <- dispute_date_bounds(disp, dd)
  pdi <- pair_dispute_intervals(dd, db, part = part)
  pairs <- bind_rows(tibble(dispnum = 1L, initiator = 10L, target = c(20L, 30L, 40L, 50L), year = 1939L), tibble(dispnum = 1L, initiator = 60L, target = c(20L, 30L), year = 1940L),
                     tibble(dispnum = c(2L, 3L), initiator = 70L, target = 80L, year = 1950L)) |>
    inner_join(pdi |> filter(spell == 1L) |> select(dispnum, initiator = ccode1, target = ccode2, ev_start_lo = start_lo, ev_start_hi = start_hi), by = c("dispnum", "initiator", "target"))
  list(dd = dd, db = db, part = part, pairs = pairs, psd = participant_spell_dates(dd, db, part))
}

test_that("coalition flag: evaluated at the GROUP onset, constant within (initiator, dispute, year), strict <= inclusive, both tie conventions at day and month resolution", {
  tc <- toy_context(part_days = TRUE)
  cc <- coalition_context(tc$pairs, tc$psd); grp <- attr(cc, "groups"); evt <- attr(cc, "events")
  expect_equal(nrow(cc), nrow(tc$pairs)); expect_false(anyNA(cc$coal_incl)); expect_false(anyNA(cc$coal_strict))
  ## the 1939 group of 10: onset = March 25 (the pair with 20); 30 / 40 / 50 enter later -> bilateral under both conventions for all four pairs
  g10 <- cc[cc$initiator == 10L, ]
  expect_true(all(!g10$coal_incl)); expect_true(all(!g10$coal_strict)); expect_true(all(g10$coal_res == "day"))
  expect_equal(grp$onset_dkey[grp$initiator == 10L], as.integer(as.Date("1939-03-25"))); expect_equal(grp$n_incl[grp$initiator == 10L], 1L)
  ## at each pair's OWN onset the flags would differ inside the group (the diagnostic): 10 -> 30 on March 31 has 20 present
  expect_equal(grp$n_event_flags_incl[grp$initiator == 10L], 2L); expect_true(evt$ev_incl[evt$initiator == 10L & evt$target == 30L]); expect_false(evt$ev_incl[evt$initiator == 10L & evt$target == 20L])
  ## 60's 1940 group: a joiner against 20 and 30 with 10, 40, 50 present since 1939 -> coalition under both conventions
  g60 <- cc[cc$initiator == 60L, ]; expect_true(all(g60$coal_incl)); expect_true(all(g60$coal_strict))
  ## the bilateral disputes 2 and 3 (70 vs 80) are bilateral; the strict flag never exceeds the inclusive one; constant within every group
  expect_true(all(!cc$coal_incl[cc$dispnum %in% 2:3]))
  expect_true(all(cc$coal_incl >= cc$coal_strict))
  chk <- cc |> group_by(dispnum, initiator, year) |> summarise(a = n_distinct(coal_incl), b = n_distinct(coal_strict), .groups = "drop"); expect_true(all(chk$a == 1L & chk$b == 1L))
  ## month resolution (no days in GML): 10's group onset is March 1939; 30 and 40 enter the SAME month -> inclusive coalition, strict bilateral
  tm <- toy_context(part_days = FALSE)
  cm <- coalition_context(tm$pairs, tm$psd); gm <- attr(cm, "groups")
  expect_true(all(cm$coal_incl[cm$initiator == 10L])); expect_true(all(!cm$coal_strict[cm$initiator == 10L])); expect_equal(gm$coal_res[gm$initiator == 10L], "month")
  expect_equal(gm$n_incl[gm$initiator == 10L], 3L); expect_equal(gm$n_strict[gm$initiator == 10L], 0L)
  ## the strict flag with a target that entered BEFORE the initiator: 60 (May 1940) against 20: 20 is the pair, the others make the coalition
  expect_true(gm$target_strict[gm$initiator == 60L]); expect_true(all(cm$coal_strict[cm$initiator == 60L]))
  ## a mismatched pair (no participant spell containing the year) is an error
  expect_error(coalition_context(tc$pairs |> mutate(year = ifelse(initiator == 70L, 1951L, year)), tc$psd), "no participant spell")
})

test_that("ambiguous onset: both states in the dispute before they were opposed (side switch) falls back to the COW year and is flagged", {
  ## 10 and 20 are in dispute 1 from 1939; they are opposed only from 1941 (20 switches sides): the pair's COW onset is 1941
  dd <- bind_rows(dd_rows(1L, 10L, 30L, 1939:1941), dd_rows(1L, 20L, 30L, 1939:1940), dd_rows(1L, 10L, 20L, 1941L, orig_b = 0L))
  dd$orig1[dd$ccode1 == 20L] <- 0L; dd$orig2[dd$ccode2 == 20L] <- 0L
  disp <- disp_row(1L, 1939L, 3L, 1000); db <- dispute_date_bounds(disp, dd)
  part <- bind_rows(part_row(1L, 10L, 1939L, 3L, 1L, 1941L, 12L, 1L), part_row(1L, 30L, 1939L, 3L, 1L, 1941L, 12L, 1L), part_row(1L, 20L, 1939L, 6L, 1L, 1940L, 12L, 1L), part_row(1L, 20L, 1941L, 1L, 1L, 1941L, 12L, 1L))
  psd <- participant_spell_dates(dd, db, part)
  pdi <- pair_dispute_intervals(dd, db, part = part)
  pairs <- tibble(dispnum = 1L, initiator = 10L, target = c(30L, 20L), year = c(1939L, 1941L)) |>
    inner_join(pdi |> filter(spell == 1L) |> select(dispnum, initiator = ccode1, target = ccode2, ev_start_lo = start_lo, ev_start_hi = start_hi), by = c("dispnum", "initiator", "target"))
  cc <- coalition_context(pairs, psd)
  expect_true(cc$coal_onset_ambiguous[cc$target == 20L]); expect_false(cc$coal_onset_ambiguous[cc$target == 30L])
  expect_equal(attr(cc, "groups")$onset_res[attr(cc, "groups")$year == 1941L], 3L)       # year resolution for the ambiguous group
  expect_true(cc$coal_incl[cc$target == 20L])                                             # 30 is present in 1941 under both conventions
  expect_true(cc$coal_strict[cc$target == 20L])
})

test_that("spec: the context switches, their suffixes (order, uniqueness, _strict only with a context, _drop3 only with coalition), the preset", {
  main <- parse_spec(character(0))
  expect_equal(main$context, "all"); expect_equal(main$coal_ties, "inclusive"); expect_equal(main$coal_drop, "none"); expect_true(spec_is_main(main))
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec("--context=bilateral")), "fmid_lag3_all_pr_M1b_bi")
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec("--context=coalition")), "fmid_lag3_all_pr_M1b_co")
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec(c("--context=coalition", "--coal_ties=strict"))), "fmid_lag3_all_pr_M1b_co_strict")
  expect_equal(cell_name("fmid", 3L, "M2b", parse_spec(c("--context=coalition", "--coal_drop=top3"))), "fmid_lag3_all_pr_M2b_K2_co_drop3")
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec(c("--context=coalition", "--coal_ties=strict", "--coal_drop=top3", "--ww2=drop", "--direction=sideA"))), "fmid_lag3_all_pr_M1b_co_strict_drop3_noww2_sideA")
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec(c("--pairs=joiners", "--context=bilateral"))), "fmid_lag3_all_pr_M1b_join_bi")
  ## coal_ties is inert under context = all: normalised, no suffix, still the main specification
  s <- parse_spec("--coal_ties=strict"); expect_equal(s$coal_ties, "inclusive"); expect_true(spec_is_main(s)); expect_equal(cell_name("mid", 3L, "M1b", s), "mid_lag3_all_pr_M1b")
  expect_error(parse_spec("--coal_drop=top3"), "needs --context=coalition"); expect_error(parse_spec(c("--coal_drop=top3", "--context=bilateral")), "needs --context=coalition")
  expect_error(parse_spec("--context=multilateral"), "--context must be one of"); expect_error(parse_spec("--coal_ties=loose"), "--coal_ties must be one of")
  ## every suffix is distinct, including the fixed ones (_K<K>, _pilot, _STALE) and the archived _samp / _dist / _noCH / _polrel
  all_suffixes <- c(unlist(SPEC_SUFFIX), "_K2", "_pilot", "_STALE", "_samp10", "_dist", "_noCH", "_polrel")
  expect_false(any(duplicated(all_suffixes)))
  expect_equal(unname(SPEC_SUFFIX$context), c("_bi", "_co")); expect_equal(unname(SPEC_SUFFIX$coal_ties), "_strict"); expect_equal(unname(SPEC_SUFFIX$coal_drop), "_drop3")
  expect_equal(spec_coal_col(parse_spec("--context=bilateral")), "coal_incl"); expect_equal(spec_coal_col(parse_spec(c("--context=bilateral", "--coal_ties=strict"))), "coal_strict")
  ## the preset: M1b, lag 3, mid and fmid, no war; three variants
  expect_true("context" %in% names(ROBUSTNESS_PRESETS)); expect_equal(CONTEXT_PRESET$cells, c("mid:3:M1b", "fmid:3:M1b"))
  expect_equal(names(CONTEXT_PRESET$variants), c("bilateral", "coalition", "coalition_drop3"))
  for (v in CONTEXT_PRESET$variants) for (nm in names(v)) expect_true(v[[nm]] %in% SPEC_CHOICES[[nm]])
  expect_equal(CONTEXT_PRESET$variants$coalition_drop3$coal_drop, "top3")
  expect_false(any(grepl("war", CONTEXT_PRESET$cells)))
})

test_that("assembly: context = all is the plain call; bilateral and coalition partition the stage-2 events under each tie convention; stage 1 and the weights are untouched", {
  td <- write_toy_cell(tempfile("toycell"))
  A  <- quiet(assemble_cell("mid", 3L, "all", data_dir = td))
  A0 <- quiet(assemble_cell("mid", 3L, "all", CONTEXT = "all", COAL_TIES = "strict", data_dir = td))
  d <- function(X) build_stan_data(X$sy, X$alts, X$ev, X$W_cols, X$X_cols, X$Z_cols, joint = 1L, x_ref = X$x_ref, allow_orphan_initiations = TRUE)
  expect_identical(d(A), d(A0)); expect_identical(A$ev, A0$ev); expect_identical(A$sy, A0$sy)
  key <- function(X) paste(X$alts$ccode1[X$ev$alt_row], X$alts$ccode2[X$ev$alt_row], X$alts$year[X$ev$alt_row])
  for (ties in c("inclusive", "strict")) {
    B <- quiet(assemble_cell("mid", 3L, "all", CONTEXT = "bilateral", COAL_TIES = ties, data_dir = td))
    C <- quiet(assemble_cell("mid", 3L, "all", CONTEXT = "coalition", COAL_TIES = ties, data_dir = td))
    expect_equal(length(intersect(key(B), key(C))), 0L); expect_setequal(c(key(B), key(C)), key(A))        # disjoint, union = all usable events
    expect_equal(nrow(B$ev) + nrow(C$ev), nrow(A$ev))
    ## stage 1 unchanged: same sy rows, same y1, same W (no recount), orphans appear
    expect_identical(B$sy, A$sy); expect_identical(C$sy, A$sy); expect_identical(B$alts, A$alts)
    expect_true(B$n_orphan > 0 || C$n_orphan > 0)
    ## the kept events carry their build weights
    ab <- A$ev |> inner_join(B$ev, by = c("sy_row", "alt_row"), suffix = c("", "_b")); expect_equal(ab$w, ab$w_b); expect_equal(nrow(ab), nrow(B$ev))
    expect_equal(B$context, "bilateral"); expect_equal(C$coal_ties, ties); expect_equal(B$n_context_dropped, nrow(C$ev)); expect_null(B$context_top3)
  }
  ## inclusive: 10 -> 20 and 30 -> 40 (stage 1 only anyway) bilateral -> B has 1 usable event; strict adds 30 -> 20
  Bi <- quiet(assemble_cell("mid", 3L, "all", CONTEXT = "bilateral", data_dir = td)); expect_equal(nrow(Bi$ev), 1L)
  Bs <- quiet(assemble_cell("mid", 3L, "all", CONTEXT = "bilateral", COAL_TIES = "strict", data_dir = td)); expect_equal(nrow(Bs$ev), 2L)
  ## the Stan data of a context cell build and validate (orphan initiations allowed, as 02_fit_model.R does)
  Ci <- quiet(assemble_cell("mid", 3L, "all", CONTEXT = "coalition", data_dir = td)); dc <- d(Ci)
  expect_equal(dc$E, 4L); expect_equal(dc$N1, 4L); expect_equal(sum(dc$y1), 3L)
  ## drop3: the 3 disputes with the most coalition events (ties: smallest dispute numbers), only with the coalition context
  expect_error(quiet(assemble_cell("mid", 3L, "all", CONTEXT = "bilateral", COAL_DROP = "top3", data_dir = td)), "needs CONTEXT = coalition")
  expect_error(quiet(assemble_cell("mid", 3L, "all", COAL_DROP = "top3", data_dir = td)), "needs CONTEXT = coalition")
  D3 <- quiet(assemble_cell("mid", 3L, "all", CONTEXT = "coalition", COAL_DROP = "top3", data_dir = td))
  ## coalition events: 10 -> 30 (4;5), 30 -> 20 (2), 30 -> 10 (6), 40 -> 20 (7): one event per dispute, top3 = 2, 4, 5 -> 10 -> 30 and 30 -> 20 leave
  expect_equal(D3$context_top3$dispnum, c(2L, 4L, 5L)); expect_equal(nrow(D3$ev), 2L); expect_equal(D3$n_context_dropped, 2L + 2L)
  expect_setequal(key(D3), c("30 10 1990", "40 20 1990"))
  ## a key whose coalition flag is NA on a non-escalation row, or strict > inclusive, is rejected
  td2 <- write_toy_cell(tempfile("toycell"), events_edit = function(e) { e$coal_incl[1] <- NA; e })
  expect_error(quiet(assemble_cell("mid", 3L, "all", data_dir = td2)))
  td3 <- write_toy_cell(tempfile("toycell"), events_edit = function(e) { e$coal_strict[1] <- TRUE; e })
  expect_error(quiet(assemble_cell("mid", 3L, "all", data_dir = td3)))
  ## a group split among single-dispute members (flag not constant within the group) is rejected
  td4 <- write_toy_cell(tempfile("toycell"), events_edit = function(e) { e$dispnum[3] <- 1L; e$dispnums[3] <- "1"; e$initiator[3] <- 10L; e$target[3] <- 40L; e })   # 10 -> 40 (dispute 1) coalition, 10 -> 20 (dispute 1) bilateral
  expect_error(quiet(assemble_cell("mid", 3L, "all", CONTEXT = "bilateral", data_dir = td4)), "not constant within the group")
  ## assemble_cell_spec() passes the switches through
  S1 <- quiet(assemble_cell_spec("mid", 3L, parse_spec(c("--context=coalition", "--coal_drop=top3")), data_dir = td)); expect_identical(S1$ev, D3$ev)
  S2 <- quiet(assemble_cell_spec("mid", 3L, parse_spec(c("--context=bilateral", "--coal_ties=strict")), data_dir = td)); expect_identical(S2$ev, Bs$ev)
})

test_that("robustness driver: the context preset is a fixed six-row table; 02_fit_model.R records the switches", {
  src <- readLines(file.path(root, "R", "05_robustness.R"))
  expect_true(any(grepl('PRESET == "context"', src, fixed = TRUE))); expect_true(any(grepl("CONTEXT_PRESET$variants", src, fixed = TRUE)))
  sb <- readLines(file.path(root, "slurm", "04_robustness.sbatch")); expect_true(any(grepl("^  context)", sb)))
  fit_src <- readLines(file.path(root, "R", "02_fit_model.R")); expect_true(any(grepl("context = spec$context, coal_ties = spec$coal_ties, coal_drop = spec$coal_drop", fit_src, fixed = TRUE)))
})

## ---- the comparison on tiny mock fits (quickly fitted; needs the compiled M1 program) --------------------------------------------------
test_that("context_compare() runs on two tiny fits (plus a drop3 fit): blocks, intervals, partition of the role split, the markdown summary", {
  sm <- skip_if_no_build("m1_latent_additive")
  set.seed(123)
  sim <- simulate_targeting("M1b", S = 6L, Tt = 16L, K = 2L, seed = 123L, target_rate = 0.6)
  ev <- sim$ev; ev$w <- rep(c(1, 0.5), length.out = nrow(ev))
  sy_keys <- data.frame(ccode = 100L + sim$sy$state_idx, year = sim$sy$year, state_idx = sim$sy$state_idx, y1 = sim$sy$y1)
  state_codes <- 100L + seq_len(6L)
  mock_fit <- function(rows, cell) {
    dat <- build_stan_data(sim$sy, sim$alts, ev[rows, ], c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"), joint = 1L, x_ref = rep(0, 4), allow_orphan_initiations = TRUE)
    rstan_options(threads_per_chain = 1)
    fit <- sampling(sm, data = dat, chains = 1, iter = 60, warmup = 30, seed = 123, refresh = 0, pars = c("z_re", "z_b", "L_Om_joint", "L_Om", "re"), include = FALSE)
    list(fit = fit, dat = dat, cell = cell, state_codes = state_codes, sy_keys = sy_keys, args = list(pilot = FALSE, chains = 1L, iter = 60L, warmup = 30L))
  }
  odd <- seq(1, nrow(ev), by = 2); even <- seq(2, nrow(ev), by = 2)
  fits <- list(bilateral = mock_fit(odd, "toy_bi"), coalition = mock_fit(even, "toy_co"), coalition_drop3 = mock_fit(even[-(1:2)], "toy_co_drop3"))
  jl <- data.frame(initiator = state_codes[sim$sy$state_idx[ev$sy_row]], target = state_codes[sim$alts$target_idx[ev$alt_row]], year = sim$sy$year[ev$sy_row],
                   is_joiner_pair = rep(c(FALSE, FALSE, TRUE), length.out = nrow(ev)))
  cc <- context_compare(fits, jl, r0 = 0.5, min_decisions = 2)
  st <- cc$stats
  expect_setequal(unique(st$block), c("context", "leverage", "role"))
  ctx <- st[st$block == "context", ]
  expect_equal(nrow(ctx), 7L)                                                                    # 3 statistics x 2 axes + type agreement
  for (ax in c("material", "relational")) {
    r <- ctx[ctx$axis == ax, ]; expect_setequal(r$statistic, c("cor_pmean", "cor_draw", "mean_shift"))
    cd <- r[r$statistic == "cor_draw", ]; expect_true(cd$lo <= cd$estimate && cd$estimate <= cd$hi); expect_true(abs(cd$estimate) <= 1)
    ms <- r[r$statistic == "mean_shift", ]; expect_true(ms$lo <= ms$estimate && ms$estimate <= ms$hi)
    expect_true(abs(r$estimate[r$statistic == "cor_pmean"]) <= 1); expect_true(is.na(r$lo[r$statistic == "cor_pmean"]))
  }
  ta <- ctx[ctx$statistic == "type_agreement", ]; expect_true(ta$estimate >= 0 && ta$estimate <= 1); expect_equal(unique(ctx$n_states), ta$n_states)
  ct <- cc$crosstab[cc$crosstab$block == "context", ]; expect_equal(sum(ct$n), ta$n_states); expect_equal(nrow(ct), 25L)
  ## the restriction: states with >= 2 weighted decisions in BOTH fits
  wb <- weighted_decisions(fits$bilateral$dat); wc <- weighted_decisions(fits$coalition$dat)
  expect_equal(ta$n_states, sum(wb >= 2 & wc >= 2)); expect_equal(sum(cc$states$in_comparison), ta$n_states)
  expect_equal(sum(wb), sum(ev$w[odd]))
  lev <- st[st$block == "leverage", ]; expect_equal(nrow(lev), 7L); expect_true(all(lev$comparison == "coalition_drop3_vs_coalition"))
  ## role block: raw weighted means per context x role, and the per-state coalition split; marked as data descriptives
  rl <- st[st$block == "role", ]
  expect_setequal(unique(rl$comparison), c("bilateral_originator", "bilateral_joiner", "coalition_originator", "coalition_joiner", "coalition_joiner_vs_originator_by_state"))
  expect_true(all(grepl("DATA DESCRIPTIVE", rl$note))); expect_true(all(is.na(rl$lo)))
  wm <- rl[rl$statistic == "weighted_mean_chosen_gap" & rl$comparison == "coalition_originator" & rl$axis == "material", ]
  fe <- fit_events(fits$coalition) |> inner_join(jl |> distinct(initiator, target, year, is_joiner_pair), by = c("initiator", "target", "year")) |> filter(!is_joiner_pair)
  expect_equal(wm$estimate, sum(fe$w * fe$x_material) / sum(fe$w))
  ## the direct per-draw correlation reproduces the table
  mu_b <- mu_state_draws(fits$bilateral$fit); mu_c <- mu_state_draws(fits$coalition$fit); R <- which(wb >= 2 & wc >= 2)
  cors <- sapply(seq_len(dim(mu_b)[1]), function(k) cor(mu_b[k, R, 1], mu_c[k, R, 1]))
  expect_equal(ctx$estimate[ctx$axis == "material" & ctx$statistic == "cor_draw"], mean(cors))
  expect_equal(ctx$estimate[ctx$axis == "relational" & ctx$statistic == "mean_shift"], mean(mu_c[, R, 2] - mu_b[, R, 2]))
  ## without the drop3 fit the leverage block is absent; the markdown summary renders in both cases
  cc2 <- context_compare(fits[c("bilateral", "coalition")], jl, r0 = 0.5, min_decisions = 2); expect_false(any(cc2$stats$block == "leverage"))
  md <- context_compare_md(cc, "toy"); expect_true(any(grepl("Coalition vs bilateral", md))); expect_true(any(grepl("Leverage: coalition without", md)))
  md2 <- context_compare_md(cc2, "toy"); expect_true(any(grepl("No coalition_drop3 fit", md2)))
  ## mismatched draw counts and a joiner lookup with gaps are errors
  short <- fits$coalition; short$fit <- fits$bilateral$fit
  expect_error(context_compare(list(bilateral = fits$bilateral, coalition = list(fit = sampling(sm, data = fits$coalition$dat, chains = 1, iter = 40, warmup = 20, seed = 1, refresh = 0, pars = c("z_re", "z_b", "L_Om_joint", "L_Om", "re"), include = FALSE),
                                                                                  dat = fits$coalition$dat, state_codes = state_codes, sy_keys = sy_keys)), jl, r0 = 0.5, min_decisions = 2), "draw dimensions")
  expect_error(context_compare(fits, jl[-1, ], r0 = 0.5, min_decisions = 2), "no is_joiner_pair")
})

## ---- the real build (when present) ------------------------------------------------------------------------------------------------------
test_that("real build (when present): coal flags constant within the D8 groups, strict <= inclusive per sample, bilateral + coalition partition the fmid lag-3 events, context = all reproduces the last dry run", {
  f <- file.path(v2_root(), "data", "events.rds"); man_path <- file.path(v2_root(), "data", "manifest.csv")
  skip_if_not(file.exists(f) && file.exists(man_path) && file.exists(file.path(v2_root(), "data", "alts_lag3_all.rds")), "real data files not built")
  ev <- readRDS(f) |> filter(!escalation)
  expect_true(all(c("coal_incl", "coal_strict", "coal_onset_ambiguous") %in% names(ev))); expect_false(anyNA(ev$coal_incl)); expect_false(anyNA(ev$coal_strict))
  expect_true(all(ev$coal_incl >= ev$coal_strict))
  for (s in c("mid", "fmid", "war")) { ff <- ev[[paste0("in_", s)]]; expect_true(sum(ev$coal_strict[ff]) <= sum(ev$coal_incl[ff])); expect_true(sum(ev$coal_incl[ff]) > 0); expect_true(sum(!ev$coal_incl[ff]) > 0) }
  ## constant within (initiator, dispute, year) among single-dispute events (a collapsed key carries the OR of its disputes)
  g <- ev |> filter(n_disputes == 1L) |> group_by(initiator, dispnum, year) |> summarise(a = n_distinct(coal_incl), b = n_distinct(coal_strict), .groups = "drop")
  expect_true(all(g$a == 1L)); expect_true(all(g$b == 1L))
  ## Germany -> Poland 1939 is bilateral under both conventions (group onset 1939-03-25, UK and France enter 1939-03-31)
  gp <- ev |> filter(initiator == 255L, target == 290L, year == 1939L); expect_equal(nrow(gp), 1L); expect_false(gp$coal_incl); expect_false(gp$coal_strict)
  ## the partition on the assembled fmid cell, and context = all against the saved dry run of the same data
  A <- quiet(assemble_cell("fmid", 3L, "all")); B <- quiet(assemble_cell("fmid", 3L, "all", CONTEXT = "bilateral")); C <- quiet(assemble_cell("fmid", 3L, "all", CONTEXT = "coalition"))
  key <- function(X) paste(X$alts$ccode1[X$ev$alt_row], X$alts$ccode2[X$ev$alt_row], X$alts$year[X$ev$alt_row])
  expect_equal(length(intersect(key(B), key(C))), 0L); expect_setequal(c(key(B), key(C)), key(A)); expect_identical(B$sy, A$sy); expect_identical(C$sy, A$sy)
  expect_true(B$n_context_dropped >= nrow(C$ev)); expect_true(C$n_context_dropped >= nrow(B$ev))      # the dropped rows include events without a stage-2 row
  expect_true(B$n_context_split_groups <= nrow(ev |> filter(in_fmid, n_disputes > 1L)))               # a split group needs a collapsed key
  sd_path <- file.path(v2_root(), "output", "standata_fmid_lag3_all_pr_M1b.rds")
  if (file.exists(sd_path)) {
    sd <- readRDS(sd_path); man <- read.csv(man_path, stringsAsFactors = FALSE)
    same_inputs <- all(sd$manifest$sha256[sd$manifest$file != "events.rds"] == man$sha256[match(sd$manifest$file[sd$manifest$file != "events.rds"], man$file)])
    if (same_inputs) {
      dat <- build_stan_data(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, x_ref = A$x_ref, allow_orphan_initiations = TRUE, D_mode = 3L, grainsize = sd$dat$grainsize)
      expect_identical(dat, sd$dat)
    } else skip("the saved dry run was made from other panel / alternatives files")
  }
})
