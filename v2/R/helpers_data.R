## ============================================================================
## helpers_data.R -- data-side helper functions of the v2 pipeline (functions only, no side effects)
## Sourced after spec.R by every script: paths and logging, assertions, the power-change table, the fatal-MID coding,
## dispute participation and the month-level escalation rule, the two choice sets, the duplicate collapse, event weights,
## the initiator-direction rule (D1), the WWII list, the gap scale, and the manifest checks. Model-side helpers (Stan data
## assembly, anchors, sampled alternatives, the brute-force likelihood, fit provenance and safety) are in helpers_model.R.
## ============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

## ---- paths -----------------------------------------------------------------
v2_root <- function() {
  root <- Sys.getenv("V2_ROOT", unset = "")
  if (root == "") {
    ## fall back: directory that contains this file's parent (v2/)
    args <- commandArgs(trailingOnly = FALSE)
    fa <- sub("^--file=", "", args[grepl("^--file=", args)])
    if (length(fa) == 1) root <- normalizePath(file.path(dirname(fa), ".."))
    else root <- normalizePath("v2", mustWork = FALSE)
  }
  root
}

## ---- logging -----------------------------------------------------------------
## .v2_log_file is the one piece of state: the path logf() appends to once log_open() has been called (NULL = console only)
.v2_log_file <- NULL
log_open <- function(name) {
  dir.create(file.path(v2_root(), "logs"), showWarnings = FALSE, recursive = TRUE)
  .v2_log_file <<- file.path(v2_root(), "logs",
                             paste0(name, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
  cat("", file = .v2_log_file)
  logf("log opened: %s", .v2_log_file)
  invisible(.v2_log_file)
}
logf <- function(fmt, ...) {
  msg <- sprintf(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", fmt), ...)
  cat(msg, "\n", sep = "")
  if (!is.null(.v2_log_file)) cat(msg, "\n", sep = "", file = .v2_log_file, append = TRUE)
  invisible(msg)
}

## ---- assertions ------------------------------------------------------------------
assert_no_dup_keys <- function(df, keys, what = "data") {
  n_dup <- nrow(df) - nrow(distinct(df[, keys, drop = FALSE]))
  if (n_dup != 0) stop(sprintf("%s: %d duplicate rows on keys (%s)", what, n_dup, paste(keys, collapse = ", ")))
  invisible(TRUE)
}
assert_no_na <- function(df, cols, what = "data") {
  bad <- cols[vapply(cols, function(cc) anyNA(df[[cc]]), logical(1))]
  if (length(bad) > 0) stop(sprintf("%s: NA present in model columns: %s", what, paste(bad, collapse = ", ")))
  invisible(TRUE)
}
assert_finite <- function(df, cols, what = "data") {
  bad <- cols[vapply(cols, function(cc) !all(is.finite(df[[cc]])), logical(1))]
  if (length(bad) > 0) stop(sprintf("%s: non-finite values (NA / NaN / Inf) in columns: %s", what, paste(bad, collapse = ", ")))
  invisible(TRUE)
}
assert_standardized <- function(df, col, by = "year", tol = 1e-6, what = "data") {
  chk <- df |> group_by(across(all_of(by))) |>
    summarise(m = mean(.data[[col]]), s = sd(.data[[col]]), n = n(), .groups = "drop") |>
    filter(n > 1)
  if (any(abs(chk$m) > tol) || any(abs(chk$s - 1) > 1e-3))
    stop(sprintf("%s: %s is not standardized within %s (max |mean| = %.2e, sd range [%.4f, %.4f])",
                 what, col, by, max(abs(chk$m)), min(chk$s), max(chk$s)))
  invisible(TRUE)
}

## ---- power changes (PI decision 2026-10-23; see README "Positions" and "Changes") -----------
## pos  : state-year levels with columns ccode, year, M (log CINC share), R_pr (log IGO PageRank)
## at_t : the (ccode, year) state-years at t for which a change is wanted (the state-year panel)
## h    : window length. The change at (ccode, t) is level(t-1) - level(t-1-h), by CALENDAR joins.
## Returns list(tbl, scale). tbl has one row per (ccode, year) in at_t with levels at both dates:
##   dM_raw_h, dRpr_raw_h : the raw change in log units (log ratio of the state's own share across the window)
##   dM_h, dRpr_h         : the raw change minus the mean raw change of all states with a change in that year,
##                          divided by ONE pooled constant per lag and variable (the SD of the year-centred change
##                          over all rows of the table). Zero means "changed like the average state that year";
##                          the pooled SD is 1 by construction.
## scale is a data.frame(h, variable, pooled_sd, n_rows, n_years) with the two constants.
## The same table serves the initiator change in stage 1 (W, X) and the target change in stage 2 (Z).
build_change_table <- function(pos, at_t, h) {
  stopifnot(all(c("ccode", "year", "M", "R_pr") %in% names(pos)), all(c("ccode", "year") %in% names(at_t)),
            length(h) == 1, h >= 1)
  h <- as.integer(h)
  assert_no_dup_keys(pos, c("ccode", "year"), "change table: positions")
  assert_finite(pos, c("M", "R_pr"), "change table: positions")
  tbl <- at_t |> distinct(ccode, year) |>
    inner_join(pos |> transmute(ccode, year = year + 1L, .M_tm1 = M, .R_tm1 = R_pr), by = c("ccode", "year")) |>
    inner_join(pos |> transmute(ccode, year = year + 1L + h, .M_w0 = M, .R_w0 = R_pr), by = c("ccode", "year")) |>
    mutate(dM_raw_h = .M_tm1 - .M_w0, dRpr_raw_h = .R_tm1 - .R_w0) |>
    group_by(year) |>
    mutate(.dM_c = dM_raw_h - mean(dM_raw_h), .dR_c = dRpr_raw_h - mean(dRpr_raw_h)) |> ungroup()
  if (nrow(tbl) < 3) stop(sprintf("change table h=%d: only %d state-years have levels at t-1 and t-1-h", h, nrow(tbl)))
  sd_M <- sd(tbl$.dM_c); sd_R <- sd(tbl$.dR_c)
  if (!is.finite(sd_M) || sd_M <= 0 || !is.finite(sd_R) || sd_R <= 0)
    stop(sprintf("change table h=%d: pooled SD of the year-centred change is not positive (M %.3g, R_pr %.3g)", h, sd_M, sd_R))
  tbl <- tbl |> mutate(dM_h = .dM_c / sd_M, dRpr_h = .dR_c / sd_R, h = h) |>
    select(ccode, year, h, dM_raw_h, dRpr_raw_h, dM_h, dRpr_h) |> arrange(year, ccode)
  assert_no_dup_keys(tbl, c("ccode", "year"), sprintf("change table h=%d", h))
  assert_finite(tbl, c("dM_raw_h", "dRpr_raw_h", "dM_h", "dRpr_h"), sprintf("change table h=%d", h))
  ## invariants: within-year mean 0 (tolerance), pooled SD exactly 1 after scaling
  yr <- tbl |> group_by(year) |> summarise(mM = mean(dM_h), mR = mean(dRpr_h), .groups = "drop")
  if (any(abs(yr$mM) > 1e-8) || any(abs(yr$mR) > 1e-8))
    stop(sprintf("change table h=%d: within-year mean of the centred change is not 0 (max |mean| M %.2e, R_pr %.2e)", h, max(abs(yr$mM)), max(abs(yr$mR))))
  if (abs(sd(tbl$dM_h) - 1) > 1e-8 || abs(sd(tbl$dRpr_h) - 1) > 1e-8)
    stop(sprintf("change table h=%d: pooled SD after scaling is not 1 (M %.6f, R_pr %.6f)", h, sd(tbl$dM_h), sd(tbl$dRpr_h)))
  list(tbl = tbl,
       scale = data.frame(h = h, variable = c("dM_h", "dRpr_h"), pooled_sd = c(sd_M, sd_R),
                          n_rows = nrow(tbl), n_years = n_distinct(tbl$year)))
}

## ---- fatal-MID coding (PI decision 2026-10-23; see README "Samples") -----------------
## fatal_cow : COW coding (TRUE / FALSE / NA when every COW fatality code is missing)
## fatal_mic : MIC 1.0 confrontation fatalmin >= 1 (NA when the dispute is not in MIC)
## war       : COW dispute hostility level 5
## fatal_gml : GML MID 2.2.1 coding of the same directed pair (NA when missing or not in GML)
## Main fmid = war (fatal by definition) OR COW fatal OR (COW missing AND MIC fatal)
##             OR (COW missing AND MIC missing AND GML fatal).
## Each source only fills cases the earlier ones leave missing; none overrides a coded value.
## Returns logical, never NA.
code_fatal_main <- function(fatal_cow, fatal_mic, war, fatal_gml) {
  n <- length(fatal_cow)
  stopifnot(length(fatal_mic) == n, length(war) == n, length(fatal_gml) == n)
  (war %in% TRUE) | (fatal_cow %in% TRUE) | (is.na(fatal_cow) & fatal_mic %in% TRUE) |
    (is.na(fatal_cow) & is.na(fatal_mic) & fatal_gml %in% TRUE)
}
MIC_CONF_SHA256 <- "41cef02d0dcb23f91658b2367864b018a517e13d5b0278394dc12775fa5a5ee0"   # mic-conf-1.0.csv

## ---- dispute participation (PI decisions 2026-10-24; README "Events", "Ongoing disputes") ----------------------
## dd: COW-MID directed dispute-years (cow_mid_dirdisps) with integer dispnum, ccode1, ccode2, year. The data hold
## both directions of every opposing pair for every year the pair is in the dispute. Every function below uses
## EVERY participant (originators and joiners, both sides), so the already-fighting rule, the stage-1 control
## ongoing_tm1 and the dyadic lagged DV dy_mid_tm1 count joiner involvement the same way the events do.
## pair_dispute_years(): one row per (dispnum, ccode1, ccode2, year) with pair_onset = the pair's first year in
## that dispute (the events' onset year uses the same definition).
pair_dispute_years <- function(dd) {
  stopifnot(all(c("dispnum", "ccode1", "ccode2", "year") %in% names(dd)))
  out <- dd |> distinct(dispnum, ccode1, ccode2, year) |> group_by(dispnum, ccode1, ccode2) |>
    mutate(pair_onset = min(year)) |> ungroup()
  assert_no_dup_keys(out, c("dispnum", "ccode1", "ccode2", "year"), "pair dispute-years")
  out
}
## already_fighting_keys(): directed (ccode1, ccode2, year, dispnum) such that the pair is a participant pair of
## dispute dispnum, the dispute is active for the pair in `year` and the pair's onset in it is BEFORE `year`.
## A stage-2 alternative (i, j, t) matching any row is "already fighting" and is dropped from the choice set.
already_fighting_keys <- function(dd) {
  pair_dispute_years(dd) |> filter(year > pair_onset) |> select(ccode1, ccode2, year, dispnum)
}
## flag_escalations(): events (initiator, target, year, dispnum) get escalation = TRUE when (initiator, target, year)
## is already fighting in a dispute OTHER than the event's own. The own dispute can never match because the event
## year is the pair's onset in it, and this is asserted.
flag_escalations <- function(events, af) {
  stopifnot(all(c("initiator", "target", "year", "dispnum") %in% names(events)), all(c("ccode1", "ccode2", "year", "dispnum") %in% names(af)))
  hit <- events |> select(initiator, target, year, dispnum) |>
    inner_join(af |> rename(initiator = ccode1, target = ccode2, af_dispnum = dispnum), by = c("initiator", "target", "year"), relationship = "many-to-many")
  if (nrow(hit) > 0 && any(hit$dispnum == hit$af_dispnum)) stop("flag_escalations: an event matched an already-fighting row of its OWN dispute (onset year is not the pair's first year)")
  esc <- hit |> distinct(initiator, target, year, dispnum) |> mutate(escalation = TRUE)
  out <- events |> left_join(esc, by = c("initiator", "target", "year", "dispnum")) |> mutate(escalation = coalesce(escalation, FALSE))
  stopifnot(nrow(out) == nrow(events), !anyNA(out$escalation))
  out
}
## active_state_years(): (ccode, year) for every state that is a participant (any role) in a dispute active that year.
## Joined at year t-1 it is the stage-1 control ongoing_tm1.
active_state_years <- function(dd) dd |> distinct(ccode = ccode1, year) |> mutate(active = 1L)
## active_pair_years(): (ccode1, ccode2, year), both directions, for every pair of opposing participants in a dispute
## active that year. Joined at year t-1 it is the dyadic lagged DV dy_mid_tm1.
active_pair_years <- function(dd) {
  a <- dd |> distinct(ccode1, ccode2, year)
  bind_rows(a, a |> rename(ccode1 = ccode2, ccode2 = ccode1)) |> distinct() |> mutate(active = 1L)
}

## ---- month-level escalation rule (PI decision 2026-10-24, second round; README "Ongoing disputes") ---------------
## The year rule above mislabels disputes that ended early in the year before a new one started (China -> USA 1950,
## Korea, was dropped because a short 1949 incident had a 1950 row) and misses within-year re-entries. The month rule
## works with COW dispute dates at month resolution. peacesciencer 1.2.0 supplies, in cow_mid_disps, the dispute
## START month and year (stmon, styear) and the dispute duration in days (mindur, maxdur) but no end month and no
## participant-level dates; cow_mid_dirdisps has only the year. So every date here is DISPUTE-level, and the end
## date is bounded, not known: the start day is unknown, so the end lies between (first day of the start month +
## mindur - 1) and (last day of the start month + maxdur), intersected with the dispute's last year in the directed
## dispute-years. Every date is an interval [lo, hi] of month indices (ym_index); exact dates have lo == hi.
## A pair that enters a dispute after the dispute's start year has an unknown entry month (interval = that whole
## year), and a pair that leaves before the dispute's last year has an unknown exit month. Comparisons between
## intervals are three-valued (TRUE / FALSE / NA); NA falls back to the year rule and is counted.
ym_index <- function(year, month) {
  stopifnot(all(month %in% 1:12))
  as.integer(year) * 12L + as.integer(month) - 1L
}
ym_year <- function(ym) ym %/% 12L
ym_label <- function(ym) sprintf("%d-%02d", ym %/% 12L, ym %% 12L + 1L)
## "a is before b" (a < b) and "a is not before b" (a >= b) for intervals; NA when the intervals overlap
interval_before <- function(a_lo, a_hi, b_lo, b_hi) case_when(a_hi < b_lo ~ TRUE, a_lo >= b_hi ~ FALSE, TRUE ~ NA)
interval_not_before <- function(a_lo, a_hi, b_lo, b_hi) case_when(a_lo >= b_hi ~ TRUE, a_hi < b_lo ~ FALSE, TRUE ~ NA)
## dispute_date_bounds(disp, dd): one row per dispute in dd with start_lo / start_hi (exact when styear equals the first
## directed dispute-year, else the whole first year with start_month_known = FALSE), end_lo / end_hi (duration bounds
## within the last year; the whole last year with end_rule = "year_only" when the bounds miss that year), y0, y1.
dispute_date_bounds <- function(disp, dd) {
  stopifnot(all(c("dispnum", "styear", "stmon", "mindur", "maxdur") %in% names(disp)), all(c("dispnum", "year") %in% names(dd)))
  yrs <- dd |> distinct(dispnum, year) |> group_by(dispnum) |> summarise(y0 = min(year), y1 = max(year), .groups = "drop")
  db <- disp |> transmute(dispnum = as.integer(dispnum), styear = as.integer(styear), stmon = as.integer(stmon),
                          mindur = as.numeric(mindur), maxdur = as.numeric(maxdur)) |> inner_join(yrs, by = "dispnum")
  if (nrow(db) != nrow(yrs)) stop(sprintf("dispute_date_bounds: %d disputes in the directed dispute-years have no dispute-level row", nrow(yrs) - nrow(db)))
  assert_no_dup_keys(db, "dispnum", "dispute date bounds")
  if (!all(db$stmon %in% 1:12)) stop(sprintf("dispute_date_bounds: %d disputes have a start month outside 1..12", sum(!db$stmon %in% 1:12)))
  if (any(db$mindur < 0 | db$maxdur < db$mindur)) stop("dispute_date_bounds: negative or inverted durations")
  first_day <- as.Date(sprintf("%d-%02d-01", db$styear, db$stmon))
  next_first <- as.Date(sprintf("%d-%02d-01", db$styear + (db$stmon == 12L), ifelse(db$stmon == 12L, 1L, db$stmon + 1L)))
  end_lo_d <- first_day + pmax(db$mindur - 1, 0); end_hi_d <- (next_first - 1) + db$maxdur
  to_ym <- function(d) ym_index(as.integer(format(d, "%Y")), as.integer(format(d, "%m")))
  db <- db |> mutate(start_month_known = styear == y0,
                     start_lo = ifelse(start_month_known, ym_index(styear, stmon), ym_index(y0, 1L)),
                     start_hi = ifelse(start_month_known, ym_index(styear, stmon), ym_index(y0, 12L)),
                     end_lo = pmax(to_ym(end_lo_d), ym_index(y1, 1L)), end_hi = pmin(to_ym(end_hi_d), ym_index(y1, 12L)),
                     end_rule = ifelse(end_lo <= end_hi, "duration", "year_only"))
  bad <- db$end_rule == "year_only"
  db$end_lo[bad] <- ym_index(db$y1[bad], 1L); db$end_hi[bad] <- ym_index(db$y1[bad], 12L)
  stopifnot(all(db$start_lo <= db$start_hi), all(db$end_lo <= db$end_hi), all(db$start_lo <= db$end_hi),
            all(ym_year(db$start_lo) == db$y0), all(ym_year(db$end_hi) == db$y1))
  db |> select(dispnum, styear, stmon, mindur, maxdur, y0, y1, start_month_known, start_lo, start_hi, end_lo, end_hi, end_rule)
}
## pair_dispute_intervals(dd, db): one row per directed (dispnum, ccode1, ccode2, spell), a spell being a run of
## consecutive years the pair is in the dispute (a pair that leaves and re-enters, e.g. France and Germany in dispute
## 258, has two spells; most pairs have one). pair_onset / pair_last are the whole pair's first and last year (the year
## rule's quantities), spell_first / spell_last the spell's. The spell's entry interval (start_lo / start_hi) is the
## dispute's start when the spell begins in the dispute's first year AND both members are originators (orig == 1 in the
## directed dispute-years; Codex review P1-2, 2026-10-24). When either member is a joiner its entry date within the first
## year is unknown, so the entry interval is [dispute start, December of y0] with start_exact = FALSE (the month rule's
## NA then falls back to the year rule). Spells beginning after the dispute's first year get the whole first year of the
## spell. The exit interval (end_lo / end_hi) is the dispute's end bounds when the spell ends in the dispute's last year,
## else the whole last year of the spell. both_orig is kept in the output for logging.
## joiner_entry_unknown = FALSE reproduces the rule before the review (joiner pairs get the exact dispute start); the build
## uses it only to log the before / after counts.
## Participant-level dates (PI decision 2026-10-25, D5). part = gml_part (peacesciencer 1.2.0, bundled): one row per
## participant SPELL with start and end year / month (dispnum, ccode, styear, stmon, endyear, endmon; the day is often -9 and
## is not used, every date here is a month index). About 2165 of the 2436 COW dispnums are in it. For a pair spell whose
## dispute is in part, the pair's entry is the LATER of the two participants' starts and its exit the EARLIER of their ends,
## both exact months, provided (a) each participant has exactly one gml_part spell overlapping the pair spell's years, (b) the
## entry month lies in the spell's first COW year and the exit month in its last COW year (GML and COW disagree on the
## participation years of some pairs; those keep the dispute-level bounds), and (c) entry <= exit. Otherwise the dispute-level
## bounds above apply; date_source records which ("participant" / "dispute") and part_reason why not.
## Clipping (D5): an unknown entry bound can exceed the feasible exit (a joiner pair entering some time in y0 when the dispute
## ended in y0 before December): start_hi is clipped to end_hi; start_clipped marks those spells.
pair_dispute_intervals <- function(dd, db, joiner_entry_unknown = TRUE, part = NULL) {
  stopifnot(all(c("orig1", "orig2") %in% names(dd)))
  orig <- dd |> group_by(dispnum, ccode1, ccode2) |>
    summarise(n_orig = n_distinct(orig1 == 1L & orig2 == 1L), both_orig = first(orig1 == 1L & orig2 == 1L), .groups = "drop")
  if (any(orig$n_orig != 1L)) stop("pair_dispute_intervals: originator status of a pair differs across its dispute-years")
  stopifnot(!anyNA(orig$both_orig))
  py <- pair_dispute_years(dd) |> arrange(dispnum, ccode1, ccode2, year) |> group_by(dispnum, ccode1, ccode2) |>
    mutate(spell = cumsum(c(1L, diff(year) != 1L)), pair_onset = min(year), pair_last = max(year)) |> ungroup()
  pdi <- py |> group_by(dispnum, ccode1, ccode2, spell, pair_onset, pair_last) |>
    summarise(spell_first = min(year), spell_last = max(year), n_years = n(), .groups = "drop") |>
    inner_join(db |> select(dispnum, y0, y1, d_start_lo = start_lo, d_start_hi = start_hi, d_end_lo = end_lo, d_end_hi = end_hi), by = "dispnum") |>
    inner_join(orig |> select(dispnum, ccode1, ccode2, both_orig), by = c("dispnum", "ccode1", "ccode2"))
  stopifnot(all(pdi$n_years == pdi$spell_last - pdi$spell_first + 1L), all(pdi$spell_first >= pdi$y0), all(pdi$spell_last <= pdi$y1),
            !anyNA(pdi$both_orig))
  exact_entry <- if (joiner_entry_unknown) pdi$both_orig else rep(TRUE, nrow(pdi))
  pdi <- pdi |> mutate(start_exact = spell_first == y0 & d_start_lo == d_start_hi & exact_entry,
                       start_lo = ifelse(spell_first == y0, d_start_lo, ym_index(spell_first, 1L)),
                       start_hi = ifelse(spell_first == y0 & exact_entry, d_start_hi, ym_index(spell_first, 12L)),
                       end_lo = ifelse(spell_last == y1, d_end_lo, ym_index(spell_last, 1L)),
                       end_hi = ifelse(spell_last == y1, d_end_hi, ym_index(spell_last, 12L))) |>
    select(dispnum, ccode1, ccode2, spell, pair_onset, pair_last, spell_first, spell_last, both_orig, start_exact, start_lo, start_hi, end_lo, end_hi) |>
    mutate(date_source = "dispute", part_reason = if (is.null(part)) "no participant data" else "dispute not in gml_part")
  stopifnot(all(pdi$start_lo <= pdi$start_hi), all(pdi$end_lo <= pdi$end_hi), all(pdi$start_lo <= pdi$end_hi), !anyNA(pdi$start_lo), !anyNA(pdi$end_hi),
            all(pdi$spell_first[pdi$spell == 1L] == pdi$pair_onset[pdi$spell == 1L]))
  ## ---- participant-level dates (D5) -----------------------------------------------------------------------------------------
  if (!is.null(part)) {
    stopifnot(all(c("dispnum", "ccode", "styear", "stmon", "endyear", "endmon") %in% names(part)))
    ps <- part |> transmute(dispnum = as.integer(dispnum), ccode = as.integer(ccode), p_y0 = as.integer(styear), p_y1 = as.integer(endyear),
                            p_start = ym_index(as.integer(styear), as.integer(stmon)), p_end = ym_index(as.integer(endyear), as.integer(endmon)))
    stopifnot(!anyNA(ps), all(ps$p_start <= ps$p_end))
    in_part <- pdi$dispnum %in% ps$dispnum
    ## participant spells of ccode1 / ccode2 overlapping the pair spell's COW years; count them per pair spell
    one_side <- function(code_col, pre) {
      pdi |> mutate(.id = row_number()) |> select(.id, dispnum, ccode = all_of(code_col), spell_first, spell_last) |>
        inner_join(ps, by = c("dispnum", "ccode"), relationship = "many-to-many") |>
        filter(p_y1 >= spell_first, p_y0 <= spell_last) |> group_by(.id) |>
        summarise(n = n(), start = p_start[1], end = p_end[1], .groups = "drop") |> rename_with(~ paste0(pre, .x), -.id)
    }
    s1 <- one_side("ccode1", "a_"); s2 <- one_side("ccode2", "b_")
    pp <- pdi |> mutate(.id = row_number()) |> left_join(s1, by = ".id") |> left_join(s2, by = ".id") |>
      mutate(a_n = coalesce(a_n, 0L), b_n = coalesce(b_n, 0L),
             p_entry_raw = pmax(a_start, b_start), p_exit_raw = pmin(a_end, b_end),
             ## GML and COW disagree on some participation years (e.g. Egypt and Israel in dispute 1035 start 1966-12 in gml_part but
             ## their COW pair-years begin in 1967). An entry before the spell's first COW year means the pair was in the dispute at
             ## January of that year at the latest: clamp to January (exact at the COW year's start); an exit after the spell's last
             ## COW year is clamped to December of it. An entry AFTER the first COW year or an exit BEFORE the last one contradicts
             ## the COW years and falls back to the dispute-level bounds.
             p_entry = pmax(p_entry_raw, ym_index(spell_first, 1L)), p_exit = pmin(p_exit_raw, ym_index(spell_last, 12L)),
             part_reason = case_when(!in_part ~ "dispute not in gml_part",
                                     a_n == 0L | b_n == 0L ~ "a participant has no gml_part spell overlapping the pair spell",
                                     a_n > 1L | b_n > 1L ~ "a participant has several gml_part spells overlapping the pair spell",
                                     ym_year(p_entry_raw) > spell_first ~ "participant entry after the spell's first COW year",
                                     ym_year(p_exit_raw) < spell_last ~ "participant exit before the spell's last COW year",
                                     p_entry > p_exit ~ "participant entry after exit",
                                     ym_year(p_entry_raw) < spell_first & ym_year(p_exit_raw) > spell_last ~ "participant dates used (entry and exit clamped to the COW years)",
                                     ym_year(p_entry_raw) < spell_first ~ "participant dates used (entry clamped to January of the first COW year)",
                                     ym_year(p_exit_raw) > spell_last ~ "participant dates used (exit clamped to December of the last COW year)",
                                     TRUE ~ "participant dates used"))
    use <- startsWith(pp$part_reason, "participant dates used")
    pdi$start_lo[use] <- pp$p_entry[use]; pdi$start_hi[use] <- pp$p_entry[use]
    pdi$end_lo[use] <- pp$p_exit[use]; pdi$end_hi[use] <- pp$p_exit[use]
    pdi$start_exact[use] <- TRUE; pdi$date_source[use] <- "participant"; pdi$part_reason <- pp$part_reason
    stopifnot(all(ym_year(pdi$start_lo) == pdi$spell_first), all(ym_year(pdi$end_hi) == pdi$spell_last), all(pdi$start_lo <= pdi$end_hi))
  }
  ## ---- clip an unknown entry bound by the feasible exit (D5) ---------------------------------------------------------------
  pdi <- pdi |> mutate(start_clipped = start_hi > end_hi, start_hi = pmin(start_hi, end_hi))
  stopifnot(all(pdi$start_lo <= pdi$start_hi), all(pdi$start_hi <= pdi$end_hi), all(pdi$end_lo <= pdi$end_hi), !any(pdi$start_clipped & pdi$date_source == "participant"))
  assert_no_dup_keys(pdi, c("dispnum", "ccode1", "ccode2", "spell"), "pair dispute intervals")
  pdi
}
## month_rule(): is the pair "already fighting" in dispute D at the reference interval [ref_lo, ref_hi]? TRUE when D's
## pair start is before the reference and its pair end is not before it; FALSE when either fails for sure; NA otherwise,
## in which case the caller falls back to the year rule.
month_rule <- function(start_lo, start_hi, end_lo, end_hi, ref_lo, ref_hi) {
  sb <- interval_before(start_lo, start_hi, ref_lo, ref_hi)
  ne <- interval_not_before(end_lo, end_hi, ref_lo, ref_hi)
  case_when(sb %in% FALSE | ne %in% FALSE ~ FALSE, sb %in% TRUE & ne %in% TRUE ~ TRUE, TRUE ~ NA)
}
## flag_escalations_month(events, pdi): events (initiator, target, year, dispnum) get
##   ev_start_lo / ev_start_hi : the event's date = the pair's entry interval in its own dispute
##   escalation                : TRUE when the pair is already fighting in ANOTHER dispute at that date (month rule, year
##                               rule where the month rule is undecided)
##   escalation_year           : the year rule alone (pair in another dispute active in `year` with pair onset before it)
##   esc_fallback              : TRUE when the decision rests on at least one year-rule fallback (no comparison was a
##                               certain TRUE and at least one was undecided)
## attr(, "comparisons") holds the per-(event, other dispute) table for logging.
flag_escalations_month <- function(events, pdi) {
  stopifnot(all(c("initiator", "target", "year", "dispnum") %in% names(events)))
  own <- pdi |> filter(spell == 1L) |> select(dispnum, initiator = ccode1, target = ccode2, pair_onset, ev_start_lo = start_lo, ev_start_hi = start_hi)
  ev <- events |> select(initiator, target, year, dispnum) |> inner_join(own, by = c("dispnum", "initiator", "target"))
  if (nrow(ev) != nrow(events)) stop("flag_escalations_month: an event has no pair-dispute interval")
  if (any(ev$pair_onset != ev$year)) stop("flag_escalations_month: an event's year is not the pair's onset year in its dispute")
  ## one comparison per (event, other dispute, spell of the pair in it); the year rule is the whole pair's onset and the
  ## spell's years, so any() over spells reproduces already_fighting_keys() exactly
  cmp <- ev |> select(-pair_onset) |>
    inner_join(pdi |> select(initiator = ccode1, target = ccode2, other_dispnum = dispnum, o_spell = spell, o_onset = pair_onset,
                             o_first = spell_first, o_last = spell_last, o_start_lo = start_lo, o_start_hi = start_hi, o_end_lo = end_lo, o_end_hi = end_hi),
               by = c("initiator", "target"), relationship = "many-to-many") |>
    filter(other_dispnum != dispnum) |>
    mutate(month_dec = month_rule(o_start_lo, o_start_hi, o_end_lo, o_end_hi, ev_start_lo, ev_start_hi),
           year_dec = o_onset < year & o_first <= year & o_last >= year,
           fallback = is.na(month_dec), decision = ifelse(fallback, year_dec, month_dec))
  stopifnot(!anyNA(cmp$decision), !anyNA(cmp$year_dec))
  agg <- cmp |> group_by(initiator, target, year, dispnum) |>
    summarise(escalation = any(decision), escalation_year = any(year_dec),
              esc_fallback = !any(month_dec %in% TRUE) & any(fallback), .groups = "drop")
  out <- events |> left_join(ev |> select(initiator, target, year, dispnum, ev_start_lo, ev_start_hi), by = c("initiator", "target", "year", "dispnum")) |>
    left_join(agg, by = c("initiator", "target", "year", "dispnum")) |>
    mutate(escalation = coalesce(escalation, FALSE), escalation_year = coalesce(escalation_year, FALSE), esc_fallback = coalesce(esc_fallback, FALSE))
  stopifnot(nrow(out) == nrow(events), !anyNA(out$escalation), !anyNA(out$ev_start_lo), !anyNA(out$ev_start_hi))
  ## the year rule here must equal the earlier already_fighting_keys() / flag_escalations() path
  attr(out, "comparisons") <- cmp
  out
}
## reference_dates(events): for every (ccode, year) that initiates (has a non-escalation event), the interval of its
## EARLIEST non-escalation event's start (lo = min of the lower bounds, hi = min of the upper bounds). State-years absent
## from this table use January of year t (the caller coalesces).
reference_dates <- function(events) {
  stopifnot(all(c("initiator", "year", "escalation", "ev_start_lo", "ev_start_hi") %in% names(events)))
  out <- events |> filter(!escalation) |> group_by(ccode = initiator, year) |>
    summarise(ref_lo = min(ev_start_lo), ref_hi = min(ev_start_hi), .groups = "drop")
  stopifnot(all(out$ref_lo <= out$ref_hi), all(ym_year(out$ref_lo) == out$year))
  out
}
## already_fighting_keys_month(pdi, ref): directed (ccode1, ccode2, year, dispnum) rows for every pair-dispute and
## every year the pair is in it, with drop = TRUE when the pair is already fighting in that dispute at the state-year's
## reference date (month rule, year fallback), drop_year = the year rule, fallback flag. The reference date of
## (ccode1, year) is its row in ref, else January of year. ref = NULL (Codex review P1-1, 2026-10-24) uses January of
## year for EVERY state-year: this is the stage-1 IV choice set, which must not depend on the state-year's own events
## (the reference-date set is the stage-2 set). No information from after the reference date is used.
already_fighting_keys_month <- function(pdi, ref = NULL) {
  if (is.null(ref)) ref <- tibble(ccode = integer(0), year = integer(0), ref_lo = integer(0), ref_hi = integer(0))
  stopifnot(all(c("ccode", "year", "ref_lo", "ref_hi") %in% names(ref)))
  keys <- pdi |> rowwise() |> mutate(year = list(seq.int(spell_first, spell_last))) |> ungroup() |> tidyr::unnest(year) |>
    left_join(ref |> rename(ccode1 = ccode), by = c("ccode1", "year")) |>
    mutate(ref_lo = coalesce(ref_lo, ym_index(year, 1L)), ref_hi = coalesce(ref_hi, ym_index(year, 1L)),
           month_dec = month_rule(start_lo, start_hi, end_lo, end_hi, ref_lo, ref_hi),
           drop_year = year > pair_onset,
           fallback = is.na(month_dec), drop = ifelse(fallback, drop_year, month_dec))
  stopifnot(!anyNA(keys$drop))
  keys |> select(ccode1, ccode2, year, dispnum, spell, ref_lo, ref_hi, drop, drop_year, fallback)
}
## choice_set_flags(alts, af_iv, af_s2) (Codex review P1-1, PI choice "A'", 2026-10-24). alts has one row per
## (ccode1, ccode2, year); af_iv / af_s2 are the distinct directed (ccode1, ccode2, year) keys already fighting at
## January (stage-1 IV set) and at the state-year's reference date (stage-2 set). Returns alts with
##   in_iv : 1 when the row is in the stage-1 IV set (not in af_iv)
##   in_s2 : 1 when the row is in the stage-2 choice set (not in af_s2)
## and WITHOUT the rows that are in neither set. One row block per state-year is the UNION of both sets; the Stan
## programs take IV[n] over the in_iv rows and the stage-2 softmax over the in_s2 rows. For a non-initiating
## state-year the two sets are identical (the reference date is January in both). in_iv never depends on af_s2.
choice_set_flags <- function(alts, af_iv, af_s2) {
  stopifnot(all(c("ccode1", "ccode2", "year") %in% names(alts)), all(c("ccode1", "ccode2", "year") %in% names(af_iv)),
            all(c("ccode1", "ccode2", "year") %in% names(af_s2)), !any(c("in_iv", "in_s2") %in% names(alts)))
  assert_no_dup_keys(alts, c("ccode1", "ccode2", "year"), "choice_set_flags: alts")
  k_alt <- paste(alts$ccode1, alts$ccode2, alts$year)
  k_iv <- unique(paste(af_iv$ccode1, af_iv$ccode2, af_iv$year)); k_s2 <- unique(paste(af_s2$ccode1, af_s2$ccode2, af_s2$year))
  out <- alts |> mutate(in_iv = as.integer(!(k_alt %in% k_iv)), in_s2 = as.integer(!(k_alt %in% k_s2)))
  n_neither <- sum(out$in_iv == 0L & out$in_s2 == 0L)
  out <- out |> filter(in_iv == 1L | in_s2 == 1L)
  stopifnot(nrow(out) == nrow(alts) - n_neither, all(out$in_iv %in% 0:1), all(out$in_s2 %in% 0:1), all(out$in_iv == 1L | out$in_s2 == 1L))
  attr(out, "n_neither") <- n_neither
  out
}
## collapse_duplicate_events(events): one event per (initiator, target, year) among the non-escalation rows (PI decision
## 2026-10-24). Sample flags and codings are ORed, is_joiner_pair is FALSE when any collapsed row is an originator pair,
## dispnum keeps the smallest dispute number and dispnums lists them all. Escalation rows are kept as they are. Every
## column must be named in one of the rules below, so a new column fails loudly instead of being collapsed by accident.
## event_weights(events) (PI decision 2026-10-25, D8; README "Coalition expansion"): per sample s, w_s = 1 / (number of the
## initiator's targets in the SAME dispute and onset year that are in sample s), so an initiator facing k Side-B members of one
## dispute contributes one decision's weight to the stage-2 composite likelihood; targets in different disputes stay separate
## decisions. 0 for rows outside the sample. Computed on the pre-collapse rows (one per dispute pair); the collapse sums the
## weights of a key's source rows and caps at 1.
event_weights <- function(events, samples = c("mid", "fmid", "fmidcow", "fmidub", "war")) {
  stopifnot(all(c("dispnum", "initiator", "year") %in% names(events)), all(paste0("in_", samples) %in% names(events)))
  out <- events
  for (s in samples) {
    f <- out[[paste0("in_", s)]]
    stopifnot(!anyNA(f), is.logical(f))
    n <- out |> mutate(.f = f) |> group_by(dispnum, initiator, year) |> mutate(.n = sum(.f)) |> ungroup() |> pull(.n)
    w <- ifelse(f, 1 / n, 0)
    stopifnot(all(is.finite(w)), all(w[f] > 0), all(w[f] <= 1), all(w[!f] == 0))
    tot <- out |> mutate(.w = w) |> group_by(dispnum, initiator, year) |> summarise(.s = sum(.w), .groups = "drop")
    stopifnot(all(abs(tot$.s - 1) < 1e-12 | tot$.s == 0))
    out[[paste0("w_", s)]] <- w
  }
  out
}
## Direction columns (2026-10-25, D1): dir_flipped is ORed, dir_source / dir_rule are pasted, war_link and n_war_match take the
## max; sidea1 (initiator on MID Side A) is TRUE only when it holds in EVERY source dispute and sidea2 (initiator on MID Side B,
## i.e. flipped) when it holds in ANY, so sidea1 + sidea2 == 1 still marks a key whose sources agree. Event weights (D8): the
## per-sample stage-2 weights w_<sample> are summed over the collapsed rows and capped at 1 (one chosen alternative can carry
## at most one decision's weight).
collapse_duplicate_events <- function(events) {
  key_cols   <- c("initiator", "target", "year")
  same_cols  <- c("target_has_pos_tm1", "fighting_at_ref")                                     # must agree within a key
  all_cols   <- c("sidea1"); any_cols <- c("sidea2")                                           # 0/1, see above
  wsum_cols  <- c("w_mid", "w_fmid", "w_fmidcow", "w_fmidub", "w_war")                         # min(1, sum)
  ## the *_noww2 flags (Codex review P2-3) are the sample flags of the NON-WWII source rows only: ORing them over a key
  ## that mixes a WWII and a non-WWII dispute keeps the event in the no-WWII sample iff a non-WWII dispute qualifies
  or_cols    <- c("war", "war_disp", "escalation_year", "esc_fallback", "is_ww2", "fatal_main", "fatal_cow", "fatal_ub",
                  "in_mid", "in_fmid", "in_fmidcow", "in_fmidub", "in_war",
                  "in_mid_noww2", "in_fmid_noww2", "in_fmidcow_noww2", "in_fmidub_noww2", "in_war_noww2", "dir_flipped")   # never NA
  tri_cols   <- c("fatal1", "fatal2", "fatal_disp", "fatal", "fatal_mic", "fatal_gml")          # TRUE if any TRUE, NA if all NA
  max_cols   <- c("orig1", "orig2", "fatality1", "fatality2", "hostlev1", "hostlev2", "fatality_disp", "hostlev_disp", "mic_fatalmin",
                  "war_link", "n_war_match")
  min_cols   <- c("dispnum", "ev_start_lo", "ev_start_hi")
  paste_cols <- c("fatal_source", "dir_source", "dir_rule")
  known <- c(key_cols, same_cols, all_cols, any_cols, wsum_cols, or_cols, tri_cols, max_cols, min_cols, paste_cols, "is_joiner_pair", "escalation")
  extra <- setdiff(names(events), known); missing <- setdiff(known, names(events))
  if (length(extra)) stop("collapse_duplicate_events: unclassified columns: ", paste(extra, collapse = ", "))
  if (length(missing)) stop("collapse_duplicate_events: expected columns absent: ", paste(missing, collapse = ", "))
  for (cc in c(or_cols, all_cols, any_cols, wsum_cols)) if (anyNA(events[[cc]])) stop("collapse_duplicate_events: NA in ", cc)
  for (cc in c(all_cols, any_cols)) if (!all(events[[cc]] %in% 0:1)) stop("collapse_duplicate_events: ", cc, " is not 0/1")
  for (cc in wsum_cols) if (any(events[[cc]] < 0 | events[[cc]] > 1)) stop("collapse_duplicate_events: ", cc, " outside [0, 1]")
  esc <- events |> filter(escalation) |> mutate(dispnums = as.character(dispnum), n_disputes = 1L)
  keep <- events |> filter(!escalation)
  tri <- function(x) if (all(is.na(x))) NA else any(x, na.rm = TRUE)
  mx <- function(x) if (all(is.na(x))) x[1] else max(x, na.rm = TRUE)
  same <- function(x) { if (length(unique(x)) != 1) stop("collapse_duplicate_events: a 'same' column differs within a key"); x[1] }
  ## dispnums and n_disputes come FIRST: inside summarise() a later expression sees the already-reduced dispnum
  out <- keep |> group_by(across(all_of(key_cols))) |>
    summarise(dispnums = paste(sort(dispnum), collapse = ";"), n_disputes = n(),
              across(all_of(same_cols), same), across(all_of(or_cols), ~ any(.x)), across(all_of(tri_cols), tri),
              across(all_of(max_cols), mx), across(all_of(min_cols), ~ min(.x)),
              across(all_of(paste_cols), ~ paste(sort(unique(.x)), collapse = ";")),
              across(all_of(all_cols), ~ as.integer(all(.x == 1L))), across(all_of(any_cols), ~ as.integer(any(.x == 1L))),
              across(all_of(wsum_cols), ~ pmin(1, sum(.x))),
              is_joiner_pair = all(is_joiner_pair), escalation = FALSE, .groups = "drop")
  stopifnot(nrow(out) == nrow(distinct(keep, initiator, target, year)), sum(out$n_disputes) == nrow(keep),
            all(out$in_mid), !any(out$escalation), all(out$w_mid > 0), all(out$w_mid <= 1),
            all(lengths(strsplit(out$dispnums, ";", fixed = TRUE)) == out$n_disputes),
            all(vapply(strsplit(out$dispnums, ";", fixed = TRUE), function(z) as.integer(z[1]), 1L) == out$dispnum),
            setequal(as.integer(unlist(strsplit(out$dispnums, ";", fixed = TRUE))), keep$dispnum),
            length(unlist(strsplit(out$dispnums, ";", fixed = TRUE))) == nrow(keep))
  ## ORed flags are exactly "any row had the flag", checked column by column against the uncollapsed rows
  for (cc in or_cols) {
    chk <- keep |> group_by(across(all_of(key_cols))) |> summarise(v = any(.data[[cc]]), .groups = "drop") |> inner_join(out |> select(all_of(key_cols), w = all_of(cc)), by = key_cols)
    stopifnot(nrow(chk) == nrow(out), all(chk$v == chk$w))
  }
  res <- bind_rows(out, esc |> select(all_of(names(out))))
  assert_no_dup_keys(res |> filter(!escalation), key_cols, "collapsed events")
  assert_no_dup_keys(res, c("dispnum", "initiator", "target"), "collapsed events (dispute keys)")
  res
}

## ---- initiator direction for war-level disputes (PI decision 2026-10-25; README "Initiator direction") -------------
## COW-MID Side A is the side that took the FIRST codable action, not the aggressor. In dispute 258 (the Second World War)
## Side A is Poland and the Allies, so the Side A rule gives Poland -> Germany and never Germany -> Poland. For disputes at
## war level (hostlev_disp == 5) that can be linked to a COW inter-state war the direction follows COW's war-INITIATING
## side instead (cow_war_inter: sidea1 / sidea2 are the war sides 1 / 2, initiator1 / initiator2 the state-level initiator
## flags; in every COW war exactly one side has initiators, asserted below).
## build_pair_table(dd, disp): one row per directed (dispnum, ccode1, ccode2) with the pair's onset year, the sides and
## originator flags of the EARLIEST year the pair appears (sidea can change when a state switches sides, e.g. France and the
## USSR in dispute 258 are opposed in 1939 and again, with the sides reversed, in 1941; the P3 fix replaces first() by the
## earliest year's value and reports the pairs where it differs), the max fatality / hostility over the pair's years and
## the dispute-level codings. attr(, "side_changes") lists the directed pairs whose sidea differs across their years.
build_pair_table <- function(dd, disp) {
  stopifnot(all(c("dispnum", "ccode1", "ccode2", "year", "sidea1", "sidea2", "orig1", "orig2", "fatality1", "fatality2", "hostlev1", "hostlev2") %in% names(dd)),
            all(c("dispnum", "fatality_disp", "hostlev_disp") %in% names(disp)))
  assert_no_dup_keys(disp, "dispnum", "dispute-level table")
  chg <- dd |> group_by(dispnum, ccode1, ccode2) |> summarise(n_side = n_distinct(sidea1), n_orig = n_distinct(orig1), .groups = "drop")
  if (any(chg$n_orig != 1L)) stop("build_pair_table: originator status of a pair differs across its dispute-years")
  pair <- dd |> arrange(dispnum, ccode1, ccode2, year) |> group_by(dispnum, ccode1, ccode2) |>
    summarise(year = min(year), sidea1 = sidea1[1], sidea2 = sidea2[1],            # arranged by year: the EARLIEST year's sides
              orig1 = orig1[1], orig2 = orig2[1],
              fatality1 = max(fatality1), fatality2 = max(fatality2),
              hostlev1 = max(hostlev1), hostlev2 = max(hostlev2), .groups = "drop") |>
    left_join(disp |> select(dispnum, fatality_disp, hostlev_disp), by = "dispnum")
  stopifnot(!anyNA(pair$hostlev_disp), all(pair$sidea1 %in% 0:1), all(pair$sidea2 %in% 0:1))
  attr(pair, "side_changes") <- chg |> filter(n_side > 1L) |> select(dispnum, ccode1, ccode2)
  pair
}
## war_dyad_years(war_inter): distinct COW inter-state war dyad-years with, per row, whether the side of ccode1 (side1_init)
## and of ccode2 (side2_init) is the war's initiating side. Exactly one side per war has initiators (asserted), so
## side2_init == !side1_init on every row.
war_dyad_years <- function(war_inter) {
  stopifnot(all(c("warnum", "ccode1", "ccode2", "year", "sidea1", "sidea2", "initiator1", "initiator2", "cowinteronset") %in% names(war_inter)))
  w <- war_inter |> transmute(warnum = as.integer(warnum), ccode1 = as.integer(ccode1), ccode2 = as.integer(ccode2), year = as.integer(year),
                              sidea1 = as.integer(sidea1), sidea2 = as.integer(sidea2), initiator1 = as.integer(initiator1),
                              initiator2 = as.integer(initiator2), cowinteronset = as.integer(cowinteronset)) |> distinct()
  stopifnot(!anyNA(w), all(w$sidea1 %in% 1:2), all(w$sidea2 %in% 1:2), all(w$sidea1 != w$sidea2), all(w$initiator1 %in% 0:1), all(w$initiator2 %in% 0:1))
  side_init <- w |> group_by(warnum, side = sidea1) |> summarise(side_has_init = any(initiator1 == 1L), .groups = "drop")
  chk <- side_init |> group_by(warnum) |> summarise(n = sum(side_has_init), .groups = "drop")
  if (any(chk$n != 1L)) stop("war_dyad_years: COW wars without exactly one initiating side: ", paste(chk$warnum[chk$n != 1L], collapse = ", "))
  out <- w |> left_join(side_init |> rename(side1_init = side_has_init), by = c("warnum", "sidea1" = "side")) |>
    left_join(side_init |> rename(side2_init = side_has_init), by = c("warnum", "sidea2" = "side"))
  stopifnot(!anyNA(out$side1_init), !anyNA(out$side2_init), all(out$side1_init != out$side2_init))
  out
}
## cow_war_direction(ev, wd, window, min_link): the direction decision for every event row of the pair table.
## ev needs dispnum, initiator, target, year, hostlev_disp; wd is war_dyad_years(). Rules (in this order):
##   pair level     : the pair (i, j) is a COW war dyad-year with year in [t, t + window] (t = the pair's MID onset year; MID
##                    onsets precede war onsets, so the window allows the war to begin in the following year). Every matched
##                    row votes: flip when j's war side initiated, keep when i's did. Unanimous -> decided, dir_source
##                    "cow_war_pair" (flipped or confirmed); mixed votes (two wars with opposite initiators) -> ambiguous,
##                    Side A kept.
##   dispute level  : for the remaining pairs of a war-level dispute the dispute is linked to the COW war most of its
##                    pair-decided events matched (strict majority, at least min_link events); the MID Side A state of every
##                    pair-decided event linked to that war is on war side sidea1, and the majority war side of MID Side A
##                    (again at least min_link events and a strict majority) is the mapping. If that side is the war's
##                    non-initiating side the unmatched pairs are flipped, else kept; dir_source "cow_war_dispute". No link,
##                    too few linked events or a split mapping -> ambiguous, Side A kept.
##   everything else: non-war disputes and unlinked war-level disputes keep Side A, dir_source "sideA".
## Returns ev with dir_flipped (logical), dir_source ("cow_war_pair" / "cow_war_dispute" / "sideA"), dir_rule (the reason),
## war_link (the COW war number behind the decision, NA otherwise), n_war_match (matched dyad-year rows), and the
## attributes "pair_matches" (every matched (event, war dyad-year) row) and "dispute_links" (one row per war-level dispute).
cow_war_direction <- function(ev, wd, window = 1L, min_link = 2L) {
  stopifnot(all(c("dispnum", "initiator", "target", "year", "hostlev_disp") %in% names(ev)), window >= 0, min_link >= 1)
  assert_no_dup_keys(ev, c("dispnum", "initiator", "target"), "cow_war_direction: events")
  key <- c("dispnum", "initiator", "target", "year")
  cand <- ev |> filter(hostlev_disp == 5L) |> select(all_of(key))
  m <- cand |> inner_join(wd |> rename(initiator = ccode1, target = ccode2, wyear = year), by = c("initiator", "target"), relationship = "many-to-many") |>
    filter(wyear >= year, wyear <= year + window) |>
    mutate(vote = ifelse(side2_init, "flip", "keep"))            # exactly one side initiated: j's side -> flip, i's side -> keep
  modal <- function(v) as.integer(names(which.max(table(v))))                # the most frequent value (ties: the smallest)
  dec <- if (nrow(m) > 0) m |> group_by(across(all_of(key))) |>
    summarise(n_war_match = n(), n_flip = sum(vote == "flip"), n_keep = sum(vote == "keep"),
              match_t = any(wyear == year), match_t1 = any(wyear > year), onset_in_window = any(cowinteronset == 1L),
              war_link = modal(warnum), n_wars = n_distinct(warnum), .groups = "drop") else
    tibble(dispnum = integer(0), initiator = integer(0), target = integer(0), year = integer(0), n_war_match = integer(0), n_flip = integer(0), n_keep = integer(0),
           match_t = logical(0), match_t1 = logical(0), onset_in_window = logical(0), war_link = integer(0), n_wars = integer(0))
  dec <- dec |> mutate(pair_dec = case_when(n_flip == n_war_match ~ "flip", n_keep == n_war_match ~ "keep", TRUE ~ "ambiguous"))
  stopifnot(all(dec$n_flip + dec$n_keep == dec$n_war_match))
  ## dispute link: modal war among the pair-decided events, strict majority and at least min_link of them
  decided <- dec |> filter(pair_dec != "ambiguous")
  if (nrow(decided) > 0) {
    link <- decided |> group_by(dispnum) |>
      summarise(n_pair_decided = n(), n_link = max(table(war_link)), war_link = modal(war_link), .groups = "drop") |>   # n_link BEFORE war_link is redefined
      mutate(link_ok = n_link >= min_link & 2L * n_link > n_pair_decided)
    ## side mapping: war side of the MID Side A state over the pair-decided events linked to war_link
    mp <- m |> semi_join(decided, by = key) |> inner_join(link |> select(dispnum, war_link), by = "dispnum") |> filter(warnum == war_link) |>
      distinct(across(all_of(key)), sidea1, side1_init) |> group_by(dispnum) |>
      summarise(n_map = n(), sideA_war_side = modal(sidea1), n_map_modal = max(table(sidea1)),
                sideA_war_side_init = any(side1_init[sidea1 == sideA_war_side]), .groups = "drop")
    stopifnot(all(mp$n_map_modal >= 1L), nrow(mp) == nrow(link))
  } else {
    link <- tibble(dispnum = integer(0), n_pair_decided = integer(0), n_link = integer(0), war_link = integer(0), link_ok = logical(0))
    mp <- tibble(dispnum = integer(0), n_map = integer(0), sideA_war_side = integer(0), n_map_modal = integer(0), sideA_war_side_init = logical(0))
  }
  link <- link |> left_join(mp, by = "dispnum") |>
    mutate(map_ok = link_ok & n_map_modal >= min_link & 2L * n_map_modal > n_map,
           dispute_dec = case_when(!link_ok ~ "no_majority_link", !map_ok ~ "split_side_mapping", !sideA_war_side_init ~ "flip", TRUE ~ "keep"))
  all_war <- cand |> distinct(dispnum) |> left_join(link, by = "dispnum") |>
    mutate(across(c(n_pair_decided, n_link, n_map, n_map_modal), ~ coalesce(.x, 0L)), link_ok = coalesce(link_ok, FALSE), map_ok = coalesce(map_ok, FALSE),
           dispute_dec = coalesce(dispute_dec, "no_pair_match"))
  out <- ev |> left_join(dec |> select(all_of(key), n_war_match, pair_dec, pair_war = war_link), by = key) |>
    left_join(all_war |> select(dispnum, dispute_dec, disp_war = war_link), by = "dispnum") |>
    mutate(n_war_match = coalesce(n_war_match, 0L),
           dir_rule = case_when(hostlev_disp != 5L ~ "non_war_dispute",
                                pair_dec %in% c("flip", "keep") ~ "pair_match",
                                pair_dec == "ambiguous" ~ "pair_ambiguous",
                                dispute_dec %in% c("flip", "keep") ~ "dispute_mapping",
                                TRUE ~ paste0("dispute_", dispute_dec)),
           dir_source = case_when(dir_rule == "pair_match" ~ "cow_war_pair", dir_rule == "dispute_mapping" ~ "cow_war_dispute", TRUE ~ "sideA"),
           dir_flipped = (dir_rule == "pair_match" & pair_dec == "flip") | (dir_rule == "dispute_mapping" & dispute_dec == "flip"),
           war_link = case_when(dir_rule == "pair_match" ~ pair_war, dir_rule == "dispute_mapping" ~ disp_war, TRUE ~ NA_integer_)) |>
    select(-pair_dec, -dispute_dec, -pair_war, -disp_war)
  stopifnot(nrow(out) == nrow(ev), !anyNA(out$dir_flipped), !anyNA(out$dir_source), all(out$dir_source %in% c("cow_war_pair", "cow_war_dispute", "sideA")),
            !any(out$dir_flipped & out$dir_source == "sideA"), all(!is.na(out$war_link) == (out$dir_source != "sideA")),
            !any(out$dir_source != "sideA" & out$hostlev_disp != 5L))
  attr(out, "pair_matches") <- m
  attr(out, "dispute_links") <- all_war
  out
}
## apply_direction(pair_ev): reverse the flipped rows of the pair table (initiator <-> target, and every per-participant
## column with it: the MID side flags, originator status, fatality and hostility), so that the originator status moves with
## the state. After this sidea1 = 1 on every unflipped row and 0 on every flipped row.
apply_direction <- function(pair_ev) {
  stopifnot(all(c("initiator", "target", "sidea1", "sidea2", "orig1", "orig2", "fatality1", "fatality2", "hostlev1", "hostlev2", "dir_flipped") %in% names(pair_ev)),
            all(pair_ev$sidea1 == 1L), all(pair_ev$sidea2 == 0L))
  f <- pair_ev$dir_flipped
  out <- pair_ev
  swap <- function(a, b) { tmp <- out[[a]][f]; out[[a]][f] <<- out[[b]][f]; out[[b]][f] <<- tmp }
  swap("initiator", "target"); swap("sidea1", "sidea2"); swap("orig1", "orig2"); swap("fatality1", "fatality2"); swap("hostlev1", "hostlev2")
  stopifnot(all(out$initiator[f] == pair_ev$target[f]), all(out$target[f] == pair_ev$initiator[f]), all(out$sidea1[f] == 0L), all(out$sidea2[f] == 1L),
            all(out$orig1[f] == pair_ev$orig2[f]), all(out$initiator[!f] == pair_ev$initiator[!f]), all(out$initiator != out$target),
            all((out$sidea1 == 0L) == out$dir_flipped))
  assert_no_dup_keys(out, c("dispnum", "initiator", "target"), "events after the direction rule")
  out
}
## COW inter-state war names for the build log (hand-entered from the COW Inter-State War Data v4.0 list; logging only,
## nothing is computed from them; a war number absent here is logged without a name)
COW_WAR_NAMES <- c(`1` = "Franco-Spanish 1823", `4` = "First Russo-Turkish 1828", `7` = "Mexican-American 1846", `10` = "Austro-Sardinian 1848",
                   `13` = "First Schleswig-Holstein 1848", `16` = "Roman Republic 1849", `19` = "La Plata 1851", `22` = "Crimean 1853",
                   `25` = "Anglo-Persian 1856", `28` = "Italian Unification 1859", `31` = "First Spanish-Moroccan 1859", `34` = "Italo-Roman 1860",
                   `37` = "Italo-Sicilian 1860", `40` = "Franco-Mexican 1862", `43` = "Ecuadorian-Colombian 1863", `46` = "Second Schleswig-Holstein 1864",
                   `49` = "Lopez (Paraguayan) 1864", `52` = "Naval War 1865", `55` = "Seven Weeks 1866", `58` = "Franco-Prussian 1870",
                   `60` = "First Central American 1876", `61` = "Second Russo-Turkish 1877", `64` = "Pacific 1879", `65` = "Anglo-Egyptian 1882",
                   `67` = "Sino-French 1884", `70` = "Second Central American 1885", `73` = "First Sino-Japanese 1894", `76` = "Greco-Turkish 1897",
                   `79` = "Spanish-American 1898", `82` = "Boxer Rebellion 1900", `83` = "Sino-Russian 1900", `85` = "Russo-Japanese 1904",
                   `88` = "Third Central American 1906", `91` = "Fourth Central American 1907", `94` = "Second Spanish-Moroccan 1909", `97` = "Italian-Turkish 1911",
                   `100` = "First Balkan 1912", `103` = "Second Balkan 1913", `106` = "World War I 1914", `107` = "Estonian Liberation 1918",
                   `108` = "Latvian Liberation 1918", `109` = "Russo-Polish 1919", `112` = "Hungarian Adversaries 1919", `115` = "Second Greco-Turkish 1919",
                   `116` = "Franco-Turkish 1919", `117` = "Lithuanian-Polish 1920", `118` = "Manchurian (Sino-Soviet) 1929", `121` = "Second Sino-Japanese 1931",
                   `124` = "Chaco 1932", `125` = "Saudi-Yemeni 1934", `127` = "Conquest of Ethiopia 1935", `130` = "Third Sino-Japanese 1937",
                   `133` = "Changkufeng 1938", `136` = "Nomonhan 1939", `139` = "World War II 1939", `142` = "Russo-Finnish 1939",
                   `145` = "Franco-Thai 1940", `147` = "First Kashmir 1947", `148` = "Arab-Israeli 1948", `151` = "Korean 1950",
                   `153` = "Off-shore Islands 1954", `155` = "Sinai 1956", `156` = "Soviet Invasion of Hungary 1956", `158` = "Ifni 1957",
                   `159` = "Taiwan Straits 1958", `160` = "Assam 1962", `163` = "Vietnam War Phase 2 1965", `166` = "Second Kashmir 1965",
                   `169` = "Six Day 1967", `170` = "Second Laotian Phase 2 1968", `172` = "War of Attrition 1969", `175` = "Football 1969",
                   `176` = "Communist Coalition 1970", `178` = "Bangladesh 1971", `181` = "Yom Kippur 1973", `184` = "Turco-Cypriot 1974",
                   `186` = "War over Angola 1975", `187` = "Second Ogaden Phase 2 1977", `189` = "Vietnamese-Cambodian 1977", `190` = "Ugandan-Tanzanian 1978",
                   `193` = "Sino-Vietnamese Punitive 1979", `199` = "Iran-Iraq 1980", `202` = "Falklands 1982", `205` = "War over Lebanon 1982",
                   `207` = "War over the Aouzou Strip 1986", `208` = "Sino-Vietnamese Border 1987", `211` = "Gulf War 1990", `215` = "Bosnian Independence 1992",
                   `216` = "Azeri-Armenian 1993", `217` = "Cenepa Valley 1995", `219` = "Badme Border 1998", `221` = "War for Kosovo 1999",
                   `223` = "Kargil 1999", `225` = "Invasion of Afghanistan 2001", `227` = "Invasion of Iraq 2003")
cow_war_name <- function(warnum) ifelse(is.na(warnum), NA_character_, coalesce(COW_WAR_NAMES[as.character(warnum)], "(name not in COW_WAR_NAMES)"))

## ---- WWII disputes (PI-approved appendix option --drop_ww2, 2026-10-24) -----------------------------------------
## COW-MID 5.0 codes the Second World War as dispute 258 (1939-03 to 1945-09, 26 states). Several other war-level
## (hostlev 5) disputes of 1939-1945 are the run-ups to a state's ENTRY into that war (Germany vs Norway, Italy vs
## Greece, Germany vs the USSR, the United States vs Japan and Germany, the USSR vs Japan in 1945, ...). COW closes
## each of them with highest action 21 ("join interstate war") and codes the same dyad as a WWII war dyad-year in
## the COW inter-state war data (warnum 139). The rule applied here, derived from the data and logged by the build:
##   WWII = 258, plus every hostlev-5 dispute that (a) starts between 1939-09 and 1945-09, (b) has dispute-level
##   highest action 21 (join interstate war) and (c) has at least one of its dyad-years in COW inter-state war 139.
## Excluded on purpose, with the reason the build logs: 157 (Third Sino-Japanese War, COW war 130, began 1937-07,
## before the window), 179 (Russo-Finnish Winter War, COW war 142), 183 (Nomonhan, COW war 136) and 613
## (Franco-Thai War, COW war 145); these begin a war of their own (highest action 20) and share no dyad-year with
## war 139. The list is asserted against WW2_DISPNUMS so a data update changes it loudly, not silently.
WW2_DISPNUMS <- c(235L, 258L, 339L, 414L, 418L, 518L, 2725L, 3701L, 3705L, 3717L, 3724L, 3813L, 3820L, 3822L, 3826L, 3827L)
WW2_HIACT_JOIN_WAR <- 21L
ww2_disputes <- function(disp, dd, war_inter) {
  stopifnot(all(c("dispnum", "styear", "stmon", "hiact", "hostlev") %in% names(disp)),
            all(c("warnum", "ccode1", "ccode2", "year") %in% names(war_inter)))
  ww_dyads <- war_inter |> filter(warnum == 139) |> transmute(ccode1 = as.integer(ccode1), ccode2 = as.integer(ccode2), year = as.integer(year)) |> distinct()
  stopifnot(nrow(ww_dyads) > 0)
  yrs <- dd |> group_by(dispnum) |> summarise(y0 = min(year), y1 = max(year), n_states = n_distinct(c(ccode1, ccode2)), .groups = "drop")
  in_ww <- dd |> distinct(dispnum, ccode1, ccode2, year) |> semi_join(ww_dyads, by = c("ccode1", "ccode2", "year")) |> count(dispnum, name = "n_dyadyears_in_war139")
  cand <- disp |> filter(hostlev == 5L) |> inner_join(yrs, by = "dispnum") |> filter(y0 <= 1945L, y1 >= 1939L) |>
    left_join(in_ww, by = "dispnum") |> mutate(n_dyadyears_in_war139 = coalesce(n_dyadyears_in_war139, 0L),
                                              start = styear * 100L + stmon,
                                              in_window = start >= 193909L & start <= 194509L,
                                              joins_war = hiact == WW2_HIACT_JOIN_WAR,
                                              is_ww2 = dispnum == 258L | (in_window & joins_war & n_dyadyears_in_war139 > 0L),
                                              reason = case_when(dispnum == 258L ~ "COW-MID dispute of the Second World War",
                                                                 is_ww2 ~ "starts 1939-09..1945-09, highest action 21 (join interstate war), dyad-years in COW war 139",
                                                                 !in_window ~ "starts before 1939-09 (a war of its own)",
                                                                 !joins_war ~ "highest action 20 (begins a war of its own), not a WWII entry",
                                                                 TRUE ~ "no dyad-year in COW inter-state war 139")) |>
    arrange(dispnum)
  got <- sort(cand$dispnum[cand$is_ww2])
  if (!identical(got, sort(WW2_DISPNUMS)))
    stop(sprintf("WWII dispute list derived from the data (%s) differs from WW2_DISPNUMS (%s); review and update helpers_data.R",
                 paste(got, collapse = ","), paste(sort(WW2_DISPNUMS), collapse = ",")))
  cand |> select(dispnum, styear, stmon, y0, y1, n_states, hiact, hostlev, n_dyadyears_in_war139, is_ww2, reason)
}

## ---- coalition context (dual typology stage a, PI-approved spec 2026-10-25; the `context` switch) ----------------------------
## Is a stage-2 event a BILATERAL or a COALITION decision? Context = the dispute has more than one state on either side at the
## onset. The unit is the D8 weight group (initiator, dispute, onset year): the flag is evaluated ONCE per group at the GROUP ONSET,
## the earliest pair onset among the initiator's pairs of that dispute-year (the date the initiator first faced an opponent in that
## dispute-year), so it is constant within the group by construction (asserted) and a group is kept or dropped whole by the
## `context` switch. Evaluating at every pair's own onset instead gives different flags inside 72 (inclusive) / 96 (strict) of the
## 2874 groups of the COW data (Germany in dispute 258, 1939: the Poland pair starts 1939-03-25 with no third state, the UK and
## France pairs 1939-03-31 with Poland present); the build logs that per-event heterogeneity as a diagnostic, it does not use it.
## Dates have three resolutions, res 1 = day, 2 = month, 3 = year, and two dates are compared at the COARSER of the two (date_cmp();
## a day-dated entry against a month-dated onset in the same month is a tie), so the tie conventions below act at whatever resolution
## the pair of dates has. Participant dates come from gml_part (D5) where the dispnum and the participant match and the GML years
## agree with the COW years (participant_spell_dates(), the participant-level analogue of pair_dispute_intervals()); otherwise the
## dispute-level bounds (dispute start month for a first-year originator, else the year). Membership is COW's (cow_mid_dirdisps).
## Two tie conventions (both stored, the `coal_ties` switch selects): INCLUSIVE counts a state entering on the same date as the
## onset as present; STRICT counts only strictly earlier entrants. A state whose spell ended before the onset for sure is absent.
## coal_incl   : at least one state besides the initiator and ONE onset-defining target is present (entry <= onset, not exited)
## coal_strict : the same with entry < onset; an onset-defining target that itself entered strictly before the onset (the initiator
##               joined an existing pair) is the pair and does not count
## date_cmp(a, b): a and b are lists / data.frames with res, y, mkey (ym_index, NA at year resolution), dkey (as.integer(Date), NA
## below day resolution); returns cmp in {-1, 0, 1} (a before / same / after b at the coarser resolution) and the resolution used.
date_keys <- function(res, y, m, d) {
  stopifnot(all(res %in% 1:3), !anyNA(y), all(!is.na(m[res <= 2L])), all(!is.na(d[res == 1L])))
  mkey <- ifelse(res <= 2L, ym_index(y, ifelse(res <= 2L, m, 1L)), NA_integer_)
  dstr <- sprintf("%d-%02d-%02d", y, ifelse(res <= 2L, m, 1L), ifelse(res == 1L, d, 1L))
  dkey <- ifelse(res == 1L, as.integer(as.Date(dstr, optional = TRUE)), NA_integer_)
  if (any(res == 1L & is.na(dkey))) stop("date_keys: invalid calendar day in ", paste(dstr[res == 1L & is.na(dkey)], collapse = ", "))
  data.frame(res = as.integer(res), y = as.integer(y), mkey = as.integer(mkey), dkey = as.integer(dkey))
}
date_cmp <- function(a, b) {
  n <- length(a$res); stopifnot(length(b$res) == n)
  r <- pmax(a$res, b$res)
  ka <- ifelse(r == 1L, a$dkey, ifelse(r == 2L, a$mkey, a$y)); kb <- ifelse(r == 1L, b$dkey, ifelse(r == 2L, b$mkey, b$y))
  stopifnot(!anyNA(ka), !anyNA(kb))
  list(cmp = as.integer(sign(ka - kb)), res = as.integer(r))
}
## date_later(a, b): the later of two dates; a tie at the coarser resolution returns that date at the coarser resolution
date_later <- function(a, b) {
  cc <- date_cmp(a, b)
  tie <- cc$cmp == 0L
  data.frame(res = ifelse(cc$cmp > 0L, a$res, ifelse(cc$cmp < 0L, b$res, cc$res)),
             y = ifelse(cc$cmp >= 0L, a$y, b$y),
             mkey = ifelse(cc$cmp > 0L, a$mkey, ifelse(cc$cmp < 0L, b$mkey, ifelse(cc$res <= 2L, a$mkey, NA_integer_))),
             dkey = ifelse(cc$cmp > 0L, a$dkey, ifelse(cc$cmp < 0L, b$dkey, ifelse(cc$res == 1L, a$dkey, NA_integer_))))
}
## date_earliest_by(d, g): the earliest date within each group g, with ties resolved at the coarsest tied resolution (the group's
## minimum year; if any date of that year is year-dated the result is that year, else the minimum month; if any date of that month is
## month-dated the result is that month, else the minimum day). One row per group, in the order of sort(unique(g)).
date_earliest_by <- function(d, g) {
  stopifnot(length(g) == length(d$res))
  tb <- tibble(g = g, res = d$res, y = d$y, mkey = d$mkey, dkey = d$dkey) |> group_by(g) |> filter(y == min(y)) |>
    mutate(.yr = any(res == 3L)) |> filter(.yr | mkey == min(mkey)) |> mutate(.mo = .yr | any(res == 2L)) |>
    filter(.mo | dkey == min(dkey)) |>
    summarise(res = if (first(.yr)) 3L else if (first(.mo)) 2L else 1L, y = first(y),
              mkey = if (first(.yr)) NA_integer_ else first(mkey), dkey = if (first(.yr) || first(.mo)) NA_integer_ else first(dkey), .groups = "drop") |>
    arrange(g)
  stopifnot(nrow(tb) == length(unique(g)), !anyNA(tb$y), all(!is.na(tb$mkey[tb$res <= 2L])), all(!is.na(tb$dkey[tb$res == 1L])))
  tb
}
## participant_spell_dates(dd, db, part): one row per (dispnum, ccode, spell), a spell being a run of consecutive COW years of the
## participant in the dispute, with the entry date (entry_res / entry_y / entry_m / entry_d and the keys) and the exit date. GML
## (part = gml_part) is used for a spell when the participant has at least one gml_part spell overlapping the COW spell's years: the entry
## is the EARLIEST GML start (several overlapping GML spells, e.g. France in dispute 258 switching sides in 1940, are one COW spell) when
## its year does not exceed the COW spell's first year (a GML start before the first COW year is clamped to January of it, at month
## resolution, as in D5); the exit is the LATEST GML end when its year is not before the last COW year (a later end is clamped to December).
## Day resolution needs a known day (stday / endday > 0). Fallback entry: the dispute start month for an originator whose spell begins
## in the dispute's first year (start month known), else the spell's first year; fallback exit: the dispute's end month when the spell
## ends in the dispute's last year and that end is an exact month, else the spell's last year. attr(, "membership") lists the
## (dispnum, ccode) where COW and GML disagree on who took part (part = NULL: no GML dates at all).
participant_spell_dates <- function(dd, db, part = NULL) {
  stopifnot(all(c("dispnum", "ccode1", "year", "orig1") %in% names(dd)), all(c("dispnum", "y0", "y1", "start_month_known", "stmon", "end_lo", "end_hi") %in% names(db)))
  sp <- dd |> distinct(dispnum, ccode = ccode1, year, orig1) |> arrange(dispnum, ccode, year) |> group_by(dispnum, ccode) |>
    mutate(spell = cumsum(c(1L, diff(year) != 1L))) |> group_by(dispnum, ccode, spell) |>
    summarise(spell_first = min(year), spell_last = max(year), orig = as.integer(orig1[1]), n_orig = n_distinct(orig1), .groups = "drop") |>
    inner_join(db |> select(dispnum, y0, y1, start_month_known, stmon, end_lo, end_hi), by = "dispnum")
  if (any(sp$n_orig != 1L)) stop("participant_spell_dates: originator status differs across a participant's dispute-years")
  g_entry <- g_exit <- NULL
  if (!is.null(part)) {
    stopifnot(all(c("dispnum", "ccode", "styear", "stmon", "stday", "endyear", "endmon", "endday") %in% names(part)))
    ps <- part |> transmute(dispnum = as.integer(dispnum), ccode = as.integer(ccode), g_sy = as.integer(styear), g_sm = as.integer(stmon), g_sd = as.integer(stday),
                            g_ey = as.integer(endyear), g_em = as.integer(endmon), g_ed = as.integer(endday))
    stopifnot(!anyNA(ps), all(ps$g_sm %in% 1:12), all(ps$g_em %in% 1:12))
    ov <- sp |> mutate(.id = row_number()) |> select(.id, dispnum, ccode, spell_first, spell_last) |>
      inner_join(ps, by = c("dispnum", "ccode"), relationship = "many-to-many") |> filter(g_ey >= spell_first, g_sy <= spell_last) |>
      group_by(.id) |> arrange(g_sy, g_sm, pmax(g_sd, 1L), .by_group = TRUE) |>
      summarise(n_gml = n(), g_sy = first(g_sy), g_sm = first(g_sm), g_sd = first(g_sd),
                g_ey = g_ey[which.max(g_ey * 400L + g_em * 32L + pmax(g_ed, 1L))], g_em = g_em[which.max(g_ey * 400L + g_em * 32L + pmax(g_ed, 1L))],
                g_ed = g_ed[which.max(g_ey * 400L + g_em * 32L + pmax(g_ed, 1L))], .groups = "drop")
    sp <- sp |> mutate(.id = row_number()) |> left_join(ov, by = ".id") |> mutate(n_gml = coalesce(n_gml, 0L), in_gml_disp = dispnum %in% ps$dispnum)
    cow_part <- sp |> distinct(dispnum, ccode); gml_part_tbl <- ps |> distinct(dispnum, ccode) |> filter(dispnum %in% cow_part$dispnum)
    membership <- bind_rows(cow_part |> filter(dispnum %in% ps$dispnum) |> anti_join(gml_part_tbl, by = c("dispnum", "ccode")) |> mutate(issue = "COW participant absent from gml_part"),
                            gml_part_tbl |> anti_join(cow_part, by = c("dispnum", "ccode")) |> mutate(issue = "gml_part participant absent from COW"))
  } else {
    sp <- sp |> mutate(.id = row_number(), n_gml = 0L, g_sy = NA_integer_, g_sm = NA_integer_, g_sd = NA_integer_, g_ey = NA_integer_, g_em = NA_integer_, g_ed = NA_integer_, in_gml_disp = FALSE)
    membership <- tibble(dispnum = integer(0), ccode = integer(0), issue = character(0))
  }
  sp <- sp |> mutate(
    entry_src = case_when(n_gml >= 1L & g_sy == spell_first ~ "participant", n_gml >= 1L & g_sy < spell_first ~ "participant (entry clamped to January)",
                          !in_gml_disp ~ "dispute (not in gml_part)", n_gml == 0L ~ "dispute (participant has no overlapping gml_part spell)",
                          TRUE ~ "dispute (GML entry after the first COW year)"),
    entry_res = case_when(entry_src == "participant" & g_sd > 0L ~ 1L, startsWith(entry_src, "participant") ~ 2L,
                          spell_first == y0 & orig == 1L & start_month_known ~ 2L, TRUE ~ 3L),
    entry_y = spell_first,
    entry_m = case_when(entry_src == "participant" ~ g_sm, startsWith(entry_src, "participant") ~ 1L, entry_res == 2L ~ stmon, TRUE ~ NA_integer_),
    entry_d = ifelse(entry_res == 1L, g_sd, NA_integer_),
    exit_src = case_when(n_gml >= 1L & g_ey == spell_last ~ "participant", n_gml >= 1L & g_ey > spell_last ~ "participant (exit clamped to December)",
                         !in_gml_disp ~ "dispute (not in gml_part)", n_gml == 0L ~ "dispute (participant has no overlapping gml_part spell)",
                         TRUE ~ "dispute (GML exit before the last COW year)"),
    exit_res = case_when(exit_src == "participant" & g_ed > 0L ~ 1L, startsWith(exit_src, "participant") ~ 2L,
                         spell_last == y1 & end_lo == end_hi ~ 2L, TRUE ~ 3L),
    exit_y = spell_last,
    exit_m = case_when(exit_src == "participant" ~ g_em, startsWith(exit_src, "participant") ~ 12L, exit_res == 2L ~ end_lo %% 12L + 1L, TRUE ~ NA_integer_),
    exit_d = ifelse(exit_res == 1L, g_ed, NA_integer_))
  ek <- date_keys(sp$entry_res, sp$entry_y, sp$entry_m, sp$entry_d); xk <- date_keys(sp$exit_res, sp$exit_y, sp$exit_m, sp$exit_d)
  sp <- sp |> mutate(entry_mkey = ek$mkey, entry_dkey = ek$dkey, exit_mkey = xk$mkey, exit_dkey = xk$dkey)
  ## the entry never lies after the exit at the coarser resolution
  cx <- date_cmp(list(res = sp$entry_res, y = sp$entry_y, mkey = sp$entry_mkey, dkey = sp$entry_dkey), list(res = sp$exit_res, y = sp$exit_y, mkey = sp$exit_mkey, dkey = sp$exit_dkey))
  bad <- cx$cmp > 0L
  if (any(bad)) {
    ## a GML spell ending before it started within the COW years cannot be used: year resolution for both
    sp$entry_res[bad] <- 3L; sp$entry_m[bad] <- NA_integer_; sp$entry_d[bad] <- NA_integer_; sp$entry_mkey[bad] <- NA_integer_; sp$entry_dkey[bad] <- NA_integer_
    sp$exit_res[bad] <- 3L; sp$exit_m[bad] <- NA_integer_; sp$exit_d[bad] <- NA_integer_; sp$exit_mkey[bad] <- NA_integer_; sp$exit_dkey[bad] <- NA_integer_
    sp$entry_src[bad] <- "dispute (GML entry after exit)"; sp$exit_src[bad] <- "dispute (GML entry after exit)"
  }
  out <- sp |> select(dispnum, ccode, spell, spell_first, spell_last, orig, entry_src, entry_res, entry_y, entry_m, entry_d, entry_mkey, entry_dkey,
                      exit_src, exit_res, exit_y, exit_m, exit_d, exit_mkey, exit_dkey)
  stopifnot(all(out$entry_y == out$spell_first), all(out$exit_y == out$spell_last), !anyNA(out$entry_res), !anyNA(out$exit_res))
  assert_no_dup_keys(out, c("dispnum", "ccode", "spell"), "participant spell dates")
  attr(out, "membership") <- membership
  attr(out, "n_entry_after_exit_fixed") <- sum(bad)
  out
}
## context_present(onsets, psd): for every onset row (dispnum, initiator, onset date keys, .g) the participants p != initiator whose
## spell is present at the onset. Returns one row per (onset row, participant spell) with cmp_entry, cmp_exit, the comparison
## resolution and the presence flags under both conventions.
context_present <- function(onsets, psd) {
  stopifnot(all(c(".g", "dispnum", "initiator", "res", "y", "mkey", "dkey") %in% names(onsets)))
  cand <- onsets |> select(.g, dispnum, initiator, o_res = res, o_y = y, o_mkey = mkey, o_dkey = dkey) |>
    inner_join(psd |> select(dispnum, ccode, spell, entry_res, entry_y, entry_mkey, entry_dkey, exit_res, exit_y, exit_mkey, exit_dkey), by = "dispnum", relationship = "many-to-many") |>
    filter(ccode != initiator)
  ce <- date_cmp(list(res = cand$entry_res, y = cand$entry_y, mkey = cand$entry_mkey, dkey = cand$entry_dkey), list(res = cand$o_res, y = cand$o_y, mkey = cand$o_mkey, dkey = cand$o_dkey))
  cx <- date_cmp(list(res = cand$exit_res, y = cand$exit_y, mkey = cand$exit_mkey, dkey = cand$exit_dkey), list(res = cand$o_res, y = cand$o_y, mkey = cand$o_mkey, dkey = cand$o_dkey))
  cand |> mutate(cmp_entry = ce$cmp, cmp_res = ce$res, cmp_exit = cx$cmp, exited = cmp_exit < 0L,
                 present_incl = cmp_entry <= 0L & !exited, present_strict = cmp_entry < 0L & !exited)
}
## context_flags(onsets, targets, psd): the flags of every onset row. onsets: .g, dispnum, initiator, res, y, mkey, dkey; targets: the
## onset-defining targets (.g, ccode), which ARE the pair and are present by construction whatever their GML exit says (two COW pairs of
## the data have the target's GML exit a few days before the initiator's GML entry in the same month; counted in n_target_exit_before_onset);
## their entry must not be after the onset (asserted). A participant with several spells is present when ANY spell is. Per onset row:
## n_incl / n_strict = the participants other than the initiator present under each convention, target_strict = an onset-defining target
## entered strictly before the onset (the initiator joined an existing pair, so that target does not make a coalition), coal_incl =
## n_incl >= 2, coal_strict = n_strict - target_strict >= 1, coal_res = the resolution of the entry comparisons (day / month / year / mixed).
context_flags <- function(onsets, targets, psd) {
  stopifnot(all(c(".g", "ccode") %in% names(targets)), !anyDuplicated(onsets$.g))
  pres <- context_present(onsets, psd)
  jt <- targets |> inner_join(pres |> group_by(.g, ccode) |> summarise(entry_ok = any(cmp_entry <= 0L), present_incl = any(present_incl), entered_before = any(cmp_entry < 0L), .groups = "drop"),
                              by = c(".g", "ccode"))
  if (nrow(jt) != nrow(targets) || !all(jt$entry_ok)) stop(sprintf("context_flags: %d onset-defining targets entered after their own onset or have no spell", nrow(targets) - sum(jt$entry_ok)))
  pres <- pres |> left_join(jt |> transmute(.g, ccode, .force = TRUE), by = c(".g", "ccode")) |>
    mutate(.force = coalesce(.force, FALSE), present_incl = present_incl | .force, present_strict = present_strict | (.force & cmp_entry < 0L))
  agg <- pres |> group_by(.g) |>
    summarise(n_incl = n_distinct(ccode[present_incl]), n_strict = n_distinct(ccode[present_strict]), res_set = paste(sort(unique(cmp_res)), collapse = ""), .groups = "drop") |>
    left_join(jt |> group_by(.g) |> summarise(target_strict = any(entered_before), n_target_exit_before_onset = sum(!present_incl), .groups = "drop"), by = ".g")
  out <- onsets |> select(.g) |> left_join(agg, by = ".g") |>
    mutate(n_incl = coalesce(n_incl, 0L), n_strict = coalesce(n_strict, 0L), target_strict = coalesce(target_strict, FALSE), n_target_exit_before_onset = coalesce(n_target_exit_before_onset, 0L),
           coal_incl = n_incl >= 2L, coal_strict = (n_strict - as.integer(target_strict)) >= 1L,
           coal_res = case_when(res_set == "1" ~ "day", res_set == "2" ~ "month", res_set == "3" ~ "year", TRUE ~ "mixed"))
  stopifnot(nrow(out) == nrow(onsets), all(out$n_incl >= 1L), !anyNA(out$coal_incl), !anyNA(out$coal_strict), all(out$coal_incl >= out$coal_strict))
  out |> select(.g, n_incl, n_strict, target_strict, n_target_exit_before_onset, coal_incl, coal_strict, coal_res)
}
## coalition_context(pairs, psd): pairs = the PRE-collapse non-escalation event rows (dispnum, initiator, target, year, ev_start_lo,
## ev_start_hi, the pair's month-level entry interval from flag_escalations_month()); psd = participant_spell_dates(). Returns pairs with
##   coal_incl, coal_strict      : the group's flag (constant within (dispnum, initiator, year), asserted)
##   coal_onset_ambiguous        : the pair onset derived from the two participants' entries (the later of the two, within their spells
##                                 containing the onset year) lies in an EARLIER year than the pair's COW onset year (both states were in the
##                                 dispute before they were opposed: a side switch), so the onset is ambiguous across a year boundary and the
##                                 COW onset year is used at year resolution; ev_start_lo / ev_start_hi (the escalation rule's month interval)
##                                 enter only the logged diagnostic onset_vs_pdi (attr "events")
##   coal_res                    : the resolution of the group's entry comparisons: "day" / "month" / "year" / "mixed"
## attr(, "groups") is the group table (onset, counts, flags, resolution, n_events, and the per-event heterogeneity diagnostic
## n_event_flags_incl / n_event_flags_strict = the number of distinct per-event flags when every pair is evaluated at its own onset).
coalition_context <- function(pairs, psd) {
  stopifnot(all(c("dispnum", "initiator", "target", "year", "ev_start_lo", "ev_start_hi") %in% names(pairs)))
  assert_no_dup_keys(pairs, c("dispnum", "initiator", "target"), "coalition_context: pairs")
  n_in <- nrow(pairs)
  pr <- pairs |> mutate(.e = row_number())
  spell_at <- function(code_col, pre) {
    pr |> select(.e, dispnum, ccode = all_of(code_col), year) |>
      inner_join(psd |> select(dispnum, ccode, spell_first, spell_last, entry_res, entry_y, entry_mkey, entry_dkey), by = c("dispnum", "ccode"), relationship = "many-to-many") |>
      filter(year >= spell_first, year <= spell_last) |> select(.e, entry_res, entry_y, entry_mkey, entry_dkey) |> rename_with(~ paste0(pre, .x), -.e)
  }
  si <- spell_at("initiator", "i_"); sj <- spell_at("target", "j_")
  if (nrow(si) != n_in || nrow(sj) != n_in) stop("coalition_context: a pair's initiator or target has no participant spell containing the onset year (or several)")
  pr <- pr |> inner_join(si, by = ".e") |> inner_join(sj, by = ".e")
  on <- date_later(list(res = pr$i_entry_res, y = pr$i_entry_y, mkey = pr$i_entry_mkey, dkey = pr$i_entry_dkey),
                   list(res = pr$j_entry_res, y = pr$j_entry_y, mkey = pr$j_entry_mkey, dkey = pr$j_entry_dkey))
  ## the participants' entries are dated within the spells containing the onset year, so the later entry is in year t or EARLIER; earlier
  ## means both were in the dispute before they were opposed (a side switch, e.g. France and the USSR in dispute 258) and the pair's
  ## onset date is ambiguous: the COW onset year at year resolution is used
  stopifnot(all(on$y <= pr$year))
  amb <- on$y != pr$year
  on$res[amb] <- 3L; on$y[amb] <- pr$year[amb]; on$mkey[amb] <- NA_integer_; on$dkey[amb] <- NA_integer_
  stopifnot(all(on$y == pr$year), all(!is.na(on$mkey[on$res <= 2L])), all(!is.na(on$dkey[on$res == 1L])))
  ## diagnostic only: does the participant-derived onset month agree with the pair's month interval of flag_escalations_month()?
  onset_vs_pdi <- case_when(amb ~ "ambiguous (year fallback)", on$res == 3L ~ "year resolution",
                            on$mkey >= pr$ev_start_lo & on$mkey <= pr$ev_start_hi ~ "inside the pair interval", TRUE ~ "outside the pair interval")
  pr <- pr |> mutate(o_res = on$res, o_y = on$y, o_mkey = on$mkey, o_dkey = on$dkey, coal_onset_ambiguous = amb, onset_vs_pdi = onset_vs_pdi)
  ## group onsets: the earliest pair onset of (dispnum, initiator, year), ties at the coarsest tied resolution
  pr <- pr |> mutate(.g = as.integer(factor(paste(dispnum, initiator, year))))
  gon <- date_earliest_by(list(res = pr$o_res, y = pr$o_y, mkey = pr$o_mkey, dkey = pr$o_dkey), pr$.g) |> rename(.g = g)
  grp <- pr |> group_by(.g, dispnum, initiator, year) |> summarise(n_events = n(), .groups = "drop") |> inner_join(gon, by = ".g")
  stopifnot(nrow(grp) == n_distinct(pr$.g), all(grp$y == grp$year))
  ## the pair onsets that define the group onset (ties at the coarser resolution)
  pg <- pr |> inner_join(grp |> select(.g, g_res = res, g_y = y, g_mkey = mkey, g_dkey = dkey), by = ".g")
  cg <- date_cmp(list(res = pg$o_res, y = pg$o_y, mkey = pg$o_mkey, dkey = pg$o_dkey), list(res = pg$g_res, y = pg$g_y, mkey = pg$g_mkey, dkey = pg$g_dkey))
  stopifnot(all(cg$cmp >= 0L))                                     # no pair starts before its group
  pr$onset_defining <- cg$cmp == 0L
  ## presence of every other participant at the group onset; the onset-defining targets are the pair (context_flags())
  fl <- context_flags(grp, pr |> filter(onset_defining) |> select(.g, ccode = target), psd)
  grp <- grp |> inner_join(fl, by = ".g")
  stopifnot(nrow(grp) == n_distinct(pr$.g), all(grp$n_incl >= 1L), all(grp$n_strict <= grp$n_incl), !anyNA(grp$coal_incl), !anyNA(grp$coal_strict), all(grp$coal_incl >= grp$coal_strict))
  ## diagnostic: the flag at every pair's OWN onset (not used); how many groups would carry different flags
  fl_e <- context_flags(pr |> transmute(.g = .e, dispnum, initiator, res = o_res, y = o_y, mkey = o_mkey, dkey = o_dkey), pr |> select(.g = .e, ccode = target), psd)
  pr <- pr |> left_join(fl_e |> select(.e = .g, ev_incl = coal_incl, ev_strict = coal_strict), by = ".e")
  stopifnot(!anyNA(pr$ev_incl), !anyNA(pr$ev_strict))
  het <- pr |> group_by(.g) |> summarise(n_event_flags_incl = n_distinct(ev_incl), n_event_flags_strict = n_distinct(ev_strict), .groups = "drop")
  grp <- grp |> left_join(het, by = ".g")
  out <- pr |> inner_join(grp |> select(.g, coal_incl, coal_strict, coal_res), by = ".g")
  stopifnot(nrow(out) == n_in, !anyNA(out$coal_incl), !anyNA(out$coal_strict))
  ## constant within every (dispnum, initiator, year) group: by construction, asserted on the output
  chk <- out |> group_by(dispnum, initiator, year) |> summarise(a = n_distinct(coal_incl), b = n_distinct(coal_strict), .groups = "drop")
  if (any(chk$a != 1L) || any(chk$b != 1L)) stop("coalition_context: the coalition flag is not constant within an (initiator, dispute, year) group")
  res <- pairs |> left_join(out |> select(dispnum, initiator, target, coal_incl, coal_strict, coal_onset_ambiguous, coal_res), by = c("dispnum", "initiator", "target"))
  stopifnot(nrow(res) == n_in, !anyNA(res$coal_incl))
  attr(res, "groups") <- grp |> select(dispnum, initiator, year, n_events, onset_res = res, onset_y = y, onset_mkey = mkey, onset_dkey = dkey, n_incl, n_strict, target_strict,
                                       n_target_exit_before_onset, coal_incl, coal_strict, coal_res, n_event_flags_incl, n_event_flags_strict)
  attr(res, "events") <- out |> select(dispnum, initiator, target, year, o_res, o_mkey, o_dkey, onset_defining, onset_vs_pdi, ev_incl, ev_strict)
  res
}

## ---- gap scale (PI decision 2026-10-23): SD of the within-year gap v_i - v_j over all ordered pairs i != j ----
## Exact, without forming the pairs. Per year with n values and SS = sum((v - mean(v))^2), sum over ordered pairs of
## (v_i - v_j)^2 = 2 n SS and there are n (n - 1) pairs; the pairs' mean is exactly 0. Years with one state add nothing.
gap_scale <- function(v, yr) {
  stopifnot(length(v) == length(yr), !anyNA(v), all(is.finite(v)))
  ss <- tapply(v, yr, function(z) sum((z - mean(z))^2)); n <- tapply(v, yr, length)
  out <- sqrt(sum(2 * n * ss) / sum(n * (n - 1)))
  stopifnot(is.finite(out), out > 0)
  out
}


## ---- manifest / digests ------------------------------------------------------------
file_digest <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
write_manifest <- function(entries, path) {
  ## entries: list of data.frames with file, nrow, ncol, sha256, note
  man <- bind_rows(entries) |> mutate(built = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  write.csv(man, path, row.names = FALSE)
  invisible(man)
}
## manifest may be the path of manifest.csv or a data.frame already read from it (a SNAPSHOT taken earlier; Codex review
## P2-5: 02_fit_model.R passes the manifest it read at load so that a rebuild between the snapshot and the load is caught)
## The file is read ONCE into memory; the digest and the returned object both come from those same bytes (Codex re-check
## P2, 2026-10-24: a separate readRDS() and file_digest() could see two different files if the file were replaced in
## between). .after_read is a test hook only: it is called after the bytes are captured, so a test can replace the file at
## that moment and confirm the returned object is the one that was hashed.
assert_matches_manifest <- function(path, manifest, .after_read = NULL) {
  man <- if (is.data.frame(manifest)) manifest else read.csv(manifest, stringsAsFactors = FALSE)
  stopifnot(all(c("file", "nrow", "sha256") %in% names(man)))
  row <- man[man$file == basename(path), , drop = FALSE]
  if (nrow(row) != 1) stop("file not in manifest: ", basename(path))
  n_bytes <- file.info(path)$size
  if (is.na(n_bytes)) stop("file not found: ", path)
  bytes <- readBin(path, what = "raw", n = n_bytes)
  if (length(bytes) != n_bytes) stop(sprintf("%s: read %d of %d bytes", basename(path), length(bytes), n_bytes))
  if (!is.null(.after_read)) .after_read()
  tmp <- tempfile(fileext = ".rds"); on.exit(unlink(tmp), add = TRUE)
  writeBin(bytes, tmp)
  obj <- readRDS(tmp)
  if (nrow(obj) != row$nrow) stop(sprintf("%s: nrow %d != manifest %d", basename(path), nrow(obj), row$nrow))
  if (digest::digest(bytes, algo = "sha256", serialize = FALSE) != row$sha256)
    stop(sprintf("%s: sha256 digest does not match the manifest%s", basename(path),
                 if (is.data.frame(manifest)) " snapshot taken at load (data rebuilt in between?)" else ""))
  obj
}
