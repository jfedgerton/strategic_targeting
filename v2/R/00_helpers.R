## ============================================================================
## 00_helpers.R -- small shared helpers for the v2 pipeline
## Sourced by every script. Keep this file small and concrete.
## ============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

V2_SEED <- 123L

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
pair_dispute_intervals <- function(dd, db, joiner_entry_unknown = TRUE) {
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
    select(dispnum, ccode1, ccode2, spell, pair_onset, pair_last, spell_first, spell_last, both_orig, start_exact, start_lo, start_hi, end_lo, end_hi)
  stopifnot(all(pdi$start_lo <= pdi$start_hi), all(pdi$end_lo <= pdi$end_hi), all(pdi$start_lo <= pdi$end_hi), !anyNA(pdi$start_lo), !anyNA(pdi$end_hi),
            all(pdi$spell_first[pdi$spell == 1L] == pdi$pair_onset[pdi$spell == 1L]))
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
collapse_duplicate_events <- function(events) {
  key_cols   <- c("initiator", "target", "year")
  same_cols  <- c("sidea1", "sidea2", "target_has_pos_tm1", "fighting_at_ref")               # must agree within a key
  ## the *_noww2 flags (Codex review P2-3) are the sample flags of the NON-WWII source rows only: ORing them over a key
  ## that mixes a WWII and a non-WWII dispute keeps the event in the no-WWII sample iff a non-WWII dispute qualifies
  or_cols    <- c("war", "war_disp", "escalation_year", "esc_fallback", "is_ww2", "fatal_main", "fatal_cow", "fatal_ub",
                  "in_mid", "in_fmid", "in_fmidcow", "in_fmidub", "in_war",
                  "in_mid_noww2", "in_fmid_noww2", "in_fmidcow_noww2", "in_fmidub_noww2", "in_war_noww2")   # never NA
  tri_cols   <- c("fatal1", "fatal2", "fatal_disp", "fatal", "fatal_mic", "fatal_gml")          # TRUE if any TRUE, NA if all NA
  max_cols   <- c("orig1", "orig2", "fatality1", "fatality2", "hostlev1", "hostlev2", "fatality_disp", "hostlev_disp", "mic_fatalmin")
  min_cols   <- c("dispnum", "ev_start_lo", "ev_start_hi")
  paste_cols <- c("fatal_source")
  known <- c(key_cols, same_cols, or_cols, tri_cols, max_cols, min_cols, paste_cols, "is_joiner_pair", "escalation")
  extra <- setdiff(names(events), known); missing <- setdiff(known, names(events))
  if (length(extra)) stop("collapse_duplicate_events: unclassified columns: ", paste(extra, collapse = ", "))
  if (length(missing)) stop("collapse_duplicate_events: expected columns absent: ", paste(missing, collapse = ", "))
  for (cc in or_cols) if (anyNA(events[[cc]])) stop("collapse_duplicate_events: NA in ", cc)
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
              is_joiner_pair = all(is_joiner_pair), escalation = FALSE, .groups = "drop")
  stopifnot(nrow(out) == nrow(distinct(keep, initiator, target, year)), sum(out$n_disputes) == nrow(keep),
            all(out$in_mid), !any(out$escalation),
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
    stop(sprintf("WWII dispute list derived from the data (%s) differs from WW2_DISPNUMS (%s); review and update 00_helpers.R",
                 paste(got, collapse = ","), paste(sort(WW2_DISPNUMS), collapse = ",")))
  cand |> select(dispnum, styear, stmon, y0, y1, n_states, hiact, hostlev, n_dyadyears_in_war139, is_ww2, reason)
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

## ---- fit provenance (audit 2026-10-23): a fit must have been made from the CURRENT data files ----------------
## obj$manifest holds the manifest rows (file, sha256) of the data the fit used; compare with the current manifest.
assert_fit_current <- function(obj, manifest_path) {
  stopifnot(!is.null(obj$manifest), nrow(obj$manifest) > 0, file.exists(manifest_path))
  cur <- read.csv(manifest_path, stringsAsFactors = FALSE)
  m <- merge(obj$manifest[, c("file", "sha256")], cur[, c("file", "sha256")], by = "file", all.x = TRUE, suffixes = c("_fit", "_now"))
  bad <- m$file[is.na(m$sha256_now) | m$sha256_fit != m$sha256_now]
  if (length(bad)) stop(sprintf("fit %s was made from data that no longer match the current build (%s); refit or rebuild", obj$cell, paste(bad, collapse = ", ")))
  invisible(TRUE)
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

## ---- Stan data assembly ----------------------------------------------------------------
## Builds the Stan data list from three data frames:
##   sy:   state-years, one row per (state_idx, year); columns y1, state_idx, W cols, X cols
##   alts: alternative rows, one row per (sy_row, target_idx); columns sy_row, target_idx, dM, dR, Z cols
##   ev:   events, columns sy_row, alt_row (index into alts after ordering by sy_row)
## alts must be sorted by sy_row; every sy row must own >= 1 alternative.
## Choice-set flags (Codex review P1-1, PI choice "A'", 2026-10-24). alts may carry two 0/1 columns, in_iv (the row is
## in the state-year's stage-1 IV set, defined at January of t) and in_s2 (the row is in the stage-2 choice set, defined
## at the state-year's reference date). A state-year's row block is the UNION of both sets; the Stan programs take
## IV[n] over the in_iv rows and the stage-2 softmax over the in_s2 rows. Absent columns mean both sets are the whole
## block (the simulator, the tests, every pre-review caller). Every state-year needs >= 1 in_iv row; every event's
## chosen row needs in_s2 = 1 and its stage-2 set >= 2 rows (a chosen-only set has no information); with sampled = 1
## the sets are resolved by sample_alternatives() and every row carries in_iv = in_s2 = 1.
## prior_mu0_sd: sd of the normal(0, .) prior on the mu intercept mu0 (the population-average ideal gap, in the
## units of x; since 2026-10-23 x is a log ratio of CINC shares / PageRank, no longer within-year SD units);
## 2 is the README value and was set when x was in SD units (flagged for PI review, not changed here).
build_stan_data <- function(sy, alts, ev, W_cols, X_cols, Z_cols, joint, sampled = 0L,
                            n_full = NULL, x_ref = NULL, K = NULL, anchor = NULL,
                            prior_beta_sd = 1, prior_mu0_sd = 2, allow_orphan_initiations = FALSE, S = NULL, prior_only = 0L,
                            D_mode = 3L, grainsize = NULL) {
  stopifnot(is.data.frame(sy), is.data.frame(alts), is.data.frame(ev))
  stopifnot(is.numeric(prior_mu0_sd), length(prior_mu0_sd) == 1, is.finite(prior_mu0_sd), prior_mu0_sd > 0)
  sampled <- as.integer(sampled)
  if (!"in_iv" %in% names(alts)) alts$in_iv <- 1L
  if (!"in_s2" %in% names(alts)) alts$in_s2 <- 1L
  alts$in_iv <- as.integer(alts$in_iv); alts$in_s2 <- as.integer(alts$in_s2)
  stopifnot(!anyNA(alts$in_iv), !anyNA(alts$in_s2), all(alts$in_iv %in% 0:1), all(alts$in_s2 %in% 0:1), all(alts$in_iv == 1L | alts$in_s2 == 1L))
  if (sampled == 1) {
    ## alts = IV rows (first n_iv, sorted by sy_row) followed by per-event blocks; ev carries ev_start/ev_end.
    ## sample_alternatives() draws the IV rows from the in_iv set and the event blocks from the in_s2 set, so the flags
    ## are resolved: every remaining row is in its own set.
    stopifnot(all(c("ev_start", "ev_end") %in% names(ev)), !is.null(n_full))
    if (!all(alts$in_iv == 1L & alts$in_s2 == 1L)) stop("build_stan_data: with sampled = 1 every row must have in_iv = in_s2 = 1 (use sample_alternatives())")
    n_iv <- min(ev$ev_start) - 1L
    iv_rows <- seq_len(n_iv)
    stopifnot(all(diff(alts$sy_row[iv_rows]) >= 0), all(ev$ev_start > n_iv), all(ev$ev_end <= nrow(alts)),
              all(ev$alt_row >= ev$ev_start & ev$alt_row <= ev$ev_end))
    if (any(ev$ev_end - ev$ev_start + 1L < 2L)) stop("build_stan_data: an event's sampled set holds only the chosen target")
  } else {
    n_iv <- nrow(alts); iv_rows <- seq_len(n_iv)
    stopifnot(all(diff(alts$sy_row) >= 0))                     # sorted by state-year
  }
  stopifnot(setequal(unique(alts$sy_row[iv_rows]), seq_len(nrow(sy))))  # every state-year has a row block
  ## every state-year has a NON-EMPTY IV set (rows with in_iv = 1)
  n_iv_sy <- tabulate(alts$sy_row[iv_rows][alts$in_iv[iv_rows] == 1L], nbins = nrow(sy))
  if (any(n_iv_sy == 0L)) stop(sprintf("build_stan_data: %d state-years have an EMPTY stage-1 IV set (rows exist but none has in_iv = 1)", sum(n_iv_sy == 0L)))
  stopifnot(all(sy$y1 %in% c(0L, 1L)))
  stopifnot(all(ev$sy_row %in% seq_len(nrow(sy))), all(ev$alt_row %in% seq_len(nrow(alts))))
  stopifnot(all(alts$sy_row[ev$alt_row] == ev$sy_row))         # chosen target belongs to the event's state-year
  stopifnot(all(sy$y1[ev$sy_row] == 1L))                        # every event sits in an initiating state-year
  ## the chosen row is in the stage-2 set, and that set is not the chosen row alone
  if (any(alts$in_s2[ev$alt_row] != 1L)) stop(sprintf("build_stan_data: %d events have a chosen row outside the stage-2 set (in_s2 = 0)", sum(alts$in_s2[ev$alt_row] != 1L)))
  if (sampled == 0) {
    n_s2_sy <- tabulate(alts$sy_row[alts$in_s2 == 1L], nbins = nrow(sy))
    bad_s2 <- unique(ev$sy_row[n_s2_sy[ev$sy_row] < 2L])
    if (length(bad_s2)) stop(sprintf("build_stan_data: %d initiating state-years have an EMPTY or CHOSEN-ONLY stage-2 set (< 2 rows with in_s2 = 1)", length(bad_s2)))
  }
  n_orphan <- sum(sy$y1[setdiff(seq_len(nrow(sy)), ev$sy_row)] == 1L)
  ## y1 == 1 without a stage-2 event happens in real data when the only target lacks positions;
  ## the simulator never produces it, so it is an error unless explicitly allowed.
  if (n_orphan > 0 && !allow_orphan_initiations) stop(n_orphan, " state-years have y1 = 1 but no event")
  S_data <- max(c(sy$state_idx, alts$target_idx))
  if (is.null(S)) S <- S_data else stopifnot(S >= S_data)   # override keeps one index space across subsets
  cs <- alts[iv_rows, , drop = FALSE] |> mutate(row = row_number()) |> group_by(sy_row) |>
    summarise(s = min(row), e = max(row), n = n(), n_iv = sum(in_iv), .groups = "drop") |> arrange(sy_row)
  stopifnot(all(cs$e - cs$s + 1 == cs$n), all(cs$n_iv >= 1L))
  if (is.null(n_full)) n_full <- cs$n_iv                        # full IV-set size (the sampling correction's numerator)
  stopifnot(length(n_full) == nrow(sy), all(n_full >= cs$n_iv))
  if (sampled == 1) { ev_start <- as.integer(ev$ev_start); ev_end <- as.integer(ev$ev_end) }
  else { ev_start <- as.integer(cs$s[ev$sy_row]); ev_end <- as.integer(cs$e[ev$sy_row]) }
  stopifnot(all(ev$alt_row >= ev_start & ev$alt_row <= ev_end))  # chosen target inside the event's own block
  if (is.null(x_ref)) x_ref <- colMeans(as.matrix(sy[, X_cols, drop = FALSE]))
  W <- as.matrix(sy[, W_cols, drop = FALSE]); X <- as.matrix(sy[, X_cols, drop = FALSE])
  Z <- as.matrix(alts[, Z_cols, drop = FALSE])
  assert_no_na(sy, c("y1", "state_idx", W_cols, X_cols), "stage-1 panel")
  assert_no_na(alts, c("sy_row", "target_idx", "dM", "dR", Z_cols), "alternatives")
  dat <- list(S = S, N1 = nrow(sy), N2 = nrow(alts), E = nrow(ev),
              P1 = ncol(W), P = ncol(X), Pz = ncol(Z),
              y1 = as.integer(sy$y1), sy_state = as.integer(sy$state_idx),
              W = W, X = X,
              cs_start = as.integer(cs$s), cs_end = as.integer(cs$e), n_full = as.integer(n_full),
              alt_target = as.integer(alts$target_idx), x = cbind(alts$dM, alts$dR), Z = Z,
              in_iv = as.integer(alts$in_iv), in_s2 = as.integer(alts$in_s2),
              ev_sy = as.integer(ev$sy_row), ev_row = as.integer(ev$alt_row), ev_start = ev_start, ev_end = ev_end,
              joint = as.integer(joint), sampled = sampled,
              x_ref = as.numeric(x_ref), prior_beta_sd = prior_beta_sd, prior_mu0_sd = prior_mu0_sd,
              prior_only = as.integer(prior_only), D_mode = as.integer(D_mode))
  ## reduce_sum_static slice size over state-years: NULL -> one slice (exact serial sum, threads = 1)
  dat$grainsize <- as.integer(if (is.null(grainsize)) dat$N1 else max(1L, grainsize))
  stopifnot(dat$D_mode %in% 1:3, dat$grainsize >= 1)
  dimnames(dat$W) <- NULL; dimnames(dat$X) <- NULL; dimnames(dat$Z) <- NULL
  if (!is.null(K)) {
    stopifnot(K >= 1, K <= S - 1)
    if (is.null(anchor)) anchor <- choose_anchors(sy, ev, K, S)
    stopifnot(length(anchor) == K, !anyDuplicated(anchor), all(anchor %in% seq_len(S)))
    dat$K <- as.integer(K); dat$anchor <- array(as.integer(anchor), dim = K)   # array so K = 1 is a 1-d array for Stan
  }
  dat
}

## Anchor states for AME identification: the K states with the most initiation events.
## Ties broken by state index (deterministic). NOT string-sorted ccodes (audit F6).
choose_anchors <- function(sy, ev, K, S = max(sy$state_idx)) {
  cnt <- tabulate(sy$state_idx[ev$sy_row], nbins = S)
  ord <- order(-cnt, seq_along(cnt))
  ord[seq_len(K)]
}

## ---- sampled alternatives --------------------------------------------------------------
## Two outcome-separated samples (fixes the Codex review items 3 and 4):
##   IV rows    : for EVERY state-year, m alternatives drawn uniformly without replacement from
##                its full set, ignoring outcomes (chosen targets are NOT forced in, so the IV
##                estimate has the same distribution for initiating and non-initiating rows).
##   event rows : EACH event gets its own set = its chosen target + m alternatives drawn
##                uniformly from the remaining alternatives, independently per event (so two
##                events of one state-year never share an outcome-dependent set).
## Returns list(alts = IV rows then event rows (sorted), ev with alt_row / ev_start / ev_end into
## the event block, n_full per state-year, n_iv = number of IV rows).
## Choice-set flags (Codex review P1-1): the IV rows are drawn from the state-year's in_iv rows (n_full counts them),
## the event blocks from its in_s2 rows; absent flags mean the whole block. The output rows carry in_iv = in_s2 = 1
## because every row is then in the set it was drawn for.
## Full choice sets (m = 0 / sampled = 0) are the main specification; this is an option only.
sample_alternatives <- function(alts, ev, m, seed = V2_SEED) {
  stopifnot(m >= 1, all(diff(alts$sy_row) >= 0))
  if (!"in_iv" %in% names(alts)) alts$in_iv <- 1L
  if (!"in_s2" %in% names(alts)) alts$in_s2 <- 1L
  stopifnot(all(alts$in_iv %in% 0:1), all(alts$in_s2 %in% 0:1), all(alts$in_iv == 1L | alts$in_s2 == 1L))
  set.seed(seed)
  alts <- alts |> mutate(.row = row_number())
  ## (i) IV rows: uniform m per state-year from the IV set, outcome-blind
  iv <- alts |> filter(in_iv == 1L) |> group_by(sy_row) |> mutate(.u = runif(n()), .rank = rank(.u, ties.method = "first")) |>
    filter(.rank <= m) |> ungroup() |> arrange(sy_row, .row) |> select(-.u, -.rank)
  n_full <- alts |> filter(in_iv == 1L) |> count(sy_row) |> arrange(sy_row) |> pull(n)
  stopifnot(length(n_full) == length(unique(alts$sy_row)), setequal(unique(iv$sy_row), unique(alts$sy_row)))
  ## (ii) event rows: per event, chosen + m others from the same state-year's stage-2 set
  if (any(alts$in_s2[ev$alt_row] != 1L)) stop("sample_alternatives: a chosen row is not in the stage-2 set (in_s2 = 0)")
  ev_blocks <- vector("list", nrow(ev))
  for (e in seq_len(nrow(ev))) {
    rows_n <- alts$.row[alts$sy_row == ev$sy_row[e] & alts$in_s2 == 1L]
    others <- setdiff(rows_n, ev$alt_row[e])
    pick <- if (length(others) > m) sample(others, m) else others
    blk <- alts[c(ev$alt_row[e], sort(pick)), , drop = FALSE]
    blk$.event <- e
    ev_blocks[[e]] <- blk
  }
  evrows <- bind_rows(ev_blocks)
  n_iv <- nrow(iv)
  out_alts <- bind_rows(iv |> mutate(.event = 0L), evrows)
  ev2 <- ev
  ev2$ev_start <- n_iv + match(seq_len(nrow(ev)), evrows$.event)                 # first row of each block (the chosen target)
  ev2$ev_end   <- n_iv + nrow(evrows) + 1L - match(seq_len(nrow(ev)), rev(evrows$.event))
  ev2$alt_row  <- ev2$ev_start                                                     # chosen target is the block's first row
  stopifnot(all(out_alts$sy_row[ev2$alt_row] == ev2$sy_row),
            all(out_alts$.row[ev2$alt_row] == alts$.row[ev$alt_row]),            # same chosen target as before
            all(ev2$ev_end >= ev2$ev_start), all(ev2$ev_start > n_iv))
  out_alts$in_iv <- 1L; out_alts$in_s2 <- 1L                                       # resolved: every row is in the set it was drawn for
  list(alts = out_alts |> select(-.row, -.event), ev = ev2, n_full = n_full, n_iv = n_iv)
}

## ---- brute-force stage-2 likelihood in R (used by tests and the Fixed_param check) -----
## pars: list(mu0 (D), beta (P x D), a_i (S x D), b_j (S), gamma (Pz), sigma_x (D), rho, U (S x K) or NULL, V)
## D = 2 when dat$D_mode == 3; D = 1 uses column dat$D_mode of x (1 = dM, 2 = dR).
## mu = mu0 + X beta + a_i; mu0 is required (no silent zero default) so a caller that forgets it fails here.
## lse[n] is over the state-year's in_iv rows, the event softmax over its in_s2 rows (sampled = 0) or the event's own
## block (sampled = 1); dat$in_iv / dat$in_s2 are required (build_stan_data() always writes them).
stage2_loglik_R <- function(dat, pars) {
  D <- if (dat$D_mode == 3L) 2L else 1L
  xcols <- if (D == 2L) 1:2 else dat$D_mode
  if (is.null(pars$mu0)) stop("stage2_loglik_R: pars$mu0 (length D intercept of mu) is missing")
  stopifnot(length(pars$mu0) == D, all(is.finite(pars$mu0)))
  if (is.null(dat$in_iv) || is.null(dat$in_s2)) stop("stage2_loglik_R: dat$in_iv / dat$in_s2 are missing (build the data with build_stan_data())")
  stopifnot(length(dat$in_iv) == dat$N2, length(dat$in_s2) == dat$N2)
  Mu <- rep(1, dat$N1) %o% as.numeric(pars$mu0) +
    dat$X %*% matrix(pars$beta, ncol = D) + matrix(pars$a_i, ncol = D)[dat$sy_state, , drop = FALSE]
  zg <- if (dat$Pz > 0) as.numeric(dat$Z %*% pars$gamma) else rep(0, dat$N2)
  util <- numeric(dat$N2); lse <- numeric(dat$N1)
  Sinv <- if (D == 2L) solve(matrix(c(pars$sigma_x[1]^2, pars$rho * prod(pars$sigma_x),
                                      pars$rho * prod(pars$sigma_x), pars$sigma_x[2]^2), 2))
          else matrix(1 / pars$sigma_x[1]^2, 1, 1)
  one_util <- function(r, n) {
    d <- dat$x[r, xcols] - Mu[n, ]
    u <- -0.5 * as.numeric(t(d) %*% Sinv %*% d) + pars$b_j[dat$alt_target[r]] + zg[r]
    if (!is.null(pars$U)) u <- u + sum(pars$U[dat$sy_state[n], ] * pars$V[dat$alt_target[r], ])
    u
  }
  for (n in seq_len(dat$N1)) {
    rows <- dat$cs_start[n]:dat$cs_end[n]
    for (r in rows) util[r] <- one_util(r, n)
    rows_iv <- rows[dat$in_iv[rows] == 1L]
    stopifnot(length(rows_iv) >= 1)
    lse[n] <- log(sum(exp(util[rows_iv])))
  }
  loglik <- numeric(dat$E)
  for (e in seq_len(dat$E)) {
    rows <- dat$ev_start[e]:dat$ev_end[e]
    if (dat$sampled == 1) for (r in rows) util[r] <- one_util(r, dat$ev_sy[e])
    else rows <- rows[dat$in_s2[rows] == 1L]
    stopifnot(dat$ev_row[e] %in% rows)
    loglik[e] <- util[dat$ev_row[e]] - log(sum(exp(util[rows])))
  }
  list(loglik = loglik, lse = lse, util = util)
}

## ---- Procrustes alignment of U/V draws across chains ---------------------------------
## U_draws: array [iter, chain, S, K]. Rotates each chain's draws (U R, V R) to align the
## chain-mean U with the reference chain's mean U via orthogonal Procrustes (R = W V' from
## the SVD of U_c' U_ref). u_i'v_j is invariant to R, so the likelihood is unchanged.
procrustes_align <- function(U_draws, V_draws, ref_chain = 1) {
  stopifnot(length(dim(U_draws)) == 4, all(dim(U_draws) == dim(V_draws)))
  n_iter <- dim(U_draws)[1]; n_chain <- dim(U_draws)[2]; K <- dim(U_draws)[4]
  U_ref <- apply(U_draws[, ref_chain, , , drop = FALSE], c(3, 4), mean)
  rotations <- vector("list", n_chain)
  for (ch in seq_len(n_chain)) {
    U_c <- apply(U_draws[, ch, , , drop = FALSE], c(3, 4), mean)
    sv <- svd(t(U_c) %*% U_ref)
    R <- sv$u %*% t(sv$v)                       # K x K orthogonal
    rotations[[ch]] <- R
    for (it in seq_len(n_iter)) {
      U_draws[it, ch, , ] <- matrix(U_draws[it, ch, , ], ncol = K) %*% R
      V_draws[it, ch, , ] <- matrix(V_draws[it, ch, , ], ncol = K) %*% R
    }
  }
  list(U = U_draws, V = V_draws, rotations = rotations)
}

## grainsize for reduce_sum_static: with 1 thread one slice (serial); otherwise about 4 slices per thread
choose_grainsize <- function(N1, threads) if (threads <= 1) N1 else max(1L, as.integer(ceiling(N1 / (4 * threads))))
