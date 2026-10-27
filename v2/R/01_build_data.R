## ============================================================================
## 01_build_data.R -- ONE script that writes EVERY data file for the v2 pipeline
##
## Usage:  Rscript v2/R/01_build_data.R
## Inputs: peacesciencer bundled data (COW state system, NMC, COW majors, Polity/V-Dem,
##         COW-MID 5.0 directed dispute-years, ATOP, COW contiguity, COW IGO v3).
##         No download_extdata() call is needed for the columns used here.
## Outputs (v2/data/, all gitignored; manifest.csv records nrow + sha256 per file):
##   sy_positions.rds              state-year positions M = log CINC share, R_pr = log IGO PageRank (both raw log
##                                 shares, NOT standardised within year; PI decision 2026-10-23), R_eig (eigenvector
##                                 centrality, still std within year; not used in the main grid)
##   events.rds                    one row per (initiator, target, year) directed event among the kept rows (initiator = COW's
##                                 war-initiating side in war-level disputes linked to a COW war, else MID Side A; D1 2026-10-25:
##                                 dir_flipped, dir_source, dir_rule, war_link; stage-2 weights w_<sample>, D8) (same-year
##                                 duplicates from different disputes collapsed; dispnum = smallest, dispnums lists all,
##                                 n_disputes), joiners included (is_joiner_pair), plus one row per escalation event
##                                 (escalation = TRUE under the month rule; excluded from every sample flag and from y1;
##                                 escalation_year = the old year rule, esc_fallback when the year rule decided),
##                                 ev_start_lo / ev_start_hi (the event's date interval, month index), fighting_at_ref
##                                 (pair already fighting at the state-year's reference date; stage 2 only), war (pair
##                                 level) and war_disp (dispute level), is_ww2 (appendix option --drop_ww2),
##                                 target_has_pos_tm1, sample flags in_mid / in_fmid / in_fmidcow / in_fmidub / in_war
##                                 (+ fatal_source), and in_*_noww2 (the sample flags of the non-WWII disputes only,
##                                 used by --drop_ww2; Codex review P2-3); coal_incl / coal_strict (coalition context of the
##                                 event's D8 group under the inclusive / strict tie convention, NA on escalation rows) and
##                                 coal_onset_ambiguous (the `context` switch; dual typology stage a, 2026-10-25)
##   changes_lag{h}.rds            ONE state-level change table per lag (build_change_table() in helpers_data.R):
##                                 raw change dM_raw_h / dRpr_raw_h = level(t-1) - level(t-1-h), and the transformed
##                                 change dM_h / dRpr_h = (raw - year mean) / one pooled SD per lag. Feeds BOTH the
##                                 initiator change in stage 1 and the target change in stage 2.
##   change_scale.csv              the pooled SD constants per lag and variable
##   stage1_panel_lag{h}.rds       state-year panel for h in LAGS: y1_{mid,fmid,fmidcow,fmidub,war} (escalations
##                                 excluded), W, X, ongoing_tm1 (participant in any MID active at t-1), the power levels M_lvl /
##                                 R_lvl at t-1 (within-year centred, gap-scaled; D2) and log_n_system_tm1 (D3); multi_mid_jan /
##                                 n_multi_jan (a multilateral MID the state had not entered is active at 1 January of t; the M1c
##                                 stage-1 control of dual typology stage b, 2026-10-27; not in the main W)
##   alts_lag{h}_{all,polrel}.rds  flattened choice sets (one row per state-year x alternative): the UNION of the
##                                 stage-1 IV set (in_iv = 1: pair not already fighting at JANUARY of t, every state-year)
##                                 and the stage-2 set (in_s2 = 1: pair not already fighting at the state-year's reference
##                                 date), rows in neither set removed (Codex review P1-1, PI choice "A'");
##                                 gaps dM / dR_pr are log ratios; target changes dM_j_h / dRpr_j_h (+ dch_j_missing
##                                 when the target has no change value); dy_mid_tm1 (pair in a MID active at t-1);
##                                 log_capdist = log1p(capital distance km) for the --distance option
##   build_log.txt                 dropped years / rows and why
##   manifest.csv
##
## Timing conventions (audit F5): event onset year = t. Positions x_ij = (M_i - M_j, R_i - R_j)
## are measured at t-1. Power changes (initiator AND target) are level(t-1) - level(t-1-h), h in LAGS.
## ALL lags are calendar-year joins on (ccode, year - h); never row lag().
##
## Initiation coding (audit F4; PI decision 2026-10-24): initiator = Side A participant (sidea == 1),
## target = Side B participant. One direction per dispute pair, onset year = the pair's first year in the
## dispute. Since 2026-10-25 (D1) war-level disputes linked to a COW inter-state war take COW's war-initiating
## side instead (Side A is the first codable action, not the aggressor); V2_DIRECTION=sideA rebuilds the old coding. For multilateral disputes every (Side A x Side B) pair is one event. Since 2026-10-24 joiners
## (orig == 0) are part of the main build (INCLUDE_JOINERS = TRUE); is_joiner_pair marks pairs with a joiner
## on either side (tabulation only, not a covariate). INCLUDE_JOINERS = FALSE reproduces the originator-only
## event base.
##
## Ongoing-dispute rule (PI decision 2026-10-24, month resolution since the second round of that day). Stage 2 is
## target SELECTION. An event is an escalation when i and j are both participants in a DIFFERENT dispute that
## started before the event's dispute start date and had not ended by it (month resolution; dispute-level COW
## dates from cow_mid_disps, end dates bounded from the duration because peacesciencer ships no end month; the
## year rule decides where the month comparison is undecided and the fallback is counted). Escalations stay in
## events.rds with escalation = TRUE but leave every sample flag, stage 2 and y1 / n_ev in stage 1. For choice
## sets each state-year has TWO sets (Codex review P1-1): the stage-1 IV set drops j when i and j are in a dispute
## active at January of t that started before it (every state-year, no event information, so IV cannot depend on
## y1); the stage-2 set uses the state-year's reference date (start of its earliest non-escalation event, else
## January of t) and drops j when i and j are in a different dispute active at that date that started before it.
## Joiner pairs entering a dispute in its first year have an unknown entry month, [dispute start, December] (P1-2).
## War is coded at the pair level (dispute war-level AND at least one participant at hostility level 5), and same-year
## events of one (initiator, target) pair from different disputes collapse into one event (helpers in helpers_data.R).
##
## Conflict-history controls (new in v2), all relative to the change window [t-1-h, t-1]:
##   (a) mon_prior_mids_h : # distinct MIDs the state was a participant in (either side, joiners included
##                          when INCLUDE_JOINERS, consistent with the events)
##                          with onset in the 10 years BEFORE the window: [t-1-h-10, t-1-h-1]
##   (b) dy_prior_mids_h  : # distinct MIDs between i and j with onset in that same period
##   (c) dy_peace_years_h : years since the last MID-ongoing year between i and j, as of t-1-h
##   plus mon_peace_years_h (monadic analogue of (c)).
## Lagged involvement at t-1 (PI decision 2026-10-24), always from every participant including joiners:
##   ongoing_tm1 (stage 1): the state was in any MID active in t-1
##   dy_mid_tm1  (stage 2): i and j were both in a MID active in t-1 (replaces dy_peace_years in the main Z)
## ============================================================================
suppressPackageStartupMessages({
  library(peacesciencer); library(igraph); library(dplyr); library(tidyr); library(digest)
})
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(this_file) == 0) this_file <- "v2/R/01_build_data.R"
source(file.path(dirname(this_file), "spec.R")); source(file.path(dirname(this_file), "helpers_data.R"))
Sys.setenv(V2_ROOT = normalizePath(file.path(dirname(this_file), "..")))
set.seed(V2_SEED)
## ---- initiator direction (PI decision 2026-10-25, D1; README "Initiator direction") -------------------------------------
## DIRECTION = "cowwar" (main): in war-level disputes linked to a COW inter-state war the initiator is the state on COW's
## war-initiating side (pair-level dyad-year match first, dispute-level side mapping for the unmatched pairs); everything else
## keeps the Side A rule. DIRECTION = "sideA" (appendix rebuild of the pre-2026-10-25 coding): Side A everywhere. Selected by the
## `direction` switch of spec.R: `Rscript v2/R/01_build_data.R --direction=sideA` (the environment variable V2_DIRECTION=sideA
## still works). The two codings write to SEPARATE directories, v2/data (cowwar) and v2/data_sideA (sideA), each with its own
## manifest and build_log, so they never overwrite each other; 02_fit_model.R --direction=sideA reads the second one and
## suffixes its cells with _sideA. The direction is the only switch the build takes; every other switch acts at assembly.
## the build takes --direction only: the archived flags and every other switch of spec.R are refused, not ignored (audit fix 14, 2026-10-26)
refuse_archived_options(commandArgs(trailingOnly = TRUE), "01_build_data.R", also = setdiff(names(SPEC_CHOICES), "direction"),
                        also_msg = "the build takes --direction only; every other switch acts at assembly, in 02_fit_model.R")
DIRECTION <- parse_spec(c(commandArgs(trailingOnly = TRUE), paste0("--direction=", Sys.getenv("V2_DIRECTION", unset = "cowwar"))))$direction
if (!DIRECTION %in% SPEC_CHOICES$direction) stop("--direction / V2_DIRECTION must be 'cowwar' (main) or 'sideA' (appendix); got ", DIRECTION)
log_open(paste0("01_build_data", if (DIRECTION == "sideA") "_sideA" else ""))
out_dir <- file.path(v2_root(), if (DIRECTION == "sideA") "data_sideA" else "data"); dir.create(out_dir, showWarnings = FALSE)
build_log <- file.path(out_dir, "build_log.txt"); cat("", file = build_log)
blog <- function(fmt, ...) { m <- logf(fmt, ...); cat(m, "\n", file = build_log, append = TRUE) }

## ---- flags (PI decisions; see README) ----------------------------------------------
LAGS            <- c(1L, 3L, 5L, 10L)
WAR_DIR_WINDOW  <- 1L             # D1: a MID pair with onset t links to a COW war dyad-year in [t, t + WAR_DIR_WINDOW] (MID onsets precede war onsets)
WAR_DIR_MIN_LINK <- 2L            # D1: dispute-level mapping needs at least this many pair-decided events behind the linked war and the side mapping
PART_DATES      <- TRUE           # D5: participant-level start / end months from peacesciencer::gml_part where the dispnum matches
INCLUDE_JOINERS <- TRUE           # main build since 2026-10-24: Side A x Side B pairs including joiners (FALSE = originators only)
IGO_FILL        <- "drop"         # "drop": years with all-missing IGO data are dropped from positions
                                  # "locf": R at year y uses the latest IGO year <= y (igo_src_year recorded)
PR_CONT_MAX     <- 5L             # politically relevant: conttype in 1..PR_CONT_MAX OR either a major power
PREWINDOW_YEARS <- 10L            # conflict-history look-back before the change window starts
POLITY_DEM_CUT  <- 6L             # democracy = polity2 >= 6 (NA -> 0, flagged in *_dem_missing)
FIRST_YEAR      <- 1816L
blog("flags: DIRECTION=%s (output dir %s) LAGS=%s INCLUDE_JOINERS=%s IGO_FILL=%s PR_CONT_MAX=%d PREWINDOW_YEARS=%d WAR_DIR_WINDOW=%d WAR_DIR_MIN_LINK=%d PART_DATES=%s",
     DIRECTION, basename(out_dir), paste(LAGS, collapse = ","), INCLUDE_JOINERS, IGO_FILL, PR_CONT_MAX, PREWINDOW_YEARS, WAR_DIR_WINDOW, WAR_DIR_MIN_LINK, PART_DATES)

## ============================================================================
## 1. State-years, CINC share, majors, democracy
## ============================================================================
mid_years <- range(cow_mid_dirdisps$year)
LAST_YEAR <- mid_years[2]
blog("COW-MID directed dispute-years cover %d-%d; panel truncated at %d", mid_years[1], mid_years[2], LAST_YEAR)
states <- create_stateyears(system = "cow", subset_years = FIRST_YEAR:LAST_YEAR) |>
  add_nmc() |> add_cow_majors() |> add_democracy() |>
  mutate(ccode = as.integer(ccode), year = as.integer(year))
assert_no_dup_keys(states, c("ccode", "year"), "state-years")
states <- states |> group_by(year) |>
  mutate(n_system = n(), cinc_share = cinc / sum(cinc, na.rm = TRUE)) |> ungroup()
blog("state-years: %d rows, %d states, cinc NA in %d rows", nrow(states), n_distinct(states$ccode), sum(is.na(states$cinc)))
entry <- states |> group_by(ccode) |> summarise(entry_year = min(year), .groups = "drop")

## ============================================================================
## 2. IGO networks -> PageRank and eigenvector centrality, standardised within year
## ============================================================================
igo <- create_dyadyears(system = "cow", directed = FALSE, mry = FALSE, subset_years = FIRST_YEAR:LAST_YEAR) |>
  add_igos() |> mutate(ccode1 = as.integer(ccode1), ccode2 = as.integer(ccode2), year = as.integer(year))
igo_avail <- igo |> group_by(year) |> summarise(n = n(), n_na = sum(is.na(dyadigos)), .groups = "drop") |>
  mutate(all_missing = n_na == n)
blog("IGO dyad-years: %d rows; years with ALL dyadigos missing: %d of %d (listed in build_log)",
     nrow(igo), sum(igo_avail$all_missing), nrow(igo_avail))
cat("IGO all-missing years: ", paste(igo_avail$year[igo_avail$all_missing], collapse = " "), "\n", file = build_log, append = TRUE)
part_missing <- igo_avail |> filter(!all_missing, n_na > 0)
if (nrow(part_missing)) blog("IGO partially missing years (NA dyads set to 0 ties): %s", paste(part_missing$year, collapse = " "))
igo_years <- igo_avail$year[!igo_avail$all_missing]

cent_rows <- vector("list", length(igo_years))
for (k in seq_along(igo_years)) {
  yy <- igo_years[k]
  members <- states$ccode[states$year == yy]
  ed <- igo |> filter(year == yy, !is.na(dyadigos), dyadigos > 0) |> select(ccode1, ccode2, weight = dyadigos)
  g <- graph_from_data_frame(ed |> mutate(ccode1 = as.character(ccode1), ccode2 = as.character(ccode2)),
                             directed = FALSE, vertices = data.frame(name = as.character(members)))
  pr <- page_rank(g, weights = E(g)$weight, directed = FALSE)$vector
  ev <- eigen_centrality(g, weights = E(g)$weight, directed = FALSE)$vector
  stopifnot(length(pr) == length(members), !anyNA(pr), !anyNA(ev))
  cent_rows[[k]] <- data.frame(ccode = as.integer(names(pr)), year = yy, pagerank = as.numeric(pr), eigen = as.numeric(ev))
}
cent <- bind_rows(cent_rows)
assert_no_dup_keys(cent, c("ccode", "year"), "centrality")
blog("centrality computed for %d IGO years (%d-%d)", length(igo_years), min(igo_years), max(igo_years))

## ============================================================================
## 3. Positions: M = log CINC share; R_pr = log IGO PageRank (raw log shares, PI decision 2026-10-23)
##    R_eig = eigenvector centrality std within year (unchanged; not in the main grid)
## ============================================================================
pos <- states |> select(ccode, year, cinc_share, n_system)
## Material position M = LOG CINC share. A zero or NA share has no log: check on the FULL state-year panel,
## BEFORE any filtering (Codex review of 36c0f7b: the old check ran after NA rows were already dropped, so a
## missing share silently removed the state from positions and choice sets). Stop and list the state-years.
bad_cinc <- pos |> filter(is.na(cinc_share) | cinc_share <= 0)
if (nrow(bad_cinc) > 0) {
  blog("cinc_share <= 0 or NA in %d state-years: %s", nrow(bad_cinc),
       paste(sprintf("%d-%d", bad_cinc$ccode, bad_cinc$year), collapse = " "))
  stop(sprintf("%d state-years have cinc_share <= 0 or NA (see build_log.txt); log(CINC share) undefined. PI decision needed.", nrow(bad_cinc)))
}
if (IGO_FILL == "drop") {
  pos <- pos |> inner_join(cent |> mutate(igo_src_year = year), by = c("ccode", "year"))
} else if (IGO_FILL == "locf") {
  ## for every state-year, take centrality from the latest IGO year <= year (calendar join)
  cand <- pos |> select(ccode, year) |>
    inner_join(cent |> rename(igo_src_year = year), by = "ccode", relationship = "many-to-many") |>
    filter(igo_src_year <= year) |> group_by(ccode, year) |> slice_max(igo_src_year, n = 1, with_ties = FALSE) |> ungroup()
  pos <- pos |> inner_join(cand, by = c("ccode", "year"))
  blog("IGO_FILL = locf: max staleness %d years", max(pos$year - pos$igo_src_year))
} else stop("unknown IGO_FILL")
## PageRank has the same requirement as the CINC share: it must be finite and > 0 for log(pagerank) to exist.
## page_rank() on a graph with every system member as a vertex gives every vertex positive mass (damping),
## so a violation means a construction error upstream. Stop and list the state-years.
bad_pr <- pos |> filter(!is.finite(pagerank) | pagerank <= 0)
if (nrow(bad_pr) > 0) {
  blog("pagerank <= 0 or non-finite in %d state-years: %s", nrow(bad_pr),
       paste(sprintf("%d-%d", bad_pr$ccode, bad_pr$year), collapse = " "))
  stop(sprintf("%d state-years have pagerank <= 0 or non-finite (see build_log.txt); log(PageRank) undefined. PI decision needed.", nrow(bad_pr)))
}
## Positions (PI decision 2026-10-23). M = log CINC share and R_pr = log PageRank, both WITHOUT within-year
## standardisation. CINC share and PageRank each sum to 1 over the system in a year, so M and R_pr are log shares
## of the system total (world capabilities, network centrality) and the stage-2 gaps M_i - M_j, R_pr_i - R_pr_j
## are log ratios between the two states. Within-year standardisation (2026-10-22) made the change variables
## move with system growth (cor with 3-year log growth in the number of states about +0.29 for M, +0.17 for R
## at h = 3) and left them with SD 0.11 / 0.25, so per-unit priors dominated. Zero/NA shares were rejected above.
## R_eig stays EXACTLY as before (eigenvector centrality standardised within year, with the sd = 0 guard): it is
## not used in the main grid (the eigenvector robustness run is unscheduled) and was not part of this decision.
stopifnot(!anyNA(pos$cinc_share), all(pos$cinc_share > 0), all(is.finite(pos$pagerank)), all(pos$pagerank > 0))
pos <- pos |> mutate(log_cinc_share = log(cinc_share), M = log_cinc_share, R_pr = log(pagerank))
pos <- pos |> group_by(year) |> filter(n() > 2) |>
  mutate(R_eig = if (sd(eigen) > 0) as.numeric(scale(eigen)) else 0) |> ungroup()
assert_no_dup_keys(pos, c("ccode", "year"), "positions")
assert_no_na(pos, c("M", "R_pr", "R_eig"), "positions")
assert_finite(pos, c("log_cinc_share", "M", "R_pr", "R_eig"), "positions")
stopifnot(all(abs(pos$M - log(pos$cinc_share)) < 1e-12), all(abs(pos$R_pr - log(pos$pagerank)) < 1e-12))
m_skew <- mean((pos$M - mean(pos$M))^3) / sd(pos$M)^3
raw <- as.numeric(unlist(tapply(pos$cinc_share, pos$year, scale)))
within_sd <- pos |> group_by(year) |> summarise(sM = sd(M), sR = sd(R_pr), .groups = "drop")
blog("M = log CINC share (raw log share, no within-year standardisation): pooled SD %.3f, within-year SD range [%.3f, %.3f], skewness %.3f (within-year z of the raw share would be %.3f)",
     sd(pos$M), min(within_sd$sM), max(within_sd$sM), m_skew, mean((raw - mean(raw))^3) / sd(raw)^3)
blog("R_pr = log IGO PageRank (raw log share, no within-year standardisation): pooled SD %.3f, within-year SD range [%.3f, %.3f]; R_eig unchanged (eigenvector, std within year)",
     sd(pos$R_pr), min(within_sd$sR), max(within_sd$sR))
dropped_years <- setdiff(unique(states$year), unique(pos$year))
blog("positions: %d state-years over %d years; %d calendar years have NO positions (IGO or CINC missing)",
     nrow(pos), n_distinct(pos$year), length(dropped_years))
cat("years without positions: ", paste(sort(dropped_years), collapse = " "), "\n", file = build_log, append = TRUE)
saveRDS(pos, file.path(out_dir, "sy_positions.rds"))

## ============================================================================
## 4. Directed initiation events from COW-MID directed dispute-years (audit F4)
## ============================================================================
dd <- cow_mid_dirdisps |> mutate(across(c(dispnum, ccode1, ccode2, year), as.integer))
## dispute-level fatality / hostility
disp <- cow_mid_disps |> transmute(dispnum = as.integer(dispnum), fatality_disp = fatality, hostlev_disp = hostlev)
assert_no_dup_keys(disp, "dispnum", "cow_mid_disps")
## one row per (dispnum, ccode1, ccode2): the onset year of the pair within the dispute, with the sides of the pair's EARLIEST
## year (P3 fix 2026-10-25: first() depended on row order; a state can switch sides within a dispute)
pair <- build_pair_table(dd, disp)
side_chg <- attr(pair, "side_changes")
blog("pair table: %d directed dispute pairs; pairs whose sidea differs across their dispute-years (earliest year's value used): %d (%s)",
     nrow(pair), nrow(side_chg), if (nrow(side_chg)) paste(sprintf("dispute %d %d-%d", side_chg$dispnum, side_chg$ccode1, side_chg$ccode2), collapse = "; ") else "none")
## A -> B direction only, originators on both sides (main) or any participants (robustness)
events <- pair |> filter(sidea1 == 1L, sidea2 == 0L)
if (!INCLUDE_JOINERS) events <- events |> filter(orig1 == 1L, orig2 == 1L)
events <- events |> rename(initiator = ccode1, target = ccode2)
## ---- initiator direction for war-level disputes (PI decision 2026-10-25, D1; helpers and the rule in helpers_data.R) -----------
## COW-MID Side A = first codable action, not the aggressor (dispute 258: Side A = Poland and the Allies). For disputes at
## hostility level 5 that link to a COW inter-state war the direction follows COW's war-initiating side: pair-level dyad-year
## match (year t or t + WAR_DIR_WINDOW), then the dispute-level side mapping for the unmatched pairs of a linked dispute; ambiguous
## cases and unlinked disputes keep Side A. Under DIRECTION = "sideA" the decisions are computed and logged but NOT applied.
war_dy <- war_dyad_years(cow_war_inter)
blog("COW inter-state wars: %d wars, %d distinct dyad-years (both directions); every war has exactly one initiating side (asserted)",
     n_distinct(war_dy$warnum), nrow(war_dy))
n_ev_dir <- nrow(events)
events <- cow_war_direction(events, war_dy, window = WAR_DIR_WINDOW, min_link = WAR_DIR_MIN_LINK)
pm <- attr(events, "pair_matches"); dl <- attr(events, "dispute_links")
stopifnot(nrow(events) == n_ev_dir)
## match-rate details (the t + 1 window is justified below; the alternatives are logged for the PI)
wl <- events |> filter(hostlev_disp == 5L) |> select(dispnum, initiator, target, year)
match_in <- function(w) nrow(wl |> semi_join(war_dy |> rename(initiator = ccode1, target = ccode2, wyear = year) |> inner_join(wl, by = c("initiator", "target"), relationship = "many-to-many") |>
                                       filter(wyear >= year, wyear <= year + w) |> distinct(dispnum, initiator, target, year), by = c("dispnum", "initiator", "target", "year")))
blog("direction (D1): war-level dispute events %d in %d disputes; pair matched to a COW war dyad-year with the window t..t+%d: %d (t only %d; t..t+1 %d; t..t+2 %d; any ongoing year >= t %d); matched events with a war ONSET in the window %d, ongoing-only %d; pair-level decisions: flip %d, keep (confirmed) %d, ambiguous (two wars, opposite initiators) %d",
     nrow(wl), n_distinct(wl$dispnum), WAR_DIR_WINDOW, sum(events$dir_rule %in% c("pair_match", "pair_ambiguous")), match_in(0L), match_in(1L), match_in(2L), match_in(500L),
     nrow(pm |> group_by(dispnum, initiator, target, year) |> summarise(o = any(cowinteronset == 1L), .groups = "drop") |> filter(o)),
     nrow(pm |> group_by(dispnum, initiator, target, year) |> summarise(o = any(cowinteronset == 1L), .groups = "drop") |> filter(!o)),
     sum(events$dir_rule == "pair_match" & events$dir_flipped), sum(events$dir_rule == "pair_match" & !events$dir_flipped), sum(events$dir_rule == "pair_ambiguous"))
blog("direction (D1): why t..t+1. The MID onset year is the pair's first year in the dispute (the crisis), the COW war dyad-year begins when the fighting reached war level; of the %d matched events %d match at t and %d only at t+1; widening to t+2 would add %d more, any later ongoing year %d more (not used)",
     match_in(WAR_DIR_WINDOW), match_in(0L), match_in(1L) - match_in(0L), match_in(2L) - match_in(1L), match_in(500L) - match_in(1L))
blog("direction (D1): dispute-level links: %d war-level disputes with pair-decided events, %d linked to one COW war by a strict majority of at least %d events, %d with an unambiguous side mapping (strict majority, at least %d events); decisions for the unmatched pairs of linked disputes: flip %d events, keep %d; unmatched pairs kept as Side A because the dispute has no pair match %d events (%d disputes), no majority link / fewer than %d events %d events (%d disputes), split side mapping %d events (%d disputes)",
     sum(dl$n_pair_decided > 0), sum(dl$link_ok), WAR_DIR_MIN_LINK, sum(dl$map_ok), WAR_DIR_MIN_LINK,
     sum(events$dir_rule == "dispute_mapping" & events$dir_flipped), sum(events$dir_rule == "dispute_mapping" & !events$dir_flipped),
     sum(events$dir_rule == "dispute_no_pair_match"), n_distinct(events$dispnum[events$dir_rule == "dispute_no_pair_match"]),
     WAR_DIR_MIN_LINK, sum(events$dir_rule == "dispute_no_majority_link"), n_distinct(events$dispnum[events$dir_rule == "dispute_no_majority_link"]),
     sum(events$dir_rule == "dispute_split_side_mapping"), n_distinct(events$dispnum[events$dir_rule == "dispute_split_side_mapping"]))
## what a minimum of 1 would do (for the PI): disputes with exactly one pair-decided event
one_link <- dl |> filter(n_pair_decided == 1L)
blog("direction (D1): with WAR_DIR_MIN_LINK = 1 the %d disputes with a single pair-decided event would map their %d unmatched pairs too (flip %d, keep %d); not applied",
     nrow(one_link), sum(events$dispnum %in% one_link$dispnum & events$dir_rule == "dispute_no_majority_link"),
     sum(events$dispnum %in% one_link$dispnum[!one_link$sideA_war_side_init %in% TRUE] & events$dir_rule == "dispute_no_majority_link"),
     sum(events$dispnum %in% one_link$dispnum[one_link$sideA_war_side_init %in% TRUE] & events$dir_rule == "dispute_no_majority_link"))
blog("direction (D1): events by rule: %s", paste(names(table(events$dir_rule)), table(events$dir_rule), sep = "=", collapse = ", "))
if (DIRECTION == "sideA") {
  blog("DIRECTION = sideA: the %d flips above are NOT applied (appendix rebuild); dir_source = sideA, dir_rule = sideA_build on every event", sum(events$dir_flipped))
  events <- events |> mutate(dir_flipped = FALSE, dir_source = "sideA", dir_rule = "sideA_build", war_link = NA_integer_)
}
## the full list of flipped events, grouped by dispute, with the COW war number and name
flipped <- events |> filter(dir_flipped) |> arrange(dispnum, year, initiator, target)
blog("direction (D1): %d events reversed in %d disputes (%s); listed below as dispute / war / year / Side A rule -> applied direction",
     nrow(flipped), n_distinct(flipped$dispnum), paste(names(table(flipped$dir_source)), table(flipped$dir_source), sep = "=", collapse = ", "))
for (d in unique(flipped$dispnum)) {
  fd <- flipped[flipped$dispnum == d, ]
  cat(sprintf("  dispute %d (COW war %d, %s; %d of %d war-level pairs reversed, %d by the pair match, %d by the dispute mapping):\n", d, fd$war_link[1], cow_war_name(fd$war_link[1]),
              nrow(fd), sum(events$dispnum == d), sum(fd$dir_source == "cow_war_pair"), sum(fd$dir_source == "cow_war_dispute")), file = build_log, append = TRUE)
  for (r in seq_len(nrow(fd))) cat(sprintf("    %d: Side A %d -> %d  becomes  %d -> %d  (%s)\n", fd$year[r], fd$initiator[r], fd$target[r], fd$target[r], fd$initiator[r], fd$dir_source[r]), file = build_log, append = TRUE)
}
amb <- events |> filter(hostlev_disp == 5L, dir_source == "sideA") |> count(dispnum, dir_rule)
blog("direction (D1): war-level dispute events kept as Side A: %d in %d disputes (listed below with the reason)", sum(amb$n), n_distinct(amb$dispnum))
for (r in seq_len(nrow(amb))) cat(sprintf("  kept as Side A: dispute %d, %d events, reason %s%s\n", amb$dispnum[r], amb$n[r], amb$dir_rule[r],
                                          if (amb$dispnum[r] %in% dl$dispnum) sprintf(" (pair-decided events %d, linked war %s, side mapping %d of %d)", dl$n_pair_decided[dl$dispnum == amb$dispnum[r]],
                                                                                        ifelse(is.na(dl$war_link[dl$dispnum == amb$dispnum[r]]), "none", as.character(dl$war_link[dl$dispnum == amb$dispnum[r]])),
                                                                                        dl$n_map_modal[dl$dispnum == amb$dispnum[r]], dl$n_map[dl$dispnum == amb$dispnum[r]]) else ""), file = build_log, append = TRUE)
events <- apply_direction(events)
stopifnot(nrow(events) == n_ev_dir, all((events$sidea1 == 0L) == events$dir_flipped), all(events$sidea1 + events$sidea2 == 1L))
if (DIRECTION == "cowwar") {
  ## the review's headline case must hold on the real data: Germany -> Poland 1939 exists, Poland -> Germany 1939 does not
  stopifnot(any(events$initiator == 255L & events$target == 290L & events$year == 1939L), !any(events$initiator == 290L & events$target == 255L & events$year == 1939L))
} else stopifnot(!any(events$dir_flipped))
rm(wl, match_in, one_link, flipped, amb)
events <- events |>
  mutate(is_joiner_pair = orig1 == 0L | orig2 == 0L,          # tabulation only; NOT a model covariate (originator status moved with the state)
         fatal1 = ifelse(fatality1 < 0, NA, fatality1 >= 1),
         fatal2 = ifelse(fatality2 < 0, NA, fatality2 >= 1),
         fatal_disp = ifelse(fatality_disp < 0, NA, fatality_disp >= 1),
         fatal = case_when(fatal1 %in% TRUE | fatal2 %in% TRUE | fatal_disp %in% TRUE ~ TRUE,
                           is.na(fatal1) & is.na(fatal2) & is.na(fatal_disp) ~ NA,
                           TRUE ~ FALSE),
         ## war at the PAIR level (PI decision 2026-10-24, revised the same day): the dispute is war-level AND at least one
         ## participant of the pair reached hostility level 5 (participant hostlev from cow_mid_dirdisps, max over the pair's
         ## years). The first version required both participants at 5 (329 war events); the PI chose the less strict rule.
         ## war_disp is the old dispute-level coding, kept for the counts only.
         war_disp = hostlev_disp == 5L,
         war = hostlev_disp == 5L & (hostlev1 == 5L | hostlev2 == 5L))
stopifnot(all(events$initiator != events$target), !anyNA(events$war), !anyNA(events$war_disp), all(events$war_disp[events$war]))
blog("war coding: dispute-level war events %d; pair-level (dispute war AND at least one participant hostlev 5) %d; dropped by the pair rule %d (joiner pairs %d, originator pairs %d)",
     sum(events$war_disp), sum(events$war), sum(events$war_disp & !events$war), sum(events$war_disp & !events$war & events$is_joiner_pair),
     sum(events$war_disp & !events$war & !events$is_joiner_pair))
## ---- ongoing-dispute rule, month resolution (PI decision 2026-10-24; helpers and the date conventions in helpers_data.R) --
## An event is an escalation when i and j are both participants in a DIFFERENT dispute that started before the event's
## dispute start date and had not ended by it, at month resolution. peacesciencer 1.2.0 gives dispute-level start month
## and year (cow_mid_disps styear / stmon) and durations (mindur / maxdur) but no end month and no participant dates
## (cow_mid_dirdisps has years only), so the end date is bounded from the duration and every date is dispute-level
## (dispute_date_bounds, pair_dispute_intervals). Where the month comparison is undecided (unknown entry or exit month
## of a late joiner / early leaver, an end bound straddling the event month) the year rule decides and the fallback is
## counted. The year rule (already_fighting_keys / flag_escalations) is kept for the comparison below.
disp_dates <- dispute_date_bounds(disp_full <- cow_mid_disps |> mutate(dispnum = as.integer(dispnum)), dd)
blog("dispute dates (dispute-level; peacesciencer has no end month and no participant dates): %d disputes; start month known %d; end month exact %d, end spans 2 months %d, more %d; end bounds outside the last dispute-year (whole year used) %d",
     nrow(disp_dates), sum(disp_dates$start_month_known), sum(disp_dates$end_lo == disp_dates$end_hi), sum(disp_dates$end_hi - disp_dates$end_lo == 1L),
     sum(disp_dates$end_hi - disp_dates$end_lo > 1L), sum(disp_dates$end_rule == "year_only"))
## Participant-level dates (PI decision 2026-10-25, D5): peacesciencer::gml_part has participant start and end year / month for
## 2165 of the 2436 COW dispnums; where the dispnum matches, the pair's entry is the later of the two starts and its exit the earlier
## of the two ends (exact months), else the dispute-level bounds stay. The pre-D5 intervals (dispute bounds only, unclipped entry)
## are rebuilt for the before / after counts only; every downstream object uses `pdi`.
part_tbl <- if (PART_DATES) gml_part else NULL
pdi_nopart <- pair_dispute_intervals(dd, disp_dates, part = NULL)
pdi <- pair_dispute_intervals(dd, disp_dates, part = part_tbl)
stopifnot(nrow(pdi) == nrow(pdi_nopart), identical(pdi$dispnum, pdi_nopart$dispnum), identical(pdi$ccode1, pdi_nopart$ccode1))
blog("pair-dispute spells (both directions): %d (pairs %d, pairs with a re-entry spell %d); entry month exact %d; entry month unknown %d; spells ending before the dispute's last year (exit month unknown) %d",
     nrow(pdi), nrow(distinct(pdi, dispnum, ccode1, ccode2)), sum(pdi$spell > 1L), sum(pdi$start_exact), sum(!pdi$start_exact),
     sum(pdi$spell_last != disp_dates$y1[match(pdi$dispnum, disp_dates$dispnum)]))
blog("participant dates (D5): gml_part covers %d of %d COW dispnums (%d participant spells); pair spells dated from participant months %d of %d (exact entry AND exit); dispute-level bounds kept for %d; reasons: %s",
     length(intersect(unique(as.integer(gml_part$dispnum)), unique(dd$dispnum))), n_distinct(dd$dispnum), nrow(gml_part), sum(pdi$date_source == "participant"), nrow(pdi),
     sum(pdi$date_source == "dispute"), paste(names(table(pdi$part_reason)), table(pdi$part_reason), sep = "=", collapse = "; "))
blog("participant dates (D5): spells whose entry bound exceeded the feasible exit (start_hi > end_hi) under the dispute-level bounds: %d of %d, now clipped to the exit (start_hi <= end_hi asserted); with participant dates %d spells still needed the clip; exact entry months before / after %d / %d, exact exit months before / after %d / %d",
     sum(pdi_nopart$start_clipped), nrow(pdi_nopart), sum(pdi$start_clipped), sum(pdi_nopart$start_exact), sum(pdi$start_exact),
     sum(pdi_nopart$end_lo == pdi_nopart$end_hi), sum(pdi$end_lo == pdi$end_hi))
## Codex review P1-2 (2026-10-24): a pair with a joiner on either side whose spell begins in the dispute's first year used
## to get the dispute's exact start month although the joiner's entry month is unknown. Its entry interval is now
## [dispute start, December of that year] (start_exact FALSE). The pre-review intervals are rebuilt here ONLY to log the
## before / after counts; every downstream object uses `pdi`.
pdi_old <- pair_dispute_intervals(dd, disp_dates, joiner_entry_unknown = FALSE, part = part_tbl)
stopifnot(nrow(pdi_old) == nrow(pdi), identical(pdi_old$start_lo, pdi$start_lo), identical(pdi_old$end_lo, pdi$end_lo), identical(pdi_old$end_hi, pdi$end_hi))
chg_spell <- pdi$start_hi != pdi_old$start_hi
stopifnot(all(!pdi$both_orig[chg_spell]), all(pdi$date_source[chg_spell] == "dispute"), all(pdi$spell_first[chg_spell] == disp_dates$y0[match(pdi$dispnum[chg_spell], disp_dates$dispnum)]),
          all(ym_year(pdi$start_hi[chg_spell]) == pdi$spell_first[chg_spell]), all(pdi$start_hi[chg_spell] %% 12L == 11L | pdi$start_clipped[chg_spell]))
blog("joiner entry month (review P1-2): %d of %d spells have a joiner on either side; %d of them begin in the dispute's first year and now carry the interval [dispute start, December] instead of the exact start month (directed, both directions; %d distinct undirected pairs); spells with an exact entry month before / after: %d / %d",
     sum(!pdi$both_orig), nrow(pdi), sum(chg_spell), nrow(pdi[chg_spell, ] |> transmute(a = pmin(ccode1, ccode2), b = pmax(ccode1, ccode2), dispnum) |> distinct()),
     sum(pdi_old$start_exact), sum(pdi$start_exact))
af_keys <- already_fighting_keys(dd)                           # the year rule, for the comparison
assert_no_dup_keys(af_keys, c("ccode1", "ccode2", "year", "dispnum"), "already-fighting keys (year rule)")
ev_year <- flag_escalations(events, af_keys)                   # year rule on the same events
ev_old <- flag_escalations_month(events, pdi_old)              # pre-review joiner intervals, for the comparison only
ev_nopart <- flag_escalations_month(events, pdi_nopart)        # pre-D5 intervals (dispute-level bounds), for the comparison only
events <- flag_escalations_month(events, pdi)
cmp <- attr(events, "comparisons"); cmp_nopart <- attr(ev_nopart, "comparisons")
stopifnot(identical(events$escalation_year, ev_year$escalation))   # the year rule inside the month helper equals the old path
stopifnot(nrow(ev_old) == nrow(events), identical(ev_old$dispnum, events$dispnum), identical(ev_old$initiator, events$initiator),
          nrow(ev_nopart) == nrow(events), identical(ev_nopart$dispnum, events$dispnum), identical(ev_nopart$escalation_year, events$escalation_year))
blog("participant dates (D5): escalations with dispute-level bounds / with participant dates %d / %d (newly escalation %d, no longer escalation %d); undecided (event, other dispute) comparisons %d / %d; events whose decision rests on a year-rule fallback %d / %d",
     sum(ev_nopart$escalation), sum(events$escalation), sum(!ev_nopart$escalation & events$escalation), sum(ev_nopart$escalation & !events$escalation),
     sum(cmp_nopart$fallback), sum(cmp$fallback), sum(ev_nopart$esc_fallback), sum(events$esc_fallback))
ex_p <- events[ev_nopart$escalation != events$escalation, ] |> head(10)
for (r in seq_len(nrow(ex_p))) cat(sprintf("  D5 change: %d -> %d %d dispute %d: escalation %s -> %s (event date %s..%s)\n", ex_p$initiator[r], ex_p$target[r], ex_p$year[r], ex_p$dispnum[r],
                                           ev_nopart$escalation[ev_nopart$dispnum == ex_p$dispnum[r] & ev_nopart$initiator == ex_p$initiator[r] & ev_nopart$target == ex_p$target[r]], ex_p$escalation[r],
                                           ym_label(ex_p$ev_start_lo[r]), ym_label(ex_p$ev_start_hi[r])), file = build_log, append = TRUE)
## same-year opposite-direction pairs among the non-escalation events (i -> j and j -> i both in year t, from different disputes)
opp_pairs <- function(ev) { k <- ev |> filter(!escalation) |> transmute(a = pmin(initiator, target), b = pmax(initiator, target), year, dir = initiator < target) |> distinct() |>
  count(a, b, year) |> filter(n == 2L); nrow(k) }
blog("participant dates (D5): same-year opposite-direction pairs among the non-escalation events (unordered (i, j, t) with both i -> j and j -> i): dispute-level bounds %d, participant dates %d",
     opp_pairs(ev_nopart), opp_pairs(events))
blog("joiner entry month (review P1-2): escalations before / after %d / %d (newly escalation %d, no longer escalation %d); events whose decision rests on a year-rule fallback before / after %d / %d; events whose own date interval widened to [start, December] %d",
     sum(ev_old$escalation), sum(events$escalation), sum(!ev_old$escalation & events$escalation), sum(ev_old$escalation & !events$escalation),
     sum(ev_old$esc_fallback), sum(events$esc_fallback), sum(ev_old$ev_start_hi != events$ev_start_hi))
ex_j <- events[ev_old$escalation != events$escalation, ] |> head(10)
for (r in seq_len(nrow(ex_j))) cat(sprintf("  joiner-rule change: %d -> %d %d dispute %d (joiner pair %s): escalation %s -> %s\n", ex_j$initiator[r], ex_j$target[r], ex_j$year[r], ex_j$dispnum[r],
                                           ex_j$is_joiner_pair[r], ev_old$escalation[ev_old$dispnum == ex_j$dispnum[r] & ev_old$initiator == ex_j$initiator[r] & ev_old$target == ex_j$target[r]], ex_j$escalation[r]), file = build_log, append = TRUE)
rm(pdi_old, ev_old)
blog("escalations: year rule %d; month rule %d (year TRUE & month FALSE %d; year FALSE & month TRUE %d); (event, other dispute) comparisons %d, undecided at month resolution (year fallback) %d; events whose decision rests on a fallback %d (of which escalation %d)",
     sum(events$escalation_year), sum(events$escalation), sum(events$escalation_year & !events$escalation), sum(!events$escalation_year & events$escalation),
     nrow(cmp), sum(cmp$fallback), sum(events$esc_fallback), sum(events$esc_fallback & events$escalation))
fb_reason <- cmp |> filter(fallback) |>
  mutate(reason = case_when(ev_start_lo != ev_start_hi ~ "event month unknown (pair entered its dispute after the dispute's start year)",
                            o_start_lo != o_start_hi ~ "other dispute's pair entry month unknown (later joiner)",
                            TRUE ~ "other dispute's end month ambiguous (duration bounds or exit month unknown)")) |> count(reason)
for (r in seq_len(nrow(fb_reason))) blog("  fallback reason: %s = %d comparisons", fb_reason$reason[r], fb_reason$n[r])
## examples for the PI (China -> USA 1950 first; then the first changes in each direction)
show_ev <- function(e) sprintf("%d -> %d %d dispute %d (start %s): year rule %s, month rule %s%s", e$initiator, e$target, e$year, e$dispnum,
                               ym_label(e$ev_start_lo), e$escalation_year, e$escalation, if (e$esc_fallback) " (fallback)" else "")
ex_cn <- events |> filter(initiator == 710L, target == 2L, year == 1950L)
stopifnot(nrow(ex_cn) == 1L, ex_cn$escalation_year, !ex_cn$escalation)   # the audit's example must flip
cat("  escalation example (audit): ", show_ev(ex_cn), "\n", file = build_log, append = TRUE)
pp <- pdi |> semi_join(cmp |> filter(initiator == 710L, target == 2L, year == 1950L, year_dec) |> transmute(dispnum = other_dispnum), by = "dispnum") |> filter(ccode1 == 710L, ccode2 == 2L)
for (r in seq_len(nrow(pp)))
  cat(sprintf("    other dispute %d between 710 and 2 (spell %d): pair years %d-%d, entry %s..%s, exit between %s and %s\n", pp$dispnum[r], pp$spell[r], pp$spell_first[r], pp$spell_last[r],
              ym_label(pp$start_lo[r]), ym_label(pp$start_hi[r]), ym_label(pp$end_lo[r]), ym_label(pp$end_hi[r])), file = build_log, append = TRUE)
ex1 <- events |> filter(escalation_year, !escalation) |> head(5); ex2 <- events |> filter(!escalation_year, escalation) |> head(5)
for (r in seq_len(nrow(ex1))) cat("  no longer an escalation: ", show_ev(ex1[r, ]), "\n", file = build_log, append = TRUE)
for (r in seq_len(nrow(ex2))) cat("  newly an escalation:     ", show_ev(ex2[r, ]), "\n", file = build_log, append = TRUE)
## ---- choice-set rules: stage-1 IV set at January, stage-2 set at the reference date (Codex review P1-1, PI choice "A'") --
## Two already-fighting rules give two sets per state-year (i, t), both at month resolution with the year fallback.
##   stage-1 IV set : alternative j leaves when i and j are in a dispute active at JANUARY of t that started before it, for
##                    EVERY state-year. The same rule for initiating and non-initiating state-years, so IV[n] cannot
##                    depend on whether (i, t) initiated (the previous rule used the event date for initiators, which leaked
##                    y1 into the stage-1 regressor).
##   stage-2 set    : reference date of (i, t) = start of its earliest non-escalation event when it initiates, else January
##                    of t; j leaves when i and j are in a different dispute active at that date that started before it.
##                    Unchanged from the previous build. Only information up to the reference date enters.
## For a non-initiating state-year the two sets are identical. The alternatives files keep ONE row block per state-year,
## the UNION of both sets, with the 0/1 flags in_iv and in_s2 (choice_set_flags() in helpers_data.R); the Stan programs take
## IV[n] over the in_iv rows and the stage-2 softmax over the in_s2 rows. The year rule's keys are kept to count the change.
## Because the stage-2 set is fixed at the reference date, a later event of the same state-year whose pair was still
## fighting at the reference date (their dispute ended between the reference date and the event) has no stage-2 row (it
## may still have an IV row with in_s2 = 0). Such events are flagged fighting_at_ref, stay in every sample flag for y1 and
## leave stage 2 only (logged below); the build stops if any other non-escalation event loses its stage-2 row.
ref_dates <- reference_dates(events)
af_s2_keys <- already_fighting_keys_month(pdi, ref_dates)      # stage-2 set: reference date of the state-year
af_iv_keys <- already_fighting_keys_month(pdi, NULL)           # stage-1 IV set: January of t for every state-year
af_s2 <- af_s2_keys |> filter(drop) |> distinct(ccode1, ccode2, year)
af_iv <- af_iv_keys |> filter(drop) |> distinct(ccode1, ccode2, year)
af_dir_year <- af_keys |> distinct(ccode1, ccode2, year)
stopifnot(identical(nrow(af_dir_year), nrow(af_s2_keys |> filter(drop_year) |> distinct(ccode1, ccode2, year))),
          nrow(af_iv_keys) == nrow(af_s2_keys), all(af_iv_keys$ref_lo == ym_index(af_iv_keys$year, 1L)), all(af_iv_keys$ref_hi == af_iv_keys$ref_lo))
## the two rules can differ only for initiating state-years (ref row present)
iv_not_s2 <- anti_join(af_iv, af_s2, by = c("ccode1", "ccode2", "year")); s2_not_iv <- anti_join(af_s2, af_iv, by = c("ccode1", "ccode2", "year"))
stopifnot(nrow(iv_not_s2 |> anti_join(ref_dates, by = c("ccode1" = "ccode", "year"))) == 0,
          nrow(s2_not_iv |> anti_join(ref_dates, by = c("ccode1" = "ccode", "year"))) == 0)
blog("already-fighting directed pair-years, stage-2 set (reference date): year rule %d distinct (i, j, t); month rule %d (in both %d; year only %d; month only %d); key-dispute comparisons %d, undecided (year fallback) %d, keys whose decision rests on a fallback %d",
     nrow(af_dir_year), nrow(af_s2), nrow(semi_join(af_s2, af_dir_year, by = c("ccode1", "ccode2", "year"))),
     nrow(anti_join(af_dir_year, af_s2, by = c("ccode1", "ccode2", "year"))), nrow(anti_join(af_s2, af_dir_year, by = c("ccode1", "ccode2", "year"))),
     nrow(af_s2_keys), sum(af_s2_keys$fallback),
     nrow(af_s2_keys |> group_by(ccode1, ccode2, year) |> summarise(fb = any(fallback) & !any(!fallback & drop), .groups = "drop") |> filter(fb)))
blog("already-fighting directed pair-years, stage-1 IV set (January of t for every state-year): %d distinct (i, j, t); dropped from the IV set but not from the stage-2 set %d, from the stage-2 set but not from the IV set %d (all in initiating state-years; %d of %d initiating state-years have different sets); undecided key comparisons (year fallback) %d",
     nrow(af_iv), nrow(iv_not_s2), nrow(s2_not_iv), nrow(bind_rows(iv_not_s2, s2_not_iv) |> distinct(ccode1, year)), nrow(ref_dates), sum(af_iv_keys$fallback))
events <- events |> left_join(af_s2 |> transmute(initiator = ccode1, target = ccode2, year, fighting_at_ref = TRUE), by = c("initiator", "target", "year")) |>
  mutate(fighting_at_ref = coalesce(fighting_at_ref, FALSE) & !escalation)
far <- events |> filter(fighting_at_ref)
if (nrow(far) > 0) {
  ## each must be a later event of a state-year whose reference date precedes it, with a different dispute of the pair
  ## that was active at the reference date; a state-year's earliest event can never be one (it would be an escalation)
  chk <- far |> inner_join(ref_dates, by = c("initiator" = "ccode", "year")) |>
    inner_join(af_s2_keys |> filter(drop) |> select(initiator = ccode1, target = ccode2, year, od = dispnum), by = c("initiator", "target", "year"), relationship = "many-to-many")
  stopifnot(nrow(chk) >= nrow(far), all(chk$od != chk$dispnum), all(chk$ev_start_hi >= chk$ref_lo), !any(chk$ev_start_hi == chk$ref_hi & chk$ev_start_lo == chk$ref_lo & chk$ev_start_lo == chk$ev_start_hi))
  blog("non-escalation events whose pair was already fighting at the state-year's reference date (fighting_at_ref; stage 2 only, y1 unchanged): %d (joiner pairs %d; war-level %d)",
       nrow(far), sum(far$is_joiner_pair), sum(far$war))
  for (r in seq_len(nrow(far))) cat(sprintf("  fighting_at_ref: %d -> %d %d dispute %d (start %s; reference date %s..%s)\n", far$initiator[r], far$target[r], far$year[r], far$dispnum[r],
                                           ym_label(far$ev_start_lo[r]), ym_label(ref_dates$ref_lo[ref_dates$ccode == far$initiator[r] & ref_dates$year == far$year[r]]),
                                           ym_label(ref_dates$ref_hi[ref_dates$ccode == far$initiator[r] & ref_dates$year == far$year[r]])), file = build_log, append = TRUE)
}
## D5 before / after on the stage-2 choice-set rule: fighting_at_ref under the dispute-level bounds
ref_np <- reference_dates(ev_nopart)
af_np <- already_fighting_keys_month(pdi_nopart, ref_np) |> filter(drop) |> distinct(ccode1, ccode2, year)
far_np <- ev_nopart |> filter(!escalation) |> semi_join(af_np |> rename(initiator = ccode1, target = ccode2), by = c("initiator", "target", "year"))
blog("participant dates (D5): fighting_at_ref events with dispute-level bounds / with participant dates %d / %d; already-fighting (i, j, t) keys at the reference date %d / %d, at January %d / %d",
     nrow(far_np), nrow(far), nrow(af_np), nrow(af_s2), nrow(already_fighting_keys_month(pdi_nopart, NULL) |> filter(drop) |> distinct(ccode1, ccode2, year)), nrow(af_iv))
rm(pdi_nopart, ev_nopart, cmp_nopart, ref_np, af_np, far_np, ex_p)
## ---- WWII disputes (switch --ww2=drop; rule and list in helpers_data.R) ------------------------------------
ww2_tbl <- ww2_disputes(disp_full, dd, cow_war_inter)
events <- events |> mutate(is_ww2 = dispnum %in% WW2_DISPNUMS)
blog("WWII disputes (is_ww2; appendix option --drop_ww2): %s", paste(sort(WW2_DISPNUMS), collapse = ", "))
for (r in seq_len(nrow(ww2_tbl)))
  cat(sprintf("  WWII candidate %5d: %d-%02d, years %d-%d, %2d states, hiact %d, dyad-years in COW war 139 %3d -> %s (%s)\n",
              ww2_tbl$dispnum[r], ww2_tbl$styear[r], ww2_tbl$stmon[r], ww2_tbl$y0[r], ww2_tbl$y1[r], ww2_tbl$n_states[r],
              ww2_tbl$hiact[r], ww2_tbl$n_dyadyears_in_war139[r], if (ww2_tbl$is_ww2[r]) "WWII" else "not WWII", ww2_tbl$reason[r]),
      file = build_log, append = TRUE)
## fatal-MID coding (PI decision 2026-10-23). Main fmid: wars are fatal by definition (COW wars need
## 1,000+ battle deaths; since 2026-10-24 "war" is the PAIR-level coding above, so a peripheral participant of a
## war-level dispute that did not itself reach hostility level 5 is not made fatal by this rule); events whose COW
## fatality codes are all missing are filled from MIC 1.0
## (Gibler et al., confrontation-level fatalmin >= 1; micnum = COW dispnum). MIC never overrides a
## COW code. fmidcow = the original COW-only coding (missing -> not fatal), kept for the appendix.
mic_path <- file.path(v2_root(), "data_raw", "mic-conf-1.0.csv")
stopifnot(file.exists(mic_path))
if (file_digest(mic_path) != MIC_CONF_SHA256) stop("data_raw/mic-conf-1.0.csv sha256 does not match MIC_CONF_SHA256 in helpers_data.R")
mic <- read.csv(mic_path) |> transmute(dispnum = as.integer(micnum), mic_fatalmin = as.integer(fatalmin))
assert_no_dup_keys(mic, "dispnum", "MIC confrontations"); assert_no_na(mic, "mic_fatalmin", "MIC confrontations")
stopifnot(nrow(mic) == 1881L, all(mic$mic_fatalmin >= 0))
## GML MID (Gibler, Miller and Little 2016; peacesciencer gml_mid_ddydisps, v2.2.1) is the third source, used only
## when COW and MIC are both missing. Same rule as COW: per-event max over dispute-years of each side's fatality
## and the dispute fatality; negative codes -> NA; fatal if any >= 1; NA if all NA. Joined on (dispnum, initiator, target).
stopifnot(exists("gml_mid_ddydisps"))
gml <- gml_mid_ddydisps |> mutate(across(c(dispnum, ccode1, ccode2), as.integer)) |>
  group_by(dispnum, initiator = ccode1, target = ccode2) |>
  summarise(g1 = max(fatality1), g2 = max(fatality2), gd = max(fatality), .groups = "drop") |>
  mutate(g1 = ifelse(g1 < 0, NA, g1 >= 1), g2 = ifelse(g2 < 0, NA, g2 >= 1), gd = ifelse(gd < 0, NA, gd >= 1),
         fatal_gml = case_when(g1 %in% TRUE | g2 %in% TRUE | gd %in% TRUE ~ TRUE,
                               is.na(g1) & is.na(g2) & is.na(gd) ~ NA,
                               TRUE ~ FALSE)) |>
  select(dispnum, initiator, target, fatal_gml)
assert_no_dup_keys(gml, c("dispnum", "initiator", "target"), "GML directed pairs")
n_ev_before <- nrow(events)
## fatal_main / fatal_cow / fatal_ub / war are the fatality codings of the event (escalations included, for the counts
## below); the sample flags in_* are these codings AND NOT escalation (PI decision 2026-10-24).
events <- events |> left_join(mic, by = "dispnum") |> left_join(gml, by = c("dispnum", "initiator", "target")) |>
  mutate(fatal_mic = ifelse(is.na(mic_fatalmin), NA, mic_fatalmin >= 1),
         fatal_main = code_fatal_main(fatal, fatal_mic, war, fatal_gml),
         fatal_cow = fatal %in% TRUE,
         fatal_ub = fatal_main | (is.na(fatal) & is.na(fatal_mic) & is.na(fatal_gml)),   # upper bound: uncoded -> fatal
         fatal_source = case_when(war %in% TRUE & !(fatal %in% TRUE) ~ "war_rule",
                                  !is.na(fatal) ~ "cow",
                                  !is.na(fatal_mic) ~ "mic_fill",
                                  !is.na(fatal_gml) ~ "gml_fill",
                                  TRUE ~ "missing_nonfatal"),
         in_mid = !escalation,
         in_fmid = fatal_main & !escalation,
         in_fmidcow = fatal_cow & !escalation,
         in_fmidub = fatal_ub & !escalation,
         in_war = (war %in% TRUE) & !escalation,
         ## no-WWII sample flags (Codex review P2-3): the sample flag of this dispute pair when it is not a WWII dispute. They
         ## are ORed in the duplicate collapse like the other flags, so an (i, j, t) that mixes a WWII and a non-WWII
         ## dispute stays in the --drop_ww2 sample only when a non-WWII dispute qualifies (its fatal / war coding comes
         ## from the non-WWII disputes alone). assemble_cell(DROP_WW2 = TRUE) uses these flags.
         in_mid_noww2 = in_mid & !is_ww2, in_fmid_noww2 = in_fmid & !is_ww2, in_fmidcow_noww2 = in_fmidcow & !is_ww2,
         in_fmidub_noww2 = in_fmidub & !is_ww2, in_war_noww2 = in_war & !is_ww2)
stopifnot(nrow(events) == n_ev_before, !anyNA(events$in_fmid), !anyNA(events$in_mid), !anyNA(events$in_war),
          all(events$in_mid_noww2 == (events$in_mid & !events$is_ww2)), all(events$in_war_noww2 <= events$in_war),
          all(events$in_fmid[events$in_war]),          # war is a subset of fmid
          all(events$in_fmid[events$in_fmidcow]),      # COW-only coding is a subset of the main coding
          all(events$in_fmidub[events$in_fmid]),       # main coding is a subset of the upper bound
          all(events$in_mid[events$in_fmidub]),        # every sample is a subset of mid (so escalations are in none)
          !any(events$in_mid & events$escalation),
          sum(events$in_fmidub & !events$in_fmid) == sum(events$fatal_source == "missing_nonfatal" & !events$escalation))
## ---- stage-2 event weights (PI decision 2026-10-25, D8; README "Coalition expansion") ---------------------------------------
## w_<sample> = 1 / (number of the initiator's targets in the same dispute and onset year that are in the sample), computed on the
## pre-collapse rows; the collapse sums a key's weights and caps at 1. Escalations carry 0 (they are in no sample).
events <- event_weights(events)
stopifnot(all(events$w_mid[events$in_mid] > 0), all(events$w_mid[!events$in_mid] == 0), all(events$w_war[!events$in_war] == 0))
for (s in c("mid", "fmid", "fmidcow", "fmidub", "war")) {
  w <- events[[paste0("w_", s)]][events[[paste0("in_", s)]]]; k <- round(1 / w)
  blog("  event weights %-7s: %d events; initiator-dispute-years %d; weight 1 (single target) %d events, 1/2 %d, 1/3..1/5 %d, 1/6..1/10 %d, below 1/10 %d; min weight %.4f (k = %d targets); sum of weights %.1f (= decisions)",
       s, length(w), nrow(events |> filter(.data[[paste0("in_", s)]]) |> distinct(dispnum, initiator, year)), sum(k == 1L), sum(k == 2L), sum(k >= 3L & k <= 5L), sum(k >= 6L & k <= 10L), sum(k > 10L),
       min(w), max(k), sum(w))
}
kfor <- events |> filter(dispnum %in% c(4186L, 4343L), year == 2000L) |> group_by(dispnum, initiator) |> summarise(n_targets = n(), w_mid = w_mid[1], .groups = "drop") |> arrange(dispnum, desc(n_targets))
blog("  event weights, KFOR / Yugoslavia 2000 (disputes 4186 and 4343): %s",
     paste(sprintf("dispute %d: %d initiators, max targets per initiator %d (initiator %d, w_mid = %.4f)", unique(kfor$dispnum),
                   tapply(kfor$initiator, kfor$dispnum, length), tapply(kfor$n_targets, kfor$dispnum, max),
                   kfor$initiator[match(unique(kfor$dispnum), kfor$dispnum)], kfor$w_mid[match(unique(kfor$dispnum), kfor$dispnum)]), collapse = "; "))
rm(kfor)
## ---- coalition context (dual typology stage a, PI-approved spec 2026-10-25; the `context` switch; helpers in helpers_data.R) -----------
## Per D8 weight group (initiator, dispute, onset year): is the dispute multilateral at the GROUP onset (the initiator's earliest pair
## onset in that dispute-year)? Participant entry / exit dates from gml_part at day resolution where the day is known, month otherwise,
## dispute-level bounds where GML and COW disagree; comparisons at the coarser resolution of the two dates. Both tie conventions are
## stored: coal_incl (same-date entrants present) and coal_strict (strictly earlier entrants only). Computed on the pre-collapse rows and
## joined to the collapsed events below (ORed over a key's source disputes).
psd <- participant_spell_dates(dd, disp_dates, part = part_tbl)
pairs_pre <- events |> filter(!escalation) |> select(dispnum, initiator, target, year, ev_start_lo, ev_start_hi)
coal <- coalition_context(pairs_pre, psd)
coal_grp <- attr(coal, "groups"); coal_ev <- attr(coal, "events"); coal_mem <- attr(psd, "membership")
stopifnot(nrow(coal) == nrow(pairs_pre), identical(coal$dispnum, pairs_pre$dispnum), identical(coal$initiator, pairs_pre$initiator), identical(coal$target, pairs_pre$target),
          !anyNA(coal$coal_incl), !anyNA(coal$coal_strict), all(coal$coal_incl >= coal$coal_strict), all(coal_grp$coal_incl >= coal_grp$coal_strict))
## the D8 groups are exactly the coalition groups, and the flag is constant within each (asserted again on the joined rows)
w_grp <- events |> filter(!escalation) |> distinct(dispnum, initiator, year)
stopifnot(nrow(w_grp) == nrow(coal_grp), nrow(anti_join(w_grp, coal_grp, by = c("dispnum", "initiator", "year"))) == 0)
chk_c <- coal |> group_by(dispnum, initiator, year) |> summarise(a = n_distinct(coal_incl), b = n_distinct(coal_strict), .groups = "drop")
stopifnot(all(chk_c$a == 1L), all(chk_c$b == 1L))
blog("coalition context: %d pre-collapse non-escalation pairs in %d (initiator, dispute, year) groups (= the D8 weight groups); groups with a coalition context: inclusive %d, strict %d; pairs: inclusive %d, strict %d; the flag is constant within every group (asserted)",
     nrow(coal), nrow(coal_grp), sum(coal_grp$coal_incl), sum(coal_grp$coal_strict), sum(coal$coal_incl), sum(coal$coal_strict))
blog("coalition context: groups with two or more onset-defining targets (several targets faced from the first moment: multilateral under BOTH conventions, audit fix 19): %d, of which %d would be strict-bilateral by the entry comparison alone",
     sum(coal_grp$n_onset_defining >= 2L), sum(coal_grp$n_onset_defining >= 2L & (coal_grp$n_strict - as.integer(coal_grp$target_strict)) < 1L))
blog("coalition context: participant spells %d (entry dated from gml_part %d, of which day resolution %d, month %d; dispute-level fallback %d: %s); exit from gml_part %d (day %d, month %d), fallback %d",
     nrow(psd), sum(startsWith(psd$entry_src, "participant")), sum(psd$entry_res == 1L), sum(psd$entry_res == 2L & startsWith(psd$entry_src, "participant")), sum(!startsWith(psd$entry_src, "participant")),
     paste(names(table(psd$entry_src[!startsWith(psd$entry_src, "participant")])), table(psd$entry_src[!startsWith(psd$entry_src, "participant")]), sep = "=", collapse = "; "),
     sum(startsWith(psd$exit_src, "participant")), sum(psd$exit_res == 1L), sum(psd$exit_res == 2L & startsWith(psd$exit_src, "participant")), sum(!startsWith(psd$exit_src, "participant")))
blog("coalition context: resolution of the entry comparisons per group: day %d, month %d, year %d, mixed (day and month, or coarser, within one group) %d; group onset itself at day / month / year resolution %d / %d / %d",
     sum(coal_grp$coal_res == "day"), sum(coal_grp$coal_res == "month"), sum(coal_grp$coal_res == "year"), sum(coal_grp$coal_res == "mixed"),
     sum(coal_grp$onset_res == 1L), sum(coal_grp$onset_res == 2L), sum(coal_grp$onset_res == 3L))
blog("coalition context: COW and GML disagree on membership in %d disputes: %d COW participants absent from gml_part (dispute present in it), %d gml_part participants absent from COW; GML dates are used for the participants both sources have, COW membership decides who is in the dispute (listed below)",
     n_distinct(coal_mem$dispnum), sum(coal_mem$issue == "COW participant absent from gml_part"), sum(coal_mem$issue == "gml_part participant absent from COW"))
for (d in sort(unique(coal_mem$dispnum)))
  cat(sprintf("  membership disagreement: dispute %d: %s\n", d, paste(sprintf("%d (%s)", coal_mem$ccode[coal_mem$dispnum == d], coal_mem$issue[coal_mem$dispnum == d]), collapse = "; ")), file = build_log, append = TRUE)
blog("coalition context: pairs whose participant-derived onset year differs from the COW pair onset year (both states in the dispute before they were opposed; onset ambiguous across a year boundary, COW onset year used at year resolution): %d in %d groups; participant-derived onset month against the pair's escalation-rule interval: %s; onset-defining targets whose GML exit precedes the initiator's GML entry in the same month (present by construction): %d",
     sum(coal$coal_onset_ambiguous), nrow(coal |> filter(coal_onset_ambiguous) |> distinct(dispnum, initiator, year)),
     paste(names(table(coal_ev$onset_vs_pdi)), table(coal_ev$onset_vs_pdi), sep = "=", collapse = "; "), sum(coal_grp$n_target_exit_before_onset))
for (r in which(coal$coal_onset_ambiguous)) cat(sprintf("  ambiguous onset: %d -> %d %d dispute %d (coal_incl %s, coal_strict %s)\n", coal$initiator[r], coal$target[r], coal$year[r], coal$dispnum[r], coal$coal_incl[r], coal$coal_strict[r]), file = build_log, append = TRUE)
## diagnostic for the PI (not used): evaluating every pair at its OWN onset instead of the group onset
blog("coalition context: DIAGNOSTIC, flag at every pair's own onset instead of the group onset: pairs coalition inclusive %d / strict %d (group-level %d / %d); groups whose pairs would carry different flags: inclusive %d, strict %d of %d (the group-level flag keeps a D8 group whole)",
     sum(coal_ev$ev_incl), sum(coal_ev$ev_strict), sum(coal$coal_incl), sum(coal$coal_strict), sum(coal_grp$n_event_flags_incl > 1L), sum(coal_grp$n_event_flags_strict > 1L), nrow(coal_grp))
het_g <- coal_grp |> filter(n_event_flags_incl > 1L | n_event_flags_strict > 1L) |> arrange(dispnum, year, initiator)
for (r in seq_len(nrow(het_g))) cat(sprintf("  per-event flags differ within group: dispute %d, initiator %d, %d (%d pairs; group flag inclusive %s / strict %s)\n", het_g$dispnum[r], het_g$initiator[r], het_g$year[r], het_g$n_events[r], het_g$coal_incl[r], het_g$coal_strict[r]), file = build_log, append = TRUE)
## the review's headline case under both conventions
gp39 <- coal |> filter(initiator == 255L, target == 290L, year == 1939L, dispnum == 258L); gg39 <- coal_grp |> filter(dispnum == 258L, initiator == 255L, year == 1939L)
if (DIRECTION == "cowwar") {
  stopifnot(nrow(gp39) == 1L, nrow(gg39) == 1L)
  blog("coalition context: Germany -> Poland 1939 (dispute 258): group onset %s (%s resolution; Germany's %d pairs of 1939), states present besides Germany: inclusive %d, strict %d -> coal_incl = %s, coal_strict = %s; at the pair's own onset the per-event flags would be inclusive %s / strict %s; the group's per-event flags differ: inclusive %s, strict %s",
       if (gg39$onset_res == 1L) as.character(as.Date(gg39$onset_dkey, origin = "1970-01-01")) else if (gg39$onset_res == 2L) ym_label(gg39$onset_mkey) else as.character(gg39$onset_y),
       c("day", "month", "year")[gg39$onset_res], gg39$n_events, gg39$n_incl, gg39$n_strict, gp39$coal_incl, gp39$coal_strict,
       coal_ev$ev_incl[coal_ev$initiator == 255L & coal_ev$target == 290L & coal_ev$year == 1939L & coal_ev$dispnum == 258L],
       coal_ev$ev_strict[coal_ev$initiator == 255L & coal_ev$target == 290L & coal_ev$year == 1939L & coal_ev$dispnum == 258L],
       gg39$n_event_flags_incl > 1L, gg39$n_event_flags_strict > 1L)
} else blog("coalition context: Germany -> Poland 1939 does not exist under DIRECTION = sideA (Poland -> Germany: coal_incl %s, coal_strict %s)",
            coal$coal_incl[coal$initiator == 290L & coal$target == 255L & coal$year == 1939L & coal$dispnum == 258L], coal$coal_strict[coal$initiator == 290L & coal$target == 255L & coal$year == 1939L & coal$dispnum == 258L])
for (s in c("mid", "fmid", "fmidcow", "fmidub", "war")) {
  f <- events[[paste0("in_", s)]][!events$escalation]
  blog("  coalition context, sample %-7s (pre-collapse): %4d events; coalition inclusive %4d (%.1f%%), strict %4d (%.1f%%); ambiguous onsets %d", s, sum(f),
       sum(coal$coal_incl[f]), 100 * mean(coal$coal_incl[f]), sum(coal$coal_strict[f]), 100 * mean(coal$coal_strict[f]), sum(coal$coal_onset_ambiguous[f]))
}
rm(w_grp, chk_c, het_g, gp39, gg39)
## target has positions at t-1 and is a system member at t: the conditions for its alternative row to exist in the
## "all" choice set (besides not already fighting); assemble_cell() asserts every such event is in its choice set
events <- events |> left_join(pos |> transmute(target = ccode, year = year + 1L, .pos = 1L), by = c("target", "year")) |>
  left_join(states |> transmute(target = ccode, year, .sys = 1L), by = c("target", "year")) |>
  mutate(target_has_pos_tm1 = !is.na(.pos) & !is.na(.sys)) |> select(-.pos, -.sys)
stopifnot(nrow(events) == n_ev_before)
cow_na <- is.na(events$fatal)
## fatality-coding diagnostics are on the CODINGS (fatal_main, fatal_cow; escalations included); the sample sizes
## after the escalation exclusion are logged in the counts table below
blog("fatal coding (joiners %s, escalations included): COW fatality missing for %d events; of these %d disputes in MIC (MIC fatal %d, MIC non-fatal %d), %d not in MIC",
     if (INCLUDE_JOINERS) "included" else "excluded", sum(cow_na), sum(cow_na & !is.na(events$fatal_mic)), sum(cow_na & events$fatal_mic %in% TRUE),
     sum(cow_na & events$fatal_mic %in% FALSE), sum(cow_na & is.na(events$fatal_mic)))
blog("fatal coding: wars added to fmid by the war rule %d; fatal_main %d vs fatal_cow (COW only) %d; fatal_source counts: %s",
     sum(events$fatal_source == "war_rule"), sum(events$fatal_main), sum(events$fatal_cow),
     paste(names(table(events$fatal_source)), table(events$fatal_source), sep = "=", collapse = ", "))
blog("fatal coding by source for joiner pairs only (%d events): %s; for originator pairs (%d events): %s",
     sum(events$is_joiner_pair), paste(names(table(events$fatal_source[events$is_joiner_pair])), table(events$fatal_source[events$is_joiner_pair]), sep = "=", collapse = ", "),
     sum(!events$is_joiner_pair), paste(names(table(events$fatal_source[!events$is_joiner_pair])), table(events$fatal_source[!events$is_joiner_pair]), sep = "=", collapse = ", "))
blog("fatal coding: COW and MIC both missing for %d non-war events; GML codes %d of them (GML fatal %d); %d stay uncoded (non-fatal)",
     sum(cow_na & is.na(events$fatal_mic) & !(events$war %in% TRUE)), sum(events$fatal_source == "gml_fill"),
     sum(events$fatal_source == "gml_fill" & events$fatal_main), sum(events$fatal_source == "missing_nonfatal"))
blog("fatal coding (info only, not used): events in MIC %d of %d; where both code fatality, COW and MIC disagree on %d events",
     sum(!is.na(events$fatal_mic)), nrow(events), sum(!is.na(events$fatal) & !is.na(events$fatal_mic) & events$fatal != events$fatal_mic))
## one direction per dispute pair: no (dispnum, j, i) when (dispnum, i, j) exists (holds with joiners too)
rev_key <- paste(events$dispnum, events$target, events$initiator)
stopifnot(!any(rev_key %in% paste(events$dispnum, events$initiator, events$target)))
assert_no_dup_keys(events, c("dispnum", "initiator", "target"), "events")
stopifnot(all(events$initiator != events$target), !anyNA(events$is_joiner_pair), !anyNA(events$escalation), !anyNA(events$is_ww2))
if (INCLUDE_JOINERS) stopifnot(sum(!events$is_joiner_pair) == 2770L)   # the originator-only event base of the builds before 2026-10-24
blog("events (Side A -> Side B, joiners %s): %d (originator pairs %d, joiner pairs %d); escalations %d (joiner pairs %d); disputes %d; multilateral pairs %d; WWII events %d (pair-level war %d, dispute-level war %d)",
     if (INCLUDE_JOINERS) "INCLUDED" else "excluded", nrow(events), sum(!events$is_joiner_pair), sum(events$is_joiner_pair),
     sum(events$escalation), sum(events$escalation & events$is_joiner_pair), n_distinct(events$dispnum), sum(duplicated(events$dispnum)),
     sum(events$is_ww2), sum(events$is_ww2 & events$war), sum(events$is_ww2 & events$war_disp))
## counts table per sample: (a) originator pairs only (the pre-2026-10-24 event base), (b) with joiners, (c) after
## removing escalations under the YEAR rule (the c137995 build), (d) after removing escalations under the month rule
## (= the sample flag before the duplicate collapse), and how many of the removed events are joiner pairs
sample_codings <- list(mid = rep(TRUE, nrow(events)), fmid = events$fatal_main, fmidcow = events$fatal_cow, fmidub = events$fatal_ub, war = events$war)
for (s in names(sample_codings)) {
  cs <- sample_codings[[s]]; flag <- events[[paste0("in_", s)]]
  stopifnot(identical(as.logical(flag), as.logical(cs & !events$escalation)))
  blog("  sample %-7s: (a) originator pairs %4d; (b) with joiners %4d; (c) year-rule escalations removed %4d; (d) month-rule escalations removed %4d; removed by the month rule %3d (joiner pairs %3d, originator pairs %3d); WWII events in (d) %3d; fighting_at_ref in (d) %3d; direction reversed in (d) %3d (pair match %3d, dispute mapping %3d)",
       s, sum(cs & !events$is_joiner_pair), sum(cs), sum(cs & !events$escalation_year), sum(flag), sum(cs & events$escalation),
       sum(cs & events$escalation & events$is_joiner_pair), sum(cs & events$escalation & !events$is_joiner_pair), sum(flag & events$is_ww2), sum(flag & events$fighting_at_ref),
       sum(flag & events$dir_flipped), sum(flag & events$dir_flipped & events$dir_source == "cow_war_pair"), sum(flag & events$dir_flipped & events$dir_source == "cow_war_dispute"))
}
blog("direction (D1): reversed events %d of %d (escalations among them %d); by source %s; by sample, see the table above",
     sum(events$dir_flipped), nrow(events), sum(events$dir_flipped & events$escalation), paste(names(table(events$dir_source[events$dir_flipped])), table(events$dir_source[events$dir_flipped]), sep = "=", collapse = ", "))
## war at the pair level versus the dispute level, after the month-rule escalation exclusion
blog("war sample: dispute-level coding would give %d events after escalations, the pair-level coding gives %d (removed %d, of which WWII %d); WWII share of the pair-level war sample %.1f%% (%d of %d); war rule in fatal coding now adds %d events",
     sum(events$war_disp & !events$escalation), sum(events$in_war), sum(events$war_disp & !events$war & !events$escalation),
     sum(events$war_disp & !events$war & !events$escalation & events$is_ww2), 100 * sum(events$in_war & events$is_ww2) / sum(events$in_war),
     sum(events$in_war & events$is_ww2), sum(events$in_war), sum(events$fatal_source == "war_rule"))
## ---- same-year duplicates (PI decision 2026-10-24): one event per (initiator, target, year) among the kept rows -----------
## Two disputes between the same pair with the same onset year, neither an escalation of the other, were two chosen
## alternatives on one row. They collapse into one event (flags ORed, dispnums listed); the likelihood then sees one chosen
## alternative per (i, j, t), asserted in assemble_cell(). Escalation rows are untouched.
n_pre <- nrow(events); pre_counts <- sapply(c("in_mid", "in_fmid", "in_fmidcow", "in_fmidub", "in_war"), function(cc) sum(events[[cc]]))
dup_keys <- events |> filter(!escalation) |> count(initiator, target, year) |> filter(n > 1)
w_cap <- events |> filter(!escalation) |> group_by(initiator, target, year) |> summarise(s = sum(w_mid), .groups = "drop") |> filter(s > 1 + 1e-12)
events <- collapse_duplicate_events(events)
blog("direction (D1) / weights (D8) in the collapse: collapsed keys whose source disputes disagree on the direction rule (one reversed, one not; sidea1 + sidea2 != 1) %d; keys whose summed mid weight exceeded 1 and was capped %d",
     sum(events$sidea1 + events$sidea2 != 1L), nrow(w_cap))
rm(w_cap)
post_counts <- sapply(c("in_mid", "in_fmid", "in_fmidcow", "in_fmidub", "in_war"), function(cc) sum(events[[cc]]))
stopifnot(nrow(events) == n_pre - sum(dup_keys$n) + nrow(dup_keys), all(events$n_disputes[events$escalation] == 1L),
          sum(events$n_disputes[!events$escalation]) == pre_counts["in_mid"], all(post_counts <= pre_counts),
          !anyNA(events$dispnums), all(events$n_disputes >= 1L))
assert_no_dup_keys(events |> filter(!escalation), c("initiator", "target", "year"), "non-escalation events after the collapse")
blog("duplicate collapse: %d non-escalation rows in %d (initiator, target, year) keys collapsed into %d events; events.rds %d -> %d rows; rows with more than one dispute %d (max %d)",
     sum(dup_keys$n), nrow(dup_keys), nrow(dup_keys), n_pre, nrow(events), sum(events$n_disputes > 1L), max(events$n_disputes))
for (s in c("mid", "fmid", "fmidcow", "fmidub", "war"))
  blog("  sample %-7s: events before the collapse %4d, after %4d (collapsed away %3d)", s, pre_counts[paste0("in_", s)], post_counts[paste0("in_", s)], pre_counts[paste0("in_", s)] - post_counts[paste0("in_", s)])
blog("is_joiner_pair after the collapse: %d of %d non-escalation events (an event is a joiner pair only when every collapsed dispute pair was one)",
     sum(events$is_joiner_pair & !events$escalation), sum(!events$escalation))
stopifnot(!anyNA(events$fighting_at_ref), !any(events$fighting_at_ref & events$escalation))
## ---- coalition context after the collapse: ORed over a key's source disputes (a key collapsed from a bilateral and a coalition dispute is a
## coalition event; such keys are the only way a D8 group can be split by the `context` switch, counted here and asserted at assembly) ----
## (.mixed_* come FIRST: inside summarise() a later expression would see the already-reduced flag)
coal_key <- coal |> group_by(initiator, target, year) |>
  summarise(.n_src = n(), .mixed_incl = n_distinct(coal_incl) > 1L, .mixed_strict = n_distinct(coal_strict) > 1L,
            coal_incl = any(coal_incl), coal_strict = any(coal_strict), coal_onset_ambiguous = any(coal_onset_ambiguous), .groups = "drop")
stopifnot(all(coal_key$.mixed_incl <= (coal_key$.n_src > 1L)), all(coal_key$coal_incl[coal_key$.mixed_incl]))
n_ev_c <- nrow(events)
events <- events |> left_join(coal_key, by = c("initiator", "target", "year")) |>
  mutate(across(c(coal_incl, coal_strict, coal_onset_ambiguous, .mixed_incl, .mixed_strict), ~ ifelse(escalation, NA, .x)), .n_src = ifelse(escalation, NA_integer_, .n_src))   # escalation rows (in no sample) carry NA
stopifnot(nrow(events) == n_ev_c, all(is.na(events$coal_incl) == events$escalation), all(is.na(events$coal_strict) == events$escalation),
          all(events$.n_src[!events$escalation] == events$n_disputes[!events$escalation]), all(events$coal_incl[!events$escalation] >= events$coal_strict[!events$escalation]))
blog("coalition context after the collapse: %d non-escalation events; coalition inclusive %d, strict %d; collapsed keys whose source disputes disagree (inclusive / strict; the ORed flag makes them coalition events): %d / %d of %d collapsed keys",
     sum(!events$escalation), sum(events$coal_incl, na.rm = TRUE), sum(events$coal_strict, na.rm = TRUE), sum(events$.mixed_incl, na.rm = TRUE), sum(events$.mixed_strict, na.rm = TRUE), sum(events$n_disputes > 1L & !events$escalation))
mixed_c <- events |> filter(!escalation, .mixed_incl | .mixed_strict)
for (r in seq_len(nrow(mixed_c))) cat(sprintf("  context differs across the collapsed disputes of %d -> %d %d (disputes %s): coal_incl %s (mixed %s), coal_strict %s (mixed %s)\n", mixed_c$initiator[r], mixed_c$target[r], mixed_c$year[r], mixed_c$dispnums[r],
                                            mixed_c$coal_incl[r], mixed_c$.mixed_incl[r], mixed_c$coal_strict[r], mixed_c$.mixed_strict[r]), file = build_log, append = TRUE)
for (s in c("mid", "fmid", "fmidcow", "fmidub", "war")) {
  f <- events[[paste0("in_", s)]]
  blog("  coalition context, sample %-7s (after the collapse, what the context switch sees): %4d events; coalition inclusive %4d (%.1f%%), strict %4d (%.1f%%); bilateral inclusive %4d, strict %4d; ambiguous onsets %d; collapsed keys with disagreeing source contexts (inclusive / strict) %d / %d",
       s, sum(f), sum(events$coal_incl[f]), 100 * mean(events$coal_incl[f]), sum(events$coal_strict[f]), 100 * mean(events$coal_strict[f]), sum(!events$coal_incl[f]), sum(!events$coal_strict[f]),
       sum(events$coal_onset_ambiguous[f]), sum(events$.mixed_incl[f]), sum(events$.mixed_strict[f]))
}
events <- events |> select(-.n_src, -.mixed_incl, -.mixed_strict)
rm(coal_key, mixed_c, n_ev_c, pairs_pre, coal, coal_grp, coal_ev, coal_mem)   # psd stays: multilateral_jan_flags() (stage b) uses it below
## no-WWII flags after the collapse (review P2-3): a subset of the sample flag; equal to in_* & !is_ww2 on every event from a
## single dispute; on mixed keys (a WWII and a non-WWII dispute) the old rule in_* & !is_ww2 dropped the event, the new flag
## keeps it when a non-WWII dispute qualifies
mixed_ww2 <- events |> filter(!escalation, n_disputes > 1L, is_ww2) |> mutate(dn = strsplit(dispnums, ";", fixed = TRUE)) |>
  filter(vapply(dn, function(z) !all(as.integer(z) %in% WW2_DISPNUMS), logical(1)))
for (s in c("mid", "fmid", "fmidcow", "fmidub", "war")) {
  f <- events[[paste0("in_", s)]]; fn <- events[[paste0("in_", s, "_noww2")]]
  stopifnot(!anyNA(fn), all(fn <= f), all(fn[events$n_disputes == 1L] == (f & !events$is_ww2)[events$n_disputes == 1L]),
            all(fn[!events$is_ww2] == f[!events$is_ww2]))
  blog("  sample %-7s --drop_ww2: old rule (in_%s & !is_ww2) keeps %4d events, the no-WWII flag keeps %4d; mixed WWII / non-WWII keys kept by the flag %d",
       s, s, sum(f & !events$is_ww2), sum(fn), sum(fn & events$is_ww2))
}
blog("duplicate collapse: %d collapsed (initiator, target, year) keys mix a WWII dispute with a non-WWII dispute (in_mid_noww2 keeps %d of them): %s",
     nrow(mixed_ww2), sum(mixed_ww2$in_mid_noww2), paste(sprintf("%d->%d %d disputes %s", mixed_ww2$initiator, mixed_ww2$target, mixed_ww2$year, mixed_ww2$dispnums), collapse = "; "))
saveRDS(events, file.path(out_dir, "events.rds"))

## ============================================================================
## 5. Conflict-history building blocks (calendar-indexed, cumulative)
## ============================================================================
## monadic: distinct MID onsets per (state, year) as an originator on either side; ongoing flag
## Each state-dispute is counted ONCE, at the state's first year in the dispute. (disponset == 1 is
## set per directed pair, so a state would otherwise get another "onset" whenever a later opponent
## joins the same dispute.) Ongoing flags are unaffected.
mon_first <- dd |> filter(orig1 == 1L | INCLUDE_JOINERS) |>
  group_by(ccode = ccode1, dispnum) |> summarise(year = min(year), .groups = "drop") |>
  count(ccode, year, name = "onsets")
mon_ongoing <- dd |> filter(orig1 == 1L | INCLUDE_JOINERS) |>
  group_by(ccode = ccode1, year) |> summarise(ongoing = as.integer(any(dispongoing == 1L)), .groups = "drop")
mon_disp <- full_join(mon_first, mon_ongoing, by = c("ccode", "year")) |>
  mutate(onsets = coalesce(onsets, 0L), ongoing = coalesce(ongoing, 0L))
assert_no_dup_keys(mon_disp, c("ccode", "year"), "monadic MID history")
mon_hist <- states |> select(ccode, year) |> left_join(mon_disp, by = c("ccode", "year")) |>
  mutate(onsets = coalesce(onsets, 0L), ongoing = coalesce(ongoing, 0L)) |>
  arrange(ccode, year) |> group_by(ccode) |>
  mutate(cum_onsets = cumsum(onsets),
         last_mid_year = cummax(ifelse(ongoing == 1L, year, min(year) - 1L)),
         ## peace years AS OF year y: y - last MID-ongoing year strictly before y (entry - 1 if none)
         mon_peace_years = year - coalesce(lag(last_mid_year), min(year) - 1L)) |> ungroup()
assert_no_dup_keys(mon_hist, c("ccode", "year"), "monadic history")
stopifnot(all(mon_hist$mon_peace_years >= 1))
## Note: lag() here is on consecutive system years of a state; a state that exits and re-enters
## has a gap, and year - last_mid_year is still a calendar difference, so no row-lag artefact.

## dyadic (undirected key, both directions written): onsets and ongoing between i and j
## dyad-dispute counted once, at the pair's first year in the dispute
dy_first <- dd |> filter((orig1 == 1L & orig2 == 1L) | INCLUDE_JOINERS) |>
  group_by(ccode1, ccode2, dispnum) |> summarise(year = min(year), .groups = "drop") |>
  count(ccode1, ccode2, year, name = "onsets")
dy_ongoing <- dd |> filter((orig1 == 1L & orig2 == 1L) | INCLUDE_JOINERS) |>
  group_by(ccode1, ccode2, year) |> summarise(ongoing = as.integer(any(dispongoing == 1L)), .groups = "drop")
dy_disp <- full_join(dy_first, dy_ongoing, by = c("ccode1", "ccode2", "year")) |>
  mutate(onsets = coalesce(onsets, 0L), ongoing = coalesce(ongoing, 0L))
assert_no_dup_keys(dy_disp, c("ccode1", "ccode2", "year"), "dyadic MID history")
blog("dyadic MID history rows: %d (joiners %s, consistent with the events)", nrow(dy_disp), if (INCLUDE_JOINERS) "included" else "excluded")
## lagged involvement at t-1 (PI decision 2026-10-24), ALWAYS from every participant including joiners
act_state <- active_state_years(dd)            # (ccode, year): in any MID active that year -> ongoing_tm1
act_pair  <- active_pair_years(dd)             # (ccode1, ccode2, year), both directions           -> dy_mid_tm1
assert_no_dup_keys(act_state, c("ccode", "year"), "active state-years"); assert_no_dup_keys(act_pair, c("ccode1", "ccode2", "year"), "active pair-years")
blog("active MID state-years %d; active MID pair-years (both directions) %d; already-fighting directed pair-years at the reference date (stage-2 set, month rule) %d, at January (stage-1 IV set) %d, year rule %d",
     nrow(act_state), nrow(act_pair), nrow(af_s2), nrow(af_iv), nrow(af_dir_year))

## cumulative onsets for a (key, year) table -> helper used below via calendar joins
## prior_mids over [a, b] = cum(b) - cum(a - 1) where cum(y) = onsets in years <= y (0 before entry)
cum_lookup <- function(cum_tbl, keys, ref_year_col, out_col, query) {
  ## query has keys + ref_year_col; returns query with out_col = cum_onsets at the LATEST table
  ## year <= ref year. The table is a per-key grid of in-system years, so:
  ##   ref in grid            -> direct calendar join
  ##   ref < first grid year  -> 0 (no history before the unit existed)
  ##   ref > last grid year   -> value at the last grid year
  ##   ref in an interior gap -> latest grid year <= ref (rare; exit and re-entry)
  q <- query |> mutate(.ref = .data[[ref_year_col]], .rowid = row_number())
  rng <- cum_tbl |> group_by(across(all_of(keys))) |>
    summarise(.first = min(year), .last = max(year), .cum_last = cum_onsets[which.max(year)], .groups = "drop")
  out <- q |> left_join(cum_tbl |> transmute(across(all_of(keys)), .ref = year, .val = cum_onsets), by = c(keys, ".ref")) |>
    left_join(rng, by = keys)
  out <- out |> mutate(.val = case_when(!is.na(.val) ~ .val,
                                        is.na(.first) | .ref < .first ~ 0L,
                                        .ref > .last ~ .cum_last,
                                        TRUE ~ NA_integer_))
  gap <- out |> filter(is.na(.val))
  if (nrow(gap) > 0) {
    fill <- gap |> select(all_of(keys), .ref, .rowid) |>
      inner_join(cum_tbl |> rename(.y = year), by = keys, relationship = "many-to-many") |>
      filter(.y <= .ref) |> group_by(.rowid) |> slice_max(.y, n = 1, with_ties = FALSE) |> ungroup() |>
      select(.rowid, .fill = cum_onsets)
    out <- out |> left_join(fill, by = ".rowid") |> mutate(.val = coalesce(.val, .fill)) |> select(-.fill)
  }
  stopifnot(!anyNA(out$.val))
  out[[out_col]] <- out$.val
  out |> select(-.ref, -.val, -.rowid, -.first, -.last, -.cum_last)
}

## ============================================================================
## 6. Per-lag builds: stage-1 panel + choice sets
## ============================================================================
dyads <- create_dyadyears(system = "cow", directed = TRUE, mry = FALSE, subset_years = FIRST_YEAR:LAST_YEAR) |>
  add_contiguity() |> add_cow_majors() |> add_atop_alliance() |> add_democracy() |>
  add_capital_distance() |>                     # bundled cow_capitals; no download_extdata() needed
  mutate(across(c(ccode1, ccode2, year), as.integer))
stopifnot(!anyNA(dyads$capdist), all(dyads$capdist >= 0))
assert_no_dup_keys(dyads, c("ccode1", "ccode2", "year"), "directed dyad-years")
dyads <- dyads |> transmute(
  ccode1, ccode2, year,
  contig_land = as.integer(coalesce(conttype, 0L) == 1L),
  contig_any  = as.integer(coalesce(conttype, 0L) %in% 1:PR_CONT_MAX),
  major1 = coalesce(cowmaj1, 0L), major2 = coalesce(cowmaj2, 0L),
  alliance_def_off = as.integer(coalesce(atop_defense, 0) == 1 | coalesce(atop_offense, 0) == 1),
  alliance_missing = as.integer(is.na(atop_defense)),
  dem1 = as.integer(coalesce(polity21 >= POLITY_DEM_CUT, FALSE)), dem2 = as.integer(coalesce(polity22 >= POLITY_DEM_CUT, FALSE)),
  dem1_missing = as.integer(is.na(polity21)), dem2_missing = as.integer(is.na(polity22)),
  log_capdist = log1p(capdist))                 # capital-to-capital great-circle distance, km (alternative Z, flag --distance)
blog("directed dyad-years: %d rows; ATOP missing %d; polity1 missing %d", nrow(dyads), sum(dyads$alliance_missing), sum(dyads$dem1_missing))
## dyadic ongoing-MID spell table on the directed dyad-year grid
dy_hist <- dyads |> select(ccode1, ccode2, year) |>
  left_join(dy_disp, by = c("ccode1", "ccode2", "year")) |>
  mutate(onsets = coalesce(onsets, 0L), ongoing = coalesce(ongoing, 0L)) |>
  arrange(ccode1, ccode2, year) |> group_by(ccode1, ccode2) |>
  mutate(cum_onsets = cumsum(onsets),
         last_mid_year = cummax(ifelse(ongoing == 1L, year, min(year) - 1L)),
         dy_peace_years = year - coalesce(lag(last_mid_year), min(year) - 1L)) |> ungroup()
stopifnot(all(dy_hist$dy_peace_years >= 1))
dy_cum <- dy_hist |> select(ccode1, ccode2, year, cum_onsets)
dy_py  <- dy_hist |> select(ccode1, ccode2, year, dy_peace_years)
mon_cum <- mon_hist |> select(ccode, year, cum_onsets)
mon_py  <- mon_hist |> select(ccode, year, mon_peace_years)
rm(dy_hist); invisible(gc())

sy_ctrl <- states |> transmute(ccode, year, cowmaj = coalesce(cowmaj, 0L),
                               dem = as.integer(coalesce(polity2 >= POLITY_DEM_CUT, FALSE)),
                               dem_missing = as.integer(is.na(polity2)), n_system)
pos_k <- pos |> select(ccode, year, M, R_pr, R_eig)

## ---- gap scale (PI decision 2026-10-23, "option 1") ---------------------------------------------------------
## The stage-2 gaps are log ratios, M_i - M_j = log(CINC_i / CINC_j) and R_pr_i - R_pr_j = log(PR_i / PR_j). The two
## axes have very different spreads, so each gap is divided by ONE constant per axis, the same for every year, lag and
## choice-set variant. The constant is the SD of the within-year gap over all ordered pairs of distinct states with
## positions in the same year (the mean of these pairs is exactly 0). Exact formula per year with n states and
## within-year sum of squares SS = sum((v - mean(v))^2), sum over ordered pairs i != j of (v_i - v_j)^2 = 2 n SS, and
## there are n (n - 1) pairs, so the pooled variance is sum_y(2 n_y SS_y) / sum_y(n_y (n_y - 1)).
## One unit of a scaled gap = GAP_SCALE log points, i.e. a power ratio of exp(GAP_SCALE).
## (gap_scale() is in helpers_data.R so it can be unit-tested.)
GAP_SCALE_M <- gap_scale(pos$M, pos$year)
GAP_SCALE_R <- gap_scale(pos$R_pr, pos$year)
stopifnot(is.finite(GAP_SCALE_M), GAP_SCALE_M > 0, is.finite(GAP_SCALE_R), GAP_SCALE_R > 0)
blog("gap scale (one constant per axis, all years): M (log CINC ratio) %.6f; R_pr (log PageRank ratio) %.6f; 1 scaled unit = power ratio %.2f (M), %.2f (R_pr)",
     GAP_SCALE_M, GAP_SCALE_R, exp(GAP_SCALE_M), exp(GAP_SCALE_R))
## D2: within-year centred, gap-scaled levels of every state with positions (joined to the panel at t-1 below)
pos_lvl <- pos_k |> group_by(year) |> mutate(M_lvl = (M - mean(M)) / GAP_SCALE_M, R_lvl = (R_pr - mean(R_pr)) / GAP_SCALE_R) |> ungroup()
lvl_chk <- pos_lvl |> group_by(year) |> summarise(mM = mean(M_lvl), mR = mean(R_lvl), .groups = "drop")
stopifnot(all(abs(lvl_chk$mM) < 1e-10), all(abs(lvl_chk$mR) < 1e-10), all(is.finite(pos_lvl$M_lvl)), all(is.finite(pos_lvl$R_lvl)),
          all(abs(pos_lvl$M_lvl * GAP_SCALE_M - (pos_lvl$M - ave(pos_lvl$M, pos_lvl$year))) < 1e-10))
blog("power levels (D2): M_lvl = (log CINC share - year mean) / %.6f, R_lvl = (log PageRank - year mean) / %.6f over all %d state-years with positions; SD %.3f / %.3f, 5/50/95 pct %.2f/%.2f/%.2f and %.2f/%.2f/%.2f",
     GAP_SCALE_M, GAP_SCALE_R, nrow(pos_lvl), sd(pos_lvl$M_lvl), sd(pos_lvl$R_lvl), quantile(pos_lvl$M_lvl, .05), median(pos_lvl$M_lvl), quantile(pos_lvl$M_lvl, .95),
     quantile(pos_lvl$R_lvl, .05), median(pos_lvl$R_lvl), quantile(pos_lvl$R_lvl, .95))
## ---- multilateral MID open at January of t (dual typology stage b, 2026-10-27; the M1c stage-1 control; helper in helpers_data.R) --------
## multi_mid_jan(i, t) = 1 when a MID with more than one state on a side (participants present by their spell dates: entered in a year
## before t or day-dated on or before 1 January of t, not exited; COW sides in year t) is active at 1 January of t and i has no spell in
## it entered before that date (audit fixes 1-2, 2026-10-27: a coarse-dated entry in t never counts, a former participant is "already
## in"). Same participant spells (psd) and date keys as the coalition flag. Computed once on every system state-year and joined to each
## lag's panel below; outcome-blind by construction (nothing dated in t enters). It enters W in the M1c model ONLY (assemble_cell(CONTEXT_MODEL = TRUE)).
multi_jan <- multilateral_jan_flags(states |> select(ccode, year), psd, dd)
multi_disp <- attr(multi_jan, "disputes")
mmj_chg <- attr(multi_jan, "changed_vs_prefix")
blog("multi_mid_jan audit fixes 1-2 (2026-10-27), against the pre-fix rule: %d spell-years were present at 1 January only through a year- or month-resolution entry in t; multilateral dispute-years %d before, %d after (%d lost, %d gained); %d (dispute, state, year) former-participant triples now count as already in; state-year flags changed by fix 1 (coarse entries) %d, by fix 2 (former participants) %d, total %d of %d; flags set %d before, %d after",
     mmj_chg[["n_spell_years_coarse_entry_in_t"]], mmj_chg[["n_multi_dispute_years_legacy"]], mmj_chg[["n_multi_dispute_years"]], mmj_chg[["n_multi_dispute_years_lost"]], mmj_chg[["n_multi_dispute_years_gained"]],
     mmj_chg[["n_former_participant_triples"]], mmj_chg[["n_flags_changed_fix1"]], mmj_chg[["n_flags_changed_fix2"]], mmj_chg[["n_flags_changed_total"]], nrow(multi_jan), mmj_chg[["n_flags_legacy"]], mmj_chg[["n_flags"]])
blog("multi_mid_jan (stage b): active multilateral dispute-years at 1 January %d (of %d active dispute-years), in %d distinct years; state-years with at least one such dispute they had not entered: %d of %d (mean %.4f); max disputes per state-year %d; participants on both sides within a dispute-year (side of t-1 used at 1 January) %d",
     nrow(multi_disp), attr(multi_jan, "n_active_dispute_years"), n_distinct(multi_disp$year), sum(multi_jan$multi_mid_jan), nrow(multi_jan), mean(multi_jan$multi_mid_jan), max(multi_jan$n_multi_jan), attr(multi_jan, "n_side_switch_years"))
blog("multi_mid_jan: the 10 largest active multilateral disputes at a 1 January (dispute, year, states present on side A / side B): %s",
     paste(with(multi_disp |> arrange(desc(n_side_a + n_side_b), year) |> slice_head(n = 10), sprintf("%d/%d %d/%d", dispnum, year, n_side_a, n_side_b)), collapse = "; "))
if (DIRECTION == "cowwar") {
  ## the Second World War (dispute 258) is open and multilateral at 1 January 1940-1945, and the USA (2) had not entered it at 1 January 1941
  stopifnot(all(1940:1945 %in% multi_disp$year[multi_disp$dispnum == 258L]), multi_jan$multi_mid_jan[multi_jan$ccode == 2L & multi_jan$year == 1941L] == 1L)
}
manifest <- list()
add_manifest <- function(obj, fname, note) {
  p <- file.path(out_dir, fname); saveRDS(obj, p)
  manifest[[fname]] <<- data.frame(file = fname, nrow = nrow(obj), ncol = ncol(obj), sha256 = file_digest(p), note = note)
  blog("wrote %s: %d x %d", fname, nrow(obj), ncol(obj))
}
add_manifest(pos, "sy_positions.rds", "state-year positions: M = log CINC share, R_pr = log PageRank (raw log shares), R_eig std within year")
add_manifest(events, "events.rds", sprintf("directed events (Side A x Side B pairs; DIRECTION=%s: war-level disputes follow COW's war-initiating side, dir_flipped / dir_source / war_link), joiners %s (is_joiner_pair); month-rule escalation flag with participant-level dates where gml_part matches (excluded from all in_* flags); pair-level war; same-year duplicates collapsed (dispnums); fighting_at_ref; is_ww2; target_has_pos_tm1; stage-2 weights w_<sample>; coalition context coal_incl / coal_strict / coal_onset_ambiguous",
                                           DIRECTION, if (INCLUDE_JOINERS) "included" else "excluded"))
## number of system states per year, for the system-growth diagnostic on the change variables (item 4 of the
## 2026-10-23 decision): g3 at t-1 = log n_system(t-1) - log n_system(t-4), constant within a year
nsys_year <- states |> distinct(year, n_system) |> arrange(year)
assert_no_dup_keys(nsys_year, "year", "system size by year")
growth3 <- function(year_tm1) log(nsys_year$n_system[match(year_tm1, nsys_year$year)]) -
  log(nsys_year$n_system[match(year_tm1 - 3L, nsys_year$year)])
scale_rows <- list()

for (h in LAGS) {
  blog("---- lag h = %d ----", h)
  ## ---- power changes: ONE state-level table per lag (PI decision 2026-10-23) ---------------------------
  ## raw change = level(t-1) - level(t-1-h) for every state in the system at t with levels at both dates;
  ## centred within year (minus the mean raw change of all states with a change that year) and divided by one
  ## pooled SD per lag and variable. Used for the initiator (stage 1, W and X) AND the target (stage 2, Z).
  chg <- build_change_table(pos_k, states |> select(ccode, year), h)
  scale_rows[[as.character(h)]] <- chg$scale
  blog("changes h=%d: %d state-years over %d years; pooled SD of the year-centred raw change: M %.5f, R_pr %.5f (the scaling constants); raw change SD M %.4f, R_pr %.4f",
       h, nrow(chg$tbl), n_distinct(chg$tbl$year), chg$scale$pooled_sd[1], chg$scale$pooled_sd[2], sd(chg$tbl$dM_raw_h), sd(chg$tbl$dRpr_raw_h))
  add_manifest(chg$tbl, sprintf("changes_lag%d.rds", h),
               sprintf("power changes h=%d: raw and (year-centred / pooled SD) transformed; pooled SD M=%.6f R_pr=%.6f", h, chg$scale$pooled_sd[1], chg$scale$pooled_sd[2]))

  ## ---- stage-1 panel: states in system at t with positions at t-1 and t-1-h ---------
  ## Power LEVELS at t-1 (PI decision 2026-10-25, D2): the within-year centred log CINC share and log PageRank, divided by the
  ## gap constants so one unit matches the stage-2 gap axes. Centred over EVERY state with positions in that year (pos, the same
  ## set the gap scale uses), not over the panel rows, so the reference (0 = the system-average level that year) does not depend
  ## on h. They enter the main stage-1 W (relative standing; W had the power CHANGES only, so a "capable states initiate more"
  ## effect could flow into mu0 through the inclusive value); they are NOT in X (open PI question).
  panel <- states |> select(ccode, year) |>
    inner_join(pos_lvl |> transmute(ccode, year = year + 1L, M_tm1 = M, Rpr_tm1 = R_pr, Reig_tm1 = R_eig, M_lvl, R_lvl), by = c("ccode", "year")) |>
    inner_join(pos_k |> transmute(ccode, year = year + 1L + h, M_w0 = M, Rpr_w0 = R_pr, Reig_w0 = R_eig), by = c("ccode", "year")) |>
    mutate(dReig_h = Reig_tm1 - Reig_w0, h = h, win_start = year - 1L - h) |>     # R_eig change unchanged (z-score difference)
    ## initiator power changes from the change table (dM_h, dRpr_h transformed; dM_raw_h, dRpr_raw_h raw)
    inner_join(chg$tbl |> select(-h), by = c("ccode", "year")) |>
    ## controls at t-1 (calendar join)
    inner_join(sy_ctrl |> transmute(ccode, year = year + 1L, cowmaj_tm1 = cowmaj, dem_tm1 = dem, dem_missing_tm1 = dem_missing,
                                    n_system_tm1 = n_system), by = c("ccode", "year"))
  ## the panel and the change table have the same domain; the raw change must equal the panel's own level difference
  stopifnot(nrow(panel) == nrow(chg$tbl),
            all(abs(panel$dM_raw_h - (panel$M_tm1 - panel$M_w0)) < 1e-10),
            all(abs(panel$dRpr_raw_h - (panel$Rpr_tm1 - panel$Rpr_w0)) < 1e-10))
  ## monadic conflict history: onsets in [win_start - 10, win_start - 1]; peace years at win_start
  panel <- cum_lookup(mon_cum, "ccode", "win_end_hi", "cum_hi", panel |> mutate(win_end_hi = win_start - 1L))
  panel <- cum_lookup(mon_cum, "ccode", "win_end_lo", "cum_lo", panel |> mutate(win_end_lo = win_start - 1L - PREWINDOW_YEARS))
  panel <- panel |> mutate(mon_prior_mids_h = cum_hi - cum_lo) |> select(-cum_hi, -cum_lo, -win_end_hi, -win_end_lo)
  panel <- panel |> left_join(mon_py |> transmute(ccode, win_start = year, mon_peace_years_h = mon_peace_years), by = c("ccode", "win_start"))
  ## a state not in the system at win_start has no spell there; use years since entry as censored value
  panel <- panel |> left_join(entry, by = "ccode") |>
    mutate(mon_peace_years_h = coalesce(mon_peace_years_h, pmax(win_start - entry_year + 1L, 1L)))
  stopifnot(all(panel$mon_prior_mids_h >= 0))
  ## ongoing involvement at t-1 (PI decision 2026-10-24): participant in any MID active in t-1, joiners included
  panel <- panel |> left_join(act_state |> transmute(ccode, year = year + 1L, ongoing_tm1 = active), by = c("ccode", "year")) |>
    mutate(ongoing_tm1 = coalesce(ongoing_tm1, 0L))
  stopifnot(all(panel$ongoing_tm1 %in% c(0L, 1L)))
  ## multilateral MID open at January of t (stage b, 2026-10-27): every panel state-year is a system state-year, so the join is exact
  n_panel <- nrow(panel)
  panel <- panel |> inner_join(multi_jan |> select(ccode, year, multi_mid_jan, n_multi_jan), by = c("ccode", "year"))
  stopifnot(nrow(panel) == n_panel, all(panel$multi_mid_jan %in% c(0L, 1L)), !anyNA(panel$n_multi_jan))
  ## initiation indicators per sample. The sample flags exclude escalations (PI decision 2026-10-24), so a state-year
  ## whose only events are escalations has y1 = 0 and n_ev = 0. n_ev_mid counts in_mid, not every event row.
  ev_sy <- events |> group_by(ccode = initiator, year) |>
    summarise(n_ev_mid = sum(in_mid), n_ev_fmid = sum(in_fmid), n_ev_fmidcow = sum(in_fmidcow), n_ev_fmidub = sum(in_fmidub),
              n_ev_war = sum(in_war), n_ev_mid_incl_esc = n(), n_ev_war_incl_esc = sum(war %in% TRUE), .groups = "drop")
  panel <- panel |> left_join(ev_sy, by = c("ccode", "year")) |>
    mutate(across(c(n_ev_mid, n_ev_fmid, n_ev_fmidcow, n_ev_fmidub, n_ev_war, n_ev_mid_incl_esc, n_ev_war_incl_esc), ~ coalesce(.x, 0L)),
           y1_mid = as.integer(n_ev_mid > 0), y1_fmid = as.integer(n_ev_fmid > 0),
           y1_fmidcow = as.integer(n_ev_fmidcow > 0),
           y1_fmidub = as.integer(n_ev_fmidub > 0), y1_war = as.integer(n_ev_war > 0),
           log_n_system_tm1 = log(n_system_tm1)) |>
    arrange(year, ccode)
  ## y1 before / after the escalation exclusion (the "before" counts are logged only, not saved as y1 columns)
  y1_incl <-sapply(list(mid = rep(TRUE, nrow(events)), fmid = events$fatal_main, fmidcow = events$fatal_cow, fmidub = events$fatal_ub, war = events$war %in% TRUE),
                    function(cd) nrow(events[cd, ] |> distinct(ccode = initiator, year) |> semi_join(panel, by = c("ccode", "year"))))
  blog("panel h=%d: stage-1 initiations incl. escalations mid/fmid/fmidcow/fmidub/war = %d/%d/%d/%d/%d; after excluding escalations = %d/%d/%d/%d/%d (state-years whose only events were escalations: mid %d, war %d)",
       h, y1_incl["mid"], y1_incl["fmid"], y1_incl["fmidcow"], y1_incl["fmidub"], y1_incl["war"],
       sum(panel$y1_mid), sum(panel$y1_fmid), sum(panel$y1_fmidcow), sum(panel$y1_fmidub), sum(panel$y1_war),
       sum(panel$y1_mid == 0L & panel$n_ev_mid_incl_esc > 0L), sum(panel$y1_war == 0L & panel$n_ev_war_incl_esc > 0L))
  panel <- panel |> select(-n_ev_mid_incl_esc, -n_ev_war_incl_esc)
  ## recount y1 from the non-escalation events (independent path through the escalation column)
  y1_re <- events |> filter(!escalation) |> distinct(ccode = initiator, year) |> mutate(.y = 1L)
  chk_y1 <- panel |> select(ccode, year, y1_mid) |> left_join(y1_re, by = c("ccode", "year")) |> mutate(.y = coalesce(.y, 0L))
  stopifnot(all(chk_y1$y1_mid == chk_y1$.y))
  blog("  ongoing_tm1 h=%d: mean %.4f (%d of %d state-years in a MID active at t-1)", h, mean(panel$ongoing_tm1), sum(panel$ongoing_tm1), nrow(panel))
  blog("  multi_mid_jan h=%d (stage b): mean %.4f (%d of %d state-years with a multilateral MID they had not entered active at 1 January of t); among initiators mid/fmid/war %.4f/%.4f/%.4f; cor with ongoing_tm1 %.3f",
       h, mean(panel$multi_mid_jan), sum(panel$multi_mid_jan), nrow(panel), mean(panel$multi_mid_jan[panel$y1_mid == 1L]), mean(panel$multi_mid_jan[panel$y1_fmid == 1L]),
       mean(panel$multi_mid_jan[panel$y1_war == 1L]), cor(panel$multi_mid_jan, panel$ongoing_tm1))
  assert_no_dup_keys(panel, c("ccode", "year"), sprintf("stage-1 panel h=%d", h))
  assert_no_na(panel, c("dM_h", "dRpr_h", "dM_raw_h", "dRpr_raw_h", "dReig_h", "M_tm1", "Rpr_tm1", "M_lvl", "R_lvl", "log_n_system_tm1", "cowmaj_tm1", "dem_tm1", "ongoing_tm1", "multi_mid_jan",
                        "mon_prior_mids_h", "mon_peace_years_h", "y1_mid", "y1_fmid", "y1_fmidcow", "y1_fmidub", "y1_war"), "stage-1 panel")
  assert_finite(panel, c("dM_h", "dRpr_h", "dM_raw_h", "dRpr_raw_h", "M_tm1", "Rpr_tm1", "M_lvl", "R_lvl", "log_n_system_tm1"), sprintf("stage-1 panel h=%d", h))
  ## D2 levels on the panel rows: they are centred over ALL states with positions at t-1, so the panel mean is near but not exactly 0
  stopifnot(all(abs(panel$M_lvl * GAP_SCALE_M - (panel$M_tm1 - ave(pos_lvl$M, pos_lvl$year)[match(panel$year - 1L, pos_lvl$year)])) < 1e-8))
  blog("  power levels h=%d (D2): M_lvl mean %.4f SD %.3f, R_lvl mean %.4f SD %.3f over the panel rows; cor(M_lvl, dM_h) %.3f, cor(R_lvl, dRpr_h) %.3f, cor(M_lvl, R_lvl) %.3f; log_n_system_tm1 (D3) mean %.3f SD %.3f",
       h, mean(panel$M_lvl), sd(panel$M_lvl), mean(panel$R_lvl), sd(panel$R_lvl), cor(panel$M_lvl, panel$dM_h), cor(panel$R_lvl, panel$dRpr_h), cor(panel$M_lvl, panel$R_lvl),
       mean(panel$log_n_system_tm1), sd(panel$log_n_system_tm1))
  ## invariants of the transformed changes on the panel rows (item 4): within-year mean 0, pooled SD 1, and no
  ## correlation with system growth (g3 is constant within a year, so the year-centring removes it by construction)
  yr_chk <- panel |> group_by(year) |> summarise(mM = mean(dM_h), mR = mean(dRpr_h), .groups = "drop")
  stopifnot(all(abs(yr_chk$mM) < 1e-8), all(abs(yr_chk$mR) < 1e-8),
            abs(sd(panel$dM_h) - 1) < 1e-8, abs(sd(panel$dRpr_h) - 1) < 1e-8)
  g3 <- growth3(panel$year - 1L)
  ok_g <- is.finite(g3)
  cor_g <- c(raw_M = cor(panel$dM_raw_h[ok_g], g3[ok_g]), raw_R = cor(panel$dRpr_raw_h[ok_g], g3[ok_g]),
             new_M = cor(panel$dM_h[ok_g], g3[ok_g]), new_R = cor(panel$dRpr_h[ok_g], g3[ok_g]))
  blog("  system-growth check h=%d (cor with log n_system(t-1) - log n_system(t-4), %d rows): raw change M %.3f, R_pr %.3f; transformed change M %.4f, R_pr %.4f (2026-10-22 z-score version was about +0.29 / +0.17 at h = 3)",
       h, sum(ok_g), cor_g["raw_M"], cor_g["raw_R"], cor_g["new_M"], cor_g["new_R"])
  if (!all(is.finite(cor_g)) || abs(cor_g["new_M"]) >= 0.02 || abs(cor_g["new_R"]) >= 0.02)
    stop(sprintf("h=%d: transformed changes still correlate with system growth (M %.3f, R_pr %.3f); year-centring failed", h, cor_g["new_M"], cor_g["new_R"]))
  blog("  changes h=%d: transformed dM_h SD %.3f, 5/50/95 pct %.3f/%.3f/%.3f; dRpr_h SD %.3f, 5/50/95 pct %.3f/%.3f/%.3f",
       h, sd(panel$dM_h), quantile(panel$dM_h, .05), median(panel$dM_h), quantile(panel$dM_h, .95),
       sd(panel$dRpr_h), quantile(panel$dRpr_h, .05), median(panel$dRpr_h), quantile(panel$dRpr_h, .95))
  blog("panel h=%d: %d state-years, years %d-%d (%d distinct), initiations mid/fmid/fmidcow/war = %d/%d/%d/%d",
       h, nrow(panel), min(panel$year), max(panel$year), n_distinct(panel$year), sum(panel$y1_mid), sum(panel$y1_fmid), sum(panel$y1_fmidcow), sum(panel$y1_war))
  ev_lost <- events |> filter(in_mid) |> anti_join(panel, by = c("initiator" = "ccode", "year"))
  blog("  non-escalation events whose initiator state-year is NOT in the panel (no positions at t-1 / t-1-h): %d of %d", nrow(ev_lost), sum(events$in_mid))
  add_manifest(panel, sprintf("stage1_panel_lag%d.rds", h),
               sprintf("stage-1 panel, h=%d; y1 excludes escalations; ongoing_tm1; multi_mid_jan (stage b); levels M_lvl / R_lvl at t-1 (within-year centred / gap constants); log_n_system_tm1; dM_h / dRpr_h year-centred and scaled by pooled SD M=%.6f R_pr=%.6f", h, chg$scale$pooled_sd[1], chg$scale$pooled_sd[2]))

  ## ---- choice sets: panel state-year x every other system member at t with positions at t-1 ----
  alts <- panel |> select(ccode1 = ccode, year, M_i = M_tm1, Rpr_i = Rpr_tm1, Reig_i = Reig_tm1, win_start) |>
    inner_join(states |> select(ccode2 = ccode, year), by = "year", relationship = "many-to-many") |>
    filter(ccode1 != ccode2) |>
    inner_join(pos_k |> transmute(ccode2 = ccode, year = year + 1L, M_j = M, Rpr_j = R_pr, Reig_j = R_eig), by = c("ccode2", "year")) |>
    ## gaps at t-1: log ratios of CINC share / PageRank, kept raw (dM_logratio, dR_pr_logratio) and, for the model,
    ## divided by the one-per-axis constants GAP_SCALE_M / GAP_SCALE_R (see "gap scale" above); dR_eig is a z-score difference
    mutate(dM_logratio = M_i - M_j, dR_pr_logratio = Rpr_i - Rpr_j,
           dM = dM_logratio / GAP_SCALE_M, dR_pr = dR_pr_logratio / GAP_SCALE_R, dR_eig = Reig_i - Reig_j) |>
    ## target power changes over the same window, from the SAME change table as the initiator (stage-2 Z columns).
    ## A target with no change value (no position at t-1-h, e.g. a recent entrant) gets 0 and dch_j_missing = 1.
    left_join(chg$tbl |> transmute(ccode2 = ccode, year, dM_j_h = dM_h, dRpr_j_h = dRpr_h), by = c("ccode2", "year"))
  stopifnot(all(is.na(alts$dM_j_h) == is.na(alts$dRpr_j_h)))                     # both come from the one table
  alts <- alts |>
    mutate(dch_j_missing = as.integer(is.na(dM_j_h)),
           dM_j_h = coalesce(dM_j_h, 0), dRpr_j_h = coalesce(dRpr_j_h, 0)) |>
    ## dyadic controls at t-1
    left_join(dyads |> mutate(year = year + 1L), by = c("ccode1", "ccode2", "year"))
  stopifnot(all(alts$dch_j_missing %in% c(0L, 1L)), all(alts$dM_j_h[alts$dch_j_missing == 1L] == 0), all(alts$dRpr_j_h[alts$dch_j_missing == 1L] == 0))
  n_miss <- sum(alts$dch_j_missing); miss_tg <- alts |> filter(dch_j_missing == 1L) |> distinct(ccode2, year)
  blog("  target change h=%d: %d of %d alternative rows have a target with NO change value (dch_j_missing = 1; changes set to 0); %d distinct target state-years, %d distinct targets",
       h, n_miss, nrow(alts), nrow(miss_tg), n_distinct(miss_tg$ccode2))
  stopifnot(all(abs(alts$dM * GAP_SCALE_M - alts$dM_logratio) < 1e-10), all(abs(alts$dR_pr * GAP_SCALE_R - alts$dR_pr_logratio) < 1e-10))
  blog("  gaps h=%d (log ratios): dM SD %.3f, 5/50/95 pct %.3f/%.3f/%.3f; dR_pr SD %.3f, 5/50/95 pct %.3f/%.3f/%.3f",
       h, sd(alts$dM_logratio), quantile(alts$dM_logratio, .05), median(alts$dM_logratio), quantile(alts$dM_logratio, .95),
       sd(alts$dR_pr_logratio), quantile(alts$dR_pr_logratio, .05), median(alts$dR_pr_logratio), quantile(alts$dR_pr_logratio, .95))
  blog("  gaps h=%d (scaled, model units): dM SD %.3f, 5/50/95 pct %.3f/%.3f/%.3f; dR_pr SD %.3f, 5/50/95 pct %.3f/%.3f/%.3f",
       h, sd(alts$dM), quantile(alts$dM, .05), median(alts$dM), quantile(alts$dM, .95),
       sd(alts$dR_pr), quantile(alts$dR_pr, .05), median(alts$dR_pr), quantile(alts$dR_pr, .95))
  n_ctrl_na <- sum(is.na(alts$contig_any))
  if (n_ctrl_na) blog("  alternatives with no dyad-year control row at t-1 (filled 0): %d", n_ctrl_na)
  stopifnot(n_ctrl_na == 0 || all(is.na(alts$log_capdist[is.na(alts$contig_any)])))
  alts <- alts |> mutate(across(c(contig_land, contig_any, major1, major2, alliance_def_off, alliance_missing,
                                  dem1, dem2, dem1_missing, dem2_missing), ~ coalesce(.x, 0L)),
                         major_count = major1 + major2, joint_dem = dem1 * dem2,
                         polrel = as.integer(contig_any == 1L | major1 == 1L | major2 == 1L))
  ## dyadic conflict history relative to the window
  alts <- cum_lookup(dy_cum, c("ccode1", "ccode2"), "win_hi", "cum_hi", alts |> mutate(win_hi = win_start - 1L))
  alts <- cum_lookup(dy_cum, c("ccode1", "ccode2"), "win_lo", "cum_lo", alts |> mutate(win_lo = win_start - 1L - PREWINDOW_YEARS))
  alts <- alts |> mutate(dy_prior_mids_h = cum_hi - cum_lo) |> select(-cum_hi, -cum_lo, -win_hi, -win_lo)
  alts <- alts |> left_join(dy_py |> rename(win_start = year, dy_peace_years_h = dy_peace_years), by = c("ccode1", "ccode2", "win_start")) |>
    left_join(entry |> rename(ccode1 = ccode, entry1 = entry_year), by = "ccode1") |>
    left_join(entry |> rename(ccode2 = ccode, entry2 = entry_year), by = "ccode2") |>
    mutate(dy_peace_years_h = coalesce(dy_peace_years_h, pmax(win_start - pmax(entry1, entry2) + 1L, 1L))) |>
    select(-entry1, -entry2)
  stopifnot(all(alts$dy_prior_mids_h >= 0), all(alts$dy_peace_years_h >= 1))
  ## ---- ongoing-dispute rules (PI decision 2026-10-24, month resolution; Codex review P1-1 "A'"): two sets per state-year.
  ## The stage-1 IV set drops (i, j, t) whose pair is in a dispute active at JANUARY of t that started before it (every
  ## state-year, no event information); the stage-2 set drops the pairs fighting at the state-year's reference date (start
  ## of its earliest non-escalation event, else January). The row block is the union with the flags in_iv / in_s2; rows in
  ## neither set leave the file (both variants are cut from this table; af_iv and af_s2 are lag-independent). Escalation
  ## events lose their stage-2 row, which is intended. A non-escalation event may lose its stage-2 row ONLY when flagged
  ## fighting_at_ref (pair still fighting at the reference date, dispute ended before the event); any other loss stops the build.
  n_alts_before <- nrow(alts)
  key_alt <- paste(alts$ccode1, alts$ccode2, alts$year)
  n_drop_iv <- sum(key_alt %in% paste(af_iv$ccode1, af_iv$ccode2, af_iv$year))
  n_drop_s2 <- sum(key_alt %in% paste(af_s2$ccode1, af_s2$ccode2, af_s2$year))
  n_drop_year <- sum(key_alt %in% paste(af_dir_year$ccode1, af_dir_year$ccode2, af_dir_year$year))
  n_drop_s2_year <- sum(key_alt %in% paste(af_s2$ccode1, af_s2$ccode2, af_s2$year) & key_alt %in% paste(af_dir_year$ccode1, af_dir_year$ccode2, af_dir_year$year))
  alts <- choice_set_flags(alts, af_iv, af_s2)
  n_neither <- attr(alts, "n_neither")
  stopifnot(nrow(alts) == n_alts_before - n_neither, sum(alts$in_iv == 0L) == n_drop_iv - n_neither, sum(alts$in_s2 == 0L) == n_drop_s2 - n_neither,
            nrow(alts |> filter(in_iv == 1L) |> semi_join(af_iv, by = c("ccode1", "ccode2", "year"))) == 0,
            nrow(alts |> filter(in_s2 == 1L) |> semi_join(af_s2, by = c("ccode1", "ccode2", "year"))) == 0)
  ## the IV flag is the January rule and nothing else: it equals "not an af_iv key" on every retained row (af_iv was built
  ## without the reference dates, so the stage-1 IV set cannot depend on the state-year's events)
  stopifnot(all((alts$in_iv == 1L) == !(paste(alts$ccode1, alts$ccode2, alts$year) %in% paste(af_iv$ccode1, af_iv$ccode2, af_iv$year))))
  ev_keep <- events |> filter(!escalation) |> transmute(ccode1 = initiator, ccode2 = target, year, dispnum, fighting_at_ref)
  ## non-escalation events whose (i, j, t) existed before the rules and is dropped from the stage-2 set (row kept with
  ## in_s2 = 0, or removed because it is in neither set)
  s2_dropped_keys <- key_alt[key_alt %in% paste(af_s2$ccode1, af_s2$ccode2, af_s2$year)]
  lost_af <- ev_keep |> filter(paste(ccode1, ccode2, year) %in% s2_dropped_keys)
  if (any(!lost_af$fighting_at_ref)) stop(sprintf("h=%d: %d non-escalation events not flagged fighting_at_ref have no stage-2 row", h, sum(!lost_af$fighting_at_ref)))
  stopifnot(nrow(ev_keep |> filter(fighting_at_ref) |> semi_join(alts |> filter(in_s2 == 1L), by = c("ccode1", "ccode2", "year"))) == 0)
  esc_keys <- events |> filter(escalation) |> transmute(ccode1 = initiator, ccode2 = target, year)
  n_esc_rows <- sum(key_alt %in% paste(af_s2$ccode1, af_s2$ccode2, af_s2$year) & key_alt %in% paste(esc_keys$ccode1, esc_keys$ccode2, esc_keys$year))
  sets_differ <- alts |> group_by(ccode1, year) |> summarise(d = any(in_iv != in_s2), .groups = "drop") |> filter(d)
  init_sy <- ref_dates |> semi_join(panel, by = c("ccode", "year"))
  stopifnot(nrow(sets_differ |> anti_join(init_sy, by = c("ccode1" = "ccode", "year"))) == 0)   # only initiating state-years can differ
  blog("  choice sets h=%d: %d alternative rows before the already-fighting rules; stage-1 IV set (January) drops %d rows, stage-2 set (reference date) drops %d (%d of them were escalation events' rows, %d fighting_at_ref events' rows); %d rows in neither set leave the file; %d rows remain in the union (in both sets %d; IV only %d; stage 2 only %d); state-years whose two sets differ %d of %d initiating state-years (%d state-years in the panel); the year rule would drop %d (stage-2 set and year rule both %d)",
       h, n_alts_before, n_drop_iv, n_drop_s2, n_esc_rows, nrow(lost_af), n_neither, nrow(alts), sum(alts$in_iv == 1L & alts$in_s2 == 1L),
       sum(alts$in_iv == 1L & alts$in_s2 == 0L), sum(alts$in_iv == 0L & alts$in_s2 == 1L), nrow(sets_differ), nrow(init_sy), nrow(panel), n_drop_year, n_drop_s2_year)
  rm(key_alt)
  ## dyadic lagged DV (PI decision 2026-10-24): i and j both in a MID active in t-1, any role, either direction.
  ## After the already-fighting cut, a remaining alternative can have dy_mid_tm1 = 1 only when every dispute of the
  ## pair active in t-1 had ended by t (or the pair's new dispute starts in t, i.e. the row is a chosen event).
  alts <- alts |> left_join(act_pair |> transmute(ccode1, ccode2, year = year + 1L, dy_mid_tm1 = active), by = c("ccode1", "ccode2", "year")) |>
    mutate(dy_mid_tm1 = coalesce(dy_mid_tm1, 0L))
  stopifnot(all(alts$dy_mid_tm1 %in% c(0L, 1L)))
  ## event flags on alternative rows: is this (i, j, t) a chosen target in each sample? A chosen target needs a row in the
  ## STAGE-2 set (in_s2 = 1); a row with in_s2 = 0 counts as absent here
  ev_flag <- events |> filter(in_mid) |> transmute(ccode1 = initiator, ccode2 = target, year, dispnum, in_mid, in_fmid, in_war, target_has_pos_tm1, fighting_at_ref)
  assert_no_dup_keys(ev_flag, c("ccode1", "ccode2", "year"), sprintf("kept events h=%d (one chosen alternative per (i, j, t))", h))
  ev_in_panel <- ev_flag |> semi_join(panel, by = c("ccode1" = "ccode", "year"))
  alts_s2 <- alts |> filter(in_s2 == 1L)
  ev_lost_target <- ev_in_panel |> anti_join(alts_s2, by = c("ccode1", "ccode2", "year"))
  blog("  events whose target is NOT in the stage-2 choice set: %d (target lacks positions at t-1 or is not a system member at t %d; pair already fighting at the reference date %d, of which %d keep an IV-only row with in_s2 = 0)",
       nrow(ev_lost_target), sum(!ev_lost_target$target_has_pos_tm1), sum(ev_lost_target$target_has_pos_tm1 & ev_lost_target$fighting_at_ref),
       nrow(ev_lost_target |> filter(fighting_at_ref) |> semi_join(alts, by = c("ccode1", "ccode2", "year"))))
  ## every non-escalation event with its initiator in the panel, a target with positions at t-1 and a pair not fighting at
  ## the reference date has its stage-2 row; nothing else is ever lost
  stopifnot(all(!ev_lost_target$target_has_pos_tm1 | ev_lost_target$fighting_at_ref),
            nrow(ev_in_panel |> filter(target_has_pos_tm1, !fighting_at_ref) |> anti_join(alts_s2, by = c("ccode1", "ccode2", "year"))) == 0)
  ev_rows <- ev_in_panel |> inner_join(alts_s2 |> select(ccode1, ccode2, year, dy_mid_tm1), by = c("ccode1", "ccode2", "year"))
  rm(alts_s2)
  blog("  dy_mid_tm1 h=%d: %d of %d alternative rows = 1 (%.3f%%); chosen events (mid) with dy_mid_tm1 = 1: %d of %d; war events: %d of %d",
       h, sum(alts$dy_mid_tm1), nrow(alts), 100 * mean(alts$dy_mid_tm1), sum(ev_rows$dy_mid_tm1), nrow(ev_rows),
       sum(ev_rows$dy_mid_tm1[ev_rows$in_war]), sum(ev_rows$in_war))
  if (length(unique(alts$dy_mid_tm1)) < 2) stop(sprintf("h=%d: dy_mid_tm1 is constant over the alternative rows", h))
  alts <- alts |> arrange(year, ccode1, ccode2) |> mutate(h = h)
  assert_no_dup_keys(alts, c("ccode1", "ccode2", "year"), sprintf("alternatives h=%d", h))
  if (anyNA(alts$log_capdist)) blog("  alternatives without a capital distance at t-1: %d (NA kept; --distance runs drop them)", sum(is.na(alts$log_capdist)))
  assert_no_na(alts, c("dM", "dR_pr", "dR_eig", "dM_j_h", "dRpr_j_h", "dch_j_missing", "contig_land", "contig_any", "alliance_def_off", "joint_dem",
                       "major_count", "dy_prior_mids_h", "dy_peace_years_h", "dy_mid_tm1", "polrel", "in_iv", "in_s2"), "alternatives")
  assert_finite(alts, c("dM", "dR_pr", "dM_j_h", "dRpr_j_h"), sprintf("alternatives h=%d", h))
  ## the target change must agree with the initiator change of the same state-year wherever the target is also a panel row
  chk_j <- alts |> select(ccode2, year, dM_j_h, dRpr_j_h, dch_j_missing) |> distinct() |>
    inner_join(panel |> select(ccode2 = ccode, year, dM_h, dRpr_h), by = c("ccode2", "year"))
  stopifnot(nrow(chk_j) > 0, all(chk_j$dch_j_missing == 0L), all(abs(chk_j$dM_j_h - chk_j$dM_h) < 1e-12), all(abs(chk_j$dRpr_j_h - chk_j$dRpr_h) < 1e-12))
  ## directedness check: wherever (j,i,t) is also a row, dM(i,j) = -dM(j,i). The reverse row
  ## exists only when j is itself a panel state-year at t (has positions at t-1-h), so coverage < 100%.
  chk <- alts |> select(ccode1, ccode2, year, dM, dR_pr) |>
    inner_join(alts |> select(ccode1 = ccode2, ccode2 = ccode1, year, dM_rev = dM, dR_rev = dR_pr), by = c("ccode1", "ccode2", "year"))
  stopifnot(nrow(chk) > 0, all(abs(chk$dM + chk$dM_rev) < 1e-10), all(abs(chk$dR_pr + chk$dR_rev) < 1e-10))
  blog("  directedness check passed on %d of %d rows that have a reverse row", nrow(chk), nrow(alts))
  for (variant in c("all", "polrel")) {
    a_v <- if (variant == "all") alts else alts |> filter(polrel == 1L)
    ## every panel state-year must keep a non-empty choice set. A state-year with no row at all (possible in polrel for an
    ## isolated non-major state) is logged and dropped by assemble_cell(); a state-year that HAS rows but none in the stage-1
    ## IV set, or an initiating state-year whose stage-2 set is empty or holds only its chosen target, stops the build
    ## (review coverage gap: empty IV set / chosen-only stage-2 set)
    empty <- panel |> anti_join(a_v, by = c("ccode" = "ccode1", "year"))
    if (nrow(empty)) blog("  variant %s: %d panel state-years have NO alternative row (assemble_cell drops them)", variant, nrow(empty))
    n_sets <- a_v |> group_by(ccode1, year) |> summarise(n_iv = sum(in_iv), n_s2 = sum(in_s2), .groups = "drop")
    if (any(n_sets$n_iv == 0L)) stop(sprintf("h=%d variant %s: %d state-years have alternative rows but an EMPTY stage-1 IV set", h, variant, sum(n_sets$n_iv == 0L)))
    ev_v <- ev_flag |> semi_join(panel, by = c("ccode1" = "ccode", "year")) |> semi_join(a_v |> filter(in_s2 == 1L), by = c("ccode1", "ccode2", "year")) |>
      distinct(ccode1, year) |> inner_join(n_sets, by = c("ccode1", "year"))
    if (any(ev_v$n_s2 < 2L)) stop(sprintf("h=%d variant %s: %d initiating state-years have a stage-2 set that is empty or holds only the chosen target", h, variant, sum(ev_v$n_s2 < 2L)))
    blog("  variant %s: %d state-years with rows; IV-set size min %d / median %d; stage-2 set size min %d among the %d initiating state-years with a stage-2 event; empty IV sets 0; chosen-only stage-2 sets 0",
         variant, nrow(n_sets), min(n_sets$n_iv), as.integer(median(n_sets$n_iv)), min(ev_v$n_s2), nrow(ev_v))
    lost_pr <- ev_flag |> semi_join(panel, by = c("ccode1" = "ccode", "year")) |> semi_join(alts |> filter(in_s2 == 1L), by = c("ccode1", "ccode2", "year")) |>
      anti_join(a_v |> filter(in_s2 == 1L), by = c("ccode1", "ccode2", "year")) |> nrow()
    if (variant == "polrel") blog("  variant polrel: %d events have a target outside the politically relevant set (dropped from stage 2)", lost_pr)
    add_manifest(a_v, sprintf("alts_lag%d_%s.rds", h, variant), sprintf("choice sets h=%d variant=%s; union of the stage-1 IV set (in_iv, January rule) and the stage-2 set (in_s2, reference-date rule), already-fighting alternatives removed (month rule); dy_mid_tm1", h, variant))
  }
  rm(alts, chk); invisible(gc())
}
## the pooled scaling constants of the change variables, one row per (lag, variable)
scale_tbl <- bind_rows(scale_rows)
stopifnot(nrow(scale_tbl) == 2L * length(LAGS), all(is.finite(scale_tbl$pooled_sd)), all(scale_tbl$pooled_sd > 0))
## the two gap constants (lag-independent) are stored in the same file with h = NA
scale_tbl <- bind_rows(scale_tbl, data.frame(h = NA_integer_, variable = c("gap_M", "gap_R_pr"), pooled_sd = c(GAP_SCALE_M, GAP_SCALE_R)))
write.csv(scale_tbl, file.path(out_dir, "change_scale.csv"), row.names = FALSE)
blog("change_scale.csv written: %s", paste(sprintf("h=%d %s=%.6f", scale_tbl$h, scale_tbl$variable, scale_tbl$pooled_sd), collapse = "; "))
write_manifest(manifest, file.path(out_dir, "manifest.csv"))
blog("manifest written with %d files; build complete", length(manifest))
