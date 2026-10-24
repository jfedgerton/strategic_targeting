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
##   events.rds                    one row per (dispnum, initiator, target) directed event
##                                 with sample flags in_mid / in_fmid / in_fmidcow / in_fmidub / in_war (+ fatal_source)
##   changes_lag{h}.rds            ONE state-level change table per lag (build_change_table() in 00_helpers.R):
##                                 raw change dM_raw_h / dRpr_raw_h = level(t-1) - level(t-1-h), and the transformed
##                                 change dM_h / dRpr_h = (raw - year mean) / one pooled SD per lag. Feeds BOTH the
##                                 initiator change in stage 1 and the target change in stage 2.
##   change_scale.csv              the pooled SD constants per lag and variable
##   stage1_panel_lag{h}.rds       state-year panel for h in LAGS: y1_{mid,fmid,fmidcow,fmidub,war}, W, X
##   alts_lag{h}_{all,polrel}.rds  flattened choice sets (one row per state-year x alternative); gaps dM / dR_pr are
##                                 log ratios; target changes dM_j_h / dRpr_j_h (+ dch_j_missing when the target has
##                                 no change value); log_capdist = log1p(capital distance km) for the --distance option
##   build_log.txt                 dropped years / rows and why
##   manifest.csv
##
## Timing conventions (audit F5): event onset year = t. Positions x_ij = (M_i - M_j, R_i - R_j)
## are measured at t-1. Power changes (initiator AND target) are level(t-1) - level(t-1-h), h in LAGS.
## ALL lags are calendar-year joins on (ccode, year - h); never row lag().
##
## Initiation coding (audit F4): initiator = Side A ORIGINATOR (sidea == 1, orig == 1);
## target = Side B ORIGINATOR. One direction per dispute pair. For multilateral disputes
## every (Side A originator x Side B originator) pair is one event. Joiners (orig == 0)
## are excluded in the main build; INCLUDE_JOINERS = TRUE adds them (robustness).
##
## Conflict-history controls (new in v2), all relative to the change window [t-1-h, t-1]:
##   (a) mon_prior_mids_h : # distinct MIDs the state was an originator in (either side)
##                          with onset in the 10 years BEFORE the window: [t-1-h-10, t-1-h-1]
##   (b) dy_prior_mids_h  : # distinct MIDs between i and j with onset in that same period
##   (c) dy_peace_years_h : years since the last MID-ongoing year between i and j, as of t-1-h
##   plus mon_peace_years_h (monadic analogue of (c)).
## ============================================================================
suppressPackageStartupMessages({
  library(peacesciencer); library(igraph); library(dplyr); library(tidyr); library(digest)
})
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(this_file) == 0) this_file <- "v2/R/01_build_data.R"
source(file.path(dirname(this_file), "00_helpers.R"))
Sys.setenv(V2_ROOT = normalizePath(file.path(dirname(this_file), "..")))
set.seed(V2_SEED)
log_open("01_build_data")
out_dir <- file.path(v2_root(), "data"); dir.create(out_dir, showWarnings = FALSE)
build_log <- file.path(out_dir, "build_log.txt"); cat("", file = build_log)
blog <- function(fmt, ...) { m <- logf(fmt, ...); cat(m, "\n", file = build_log, append = TRUE) }

## ---- flags (PI decisions; see README) ----------------------------------------------
LAGS            <- c(1L, 3L, 5L, 10L)
INCLUDE_JOINERS <- FALSE          # TRUE = robustness build that also counts joiners as initiators/targets
IGO_FILL        <- "drop"         # "drop": years with all-missing IGO data are dropped from positions
                                  # "locf": R at year y uses the latest IGO year <= y (igo_src_year recorded)
PR_CONT_MAX     <- 5L             # politically relevant: conttype in 1..PR_CONT_MAX OR either a major power
PREWINDOW_YEARS <- 10L            # conflict-history look-back before the change window starts
POLITY_DEM_CUT  <- 6L             # democracy = polity2 >= 6 (NA -> 0, flagged in *_dem_missing)
FIRST_YEAR      <- 1816L
blog("flags: LAGS=%s INCLUDE_JOINERS=%s IGO_FILL=%s PR_CONT_MAX=%d PREWINDOW_YEARS=%d",
     paste(LAGS, collapse = ","), INCLUDE_JOINERS, IGO_FILL, PR_CONT_MAX, PREWINDOW_YEARS)

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
## one row per (dispnum, ccode1, ccode2): the onset year of the pair within the dispute
pair <- dd |> group_by(dispnum, ccode1, ccode2) |>
  summarise(year = min(year), sidea1 = first(sidea1), sidea2 = first(sidea2),
            orig1 = first(orig1), orig2 = first(orig2),
            fatality1 = max(fatality1), fatality2 = max(fatality2),
            hostlev1 = max(hostlev1), hostlev2 = max(hostlev2), .groups = "drop") |>
  left_join(disp, by = "dispnum")
## A -> B direction only, originators on both sides (main) or any participants (robustness)
events <- pair |> filter(sidea1 == 1L, sidea2 == 0L)
if (!INCLUDE_JOINERS) events <- events |> filter(orig1 == 1L, orig2 == 1L)
events <- events |> rename(initiator = ccode1, target = ccode2) |>
  mutate(fatal1 = ifelse(fatality1 < 0, NA, fatality1 >= 1),
         fatal2 = ifelse(fatality2 < 0, NA, fatality2 >= 1),
         fatal_disp = ifelse(fatality_disp < 0, NA, fatality_disp >= 1),
         fatal = case_when(fatal1 %in% TRUE | fatal2 %in% TRUE | fatal_disp %in% TRUE ~ TRUE,
                           is.na(fatal1) & is.na(fatal2) & is.na(fatal_disp) ~ NA,
                           TRUE ~ FALSE),
         war = hostlev_disp == 5L,
         in_mid = TRUE, in_war = war %in% TRUE)
## fatal-MID coding (PI decision 2026-10-23). Main fmid: wars are fatal by definition (COW wars need
## 1,000+ battle deaths); events whose COW fatality codes are all missing are filled from MIC 1.0
## (Gibler et al., confrontation-level fatalmin >= 1; micnum = COW dispnum). MIC never overrides a
## COW code. fmidcow = the original COW-only coding (missing -> not fatal), kept for the appendix.
mic_path <- file.path(v2_root(), "data_raw", "mic-conf-1.0.csv")
stopifnot(file.exists(mic_path))
if (file_digest(mic_path) != MIC_CONF_SHA256) stop("data_raw/mic-conf-1.0.csv sha256 does not match MIC_CONF_SHA256 in 00_helpers.R")
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
events <- events |> left_join(mic, by = "dispnum") |> left_join(gml, by = c("dispnum", "initiator", "target")) |>
  mutate(fatal_mic = ifelse(is.na(mic_fatalmin), NA, mic_fatalmin >= 1),
         in_fmid = code_fatal_main(fatal, fatal_mic, war, fatal_gml),
         in_fmidcow = fatal %in% TRUE,
         in_fmidub = in_fmid | (is.na(fatal) & is.na(fatal_mic) & is.na(fatal_gml)),   # upper bound: uncoded -> fatal
         fatal_source = case_when(war %in% TRUE & !(fatal %in% TRUE) ~ "war_rule",
                                  !is.na(fatal) ~ "cow",
                                  !is.na(fatal_mic) ~ "mic_fill",
                                  !is.na(fatal_gml) ~ "gml_fill",
                                  TRUE ~ "missing_nonfatal"))
stopifnot(nrow(events) == n_ev_before, !anyNA(events$in_fmid),
          all(events$in_fmid[events$in_war]),          # war is a subset of fmid
          all(events$in_fmid[events$in_fmidcow]),      # COW-only coding is a subset of the main coding
          all(events$in_fmidub[events$in_fmid]),       # main coding is a subset of the upper bound
          sum(events$in_fmidub & !events$in_fmid) == sum(events$fatal_source == "missing_nonfatal"))
cow_na <- is.na(events$fatal)
blog("fatal coding: COW fatality missing for %d events; of these %d disputes in MIC (MIC fatal %d, MIC non-fatal %d), %d not in MIC",
     sum(cow_na), sum(cow_na & !is.na(events$fatal_mic)), sum(cow_na & events$fatal_mic %in% TRUE),
     sum(cow_na & events$fatal_mic %in% FALSE), sum(cow_na & is.na(events$fatal_mic)))
blog("fatal coding: wars added to fmid by the war rule %d; fmid (main) %d vs fmidcow (COW only) %d; fatal_source counts: %s",
     sum(events$fatal_source == "war_rule"), sum(events$in_fmid), sum(events$in_fmidcow),
     paste(names(table(events$fatal_source)), table(events$fatal_source), sep = "=", collapse = ", "))
blog("fatal coding: COW and MIC both missing for %d non-war events; GML codes %d of them (GML fatal %d); %d stay uncoded (non-fatal)",
     sum(cow_na & is.na(events$fatal_mic) & !events$in_war), sum(events$fatal_source == "gml_fill"),
     sum(events$fatal_source == "gml_fill" & events$in_fmid), sum(events$fatal_source == "missing_nonfatal"))
blog("fatal coding (info only, not used): events in MIC %d of %d; where both code fatality, COW and MIC disagree on %d events",
     sum(!is.na(events$fatal_mic)), nrow(events), sum(!is.na(events$fatal) & !is.na(events$fatal_mic) & events$fatal != events$fatal_mic))
## one direction per dispute pair: no (dispnum, j, i) when (dispnum, i, j) exists
rev_key <- paste(events$dispnum, events$target, events$initiator)
stopifnot(!any(rev_key %in% paste(events$dispnum, events$initiator, events$target)))
assert_no_dup_keys(events, c("dispnum", "initiator", "target"), "events")
stopifnot(all(events$initiator != events$target))
blog("events (Side A originator -> Side B originator, joiners %s): %d; fatal %d (fatal NA %d); war %d; disputes %d; multilateral pairs %d",
     if (INCLUDE_JOINERS) "INCLUDED" else "excluded", nrow(events), sum(events$in_fmid), sum(is.na(events$fatal)),
     sum(events$in_war), n_distinct(events$dispnum), sum(duplicated(events$dispnum)))
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
blog("dyadic MID history rows: %d", nrow(dy_disp))

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
## (gap_scale() is in 00_helpers.R so it can be unit-tested.)
GAP_SCALE_M <- gap_scale(pos$M, pos$year)
GAP_SCALE_R <- gap_scale(pos$R_pr, pos$year)
stopifnot(is.finite(GAP_SCALE_M), GAP_SCALE_M > 0, is.finite(GAP_SCALE_R), GAP_SCALE_R > 0)
blog("gap scale (one constant per axis, all years): M (log CINC ratio) %.6f; R_pr (log PageRank ratio) %.6f; 1 scaled unit = power ratio %.2f (M), %.2f (R_pr)",
     GAP_SCALE_M, GAP_SCALE_R, exp(GAP_SCALE_M), exp(GAP_SCALE_R))
manifest <- list()
add_manifest <- function(obj, fname, note) {
  p <- file.path(out_dir, fname); saveRDS(obj, p)
  manifest[[fname]] <<- data.frame(file = fname, nrow = nrow(obj), ncol = ncol(obj), sha256 = file_digest(p), note = note)
  blog("wrote %s: %d x %d", fname, nrow(obj), ncol(obj))
}
add_manifest(pos, "sy_positions.rds", "state-year positions: M = log CINC share, R_pr = log PageRank (raw log shares), R_eig std within year")
add_manifest(events, "events.rds", sprintf("directed events, joiners %s", if (INCLUDE_JOINERS) "included" else "excluded"))
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
  panel <- states |> select(ccode, year) |>
    inner_join(pos_k |> transmute(ccode, year = year + 1L, M_tm1 = M, Rpr_tm1 = R_pr, Reig_tm1 = R_eig), by = c("ccode", "year")) |>
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
  ## initiation indicators per sample
  ev_sy <- events |> group_by(ccode = initiator, year) |>
    summarise(n_ev_mid = n(), n_ev_fmid = sum(in_fmid), n_ev_fmidcow = sum(in_fmidcow), n_ev_fmidub = sum(in_fmidub),
              n_ev_war = sum(in_war), .groups = "drop")
  panel <- panel |> left_join(ev_sy, by = c("ccode", "year")) |>
    mutate(across(c(n_ev_mid, n_ev_fmid, n_ev_fmidcow, n_ev_fmidub, n_ev_war), ~ coalesce(.x, 0L)),
           y1_mid = as.integer(n_ev_mid > 0), y1_fmid = as.integer(n_ev_fmid > 0),
           y1_fmidcow = as.integer(n_ev_fmidcow > 0),
           y1_fmidub = as.integer(n_ev_fmidub > 0), y1_war = as.integer(n_ev_war > 0),
           log_n_system_tm1 = log(n_system_tm1)) |>
    arrange(year, ccode)
  assert_no_dup_keys(panel, c("ccode", "year"), sprintf("stage-1 panel h=%d", h))
  assert_no_na(panel, c("dM_h", "dRpr_h", "dM_raw_h", "dRpr_raw_h", "dReig_h", "M_tm1", "Rpr_tm1", "cowmaj_tm1", "dem_tm1",
                        "mon_prior_mids_h", "mon_peace_years_h", "y1_mid", "y1_fmid", "y1_fmidcow", "y1_fmidub", "y1_war"), "stage-1 panel")
  assert_finite(panel, c("dM_h", "dRpr_h", "dM_raw_h", "dRpr_raw_h", "M_tm1", "Rpr_tm1"), sprintf("stage-1 panel h=%d", h))
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
  ev_lost <- events |> anti_join(panel, by = c("initiator" = "ccode", "year"))
  blog("  events whose initiator state-year is NOT in the panel (no positions at t-1 / t-1-h): %d of %d", nrow(ev_lost), nrow(events))
  add_manifest(panel, sprintf("stage1_panel_lag%d.rds", h),
               sprintf("stage-1 panel, h=%d; dM_h / dRpr_h year-centred and scaled by pooled SD M=%.6f R_pr=%.6f", h, chg$scale$pooled_sd[1], chg$scale$pooled_sd[2]))

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
  ## event flags on alternative rows: is this (i, j, t) a chosen target in each sample?
  ev_flag <- events |> transmute(ccode1 = initiator, ccode2 = target, year, dispnum, in_mid, in_fmid, in_war)
  n_ev_target_lost <- ev_flag |> semi_join(panel, by = c("ccode1" = "ccode", "year")) |>
    anti_join(alts, by = c("ccode1", "ccode2", "year")) |> nrow()
  blog("  events whose target is NOT in the choice set (target lacks positions at t-1): %d", n_ev_target_lost)
  alts <- alts |> arrange(year, ccode1, ccode2) |> mutate(h = h)
  assert_no_dup_keys(alts, c("ccode1", "ccode2", "year"), sprintf("alternatives h=%d", h))
  if (anyNA(alts$log_capdist)) blog("  alternatives without a capital distance at t-1: %d (NA kept; --distance runs drop them)", sum(is.na(alts$log_capdist)))
  assert_no_na(alts, c("dM", "dR_pr", "dR_eig", "dM_j_h", "dRpr_j_h", "dch_j_missing", "contig_land", "contig_any", "alliance_def_off", "joint_dem",
                       "major_count", "dy_prior_mids_h", "dy_peace_years_h", "polrel"), "alternatives")
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
    ## every panel state-year must keep a non-empty choice set
    empty <- panel |> anti_join(a_v, by = c("ccode" = "ccode1", "year"))
    if (nrow(empty)) blog("  variant %s: %d panel state-years have an EMPTY choice set (02_fit_model drops them)", variant, nrow(empty))
    lost_pr <- ev_flag |> semi_join(panel, by = c("ccode1" = "ccode", "year")) |> semi_join(alts, by = c("ccode1", "ccode2", "year")) |>
      anti_join(a_v, by = c("ccode1", "ccode2", "year")) |> nrow()
    if (variant == "polrel") blog("  variant polrel: %d events have a target outside the politically relevant set (dropped from stage 2)", lost_pr)
    add_manifest(a_v, sprintf("alts_lag%d_%s.rds", h, variant), sprintf("choice sets h=%d variant=%s", h, variant))
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
