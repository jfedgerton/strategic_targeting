## ============================================================================
## assemble_cell.R -- assemble the sy / alts / ev frames for one real-data cell
## Used by 02_fit_model.R (fits), run_measurement_real.R (measurement tests) so that every
## consumer sees the SAME rows (audit F3/F8: one assembly path, manifest-checked inputs).
## Returns list(sy, alts, ev, W_cols, X_cols, Z_cols, state_codes, x_ref, n_full, sampled, log lines)
## ============================================================================
## Control sets (PI decision 2026-10-24, minimal controls). Major-power status and alliance are components of material
## and relational power, so they are bad controls and left the main specification; the main sets are as small as the
## design allows.
##   W (stage 1), main : dM_i, dR_i, ongoing_tm1, log1p(monadic prior MIDs) [CONFLICT_HISTORY], dem_tm1
##   X (mu shifters), main : dM_i, dR_i only (no monadic controls)
##   Z (dyadic), main : land contiguity [CONTIGUITY], joint democracy, log1p(dyadic prior MIDs) [CONFLICT_HISTORY],
##                      dy_mid_tm1, [log1p dyadic peace years, DY_PEACE_YEARS], [log1p capital distance, DISTANCE],
##                      the target's power changes dM_j_h, dRpr_j_h, [dch_j_missing when non-constant, TARGET_MISSING]
## FULL_CONTROLS = TRUE (fit option --full_controls, cell suffix _fullctrl) reproduces EXACTLY the base specification of
## commit 8bf0ce1 (`git show 8bf0ce1:v2/R/assemble_cell.R`; Codex review P2-4, 2026-10-24), same columns in the same order:
##   W = X = dM_i, dR_i, cowmaj_tm1, dem_tm1, [log_prior_mids, log_peace_years : CONFLICT_HISTORY]        (P1 = P = 6)
##   Z     = [contig_land : CONTIGUITY], alliance_def_off, joint_dem, major_count,
##           [log_dy_prior_mids, log_dy_peace_years : CONFLICT_HISTORY], [log_capdist : DISTANCE],
##           dM_j_h, dRpr_j_h, [dch_j_missing : TARGET_MISSING]                                           (Pz = 8 or 9)
## The base had no ongoing_tm1 in W / X and no dy_mid_tm1 in Z (both were added on 2026-10-24, after 8bf0ce1), and its
## dyadic peace years followed CONFLICT_HISTORY. Under FULL_CONTROLS, DY_PEACE_YEARS therefore adds nothing when
## CONFLICT_HISTORY is TRUE (the column is already there) and restores log_dy_peace_years at its base position when
## CONFLICT_HISTORY is FALSE. The lists below are hard-coded per set (no shared abstraction) so that each one can be
## read off and tested against the base commit. w_controls() / x_controls() return the monadic controls (the columns
## after dM_i, dR_i) so that assemble_cell(), benchmarks.R and the tests share one list.
w_controls <- function(CONFLICT_HISTORY = TRUE, FULL_CONTROLS = FALSE) {
  if (FULL_CONTROLS) c("cowmaj_tm1", "dem_tm1", if (CONFLICT_HISTORY) c("log_prior_mids", "log_peace_years"))
  else c("ongoing_tm1", if (CONFLICT_HISTORY) "log_prior_mids", "dem_tm1")
}
x_controls <- function(CONFLICT_HISTORY = TRUE, FULL_CONTROLS = FALSE) {
  if (FULL_CONTROLS) c("cowmaj_tm1", "dem_tm1", if (CONFLICT_HISTORY) c("log_prior_mids", "log_peace_years"))
  else character(0)
}
z_columns <- function(CONFLICT_HISTORY = TRUE, DISTANCE = FALSE, CONTIGUITY = TRUE, TARGET_MISSING = FALSE, DY_PEACE_YEARS = FALSE,
                      FULL_CONTROLS = FALSE) {
  if (FULL_CONTROLS) {
    ## the base Z of 8bf0ce1, verbatim order
    c(if (CONTIGUITY) "contig_land", "alliance_def_off", "joint_dem", "major_count",
      if (CONFLICT_HISTORY) c("log_dy_prior_mids", "log_dy_peace_years") else if (DY_PEACE_YEARS) "log_dy_peace_years",
      if (DISTANCE) "log_capdist",
      "dM_j_h", "dRpr_j_h",
      if (TARGET_MISSING) "dch_j_missing")
  } else {
    c(if (CONTIGUITY) "contig_land", "joint_dem",
      if (CONFLICT_HISTORY) "log_dy_prior_mids",
      "dy_mid_tm1",
      if (DY_PEACE_YEARS) "log_dy_peace_years",
      if (DISTANCE) "log_capdist",
      "dM_j_h", "dRpr_j_h",
      if (TARGET_MISSING) "dch_j_missing")
  }
}
## continuous Z columns that assemble_cell() centres over the assembled alternative rows (the 0/1 indicators, dy_mid_tm1
## included, stay as is)
Z_CENTRED <- c("log_dy_prior_mids", "log_dy_peace_years", "log_capdist", "dM_j_h", "dRpr_j_h")
## Monadic 0/1 indicators available in the panel (major power, democracy, ongoing_tm1 = participant in any MID active at
## t-1). Those that enter W or X are centred like the conflict-history counts; the main W uses dem_tm1 and ongoing_tm1.
SY_INDICATORS <- c("cowmaj_tm1", "dem_tm1", "ongoing_tm1")
## manifest: NULL reads data_dir/manifest.csv now; a data.frame is a SNAPSHOT read earlier by the caller (02_fit_model.R passes
## the manifest it read at load, so a rebuild between that read and the loading of the inputs stops the run; review P2-5).
assemble_cell <- function(SAMPLE, LAG, VARIANT, CENT = "pr", CONFLICT_HISTORY = TRUE, SAMPLED_M = 0L,
                          DISTANCE = FALSE, CONTIGUITY = TRUE, DY_PEACE_YEARS = FALSE, DROP_WW2 = FALSE, FULL_CONTROLS = FALSE,
                          data_dir = file.path(v2_root(), "data"), manifest = NULL) {
  stopifnot(SAMPLE %in% c("mid", "fmid", "fmidcow", "fmidub", "war"), LAG %in% c(1L, 3L, 5L, 10L), VARIANT %in% c("all", "polrel"), CENT %in% c("pr", "eig"))
  manifest_path <- file.path(data_dir, "manifest.csv")
  stopifnot(file.exists(manifest_path))
  if (is.null(manifest)) { manifest <- read.csv(manifest_path, stringsAsFactors = FALSE); manifest_source <- "read now" }
  else { stopifnot(is.data.frame(manifest), all(c("file", "nrow", "sha256") %in% names(manifest))); manifest_source <- "snapshot passed by the caller" }
  ## ---- load inputs, each asserted against the manifest (snapshot) -------------------------------------
  panel  <- assert_matches_manifest(file.path(data_dir, sprintf("stage1_panel_lag%d.rds", LAG)), manifest)
  alts   <- assert_matches_manifest(file.path(data_dir, sprintf("alts_lag%d_%s.rds", LAG, VARIANT)), manifest)
  events <- assert_matches_manifest(file.path(data_dir, "events.rds"), manifest)
  stopifnot(all(panel$h == LAG), all(alts$h == LAG))
  stopifnot(all(c("escalation", "is_joiner_pair", "is_ww2", "target_has_pos_tm1", "fighting_at_ref", "dispnums", "n_disputes", "war_disp",
                  "in_mid_noww2", "in_fmid_noww2", "in_fmidcow_noww2", "in_fmidub_noww2", "in_war_noww2") %in% names(events)),
            all(c("ongoing_tm1") %in% names(panel)), all(c("dy_mid_tm1", "in_iv", "in_s2") %in% names(alts)))
  stopifnot(!anyNA(alts$in_iv), !anyNA(alts$in_s2), all(alts$in_iv %in% 0:1), all(alts$in_s2 %in% 0:1), all(alts$in_iv == 1L | alts$in_s2 == 1L))
  assert_no_dup_keys(events |> filter(!escalation), c("initiator", "target", "year"), "events (one per initiator, target, year after the collapse)")
  logf("inputs verified against the manifest (%s): panel %d rows, alts %d rows (in both sets %d, IV only %d, stage 2 only %d), events %d rows (%d escalations, %d joiner pairs, %d WWII, %d collapsed from more than one dispute, %d fighting at the reference date)",
       manifest_source, nrow(panel), nrow(alts), sum(alts$in_iv == 1L & alts$in_s2 == 1L), sum(alts$in_iv == 1L & alts$in_s2 == 0L), sum(alts$in_iv == 0L & alts$in_s2 == 1L),
       nrow(events), sum(events$escalation), sum(events$is_joiner_pair), sum(events$is_ww2), sum(events$n_disputes > 1L), sum(events$fighting_at_ref))

  ## ---- stage-1 panel -> sy ------------------------------------------------------------------
  y1_col <- paste0("y1_", SAMPLE); ev_flag <- paste0("in_", SAMPLE)
  dR_col <- if (CENT == "pr") "dRpr_h" else "dReig_h"
  dRalt_col <- if (CENT == "pr") "dR_pr" else "dR_eig"
  ## drop panel state-years with no alternative row at all in this variant (logged by the build). A state-year that has
  ## rows but none in the stage-1 IV set is an error, not a drop.
  has_cs <- panel |> semi_join(alts, by = c("ccode" = "ccode1", "year"))
  logf("panel state-years with a non-empty choice set: %d of %d", nrow(has_cs), nrow(panel))
  n_iv_sy <- alts |> group_by(ccode1, year) |> summarise(n_iv = sum(in_iv), .groups = "drop")
  if (any(n_iv_sy$n_iv == 0L)) stop(sprintf("%d state-years have alternative rows but an EMPTY stage-1 IV set (no row with in_iv = 1); rebuild the data", sum(n_iv_sy$n_iv == 0L)))
  sy <- has_cs |> arrange(year, ccode) |>
    mutate(y1 = .data[[y1_col]],
           log_prior_mids = log1p(mon_prior_mids_h), log_peace_years = log1p(mon_peace_years_h),
           dM_i = dM_h, dR_i = .data[[dR_col]])
  ## state index: union of initiators (panel) and targets (alts), sorted by ccode (index only; anchors are NOT taken from this order)
  state_codes <- sort(unique(c(sy$ccode, alts$ccode2)))
  sy$state_idx <- match(sy$ccode, state_codes)
  ctrl_w <- w_controls(CONFLICT_HISTORY, FULL_CONTROLS); ctrl_x <- x_controls(CONFLICT_HISTORY, FULL_CONTROLS)
  ctrl_sy <- union(ctrl_w, ctrl_x)
  ## centre every stage-1 / mu control that enters W or X (the 0/1 indicators as well as the conflict-history counts; PI
  ## decisions 2026-10-23 / 2026-10-24) so that alpha0 in stage 1 is the log-odds at zero power change and average
  ## controls, and (under --full_controls) the mu intercept mu0 is the ideal gap at those averages. In the main
  ## specification X has no controls, so mu0 is the ideal gap of a state with zero power change. dM_i and dR_i are NOT
  ## centred (zero power change is the reference).
  stopifnot(all(ctrl_sy %in% names(sy)), all(ctrl_x %in% ctrl_w), setequal(ctrl_x, if (FULL_CONTROLS) ctrl_w else character(0)))
  for (cc in SY_INDICATORS) stopifnot(!anyNA(sy[[cc]]), all(sy[[cc]] %in% c(0, 1)))
  logf("ongoing_tm1: mean %.4f over %d state-years (%d in a MID active at t-1)", mean(sy$ongoing_tm1), nrow(sy), sum(sy$ongoing_tm1))
  for (cc in ctrl_sy) sy[[cc]] <- sy[[cc]] - mean(sy[[cc]])
  stopifnot(all(abs(colMeans(as.matrix(sy[, ctrl_sy, drop = FALSE]))) < 1e-8))
  W_cols <- c("dM_i", "dR_i", ctrl_w)      # PI decision 2026-10-22: initiator power changes also enter stage 1 directly
  X_cols <- c("dM_i", "dR_i", ctrl_x)      # main: power changes only; --full_controls: the same controls as W
  logf("control set: %s; W = %s; X = %s", if (FULL_CONTROLS) "FULL (base specification of commit 8bf0ce1)" else "MINIMAL (main)", paste(W_cols, collapse = ","), paste(X_cols, collapse = ","))
  
  ## ---- alternatives -----------------------------------------------------------------------------
  sy_key <- sy |> transmute(ccode1 = ccode, year, sy_row = row_number())
  alts <- alts |> inner_join(sy_key, by = c("ccode1", "year")) |> arrange(sy_row, ccode2) |>
    mutate(state_idx = match(ccode1, state_codes), target_idx = match(ccode2, state_codes), dR = .data[[dRalt_col]],
           log_dy_prior_mids = log1p(dy_prior_mids_h), log_dy_peace_years = log1p(dy_peace_years_h))
  stopifnot(!anyNA(alts$target_idx), !anyNA(alts$sy_row))
  if (DISTANCE) {
    n_na <- sum(is.na(alts$log_capdist))
    if (n_na > 0) { logf("--distance: dropping %d alternatives without a capital distance", n_na); alts <- alts |> filter(!is.na(log_capdist)) }
  }
  ## target power changes (PI decision 2026-10-23): dM_j_h / dRpr_j_h are always in Z. The missing-target indicator
  ## enters only when it is non-constant in THIS assembled cell (a constant column has no coefficient to estimate).
  stopifnot(all(c("dM_j_h", "dRpr_j_h", "dch_j_missing") %in% names(alts)), all(alts$dch_j_missing %in% c(0L, 1L)),
            all(alts$dM_j_h[alts$dch_j_missing == 1L] == 0), all(alts$dRpr_j_h[alts$dch_j_missing == 1L] == 0))
  n_miss <- sum(alts$dch_j_missing); target_missing <- n_miss > 0 && n_miss < nrow(alts)
  logf("target changes: %d of %d alternative rows have a target without a change value (set to 0); dch_j_missing %s Z (%s)",
       n_miss, nrow(alts), if (target_missing) "enters" else "is NOT in",
       if (target_missing) "non-constant in this cell" else "constant in this cell")
  Z_cols <- z_columns(CONFLICT_HISTORY, DISTANCE, CONTIGUITY, TARGET_MISSING = target_missing, DY_PEACE_YEARS = DY_PEACE_YEARS, FULL_CONTROLS = FULL_CONTROLS)
  stopifnot(all(Z_cols %in% names(alts)))
  ## centre the continuous Z columns (conflict history, distance, target changes). The centring constants are the means over
  ## the stage-1 IV rows (in_iv = 1) ONLY (Codex re-check P2, 2026-10-24): the IV set is defined at January of t and does
  ## not depend on year-t events, whereas the stage-2-only rows do. Centring over the whole union block let a stage-2-only
  ## row move the means, which shifts every utility by a common constant: harmless in the stage-2 softmax but it enters
  ## stage 1 through tau * IV. The same constant is subtracted from every row (IV and stage-2 alike), so stage-2
  ## likelihoods are unchanged.
  z_cent <- intersect(Z_cols, Z_CENTRED)
  iv_rows <- alts$in_iv == 1L
  stopifnot(any(iv_rows))
  for (cc in z_cent) alts[[cc]] <- alts[[cc]] - mean(alts[[cc]][iv_rows])
  stopifnot(all(abs(colMeans(as.matrix(alts[iv_rows, z_cent, drop = FALSE]))) < 1e-8))
  logf("Z centring: %d continuous Z columns centred over the %d stage-1 IV rows (of %d assembled rows)", length(z_cent), sum(iv_rows), nrow(alts))
  ## dyadic lagged DV (PI decision 2026-10-24): 0/1, not centred, must vary within the assembled cell when it is in Z (it is
  ## not under --full_controls, whose Z is the base specification of 8bf0ce1)
  stopifnot(!anyNA(alts$dy_mid_tm1), all(alts$dy_mid_tm1 %in% c(0L, 1L)))
  if ("dy_mid_tm1" %in% Z_cols && length(unique(alts$dy_mid_tm1)) < 2)
    stop(sprintf("dy_mid_tm1 is constant (%d) over the %d assembled alternative rows of this cell (sample %s, lag %d, variant %s); it cannot enter Z",
                 alts$dy_mid_tm1[1], nrow(alts), SAMPLE, LAG, VARIANT))
  logf("choice-set flags: %d alternative rows; in the stage-1 IV set %d, in the stage-2 set %d, in both %d; state-years whose two sets differ %d of %d",
       nrow(alts), sum(alts$in_iv), sum(alts$in_s2), sum(alts$in_iv == 1L & alts$in_s2 == 1L),
       nrow(alts |> group_by(sy_row) |> summarise(d = any(in_iv != in_s2), .groups = "drop") |> filter(d)), nrow(sy))

  ## ---- events in this sample -------------------------------------------------------------------------
  ## The sample flag in_<sample> already excludes escalations (events whose pair was fighting in another dispute at t);
  ## asserted here, not assumed. --drop_ww2 (appendix) uses the build's in_<sample>_noww2 flag (review P2-3): the sample flag
  ## of the non-WWII disputes of the (i, j, t) alone, so an event collapsed from a WWII and a non-WWII dispute stays when the
  ## non-WWII dispute qualifies (the old rule in_<sample> & !is_ww2 dropped it). Their stage-1 initiations go with them.
  ev_all <- events |> filter(.data[[ev_flag]])
  stopifnot(!any(ev_all$escalation))
  n_ww2 <- sum(ev_all$is_ww2)
  if (DROP_WW2) {
    ev_flag_noww2 <- paste0(ev_flag, "_noww2")
    n_before <- nrow(ev_all)
    ev_all <- events |> filter(.data[[ev_flag_noww2]])
    stopifnot(!any(ev_all$escalation), all(events[[ev_flag]][events[[ev_flag_noww2]]]), !any(ev_all$is_ww2 & ev_all$n_disputes == 1L))
    logf("--drop_ww2: sample %s %d -> %d events using %s (%d events flagged is_ww2 removed; %d events collapsed from a WWII and a non-WWII dispute kept because a non-WWII dispute qualifies, whereas in_%s & !is_ww2 would keep %d events); WWII disputes %s",
         SAMPLE, n_before, nrow(ev_all), ev_flag_noww2, n_before - nrow(ev_all), sum(ev_all$is_ww2), SAMPLE, sum(events[[ev_flag]] & !events$is_ww2),
         paste(sort(unique(events$dispnum[events$is_ww2])), collapse = ","))
  } else logf("WWII events in sample %s: %d of %d (kept; --drop_ww2=TRUE removes them)", SAMPLE, n_ww2, nrow(ev_all))
  ## stage-1 outcome recounted from the kept events (independent path through the escalation and WWII columns). Without
  ## --drop_ww2 it must equal the panel's y1_<sample>; with it the recount IS y1, so WWII-only state-years get y1 = 0.
  y1_re <- ev_all |> distinct(ccode = initiator, year) |> mutate(.y1 = 1L)
  sy <- sy |> left_join(y1_re, by = c("ccode", "year")) |> mutate(.y1 = coalesce(.y1, 0L))
  if (!DROP_WW2) {
    if (!all(sy$y1 == sy$.y1)) stop(sprintf("y1_%s disagrees with the recount from non-escalation events on %d state-years", SAMPLE, sum(sy$y1 != sy$.y1)))
  } else {
    logf("--drop_ww2: stage-1 initiations %d -> %d (%d state-years whose only %s events were WWII now have y1 = 0)", sum(sy$y1), sum(sy$.y1), sum(sy$y1 == 1L & sy$.y1 == 0L), SAMPLE)
    stopifnot(all(sy$.y1 <= sy$y1))
    sy$y1 <- sy$.y1
  }
  sy <- sy |> select(-.y1)
  ## one chosen alternative per (i, j, t): same-year duplicates were collapsed by the build (PI decision 2026-10-24)
  ev_keys <- ev_all |> transmute(ccode1 = initiator, ccode2 = target, year, dispnum, dispnums, target_has_pos_tm1, fighting_at_ref)
  assert_no_dup_keys(ev_keys, c("ccode1", "ccode2", "year"), sprintf("events of sample %s", SAMPLE))
  ## an event whose pair was already fighting at the state-year's reference date (its earlier dispute ended before the
  ## event) has no STAGE-2 row by construction and enters stage 1 only; its (i, j, t) may keep an IV-only row (in_s2 = 0)
  far_keys <- ev_keys |> filter(fighting_at_ref)
  if (nrow(far_keys |> semi_join(alts |> filter(in_s2 == 1L), by = c("ccode1", "ccode2", "year"))) > 0) stop("fighting_at_ref events have stage-2 rows (in_s2 = 1); rebuild the data")
  logf("events in sample %s flagged fighting_at_ref (pair already fighting at the reference date; stage 1 only): %d of %d (%d keep an IV-only row)", SAMPLE, nrow(far_keys), nrow(ev_keys),
       nrow(far_keys |> semi_join(alts, by = c("ccode1", "ccode2", "year"))))
  ## the chosen alternative must be a STAGE-2 row of the state-year's block (alt_row indexes the full block)
  ev <- ev_keys |> inner_join(sy_key, by = c("ccode1", "year")) |>
    inner_join(alts |> mutate(alt_row = row_number()) |> filter(in_s2 == 1L) |> select(ccode1, ccode2, year, alt_row, dy_mid_tm1), by = c("ccode1", "ccode2", "year"))
  logf("events in sample %s: %d (from %d dispute pairs); usable (initiator in panel AND target in the stage-2 choice set): %d", SAMPLE, nrow(ev_keys),
       sum(lengths(strsplit(ev_keys$dispnums, ";", fixed = TRUE))), nrow(ev))
  stopifnot(nrow(ev) > 0, all(sy$y1[ev$sy_row] == 1L), all(alts$in_s2[ev$alt_row] == 1L), all(alts$sy_row[ev$alt_row] == ev$sy_row))
  assert_no_dup_keys(ev, c("sy_row", "alt_row"), "events (one chosen alternative per state-year and target)")
  ## every event whose chosen alternative can exist (initiator in sy, target with positions at t-1 and in the system at t,
  ## pair not fighting at the reference date) has it in its stage-2 set; in the "all" variant without the distance filter
  ## this is exact, elsewhere a lower bound
  expected <- ev_keys |> semi_join(sy_key, by = c("ccode1", "year")) |> filter(target_has_pos_tm1, !fighting_at_ref)
  if (VARIANT == "all" && !DISTANCE) {
    if (nrow(ev) != nrow(expected)) stop(sprintf("%d events with initiator in the panel, a target with positions at t-1 and a pair not fighting at the reference date are not in their stage-2 choice set (expected %d usable, got %d)", abs(nrow(expected) - nrow(ev)), nrow(expected), nrow(ev)))
  } else stopifnot(nrow(ev) <= nrow(expected))
  stopifnot(nrow(ev_keys |> semi_join(sy_key, by = c("ccode1", "year")) |> filter(!target_has_pos_tm1) |> semi_join(alts, by = c("ccode1", "ccode2", "year"))) == 0)
  ## no initiating state-year with a usable event may have an empty or chosen-only stage-2 set (review coverage gap); the
  ## empty-IV-set check is above, so a state-year here always has >= 1 IV row as well
  n_s2_sy <- alts |> group_by(sy_row) |> summarise(n_s2 = sum(in_s2), .groups = "drop")
  s2_small <- n_s2_sy |> filter(sy_row %in% ev$sy_row, n_s2 < 2L)
  if (nrow(s2_small) > 0) stop(sprintf("%d initiating state-years have a stage-2 set that is empty or holds only the chosen target (sample %s, lag %d, variant %s)", nrow(s2_small), SAMPLE, LAG, VARIANT))
  logf("stage-2 sets of the %d initiating state-years with a usable event: size min %d, median %d; IV sets of all %d state-years: size min %d, median %d; empty IV sets 0, chosen-only stage-2 sets 0",
       n_distinct(ev$sy_row), min(n_s2_sy$n_s2[n_s2_sy$sy_row %in% ev$sy_row]), as.integer(median(n_s2_sy$n_s2[n_s2_sy$sy_row %in% ev$sy_row])),
       nrow(sy), min(n_iv_sy$n_iv), as.integer(median(n_iv_sy$n_iv)))
  logf("dy_mid_tm1: %d of %d alternative rows = 1; chosen events with dy_mid_tm1 = 1: %d of %d", sum(alts$dy_mid_tm1), nrow(alts), sum(ev$dy_mid_tm1), nrow(ev))
  n_orphan <- sum(sy$y1 == 1L & !(seq_len(nrow(sy)) %in% ev$sy_row))
  logf("state-years with y1 = 1 but no usable stage-2 event (kept in stage 1 only): %d", n_orphan)
  ev <- ev |> select(sy_row, alt_row)
  
  ## ---- optional sampled alternatives ------------------------------------------------------------------
  n_full <- NULL; sampled <- 0L
  if (SAMPLED_M > 0) {
    sa <- sample_alternatives(alts, ev, m = SAMPLED_M, seed = V2_SEED)
    alts <- sa$alts; ev <- sa$ev; n_full <- sa$n_full; sampled <- 1L
    logf("sampled alternatives: m = %d -> %d alternative rows", SAMPLED_M, nrow(alts))
    if ("dch_j_missing" %in% Z_cols && length(unique(alts$dch_j_missing)) < 2)
      stop("dch_j_missing is in Z but became constant after sampling alternatives; rerun with a larger --sampled_m or full choice sets")
  }
  
  ## ---- reference X for mu_state: mean X over INITIATING state-years in this sample -----------------------
  x_ref <- colMeans(as.matrix(sy[sy$y1 == 1L, X_cols, drop = FALSE]))
  logf("x_ref (mean X over initiating state-years): %s", paste(sprintf("%s=%.3f", X_cols, x_ref), collapse = ", "))
  list(sy = sy, alts = alts, ev = ev, W_cols = W_cols, X_cols = X_cols, Z_cols = Z_cols, state_codes = state_codes,
       x_ref = x_ref, n_full = n_full, sampled = sampled, n_orphan = n_orphan, manifest_path = manifest_path,
       n_ww2 = n_ww2, drop_ww2 = DROP_WW2, dy_peace_years = DY_PEACE_YEARS, full_controls = FULL_CONTROLS, n_fighting_at_ref = nrow(far_keys))
}
