## ============================================================================
## assemble_cell.R -- assemble the sy / alts / ev frames for one real-data cell
## Used by 02_fit_model.R (fits), run_measurement_real.R (measurement tests) so that every
## consumer sees the SAME rows (audit F3/F8: one assembly path, manifest-checked inputs).
## Returns list(sy, alts, ev, W_cols, X_cols, Z_cols, state_codes, x_ref, n_full, sampled, log lines)
## ============================================================================
## Control sets (PI decision 2026-10-24, minimal controls). Major-power status and alliance are components of material
## and relational power, so they are bad controls and left the main specification; the main sets are as small as the
## design allows. FULL_CONTROLS = TRUE (fit option --full_controls, cell suffix _fullctrl) restores the previous main
## specification exactly (W and X with cowmaj_tm1, dem_tm1, ongoing_tm1, log prior MIDs and log monadic peace years; Z with
## alliance and the major-power count as well).
##   W (stage 1), main : dM_i, dR_i, ongoing_tm1, log1p(monadic prior MIDs) [CONFLICT_HISTORY], dem_tm1
##   X (mu shifters), main : dM_i, dR_i only (no monadic controls)
##   Z (dyadic), main : land contiguity [CONTIGUITY], joint democracy, log1p(dyadic prior MIDs) [CONFLICT_HISTORY],
##                      dy_mid_tm1, [log1p dyadic peace years, DY_PEACE_YEARS], [log1p capital distance, DISTANCE],
##                      the target's power changes dM_j_h, dRpr_j_h, [dch_j_missing when non-constant, TARGET_MISSING]
## The option flags combine with either control set the same way. --conflict_history=FALSE removes the prior-MID counts
## (and, under FULL_CONTROLS, the monadic peace years) from W, X and Z; --dy_peace_years, --distance and --contiguity act
## on Z only. w_controls() / x_controls() return the monadic controls (the columns after dM_i, dR_i) so that
## assemble_cell(), benchmarks.R and the tests share one list.
w_controls <- function(CONFLICT_HISTORY = TRUE, FULL_CONTROLS = FALSE) {
  if (FULL_CONTROLS) c(SY_INDICATORS, if (CONFLICT_HISTORY) c("log_prior_mids", "log_peace_years"))
  else c("ongoing_tm1", if (CONFLICT_HISTORY) "log_prior_mids", "dem_tm1")
}
x_controls <- function(CONFLICT_HISTORY = TRUE, FULL_CONTROLS = FALSE) {
  if (FULL_CONTROLS) w_controls(CONFLICT_HISTORY, TRUE) else character(0)
}
z_columns <- function(CONFLICT_HISTORY = TRUE, DISTANCE = FALSE, CONTIGUITY = TRUE, TARGET_MISSING = FALSE, DY_PEACE_YEARS = FALSE,
                      FULL_CONTROLS = FALSE) {
  c(if (CONTIGUITY) "contig_land", if (FULL_CONTROLS) "alliance_def_off", "joint_dem", if (FULL_CONTROLS) "major_count",
    if (CONFLICT_HISTORY) "log_dy_prior_mids",
    "dy_mid_tm1",
    if (DY_PEACE_YEARS) "log_dy_peace_years",
    if (DISTANCE) "log_capdist",
    "dM_j_h", "dRpr_j_h",
    if (TARGET_MISSING) "dch_j_missing")
}
## continuous Z columns that assemble_cell() centres over the assembled alternative rows (the 0/1 indicators, dy_mid_tm1
## included, stay as is)
Z_CENTRED <- c("log_dy_prior_mids", "log_dy_peace_years", "log_capdist", "dM_j_h", "dRpr_j_h")
## Monadic 0/1 indicators available in the panel (major power, democracy, ongoing_tm1 = participant in any MID active at
## t-1). Those that enter W or X are centred like the conflict-history counts; the main W uses dem_tm1 and ongoing_tm1.
SY_INDICATORS <- c("cowmaj_tm1", "dem_tm1", "ongoing_tm1")
assemble_cell <- function(SAMPLE, LAG, VARIANT, CENT = "pr", CONFLICT_HISTORY = TRUE, SAMPLED_M = 0L,
                          DISTANCE = FALSE, CONTIGUITY = TRUE, DY_PEACE_YEARS = FALSE, DROP_WW2 = FALSE, FULL_CONTROLS = FALSE,
                          data_dir = file.path(v2_root(), "data")) {
  stopifnot(SAMPLE %in% c("mid", "fmid", "fmidcow", "fmidub", "war"), LAG %in% c(1L, 3L, 5L, 10L), VARIANT %in% c("all", "polrel"), CENT %in% c("pr", "eig"))
  manifest_path <- file.path(data_dir, "manifest.csv")
  stopifnot(file.exists(manifest_path))
  ## ---- load inputs, each asserted against the manifest -------------------------------------
  panel  <- assert_matches_manifest(file.path(data_dir, sprintf("stage1_panel_lag%d.rds", LAG)), manifest_path)
  alts   <- assert_matches_manifest(file.path(data_dir, sprintf("alts_lag%d_%s.rds", LAG, VARIANT)), manifest_path)
  events <- assert_matches_manifest(file.path(data_dir, "events.rds"), manifest_path)
  stopifnot(all(panel$h == LAG), all(alts$h == LAG))
  stopifnot(all(c("escalation", "is_joiner_pair", "is_ww2", "target_has_pos_tm1", "fighting_at_ref", "dispnums", "n_disputes", "war_disp") %in% names(events)),
            all(c("ongoing_tm1") %in% names(panel)), all(c("dy_mid_tm1") %in% names(alts)))
  assert_no_dup_keys(events |> filter(!escalation), c("initiator", "target", "year"), "events (one per initiator, target, year after the collapse)")
  logf("inputs verified: panel %d rows, alts %d rows, events %d rows (%d escalations, %d joiner pairs, %d WWII, %d collapsed from more than one dispute, %d fighting at the reference date)",
       nrow(panel), nrow(alts), nrow(events), sum(events$escalation), sum(events$is_joiner_pair), sum(events$is_ww2), sum(events$n_disputes > 1L), sum(events$fighting_at_ref))
  
  ## ---- stage-1 panel -> sy ------------------------------------------------------------------
  y1_col <- paste0("y1_", SAMPLE); ev_flag <- paste0("in_", SAMPLE)
  dR_col <- if (CENT == "pr") "dRpr_h" else "dReig_h"
  dRalt_col <- if (CENT == "pr") "dR_pr" else "dR_eig"
  ## drop panel state-years with an empty choice set in this variant (logged by the build)
  has_cs <- panel |> semi_join(alts, by = c("ccode" = "ccode1", "year"))
  logf("panel state-years with a non-empty choice set: %d of %d", nrow(has_cs), nrow(panel))
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
  stopifnot(all(ctrl_sy %in% names(sy)), all(ctrl_x %in% ctrl_w))
  for (cc in SY_INDICATORS) stopifnot(!anyNA(sy[[cc]]), all(sy[[cc]] %in% c(0, 1)))
  logf("ongoing_tm1: mean %.4f over %d state-years (%d in a MID active at t-1)", mean(sy$ongoing_tm1), nrow(sy), sum(sy$ongoing_tm1))
  for (cc in ctrl_sy) sy[[cc]] <- sy[[cc]] - mean(sy[[cc]])
  stopifnot(all(abs(colMeans(as.matrix(sy[, ctrl_sy, drop = FALSE]))) < 1e-8))
  W_cols <- c("dM_i", "dR_i", ctrl_w)      # PI decision 2026-10-22: initiator power changes also enter stage 1 directly
  X_cols <- c("dM_i", "dR_i", ctrl_x)      # main: power changes only; --full_controls: the same controls as W
  logf("control set: %s; W = %s; X = %s", if (FULL_CONTROLS) "FULL (previous main specification)" else "MINIMAL (main)", paste(W_cols, collapse = ","), paste(X_cols, collapse = ","))
  
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
  ## centre the continuous Z columns over the assembled alternative rows (as before for the conflict-history and distance
  ## columns; the target changes are treated the same way). A constant shift of a Z column is identical for every
  ## alternative in every choice set and cancels in the softmax, so this changes no likelihood value.
  for (cc in intersect(Z_cols, Z_CENTRED)) alts[[cc]] <- alts[[cc]] - mean(alts[[cc]])
  stopifnot(all(abs(colMeans(as.matrix(alts[, intersect(Z_cols, Z_CENTRED), drop = FALSE]))) < 1e-8))
  ## dyadic lagged DV (PI decision 2026-10-24): 0/1, not centred, must vary within the assembled cell
  stopifnot(!anyNA(alts$dy_mid_tm1), all(alts$dy_mid_tm1 %in% c(0L, 1L)))
  if (length(unique(alts$dy_mid_tm1)) < 2)
    stop(sprintf("dy_mid_tm1 is constant (%d) over the %d assembled alternative rows of this cell (sample %s, lag %d, variant %s); it cannot enter Z",
                 alts$dy_mid_tm1[1], nrow(alts), SAMPLE, LAG, VARIANT))
  
  ## ---- events in this sample -------------------------------------------------------------------------
  ## The sample flag in_<sample> already excludes escalations (events whose pair was fighting in another dispute at t);
  ## asserted here, not assumed. --drop_ww2 (appendix) removes the WWII events and their stage-1 initiations as well.
  ev_all <- events |> filter(.data[[ev_flag]])
  stopifnot(!any(ev_all$escalation))
  n_ww2 <- sum(ev_all$is_ww2)
  if (DROP_WW2) {
    ev_all <- ev_all |> filter(!is_ww2)
    logf("--drop_ww2: %d WWII events (disputes %s) removed from sample %s before stage 1 and stage 2", n_ww2, paste(sort(unique(events$dispnum[events$is_ww2])), collapse = ","), SAMPLE)
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
  ## event) has no alternative row by construction and enters stage 1 only; its (i, j, t) must be absent from alts
  far_keys <- ev_keys |> filter(fighting_at_ref)
  if (nrow(far_keys |> semi_join(alts, by = c("ccode1", "ccode2", "year"))) > 0) stop("fighting_at_ref events have alternative rows; rebuild the data")
  logf("events in sample %s flagged fighting_at_ref (pair already fighting at the reference date; stage 1 only): %d of %d", SAMPLE, nrow(far_keys), nrow(ev_keys))
  ev <- ev_keys |> inner_join(sy_key, by = c("ccode1", "year")) |>
    inner_join(alts |> transmute(ccode1, ccode2, year, alt_row = row_number(), dy_mid_tm1), by = c("ccode1", "ccode2", "year"))
  logf("events in sample %s: %d (from %d dispute pairs); usable (initiator in panel AND target in choice set): %d", SAMPLE, nrow(ev_keys),
       sum(lengths(strsplit(ev_keys$dispnums, ";", fixed = TRUE))), nrow(ev))
  stopifnot(nrow(ev) > 0, all(sy$y1[ev$sy_row] == 1L))
  assert_no_dup_keys(ev, c("sy_row", "alt_row"), "events (one chosen alternative per state-year and target)")
  ## every event whose chosen alternative can exist (initiator in sy, target with positions at t-1 and in the system at t,
  ## pair not fighting at the reference date) has it in its choice set; in the "all" variant without the distance filter
  ## this is exact, elsewhere a lower bound
  expected <- ev_keys |> semi_join(sy_key, by = c("ccode1", "year")) |> filter(target_has_pos_tm1, !fighting_at_ref)
  if (VARIANT == "all" && !DISTANCE) {
    if (nrow(ev) != nrow(expected)) stop(sprintf("%d events with initiator in the panel, a target with positions at t-1 and a pair not fighting at the reference date are not in their choice set (expected %d usable, got %d)", abs(nrow(expected) - nrow(ev)), nrow(expected), nrow(ev)))
  } else stopifnot(nrow(ev) <= nrow(expected))
  stopifnot(nrow(ev_keys |> semi_join(sy_key, by = c("ccode1", "year")) |> filter(!target_has_pos_tm1) |> semi_join(alts, by = c("ccode1", "ccode2", "year"))) == 0)
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
