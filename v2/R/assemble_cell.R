## ============================================================================
## assemble_cell.R -- assemble the sy / alts / ev frames for one real-data cell
## Used by 02_fit_model.R (fits), run_measurement_real.R (measurement tests) so that every
## consumer sees the SAME rows (audit F3/F8: one assembly path, manifest-checked inputs).
## Returns list(sy, alts, ev, W_cols, X_cols, Z_cols, state_codes, x_ref, n_full, sampled, log lines)
## ============================================================================
## Dyadic control columns. Main Z = land contiguity, alliance, joint democracy, major-power count
## (+ conflict history). DISTANCE = TRUE adds log1p(capital distance km) (alternative specification,
## not the main spec); CONTIGUITY = FALSE drops contig_land (robustness option).
z_columns <- function(CONFLICT_HISTORY = TRUE, DISTANCE = FALSE, CONTIGUITY = TRUE) {
  c(if (CONTIGUITY) "contig_land", "alliance_def_off", "joint_dem", "major_count",
    if (CONFLICT_HISTORY) c("log_dy_prior_mids", "log_dy_peace_years"),
    if (DISTANCE) "log_capdist")
}
assemble_cell <- function(SAMPLE, LAG, VARIANT, CENT = "pr", CONFLICT_HISTORY = TRUE, SAMPLED_M = 0L,
                          DISTANCE = FALSE, CONTIGUITY = TRUE, data_dir = file.path(v2_root(), "data")) {
  stopifnot(SAMPLE %in% c("mid", "fmid", "fmidcow", "fmidub", "war"), LAG %in% c(1L, 3L, 5L, 10L), VARIANT %in% c("all", "polrel"), CENT %in% c("pr", "eig"))
  manifest_path <- file.path(data_dir, "manifest.csv")
  stopifnot(file.exists(manifest_path))
  ## ---- load inputs, each asserted against the manifest -------------------------------------
  panel  <- assert_matches_manifest(file.path(data_dir, sprintf("stage1_panel_lag%d.rds", LAG)), manifest_path)
  alts   <- assert_matches_manifest(file.path(data_dir, sprintf("alts_lag%d_%s.rds", LAG, VARIANT)), manifest_path)
  events <- assert_matches_manifest(file.path(data_dir, "events.rds"), manifest_path)
  stopifnot(all(panel$h == LAG), all(alts$h == LAG))
  logf("inputs verified: panel %d rows, alts %d rows, events %d rows", nrow(panel), nrow(alts), nrow(events))
  
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
  ctrl_sy <- c("cowmaj_tm1", "dem_tm1", if (CONFLICT_HISTORY) c("log_prior_mids", "log_peace_years"))
  ## centre the continuous controls so x_ref (below) and the intercept are interpretable
  for (cc in intersect(ctrl_sy, c("log_prior_mids", "log_peace_years"))) sy[[cc]] <- sy[[cc]] - mean(sy[[cc]])
  W_cols <- c("dM_i", "dR_i", ctrl_sy)      # PI decision 2026-10-22: initiator power changes also enter stage 1 directly
  X_cols <- c("dM_i", "dR_i", ctrl_sy)      # power changes + the same controls (documented in README)
  
  ## ---- alternatives -----------------------------------------------------------------------------
  sy_key <- sy |> transmute(ccode1 = ccode, year, sy_row = row_number())
  alts <- alts |> inner_join(sy_key, by = c("ccode1", "year")) |> arrange(sy_row, ccode2) |>
    mutate(state_idx = match(ccode1, state_codes), target_idx = match(ccode2, state_codes), dR = .data[[dRalt_col]],
           log_dy_prior_mids = log1p(dy_prior_mids_h), log_dy_peace_years = log1p(dy_peace_years_h))
  stopifnot(!anyNA(alts$target_idx), !anyNA(alts$sy_row))
  Z_cols <- z_columns(CONFLICT_HISTORY, DISTANCE, CONTIGUITY)
  if (DISTANCE) {
    n_na <- sum(is.na(alts$log_capdist))
    if (n_na > 0) { logf("--distance: dropping %d alternatives without a capital distance", n_na); alts <- alts |> filter(!is.na(log_capdist)) }
  }
  for (cc in intersect(Z_cols, c("log_dy_prior_mids", "log_dy_peace_years", "log_capdist"))) alts[[cc]] <- alts[[cc]] - mean(alts[[cc]])
  
  ## ---- events in this sample -------------------------------------------------------------------------
  ev_all <- events |> filter(.data[[ev_flag]]) |> transmute(ccode1 = initiator, ccode2 = target, year, dispnum)
  ev <- ev_all |> inner_join(sy_key, by = c("ccode1", "year")) |>
    inner_join(alts |> transmute(ccode1, ccode2, year, alt_row = row_number()), by = c("ccode1", "ccode2", "year"))
  logf("events in sample %s: %d; usable (initiator in panel AND target in choice set): %d", SAMPLE, nrow(ev_all), nrow(ev))
  stopifnot(nrow(ev) > 0, all(sy$y1[ev$sy_row] == 1L))
  n_orphan <- sum(sy$y1 == 1L & !(seq_len(nrow(sy)) %in% ev$sy_row))
  logf("state-years with y1 = 1 but no usable stage-2 event (kept in stage 1 only): %d", n_orphan)
  ev <- ev |> select(sy_row, alt_row)
  
  ## ---- optional sampled alternatives ------------------------------------------------------------------
  n_full <- NULL; sampled <- 0L
  if (SAMPLED_M > 0) {
    sa <- sample_alternatives(alts, ev, m = SAMPLED_M, seed = V2_SEED)
    alts <- sa$alts; ev <- sa$ev; n_full <- sa$n_full; sampled <- 1L
    logf("sampled alternatives: m = %d -> %d alternative rows", SAMPLED_M, nrow(alts))
  }
  
  ## ---- reference X for mu_state: mean X over INITIATING state-years in this sample -----------------------
  x_ref <- colMeans(as.matrix(sy[sy$y1 == 1L, X_cols, drop = FALSE]))
  logf("x_ref (mean X over initiating state-years): %s", paste(sprintf("%s=%.3f", X_cols, x_ref), collapse = ", "))
  list(sy = sy, alts = alts, ev = ev, W_cols = W_cols, X_cols = X_cols, Z_cols = Z_cols, state_codes = state_codes,
       x_ref = x_ref, n_full = n_full, sampled = sampled, n_orphan = n_orphan, manifest_path = manifest_path)
}
