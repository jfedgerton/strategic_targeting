## ============================================================================
## assemble_cell.R -- assemble the sy / alts / ev frames for one real-data cell
## Used by 02_fit_model.R (fits) and 05_robustness.R (refit presets, measurement tests) so that every consumer sees the SAME
## rows (audit F3/F8: one assembly path, manifest-checked inputs). assemble_cell_spec() maps a spec (spec.R) onto the arguments.
## Returns list(sy, alts, ev, W_cols, X_cols, Z_cols, state_codes, x_ref, n_full, sampled, log lines)
## ============================================================================
## Switches of spec.R handled here (every other switch is a constant of the main specification since 2026-10-25):
##   PAIRS         = "all" (main) | "originators" (events with is_joiner_pair FALSE) | "joiners" (TRUE). The events of the sample
##                   are filtered and y1 is recomputed from the kept events exactly as DROP_WW2 does, so a state-year whose
##                   only events were in the other group gets y1 = 0. Joiners-only is a STAGE-2 comparison: its stage 1 counts a
##                   state as initiating only when it joined (is_joiner_pair is FALSE for a collapsed key with any originator pair).
##                   D8 weights (PI decision 2026-10-26): in every (initiator, dispute, year) group the filter splits, the kept events'
##                   weights are divided by their kept sum so the group sums to 1 again (EVENT_WEIGHTS = "dispute" only); groups kept
##                   whole keep the build's weights; a collapsed key is grouped by its dispnum column (its smallest dispute).
##   CHANGES      = "raw" (main) | "winsor": the initiator changes dM_i, dR_i and the target changes dM_j_h, dRpr_j_h are clipped
##                   at the WINSOR_PROBS (1st / 99th percentile) quantiles of the stage-1 panel's dM_h / dRpr_h for the cell's lag
##                   (the target change is the same state-year change table, so one pair of cut points per axis), BEFORE the Z
##                   centring; cut points and the number of clipped values are logged. A target without a change value stays 0.
##   CONTEXT       = "all" (main) | "bilateral" | "coalition" (dual typology stage a, 2026-10-25): keeps the STAGE-2 events whose D8 group
##                   (initiator, dispute, onset year) has the given coalition context, read from the build's coal_<COAL_TIES> column
##                   (COAL_TIES = "inclusive" -> coal_incl, "strict" -> coal_strict). Stage 1 is NOT touched: y1 stays the sample's, so a
##                   state-year whose events are all in the other context keeps y1 = 1 without a stage-2 event (counted as an orphan).
##                   A D8 group is kept or dropped whole; the only exception is a key the build collapsed from a bilateral and a coalition
##                   dispute (ORed to coalition), asserted and logged. COAL_DROP = "top3" (CONTEXT = "coalition" only) also removes the
##                   events of the 3 disputes contributing the most coalition events of the sample (a collapsed key counts for each of its
##                   disputes and leaves when any of them is in the top 3; ties at the third place are logged).
##   MAJOR         = "none" (main) | "joint" (PI-approved 2026-10-27, cell suffix _major): Z gains two 0/1 indicators built from the
##                   alternatives' major_count (COW major-power status of both states at t-1, major1 + major2 from the build), one_major =
##                   (major_count == 1) and both_major = (major_count == 2), zero majors the reference; appended after the main Z columns,
##                   NOT centred (the convention for every 0/1 Z column, contig_land, joint_dem, dy_mid_tm1, dch_j_missing; only the
##                   continuous Z columns in Z_CENTRED are centred over the IV rows). They are identified only through within-choice-set
##                   variation (b_j absorbs time-invariant target attractiveness), so the assembly stops unless each varies within some
##                   IV set and within some stage-2 set of an initiating state-year, and logs how many sets vary and the shares of
##                   alternatives and of stage-2 events in each category. W gains the initiator's own indicator cowmaj_tm1 (the panel's
##                   major-power status at t-1, the same source and lag as major1 in major_count), appended after the minimal W controls
##                   and grand-centred like the other 0/1 W controls (dem_tm1, ongoing_tm1; SY_INDICATORS); under FULL_CONTROLS it is
##                   already in W and is not added twice. X is unchanged. Its mean and its correlation with M_lvl are logged.
##   INTERACT      = "none" (main) | "prefixed" (PI-approved 2026-10-27, cell suffix _ix): Z gains four products of the change columns
##                   already in the data, in the scaled change units of X / Z: ix_dMj_dRj = dM_j_h x dRpr_j_h, ix_dMi_dMj = dM_i x dM_j_h,
##                   ix_dRi_dRj = dR_i x dRpr_j_h, ix_dMi_dRj = dM_i x dRpr_j_h, with dM_i / dR_i the initiator changes of the row's
##                   state-year exactly as they enter X and dM_j_h / dRpr_j_h the target changes exactly as they enter Z (both after the
##                   winsor clipping when CHANGES = "winsor"). The products are formed from the UNCENTRED scaled values (before the Z
##                   centring of dM_j_h / dRpr_j_h) and then centred over the IV rows like every continuous Z column (Z_CENTRED); forming
##                   them from the centred target changes instead would be the same model up to a reparameterisation absorbed by the Z
##                   main effects, alpha0 and alpha (dM_i, dR_i are in W). Appended after the major indicators when both switches are on.
##   DROP_WW2, FULL_CONTROLS, EVENT_WEIGHTS: the ww2, controls and weights switches (below).
##   CONTEXT_MODEL = FALSE (every model of MODELS) | TRUE (the M1c model of dual typology stage b, 2026-10-27): the build's state-year
##                   column multi_mid_jan (a multilateral MID the state had not entered is active at 1 January of t) is appended to W as
##                   its last column, grand-centred like the other 0/1 controls. Nothing else changes, so the main specification's W is
##                   untouched. The three stage-2 event groups (bilateral / coalition leader / joiner; stage2_groups() in helpers_data.R,
##                   from coal_incl and is_joiner_pair; COAL_TIES applies to the context filter only) are computed, asserted and logged for EVERY model and returned as ev$grp;
##                   only the M1c model passes them to Stan (02_fit_model.R), so the Stan data of the existing models are unchanged.
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
## Main W since 2026-10-25 (PI decisions D2 and D3): the initiator's power LEVELS at t-1, M_lvl and R_lvl (within-year centred log
## CINC share / log PageRank divided by the gap constants, built by 01_build_data.R; 0 = the system-average level that year, so
## they are NOT grand-centred here), and the centred log system size at t-1, log_n_system_tm1 (the IV grows with the number of
## alternatives, so tau * IV carried a system-size trend). Neither enters X. FULL_CONTROLS stays exactly 8bf0ce1.
W_LEVELS <- c("M_lvl", "R_lvl")
## MAJOR = "joint" (2026-10-27) appends the initiator's major-power indicator cowmaj_tm1 to the minimal W; the full set of 8bf0ce1 holds it already
w_controls <- function(CONFLICT_HISTORY = TRUE, FULL_CONTROLS = FALSE, MAJOR = "none") {
  stopifnot(MAJOR %in% c("none", "joint"))
  if (FULL_CONTROLS) c("cowmaj_tm1", "dem_tm1", if (CONFLICT_HISTORY) c("log_prior_mids", "log_peace_years"))
  else c("M_lvl", "R_lvl", "ongoing_tm1", if (CONFLICT_HISTORY) "log_prior_mids", "dem_tm1", "log_n_system_tm1", if (MAJOR == "joint") "cowmaj_tm1")
}
x_controls <- function(CONFLICT_HISTORY = TRUE, FULL_CONTROLS = FALSE) {
  if (FULL_CONTROLS) c("cowmaj_tm1", "dem_tm1", if (CONFLICT_HISTORY) c("log_prior_mids", "log_peace_years"))
  else character(0)
}
## the Z columns the 2026-10-27 switches append (after the base list, major indicators first): Z_MAJOR under MAJOR = "joint" (0/1, not centred),
## Z_INTERACT under INTERACT = "prefixed" (continuous, centred over the IV rows after the products are formed)
Z_MAJOR    <- c("one_major", "both_major")
Z_INTERACT <- c("ix_dMj_dRj", "ix_dMi_dMj", "ix_dRi_dRj", "ix_dMi_dRj")
z_columns <- function(CONFLICT_HISTORY = TRUE, DISTANCE = FALSE, CONTIGUITY = TRUE, TARGET_MISSING = FALSE, DY_PEACE_YEARS = FALSE,
                      FULL_CONTROLS = FALSE, MAJOR = "none", INTERACT = "none") {
  stopifnot(MAJOR %in% c("none", "joint"), INTERACT %in% c("none", "prefixed"))
  base <- if (FULL_CONTROLS) {
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
  c(base, if (MAJOR == "joint") Z_MAJOR, if (INTERACT == "prefixed") Z_INTERACT)
}
## continuous Z columns that assemble_cell() centres over the assembled alternative rows (the 0/1 indicators, dy_mid_tm1
## included, stay as is); the interaction products (2026-10-27) are continuous and centred like the target changes they are built from
Z_CENTRED <- c("log_dy_prior_mids", "log_dy_peace_years", "log_capdist", "dM_j_h", "dRpr_j_h", Z_INTERACT)
## Monadic 0/1 indicators available in the panel (major power, democracy, ongoing_tm1 = participant in any MID active at
## t-1). Those that enter W or X are centred like the conflict-history counts; the main W uses dem_tm1 and ongoing_tm1.
SY_INDICATORS <- c("cowmaj_tm1", "dem_tm1", "ongoing_tm1")
## manifest: NULL reads data_dir/manifest.csv now; a data.frame is a SNAPSHOT read earlier by the caller (02_fit_model.R passes
## the manifest it read at load, so a rebuild between that read and the loading of the inputs stops the run; review P2-5).
## EVENT_WEIGHTS (PI decision 2026-10-25, D8): "dispute" (main) gives every stage-2 event the build's weight w_<sample> = 1 / (the
## initiator's targets in the same dispute and onset year), a composite likelihood with dispute-level weighting; "none" gives weight 1
## (the robustness comparison, cell suffix _w1). The weights multiply the stage-2 log-lik terms in the Stan programs; log_lik_stage2
## stays unweighted per event.
assemble_cell <- function(SAMPLE, LAG, VARIANT, CENT = "pr", CONFLICT_HISTORY = TRUE, SAMPLED_M = 0L,
                          DISTANCE = FALSE, CONTIGUITY = TRUE, DY_PEACE_YEARS = FALSE, DROP_WW2 = FALSE, FULL_CONTROLS = FALSE,
                          EVENT_WEIGHTS = "dispute", PAIRS = "all", CHANGES = "raw", CONTEXT = "all", COAL_TIES = "inclusive", COAL_DROP = "none",
                          MAJOR = "none", INTERACT = "none",
                          CONTEXT_MODEL = FALSE, data_dir = file.path(v2_root(), "data"), manifest = NULL) {
  stopifnot(SAMPLE %in% c("mid", "fmid", "fmidcow", "fmidub", "war"), LAG %in% c(1L, 3L, 5L, 10L), VARIANT %in% c("all", "polrel"), CENT %in% c("pr", "eig"),
            EVENT_WEIGHTS %in% c("dispute", "none"), PAIRS %in% c("all", "originators", "joiners"), CHANGES %in% c("raw", "winsor"),
            CONTEXT %in% c("all", "bilateral", "coalition"), COAL_TIES %in% c("inclusive", "strict"), COAL_DROP %in% c("none", "top3"),
            MAJOR %in% c("none", "joint"), INTERACT %in% c("none", "prefixed"))
  if (COAL_DROP != "none" && CONTEXT != "coalition") stop("COAL_DROP = top3 needs CONTEXT = coalition")
  stopifnot(is.logical(CONTEXT_MODEL), length(CONTEXT_MODEL) == 1L, !is.na(CONTEXT_MODEL))
  if (CONTEXT_MODEL && FULL_CONTROLS) stop("CONTEXT_MODEL (M1c) is defined on the minimal control set only (FULL_CONTROLS = TRUE not supported)")
  ## merge 2026-10-27: the M1c model and the major / interact switches were developed independently; their combination is not specified
  ## (W would carry cowmaj_tm1 before multi_mid_jan, Z the extra columns), so it is refused here until the PI decides
  if (CONTEXT_MODEL && (MAJOR != "none" || INTERACT != "none")) stop("CONTEXT_MODEL (M1c) with MAJOR != none or INTERACT != none is not specified; run M1c with MAJOR = none, INTERACT = none")
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
            all(c("ongoing_tm1", "M_lvl", "R_lvl", "log_n_system_tm1") %in% names(panel)), all(c("dy_mid_tm1", "in_iv", "in_s2") %in% names(alts)),
            all(c("dir_flipped", "dir_source", paste0("w_", SAMPLE), "coal_incl", "coal_strict") %in% names(events)))
  stopifnot(!anyNA(events$coal_incl[!events$escalation]), !anyNA(events$coal_strict[!events$escalation]), all(events$coal_incl[!events$escalation] >= events$coal_strict[!events$escalation]))
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
  no_cs <- panel |> anti_join(alts, by = c("ccode" = "ccode1", "year"))
  logf("panel state-years with a non-empty choice set: %d of %d", nrow(has_cs), nrow(panel))
  if (nrow(no_cs) > 0) {
    ## P3 (2026-10-25): in the all-states variant every panel state-year has alternatives by construction, so a missing block is an
    ## error; in polrel an isolated non-major state can have none and is dropped, listed with its ccodes
    msg <- sprintf("%d panel state-years have NO alternative row in variant %s: %s", nrow(no_cs), VARIANT, paste(sprintf("%d-%d", no_cs$ccode, no_cs$year), collapse = " "))
    if (VARIANT == "all") stop(msg) else logf("dropped from stage 1: %s", msg)
  }
  n_iv_sy <- alts |> group_by(ccode1, year) |> summarise(n_iv = sum(in_iv), .groups = "drop")
  if (any(n_iv_sy$n_iv == 0L)) stop(sprintf("%d state-years have alternative rows but an EMPTY stage-1 IV set (no row with in_iv = 1); rebuild the data", sum(n_iv_sy$n_iv == 0L)))
  sy <- has_cs |> arrange(year, ccode) |>
    mutate(y1 = .data[[y1_col]],
           log_prior_mids = log1p(mon_prior_mids_h), log_peace_years = log1p(mon_peace_years_h),
           dM_i = dM_h, dR_i = .data[[dR_col]])
  ## changes = winsor: cut points from the stage-1 PANEL of this lag (every state-year the build wrote, before the choice-set join),
  ## one pair per axis, applied to the initiator changes here and to the target changes below (same change table)
  winsor <- NULL
  if (CHANGES == "winsor") {
    winsor <- list(M = unname(quantile(panel$dM_h, WINSOR_PROBS, type = 7)), R = unname(quantile(panel[[dR_col]], WINSOR_PROBS, type = 7)))
    stopifnot(all(is.finite(unlist(winsor))), winsor$M[1] < 0, winsor$M[2] > 0, winsor$R[1] < 0, winsor$R[2] > 0)   # 0 stays inside the caps
    n_clip_M <- sum(sy$dM_i < winsor$M[1] | sy$dM_i > winsor$M[2]); n_clip_R <- sum(sy$dR_i < winsor$R[1] | sy$dR_i > winsor$R[2])
    sy <- sy |> mutate(dM_i = pmin(pmax(dM_i, winsor$M[1]), winsor$M[2]), dR_i = pmin(pmax(dR_i, winsor$R[1]), winsor$R[2]))
    logf("changes = winsor: cut points at the %s / %s percentiles of the lag-%d panel (%d state-years): dM_h [%.4f, %.4f], %s [%.4f, %.4f]; initiator values clipped: dM_i %d, dR_i %d of %d",
         WINSOR_PROBS[1] * 100, WINSOR_PROBS[2] * 100, LAG, nrow(panel), winsor$M[1], winsor$M[2], dR_col, winsor$R[1], winsor$R[2], n_clip_M, n_clip_R, nrow(sy))
    stopifnot(all(sy$dM_i >= winsor$M[1] & sy$dM_i <= winsor$M[2]), all(sy$dR_i >= winsor$R[1] & sy$dR_i <= winsor$R[2]))
  }
  ## state index: union of initiators (panel) and targets (alts), sorted by ccode (index only; anchors are NOT taken from this order)
  state_codes <- sort(unique(c(sy$ccode, alts$ccode2)))
  sy$state_idx <- match(sy$ccode, state_codes)
  ctrl_w <- w_controls(CONFLICT_HISTORY, FULL_CONTROLS, MAJOR = MAJOR); ctrl_x <- x_controls(CONFLICT_HISTORY, FULL_CONTROLS)
  stopifnot(!anyDuplicated(ctrl_w))
  if (CONTEXT_MODEL) {
    ## stage b: the outcome-blind multilateral-opportunity control, LAST in W (0/1 from the build; grand-centred below like the other indicators)
    if (!"multi_mid_jan" %in% names(sy)) stop("CONTEXT_MODEL (M1c) needs the panel column multi_mid_jan (rebuild the data with 01_build_data.R of 2026-10-27 or later)")
    stopifnot(!anyNA(sy$multi_mid_jan), all(sy$multi_mid_jan %in% c(0, 1)))
    logf("multi_mid_jan (stage b, M1c): mean %.4f over %d state-years (%d with a multilateral MID they had not entered active at 1 January of t); among initiators (y1 = 1) %.4f", mean(sy$multi_mid_jan), nrow(sy), sum(sy$multi_mid_jan), mean(sy$multi_mid_jan[sy$y1 == 1L]))
    ctrl_w <- c(ctrl_w, "multi_mid_jan")
    stopifnot(!anyDuplicated(ctrl_w))
  }
  ctrl_sy <- union(ctrl_w, ctrl_x)
  ## major = joint (2026-10-27): the initiator's own major-power indicator at t-1 enters W (grand-centred below like dem_tm1 and ongoing_tm1);
  ## logged BEFORE the centring with its correlation with the material level M_lvl, which it is bound to track
  if (MAJOR == "joint") {
    stopifnot("cowmaj_tm1" %in% ctrl_w, "cowmaj_tm1" %in% names(sy), all(sy$cowmaj_tm1 %in% c(0L, 1L)))
    if (length(unique(sy$cowmaj_tm1)) < 2) stop("major = joint: cowmaj_tm1 is constant over the stage-1 panel of this cell; it cannot enter W")
    logf("major = joint: cowmaj_tm1 in W (%s): mean %.4f over %d state-years (%d major-power state-years); cor(cowmaj_tm1, M_lvl) = %.3f; grand-centred like the other 0/1 W controls",
         if (FULL_CONTROLS) "already in the full control set of 8bf0ce1, not added twice" else "appended to the minimal controls",
         mean(sy$cowmaj_tm1), nrow(sy), sum(sy$cowmaj_tm1), cor(sy$cowmaj_tm1, sy$M_lvl))
  }
  ## centre every stage-1 / mu control that enters W or X (the 0/1 indicators as well as the conflict-history counts; PI
  ## decisions 2026-10-23 / 2026-10-24) so that alpha0 in stage 1 is the log-odds at zero power change and average
  ## controls, and (under --full_controls) the mu intercept mu0 is the ideal gap at those averages. In the main
  ## specification X has no controls, so mu0 is the ideal gap of a state with zero power change. dM_i and dR_i are NOT
  ## centred (zero power change is the reference).
  stopifnot(all(ctrl_sy %in% names(sy)), all(ctrl_x %in% ctrl_w), setequal(ctrl_x, if (FULL_CONTROLS) ctrl_w else character(0)))
  for (cc in SY_INDICATORS) stopifnot(!anyNA(sy[[cc]]), all(sy[[cc]] %in% c(0, 1)))
  logf("ongoing_tm1: mean %.4f over %d state-years (%d in a MID active at t-1)", mean(sy$ongoing_tm1), nrow(sy), sum(sy$ongoing_tm1))
  ## the power levels (D2) are already centred within year over every state with positions (0 = the system average that year) and
  ## are left as they are, like dM_i and dR_i; every other control is grand-centred
  ctrl_cent <- setdiff(ctrl_sy, W_LEVELS)
  for (cc in ctrl_cent) sy[[cc]] <- sy[[cc]] - mean(sy[[cc]])
  stopifnot(all(abs(colMeans(as.matrix(sy[, ctrl_cent, drop = FALSE]))) < 1e-8), all(is.finite(as.matrix(sy[, ctrl_sy, drop = FALSE]))))
  if (any(W_LEVELS %in% ctrl_sy)) logf("power levels in W (D2, within-year centred over all states with positions, gap-scaled, NOT grand-centred): M_lvl mean %.4f SD %.3f, R_lvl mean %.4f SD %.3f",
                                        mean(sy$M_lvl), sd(sy$M_lvl), mean(sy$R_lvl), sd(sy$R_lvl))
  W_cols <- c("dM_i", "dR_i", ctrl_w)      # PI decision 2026-10-22: initiator power changes also enter stage 1 directly
  X_cols <- c("dM_i", "dR_i", ctrl_x)      # main: power changes only; --full_controls: the same controls as W
  logf("control set: %s%s; W = %s; X = %s", if (FULL_CONTROLS) "FULL (base specification of commit 8bf0ce1)" else "MINIMAL (main)", if (CONTEXT_MODEL) " + multi_mid_jan (M1c)" else "", paste(W_cols, collapse = ","), paste(X_cols, collapse = ","))
  
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
  if (CHANGES == "winsor") {
    ## the target change is the same state-year change as the initiator's (asserted by the build), so the SAME cut points apply;
    ## the relational target change follows the PageRank table (dRpr_j_h) in both centralities, as in z_columns()
    n_clip_Mj <- sum(alts$dM_j_h < winsor$M[1] | alts$dM_j_h > winsor$M[2]); n_clip_Rj <- sum(alts$dRpr_j_h < winsor$R[1] | alts$dRpr_j_h > winsor$R[2])
    alts <- alts |> mutate(dM_j_h = pmin(pmax(dM_j_h, winsor$M[1]), winsor$M[2]), dRpr_j_h = pmin(pmax(dRpr_j_h, winsor$R[1]), winsor$R[2]))
    stopifnot(all(alts$dM_j_h[alts$dch_j_missing == 1L] == 0), all(alts$dRpr_j_h[alts$dch_j_missing == 1L] == 0))
    logf("changes = winsor: target values clipped before the Z centring: dM_j_h %d, dRpr_j_h %d of %d alternative rows", n_clip_Mj, n_clip_Rj, nrow(alts))
  }
  logf("target changes: %d of %d alternative rows have a target without a change value (set to 0); dch_j_missing %s Z (%s)",
       n_miss, nrow(alts), if (target_missing) "enters" else "is NOT in",
       if (target_missing) "non-constant in this cell" else "constant in this cell")
  ## major = joint (2026-10-27): the dyad's number of major powers at t-1 (major_count = major1 + major2 from the build, COW majors lagged like
  ## every dyadic control) as two 0/1 indicators with zero majors as the reference. Shares over the alternative rows are logged here; the
  ## within-choice-set variation (the only source of identification next to b_j) and the event shares are checked below, once ev exists.
  major_stats <- NULL
  if (MAJOR == "joint") {
    stopifnot("major_count" %in% names(alts), !anyNA(alts$major_count), all(alts$major_count %in% 0:2))
    alts <- alts |> mutate(one_major = as.integer(major_count == 1L), both_major = as.integer(major_count == 2L))
    stopifnot(all(alts$one_major + alts$both_major <= 1L), all((alts$one_major == 1L) == (alts$major_count == 1L)), all((alts$both_major == 1L) == (alts$major_count == 2L)))
    sh_alt <- tabulate(alts$major_count + 1L, nbins = 3L) / nrow(alts)
    logf("major = joint: alternative rows by the dyad's major-power count at t-1: 0 majors %.4f, 1 major (one_major) %.4f, 2 majors (both_major) %.4f of %d rows; the indicators enter Z uncentred (0/1 convention)",
         sh_alt[1], sh_alt[2], sh_alt[3], nrow(alts))
    major_stats <- list(share_alts = setNames(sh_alt, c("zero", "one", "both")))
  }
  ## interact = prefixed (2026-10-27): four products of the change columns, formed from the UNCENTRED scaled values exactly as they stand after
  ## the winsor step (dM_i / dR_i of the row's state-year as in X, dM_j_h / dRpr_j_h of the target as in Z before its centring), then centred
  ## over the IV rows with the other continuous Z columns below (Z_INTERACT is in Z_CENTRED)
  if (INTERACT == "prefixed") {
    stopifnot(!anyNA(alts$sy_row), all(alts$sy_row %in% seq_len(nrow(sy))))
    dMi <- sy$dM_i[alts$sy_row]; dRi <- sy$dR_i[alts$sy_row]
    stopifnot(all(is.finite(dMi)), all(is.finite(dRi)), all(is.finite(alts$dM_j_h)), all(is.finite(alts$dRpr_j_h)))
    alts <- alts |> mutate(ix_dMj_dRj = dM_j_h * dRpr_j_h, ix_dMi_dMj = .env$dMi * dM_j_h, ix_dRi_dRj = .env$dRi * dRpr_j_h, ix_dMi_dRj = .env$dMi * dRpr_j_h)
    stopifnot(all(is.finite(as.matrix(alts[, Z_INTERACT]))))
    logf("interact = prefixed: four change products in Z, formed before the Z centring (%s): SD before centring ix_dMj_dRj %.4f, ix_dMi_dMj %.4f, ix_dRi_dRj %.4f, ix_dMi_dRj %.4f over %d rows; a target without a change value contributes 0 to every product",
         if (CHANGES == "winsor") "winsorised changes" else "raw changes", sd(alts$ix_dMj_dRj), sd(alts$ix_dMi_dMj), sd(alts$ix_dRi_dRj), sd(alts$ix_dMi_dRj), nrow(alts))
  }
  Z_cols <- z_columns(CONFLICT_HISTORY, DISTANCE, CONTIGUITY, TARGET_MISSING = target_missing, DY_PEACE_YEARS = DY_PEACE_YEARS, FULL_CONTROLS = FULL_CONTROLS,
                      MAJOR = MAJOR, INTERACT = INTERACT)
  stopifnot(all(Z_cols %in% names(alts)), !anyDuplicated(Z_cols))
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
  ## ---- D3 diagnostic: the IV at a null point against the system size -----------------------------------------------------------
  ## IV_null[n] = log sum over the IV rows of exp(-1/2 (dM^2 + dR^2)) (mu = 0, Sigma = I, b = gamma = 0, no u'v): it grows with the
  ## number of alternatives, so without log_n_system_tm1 in W the trend in the system's size entered stage 1 through tau * IV.
  ## Logged: cor(IV_null, log_n_system_tm1), and the share of IV_null's variance explained by W BEFORE (W without log_n_system_tm1)
  ## and AFTER (with it) adding the column.
  iv_null <- as.numeric(rowsum(exp(-0.5 * (alts$dM^2 + alts$dR^2)) * (alts$in_iv == 1L), alts$sy_row, reorder = TRUE))
  stopifnot(length(iv_null) == nrow(sy), all(iv_null > 0)); iv_null <- log(iv_null)
  r2 <- function(cols) if (length(cols)) summary(lm(iv_null ~ as.matrix(sy[, cols, drop = FALSE])))$r.squared else 0
  w_before <- setdiff(c("dM_i", "dR_i", ctrl_w), "log_n_system_tm1"); w_after <- c("dM_i", "dR_i", ctrl_w)
  logf("IV at a null point (D3): cor(IV_null, log_n_system_tm1) = %.3f; R^2 of IV_null on W without / with log_n_system_tm1 = %.3f / %.3f (log_n_system_tm1 %s W)",
       cor(iv_null, sy$log_n_system_tm1), r2(w_before), r2(w_after), if ("log_n_system_tm1" %in% ctrl_w) "in" else "NOT in")

  ## ---- events in this sample -------------------------------------------------------------------------
  ## The sample flag in_<sample> already excludes escalations (events whose pair was fighting in another dispute at t);
  ## asserted here, not assumed. --drop_ww2 (appendix) uses the build's in_<sample>_noww2 flag (review P2-3): the sample flag
  ## of the non-WWII disputes of the (i, j, t) alone, so an event collapsed from a WWII and a non-WWII dispute stays when the
  ## non-WWII dispute qualifies (the old rule in_<sample> & !is_ww2 dropped it). Their stage-1 initiations go with them.
  ev_all <- events |> filter(.data[[ev_flag]])
  stopifnot(!any(ev_all$escalation))
  n_ww2 <- sum(ev_all$is_ww2)
  ## D8 weight column of the sample; every in-sample event carries a weight in (0, 1] from the build (the pairs rule below may rescale it)
  w_col <- paste0("w_", SAMPLE)
  stopifnot(!anyNA(ev_all[[w_col]]), all(ev_all[[w_col]] > 0), all(ev_all[[w_col]] <= 1))
  if (DROP_WW2) {
    ev_flag_noww2 <- paste0(ev_flag, "_noww2")
    n_before <- nrow(ev_all); ev_before <- ev_all
    ev_all <- events |> filter(.data[[ev_flag_noww2]])
    stopifnot(!any(ev_all$escalation), all(events[[ev_flag]][events[[ev_flag_noww2]]]), !any(ev_all$is_ww2 & ev_all$n_disputes == 1L))
    logf("--drop_ww2: sample %s %d -> %d events using %s (%d events flagged is_ww2 removed; %d events collapsed from a WWII and a non-WWII dispute kept because a non-WWII dispute qualifies, whereas in_%s & !is_ww2 would keep %d events); WWII disputes %s",
         SAMPLE, n_before, nrow(ev_all), ev_flag_noww2, n_before - nrow(ev_all), sum(ev_all$is_ww2), SAMPLE, sum(events[[ev_flag]] & !events$is_ww2),
         paste(sort(unique(events$dispnum[events$is_ww2])), collapse = ","))
    ## D8 group integrity (2026-10-26): whole disputes leave, so every (initiator, dispute, year) group is untouched or dropped whole; a split
    ## among single-dispute members stops the run; a split through a collapsed key kept by its non-WWII dispute is logged, NOT rescaled
    ww <- d8_assert_groups_whole(ev_before, ev_all, w_col, "--ww2=drop")
    logf("--ww2=drop D8 groups (initiator, dispute, year) of sample %s: %d kept whole, %d dropped whole, %d split through a collapsed key kept by its non-WWII dispute%s; weights %.1f -> %.1f, not rescaled",
         SAMPLE, ww$n_whole_kept, ww$n_whole_dropped, ww$n_split_collapsed, if (ww$n_split_collapsed > 0) " (WARNING: these groups no longer sum to one decision)" else "",
         sum(ev_before[[w_col]]), sum(ev_all[[w_col]]))
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
  ## pairs = originators | joiners: keep one group of event pairs and recompute y1 from the kept events (as --ww2=drop does), so the
  ## stage-1 outcome and the stage-2 events describe the same group. is_joiner_pair is FALSE for a collapsed key with any originator
  ## pair. Joiners-only is a stage-2 comparison (its stage 1 counts a state as initiating only when it joined).
  ## D8 weights under pairs != all (PI decision 2026-10-26): the filter can keep part of an (initiator, dispute, year) group (a dispute
  ## where the initiator faced originators AND joiners), so the kept members' weights are divided by their kept sum and every split group
  ## sums to 1 again; groups kept whole keep the build's weights exactly. A collapsed key is grouped by its dispnum column (the smallest
  ## of its disputes), the only group it can belong to after the collapse (d8_rescale_split_groups, helpers_data.R). Only this filter
  ## rescales: ww2 = drop and coal_drop = top3 remove whole disputes (asserted); context keeps groups whole (asserted above / below); the
  ## events that leave stage 2 through the choice set are logged, not rescaled. With weights = none every event has weight 1 anyway.
  pairs_w <- list(n_groups = 0L, n_groups_split = 0L, n_groups_rescaled = 0L, n_events_rescaled = 0L, w_before = NA_real_, w_after = NA_real_)
  if (PAIRS != "all") {
    n_before <- nrow(ev_all); n_join_before <- sum(ev_all$is_joiner_pair); y1_before <- sum(sy$y1); ev_before <- ev_all
    ev_all <- if (PAIRS == "originators") ev_all |> filter(!is_joiner_pair) else ev_all |> filter(is_joiner_pair)
    if (nrow(ev_all) == 0) stop(sprintf("pairs = %s leaves no event in sample %s", PAIRS, SAMPLE))
    y1_pairs <- ev_all |> distinct(ccode = initiator, year) |> mutate(.y1 = 1L)
    sy <- sy |> left_join(y1_pairs, by = c("ccode", "year")) |> mutate(.y1 = coalesce(.y1, 0L))
    stopifnot(all(sy$.y1 <= sy$y1), all(sy$.y1 %in% 0:1))
    sy$y1 <- sy$.y1; sy <- sy |> select(-.y1)
    logf("pairs = %s: sample %s %d -> %d events (%d originator pairs and %d joiner pairs before the filter); stage-1 initiations %d -> %d (%d state-years whose %s events were all in the other group now have y1 = 0)",
         PAIRS, SAMPLE, n_before, nrow(ev_all), n_before - n_join_before, n_join_before, y1_before, sum(sy$y1), y1_before - sum(sy$y1), SAMPLE)
    if (EVENT_WEIGHTS == "dispute") {
      ev_all <- d8_rescale_split_groups(ev_before, ev_all, w_col)
      g <- attr(ev_all, "groups")
      pairs_w <- list(n_groups = nrow(g), n_groups_split = attr(ev_all, "n_groups_split"), n_groups_rescaled = attr(ev_all, "n_groups_rescaled"),
                      n_events_rescaled = attr(ev_all, "n_events_rescaled"), w_before = attr(ev_all, "w_before"), w_after = attr(ev_all, "w_after"))
      stopifnot(abs(pairs_w$w_before - sum(ev_before[[w_col]][ev_before$is_joiner_pair == (PAIRS == "joiners")])) < 1e-9, nrow(ev_all) == sum(ev_before$is_joiner_pair == (PAIRS == "joiners")))
      logf("pairs = %s, D8 weights rescaled within the (initiator, dispute, year) groups the filter split (PI decision 2026-10-26; a collapsed key belongs to the group of its smallest dispute number): %d groups in sample %s, %d kept whole (build weights unchanged; %d of them do not sum to 1 through collapsed keys, as under pairs = all), %d dropped whole, %d split (%d holding a collapsed key); %d split groups rescaled to sum 1 (%d had a kept sum of 1 already), %d of %d kept events rescaled; sum of weights over the kept events %.1f -> %.1f",
           PAIRS, pairs_w$n_groups, SAMPLE, sum(g$n_kept == g$n), sum(g$n_kept == g$n & abs(g$w_before - 1) > 1e-9), sum(g$n_kept == 0L),
           pairs_w$n_groups_split, attr(ev_all, "n_split_with_collapsed"), pairs_w$n_groups_rescaled, pairs_w$n_groups_split - pairs_w$n_groups_rescaled,
           pairs_w$n_events_rescaled, nrow(ev_all), pairs_w$w_before, pairs_w$w_after)
      for (a in c("n_groups_split", "n_groups_rescaled", "n_events_rescaled", "n_events_in_split", "w_before", "w_after", "n_split_with_collapsed", "groups"))
        attr(ev_all, a) <- NULL                                       # the helper's attributes do not travel on
    } else logf("pairs = %s with weights = none: weight 1 on every event, no D8 rescale", PAIRS)
  } else logf("pairs = all: %d events in sample %s (%d joiner pairs); D8 weights as built", nrow(ev_all), SAMPLE, sum(ev_all$is_joiner_pair))
  ## context = bilateral | coalition (dual typology stage a, 2026-10-25): STAGE-2 events of one coalition context only. y1 is NOT recomputed
  ## (stage 1 is unchanged by design; a state-year whose events are all in the other context keeps y1 = 1 and is counted as an orphan below).
  coal_col <- c(inclusive = "coal_incl", strict = "coal_strict")[[COAL_TIES]]
  context_top3 <- NULL; n_context_dropped <- 0L; n_split_groups <- 0L
  if (CONTEXT != "all") {
    n_before <- nrow(ev_all); w_before <- sum(ev_all[[w_col]])
    keep <- if (CONTEXT == "coalition") ev_all[[coal_col]] else !ev_all[[coal_col]]
    stopifnot(!anyNA(keep))
    dropped <- ev_all[!keep, ]; ev_all <- ev_all[keep, ]
    if (nrow(ev_all) == 0) stop(sprintf("context = %s (%s) leaves no event in sample %s", CONTEXT, COAL_TIES, SAMPLE))
    n_context_dropped <- nrow(dropped)
    ## D8 group integrity: every (initiator, dispute, year) group of the sample is kept or dropped WHOLE, except through a key the build
    ## collapsed from a bilateral and a coalition dispute (ORed to coalition). Per group: the single-dispute members must agree, and any
    ## member with the other status must be a collapsed key. Anything else is an error.
    gs <- bind_rows(ev_all |> mutate(.kept = TRUE), dropped |> mutate(.kept = FALSE)) |> select(initiator, year, dispnums, n_disputes, .kept) |>
      mutate(.disp = strsplit(dispnums, ";", fixed = TRUE)) |> tidyr::unnest(.disp) |> mutate(dispnum = as.integer(.disp)) |>
      group_by(initiator, dispnum, year) |>
      summarise(n = n(), n_kept = sum(.kept), n_single_kept = sum(.kept & n_disputes == 1L), n_single = sum(n_disputes == 1L),
                n_coll_kept = sum(.kept & n_disputes > 1L), n_coll = sum(n_disputes > 1L), .groups = "drop")
    split_g <- gs |> filter(n_kept > 0L, n_kept < n)
    bad_single <- split_g |> filter(n_single_kept > 0L, n_single_kept < n_single)
    if (nrow(bad_single) > 0) stop(sprintf("context = %s: %d (initiator, dispute, year) groups are split between kept and dropped events among their single-dispute members; the coalition flag is not constant within the group (rebuild the data)", CONTEXT, nrow(bad_single)))
    ## in a split group the minority members are collapsed keys (given the check above, by construction), and the split can go ONE way
    ## only, because the collapse ORs the flag over a key's disputes: a collapsed key in a COALITION group is itself coalition, so under
    ## context = bilateral the single-dispute members of a split group must all be kept (bilateral group, collapsed keys dropped through
    ## their other dispute) and under context = coalition they must all be dropped (bilateral group, collapsed keys kept through the other
    ## dispute). The reverse split means a collapsed key's flag contradicts its own group's flag (audit fix 21, 2026-10-26: the former
    ## bad_mix condition was unreachable after the single-member check; a group of collapsed keys only cannot be checked this way).
    bad_dir <- split_g |> filter(n_single > 0L, if (CONTEXT == "bilateral") n_single_kept < n_single else n_single_kept > 0L)
    if (nrow(bad_dir) > 0) stop(sprintf("context = %s: %d split groups have a collapsed key whose flag contradicts the group's own flag (a key collapsed from a %s dispute must carry the coalition flag of the OR); rebuild the data",
                                        CONTEXT, nrow(bad_dir), "coalition"))
    stopifnot(all(split_g$n_kept == split_g$n_single_kept + split_g$n_coll_kept))
    n_split_groups <- nrow(split_g)
    logf("context = %s (ties %s, column %s): sample %s %d -> %d stage-2 events (%d removed; weights %.1f -> %.1f); stage 1 unchanged (y1 stays %d initiations); (initiator, dispute, year) groups of the sample %d, kept whole %d, dropped whole %d, split through collapsed keys mixing a bilateral and a coalition dispute %d (their weights no longer sum to one decision; logged, not rescaled)",
         CONTEXT, COAL_TIES, coal_col, SAMPLE, n_before, nrow(ev_all), n_before - nrow(ev_all), w_before, sum(ev_all[[w_col]]), sum(sy$y1),
         nrow(gs), sum(gs$n_kept == gs$n), sum(gs$n_kept == 0L), n_split_groups)
    if (COAL_DROP == "top3") {
      ## the 3 disputes contributing the most coalition events of this sample (a collapsed key counts for each of its disputes)
      cnt <- ev_all |> select(dispnums) |> mutate(.disp = strsplit(dispnums, ";", fixed = TRUE)) |> tidyr::unnest(.disp) |> count(dispnum = as.integer(.disp), name = "n_events") |> arrange(desc(n_events), dispnum)
      if (nrow(cnt) < 4L) stop(sprintf("coal_drop = top3: only %d disputes contribute coalition events in sample %s", nrow(cnt), SAMPLE))
      context_top3 <- cnt[1:3, ]
      tie <- cnt$n_events[4] == cnt$n_events[3]
      in_top3 <- vapply(strsplit(ev_all$dispnums, ";", fixed = TRUE), function(z) any(as.integer(z) %in% context_top3$dispnum), logical(1))
      n_b <- nrow(ev_all); w_b <- sum(ev_all[[w_col]]); ev_b <- ev_all
      ev_all <- ev_all[!in_top3, ]
      if (nrow(ev_all) == 0) stop("coal_drop = top3 leaves no event")
      n_context_dropped <- n_context_dropped + sum(in_top3)
      logf("coal_drop = top3: the 3 disputes with the most coalition events in sample %s are %s; %d -> %d events (%d removed, weights %.1f -> %.1f)%s",
           SAMPLE, paste(sprintf("%d (%d events)", context_top3$dispnum, context_top3$n_events), collapse = ", "), n_b, nrow(ev_all), sum(in_top3), w_b, sum(ev_all[[w_col]]),
           if (tie) sprintf("; WARNING: dispute %d ties with the third place (%d events) and is kept (smallest dispute numbers win ties)", cnt$dispnum[4], cnt$n_events[4]) else "")
      ## D8 group integrity (2026-10-26): whole disputes leave, so groups are untouched or dropped whole; the only possible split is a collapsed
      ## key leaving through its OTHER (top-3) dispute while its own group stays (none on the current data); logged, NOT rescaled
      cd <- d8_assert_groups_whole(ev_b, ev_all, w_col, "coal_drop = top3")
      logf("coal_drop = top3 D8 groups (initiator, dispute, year) of the coalition events: %d kept whole, %d dropped whole, %d split through a collapsed key%s; weights not rescaled",
           cd$n_whole_kept, cd$n_whole_dropped, cd$n_split_collapsed, if (cd$n_split_collapsed > 0) " (WARNING: these groups no longer sum to one decision)" else "")
    }
  } else logf("context = all: %d events in sample %s (coalition context: %d inclusive, %d strict)", nrow(ev_all), SAMPLE, sum(ev_all$coal_incl), sum(ev_all$coal_strict))
  ## one chosen alternative per (i, j, t): same-year duplicates were collapsed by the build (PI decision 2026-10-24)
  ## stage-2 event groups (dual typology stage b, 2026-10-27): bilateral / coalition leader / joiner from coal_incl (the spec's column; the
  ## stage-a tie convention COAL_TIES applies to the context FILTER only) and is_joiner_pair (stage2_groups(): exactly one group per event,
  ## asserted); counted with their D8 weight sums here and again on the usable events below
  grp_all <- stage2_groups(ev_all, "coal_incl")
  ev_keys <- ev_all |> transmute(ccode1 = initiator, ccode2 = target, year, dispnum, dispnums, target_has_pos_tm1, fighting_at_ref, dir_flipped,
                                 w = if (EVENT_WEIGHTS == "dispute") .data[[w_col]] else 1, grp = grp_all)
  assert_no_dup_keys(ev_keys, c("ccode1", "ccode2", "year"), sprintf("events of sample %s", SAMPLE))
  stopifnot(all(ev_keys$grp %in% 0:2), sum(ev_keys$grp == 0L) + sum(ev_keys$grp == 1L) + sum(ev_keys$grp == 2L) == nrow(ev_keys))
  logf("stage-2 event groups (stage b; coalition column %s, joiner flag is_joiner_pair) over the %d sample events: bilateral %d (weights %.1f), coalition leader %d (%.1f), joiner %d (%.1f); joiner pairs whose initiator is an originator (orig1 = 1, target a joiner) %d",
       "coal_incl", nrow(ev_keys), sum(ev_keys$grp == 0L), sum(ev_keys$w[ev_keys$grp == 0L]), sum(ev_keys$grp == 1L), sum(ev_keys$w[ev_keys$grp == 1L]), sum(ev_keys$grp == 2L), sum(ev_keys$w[ev_keys$grp == 2L]),
       if ("orig1" %in% names(ev_all)) sum(ev_all$is_joiner_pair & ev_all$orig1 == 1L) else NA_integer_)
  stopifnot(all(is.finite(ev_keys$w)), all(ev_keys$w > 0), all(ev_keys$w <= 1))
  logf("direction (D1): %d of %d events in sample %s have the reversed direction (COW war initiator)", sum(ev_keys$dir_flipped), nrow(ev_keys), SAMPLE)
  ## an event whose pair was already fighting at the state-year's reference date (its earlier dispute ended before the
  ## event) has no STAGE-2 row by construction and enters stage 1 only; its (i, j, t) may keep an IV-only row (in_s2 = 0)
  far_keys <- ev_keys |> filter(fighting_at_ref)
  if (nrow(far_keys |> semi_join(alts |> filter(in_s2 == 1L), by = c("ccode1", "ccode2", "year"))) > 0) stop("fighting_at_ref events have stage-2 rows (in_s2 = 1); rebuild the data")
  logf("events in sample %s flagged fighting_at_ref (pair already fighting at the reference date; stage 1 only): %d of %d (%d keep an IV-only row)", SAMPLE, nrow(far_keys), nrow(ev_keys),
       nrow(far_keys |> semi_join(alts, by = c("ccode1", "ccode2", "year"))))
  ## the chosen alternative must be a STAGE-2 row of the state-year's block (alt_row indexes the full block)
  ev <- ev_keys |> inner_join(sy_key, by = c("ccode1", "year")) |>
    inner_join(alts |> mutate(alt_row = row_number()) |> filter(in_s2 == 1L) |> select(ccode1, ccode2, year, alt_row, dy_mid_tm1), by = c("ccode1", "ccode2", "year"))
  stopifnot(all(ev$grp %in% 0:2))
  group_table <- data.frame(group = STAGE2_GROUPS, grp = 0:2,
                            n_sample = as.integer(tabulate(ev_keys$grp + 1L, nbins = 3L)), w_sample = vapply(0:2, function(g) sum(ev_keys$w[ev_keys$grp == g]), numeric(1)),
                            n_usable = as.integer(tabulate(ev$grp + 1L, nbins = 3L)), w_usable = vapply(0:2, function(g) sum(ev$w[ev$grp == g]), numeric(1)),   # (rowsum() would drop an empty group)
                            n_states = vapply(0:2, function(g) n_distinct(ev$ccode1[ev$grp == g]), 1L), stringsAsFactors = FALSE)
  stopifnot(sum(group_table$n_usable) == nrow(ev), sum(group_table$n_sample) == nrow(ev_keys), abs(sum(group_table$w_usable) - sum(ev$w)) < 1e-9)
  logf("stage-2 event groups over the %d usable events (sample %s): bilateral %d (weights %.1f, %d initiating states), coalition leader %d (%.1f, %d), joiner %d (%.1f, %d)",
       nrow(ev), SAMPLE, group_table$n_usable[1], group_table$w_usable[1], group_table$n_states[1], group_table$n_usable[2], group_table$w_usable[2], group_table$n_states[2],
       group_table$n_usable[3], group_table$w_usable[3], group_table$n_states[3])
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
  ## major = joint: the event shares per category and the within-choice-set variation of each indicator. one_major / both_major are identified
  ## only by variation INSIDE a choice set (b_j absorbs anything constant per target; a constant per state-year drops from the softmax), so
  ## each must vary within some IV set (it enters the IV of stage 1) and within some stage-2 set of an initiating state-year with a usable
  ## event (the sets the stage-2 likelihood runs over); otherwise the run stops. For a major initiator one_major marks the non-major targets
  ## and both_major the major ones; for a non-major initiator one_major marks the major targets and both_major is 0 throughout its set.
  if (MAJOR == "joint") {
    mc_ev <- alts$major_count[ev$alt_row]
    sh_ev <- tabulate(mc_ev + 1L, nbins = 3L) / nrow(ev)
    var_in <- function(rows) {
      g <- alts[rows, c("sy_row", "one_major", "both_major")] |> group_by(sy_row) |>
        summarise(v1 = n_distinct(one_major) > 1L, v2 = n_distinct(both_major) > 1L, .groups = "drop")
      c(n_sets = nrow(g), one_major = sum(g$v1), both_major = sum(g$v2))
    }
    v_iv <- var_in(alts$in_iv == 1L)
    v_s2 <- var_in(alts$in_s2 == 1L & alts$sy_row %in% ev$sy_row)
    logf("major = joint: stage-2 events by the dyad's major-power count: 0 majors %.4f (%d), 1 major %.4f (%d), 2 majors %.4f (%d) of %d usable events",
         sh_ev[1], sum(mc_ev == 0L), sh_ev[2], sum(mc_ev == 1L), sh_ev[3], sum(mc_ev == 2L), nrow(ev))
    logf("major = joint: within-choice-set variation: IV sets (all %d state-years) with variation in one_major %d, in both_major %d; stage-2 sets of the %d initiating state-years with a usable event: one_major %d, both_major %d (both_major varies only in a major initiator's set)",
         v_iv[["n_sets"]], v_iv[["one_major"]], v_iv[["both_major"]], v_s2[["n_sets"]], v_s2[["one_major"]], v_s2[["both_major"]])
    if (v_iv[["one_major"]] == 0L || v_iv[["both_major"]] == 0L || v_s2[["one_major"]] == 0L || v_s2[["both_major"]] == 0L)
      stop(sprintf("major = joint: an indicator does not vary within any choice set (IV sets: one_major %d, both_major %d; event-bearing stage-2 sets: one_major %d, both_major %d); it is not identified in this cell",
                   v_iv[["one_major"]], v_iv[["both_major"]], v_s2[["one_major"]], v_s2[["both_major"]]))
    major_stats$share_events <- setNames(sh_ev, c("zero", "one", "both"))
    major_stats$n_events <- setNames(c(sum(mc_ev == 0L), sum(mc_ev == 1L), sum(mc_ev == 2L)), c("zero", "one", "both"))
    major_stats$var_iv_sets <- v_iv; major_stats$var_s2_sets <- v_s2
  }
  n_orphan <- sum(sy$y1 == 1L & !(seq_len(nrow(sy)) %in% ev$sy_row))
  logf("state-years with y1 = 1 but no usable stage-2 event (kept in stage 1 only): %d", n_orphan)
  ## stage-2 event weights (D8): w_<sample> from the build (or 1 under EVENT_WEIGHTS = "none"); build_stan_data() passes them as ev_w
  logf("event weights (D8, %s): %d usable events, sum of weights %.1f; weight 1 on %d events, min %.4f; weights of the events in the panel's (initiator, dispute, year) groups whose other targets left stage 2 through the choice set no longer sum to 1 (logged, not rescaled; only the pairs filter rescales, see above)",
       EVENT_WEIGHTS, nrow(ev), sum(ev$w), sum(ev$w == 1), min(ev$w))
  ev <- ev |> select(sy_row, alt_row, w, grp)        # grp: the stage-b group (0 bilateral, 1 leader, 2 joiner); build_stan_data() ignores it, 02_fit_model.R passes it for M1c only
  
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
       n_ww2 = n_ww2, drop_ww2 = DROP_WW2, dy_peace_years = DY_PEACE_YEARS, full_controls = FULL_CONTROLS, n_fighting_at_ref = nrow(far_keys),
       event_weights = EVENT_WEIGHTS, n_flipped = sum(ev_keys$dir_flipped), iv_null_cor_logn = cor(iv_null, sy$log_n_system_tm1),
       pairs = PAIRS, pairs_weights = pairs_w, changes = CHANGES, winsor = winsor, context = CONTEXT, coal_ties = COAL_TIES, coal_drop = COAL_DROP, context_top3 = context_top3,
       n_context_dropped = n_context_dropped, n_context_split_groups = n_split_groups, context_model = CONTEXT_MODEL, group_table = group_table,
       major = MAJOR, major_stats = major_stats, interact = INTERACT)
}
## assemble_cell_spec(SAMPLE, LAG, spec, ...): the main-specification constants of spec.R plus the switches of `spec`
## model = NULL or a model code: the M1c model (CONTEXT_MODELS) switches CONTEXT_MODEL on; every other model leaves it off
assemble_cell_spec <- function(SAMPLE, LAG, spec, data_dir = spec_data_dir(spec, v2_root()), manifest = NULL, model = NULL) {
  stopifnot(is.null(model) || model %in% ALL_MODELS)
  assemble_cell(SAMPLE, LAG, VARIANT, CENTRALITY, CONFLICT_HISTORY, SAMPLED_M, DISTANCE = DISTANCE, CONTIGUITY = CONTIGUITY,
                DY_PEACE_YEARS = DY_PEACE_YEARS, DROP_WW2 = spec$ww2 == "drop", FULL_CONTROLS = spec$controls == "full",
                EVENT_WEIGHTS = spec$weights, PAIRS = spec$pairs, CHANGES = spec$changes, CONTEXT = spec$context, COAL_TIES = spec$coal_ties,
                COAL_DROP = spec$coal_drop, MAJOR = spec$major, INTERACT = spec$interact,
                CONTEXT_MODEL = !is.null(model) && model %in% CONTEXT_MODELS, data_dir = data_dir, manifest = manifest)
}
