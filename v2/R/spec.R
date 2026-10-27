## ============================================================================
## spec.R -- THE main specification and its switches. Every v2 script sources this file first.
##
## One block defines (a) the main specification (the defaults below), (b) the named switches a run can flip from the
## command line as --name=value, (c) the cell-name suffix every non-default value adds, so a variant NEVER overwrites a
## main-spec output, and (d) the sampler settings per model family. Nothing here has a side effect; the scripts call
## parse_spec() on their command-line arguments and read the constants.
##
## Switches added 2026-10-25 (dual typology stage a): context (all | bilateral | coalition, STAGE-2 events only), coal_ties
## (inclusive | strict), coal_drop (none | top3, context = coalition only); suffixes _bi / _co, _strict, _drop3.
## Model code added 2026-10-27 (dual typology stage b): M1c (CONTEXT_MODELS, stan/m1_context.stan), outside MODELS and the main grid;
## the context_b preset runs it at lag 3 on mid and fmid (CONTEXT_B_PRESET).
## Switches added 2026-10-27 (PI-approved stage-2 extensions): major (none | joint: one_major / both_major in Z, cowmaj_tm1 in W) and
## interact (none | prefixed: four change products in Z); suffixes _major, _ix; preset major_ix (MAJOR_IX_PRESET).
## Main specification (PI decisions through 2026-10-25; README "The main specification"):
##   sample in {mid, fmid, war} x lag in {1, 3, 5, 10} x model in {M1m, M1r, M1a, M1b, M2m, M2r, M2a, M2b}
##   choice set = every other system member at t with positions at t-1 (variant "all"), full sets (no sampled alternatives)
##   relational position = IGO PageRank; prior-MID counts in W and Z; land contiguity in Z; no capital distance; no dyadic
##   peace years (these were CLI options before 2026-10-25 and are now fixed constants; see README "Archived options")
## ============================================================================
V2_SEED <- 123L

## ---- switches: name -> allowed values; the FIRST value is the main specification ---------------------------------------------
SPEC_CHOICES <- list(
  pairs     = c("all", "originators", "joiners"),   # event pairs: all | originator pairs only | joiner pairs only (stage-2 comparison)
  context   = c("all", "bilateral", "coalition"),   # STAGE-2 events by coalition context of their D8 group (coal_<ties> in events.rds): all | bilateral only | coalition only; stage 1 unchanged
  coal_ties = c("inclusive", "strict"),             # tie convention of the context flag: same-date entrants present (coal_incl) | strictly earlier entrants only (coal_strict); inert with context = all
  coal_drop = c("none", "top3"),                    # context = coalition only: drop the 3 disputes contributing the most coalition events of the sample (leverage check)
  direction = c("cowwar", "sideA"),                 # initiator direction: COW war initiator in war-level disputes (v2/data) | MID Side A (v2/data_sideA)
  ww2       = c("include", "drop"),                 # the Second World War events: kept | removed from stage 2 and stage 1
  weights   = c("dispute", "none"),                 # stage-2 event weights: 1 / (initiator's targets in the dispute-year) | 1 everywhere
  controls  = c("minimal", "full"),                 # control sets: minimal (main) | the full base specification of commit 8bf0ce1
  major     = c("none", "joint"),                   # major-power status (2026-10-27): none (main) | joint = Z gains one_major / both_major (dyad's major count, 0 = reference) and W gains cowmaj_tm1
  interact  = c("none", "prefixed"),                # change interactions in Z (2026-10-27): none (main) | prefixed = the four products dM_j x dRpr_j, dM_i x dM_j, dR_i x dRpr_j, dM_i x dRpr_j
  changes   = c("raw", "winsor")                    # initiator and target power changes: as built | winsorised at WINSOR_PROBS of the panel
)
SPEC_DEFAULTS <- lapply(SPEC_CHOICES, function(v) v[1])
SPEC_DEFAULTS$ame_K <- 2L                            # latent dimensions of the M2 (AME) family; "_K<K>" is always in an M2 cell name
## cell-name suffix of every NON-default value, in the order they are appended. The four that existed before 2026-10-25
## (_noww2, _fullctrl, _w1, _sideA) keep their names and relative order, so earlier output files still match.
## The context switches (2026-10-25, dual typology stage a) follow pairs: _bi / _co, then _strict (only with context != all; coal_ties is
## normalised to inclusive under context = all, so the main cell name never carries it), then _drop3 (context = coalition only).
## The stage-2 switches of 2026-10-27 follow controls (they extend the control sets): _major, then _ix, so both together give _major_ix.
SPEC_SUFFIX <- list(
  pairs     = c(originators = "_orig", joiners = "_join"),
  context   = c(bilateral = "_bi", coalition = "_co"),
  coal_ties = c(strict = "_strict"),
  coal_drop = c(top3 = "_drop3"),
  ww2       = c(drop = "_noww2"),
  controls  = c(full = "_fullctrl"),
  major     = c(joint = "_major"),
  interact  = c(prefixed = "_ix"),
  weights   = c(none = "_w1"),
  changes   = c(winsor = "_winsor"),
  direction = c(sideA = "_sideA")
)
## winsorisation caps (changes = winsor): the 1st and 99th percentiles of the stage-1 panel's transformed change (dM_h, dRpr_h)
## for the cell's lag, applied to dM_i, dR_i (stage 1 / mu shifters) and dM_j_h, dRpr_j_h (target changes in Z) at assembly
WINSOR_PROBS <- c(0.01, 0.99)

## ---- fixed design constants (not switches) -----------------------------------------------------------------------------------
MAIN_SAMPLES <- c("mid", "fmid", "war")              # the main grid; fmidcow and fmidub are appendix samples a fit can still take
ALL_SAMPLES  <- c("mid", "fmid", "fmidcow", "fmidub", "war")
LAGS         <- c(1L, 3L, 5L, 10L)
MODELS       <- c("M1m", "M1r", "M1a", "M1b", "M2m", "M2r", "M2a", "M2b")
## dual typology stage (b), 2026-10-27: M1c = the M1b joint model with group-specific ideal-point shifts (delta_lead, delta_join) and
## group-specific tolerance scales (sigma_g, shared correlation rho) over the three stage-2 event groups bilateral / coalition leader /
## joiner, plus one extra stage-1 control (multi_mid_jan). Its own Stan program (stan/m1_context.stan). NOT part of MODELS (the main
## grid's 0-95 index mapping and 04_compare.R's eight-model groups stay as they are); 02_fit_model.R and 03_extract.R accept it.
CONTEXT_MODELS <- c("M1c")
ALL_MODELS     <- c(MODELS, CONTEXT_MODELS)
## stan_program(model): the Stan source / build name of a model code (m1_latent_additive | m2_ame | m1_context)
stan_program <- function(model) { stopifnot(model %in% ALL_MODELS); if (model %in% CONTEXT_MODELS) "m1_context" else if (substr(model, 1, 2) == "M1") "m1_latent_additive" else "m2_ame" }
STAN_PROGRAMS  <- c("m1_latent_additive", "m2_ame", "m1_context")
VARIANT      <- "all"                                # choice set: every other system member (polrel is an archived option)
CENTRALITY   <- "pr"                                 # IGO PageRank (eigenvector centrality is an archived option)
SAMPLED_M    <- 0L                                   # full choice sets (sampled alternatives are an archived option)
CONFLICT_HISTORY <- TRUE                             # prior-MID counts in W and Z (switching them off is an archived option)
DISTANCE     <- FALSE                                # capital distance in Z (archived option)
CONTIGUITY   <- TRUE                                 # land contiguity in Z (switching it off is an archived option)
DY_PEACE_YEARS <- FALSE                              # dyadic peace years in Z (archived option)
R0           <- 0.5                                  # typology peer radius in scaled gap units (03_extract.R default)

## ---- sampler settings per model family (02_fit_model.R defaults; --chains / --iter / --warmup / --adapt_delta / --max_treedepth
## override them). PI decision 2026-10-25 (M2 convergence): the M2 family runs 4000 iterations with 2000 warmup and adapt_delta 0.95;
## M1 and max_treedepth unchanged. The main-spec M2 cell names stay the same. ----------------------------------------------------
SAMPLER <- list(
  M1 = list(chains = 4L, iter = 2000L, warmup = 1000L, adapt_delta = 0.90, max_treedepth = 10L),
  M2 = list(chains = 4L, iter = 4000L, warmup = 2000L, adapt_delta = 0.95, max_treedepth = 10L)
)
## M2 initial values (PI decision 2026-10-25): every M2 chain starts at the posterior means of the matching M1 fit (same sample, lag
## and switches; M2x from M1x), each chain jittered by N(0, INIT_JITTER_SD) on the unconstrained parameters (exp(N(0, INIT_JITTER_SD))
## on the positive scales, N(0, INIT_JITTER_SD / 5) on the correlations) and the AME factors drawn uniformly from INIT_UV_RANGE
## (u_diag on its positive upper half), seed V2_SEED + chain. 02_fit_model.R --init_from=<fit file> names the M1 fit; by default it is
## located in v2/output (fit_<cell with M1x>.rds); when none is found the fit falls back to Stan's random inits and says so loudly in the
## log, the fit object (obj$init) and the summary (init_source).
INIT_JITTER_SD <- 0.05
INIT_UV_RANGE  <- c(-0.1, 0.1)

## ---- robustness presets (05_robustness.R): preset -> the switch it flips ------------------------------------------------------
ROBUSTNESS_PRESETS <- list(
  noww2        = list(ww2 = "drop"),
  sideA        = list(direction = "sideA"),
  noweights    = list(weights = "none"),
  fullcontrols = list(controls = "full"),
  winsor       = list(changes = "winsor"),
  perm_error   = list(),                             # the permutation and measurement-error tests on the main cell (no switch)
  context      = list(),                             # the context split (CONTEXT_PRESET below): three switch combinations x two cells
  context_b    = list(),                             # dual typology stage (b): the M1c model on the main-spec cells of CONTEXT_B_PRESET (no switch)
  major_ix     = list()                              # the major-power and interaction switches (MAJOR_IX_PRESET below): three combinations x two cells
)
## context_b preset (dual typology stage b, 2026-10-27): the M1c model (one model, group-specific shifts) at lag 3 on mid and fmid with every
## switch at its default; 05_robustness.R maps --cell_index = 0..1 onto the cell list (0 mid, 1 fmid). M1b family only, no war sample.
CONTEXT_B_PRESET <- list(cells = c("mid:3:M1c", "fmid:3:M1c"))
## context preset (dual typology stage a, 2026-10-25): M1b, lag 3, samples mid and fmid, default tie convention; the three variants are
## context = bilateral, context = coalition, and context = coalition without the top-3 coalition disputes (coal_drop = top3). 05_robustness.R
## maps --cell_index = 0..5 onto variant (slowest) x cell: 0-1 bilateral, 2-3 coalition, 4-5 coalition_drop3. No war sample.
CONTEXT_PRESET <- list(cells = c("mid:3:M1b", "fmid:3:M1b"),
                       variants = list(bilateral = list(context = "bilateral"), coalition = list(context = "coalition"),
                                       coalition_drop3 = list(context = "coalition", coal_drop = "top3")))
## major_ix preset (PI-approved stage-2 switches, 2026-10-27): M1b, lag 3, samples fmid and mid; the three variants are major = joint alone
## (_major), interact = prefixed alone (_ix), and both (_major_ix). 05_robustness.R maps --cell_index = 0..5 onto variant (slowest) x cell:
## 0-1 major (fmid, mid), 2-3 interact, 4-5 both. M1b only, so no M1 -> M2 dependency; one Slurm array.
MAJOR_IX_PRESET <- list(cells = c("fmid:3:M1b", "mid:3:M1b"),
                        variants = list(major = list(major = "joint"), interact = list(interact = "prefixed"),
                                        both = list(major = "joint", interact = "prefixed")))
## default cell list of the refit presets, "sample:lag:model"; override with --cells=
ROBUSTNESS_CELLS <- c("mid:3:M1b", "fmid:3:M1b", "war:3:M1b", "mid:3:M2b", "fmid:3:M2b", "war:3:M2b")
ROBUSTNESS_MAIN_CELL <- list(sample = "mid", lag = 5L)   # the perm_error preset's cell (the a / b fits of both families)

## ---- archived CLI options (README "Archived options"; audit fix 14, 2026-10-26) ---------------------------------------------------
## Every script refuses the SAME list instead of ignoring an archived flag silently: the former --variant / --sampled_m / --centrality /
## --conflict_history / --distance / --contiguity / --dy_peace_years options and the pre-spec.R spellings of today's switches
## (--drop_ww2 -> --ww2=drop, --full_controls -> --controls=full, --event_weights -> --weights, --K -> --ame_K).
ARCHIVED_OPTIONS <- c("variant", "sampled_m", "centrality", "conflict_history", "distance", "contiguity", "dy_peace_years",
                      "drop_ww2", "full_controls", "event_weights", "K")
## refuse_archived_options(args, script, also, also_msg): stop() when args hold --<archived>=... (or the bare flag); `also` adds
## script-specific refusals with their own message (the build takes --direction only, so it refuses every other switch of SPEC_CHOICES)
refuse_archived_options <- function(args, script, also = character(0), also_msg = "it is not an option of this script") {
  for (old in ARCHIVED_OPTIONS)
    if (length(grep(paste0("^--", old, "(=|$)"), args)))
      stop(sprintf("--%s is not an option of %s (archived 2026-10-25); use the switches of spec.R (--ww2, --controls, --weights, --ame_K, ...) or see README 'Archived options'", old, script))
  for (old in also)
    if (length(grep(paste0("^--", old, "(=|$)"), args))) stop(sprintf("--%s is not an option of %s (%s)", old, script, also_msg))
  invisible(TRUE)
}

## ---- parsing -----------------------------------------------------------------------------------------------------------------
## get_arg(args, name, default): the value of --name=value in args, typed like default (logical / numeric / character)
get_arg <- function(args, name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  val <- sub(paste0("^--", name, "="), "", hit[1])
  if (is.logical(default)) { out <- as.logical(val); if (is.na(out)) stop(sprintf("--%s must be TRUE or FALSE, got '%s'", name, val)); return(out) }
  if (is.numeric(default)) { out <- suppressWarnings(as.numeric(val)); if (is.na(out)) stop(sprintf("--%s must be numeric, got '%s'", name, val)); return(out) }
  val
}
## parse_spec(args): the switches from the command line, defaults = the main specification; every value is validated
parse_spec <- function(args = character(0)) {
  spec <- SPEC_DEFAULTS
  for (nm in names(SPEC_CHOICES)) {
    spec[[nm]] <- get_arg(args, nm, SPEC_DEFAULTS[[nm]])
    if (!spec[[nm]] %in% SPEC_CHOICES[[nm]])
      stop(sprintf("--%s must be one of %s, got '%s'", nm, paste(SPEC_CHOICES[[nm]], collapse = " | "), spec[[nm]]))
  }
  spec$ame_K <- as.integer(get_arg(args, "ame_K", SPEC_DEFAULTS$ame_K))
  if (is.na(spec$ame_K) || spec$ame_K < 1L) stop("--ame_K must be a positive integer")
  ## context switches: coal_drop needs the coalition context; coal_ties has no effect without a context filter and is normalised to its
  ## default there, so context = all is byte-identical to the main specification whatever --coal_ties says
  if (spec$coal_drop != "none" && spec$context != "coalition") stop(sprintf("--coal_drop=%s needs --context=coalition (got context=%s)", spec$coal_drop, spec$context))
  if (spec$context == "all") spec$coal_ties <- SPEC_DEFAULTS$coal_ties
  spec
}
## spec_suffix(spec): the cell-name suffix of the non-default switch values, in SPEC_SUFFIX order ("" for the main specification)
spec_suffix <- function(spec) {
  out <- ""
  for (nm in names(SPEC_SUFFIX)) {
    v <- spec[[nm]]
    if (v != SPEC_DEFAULTS[[nm]]) { s <- SPEC_SUFFIX[[nm]][[v]]; stopifnot(!is.null(s)); out <- paste0(out, s) }
  }
  out
}
## spec_is_main(spec): TRUE when every switch is at its default
spec_is_main <- function(spec) all(vapply(names(SPEC_CHOICES), function(nm) identical(spec[[nm]], SPEC_DEFAULTS[[nm]]), logical(1))) && spec$ame_K == SPEC_DEFAULTS$ame_K
## spec_data_dir(spec, root): v2/data for the main direction, v2/data_sideA for direction = sideA
spec_data_dir <- function(spec, root) file.path(root, if (spec$direction == "sideA") "data_sideA" else "data")
## spec_string(spec): one line for the logs
spec_string <- function(spec) paste(sprintf("%s=%s", names(spec), unlist(spec)), collapse = " ")
## spec_coal_col(spec): the events.rds column the context switch filters on ("coal_incl" / "coal_strict")
spec_coal_col <- function(spec) c(inclusive = "coal_incl", strict = "coal_strict")[[spec$coal_ties]]
## cell_name(sample, lag, model, spec): <sample>_lag<h>_all_pr_<model>[_K<K>]<suffix>, e.g. fmid_lag3_all_pr_M2b_K2_noww2
## An M1 cell name carries no K, so ame_K != 2 with an M1 model would be dropped silently: refused (audit fix 16, 2026-10-26).
## M1c (stage b, 2026-10-27) is named like an M1 model, <sample>_lag<h>_all_pr_M1c<suffix>: no existing cell name contains "M1c".
cell_name <- function(sample, lag, model, spec) {
  stopifnot(sample %in% ALL_SAMPLES, lag %in% LAGS, model %in% ALL_MODELS)
  if (model_family(model) == "M1" && !identical(as.integer(spec$ame_K), SPEC_DEFAULTS$ame_K))
    stop(sprintf("cell_name: --ame_K=%s with the M1 model %s (K applies to the M2 family only and would not appear in the cell name)", spec$ame_K, model))
  sprintf("%s_lag%d_%s_%s_%s%s%s", sample, as.integer(lag), VARIANT, CENTRALITY, model,
          if (grepl("^M2", model)) sprintf("_K%d", spec$ame_K) else "", spec_suffix(spec))
}
## compare_group_name(sample, lag, spec): the <group> of 04_compare.R, <sample>_lag<h>_all_pr[_K<K>]<suffix>; "_K<K>" only for a non-default
## ame_K, so the documented main-grid names (mid_lag5_all_pr, fmid_lag3_all_pr_noww2, ...) are unchanged and a K=3 comparison never
## overwrites the K=2 one (review 2026-10-26, P2)
compare_group_name <- function(sample, lag, spec) {
  stopifnot(sample %in% ALL_SAMPLES, lag %in% LAGS)
  sprintf("%s_lag%d_%s_%s%s%s", sample, as.integer(lag), VARIANT, CENTRALITY,
          if (!identical(as.integer(spec$ame_K), SPEC_DEFAULTS$ame_K)) sprintf("_K%d", spec$ame_K) else "", spec_suffix(spec))
}
## measurement_cell_name(sample, lag, family, spec): the <cell> of the 05_robustness.R perm_error outputs (measurement_<cell>_*.csv),
## <sample>_lag<h>_all_pr_<M1|M2>[_K<K>]<suffix>. The M2 name ALWAYS carries its K, as cell_name() does for M2 fits (review 2026-10-26,
## P2: K=2 and K=3 wrote the same mid_lag5_all_pr_M2 files); the M1 name never does.
measurement_cell_name <- function(sample, lag, family, spec) {
  stopifnot(sample %in% ALL_SAMPLES, lag %in% LAGS, family %in% c("m1", "m2"))
  sprintf("%s_lag%d_%s_%s_%s%s%s", sample, as.integer(lag), VARIANT, CENTRALITY, toupper(family),
          if (family == "m2") sprintf("_K%d", spec$ame_K) else "", spec_suffix(spec))
}
## m1_spec(spec): the spec of the M1 cell that goes with an M2 cell, i.e. the same switches with ame_K reset to the default (K belongs to
## the M2 family only; cell_name() refuses a non-default K with an M1 model). The M1 init lookup of an M2 fit and every M1 name in the
## comparison / measurement tools derive their names from this (review 2026-10-26, P2: `--model=M2b --ame_K=3` threw before checking
## whether the M1 file existed).
m1_spec <- function(spec) { spec$ame_K <- SPEC_DEFAULTS$ame_K; spec }
## cell_spec(model, spec): the spec whose cell name belongs to `model` under `spec`: m1_spec(spec) for the M1 family, spec for M2
cell_spec <- function(model, spec) if (model_family(model) == "M1") m1_spec(spec) else spec
## model_family(model): "M1" or "M2"; model_D_mode(model): 1 (material), 2 (relational), 3 (both); model_joint(model): 1 for the b models
## M1c (stage b): family M1 (sampler settings, name without a K), both axes, joint (full Sigma per group, full Omega)
model_family <- function(model) { stopifnot(model %in% ALL_MODELS); substr(model, 1, 2) }
model_D_mode <- function(model) { stopifnot(model %in% ALL_MODELS); switch(substr(model, 3, 3), m = 1L, r = 2L, a = 3L, b = 3L, c = 3L) }
model_joint  <- function(model) { stopifnot(model %in% ALL_MODELS); as.integer(grepl("[bc]$", model)) }
## sampler_settings(model, args): the family defaults, each overridable from the command line
sampler_settings <- function(model, args = character(0)) {
  s <- SAMPLER[[model_family(model)]]
  list(chains = as.integer(get_arg(args, "chains", s$chains)), iter = as.integer(get_arg(args, "iter", s$iter)),
       warmup = as.integer(get_arg(args, "warmup", s$warmup)), adapt_delta = get_arg(args, "adapt_delta", s$adapt_delta),
       max_treedepth = as.integer(get_arg(args, "max_treedepth", s$max_treedepth)))
}
