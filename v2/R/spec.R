## ============================================================================
## spec.R -- THE main specification and its switches. Every v2 script sources this file first.
##
## One block defines (a) the main specification (the defaults below), (b) the named switches a run can flip from the
## command line as --name=value, (c) the cell-name suffix every non-default value adds, so a variant NEVER overwrites a
## main-spec output, and (d) the sampler settings per model family. Nothing here has a side effect; the scripts call
## parse_spec() on their command-line arguments and read the constants.
##
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
  direction = c("cowwar", "sideA"),                 # initiator direction: COW war initiator in war-level disputes (v2/data) | MID Side A (v2/data_sideA)
  ww2       = c("include", "drop"),                 # the Second World War events: kept | removed from stage 2 and stage 1
  weights   = c("dispute", "none"),                 # stage-2 event weights: 1 / (initiator's targets in the dispute-year) | 1 everywhere
  controls  = c("minimal", "full"),                 # control sets: minimal (main) | the full base specification of commit 8bf0ce1
  changes   = c("raw", "winsor")                    # initiator and target power changes: as built | winsorised at WINSOR_PROBS of the panel
)
SPEC_DEFAULTS <- lapply(SPEC_CHOICES, function(v) v[1])
SPEC_DEFAULTS$ame_K <- 2L                            # latent dimensions of the M2 (AME) family; "_K<K>" is always in an M2 cell name
## cell-name suffix of every NON-default value, in the order they are appended. The four that existed before 2026-10-25
## (_noww2, _fullctrl, _w1, _sideA) keep their names and relative order, so earlier output files still match.
SPEC_SUFFIX <- list(
  pairs     = c(originators = "_orig", joiners = "_join"),
  ww2       = c(drop = "_noww2"),
  controls  = c(full = "_fullctrl"),
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
VARIANT      <- "all"                                # choice set: every other system member (polrel is an archived option)
CENTRALITY   <- "pr"                                 # IGO PageRank (eigenvector centrality is an archived option)
SAMPLED_M    <- 0L                                   # full choice sets (sampled alternatives are an archived option)
CONFLICT_HISTORY <- TRUE                             # prior-MID counts in W and Z (switching them off is an archived option)
DISTANCE     <- FALSE                                # capital distance in Z (archived option)
CONTIGUITY   <- TRUE                                 # land contiguity in Z (switching it off is an archived option)
DY_PEACE_YEARS <- FALSE                              # dyadic peace years in Z (archived option)
R0           <- 0.5                                  # typology peer radius in scaled gap units (03_extract.R default)

## ---- sampler settings per model family (02_fit_model.R defaults; --chains / --iter / --warmup / --adapt_delta / --max_treedepth
## override them) ---------------------------------------------------------------------------------------------------------------
SAMPLER <- list(
  M1 = list(chains = 4L, iter = 2000L, warmup = 1000L, adapt_delta = 0.90, max_treedepth = 10L),
  M2 = list(chains = 4L, iter = 2000L, warmup = 1000L, adapt_delta = 0.90, max_treedepth = 10L)
)

## ---- robustness presets (05_robustness.R): preset -> the switch it flips ------------------------------------------------------
ROBUSTNESS_PRESETS <- list(
  noww2        = list(ww2 = "drop"),
  sideA        = list(direction = "sideA"),
  noweights    = list(weights = "none"),
  fullcontrols = list(controls = "full"),
  winsor       = list(changes = "winsor"),
  perm_error   = list()                              # the permutation and measurement-error tests on the main cell (no switch)
)
## default cell list of the refit presets, "sample:lag:model"; override with --cells=
ROBUSTNESS_CELLS <- c("mid:3:M1b", "fmid:3:M1b", "war:3:M1b", "mid:3:M2b", "fmid:3:M2b", "war:3:M2b")
ROBUSTNESS_MAIN_CELL <- list(sample = "mid", lag = 5L)   # the perm_error preset's cell (the a / b fits of both families)

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
## cell_name(sample, lag, model, spec): <sample>_lag<h>_all_pr_<model>[_K<K>]<suffix>, e.g. fmid_lag3_all_pr_M2b_K2_noww2
cell_name <- function(sample, lag, model, spec) {
  stopifnot(sample %in% ALL_SAMPLES, lag %in% LAGS, model %in% MODELS)
  sprintf("%s_lag%d_%s_%s_%s%s%s", sample, as.integer(lag), VARIANT, CENTRALITY, model,
          if (grepl("^M2", model)) sprintf("_K%d", spec$ame_K) else "", spec_suffix(spec))
}
## model_family(model): "M1" or "M2"; model_D_mode(model): 1 (material), 2 (relational), 3 (both); model_joint(model): 1 for the b models
model_family <- function(model) { stopifnot(model %in% MODELS); substr(model, 1, 2) }
model_D_mode <- function(model) { stopifnot(model %in% MODELS); switch(substr(model, 3, 3), m = 1L, r = 2L, a = 3L, b = 3L) }
model_joint  <- function(model) { stopifnot(model %in% MODELS); as.integer(grepl("b$", model)) }
## sampler_settings(model, args): the family defaults, each overridable from the command line
sampler_settings <- function(model, args = character(0)) {
  s <- SAMPLER[[model_family(model)]]
  list(chains = as.integer(get_arg(args, "chains", s$chains)), iter = as.integer(get_arg(args, "iter", s$iter)),
       warmup = as.integer(get_arg(args, "warmup", s$warmup)), adapt_delta = get_arg(args, "adapt_delta", s$adapt_delta),
       max_treedepth = as.integer(get_arg(args, "max_treedepth", s$max_treedepth)))
}
