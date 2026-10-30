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
## Switch added 2026-10-28 (PI request): distance (none | capital: log1p capital distance in Z, alternatives without one dropped); suffix
## _dist; preset dist (DIST_PRESET: distance alone and distance + major). Contiguity stays in Z under either value.
## Switch added 2026-10-28 (PI request, same day): contiguity (land | none: contig_land in Z | dropped from Z); none is refused unless
## distance = capital (distance REPLACES contiguity); suffix _nocontig right after _dist; preset dist_nocontig (DIST_NOCONTIG_PRESET: the
## four dist cells without contiguity, so the two presets can be compared cell by cell).
## Switches added 2026-10-30 (PI-approved prior-sensitivity checks on the jointness result Sigma[1,2] / rho of M1b): rho_prior (beta25 | flat:
## (rho + 1) / 2 ~ Beta(2.5, 2.5) | Beta(1, 1) = uniform on (-1, 1)) and sigma_prior (ln1 | ln05: sigma_x ~ lognormal(0, 1) | lognormal(0, 0.5));
## suffixes _rhoflat, _sigln05, appended AFTER every data switch. They change no data: a non-default value routes the M1 family to the copy
## stan/m1_latent_prior.stan, which takes both hyperparameters as data (prior_hyper()); the main program and every existing fit are untouched.
## Run tag added 2026-10-30 (same request): --tag=<alnum> (spec$tag, default "") appends _t<tag> LAST to the cell name, so a refit with other
## sampler settings (e.g. --iter=6000 --warmup=3000 --tag=long) never overwrites the fit file a running job reads. Not a switch (no model change).
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
  controls  = c("minimal", "full"),                 # control sets: minimal (main) | full = main + log_peace_years (W), alliance_def_off + log_dy_peace_years (Z); no major-power status (PI 2026-10-27)
  major     = c("none", "joint"),                   # major-power status (2026-10-27): none (main) | joint = Z gains one_major / both_major (dyad's major count, 0 = reference) and W gains cowmaj_tm1
  interact  = c("none", "prefixed"),                # change interactions in Z (2026-10-27): none (main) | prefixed = the four products dM_j x dRpr_j, dM_i x dM_j, dR_i x dRpr_j, dM_i x dRpr_j
  distance  = c("none", "capital"),                 # capital distance (2026-10-28): none (main) | capital = Z gains log_capdist (log1p km at t-1, centred over the IV rows); alternatives without one are dropped from both sets
  contiguity = c("land", "none"),                   # land contiguity in Z (2026-10-28): land (main) = contig_land in Z | none = dropped from Z; none needs distance = capital (parse_spec() refuses it otherwise)
  changes   = c("raw", "winsor"),                   # initiator and target power changes: as built | winsorised at WINSOR_PROBS of the panel
  rho_prior = c("beta25", "flat"),                  # prior on Sigma's correlation rho (M1 joint models only, 2026-10-30): (rho + 1) / 2 ~ Beta(2.5, 2.5) (main) | Beta(1, 1), uniform on (-1, 1)
  sigma_prior = c("ln1", "ln05")                    # prior on the tolerance scales sigma_x (M1 family only, 2026-10-30): lognormal(0, 1) (main) | lognormal(0, 0.5)
)
SPEC_DEFAULTS <- lapply(SPEC_CHOICES, function(v) v[1])
SPEC_DEFAULTS$ame_K <- 2L                            # latent dimensions of the M2 (AME) family; "_K<K>" is always in an M2 cell name
SPEC_DEFAULTS$tag   <- ""                            # run tag (2026-10-30): "" = no tag (names unchanged); "<alnum>" appends _t<alnum> last
## the hyperparameters behind the two prior switches, passed to stan/m1_latent_prior.stan as data (prior_hyper()); the first value of each
## is the main program's literal (m1_latent_additive.stan: beta_lpdf(... | 2.5, 2.5), lognormal(0, 1))
PRIOR_RHO_SHAPE  <- c(beta25 = 2.5, flat = 1)
PRIOR_SIGMA_X_SD <- c(ln1 = 1, ln05 = 0.5)
PRIOR_SWITCHES   <- c("rho_prior", "sigma_prior")
## cell-name suffix of every NON-default value, in the order they are appended. The four that existed before 2026-10-25
## (_noww2, _fullctrl, _w1, _sideA) keep their names and relative order, so earlier output files still match.
## The context switches (2026-10-25, dual typology stage a) follow pairs: _bi / _co, then _strict (only with context != all; coal_ties is
## normalised to inclusive under context = all, so the main cell name never carries it), then _drop3 (context = coalition only).
## The stage-2 switches of 2026-10-27 follow controls (they extend the control sets): _major, then _ix, so both together give _major_ix.
## The distance switch of 2026-10-28 goes FIRST among the control extensions, right after controls and before _major / _ix: it is the only
## one that also changes the alternative rows (NA distances are dropped), and the existing _major / _ix / _major_ix names stay as they are.
## distance + major therefore gives _dist_major (DIST_PRESET).
## The contiguity switch (2026-10-28, same day) follows distance immediately: _nocontig can only appear next to _dist (contiguity = none needs
## distance = capital), so the dist_nocontig cells are _dist_nocontig and _dist_nocontig_major (DIST_NOCONTIG_PRESET); no existing name changes.
SPEC_SUFFIX <- list(
  pairs     = c(originators = "_orig", joiners = "_join"),
  context   = c(bilateral = "_bi", coalition = "_co"),
  coal_ties = c(strict = "_strict"),
  coal_drop = c(top3 = "_drop3"),
  ww2       = c(drop = "_noww2"),
  controls  = c(full = "_fullctrl"),
  distance  = c(capital = "_dist"),
  contiguity = c(none = "_nocontig"),
  major     = c(joint = "_major"),
  interact  = c(prefixed = "_ix"),
  weights   = c(none = "_w1"),
  changes   = c(winsor = "_winsor"),
  direction = c(sideA = "_sideA"),
  ## the prior switches (2026-10-30) come AFTER every data switch: they change the program's prior, not the cell's data, so every existing
  ## name is unchanged and a prior-sensitivity fit never collides with a main or preset fit (_rhoflat, _sigln05, both: _rhoflat_sigln05)
  rho_prior = c(flat = "_rhoflat"),
  sigma_prior = c(ln05 = "_sigln05")
)
## the run tag (2026-10-30) is appended by spec_suffix() after every switch suffix as _t<tag> (TAG_PREFIX), e.g. fmid_lag3_all_pr_M1b_dist_tlong;
## no switch suffix starts with "_t" followed by an alphanumeric run, so a tagged name can never be read as another switch combination
TAG_PREFIX <- "_t"
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
## stan_program(model, spec): the Stan source / build name of a model code (m1_latent_additive | m2_ame | m1_context), and, since 2026-10-30,
## m1_latent_prior for an M1 m / r / a / b model whose spec has a prior switch off its default (the copy of m1_latent_additive.stan with the
## hyperparameters as data). spec = NULL (or a spec saved before the switches existed, which lacks the fields) means the main priors.
stan_program <- function(model, spec = NULL) {
  stopifnot(model %in% ALL_MODELS)
  if (model %in% CONTEXT_MODELS) return("m1_context")
  if (substr(model, 1, 2) == "M2") return("m2_ame")
  if (spec_prior_is_main(spec)) "m1_latent_additive" else "m1_latent_prior"
}
STAN_PROGRAMS  <- c("m1_latent_additive", "m2_ame", "m1_context", "m1_latent_prior")
## spec_prior_is_main(spec): TRUE when both prior switches are at their default (a NULL spec, or one saved before 2026-10-30 without the
## fields, counts as the main priors: every fit before that date sampled the main program)
spec_prior_is_main <- function(spec) {
  if (is.null(spec)) return(TRUE)
  all(vapply(PRIOR_SWITCHES, function(nm) is.null(spec[[nm]]) || identical(spec[[nm]], SPEC_DEFAULTS[[nm]]), logical(1)))
}
## prior_hyper(spec): the data block of the prior program (prior_rho_shape, prior_sigma_x_sd) for a spec; the defaults give 2.5 and 1
prior_hyper <- function(spec) {
  rp <- if (is.null(spec$rho_prior)) SPEC_DEFAULTS$rho_prior else spec$rho_prior
  sp <- if (is.null(spec$sigma_prior)) SPEC_DEFAULTS$sigma_prior else spec$sigma_prior
  stopifnot(rp %in% names(PRIOR_RHO_SHAPE), sp %in% names(PRIOR_SIGMA_X_SD))
  list(prior_rho_shape = unname(PRIOR_RHO_SHAPE[[rp]]), prior_sigma_x_sd = unname(PRIOR_SIGMA_X_SD[[sp]]))
}
VARIANT      <- "all"                                # choice set: every other system member (polrel is an archived option)
CENTRALITY   <- "pr"                                 # IGO PageRank (eigenvector centrality is an archived option)
SAMPLED_M    <- 0L                                   # full choice sets (sampled alternatives are an archived option)
CONFLICT_HISTORY <- TRUE                             # prior-MID counts in W and Z (switching them off is an archived option)
## capital distance in Z is the `distance` switch since 2026-10-28 (spec$distance, "none" | "capital"); the former global DISTANCE constant
## is gone so that no script can set it apart from the switch (assemble_cell_spec() passes spec$distance; the cell name carries _dist)
## land contiguity in Z is the `contiguity` switch since 2026-10-28 (spec$contiguity, "land" | "none"); the former global CONTIGUITY constant
## is gone for the same reason (assemble_cell_spec() passes spec$contiguity; the cell name carries _nocontig)
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
  major_ix     = list(),                             # the major-power and interaction switches (MAJOR_IX_PRESET below): three combinations x two cells
  dist         = list(),                             # the capital-distance switch (DIST_PRESET below): distance alone and distance + major, x two cells
  dist_nocontig = list()                             # the same four cells with contiguity = none (DIST_NOCONTIG_PRESET below): distance replaces contiguity
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
## dist preset (PI request 2026-10-28: distance and distance + major powers; contiguity stays in): M1b, lag 3, samples fmid and mid; the two
## variants are distance = capital alone (_dist) and distance = capital with major = joint (_dist_major). No interactions. 05_robustness.R maps
## --cell_index = 0..3 onto variant (slowest) x cell: 0-1 distance (fmid, mid), 2-3 distance_major. M1b only, so one Slurm array, no dependency.
DIST_PRESET <- list(cells = c("fmid:3:M1b", "mid:3:M1b"),
                    variants = list(distance = list(distance = "capital"), distance_major = list(distance = "capital", major = "joint")))
## dist_nocontig preset (PI request 2026-10-28, same day: the dist cells WITHOUT contiguity, so the two can be compared): the four DIST_PRESET
## cells with contiguity = none on top (distance replaces contiguity); same cells, same variant names, same --cell_index mapping 0-1 distance
## (fmid, mid), 2-3 distance_major. Cells <sample>_lag3_all_pr_M1b_dist_nocontig and _dist_nocontig_major. M1b only, one Slurm array, no dependency.
DIST_NOCONTIG_PRESET <- list(cells = c("fmid:3:M1b", "mid:3:M1b"),
                             variants = list(distance = list(distance = "capital", contiguity = "none"),
                                             distance_major = list(distance = "capital", contiguity = "none", major = "joint")))
## default cell list of the refit presets, "sample:lag:model"; override with --cells=
ROBUSTNESS_CELLS <- c("mid:3:M1b", "fmid:3:M1b", "war:3:M1b", "mid:3:M2b", "fmid:3:M2b", "war:3:M2b")
ROBUSTNESS_MAIN_CELL <- list(sample = "mid", lag = 5L)   # the perm_error preset's cell (the a / b fits of both families)

## ---- archived CLI options (README "Archived options"; audit fix 14, 2026-10-26) ---------------------------------------------------
## Every script refuses the SAME list instead of ignoring an archived flag silently: the former --variant / --sampled_m / --centrality /
## --conflict_history / --contiguity / --dy_peace_years options and the pre-spec.R spellings of today's switches
## (--drop_ww2 -> --ww2=drop, --full_controls -> --controls=full, --event_weights -> --weights, --K -> --ame_K).
## --distance left this list on 2026-10-28: it is a switch again (none | capital), so --distance=TRUE (the archived spelling) is refused by
## parse_spec() as a value outside SPEC_CHOICES$distance, not here. --contiguity left it the same day (switch land | none): --contiguity=TRUE
## and --contiguity=FALSE are refused by parse_spec() as bad values, the bare --contiguity as a flag without a value.
ARCHIVED_OPTIONS <- c("variant", "sampled_m", "centrality", "conflict_history", "dy_peace_years",
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
  ## a bare --<switch> without =value is refused, never ignored (audit of ec9d8ae, P3 1: `--distance` alone ran the MAIN cell silently)
  for (nm in c(names(SPEC_CHOICES), "ame_K", "tag"))
    if (length(grep(paste0("^--", nm, "$"), args))) stop(sprintf("--%s needs a value (--%s=%s); a bare flag is refused", nm, nm, if (nm == "ame_K") "<K>" else if (nm == "tag") "<alnum>" else paste(SPEC_CHOICES[[nm]], collapse = " | ")))
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
  ## contiguity = none is only meaningful when the distance replaces it (PI request 2026-10-28): a Z without contiguity AND without distance is refused
  if (spec$contiguity != "land" && spec$distance != "capital") stop(sprintf("--contiguity=%s needs --distance=capital (distance replaces contiguity; got distance=%s)", spec$contiguity, spec$distance))
  ## run tag (2026-10-30): letters and digits only, so the suffix _t<tag> is one unambiguous token of the file name; "" = no tag
  spec$tag <- get_arg(args, "tag", SPEC_DEFAULTS$tag)
  if (!is.character(spec$tag) || length(spec$tag) != 1 || is.na(spec$tag)) stop("--tag must be one string")
  if (length(grep("^--tag=$", args))) stop("--tag= with an empty value is refused (leave the option out for an untagged run)")
  if (nzchar(spec$tag) && !grepl("^[A-Za-z0-9]+$", spec$tag)) stop(sprintf("--tag=%s must be letters and digits only (it becomes the cell-name suffix %s%s)", spec$tag, TAG_PREFIX, spec$tag))
  spec
}
## spec_suffix(spec): the cell-name suffix of the non-default switch values, in SPEC_SUFFIX order ("" for the main specification), then the
## run tag as _t<tag> when one is set (a spec saved before 2026-10-30 has no tag field: no tag)
spec_suffix <- function(spec) {
  out <- ""
  for (nm in names(SPEC_SUFFIX)) {
    v <- spec[[nm]]
    if (is.null(v)) next                                         # a spec saved before the switch existed is at that switch's default
    if (v != SPEC_DEFAULTS[[nm]]) { s <- SPEC_SUFFIX[[nm]][[v]]; stopifnot(!is.null(s)); out <- paste0(out, s) }
  }
  if (!is.null(spec$tag) && nzchar(spec$tag)) out <- paste0(out, TAG_PREFIX, spec$tag)
  out
}
## spec_tag(spec): the run tag ("" when none, also for a spec saved before 2026-10-30)
spec_tag <- function(spec) if (is.null(spec$tag)) "" else spec$tag
## spec_is_main(spec): TRUE when every switch is at its default. The run tag is NOT a switch (it changes no model and no data): a tagged
## main-spec run is still the main specification, under another file name; callers that must not see a tag check spec_tag() themselves.
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
## The prior switches (2026-10-30) apply to the M1 m / r / a / b programs only: with M2 or M1c (other programs, no hyperparameter data) a
## non-default value is refused here, so no name of those families ever carries _rhoflat / _sigln05, and --rho_prior=flat with a non-joint
## M1 model (M1m, M1r, M1a: no rho_x) is refused as well, since the fit would be the main model under a name that claims another prior.
cell_name <- function(sample, lag, model, spec) {
  stopifnot(sample %in% ALL_SAMPLES, lag %in% LAGS, model %in% ALL_MODELS)
  if (model_family(model) == "M1" && !identical(as.integer(spec$ame_K), SPEC_DEFAULTS$ame_K))
    stop(sprintf("cell_name: --ame_K=%s with the M1 model %s (K applies to the M2 family only and would not appear in the cell name)", spec$ame_K, model))
  if (!spec_prior_is_main(spec)) {
    off <- PRIOR_SWITCHES[vapply(PRIOR_SWITCHES, function(nm) !is.null(spec[[nm]]) && !identical(spec[[nm]], SPEC_DEFAULTS[[nm]]), logical(1))]
    if (stan_program(model) != "m1_latent_additive")
      stop(sprintf("cell_name: --%s=%s with the model %s (the prior switches apply to the M1 m / r / a / b programs only; %s has no such hyperparameter)", off[1], spec[[off[1]]], model, stan_program(model)))
    if ("rho_prior" %in% off && model_joint(model) != 1L)
      stop(sprintf("cell_name: --rho_prior=%s with the model %s, which has no Sigma correlation rho (the rho prior applies to the joint model M1b only)", spec$rho_prior, model))
  }
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
