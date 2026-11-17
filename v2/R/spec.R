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
  sigma_prior = c("ln1", "ln05"),                   # prior on the tolerance scales sigma_x (M1 family only, 2026-10-30): lognormal(0, 1) (main) | lognormal(0, 0.5)
  comparison = c("none", "conventional", "directed"), # the conventional comparison specification (PI-approved 2026-10-31): none (main) | conventional = W, Z, X fixed to W_CONVENTIONAL / Z_CONVENTIONAL / no X (P = 0) | directed (PI-approved 2026-11-01) = conventional with the directed Z_DIRECTED; see below
  roles     = c("stageb", "aprime", "bc", "bcj"),                # stage-2 role groups of the M1c / M1cr models (PI decision 2026-10-31): stageb (main) = the 3 stage-b groups (bilateral / leader / joiner) | aprime = 4 groups (alone originator / co-originator / joiner / target joined late; roles_aprime.rds, aprime_roles())
  mu_shift  = c("none", "traits"),                  # ideal-point shifters under the directed comparison (PI-approved 2026-11-02): none (main; X = none, P = 0) | traits = X = the initiator's M_lvl, R_lvl, cowmaj_tm1, dem_tm1 at t-1 as they enter W (X_MU_TRAITS, assemble_cell.R; P = 4); needs --comparison=directed
  joiner_set = c("all", "opposing"),                 # joiner choice set (PI-approved 2026-11-04): all (main; every stage-2 alternative) | opposing = a joiner event (aprime group 2) is scored only against the states on the OTHER side of the dispute it joins in the onset year, intersected with the state-year's stage-2 set (joiner_sets_aprime.rds); needs --roles=aprime
  entry     = c("keep", "later"),                    # joiner entry direction (PI decision 2026-11-04): keep (main; v2/data) | later = a joiner initiates against the side it joins, two joiners: the later entrant (v2/data_entry, built by 01_build_data.R --entry=later); suffix _ent
  events    = c("pairs", "dymid"),                   # event foundation (PI-approved 2026-11-04): pairs (main; every Side A x Side B pair of a MID, v2/data) | dymid = the interacting dyads of Maoz's Dyadic MID 4.03 + Directed Dyadic War data, with the entry rule built in (v2/data_dymid, built by 01_build_data.R --events=dymid; README section 6.9); suffix _dy
  variant   = c("all", "polrel", "mid3"),            # choice set (un-archived 2026-11-04 for the role / AME models only): all (main) | polrel = the politically relevant alternatives (contiguous, conttype 1..5, or either state a major power; alts_lag<h>_polrel.rds of the build, data or data_entry); no suffix: it takes the "all" slot of the cell name
                                                    # mid3 (2026-11-10): the lag-3 MID-history choice set, built at assembly from alts_lag<h>_all.rds and events.rds (below)
  reach     = c("none", "capdist"),                 # run 0a (2026-11-10): none (main) | capdist = the directed Z plus log_capdist^2 and log_capdist x the initiator's M_lvl (Z_REACH); needs --comparison=directed, model M1c; suffix _reach
  subsample = c("none", "prtarget", "bilateral", "majinit", "mininit"),  # run 0b-0d (2026-11-10): none (main) | prtarget = stage-2 events whose chosen target is politically relevant (full choice set; stage 1 unchanged) | bilateral = stage-2 events of bcj role bilateral (M1b; stage 1 unchanged) | majinit / mininit = BOTH stages restricted to the state-years of initiators that are / are not COW major powers at t-1; needs --comparison=directed
  opp_weight = c("none", "0", "0.5", "1"),          # run 1 (2026-11-10): none (main) | lambda = 0 | 0.5 | 1: the fixed offset lambda x ln q_ij of the MID-propensity opportunity weight in the stage-2 utility AND inside the stage-1 IV (M1c, stan/m1_context_offset.stan); suffixes _opw0 / _opw05 / _opw1
  opp_outcome = c("mid", "host3"),                  # run 1: the auxiliary outcome behind q: mid = any Dyadic MID onset (main) | host3 = onsets with dyadic hostility >= 3; suffix _oh3; needs an opp_weight
  opp_draw  = c("0", "1", "2", "3", "4", "5")       # run 1: 0 = q at the auxiliary MLE (main) | d = the d-th coherent full-vector draw of the auxiliary coefficients; suffix _d<d>; needs --opp_weight=1
)
## ---- the conventional comparison specification (PI-approved 2026-10-31; switch `comparison`, suffix _conv) --------------------------
## One switch fixes the three design matrices at once; it is NOT a combination of the other switches, although it carries the distance
## (log_capdist in Z, alternatives without one dropped, as distance = capital), the major-power indicators (one_major / both_major in Z and
## cowmaj_tm1 in W, as major = joint), land contiguity and the defense-alliance control of the full set; no peace years in either stage (PI decision
## 2026-10-31; lagged conflict and the prior-MID counts stay). The power CHANGES leave every matrix: no dM_i /
## dR_i in W, no X at all (P = 0: mu = mu0 + a_i, or + delta[g] for M1c / M1cr), no dM_j_h / dRpr_j_h / dch_j_missing in Z.
## multi_mid_jan is in W for EVERY model of the spec (M1a, M1b, M1c, M1cr), not for M1c alone as in the main specification. The column lists are
## W_CONVENTIONAL / Z_CONVENTIONAL in assemble_cell.R (asserted exactly at assembly and in 02_fit_model.R). Models: CONVENTIONAL_MODELS only
## (cell_name() refuses the others). The switches the spec fixes (CONVENTIONAL_FIXED) are refused on the command line in ANY value, the
## default included (--controls=minimal with --comparison=conventional is a conflict, not a no-op): the W / Z sets are the spec's, not theirs.
## The spec fields of those switches stay at their defaults (spec$distance = "none" etc.); assemble_cell_spec() maps the conventional
## spec onto DISTANCE = "capital", MAJOR = "joint" and CONVENTIONAL = TRUE, so the existing distance / major machinery (NA-distance drop,
## event reconciliation, within-set variation checks) runs unchanged and the fit records the assembled values (args$distance = "capital").
## The data switches that do not touch W / Z / X (pairs, context, coal_*, ww2, weights, direction) and the prior switches combine with it as
## with any cell; the suffix _conv sits after _ix (the last control extension) and before _w1, so e.g. fmid_lag3_all_pr_M1b_conv,
## fmid_lag3_all_pr_M1b_noww2_conv, fmid_lag3_all_pr_M1b_conv_w1. 05_robustness.R refuses the switch for every preset (single cells).
CONVENTIONAL_MODELS <- c("M1a", "M1b", "M1c", "M1cr")       # M1cr (role-specific rho, 2026-10-31) added with its program
CONVENTIONAL_FIXED  <- c("controls", "distance", "contiguity", "major", "interact", "changes")
## ---- the directed comparison specification (PI-approved 2026-11-01; comparison = directed, suffix _cdir) ------------------------------
## Identical to the conventional comparison (same W = W_CONVENTIONAL, no X / P = 0, distance = capital and major = joint fixed internally, lag 3
## and variant all enforced, the same switches refused in any value: CONVENTIONAL_FIXED) EXCEPT stage-2 Z = Z_DIRECTED (assemble_cell.R):
## contig_land, log_capdist, tgt_dem, tgt_dem_x_init_dem, alliance_def_off, tgt_major, tgt_major_x_init_major, dy_mid_tm1, log_dy_prior_mids
## (Pz = 9). In a conditional logit over one initiator's targets an initiator-only attribute cancels, so the conventional joint_dem and
## one_major / both_major are replaced by the target's own status and its product with the initiator's: tgt_dem = dem2, tgt_dem_x_init_dem =
## dem2 * dem1 (= joint_dem, asserted), tgt_major = major2, tgt_major_x_init_major = major2 * major1 (= both_major; one_major = tgt + init -
## 2 tgt init, asserted). Models: DIRECTED_MODELS, i.e. the four M1 models of the conventional spec plus the AME model M2b (no other M2 model);
## M2b starts from the M1b fit of the SAME directed cell (m1_spec(); 02_fit_model.R asserts it). The suffix _cdir takes the place of _conv in
## the suffix order (the two values are exclusive), e.g. fmid_lag3_all_pr_M1b_cdir, fmid_lag3_all_pr_M2b_K1_cdir, mid_lag3_all_pr_M1c_cdir_aprime.
DIRECTED_MODELS <- c("M1a", "M1b", "M1c", "M1cr", "M2b")
## ---- ideal-point shifters and the rho regression under the directed comparison (PI-approved 2026-11-02) ------------------------------
## mu_shift = traits (suffix _mst, right after _cdir): X = X_MU_TRAITS (assemble_cell.R) = the initiator's M_lvl, R_lvl, cowmaj_tm1, dem_tm1 at t-1,
## exactly the W columns (same values, same centring: the levels within-year centred, the 0/1 indicators grand-centred), P = 4, through the
## existing X / beta machinery (beta ~ normal(0, prior_beta_sd)); x_ref stays the mean X over the initiating state-years of the cell. Only with
## comparison = directed (parse_spec() refuses it otherwise); models MU_SHIFT_MODELS. A directed M2b with the switch starts from the M1b fit of the
## same directed cell WITH the switch (m1_spec() keeps it; 02_fit_model.R asserts the X columns).
## M1bz (model code; stan/m1_latent_rhoreg.stan): the joint M1b model with a state-year specific stage-2 correlation rho_n = tanh(g0 + Zr[n] g),
## Zr = the same four columns (N1 x Q, Q = 4, data). Only with comparison = directed, with or without mu_shift = traits. Not in MODELS (main grid).
RHOREG_MODELS   <- c("M1bz")
## ---- the AME model with initiation types (PI-approved 2026-11-04): M2c = M2b (m2_ame.stan, joint b variant) plus the M1c role structure (mu_g =
## mu0 + delta[g - 1], sigma_g[g], one shared rho), stan/m2_context.stan. Defined ONLY with --comparison=directed and --roles=aprime (cell_name()
## refuses anything else; no mu_shift, no prior switch); K = 1 | 2; every chain starts from the matching M1c fit (same sample, lag, comparison, roles,
## joiner_set, tag; 02_fit_model.R refuses to run without it). Family M2 (sampler settings, _K<K> in the name), not in MODELS (the main grid).
AME_CONTEXT_MODELS <- c("M2c")
## the role codings built on the aprime file (2026-11-04): those with a joiner group (joiner_set), the ones M2c takes; bc / bcj need entry = later
JOINER_ROLES <- c("aprime", "bc", "bcj")
ENTRY_ROLES  <- c("bc", "bcj")
## ---- the joiner choice-set restriction (PI-approved 2026-11-04; switch joiner_set, suffix _jset right after _aprime) ---------------------------
## joiner_set = opposing: an event of aprime group 2 (joiner) is scored against the states on the side of its dispute OPPOSITE to the joiner's own side
## (COW sides after the direction rule, as in aprime_roles()) that take part in the dispute in the event's onset year, intersected with the state-year's
## stage-2 set; every other event keeps the full stage-2 set; stage 1 / IV and the weights are unchanged. Needs roles = aprime (parse_spec() refuses it
## otherwise), so only the role models (ROLE_MODELS) can take it.
MU_SHIFT_MODELS <- c(DIRECTED_MODELS, RHOREG_MODELS)
## ---- run 0 and the lag-3 MID-history choice set (2026-11-10; README section 13) -------------------------------------------------------------------
## reach = capdist (run 0a, suffix _reach right after _cdir / _mst): Z = Z_DIRECTED + Z_REACH (assemble_cell.R): log_capdist_sq = (centred log_capdist)^2 and
## log_capdist_x_init_M = (centred log_capdist) x the initiator's M_lvl (the stage-1 level, as in W), both then centred over the IV rows. Contiguity stays.
## Only with comparison = directed, model M1c (REACH_MODELS), variant all, no mu_shift, no subsample.
## subsample (runs 0b-0d, suffixes _sprt / _sbil / _smaj / _smin right after the role / joiner-set suffixes): only with comparison = directed, variant all,
## pairs = all, context = all, ww2 = include, joiner_set = all, reach = none.
##   prtarget  (0b): the stage-2 events whose chosen target is NOT politically relevant (the alternatives' polrel column: contiguous 1..5 or a major in the
##                   dyad) are removed; the choice set stays the FULL set; stage 1 is unchanged (y1 not recomputed). Model M1c.
##   bilateral (0c): the stage-2 events of bcj role bilateral (aprime alone originator, roles_aprime.rds) only; stage 1 unchanged. One role, so the pooled
##                   model M1b (M1c stops on an empty group); needs --entry=later or --events=dymid (the bcj coding's data) and no --roles.
##   majinit / mininit (0d): BOTH stages restricted to the state-years whose initiator is / is not a COW major power at t-1 (the panel's cowmaj_tm1, the
##                   stage-1 W column). cowmaj_tm1 (constant in the subset) leaves W and tgt_major_x_init_major (= tgt_major among majors, = 0 among
##                   minors) leaves Z. Model M1c.
## After an event filter (prtarget, bilateral) the D8 weights are rescaled within the (initiator, dispute, year) groups the filter split, exactly as under
## pairs = originators | joiners (d8_rescale_split_groups()).
## SUPERSEDED the same day (audit of d60da78, PI decision 2026-11-10): the event subsamples do NOT rescale the D8 weights. The polrel cell does not rescale
## either (events that leave through the choice set are logged, not rescaled), so 0b carries exactly the polrel cell's events, roles and weights and differs
## from it in the choice set only.
## variant = mid3 (the lag-3 MID-history choice set; takes the variant slot of the cell name like polrel): alternatives j of initiator i in year t with a
## MID onset between i and j (either direction; every row of the build's events.rds, all hostility levels, escalations included) in [t - MID3_WINDOW, t - 1],
## OR dy_mid_tm1 == 1. Year t is never used. Built at assembly from alts_lag<h>_all.rds (no file of its own). Only with comparison = directed, models
## MID3_MODELS, pairs = all, context = all, ww2 = include, joiner_set = all, subsample = none. REFUSED since the audit of 2026-11-10 (MID3_REFUSAL below).
REACH_MODELS      <- c("M1c")
## ---- run 1 (2026-11-10; README section 14): the MID-propensity opportunity-weight SENSITIVITY ------------------------------------------------------
## opp_weight = 0 | 0.5 | 1 (lambda): U_ij = spatial + b_j + Z_ij' gamma + lambda x ln q_ij, the same term inside the stage-1 IV (reference-group utility as
## before); q_ij = the opportunity weight of 08_opportunity_weights.R (output/oppw_<foundation>_<opp_outcome>_lag3.rds). Fixed data-side offset, no new
## parameter. Model M1c only (OPP_MODELS), through the copy stan/m1_context_offset.stan (m1_context.stan and every existing fit are untouched); lambda = 0
## is the offset program with a zero offset: the Stan data of the cell without the switch plus opp_off = 0, and the same log density.
## Only with comparison = directed, roles = bcj, variant = all, and pairs / context / ww2 / joiner_set / mu_shift / reach / subsample at their defaults
## (WWII stays in). opp_draw = 1..5 needs opp_weight = 1; opp_outcome = host3 needs an opp_weight. Suffixes _opw0 / _opw05 / _opw1, then _oh3, then _d<d>,
## after the data-variant suffixes (_ent / _dy) and before the prior switches and the tag, e.g. fmid_lag3_all_pr_M1c_cdir_bcj_dy_opw1_d3_tlong.
OPP_MODELS        <- c("M1c")
OPP_LAMBDA        <- c("0" = 0, "0.5" = 0.5, "1" = 1)
## spec_opp_weight(spec) / spec_opp_outcome(spec) / spec_opp_draw(spec): the values of the three run-1 switches ("none" / "mid" / "0" for a spec saved
## before 2026-11-10, which has none of the fields); spec_opp_on(spec): TRUE when an offset is requested (lambda = 0 included)
spec_opp_weight  <- function(spec) if (is.null(spec$opp_weight)) "none" else spec$opp_weight
spec_opp_outcome <- function(spec) if (is.null(spec$opp_outcome)) "mid" else spec$opp_outcome
spec_opp_draw    <- function(spec) if (is.null(spec$opp_draw)) "0" else spec$opp_draw
spec_opp_on      <- function(spec) !identical(spec_opp_weight(spec), "none")
SUBSAMPLE_MODELS  <- list(prtarget = "M1c", bilateral = "M1b", majinit = "M1c", mininit = "M1c")
INIT_SUBSAMPLES   <- c("majinit", "mininit")          # the subsamples that restrict the stage-1 population as well
EVENT_SUBSAMPLES  <- c("prtarget", "bilateral")       # the subsamples that filter stage-2 events only (stage 1 unchanged)
MID3_MODELS       <- c("M1c", "M2c")
MID3_WINDOW       <- 3L                               # years of MID history before t: onsets in [t - 3, t - 1]
## REFUSED since the audit of 2026-11-10 (P1-1; PI decision the same day): --variant=mid3 stops in parse_spec() and cell_name(), so no script can run a mid3
## cell and no mid3 cell name exists. The construction (mid3_history(), the VARIANT = "mid3" block of assemble_cell()) stays in the code, reachable only by
## calling the assembly directly, and is exercised by unit tests alone.
MID3_REFUSAL <- paste("--variant=mid3 is refused (audit P1-1, PI decision 2026-11-10): with a MID-history window of t-3..t-1 the chosen target must be kept in the set, and a chosen",
                      "target without history is then the ONLY non-history member of its set (60% of the usable fatal pairs events, 53% of the jointly fatal ones), so its stage-2",
                      "term is not a valid conditional probability; widening the window to include year t would define the set by the outcome instead. Use --variant=all or",
                      "--variant=polrel. The construction code is kept for reference (README section 13.2).")
## spec_reach(spec) / spec_subsample(spec): the values of the two run-0 switches ("none" for a spec saved before 2026-11-10, which has neither field)
spec_reach     <- function(spec) if (is.null(spec$reach)) "none" else spec$reach
spec_subsample <- function(spec) if (is.null(spec$subsample)) "none" else spec$subsample
## spec_alts_variant(spec): the variant of the alternatives FILE a spec reads (alts_lag<h>_<this>.rds): mid3 is cut from the all-states file at assembly
spec_alts_variant <- function(spec) { v <- spec_variant(spec); if (identical(v, "mid3")) "all" else v }
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
  ## the conventional comparison (2026-10-31) follows the control extensions it replaces (it is refused with any of them) and precedes _w1
  comparison = c(conventional = "_conv", directed = "_cdir"),   # directed (2026-11-01): same position, exclusive with _conv
  ## the ideal-point shifters (2026-11-02) follow the comparison switch they need: e.g. fmid_lag3_all_pr_M1b_cdir_mst, mid_lag3_all_pr_M1c_cdir_mst_aprime
  mu_shift  = c(traits = "_mst"),
  ## run 0a (2026-11-10) follows the comparison / mu_shift suffixes: e.g. fmid_lag3_all_pr_M1c_cdir_reach_bcj_ent
  reach     = c(capdist = "_reach"),
  ## the role coding (2026-10-31) follows the comparison switch: e.g. fmid_lag3_all_pr_M1c_conv_aprime
  roles     = c(aprime = "_aprime", bc = "_bc", bcj = "_bcj"),   # bc / bcj (2026-11-04): need entry = later
  ## the joiner choice-set restriction (2026-11-04) follows the role coding it needs: e.g. fmid_lag3_all_pr_M1c_cdir_aprime_jset
  joiner_set = c(opposing = "_jset"),
  ## runs 0b-0d (2026-11-10) follow the role / joiner-set suffixes: e.g. fmid_lag3_all_pr_M1c_cdir_bcj_sprt_ent, fmid_lag3_all_pr_M1b_cdir_sbil_dy
  subsample = c(prtarget = "_sprt", bilateral = "_sbil", majinit = "_smaj", mininit = "_smin"),
  weights   = c(none = "_w1"),
  changes   = c(winsor = "_winsor"),
  direction = c(sideA = "_sideA"),
  entry     = c(later = "_ent"),                     # joiner entry direction (2026-11-04), a data variant like direction
  events    = c(dymid = "_dy"),                      # event foundation (2026-11-04), a data variant like entry: e.g. fmid_lag3_all_pr_M2c_K2_cdir_bcj_jset_dy
  ## run 1 (2026-11-10) follows the data variants: the offset weight, then the auxiliary outcome, then the auxiliary draw, e.g. fmid_lag3_all_pr_M1c_cdir_bcj_dy_opw1_d3
  opp_weight = c("0" = "_opw0", "0.5" = "_opw05", "1" = "_opw1"),
  opp_outcome = c(host3 = "_oh3"),
  opp_draw  = c("1" = "_d1", "2" = "_d2", "3" = "_d3", "4" = "_d4", "5" = "_d5"),
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
ALL_SAMPLES  <- c("mid", "fmid", "fmidcow", "fmidub", "war", "fmidany")
## the samples each event foundation's build carries (2026-11-04): the pairs build has the COW / MIC / GML fatality codings (fmidcow), the dymid build the
## dyadic ones (fmid = both states suffered deaths in the dyad's fighting or a war dyad; fmidany = the dyadic dispute was fatal or a war dyad; fmidub =
## fmid with missing side fatalities counted fatal); a sample outside the cell's foundation is refused (sample_in_events())
EVENTS_SAMPLES <- list(pairs = c("mid", "fmid", "fmidcow", "fmidub", "war"), dymid = c("mid", "fmid", "fmidany", "fmidub", "war"))
## spec_events(spec): the event foundation of a spec ("pairs" for a spec saved before 2026-11-04, which has no events field)
spec_events <- function(spec) if (is.null(spec$events)) "pairs" else spec$events
## sample_in_events(sample, spec): stop unless the sample exists in the build of the spec's event foundation
sample_in_events <- function(sample, spec) {
  ev <- spec_events(spec)
  if (!sample %in% EVENTS_SAMPLES[[ev]]) stop(sprintf("--sample=%s does not exist under --events=%s (its samples: %s)", sample, ev, paste(EVENTS_SAMPLES[[ev]], collapse = ", ")))
  invisible(TRUE)
}
LAGS         <- c(1L, 3L, 5L, 10L)
MODELS       <- c("M1m", "M1r", "M1a", "M1b", "M2m", "M2r", "M2a", "M2b")
## dual typology stage (b), 2026-10-27: M1c = the M1b joint model with group-specific ideal-point shifts (delta_lead, delta_join) and
## group-specific tolerance scales (sigma_g, shared correlation rho) over the three stage-2 event groups bilateral / coalition leader /
## joiner, plus one extra stage-1 control (multi_mid_jan). Its own Stan program (stan/m1_context.stan). NOT part of MODELS (the main
## grid's 0-95 index mapping and 04_compare.R's eight-model groups stay as they are); 02_fit_model.R and 03_extract.R accept it.
## M1cr (PI request 2026-10-31): M1c with a ROLE-SPECIFIC stage-2 correlation rho_g[G] (stan/m1_context_rho.stan); same data, same groups
CONTEXT_MODELS <- c("M1c", "M1cr")
ALL_MODELS     <- c(MODELS, CONTEXT_MODELS, RHOREG_MODELS, AME_CONTEXT_MODELS)   # RHOREG_MODELS (M1bz, 2026-11-02), AME_CONTEXT_MODELS (M2c, 2026-11-04) defined above
## ROLE_MODELS (2026-11-04): every model that takes the stage-2 groups (ev_grp, G), the role switch and the joiner choice-set restriction
ROLE_MODELS    <- c(CONTEXT_MODELS, AME_CONTEXT_MODELS)
## stan_program(model, spec): the Stan source / build name of a model code (m1_latent_additive | m2_ame | m1_context), and, since 2026-10-30,
## m1_latent_prior for an M1 m / r / a / b model whose spec has a prior switch off its default (the copy of m1_latent_additive.stan with the
## hyperparameters as data). spec = NULL (or a spec saved before the switches existed, which lacks the fields) means the main priors.
stan_program <- function(model, spec = NULL) {
  stopifnot(model %in% ALL_MODELS)
  if (model %in% RHOREG_MODELS) return("m1_latent_rhoreg")      # M1bz (2026-11-02): the rho-regression program, no prior-switch copy
  if (model %in% AME_CONTEXT_MODELS) return("m2_context")        # M2c (2026-11-04): AME with initiation types
  if (model == "M1cr") return("m1_context_rho")
  if (model %in% OPP_MODELS && !is.null(spec) && spec_opp_on(spec)) return("m1_context_offset")   # run 1 (2026-11-10): M1c with the opportunity-weight offset (lambda = 0 included)
  if (model %in% CONTEXT_MODELS) return("m1_context")
  if (substr(model, 1, 2) == "M2") return("m2_ame")
  if (spec_prior_is_main(spec)) "m1_latent_additive" else "m1_latent_prior"
}
STAN_PROGRAMS  <- c("m1_latent_additive", "m2_ame", "m1_context", "m1_latent_prior", "m1_context_rho", "m1_latent_rhoreg", "m2_context", "m1_context_offset")   # m1_context_offset: run 1 (2026-11-10)
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
VARIANT      <- "all"                                # choice set of the main specification: every other system member (polrel: the `variant` switch, 2026-11-04)
## spec_variant(spec): the choice-set variant of a spec ("all" for a spec saved before 2026-11-04, which has no variant field)
spec_variant <- function(spec) if (is.null(spec$variant)) VARIANT else spec$variant
POLREL_MODELS <- c("M1c", "M1cr", "M2c")             # the models the polrel variant is specified for (PI request 2026-11-04)
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
ARCHIVED_OPTIONS <- c("sampled_m", "centrality", "conflict_history", "dy_peace_years",
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
  ## the conventional comparison (2026-10-31) fixes W, Z and X itself (checked FIRST, so a conflicting switch is reported as a conflict): the switches of CONVENTIONAL_FIXED are refused on the command line in
  ## ANY value (the default included), since passing one claims a control set the spec does not take; their spec fields stay at the defaults
  if (spec$comparison == "conventional") {
    for (nm in CONVENTIONAL_FIXED) if (length(grep(paste0("^--", nm, "="), args)))
      stop(sprintf("--%s=%s conflicts with --comparison=conventional, which fixes W, Z and X itself (W_CONVENTIONAL / Z_CONVENTIONAL, no X: distance, majors, contiguity and alliance are part of the spec); leave --%s out", nm, spec[[nm]], nm))
    stopifnot(all(vapply(CONVENTIONAL_FIXED, function(nm) identical(spec[[nm]], SPEC_DEFAULTS[[nm]]), logical(1))))
  }
  ## the directed comparison (2026-11-01) fixes the same switches as the conventional one, refused the same way (any value, the default included)
  if (spec$comparison == "directed") {
    for (nm in CONVENTIONAL_FIXED) if (length(grep(paste0("^--", nm, "="), args)))
      stop(sprintf("--%s=%s conflicts with --comparison=directed, which fixes W, Z and X itself (W_CONVENTIONAL / Z_DIRECTED, no X: distance, majors, contiguity and alliance are part of the spec); leave --%s out", nm, spec[[nm]], nm))
    stopifnot(all(vapply(CONVENTIONAL_FIXED, function(nm) identical(spec[[nm]], SPEC_DEFAULTS[[nm]]), logical(1))))
  }
  ## the ideal-point shifters (2026-11-02) are defined under the directed comparison only (X_MU_TRAITS are its W columns)
  ## the joiner choice-set restriction (2026-11-04) is defined on the aprime joiner group only
  if (spec$joiner_set != SPEC_DEFAULTS$joiner_set && !spec$roles %in% JOINER_ROLES)
    stop(sprintf("--joiner_set=%s needs --roles=aprime, bc or bcj (the restriction applies to the joiner events; got roles=%s)", spec$joiner_set, spec$roles))
  ## the bc / bcj codings (2026-11-04) have no late-target group, so they need the entry = later events (keep still has late-target events) or the dymid
  ## events, whose direction rule has the entry rule built in
  if (spec$roles %in% ENTRY_ROLES && spec$entry != "later" && spec$events != "dymid") stop(sprintf("--roles=%s needs --entry=later or --events=dymid (under --entry=keep the late-target events have no group)", spec$roles))
  ## the interacting-dyads foundation (2026-11-04): built on the D1 direction with the entry rule as its rule (a), so --direction=sideA and --entry=later are
  ## refused with it (the second would be redundant: the dymid build applies the entry rule itself); roles = aprime is refused (its group 3, target joined
  ## late, cannot occur under the entry rule, so the four-group coding would carry an empty group; use bcj or bc)
  if (spec$events == "dymid") {
    if (spec$direction != SPEC_DEFAULTS$direction) stop("--events=dymid is built on the cowwar direction only (not with --direction=sideA)")
    if (spec$entry != SPEC_DEFAULTS$entry) stop("--events=dymid has the entry rule built in (its direction rule (a)); --entry=later with it is redundant and refused; leave --entry out")
    if (spec$roles == "aprime") stop("--roles=aprime is refused under --events=dymid (no target-joined-late event exists under the entry rule, so group 3 would be empty); use --roles=bcj (main) or --roles=bc")
  }
  if (spec$entry != SPEC_DEFAULTS$entry && spec$direction != SPEC_DEFAULTS$direction) stop("--entry=later is built on the cowwar direction only (not with --direction=sideA)")
  if (spec$mu_shift != SPEC_DEFAULTS$mu_shift && spec$comparison != "directed")
    stop(sprintf("--mu_shift=%s needs --comparison=directed (the shifters are the directed comparison's W columns M_lvl, R_lvl, cowmaj_tm1, dem_tm1; got comparison=%s)", spec$mu_shift, spec$comparison))
  ## contiguity = none is only meaningful when the distance replaces it (PI request 2026-10-28): a Z without contiguity AND without distance is refused
  if (spec$contiguity != "land" && spec$distance != "capital") stop(sprintf("--contiguity=%s needs --distance=capital (distance replaces contiguity; got distance=%s)", spec$contiguity, spec$distance))
  ## run 0 and the mid3 choice set (2026-11-10): invalid combinations are refused, never ignored
  if (spec$reach != SPEC_DEFAULTS$reach) {
    if (spec$comparison != "directed") stop(sprintf("--reach=%s needs --comparison=directed (the terms extend the directed Z; got comparison=%s)", spec$reach, spec$comparison))
    if (spec$mu_shift != SPEC_DEFAULTS$mu_shift) stop(sprintf("--reach=%s with --mu_shift=%s is not specified (run 0a is the directed M1c cell without ideal-point shifters)", spec$reach, spec$mu_shift))
    if (spec$variant == "mid3") stop(MID3_REFUSAL)
    if (spec$variant != "all") stop(sprintf("--reach=%s is defined on the full choice set (got --variant=%s)", spec$reach, spec$variant))
    if (spec$subsample != SPEC_DEFAULTS$subsample) stop(sprintf("--reach=%s with --subsample=%s is not specified (the run-0 cells flip one switch each)", spec$reach, spec$subsample))
  }
  if (spec$subsample != SPEC_DEFAULTS$subsample) {
    if (spec$comparison != "directed") stop(sprintf("--subsample=%s needs --comparison=directed (run 0 is defined on the directed specification; got comparison=%s)", spec$subsample, spec$comparison))
    if (spec$variant == "mid3") stop(MID3_REFUSAL)
    if (spec$variant != "all") stop(sprintf("--subsample=%s is defined on the full choice set (got --variant=%s)", spec$subsample, spec$variant))
    for (nm in c("pairs", "context", "ww2", "joiner_set", "mu_shift")) if (spec[[nm]] != SPEC_DEFAULTS[[nm]])
      stop(sprintf("--subsample=%s with --%s=%s is not specified (one event / sample filter at a time)", spec$subsample, nm, spec[[nm]]))
    if (spec$subsample == "bilateral") {
      if (spec$roles != SPEC_DEFAULTS$roles) stop(sprintf("--subsample=bilateral keeps one role (bcj bilateral), so it takes the pooled model M1b and no --roles (got roles=%s)", spec$roles))
      if (spec$entry != "later" && spec$events != "dymid") stop("--subsample=bilateral needs --entry=later or --events=dymid (the bcj role coding is defined on those events only)")
    }
  }
  if (spec$variant == "mid3") stop(MID3_REFUSAL)        # audit P1-1 (2026-11-10): refused on every command line
  ## run 1 (2026-11-10): the opportunity-weight offset is specified for the directed bcj cells on the full choice set only; WWII stays in
  if (spec$opp_weight != SPEC_DEFAULTS$opp_weight) {
    if (spec$comparison != "directed") stop(sprintf("--opp_weight=%s needs --comparison=directed (run 1 is defined on the directed specification; got comparison=%s)", spec$opp_weight, spec$comparison))
    if (spec$roles != "bcj") stop(sprintf("--opp_weight=%s needs --roles=bcj (run 1 is the directed bcj M1c cell; got roles=%s)", spec$opp_weight, spec$roles))
    if (spec$variant != "all") stop(sprintf("--opp_weight=%s is defined on the full choice set (got --variant=%s)", spec$opp_weight, spec$variant))
    if (spec$direction != SPEC_DEFAULTS$direction) stop(sprintf("--opp_weight=%s is not specified with --direction=%s", spec$opp_weight, spec$direction))
    for (nm in c("pairs", "context", "ww2", "joiner_set", "mu_shift", "reach", "subsample", "weights")) if (spec[[nm]] != SPEC_DEFAULTS[[nm]])
      stop(sprintf("--opp_weight=%s with --%s=%s is not specified (run 1 changes the offset only; the Second World War stays in)", spec$opp_weight, nm, spec[[nm]]))
  }
  if (spec$opp_outcome != SPEC_DEFAULTS$opp_outcome && spec$opp_weight == SPEC_DEFAULTS$opp_weight) stop(sprintf("--opp_outcome=%s needs --opp_weight=0 | 0.5 | 1 (it names the auxiliary outcome behind the offset)", spec$opp_outcome))
  if (spec$opp_draw != SPEC_DEFAULTS$opp_draw && spec$opp_weight != "1") stop(sprintf("--opp_draw=%s needs --opp_weight=1 (the auxiliary draws are run at lambda = 1 only; got opp_weight=%s)", spec$opp_draw, spec$opp_weight))
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
## entry = later (2026-11-04): v2/data_entry (a spec saved before has no entry field: keep)
## events = dymid (2026-11-04): v2/data_dymid (a spec saved before has no events field: pairs)
spec_data_dir <- function(spec, root) file.path(root, if (spec$direction == "sideA") "data_sideA" else if (!is.null(spec$entry) && identical(spec$entry, "later")) "data_entry" else
                                                  if (identical(spec_events(spec), "dymid")) "data_dymid" else "data")
## spec_string(spec): one line for the logs
## mu_shift (2026-11-02) and joiner_set (2026-11-04) are left out at their defaults, so the spec string of every existing cell is unchanged
spec_string <- function(spec) {
  if (!is.null(spec$mu_shift) && identical(spec$mu_shift, SPEC_DEFAULTS$mu_shift)) spec <- spec[setdiff(names(spec), "mu_shift")]
  if (!is.null(spec$joiner_set) && identical(spec$joiner_set, SPEC_DEFAULTS$joiner_set)) spec <- spec[setdiff(names(spec), "joiner_set")]
  if (!is.null(spec$entry) && identical(spec$entry, SPEC_DEFAULTS$entry)) spec <- spec[setdiff(names(spec), "entry")]
  if (!is.null(spec$events) && identical(spec$events, SPEC_DEFAULTS$events)) spec <- spec[setdiff(names(spec), "events")]
  if (!is.null(spec$variant) && identical(spec$variant, SPEC_DEFAULTS$variant)) spec <- spec[setdiff(names(spec), "variant")]
  for (nm in c("reach", "subsample", "opp_weight", "opp_outcome", "opp_draw")) if (!is.null(spec[[nm]]) && identical(spec[[nm]], SPEC_DEFAULTS[[nm]])) spec <- spec[setdiff(names(spec), nm)]   # 2026-11-10: left out at their defaults
  paste(sprintf("%s=%s", names(spec), unlist(spec)), collapse = " ")
}
## spec_is_conventional(spec): TRUE under --comparison=conventional (a spec saved before 2026-10-31 has no field: not conventional)
spec_is_conventional <- function(spec) !is.null(spec$comparison) && identical(spec$comparison, "conventional")
## spec_is_directed(spec): TRUE under --comparison=directed (2026-11-01; a spec saved before has no such value: not directed)
spec_is_directed <- function(spec) !is.null(spec$comparison) && identical(spec$comparison, "directed")
## perm_part_spec_string(spec): the spec string written into the permutation refit parts of 05_robustness.R and checked by its combiner
## (audit R1 of b1f0988: the 7c4645e form, without the prior switches and the tag, so parts written by that commit still combine). The
## comparison switch (2026-10-31) is dropped at its default for the same reason and KEPT when non-default: it changes W, Z and X, so a
## conventional cell's parts must never combine with another cell's (the driver refuses the switch anyway, see 05_robustness.R).
perm_part_spec_string <- function(spec) {
  drop <- c(PRIOR_SWITCHES, "tag", if (!spec_is_conventional(spec) && !spec_is_directed(spec)) "comparison", if (is.null(spec$roles) || identical(spec$roles, SPEC_DEFAULTS$roles)) "roles",
            if (is.null(spec$mu_shift) || identical(spec$mu_shift, SPEC_DEFAULTS$mu_shift)) "mu_shift",
            if (is.null(spec$joiner_set) || identical(spec$joiner_set, SPEC_DEFAULTS$joiner_set)) "joiner_set",
            if (is.null(spec$entry) || identical(spec$entry, SPEC_DEFAULTS$entry)) "entry",
            if (is.null(spec$events) || identical(spec$events, SPEC_DEFAULTS$events)) "events",
            if (is.null(spec$variant) || identical(spec$variant, SPEC_DEFAULTS$variant)) "variant",
            if (is.null(spec$reach) || identical(spec$reach, SPEC_DEFAULTS$reach)) "reach", if (is.null(spec$subsample) || identical(spec$subsample, SPEC_DEFAULTS$subsample)) "subsample",
            unlist(lapply(c("opp_weight", "opp_outcome", "opp_draw"), function(nm) if (is.null(spec[[nm]]) || identical(spec[[nm]], SPEC_DEFAULTS[[nm]])) nm)))   # joiner_set (2026-11-04): likewise   # mu_shift (2026-11-02): dropped at its default, kept when on
  spec_string(spec[setdiff(names(spec), drop)])
}
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
  if (identical(spec_variant(spec), "mid3")) stop(paste("cell_name:", MID3_REFUSAL))      # audit P1-1 (2026-11-10): no mid3 cell name exists
  sample_in_events(sample, spec)                    # 2026-11-04: fmidany exists under events = dymid only, fmidcow under pairs only
  if (model_family(model) == "M1" && !identical(as.integer(spec$ame_K), SPEC_DEFAULTS$ame_K))
    stop(sprintf("cell_name: --ame_K=%s with the M1 model %s (K applies to the M2 family only and would not appear in the cell name)", spec$ame_K, model))
  ## the conventional comparison (2026-10-31) is specified for M1a, M1b and M1c only
  if (spec_is_conventional(spec) && !model %in% CONVENTIONAL_MODELS)
    stop(sprintf("cell_name: --comparison=conventional with the model %s (the conventional comparison is specified for %s only)", model, paste(CONVENTIONAL_MODELS, collapse = ", ")))
  ## the directed comparison (2026-11-01): the conventional models plus the AME model M2b; the rho-regression model M1bz (2026-11-02) as well
  if (spec_is_directed(spec) && !model %in% c(DIRECTED_MODELS, RHOREG_MODELS, AME_CONTEXT_MODELS))
    stop(sprintf("cell_name: --comparison=directed with the model %s (the directed comparison is specified for %s only; and the rho regression %s)", model, paste(DIRECTED_MODELS, collapse = ", "), paste(RHOREG_MODELS, collapse = ", ")))
  ## M1bz (2026-11-02) is defined under the directed comparison only
  if (model %in% RHOREG_MODELS && !spec_is_directed(spec))
    stop(sprintf("cell_name: the model %s needs --comparison=directed (the rho regression is specified under the directed comparison only; got comparison=%s)", model, if (is.null(spec$comparison)) "none" else spec$comparison))
  ## mu_shift = traits (2026-11-02): directed only (parse_spec() refuses it otherwise) and the models MU_SHIFT_MODELS
  if (!is.null(spec$mu_shift) && !identical(spec$mu_shift, SPEC_DEFAULTS$mu_shift)) {
    if (!spec_is_directed(spec)) stop(sprintf("cell_name: --mu_shift=%s needs --comparison=directed", spec$mu_shift))
    if (!model %in% MU_SHIFT_MODELS) stop(sprintf("cell_name: --mu_shift=%s with the model %s (specified for %s only)", spec$mu_shift, model, paste(MU_SHIFT_MODELS, collapse = ", ")))
  }
  ## the role coding (2026-10-31) applies to the stage-2 groups of the context models only (no other model sees the groups); M2c (2026-11-04) as well
  if (!is.null(spec$roles) && !identical(spec$roles, SPEC_DEFAULTS$roles) && !model %in% ROLE_MODELS)
    stop(sprintf("cell_name: --roles=%s with the model %s (the role coding applies to %s only)", spec$roles, model, paste(ROLE_MODELS, collapse = ", ")))
  ## the joiner choice-set restriction (2026-11-04): roles = aprime only (parse_spec() refuses it otherwise), hence the role models only
  if (!is.null(spec$joiner_set) && !identical(spec$joiner_set, SPEC_DEFAULTS$joiner_set)) {
    if (!isTRUE(spec$roles %in% JOINER_ROLES)) stop(sprintf("cell_name: --joiner_set=%s needs --roles=aprime, bc or bcj", spec$joiner_set))
    if (!model %in% ROLE_MODELS) stop(sprintf("cell_name: --joiner_set=%s with the model %s (the restriction applies to %s only)", spec$joiner_set, model, paste(ROLE_MODELS, collapse = ", ")))
  }
  ## M2c (2026-11-04): the AME model with initiation types is specified under --comparison=directed with --roles=aprime only
  if (model %in% AME_CONTEXT_MODELS && (!spec_is_directed(spec) || !isTRUE(spec$roles %in% JOINER_ROLES)))
    stop(sprintf("cell_name: the model %s needs --comparison=directed and --roles=aprime, bc or bcj (got comparison=%s, roles=%s)", model, if (is.null(spec$comparison)) "none" else spec$comparison, if (is.null(spec$roles)) "stageb" else spec$roles))
  if (!spec_prior_is_main(spec)) {
    off <- PRIOR_SWITCHES[vapply(PRIOR_SWITCHES, function(nm) !is.null(spec[[nm]]) && !identical(spec[[nm]], SPEC_DEFAULTS[[nm]]), logical(1))]
    if (stan_program(model) != "m1_latent_additive")
      stop(sprintf("cell_name: --%s=%s with the model %s (the prior switches apply to the M1 m / r / a / b programs only; %s has no such hyperparameter)", off[1], spec[[off[1]]], model, stan_program(model)))
    if ("rho_prior" %in% off && model_joint(model) != 1L)
      stop(sprintf("cell_name: --rho_prior=%s with the model %s, which has no Sigma correlation rho (the rho prior applies to the joint model M1b only)", spec$rho_prior, model))
  }
  ## the politically relevant choice set (2026-11-04): the role / AME models only; its name takes the variant slot (<sample>_lag<h>_polrel_pr_...)
  ## run 0 and mid3 (2026-11-10): the models each switch is specified for
  if (spec_reach(spec) != "none" && !model %in% REACH_MODELS)
    stop(sprintf("cell_name: --reach=%s with the model %s (run 0a is specified for %s only)", spec_reach(spec), model, paste(REACH_MODELS, collapse = ", ")))
  if (spec_subsample(spec) != "none" && !identical(model, SUBSAMPLE_MODELS[[spec_subsample(spec)]]))
    stop(sprintf("cell_name: --subsample=%s with the model %s (specified for %s only%s)", spec_subsample(spec), model, SUBSAMPLE_MODELS[[spec_subsample(spec)]],
                 if (spec_subsample(spec) == "bilateral") ": one role, so the pooled model; M1c stops on an empty group" else ""))
  if (spec_opp_on(spec) && !model %in% OPP_MODELS)
    stop(sprintf("cell_name: --opp_weight=%s with the model %s (run 1 is specified for %s only; no AME model takes the offset)", spec_opp_weight(spec), model, paste(OPP_MODELS, collapse = ", ")))
  if (identical(spec_variant(spec), "polrel") && !model %in% POLREL_MODELS)
    stop(sprintf("cell_name: --variant=%s with the model %s (the politically relevant choice set is specified for %s only)", spec_variant(spec), model, paste(POLREL_MODELS, collapse = ", ")))
  sprintf("%s_lag%d_%s_%s_%s%s%s", sample, as.integer(lag), spec_variant(spec), CENTRALITY, model,
          if (grepl("^M2", model)) sprintf("_K%d", spec$ame_K) else "", spec_suffix(spec))
}
## compare_group_name(sample, lag, spec): the <group> of 04_compare.R, <sample>_lag<h>_all_pr[_K<K>]<suffix>; "_K<K>" only for a non-default
## ame_K, so the documented main-grid names (mid_lag5_all_pr, fmid_lag3_all_pr_noww2, ...) are unchanged and a K=3 comparison never
## overwrites the K=2 one (review 2026-10-26, P2)
compare_group_name <- function(sample, lag, spec) {
  stopifnot(sample %in% ALL_SAMPLES, lag %in% LAGS)
  sprintf("%s_lag%d_%s_%s%s%s", sample, as.integer(lag), spec_variant(spec), CENTRALITY,
          if (!identical(as.integer(spec$ame_K), SPEC_DEFAULTS$ame_K)) sprintf("_K%d", spec$ame_K) else "", spec_suffix(spec))
}
## measurement_cell_name(sample, lag, family, spec): the <cell> of the 05_robustness.R perm_error outputs (measurement_<cell>_*.csv),
## <sample>_lag<h>_all_pr_<M1|M2>[_K<K>]<suffix>. The M2 name ALWAYS carries its K, as cell_name() does for M2 fits (review 2026-10-26,
## P2: K=2 and K=3 wrote the same mid_lag5_all_pr_M2 files); the M1 name never does.
measurement_cell_name <- function(sample, lag, family, spec) {
  stopifnot(sample %in% ALL_SAMPLES, lag %in% LAGS, family %in% c("m1", "m2"))
  sprintf("%s_lag%d_%s_%s_%s%s%s", sample, as.integer(lag), spec_variant(spec), CENTRALITY, toupper(family),
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
## M1cr (2026-10-31): the third character decides (M1cr is a "c" model; a regex on the LAST character would read its "r" as relational-only);
## M1bz (2026-11-02) reads as "b" (both axes, joint): it is the M1b model with rho regressed on the initiator's traits
model_D_mode <- function(model) { stopifnot(model %in% ALL_MODELS); switch(substr(model, 3, 3), m = 1L, r = 2L, a = 3L, b = 3L, c = 3L) }
model_joint  <- function(model) { stopifnot(model %in% ALL_MODELS); as.integer(substr(model, 3, 3) %in% c("b", "c")) }
## sampler_settings(model, args): the family defaults, each overridable from the command line
sampler_settings <- function(model, args = character(0)) {
  s <- SAMPLER[[model_family(model)]]
  list(chains = as.integer(get_arg(args, "chains", s$chains)), iter = as.integer(get_arg(args, "iter", s$iter)),
       warmup = as.integer(get_arg(args, "warmup", s$warmup)), adapt_delta = get_arg(args, "adapt_delta", s$adapt_delta),
       max_treedepth = as.integer(get_arg(args, "max_treedepth", s$max_treedepth)))
}
