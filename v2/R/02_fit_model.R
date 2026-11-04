## ============================================================================
## 02_fit_model.R -- assemble the Stan data for ONE cell (sample, lag, model) and fit it
##
## Usage:
##   Rscript v2/R/02_fit_model.R --sample=mid --lag=3 --model=M1b [switches] [sampler] [--pilot=FALSE] [--dry_run=FALSE]
##          [--threads_per_chain=1]
##
##   sample : mid (all MIDs) | fmid (fatal MIDs, main coding) | war (pair-level war); appendix: fmidcow | fmidub
##   lag    : 1 | 3 | 5 | 10   change-window length h
##   model  : M1m | M1r | M1a | M1b | M2m | M2r | M2a | M2b | M1c | M1cr (M1cr, 2026-10-31: M1c with a role-specific stage-2 correlation rho_g, stan/m1_context_rho.stan)
##            | M2c (PI-approved 2026-11-04; stan/m2_context.stan): AME with initiation types = M2b plus the M1c role structure (mu_g = mu0 + delta[g - 1],
##            sigma_g[g], shared rho); --comparison=directed and --roles=aprime only, --ame_K=1|2, with or without --joiner_set; every chain starts from
##            the matching M1c fit (same sample, lag, switches, tag; K dropped), and the run REFUSES to sample when that fit is missing
##            m = material only, r = relational only, a = both independent, b = both joint; M2 = M1 + AME factors (ame_K)
##            M1bz (PI-approved 2026-11-02; stan/m1_latent_rhoreg.stan; --comparison=directed only, with or without --mu_shift=traits): M1b with
##            a state-year specific stage-2 correlation rho_n = tanh(g0 + Zr[n] g), Zr = the initiator's M_lvl, R_lvl, cowmaj_tm1, dem_tm1 as in W
##            (rhoreg_Zr(), assemble_cell.R; Q = 4), used for that state-year's stage-2 softmax and its IV; g0 ~ N(0, 1), g ~ N(0, 0.5)
##            M1c (dual typology stage b, 2026-10-27; stan/m1_context.stan): M1b with group-specific ideal-point shifts delta_lead /
##            delta_join and tolerance scales sigma_g over the stage-2 groups bilateral / coalition leader / joiner (shared rho), and
##            multi_mid_jan appended to W. Cell <sample>_lag<h>_all_pr_M1c<suffix>; M1 sampler settings; not in the main grid
##   switches (spec.R; every non-default value adds a cell-name suffix, so variants never overwrite main outputs):
##     --pairs=all|originators|joiners  --direction=cowwar|sideA  --ww2=include|drop  --weights=dispute|none
##     --controls=minimal|full          --changes=raw|winsor      --ame_K=2
##     --context=all|bilateral|coalition (STAGE-2 events of one coalition context; stage 1 unchanged)  --coal_ties=inclusive|strict
##     --coal_drop=none|top3 (context=coalition only: without the 3 disputes contributing the most coalition events)
##     --major=none|joint (one_major / both_major in Z, cowmaj_tm1 in W)  --interact=none|prefixed (four change products in Z)
##     --distance=none|capital (2026-10-28: log_capdist in Z, alternatives without a capital distance dropped from both sets; suffix _dist)
##     --contiguity=land|none (2026-10-28: contig_land dropped from Z; needs --distance=capital, distance replaces contiguity; suffix _nocontig)
##     --comparison=none|conventional (PI-approved 2026-10-31, suffix _conv; M1a / M1b / M1c / M1cr only): W, Z, X fixed to the conventional comparison
##            (W_CONVENTIONAL / Z_CONVENTIONAL of assemble_cell.R, no X: P = 0); carries the distance, the major indicators, contiguity, the
##            defense alliance itself (no peace years, PI 2026-10-31), so --controls, --distance, --contiguity, --major, --interact and --changes
##            are refused with it in any value (parse_spec()); multi_mid_jan is in W for all four models
##     --comparison=directed (PI-approved 2026-11-01, suffix _cdir; M1a / M1b / M1c / M1cr / M2b only): the conventional comparison with the
##            directed stage-2 Z (Z_DIRECTED of assemble_cell.R: tgt_dem, tgt_dem_x_init_dem, tgt_major, tgt_major_x_init_major in place of
##            joint_dem, one_major, both_major; Pz = 9); same W, no X (P = 0, also in m2_ame.stan since 2026-11-01), same refusals; M2b starts
##            from the M1b fit of the same directed cell (asserted: never a _conv or default M1b fit)
##     --mu_shift=none|traits (PI-approved 2026-11-02, suffix _mst after _cdir; --comparison=directed only; M1a / M1b / M1c / M1cr / M2b / M1bz):
##            X = M_lvl, R_lvl, cowmaj_tm1, dem_tm1 (X_MU_TRAITS, the W columns centred as in W; P = 4), mu = mu0 (+ delta[g]) + X beta + a_i, x_ref the
##            mean X over the initiating state-years; a directed M2b with the switch starts from the directed M1b fit WITH the switch (asserted)
##     --roles=stageb|aprime (PI decision 2026-10-31, suffix _aprime; M1c / M1cr only): the stage-2 groups; aprime = 0 alone originator, 1 co-originator,
##            2 joiner, 3 target joined late (G = 4; roles_aprime.rds of the build); stageb (default) = bilateral / coalition leader / joiner
##     --joiner_set=all|opposing (PI-approved 2026-11-04, suffix _jset after _aprime; --roles=aprime only, M1c / M1cr / M2c): opposing scores every joiner
##            event (aprime group 2) only against the states on the opposite side of its dispute in the onset year (joiner_sets_aprime.rds of the build),
##            intersected with the state-year's stage-2 set; other events, stage 1 / IV and the weights unchanged; all = byte-identical to before
##     --rho_prior=beta25|flat  --sigma_prior=ln1|ln05 (2026-10-30, prior sensitivity of the jointness result; M1 m / r / a / b only, the rho
##            prior with M1b only): a non-default value fits the SAME model with stan/m1_latent_prior.stan, the copy of the main program that
##            takes the Beta shape on (rho + 1) / 2 (2.5 | 1 = uniform on (-1, 1)) and the lognormal sd of sigma_x (1 | 0.5) as data
##            (prior_hyper(), spec.R); suffixes _rhoflat / _sigln05 after every data switch; refused with M2 and M1c (cell_name())
##   tag    : --tag=<alnum> (2026-10-30) appends _t<tag> LAST to the cell name, so the fit, summary, standata and log files of this run never
##            overwrite an existing fit (e.g. --iter=6000 --warmup=3000 --tag=long for a longer-chain refit while jobs read the untagged fit).
##            Default "" = every name unchanged. The tag reaches every name built from the spec: an M2 run with a tag looks for the
##            TAGGED M1 fit (and falls back loudly to random inits when there is none; it never takes the untagged file)
##   sampler: --chains --iter --warmup --adapt_delta --max_treedepth (defaults per family in spec.R, SAMPLER: M2 runs 4000 / 2000
##            with adapt_delta 0.95 since the PI decision of 2026-10-25)
##   init_from: (M2 only) the M1 fit whose posterior means start every chain (jittered; AME factors near 0; spec.R INIT_*). Default:
##            the matching M1 fit of the same cell in v2/output (M2x from M1x, same sample, lag and switches). When no such FILE exists
##            the run falls back to Stan's random inits and says so in the log, the fit object and the summary. A fit file that exists
##            but fails the checks (readable, production, current data, same cell, same D_mode / joint) is an ERROR, whether auto-located
##            or named with --init_from (audit fix 1, 2026-10-26).
##   pilot  : TRUE = 1 chain x 60 iterations (timing only); saved as fit_<cell>_pilot.rds, which every production tool refuses
##   dry_run: TRUE = assemble and validate the Stan data, write output/standata_<cell>.rds, do NOT sample
##   threads_per_chain: within-chain parallelism via reduce_sum over state-years (1 = serial, exact); Slurm cpus = chains x threads
##
## The choice set (every other system member, full sets), the PageRank centrality and the control-set constants are fixed in
## spec.R (the former --variant / --sampled_m / --centrality / --conflict_history / --dy_peace_years options are
## archived; see README "Archived options"; --distance returned as the switch none | capital on 2026-10-28, --contiguity as land | none).
##
## Every input file is checked against <data dir>/manifest.csv (nrow + sha256) before use (audit F3/F8/F9). Output:
##   v2/output/fit_<cell>.rds    fit, Stan data, cell description, manifest rows, sampling configuration, compiled-model hashes
##                               (written atomically: fit_<cell>.rds.tmp then renamed; 2026-10-25 C2); fit_<cell>_STALE.rds
##                               when the data manifest changed while sampling (downstream tools need --allow_stale). Written ONLY
##                               when the stanfit holds every requested chain with all its post-warmup draws (assert_fit_chains();
##                               review 2026-10-26 P1: rstan drops a failed parallel chain and returns the rest with mode 0)
##   v2/output/summary_<cell>.csv posterior summaries of the structural parameters, labelled by axis and covariate (written AFTER
##                               the fit file)
## Cell name: <sample>_lag<h>_all_pr_<model>[_K<K>][_orig|_join][_bi|_co][_strict][_drop3][_noww2][_fullctrl][_dist][_nocontig][_major][_ix][_conv|_cdir][_mst][_aprime][_jset][_w1][_winsor][_sideA][_rhoflat][_sigln05][_t<tag>][_pilot]
##            (model M1c: <sample>_lag<h>_all_pr_M1c<suffix>; its fit records args$stan_program = "m1_context" and args$context_model = TRUE)
## Seed 123 everywhere (audit F11).
## ============================================================================
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tidyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
R_DIR <- dirname(this_file)
source(file.path(R_DIR, "spec.R")); source(file.path(R_DIR, "helpers_data.R")); source(file.path(R_DIR, "helpers_model.R"))
source(file.path(R_DIR, "assemble_cell.R"))
Sys.setenv(V2_ROOT = normalizePath(file.path(R_DIR, "..")))

## ---- arguments -------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
SAMPLE <- get_arg(args, "sample", "mid"); LAG <- as.integer(get_arg(args, "lag", 5)); MODEL <- get_arg(args, "model", "M1b")
stopifnot(SAMPLE %in% ALL_SAMPLES, LAG %in% LAGS, MODEL %in% ALL_MODELS)
spec <- parse_spec(args)
smp <- sampler_settings(MODEL, args)
CHAINS <- smp$chains; ITER <- smp$iter; WARMUP <- smp$warmup; ADAPT_DELTA <- smp$adapt_delta; MAX_TD <- smp$max_treedepth
PILOT <- get_arg(args, "pilot", FALSE); DRY_RUN <- get_arg(args, "dry_run", FALSE)
THREADS <- as.integer(get_arg(args, "threads_per_chain", 1))
INIT_FROM <- get_arg(args, "init_from", "")
refuse_archived_options(args, "02_fit_model.R")          # archived options are refused instead of being ignored silently (ARCHIVED_OPTIONS, spec.R)
stopifnot(CHAINS >= 1, ITER > WARMUP, WARMUP >= 0, THREADS >= 1, is.finite(ADAPT_DELTA), ADAPT_DELTA > 0, ADAPT_DELTA < 1, MAX_TD >= 1)

D_MODE <- model_D_mode(MODEL); joint <- model_joint(MODEL); family <- stan_program(MODEL, spec)   # m1_latent_additive | m2_ame | m1_context (M1c) | m1_latent_prior (a prior switch off its default)
PRIOR_PROGRAM <- family == "m1_latent_prior"; prior <- prior_hyper(spec)   # the hyperparameters the prior program takes as data (2.5 / 1 = the main program's literals)
CONTEXT_MODEL <- MODEL %in% ROLE_MODELS                  # M1c / M1cr, and M2c (2026-11-04): the models that take the stage-2 groups
AME_FAMILY <- family %in% c("m2_ame", "m2_context")      # the AME programs (K, anchors, U / V; M2c 2026-11-04)
JOINER_SET <- if (is.null(spec$joiner_set)) "all" else spec$joiner_set
RHOREG <- MODEL %in% RHOREG_MODELS                       # M1bz (2026-11-02): the rho-regression program m1_latent_rhoreg
MU_SHIFT <- if (is.null(spec$mu_shift)) "none" else spec$mu_shift   # 2026-11-02
if (RHOREG) stopifnot(family == "m1_latent_rhoreg", D_MODE == 3L, joint == 1L, spec_is_directed(spec))
if (CONTEXT_MODEL && (spec$controls != "minimal" || spec$context != "all")) stop("M1c is defined on the minimal control set with every stage-2 event (--controls=minimal, --context=all); it models the contexts itself")
## merge 2026-10-27 (stage b x major / interact): the two features were developed independently and M1c under the major or interact switches
## has not been specified by the PI; refused here until it is, rather than silently combining the W / Z extensions with the M1c program
if (CONTEXT_MODEL && (spec$distance != "none" || spec$contiguity != "land")) stop("M1c with --distance=capital or --contiguity=none is not specified (the stage-b model is defined on the main Z); run M1c with --distance=none --contiguity=land")
if (CONTEXT_MODEL && (spec$major != "none" || spec$interact != "none")) stop("M1c with --major=joint or --interact=prefixed is not specified (the stage-b model and the 2026-10-27 stage-2 switches were developed independently); run M1c with --major=none --interact=none")
K <- spec$ame_K
cell <- cell_name(SAMPLE, LAG, MODEL, spec)
if (PILOT) cell <- paste0(cell, "_pilot")          # C1: pilots never share a filename with production fits
log_open(paste0("02_fit_", cell))
logf("cell = %s; spec: %s (%s)", cell, spec_string(spec), if (spec_is_main(spec)) "MAIN specification" else paste0("variant, suffix ", spec_suffix(spec)))
if (nzchar(spec_tag(spec))) logf("run tag %s (2026-10-30): every output of this run carries the suffix %s%s (fit_%s.rds, summary_%s.csv, standata_%s.rds, this log); the untagged files are not touched", spec_tag(spec), TAG_PREFIX, spec_tag(spec), cell, cell, cell)
logf("Stan program: %s%s", family, if (PRIOR_PROGRAM) sprintf(" (prior sensitivity, 2026-10-30: rho_prior = %s -> Beta(%.2f, %.2f) on (rho + 1) / 2; sigma_prior = %s -> sigma_x ~ lognormal(0, %.2f); the main program has 2.5 and 1)", spec$rho_prior, prior$prior_rho_shape, prior$prior_rho_shape, spec$sigma_prior, prior$prior_sigma_x_sd) else " (the main priors)")
logf("sampler (%s family): chains %d, iter %d, warmup %d, adapt_delta %.3f, max_treedepth %d", model_family(MODEL), CHAINS, ITER, WARMUP, ADAPT_DELTA, MAX_TD)
set.seed(V2_SEED)
rstan_options(auto_write = FALSE)
rstan_options(threads_per_chain = THREADS)                 # needs the model compiled with -DSTAN_THREADS (rstan 2.32 default)
options(mc.cores = min(CHAINS, max(1L, n_cores() %/% THREADS)))        # n_cores(): detectCores() with a fallback when it returns NA
logf("threads_per_chain = %d; mc.cores = %d", THREADS, getOption("mc.cores"))
data_dir <- spec_data_dir(spec, v2_root()); out_dir <- file.path(v2_root(), "output"); dir.create(out_dir, showWarnings = FALSE)
if (!dir.exists(data_dir)) stop("data directory ", data_dir, " does not exist (for --direction=sideA run `Rscript v2/R/01_build_data.R --direction=sideA` first)")
logf("direction = %s (data directory %s); event weights = %s", spec$direction, basename(data_dir), spec$weights)
manifest_path <- file.path(data_dir, "manifest.csv")
## ---- provenance at LOAD time (audit 2026-10-24: a fit could be stamped with data that changed while it sampled) -------
## The manifest digests and the git commit are read BEFORE the data are assembled and sampling starts; the manifest is
## read again just before the fit is saved and the run stops if any digest changed in between. Slurm ids are recorded.
stopifnot(file.exists(manifest_path))
manifest_at_load <- read.csv(manifest_path, stringsAsFactors = FALSE)
stopifnot(all(c("file", "nrow", "sha256") %in% names(manifest_at_load)), nrow(manifest_at_load) > 0, !anyNA(manifest_at_load$sha256))
git_commit_at_load <- tryCatch(system2("git", c("-C", v2_root(), "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[1], error = function(e) NA_character_)
if (is.null(git_commit_at_load) || length(git_commit_at_load) != 1 || is.na(git_commit_at_load) || !grepl("^[0-9a-f]{40}$", git_commit_at_load)) {
  logf("git rev-parse HEAD unavailable at load (recorded as NA)"); git_commit_at_load <- NA_character_
}
git_dirty_at_load <- if (is.na(git_commit_at_load)) NA else
  tryCatch(length(system2("git", c("-C", v2_root(), "status", "--porcelain", "--", "R", "stan"), stdout = TRUE, stderr = FALSE)) > 0, error = function(e) NA)
slurm_job_id <- Sys.getenv("SLURM_JOB_ID", unset = NA_character_); slurm_task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA_character_)
if (!is.na(slurm_job_id) && slurm_job_id == "") slurm_job_id <- NA_character_
if (!is.na(slurm_task_id) && slurm_task_id == "") slurm_task_id <- NA_character_
logf("provenance at load: git commit %s (uncommitted changes under v2/R or v2/stan: %s); manifest built %s with %d files; SLURM_JOB_ID %s; SLURM_ARRAY_TASK_ID %s",
     git_commit_at_load, git_dirty_at_load, if ("built" %in% names(manifest_at_load)) manifest_at_load$built[1] else "NA", nrow(manifest_at_load), slurm_job_id, slurm_task_id)
cell_files <- c(sprintf("stage1_panel_lag%d.rds", LAG), sprintf("alts_lag%d_%s.rds", LAG, VARIANT), "events.rds",
                if (identical(spec$roles, "aprime")) "roles_aprime.rds",   # roles = aprime (2026-10-31) reads the build's role file as well
                if (JOINER_SET == "opposing") "joiner_sets_aprime.rds")   # joiner_set = opposing (2026-11-04) reads the joiner sets as well
stopifnot(all(cell_files %in% manifest_at_load$file))
manifest_cell <- manifest_at_load[manifest_at_load$file %in% cell_files, ]

## the manifest SNAPSHOT read above is passed in (review P2-5): every input file's nrow and sha256 are checked against it
## at load, so a rebuild between the snapshot and the load stops the run instead of being read silently from disk
A <- assemble_cell_spec(SAMPLE, LAG, spec, data_dir = data_dir, manifest = manifest_at_load, model = MODEL)   # model = M1c appends multi_mid_jan to W
sy <- A$sy; alts <- A$alts; ev <- A$ev; W_cols <- A$W_cols; X_cols <- A$X_cols; Z_cols <- A$Z_cols
state_codes <- A$state_codes; x_ref <- A$x_ref; n_full <- A$n_full; sampled <- A$sampled

## ---- Stan data ---------------------------------------------------------------------------------------------
dat <- build_stan_data(sy, alts, ev, W_cols, X_cols, Z_cols, joint = joint, sampled = sampled, n_full = n_full,
                       x_ref = x_ref, K = if (AME_FAMILY) K else NULL, allow_orphan_initiations = TRUE,
                       D_mode = D_MODE, grainsize = choose_grainsize(nrow(sy), THREADS))
if (CONTEXT_MODEL) {
  ## stage b: the stage-2 group of every usable event goes to Stan as ev_grp (0..G-1; 0 = the reference group) with the group count G
  ## (data since 2026-10-31; G = 3 with the stage-b coding: 0 bilateral, 1 coalition leader, 2 joiner); the existing programs never see either
  ## (they are not added to their data, so their Stan data stay identical to the main grid's)
  gt <- A$group_table
  G_GRP <- nrow(gt)
  stopifnot(G_GRP >= 2L, identical(gt$grp, 0:(G_GRP - 1L)), identical(gt$group, stage2_group_labels(spec$roles)), identical(A$roles, spec$roles),
            "grp" %in% names(ev), length(ev$grp) == dat$E, all(ev$grp %in% 0:(G_GRP - 1L)), "multi_mid_jan" %in% W_cols, identical(as.integer(dat$ev_sy), as.integer(ev$sy_row)))
  dat$G <- G_GRP
  dat$ev_grp <- as.integer(ev$grp)
  if (any(gt$n_usable == 0L)) stop(sprintf("%s needs events in every stage-2 group; usable events per group: %s", MODEL, paste(sprintf("%s %d", gt$group, gt$n_usable), collapse = ", ")))
  logf("%s stage-2 groups (ev_grp, G = %d; usable events, weight sum, distinct initiating states; sample events before the choice-set join): %s", MODEL, G_GRP,
       paste(sprintf("%d %s %d (%.1f, %d; %d)", gt$grp, gt$group, gt$n_usable, gt$w_usable, gt$n_states, gt$n_sample), collapse = "; "))
  ## event-level allowed stage-2 rows (2026-11-04): the whole stage-2 set of the state-year for every event under joiner_set = all (the programs then
  ## compute exactly what they computed before), the joiner events' opposing-side rows under joiner_set = opposing
  stopifnot(identical(A$joiner_set, JOINER_SET), is.null(A$ev_sets) == (JOINER_SET == "all"))
  dat <- event_alt_sets(dat, A$ev_sets)
  s2_n <- tabulate(findInterval(which(dat$in_s2 == 1L), dat$cs_start), nbins = dat$N1)      # stage-2 set size per state-year (cs_start increasing)
  stopifnot(identical(dat$cs_start, sort(dat$cs_start)), all(dat$ev_alt_n <= s2_n[dat$ev_sy]))
  is_restr <- dat$ev_alt_n < s2_n[dat$ev_sy]; n_restr <- sum(is_restr)
  if (JOINER_SET == "opposing") {
    js <- A$jset_stats
    stopifnot(!is.null(js), all(dat$ev_grp[is_restr] == 2L), n_restr == js$n_joiner_events - js$n_unrestricted)
    logf("joiner_set = opposing (2026-11-04): %d joiner events (aprime group 2); restricted set size min %d / median %.1f / max %d; size 1: %d (log-prob 0 in stage 2, kept; weights %.1f); %d events actually restricted in Stan (ev_restr = 1); alternatives removed %d; n_ev_alt = %d",
         js$n_joiner_events, js$size_min, js$size_median, js$size_max, js$n_singleton, js$w_singleton, n_restr, js$n_removed, dat$n_ev_alt)
  } else {
    stopifnot(n_restr == 0L)
    logf("joiner_set = all: every event is scored against its state-year's whole stage-2 set (n_ev_alt = %d)", dat$n_ev_alt)
  }
}
if (RHOREG) {
  ## M1bz (2026-11-02): the traits of the rho regression, N1 x Q, exactly the W columns M_lvl, R_lvl, cowmaj_tm1, dem_tm1 (rhoreg_Zr() asserts them)
  Zr <- rhoreg_Zr(A)
  stopifnot(nrow(Zr) == dat$N1, ncol(Zr) == 4L, identical(unname(Zr), unname(dat$W[, match(X_MU_TRAITS, W_cols), drop = FALSE])), all(is.finite(Zr)))
  dat$Q <- ncol(Zr); dat$Zr <- Zr
  logf("M1bz rho regression (2026-11-02): rho_n = tanh(g0 + Zr[n] g), Zr = %s (Q = %d, the W columns as they enter W); column means over the %d state-years %s; g0 ~ N(0, 1), g ~ N(0, 0.5)",
       paste(X_MU_TRAITS, collapse = ","), dat$Q, dat$N1, paste(sprintf("%.4f", colMeans(Zr)), collapse = " / "))
}
CONVENTIONAL <- spec_is_conventional(spec)
DIRECTED <- spec_is_directed(spec)                           # 2026-11-01
if (CONVENTIONAL) {
  ## the conventional comparison (PI-approved 2026-10-31): the three matrices are exactly the spec's lists, no X (P = 0; the Stan programs take
  ## P = 0 since this date: beta is 0 x D, X is N1 x 0, x_ref is empty), no change column anywhere, multi_mid_jan in W for every model
  stopifnot(isTRUE(A$conventional), A$distance == "capital", A$major == "joint", identical(W_cols, W_CONVENTIONAL), identical(Z_cols, Z_CONVENTIONAL), length(X_cols) == 0L,
            dat$P == 0L, ncol(dat$X) == 0L, nrow(dat$X) == dat$N1, length(dat$x_ref) == 0L, dat$P1 == length(W_CONVENTIONAL), dat$Pz == length(Z_CONVENTIONAL),
            !any(CHANGE_COLUMNS %in% c(W_cols, X_cols, Z_cols)), !any(PEACE_COLUMNS %in% c(W_cols, Z_cols)), !anyNA(dat$W), !anyNA(dat$Z), "multi_mid_jan" %in% W_cols, MODEL %in% CONVENTIONAL_MODELS)
  logf("conventional comparison (2026-10-31, suffix _conv): W = %s (P1 = %d); X = none (P = 0, mu = mu0 + a_i%s); Z = %s (Pz = %d); distance = capital and major = joint are part of the spec (not of the switches)",
       paste(W_cols, collapse = ","), dat$P1, if (CONTEXT_MODEL) " + delta[ev_grp]" else "", paste(Z_cols, collapse = ","), dat$Pz)
} else if (DIRECTED) {
  ## the directed comparison (PI-approved 2026-11-01): the conventional W, no X (P = 0, in the M1 programs and m2_ame.stan), Z exactly Z_DIRECTED;
  ## under mu_shift = traits (2026-11-02) X is exactly X_MU_TRAITS (P = 4), the same values as those W columns
  X_EXPECTED <- if (MU_SHIFT == "traits") X_MU_TRAITS else character(0)
  stopifnot(MU_SHIFT %in% c("none", "traits"), MU_SHIFT == "none" || MODEL %in% MU_SHIFT_MODELS)
  if (MU_SHIFT == "traits") stopifnot(identical(X_cols, c("M_lvl", "R_lvl", "cowmaj_tm1", "dem_tm1")), dat$P == 4L, ncol(dat$X) == 4L, length(dat$x_ref) == 4L, all(is.finite(dat$x_ref)),
                                      identical(unname(dat$X), unname(dat$W[, match(X_MU_TRAITS, W_cols), drop = FALSE])))
  stopifnot(isTRUE(A$directed), !isTRUE(A$conventional), A$distance == "capital", A$major == "joint", identical(W_cols, W_CONVENTIONAL), identical(Z_cols, Z_DIRECTED), identical(X_cols, X_EXPECTED),
            dat$P == length(X_EXPECTED), ncol(dat$X) == length(X_EXPECTED), nrow(dat$X) == dat$N1, length(dat$x_ref) == length(X_EXPECTED), dat$P1 == length(W_CONVENTIONAL), dat$Pz == length(Z_DIRECTED), dat$Pz == 9L,
            identical(colnames(dat$Z), Z_DIRECTED) || is.null(colnames(dat$Z)),
            !any(CHANGE_COLUMNS %in% c(W_cols, X_cols, Z_cols)), !any(PEACE_COLUMNS %in% c(W_cols, Z_cols)), !any(c("joint_dem", Z_MAJOR) %in% Z_cols), !anyNA(dat$W), !anyNA(dat$Z),
            "multi_mid_jan" %in% W_cols, MODEL %in% c(DIRECTED_MODELS, RHOREG_MODELS, AME_CONTEXT_MODELS))
  logf("directed comparison (2026-11-01, suffix _cdir): W = %s (P1 = %d); X = %s (P = %d, mu = mu0 + %sa_i%s); Z = %s (Pz = %d); distance = capital and major = joint are part of the spec (not of the switches)",
       paste(W_cols, collapse = ","), dat$P1, if (dat$P > 0L) paste(X_cols, collapse = ",") else "none", dat$P, if (dat$P > 0L) "X beta + " else "", if (CONTEXT_MODEL) " + delta[ev_grp]" else "", paste(Z_cols, collapse = ","), dat$Pz)
  if (MU_SHIFT == "traits") logf("mu_shift = traits (2026-11-02, suffix _mst): X = %s (P = %d; the W columns, centred as in W); x_ref (mean X over the initiating state-years) = %s",
                                 paste(X_cols, collapse = ","), dat$P, paste(sprintf("%s=%.4f", X_cols, dat$x_ref), collapse = ", "))
} else stopifnot(!isTRUE(A$conventional), !isTRUE(A$directed), dat$P >= 1L, MU_SHIFT == "none", !RHOREG)
if (PRIOR_PROGRAM) {
  ## prior sensitivity (2026-10-30): the hyperparameters go to Stan as data; the main program's data stay exactly as before (no field added)
  stopifnot(is.finite(prior$prior_rho_shape), prior$prior_rho_shape > 0, is.finite(prior$prior_sigma_x_sd), prior$prior_sigma_x_sd > 0)
  dat$prior_rho_shape <- prior$prior_rho_shape; dat$prior_sigma_x_sd <- prior$prior_sigma_x_sd
}
logf("reduce_sum grainsize = %d (N1 = %d)", dat$grainsize, dat$N1)
logf("D_mode = %d (%s)", D_MODE, c("material only", "relational only", "both dimensions")[D_MODE])
if (AME_FAMILY) logf("AME anchors (most initiations first): ccodes %s", paste(state_codes[dat$anchor], collapse = ", "))
logf("Stan data: S=%d N1=%d N2=%d E=%d P1=%d P=%d Pz=%d joint=%d sampled=%d D_mode=%d; sum(y1)=%d; rows in the IV set %d, in the stage-2 set %d; event weights sum %.1f (min %.4f, %d events with weight 1); W = %s; X = %s; Z = %s",
     dat$S, dat$N1, dat$N2, dat$E, dat$P1, dat$P, dat$Pz, dat$joint, dat$sampled, dat$D_mode, sum(dat$y1), sum(dat$in_iv), sum(dat$in_s2),
     sum(dat$ev_w), min(dat$ev_w), sum(dat$ev_w == 1), paste(W_cols, collapse = ","), paste(X_cols, collapse = ","), paste(Z_cols, collapse = ","))
logf("D2 levels in W: %s; D3 log system size in W: %s; D1 reversed events in the sample: %d", paste(intersect(W_cols, W_LEVELS), collapse = ","),
     "log_n_system_tm1" %in% W_cols, A$n_flipped)
if (spec$pairs != "all" && spec$weights == "dispute") logf("pairs = %s D8 rescale (PI decision 2026-10-26): %d of %d (initiator, dispute, year) groups split by the filter, %d rescaled to sum 1, %d events rescaled; sum of weights over the kept events %.1f -> %.1f",
                                                            spec$pairs, A$pairs_weights$n_groups_split, A$pairs_weights$n_groups, A$pairs_weights$n_groups_rescaled, A$pairs_weights$n_events_rescaled,
                                                            A$pairs_weights$w_before, A$pairs_weights$w_after)
if (spec$context != "all") logf("context = %s (%s%s): %d stage-2 events removed, %d D8 groups split through collapsed keys; stage 1 unchanged, %d state-years with y1 = 1 and no stage-2 event",
                                spec$context, spec$coal_ties, if (spec$coal_drop == "top3") sprintf(", without disputes %s", paste(A$context_top3$dispnum, collapse = ", ")) else "",
                                A$n_context_dropped, A$n_context_split_groups, A$n_orphan)
if (A$major == "joint") logf("major = joint (2026-10-27%s): cowmaj_tm1 in W; one_major / both_major in Z; alternatives by major count (0 / 1 / 2) %.4f / %.4f / %.4f, stage-2 events %.4f / %.4f / %.4f (%d / %d / %d of %d); choice sets with within-set variation: IV sets one_major %d, both_major %d of %d; event-bearing stage-2 sets one_major %d, both_major %d of %d",
                                if (CONVENTIONAL) ", part of the conventional comparison" else if (DIRECTED) ", built for the identity checks of the directed comparison; one_major / both_major NOT in Z" else "", A$major_stats$share_alts[["zero"]], A$major_stats$share_alts[["one"]], A$major_stats$share_alts[["both"]],
                                A$major_stats$share_events[["zero"]], A$major_stats$share_events[["one"]], A$major_stats$share_events[["both"]],
                                A$major_stats$n_events[["zero"]], A$major_stats$n_events[["one"]], A$major_stats$n_events[["both"]], dat$E,
                                A$major_stats$var_iv_sets[["one_major"]], A$major_stats$var_iv_sets[["both_major"]], A$major_stats$var_iv_sets[["n_sets"]],
                                A$major_stats$var_s2_sets[["one_major"]], A$major_stats$var_s2_sets[["both_major"]], A$major_stats$var_s2_sets[["n_sets"]])
if (A$distance == "capital") logf("distance = capital (2026-10-28%s): log_capdist in Z (centred over the IV rows); %d alternative rows without a capital distance dropped from both sets (%d IV, %d stage 2; %d state-years), %d expected events lost with them (counted and reconciled); IV-row mean %.4f, SD %.4f, cor with contig_land %.3f; mean uncentred log_capdist chosen %.4f vs all stage-2 alternatives %.4f",
                                     if (CONVENTIONAL) ", part of the conventional comparison" else if (DIRECTED) ", part of the directed comparison" else "", A$distance_stats$n_rows_dropped, A$distance_stats$n_rows_dropped_iv, A$distance_stats$n_rows_dropped_s2, A$distance_stats$n_sy_affected, A$distance_stats$n_events_dropped,
                                     A$distance_stats$mean_iv, A$distance_stats$sd_iv, A$distance_stats$cor_contig_iv, A$distance_stats$mean_chosen, A$distance_stats$mean_s2)
if (spec$contiguity == "none") logf("contiguity = none (2026-10-28): contig_land is NOT in Z; the distance replaces it (Z = the --distance=capital Z without its contig_land column; same rows, same events)")
if (spec$interact == "prefixed") logf("interact = prefixed (2026-10-27): Z carries %s (products of the uncentred %s changes, centred over the IV rows)", paste(Z_INTERACT, collapse = ","), spec$changes)
if (DRY_RUN) {
  saveRDS(list(dat = dat, cell = cell, state_codes = state_codes, W_cols = W_cols, X_cols = X_cols, Z_cols = Z_cols,
               manifest = manifest_cell, git_commit = git_commit_at_load, spec = spec, winsor = A$winsor, model = MODEL,
               full_controls = spec$controls == "full", event_weights = spec$weights, direction = spec$direction, group_table = A$group_table,
               stan_program = family, prior = c(spec[PRIOR_SWITCHES], prior), tag = spec_tag(spec), conventional = CONVENTIONAL, directed = DIRECTED,
               mu_shift = MU_SHIFT, Zr_cols = if (RHOREG) X_MU_TRAITS else NULL, joiner_set = JOINER_SET, jset_stats = A$jset_stats),
          file.path(out_dir, paste0("standata_", cell, ".rds")))
  logf("DRY RUN: Stan data assembled and saved; no sampling"); quit(save = "no", status = 0)
}
if (PILOT) { ITER <- 60L; WARMUP <- 30L; CHAINS <- 1L; logf("PILOT run: chains=1 iter=60 (timing only); saved with the _pilot suffix") }

## ---- fit ------------------------------------------------------------------------------------------------------------
sm <- readRDS(file.path(v2_root(), "stan", "build", paste0(family, ".rds")))
## compiled-model provenance (2026-10-25 C3): the compiled model_code must be the checked-out .stan source, and the sidecar hash
## written by compile_models.R must match it; both hashes are saved with the fit. (rstan silently ignores unused data, so an old
## binary would otherwise run on new data.)
stan_src <- file.path(v2_root(), "stan", paste0(family, ".stan"))
stan_prov <- assert_compiled_current(sm, stan_src, file.path(v2_root(), "stan", "build", paste0(family, ".sha256")))
logf("compiled %s matches %s (source sha256 %s)", family, basename(stan_src), stan_prov$stan_source_sha256)
## ---- initial values (PI decision 2026-10-25): an M2 fit starts at the posterior means of the matching M1 fit ------------------------
## M1 fits use Stan's random inits. For M2 the matching M1 fit (M2x -> M1x, same cell otherwise) is located in out_dir, or named with
## --init_from; it must be a production fit of the CURRENT data that scores exactly this cell (C1 / C4) with the same D_mode and joint
## flag. m1_inits_for_m2() (helpers_model.R) copies every shared parameter, reconstructs z_re / z_b, jitters each chain and starts the
## AME factors near 0. ONLY "no M1 fit file found" falls back to random inits (loudly, recorded everywhere); an M1 file that exists
## but is unusable (unreadable, pilot, stale, other data, other cell, other D_mode / joint, a failure inside the init construction) is an
## ERROR: the run stops (audit fix 1, 2026-10-26; before, every rejection fell back silently to random inits with exit 0).
init_arg <- "random"
init_info <- list(source = "stan_random", m1_fit = NA_character_, m1_cell = NA_character_, reason = "M1 family: Stan's random inits",
                  jitter_sd = NA_real_, uv_range = NA_real_, n_copied = NA_integer_)
if (family == "m2_ame") {
  explicit <- nzchar(INIT_FROM)
  m1_path <- if (explicit) INIT_FROM else locate_m1_fit(SAMPLE, LAG, MODEL, spec, out_dir)
  if (is.na(m1_path) || !file.exists(m1_path)) {
    init_info$reason <- sprintf("no matching M1 fit (%s) in %s", paste0("fit_", m1_cell_for_m2(SAMPLE, LAG, MODEL, spec), ".rds"), out_dir)
    if (explicit) stop("--init_from=", INIT_FROM, " does not exist")
    logf("WARNING: %s; the M2 fit starts from Stan's RANDOM inits (run the M1 cell first to start from its posterior means)", init_info$reason)
  } else {
    o1 <- tryCatch(readRDS(m1_path), error = function(e) NULL)
    if (is.null(o1) || is.null(o1$fit) || is.null(o1$dat)) stop("the M1 fit ", m1_path, " exists but is unreadable or not a fit object; the M2 fit does NOT fall back to random inits (delete or refit the M1 file)")
    if (!PILOT) assert_production_fit(o1, allow_stale = FALSE, what = "init_from M1 fit")
    assert_fit_current(o1, manifest_path)
    assert_fit_matches_cell(o1, dat, W_cols, X_cols, Z_cols, what = "init_from M1 fit")
    if (o1$dat$D_mode != D_MODE || o1$dat$joint != joint) stop(sprintf("D_mode / joint of the M1 fit (%d / %d) differ from this cell (%d / %d)", o1$dat$D_mode, o1$dat$joint, D_MODE, joint))
    if (is.null(o1$args$model) || model_family(o1$args$model) != "M1") stop("the init fit is not an M1 fit (", o1$args$model, ")")
    ## directed comparison (2026-11-01): the init fit must be the M1 fit of the SAME directed cell (never a _conv or a default-spec M1 fit)
    if (DIRECTED) {
      want <- m1_cell_for_m2(SAMPLE, LAG, MODEL, spec)
      ## mu_shift (2026-11-02): the init fit must carry the same switch and the same X columns (P = 0 without it, X_MU_TRAITS with it)
      o1_mu <- if (is.null(o1$spec$mu_shift)) "none" else o1$spec$mu_shift
      if (!identical(o1$cell, want) || !spec_is_directed(o1$spec) || !isTRUE(o1$args$directed) || !identical(o1$Z_cols, Z_DIRECTED) || !identical(o1_mu, MU_SHIFT) ||
          !identical(o1$X_cols, X_cols) || o1$dat$P != length(X_cols) || !identical(o1$args$model, "M1b"))
        stop(sprintf("--comparison=directed: the init fit %s (cell %s, comparison %s, mu_shift %s) is not the M1 fit of this directed cell (%s)", basename(m1_path), o1$cell, if (is.null(o1$spec$comparison)) "none" else o1$spec$comparison, o1_mu, want))
      logf("directed comparison: the init fit is the M1 fit of the same directed cell (%s)", o1$cell)
    }
    res <- m1_inits_for_m2(o1, dat, CHAINS, seed = V2_SEED, jitter_sd = INIT_JITTER_SD, uv_range = INIT_UV_RANGE)
    stopifnot(is.list(res), length(res) == CHAINS, !is.null(attr(res, "n_copied")))
    init_arg <- res
    init_info <- list(source = "m1_posterior_means", m1_fit = basename(m1_path), m1_cell = o1$cell, reason = if (explicit) "--init_from" else "auto-located in output/",
                      jitter_sd = INIT_JITTER_SD, uv_range = INIT_UV_RANGE, n_copied = attr(res, "n_copied"),
                      m1_git_commit = o1$args$git_commit, m1_saved = o1$args$saved)
    logf("inits: %d chains start at the posterior means of %s (%s; %d parameter blocks copied, jitter sd %.3f, AME factors uniform on [%.2f, %.2f])",
         CHAINS, basename(m1_path), init_info$reason, init_info$n_copied, INIT_JITTER_SD, INIT_UV_RANGE[1], INIT_UV_RANGE[2])
  }
}
## ---- M2c (PI-approved 2026-11-04): every chain starts from the matching M1c fit (same sample, lag, comparison, roles, joiner_set, tag; K dropped:
## m1_cell_for_m2() maps M2c -> M1c). Unlike M2b there is NO random-init fallback: a missing M1c fit stops the run (PI spec). The M1c fit must be a
## production fit of the current data scoring exactly this cell (events, weights, W / X / Z, groups, allowed sets).
if (family == "m2_context") {
  explicit <- nzchar(INIT_FROM)
  want <- m1_cell_for_m2(SAMPLE, LAG, MODEL, spec)
  m1_path <- if (explicit) INIT_FROM else locate_m1_fit(SAMPLE, LAG, MODEL, spec, out_dir)
  if (is.na(m1_path) || !file.exists(m1_path))
    stop(sprintf("M2c needs the matching M1c fit fit_%s.rds in %s (or --init_from=); it is missing, so the run stops (no random-init fallback for M2c)", want, out_dir))
  o1 <- tryCatch(readRDS(m1_path), error = function(e) NULL)
  if (is.null(o1) || is.null(o1$fit) || is.null(o1$dat)) stop("the M1c fit ", m1_path, " exists but is unreadable or not a fit object")
  if (!PILOT) assert_production_fit(o1, allow_stale = FALSE, what = "init_from M1c fit")
  assert_fit_current(o1, manifest_path)
  assert_fit_matches_cell(o1, dat, W_cols, X_cols, Z_cols, what = "init_from M1c fit")
  o1_js <- if (is.null(o1$spec$joiner_set)) "all" else o1$spec$joiner_set
  if (!identical(o1$args$model, "M1c") || !identical(o1$cell, want) || !spec_is_directed(o1$spec) || !identical(o1$spec$roles, "aprime") || !identical(o1_js, JOINER_SET) ||
      !identical(spec_tag(o1$spec), spec_tag(spec)) || o1$dat$D_mode != 3L || o1$dat$joint != 1L || !identical(as.integer(o1$dat$G), as.integer(dat$G)) ||
      !identical(as.integer(o1$dat$ev_grp), as.integer(dat$ev_grp)) || !identical(o1$group_labels, A$group_table$group))
    stop(sprintf("the init fit %s (cell %s, model %s, roles %s, joiner_set %s) is not the M1c fit of this M2c cell (%s)", basename(m1_path), o1$cell, o1$args$model,
                 if (is.null(o1$spec$roles)) "stageb" else o1$spec$roles, o1_js, want))
  ## the allowed sets: identical when the M1c fit carries them (fits saved before 2026-11-04 have none and are joiner_set = all, i.e. the full sets)
  if (!is.null(o1$dat$ev_alt_n)) {
    if (!identical(as.integer(o1$dat$ev_alt_n), as.integer(dat$ev_alt_n)) || !identical(as.integer(o1$dat$ev_alt_loc), as.integer(dat$ev_alt_loc)))
      stop(sprintf("the init fit %s scores other event-level choice sets than this M2c cell", basename(m1_path)))
  } else if (JOINER_SET != "all") stop(sprintf("the init fit %s carries no event-level choice sets but this cell has --joiner_set=%s", basename(m1_path), JOINER_SET))
  res <- m1c_inits_for_m2c(o1, dat, CHAINS, seed = V2_SEED, jitter_sd = INIT_JITTER_SD, uv_range = INIT_UV_RANGE)
  stopifnot(is.list(res), length(res) == CHAINS, !is.null(attr(res, "n_copied")))
  init_arg <- res
  init_info <- list(source = "m1_posterior_means", m1_fit = basename(m1_path), m1_cell = o1$cell, reason = if (explicit) "--init_from" else "auto-located in output/",
                    jitter_sd = INIT_JITTER_SD, uv_range = INIT_UV_RANGE, n_copied = attr(res, "n_copied"), m1_git_commit = o1$args$git_commit, m1_saved = o1$args$saved)
  logf("inits (M2c): %d chains start at the posterior means of the M1c fit %s (%s; %d parameter blocks copied, jitter sd %.3f, AME factors uniform on [%.2f, %.2f])",
       CHAINS, basename(m1_path), init_info$reason, init_info$n_copied, INIT_JITTER_SD, INIT_UV_RANGE[1], INIT_UV_RANGE[2])
}
stopifnot(identical(init_info$source, if (AME_FAMILY && !is.na(init_info$m1_fit)) "m1_posterior_means" else "stan_random"))
if (nzchar(INIT_FROM) && !AME_FAMILY) stop("--init_from applies to the M2 family only")
## test hook (tests/testthat/test_fit_smoke.R; review 2026-10-26, P1): V2_TEST_REJECT_CHAIN=<k> gives chain k an init Stan must reject
## (sigma_x < 0) so the test can check that a run with ONE failed parallel chain saves no production fit and exits non-zero. Never set it
## outside the tests; the log says loudly when it is.
reject_chain <- suppressWarnings(as.integer(Sys.getenv("V2_TEST_REJECT_CHAIN", unset = "")))
if (length(reject_chain) == 1 && !is.na(reject_chain)) {
  stopifnot(reject_chain >= 1L, reject_chain <= CHAINS)
  logf("WARNING: TEST HOOK V2_TEST_REJECT_CHAIN=%d: chain %d gets an init Stan rejects (sigma_x < 0); this run must NOT produce a production fit", reject_chain, reject_chain)
  if (identical(init_arg, "random")) init_arg <- lapply(seq_len(CHAINS), function(ch) list(alpha0 = 0))   # partial named inits: the rest stays random
  init_arg[[reject_chain]]$sigma_x <- array(-1, dim = if (D_MODE == 3L) 2L else 1L)
}
t0 <- Sys.time()
fit <- sampling(sm, data = dat, chains = CHAINS, iter = ITER, warmup = WARMUP, seed = V2_SEED, init = init_arg,
                refresh = max(1L, ITER %/% 20L), control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TD),
                pars = c("z_re", "z_b", "V_raw", "U_free", "L_Om_joint", "L_Om", "re"), include = FALSE)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
if (fit@mode != 0L) stop(sprintf("sampling did not run (stanfit mode %d: rstan reported an error above, e.g. inits rejected); nothing saved", fit@mode))
## every requested chain must be in the stanfit with all its post-warmup draws (review 2026-10-26, P1): with parallel chains rstan drops a
## failed chain and returns the rest with mode 0, which would otherwise be saved as a production fit and release the afterok dependents
assert_fit_chains(fit, chains = CHAINS, iter = ITER, warmup = WARMUP, what = "sampled fit", cell = cell)
logf("sampling finished in %.1f min; %d chains x %d post-warmup draws, as requested", elapsed, fit_chain_summary(fit)$chains, fit_chain_summary(fit)$draws[1])
sp <- get_sampler_params(fit, inc_warmup = FALSE)
logf("divergences: %d; max-treedepth hits: %d", sum(sapply(sp, function(m) sum(m[, "divergent__"]))),
     sum(sapply(sp, function(m) sum(m[, "treedepth__"] >= MAX_TD))))
pars_struct <- c("alpha0", "alpha", "tau", "mu0", if (dat$P > 0L) "beta", if (RHOREG) c("g0", "g"), if (CONTEXT_MODEL) c("delta", "sigma_g", if (family == "m1_context_rho") "rho_g" else "rho_x", "mu0_g") else "sigma_x",
                 "sd_re", "sigma_b", "gamma", "Sigma", "Omega", "R_Om", if (AME_FAMILY) "u_diag")
sm_tab <- summary(fit, pars = pars_struct)$summary
logf("worst Rhat among structural parameters: %.3f; min n_eff: %.0f",
     max(sm_tab[, "Rhat"], na.rm = TRUE), min(sm_tab[, "n_eff"], na.rm = TRUE))
## ---- provenance just before saving: the data files must not have changed while sampling ---------------------------
man_now <- read.csv(manifest_path, stringsAsFactors = FALSE)
chg <- merge(manifest_at_load[, c("file", "sha256")], man_now[, c("file", "sha256")], by = "file", all = TRUE, suffixes = c("_load", "_now"))
chg_bad <- chg$file[is.na(chg$sha256_load) | is.na(chg$sha256_now) | chg$sha256_load != chg$sha256_now]
disk_bad <- cell_files[vapply(cell_files, function(f) file_digest(file.path(data_dir, f)) != manifest_cell$sha256[manifest_cell$file == f], logical(1))]
stale <- NULL
if (length(chg_bad) > 0 || length(disk_bad) > 0) {
  ## C2 (2026-10-25): the fit is NOT discarded; it is saved as fit_<cell>_STALE.rds with the mismatch recorded, and every downstream
  ## tool refuses it unless --allow_stale is given
  stale <- list(reason = "data manifest or cell files changed while sampling", changed_files = union(chg_bad, disk_bad),
                manifest_at_load = manifest_at_load, manifest_at_save = man_now, detected = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  logf("WARNING: manifest.csv / cell files changed while sampling (%s); the fit is saved as STALE", paste(stale$changed_files, collapse = ", "))
}
git_commit_now <- tryCatch(system2("git", c("-C", v2_root(), "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[1], error = function(e) NA_character_)
if (!identical(git_commit_now, git_commit_at_load)) logf("WARNING: git HEAD moved while sampling (%s -> %s); the fit records the commit at load", git_commit_at_load, git_commit_now)
## the fit object records the sampling configuration (C1), the switches, the compiled-model hashes (C3) and, when the manifest changed
## while sampling, the mismatch (C2)
obj <- list(fit = fit, dat = dat, cell = cell, state_codes = state_codes, W_cols = W_cols, X_cols = X_cols, Z_cols = Z_cols,
            sy_keys = sy |> select(ccode, year, state_idx, y1), ev = ev, x_ref = x_ref, elapsed_min = elapsed,
            manifest = manifest_cell, data_dir = basename(data_dir), stale = stale, spec = spec, winsor = A$winsor, init = init_info, group_table = A$group_table,
            prior = c(spec[PRIOR_SWITCHES], prior, list(stan_program = family)),   # the prior switches and the hyperparameters the program saw (2026-10-30)
            param_labels = label_parameters(rownames(sm_tab), W_cols, X_cols, Z_cols, D_MODE, groups = A$group_table$group, Zr_cols = if (RHOREG) X_MU_TRAITS else NULL),
            Zr_cols = if (RHOREG) X_MU_TRAITS else NULL,                     # M1bz (2026-11-02): the traits of rho_n, dat$Zr
            jset_stats = A$jset_stats,                                        # joiner_set = opposing (2026-11-04): joiner events, set sizes, singletons, removed rows
            group_labels = if (CONTEXT_MODEL) A$group_table$group else NULL,   # M1c / M1cr: the label of ev_grp g is group_labels[g + 1] (2026-10-31)
            args = list(sample = SAMPLE, lag = LAG, variant = VARIANT, model = MODEL, centrality = CENTRALITY, K = K, chains = CHAINS,
                        iter = ITER, warmup = WARMUP, pilot = PILOT, sampled_m = SAMPLED_M, conflict_history = CONFLICT_HISTORY,
                        distance = A$distance, distance_stats = A$distance_stats, contiguity = spec$contiguity, dy_peace_years = DY_PEACE_YEARS, drop_ww2 = spec$ww2 == "drop",   # A$distance / A$major: the assembled values (= the switch, or "capital" / "joint" under the conventional comparison)
                        full_controls = spec$controls == "full", event_weights = spec$weights, direction = spec$direction,
                        pairs = spec$pairs, changes = spec$changes, context = spec$context, coal_ties = spec$coal_ties, coal_drop = spec$coal_drop,
                        context_top3 = A$context_top3, major = A$major, major_stats = A$major_stats, interact = spec$interact,
                        comparison = spec$comparison, conventional = CONVENTIONAL, directed = DIRECTED, directed_stats = A$directed_stats,
                        mu_shift = MU_SHIFT, rhoreg = RHOREG,                  # 2026-11-02
                        joiner_set = JOINER_SET,                               # 2026-11-04
                        D_mode = D_MODE, threads_per_chain = THREADS, context_model = CONTEXT_MODEL, stan_program = family,
                        rho_prior = spec$rho_prior, sigma_prior = spec$sigma_prior, prior_rho_shape = prior$prior_rho_shape, prior_sigma_x_sd = prior$prior_sigma_x_sd, tag = spec_tag(spec),
                        init_source = init_info$source, init_from = init_info$m1_fit,
                        grainsize = dat$grainsize, seed = V2_SEED, adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TD,
                        stan_source_sha256 = stan_prov$stan_source_sha256, compiled_model_code_sha256 = stan_prov$compiled_model_code_sha256,
                        git_commit = git_commit_at_load, git_commit_at_save = git_commit_now, git_dirty_at_load = git_dirty_at_load,
                        manifest_built = if ("built" %in% names(manifest_at_load)) manifest_at_load$built[1] else NA_character_,
                        slurm_job_id = slurm_job_id, slurm_array_task_id = slurm_task_id,
                        started = format(t0, "%Y-%m-%d %H:%M:%S"), saved = format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
## atomic save (C2): write to a temporary name, then rename, so a killed job never leaves a truncated fit_<cell>.rds
fit_name <- paste0("fit_", cell, if (!is.null(stale)) "_STALE" else "", ".rds")
fit_path <- file.path(out_dir, fit_name); tmp_path <- paste0(fit_path, ".tmp")
saveRDS(obj, tmp_path)
chk <- tryCatch(readRDS(tmp_path), error = function(e) NULL)
if (is.null(chk) || !identical(chk$cell, cell)) stop("the saved fit could not be read back from ", tmp_path, "; not renamed")
if (!file.rename(tmp_path, fit_path)) stop("rename of ", tmp_path, " to ", fit_path, " failed")
rm(chk)
logf("saved %s%s", fit_path, if (!is.null(stale)) " (STALE: the manifest changed while sampling; downstream tools need --allow_stale)" else "")
## the summary CSV AFTER the fit file (C2), labelled by axis and covariate (B2)
write.csv(data.frame(obj$param_labels, sm_tab, init_source = init_info$source, init_from = init_info$m1_fit,
                     rho_prior = spec$rho_prior, sigma_prior = spec$sigma_prior, stan_program = family, tag = spec_tag(spec), check.names = FALSE),   # the priors and the tag of the run (2026-10-30)
          file.path(out_dir, paste0("summary_", cell, if (!is.null(stale)) "_STALE" else "", ".csv")), row.names = FALSE)
logf("summary written: %s", paste0("summary_", cell, if (!is.null(stale)) "_STALE" else "", ".csv"))
