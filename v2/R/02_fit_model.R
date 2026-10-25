## ============================================================================
## 02_fit_model.R -- assemble Stan data for one (sample, lag, variant, model) cell and fit
##
## Usage:
##   Rscript v2/R/02_fit_model.R --sample=mid --lag=5 --variant=all --model=M1b \
##          [--centrality=pr|eig] [--K=2] [--chains=4] [--iter=2000] [--warmup=1000] \
##          [--sampled_m=0] [--conflict_history=TRUE] [--distance=FALSE] [--contiguity=TRUE]
##          [--dy_peace_years=FALSE] [--drop_ww2=FALSE] [--full_controls=FALSE] [--event_weights=dispute|none] [--direction=cowwar|sideA]
##          [--adapt_delta=0.9] [--max_treedepth=10] [--pilot=FALSE] [--dry_run=FALSE] [--threads_per_chain=1]
##
##   sample    : mid (all MIDs) | fmid (fatal MIDs, main coding) | fmidcow (fatal MIDs, COW-only coding; appendix) | fmidub (upper bound: uncoded fatality -> fatal; appendix) | war (hostility level 5)
##   lag       : 1 | 3 | 5 | 10         change-window length h
##   variant   : all | polrel            choice set: all system members | politically relevant
##   model     : M1m | M1r | M1a | M1b | M2m | M2r | M2a | M2b
##               m = material only, r = relational only, a = both independent, b = both joint; M2 = AME with K factors
##   centrality: pr (IGO PageRank, main) | eig (eigenvector, robustness)
##   sampled_m : 0 = full choice sets; m > 0 = chosen target(s) plus m uniformly sampled alternatives
##   conflict_history: TRUE = include the prior-MID counts in W and Z (main spec; under --full_controls also in X, with
##               the monadic peace years); FALSE = drop them; cell name gets "_noCH"
##   full_controls: FALSE (main, minimal controls since 2026-10-24: W = dM_i, dR_i, ongoing_tm1, log prior MIDs, dem_tm1;
##               X = dM_i, dR_i; Z = contiguity, joint democracy, log dyadic prior MIDs, dy_mid_tm1, target changes) | TRUE
##               reproduces exactly the base specification of commit 8bf0ce1 (W = X = dM_i, dR_i, cowmaj_tm1, dem_tm1,
##               log prior MIDs, log monadic peace years; Z = contiguity, alliance, joint democracy, major-power count,
##               log dyadic prior MIDs, log dyadic peace years, target changes; no ongoing_tm1, no dy_mid_tm1); cell
##               name gets "_fullctrl"
##   distance  : FALSE (main) | TRUE adds log1p(capital-to-capital km) to Z; cell name gets "_dist"
##   contiguity: TRUE (main) | FALSE drops land contiguity from Z; cell name gets "_nocontig"
##   dy_peace_years: FALSE (main) | TRUE restores log1p(dyadic peace years) to Z next to dy_mid_tm1; cell name gets "_dypy"
##   drop_ww2  : FALSE (main) | TRUE removes the WWII events (disputes in WW2_DISPNUMS, 00_helpers.R) from stage 2 and
##               their initiations from stage 1 (appendix) using the build's in_<sample>_noww2 flags (an event collapsed
##               from a WWII and a non-WWII dispute stays when the non-WWII dispute qualifies); cell name gets "_noww2"
##   event_weights: dispute (main; PI decision 2026-10-25, D8) weights every stage-2 event by 1 / (the initiator's targets in the
##               same dispute and onset year), a composite likelihood with dispute-level weighting; none = weight 1 (robustness,
##               cell suffix "_w1")
##   direction : cowwar (main; PI decision 2026-10-25, D1) reads v2/data, built with COW's war-initiating side as the initiator of
##               war-level disputes; sideA reads v2/data_sideA (the appendix rebuild under V2_DIRECTION=sideA, Side A everywhere),
##               cell suffix "_sideA"
##   pilot     : TRUE = 1 chain x 60 iterations (timing only); the fit is saved as fit_<cell>_pilot.rds and every production tool
##               refuses it (assert_production_fit(), 2026-10-25 C1)
##   threads_per_chain: within-chain parallelism via reduce_sum over state-years (default 1 = serial, exact).
##               Request cpus-per-task = chains x threads_per_chain in Slurm.
##   Main specification: --variant=all --sampled_m=0 (full all-states choice sets), the defaults.
##
## Every input file is checked against v2/data/manifest.csv (nrow + sha256) before use
## (audit F3/F8/F9: no silent dependence on which file happens to be on disk).
## Output: v2/output/fit_<cell>.rds containing fit, Stan data, cell description, manifest rows, the sampling configuration and
##         the compiled-model hashes (written atomically: fit_<cell>.rds.tmp then renamed; 2026-10-25 C2), plus
##         v2/output/summary_<cell>.csv (posterior summaries of the structural parameters, labelled by axis and covariate; written
##         AFTER the fit file). If the data manifest changed while sampling the fit is saved as fit_<cell>_STALE.rds with the
##         mismatch recorded in obj$stale instead of being discarded; downstream tools need --allow_stale to read it.
## Seed 123 everywhere (audit F11).
## ============================================================================
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tidyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(this_file), "00_helpers.R"))
Sys.setenv(V2_ROOT = normalizePath(file.path(dirname(this_file), "..")))

## ---- arguments -------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  val <- sub(paste0("^--", name, "="), "", hit[1])
  if (is.logical(default)) return(as.logical(val))
  if (is.numeric(default)) return(as.numeric(val))
  val
}
SAMPLE <- get_arg("sample", "mid"); LAG <- as.integer(get_arg("lag", 5)); VARIANT <- get_arg("variant", "all")
MODEL <- get_arg("model", "M1b"); CENT <- get_arg("centrality", "pr"); K <- as.integer(get_arg("K", 2))
CHAINS <- as.integer(get_arg("chains", 4)); ITER <- as.integer(get_arg("iter", 2000)); WARMUP <- as.integer(get_arg("warmup", 1000))
SAMPLED_M <- as.integer(get_arg("sampled_m", 0)); CONFLICT_HISTORY <- get_arg("conflict_history", TRUE)
ADAPT_DELTA <- get_arg("adapt_delta", 0.9); MAX_TD <- as.integer(get_arg("max_treedepth", 10)); PILOT <- get_arg("pilot", FALSE)
DRY_RUN <- get_arg("dry_run", FALSE)      # TRUE = assemble and validate the Stan data, write it, do NOT sample
THREADS <- as.integer(get_arg("threads_per_chain", 1))   # within-chain parallelism (reduce_sum); cpus = chains x threads
DISTANCE <- get_arg("distance", FALSE)    # TRUE = add log1p(capital distance km) to Z (alternative spec, not main)
CONTIGUITY <- get_arg("contiguity", TRUE) # FALSE = drop land contiguity from Z (robustness option)
DY_PEACE_YEARS <- get_arg("dy_peace_years", FALSE)  # TRUE = put log1p(dyadic peace years) back into Z (left the main Z 2026-10-24)
DROP_WW2 <- get_arg("drop_ww2", FALSE)    # TRUE = appendix run without the WWII events / initiations
FULL_CONTROLS <- get_arg("full_controls", FALSE)  # TRUE = the pre-2026-10-24 control sets (major power, alliance, peace years as well)
EVENT_WEIGHTS <- get_arg("event_weights", "dispute")   # dispute (main) | none (weight 1 on every stage-2 event; "_w1")
DIRECTION <- get_arg("direction", "cowwar")       # cowwar (main, v2/data) | sideA (appendix rebuild in v2/data_sideA; "_sideA")
stopifnot(SAMPLE %in% c("mid", "fmid", "fmidcow", "fmidub", "war"), LAG %in% c(1L, 3L, 5L, 10L), VARIANT %in% c("all", "polrel"),
          MODEL %in% c("M1m", "M1r", "M1a", "M1b", "M2m", "M2r", "M2a", "M2b"), CENT %in% c("pr", "eig"), K >= 1,
          EVENT_WEIGHTS %in% c("dispute", "none"), DIRECTION %in% c("cowwar", "sideA"))
## model code: M1 / M2 = family; m = material only (D_mode 1), r = relational only (D_mode 2),
## a = both dimensions independent (D_mode 3, joint 0), b = both dimensions joint (D_mode 3, joint 1)
D_MODE <- switch(substr(MODEL, 3, 3), m = 1L, r = 2L, a = 3L, b = 3L)
cell <- sprintf("%s_lag%d_%s_%s_%s%s%s%s", SAMPLE, LAG, VARIANT, CENT, MODEL,
                if (grepl("^M2", MODEL)) sprintf("_K%d", K) else "",
                if (!CONFLICT_HISTORY) "_noCH" else "",
                if (SAMPLED_M > 0) sprintf("_samp%d", SAMPLED_M) else "")
cell <- paste0(cell, if (DISTANCE) "_dist" else "", if (!CONTIGUITY) "_nocontig" else "",
               if (DY_PEACE_YEARS) "_dypy" else "", if (DROP_WW2) "_noww2" else "", if (FULL_CONTROLS) "_fullctrl" else "",
               if (EVENT_WEIGHTS == "none") "_w1" else "", if (DIRECTION == "sideA") "_sideA" else "")
if (PILOT) cell <- paste0(cell, "_pilot")          # C1: pilots never share a filename with production fits
log_open(paste0("02_fit_", cell))
logf("cell = %s", cell)
set.seed(V2_SEED)
rstan_options(auto_write = FALSE)
rstan_options(threads_per_chain = THREADS)                 # needs the model compiled with -DSTAN_THREADS (rstan 2.32 default)
options(mc.cores = min(CHAINS, max(1L, parallel::detectCores() %/% THREADS)))
logf("threads_per_chain = %d; mc.cores = %d", THREADS, getOption("mc.cores"))
data_dir <- file.path(v2_root(), if (DIRECTION == "sideA") "data_sideA" else "data"); out_dir <- file.path(v2_root(), "output"); dir.create(out_dir, showWarnings = FALSE)
if (!dir.exists(data_dir)) stop("data directory ", data_dir, " does not exist (for --direction=sideA run `V2_DIRECTION=sideA Rscript v2/R/01_build_data.R` first)")
logf("direction = %s (data directory %s); event weights = %s", DIRECTION, basename(data_dir), EVENT_WEIGHTS)
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
cell_files <- c(sprintf("stage1_panel_lag%d.rds", LAG), sprintf("alts_lag%d_%s.rds", LAG, VARIANT), "events.rds")
stopifnot(all(cell_files %in% manifest_at_load$file))
manifest_cell <- manifest_at_load[manifest_at_load$file %in% cell_files, ]

source(file.path(dirname(this_file), "assemble_cell.R"))
## the manifest SNAPSHOT read above is passed in (review P2-5): every input file's nrow and sha256 are checked against it
## at load, so a rebuild between the snapshot and the load stops the run instead of being read silently from disk
A <- assemble_cell(SAMPLE, LAG, VARIANT, CENT, CONFLICT_HISTORY, SAMPLED_M, DISTANCE = DISTANCE, CONTIGUITY = CONTIGUITY,
                   DY_PEACE_YEARS = DY_PEACE_YEARS, DROP_WW2 = DROP_WW2, FULL_CONTROLS = FULL_CONTROLS, EVENT_WEIGHTS = EVENT_WEIGHTS,
                   data_dir = data_dir, manifest = manifest_at_load)
sy <- A$sy; alts <- A$alts; ev <- A$ev; W_cols <- A$W_cols; X_cols <- A$X_cols; Z_cols <- A$Z_cols
state_codes <- A$state_codes; x_ref <- A$x_ref; n_full <- A$n_full; sampled <- A$sampled

## ---- Stan data ---------------------------------------------------------------------------------------------
joint <- as.integer(grepl("b$", MODEL)); family <- if (grepl("^M1", MODEL)) "m1_latent_additive" else "m2_ame"
dat <- build_stan_data(sy, alts, ev, W_cols, X_cols, Z_cols, joint = joint, sampled = sampled, n_full = n_full,
                       x_ref = x_ref, K = if (family == "m2_ame") K else NULL, allow_orphan_initiations = TRUE,
                       D_mode = D_MODE, grainsize = choose_grainsize(nrow(sy), THREADS))
logf("reduce_sum grainsize = %d (N1 = %d)", dat$grainsize, dat$N1)
logf("D_mode = %d (%s)", D_MODE, c("material only", "relational only", "both dimensions")[D_MODE])
if (family == "m2_ame") logf("AME anchors (most initiations first): ccodes %s", paste(state_codes[dat$anchor], collapse = ", "))
logf("Stan data: S=%d N1=%d N2=%d E=%d P1=%d P=%d Pz=%d joint=%d sampled=%d D_mode=%d; sum(y1)=%d; rows in the IV set %d, in the stage-2 set %d; event weights sum %.1f (min %.4f, %d events with weight 1); W = %s; X = %s; Z = %s",
     dat$S, dat$N1, dat$N2, dat$E, dat$P1, dat$P, dat$Pz, dat$joint, dat$sampled, dat$D_mode, sum(dat$y1), sum(dat$in_iv), sum(dat$in_s2),
     sum(dat$ev_w), min(dat$ev_w), sum(dat$ev_w == 1), paste(W_cols, collapse = ","), paste(X_cols, collapse = ","), paste(Z_cols, collapse = ","))
logf("D2 levels in W: %s; D3 log system size in W: %s; D1 reversed events in the sample: %d", paste(intersect(W_cols, W_LEVELS), collapse = ","),
     "log_n_system_tm1" %in% W_cols, A$n_flipped)
if (DRY_RUN) {
  saveRDS(list(dat = dat, cell = cell, state_codes = state_codes, W_cols = W_cols, X_cols = X_cols, Z_cols = Z_cols,
               manifest = manifest_cell, git_commit = git_commit_at_load, full_controls = FULL_CONTROLS, event_weights = EVENT_WEIGHTS, direction = DIRECTION),
          file.path(out_dir, paste0("standata_", cell, ".rds")))
  logf("DRY RUN: Stan data assembled and saved; no sampling"); quit(save = "no", status = 0)
}
if (PILOT) { ITER <- 60L; WARMUP <- 30L; CHAINS <- 1L; logf("PILOT run: chains=1 iter=60 (timing only); saved with the _pilot suffix") }

## ---- fit ------------------------------------------------------------------------------------------------------------
sm <- readRDS(file.path(v2_root(), "stan", "build", paste0(family, ".rds")))
## compiled-model provenance (2026-10-25 C3): the compiled model_code must be the checked-out .stan source, and the sidecar hash
## written by compile_models.R must match it; both hashes are saved with the fit. (rstan silently ignores unused data, so an old
## binary would otherwise run on new data; this replaces the mu0 / in_iv marker checks of 2026-10-23/24.)
stan_src <- file.path(v2_root(), "stan", paste0(family, ".stan"))
stan_prov <- assert_compiled_current(sm, stan_src, file.path(v2_root(), "stan", "build", paste0(family, ".sha256")))
logf("compiled %s matches %s (source sha256 %s)", family, basename(stan_src), stan_prov$stan_source_sha256)
t0 <- Sys.time()
fit <- sampling(sm, data = dat, chains = CHAINS, iter = ITER, warmup = WARMUP, seed = V2_SEED,
                refresh = max(1L, ITER %/% 20L), control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TD),
                pars = c("z_re", "z_b", "V_raw", "U_free", "L_Om_joint", "L_Om", "re"), include = FALSE)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
logf("sampling finished in %.1f min", elapsed)
sp <- get_sampler_params(fit, inc_warmup = FALSE)
logf("divergences: %d; max-treedepth hits: %d", sum(sapply(sp, function(m) sum(m[, "divergent__"]))),
     sum(sapply(sp, function(m) sum(m[, "treedepth__"] >= MAX_TD))))
pars_struct <- c("alpha0", "alpha", "tau", "mu0", "beta", "sigma_x", "sd_re", "sigma_b", "gamma", "Sigma", "Omega", "R_Om",
                 if (family == "m2_ame") "u_diag")
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
## the fit object (saved below) records the sampling configuration (C1), the compiled-model hashes (C3) and, when the manifest changed
## while sampling, the mismatch (C2: saved as fit_<cell>_STALE.rds instead of being discarded)
obj <- list(fit = fit, dat = dat, cell = cell, state_codes = state_codes, W_cols = W_cols, X_cols = X_cols, Z_cols = Z_cols,
            sy_keys = sy |> select(ccode, year, state_idx, y1), ev = ev, x_ref = x_ref, elapsed_min = elapsed,
            manifest = manifest_cell, data_dir = basename(data_dir), stale = stale,
            param_labels = label_parameters(rownames(sm_tab), W_cols, X_cols, Z_cols, D_MODE),
            args = list(sample = SAMPLE, lag = LAG, variant = VARIANT, model = MODEL, centrality = CENT, K = K, chains = CHAINS,
                        iter = ITER, warmup = WARMUP, pilot = PILOT, sampled_m = SAMPLED_M, conflict_history = CONFLICT_HISTORY,
                        distance = DISTANCE, contiguity = CONTIGUITY, dy_peace_years = DY_PEACE_YEARS, drop_ww2 = DROP_WW2,
                        full_controls = FULL_CONTROLS, event_weights = EVENT_WEIGHTS, direction = DIRECTION, D_mode = D_MODE, threads_per_chain = THREADS,
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
write.csv(data.frame(obj$param_labels, sm_tab, check.names = FALSE), file.path(out_dir, paste0("summary_", cell, if (!is.null(stale)) "_STALE" else "", ".csv")), row.names = FALSE)
logf("summary written: %s", paste0("summary_", cell, if (!is.null(stale)) "_STALE" else "", ".csv"))
