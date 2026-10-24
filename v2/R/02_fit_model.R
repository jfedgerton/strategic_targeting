## ============================================================================
## 02_fit_model.R -- assemble Stan data for one (sample, lag, variant, model) cell and fit
##
## Usage:
##   Rscript v2/R/02_fit_model.R --sample=mid --lag=5 --variant=all --model=M1b \
##          [--centrality=pr|eig] [--K=2] [--chains=4] [--iter=2000] [--warmup=1000] \
##          [--sampled_m=0] [--conflict_history=TRUE] [--distance=FALSE] [--contiguity=TRUE]
##          [--dy_peace_years=FALSE] [--drop_ww2=FALSE]
##          [--adapt_delta=0.9] [--max_treedepth=10] [--pilot=FALSE] [--dry_run=FALSE] [--threads_per_chain=1]
##
##   sample    : mid (all MIDs) | fmid (fatal MIDs, main coding) | fmidcow (fatal MIDs, COW-only coding; appendix) | fmidub (upper bound: uncoded fatality -> fatal; appendix) | war (hostility level 5)
##   lag       : 1 | 3 | 5 | 10         change-window length h
##   variant   : all | polrel            choice set: all system members | politically relevant
##   model     : M1m | M1r | M1a | M1b | M2m | M2r | M2a | M2b
##               m = material only, r = relational only, a = both independent, b = both joint; M2 = AME with K factors
##   centrality: pr (IGO PageRank, main) | eig (eigenvector, robustness)
##   sampled_m : 0 = full choice sets; m > 0 = chosen target(s) plus m uniformly sampled alternatives
##   conflict_history: TRUE = include conflict-history controls in W, X and Z (main spec); FALSE = drop them
##   distance  : FALSE (main) | TRUE adds log1p(capital-to-capital km) to Z; cell name gets "_dist"
##   contiguity: TRUE (main) | FALSE drops land contiguity from Z; cell name gets "_nocontig"
##   dy_peace_years: FALSE (main) | TRUE restores log1p(dyadic peace years) to Z next to dy_mid_tm1; cell name gets "_dypy"
##   drop_ww2  : FALSE (main) | TRUE removes the WWII events (disputes in WW2_DISPNUMS, 00_helpers.R) from stage 2 and
##               their initiations from stage 1 (appendix); cell name gets "_noww2"
##   threads_per_chain: within-chain parallelism via reduce_sum over state-years (default 1 = serial, exact).
##               Request cpus-per-task = chains x threads_per_chain in Slurm.
##   Main specification: --variant=all --sampled_m=0 (full all-states choice sets), the defaults.
##
## Every input file is checked against v2/data/manifest.csv (nrow + sha256) before use
## (audit F3/F8/F9: no silent dependence on which file happens to be on disk).
## Output: v2/output/fit_<cell>.rds containing fit, Stan data, cell description, manifest rows,
##         plus v2/output/summary_<cell>.csv (posterior summaries of the structural parameters).
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
stopifnot(SAMPLE %in% c("mid", "fmid", "fmidcow", "fmidub", "war"), LAG %in% c(1L, 3L, 5L, 10L), VARIANT %in% c("all", "polrel"),
          MODEL %in% c("M1m", "M1r", "M1a", "M1b", "M2m", "M2r", "M2a", "M2b"), CENT %in% c("pr", "eig"), K >= 1)
## model code: M1 / M2 = family; m = material only (D_mode 1), r = relational only (D_mode 2),
## a = both dimensions independent (D_mode 3, joint 0), b = both dimensions joint (D_mode 3, joint 1)
D_MODE <- switch(substr(MODEL, 3, 3), m = 1L, r = 2L, a = 3L, b = 3L)
cell <- sprintf("%s_lag%d_%s_%s_%s%s%s%s", SAMPLE, LAG, VARIANT, CENT, MODEL,
                if (grepl("^M2", MODEL)) sprintf("_K%d", K) else "",
                if (!CONFLICT_HISTORY) "_noCH" else "",
                if (SAMPLED_M > 0) sprintf("_samp%d", SAMPLED_M) else "")
cell <- paste0(cell, if (DISTANCE) "_dist" else "", if (!CONTIGUITY) "_nocontig" else "",
               if (DY_PEACE_YEARS) "_dypy" else "", if (DROP_WW2) "_noww2" else "")
log_open(paste0("02_fit_", cell))
logf("cell = %s", cell)
set.seed(V2_SEED)
rstan_options(auto_write = FALSE)
rstan_options(threads_per_chain = THREADS)                 # needs the model compiled with -DSTAN_THREADS (rstan 2.32 default)
options(mc.cores = min(CHAINS, max(1L, parallel::detectCores() %/% THREADS)))
logf("threads_per_chain = %d; mc.cores = %d", THREADS, getOption("mc.cores"))
data_dir <- file.path(v2_root(), "data"); out_dir <- file.path(v2_root(), "output"); dir.create(out_dir, showWarnings = FALSE)
manifest_path <- file.path(data_dir, "manifest.csv")

source(file.path(dirname(this_file), "assemble_cell.R"))
A <- assemble_cell(SAMPLE, LAG, VARIANT, CENT, CONFLICT_HISTORY, SAMPLED_M, DISTANCE = DISTANCE, CONTIGUITY = CONTIGUITY,
                   DY_PEACE_YEARS = DY_PEACE_YEARS, DROP_WW2 = DROP_WW2, data_dir = data_dir)
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
logf("Stan data: S=%d N1=%d N2=%d E=%d P1=%d P=%d Pz=%d joint=%d sampled=%d D_mode=%d; sum(y1)=%d; W = %s; Z = %s",
     dat$S, dat$N1, dat$N2, dat$E, dat$P1, dat$P, dat$Pz, dat$joint, dat$sampled, dat$D_mode, sum(dat$y1),
     paste(W_cols, collapse = ","), paste(Z_cols, collapse = ","))
if (DRY_RUN) {
  saveRDS(list(dat = dat, cell = cell, state_codes = state_codes, W_cols = W_cols, X_cols = X_cols, Z_cols = Z_cols),
          file.path(out_dir, paste0("standata_", cell, ".rds")))
  logf("DRY RUN: Stan data assembled and saved; no sampling"); quit(save = "no", status = 0)
}
if (PILOT) { ITER <- 60L; WARMUP <- 30L; CHAINS <- 1L; logf("PILOT run: chains=1 iter=60 (timing only)") }

## ---- fit ------------------------------------------------------------------------------------------------------------
sm <- readRDS(file.path(v2_root(), "stan", "build", paste0(family, ".rds")))
## audit 2026-10-23: rstan silently ignores unused data, so an old compiled program would run on new data without mu0.
if (!grepl("row_vector[D] mu0;", sm@model_code, fixed = TRUE)) stop(sprintf("compiled %s.rds has no mu0 parameter: stale binary; rerun v2/R/compile_models.R", family))
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
write.csv(data.frame(parameter = rownames(sm_tab), sm_tab), file.path(out_dir, paste0("summary_", cell, ".csv")), row.names = FALSE)
logf("worst Rhat among structural parameters: %.3f; min n_eff: %.0f",
     max(sm_tab[, "Rhat"], na.rm = TRUE), min(sm_tab[, "n_eff"], na.rm = TRUE))
man <- read.csv(manifest_path, stringsAsFactors = FALSE)
saveRDS(list(fit = fit, dat = dat, cell = cell, state_codes = state_codes, W_cols = W_cols, X_cols = X_cols, Z_cols = Z_cols,
             sy_keys = sy |> select(ccode, year, state_idx, y1), ev = ev, x_ref = x_ref, elapsed_min = elapsed,
             manifest = man[man$file %in% c(sprintf("stage1_panel_lag%d.rds", LAG), sprintf("alts_lag%d_%s.rds", LAG, VARIANT), "events.rds"), ],
             args = list(sample = SAMPLE, lag = LAG, variant = VARIANT, model = MODEL, centrality = CENT, K = K, chains = CHAINS,
                         iter = ITER, warmup = WARMUP, sampled_m = SAMPLED_M, conflict_history = CONFLICT_HISTORY,
                         distance = DISTANCE, contiguity = CONTIGUITY, dy_peace_years = DY_PEACE_YEARS, drop_ww2 = DROP_WW2,
                         D_mode = D_MODE, threads_per_chain = THREADS,
                         grainsize = dat$grainsize, seed = V2_SEED,
                         git_commit = tryCatch(system2("git", c("-C", v2_root(), "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[1],
                                               error = function(e) NA_character_))),
        file.path(out_dir, paste0("fit_", cell, ".rds")))
logf("saved %s", file.path(out_dir, paste0("fit_", cell, ".rds")))
