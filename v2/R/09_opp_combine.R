## ============================================================================
## 09_opp_combine.R -- run 1 (built 2026-11-10; README section 14): the equally weighted MIXTURE of the draw-conditional posteriors at lambda = 1
##
## Usage: Rscript v2/R/09_opp_combine.R --sample=fmid --comparison=directed --roles=bcj (--entry=later | --events=dymid) --opp_weight=1 [--opp_outcome=mid]
##          [--tag=long] [--allow_stale=FALSE]
##
## The switches name the POINT cell (q at the auxiliary MLE; --opp_draw must be left out). The script reads its fit and the five fits of the same cell with
## --opp_draw=1..5 (fit_<cell>_d<d>...rds: each conditional on one coherent draw of the full auxiliary coefficient vector), checks that the six are
## production fits of the current data and Stan program, made with the SAME weights file, the same Stan data apart from the offset, the same lambda and
## different offsets, and stacks the posterior draws of the five conditional fits (same number of draws each, so the mixture is equally weighted).
## Reported for mu (mu0_g: the population ideal point of every role), sigma (sigma_g) and rho (rho_x): the point fit's mean and SD, the mixture's mean,
## SD and 5 / 50 / 95% quantiles, the SD of the conditional means across the auxiliary draws, the ratio mixture SD / point SD and whether it exceeds 1.10.
## Output: v2/output/oppw_mixture_<point cell>.csv and oppw_mixture_<point cell>.rds (the table and the stacked draws). M1c only.
## ============================================================================
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(digest) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
R_DIR <- dirname(this_file)
source(file.path(R_DIR, "spec.R")); source(file.path(R_DIR, "helpers_data.R")); source(file.path(R_DIR, "helpers_model.R")); source(file.path(R_DIR, "opp_weights.R"))
Sys.setenv(V2_ROOT = normalizePath(file.path(R_DIR, "..")))
args <- commandArgs(trailingOnly = TRUE)
SAMPLE <- get_arg(args, "sample", "fmid"); MODEL <- get_arg(args, "model", "M1c"); ALLOW_STALE <- get_arg(args, "allow_stale", FALSE)
refuse_archived_options(args, "09_opp_combine.R")
spec <- parse_spec(args)
if (!MODEL %in% OPP_MODELS) stop("09_opp_combine.R: the mixture is specified for ", paste(OPP_MODELS, collapse = ", "), " only")
if (spec_opp_weight(spec) != "1") stop("09_opp_combine.R needs --opp_weight=1 (the auxiliary draws are run at lambda = 1 only)")
if (spec_opp_draw(spec) != "0") stop("09_opp_combine.R takes the POINT cell: leave --opp_draw out (the script reads the draws 1..", OPPW_N_DRAWS, " itself)")
point_cell <- cell_name(SAMPLE, OPPW_LAG, MODEL, spec)
log_open(paste0("09_oppmix_", point_cell))
out_dir <- file.path(v2_root(), "output"); manifest_path <- file.path(spec_data_dir(spec, v2_root()), "manifest.csv")
PARS <- c("mu0_g", "sigma_g", "rho_x")                       # mu, sigma, rho of the M1c role model
## read_one(sp, draw): the checked fit of one cell; returns the posterior draws (iterations x parameters) and what identifies its data
read_one <- function(sp, draw) {
  cell <- cell_name(SAMPLE, OPPW_LAG, MODEL, sp); path <- file.path(out_dir, paste0("fit_", cell, if (ALLOW_STALE && !file.exists(file.path(out_dir, paste0("fit_", cell, ".rds")))) "_STALE" else "", ".rds"))
  if (!file.exists(path)) stop("missing fit: ", path)
  o <- readRDS(path)
  assert_production_fit(o, allow_stale = ALLOW_STALE, what = basename(path)); assert_fit_current(o, manifest_path); assert_fit_stan_current(o, what = basename(path))
  if (!identical(o$cell, cell) || !identical(o$args$model, MODEL) || !identical(o$args$stan_program, "m1_context_offset") || is.null(o$opp) || o$opp$lambda != 1 || o$opp$draw != draw ||
      !identical(o$opp$outcome, spec_opp_outcome(spec))) stop("the fit ", basename(path), " is not the lambda = 1, draw ", draw, " fit of this cell")
  m <- as.matrix(o$fit, pars = PARS)
  stopifnot(all(is.finite(m)))
  res <- list(cell = cell, draws = m, dat_sha = digest::digest(o$dat[setdiff(names(o$dat), c("opp_off", "grainsize"))], algo = "sha256"), off_sha = digest::digest(o$dat$opp_off, algo = "sha256"),
              weights_sha256 = o$opp$weights_sha256, chains = o$args$chains, iter = o$args$iter, warmup = o$args$warmup, group_labels = o$group_labels)
  rm(o); invisible(gc()); res
}
pt <- read_one(spec, 0L)
logf("point fit %s: %d posterior draws (%d chains, iter %d, warmup %d); weights sha256 %s", pt$cell, nrow(pt$draws), pt$chains, pt$iter, pt$warmup, pt$weights_sha256)
dl <- lapply(seq_len(OPPW_N_DRAWS), function(d) { sp <- spec; sp$opp_draw <- as.character(d); r <- read_one(sp, d)
  if (!identical(r$weights_sha256, pt$weights_sha256)) stop("draw ", d, " was fitted with another weights file than the point fit")
  if (!identical(r$dat_sha, pt$dat_sha)) stop("draw ", d, ": its Stan data differ from the point fit's in more than the offset")
  if (!identical(c(r$chains, r$iter, r$warmup), c(pt$chains, pt$iter, pt$warmup))) stop("draw ", d, ": other sampler settings than the point fit")
  logf("draw %d: %s, %d posterior draws", d, r$cell, nrow(r$draws)); r })
offs <- c(pt$off_sha, vapply(dl, function(r) r$off_sha, character(1)))
if (anyDuplicated(offs)) stop("two of the six fits carry the same offset vector (the auxiliary draws must differ from each other and from the MLE)")
mx <- oppw_mixture(lapply(dl, function(r) r$draws), pt$draws)
tab <- mx$table
tab$family <- ifelse(grepl("^mu0_g", tab$parameter), "mu", ifelse(grepl("^sigma_g", tab$parameter), "sigma", "rho"))
## labels: mu0_g[g,d] / sigma_g[g,d] -> role and axis
idx <- regmatches(tab$parameter, regexec("\\[(\\d+),(\\d+)\\]", tab$parameter))
tab$role <- vapply(idx, function(v) if (length(v) == 3L) pt$group_labels[as.integer(v[2])] else "shared", character(1))
tab$axis <- vapply(idx, function(v) if (length(v) == 3L) c("material", "relational")[as.integer(v[3])] else NA_character_, character(1))
tab <- cbind(cell = point_cell, tab[, c("family", "parameter", "role", "axis", setdiff(names(tab), c("family", "parameter", "role", "axis")))])
write.csv(tab, file.path(out_dir, paste0("oppw_mixture_", point_cell, ".csv")), row.names = FALSE)
saveRDS(list(cell = point_cell, table = tab, mix = mx$mix, draw_cells = vapply(dl, function(r) r$cell, character(1)), weights_sha256 = pt$weights_sha256,
             saved = format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file.path(out_dir, paste0("oppw_mixture_", point_cell, ".rds")))
for (fam in c("mu", "sigma", "rho")) { t <- tab[tab$family == fam, ]
  logf("%s: mixture SD / point SD min %.3f, max %.3f; exceeds 1.10 in %d of %d parameters%s", fam, min(t$sd_ratio), max(t$sd_ratio), sum(t$exceeds_10pct), nrow(t),
       if (any(t$exceeds_10pct)) paste0(": ", paste(sprintf("%s (%.2f)", t$parameter[t$exceeds_10pct], t$sd_ratio[t$exceeds_10pct]), collapse = ", ")) else "") }
logf("VERDICT: the mixture SDs %s the point-fit SDs by more than 10%% (%d of %d parameters of mu, sigma, rho); written oppw_mixture_%s.csv / .rds",
     if (any(tab$exceeds_10pct)) "EXCEED" else "do not exceed", sum(tab$exceeds_10pct), nrow(tab), point_cell)
