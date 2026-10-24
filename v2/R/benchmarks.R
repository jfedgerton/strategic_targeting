## ============================================================================
## benchmarks.R -- two comparison models, built on the SAME v2 data files and timing
##
## (i)  Directed-dyad logit of initiation (i -> j in year t) on the CINC ratio, the centrality
##      difference dR, and the dyadic controls. Rows = every alternative row of the chosen choice
##      set file (state-year x alternative), outcome = 1 if (i, j, t) is an event in the sample.
##      Estimated with brms bernoulli (seed 123). This is the "standard" benchmark.
## (ii) The OLD gap regression: bivariate brms mvbind(dM, dR) on initiator power changes and
##      controls, restricted to EVENT rows (the realised target's position relative to the
##      initiator), with initiator and target random intercepts. Rebuilt with the fixed data:
##      positions at t-1 and change windows [t-1-h, t-1] (audit F5), Side A originators (F4).
##      It is kept only as a benchmark; the v2 model replaces it.
##
## Usage: Rscript v2/R/benchmarks.R --sample=mid --lag=5 --variant=polrel [--centrality=pr]
##        [--chains=4] [--iter=2000] [--warmup=1000] [--max_rows=0] [--dry_run=FALSE]
## The directed-dyad logit on the "all" variant has ~1.7M rows; use --variant=polrel (default
## here) or --max_rows to subsample non-events (then the intercept is offset-corrected: see code).
## Outputs: v2/output/bench_dyadlogit_<cell>.rds, bench_gap_<cell>.rds and summary CSVs.
## Tested on simulated data by v2/sim/run_measurement_tests.R (dry assembly + short brms fit).
## ============================================================================
suppressPackageStartupMessages({ library(brms); library(dplyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
## when sourced from another script (v2/sim/run_measurement_tests.R) the helpers are already loaded
if (!exists("build_stan_data")) {
  source(file.path(dirname(this_file), "00_helpers.R"))
  Sys.setenv(V2_ROOT = normalizePath(file.path(dirname(this_file), "..")))
}

## ---- assembly functions (also used by the sim test) -------------------------------------------------
## alts: alternative rows with dM, dR, Z columns, ccode1/ccode2/year (or state_idx/target_idx/year);
## ev_keys: data.frame of event (ccode1, ccode2, year) rows in the sample
assemble_dyad_logit <- function(alts, ev_keys, id_cols = c("ccode1", "ccode2", "year"),
                                z_cols = c("contig_land", "alliance_def_off", "joint_dem", "major_count"),
                                cinc_cols = NULL) {
  stopifnot(all(id_cols %in% names(alts)), all(id_cols %in% names(ev_keys)), all(z_cols %in% names(alts)))
  ev_keys$y <- 1L
  d <- alts |> left_join(ev_keys |> distinct(across(all_of(id_cols)), y), by = id_cols) |> mutate(y = coalesce(y, 0L))
  stopifnot(sum(d$y) <= nrow(ev_keys), sum(d$y) > 0)
  ## CINC ratio: the model has dM = M_i - M_j, a log ratio of CINC shares divided by the gap constant GAP_SCALE_M
  ## (change_scale.csv, h = NA); the "ratio" benchmark uses the raw share when available
  ## (cinc_i / (cinc_i + cinc_j)); otherwise dM stands in and this is logged.
  if (!is.null(cinc_cols) && all(cinc_cols %in% names(d))) d$cinc_ratio <- d[[cinc_cols[1]]] / (d[[cinc_cols[1]]] + d[[cinc_cols[2]]])
  else d$cinc_ratio <- d$dM
  d
}
assemble_gap_regression <- function(sy, alts, ev, X_cols) {
  ## event rows: initiator state-year covariates X and the realised target's (dM, dR)
  stopifnot(all(ev$sy_row %in% seq_len(nrow(sy))), all(ev$alt_row %in% seq_len(nrow(alts))))
  d <- sy[ev$sy_row, c("state_idx", "year", X_cols), drop = FALSE]
  d$target_idx <- alts$target_idx[ev$alt_row]; d$dM <- alts$dM[ev$alt_row]; d$dR <- alts$dR[ev$alt_row]
  stopifnot(nrow(d) == nrow(ev), !anyNA(d$dM), !anyNA(d$dR))
  d$initiator <- factor(d$state_idx); d$target <- factor(d$target_idx)
  d
}
gap_formula <- function(X_cols) {
  rhs <- paste(c(X_cols, "(1 | p | initiator)", "(1 | q | target)"), collapse = " + ")
  bf(as.formula(paste("mvbind(dM, dR) ~", rhs))) + set_rescor(TRUE)
}

## ---- CLI -------------------------------------------------------------------------------------------------
if (sys.nframe() == 0 && length(grep("^--sample=", commandArgs(TRUE))) == 1) {
  args <- commandArgs(TRUE)
  source(file.path(dirname(this_file), "assemble_cell.R"))   # z_columns() and SY_INDICATORS: the model's own control lists
  get_arg <- function(name, default) { h <- grep(paste0("^--", name, "="), args, value = TRUE); if (!length(h)) return(default)
    v <- sub(paste0("^--", name, "="), "", h[1]); if (is.logical(default)) as.logical(v) else if (is.numeric(default)) as.numeric(v) else v }
  SAMPLE <- get_arg("sample", "mid"); LAG <- as.integer(get_arg("lag", 5)); VARIANT <- get_arg("variant", "polrel"); CENT <- get_arg("centrality", "pr")
  CHAINS <- as.integer(get_arg("chains", 4)); ITER <- as.integer(get_arg("iter", 2000)); WARMUP <- as.integer(get_arg("warmup", 1000))
  MAX_ROWS <- as.integer(get_arg("max_rows", 0)); DRY <- get_arg("dry_run", FALSE)
  cell <- sprintf("%s_lag%d_%s_%s", SAMPLE, LAG, VARIANT, CENT)
  log_open(paste0("benchmarks_", cell)); set.seed(V2_SEED)
  data_dir <- file.path(v2_root(), "data"); out_dir <- file.path(v2_root(), "output"); dir.create(out_dir, showWarnings = FALSE)
  man <- file.path(data_dir, "manifest.csv")
  panel <- assert_matches_manifest(file.path(data_dir, sprintf("stage1_panel_lag%d.rds", LAG)), man)
  alts <- assert_matches_manifest(file.path(data_dir, sprintf("alts_lag%d_%s.rds", LAG, VARIANT)), man)
  events <- assert_matches_manifest(file.path(data_dir, "events.rds"), man)
  dR_alt <- if (CENT == "pr") "dR_pr" else "dR_eig"; dR_h <- if (CENT == "pr") "dRpr_h" else "dReig_h"
  alts <- alts |> mutate(dR = .data[[dR_alt]], log_dy_prior_mids = log1p(dy_prior_mids_h), log_dy_peace_years = log1p(dy_peace_years_h))
  ev_keys <- events |> filter(.data[[paste0("in_", SAMPLE)]]) |> transmute(ccode1 = initiator, ccode2 = target, year)
  z_cols <- z_columns()                      # the main-spec Z (audit 2026-10-24: no separate hard-coded list here)
  stopifnot(all(z_cols %in% names(alts)))
  ## (i) directed-dyad logit
  d1 <- assemble_dyad_logit(alts, ev_keys, z_cols = z_cols)
  offset_term <- 0
  if (MAX_ROWS > 0 && nrow(d1) > MAX_ROWS) {
    ## keep all events, subsample non-events; the intercept is corrected by log(sampling fraction)
    ne <- d1 |> filter(y == 1L); nn <- d1 |> filter(y == 0L)
    frac <- (MAX_ROWS - nrow(ne)) / nrow(nn); nn <- nn |> slice_sample(prop = frac)
    d1 <- bind_rows(ne, nn); offset_term <- -log(frac)
    logf("non-events subsampled at fraction %.4f; intercept offset %.3f applied", frac, offset_term)
  }
  d1$off <- offset_term
  logf("dyad logit rows %d, events %d", nrow(d1), sum(d1$y))
  f1 <- bf(as.formula(paste("y ~ cinc_ratio + dR +", paste(z_cols, collapse = " + "), "+ offset(off) + (1 | ccode1) + (1 | ccode2)")), family = bernoulli())
  ## (ii) old gap regression on event rows (initiator X from the panel)
  sy <- panel |> mutate(state_idx = ccode, dM_i = dM_h, dR_i = .data[[dR_h]], log_prior_mids = log1p(mon_prior_mids_h), log_peace_years = log1p(mon_peace_years_h))
  X_cols <- c("dM_i", "dR_i", SY_INDICATORS, "log_prior_mids", "log_peace_years")   # same monadic controls as the model (ongoing_tm1 included)
  sy_key <- sy |> transmute(ccode1 = ccode, year, sy_row = row_number())
  alts2 <- alts |> inner_join(sy_key, by = c("ccode1", "year")) |> arrange(sy_row, ccode2) |> mutate(target_idx = ccode2)
  ev <- ev_keys |> inner_join(sy_key, by = c("ccode1", "year")) |>
    inner_join(alts2 |> transmute(ccode1, ccode2, year, alt_row = row_number()), by = c("ccode1", "ccode2", "year")) |> select(sy_row, alt_row)
  d2 <- assemble_gap_regression(sy, alts2, ev, X_cols)
  logf("gap regression rows %d", nrow(d2))
  if (DRY) { logf("DRY RUN: benchmark data assembled (%d, %d rows); no fits", nrow(d1), nrow(d2)); quit(save = "no") }
  fit1 <- brm(f1, data = d1, prior = c(set_prior("normal(0, 1)", class = "b"), set_prior("student_t(3, 0, 1.5)", class = "sd")),
              chains = CHAINS, cores = CHAINS, iter = ITER, warmup = WARMUP, seed = V2_SEED, refresh = 0)
  saveRDS(list(fit = fit1, cell = cell, n = nrow(d1), offset = offset_term), file.path(out_dir, paste0("bench_dyadlogit_", cell, ".rds")))
  write.csv(as.data.frame(fixef(fit1)), file.path(out_dir, paste0("bench_dyadlogit_", cell, ".csv")))
  fit2 <- brm(gap_formula(X_cols), data = d2, prior = set_prior("normal(0, 1)", class = "b"),
              chains = CHAINS, cores = CHAINS, iter = ITER, warmup = WARMUP, seed = V2_SEED, refresh = 0)
  saveRDS(list(fit = fit2, cell = cell, n = nrow(d2)), file.path(out_dir, paste0("bench_gap_", cell, ".rds")))
  write.csv(as.data.frame(fixef(fit2)), file.path(out_dir, paste0("bench_gap_", cell, ".csv")))
  logf("benchmarks done")
}
