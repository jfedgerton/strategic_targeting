## ============================================================================
## 08_opportunity_weights.R -- run 1 (built 2026-11-10; README section 14): the auxiliary MID-propensity model and the OPPORTUNITY WEIGHTS q_ij
##
## Usage: Rscript v2/R/08_opportunity_weights.R (--entry=later | --events=dymid) [--opp_outcome=mid|host3] [--glm_check=TRUE]
##
##   --entry=later | --events=dymid   the FOUNDATION of the cells the weights are for (v2/data_entry | v2/data_dymid): its lag-3 alternatives file gives
##                                    the rows, its panel the initiator's level. The OUTCOME and the MID history always come from the Dyadic MID build
##                                    (v2/data_dymid/events.rds), whatever the foundation. Exactly one of the two is required (bcj cells exist on those only).
##   --opp_outcome=mid | host3        mid (main): any directly interacting Dyadic MID onset between i and j in year t, undirected; host3: onsets whose
##                                    dyadic hostility max(hostlev1, hostlev2) is >= 3
##
## Rows: every directed dyad-year (i, j, t) of alts_lag3_all.rds. The model is FITTED on the January stage-1 IV rows (in_iv = 1: already-fighting at January
## removed; the rows the IV sums over) and predicts every row, the stage-2-only rows included (they need an offset in stage 2).
## Predictors A (opp_weights.R, OPPW_A_COLS; all dated t - 1 or 1 January of t): land contiguity, sea contiguity (conttype 2..5), log capital distance,
## log distance x the initiator's log CINC level (M_lvl; no main effect of the level), an indicator of a Dyadic MID onset in [t - 3, t - 1], the joinable-
## dispute indicator A4 at 1 January of t, decade effects; centred and scaled per training window. NO target-strength variable (Choice A: target selection).
## Fit: logit, EXPANDING WINDOW: year t is predicted from a fit on the years < t; the years <= 1865 (burn-in) are fitted in sample by one fit, which also
## predicts 1866. Saved: ln q at the auxiliary MLE (d0) and for five coherent draws of the FULL coefficient vector (d1..d5; seed 123; oppw_expanding()).
##
## Output (v2/output, gitignored):
##   oppw_<foundation>_<outcome>_lag3.rds         keys (ccode1, ccode2, year), in_iv, y, lq (rows x d0..d5), windows, coefs, z, input hashes, provenance
##   oppw_diag_<foundation>_<outcome>_lag3.csv    every diagnostic table, long (table, by, group, ...)
##   oppw_diag_<foundation>_<outcome>_lag3.md     the same as a report: held-out log loss against the contiguity + distance baseline, precision-recall
##                                                summary, calibration by era, initiator role, initiator region, initiator capability tercile and q bin
##   oppw_coefs_<foundation>_<outcome>_lag3.csv   the MLE coefficients (standardised scale) of every window
## AUC is not reported (uninformative at a positive rate near 0.2%). Nothing here reads a fit or a cluster.
## ============================================================================
suppressPackageStartupMessages({ library(peacesciencer); library(dplyr); library(tidyr); library(digest) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
R_DIR <- dirname(this_file)
source(file.path(R_DIR, "spec.R")); source(file.path(R_DIR, "helpers_data.R")); source(file.path(R_DIR, "opp_weights.R"))
Sys.setenv(V2_ROOT = normalizePath(file.path(R_DIR, "..")))
args <- commandArgs(trailingOnly = TRUE)
refuse_archived_options(args, "08_opportunity_weights.R", also = setdiff(names(SPEC_CHOICES), c("entry", "events", "opp_outcome")),
                        also_msg = "the auxiliary model takes --entry / --events (the foundation) and --opp_outcome only")
OUTCOME <- get_arg(args, "opp_outcome", SPEC_DEFAULTS$opp_outcome)
if (!OUTCOME %in% SPEC_CHOICES$opp_outcome) stop(sprintf("--opp_outcome must be one of %s, got '%s'", paste(SPEC_CHOICES$opp_outcome, collapse = " | "), OUTCOME))
spec <- parse_spec(grep("^--opp_outcome=", args, value = TRUE, invert = TRUE)); spec$opp_outcome <- OUTCOME
if (spec$entry != "later" && spec$events != "dymid") stop("08_opportunity_weights.R needs --entry=later or --events=dymid (the foundations of the bcj cells)")
BURN_IN <- OPPW_BURN_IN; GLM_CHECK <- get_arg(args, "glm_check", TRUE)
FOUNDATION <- oppw_foundation(spec); stopifnot(FOUNDATION %in% c("data_entry", "data_dymid"))
tagname <- sprintf("%s_%s_lag%d", FOUNDATION, OUTCOME, OPPW_LAG)
log_open(paste0("08_oppw_", tagname))
set.seed(V2_SEED)
data_dir <- spec_data_dir(spec, v2_root()); dymid_dir <- file.path(v2_root(), "data_dymid"); out_dir <- file.path(v2_root(), "output"); dir.create(out_dir, showWarnings = FALSE)
stopifnot(dir.exists(data_dir), dir.exists(dymid_dir))
logf("opportunity weights (run 1): foundation %s, outcome %s (%s), lag %d, burn-in through %d, %d full-vector draws, seed %d; R %s, peacesciencer %s",
     FOUNDATION, OUTCOME, OPPW_OUTCOMES[[OUTCOME]], OPPW_LAG, BURN_IN, OPPW_N_DRAWS, V2_SEED, getRversion(), as.character(packageVersion("peacesciencer")))

## ---- inputs, manifest-checked ------------------------------------------------------------------------------------------------------------------------
man <- read.csv(file.path(data_dir, "manifest.csv"), stringsAsFactors = FALSE); man_dy <- read.csv(file.path(dymid_dir, "manifest.csv"), stringsAsFactors = FALSE)
f_alts <- sprintf("alts_lag%d_all.rds", OPPW_LAG); f_panel <- sprintf("stage1_panel_lag%d.rds", OPPW_LAG)
alts <- assert_matches_manifest(file.path(data_dir, f_alts), man)
panel <- assert_matches_manifest(file.path(data_dir, f_panel), man)
ev_fdn <- assert_matches_manifest(file.path(data_dir, "events.rds"), man)
roles_ap <- assert_matches_manifest(file.path(data_dir, "roles_aprime.rds"), man)
ev_dy <- if (FOUNDATION == "data_dymid") ev_fdn else assert_matches_manifest(file.path(dymid_dir, "events.rds"), man_dy)
input <- rbind(data.frame(file = c(f_alts, f_panel), sha256 = man$sha256[match(c(f_alts, f_panel), man$file)], role = "foundation"),
               data.frame(file = "events.rds", sha256 = man_dy$sha256[man_dy$file == "events.rds"], role = "dymid_outcome"))
stopifnot(!anyNA(input$sha256))
logf("inputs: %s (%d rows), %s (%d state-years) of %s; Dyadic MID events.rds of data_dymid (%d rows, %d-%d; hostility >= 3: %d)", f_alts, nrow(alts), f_panel, nrow(panel), FOUNDATION,
     nrow(ev_dy), min(ev_dy$year), max(ev_dy$year), sum(pmax(ev_dy$hostlev1, ev_dy$hostlev2) >= 3L))

## ---- A4 from COW's participant spells (the source of multi_mid_jan under both foundations; 01_build_data.R, psd_hist) -------------------------------
dd <- cow_mid_dirdisps |> mutate(across(c(dispnum, ccode1, ccode2, year), as.integer))
psd <- participant_spell_dates(dd, dispute_date_bounds(cow_mid_disps |> mutate(dispnum = as.integer(dispnum)), dd), part = gml_part)
jp <- oppw_jan_presence(psd, dd, panel |> select(ccode, year))
chk <- panel |> select(ccode, year, multi_mid_jan) |> inner_join(jp$flags, by = c("ccode", "year"), suffix = c("_panel", "_here"))
if (nrow(chk) != nrow(panel) || any(chk$multi_mid_jan_panel != chk$multi_mid_jan_here)) stop("the January presence recomputed here does not reproduce the panel's multi_mid_jan (other peacesciencer data than the build's?)")
logf("A4: %d dispute-years multilateral at 1 January (>= 3 present participants, both sides), %d present participant rows; the panel's multi_mid_jan is reproduced on all %d state-years",
     nrow(jp$multi), nrow(jp$present), nrow(panel))

## ---- the frame ---------------------------------------------------------------------------------------------------------------------------------------
fr <- oppw_frame(alts, panel, ev_dy, jp); rm(alts); invisible(gc())
ycol <- OPPW_OUTCOMES[[OUTCOME]]; iv <- fr$in_iv == 1L
logf("frame: %d directed dyad-years %d-%d (%d without a capital distance dropped); IV rows %d (fitted), stage-2-only rows %d (predicted only); onsets on the IV rows: y_mid %d (%.4f%%), y_host3 %d (%.4f%%)",
     nrow(fr), min(fr$year), max(fr$year), attr(fr, "n_nodist"), sum(iv), sum(!iv), sum(fr$y_mid[iv]), 100 * mean(fr$y_mid[iv]), sum(fr$y_host3[iv]), 100 * mean(fr$y_host3[iv]))
logf("predictor means on the IV rows: %s", paste(sprintf("%s %.4f", c("contig_land", "contig_sea", "dymid_hist3", "a4_joinable"), colMeans(fr[iv, c("contig_land", "contig_sea", "dymid_hist3", "a4_joinable")])), collapse = "; "))

## ---- the expanding-window fits: the model (with the draws) and the contiguity + distance baseline (MLE only) -----------------------------------------
t0 <- Sys.time()
fit <- oppw_expanding(fr, ycol, OPPW_A_COLS, burn_in = BURN_IN, n_draws = OPPW_N_DRAWS, seed = V2_SEED, verbose = TRUE, glm_check = GLM_CHECK)
logf("model: %d windows in %.1f min; coefficients %d (first window) to %d (last); glm.fit check on the first / last window: max coefficient difference %s",
     nrow(fit$windows), as.numeric(difftime(Sys.time(), t0, units = "mins")), fit$windows$p[1], fit$windows$p[nrow(fit$windows)],
     if (GLM_CHECK) paste(sprintf("%.2e", fit$windows$glm_max_diff[c(1, nrow(fit$windows))]), collapse = " / ") else "skipped")
base <- oppw_expanding(fr, ycol, OPPW_A_BASELINE, burn_in = BURN_IN, n_draws = 0L, seed = V2_SEED, verbose = FALSE, glm_check = FALSE)
stopifnot(all(fit$windows$train_year_max[-1] < fit$windows$pred_from[-1]), identical(fit$windows$train_max, base$windows$train_max))

## ---- diagnostics on the IV rows ----------------------------------------------------------------------------------------------------------------------
y <- fr[[ycol]]; lq0 <- fit$lq[, "d0"]; lqb <- base$lq[, "d0"]; ho <- iv & fit$held_out; ins <- iv & !fit$held_out
ll <- rbind(data.frame(rows = "held_out (years > burn-in)", n = sum(ho), n_y = sum(y[ho]), logloss_model = oppw_logloss(y[ho], lq0[ho]), logloss_baseline = oppw_logloss(y[ho], lqb[ho]),
                       logloss_constant = oppw_logloss(y[ho], rep(log(mean(y[ins])), sum(ho)))),
            data.frame(rows = "burn_in (in sample)", n = sum(ins), n_y = sum(y[ins]), logloss_model = oppw_logloss(y[ins], lq0[ins]), logloss_baseline = oppw_logloss(y[ins], lqb[ins]),
                       logloss_constant = oppw_logloss(y[ins], rep(log(mean(y[ins])), sum(ins)))))
ll$model_over_baseline <- ll$logloss_model / ll$logloss_baseline
ll_draws <- data.frame(draw = colnames(fit$lq), logloss_held_out = apply(fit$lq, 2, function(v) oppw_logloss(y[ho], v[ho])), mean_q_held_out = apply(fit$lq, 2, function(v) mean(exp(v[ho]))))
pr <- rbind(cbind(model = "A (run 1)", oppw_pr(y[ho], lq0[ho])), cbind(model = "contiguity + distance baseline", oppw_pr(y[ho], lqb[ho])))
## groups: era; initiator role (the foundation's own mid events, bcj coding; a state-year with two roles counts in both); initiator region; initiator
## capability tercile (terciles of M_lvl over the held-out state-years); q bin
era <- as.character(cut(fr$year, breaks = c(-Inf, BURN_IN, 1913, 1945, 1989, Inf), labels = c(sprintf("<=%d (in sample)", BURN_IN), sprintf("%d-1913", BURN_IN + 1L), "1914-1945", "1946-1989", "1990-")))
rl <- ev_fdn |> filter(in_mid, !escalation) |> select(initiator, target, year) |> left_join(roles_ap |> select(initiator, target, year, grp_aprime), by = c("initiator", "target", "year"))
if (anyNA(rl$grp_aprime) || any(rl$grp_aprime == 3L)) stop("initiator roles: a mid event without a bcj role in roles_aprime.rds")
rl <- distinct(data.frame(ccode1 = rl$initiator, year = rl$year, role = BCJ_GROUPS[ROLE_FROM_APRIME$bcj[rl$grp_aprime + 1L] + 1L]))
k_sy <- oppw_key(fr$ccode1, 1L, fr$year)
is_init <- k_sy %in% oppw_key(rl$ccode1, 1L, rl$year)
sy_ho <- distinct(fr[ho, c("ccode1", "year", "M_lvl_i")]); cuts <- stats::quantile(sy_ho$M_lvl_i, c(1 / 3, 2 / 3))
terc <- as.character(cut(fr$M_lvl_i, breaks = c(-Inf, cuts, Inf), labels = c("low", "middle", "high")))
qbin <- as.character(cut(exp(lq0), breaks = c(0, 1e-4, 1e-3, 1e-2, 1e-1, 1), labels = c("q<1e-4", "1e-4..1e-3", "1e-3..1e-2", "1e-2..0.1", "q>=0.1"), right = FALSE, include.lowest = TRUE))
cal <- rbind(oppw_calibration(y[iv], lq0[iv], era[iv], "era"),
             do.call(rbind, lapply(BCJ_GROUPS, function(r) { i <- ho & k_sy %in% oppw_key(rl$ccode1[rl$role == r], 1L, rl$year[rl$role == r]); oppw_calibration(y[i], lq0[i], rep(r, sum(i)), "initiator_role (held out)") })),
             oppw_calibration(y[ho & !is_init], lq0[ho & !is_init], rep("non_initiator", sum(ho & !is_init)), "initiator_role (held out)"),
             oppw_calibration(y[ho], lq0[ho], oppw_region(fr$ccode1[ho]), "initiator_region (held out)"),
             oppw_calibration(y[ho], lq0[ho], terc[ho], "initiator_capability_tercile (held out)"),
             oppw_calibration(y[ho], lq0[ho], qbin[ho], "q_bin (held out)"))
last <- fit$coefs[fit$coefs$window == max(fit$coefs$window), ]
for (r in seq_len(nrow(ll))) logf("log loss, %s: model %.6f, baseline %.6f (ratio %.4f), constant rate %.6f; n = %d, onsets %d", ll$rows[r], ll$logloss_model[r], ll$logloss_baseline[r], ll$model_over_baseline[r], ll$logloss_constant[r], ll$n[r], ll$n_y[r])
logf("precision-recall (held out): prevalence %.5f; average precision model %.4f, baseline %.4f; precision at recall 0.25 / 0.5 / 0.75: model %.4f / %.4f / %.4f", pr$prevalence[1], pr$avg_precision[1], pr$avg_precision[2], pr$prec_at_rec25[1], pr$prec_at_rec50[1], pr$prec_at_rec75[1])
logf("last window (years <= %d), standardised coefficients: %s", last$train_max[1], paste(sprintf("%s %.3f (%.3f)", last$term[1:(1 + length(OPPW_A_COLS))], last$coef[1:(1 + length(OPPW_A_COLS))], last$se[1:(1 + length(OPPW_A_COLS))]), collapse = "; "))
if (ll$logloss_model[1] >= ll$logloss_baseline[1]) logf("WARNING: the held-out log loss of the model is NOT below the contiguity + distance baseline")

## ---- save --------------------------------------------------------------------------------------------------------------------------------------------
git_commit <- tryCatch(system2("git", c("-C", v2_root(), "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[1], error = function(e) NA_character_)
obj <- list(keys = fr[, c("ccode1", "ccode2", "year")], in_iv = fr$in_iv, y = y, lq = fit$lq, windows = fit$windows, coefs = fit$coefs, z = fit$z, input = input,
            a_cols = OPPW_A_COLS, outcome = OUTCOME, foundation = FOUNDATION, lag = OPPW_LAG, burn_in = BURN_IN, n_draws = OPPW_N_DRAWS, seed = V2_SEED,
            diag = list(logloss = ll, logloss_draws = ll_draws, pr = pr, calibration = cal), git_commit = git_commit, built = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            r_version = as.character(getRversion()))
stopifnot(nrow(obj$keys) == nrow(obj$lq), !anyNA(obj$lq), ncol(obj$lq) == 1L + OPPW_N_DRAWS)
out_path <- file.path(out_dir, paste0("oppw_", tagname, ".rds")); tmp <- paste0(out_path, ".tmp")
saveRDS(obj, tmp); if (!file.rename(tmp, out_path)) stop("rename of ", tmp, " failed")
write.csv(fit$coefs, file.path(out_dir, paste0("oppw_coefs_", tagname, ".csv")), row.names = FALSE)
long <- function(tab, name) { tab <- as.data.frame(tab); data.frame(table = name, row = seq_len(nrow(tab)), tab, check.names = FALSE) }
tabs <- list(logloss = ll, logloss_draws = ll_draws, precision_recall = pr, calibration = cal)
all_cols <- unique(unlist(lapply(tabs, names)))
write.csv(do.call(rbind, lapply(names(tabs), function(nm) { t <- long(tabs[[nm]], nm); for (cc in setdiff(all_cols, names(t))) t[[cc]] <- NA; t[, c("table", "row", all_cols)] })),
          file.path(out_dir, paste0("oppw_diag_", tagname, ".csv")), row.names = FALSE)
md_table <- function(tab) { f <- function(v) if (is.numeric(v)) ifelse(abs(v) >= 1000 | v == round(v), format(v, big.mark = ",", scientific = FALSE, trim = TRUE), formatC(v, digits = 4, format = "g")) else as.character(v)
  c(paste0("| ", paste(names(tab), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(tab)), collapse = "|"), "|"), apply(as.data.frame(lapply(tab, f)), 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))) }
writeLines(c(sprintf("# Opportunity weights, auxiliary model diagnostics (%s)", tagname), "",
             sprintf("Built %s, commit %s. Foundation %s, outcome %s, lag %d, burn-in through %d (in sample), expanding window after it, %d full-vector draws, seed %d. Rows are the January IV rows. No AUC (positive rate %.3f%%).",
                     obj$built, git_commit, FOUNDATION, OUTCOME, OPPW_LAG, BURN_IN, OPPW_N_DRAWS, V2_SEED, 100 * mean(y[iv])), "",
             "## Log loss", "", md_table(ll), "", "## Log loss of the draws (held out)", "", md_table(ll_draws), "", "## Precision-recall (held out)", "", md_table(pr), "",
             "## Calibration", "", md_table(cal), "", sprintf("## Coefficients of the last window (years <= %d, standardised predictors)", last$train_max[1]), "", md_table(last[, c("term", "coef", "se")])),
           file.path(out_dir, paste0("oppw_diag_", tagname, ".md")))
logf("saved %s (%d rows x %d ln q columns), oppw_diag_%s.csv / .md, oppw_coefs_%s.csv", out_path, nrow(obj$lq), ncol(obj$lq), tagname, tagname)
