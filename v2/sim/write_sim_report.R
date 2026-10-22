## ============================================================================
## write_sim_report.R -- assemble v2/sim/SIM_REPORT.md from the result CSVs written by
## run_sim.R, run_sim_1d.R and run_measurement_tests.R. Numbers are read, never typed.
## Usage: Rscript v2/sim/write_sim_report.R
## ============================================================================
suppressPackageStartupMessages(library(dplyr))
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
sim_dir <- normalizePath(dirname(this_file)); res <- file.path(sim_dir, "results")
rd <- function(f) { p <- file.path(res, f); if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else NULL }
f3 <- function(x) formatC(x, digits = 3, format = "f"); f2 <- function(x) formatC(x, digits = 2, format = "f"); f1 <- function(x) formatC(x, digits = 1, format = "f")
md_table <- function(df, digits = 3) {
  num <- vapply(df, is.numeric, logical(1)); df[num] <- lapply(df[num], function(v) formatC(v, digits = digits, format = "f"))
  c(paste0("| ", paste(names(df), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
L <- c()
add <- function(...) L <<- c(L, ...)
conv <- rd("convergence.csv"); rec <- rd("recovery.csv"); cls <- rd("classification_accuracy.csv"); loo <- rd("loo.csv")
cmp <- rd("loo_compare.csv"); mis <- rd("misclassification_a_vs_b.csv"); ame <- rd("ame_recovery.csv")
stopifnot(!is.null(conv), !is.null(rec), !is.null(cls), !is.null(loo))

add("# Simulation report (v2 pipeline)", "",
    sprintf("Generated %s by `v2/sim/write_sim_report.R` from `v2/sim/results/*.csv`. Seed 123 throughout.", format(Sys.time(), "%Y-%m-%d %H:%M")), "",
    "## Design", "",
    "Synthetic world from `v2/sim/simulate.R`: 40 states, 25 usable years (1000 state-years), full choice sets of 39 alternatives (39,000 alternative rows), one target per initiating state-year. Planted initiator types as offsets a_i in the (dM, dR) plane (bully (+1.6, +1.6), peer (0, 0), challenger (−1.6, −1.6), hybrids (±1.6, ∓1.6), noise SD 0.35); c_i correlated with a_i so Ω has non-zero correlations; Σ with SDs 0.8 and ρ = 0.5; β with non-zero cross terms; τ = 0.5; γ = (0.8, −0.5). Truth M2b adds a K = 2 u'v structure (entries N(0, 0.7)). Because a_i is a type mixture, the 'true Ω' is the realised covariance of (c_i, a_i1, a_i2).", "",
    "Fits: M1a, M1b, M2a, M2b (K = 2) on data from M1b and from M2b, with the FIXED programs (post-review commits 004317c..9374e8d). Container limits (2 cores, 8 GB): 2 chains × 600 iterations (300 warmup) per fit instead of the 4 × 2000 planned for Roar, so ESS is small and intervals are approximate. The one fit that no fix touched (M1b on M1b data, full choice sets, joint M1 program) was kept from the first run; every other fit was rerun after the fixes. The `reused_saved_fit` column below records this. Simulations used full choice sets only (sampled = 0), so the sampled-alternatives fixes (items 3 and 4) do not bear on them.", "")

add("## Convergence", "", md_table(conv |> select(truth, model, reused_saved_fit, elapsed_min, divergences, max_treedepth_hits, worst_rhat, min_ess_bulk, frac_rhat_gt_1.05)), "")

add("## Parameter recovery and 90% interval coverage", "")
cov <- rec |> filter(!fixed_zero) |> group_by(truth_family, model) |> summarise(n_params = n(), coverage90 = mean(covered), mean_abs_err = mean(abs_err), .groups = "drop")
add(md_table(cov), "")
add("Coverage is the share of free structural parameters (α0, α, τ, β (8), γ, Σ, Ω) whose truth lies in the 90% interval; parameters fixed at zero by construction in the 'a' models (Σ[1,2], Ω[2,3]) are excluded from their coverage count. Per-parameter detail follows for the two matched fits (truth = fitted family).", "")
for (tf in c("M1b", "M2b")) {
  d <- rec |> filter(truth_family == tf, model == tf) |> select(param, truth, mean, lo, hi, covered)
  if (nrow(d)) add(sprintf("### %s fitted to %s data", tf, tf), "", md_table(d), "")
}
add("### Independent models on joint data: Σ and Ω cross terms", "")
d <- rec |> filter(param %in% c("Sigma[1,2]", "Omega[2,3]", "Omega[1,2]", "Omega[1,3]")) |> select(truth_family, model, param, truth, mean, lo, hi, covered)
add(md_table(d), "")

add("## Type classification from posterior μ_i (x_ref = 0, r0 = 0.8)", "", md_table(cls |> select(truth, model, accuracy, mean_abs_err_mu)), "")
if (!is.null(mis)) add("Disagreement between the independent (a) and joint (b) classifications, and their accuracies:", "", md_table(mis), "")

add("## LOO (stage 1 and stage 2 kept separate)", "", md_table(loo |> select(truth, model, elpd_stage2, se_stage2, bad_k_stage2, elpd_stage1, se_stage1, bad_k_stage1), 1), "")
if (!is.null(cmp)) add("Pairwise comparisons (positive = second model better):", "", md_table(cmp, 1), "")

if (!is.null(ame)) add("## AME: recovery of u'v (doubly centred truth; U, V Procrustes-aligned across chains)", "", md_table(ame), "",
                       "`cor_UV_offdiag` is the correlation between the posterior-mean u_i'v_j and the true doubly-centred u_i'v_j over off-diagonal pairs; `coverage90_UV_offdiag` the share of pairs whose truth lies in the 90% interval; for M2 fitted to M1b data (no true u'v) `sd_UV_posterior_mean` shows how much structure the model invents. `max_procrustes_rotation_dev` near 0 means the anchor constraints already fixed the rotation.", "")

s1 <- rd("sim1d_recovery.csv"); c1 <- rd("sim1d_classification.csv"); l1 <- rd("sim1d_loo.csv"); v1 <- rd("sim1d_convergence.csv")
if (!is.null(s1)) {
  add("## One-dimensional model (M1m) recovery", "",
      "Smaller world (30 states, 12 years) generated from the material-only truth; M1m (right axis), M1r (wrong axis) and M1b (both axes) fitted to the same data.", "",
      md_table(v1), "", md_table(s1 |> filter(model == "M1m") |> select(param, truth, mean, lo, hi, covered)), "",
      md_table(c1), "", md_table(l1, 1), "")
}
ml <- rd("measurement_loo_ab.csv"); mm <- rd("measurement_misclassification.csv"); th <- rd("measurement_temporal_holdout.csv")
pt <- rd("measurement_permutation.csv"); me <- rd("measurement_error_inflation.csv"); bm <- rd("benchmarks_sim.csv")
if (!is.null(ml)) {
  add("## Measurement tests (small world: 30 states, 12 years, truth M1b)", "",
      "### LOO a vs b", "", md_table(ml, 1), "",
      "### Misclassification a vs b", "", md_table(mm), "",
      "### Temporal holdout (train on the first 9 usable years, score the last 3)", "", md_table(th), "",
      "### Permutation test (dR shuffled within each choice set; rescore mode, 20 permutations)", "", md_table(pt, 1), "",
      "### Systematic measurement error (R inflated by delta for states that initiated in t−1; joint model refit)", "", md_table(me), "")
  if (!is.null(bm)) add("### Benchmarks on the simulated frames (brms)", "", md_table(bm |> select(model, term, Estimate, Q2.5, Q97.5)), "")
}
add("## What failed or is fragile", "",
    "* Chains and iterations are far below the Roar settings; R-hat and ESS in the convergence table are the honest state of these short runs, not of the production fits.",
    "* Ω is recovered as the covariance of a type mixture; its coverage reflects the MVN approximation to planted clusters.",
    "* The AME u'v is identified only up to the U/V scale fixed by the N(0,1) priors and the rotation fixed by the anchors; read the u'v statistics, not U or V.",
    "* The 'a' model on joint data: whether it loses on LOO and misclassifies depends on how far the planted correlations push the data; the numbers above are from one seed.", "")
writeLines(L, file.path(sim_dir, "SIM_REPORT.md"))
cat("wrote", file.path(sim_dir, "SIM_REPORT.md"), "\n")
