## spec.R: the main specification, the switches, their cell-name suffixes and the sampler defaults (2026-10-25)
suppressPackageStartupMessages({ library(testthat) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "spec.R"))

test_that("the defaults are the main specification and every switch value is validated", {
  main <- parse_spec(character(0))
  expect_equal(main[names(SPEC_CHOICES)], list(pairs = "all", context = "all", coal_ties = "inclusive", coal_drop = "none", direction = "cowwar", ww2 = "include", weights = "dispute", controls = "minimal", changes = "raw"))
  expect_equal(main$ame_K, 2L); expect_true(spec_is_main(main)); expect_equal(spec_suffix(main), "")
  expect_equal(parse_spec("--pairs=joiners")$pairs, "joiners"); expect_false(spec_is_main(parse_spec("--pairs=joiners")))
  expect_error(parse_spec("--pairs=both"), "--pairs must be one of")
  expect_error(parse_spec("--ww2=TRUE"), "--ww2 must be one of")
  expect_error(parse_spec("--ame_K=0"), "positive integer")
  expect_false(spec_is_main(parse_spec("--ame_K=3")))
  ## get_arg types its result like the default and fails on garbage
  expect_equal(get_arg("--lag=3", "lag", 5), 3); expect_true(get_arg("--dry_run=TRUE", "dry_run", FALSE)); expect_equal(get_arg(character(0), "x", "y"), "y")
  expect_error(get_arg("--lag=three", "lag", 5), "numeric"); expect_error(get_arg("--dry_run=maybe", "dry_run", FALSE), "TRUE or FALSE")
})

test_that("main-spec cell names are the pre-2026-10-25 names, and every non-default switch adds its own suffix", {
  main <- parse_spec(character(0))
  expect_equal(cell_name("fmid", 3L, "M1b", main), "fmid_lag3_all_pr_M1b")
  expect_equal(cell_name("war", 3L, "M2b", main), "war_lag3_all_pr_M2b_K2")
  expect_equal(cell_name("mid", 10L, "M1m", main), "mid_lag10_all_pr_M1m")
  ## the four pre-existing options keep their suffixes
  expect_equal(cell_name("war", 3L, "M1b", parse_spec("--ww2=drop")), "war_lag3_all_pr_M1b_noww2")
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec("--controls=full")), "fmid_lag3_all_pr_M1b_fullctrl")
  expect_equal(cell_name("war", 3L, "M1b", parse_spec("--weights=none")), "war_lag3_all_pr_M1b_w1")
  expect_equal(cell_name("war", 3L, "M2b", parse_spec("--direction=sideA")), "war_lag3_all_pr_M2b_K2_sideA")
  expect_equal(cell_name("war", 3L, "M2b", parse_spec(c("--ww2=drop", "--direction=sideA"))), "war_lag3_all_pr_M2b_K2_noww2_sideA")   # the old order
  ## the new switches
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec("--pairs=originators")), "fmid_lag3_all_pr_M1b_orig")
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec("--pairs=joiners")), "fmid_lag3_all_pr_M1b_join")
  expect_equal(cell_name("fmid", 3L, "M1b", parse_spec("--changes=winsor")), "fmid_lag3_all_pr_M1b_winsor")
  expect_equal(cell_name("fmid", 3L, "M2a", parse_spec("--ame_K=3")), "fmid_lag3_all_pr_M2a_K3")
  all_on <- parse_spec(c("--pairs=originators", "--direction=sideA", "--ww2=drop", "--weights=none", "--controls=full", "--changes=winsor"))
  expect_equal(spec_suffix(all_on), "_orig_noww2_fullctrl_w1_winsor_sideA")
  ## every switch value has a suffix entry, and the suffixes are distinct
  for (nm in names(SPEC_CHOICES)) for (v in SPEC_CHOICES[[nm]][-1]) expect_true(!is.null(SPEC_SUFFIX[[nm]][[v]]))
  expect_false(any(duplicated(unlist(SPEC_SUFFIX))))
  expect_error(cell_name("mid", 7L, "M1b", main)); expect_error(cell_name("mid", 3L, "M3b", main))
  expect_equal(basename(spec_data_dir(main, "/r")), "data"); expect_equal(basename(spec_data_dir(parse_spec("--direction=sideA"), "/r")), "data_sideA")
})

test_that("model codes and sampler settings", {
  expect_equal(model_D_mode("M1m"), 1L); expect_equal(model_D_mode("M2r"), 2L); expect_equal(model_D_mode("M1a"), 3L); expect_equal(model_D_mode("M2b"), 3L)
  expect_equal(model_joint("M1b"), 1L); expect_equal(model_joint("M1a"), 0L); expect_equal(model_family("M2m"), "M2")
  s1 <- sampler_settings("M1b"); s2 <- sampler_settings("M2b")
  expect_equal(s1, list(chains = 4L, iter = 2000L, warmup = 1000L, adapt_delta = 0.9, max_treedepth = 10L))
  expect_equal(s2, SAMPLER$M2)
  o <- sampler_settings("M2b", c("--iter=300", "--warmup=100", "--adapt_delta=0.99"))
  expect_equal(o$iter, 300L); expect_equal(o$warmup, 100L); expect_equal(o$adapt_delta, 0.99); expect_equal(o$chains, SAMPLER$M2$chains)
  ## the robustness presets name real switches and values
  for (p in names(ROBUSTNESS_PRESETS)) for (nm in names(ROBUSTNESS_PRESETS[[p]])) expect_true(ROBUSTNESS_PRESETS[[p]][[nm]] %in% SPEC_CHOICES[[nm]])
  expect_true(all(vapply(strsplit(ROBUSTNESS_CELLS, ":"), length, 1L) == 3L))
  ## submit_robustness.sh hard-codes the default cell list (no R on the login node); it must be ROBUSTNESS_CELLS (audit fix 3)
  sh <- readLines(file.path(root, "slurm", "submit_robustness.sh"))
  dflt <- sub('^DEFAULT_CELLS="([^"]*)".*$', "\\1", grep("^DEFAULT_CELLS=", sh, value = TRUE))
  expect_equal(strsplit(dflt, ",")[[1]], ROBUSTNESS_CELLS)
  ## the main-grid submit script's default index lists are the M1 / M2 halves of the 0-95 mapping (lag slowest, model, sample)
  sh <- readLines(file.path(root, "slurm", "submit_main_grid.sh"))
  idx <- function(nm) sub(sprintf('^%s="\\$\\{%s:-([^}]*)\\}".*$', nm, nm), "\\1", grep(sprintf("^%s=", nm), sh, value = TRUE))
  expand <- function(s) unlist(lapply(strsplit(s, ",")[[1]], function(r) { ab <- as.integer(strsplit(r, "-")[[1]]); ab[1]:ab[2] }))
  grid <- expand.grid(sample = seq_along(MAIN_SAMPLES), model = seq_along(MODELS), lag = seq_along(LAGS))      # sample fastest, lag slowest = index order
  expect_equal(expand(idx("M1_IDX")), which(model_family(MODELS[grid$model]) == "M1") - 1L)
  expect_equal(expand(idx("M2_IDX")), which(model_family(MODELS[grid$model]) == "M2") - 1L)
})

test_that("audit fixes 14 and 16: archived flags are refused by one shared helper; ame_K != 2 with an M1 model is an error", {
  expect_true(refuse_archived_options(c("--sample=mid", "--ww2=drop", "--ame_K=2"), "x"))
  for (old in ARCHIVED_OPTIONS) expect_error(refuse_archived_options(paste0("--", old, "=1"), "02_fit_model.R"), paste0("--", old, " is not an option"), fixed = TRUE)
  expect_error(refuse_archived_options("--drop_ww2", "x"), "--drop_ww2")                                     # the bare flag form too
  expect_error(refuse_archived_options("--ww2=drop", "01_build_data.R", also = setdiff(names(SPEC_CHOICES), "direction")), "--ww2 is not an option")
  expect_true(refuse_archived_options("--direction=sideA", "01_build_data.R", also = setdiff(names(SPEC_CHOICES), "direction")))
  ## every script that parses the switches calls the helper
  for (f in c("01_build_data.R", "02_fit_model.R", "04_compare.R", "05_robustness.R", "06_context_compare.R"))
    expect_true(any(grepl("refuse_archived_options(", readLines(file.path(root, "R", f)), fixed = TRUE)), info = f)
  k3 <- parse_spec("--ame_K=3")
  expect_equal(cell_name("mid", 3L, "M2b", k3), "mid_lag3_all_pr_M2b_K3")
  expect_error(cell_name("mid", 3L, "M1b", k3), "K applies to the M2 family only")
  expect_equal(cell_name("mid", 3L, "M1b", parse_spec("--ame_K=2")), "mid_lag3_all_pr_M1b")
})

test_that("review P2 (2026-10-26): comparison and measurement names carry a non-default AME K; the M1 names never do", {
  main <- parse_spec(character(0)); k1 <- parse_spec("--ame_K=1"); k3 <- parse_spec(c("--ame_K=3", "--ww2=drop"))
  ## 04_compare.R group: the documented main-grid name for K = 2, _K<K> before the switch suffix otherwise
  expect_equal(compare_group_name("mid", 5L, main), "mid_lag5_all_pr")
  expect_equal(compare_group_name("fmid", 3L, parse_spec("--ww2=drop")), "fmid_lag3_all_pr_noww2")
  expect_equal(compare_group_name("mid", 3L, k1), "mid_lag3_all_pr_K1")
  expect_equal(compare_group_name("mid", 3L, k3), "mid_lag3_all_pr_K3_noww2")
  ## the per-model fit names of a comparison group: M1 fits are K-free, M2 fits carry the group's K
  names_for <- function(spec) vapply(MODELS, function(m) cell_name("mid", 3L, m, cell_spec(m, spec)), "")
  n1 <- names_for(k1); n2 <- names_for(main); n3 <- names_for(parse_spec("--ame_K=3"))
  expect_equal(unname(n1[c("M1b", "M2b")]), c("mid_lag3_all_pr_M1b", "mid_lag3_all_pr_M2b_K1"))
  expect_equal(unname(n3[c("M1a", "M2a")]), c("mid_lag3_all_pr_M1a", "mid_lag3_all_pr_M2a_K3"))
  expect_equal(n1[startsWith(MODELS, "M1")], n2[startsWith(MODELS, "M1")]); expect_equal(n1[startsWith(MODELS, "M1")], n3[startsWith(MODELS, "M1")])
  expect_false(any(n1[startsWith(MODELS, "M2")] %in% n2)); expect_false(any(n3[startsWith(MODELS, "M2")] %in% n2))
  ## 05_robustness.R perm_error outputs: the M2 name always has its K (cell_name() convention), so K = 2 and K = 3 never share a file name
  expect_equal(measurement_cell_name("mid", 5L, "m1", main), "mid_lag5_all_pr_M1")
  expect_equal(measurement_cell_name("mid", 5L, "m2", main), "mid_lag5_all_pr_M2_K2")
  expect_equal(measurement_cell_name("mid", 5L, "m2", parse_spec("--ame_K=3")), "mid_lag5_all_pr_M2_K3")
  expect_equal(measurement_cell_name("mid", 5L, "m2", k3), "mid_lag5_all_pr_M2_K3_noww2")
  expect_equal(measurement_cell_name("mid", 5L, "m1", k3), "mid_lag5_all_pr_M1_noww2")                  # the M1 tests ignore K
  ks <- vapply(1:4, function(k) measurement_cell_name("mid", 5L, "m2", parse_spec(sprintf("--ame_K=%d", k))), "")
  expect_equal(anyDuplicated(ks), 0L)
  expect_error(measurement_cell_name("mid", 5L, "M2", main)); expect_error(compare_group_name("mid", 7L, main))
})

test_that("PI request 2026-10-26: the lag-3 K = 1 M2 array is submitted on the same M1 array as the K = 2 one, with its own switches and distinct names", {
  ## submit_main_grid.sh against a stub sbatch that records its arguments and the SWITCHES it would export, and returns job ids 1001, 1002, ...
  tmp <- tempfile("sub_"); dir.create(file.path(tmp, "v2"), recursive = TRUE); dir.create(file.path(tmp, "bin"))
  file.copy(file.path(root, "slurm"), file.path(tmp, "v2"), recursive = TRUE)
  rec <- file.path(tmp, "sbatch_calls.txt"); ctr <- file.path(tmp, "counter")
  writeLines(c("#!/bin/bash", sprintf("n=$(cat %s 2>/dev/null || echo 1000); n=$((n + 1)); echo $n > %s", ctr, ctr),
              sprintf("echo \"job=$n SWITCHES='${SWITCHES-<unset>}' args: $*\" >> %s", rec), "echo \"$n;cluster\""), file.path(tmp, "bin", "sbatch"))
  Sys.chmod(file.path(tmp, "bin", "sbatch"), "0755")
  run <- function(env, ...) { out <- suppressWarnings(system2("bash", c("v2/slurm/submit_main_grid.sh", ...), stdout = TRUE, stderr = TRUE,
                                                              env = c(sprintf("PATH=%s:%s", file.path(tmp, "bin"), Sys.getenv("PATH")), env)))
    list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  old <- setwd(tmp); on.exit(setwd(old), add = TRUE)
  ## (a) + (b): the lag-3 M1 array and the K = 2 M2 array depending on it (the existing mechanism)
  r <- run(c("M1_IDX=24-35", "M2_IDX=36-47")); expect_equal(r$status, 0L, info = paste(r$out, collapse = "\n"))
  calls <- readLines(rec); expect_length(calls, 2L)
  expect_match(calls[1], "^job=1001 SWITCHES='<unset>' args: --parsable --array=24-35 v2/slurm/02_main_grid.sbatch$")
  expect_match(calls[2], "^job=1002 SWITCHES='' args: --parsable --array=36-47 --dependency=afterok:1001 --kill-on-invalid-dep=yes v2/slurm/02_main_grid.sbatch$")
  expect_true(any(grepl("M1 array submitted: job 1001", r$out))); expect_true(any(grepl("M2 array submitted: job 1002 \\(indices 36-47\\), afterok:1001", r$out)))
  ## (c): the K = 1 M2 array on the SAME M1 job, no new M1 array, SWITCHES=--ame_K=1 exported to the M2 tasks only
  r <- run(c("M1_JOB=1001", "M2_IDX=36-47", "M2_SWITCHES=--ame_K=1"), "--account=open"); expect_equal(r$status, 0L, info = paste(r$out, collapse = "\n"))
  calls <- readLines(rec); expect_length(calls, 3L)
  expect_match(calls[3], "^job=1003 SWITCHES='--ame_K=1' args: --parsable --array=36-47 --dependency=afterok:1001 --kill-on-invalid-dep=yes --account=open v2/slurm/02_main_grid.sbatch$")
  expect_true(any(grepl("M1 array NOT submitted: the M2 array depends on the existing job 1001", r$out)))
  expect_true(any(grepl("M2 array submitted: job 1003 \\(indices 36-47, switches '--ame_K=1'\\), afterok:1001, kill-on-invalid-dep=yes", r$out)))
  ## refusals: a bad M1_JOB, --array in the arguments, K = 2 as M2_SWITCHES (that is the main grid)
  expect_true(run(c("M1_JOB=abc", "M2_IDX=36-47"))$status != 0L); expect_true(run(character(), "--array=0-5")$status != 0L)
  expect_true(run(c("M1_JOB=1001", "M2_SWITCHES=--ame_K=2"))$status != 0L); expect_length(readLines(rec), 3L)
  ## the sbatch appends SWITCHES to the fit call, so the K = 1 tasks produce the _K1 cells from the K-free M1 fits
  sb <- readLines(file.path(root, "slurm", "02_main_grid.sbatch"))
  expect_true(any(grepl("^SW=\\$\\{SWITCHES:-\\}", sb))); expect_true(any(grepl("^Rscript v2/R/02_fit_model.R .* \\$SW$", sb)))
  k1 <- parse_spec("--ame_K=1")
  for (s in MAIN_SAMPLES) for (m in c("M2m", "M2r", "M2a", "M2b")) {
    expect_equal(cell_name(s, 3L, m, k1), sprintf("%s_lag3_all_pr_%s_K1", s, m)); expect_false(cell_name(s, 3L, m, k1) == cell_name(s, 3L, m, parse_spec(character(0))))
    expect_equal(m1_cell_for_m2(s, 3L, m, k1), m1_cell_for_m2(s, 3L, m, parse_spec(character(0))))      # the same M1 fit for K = 1 and K = 2
  }
  ## 03_post.sbatch names the group's eight fits through cell_spec(), so SWITCHES="--ame_K=1" gives the four M1 fits and the four _K1 M2 fits
  post <- readLines(file.path(root, "slurm", "03_post.sbatch")); expect_true(any(grepl("cell_name(s, l, m, cell_spec(m, spec))", post, fixed = TRUE)))
  cells <- vapply(MODELS, function(m) cell_name("mid", 3L, m, cell_spec(m, k1)), "")
  expect_equal(unname(cells), c("mid_lag3_all_pr_M1m", "mid_lag3_all_pr_M1r", "mid_lag3_all_pr_M1a", "mid_lag3_all_pr_M1b", "mid_lag3_all_pr_M2m_K1", "mid_lag3_all_pr_M2r_K1", "mid_lag3_all_pr_M2a_K1", "mid_lag3_all_pr_M2b_K1"))
  expect_equal(compare_group_name("mid", 3L, k1), "mid_lag3_all_pr_K1")
})
