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
})
