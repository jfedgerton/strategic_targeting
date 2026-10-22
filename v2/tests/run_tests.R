## Run all v2 unit tests:  Rscript v2/tests/run_tests.R
suppressPackageStartupMessages(library(testthat))
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(this_file), ".."))
Sys.setenv(V2_ROOT = root)
res <- test_dir(file.path(root, "tests", "testthat"), reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("\n%d tests: %d passed, %d failed, %d skipped\n", nrow(df), sum(df$passed), sum(df$failed), sum(df$skipped)))
if (sum(df$failed) > 0) quit(status = 1)
