## ============================================================================
## compile_models.R -- compile both Stan programs once and cache them as .rds
## Usage: Rscript v2/R/compile_models.R
## Output: v2/stan/build/{m1_latent_additive,m2_ame}.rds  (gitignored)
## ============================================================================
suppressPackageStartupMessages(library(rstan))
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "00_helpers.R"))
root <- v2_root()
dir.create(file.path(root, "stan", "build"), showWarnings = FALSE, recursive = TRUE)
for (m in c("m1_latent_additive", "m2_ame")) {
  src <- file.path(root, "stan", paste0(m, ".stan"))
  out <- file.path(root, "stan", "build", paste0(m, ".rds"))
  if (file.exists(out) && file.mtime(out) > file.mtime(src)) { cat(m, ": cached\n"); next }
  t0 <- Sys.time()
  sm <- stan_model(src, model_name = m, auto_write = FALSE)
  saveRDS(sm, out)
  cat(sprintf("%s compiled in %.0f s -> %s\n", m, as.numeric(difftime(Sys.time(), t0, units = "secs")), out))
}
## get_model(): load a compiled model, compiling if missing
