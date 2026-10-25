## ============================================================================
## compile_models.R -- compile both Stan programs once and cache them as .rds
## Usage: Rscript v2/R/compile_models.R
## Output: v2/stan/build/{m1_latent_additive,m2_ame}.rds  (gitignored) plus <model>.sha256, the sha256 of the .stan source
##         compiled into it (provenance, 2026-10-25 C3); 02_fit_model.R refuses a binary whose model_code is not the checked-out source
## ============================================================================
suppressPackageStartupMessages(library(rstan))
## Threading: rstan 2.32 puts -DSTAN_THREADS in its own CPPFLAGS (rstan:::PKG_CPPFLAGS_env_fun), so every
## model compiled here supports reduce_sum multi-threading via rstan_options(threads_per_chain = k).
## Assert it, so a different rstan build on Roar fails loudly instead of silently running serial.
cppflags <- tryCatch(rstan:::PKG_CPPFLAGS_env_fun(), error = function(e) "")
if (!grepl("STAN_THREADS", cppflags) && !grepl("STAN_THREADS", Sys.getenv("PKG_CPPFLAGS")))
  stop("rstan's CPPFLAGS lack -DSTAN_THREADS; export PKG_CPPFLAGS=-DSTAN_THREADS (see v2/slurm/00_env.sh) and recompile")
cat("compiling with -DSTAN_THREADS (rstan ", as.character(packageVersion("rstan")), ")\n", sep = "")
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "00_helpers.R"))
root <- v2_root()
dir.create(file.path(root, "stan", "build"), showWarnings = FALSE, recursive = TRUE)
for (m in c("m1_latent_additive", "m2_ame")) {
  src <- file.path(root, "stan", paste0(m, ".stan"))
  out <- file.path(root, "stan", "build", paste0(m, ".rds"))
  side <- file.path(root, "stan", "build", paste0(m, ".sha256"))
  src_sha <- stan_source_sha256(src)
  ## cached only when the binary is newer than the source AND its sidecar hash equals the source hash (provenance, 2026-10-25 C3)
  if (file.exists(out) && file.mtime(out) > file.mtime(src) && file.exists(side) && identical(trimws(readLines(side, warn = FALSE)[1]), src_sha)) {
    assert_compiled_current(readRDS(out), src, side); cat(m, ": cached (source sha256 ", src_sha, ")\n", sep = ""); next
  }
  t0 <- Sys.time()
  sm <- stan_model(src, model_name = m, auto_write = FALSE)
  saveRDS(sm, out)
  writeLines(src_sha, side)
  prov <- assert_compiled_current(sm, src, side)         # the compiled model_code must be the source just read
  cat(sprintf("%s compiled in %.0f s -> %s (source sha256 %s; sidecar %s)\n", m, as.numeric(difftime(Sys.time(), t0, units = "secs")), out, prov$stan_source_sha256, basename(side)))
}
## get_model(): load a compiled model, compiling if missing
