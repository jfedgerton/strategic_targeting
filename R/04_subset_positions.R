## ============================================================================
## 26_rel_subset_positions.R  --  build a positions file for a named subset of measures (PI request 2026-11-18, test 1:
## does the target-defensibility control absorb the chooser relational effect only because ATOP layers are inside the measure?)
## Reads the sweep's positions_*.rds files and writes ONE positions file in the format 22_rel_sweep.R reads (key, R, specs).
##   --set=noatop  every plausible measure with no ATOP layer (def/off/neu/nonagg/consul), plus the no-ATOP analogue of the pick
##                 (W2contig_ptaall|igo+dca+un|pr) first, kept even if it fails plausibility (flag in specs$plausible).
## Usage: Rscript v2/R/26_rel_subset_positions.R <positions_dir> <selection.rds> <out.rds> --set=noatop|ids [--ids=<file>]
## ============================================================================
suppressMessages(library(data.table))
set.seed(123)
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) >= 3)
PDIR <- args[1]; SELF <- args[2]; OUT <- args[3]
opt <- function(name, default) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
SET <- opt("set", ""); stopifnot(SET %in% c("noatop", "ids"))
logf <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")

## SELF: a selection .rds (list with P$spec_id = plausible ids) or the measure_validation.csv of 23_rel_validate.R (column plausible)
plaus <- if (grepl("\\.csv$", SELF)) { v <- read.csv(SELF, stringsAsFactors = FALSE); v$spec_id[v$plausible %in% c(TRUE, "TRUE")] } else readRDS(SELF)$P$spec_id
stopifnot(length(plaus) > 0)
layers_of <- function(s) sub("^[^|]*\\|([^|]*)\\|.*$", "\\1", s)
ATOP <- "(^|\\+)(def|off|neu|nonagg|consul)(\\+|$)"
ANALOGUE <- "W2contig_ptaall|igo+dca+un|pr"
if (SET == "noatop") {
  want <- c(ANALOGUE, sort(plaus[!grepl(ATOP, layers_of(plaus))]))
  want <- unique(want)
  stopifnot(!any(grepl(ATOP, layers_of(want))))
  logf("set %s: %d measures (%d plausible; analogue plausible = %s)", SET, length(want), sum(want %in% plaus), ANALOGUE %in% plaus)
} else {
  ## --set=ids --ids=<file>: exactly the listed spec ids, in file order (e.g. the Phase A pick + 200 sample rebuilt from new positions files)
  want <- readLines(opt("ids", "")); want <- want[nzchar(want)]; stopifnot(length(want) > 0, !anyDuplicated(want))
  logf("set ids: %d measures (%d plausible under %s)", length(want), sum(want %in% plaus), basename(SELF))
}

files <- list.files(PDIR, pattern = "^positions_.*\\.rds$", full.names = TRUE); stopifnot(length(files) > 0)
key <- NULL; cols <- list(); specs <- list()
for (f in files) {
  P <- readRDS(f); hit <- intersect(want, colnames(P$R)); if (!length(hit)) next
  k <- as.data.table(P$key)[, .(ccode, year)]
  if (is.null(key)) key <- k else stopifnot(identical(key, k))          # all positions files share one key (same country-year panel)
  for (s in hit) { if (!is.null(cols[[s]])) next; cols[[s]] <- P$R[, s]; specs[[s]] <- P$specs[P$specs$spec_id == s, , drop = FALSE] }
}
miss <- setdiff(want, names(cols)); if (length(miss)) stop(sprintf("measures not found in %s: %s", PDIR, paste(miss, collapse = ", ")))
R <- do.call(cbind, cols[want]); colnames(R) <- want
sp <- rbindlist(specs[want], fill = TRUE); sp$plausible <- sp$spec_id %in% plaus
stopifnot(identical(sp$spec_id, colnames(R)), nrow(R) == nrow(key), !anyNA(R))
saveRDS(list(key = as.data.frame(key), R = R, specs = as.data.frame(sp)), OUT)
logf("wrote %s: %d measures x %d country-years", OUT, ncol(R), nrow(R))
