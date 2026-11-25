## R/grid_build.R -- positions (log PageRank) for the 236 grid measures (R/lib/grid.R), one positions file per block, in the format the
## stage-1/2 cell code reads (list(R = state-year x spec matrix, key, specs)). Edges are not stored: s2_spec_edges rebuilds them from
## grid_layers.rds with the same function (grid_edges), so positions and stage-2 ties always use identical edges.
## Check: the raw-block spec with the main layers, count IGO and book PTA must equal the main measure (BOOKraw_ptadesta) to 1e-10.
##        With --c12=FILE (the weighted book measure positions) the sry contig12 main-layer spec must equal BOOKc12_ptadesta (prpos) to 1e-8.
## Usage: Rscript R/grid_build.R GRID_LAYERS_RDS KEY_POSITIONS_RDS MAIN_POSITIONS_RDS OUT_DIR [--cores=N] [--c12=FILE]
suppressPackageStartupMessages({ library(data.table); library(igraph); library(parallel) })
set.seed(123)
source("R/lib/grid.R")
a <- commandArgs(trailingOnly = TRUE); stopifnot(length(a) >= 4)
CORES <- as.integer(sub("--cores=", "", grep("^--cores=", a, value = TRUE)[1])); if (is.na(CORES)) CORES <- 2L
GLF <- a[1]; GL <- readRDS(GLF); setDT(GL); setorder(GL, year, ccode1, ccode2)
KEY <- as.data.table(readRDS(a[2])$key)[, .(ccode, year)]; OUTD <- a[4]; dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
SP <- grid_specs()
yr_ix <- split(seq_len(nrow(GL)), GL$year)
pr_year <- function(w) { out <- vector("list", length(yr_ix)); names(out) <- names(yr_ix)
  for (y in names(yr_ix)) { ix <- yr_ix[[y]]; vs <- KEY[year == as.integer(y), ccode]; k <- ix[w[ix] > 0]
    g <- graph_from_data_frame(data.frame(as.character(GL$ccode1[k]), as.character(GL$ccode2[k]), weight = w[k]), directed = FALSE,
                               vertices = data.frame(name = as.character(vs)))
    out[[y]] <- data.table(ccode = vs, year = as.integer(y), R = log(page_rank(g, directed = FALSE)$vector[as.character(vs)])) }
  P <- rbindlist(out); v <- P$R[match(paste(KEY$ccode, KEY$year), paste(P$ccode, P$year))]; stopifnot(!anyNA(v), all(is.finite(v))); v }
t0 <- Sys.time()
RV <- mclapply(SP$spec_id, function(sp) pr_year(grid_weight(GL, sp)), mc.cores = CORES, mc.preschedule = FALSE)
bad <- which(!vapply(RV, is.numeric, TRUE)); if (length(bad)) stop("PageRank failed for ", paste(SP$spec_id[bad], collapse = ", "), ": ", as.character(RV[[bad[1]]]))
R <- do.call(cbind, RV); colnames(R) <- SP$spec_id
cat(sprintf("%d specs in %.1f min\n", ncol(R), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
## equivalence with the main measure
M <- readRDS(a[3]); stopifnot(identical(as.integer(M$key$ccode), as.integer(KEY$ccode)), identical(as.integer(M$key$year), as.integer(KEY$year)))
dmax <- max(abs(R[, "Graw_count_book|igo+def+off+dca+pta|pr"] - M$R[, "BOOKraw_ptadesta|igo+def+off+pta+dca|pr"]))
cat(sprintf("grid raw count/book main layers vs BOOKraw_ptadesta: max |diff| %.2e\n", dmax)); stopifnot(dmax < 1e-10)
if (length(grep("^--c12=", a))) {   # run-stopping check listed in docs/PAPER_PIPELINE.md (audit 2026-11-24)
  M2 <- readRDS(sub("^--c12=", "", grep("^--c12=", a, value = TRUE)[1])); stopifnot(identical(as.integer(M2$key$ccode), as.integer(KEY$ccode)), identical(as.integer(M2$key$year), as.integer(KEY$year)))
  d2 <- max(abs(R[, "Gsryc12_count_book|igo+def+off+dca+pta|pr"] - M2$R[, "BOOKc12_ptadesta|igo+def+off+pta+dca|prpos"]))
  cat(sprintf("grid sry contig12 main layers vs BOOKc12_ptadesta prpos: max |diff| %.2e\n", d2)); stopifnot(d2 < 1e-8) }
## constant positions in a year (would be refused by the cells) are reported here
cst <- vapply(seq_len(ncol(R)), function(j) sum(tapply(R[, j], KEY$year, function(z) sd(z) < 1e-8), na.rm = TRUE), 0)
if (any(cst > 0)) cat("specs with constant-position years:", paste(sprintf("%s (%d)", SP$spec_id[cst > 0], cst[cst > 0]), collapse = "; "), "\n")
for (b in unique(SP$block)) { S <- SP[block == b]
  saveRDS(list(R = R[, S$spec_id, drop = FALSE], key = as.data.frame(KEY), specs = as.data.frame(S), built = format(Sys.time()),
               provenance = grid_provenance(GLF)), file.path(OUTD, sprintf("positions_core_G%s.rds", b))) }   # GPT P2-5: positions bound to the code and layers that built them
fwrite(SP, file.path(OUTD, "grid_specs.csv"))
cat(sprintf("GRID_BUILD_DONE %s: %s\n", OUTD, paste(sprintf("%s %d", names(table(SP$block)), table(SP$block)), collapse = ", ")))
