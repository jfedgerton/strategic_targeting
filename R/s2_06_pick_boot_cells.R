## R/s2_06_pick_boot_cells.R -- the bootstrapped sub-multiverse of the stage-2 appendix (PI 2026-11-22, layer 2): a fixed, rule-based list of
## cells chosen before any result is seen, each refitted by the two-stage dispute bootstrap (R/s2_07_boot_cells.R).
##   (1) the main measure at every window x sample (fmid, mid)
##   (2) main window, fmid: the main measure's layer set and centrality in every other positions file (weighting x PTA rule x arm); where a
##       file lacks that exact layer set, the same-centrality spec with the closest layer set (Jaccard; ties broken alphabetically by spec id); the IGO arm: its specs with the main centrality
##   (3) main window, fmid, main file: the main layer set under every other centrality
##   (4) main window, fmid, main file, main centrality: the main layer set minus one layer (each layer in turn)
##   (5) main window, fmid: main layers + pta under each PTA coding (all / inforce / any), main weighting and centrality
## Centrality is crossed only in the main file (rule 3); other files and the IGO arm use the main centrality. A selected sub-multiverse, not a factorial sample.
## Draws: B_MAIN for the main cell, B_CELL for the rest.
## Usage: Rscript R/s2_06_pick_boot_cells.R --posdir=DIR --main_spec=SPEC --main_window=1 --windows="0 1 2" --b_main=2000 --b_cell=500 --out=cells.csv
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
PD <- opt("posdir"); MAIN <- opt("main_spec"); WM <- as.integer(opt("main_window", "1")); WS <- as.integer(strsplit(opt("windows", "0 1 2 3 4 5 6"), " +")[[1]])
BM <- as.integer(opt("b_main", "2000")); BC <- as.integer(opt("b_cell", "500")); OUTF <- opt("out")
stopifnot(!is.null(PD), !is.null(MAIN), !is.null(OUTF))
pf <- sort(list.files(PD, pattern = "^positions_.*\\.rds$", full.names = TRUE)); pf <- pf[!grepl("^positions_cmp", basename(pf))]
SP <- rbindlist(lapply(pf, function(f) data.table(file = basename(f), spec_id = colnames(readRDS(f)$R))))
SP[, `:=`(p1 = sub("\\|.*$", "", spec_id), lay = sub("^[^|]*\\|([^|]*)\\|.*$", "\\1", spec_id), cent = sub("^.*\\|", "", spec_id))]
mp <- strsplit(MAIN, "|", fixed = TRUE)[[1]]; ML <- strsplit(mp[2], "+", fixed = TRUE)[[1]]; MC <- mp[3]
MF <- SP[spec_id == MAIN & grepl("^positions_core_", file), file]; if (length(MF) != 1) stop("main spec not found exactly once among core positions files")
jac <- function(l) { s <- strsplit(l, "+", fixed = TRUE)[[1]]; length(intersect(s, ML)) / length(union(s, ML)) }
C <- list()
add <- function(file, spec, W, smp, rule) C[[length(C) + 1L]] <<- data.table(file = file, spec_id = spec, window = W, sample = smp, rule = rule)
for (W in WS) for (s in c("fmid", "mid")) add(MF, MAIN, W, s, "1 main measure, window x sample")
for (f in setdiff(unique(SP$file), MF)) { S <- SP[file == f]
  if (grepl("igo_arm", f)) { for (sp in S[cent == MC, spec_id]) add(f, sp, WM, "fmid", "2 IGO arm, main centrality"); next }
  MCf <- if (MC %in% S$cent) MC else "pr"   # the book measure's centralities (prpos / prsigned) exist only in its own file; others use PageRank
  S <- S[cent == MCf]; if (!nrow(S)) next
  S[, j := vapply(lay, jac, 0)]; add(f, S[order(-j, spec_id)][1, spec_id], WM, "fmid", if (max(S$j) == 1) "2 main layers, other file" else "2 closest layers, other file") }
for (sp in SP[file == MF & cent != MC & vapply(lay, function(x) setequal(strsplit(x, "+", fixed = TRUE)[[1]], ML), TRUE), spec_id]) add(MF, sp, WM, "fmid", "3 main layers, other centrality")
MW0 <- sub("_pta[a-z]+$", "", mp[1])
for (l in ML) { keep <- setdiff(ML, l)
  ## a set without PTA lives only in the weighting's ptaall file (in force / any files hold PTA combinations only)
  hit <- function(f) SP[file == f & cent == MC & vapply(lay, function(x) setequal(strsplit(x, "+", fixed = TRUE)[[1]], keep), TRUE), spec_id]
  f4 <- MF; sp <- hit(f4)   # the main file first (the book measure file holds its own drop-one sets)
  if (!length(sp) && !("pta" %in% keep)) { f4 <- sprintf("positions_core_%s_ptaall.rds", MW0); sp <- hit(f4) }
  if (length(sp)) add(f4, sp[1], WM, "fmid", paste0("4 drop layer ", l)) else cat(sprintf("rule 4: no spec for main layers minus %s in %s\n", l, f4)) }
## (5) GPT P2-6, matched PTA comparison: the main layers + pta, main weighting and centrality, under each PTA coding (all / inforce / any), so
##     the PTA comparison changes only the PTA coding, not the layer set
MW <- sub("_pta[a-z]+$", "", mp[1])
for (pt in c("ptaall", "ptainforce", "ptaany")) { f <- sprintf("positions_core_%s_%s.rds", MW, pt)
  sp <- SP[file == f & cent == MC & vapply(lay, function(x) setequal(strsplit(x, "+", fixed = TRUE)[[1]], c(ML, "pta")), TRUE), spec_id]
  if (length(sp)) add(f, sp[1], WM, "fmid", paste0("5 matched PTA: main layers + pta, ", pt)) else cat(sprintf("rule 5: no main layers + pta spec in %s\n", f)) }
C <- unique(rbindlist(C), by = c("file", "spec_id", "window", "sample"))
C[, n_boot := ifelse(spec_id == MAIN & file == MF & window == WM & sample == "fmid", BM, BC)]
C[, cell_id := sprintf("b%03d", seq_len(.N))]
fwrite(C, OUTF)
cat(sprintf("wrote %s: %d cells (%s); total draws %d\n", OUTF, nrow(C), paste(sprintf("%s %d", names(table(substr(C$rule, 1, 1))), table(substr(C$rule, 1, 1))), collapse = ", "), sum(C$n_boot)))
