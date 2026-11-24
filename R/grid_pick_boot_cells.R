## R/grid_pick_boot_cells.R -- the bootstrapped cells of the network-measure grid (PI 2026-11-24), a fixed rule-based list chosen before any
## result is seen; each is refitted by the two-stage dispute bootstrap (R/s2_07_boot_cells.R).
##   (1) the main measure (BOOKraw_ptadesta, unweighted book measure) at every window x sample (fmid, mid) except the main window (main table cells m_fmid / m_mid)
##   (2) main window, fmid: main layers under the other IGO codings (cos, hyper; book PTA) and under PTA depth (count IGO)
##   (3) main window, fmid: the main layers minus one layer (raw, count IGO, book PTA where present), each layer in turn
##   (4) main window, fmid: the main layers under each weighting scheme (sry contig12, sry distance, sy contig12, sy distance)
##   (5) main window, fmid: main layers with DCA v1; main layers + the other-alliances layer
## Draws: B_CELL for every cell here (the main cells at the main window get B_MAIN in the main table).
## Usage: Rscript R/grid_pick_boot_cells.R --posdir=DIR --main_spec=SPEC --main_window=1 --windows="0 1 2 3 4 5 6" --b_main=2000 --b_cell=500 --out=cells.csv
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
PD <- opt("posdir"); MAIN <- opt("main_spec"); WM <- as.integer(opt("main_window", "1")); WS <- as.integer(strsplit(opt("windows", "0 1 2 3 4 5 6"), " +")[[1]])
BM <- as.integer(opt("b_main", "2000")); BC <- as.integer(opt("b_cell", "500")); OUTF <- opt("out"); stopifnot(!is.null(PD), !is.null(MAIN), !is.null(OUTF), WM %in% WS)
MF <- sprintf("positions_core_%s.rds", sub("\\|.*$", "", MAIN)); stopifnot(file.exists(file.path(PD, MF)), MAIN %in% colnames(readRDS(file.path(PD, MF))$R))
have <- function(f, sp) { p <- file.path(PD, f); if (!file.exists(p)) stop("missing ", f); if (!sp %in% colnames(readRDS(p)$R)) stop("missing spec ", sp, " in ", f); TRUE }
C <- list(); add <- function(f, sp, W, s, rule) { have(f, sp); C[[length(C) + 1L]] <<- data.table(file = f, spec_id = sp, window = W, sample = s, rule = rule) }
for (W in setdiff(WS, WM)) for (s in c("fmid", "mid")) add(MF, MAIN, W, s, "1 main measure, window x sample")   # main window x {fmid, mid} is bootstrapped once, by the main table (step 3, B_MAIN), not repeated here (audit 2026-11-24)
ML <- "igo+def+off+dca+pta"
add("positions_core_Graw.rds", sprintf("Graw_cos_book|%s|pr", ML), WM, "fmid", "2 IGO cosine")
add("positions_core_Graw.rds", sprintf("Graw_hyper_book|%s|pr", ML), WM, "fmid", "2 IGO hyperedge")
add("positions_core_Graw.rds", sprintf("Graw_count_depth|%s|pr", ML), WM, "fmid", "2 PTA depth")
for (l in c("igo", "def", "off", "dca", "pta")) { L <- setdiff(c("igo", "def", "off", "dca", "pta"), l)
  sp <- sprintf("Graw_%s_%s|%s|pr", if ("igo" %in% L) "count" else "na", if ("pta" %in% L) "book" else "na", paste(L, collapse = "+"))
  add("positions_core_Graw.rds", sp, WM, "fmid", paste0("3 drop layer ", l)) }
for (b in c("sryc12", "srydist", "syc12", "sydist")) add(sprintf("positions_core_G%s.rds", b), sprintf("G%s_count_book|%s|pr", b, ML), WM, "fmid", paste0("4 weighting ", b))
add("positions_core_Gdca1.rds", sprintf("Gdca1_count_book|%s|pr", ML), WM, "fmid", "5 DCA v1")
add("positions_core_Gally.rds", sprintf("Gally_count_book|%s+ally|pr", ML), WM, "fmid", "5 + other alliances")
C <- rbindlist(C); stopifnot(!anyDuplicated(C[, .(file, spec_id, window, sample)]))
C[, n_boot := BC][, cell_id := sprintf("g%03d", seq_len(.N))]
fwrite(C, OUTF)
cat(sprintf("wrote %s: %d cells; total draws %d\n", OUTF, nrow(C), sum(C$n_boot)))
