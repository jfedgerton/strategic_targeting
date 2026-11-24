## R/book_measure_raw.R -- the book's UNWEIGHTED measure, rebuilt independently from the raw sources for replication (PI 2026-11-24).
## Book reference: code/04_multiverse/01_build_measure.R, LAYERS=custom LAYER_SET=igo+def+off+dca+pta IGO=dyadic WEIGHT=raw.
##   Layers (5), each in [0,1]:
##     IGO  shared full IGO memberships (COW IGO v3 count), min-max within year (raw count, NOT log1p as in W1raw)
##     def  ATOP defense (binary, symmetric); off ATOP offense (binary; ties end 1990); dca Kinne DCAD dcaAnyV2 binary (2010 carried to 2011-14)
##     pta  the book layer layer_pta_desta (R/book_layers.R): all members of an in-force DESTA treaty linked (ACP-ACP via Lome/Cotonou),
##          actual withdrawals only, predecessor mapping (Czechoslovakia, North Yemen)
##   True zeros everywhere (no NA). No expected-tie model, no residualization, no negatives.
##   Composite = mean of the 5 layers; symmetric by construction; tie_pos = composite / its within-year maximum.
##   Position = log PageRank (igraph, damping 0.85, undirected, weighted by tie_pos; every state in the year's key is a vertex).
## Spec ids: BOOKraw_ptadesta|<layers>|pr for the 5-layer set and each drop-one set.
## PI 2026-11-24 (book alignment): PTA = book layer_pta_desta and DCA = Kinne dcaAnyV2, both from R/book_layers.R (a[2]); IGO raw count min-max.
## Usage: Rscript R/book_measure_raw.R LAYERS_RDS BOOK_LAYERS_RDS KEY_POSITIONS_RDS OUT_POSITIONS_RDS OUT_EDGES_RDS
suppressPackageStartupMessages({ library(data.table); library(igraph) })
set.seed(123)
a <- commandArgs(trailingOnly = TRUE); stopifnot(length(a) == 5)
DY <- as.data.table(readRDS(a[1])); KEY <- as.data.table(readRDS(a[3])$key)[, .(ccode, year)]
setorder(DY, year, ccode1, ccode2); stopifnot(!anyDuplicated(DY[, .(ccode1, ccode2, year)]), all(DY$ccode1 < DY$ccode2))
stopifnot(nrow(DY) == 918315L)
PT <- as.data.table(readRDS(a[2]))[, .(ccode1, ccode2, year, pta_desta, dca_v2)]
DY <- merge(DY, PT, by = c("ccode1", "ccode2", "year"), all.x = TRUE); stopifnot(!anyNA(DY$pta_desta), !anyNA(DY$dca_v2))
setorder(DY, year, ccode1, ccode2)
cat(sprintf("PTA layer %d linked dyad-years; DCA v2 %d\n", sum(DY$pta_desta), sum(DY$dca_v2)))

## ---- other layers ----
DY[, igo_mm := { r <- range(igo); if (r[2] > r[1]) (igo - r[1]) / (r[2] - r[1]) else 0 }, by = year]
DY[, dca_bin := dca_v2]
LAY <- c(igo = "igo_mm", def = "def", off = "off", pta = "pta_desta", dca = "dca_bin")
stopifnot(!anyNA(DY[, LAY, with = FALSE]), all(DY$def %in% 0:1), all(DY$off %in% 0:1), all(DY$igo_mm >= 0 & DY$igo_mm <= 1))
stopifnot(max(DY$year[DY$off > 0]) == 1990L, min(DY$year[DY$dca_bin > 0]) == 1980L,
          all(DY[year %in% 2011:2014, sum(dca_bin), by = year]$V1 == DY[year == 2010, sum(dca_bin)]))

## ---- composite, rescale, PageRank ----
n <- nrow(DY); yr_ix <- split(seq_len(n), DY$year)
pr_year <- function(w) { out <- vector("list", length(yr_ix)); names(out) <- names(yr_ix)
  for (y in names(yr_ix)) { ix <- yr_ix[[y]]; vs <- KEY[year == as.integer(y), ccode]; k <- ix[w[ix] > 0]
    g <- graph_from_data_frame(data.frame(as.character(DY$ccode1[k]), as.character(DY$ccode2[k]), weight = w[k]), directed = FALSE,
                               vertices = data.frame(name = as.character(vs)))
    out[[y]] <- data.table(ccode = vs, year = as.integer(y), R = log(page_rank(g, directed = FALSE)$vector[as.character(vs)])) }
  P <- rbindlist(out); v <- P$R[match(paste(KEY$ccode, KEY$year), paste(P$ccode, P$year))]; stopifnot(!anyNA(v)); v }
L5 <- c("igo", "def", "off", "pta", "dca")
LAYSETS <- c(list(L5), lapply(L5, function(l) setdiff(L5, l)))
RL <- list(); EL <- list(); TIES <- NULL
for (L in LAYSETS) {
  comp <- rowMeans(as.matrix(DY[, LAY[L], with = FALSE]))
  tp <- numeric(n); for (ix in yr_ix) { m <- max(comp[ix]); if (m > 0) tp[ix] <- comp[ix] / m }
  stopifnot(all(tp >= 0 & tp <= 1))
  sp <- sprintf("BOOKraw_ptadesta|%s|pr", paste(L, collapse = "+"))
  RL[[sp]] <- pr_year(tp); EL[[sp]] <- data.table(a = DY$ccode1, b = DY$ccode2, year = DY$year, w = tp)[w > 0]
  if (length(L) == 5) { TIES <- DY[, .(ccode1, ccode2, year, tie_signed = comp, tie_pos = tp)]
    cat(sprintf("main set: share of dyad-years with a positive tie, pooled %.3f (book build 0.952); mean of yearly shares %.3f (book 0.701)\n",
                mean(tp > 0), mean(tapply(tp > 0, DY$year, mean))))
    stopifnot(abs(mean(tp > 0) - 0.952) < 0.001) }
}
R <- do.call(cbind, RL); colnames(R) <- names(RL)
saveRDS(list(R = R, key = as.data.frame(KEY), ties = TIES), a[4])
saveRDS(list(edges = EL), a[5])
cat(sprintf("BOOK_RAW_DONE %s (%d specs) %s\n", a[4], ncol(R), a[5]))
