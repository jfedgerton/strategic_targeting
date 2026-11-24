## R/book_measure_raw.R -- the book's UNWEIGHTED measure, rebuilt independently from the raw sources for replication (PI 2026-11-24).
## Book reference: code/04_multiverse/01_build_measure.R, LAYERS=custom LAYER_SET=igo+def+off+dca+pta IGO=dyadic WEIGHT=raw.
##   Layers (5), each in [0,1]:
##     IGO  shared full IGO memberships (COW IGO v3 count), min-max within year (raw count, NOT log1p as in W1raw)
##     def  ATOP defense (binary, symmetric); off ATOP offense (binary; ties end 1990); dca Kinne DCAD binary (2010 carried to 2011-14)
##     pta  binary joint membership in an IN-FORCE PTA, DESTA, with every plurilateral treaty treated as a CLIQUE over its members, so
##          ACP-ACP pairs are linked through Yaounde / Lome / Cotonou (DESTA's dyadic file codes those as EC-ACP pairs only).
##          A member's join year = its earliest entry-into-force year in that base treaty (accessions included); a pair's start = the
##          later of the two join years. Treaties without an entry-into-force year are not counted. No exit years (as in the layer data).
##   True zeros everywhere (no NA). No expected-tie model, no residualization, no negatives.
##   Composite = mean of the 5 layers; symmetric by construction; tie_pos = composite / its within-year maximum.
##   Position = log PageRank (igraph, damping 0.85, undirected, weighted by tie_pos; every state in the year's key is a vertex).
## Spec ids: BOOKraw_ptaclq|<layers>|pr for the 5-layer set and each drop-one set.
## Usage: Rscript R/book_measure_raw.R LAYERS_RDS DESTA_RDA KEY_POSITIONS_RDS OUT_POSITIONS_RDS OUT_EDGES_RDS
suppressPackageStartupMessages({ library(data.table); library(igraph); library(countrycode) })
set.seed(123)
a <- commandArgs(trailingOnly = TRUE); stopifnot(length(a) == 5)
DY <- as.data.table(readRDS(a[1])); KEY <- as.data.table(readRDS(a[3])$key)[, .(ccode, year)]
setorder(DY, year, ccode1, ccode2); stopifnot(!anyDuplicated(DY[, .(ccode1, ccode2, year)]), all(DY$ccode1 < DY$ccode2))
stopifnot(nrow(DY) == 918315L)

## ---- PTA: clique over members of each in-force DESTA base treaty ----
e <- new.env(); suppressWarnings(load(a[2], envir = e)); ds <- as.data.table(e$desta)
ds <- ds[!is.na(entryforceyear)]
cw <- function(iso, nm) { x <- countrycode(iso, "iso3n", "cown", warn = FALSE); x[nm == "Serbia"] <- 345L; x[nm == "Kosovo"] <- 347L; as.integer(x) }
MEM <- rbind(ds[, .(bt = base_treaty, c = cw(iso1, country1), nm = country1, start = entryforceyear)],
             ds[, .(bt = base_treaty, c = cw(iso2, country2), nm = country2, start = entryforceyear)])
cat(sprintf("DESTA in-force rows %d; member rows with no COW code (dropped) %d: %s\n", nrow(ds), sum(is.na(MEM$c)),
            paste(sort(unique(MEM[is.na(c), nm])), collapse = "; ")))
MEM <- MEM[!is.na(c)]
succ <- c(`255` = 260L, `316` = 315L, `679` = 678L)            # unified Germany / Czechia / Yemen -> COW predecessors (as R/01_build_layers.R)
MEM <- rbind(MEM, MEM[as.character(c) %in% names(succ)][, c := succ[as.character(c)]])
MEM <- MEM[, .(join = min(start)), by = .(bt, c)]
PR <- MEM[MEM, on = "bt", allow.cartesian = TRUE, nomatch = 0L][c < i.c, .(ccode1 = c, ccode2 = i.c, start = pmax(join, i.join))]
PR <- PR[, .(start = min(start)), by = .(ccode1, ccode2)]
DY <- PR[DY, on = c("ccode1", "ccode2")]
DY[, pta_clq := as.integer(!is.na(start) & start <= year)][, start := NULL]
setorder(DY, year, ccode1, ccode2)
## the dyadic in-force coding must be a subset of the clique coding (every dyadic treaty pair is a clique pair)
stopifnot(all(DY$pta_clq[DY$pta_inforce > 0] == 1L))
acp <- c(51L, 452L, 433L, 950L, 501L)                            # Jamaica, Ghana, Senegal, Fiji, Kenya
chk <- DY[year == 1995 & ccode1 %in% acp & ccode2 %in% acp, .(n = .N, dyadic = sum(pta_inforce > 0), clique = sum(pta_clq))]
cat(sprintf("ACP check 1995 (Jamaica, Ghana, Senegal, Fiji, Kenya; %d pairs): linked dyadic %d, clique %d\n", chk$n, chk$dyadic, chk$clique))
stopifnot(chk$clique == chk$n)
cat(sprintf("PTA dyad-years linked: dyadic in force %d, clique %d (+%d); clique years %d-%d\n", sum(DY$pta_inforce > 0), sum(DY$pta_clq),
            sum(DY$pta_clq) - sum(DY$pta_inforce > 0), min(DY$year[DY$pta_clq > 0]), max(DY$year[DY$pta_clq > 0])))

## ---- other layers ----
DY[, igo_mm := { r <- range(igo); if (r[2] > r[1]) (igo - r[1]) / (r[2] - r[1]) else 0 }, by = year]
DY[, dca_bin := as.integer(dca > 0)]
LAY <- c(igo = "igo_mm", def = "def", off = "off", pta = "pta_clq", dca = "dca_bin")
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
  sp <- sprintf("BOOKraw_ptaclq|%s|pr", paste(L, collapse = "+"))
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
