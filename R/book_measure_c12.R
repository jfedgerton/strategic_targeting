## R/book_measure_c12.R -- the book / JOP 2024 main measure built exactly to the PI's specification (2026-11-23), standalone (not yet a grid
## weighting). Output is in positions format (list(R = matrix, key)) plus edge tables, so the stage-1/stage-2 cell code can use it directly.
##   Layers (5): IGO (shared full memberships, min-max scaled to [0,1] within year, continuous), ATOP defense, ATOP offense, DCA (2010 carried
##   to 2011-14 as in the layer data), PTA (binary: book layer layer_pta_desta from a[5], see below). Years with no tie anywhere are structural zeros (score 0, via the year FE).
##   Expected tie per layer on ALL directed dyads (i,j) and (j,i): logit (binary) / Poisson PML (IGO), y ~ contig12 | sender + receiver + year,
##   contig12 = COW conttype 1 or 2 (land, or water <= 12 miles). Rows dropped as perfectly predicted by the FE must have y = 0; expected 0.
##   Layer score = observed - expected (no per-layer rescaling). Directed composite = mean of the 5 scores. Symmetrize: (c_ij + c_ji) / 2.
##   Within year: tie_norm = s / max|s| (signed, [-1,1]); tie_pos = max(s,0) / max positive s ([0,1]).
## Two measures (log PageRank, igraph damping 0.85, undirected, weighted; every state in the year's key is a vertex):
##   prpos   : PageRank on tie_pos (unsigned: below-expected ties carry no weight)
##   prsigned: PageRank on tie_norm min-max rescaled to [0,1] within year (signed information kept: below-expected ties are weak but
##                   ordered, the year's most negative pair gets 0). PageRank cannot take negative weights, so this is the signed variant.
## Spec ids: BOOKc12_ptadesta|<layers>|prpos (unsigned; the main measure, PI 2026-11-23) and |prsigned (signed; appendix), for the 5-layer set and
## each drop-one set. Outputs: a positions file (R, key) and an edges file (list(edges = per-spec edge tables)) read by s2_spec_edges.
## Usage: Rscript R/book_measure_c12.R LAYERS_RDS KEY_POSITIONS_RDS OUT_POSITIONS_RDS OUT_EDGES_RDS BOOK_LAYERS_RDS
suppressPackageStartupMessages({ library(data.table); library(fixest); library(igraph) })
a <- commandArgs(trailingOnly = TRUE); stopifnot(length(a) == 5)
DY <- as.data.table(readRDS(a[1])); KEY <- as.data.table(readRDS(a[2])$key)[, .(ccode, year)]
## PTA and DCA (PI 2026-11-24, book alignment): a[5] = book_layers.rds from R/book_layers.R. PTA = book layer_pta_desta (all members of an
## in-force DESTA treaty linked, actual withdrawals, predecessor mapping), replacing the dyadic "any PTA" coding; DCA = Kinne dcaAnyV2
## (2010 carried to 2011-14), replacing layers_dy dcaAnyV1.
PT <- as.data.table(readRDS(a[5]))[, .(ccode1, ccode2, year, pta_desta, dca_v2)]
DY[, dca := NULL]; setnames(PT, "dca_v2", "dca")
DY <- merge(DY, PT, by = c("ccode1", "ccode2", "year"), all.x = TRUE); stopifnot(!anyNA(DY$pta_desta), !anyNA(DY$dca))
cat(sprintf("PTA layer: %d linked dyad-years\n", sum(DY$pta_desta)))
setorder(DY, year, ccode1, ccode2); stopifnot(!anyDuplicated(DY[, .(ccode1, ccode2, year)]), all(DY$ccode1 < DY$ccode2))
DY[, `:=`(contig12 = as.integer(!is.na(conttype) & conttype %in% 1:2), pta_bin = as.integer(pta_desta > 0), dca_bin = as.integer(dca > 0))]
DY[, igo_mm := { r <- range(igo); if (r[2] > r[1]) (igo - r[1]) / (r[2] - r[1]) else 0 }, by = year]
stopifnot(!anyNA(DY[, .(igo_mm, def, off, dca_bin, pta_bin, contig12)]), all(DY$def %in% 0:1), all(DY$off %in% 0:1))
cat(sprintf("dyad-years %d; contig12 share %.4f; DCA years %d-%d; PTA years %d-%d\n", nrow(DY), mean(DY$contig12),
            min(DY$year[DY$dca_bin > 0]), max(DY$year[DY$dca_bin > 0]), min(DY$year[DY$pta_bin > 0]), max(DY$year[DY$pta_bin > 0])))
## directed rows: original orientation (i = ccode1) first, reversed second
n <- nrow(DY)
DD <- rbind(DY[, .(a = ccode1, b = ccode2, year, contig12, igo_mm, def, off, dca_bin, pta_bin)],
            DY[, .(a = ccode2, b = ccode1, year, contig12, igo_mm, def, off, dca_bin, pta_bin)])
LAY <- c(igo = "igo_mm", def = "def", off = "off", dca = "dca_bin", pta = "pta_bin")
SC <- sapply(names(LAY), function(l) {
  y <- DD[[LAY[[l]]]]; d <- data.frame(y = y, contig12 = DD$contig12, a = DD$a, b = DD$b, year = DD$year)
  w <- character(0)
  m <- withCallingHandlers(if (l == "igo") fepois(y ~ contig12 | a + b + year, data = d, notes = FALSE, glm.iter = 50)
                           else feglm(y ~ contig12 | a + b + year, data = d, family = binomial("logit"), notes = FALSE, glm.iter = 50),
                           warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
  if (length(w)) stop(l, ": fixest warning(s): ", paste(unique(w), collapse = " | "))
  if (!isTRUE(m$convStatus)) stop(l, ": GLM did not converge")
  ex <- rep(NA_real_, nrow(d)); ex[obs(m)] <- fitted(m); om <- is.na(ex)
  if (any(y[om] != 0)) stop(l, ": ", sum(y[om] != 0), " dropped rows have y != 0 (spec requires y = 0 for perfectly predicted rows)")
  ex[om] <- 0
  cat(sprintf("%-4s %-6s n %d | contig12 %.3f (%.3f) | dropped as perfectly predicted %d (all y = 0) | iterations %d | mean score %.4f\n",
              l, if (l == "igo") "PPML" else "logit", nrow(d), coef(m)[["contig12"]], se(m)[["contig12"]], sum(om), m$iterations, mean(y - ex)))
  y - ex })
## measures: the main 5-layer set and each drop-one set (for the bootstrapped cells), unsigned and signed
LAYSETS <- c(list(c("igo", "def", "off", "pta", "dca")), lapply(c("igo", "def", "off", "pta", "dca"), function(l) setdiff(c("igo", "def", "off", "pta", "dca"), l)))
yr_ix <- split(seq_len(n), DY$year)
pr_year <- function(w) { out <- vector("list", length(yr_ix)); names(out) <- names(yr_ix)
  for (y in names(yr_ix)) { ix <- yr_ix[[y]]; vs <- KEY[year == as.integer(y), ccode]; k <- ix[w[ix] > 0]
    g <- graph_from_data_frame(data.frame(as.character(DY$ccode1[k]), as.character(DY$ccode2[k]), weight = w[k]), directed = FALSE,
                               vertices = data.frame(name = as.character(vs)))
    out[[y]] <- data.table(ccode = vs, year = as.integer(y), R = log(page_rank(g, directed = FALSE)$vector[as.character(vs)])) }
  P <- rbindlist(out); v <- P$R[match(paste(KEY$ccode, KEY$year), paste(P$ccode, P$year))]; stopifnot(!anyNA(v)); v }
RL <- list(); EL <- list(); TIES <- NULL
for (L in LAYSETS) {
  comp <- rowMeans(SC[, L, drop = FALSE])                         # directed composite, mean of the layer scores
  s <- (comp[seq_len(n)] + comp[n + seq_len(n)]) / 2              # symmetrize
  tn <- tp <- tm <- numeric(n)
  for (ix in yr_ix) { m <- max(abs(s[ix])); if (m > 0) tn[ix] <- s[ix] / m
    p <- pmax(s[ix], 0); mp <- max(p); if (mp > 0) tp[ix] <- p / mp
    r <- range(tn[ix]); if (r[2] > r[1]) tm[ix] <- (tn[ix] - r[1]) / (r[2] - r[1]) }
  stopifnot(all(abs(tn) <= 1 + 1e-12), all(tp >= 0 & tp <= 1), all(tm >= 0 & tm <= 1))
  lab <- paste(L, collapse = "+")
  for (v in c("pos", "signed")) { sp <- sprintf("BOOKc12_ptadesta|%s|pr%s", lab, v); w <- if (v == "pos") tp else tm
    RL[[sp]] <- pr_year(w); EL[[sp]] <- data.table(a = DY$ccode1, b = DY$ccode2, year = DY$year, w = w)[w > 0] }
  if (length(L) == 5) { TIES <- DY[, .(ccode1, ccode2, year)][, `:=`(s = s, tie_norm = tn, tie_pos = tp, tie_mm = tm)]
    cat(sprintf("main set: share of dyad-years with signed score > 0 %.3f\n", mean(s > 0))) }
}
R <- do.call(cbind, RL); colnames(R) <- names(RL)
saveRDS(list(R = R, key = as.data.frame(KEY), ties = TIES), a[3])
saveRDS(list(edges = EL), a[4])
cat(sprintf("BOOK_MEASURE_DONE %s (%d specs) %s\n", a[3], ncol(R), a[4]))
