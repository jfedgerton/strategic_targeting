## R/grid_layers.R -- every layer and every weighted layer score the network-measure grid needs (PI 2026-11-24), on the undirected frame.
## The grid itself (which layers, codings and weightings combine into each measure) is built by R/grid_build.R from this file.
##
## Raw layers, each in [0,1] (true zeros, no NA):
##   igo_count, igo_cos, igo_hyper : COW IGO v3 full memberships (code 1), IGOs with >= 2 full members (book 01_data/06_igo_variants.R)
##        count = the book / main dyadic layer (layers_dy igo); cos = count / sqrt(d_i d_j); hyper = sum_k B_ik B_jk / (n_k - 1); each min-max scaled within year
##   def, off      : ATOP defense / offense pledge, either direction (layers_dy; equals the book atop_un_layers)
##   othally       : 1 if a neutrality, non-aggression or consultation pledge links the pair (max of the three ATOP types; layers_dy)
##   dca2, dca1    : Kinne DCAD dcaAnyV2 (book main; book_layers.rds) and dcaAnyV1 (layers_dy), 2010 carried to 2011-14
##   pta_book      : book layer_pta_desta (book_layers.rds; 196,196 linked)
##   pta_depth     : max DESTA 2.3 depth_index over in-force agreements / 7 (book 01_data/04_pta_desta_costly.R: dyadic DESTA 2.3,
##                   in force from entryforceyear until the year before a recorded withdrawal (withdrawal year = its entry-into-force year, else its signature year); accessions take the base treaty's depth),
##                   with the predecessor mapping used for pta_book (Germany 260 through 1990, Czechoslovakia 315 through 1992,
##                   North Yemen 678 through 1989)
## Weighted layer scores (appendix weighting schemes; book 04_multiverse/01_build_measure.R), for igo (count), def, off, dca2, pta_book:
##   sry_c12 / sry_dist : y - yhat, yhat from feglm logit (binary) / fepois (IGO) on ALL directed dyads, y ~ GEO | ccode1 + ccode2 + year
##   sy_c12  / sy_dist  : same with FE ccode1^year + ccode2^year
##   GEO = contig12 (COW conttype 1 or 2) or ldist = log(1 + capital distance). Each layer is fitted on the years in which it varies;
##   other years score 0. Observations dropped as perfectly predicted: sry -> yhat = 0 (must have y = 0, else stop); sy -> yhat = y (score 0).
##   Stored per layer as the symmetrized score (s_ij + s_ji) / 2 (the composite is a mean over layers, so symmetrizing per layer first is identical).
## Usage: Rscript R/grid_layers.R LAYERS_DY_RDS BOOK_LAYERS_RDS DATA_BOOK_DIR OUT_RDS [--cores=N]
suppressPackageStartupMessages({ library(data.table); library(fixest); library(countrycode) })
set.seed(123)
a <- commandArgs(trailingOnly = TRUE); stopifnot(length(a) >= 4)
CORES <- as.integer(sub("--cores=", "", grep("^--cores=", a, value = TRUE)[1])); if (is.na(CORES)) CORES <- 2L; setFixest_nthreads(CORES)
DY <- as.data.table(readRDS(a[1]))[, .(ccode1, ccode2, year, def, off, neu, nonagg, consul, dca1 = as.integer(dca > 0), conttype, capdist)]
setorder(DY, year, ccode1, ccode2); stopifnot(nrow(DY) == 918315L, !anyDuplicated(DY[, .(ccode1, ccode2, year)]), all(DY$ccode1 < DY$ccode2))
BL <- as.data.table(readRDS(a[2]))[, .(ccode1, ccode2, year, pta_book = pta_desta, dca2 = dca_v2)]
DY <- merge(DY, BL, by = c("ccode1", "ccode2", "year")); stopifnot(nrow(DY) == 918315L)
DB <- a[3]; OUT <- a[4]
DY[, othally := as.integer(neu > 0 | nonagg > 0 | consul > 0)][, c("neu", "nonagg", "consul") := NULL]
DY[, `:=`(contig12 = as.integer(!is.na(conttype) & conttype %in% 1:2), ldist = log1p(capdist))]
stopifnot(!anyNA(DY[, .(def, off, othally, dca1, dca2, pta_book, contig12, ldist)]))
mm <- function(x, yr) { o <- numeric(length(x)); for (ix in split(seq_along(x), yr)) { r <- range(x[ix]); if (r[2] > r[1]) o[ix] <- (x[ix] - r[1]) / (r[2] - r[1]) }; o }

## ---- IGO codings (book 06_igo_variants.R) ----
w <- fread(file.path(DB, "cow_igo_state_year_formatv3.csv.gz")); cols <- setdiff(names(w), c("ccode", "year", "state"))
ig <- list()
for (yr in sort(unique(DY$year))) { x <- w[year == yr]; if (!nrow(x)) next
  B <- as.matrix(x[, ..cols]) == 1; B[is.na(B)] <- FALSE; storage.mode(B) <- "double"; rownames(B) <- x$ccode
  B <- B[, colSums(B) >= 2, drop = FALSE]; if (!ncol(B)) next
  sh <- tcrossprod(B); dg <- rowSums(B); nk <- colSums(B); hy <- B %*% (t(B) / (nk - 1))
  cs <- sh / sqrt(outer(dg, dg)); cs[!is.finite(cs)] <- 0
  ix <- which(upper.tri(sh), arr.ind = TRUE); c1 <- as.integer(rownames(B)[ix[, 1]]); c2 <- as.integer(rownames(B)[ix[, 2]])
  ig[[length(ig) + 1]] <- data.table(ccode1 = pmin(c1, c2), ccode2 = pmax(c1, c2), year = yr, cnt = sh[ix], cs = cs[ix], hy = hy[ix]) }
ig <- rbindlist(ig); stopifnot(!anyDuplicated(ig[, .(ccode1, ccode2, year)]))
DY <- merge(DY, ig, by = c("ccode1", "ccode2", "year"), all.x = TRUE); setorder(DY, year, ccode1, ccode2)
for (v in c("cnt", "cs", "hy")) set(DY, which(is.na(DY[[v]])), v, 0)
## igo_count is the book / main layer itself (layers_dy igo = peacesciencer cow_igo_ndy dyadic counts, identical to the book layer_igo), so the
## grid's count coding equals the main measure. The state-year file's count (used for cos and hyper, as in the book) differs by at most a few
## IGOs in a small share of dyad-years (COW dyadic vs state-year release); reported, not used.
lay0 <- as.data.table(readRDS(a[1]))[, .(ccode1, ccode2, year, igo)]; DY <- merge(DY, lay0, by = c("ccode1", "ccode2", "year")); setorder(DY, year, ccode1, ccode2)
cat(sprintf("IGO count, state-year file vs dyadic layer: dyad-years differing %d of %d (max |diff| %g); within-year correlation min %.4f\n",
            DY[, sum(cnt != igo)], nrow(DY), DY[, max(abs(cnt - igo))], DY[, .(r = suppressWarnings(cor(cnt, igo))), by = year][is.finite(r), min(r)]))
stopifnot(DY[, mean(cnt != igo)] < 0.05, DY[, max(abs(cnt - igo))] <= 5)
DY[, `:=`(igo_count = mm(igo, year), igo_cos = mm(cs, year), igo_hyper = mm(hy, year))][, c("cnt", "cs", "hy", "igo") := NULL]

## ---- PTA depth (book 04_pta_desta_costly.R + predecessor mapping) ----
dd <- fread(file.path(DB, "desta23_dyads.csv"), encoding = "Latin-1")[, .(iso1, iso2, number = as.character(number), base_treaty, eif = entryforceyear)]
wd <- fread(file.path(DB, "desta23_dyadic_withdrawal.csv"), encoding = "Latin-1")[, .(a = pmin(iso1, iso2), b = pmax(iso1, iso2), base_treaty, wy = fcoalesce(entryforceyear, year))][, .(wy = min(wy, na.rm = TRUE)), by = .(a, b, base_treaty)]   # audit 2026-11-24: 470 withdrawal rows lack an entry-into-force year; use the signature year, as pta_book does (book 04 script ignores them)
ixd <- fread(file.path(DB, "desta23_indices.csv"), encoding = "Latin-1")[, .(number = as.character(number), base_treaty, depth = depth_index)]
bd <- ixd[number == as.character(base_treaty), .(base_treaty, bdepth = depth)]
dd <- merge(dd, ixd[, .(number, depth)], by = "number", all.x = TRUE); dd <- merge(dd, bd, by = "base_treaty", all.x = TRUE)
dd[is.na(depth), depth := bdepth][is.na(depth), depth := 0]; dd <- dd[!is.na(eif)]
dd[, `:=`(a = pmin(iso1, iso2), b = pmax(iso1, iso2))]; dd <- merge(dd, wd, by = c("a", "b", "base_treaty"), all.x = TRUE)[is.na(wy), wy := 9999L]
yrs <- 1948:2014; ex <- dd[, .(year = yrs[yrs >= eif & yrs < wy]), by = .(a, b, number, depth)]
ex <- ex[, .(maxdepth = max(depth)), by = .(a, b, year)]
map <- function(iso, yr) { cc <- suppressWarnings(countrycode(iso, "iso3n", "cown")); cc[iso == 688] <- 345L; cc[iso == 900] <- 347L
  cc[!is.na(cc) & cc == 255 & yr <= 1990] <- 260L; cc[!is.na(cc) & cc == 316 & yr <= 1992] <- 315L; cc[!is.na(cc) & cc == 679 & yr <= 1989] <- 678L; as.integer(cc) }
ex[, `:=`(c1 = map(a, year), c2 = map(b, year))]
ex <- ex[!is.na(c1) & !is.na(c2) & c1 != c2, .(maxdepth = max(maxdepth)), by = .(ccode1 = pmin(c1, c2), ccode2 = pmax(c1, c2), year)]
DY <- merge(DY, ex, by = c("ccode1", "ccode2", "year"), all.x = TRUE); setorder(DY, year, ccode1, ccode2)
DY[is.na(maxdepth), maxdepth := 0]; DY[, pta_depth := maxdepth / 7][, maxdepth := NULL]
stopifnot(all(DY$pta_depth >= 0 & DY$pta_depth <= 1))
cat(sprintf("PTA depth: dyad-years with an in-force DESTA 2.3 agreement %d; mean depth/7 among them %.3f; pta_book linked %d\n",
            sum(DY$pta_depth > 0), mean(DY$pta_depth[DY$pta_depth > 0]), sum(DY$pta_book)))

## ---- weighted layer scores ----
n <- nrow(DY)
DD <- rbind(DY[, .(ccode1, ccode2, year, contig12, ldist, igo = igo_count, def, off, dca = dca2, pta = pta_book)],
            DY[, .(ccode1 = ccode2, ccode2 = ccode1, year, contig12, ldist, igo = igo_count, def, off, dca = dca2, pta = pta_book)])
for (sch in c("sry_c12", "sry_dist", "sy_c12", "sy_dist")) {
  GEO <- if (grepl("c12", sch)) "contig12" else "ldist"; FE <- if (startsWith(sch, "sry")) "ccode1 + ccode2 + year" else "ccode1^year + ccode2^year"
  for (v in c("igo", "def", "off", "dca", "pta")) {
    y <- DD[[v]]; vy <- DD[, .(nu = uniqueN(get(v))), by = year][nu > 1, year]; idx <- which(DD$year %in% vy)
    d <- DD[idx]; fml <- as.formula(paste0(v, " ~ ", GEO, " | ", FE)); wrn <- character(0)
    fit <- withCallingHandlers(if (v == "igo") fepois(fml, data = d, notes = FALSE, glm.iter = 50) else feglm(fml, data = d, family = binomial("logit"), notes = FALSE, glm.iter = 50),
                               warning = function(e) { wrn <<- c(wrn, conditionMessage(e)); invokeRestart("muffleWarning") })
    if (!isTRUE(fit$convStatus)) stop(sch, " ", v, ": GLM did not converge")
    if (length(wrn)) stop(sprintf("%s %s: fixest warnings (run-stopping per docs/PAPER_PIPELINE.md): %s", sch, v, paste(unique(wrn), collapse = " | ")))
    pred <- rep(NA_real_, length(idx)); pred[obs(fit)] <- fitted(fit); dr <- is.na(pred)
    if (startsWith(sch, "sry")) { if (any(d[[v]][dr] != 0)) stop(sch, " ", v, ": dropped rows with y != 0"); pred[dr] <- 0 } else pred[dr] <- d[[v]][dr]
    s <- numeric(2 * n); s[idx] <- d[[v]] - pred
    sym <- (s[seq_len(n)] + s[n + seq_len(n)]) / 2
    set(DY, j = paste0("s_", sch, "_", v), value = sym)
    cat(sprintf("  %-8s %-3s years %d-%d | %s %.3f | dropped %d | mean score %.4f\n", sch, v, min(vy), max(vy), GEO, coef(fit)[[GEO]], sum(dr), mean(sym)))
  }
}
DY[, c("conttype", "capdist", "ldist") := NULL]
saveRDS(DY, OUT)
cat("GRID_LAYERS_DONE", OUT, ncol(DY), "columns\n")
