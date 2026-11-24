###############################################################################
# 06_igo_variants.R   (Jared 2026-11-23) IGO layer constructions from COW IGO v3 full memberships (code 1), all IGOs.
# B_t = state x IGO full-membership matrix in year t; d_i = number of IGOs of state i; n_k = full members of IGO k.
#   igo_count : shared_ij = sum_k B_ik B_jk                       (original JOP layer, before scaling)
#   igo_log   : log(1 + shared_ij)
#   igo_cos   : shared_ij / sqrt(d_i d_j)                         (bipartite projection normalized by state degree)
#   igo_hyper : sum_k B_ik B_jk / (n_k - 1), IGOs with n_k >= 2   (hyperedge / Newman projection)
# Output: data/processed/igo_variants.rds (ccode1 < ccode2, year, four raw columns; builder min-max scales within year)
###############################################################################
suppressMessages(library(data.table)); set.seed(123)
BOOK <- Sys.getenv("BOOK"); stopifnot(nzchar(BOOK))
w <- fread(file.path(BOOK, "data", "raw", "cow_igo", "state_year_formatv3.csv"))
cols <- setdiff(names(w), c("ccode", "year", "state"))
fr <- readRDS(file.path(BOOK, "data", "processed", "layers_dyad_year.rds")); setDT(fr); fr <- fr[ccode1 < ccode2, .(ccode1, ccode2, year)]
out <- list()
for (yr in sort(unique(fr$year))) {
  x <- w[year == yr]; if (!nrow(x)) next
  B <- as.matrix(x[, ..cols]) == 1; B[is.na(B)] <- FALSE; storage.mode(B) <- "double"; rownames(B) <- x$ccode
  B <- B[, colSums(B) >= 2, drop = FALSE]
  if (!ncol(B)) next
  sh <- tcrossprod(B); dg <- rowSums(B); nk <- colSums(B)
  hy <- B %*% (t(B) / (nk - 1))
  cs <- sh / sqrt(outer(dg, dg)); cs[!is.finite(cs)] <- 0
  ix <- which(upper.tri(sh), arr.ind = TRUE)
  c1 <- as.integer(rownames(B)[ix[, 1]]); c2 <- as.integer(rownames(B)[ix[, 2]])
  out[[length(out) + 1]] <- data.table(ccode1 = pmin(c1, c2), ccode2 = pmax(c1, c2), year = yr, igo_count = sh[ix], igo_log = log1p(sh[ix]), igo_cos = cs[ix], igo_hyper = hy[ix])
}
D <- rbindlist(out)
D <- merge(fr, D, by = c("ccode1", "ccode2", "year"), all.x = TRUE)
for (v in c("igo_count", "igo_log", "igo_cos", "igo_hyper")) set(D, which(is.na(D[[v]])), v, 0)
stopifnot(!anyDuplicated(D[, .(ccode1, ccode2, year)]), all(D$igo_cos <= 1 + 1e-9), all(D$igo_hyper <= D$igo_count + 1e-9))
# check: count matches the original book layer up to within-year scaling
L <- readRDS(file.path(BOOK, "data", "processed", "layers_dyad_year.rds")); setDT(L)
k <- merge(D, L[ccode1 < ccode2, .(ccode1, ccode2, year, layer_igo)], by = c("ccode1", "ccode2", "year"))
cat("correlation with original layer_igo (pooled):", round(cor(k$igo_count, k$layer_igo), 3), "; within-year mean:", round(k[, cor(igo_count, layer_igo), by = year][, mean(V1, na.rm = TRUE)], 3), "\n")
print(D[, .(share_pos = round(mean(igo_count > 0), 3), r_count_log = round(cor(igo_count, igo_log), 2), r_count_cos = round(cor(igo_count, igo_cos), 2), r_count_hyper = round(cor(igo_count, igo_hyper), 2)),
        by = .(era = cut(year, c(1815, 1913, 1945, 1991, 2014)))][order(era)])
saveRDS(D, file.path(BOOK, "data", "processed", "igo_variants.rds")); cat("saved", nrow(D), "\n")
