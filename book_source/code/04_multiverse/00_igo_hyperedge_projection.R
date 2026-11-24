###############################################################################
# 00_igo_hyperedge_projection.R   (Jared, 2026-11-17)
# IGO layer coded as hyperedges (each IGO = one hyperedge over its full members), projected
# to state pairs with size weighting, as in dynamic_multiplex.layers_from_membership(weighting="size"):
#   w_ik = sum_j B_ij B_kj / (s_j - 1),  s_j = number of full members of IGO j in year t.
# A universal IGO adds little to any one pair; a small club adds a lot.
# Members = COW IGO v3 code 1 (full); IGOs with < 2 members in year t are dropped.
# Check: layer_igo in layers_dyad_year.rds is shared COW IGO memberships min-max scaled within year,
#   so the unweighted projection sum_j B_ij B_kj, min-max scaled within year, must match it;
#   mismatches (|diff| > .01) are reported and fail loudly if > 0.5% of dyad-years.
# Output: data/processed/igo_hyper_size.rds (ccode1 < ccode2, year, igo_hyper, igo_count_check)
###############################################################################
suppressMessages({library(data.table); library(Matrix)})
set.seed(123)
BOOK <- Sys.getenv("BOOK_ROOT", getwd())
sy <- fread(file.path(BOOK, "data", "raw", "cow_igo", "state_year_formatv3.csv"))
igo_cols <- setdiff(names(sy), c("ccode", "year", "state"))
long <- melt(sy[, c("ccode", "year", igo_cols), with = FALSE], id.vars = c("ccode", "year"), variable.name = "igo", value.name = "code")
long <- long[code == 1, .(ccode, year, igo)]
cat(sprintf("full memberships: %d state-IGO-years, %d IGOs\n", nrow(long), uniqueN(long$igo)))
L <- readRDS(file.path(BOOK, "data", "processed", "layers_dyad_year.rds")); setDT(L)
L <- L[ccode1 < ccode2, .(ccode1, ccode2, year, layer_igo)]
out <- list()
for (y in sort(unique(L$year))) {
  m <- long[year == y]
  sz <- m[, .N, by = igo][N >= 2]
  m <- m[igo %in% sz$igo]
  if (!nrow(m)) next
  st <- sort(unique(m$ccode)); ig <- sz$igo
  B <- sparseMatrix(i = match(m$ccode, st), j = match(m$igo, ig), x = 1, dims = c(length(st), length(ig)))
  s <- colSums(B)
  W <- B %*% Diagonal(x = 1 / (s - 1)) %*% t(B)
  C <- B %*% t(B)
  W <- as(triu(W, 1), "TsparseMatrix"); C <- as(triu(C, 1), "TsparseMatrix")
  a <- data.table(ccode1 = st[W@i + 1], ccode2 = st[W@j + 1], igo_hyper = W@x)
  b <- data.table(ccode1 = st[C@i + 1], ccode2 = st[C@j + 1], igo_count_check = C@x)
  out[[length(out) + 1]] <- merge(a, b, by = c("ccode1", "ccode2"), all = TRUE)[, year := y]
}
H <- rbindlist(out)
stopifnot(all(H$ccode1 < H$ccode2), !anyDuplicated(H[, .(ccode1, ccode2, year)]))
H <- merge(L, H, by = c("ccode1", "ccode2", "year"), all.x = TRUE)
H[is.na(igo_hyper), igo_hyper := 0]; H[is.na(igo_count_check), igo_count_check := 0]
H[, cnt01 := (igo_count_check - min(igo_count_check)) / (max(igo_count_check) - min(igo_count_check)), by = year]
mis <- H[, mean(abs(cnt01 - layer_igo) > .01)]
cat(sprintf("panel dyad-years %d; share where projected count != layer_igo: %.4f; cor(hyper, count) %.3f\n",
            nrow(H), mis, H[, cor(igo_hyper, layer_igo)]))
print(H[abs(cnt01 - layer_igo) > .01, .N, by = .(decade = 10 * (year %/% 10))][order(decade)][1:10])
stopifnot(mis < 0.005)
saveRDS(H[, .(ccode1, ccode2, year, igo_hyper, igo_count_check)], file.path(BOOK, "data", "processed", "igo_hyper_size.rds"))
