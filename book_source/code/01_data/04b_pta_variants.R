###############################################################################
# 04b_pta_variants.R   (Jared 2026-11-23; PTA fixes E and F, validation only)
# Same in-force and depth coding as 04_pta_desta_costly.R. Each variant replaces layer_ptad (0-1).
#   e2  : depth weight only for agreements with depth >= 2 (drops depth-1 deals such as Yaounde, Bangkok)
#   e4  : depth weight only for agreements with depth >= 4 (DESTA third quartile; deep agreements only)
#   f1  : drop non-reciprocal / preference schemes: EEC-ACP chain (Yaounde I-II, Arusha I-II, Lome I-IV,
#         Cotonou + accessions) and developing-country preference schemes (GSTP, PTN). base_treaty ids hard-coded.
#   f2  : size discount. Dyad value = sum over in-force base treaties of (depth/7) / (n_members_t - 1),
#         n_members_t = distinct members of that base treaty in force in year t; divided by the yearly max.
# Output: data/processed/pta_variants.rds (ccode1, ccode2, year, ptad_e2, ptad_e4, ptad_f1, ptad_f2; both dyad orders)
###############################################################################
suppressMessages({library(data.table); library(countrycode)})
set.seed(123)
BOOK <- Sys.getenv("BOOK"); stopifnot(nzchar(BOOK)); R <- file.path(BOOK, "data", "raw", "desta")
F1_DROP <- c(61L, 253L, 461L, 584L, 585L, 586L, 587L, 647L, 689L, 690L, 700L)
dy <- fread(file.path(R, "dyads.csv"), encoding = "Latin-1")[, .(iso1, iso2, number, base_treaty, eif = entryforceyear)]
wd <- fread(file.path(R, "withdrawal.csv"), encoding = "Latin-1")[, .(iso1, iso2, base_treaty, wyear = entryforceyear)]
ix <- fread(file.path(R, "indices.csv"), encoding = "Latin-1")[, .(number = as.character(number), base_treaty, depth = depth_index)]
tr <- fread(file.path(R, "treaties.csv"), encoding = "Latin-1")[, .(base_treaty, name)]
stopifnot(all(F1_DROP %in% tr$base_treaty))
cat("F1 drops:\n"); print(unique(tr[base_treaty %in% F1_DROP & !grepl("accession|withdrawal", name)]))
dy[, number := as.character(number)]
bd <- ix[number == as.character(base_treaty), .(base_treaty, bdepth = depth)]
dy <- merge(dy, ix[, .(number, depth)], by = "number", all.x = TRUE)
dy <- merge(dy, bd, by = "base_treaty", all.x = TRUE)
dy[is.na(depth), depth := bdepth][is.na(depth), depth := 0]
dy <- dy[!is.na(eif)]
wd <- wd[, .(wyear = min(wyear, na.rm = TRUE)), by = .(a = pmin(iso1, iso2), b = pmax(iso1, iso2), base_treaty)]
dy[, `:=`(a = pmin(iso1, iso2), b = pmax(iso1, iso2))]
dy <- merge(dy, wd, by = c("a", "b", "base_treaty"), all.x = TRUE)
dy[is.na(wyear), wyear := 9999L]
yrs <- 1948:2014
# dyad x base treaty x year, max depth over entries of that base treaty
dy[, rid := .I]
ex <- dy[, .(year = yrs[yrs >= eif & yrs < wyear]), by = .(rid, a, b, base_treaty, depth)]
ex <- ex[, .(depth = max(depth)), by = .(a, b, base_treaty, year)]
nm <- ex[, .(n = uniqueN(c(a, b))), by = .(base_treaty, year)]
ex <- merge(ex, nm, by = c("base_treaty", "year"))
stopifnot(all(ex$n >= 2))
dd <- ex[, .(e2 = max(fifelse(depth >= 2, depth, 0L)) / 7, e4 = max(fifelse(depth >= 4, depth, 0L)) / 7,
             f1 = { z <- depth[!base_treaty %in% F1_DROP]; if (length(z)) max(z) / 7 else 0 },
             f2 = sum(depth / 7 / (n - 1)), base = max(depth) / 7), by = .(a, b, year)]
map <- function(iso, yr) { cc <- suppressWarnings(countrycode(iso, "iso3n", "cown")); cc[!is.na(cc) & cc == 255 & yr <= 1990] <- 260L; cc }
dd[, `:=`(c1 = map(a, year), c2 = map(b, year))]
dd <- dd[!is.na(c1) & !is.na(c2) & c1 != c2, .(e2 = max(e2), e4 = max(e4), f1 = max(f1), f2 = max(f2), base = max(base)),
         by = .(ccode1 = pmin(c1, c2), ccode2 = pmax(c1, c2), year)]
L <- readRDS(file.path(BOOK, "data", "processed", "layers_dyad_year.rds"))[, .(ccode1, ccode2, year)]; setDT(L)
u <- merge(L[ccode1 < ccode2], dd, by = c("ccode1", "ccode2", "year"), all.x = TRUE)
for (v in c("e2", "e4", "f1", "f2", "base")) set(u, which(is.na(u[[v]])), v, 0)
u[, f2 := { m <- max(f2); if (m > 0) f2 / m else 0 }, by = year]
# check: base must reproduce pta_costly layer_pta_depth exactly
PC <- readRDS(file.path(BOOK, "data", "processed", "pta_costly.rds")); setDT(PC)
k <- merge(u, PC[ccode1 < ccode2, .(ccode1, ccode2, year, old = layer_pta_depth)], by = c("ccode1", "ccode2", "year"))
k[, dif := base - old]; cat("mismatch rows", k[abs(dif) > 1e-12, .N], "\n"); print(head(k[abs(dif) > 1e-12], 10)); print(warnings()[1:3]); stopifnot(nrow(k) == nrow(u), max(abs(k$dif)) < 1e-12)
stopifnot(all(u$e4 <= u$e2), all(u$e2 <= u$base), all(u$f1 <= u$base), all(u$f2 >= 0 & u$f2 <= 1))
print(u[, .(base = mean(base > 0), e2 = mean(e2 > 0), e4 = mean(e4 > 0), f1 = mean(f1 > 0), f2_mean = mean(f2), f2_mean_if_pos = mean(f2[f2 > 0])),
        by = .(era = cut(year, c(1947, 1970, 1990, 2005, 2014)))][order(era)], digits = 3)
u <- u[, .(ccode1, ccode2, year, ptad_e2 = e2, ptad_e4 = e4, ptad_f1 = f1, ptad_f2 = f2)]
out <- rbind(u, u[, .(ccode1 = ccode2, ccode2 = ccode1, year, ptad_e2, ptad_e4, ptad_f1, ptad_f2)])
stopifnot(nrow(out) == nrow(L), !anyDuplicated(out[, .(ccode1, ccode2, year)]))
saveRDS(out, file.path(BOOK, "data", "processed", "pta_variants.rds")); cat("saved", nrow(out), "rows\n")
