###############################################################################
# 04_pta_desta_costly.R   (Jared 2026-11-23; PTA plan steps 1-4)
# Costly PTA layer from DESTA 2.3 (designoftradeagreements.org, downloaded 2026-11-23):
#   desta_list_of_treaties_02_03_dyads.csv (base treaties, accessions, consolidated versions),
#   desta_dyadic_withdrawal_02_03.csv (exits), desta_indices_version_02_03.csv (depth_index 0-7).
# Coding (fixed before any outcome model):
#   in force   : a dyad's agreement counts from its entry-into-force year (entryforceyear) until the year
#                before a recorded withdrawal of either member from that base treaty; agreements with no
#                entry-into-force year are treated as not in force (reported).
#   depth      : DESTA depth_index of the agreement entry (consolidated versions use their own coding; accessions
#                take the base treaty's depth); dyad-year value = max depth over in-force agreements.
#   layers     : layer_pta_inforce (0/1), layer_pta_depth (max depth / 7, 0-1), layer_pta_deep (max depth >= 4;
#                threshold = DESTA third quartile, set before outcomes).
# Country codes: ISO numeric -> COW (countrycode 1.9.0); Germany 276 -> 260 through 1990, 255 after.
# Output: data/processed/pta_costly.rds (ccode1, ccode2, year, three layers; both dyad orders, 1816-2014)
###############################################################################
suppressMessages({library(data.table); library(countrycode)})
set.seed(123)
BOOK <- Sys.getenv("BOOK"); stopifnot(nzchar(BOOK)); R <- file.path(BOOK, "data", "raw", "desta")
cat("countrycode", as.character(packageVersion("countrycode")), "\n")
dy <- fread(file.path(R, "dyads.csv"), encoding = "Latin-1")[, .(iso1, iso2, number, base_treaty, entry_type, year, eif = entryforceyear)]
wd <- fread(file.path(R, "withdrawal.csv"), encoding = "Latin-1")[, .(iso1, iso2, base_treaty, wyear = entryforceyear)]
ix <- fread(file.path(R, "indices.csv"), encoding = "Latin-1")[, .(number = as.character(number), base_treaty, depth = depth_index)]
dy[, number := as.character(number)]
cat(sprintf("DESTA dyad rows %d; entries without entry-into-force year %.3f (treated as not in force)\n", nrow(dy), mean(is.na(dy$eif))))
# depth: own coding if present, else base treaty's
bd <- ix[number == as.character(base_treaty), .(base_treaty, bdepth = depth)]
dy <- merge(dy, ix[, .(number, depth)], by = "number", all.x = TRUE)
dy <- merge(dy, bd, by = "base_treaty", all.x = TRUE)
dy[is.na(depth), depth := bdepth]
cat(sprintf("entries with no depth coding (own or base): %.3f -> depth set to 0 but in force\n", mean(is.na(dy$depth))))
dy[is.na(depth), depth := 0]
dy <- dy[!is.na(eif)]
# withdrawals: first withdrawal year per dyad x base treaty ends that agreement for the dyad
wd <- wd[, .(wyear = min(wyear, na.rm = TRUE)), by = .(a = pmin(iso1, iso2), b = pmax(iso1, iso2), base_treaty)]
dy[, `:=`(a = pmin(iso1, iso2), b = pmax(iso1, iso2))]
dy <- merge(dy, wd, by = c("a", "b", "base_treaty"), all.x = TRUE)
dy[is.na(wyear), wyear := 9999L]
cat(sprintf("dyad-agreement rows ended by withdrawal: %d\n", dy[wyear < 9999, .N]))
# expand to dyad-years 1948-2014
yrs <- 1948:2014
ex <- dy[, .(year = yrs[yrs >= eif & yrs < wyear]), by = .(a, b, number, depth)]
dd <- ex[, .(inforce = 1L, maxdepth = max(depth)), by = .(a, b, year)]
# ISO numeric -> COW
map <- function(iso, yr) { cc <- suppressWarnings(countrycode(iso, "iso3n", "cown")); cc[!is.na(cc) & cc == 255 & yr <= 1990] <- 260L; cc }
dd[, `:=`(c1 = map(a, year), c2 = map(b, year))]
cat(sprintf("dyad-years unmapped to COW: %.4f (dropped)\n", dd[, mean(is.na(c1) | is.na(c2))]))
dd <- dd[!is.na(c1) & !is.na(c2) & c1 != c2, .(inforce = max(inforce), maxdepth = max(maxdepth)), by = .(ccode1 = pmin(c1, c2), ccode2 = pmax(c1, c2), year)]
# frame: book dyad-years
L <- readRDS(file.path(BOOK, "data", "processed", "layers_dyad_year.rds"))[, .(ccode1, ccode2, year)]; setDT(L)
u <- L[ccode1 < ccode2]
u <- merge(u, dd, by = c("ccode1", "ccode2", "year"), all.x = TRUE)
cat(sprintf("DESTA dyad-years matched to the book frame: %.3f\n", nrow(merge(dd, L[ccode1 < ccode2], by = c("ccode1", "ccode2", "year"))) / nrow(dd[year <= 2014])))
u[is.na(inforce), `:=`(inforce = 0L, maxdepth = 0L)]
u[, `:=`(layer_pta_inforce = inforce, layer_pta_depth = maxdepth / 7, layer_pta_deep = as.integer(maxdepth >= 4))]
stopifnot(all(u$layer_pta_depth >= 0 & u$layer_pta_depth <= 1), !anyNA(u$layer_pta_inforce))
# --- validation ---
P <- readRDS(file.path(BOOK, "data", "processed", "pta_layers.rds")); setDT(P)
v <- merge(u, P[ccode1 < ccode2, .(ccode1, ccode2, year, old = layer_pta_desta)], by = c("ccode1", "ccode2", "year"))
cat(sprintf("\nshare of dyad-years in force: new %.3f vs old book layer %.3f (1950-2014: %.3f vs %.3f); agreement %.3f\n",
  mean(v$layer_pta_inforce), mean(v$old), v[year >= 1950, mean(layer_pta_inforce)], v[year >= 1950, mean(old)], mean(v$layer_pta_inforce == v$old)))
chk <- function(a, b, yr) u[ccode1 == min(a, b) & ccode2 == max(a, b) & year == yr, sprintf("%d-%d %d: inforce %d, depth %d/7", a, b, yr, layer_pta_inforce, round(layer_pta_depth * 7))]
cat("\nspot checks\n"); for (x in list(c(220, 260, 1960), c(220, 260, 1957), c(2, 70, 1993), c(2, 70, 1995), c(2, 20, 1990), c(140, 160, 1992), c(140, 160, 1996),
  c(200, 220, 1972), c(200, 220, 1974), c(820, 850, 1995), c(2, 666, 2002), c(365, 369, 2000), c(710, 732, 2014))) print(chk(x[1], x[2], x[3]))
print(u[, .(inforce = round(mean(layer_pta_inforce), 3), deep = round(mean(layer_pta_deep), 3), mean_depth_if_inforce = round(7 * mean(layer_pta_depth[layer_pta_inforce == 1]), 2)),
        by = .(era = cut(year, c(1815, 1945, 1970, 1990, 2005, 2014)))][order(era)])
out <- rbind(u[, .(ccode1, ccode2, year, layer_pta_inforce, layer_pta_depth, layer_pta_deep)],
             u[, .(ccode1 = ccode2, ccode2 = ccode1, year, layer_pta_inforce, layer_pta_depth, layer_pta_deep)])
stopifnot(nrow(out) == nrow(L), !anyDuplicated(out[, .(ccode1, ccode2, year)]))
saveRDS(out, file.path(BOOK, "data", "processed", "pta_costly.rds")); cat("saved", nrow(out), "rows\n")
