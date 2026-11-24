###############################################################################
# 03_atop_un_layers.R   (Jared 2026-11-21; nine-layer grid from relational_measures/01_build_layers.R)
# Adds to the book's dyad-year frame:
#   layer_def, layer_off, layer_neu, layer_nonagg, layer_consul  ATOP obligation types (0/1), from
#     peacesciencer::atop_alliance; an obligation recorded in either direction counts for the pair
#     (same rule as relational_measures). Same coding as the targeting project.
#   igo_count  shared COW IGO memberships (peacesciencer cow_igo_ndy, dyadigos); within IGO-data years a
#              missing dyad is 0; years with no IGO data NA (as relational_measures/01_build_layers.R)
#   layer_un   UN voting agreement, kappa (Hage FPSIM, kappavv; peacesciencer extdata), 1946+.
#              Pairs with no kappa in a UN year (non-members) are NA here; the measure builder
#              treats them as 0 surprise. Years before the first UN year are NA (structural).
# Run with an R that has peacesciencer 1.2.0 and its extdata (download_extdata()).
# Input : BOOK/data/processed/layers_dyad_year.rds (directed, both orders)
# Output: BOOK/data/processed/atop_un_layers.rds (ccode1, ccode2, year, six layers; both orders)
###############################################################################
suppressMessages({library(data.table); library(peacesciencer)})
set.seed(123)
BOOK <- Sys.getenv("BOOK"); stopifnot(nzchar(BOOK))
cat("peacesciencer", as.character(packageVersion("peacesciencer")), "\n")
d <- readRDS(file.path(BOOK, "data", "processed", "layers_dyad_year.rds"))[, .(ccode1, ccode2, year, layer_defensive, layer_offensive)]
setDT(d); stopifnot(!anyDuplicated(d[, .(ccode1, ccode2, year)]))

a <- as.data.table(peacesciencer::atop_alliance)
a[, `:=`(c1 = as.integer(pmin(ccode1, ccode2)), c2 = as.integer(pmax(ccode1, ccode2)), year = as.integer(year))]
a <- a[, .(def = as.integer(max(atop_defense) > 0), off = as.integer(max(atop_offense) > 0), neu = as.integer(max(atop_neutral) > 0),
           nonagg = as.integer(max(atop_nonagg) > 0), consul = as.integer(max(atop_consul) > 0)), by = .(c1, c2, year)]
stopifnot(a[c1 == 211 & c2 == 220 & year == 1832, def] == 1L)   # one-direction obligation survives (France-Belgium 1832)

fp <- as.data.table(readRDS(file.path(system.file("extdata", package = "peacesciencer"), "dyadic_fp_similarity.rds")))
fp <- fp[!is.na(kappavv), .(un = mean(kappavv)), by = .(c1 = as.integer(pmin(ccode1, ccode2)), c2 = as.integer(pmax(ccode1, ccode2)), year = as.integer(year))]
un_first <- min(fp$year)
stopifnot(un_first >= 1945L, all(fp$un >= -1 - 1e-9 & fp$un <= 1 + 1e-9))

ig <- as.data.table(peacesciencer::cow_igo_ndy)
ig <- unique(ig[, .(c1 = as.integer(pmin(ccode1, ccode2)), c2 = as.integer(pmax(ccode1, ccode2)), year = as.integer(year),
                    igo_count = fifelse(dyadigos < 0, NA_real_, as.numeric(dyadigos)))], by = c("c1", "c2", "year"))
d[, `:=`(c1 = pmin(ccode1, ccode2), c2 = pmax(ccode1, ccode2))]
n0 <- nrow(d)
d <- merge(d, a, by = c("c1", "c2", "year"), all.x = TRUE)
d <- merge(d, fp, by = c("c1", "c2", "year"), all.x = TRUE)
d <- merge(d, ig, by = c("c1", "c2", "year"), all.x = TRUE)
igy <- d[, .(ok = any(!is.na(igo_count))), by = year]
d <- merge(d, igy, by = "year"); d[ok == TRUE & is.na(igo_count), igo_count := 0][, ok := NULL]
cat(sprintf("IGO count: years without data %s\n", paste(igy[ok == FALSE, year], collapse = " ")))
stopifnot(nrow(d) == n0)
for (v in c("def", "off", "neu", "nonagg", "consul")) set(d, which(is.na(d[[v]])), v, 0L)
d[year < un_first, un := NA_real_]
cat(sprintf("UN kappa first year %d; NA share in UN years %.3f\n", un_first, d[year >= un_first, mean(is.na(un))]))
cat(sprintf("agreement with book layers: def %.4f, off %.4f\n", mean(d$def == d$layer_defensive), mean(d$off == d$layer_offensive)))
print(d[, .(def = mean(def), off = mean(off), neu = mean(neu), nonagg = mean(nonagg), consul = mean(consul), un = mean(un, na.rm = TRUE)),
        by = .(era = cut(year, c(1815, 1870, 1913, 1945, 1991, 2014)))][order(era)])
out <- d[, .(ccode1, ccode2, year, layer_def = def, layer_off = off, layer_neu = neu, layer_nonagg = nonagg, layer_consul = consul, layer_un = un, igo_count)]
setorder(out, year, ccode1, ccode2)
saveRDS(out, file.path(BOOK, "data", "processed", "atop_un_layers.rds"))
cat("saved", nrow(out), "rows\n")
