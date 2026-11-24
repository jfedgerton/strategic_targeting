## R/book_layers.R -- the book's PTA and DCA layers. PTA = the book layer (layer_pta_desta), ported from book code/01_data/03_trade_agreements.R (withdrawal fix
## 2026-11-24) plus the predecessor-state mapping agreed 2026-11-24 (PI): ISO 203 -> COW 315 (Czechoslovakia) through 1992 and ISO 887 ->
## COW 678 (North Yemen) through 1989. Only layer_pta_desta is built (the Mansfield-Pevehouse layer of the book script is not used here).
##   - base_treaty and accession rows: listed members join in entryforceyear (rows without one are not in force)
##   - withdrawal rows: only the states that actually withdrew leave (DESTA 2.3 dyadic withdrawal file, same event `number`: the states
##     in the most withdrawal pairs; a one-pair event withdraws both), from the withdrawal year on
##   - every member pair of a base treaty in force in year t is linked (binary); 1816-2014
##   - ISO numeric -> COW (countrycode); Germany 276 -> 260 for 1955-1990; Serbia 688 -> 345; Kosovo 900 -> 347
## Check: 196,693 linked undirected dyad-years on the layer frame (book file without the predecessor mapping: 196,669); stops otherwise.
## DCA (PI 2026-11-24, book alignment): Kinne DCAD dcaAnyV2 (the book / JOP main; layers_dy holds dcaAnyV1), binary, true zeros before 1980,
##   each dyad's 2010 value carried to 2011-14. Check: equals the book layer_dca for 1980-2010.
## Usage: Rscript R/book_layers.R LAYERS_RDS DESTA_02_02_TXT DESTA23_WITHDRAWAL_CSV DCA_DYAD_CSV OUT_RDS   (out: ccode1, ccode2, year, pta_desta, dca_v2)
suppressPackageStartupMessages({ library(data.table); library(countrycode) })
set.seed(123)
a <- commandArgs(trailingOnly = TRUE); stopifnot(length(a) == 5)
FR <- as.data.table(readRDS(a[1]))[, .(ccode1, ccode2, year)]; stopifnot(all(FR$ccode1 < FR$ccode2), !anyDuplicated(FR))
x <- fread(a[2])
mcols <- grep("^c[0-9]+$", names(x), value = TRUE)
long <- melt(x[, c("number", "base_treaty", "entry_type", "year", "entryforceyear", mcols), with = FALSE],
             id.vars = c("number", "base_treaty", "entry_type", "year", "entryforceyear"), value.name = "iso", na.rm = TRUE)
long <- long[!is.na(iso) & iso != ""][, iso := as.integer(iso)]
join <- long[entry_type %in% c("base_treaty", "accession") & !is.na(entryforceyear), .(base_treaty, iso, start = entryforceyear)]
w23 <- fread(a[3], encoding = "Latin-1")[, .(number, base_treaty, iso1, iso2)]
wdr <- w23[, { cnt <- table(c(iso1, iso2)); .(iso = as.integer(names(cnt)[cnt == max(cnt)])) }, by = .(number, base_treaty)]
ev22 <- unique(long[entry_type == "withdrawal", .(number, base_treaty)])
stopifnot(nrow(merge(ev22, unique(wdr[, .(number, base_treaty)]), by = c("number", "base_treaty"))) == nrow(ev22))
wd <- long[entry_type == "withdrawal", .(number, base_treaty, iso, stop = fifelse(is.na(entryforceyear), year, entryforceyear))]
n_listed <- nrow(wd)
wd <- merge(wd, wdr, by = c("number", "base_treaty", "iso"))[, .(base_treaty, iso, stop)]
cat(sprintf("withdrawal rows: %d listed members, %d actual withdrawers\n", n_listed, nrow(wd)))
stopifnot(nrow(wd) == 93L)
mem <- join[, .(start = min(start)), by = .(base_treaty, iso)]
mem <- merge(mem, wd[, .(stop = min(stop)), by = .(base_treaty, iso)], by = c("base_treaty", "iso"), all.x = TRUE)
mem[is.na(stop), stop := 9999L]; mem <- mem[start < stop]
my <- mem[, .(year = seq(max(start, 1816L), min(stop - 1L, 2014L))), by = .(base_treaty, iso)][!is.na(year)]
map <- data.table(iso = unique(my$iso)); map[, ccode := suppressWarnings(countrycode(iso, "iso3n", "cown"))]
map[iso == 688, ccode := 345L]; map[iso == 900, ccode := 347L]
cat("ISO codes without a COW code (non-COW territories):", paste(sort(map[is.na(ccode), iso]), collapse = ", "), "\n")
my <- merge(my, map, by = "iso")[!is.na(ccode)]
my[iso == 276 & year %between% c(1955, 1990), ccode := 260L]
my[iso == 203 & year <= 1992, ccode := 315L]   # Czechoslovakia (DESTA codes it as the Czech Republic)
my[iso == 887 & year <= 1989, ccode := 678L]   # North Yemen (no PTA before 1990 in DESTA; for consistency)
my <- unique(my[, .(base_treaty, ccode, year)])
pr <- merge(my, my, by = c("base_treaty", "year"), allow.cartesian = TRUE)[ccode.x < ccode.y]
P <- unique(pr[, .(ccode1 = ccode.x, ccode2 = ccode.y, year)])[, pta_desta := 1L]
OUT <- merge(FR, P, by = c("ccode1", "ccode2", "year"), all.x = TRUE)[is.na(pta_desta), pta_desta := 0L]
stopifnot(nrow(OUT) == nrow(FR))
n <- sum(OUT$pta_desta); cat(sprintf("linked undirected dyad-years %d (years %d-%d); Czechoslovakia (315) %d\n", n,
  min(OUT[pta_desta == 1, year]), max(OUT[pta_desta == 1, year]), OUT[pta_desta == 1 & (ccode1 == 315 | ccode2 == 315), .N]))
stopifnot(n == 196693L)
## DCA v2
dc <- fread(a[4])[, .(ccode1 = pmin(ccode1, ccode2), ccode2 = pmax(ccode1, ccode2), year, v = as.integer(dcaAnyV2 > 0))][, .(v = max(v)), by = .(ccode1, ccode2, year)]
stopifnot(min(dc$year) == 1980L, max(dc$year) == 2010L)
OUT <- merge(OUT, dc, by = c("ccode1", "ccode2", "year"), all.x = TRUE)
d10 <- OUT[year == 2010L, .(ccode1, ccode2, v10 = fcoalesce(v, 0L))]
OUT <- merge(OUT, d10, by = c("ccode1", "ccode2"), all.x = TRUE)
OUT[, dca_v2 := fifelse(year < 1980L, 0L, fifelse(year <= 2010L, fcoalesce(v, 0L), fcoalesce(v10, 0L)))][, c("v", "v10") := NULL]
setorder(OUT, year, ccode1, ccode2); stopifnot(nrow(OUT) == nrow(FR), !anyNA(OUT$dca_v2))
cat(sprintf("DCA v2: linked dyad-years 1980-2010 %d (book 11,031); 2011-14 carried %d\n", OUT[year <= 2010L, sum(dca_v2)], OUT[year > 2010L, sum(dca_v2)]))
stopifnot(OUT[year <= 2010L, sum(dca_v2)] == 11031L)
saveRDS(OUT, a[5])
cat("BOOK_LAYERS_DONE", a[5], "\n")
