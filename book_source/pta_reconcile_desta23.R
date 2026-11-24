## Diagnostic (2026-11-24): independent DESTA 2.3 clique rebuild vs the book prebuilt PTA layer (layer_pta_desta). Paths are sandbox paths; run from a folder holding data/raw/ as in data_entry/book.
suppressPackageStartupMessages({ library(data.table); library(countrycode) })
B <- as.data.table(readRDS("data/raw/prebuilt/pta_layers.rds"))[ccode1 < ccode2 & layer_pta_desta > 0, .(ccode1, ccode2, year, book = 1L)]
FR <- as.data.table(readRDS("/mnt/user-data/uploads/targeting_multiverse/output/layers_dy.rds"))[, .(ccode1, ccode2, year)]
cw <- function(iso, yr) { x <- suppressWarnings(countrycode(iso, "iso3n", "cown")); x[!is.na(x) & x == 255 & yr <= 1990] <- 260L; x[iso == 688] <- 345L; as.integer(x) }
dy <- fread("data/raw/desta/dyads.csv", encoding = "Latin-1")
wd <- fread("data/raw/desta/withdrawal.csv", encoding = "Latin-1")[, .(iso = c(iso1, iso2), base_treaty, wyear = entryforceyear)]
cmp <- function(lab, PR) {           # PR: ccode1 ccode2 year linked
  x <- merge(FR, PR[, .(ccode1, ccode2, year, mine = 1L)], by = c("ccode1","ccode2","year"))
  m <- merge(B, x, all = TRUE); m[is.na(book), book := 0L]; m[is.na(mine), mine := 0L]
  cat(sprintf("%-38s total %7d | both %7d | book-only %6d | mine-only %6d\n", lab, sum(m$mine), sum(m$book & m$mine), sum(m$book & !m$mine), sum(!m$book & m$mine))) }
yrs <- 1948:2014
## (1) DESTA 2.3 dyadic in force with withdrawals (book's 04 rule)
w2 <- fread("data/raw/desta/withdrawal.csv", encoding = "Latin-1")[, .(a = pmin(iso1, iso2), b = pmax(iso1, iso2), base_treaty, wy = entryforceyear)][, .(wy = min(wy)), by = .(a, b, base_treaty)]
d1 <- dy[!is.na(entryforceyear), .(a = pmin(iso1, iso2), b = pmax(iso1, iso2), base_treaty, eif = entryforceyear)]
d1 <- merge(d1, w2, by = c("a","b","base_treaty"), all.x = TRUE)[is.na(wy), wy := 9999L]
e1 <- d1[, .(year = yrs[yrs >= eif & yrs < wy]), by = .(a, b)][, unique(.SD)]
e1[, `:=`(c1 = cw(a, year), c2 = cw(b, year))]; e1 <- e1[!is.na(c1) & !is.na(c2) & c1 != c2, .(ccode1 = pmin(c1, c2), ccode2 = pmax(c1, c2), year)][, unique(.SD)]
cmp("2.3 dyadic in force, withdrawals", e1)
## (2) DESTA 2.3 clique in force, member join = min eif in base treaty, member exit = withdrawal year
mem <- rbind(dy[, .(bt = base_treaty, iso = iso1, eif = entryforceyear)], dy[, .(bt = base_treaty, iso = iso2, eif = entryforceyear)])[!is.na(eif)][, .(join = min(eif)), by = .(bt, iso)]
ww <- wd[, .(exit = min(wyear)), by = .(bt = base_treaty, iso)]
mem <- merge(mem, ww, by = c("bt","iso"), all.x = TRUE)[is.na(exit), exit := 9999L]
pr <- mem[mem, on = "bt", allow.cartesian = TRUE][iso < i.iso, .(a = iso, b = i.iso, s = pmax(join, i.join), e = pmin(exit, i.exit))][s < e]
e2 <- pr[, .(year = yrs[yrs >= s & yrs < e]), by = .(a, b, s, e)][, .(a, b, year)][, unique(.SD)]
e2[, `:=`(c1 = cw(a, year), c2 = cw(b, year))]; e2 <- e2[!is.na(c1) & !is.na(c2) & c1 != c2, .(ccode1 = pmin(c1, c2), ccode2 = pmax(c1, c2), year)][, unique(.SD)]
cmp("2.3 clique in force, withdrawals", e2)
## (3) clique without withdrawals
pr3 <- mem[mem, on = "bt", allow.cartesian = TRUE][iso < i.iso, .(a = iso, b = i.iso, s = pmax(join, i.join))]
e3 <- pr3[, .(year = yrs[yrs >= s]), by = .(a, b, s)][, .(a, b, year)][, unique(.SD)]
e3[, `:=`(c1 = cw(a, year), c2 = cw(b, year))]; e3 <- e3[!is.na(c1) & !is.na(c2) & c1 != c2, .(ccode1 = pmin(c1, c2), ccode2 = pmax(c1, c2), year)][, unique(.SD)]
cmp("2.3 clique in force, no withdrawals", e3)
d <- readRDS("/home/claude/pta_diag_dy.rds")[A == 1L, .(ccode1, ccode2, year)]; cmp("old destaRaw clique in force (pta_clq)", d)
x <- merge(FR, e2[, .(ccode1, ccode2, year, mine = 1L)], by = c("ccode1","ccode2","year")); m <- merge(B, x, all = TRUE); m[is.na(book), book := 0L]; m[is.na(mine), mine := 0L]
mo <- m[mine == 1 & book == 0]; cat("\nmine-only by year (selected):\n"); print(mo[, .N, by = year][order(year)][N > 100])
cat("\nmine-only top countries:\n"); print(head(sort(table(c(mo$ccode1, mo$ccode2)), decreasing = TRUE), 25))
## which base treaties generate the mine-only pairs (2007 sample)
s <- mo[year == 2008][1:5]; print(s)
bo <- m[mine == 0 & book == 1]; cat("\nbook-only top pairs:\n"); print(head(bo[, .N, by = .(ccode1, ccode2)][order(-N)], 10))
saveRDS(list(mem = mem, e2 = e2), "/home/claude/bookrep/pta23_clique.rds")
mem2 <- copy(mem); mem2[bt == 647, exit := pmin(exit, 1981L)]; mem2[bt == 461, exit := pmin(exit, 2007L)]
pr4 <- mem2[mem2, on = "bt", allow.cartesian = TRUE][iso < i.iso, .(a = iso, b = i.iso, s = pmax(join, i.join), e = pmin(exit, i.exit))][s < e]
e4 <- pr4[, .(year = yrs[yrs >= s & yrs < e]), by = .(a, b, s, e)][, .(a, b, year)][, unique(.SD)]
e4[, `:=`(c1 = cw(a, year), c2 = cw(b, year))]; e4 <- e4[!is.na(c1) & !is.na(c2) & c1 != c2, .(ccode1 = pmin(c1, c2), ccode2 = pmax(c1, c2), year)][, unique(.SD)]
cmp("2.3 clique, PTN to 1980, GSTP to 2006", e4)
x <- merge(FR, e4[, .(ccode1, ccode2, year, mine = 1L)], by = c("ccode1","ccode2","year")); m <- merge(B, x, all = TRUE); m[is.na(book), book := 0L]; m[is.na(mine), mine := 0L]
cat("remaining mine-only by year:\n"); print(m[mine==1&book==0, .N, by=year][N>20][order(year)]); print(head(m[mine==1&book==0, .N, by=.(ccode1,ccode2)][order(-N)],10))
