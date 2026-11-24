###############################################################################
# 03_trade_agreements.R
# Dyad-year preferential trade agreement layers (withdrawal bug fixed 2026-11-24, see below; ported from ottoman_byzantine,
# replication-sweep raw data). See data/raw/trade_agreements/SOURCES.md.
#
# layer_pta_desta    1 if i and j are both members of at least one DESTA agreement
#                    in force in year t.
#   - base_treaty and accession rows: listed members join in entryforceyear
#     (rows without an entry-into-force year are treated as not in force)
#   - withdrawal rows: the states that actually withdrew (identified from the DESTA 2.3
#     dyadic withdrawal file for the same event) leave that base treaty from the
#     withdrawal year on (fixed 2026-11-24; previously every listed member left)
#   - negotiation, consolidated, protocol/amendment rows ignored
#   - members are ISO numeric (c1..cN) -> COW via countrycode; Germany (276)
#     maps to FRG (260) for 1955-1990 and to GMY (255) otherwise; Serbia 688 -> 345,
#     Kosovo 900 -> 347
# layer_pta_mansfield 1 if PTASame > 0 in Mansfield & Pevehouse (2000) PTAMID,
#                    1947-2000 (0 outside, so the layer drops out after 2000
#                    under the varying-years rule)
#
# Output: data/processed/pta_layers.rds  (directed: (i,j) and (j,i), like the panel)
#         results/main/pta_coverage.csv
###############################################################################
suppressMessages({library(data.table); library(countrycode); library(haven)})
set.seed(123)
BOOK <- Sys.getenv("BOOK_ROOT", getwd())
RAW <- file.path(BOOK, "data", "raw", "trade_agreements")
panel <- readRDS(file.path(BOOK, "data", "processed", "layers_dyad_year.rds"))[, .(ccode1, ccode2, year)]

# ---- DESTA ------------------------------------------------------------------
x <- fread(file.path(RAW, "desta_list_of_treaties_02_02.txt"))
mcols <- grep("^c[0-9]+$", names(x), value = TRUE)
long <- melt(x[, c("number", "base_treaty", "entry_type", "year", "entryforceyear", mcols), with = FALSE],
             id.vars = c("number", "base_treaty", "entry_type", "year", "entryforceyear"),
             value.name = "iso", na.rm = TRUE)
long <- long[!is.na(iso) & iso != ""]
long[, iso := as.integer(iso)]

join <- long[entry_type %in% c("base_treaty", "accession")]
n_no_eif <- uniqueN(join[is.na(entryforceyear), number])
join <- join[!is.na(entryforceyear), .(base_treaty, iso, start = entryforceyear)]
# BUG FIX (Jared 2026-11-24): DESTA v02_02 withdrawal rows list the whole membership, not only the states that left
#   (e.g. PTN 1981 lists all 16 members, GSTP 2007 all 48), so withdrawing everyone listed ended PTN after 1980 and GSTP
#   after 2006 for all members. The withdrawing states are now taken from the DESTA 2.3 dyadic withdrawal file
#   (data/raw/desta/withdrawal.csv), which codes the same 72 events (matched by `number`): in each event the withdrawers
#   are the states that appear in the most withdrawal pairs (a state leaving pairs with every remaining member). A
#   one-pair event (bilateral treaty) withdraws both. Only listed members who are withdrawers stop.
w23 <- fread(file.path(BOOK, "data", "raw", "desta", "withdrawal.csv"), encoding = "Latin-1")[, .(number, base_treaty, iso1, iso2)]
wdr <- w23[, { cnt <- table(c(iso1, iso2)); .(iso = as.integer(names(cnt)[cnt == max(cnt)])) }, by = .(number, base_treaty)]
ev22 <- unique(long[entry_type == "withdrawal", .(number, base_treaty)])
stopifnot(nrow(merge(ev22, unique(wdr[, .(number, base_treaty)]), by = c("number", "base_treaty"))) == nrow(ev22))   # every v02_02 event found in 2.3
wd <- long[entry_type == "withdrawal", .(number, base_treaty, iso, stop = fifelse(is.na(entryforceyear), year, entryforceyear))]
n_listed <- nrow(wd)
wd <- merge(wd, wdr, by = c("number", "base_treaty", "iso"))[, .(base_treaty, iso, stop)]
cat(sprintf("withdrawal rows: %d listed members, %d actual withdrawers kept (DESTA 2.3)\n", n_listed, nrow(wd)))
stopifnot(nrow(wd) > 0, nrow(wd) < n_listed)
cat(sprintf("DESTA: %d join rows (dropped %d treaty/accession entries with no entry into force), %d withdrawal rows\n",
            nrow(join), n_no_eif, nrow(wd)))

mem <- join[, .(start = min(start)), by = .(base_treaty, iso)]
mem <- merge(mem, wd[, .(stop = min(stop)), by = .(base_treaty, iso)], by = c("base_treaty", "iso"), all.x = TRUE)
mem[is.na(stop), stop := 9999L]
mem <- mem[start < stop]

# Expand to member-years within 1816-2014 and map ISO -> COW
my <- mem[, .(year = seq(max(start, 1816L), min(stop - 1L, 2014L))), by = .(base_treaty, iso)][!is.na(year)]
iso_u <- unique(my$iso)
cow <- suppressWarnings(countrycode(iso_u, "iso3n", "cown"))
map <- data.table(iso = iso_u, ccode = cow)
map[iso == 688, ccode := 345L]   # Serbia -> COW 345 (Yugoslavia/Serbia)
map[iso == 900, ccode := 347L]   # Kosovo (DESTA code 900) -> COW 347
# Remaining unmapped codes are non-COW territories (HK, Macao, Palestine, Greenland, etc.)
cat("ISO codes without a COW code:", paste(map[is.na(ccode), iso], collapse = ", "), "\n")
my <- merge(my, map, by = "iso")[!is.na(ccode)]
my[iso == 276 & year %between% c(1955, 1990), ccode := 260L]
my <- unique(my[, .(base_treaty, ccode, year)])

# Dyad-years sharing at least one agreement
pairs <- merge(my, my, by = c("base_treaty", "year"), allow.cartesian = TRUE)[ccode.x != ccode.y]
desta <- unique(pairs[, .(ccode1 = ccode.x, ccode2 = ccode.y, year)])
desta[, layer_pta_desta := 1]

# ---- Mansfield & Pevehouse PTAMID -------------------------------------------
pm <- as.data.table(read_dta(file.path(RAW, "PTAMIDTables.dta"), col_select = c("year", "dyadid", "PTASame")))
pm[, `:=`(c1 = as.integer(dyadid %/% 1000), c2 = as.integer(dyadid %% 1000))]
pm <- pm[!is.na(PTASame), .(pta = as.integer(max(PTASame) > 0)), by = .(c1, c2, year)]
pm <- unique(rbind(pm[, .(ccode1 = c1, ccode2 = c2, year, pta)], pm[, .(ccode1 = c2, ccode2 = c1, year, pta)]))
pm <- pm[, .(layer_pta_mansfield = max(pta)), by = .(ccode1, ccode2, year)]

# ---- Merge onto the panel ---------------------------------------------------
out <- merge(panel, desta, by = c("ccode1", "ccode2", "year"), all.x = TRUE)
out <- merge(out, pm, by = c("ccode1", "ccode2", "year"), all.x = TRUE)
cat(sprintf("DESTA dyad-years matched to panel: %d of %d\n", out[!is.na(layer_pta_desta), .N], nrow(desta)))
out[is.na(layer_pta_desta), layer_pta_desta := 0]
out[is.na(layer_pta_mansfield), layer_pta_mansfield := 0]
stopifnot(nrow(out) == nrow(panel), !anyDuplicated(out[, .(ccode1, ccode2, year)]))

# symmetry check (undirected by construction)
r <- out[, .(ccode1 = ccode2, ccode2 = ccode1, year, d_r = layer_pta_desta, m_r = layer_pta_mansfield)]
chk <- merge(out, r, by = c("ccode1", "ccode2", "year"))
stopifnot(all(chk$layer_pta_desta == chk$d_r), all(chk$layer_pta_mansfield == chk$m_r))

cov <- rbind(
  out[, .(share = mean(layer_pta_desta)), by = year][, layer := "pta_desta"],
  out[, .(share = mean(layer_pta_mansfield)), by = year][, layer := "pta_mansfield"])
fwrite(dcast(cov, year ~ layer, value.var = "share"), file.path(BOOK, "results", "main", "pta_coverage.csv"))
print(cov[year %in% c(1950, 1960, 1980, 2000, 2010)][order(layer, year)])

# Cross-check against PTAMID where both exist (1949-2000)
ov <- out[year %between% c(1949, 2000) & ccode1 < ccode2]
cat(sprintf("DESTA vs PTAMID agreement 1949-2000: %.3f; corr %.3f\n",
            mean(ov$layer_pta_desta == ov$layer_pta_mansfield), cor(ov$layer_pta_desta, ov$layer_pta_mansfield)))
saveRDS(out, file.path(BOOK, "data", "processed", "pta_layers.rds"), compress = "xz")
