## ============================================================================
## 20_build_layers.R  --  Part 1 of the relational-measure sweep
## Builds ONE undirected dyad-year panel (ccode1 < ccode2, COW system members, 1816-2014) with one column
## per cooperative layer, plus contiguity and capital distance for later residualisation.
##   igo   : shared COW IGO memberships (peacesciencer cow_igo_ndy). Raw; years with no IGO data stay NA.
##   ally  : any ATOP obligation (defense/offense/neutral/nonagg/consul), binary (peacesciencer atop_alliance; descriptive only)
##   def, off, neu, nonagg, consul : the five ATOP obligation types, each binary (PI 2026-11-18: nine-layer grid)
##   un    : UN General Assembly voting agreement, Cohen's kappa (Hage 2011 via peacesciencer dyadic_fp_similarity, kappavv),
##           1946+; before 1946 a STRUCTURAL zero (no UN); in UN years a dyad without a kappa (non-member, not voting) -> NA
##   dipl  : COW diplomatic exchange DE (both sides represented), binary; 5-year waves carried forward
##           to the year before the next wave (1816 takes the 1817 wave)
##   pta   : count of DESTA PTAs in force (entry-into-force year, else signing year; no exits recorded)
##   dca   : Kinne DCAD dcaAnyV1 (any DCA in force), binary; STRUCTURAL zeros before 1980; 2010 carried forward to 2014
##   trade : COW trade, built but NOT in the grid (PI 2026-11-18: skipped for missingness)
## *_src columns record the source year used, so staleness is auditable.
## Usage: Rscript v2/R/20_build_layers.R <netdata_dir> <out_rds>
## ============================================================================
suppressMessages({ library(dplyr); library(tidyr); library(peacesciencer); library(countrycode) })
set.seed(123)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
NETDIR <- args[1]; OUT <- args[2]
FIRST_YEAR <- 1816L; LAST_YEAR <- 2014L
logf <- function(...) cat(sprintf(...), "\n")
logf("R %s | peacesciencer %s | countrycode %s", getRversion(), packageVersion("peacesciencer"), packageVersion("countrycode"))

## ---- 0. dyad-year frame with contiguity and capital distance ----------------
dy <- create_dyadyears(system = "cow", directed = FALSE, mry = FALSE, subset_years = FIRST_YEAR:LAST_YEAR) |>
  add_contiguity() |> add_capital_distance() |>
  transmute(ccode1 = as.integer(ccode1), ccode2 = as.integer(ccode2), year = as.integer(year),
            contig_land = as.integer(conttype == 1), conttype = as.integer(conttype), capdist = capdist)
stopifnot(all(dy$ccode1 < dy$ccode2), !anyDuplicated(dy[, c("ccode1", "ccode2", "year")]))
logf("frame: %d dyad-years, %d-%d; capdist NA %d; conttype NA %d", nrow(dy), min(dy$year), max(dy$year),
     sum(is.na(dy$capdist)), sum(is.na(dy$conttype)))

## ---- 1. IGO -------------------------------------------------------------------
igo <- peacesciencer::cow_igo_ndy |>
  transmute(ccode1 = as.integer(pmin(ccode1, ccode2)), ccode2 = as.integer(pmax(ccode1, ccode2)), year = as.integer(year),
            igo = ifelse(dyadigos < 0, NA_real_, as.numeric(dyadigos))) |>
  distinct(ccode1, ccode2, year, .keep_all = TRUE)
dy <- dy |> left_join(igo, by = c("ccode1", "ccode2", "year"))
igo_years <- dy |> group_by(year) |> summarise(ok = any(!is.na(igo)), .groups = "drop")
logf("IGO: years with any data %d of %d; all-missing years: %s", sum(igo_years$ok), nrow(igo_years),
     paste(igo_years$year[!igo_years$ok], collapse = " "))
## within an IGO-data year, a missing dyad has no recorded shared membership -> 0 (matches 01_build_data.R)
dy <- dy |> left_join(igo_years, by = "year") |> mutate(igo = ifelse(ok & is.na(igo), 0, igo)) |> select(-ok)

## ---- 2. ATOP ------------------------------------------------------------------
atop <- peacesciencer::atop_alliance |>
  transmute(ccode1 = as.integer(pmin(ccode1, ccode2)), ccode2 = as.integer(pmax(ccode1, ccode2)), year = as.integer(year),
            ally = as.integer((atop_defense + atop_offense + atop_neutral + atop_nonagg + atop_consul) > 0),
            def = as.integer(atop_defense > 0), off = as.integer(atop_offense > 0), neu = as.integer(atop_neutral > 0),
            nonagg = as.integer(atop_nonagg > 0), consul = as.integer(atop_consul > 0)) |>
  group_by(ccode1, ccode2, year) |>
  summarise(across(c(ally, def, off, neu, nonagg, consul), max), .groups = "drop")
dy <- dy |> left_join(atop, by = c("ccode1", "ccode2", "year")) |>
  mutate(across(c(ally, def, off, neu, nonagg, consul), ~ coalesce(.x, 0L)))
stopifnot(all(dy$ally == pmax(dy$def, dy$off, dy$neu, dy$nonagg, dy$consul)))

## ---- 2b. UN voting agreement (kappa), 1946+ ---------------------------------------
fp <- readRDS(file.path(system.file("extdata", package = "peacesciencer"), "dyadic_fp_similarity.rds")) |>
  transmute(ccode1 = as.integer(pmin(ccode1, ccode2)), ccode2 = as.integer(pmax(ccode1, ccode2)), year = as.integer(year), un = kappavv) |>
  filter(!is.na(un)) |> group_by(ccode1, ccode2, year) |> summarise(un = mean(un), .groups = "drop")
un_first <- min(fp$year)
stopifnot(un_first >= 1945L, all(fp$un >= -1 - 1e-9 & fp$un <= 1 + 1e-9))
dy <- dy |> left_join(fp, by = c("ccode1", "ccode2", "year")) |> mutate(un = ifelse(year < un_first, 0, un))
logf("UN kappa: first year %d (structural 0 before); NA share in UN years %.3f", un_first, mean(is.na(dy$un[dy$year >= un_first])))

## ---- 3. Diplomatic exchange (5-year waves, carried forward) -------------------
dip <- read.csv(file.path(NETDIR, "dipl_dyad.csv")) |>
  transmute(ccode1 = as.integer(pmin(ccode1, ccode2)), ccode2 = as.integer(pmax(ccode1, ccode2)), year = as.integer(year),
            dipl = as.integer(DE == 1)) |>
  group_by(ccode1, ccode2, year) |> summarise(dipl = max(dipl), .groups = "drop")
waves <- sort(unique(dip$year))
stopifnot(all(dip$dipl %in% 0:1))
wave_of <- function(y) waves[pmax(1L, findInterval(y, waves))]   # 1816 -> 1817 wave (index 0 would be silently dropped)
dy <- dy |> mutate(dipl_src = wave_of(year))
stopifnot(!anyNA(dy$dipl_src), all(dy$dipl_src <= pmax(dy$year, min(waves))), all(dy$year - dy$dipl_src < 10))
dy <- dy |>
  left_join(dip |> rename(dipl_src = year), by = c("ccode1", "ccode2", "dipl_src")) |>
  mutate(dipl = coalesce(dipl, 0L))
logf("dipl: waves %s; max staleness %d years", paste(range(waves), collapse = "-"), max(dy$year - dy$dipl_src))

## ---- 4. DESTA PTAs in force ---------------------------------------------------
load(file.path(NETDIR, "destaRaw.rda"))   # object `desta`
pt <- desta |> mutate(c1 = countrycode(iso1, "iso3n", "cown", warn = FALSE),
                      c2 = countrycode(iso2, "iso3n", "cown", warn = FALSE),
                      start = as.integer(ifelse(is.na(entryforceyear), year, entryforceyear))) |>
  mutate(c1 = ifelse(country1 == "Serbia", 345L, ifelse(country1 == "Kosovo", 347L, c1)),     # COW Yugoslavia/Serbia, Kosovo
         c2 = ifelse(country2 == "Serbia", 345L, ifelse(country2 == "Kosovo", 347L, c2)))
unm <- pt |> filter(is.na(c1) | is.na(c2))
logf("DESTA rows %d; rows with an unmatched country (dropped) %d; unmatched names: %s", nrow(pt), nrow(unm),
     paste(sort(unique(c(unm$country1[is.na(unm$c1)], unm$country2[is.na(unm$c2)]))), collapse = "; "))
pt <- pt |> filter(!is.na(c1), !is.na(c2), c1 != c2) |>
  transmute(ccode1 = as.integer(pmin(c1, c2)), ccode2 = as.integer(pmax(c1, c2)), number, start) |>
  group_by(ccode1, ccode2, number) |> summarise(start = min(start), .groups = "drop")   # one entry per treaty-dyad
pta <- dy |> distinct(ccode1, ccode2, year) |> inner_join(pt, by = c("ccode1", "ccode2"), relationship = "many-to-many") |>
  filter(start <= year) |> count(ccode1, ccode2, year, name = "pta")
dy <- dy |> left_join(pta, by = c("ccode1", "ccode2", "year")) |> mutate(pta = coalesce(pta, 0L))

## ---- 5. DCA ---------------------------------------------------------------------
dca <- read.csv(file.path(NETDIR, "dca_dyad.csv")) |>
  transmute(ccode1 = as.integer(pmin(ccode1, ccode2)), ccode2 = as.integer(pmax(ccode1, ccode2)), year = as.integer(year),
            dca = as.integer(dcaAnyV1 > 0)) |>
  group_by(ccode1, ccode2, year) |> summarise(dca = max(dca), .groups = "drop")
dca_rng <- range(dca$year)
stopifnot(dca_rng[1] == 1980L, dca_rng[2] == 2010L)
## PI decision 2026-11-18: pre-1980 zeros are STRUCTURAL (no DCAs existed); 2011-2014 carry the 2010 value forward
dca_last <- dca |> filter(year == dca_rng[2]) |> select(ccode1, ccode2, dca_last = dca)
dy <- dy |> left_join(dca, by = c("ccode1", "ccode2", "year")) |>
  left_join(dca_last, by = c("ccode1", "ccode2")) |>
  mutate(dca_src = ifelse(year < dca_rng[1], NA_integer_, pmin(year, dca_rng[2])),
         dca = ifelse(year < dca_rng[1], 0L, ifelse(year <= dca_rng[2], coalesce(dca, 0L), coalesce(dca_last, 0L)))) |>
  select(-dca_last)
logf("DCA: structural zeros before %d; %d carried forward to %d", dca_rng[1], dca_rng[2], LAST_YEAR)

## ---- 5b. COW bilateral trade (1870-2014; NA before 1870 = true missingness) -------------------
## log(1 + smoothflow1 + smoothflow2), millions current USD. Missing flows (NA in COW, or dyad absent) stay NA.
## COW itself records some implausible zeros (e.g. USSR-Poland 1960); those are kept as reported.
## PI decision 2026-11-18: any spec that uses the trade layer is subset to 1870+.
trd <- readRDS(file.path(system.file("extdata", package = "peacesciencer"), "cow_trade_ddy.rds")) |>
  transmute(ccode1 = as.integer(pmin(ccode1, ccode2)), ccode2 = as.integer(pmax(ccode1, ccode2)), year = as.integer(year),
            flow = ifelse(is.na(smoothflow1) & is.na(smoothflow2), NA_real_,
                          pmax(coalesce(smoothflow1, 0), 0) + pmax(coalesce(smoothflow2, 0), 0))) |>
  group_by(ccode1, ccode2, year) |> summarise(flow = max(flow), .groups = "drop")
trd_rng <- range(trd$year)
stopifnot(trd_rng[1] == 1870L)
dy <- dy |> left_join(trd, by = c("ccode1", "ccode2", "year")) |>
  mutate(trade = ifelse(year < trd_rng[1], NA_real_, log1p(flow))) |> select(-flow)   # missing flows stay NA (no fabricated zeros)
logf("trade: 1870+ dyad-years with NA trade (missing in COW or absent): %d (%.3f)", sum(is.na(dy$trade) & dy$year >= trd_rng[1]),
     mean(is.na(dy$trade[dy$year >= trd_rng[1]])))
logf("trade: %d-%d, NA before %d; share of 1870+ dyad-years with positive trade %.3f", trd_rng[1], trd_rng[2], trd_rng[1],
     mean(dy$trade[dy$year >= trd_rng[1]] > 0, na.rm = TRUE))

## ---- 6. checks -------------------------------------------------------------------
stopifnot(!anyDuplicated(dy[, c("ccode1", "ccode2", "year")]))
era <- cut(dy$year, c(1815, 1870, 1914, 1945, 1979, 1990, 2014), labels = c("1816-70", "1871-1914", "1915-45", "1946-79", "1980-90", "1991-2014"))
dens <- dy |> mutate(era = era) |> group_by(era) |>
  summarise(dyad_years = n(), igo_mean = mean(igo, na.rm = TRUE), igo_na = mean(is.na(igo)), ally = mean(ally), def = mean(def),
            off = mean(off), neu = mean(neu), nonagg = mean(nonagg), consul = mean(consul), un_mean = mean(un, na.rm = TRUE), un_na = mean(is.na(un)),
            dipl = mean(dipl), pta_any = mean(pta > 0), dca = mean(dca), trade_pos = mean(trade > 0, na.rm = TRUE), .groups = "drop")
print(as.data.frame(dens), digits = 3, width = 200)
chk <- function(a, b, y) dy |> filter(ccode1 == min(a, b), ccode2 == max(a, b), year == y) |>
  transmute(pair = sprintf("%d-%d", a, b), year, igo, def, off, neu, nonagg, consul, un = round(un, 2), dipl, pta, dca)
face <- bind_rows(chk(2, 200, 1955), chk(2, 365, 1925), chk(2, 365, 1960), chk(365, 290, 1960), chk(200, 255, 1960),
                  chk(220, 255, 1995), chk(2, 840, 1995), chk(2, 70, 2000), chk(2, 365, 1935), chk(200, 220, 1880), chk(2, 840, 2012), chk(2, 20, 1990))
print(as.data.frame(face))
saveRDS(dy, OUT)
logf("saved %s (%d rows, %d cols)", OUT, nrow(dy), ncol(dy))
