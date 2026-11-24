###############################################################################
# 05b_igo_mia_fixed.R   (Jared 2026-11-23; approved IGO scheme on top of the F2 costly-ties network)
# Same MIA v4 x COW IGO v3 projection as 05_igo_mia_backfill.R, with the IGO audit fixes and scheme changes:
#   (1) Predecessor bodies inherit their successor's MIA score: EU -> ECSC, EEC, EURATOM, EU; OECD -> OEEC, OECD;
#       LAIA -> LAFTA, LAIA; CEMAC -> UDEAC, CEMAC; SADC -> SADCC, SADC; SICA -> OCAS, SICA; EAC1 -> EACM (was EAC);
#       OAU/AU -> OAU, AU. Membership of an MIA IO in year t = full member of ANY of its mapped COW bodies
#       (union, so overlapping bodies such as ECSC/EEC or OAU/AU in 2002 are counted once).
#   (2) COW -9 (missing) cells: last observed status of that state in that body carried forward; else 0.
#   (3) Germany 1990: COW codes the unified state 255 in 1990; the book frame also holds 260 and 265 in 1990.
#       260 in 1990 takes 255's 1990 memberships; 265 in 1990 takes its own 1989 memberships.
#   (4) Backfill cap: first MIA score carried back at most 10 years before the IO's first MIA year
#       (MIA authority stickiness evidence is r = 0.90 at ten years); earlier member-years get authority 0.
#   (5) Universal IGOs dropped: IO-years with full members > 2/3 of system states (book frame) in that year.
# Two tie weights (raw sums; the builder divides by the yearly max):
#   disc   = sum_j auth_jt / (n_jt - 1)   (size discount, as before)
#   nodisc = sum_j auth_jt                 (no size discount)
# Output: data/processed/igo_mia2.rds (ccode1 < ccode2, year, disc, nodisc); data/processed/igo_mia2_dropped_universal.csv
###############################################################################
suppressMessages({library(data.table); library(haven)})
set.seed(123)
BOOK <- Sys.getenv("BOOK"); stopifnot(nzchar(BOOK))
CAP <- as.integer(Sys.getenv("CAP", "10")); UNIV <- 2 / 3   # CAP=999 -> full backfill to first COW membership year
NOREG <- Sys.getenv("NOREG", "0") == "1"   # (6) Jared 2026-11-23: drop regional IOs (modal COW-region share of members >= 0.80, averaged over years) + Arab League
OUTF <- paste0(if (CAP == 10L) "igo_mia2" else sprintf("igo_mia2_cap%d", CAP), if (NOREG) "_noreg" else "", ".rds")
m <- as.data.table(read_dta(file.path(BOOK, "data", "raw", "mia", "MIA_v4_authority.dta")))
m[, auth := (delegation + pooling) / 2]
stopifnot(all(m$auth >= 0 & m$auth <= 1, na.rm = TRUE))
w <- fread(file.path(BOOK, "data", "raw", "cow_igo", "state_year_formatv3.csv"))
cow <- setdiff(names(w), c("ccode", "year", "state"))
norm <- function(x) gsub("[^A-Z]", "", toupper(x))
io <- unique(m[, .(ionumber, ioname, acronym, inception)])
io[, cow := cow[match(norm(acronym), norm(cow))]]
hand <- c("LAIA/ALADI" = "LAFTA|LAIA", "OAU/ AU" = "OAU|AU", "Andean /CAN" = "Andean", "CIS/SNG" = "CIS", "COMECON" = "CMEA", "ECCAS-CEEAS" = "ECCAS",
          "GCC/CCASG" = "GCC", "IOM" = "IOMig", "ISA/ISBA" = "ISA", "Iwhale/IWC" = "IWhale", "OIF/ACCT" = "ACCT", "SICA/CAIS" = "OCAS|SICA", "UNWTO" = "WTOURO",
          "EU" = "ECSC|EEC|EURATOM|EU", "OECD" = "OEEC|OECD", "CEMAC" = "UDEAC|CEMAC", "SADC" = "SADCC|SADC", "EAC1" = "EACM")
stopifnot(all(names(hand) %in% io$acronym))
io[acronym %in% names(hand), cow := hand[acronym]]
map <- io[!is.na(cow), .(cow = unlist(strsplit(cow, "\\|"))), by = .(ionumber, acronym)]
stopifnot(all(map$cow %in% cow), !anyDuplicated(map$cow))
cat(sprintf("MIA IOs mapped: %d of %d; COW bodies used: %d\n", uniqueN(map$ionumber), nrow(io), nrow(map)))
# (3) Germany 1990
if (!nrow(w[ccode == 260 & year == 1990])) w <- rbind(w, w[ccode == 255 & year == 1990][, ccode := 260L], fill = TRUE)
if (!nrow(w[ccode == 265 & year == 1990])) w <- rbind(w, w[ccode == 265 & year == 1989][, year := 1990L], fill = TRUE)
stopifnot(!anyDuplicated(w[, .(ccode, year)]))
# (2) -9 -> last observed status, else 0
setorder(w, ccode, year)
n9 <- 0L
for (v in unique(map$cow)) {
  x <- w[[v]]; n9 <- n9 + sum(x == -9, na.rm = TRUE); x[x == -9] <- NA
  set(w, j = v, value = x)
  w[, (v) := nafill(get(v), type = "locf"), by = ccode]
  set(w, which(is.na(w[[v]])), v, 0L)
}
cat(sprintf("-9 cells in mapped bodies replaced by last observed status: %d\n", n9))
# (1) membership by MIA IO = union over mapped COW bodies
L <- unique(rbindlist(lapply(seq_len(nrow(map)), function(i) w[get(map$cow[i]) == 1, .(ccode, year, ionumber = map$ionumber[i])])))
A <- m[, .(ionumber, year, auth)]
first <- A[!is.na(auth), .SD[which.min(year)], by = ionumber][, .(ionumber, y0 = year, a0 = auth)]
L <- merge(L, A, by = c("ionumber", "year"), all.x = TRUE)
L <- merge(L, first, by = "ionumber", all.x = TRUE)
# (4) capped backfill
L[is.na(auth) & year < y0 & year >= y0 - CAP, auth := a0]
cat(sprintf("member-years backfilled (cap %d years): %d; member-years before the cap set to 0: %d\n", CAP, L[year < y0 & year >= y0 - CAP, .N], L[year < y0 - CAP, .N]))
L[is.na(auth), auth := 0]
L <- L[auth > 0]
# (6) regional IOs (before the universal filter; uses union membership by MIA IO)
if (NOREG) {
  reg <- function(cc) fifelse(cc < 200, "Americas", fifelse(cc < 400, "Europe", fifelse(cc < 627, "Africa", fifelse(cc < 700, "MidEast", "Asia-Oceania"))))
  rs <- L[, .(share = max(table(reg(ccode))) / .N), by = .(ionumber, year)][, .(modal_share = mean(share)), by = ionumber]
  rs <- merge(rs, io[, .(ionumber, acronym)], by = "ionumber")
  regional <- rs[modal_share >= 0.80 | acronym == "LOAS"]
  cat(sprintf("regional IOs dropped (%d): %s\n", nrow(regional), paste(sort(regional$acronym), collapse = ", ")))
  stopifnot("LOAS" %in% regional$acronym, "EU" %in% regional$acronym, !("OECD" %in% regional$acronym))
  fwrite(rs[order(-modal_share)], file.path(BOOK, "data", "processed", "igo_mia2_regional_share.csv"))
  L <- L[!ionumber %in% regional$ionumber]
}
# (5) universal IOs
fr <- readRDS(file.path(BOOK, "data", "processed", "layers_dyad_year.rds")); setDT(fr)
Ns <- fr[, .(N = uniqueN(ccode1)), by = year]
n <- L[, .(n = uniqueN(ccode)), by = .(ionumber, year)]
n <- merge(n, Ns, by = "year")
drop <- n[n > UNIV * N]
dd <- merge(drop, io[, .(ionumber, acronym)], by = "ionumber")[, .(years = paste(range(year), collapse = "-"), n_years = .N), by = acronym][order(-n_years)]
cat("universal IO-years dropped (members > 2/3 of states):\n"); print(dd)
fwrite(dd, file.path(BOOK, "data", "processed", "igo_mia2_dropped_universal.csv"))
L <- merge(L, n[n <= UNIV * N & n >= 2, .(ionumber, year, n)], by = c("ionumber", "year"))
# dyadic projection
P <- L[, { cc <- sort(unique(ccode)); pr <- t(combn(cc, 2)); list(ccode1 = pr[, 1], ccode2 = pr[, 2], a = auth[1], n = n[1]) }, by = .(ionumber, year)]
D <- P[, .(disc = sum(a / (n - 1)), nodisc = sum(a)), by = .(ccode1, ccode2, year)]
stopifnot(all(D$ccode1 < D$ccode2), all(D$disc > 0), all(D$nodisc >= D$disc), !anyDuplicated(D[, .(ccode1, ccode2, year)]))
print(D[, .(dyads = .N, disc = round(mean(disc), 4), nodisc = round(mean(nodisc), 3)), by = .(era = cut(year, c(1815, 1914, 1945, 1990, 2014)))][order(era)])
chk <- function(a, b, yr) { x <- D[ccode1 == min(a, b) & ccode2 == max(a, b) & year == yr]; sprintf("%d-%d %d: disc %.4f nodisc %.3f", a, b, yr, if (nrow(x)) x$disc else 0, if (nrow(x)) x$nodisc else 0) }
for (x in list(c(220, 325, 1970), c(220, 325, 1992), c(220, 325, 1993), c(220, 260, 1990), c(2, 200, 1960), c(2, 365, 1960), c(365, 290, 1960), c(2, 70, 1960)))
  print(chk(x[1], x[2], x[3]))
saveRDS(D, file.path(BOOK, "data", "processed", OUTF)); cat("saved", nrow(D), "\n")
