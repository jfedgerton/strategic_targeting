###############################################################################
# 05_igo_mia_backfill.R   (Jared 2026-11-23; approved IGO robustness layer)
# Authority-weighted IGO ties from the Measure of International Authority (MIA v4, Hooghe, Marks et al.; 76 IOs, 1950-2019)
# and COW IGO v3 full memberships (state-year format; 1 = full member).
#   authority_jt = (delegation + pooling) / 2 for MIA IOs (0-1).
#   Backfill: MIA scores start in 1950; for IOs founded earlier the first MIA score is carried back to the first COW
#     membership year (MIA authority is sticky: 79% of IO-decades change < 0.05; r = 0.90 at ten years).
#   Pre-1950 IGOs absent from MIA (e.g., League of Nations): hand-coded drafts pending Jared's review are NOT used here
#     (authority 0) until approved; see data/processed/igo_mia_handcode_drafts.csv.
#   All other COW IGOs: authority 0 (technical/functional bodies).
#   Tie: w_ik,t = sum_j authority_jt * B_ij B_kj / (n_jt - 1), n_jt = full members of j in t.
# Output: data/processed/igo_mia.rds (ccode1 < ccode2, year, igo_mia raw sum), data/processed/igo_mia_map.csv
###############################################################################
suppressMessages({library(data.table); library(haven)})
set.seed(123)
BOOK <- Sys.getenv("BOOK"); stopifnot(nzchar(BOOK))
m <- as.data.table(read_dta(file.path(BOOK, "data", "raw", "mia", "MIA_v4_authority.dta")))
m[, auth := (delegation + pooling) / 2]
stopifnot(all(m$auth >= 0 & m$auth <= 1, na.rm = TRUE))
w <- fread(file.path(BOOK, "data", "raw", "cow_igo", "state_year_formatv3.csv"))
cow <- setdiff(names(w), c("ccode", "year", "state"))
# MIA acronym -> COW column: automatic normalized match, plus hand map for the 14 that differ
norm <- function(x) gsub("[^A-Z]", "", toupper(x))
io <- unique(m[, .(ionumber, ioname, acronym, inception)])
io[, cow := cow[match(norm(acronym), norm(cow))]]
hand <- c("LAIA/ALADI" = "LAIA", "OAU/ AU" = "OAU|AU", "Andean /CAN" = "Andean", "CIS/SNG" = "CIS", "COMECON" = "CMEA", "ECCAS-CEEAS" = "ECCAS",
          "GCC/CCASG" = "GCC", "IOM" = "IOMig", "ISA/ISBA" = "ISA", "Iwhale/IWC" = "IWhale", "OIF/ACCT" = "ACCT", "SICA/CAIS" = "SICA", "UNWTO" = "WTOURO")
io[acronym %in% names(hand), cow := hand[acronym]]
cat(sprintf("MIA IOs mapped to COW: %d of %d; unmapped: %s\n", io[!is.na(cow), .N], nrow(io), paste(io[is.na(cow), acronym], collapse = ", ")))
map <- io[!is.na(cow), .(cow = unlist(strsplit(cow, "\\|"))), by = .(ionumber, ioname, acronym, inception)]
stopifnot(all(map$cow %in% cow))
fwrite(map, file.path(BOOK, "data", "processed", "igo_mia_map.csv"))
# long membership (full members only)
L <- rbindlist(lapply(map$cow, function(v) w[get(v) == 1, .(ccode, year, igo = v)]))
L <- unique(L)
A <- merge(map[, .(ionumber, igo = cow)], m[, .(ionumber, year, auth)], by = "ionumber", allow.cartesian = TRUE)
A <- A[, .(auth = mean(auth, na.rm = TRUE)), by = .(igo, year)]
# backfill: years before the IO's first MIA year take its first MIA score
first <- A[!is.na(auth), .SD[which.min(year)], by = igo][, .(igo, y0 = year, a0 = auth)]
L <- merge(L, A, by = c("igo", "year"), all.x = TRUE)
L <- merge(L, first, by = "igo", all.x = TRUE)
L[is.na(auth) & year < y0, auth := a0]
cat(sprintf("membership-years backfilled before each IO's first MIA year: %d (IOs: %s)\n", L[year < y0, .N], paste(sort(unique(L[year < y0, igo])), collapse = " ")))
L[is.na(auth), auth := 0]   # e.g., years after MIA coverage ends (2020+ not in panel) or gaps
n <- L[, .(n = uniqueN(ccode)), by = .(igo, year)]
L <- merge(L, n, by = c("igo", "year"))[n >= 2]
# dyadic projection
P <- L[, { cc <- sort(unique(ccode)); if (length(cc) >= 2) { pr <- t(combn(cc, 2)); list(ccode1 = pr[, 1], ccode2 = pr[, 2], wt = auth[1] / (n[1] - 1)) } }, by = .(igo, year)]
D <- P[, .(igo_mia = sum(wt)), by = .(ccode1, ccode2, year)]
stopifnot(all(D$ccode1 < D$ccode2), all(D$igo_mia >= 0))
print(D[, .(dyads = .N, mean = round(mean(igo_mia), 4)), by = .(era = cut(year, c(1815, 1870, 1914, 1945, 1990, 2014)))][order(era)])
chk <- function(a, b, yr) { x <- D[ccode1 == min(a, b) & ccode2 == max(a, b) & year == yr, igo_mia]; sprintf("%d-%d %d: %.4f", a, b, yr, if (length(x)) x else 0) }
for (x in list(c(220, 260, 1960), c(2, 200, 1960), c(2, 365, 1960), c(365, 290, 1960), c(220, 255, 1900), c(2, 20, 1990), c(220, 260, 2000))) print(chk(x[1], x[2], x[3]))
saveRDS(D, file.path(BOOK, "data", "processed", "igo_mia.rds")); cat("saved", nrow(D), "\n")
# drafts for hand coding (not used): politically important pre-1950 COW IGOs absent from MIA
dr <- data.table(cow = c("LoN", "RepCom", "Danube", "IIA", "IACI"),
  name = c("League of Nations", "Reparation Commission", "Danube Commission", "International Institute of Agriculture", "COW IACI (Inter-American body; verify identity)"),
  proposed_delegation = NA_real_, proposed_pooling = NA_real_, status = "pending Jared review; coded with MIA codebook rules")
fwrite(dr, file.path(BOOK, "data", "processed", "igo_mia_handcode_drafts.csv"))
