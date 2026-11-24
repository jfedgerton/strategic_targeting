###############################################################################
# 01_build_measure.R   (multiverse, Jared 2026-11-17)
# One tie measure per combination of
#   LAYERS : alliance  (defensive OR offensive pledge, binary)
#            igo       (shared IGO memberships)
#            allyigo   (alliance + igo; current spine)
#            main      (defensive, offensive, igo, trade [from 1870])
#            trade     (trade only; panel starts 1870)                     [added 2026-11-17]
#            allytrade (alliance + trade)                                   [added 2026-11-17]
#            allyigotrade (alliance + igo + trade)                          [added 2026-11-17]
#            LIO-era layers (Jared 2026-11-17): DCA (Kinne DCAD, 1980-2010) and PTA (DESTA, 1950s on)
#            treated like offensive alliances -- years before the instrument existed are true zeros,
#            not missing. DCAD ends in 2010: each dyad's 2010 DCA status is carried forward to 2011-14
#            (assumption; DCAs are multi-year agreements). Sets:
#            dca, pta, allydca (alliance + DCA), tradepta (trade + PTA),
#            allyigodcapta (alliance + igo + DCA + PTA), full (alliance + DCA + igo + trade + PTA)
#   IGO    : dyadic    (layer_igo = shared full memberships, min-max scaled within year)
#            hyper     (IGOs as hyperedges, size-weighted projection, min-max scaled within year;
#                       00_igo_hyperedge_projection.R). Not used for LAYERS = alliance.
#   WEIGHT : sry       surprise y - y_hat, expected model with sender + receiver + year FE
#            sy        surprise y - y_hat, sender x year + receiver x year FE
#            hybrid    sry FE; binary: present 1, absent -p_hat; continuous: y if y >= y_hat else y - y_hat
#            raw       observed layer values, no expected-tie model (no negative ties)
#   LAYERS = custom (Jared 2026-11-21, nine-layer grid): LAYER_SET = "+"-separated subset of
#            igo, def, off, neu, nonagg, consul (ATOP types, 0/1), dca, pta, un (UN kappa, 1946+), trade.
#            ATOP types and UN from data/processed/atop_un_layers.rds (01_data/03_atop_un_layers.R).
#            UN is continuous and signed: expected tie by OLS (same FE and log distance), score y - y_hat;
#            years before 1946 and pairs without a kappa (non-members) get score 0. Neutrality,
#            non-aggression and consultation: years with no such tie anywhere are true zeros (score 0).
#            WEIGHT = w2contig (Jared 2026-11-21): the targeting project's W2contig, copied from
#            relational_measures/02_build_positions.R. igo and pta as log(1 + count); per layer, residual of
#            v ~ contig_land + year (pooled OLS over existence years with a nonzero tie), clipped at 0, divided by
#            the layer-year maximum; combo = sqrt(sum of squares) rescaled to [0, 1] within year. No negative ties.
#            UN without kappa = no tie. PTA here is the book's binary DESTA layer (no counts).
#            Weightings sry, sy, raw (raw: observed layer values; raw UN kappa keeps its sign).
# Expected-tie model: logit (binary layers) or Poisson pseudo-likelihood (continuous), log distance.
# GEO (Jared 2026-11-23): geography term in the expected-tie model. dist (default) = log(1 + capital
#   distance); contig12 = binary contiguity, COW conttype 1 or 2 (land, or water up to 12 miles).
# Composite, symmetrisation, tie_norm and tie_pos exactly as in 01_main_measure.R.
# Output (in OUT_ROOT, a per-specification book root):
#   data/processed/main_measure_allyigo.rds and main_measure_edges_allyigo.csv.gz
#   (the "allyigo" file names are kept so the chapter scripts run unchanged in OUT_ROOT)
###############################################################################
suppressMessages({library(data.table); library(fixest)})
set.seed(123)
setFixest_nthreads(1)
BOOK <- Sys.getenv("BOOK_ROOT", getwd())
OUT <- Sys.getenv("OUT_ROOT")
LAYERS <- Sys.getenv("LAYERS"); IGO <- Sys.getenv("IGO"); WEIGHT <- Sys.getenv("WEIGHT"); LAYER_SET <- Sys.getenv("LAYER_SET")
stopifnot(nzchar(OUT), LAYERS %in% c("alliance", "igo", "allyigo", "main", "trade", "allytrade", "allyigotrade",
                                     "dca", "pta", "allydca", "tradepta", "allyigodcapta", "full", "custom"),
          IGO %in% c("dyadic", "hyper", "mia", "mia2", "mia2nd", "var", "none"), WEIGHT %in% c("sry", "sy", "hybrid", "raw", "w2contig", "igofe"))
if (LAYERS == "custom") {
  CS <- strsplit(LAYER_SET, "+", fixed = TRUE)[[1]]
  stopifnot(length(CS) >= 1, !anyDuplicated(CS), all(CS %in% c("igo", "def", "off", "neu", "nonagg", "consul", "dca", "pta", "ptai", "ptad", "ptadeep", "un", "trade")),
            any(c("igo", "def") %in% CS), ("igo" %in% CS) == (IGO != "none"), WEIGHT %in% c("sry", "sy", "raw", "w2contig", "igofe"))
} else stopifnot((LAYERS %in% c("alliance", "trade", "allytrade", "dca", "pta", "allydca", "tradepta")) == (IGO == "none"))
cat(sprintf("== LAYERS=%s IGO=%s WEIGHT=%s LAYER_SET=%s\n", LAYERS, IGO, WEIGHT, LAYER_SET))

d <- readRDS(file.path(BOOK, "data", "processed", "layers_dyad_year.rds")); setDT(d)
stopifnot(!anyDuplicated(d[, .(ccode1, ccode2, year)]))
d[, layer_alliance := pmax(layer_defensive, layer_offensive)]
# DCA: DCAD ends in 2010. DCA_MODE = "carry" (default): each dyad's 2010 status carried forward
#   to 2011-2014 (used for communities, conflict, diffusion, convergence). DCA_MODE = "trunc2010":
#   no carry-forward and the panel ends in 2010 (used for churn, so carried ties cannot create
#   artificial zero churn) -- Jared 2026-11-17.
DCA_MODE <- Sys.getenv("DCA_MODE", "carry")
stopifnot(DCA_MODE %in% c("carry", "trunc2010"))
stopifnot(d[year > 2010, all(layer_dca == 0)], all(d$layer_dca %in% c(0, 1)))
if (DCA_MODE == "carry") {
  d10 <- d[year == 2010, .(ccode1, ccode2, dca10 = layer_dca)]
  d <- merge(d, d10, by = c("ccode1", "ccode2"), all.x = TRUE)
  d[year > 2010, layer_dca := fifelse(is.na(dca10), 0, dca10)][, dca10 := NULL]
  cat(sprintf("  DCA dyad-years carried forward to 2011-14: %d\n", d[year > 2010, sum(layer_dca)]))
} else {
  d <- d[year <= 2010]
  cat("  DCA_MODE trunc2010: panel ends in 2010\n")
}
P <- readRDS(file.path(BOOK, "data", "processed", "pta_layers.rds")); setDT(P)
n0 <- nrow(d)
d <- merge(d, P[, .(ccode1, ccode2, year, layer_pta = layer_pta_desta)], by = c("ccode1", "ccode2", "year"), all.x = TRUE)
stopifnot(nrow(d) == n0, !anyNA(d$layer_pta), all(d$layer_pta %in% c(0, 1)))
# Costly PTA layers from DESTA 2.3 (01_data/04_pta_desta_costly.R; Jared 2026-11-23): ptai = in force (0/1),
#   ptad = max depth of in-force agreements / 7 (continuous 0-1), ptadeep = max depth >= 4 (0/1)
if (file.exists(file.path(BOOK, "data", "processed", "pta_costly.rds"))) {
  PC <- readRDS(file.path(BOOK, "data", "processed", "pta_costly.rds")); setDT(PC)
  n0 <- nrow(d)
  d <- merge(d, PC[, .(ccode1, ccode2, year, layer_ptai = layer_pta_inforce, layer_ptad = layer_pta_depth, layer_ptadeep = layer_pta_deep)], by = c("ccode1", "ccode2", "year"), all.x = TRUE)
  stopifnot(nrow(d) == n0)
  for (v in c("layer_ptai", "layer_ptad", "layer_ptadeep")) set(d, which(is.na(d[[v]])), v, 0)
}
# PTA fix variants (01_data/04b_pta_variants.R; Jared 2026-11-23): PTA_VARIANT in e2, e4, f1, f2 replaces layer_ptad
PTA_VARIANT <- Sys.getenv("PTA_VARIANT")
if (nzchar(PTA_VARIANT)) {
  stopifnot(PTA_VARIANT %in% c("e2", "e4", "f1", "f2"))
  PV <- readRDS(file.path(BOOK, "data", "processed", "pta_variants.rds")); setDT(PV)
  n0 <- nrow(d); d[, layer_ptad := NULL]
  d <- merge(d, PV[, .(ccode1, ccode2, year, layer_ptad = get(paste0("ptad_", PTA_VARIANT)))], by = c("ccode1", "ccode2", "year"), all.x = TRUE)
  stopifnot(nrow(d) == n0); set(d, which(is.na(d$layer_ptad)), "layer_ptad", 0)
  stopifnot(all(d$layer_ptad >= 0 & d$layer_ptad <= 1))
  cat("PTA_VARIANT", PTA_VARIANT, ": share of dyad-years with ptad > 0 =", round(mean(d$layer_ptad > 0), 4), "\n")
}
if (LAYERS == "custom") {
  AU <- readRDS(file.path(BOOK, "data", "processed", "atop_un_layers.rds")); setDT(AU)
  n0 <- nrow(d)
  d <- merge(d, AU, by = c("ccode1", "ccode2", "year"), all.x = TRUE)
  stopifnot(nrow(d) == n0, !anyNA(d[, .(layer_def, layer_off, layer_neu, layer_nonagg, layer_consul)]),
            all(d$layer_def == d$layer_defensive), all(d$layer_off == d$layer_offensive))
}
if (IGO == "mia") {   # MIA authority-weighted IGO ties with backfill (01_data/05_igo_mia_backfill.R; Jared 2026-11-23)
  H <- readRDS(file.path(BOOK, "data", "processed", "igo_mia.rds")); setDT(H)
  H[, igo01 := { r <- max(igo_mia); if (r > 0) igo_mia / r else 0 }, by = year]
  H2 <- rbind(H[, .(ccode1, ccode2, year, igo01)], H[, .(ccode1 = ccode2, ccode2 = ccode1, year, igo01)])
  n0 <- nrow(d)
  d <- merge(d, H2, by = c("ccode1", "ccode2", "year"), all.x = TRUE)
  stopifnot(nrow(d) == n0)
  d[is.na(igo01), igo01 := 0]
  d[, layer_igo := igo01][, igo01 := NULL]
}
if (IGO %in% c("mia2", "mia2nd")) {   # fixed MIA layer (01_data/05b_igo_mia_fixed.R; Jared 2026-11-23): mia2 = size-discounted, mia2nd = no discount
  MIA2_CAP <- Sys.getenv("MIA2_CAP", "10"); MIA2_NOREG <- Sys.getenv("MIA2_NOREG", "0") == "1"
  H <- readRDS(file.path(BOOK, "data", "processed", paste0(if (MIA2_CAP == "10") "igo_mia2" else paste0("igo_mia2_cap", MIA2_CAP), if (MIA2_NOREG) "_noreg" else "", ".rds"))); setDT(H)
  cat("MIA2 file: cap", MIA2_CAP, "noreg", MIA2_NOREG, "\n")
  H[, raw := if (IGO == "mia2") disc else nodisc]
  H[, igo01 := { r <- max(raw); if (r > 0) raw / r else 0 }, by = year]
  H2 <- rbind(H[, .(ccode1, ccode2, year, igo01)], H[, .(ccode1 = ccode2, ccode2 = ccode1, year, igo01)])
  n0 <- nrow(d)
  d <- merge(d, H2, by = c("ccode1", "ccode2", "year"), all.x = TRUE)
  stopifnot(nrow(d) == n0)
  d[is.na(igo01), igo01 := 0]
  d[, layer_igo := igo01][, igo01 := NULL]
  IGO_WT <- as.numeric(Sys.getenv("IGO_WT", "1")); stopifnot(IGO_WT > 0, IGO_WT <= 1)
  d[, layer_igo := layer_igo * IGO_WT]   # down-weighted IGO layer in the composite mean
  cat("IGO", IGO, "IGO_WT", IGO_WT, ": share of dyad-years with igo > 0 =", round(mean(d$layer_igo > 0), 4), "\n")
}
if (IGO == "var") {   # IGO construction variants (01_data/06_igo_variants.R; Jared 2026-11-23): IGO_VAR in count, log, cos, hyper; min-max within year
  IGO_VAR <- Sys.getenv("IGO_VAR"); stopifnot(IGO_VAR %in% c("count", "log", "cos", "hyper"))
  H <- readRDS(file.path(BOOK, "data", "processed", "igo_variants.rds")); setDT(H)
  H[, raw := get(paste0("igo_", IGO_VAR))]
  H[, igo01 := { r <- max(raw) - min(raw); if (r > 0) (raw - min(raw)) / r else 0 }, by = year]
  H2 <- rbind(H[, .(ccode1, ccode2, year, igo01)], H[, .(ccode1 = ccode2, ccode2 = ccode1, year, igo01)])
  n0 <- nrow(d); d[, layer_igo := NULL]
  d <- merge(d, H2, by = c("ccode1", "ccode2", "year"), all.x = TRUE)
  stopifnot(nrow(d) == n0, !anyNA(d$igo01))
  d[, layer_igo := igo01][, igo01 := NULL]
  cat("IGO var", IGO_VAR, ": share > 0 =", round(mean(d$layer_igo > 0), 4), "\n")
}
if (IGO == "hyper") {
  H <- readRDS(file.path(BOOK, "data", "processed", "igo_hyper_size.rds"))
  H[, igo01 := { r <- max(igo_hyper) - min(igo_hyper); if (r > 0) (igo_hyper - min(igo_hyper)) / r else 0 }, by = year]
  H2 <- rbind(H[, .(ccode1, ccode2, year, igo01)], H[, .(ccode1 = ccode2, ccode2 = ccode1, year, igo01)])
  n0 <- nrow(d)
  d <- merge(d, H2, by = c("ccode1", "ccode2", "year"), all.x = TRUE)
  stopifnot(nrow(d) == n0, !anyNA(d$igo01))
  d[, layer_igo := igo01][, igo01 := NULL]
}
layers <- switch(LAYERS,
  alliance = "layer_alliance",
  igo      = "layer_igo",
  allyigo  = c("layer_alliance", "layer_igo"),
  main     = c("layer_defensive", "layer_offensive", "layer_igo", "layer_trade"),
  trade    = "layer_trade",
  allytrade = c("layer_alliance", "layer_trade"),
  allyigotrade = c("layer_alliance", "layer_igo", "layer_trade"),
  dca      = "layer_dca",
  pta      = "layer_pta",
  allydca  = c("layer_alliance", "layer_dca"),
  tradepta = c("layer_trade", "layer_pta"),
  allyigodcapta = c("layer_alliance", "layer_igo", "layer_dca", "layer_pta"),
  full     = c("layer_alliance", "layer_dca", "layer_igo", "layer_trade", "layer_pta"),
  custom   = paste0("layer_", CS))
binary <- intersect(layers, c("layer_alliance", "layer_defensive", "layer_offensive", "layer_dca", "layer_pta", "layer_ptai", "layer_ptadeep",
                              "layer_def", "layer_off", "layer_neu", "layer_nonagg", "layer_consul"))
for (v in binary) stopifnot(all(d[[v]] %in% c(0, 1)))
d[, ldist := log1p(capdist)]
stopifnot(!anyNA(d$ldist))
GEO <- Sys.getenv("GEO", "dist"); stopifnot(GEO %in% c("dist", "contig12"))
d[, contig12 := as.integer(!is.na(conttype) & conttype %in% 1:2)]
GEOVAR <- if (GEO == "dist") "ldist" else "contig12"
cat(sprintf("  geography term in expected-tie model: %s (contig12 share %.4f)\n", GEOVAR, mean(d$contig12)))
FE_TXT <- if (WEIGHT == "sy") "ccode1^year + ccode2^year" else "ccode1 + ccode2 + year"

if (WEIGHT == "w2contig") {
  d[, contig_land := as.integer(!is.na(conttype) & conttype == 1L)]
  range01 <- function(x) { r <- range(x); if (r[2] - r[1] < 1e-12) rep(0, length(x)) else (x - r[1]) / (r[2] - r[1]) }
  yr_split <- split(seq_len(nrow(d)), d$year)
  EL <- matrix(0, nrow(d), length(layers), dimnames = list(NULL, layers))
  for (v in layers) {
    x0 <- if (v == "layer_igo") log1p(d$igo_count) else as.numeric(d[[v]])
    if (v == "layer_igo") stopifnot(!anyNA(x0))
    if (v != "layer_un") stopifnot(!anyNA(x0))
    first_year <- min(d$year[!is.na(x0) & x0 != 0])
    yr_has <- d[, .(h = any(!is.na(x0[.I]) & x0[.I] != 0)), by = year]
    ok <- d$year >= first_year & d$year %in% yr_has[h == TRUE, year] & !is.na(x0)
    r <- resid(lm(v ~ contig_land + year, data = data.frame(v = x0[ok], contig_land = d$contig_land[ok], year = d$year[ok])))
    x <- numeric(nrow(d)); x[ok] <- r
    tol <- 1e-8 * max(1, max(abs(x0[ok]))); x[abs(x) < tol] <- 0; x <- pmax(x, 0)
    e <- numeric(nrow(d))
    for (ix in yr_split) { j <- ix[ok[ix]]; if (length(j)) { m <- max(x[j]); if (m > 0) e[j] <- x[j] / m } }
    EL[, v] <- e
    cat(sprintf("  %-16s W2contig: first year %d, mean weight %.3f\n", v, first_year, mean(e)))
  }
  stopifnot(!anyNA(EL), all(EL >= 0), all(EL <= 1 + 1e-12))
  S <- sqrt(rowSums(EL^2)); E <- numeric(length(S)); for (ix in yr_split) E[ix] <- range01(S[ix])
  d[, tie_dir := E]; d[, n_layers := length(layers)]
  score_cols <- character(0)
} else {
score_cols <- sub("layer_", "score_", layers)
for (v in layers) {
  sname <- sub("layer_", "score_", v)
  d[, (sname) := NA_real_]
  vy  <- d[!is.na(get(v)), .(nu = uniqueN(get(v))), by = year][nu > 1, year]
  idx <- which(d$year %in% vy & !is.na(d[[v]]))
  dd  <- d[idx]
  # igofe (Jared 2026-11-23): only the IGO layer is residualized on sender, receiver and year FE + geography (GEO),
  #   residuals below zero set to 0, then scaled 0-1 by the yearly max; all other layers raw
  if (WEIGHT == "igofe" && v == "layer_igo") {
    fit <- fepois(as.formula(paste0(v, " ~ ", GEOVAR, " | ccode1 + ccode2 + year")), data = dd, warn = FALSE, notes = FALSE)
    pred <- predict(fit, newdata = dd, type = "response"); dropped <- is.na(pred)
    if (any(dd[[v]][dropped] != 0)) stop(v, ": dropped observations with nonzero ties"); pred[dropped] <- 0
    r <- pmax(dd[[v]] - pred, 0)
    r <- ave(r, dd$year, FUN = function(z) { m <- max(z); if (m > 0) z / m else z })
    stopifnot(!anyNA(r), all(r >= 0 & r <= 1 + 1e-12))
    set(d, i = idx, j = sname, value = r)
    cat(sprintf("  %-16s igofe: %s b = %.3f, share > 0 after censoring %.3f (raw %.3f)\n", v, GEOVAR, coef(fit)[GEOVAR], mean(r > 0), mean(dd[[v]] > 0)))
    next
  }
  if (WEIGHT %in% c("raw", "igofe")) {
    set(d, i = idx, j = sname, value = dd[[v]])
    cat(sprintf("  %-16s raw, years %d-%d\n", v, min(vy), max(vy)))
    next
  }
  fml <- as.formula(paste0(v, " ~ ", GEOVAR, " | ", FE_TXT))
  fit <- if (v %in% binary) feglm(fml, data = dd, family = binomial("logit"), warn = FALSE, notes = FALSE) else
           if (v == "layer_un") feols(fml, data = dd, warn = FALSE, notes = FALSE) else
                            fepois(fml, data = dd, warn = FALSE, notes = FALSE)
  pred <- if (v == "layer_un") predict(fit, newdata = dd) else predict(fit, newdata = dd, type = "response")
  dropped <- is.na(pred)
  if (WEIGHT == "sy" || v == "layer_un") {
    pred[dropped] <- dd[[v]][dropped]           # perfectly predicted state-years carry no information
  } else {
    if (any(dd[[v]][dropped] != 0)) stop(v, ": dropped observations with nonzero ties")
    pred[dropped] <- 0
  }
  stopifnot(!anyNA(pred), v == "layer_un" || min(pred) >= 0)
  y <- dd[[v]]
  sc <- if (WEIGHT %in% c("sry", "sy")) y - pred else
          if (v %in% binary) ifelse(y == 1, 1, -pred) else ifelse(y >= pred, y, y - pred)
  set(d, i = idx, j = sname, value = sc)
  cat(sprintf("  %-16s years %d-%d  %s b = %.3f  dropped = %d\n", v, min(vy), max(vy), GEOVAR, coef(fit)[GEOVAR], sum(dropped)))
}
# Offensive alliances, DCAs, PTAs: years with no such tie anywhere are true zeros -> score 0
if ("score_un" %in% score_cols) {   # before 1946 and non-member pairs: no information, score 0
  cat(sprintf("  score_un: %d dyad-years without kappa set to 0\n", sum(is.na(d$score_un))))
  set(d, i = which(is.na(d$score_un)), j = "score_un", value = 0)
}
for (sc in intersect(score_cols, c("score_offensive", "score_dca", "score_pta", "score_ptai", "score_ptad", "score_ptadeep", "score_off", "score_neu", "score_nonagg", "score_consul", "score_igo"))) {   # score_igo: fixed MIA layer is empty before 1940 (true zeros)
  lv <- sub("score_", "layer_", sc)
  stopifnot(all(d[is.na(get(sc)), get(lv)] == 0))
  cat(sprintf("  %s: %d dyad-years in years with no tie set to 0\n", sc, sum(is.na(d[[sc]]))))
  set(d, i = which(is.na(d[[sc]])), j = sc, value = 0)
}
for (sc in setdiff(score_cols, "score_trade")) stopifnot(!anyNA(d[[sc]]))
if ("score_trade" %in% score_cols) stopifnot(all(is.na(d$score_trade) == (d$year < 1870)))

d[, n_layers := rowSums(!is.na(.SD)), .SDcols = score_cols]
if (LAYERS == "trade") {          # trade starts in 1870: drop earlier dyad-years (no layer at all)
  stopifnot(all(d[n_layers == 0, year] < 1870), all(d[year >= 1870, n_layers] == 1))
  d <- d[n_layers >= 1]
  cat(sprintf("  trade only: panel restricted to %d-%d\n", min(d$year), max(d$year)))
}
if (LAYERS %in% c("dca", "pta")) {  # single LIO-era layer: keep only years in which the instrument exists
  lv <- layers[1]; yrs <- d[, .(any_tie = any(get(lv) == 1)), by = year][any_tie == TRUE, year]
  d <- d[year %in% yrs]
  cat(sprintf("  %s only: panel restricted to %d-%d\n", LAYERS, min(d$year), max(d$year)))
}
stopifnot(all(d$n_layers >= 1))
d[, tie_dir := rowSums(.SD, na.rm = TRUE) / n_layers, .SDcols = score_cols]
}
rev_d <- d[, .(ccode1 = ccode2, ccode2 = ccode1, year, tie_dir_rev = tie_dir)]
d <- merge(d, rev_d, by = c("ccode1", "ccode2", "year"), all.x = TRUE)
stopifnot(!anyNA(d$tie_dir_rev))
cat(sprintf("  max |tie(i,j) - tie(j,i)| before averaging: %.2e\n", max(abs(d$tie_dir - d$tie_dir_rev))))
d[, tie_signed := (tie_dir + tie_dir_rev) / 2]
und <- d[ccode1 < ccode2]
stopifnot(nrow(und) * 2 == nrow(d))
und[, tie_norm := { m <- max(abs(tie_signed)); if (m > 0) tie_signed / m else 0 }, by = year]
und[, tie_pos := pmax(tie_signed, 0)]
und[, tie_pos := if (max(tie_pos) > 0) tie_pos / max(tie_pos) else 0, by = year]
stopifnot(!anyNA(und[, .(tie_signed, tie_norm, tie_pos)]), min(und$tie_pos) >= 0, max(und$tie_pos) <= 1)
if (WEIGHT %in% c("raw", "w2contig", "igofe") && !(WEIGHT == "raw" && "layer_un" %in% layers)) stopifnot(min(und$tie_signed) >= 0)   # raw UN kappa can be negative

keep <- c("ccode1", "ccode2", "year", layers, "n_layers", "tie_signed", "tie_norm", "tie_pos",
          "cowmaj1", "cowmaj2", "conttype", "capdist", "cinc1", "cinc2", "cowmidonset", "cowmidongoing", "cowmidfatal")
und <- und[, ..keep]
setorder(und, year, ccode1, ccode2)
dir.create(file.path(OUT, "data", "processed"), recursive = TRUE, showWarnings = FALSE)
diag <- und[, .(share_negative = mean(tie_signed < 0), share_positive = mean(tie_signed > 0)), by = year]
cat(sprintf("  dyad-years %d; mean share negative %.2f; mean share positive %.2f\n",
            nrow(und), mean(diag$share_negative), mean(diag$share_positive)))
saveRDS(und, file.path(OUT, "data", "processed", "main_measure_allyigo.rds"))
fwrite(und[, .(ccode1, ccode2, year, tie_signed, tie_norm, tie_pos)],
       file.path(OUT, "data", "processed", "main_measure_edges_allyigo.csv.gz"))
cat("  saved\n")
