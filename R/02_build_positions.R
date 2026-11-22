## ============================================================================
## 21_rel_measures.R  --  Part 2b of the relational-measure sweep (PI design 2026-11-18)
## From the dyad-year layer panel (20_build_layers.R) build every candidate relational POSITION (state-year), so that
## 22_rel_sweep.R can swap it into the target-choice model in place of log IGO PageRank.
##
## Grid (fixed in advance; every spec is reported):
##   layers   : 9 core layers  igo, def, off, neu, nonagg, consul, pta, dca, un
##              combos = every non-empty subset that contains igo OR def (2^9 - 1 - (2^7 - 1) = 384)
##              --arm=dipl adds diplomatic exchange to every core combo (384 more; separate arm, 5-year waves carried forward)
##   weighting: W1 raw   : each layer rescaled to [0,1] within year; combo edge = layer sum / within-year max of the sum
##              W2 contig: per layer, residual of v ~ contig_land + year (pooled OLS over the layer's existence years)
##              W3 dist  : per layer, residual of v ~ log(1 + capdist) + year
##              W4 fe    : per layer, residual of v ~ contig_land + log(1 + capdist) | state FE (both ends) + year FE
##              (W5, state x year FE, dropped by the PI 2026-11-18: it removes each state's tie level, so centrality is unidentified)
##              W2-W4: residuals clipped at 0 and divided by the layer-year maximum; combo edge = sqrt(sum of squares) rescaled to
##              [0,1] within year (IO-9933 Euclidean root of clipped positive residuals; PI 2026-11-18)
##   PTAs: --pta=all (signed or in force) | inforce (only treaties with an entry-into-force year); both are run (PI 2026-11-18)
##         | any (PI 2026-11-22): unweighted, 1 if the dyad has any PTA (signed or in force), not logged; like inforce, only combos with PTA
##   Book-aligned weightings (PI 2026-11-22; --weightings=W6srybook,W7sryyear, core arm only), the book project's 'sry' tie model:
##              per layer an expected-tie model on the doubled data -- logit for binary layers (ATOP types, DCA, PTA any), Poisson for counts
##              (IGO, PTA all/inforce; raw counts), OLS for UN kappa -- with W6srybook: log(1+capdist) | sender + receiver + year FE (the book),
##              W7sryyear: log(1+capdist) + contig_land | year FE (no state FE). Tie score = observed - expected (symmetric), divided by the
##              layer-year max |score|; combo score = mean over the combo's layers that apply to the dyad-year (signed); edge = positive part,
##              rescaled to [0,1] within year. Validated against a standalone build on 2026-11-21 (docs/measure_alignment_2026-11-21.md).
##   centrality: pr  = log PageRank (igraph default damping 0.85, undirected, weighted; isolates keep the teleport share)
##               deg = weighted degree / (n_t - 1)  (no log: isolates have 0)
##   IGO arm (--arm=igo): bipartite IGO memberships (igoR, full members), Newman weights 1/(k_g - 1), and cosine-normalised
##               shared-IGO ties, plus the raw shared-IGO count network of 01_build_data.R as an anchor; pr and deg of each (6 specs).
##
## Conventions
##   * Layer values: igo and pta enter as log(1 + count); the binary ATOP / DCA layers and the UN kappa enter as they are.
##   * A layer contributes nothing (edge weight 0) in years before it exists (first year with any non-zero value: DCA 1980,
##     UN 1946, PTA first in-force year): STRUCTURAL zeros (PI 2026-11-18). Residual models are fit on existence years only.
##   * UN kappa missing within UN years (non-member / not voting) -> layer weight 0 for that dyad (no tie).
##   * Undirected dyads: the FE models (W4, W5) are fit on the doubled (both orientations) data so the state effects are symmetric;
##     the residual of the original orientation is kept (identical in both orientations by symmetry; asserted).
##   * Vertices of year t = all COW system members in t (create_stateyears); edges with weight 0 are dropped.
##   * A spec whose position is constant within a year is marked (const_years); 22/23 refuse it if an event uses that year.
## Output: <out_dir>/positions_<arm>_<W>.rds  list(key = data.frame(ccode, year), R = matrix [state-year x spec], specs = data.frame)
##         <out_dir>/specs_<arm>.csv
## Usage: Rscript v2/R/21_rel_measures.R <layers_dy.rds> <out_dir> [--arm=core|dipl|igo] [--pta=all|inforce] [--cores=N]
## ============================================================================
suppressMessages({ library(dplyr); library(fixest); library(igraph); library(peacesciencer); library(parallel) })
set.seed(123)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2)
LAYER_RDS <- args[1]; OUT_DIR <- args[2]
opt <- function(name, default) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
ARM <- opt("arm", "core"); CORES <- as.integer(opt("cores", "1"))
PTA <- opt("pta", "all"); stopifnot(PTA %in% c("all", "inforce", "any"))   # PI 2026-11-18: PTAs both ways; PI 2026-11-22: + unweighted (any)
TAG <- paste0("pta", PTA)
ONLY_W <- opt("weightings", ""); MAX_COMBOS <- as.integer(opt("max_combos", "0"))   # smoke-test options; a full run uses neither
stopifnot(ARM %in% c("core", "dipl", "igo"), CORES >= 1)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
logf <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
logf("R %s | fixest %s | igraph %s | peacesciencer %s | arm %s | cores %d", getRversion(), packageVersion("fixest"), packageVersion("igraph"),
     packageVersion("peacesciencer"), ARM, CORES)

range01 <- function(x) { r <- range(x); if (r[2] - r[1] < 1e-12) rep(0, length(x)) else (x - r[1]) / (r[2] - r[1]) }

## ---- 0. state-years (vertices) and dyad-years (edges) ---------------------------------------------------------
sy <- create_stateyears(system = "cow", subset_years = 1816:2014) |> transmute(ccode = as.integer(ccode), year = as.integer(year)) |>
  arrange(year, ccode)
stopifnot(!anyDuplicated(sy))
dy <- readRDS(LAYER_RDS) |> arrange(year, ccode1, ccode2)
stopifnot(!anyDuplicated(dy[, c("ccode1", "ccode2", "year")]), all(dy$ccode1 < dy$ccode2))
## every dyad endpoint must be a vertex of its year
k_sy <- paste(sy$ccode, sy$year); stopifnot(all(paste(dy$ccode1, dy$year) %in% k_sy), all(paste(dy$ccode2, dy$year) %in% k_sy))
years <- sort(unique(sy$year))
logf("vertices %d state-years (%d-%d); edges %d dyad-years", nrow(sy), min(years), max(years), nrow(dy))

## ---- 1. layer values and existence years --------------------------------------------------------------------------
CORE <- c("igo", "def", "off", "neu", "nonagg", "consul", "pta", "dca", "un")
LAYERS <- if (ARM == "dipl") c(CORE, "dipl") else CORE
V <- sapply(LAYERS, function(l) {
  v <- as.numeric(dy[[if (l == "pta" && PTA == "inforce") "pta_inforce" else l]])
  if (l == "pta" && PTA == "any") v <- as.numeric(v > 0)          # PI 2026-11-22: unweighted PTA (not logged)
  else if (l %in% c("igo", "pta")) v <- log1p(v)
  if (l == "un") v[is.na(v)] <- NA_real_                        # handled below: no tie
  v })
stopifnot(!anyNA(V[, setdiff(LAYERS, "un")]))
first_year <- sapply(LAYERS, function(l) { v <- V[, l]; min(dy$year[!is.na(v) & v != 0]) })
logf("layer existence (first non-zero year): %s", paste(sprintf("%s %d", LAYERS, first_year), collapse = ", "))
dens <- sapply(LAYERS, function(l) mean(V[dy$year >= first_year[l], l] != 0, na.rm = TRUE))
logf("share of non-zero dyad-years within existence years: %s", paste(sprintf("%s %.3f", LAYERS, dens), collapse = ", "))

## ---- 2. per-layer edge weights under each weighting ----------------------------------------------------------------
WEIGHTINGS <- c("W1raw", "W2contig", "W3dist", "W4fe")      # PI 2026-11-18: W5 (state x year FE) dropped; default set unchanged
SRY <- c("W6srybook", "W7sryyear")                           # PI 2026-11-22: book-aligned weightings, requested explicitly with --weightings
if (nzchar(ONLY_W)) { WEIGHTINGS <- strsplit(ONLY_W, ",")[[1]]; stopifnot(all(WEIGHTINGS %in% c("W1raw", "W2contig", "W3dist", "W4fe", SRY))) }
if (any(WEIGHTINGS %in% SRY)) stopifnot(ARM == "core")        # no tie model is specified for the diplomacy layer
## layer-years with no non-zero value (e.g. offensive alliances in years nobody holds one) carry ZERO weight, like pre-existence
## years (PI 2026-11-18, review P1-2): otherwise the residual of an all-zero year is a pure contiguity / distance surface
yr_has <- sapply(LAYERS, function(l) { v <- V[, l]; a <- tapply(!is.na(v) & v != 0, dy$year, any); a[as.character(dy$year)] })
logf("empty layer-years after existence (zero weight): %s", paste(sprintf("%s %d", LAYERS, sapply(LAYERS, function(l)
  sum(tapply(yr_has[, l], dy$year, `[`, 1) == FALSE & as.integer(names(tapply(yr_has[, l], dy$year, `[`, 1))) >= first_year[[l]]))), collapse = ", "))
dy$lcap <- log1p(dy$capdist)
yr_split <- split(seq_len(nrow(dy)), dy$year)
within_year_range01 <- function(x) { out <- numeric(length(x)); for (ix in yr_split) out[ix] <- range01(x[ix]); out }
## doubled data for symmetric state effects
dbl_idx <- c(seq_len(nrow(dy)), seq_len(nrow(dy)))
dbl <- data.frame(a = c(dy$ccode1, dy$ccode2), b = c(dy$ccode2, dy$ccode1), year = c(dy$year, dy$year),
                  contig_land = c(dy$contig_land, dy$contig_land), lcap = c(dy$lcap, dy$lcap))

edge_w <- function(W, l) {
  v <- V[, l]; ex <- dy$year >= first_year[[l]]
  ex <- ex & yr_has[, l]
  if (l == "un") { ok <- ex & !is.na(v) } else ok <- ex
  e <- numeric(nrow(dy))
  if (W == "W1raw") {
    x <- ifelse(ok, v, 0); x[!ex] <- 0
    e <- within_year_range01(x); e[!ok] <- 0
    return(e)
  }
  d <- data.frame(v = v, year = dy$year, contig_land = dy$contig_land, lcap = dy$lcap)[ok, ]
  if (W == "W2contig") r <- resid(lm(v ~ contig_land + year, data = d))
  if (W == "W3dist")   r <- resid(lm(v ~ lcap + year, data = d))
  if (W == "W4fe") {
    okd <- c(ok, ok); dd <- dbl[okd, ]; dd$v <- c(v, v)[okd]
    m <- feols(v ~ contig_land + lcap | a + b + year, data = dd, notes = FALSE, warn = FALSE, fixef.rm = "none")
    stopifnot(nobs(m) == nrow(dd))                                 # no singleton removal: residuals align with rows
    rr <- resid(m); n1 <- sum(ok)
    stopifnot(max(abs(rr[seq_len(n1)] - rr[n1 + seq_len(n1)])) < 1e-6)   # symmetric by construction
    r <- rr[seq_len(n1)]
  }
  stopifnot(length(r) == sum(ok), !anyNA(r))
  x <- numeric(nrow(dy)); x[ok] <- r
  ## PI 2026-11-18: residuals CLIPPED at zero (only ties stronger than expected count; IO-9933 Euclidean version), then divided by the
  ## layer-year maximum so each layer lies in [0, 1] before combining; a dyad at or below expectation has weight 0
  ## GPT review P1-2: residuals within numerical tolerance of zero are zero, so solver noise in a layer the covariates fit exactly
  ## cannot be stretched into weight-1 ties by the year-maximum division; tolerance relative to the layer's own scale
  tol <- 1e-8 * max(1, max(abs(v[ok])))
  x[abs(x) < tol] <- 0
  x <- pmax(x, 0)
  n_num_empty <- 0L
  for (ix in yr_split) { j <- ix[ok[ix]]; if (length(j)) { m <- max(x[j]); if (m > 0) e[j] <- x[j] / m else n_num_empty <- n_num_empty + 1L } }
  e[!ok] <- 0
  attr(e, "n_num_empty") <- n_num_empty
  e
}

## PI 2026-11-22: book-aligned "sry" tie score per layer (signed, NA where the layer does not apply); see the header
SRY_FAM <- c(igo = "pois", def = "logit", off = "logit", neu = "logit", nonagg = "logit", consul = "logit", dca = "logit", un = "ols",
             pta = if (PTA == "any") "logit" else "pois")
sry_score <- function(W, l) {
  v <- V[, l]; ok <- dy$year >= first_year[[l]] & yr_has[, l] & !is.na(v)
  if (SRY_FAM[[l]] == "pois") v <- round(expm1(v))                # Poisson on the raw counts (V holds log1p for igo / pta counts)
  stopifnot(all(abs(v[ok] - round(v[ok])) < 1e-6) || SRY_FAM[[l]] != "pois", all(v[ok] %in% c(0, 1)) || SRY_FAM[[l]] != "logit")
  n1 <- sum(ok); okd <- c(ok, ok); dd <- dbl[okd, ]; dd$v <- c(v, v)[okd]
  fml <- if (W == "W6srybook") v ~ lcap | a + b + year else v ~ lcap + contig_land | year
  m <- switch(SRY_FAM[[l]], ols = feols(fml, data = dd, notes = FALSE, warn = FALSE),
              logit = feglm(fml, data = dd, family = binomial("logit"), notes = FALSE, warn = FALSE),
              pois = fepois(fml, data = dd, notes = FALSE, warn = FALSE))
  ex <- rep(NA_real_, nrow(dd)); ex[obs(m)] <- fitted(m)
  ndrop <- sum(is.na(ex)); ex[is.na(ex)] <- dd$v[is.na(ex)]       # rows dropped as perfectly predicted (all-0 or all-1 FE groups): residual 0
  r2 <- dd$v - ex; asym <- max(abs(r2[seq_len(n1)] - r2[n1 + seq_len(n1)])); stopifnot(asym < 1e-4)
  x <- rep(NA_real_, nrow(dy)); x[ok] <- (r2[seq_len(n1)] + r2[n1 + seq_len(n1)]) / 2
  for (ix in yr_split) { j <- ix[ok[ix]]; if (length(j)) { mx <- max(abs(x[j])); x[j] <- if (mx > 0) x[j] / mx else 0 } }
  attr(x, "info") <- sprintf("%s %s: n %d, perfectly predicted %d, share > 0 %.3f", l, SRY_FAM[[l]], n1, ndrop, mean(x[ok] > 0))
  x
}

## ---- 3. combos -------------------------------------------------------------------------------------------------------
if (ARM %in% c("core", "dipl")) {
  combos <- list()
  for (mask in 1:(2^length(CORE) - 1)) {
    inc <- CORE[bitwAnd(mask, 2^(seq_along(CORE) - 1)) > 0]
    if (!("igo" %in% inc || "def" %in% inc)) next
    combos[[length(combos) + 1]] <- if (ARM == "dipl") c(inc, "dipl") else inc
  }
  stopifnot(length(combos) == 384L)
  ## GPT review P2-8: under --pta=inforce, combos without the PTA layer are identical to the --pta=all ones; build only combos with PTA
  if (PTA %in% c("inforce", "any")) { combos <- combos[vapply(combos, function(x) "pta" %in% x, TRUE)]; stopifnot(length(combos) == 192L)
    logf("pta = %s: only the %d combos containing PTA are built (the rest duplicate --pta=all)", PTA, length(combos)) }
  if (MAX_COMBOS > 0) { combos <- combos[c(1:min(MAX_COMBOS, 384))]; logf("SMOKE TEST: first %d combos only", length(combos)) }
  combo_id <- vapply(combos, function(x) paste(x, collapse = "+"), "")
  logf("%d layer combos (each contains igo or def)%s", length(combos), if (ARM == "dipl") " + dipl" else "")
}

## ---- 4. centrality per year ---------------------------------------------------------------------------------------
sy_year_idx <- split(seq_len(nrow(sy)), sy$year)
## vertex index of each dyad endpoint within its year
pos_in_year <- integer(nrow(sy)); for (ix in sy_year_idx) pos_in_year[ix] <- seq_along(ix)
key_sy <- setNames(seq_len(nrow(sy)), k_sy)
u1 <- key_sy[paste(dy$ccode1, dy$year)]; u2 <- key_sy[paste(dy$ccode2, dy$year)]
stopifnot(!anyNA(u1), !anyNA(u2))
cent <- function(E) {                                           # E: edge weight per dyad-year row of dy (>= 0)
  pr <- numeric(nrow(sy)); dg <- numeric(nrow(sy))
  for (y in names(sy_year_idx)) {
    vi <- sy_year_idx[[y]]; n <- length(vi); ei <- yr_split[[y]]
    w <- E[ei]; keep <- w > 0
    a <- pos_in_year[u1[ei][keep]]; b <- pos_in_year[u2[ei][keep]]; w <- w[keep]
    s <- numeric(n)
    if (length(w)) { agg <- rowsum(c(w, w), c(a, b)); s[as.integer(rownames(agg))] <- agg[, 1] }
    g <- make_empty_graph(n = n, directed = FALSE)
    if (length(w)) g <- add_edges(g, rbind(a, b), attr = list(weight = w))
    p <- page_rank(g, directed = FALSE, weights = if (length(w)) E(g)$weight else NULL)$vector
    stopifnot(length(p) == n, all(p > 0), abs(sum(p) - 1) < 1e-8)
    pr[vi] <- log(p); dg[vi] <- s / (n - 1)
  }
  list(pr = pr, deg = dg)
}

ident_check <- function(R) {                                   # share of years in which the position is (near) constant
  sds <- vapply(sy_year_idx, function(ix) sd(R[ix]), 0)
  mean(sds < 1e-8)
}
## review fix (P1-1): the years in which the position is constant are recorded per spec; 22/23 refuse a spec whose constant years
## are used by any event (as t for the year scale or t-1 for positions), whatever the share
const_years <- function(R) { sds <- vapply(sy_year_idx, function(ix) sd(R[ix]), 0); paste(names(sds)[sds < 1e-8], collapse = " ") }

if (ARM %in% c("core", "dipl")) {
  spec_rows <- list()
  for (W in WEIGHTINGS) {
    t0 <- Sys.time()
    if (W %in% SRY) {
      EW <- lapply(LAYERS, function(l) sry_score(W, l)); names(EW) <- LAYERS
      for (l in LAYERS) logf("%s: %s", W, attr(EW[[l]], "info"))
      EL <- sapply(EW, as.numeric); colnames(EL) <- LAYERS; rm(EW)
      stopifnot(all(is.na(EL) | (EL >= -1 - 1e-12 & EL <= 1 + 1e-12)))
    } else {
    EW <- lapply(LAYERS, function(l) edge_w(W, l)); names(EW) <- LAYERS
    if (W != "W1raw") logf("%s: layer-years with no positive residual (zero weight) %s", W,
                           paste(sprintf("%s %d", LAYERS, vapply(EW, function(e) attr(e, "n_num_empty"), 0L)), collapse = ", "))
    EL <- sapply(EW, as.numeric); colnames(EL) <- LAYERS; rm(EW)   # dyad-year x layer
    stopifnot(!anyNA(EL), all(EL >= 0), all(EL <= 1 + 1e-12))
    logf("%s: per-layer edge weights built (mean weight by layer: %s)", W, paste(sprintf("%s %.3f", LAYERS, colMeans(EL)), collapse = ", "))
    }
    one <- function(k) {
      L <- combos[[k]]
      if (W %in% SRY) { M <- EL[, L, drop = FALSE]; k <- rowSums(!is.na(M)); sc <- ifelse(k > 0, rowSums(M, na.rm = TRUE) / pmax(k, 1), 0)
        S <- pmax(sc, 0); E <- numeric(length(S)); for (ix in yr_split) { m <- max(S[ix]); if (m > 0) E[ix] <- S[ix] / m } }
      else if (W == "W1raw") { S <- rowSums(EL[, L, drop = FALSE]); E <- numeric(length(S)); for (ix in yr_split) { m <- max(S[ix]); if (m > 0) E[ix] <- S[ix] / m } }
      else { S <- sqrt(rowSums(EL[, L, drop = FALSE]^2)); E <- within_year_range01(S) }
      cent(E)
    }
    res <- if (CORES > 1) mclapply(seq_along(combos), one, mc.cores = CORES, mc.set.seed = FALSE) else lapply(seq_along(combos), one)
    bad <- vapply(res, function(z) inherits(z, "try-error") || is.null(z$pr), TRUE); if (any(bad)) stop(sprintf("%s: %d combos failed", W, sum(bad)))
    Rm <- cbind(sapply(res, `[[`, "pr"), sapply(res, `[[`, "deg"))
    sid <- c(paste(paste0(W, "_", TAG), combo_id, "pr", sep = "|"), paste(paste0(W, "_", TAG), combo_id, "deg", sep = "|"))
    colnames(Rm) <- sid
    specs <- data.frame(spec_id = sid, arm = ARM, weighting = W, pta = PTA, layers = rep(combo_id, 2), n_layers = rep(lengths(combos), 2),
                        centrality = rep(c("pr", "deg"), each = length(combos)), stringsAsFactors = FALSE)
    specs$const_year_share <- apply(Rm, 2, ident_check)
    specs$const_years <- apply(Rm, 2, const_years)
    specs$identified <- specs$const_year_share < 0.05
    stopifnot(!anyNA(Rm), all(is.finite(Rm)))
    saveRDS(list(key = sy, R = Rm, specs = specs, built = format(Sys.time()), layer_rds = normalizePath(LAYER_RDS)), file.path(OUT_DIR, sprintf("positions_%s_%s_%s.rds", ARM, W, TAG)))
    spec_rows[[W]] <- specs
    logf("%s: %d specs in %.1f min; not identified (constant in >= 5%% of years): %d", W, ncol(Rm), as.numeric(difftime(Sys.time(), t0, units = "mins")), sum(!specs$identified))
  }
  write.csv(bind_rows(spec_rows), file.path(OUT_DIR, sprintf("specs_%s_%s%s.csv", ARM, TAG, if (all(WEIGHTINGS %in% SRY)) "_sry" else "")), row.names = FALSE)
}

## ---- 5. IGO arm: bipartite (Newman) and cosine-normalised shared-IGO ties --------------------------------------------
if (ARM == "igo") {
  d <- igoR::state_year_format3
  igo_cols <- setdiff(names(d), c("ccode", "year", "state"))
  Enew <- numeric(nrow(dy)); Ecos <- numeric(nrow(dy))
  for (y in names(yr_split)) {
    dd <- d[d$year == as.integer(y), ]
    B <- (as.matrix(dd[, igo_cols]) == 1) * 1; rownames(B) <- dd$ccode
    Bfull <- B                                                   # GPT review P2-4: cosine norms use the FULL membership vectors
    B <- B[, colSums(B) >= 2, drop = FALSE]                      # Newman needs k_g >= 2
    ei <- yr_split[[y]]
    if (ncol(B) == 0) next
    kg <- colSums(B)
    N <- B %*% (t(B) / (kg - 1)); diag(N) <- 0                  # Newman: each shared IGO weighted 1 / (k_g - 1)
    A <- Bfull %*% t(Bfull); ks <- sqrt(diag(A)); C <- A / outer(pmax(ks, 1e-12), pmax(ks, 1e-12)); diag(C) <- 0   # cosine of full membership vectors
    i1 <- match(as.character(dy$ccode1[ei]), rownames(B)); i2 <- match(as.character(dy$ccode2[ei]), rownames(B))
    okp <- !is.na(i1) & !is.na(i2)
    Enew[ei[okp]] <- N[cbind(i1[okp], i2[okp])]; Ecos[ei[okp]] <- C[cbind(i1[okp], i2[okp])]
  }
  logf("IGO arm: dyad-years with a positive Newman tie %.3f, cosine tie %.3f", mean(Enew > 0), mean(Ecos > 0))
  cN <- cent(Enew); cC <- cent(Ecos)
  cR <- cent(as.numeric(dy$igo))                                  # anchor: raw shared-IGO counts, as 01_build_data.R (the current measure)
  Rm <- cbind(cN$pr, cN$deg, cC$pr, cC$deg, cR$pr, cR$deg)
  sid <- c("IGOarm|bip_newman|pr", "IGOarm|bip_newman|deg", "IGOarm|cosine|pr", "IGOarm|cosine|deg", "IGOarm|raw_count|pr", "IGOarm|raw_count|deg"); colnames(Rm) <- sid
  specs <- data.frame(spec_id = sid, arm = "igo", weighting = c("newman", "newman", "cosine", "cosine", "raw_count", "raw_count"), pta = NA, layers = "igo", n_layers = 1L,
                      centrality = rep(c("pr", "deg"), 3), stringsAsFactors = FALSE)
  ## note (review P3-16): Newman-weighted degree equals each state's IGO membership count (sum over its IGOs of k_g-1 partners x 1/(k_g-1))
  specs$const_year_share <- apply(Rm, 2, ident_check); specs$const_years <- apply(Rm, 2, const_years); specs$identified <- specs$const_year_share < 0.05
  saveRDS(list(key = sy, R = Rm, specs = specs, built = format(Sys.time()), layer_rds = normalizePath(LAYER_RDS)), file.path(OUT_DIR, "positions_igo_arm.rds"))
  write.csv(specs, file.path(OUT_DIR, "specs_igo.csv"), row.names = FALSE)
  logf("IGO arm: 6 specs written (raw_count|pr is the current measure rebuilt here; 23 reports its correlation with the build)")
}
logf("done")
