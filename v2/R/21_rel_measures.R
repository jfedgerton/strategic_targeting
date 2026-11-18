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
##              W5 fesy  : per layer, residual of v ~ 1 | state x year FE (both ends)          (PI 2026-11-18; FE only)
##              W2-W5: residuals rescaled to [0,1] within layer-year; combo edge = sqrt(sum of squares) rescaled to [0,1] within year
##              (the IO-9933 / 2_downweight_data.R construction: Euclidean root, range01)
##   centrality: pr  = log PageRank (igraph default damping 0.85, undirected, weighted; isolates keep the teleport share)
##               deg = weighted degree / (n_t - 1)  (no log: isolates have 0)
##   IGO arm (--arm=igo): bipartite IGO memberships (igoR, full members), Newman weights 1/(k_g - 1), and cosine-normalised
##               shared-IGO ties; pr and deg of each (4 specs).
##
## Conventions
##   * Layer values: igo and pta enter as log(1 + count); the binary ATOP / DCA layers and the UN kappa enter as they are.
##   * A layer contributes nothing (edge weight 0) in years before it exists (first year with any non-zero value: DCA 1980,
##     UN 1946, PTA first in-force year): STRUCTURAL zeros (PI 2026-11-18). Residual models are fit on existence years only.
##   * UN kappa missing within UN years (non-member / not voting) -> layer weight 0 for that dyad (no tie).
##   * Undirected dyads: the FE models (W4, W5) are fit on the doubled (both orientations) data so the state effects are symmetric;
##     the residual of the original orientation is kept (identical in both orientations by symmetry; asserted).
##   * Vertices of year t = all COW system members in t (create_stateyears); edges with weight 0 are dropped.
##   * KNOWN DEGENERACY (flagged, not fixed): under W5 every state-year's residuals sum to zero, so after the affine within-year
##     rescaling every node of a single-layer network has the SAME strength; deg is then constant within year and PageRank of an
##     undirected graph is close to strength-proportional. Such specs are detected (within-year SD ~ 0) and marked not identified.
## Output: <out_dir>/positions_<arm>_<W>.rds  list(key = data.frame(ccode, year), R = matrix [state-year x spec], specs = data.frame)
##         <out_dir>/specs_<arm>.csv
## Usage: Rscript v2/R/21_rel_measures.R <layers_dy.rds> <out_dir> [--arm=core|dipl|igo] [--cores=N]
## ============================================================================
suppressMessages({ library(dplyr); library(fixest); library(igraph); library(peacesciencer); library(parallel) })
set.seed(123)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2)
LAYER_RDS <- args[1]; OUT_DIR <- args[2]
opt <- function(name, default) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
ARM <- opt("arm", "core"); CORES <- as.integer(opt("cores", "1"))
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
  v <- as.numeric(dy[[l]])
  if (l %in% c("igo", "pta")) v <- log1p(v)
  if (l == "un") v[is.na(v)] <- NA_real_                        # handled below: no tie
  v })
stopifnot(!anyNA(V[, setdiff(LAYERS, "un")]))
first_year <- sapply(LAYERS, function(l) { v <- V[, l]; min(dy$year[!is.na(v) & v != 0]) })
logf("layer existence (first non-zero year): %s", paste(sprintf("%s %d", LAYERS, first_year), collapse = ", "))
dens <- sapply(LAYERS, function(l) mean(V[dy$year >= first_year[l], l] != 0, na.rm = TRUE))
logf("share of non-zero dyad-years within existence years: %s", paste(sprintf("%s %.3f", LAYERS, dens), collapse = ", "))

## ---- 2. per-layer edge weights under each weighting ----------------------------------------------------------------
WEIGHTINGS <- c("W1raw", "W2contig", "W3dist", "W4fe", "W5fesy")
if (nzchar(ONLY_W)) { WEIGHTINGS <- strsplit(ONLY_W, ",")[[1]]; stopifnot(all(WEIGHTINGS %in% c("W1raw", "W2contig", "W3dist", "W4fe", "W5fesy"))) }
dy$lcap <- log1p(dy$capdist)
yr_split <- split(seq_len(nrow(dy)), dy$year)
within_year_range01 <- function(x) { out <- numeric(length(x)); for (ix in yr_split) out[ix] <- range01(x[ix]); out }
## doubled data for symmetric state effects
dbl_idx <- c(seq_len(nrow(dy)), seq_len(nrow(dy)))
dbl <- data.frame(a = c(dy$ccode1, dy$ccode2), b = c(dy$ccode2, dy$ccode1), year = c(dy$year, dy$year),
                  contig_land = c(dy$contig_land, dy$contig_land), lcap = c(dy$lcap, dy$lcap))

edge_w <- function(W, l) {
  v <- V[, l]; ex <- dy$year >= first_year[[l]]
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
  if (W %in% c("W4fe", "W5fesy")) {
    okd <- c(ok, ok); dd <- dbl[okd, ]; dd$v <- c(v, v)[okd]
    m <- if (W == "W4fe") feols(v ~ contig_land + lcap | a + b + year, data = dd, notes = FALSE, warn = FALSE, fixef.rm = "none")
         else feols(v ~ 1 | a^year + b^year, data = dd, notes = FALSE, warn = FALSE, fixef.rm = "none")
    stopifnot(nobs(m) == nrow(dd))                                 # no singleton removal: residuals align with rows
    rr <- resid(m); n1 <- sum(ok)
    stopifnot(max(abs(rr[seq_len(n1)] - rr[n1 + seq_len(n1)])) < 1e-6)   # symmetric by construction
    r <- rr[seq_len(n1)]
  }
  stopifnot(length(r) == sum(ok), !anyNA(r))
  x <- numeric(nrow(dy)); x[ok] <- r
  ## within layer-year rescaling over the dyads that carry the layer; dyads without it (pre-existence, UN NA) stay 0
  for (ix in yr_split) { j <- ix[ok[ix]]; if (length(j)) e[j] <- range01(x[j]) }
  e[!ok] <- 0
  e
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
    EL <- sapply(LAYERS, function(l) edge_w(W, l))          # dyad-year x layer
    stopifnot(!anyNA(EL), all(EL >= 0), all(EL <= 1 + 1e-12))
    logf("%s: per-layer edge weights built (mean weight by layer: %s)", W, paste(sprintf("%s %.3f", LAYERS, colMeans(EL)), collapse = ", "))
    one <- function(k) {
      L <- combos[[k]]
      if (W == "W1raw") { S <- rowSums(EL[, L, drop = FALSE]); E <- numeric(length(S)); for (ix in yr_split) { m <- max(S[ix]); if (m > 0) E[ix] <- S[ix] / m } }
      else { S <- sqrt(rowSums(EL[, L, drop = FALSE]^2)); E <- within_year_range01(S) }
      cent(E)
    }
    res <- if (CORES > 1) mclapply(seq_along(combos), one, mc.cores = CORES, mc.set.seed = FALSE) else lapply(seq_along(combos), one)
    bad <- vapply(res, function(z) inherits(z, "try-error") || is.null(z$pr), TRUE); if (any(bad)) stop(sprintf("%s: %d combos failed", W, sum(bad)))
    Rm <- cbind(sapply(res, `[[`, "pr"), sapply(res, `[[`, "deg"))
    sid <- c(paste(W, combo_id, "pr", sep = "|"), paste(W, combo_id, "deg", sep = "|"))
    colnames(Rm) <- sid
    specs <- data.frame(spec_id = sid, arm = ARM, weighting = W, layers = rep(combo_id, 2), n_layers = rep(lengths(combos), 2),
                        centrality = rep(c("pr", "deg"), each = length(combos)), stringsAsFactors = FALSE)
    specs$const_year_share <- apply(Rm, 2, ident_check)
    specs$const_years <- apply(Rm, 2, const_years)
    specs$identified <- specs$const_year_share < 0.05
    stopifnot(!anyNA(Rm), all(is.finite(Rm)))
    saveRDS(list(key = sy, R = Rm, specs = specs, built = format(Sys.time()), layer_rds = normalizePath(LAYER_RDS)), file.path(OUT_DIR, sprintf("positions_%s_%s.rds", ARM, W)))
    spec_rows[[W]] <- specs
    logf("%s: %d specs in %.1f min; not identified (constant in >= 5%% of years): %d", W, ncol(Rm), as.numeric(difftime(Sys.time(), t0, units = "mins")), sum(!specs$identified))
  }
  write.csv(bind_rows(spec_rows), file.path(OUT_DIR, sprintf("specs_%s.csv", ARM)), row.names = FALSE)
}

## ---- 5. IGO arm: bipartite (Newman) and cosine-normalised shared-IGO ties --------------------------------------------
if (ARM == "igo") {
  d <- igoR::state_year_format3
  igo_cols <- setdiff(names(d), c("ccode", "year", "state"))
  Enew <- numeric(nrow(dy)); Ecos <- numeric(nrow(dy))
  for (y in names(yr_split)) {
    dd <- d[d$year == as.integer(y), ]
    B <- (as.matrix(dd[, igo_cols]) == 1) * 1; rownames(B) <- dd$ccode
    B <- B[, colSums(B) >= 2, drop = FALSE]
    ei <- yr_split[[y]]
    if (ncol(B) == 0) next
    kg <- colSums(B)
    N <- B %*% (t(B) / (kg - 1)); diag(N) <- 0                  # Newman: each shared IGO weighted 1 / (k_g - 1)
    A <- B %*% t(B); ks <- sqrt(diag(A)); C <- A / outer(pmax(ks, 1e-12), pmax(ks, 1e-12)); diag(C) <- 0   # cosine of membership vectors
    i1 <- match(as.character(dy$ccode1[ei]), rownames(B)); i2 <- match(as.character(dy$ccode2[ei]), rownames(B))
    okp <- !is.na(i1) & !is.na(i2)
    Enew[ei[okp]] <- N[cbind(i1[okp], i2[okp])]; Ecos[ei[okp]] <- C[cbind(i1[okp], i2[okp])]
  }
  logf("IGO arm: dyad-years with a positive Newman tie %.3f, cosine tie %.3f", mean(Enew > 0), mean(Ecos > 0))
  cN <- cent(Enew); cC <- cent(Ecos)
  Rm <- cbind(cN$pr, cN$deg, cC$pr, cC$deg)
  sid <- c("IGOarm|bip_newman|pr", "IGOarm|bip_newman|deg", "IGOarm|cosine|pr", "IGOarm|cosine|deg"); colnames(Rm) <- sid
  specs <- data.frame(spec_id = sid, arm = "igo", weighting = c("newman", "newman", "cosine", "cosine"), layers = "igo", n_layers = 1L,
                      centrality = c("pr", "deg", "pr", "deg"), stringsAsFactors = FALSE)
  ## note (review P3-16): Newman-weighted degree equals each state's IGO membership count (sum over its IGOs of k_g-1 partners x 1/(k_g-1))
  specs$const_year_share <- apply(Rm, 2, ident_check); specs$const_years <- apply(Rm, 2, const_years); specs$identified <- specs$const_year_share < 0.05
  saveRDS(list(key = sy, R = Rm, specs = specs, built = format(Sys.time()), layer_rds = normalizePath(LAYER_RDS)), file.path(OUT_DIR, "positions_igo_arm.rds"))
  write.csv(specs, file.path(OUT_DIR, "specs_igo.csv"), row.names = FALSE)
  logf("IGO arm: 4 specs written")
}
logf("done")
