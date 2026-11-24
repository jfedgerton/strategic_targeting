###############################################################################
# 01_import_layers.R
# Import the dyad-year layer panel built for Tremors and keep only what the
# book's main measure needs.
#
# Input : data/raw/tremors_panel/dyad_year_network_layers.rds
#         (built by data/raw/tremors_panel/SOURCE_01_build_network_layers.R;
#          peacesciencer base + ATOP 5.1, COW IGO, COW trade, Kinne DCA v1.0)
# Output: data/processed/layers_dyad_year.rds
#
# Layer definitions inherited from the source script
#   layer_defensive  ATOP defense pledge (binary)
#   layer_offensive  ATOP offense pledge (binary)
#   layer_dca        Kinne DCA, dcaAnyV2 (binary); 0 outside 1980-2010
#   layer_igo        shared COW IGO memberships, min-max scaled within year
#   layer_trade      smoothflow1 + smoothflow2, min-max scaled within year
# The panel is directed: every undirected dyad appears as (i,j) and (j,i).
###############################################################################

suppressMessages(library(data.table))
set.seed(123)

BOOK <- Sys.getenv("BOOK_ROOT", getwd())
in_file  <- file.path(BOOK, "data", "raw", "tremors_panel", "dyad_year_network_layers.rds")
out_file <- file.path(BOOK, "data", "processed", "layers_dyad_year.rds")

d <- readRDS(in_file)
setDT(d)
cat(sprintf("Read %s rows, %d columns\n", format(nrow(d), big.mark = ","), ncol(d)))

layers <- c("layer_defensive", "layer_offensive", "layer_dca", "layer_igo", "layer_trade")

# Logged alternatives for the continuous layers (multiverse; Jared 2026-10-11):
#   layer_igo_log   = within-year 0-1 scaling of log(1 + shared IGOs)
#   layer_trade_log = within-year 0-1 scaling of log(1 + smoothflow1 + smoothflow2)
range01 <- function(x) { r <- max(x) - min(x); if (r == 0) rep(0, length(x)) else (x - min(x)) / r }
d[, igo_raw := fifelse(is.na(dyadigos), 0, dyadigos)]
d[, trade_raw := fifelse(is.na(smoothflow1 + smoothflow2) | (smoothflow1 + smoothflow2) < 0, 0, smoothflow1 + smoothflow2)]
d[, layer_igo_log := range01(log1p(igo_raw)), by = year]
d[, layer_trade_log := range01(log1p(trade_raw)), by = year]

d <- d[, c("ccode1", "ccode2", "year", layers, "layer_igo_log", "layer_trade_log",
           "cowmaj1", "cowmaj2", "conttype", "capdist", "cinc1", "cinc2",
           "cowmidonset", "cowmidongoing", "cowmidfatal"), with = FALSE]

# ---- Assertions ------------------------------------------------------------
stopifnot(!anyDuplicated(d[, .(ccode1, ccode2, year)]))
stopifnot(all(d$ccode1 != d$ccode2))
stopifnot(range(d$year)[1] == 1816, range(d$year)[2] == 2014)
for (v in layers) {
  if (anyNA(d[[v]])) stop("NA in ", v)
  if (min(d[[v]]) < 0 || max(d[[v]]) > 1) stop(v, " outside [0,1]")
}
# Directed panel: ATOP pledges can be one-sided (e.g., FRN->BEL 1832,
# RUS->TUR 1833). The book treats ties as undirected, so a dyad is allied if
# either side pledges: symmetrize alliance layers with the max of (i,j),(j,i).
# DCA, IGO, and trade must already be symmetric.
rev_d <- d[, c("ccode2", "ccode1", "year", layers), with = FALSE]
setnames(rev_d, c("ccode1", "ccode2", "year", paste0(layers, "_rev")))
d <- merge(d, rev_d, by = c("ccode1", "ccode2", "year"), all.x = TRUE)
stopifnot(!anyNA(d[, paste0(layers, "_rev"), with = FALSE]))
for (v in c("layer_defensive", "layer_offensive")) {
  n_asym <- sum(d[[v]] != d[[paste0(v, "_rev")]])
  cat(sprintf("  %s: %d directed dyad-years one-sided; symmetrized with max\n", v, n_asym))
  set(d, j = v, value = pmax(d[[v]], d[[paste0(v, "_rev")]]))
}
for (v in layers) {
  gap <- max(abs(d[[v]] - pmax(d[[v]], d[[paste0(v, "_rev")]])))
  gap2 <- max(abs(d[[v]] - d[[paste0(v, "_rev")]]))
  if (!(v %in% c("layer_defensive", "layer_offensive")) && gap2 > 1e-9)
    stop(v, " is not symmetric across (i,j)/(j,i); max gap ", gap2)
}
d[, (paste0(layers, "_rev")) := NULL]
rm(rev_d)

# ---- Coverage table (years in which each layer varies) ----------------------
cov <- rbindlist(lapply(layers, function(v) {
  vy <- d[, .(nu = uniqueN(get(v))), by = year][nu > 1, year]
  data.table(layer = v, first_year = min(vy), last_year = max(vy),
             n_years_varying = length(vy),
             share_dyadyears_positive = mean(d[[v]] > 0))
}))
print(cov)
fwrite(cov, file.path(BOOK, "results", "main", "layer_coverage.csv"))

saveRDS(d, out_file, compress = "xz")
cat("Saved", out_file, "\n")
