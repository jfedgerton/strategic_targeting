## R/lib/grid.R -- the network-measure grid (PI 2026-11-24): spec definitions and edge weights, from grid_layers.rds (R/grid_layers.R).
## One function makes the edges of a spec, so the positions (R/grid_build.R) and the stage-2 ties (s2_spec_edges) use identical edges.
## Spec id: G<block>_<igo>_<pta>|<layers joined by +>|pr
##   block raw    : raw layers (main sweep)             igo in count / cos / hyper (na if no igo); pta = book (na if no pta). PTA depth dropped from the sweep (PI 2026-11-24)
##         dca1   : raw, DCA = dcaAnyV1 (combos with dca only)
##         ally   : raw, + the other-alliances layer (neutrality / non-aggression / consultation as one unit)
##         sryc12 / srydist / syc12 / sydist : weighted layer scores (igo = count, pta = book only)
##   layers: subset of igo, def, off, dca, pta (+ ally in block ally) containing igo or def
## Edge weight: raw blocks = mean of the selected layers / its within-year maximum; weighted blocks = mean of the symmetrized layer scores,
## positive part / its within-year maximum (book tie_pos). Only w > 0 rows are returned (a, b = ccode1 < ccode2).
GRID_BASE <- c("igo", "def", "off", "dca", "pta")
grid_combos <- function() { out <- list()
  for (k in 1:5) for (cb in combn(GRID_BASE, k, simplify = FALSE)) if (any(c("igo", "def") %in% cb)) out[[length(out) + 1]] <- cb
  stopifnot(length(out) == 24L); out }
grid_specs <- function() {
  S <- list(); add <- function(block, L, igo, pta) S[[length(S) + 1]] <<- data.table(block = block, layers = paste(L, collapse = "+"), igo = igo, pta = pta,
    spec_id = sprintf("G%s_%s_%s|%s|pr", block, igo, pta, paste(L, collapse = "+")))
  for (L in grid_combos()) {
    IG <- if ("igo" %in% L) c("count", "cos", "hyper") else "na"; PT <- if ("pta" %in% L) "book" else "na"
    for (i in IG) for (p in PT) { add("raw", L, i, p); if ("dca" %in% L) add("dca1", L, i, p); add("ally", c(L, "ally"), i, p) }
    for (b in c("sryc12", "srydist", "syc12", "sydist")) add(b, L, if ("igo" %in% L) "count" else "na", if ("pta" %in% L) "book" else "na")
  }
  S <- rbindlist(S); stopifnot(nrow(S) == 236L, !anyDuplicated(S$spec_id), all(table(S$block)[c("raw", "dca1", "ally", "sryc12", "srydist", "syc12", "sydist")] == c(56L, 28L, 56L, 24L, 24L, 24L, 24L))); S }
grid_parse <- function(spec_id) { p <- strsplit(spec_id, "|", fixed = TRUE)[[1]]; stopifnot(length(p) == 3, p[3] == "pr", startsWith(p[1], "G"))
  t <- strsplit(sub("^G", "", p[1]), "_", fixed = TRUE)[[1]]; stopifnot(length(t) == 3)
  list(block = t[1], igo = t[2], pta = t[3], layers = strsplit(p[2], "+", fixed = TRUE)[[1]]) }
grid_weight <- function(GL, spec_id) {   # GL: data.table from grid_layers.rds, sorted by year, ccode1, ccode2
  s <- grid_parse(spec_id); L <- s$layers; wtd <- s$block %in% c("sryc12", "srydist", "syc12", "sydist")
  stopifnot(any(c("igo", "def") %in% L), ("igo" %in% L) == (s$igo != "na"), ("pta" %in% L) == (s$pta != "na"), ("ally" %in% L) == (s$block == "ally"))
  col <- function(l) {
    if (wtd) { sch <- c(sryc12 = "sry_c12", srydist = "sry_dist", syc12 = "sy_c12", sydist = "sy_dist")[[s$block]]
      stopifnot(l != "ally", l != "igo" || s$igo == "count", l != "pta" || s$pta == "book"); return(paste0("s_", sch, "_", l)) }
    switch(l, igo = paste0("igo_", s$igo), def = "def", off = "off", dca = if (s$block == "dca1") "dca1" else "dca2",
           pta = { stopifnot(s$pta == "book"); "pta_book" }, ally = "othally") }
  cols <- vapply(L, col, ""); stopifnot(all(cols %in% names(GL)))
  comp <- rowMeans(as.matrix(GL[, cols, with = FALSE]))
  if (wtd) comp <- pmax(comp, 0) else stopifnot(all(comp >= 0))
  w <- numeric(length(comp)); for (ix in split(seq_along(comp), GL$year)) { m <- max(comp[ix]); if (m > 0) w[ix] <- comp[ix] / m }
  stopifnot(all(w >= 0 & w <= 1 + 1e-12)); pmin(w, 1) }
grid_edges <- function(GL, spec_id) { w <- grid_weight(GL, spec_id); k <- w > 0; data.table(a = GL$ccode1[k], b = GL$ccode2[k], year = GL$year[k], w = w[k]) }

## GPT P2-5: hashes of the code and layer file behind grid positions; written into every grid positions file by grid_build.R and verified by the
## sweep and the bootstrap before fitting, so positions can never pair with ties rebuilt from a different grid.R or grid_layers.rds
grid_provenance <- function(grid_layers_file, grid_code = file.path("R", "lib", "grid.R")) {
  stopifnot(file.exists(grid_layers_file), file.exists(grid_code))
  c(grid_R = unname(tools::md5sum(grid_code)), grid_layers = unname(tools::md5sum(grid_layers_file))) }
grid_check_provenance <- function(pz, grid_layers_file) {
  if (is.null(pz$provenance)) stop("grid positions file has no provenance (rebuild with R/grid_build.R)")
  now <- grid_provenance(grid_layers_file)
  if (!identical(pz$provenance[names(now)], now)) stop("grid positions were built from a different grid.R or grid_layers.rds than the current ones; rebuild step 1")
  invisible(TRUE) }
