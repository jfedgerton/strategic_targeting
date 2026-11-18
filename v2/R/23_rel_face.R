## ============================================================================
## 23_rel_face.R  --  Part 2c of the relational-measure sweep: plausibility of each candidate position, NO models
## (PI 2026-11-18: "genuine and real differences in relational embeddedness that seem plausible", extracted separately,
##  read before any coefficient is seen).
## Per spec (and the build's log IGO PageRank as BASELINE):
##   cor_baseline   mean within-year Spearman with the build's log IGO PageRank
##   cor_cinc       mean within-year Spearman with log CINC share (high = the measure re-measures material power)
##   stability      mean Spearman of a state's position between t-1 and t (states present in both years)
##   sd_within      median within-year SD of the position, divided by the median within-year SD over all specs of the same centrality
##   anchors        mean within-year percentile (1 = most central) for UK and France 1816-1913, USA 1946-2014, USSR 1946-1990,
##                  Germany 1991-2014 (expected high); Cuba 1962-1990, North Korea 1991-2014, Iran 1980-2014 (expected low)
##   target_gap_<type>  weighted mean over events of (chosen target's dR - mean dR of the event's choice set), rebuilt per measure
##                  with the dump's gap convention; > 0 = targets are LESS central than the alternatives
##   top10_by_era.csv   the ten states with the highest mean within-year percentile, per era
## Usage: Rscript v2/R/23_rel_face.R <long_dump_bc.rds> <out_dir> <positions_1.rds> [positions_2.rds ...]
## ============================================================================
suppressMessages({ library(dplyr) })
set.seed(123)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
DUMP <- args[1]; OUT_DIR <- args[2]; POSF <- args[-(1:2)]
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
logf <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
gap_scale <- function(v, yr) {   # verbatim from helpers_data.R
  stopifnot(length(v) == length(yr), !anyNA(v), all(is.finite(v)))
  ss <- tapply(v, yr, function(z) sum((z - mean(z))^2)); n <- tapply(v, yr, length)
  out <- sqrt(sum(2 * n * ss) / sum(n * (n - 1))); stopifnot(is.finite(out), out > 0); out
}

dump <- readRDS(DUMP); D <- dump$D; ROLES <- dump$roles
posy <- readRDS(file.path(dump$data_dir, "sy_positions.rds"))
P1 <- readRDS(POSF[1]); key <- P1$key
kk <- paste(key$ccode, key$year)
Rbase <- posy$R_pr[match(kk, paste(posy$ccode, posy$year))]; Mcinc <- posy$M[match(kk, paste(posy$ccode, posy$year))]
stopifnot(!anyNA(Rbase), !anyNA(Mcinc))
yrs <- split(seq_len(nrow(key)), key$year)
kidx <- setNames(seq_len(nrow(key)), kk)
i_row <- kidx[paste(D$c1, D$yr - 1L)]; j_row <- kidx[paste(D$c2, D$yr - 1L)]; stopifnot(!anyNA(i_row), !anyNA(j_row))
yr_ev <- as.character(D$yr)
prev <- kidx[paste(key$ccode, key$year - 1L)]                   # same state's previous-year row (NA in its first year)

ANCH <- list(UK_1816_1913 = c(200, 1816, 1913), FRA_1816_1913 = c(220, 1816, 1913), USA_1946_2014 = c(2, 1946, 2014), USSR_1946_1990 = c(365, 1946, 1990),
             GER_1991_2014 = c(255, 1991, 2014), CUB_1962_1990 = c(40, 1962, 1990), PRK_1991_2014 = c(731, 1991, 2014), IRN_1980_2014 = c(630, 1980, 2014))
ERAS <- list(`1816-1870` = c(1816, 1870), `1871-1913` = c(1871, 1913), `1914-1945` = c(1914, 1945), `1946-1989` = c(1946, 1989), `1990-2014` = c(1990, 2014))

face_one <- function(R, spec) {
  pct <- numeric(length(R)); for (ix in yrs) pct[ix] <- rank(R[ix]) / length(ix)
  cb <- mean(vapply(yrs, function(ix) suppressWarnings(cor(R[ix], Rbase[ix], method = "spearman")), 0), na.rm = TRUE)
  cc <- mean(vapply(yrs, function(ix) suppressWarnings(cor(R[ix], Mcinc[ix], method = "spearman")), 0), na.rm = TRUE)
  st <- { ok <- !is.na(prev); d <- data.frame(y = key$year[ok], a = R[ok], b = R[prev[ok]])
          mean(vapply(split(d, d$y), function(z) suppressWarnings(cor(z$a, z$b, method = "spearman")), 0), na.rm = TRUE) }
  sdw <- median(vapply(yrs, function(ix) sd(R[ix]), 0))
  an <- vapply(ANCH, function(a) { s <- key$ccode == a[1] & key$year >= a[2] & key$year <= a[3]; if (any(s)) mean(pct[s]) else NA_real_ }, 0)
  ys <- vapply(yrs, function(ix) gap_scale(R[ix], rep(1L, length(ix))), 0)
  dR <- if (dump$rel_scale == "year") (R[i_row] - R[j_row]) / ys[yr_ev] else (R[i_row] - R[j_row]) / gap_scale(R, key$year)
  tg <- tapply(seq_len(nrow(D)), D$ev_id, function(r) dR[r][D$chosen[r] == 1L] - mean(dR[r]))
  ew <- tapply(D$w, D$ev_id, `[`, 1); et <- tapply(D$type, D$ev_id, `[`, 1)
  tgap <- vapply(seq_along(ROLES) - 1L, function(g) { s <- et == g; sum(tg[s] * ew[s]) / sum(ew[s]) }, 0); names(tgap) <- paste0("target_gap_", ROLES)
  top <- do.call(rbind, lapply(names(ERAS), function(e) { s <- key$year >= ERAS[[e]][1] & key$year <= ERAS[[e]][2]
    m <- tapply(pct[s], key$ccode[s], mean); n <- tapply(pct[s], key$ccode[s], length); m <- m[n >= 5]
    o <- order(-m)[1:min(10, length(m))]; data.frame(spec_id = spec, era = e, rank = seq_along(o), ccode = as.integer(names(m)[o]), mean_pct = round(m[o], 3)) }))
  list(row = data.frame(spec_id = spec, cor_baseline = cb, cor_cinc = cc, stability = st, sd_within = sdw, t(an), t(tgap), stringsAsFactors = FALSE), top = top)
}

rows <- list(); tops <- list()
b <- face_one(Rbase, "BASELINE|igo_pagerank_build|pr"); rows[[1]] <- b$row; tops[[1]] <- b$top
for (f in POSF) {
  P <- readRDS(f); stopifnot(identical(paste(P$key$ccode, P$key$year), kk))
  t0 <- Sys.time()
  used_years <- sort(unique(c(D$yr, D$yr - 1L)))
  for (k in seq_len(ncol(P$R))) {
    cy <- as.integer(strsplit(P$specs$const_years[k], " ")[[1]])
    bad_y <- intersect(cy, used_years)
    if (!P$specs$identified[k] || length(bad_y)) {   # review P1-1 / P2-13: listed, not silently dropped
      rows[[length(rows) + 1]] <- data.frame(spec_id = colnames(P$R)[k], skipped = sprintf("constant position in %d years (%d used by events)", length(cy), length(bad_y)), stringsAsFactors = FALSE)
      next }
    o <- face_one(P$R[, k], colnames(P$R)[k]); rows[[length(rows) + 1]] <- o$row; tops[[length(tops) + 1]] <- o$top }
  logf("%s: %d specs in %.1f min", basename(f), ncol(P$R), as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
fs <- bind_rows(rows)
fs$centrality <- sub(".*\\|", "", fs$spec_id)
if (!"skipped" %in% names(fs)) fs$skipped <- NA_character_
## review P2-13: sd_within relative to a FIXED reference (the W1raw|igo spec of the same centrality), not the median of the files passed
ref_sd <- c(pr = fs$sd_within[fs$spec_id == "W1raw|igo|pr"], deg = fs$sd_within[fs$spec_id == "W1raw|igo|deg"])
if (length(ref_sd) == 2) fs$sd_within_rel <- fs$sd_within / ref_sd[fs$centrality] else fs$sd_within_rel <- NA_real_
write.csv(fs, file.path(OUT_DIR, "face_summary.csv"), row.names = FALSE)
write.csv(bind_rows(tops), file.path(OUT_DIR, "top10_by_era.csv"), row.names = FALSE)
logf("wrote face_summary.csv (%d specs, %d listed as skipped) and top10_by_era.csv", nrow(fs), sum(!is.na(fs$skipped)))
