## ============================================================================
## 23_rel_validate.R  --  measurement validation of every relational position (2026-11-19; replaces 23_rel_face.R in the multiverse loop)
## OUTCOME-BLIND: uses positions, COW major-power status and CINC only; no event, target or model output is read.
## Part A, the plausibility rule set (PI 2026-11-18; the rule behind plausibility_rules_20261118.csv, reproduced exactly):
##   percentile = rank within year / number of states that year (ties averaged)
##   USA_cw, USSR_cw = mean percentile 1950-1989 (US 2, USSR/Russia 365)          rule: both >= 0.5
##   CHN_rise        = mean percentile 1990-2014 minus 1950-1975 (China 710)        rule: > 0
##   WEU             = mean of the UK (200), France (220), Germany (255) means 1816-1913   rule: >= 0.5
##   PRK_post        = mean percentile 1991-2014 (North Korea 731)                  rule: < 0.9
##   (IND_rise is computed and reported but is not part of the rule; the India rule was dropped on 2026-11-18)
## Part B, "does the measure put the powerful states where they belong?" (PI 2026-11-19), by era (1816-1913, 1919-1938, 1946-1989, 1990-2014):
##   rho_cinc  = mean over years of the within-year Spearman correlation between the position and CINC
##   maj_pct   = mean within-year percentile of COW major powers
##   maj_auc   = mean within-year probability that a random major power outranks a random non-major (AUC)
##   stab      = mean Spearman correlation of the position between consecutive years (states present in both)
## Part C, top-10 states by mean within-year percentile per era, for the measures listed with --top (default: the first measure of each file).
## Usage: Rscript v2/R/23_rel_validate.R <out_dir> <positions_1.rds> [<positions_2.rds> ...] [--top=<file with spec ids>] [--cores=N]
## ============================================================================
suppressMessages({ library(data.table); library(parallel); library(peacesciencer) })
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
fl <- sub("^--([A-Za-z_]+)=.*$", "\\1", grep("^--", args, value = TRUE)); if (any(!fl %in% c("top", "cores"))) stop("unknown option(s): ", paste(setdiff(fl, c("top", "cores")), collapse = ", "))
pos <- args[!grepl("^--", args)]; stopifnot(length(pos) >= 2)
OUT <- pos[1]; FILES <- pos[-1]; stopifnot(all(file.exists(FILES)))
CORES <- as.integer(opt("cores", "1")); TOPF <- opt("top", "")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
logf <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
RULE <- "US & USSR 1950-89 >= .5; China rises 1950-75 -> 1990-2014; mean(UK,FRA,GER) 1816-1913 >= .5; PRK 1991-2014 < .9 (PI 2026-11-18, India rule dropped)"
ERAS <- list(e1816_1913 = 1816:1913, e1919_1938 = 1919:1938, e1946_1989 = 1946:1989, e1990_2014 = 1990:2014)

sy <- as.data.table(peacesciencer::add_cow_majors(peacesciencer::add_nmc(peacesciencer::create_stateyears(system = "cow", subset_years = 1816:2014))))
sy <- sy[, .(ccode = as.integer(ccode), year = as.integer(year), cinc, maj = as.integer(cowmaj))]
stopifnot(!anyDuplicated(sy[, .(ccode, year)]))

one_measure <- function(R, key) {
  dt <- data.table(ccode = key$ccode, year = key$year, R = R)
  dt[, pct := frank(R, ties.method = "average") / .N, by = year]
  m <- function(cc, y) dt[ccode == cc & year %in% y, mean(pct)]
  a <- c(USA_cw = m(2, 1950:1989), USSR_cw = m(365, 1950:1989), IND_rise = m(750, 1990:2014) - m(750, 1950:1975),
         CHN_rise = m(710, 1990:2014) - m(710, 1950:1975), UK_early = m(200, 1816:1913), FRA_early = m(220, 1816:1913), GER_early = m(255, 1816:1913),
         PRK_post = m(731, 1991:2014))
  a["WEU"] <- mean(a[c("UK_early", "FRA_early", "GER_early")])
  x <- merge(dt, sy, by = c("ccode", "year"), all.x = TRUE)
  yr <- x[, .(rho = if (sum(!is.na(cinc)) > 3 && sd(R) > 0) suppressWarnings(cor(R, cinc, method = "spearman", use = "complete.obs")) else NA_real_,
              mp = if (any(maj %in% 1L)) mean(pct[maj %in% 1L]) else NA_real_,
              auc = { k <- maj %in% 1L; if (any(k) && any(!k)) { r <- frank(R, ties.method = "average"); (sum(r[k]) - sum(k) * (sum(k) + 1) / 2) / (sum(k) * sum(!k)) } else NA_real_ }), by = year]
  setkey(dt, ccode, year)
  st <- dt[, .(year, ccode, R)]; st2 <- copy(st)[, year := year + 1L]
  sb <- merge(st, st2, by = c("ccode", "year"), suffixes = c("", "_prev"))[, .(s = if (.N > 3) suppressWarnings(cor(R, R_prev, method = "spearman")) else NA_real_), by = year]
  b <- unlist(lapply(names(ERAS), function(e) { y <- ERAS[[e]]
    c(setNames(mean(yr[year %in% y, rho], na.rm = TRUE), paste0("rho_cinc_", e)), setNames(mean(yr[year %in% y, mp], na.rm = TRUE), paste0("maj_pct_", e)),
      setNames(mean(yr[year %in% y, auc], na.rm = TRUE), paste0("maj_auc_", e)), setNames(mean(sb[year %in% y, s], na.rm = TRUE), paste0("stab_", e))) }))
  c(a, b)
}
top10 <- function(R, key, spec) {
  dt <- data.table(ccode = key$ccode, year = key$year, R = R)[, pct := frank(R, ties.method = "average") / .N, by = year]
  rbindlist(lapply(names(ERAS), function(e) dt[year %in% ERAS[[e]], .(mean_pct = mean(pct), years = .N), by = ccode][order(-mean_pct)][1:10][, `:=`(era = e, rank = 1:10, spec_id = spec)]))
}
want_top <- if (nzchar(TOPF)) readLines(TOPF) else character(0)

res <- list(); tops <- list()
for (f in FILES) {
  P <- readRDS(f); key <- as.data.frame(P$key); stopifnot(nrow(key) == nrow(P$R), !anyNA(P$R))
  ids <- colnames(P$R); t0 <- Sys.time()
  rows <- mclapply(seq_along(ids), function(k) one_measure(P$R[, k], key), mc.cores = CORES)
  bad <- !vapply(rows, is.numeric, TRUE); if (any(bad)) stop(sprintf("%s: %d measures failed", basename(f), sum(bad)))
  r <- as.data.table(do.call(rbind, rows)); r[, `:=`(spec_id = ids, file = basename(f))]
  res[[f]] <- r
  tk <- if (length(want_top)) intersect(want_top, ids) else ids[1]
  for (s in tk) tops[[paste(f, s)]] <- top10(P$R[, s], key, s)
  logf("%s: %d measures in %.1f min", basename(f), length(ids), as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
A <- rbindlist(res)
A[, plausible := USA_cw >= 0.5 & USSR_cw >= 0.5 & CHN_rise > 0 & WEU >= 0.5 & PRK_post < 0.9]; A[, rule := RULE]
stopifnot(!anyDuplicated(A$spec_id), !anyNA(A$plausible))
cc <- countrycode::countrycode
TP <- rbindlist(tops)[, state := cc(ccode, "cown", "country.name", warn = FALSE)]
fwrite(TP, file.path(OUT, "top10_by_era.csv")); fwrite(A, file.path(OUT, "measure_validation.csv.tmp")); stopifnot(file.rename(file.path(OUT, "measure_validation.csv.tmp"), file.path(OUT, "measure_validation.csv")))   # atomic
logf("wrote %s: %d measures, %d plausible (%.1f%%)", file.path(OUT, "measure_validation.csv"), nrow(A), sum(A$plausible), 100 * mean(A$plausible))
