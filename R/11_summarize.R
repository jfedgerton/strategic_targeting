## ============================================================================
## 30_summarize_loop.R  --  appendix tables for the multiverse loop (2026-11-19; rewritten after the loop review P2-13 and C)
## Reads ONLY the sweep outputs listed in <manifest.tsv> (written by run_all.sh: group, path, n_requested), rejects duplicate keys and
## unexpected or missing spec ids, and writes into <out_dir>:
##   census.csv          per run x model: requested / valid (estimate rows) / invalid / skipped / errored measures
##   estimates_long.csv  per run x model x quantity x validity set (all / plausible / failing): main-measure estimate and CI, n, median, IQR,
##                       % of measures whose CI is above / below 0 (row-weighted GRID summary), and the family-balanced share (each family cell
##                       = weighting x centrality x diplomacy arm counts equally)
##   by_family.csv       the same per family cell (weighting, centrality, diplomacy, ATOP in measure, PTA coding)
##   wald.csv            Wald tests: main-measure statistic, % p < .05 among computable tests, n computable / n requested
##   summary_main.md     readable tables (estimates first; percentages are labelled as grid summaries)
##   (delta review N6) the manifest carries, per output, the file of EXPECTED spec ids; the ids seen in each group must equal the union of its
##   expected ids exactly (not just in number), and every id must be in the validation file. Wald tests report requested (the expected ids of
##   the run) / fit-valid (the model's fit is valid) / computable (finite p) separately. M0t has its own census columns. Outputs written atomically.
## Usage: Rscript v2/R/30_summarize_loop.R <manifest.tsv> <measure_validation.csv> <main_spec_id> <out_dir>
## ============================================================================
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) == 4)
MAN <- args[1]; VAL <- args[2]; MAIN <- args[3]; OUT <- args[4]; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
V <- fread(VAL)[, .(spec_id, plausible)]; stopifnot(!anyDuplicated(V$spec_id))
M <- fread(MAN, sep = "\t", header = TRUE); stopifnot(all(c("group", "path", "n_requested", "ids_file") %in% names(M)), !anyDuplicated(M$path))
miss <- M$path[!file.exists(M$path)]; if (length(miss)) stop(sprintf("%d manifest outputs are missing, e.g. %s", length(miss), miss[1]))
run_of <- function(f) sub("_(fmid|mid)_.*$", "_\\1", sub("\\.csv$", "", sub("^sweep_", "", basename(f))))
M[, run := run_of(path)]
EXP <- rbindlist(lapply(seq_len(nrow(M)), function(i) { ids <- readLines(M$ids_file[i]); stopifnot(length(ids) == M$n_requested[i], !anyDuplicated(ids))
  data.table(group = M$group[i], run = M$run[i], path = M$path[i], spec_id = ids) }))
unk <- setdiff(unique(EXP$spec_id), V$spec_id); if (length(unk)) stop(sprintf("%d requested ids are not in the validation file, e.g. %s", length(unk), unk[1]))
X <- rbindlist(lapply(seq_len(nrow(M)), function(i) { d <- fread(M$path[i]); d[, `:=`(group = M$group[i], run = M$run[i], file = M$path[i])]; d }), fill = TRUE)
X <- X[!grepl("^BASELINE", spec_id)]
## exact id sets per output file
for (f in unique(M$path)) { s <- unique(X$spec_id[X$file == f]); e <- EXP$spec_id[EXP$path == f]
  if (!setequal(s, e)) stop(sprintf("%s: analysed ids differ from the expected ids (%d unexpected, %d missing)", f, length(setdiff(s, e)), length(setdiff(e, s)))) }
## family of each measure (from its id): weighting, PTA coding, layer set, centrality
fam <- function(id) { p <- tstrsplit(id, "|", fixed = TRUE); w <- sub("_pta.*$", "", p[[1]]); pta <- ifelse(grepl("_pta", p[[1]]), sub("^.*_", "", p[[1]]), "none")
  data.table(spec_id = id, weighting = w, pta = pta, centrality = p[[3]], dipl = grepl("(^|\\+)dipl(\\+|$)", p[[2]]),
             atop = grepl("(^|\\+)(def|off|neu|nonagg|consul)(\\+|$)", p[[2]])) }
FM <- fam(unique(X$spec_id)); FM[, cell := paste(weighting, centrality, ifelse(dipl, "dipl", "nodipl"), sep = "/")]
FM <- merge(FM, V, by = "spec_id", all.x = TRUE); stopifnot(!anyNA(FM$plausible))
## valid fit of model m = the spec has an estimate row of m other than the rows written for every spec (iter, max_abs_dR, merror_rank_cor)
NONEST <- c("iter", "max_abs_dR", "merror_rank_cor")
nvalid <- function(d, m) uniqueN(d$spec_id[d$model %in% m & !d$quantity %in% NONEST])
CEN <- X[, .(n_valid_M0 = nvalid(.SD, c("M0", "B0")), n_valid_M1 = nvalid(.SD, c("M1", "B1")), n_valid_M0t = nvalid(.SD, "M0t"),
             n_invalid_M0 = uniqueN(spec_id[model == "INVALID" & grepl("^(M0|B0) fit invalid", quantity)]),
             n_invalid_M1 = uniqueN(spec_id[model == "INVALID" & grepl("^(M1|B1) fit invalid", quantity)]),
             n_invalid_M0t = uniqueN(spec_id[model == "INVALID" & grepl("^M0t fit invalid", quantity)]),
             n_skipped = uniqueN(spec_id[model == "SKIPPED"]), n_error = uniqueN(spec_id[model == "ERROR"]), n_seen = uniqueN(spec_id)), by = group]
CEN <- merge(CEN, EXP[, .(n_requested = .N), by = group], by = "group")
stopifnot(all(CEN$n_seen == CEN$n_requested))
E <- X[model %in% c("M0", "M0t", "M1", "B0", "B1") & !is.na(se) & (is.na(df) | df == 1L)]
dup <- E[, .N, by = .(group, spec_id, model, quantity)][N > 1L]; if (nrow(dup)) stop(sprintf("%d duplicate estimate keys, e.g. %s %s %s", nrow(dup), dup$group[1], dup$spec_id[1], dup$quantity[1]))
E <- merge(E, FM, by = "spec_id", all.x = TRUE)
summ <- function(d, by) d[, .(n = .N, main = if (any(spec_id == MAIN)) sprintf("%.2f [%.2f, %.2f]", est[spec_id == MAIN][1], lo[spec_id == MAIN][1], hi[spec_id == MAIN][1]) else NA_character_,
                             median = round(median(est), 3), q25 = round(quantile(est, .25), 3), q75 = round(quantile(est, .75), 3),
                             grid_pct_pos = round(100 * mean(lo > 0), 1), grid_pct_neg = round(100 * mean(hi < 0), 1),
                             family_bal_pct_pos = round(100 * mean(tapply(lo > 0, cell, mean)), 1), n_cells = uniqueN(cell)), by = by]
B <- c("group", "run", "model", "quantity")
EL <- rbind(summ(E, B)[, set := "all"], summ(E[plausible == TRUE], B)[, set := "plausible"], summ(E[plausible == FALSE], B)[, set := "failing"], fill = TRUE)
BF <- E[, .(n = .N, median = round(median(est), 3), grid_pct_pos = round(100 * mean(lo > 0), 1), grid_pct_neg = round(100 * mean(hi < 0), 1)),
        by = c(B, "weighting", "centrality", "dipl", "atop", "pta")]
W <- X[model %in% c("M1", "B1") & !is.na(df) & df > 1L]
## denominators taken BEFORE dropping invalid fits: requested = expected ids of the run; fit_valid = valid fits of that model in the run
REQ <- EXP[, .(n_requested = .N), by = .(group, run)]
FV <- X[model %in% c("M1", "B1") & !quantity %in% NONEST, .(n_fit_valid = uniqueN(spec_id)), by = .(group, run, model)]
WS <- W[, .(n_computable = sum(!is.na(p)), main = if (any(spec_id == MAIN)) sprintf("chi2 %.1f (%d), p %.3g", est[spec_id == MAIN][1], df[spec_id == MAIN][1], p[spec_id == MAIN][1]) else NA_character_,
            pct_p05_of_computable = round(100 * mean(p < 0.05, na.rm = TRUE), 1)), by = B]
WS <- merge(merge(WS, REQ, by = c("group", "run")), FV, by = c("group", "run", "model"), all.x = TRUE); WS[is.na(n_fit_valid), n_fit_valid := 0L]
stopifnot(all(WS$n_computable <= WS$n_fit_valid), all(WS$n_fit_valid <= WS$n_requested))
afw <- function(x, f) { fwrite(x, paste0(f, ".tmp")); stopifnot(file.rename(paste0(f, ".tmp"), f)) }   # atomic (delta review N5)
afw(CEN, file.path(OUT, "census.csv")); afw(EL, file.path(OUT, "estimates_long.csv")); afw(BF, file.path(OUT, "by_family.csv")); afw(WS, file.path(OUT, "wald.csv"))
## readable tables: main-measure estimate first; grid and family-balanced percentages as secondary columns
key_q <- "u_R_less|u_M_stronger|u_M_weaker"
hd <- EL[set == "all" & grepl(paste0("(^|:)(", key_q, ")$"), quantity)]
con <- file(file.path(OUT, "summary_main.md.tmp"), "w")
writeLines(sprintf("# Multiverse loop summary\n\nMain measure `%s`. Columns: main-measure estimate [95%% CI]; median and IQR across measures; GRID %% of measures with CI > 0 / < 0 (correlated constructions, a descriptive summary, not replications); FAMILY-BALANCED %% (each weighting x centrality x diplomacy cell counts equally).\n", MAIN), con)
for (g in unique(hd$group)) { d <- hd[group == g]
  writeLines(sprintf("\n## %s\n\n| model | quantity | main | median [IQR] | grid %%>0 / %%<0 | family-bal %%>0 | n |\n|---|---|---|---|---|---|---|", g), con)
  for (i in seq_len(nrow(d))) writeLines(sprintf("| %s | %s | %s | %.2f [%.2f, %.2f] | %s / %s | %s | %d |", d$model[i], d$quantity[i], d$main[i], d$median[i], d$q25[i], d$q75[i], d$grid_pct_pos[i], d$grid_pct_neg[i], d$family_bal_pct_pos[i], d$n[i]), con)
  w <- WS[group == g]; if (nrow(w)) { writeLines("\n| Wald test | main | % p<.05 (computable) | computable / fit-valid / requested |\n|---|---|---|---|", con)
    for (i in seq_len(nrow(w))) writeLines(sprintf("| %s %s | %s | %s | %d / %d / %d |", w$model[i], w$quantity[i], w$main[i], w$pct_p05_of_computable[i], w$n_computable[i], w$n_fit_valid[i], w$n_requested[i]), con) }
}
close(con); stopifnot(file.rename(file.path(OUT, "summary_main.md.tmp"), file.path(OUT, "summary_main.md")))
cat(sprintf("summarized %d manifest outputs (%d groups) into %s\n", nrow(M), uniqueN(M$group), OUT))
