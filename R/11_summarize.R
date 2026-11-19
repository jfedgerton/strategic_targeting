## ============================================================================
## 30_summarize_loop.R  --  appendix tables from every sweep output of the multiverse loop (2026-11-19)
## Reads every sweep_*.csv under <sweep_dir> (recursively) and the measure validation table, and writes into <out_dir>:
##   summary_long.csv   one row per run x quantity: main-measure estimate and CI, number of valid measures, % of measures whose 95% CI is
##                      above / below 0, median and IQR across measures -- for ALL measures and for the PLAUSIBLE ones separately
##   summary_main.md    the headline quantities per run in readable tables (no interpretation)
## A run = one sweep file (scheme tag x positions file); runs of the same scheme tag on different positions files are pooled (multiverse).
## Usage: Rscript v2/R/30_summarize_loop.R <sweep_dir> <measure_validation.csv> <main_spec_id> <out_dir>
## ============================================================================
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) == 4)
SW <- args[1]; VAL <- args[2]; MAIN <- args[3]; OUT <- args[4]; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
V <- fread(VAL)[, .(spec_id, plausible)]
files <- list.files(SW, pattern = "^sweep_.*\\.csv$", recursive = TRUE, full.names = TRUE); stopifnot(length(files) > 0)
tag_of <- function(f) sub("_(positions|phaseA|pick|noatop).*$", "", sub("^sweep_", "", basename(f)))
grp_of <- function(f) basename(dirname(f))
X <- rbindlist(lapply(files, function(f) { d <- fread(f); d[, `:=`(run = tag_of(f), group = grp_of(f))]; d }), fill = TRUE)
X <- X[!grepl("^BASELINE", spec_id)]
## census per run
CEN <- X[, .(measures = uniqueN(spec_id), invalid = uniqueN(spec_id[model == "INVALID"]), skipped = uniqueN(spec_id[model == "SKIPPED"]), errors = uniqueN(spec_id[model == "ERROR"])), by = .(group, run)]
E <- X[model %in% c("M0", "M0t", "M1", "B0", "B1") & !is.na(se) & df == 1L]
E <- merge(E, V, by = "spec_id", all.x = TRUE)
summ <- function(d) d[, .(n = .N, main = if (any(spec_id == MAIN)) sprintf("%.2f [%.2f, %.2f]", est[spec_id == MAIN][1], lo[spec_id == MAIN][1], hi[spec_id == MAIN][1]) else NA_character_,
                         pct_pos = round(100 * mean(lo > 0), 1), pct_neg = round(100 * mean(hi < 0), 1), median = round(median(est), 3),
                         q25 = round(quantile(est, .25), 3), q75 = round(quantile(est, .75), 3)), by = .(group, run, model, quantity)]
S_all <- summ(E)[, set := "all"]
S_pl <- summ(E[plausible == TRUE])[, set := "plausible"]
W <- X[model %in% c("M1", "B1") & df > 1L, .(n = .N, main = if (any(spec_id == MAIN)) sprintf("chi2 %.1f (%d), p %.3g", est[spec_id == MAIN][1], df[spec_id == MAIN][1], p[spec_id == MAIN][1]) else NA_character_,
                                             pct_p05 = round(100 * mean(p < 0.05, na.rm = TRUE), 1)), by = .(group, run, model, quantity)]
OUTL <- rbind(S_all, S_pl, fill = TRUE)
fwrite(OUTL, file.path(OUT, "summary_long.csv")); fwrite(W, file.path(OUT, "summary_wald.csv")); fwrite(CEN, file.path(OUT, "census.csv"))
## readable headline tables
key_q <- c("u_R_less", "u_M_stronger", "u_M_weaker")
hd <- OUTL[set == "all" & ((model %in% c("M0", "B0") & quantity %in% key_q) | (model %in% c("M1", "B1") & grepl(paste0(":(", paste(key_q, collapse = "|"), ")$"), quantity)))]
con <- file(file.path(OUT, "summary_main.md"), "w")
writeLines(sprintf("# Multiverse loop summary\n\nMain measure: `%s`. Cells: main-measure estimate [95%% CI] | %% of measures with CI > 0 / < 0 | n measures.\n", MAIN), con)
for (g in unique(hd$group)) for (r in unique(hd[group == g, run])) {
  d <- hd[group == g & run == r]; writeLines(sprintf("\n## %s / %s\n\n| model | quantity | main | %%>0 | %%<0 | n |\n|---|---|---|---|---|---|", g, r), con)
  for (i in seq_len(nrow(d))) writeLines(sprintf("| %s | %s | %s | %s | %s | %d |", d$model[i], d$quantity[i], d$main[i], d$pct_pos[i], d$pct_neg[i], d$n[i]), con)
  w <- W[group == g & run == r]; if (nrow(w)) { writeLines("\n| Wald test | main | % p<.05 |\n|---|---|---|", con)
    for (i in seq_len(nrow(w))) writeLines(sprintf("| %s %s | %s | %s |", w$model[i], w$quantity[i], w$main[i], w$pct_p05[i]), con) }
}
close(con)
cat(sprintf("summarized %d sweep files (%d runs) into %s\n", length(files), nrow(CEN), OUT))
