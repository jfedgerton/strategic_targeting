## R/grid_summary.R -- appendix summaries of the network-measure grid sweep (PI 2026-11-24; Tremors-style decision regressions).
## Reads the R/s2_05_sweep.R outputs (one folder per sample), keeps status OK cells, and for each headline quantity:
##   (a) distribution across measures (share with a 95% interval excluding 0, by sign), with the main measure marked
##   (b) decision regressions: the estimate and the two-sided p-value regressed on the measurement choices (layer included: igo, def, off,
##       dca, pta, ally; IGO coding: cos, hyper vs count; PTA coding: depth vs book; DCA v1; weighting scheme vs raw; sample mid vs fmid),
##       OLS with HC1 standard errors; coefficient plots (one panel per quantity)
##   (c) specification curve of the estimate (sorted), with the choices below
## Headline quantities: S1NF coal R +1 (coalition, less central target, between targets), S1NF bil R +1, S1NF coal M +1 / -1, theta A, theta B,
##   S2 B M -1 (joining the target vs a stronger attacker), S2 B R +1.
## Usage: Rscript R/grid_summary.R --sweeps=DIR1,DIR2 --main_spec=SPEC --out=DIR
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(sandwich) })
source("R/lib/grid.R")
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
SW <- strsplit(opt("sweeps"), ",")[[1]]; MAIN <- opt("main_spec"); OUT <- opt("out"); stopifnot(length(SW) > 0, !is.null(MAIN), !is.null(OUT)); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
fs <- unlist(lapply(SW, function(d) list.files(d, pattern = "\\.csv$", full.names = TRUE))); fs <- fs[!grepl("census|log", basename(fs))]
A <- rbindlist(lapply(fs, fread), fill = TRUE); stopifnot(all(c("spec_id", "status", "sample", "stage", "alt", "axis", "gap", "est", "se") %in% names(A)))
cen <- unique(A[, .(spec_id, sample, status)]); fwrite(cen[, .N, by = .(sample, status)], file.path(OUT, "grid_census.csv"))
A <- A[status == "OK"]
Q <- data.table(stage = c("S1NF", "S1NF", "S1NF", "S1NF", "S2", "S2", "S2", "S2"), alt = c("coal", "bil", "coal", "coal", "A", "B", "B", "B"),
                axis = c("R", "R", "M", "M", "theta", "theta", "M", "R"), gap = c(1, 1, 1, -1, NA, NA, -1, 1),
                label = c("Coalition: less central target", "Bilateral: less central target", "Coalition: weaker target", "Coalition: stronger target",
                          "theta A (join attacker)", "theta B (join target)", "Join target vs stronger attacker", "Join target vs less central attacker"))
K <- merge(A, Q, by = c("stage", "alt", "axis", "gap"))
K[, `:=`(lo = est - 1.96 * se, hi = est + 1.96 * se, p = 2 * pnorm(-abs(est / se)))]
## decision indicators
dec <- function(sp) { if (startsWith(sp, "G")) { g <- grid_parse(sp); L <- g$layers; b <- g$block } else { L <- strsplit(strsplit(sp, "|", fixed = TRUE)[[1]][2], "+", fixed = TRUE)[[1]]; b <- if (startsWith(sp, "BOOKraw")) "raw" else "other"; g <- list(igo = "count", pta = "book") }
  data.table(spec_id = sp, block = b, inc_igo = "igo" %in% L, inc_def = "def" %in% L, inc_off = "off" %in% L, inc_dca = "dca" %in% L, inc_pta = "pta" %in% L, inc_ally = "ally" %in% L,
             igo_cos = g$igo == "cos", igo_hyper = g$igo == "hyper", pta_depth = g$pta == "depth", dca_v1 = b == "dca1",
             w_sryc12 = b == "sryc12", w_srydist = b == "srydist", w_syc12 = b == "syc12", w_sydist = b == "sydist") }
DEC <- rbindlist(lapply(unique(K$spec_id), dec)); K <- merge(K, DEC, by = "spec_id"); K[, mid := sample == "mid"][, is_main := spec_id == MAIN]
## (a) distributions
S <- K[, .(measures = .N, median = median(est), share_ci_pos = mean(lo > 0), share_ci_neg = mean(hi < 0), main_est = est[is_main][1]), by = .(label, sample)]
fwrite(S, file.path(OUT, "grid_distribution.csv"))
## (b) decision regressions
X <- c("inc_igo", "inc_def", "inc_off", "inc_dca", "inc_pta", "inc_ally", "igo_cos", "igo_hyper", "pta_depth", "dca_v1", "w_sryc12", "w_srydist", "w_syc12", "w_sydist", "mid")
R <- list()
for (lb in unique(K$label)) for (y in c("est", "p")) { d <- K[label == lb]; xs <- X[vapply(X, function(v) uniqueN(d[[v]]) > 1, TRUE)]
  m <- lm(as.formula(paste(y, "~", paste(xs, collapse = " + "))), data = d); V <- vcovHC(m, type = "HC1"); b <- coef(m)[-1]; s <- sqrt(diag(V))[-1]
  R[[length(R) + 1]] <- data.table(label = lb, outcome = y, term = sub("TRUE$", "", names(b)), b = b, se = s, n = nobs(m)) }
R <- rbindlist(R); fwrite(R, file.path(OUT, "grid_decision_regressions.csv"))
nice <- c(inc_igo = "IGO layer", inc_def = "Defense alliances", inc_off = "Offense alliances", inc_dca = "DCA", inc_pta = "PTA", inc_ally = "Other alliances",
          igo_cos = "IGO cosine", igo_hyper = "IGO hyperedge", pta_depth = "PTA depth", dca_v1 = "DCA v1", w_sryc12 = "Weighted sry contig", w_srydist = "Weighted sry distance",
          w_syc12 = "Weighted sy contig", w_sydist = "Weighted sy distance", mid = "All MIDs (vs fatal)")
R[, term_lab := factor(nice[term], levels = rev(nice))]
th <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), strip.placement = "outside")
for (y in c("est", "p")) { p <- ggplot(R[outcome == y], aes(b, term_lab)) + geom_vline(xintercept = 0, linetype = 2) +
    geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0) + geom_point(size = 1.3) + facet_wrap(~ label, ncol = 4, scales = "free_x") +
    labs(x = if (y == "est") "Change in the estimate" else "Change in the p-value", y = NULL) + th
  ggsave(file.path(OUT, sprintf("grid_decisions_%s.pdf", y)), p, width = 11, height = 6.5); ggsave(file.path(OUT, sprintf("grid_decisions_%s.png", y)), p, width = 11, height = 6.5, dpi = 200) }
## (c) specification curves
K[, rk := frank(est, ties.method = "first"), by = .(label, sample)]
p <- ggplot(K[sample == "fmid"], aes(rk, est)) + geom_hline(yintercept = 0, linetype = 2) + geom_linerange(aes(ymin = lo, ymax = hi), colour = "grey70") +
  geom_point(aes(colour = is_main), size = .7) + scale_colour_manual(values = c(`FALSE` = "grey20", `TRUE` = "firebrick"), labels = c("Grid measure", "Main measure"), name = NULL) +
  facet_wrap(~ label, ncol = 4, scales = "free") + labs(x = "Measures, ranked by estimate (fatal MIDs)", y = "Estimate with 95% interval") + th + theme(legend.position = "bottom")
ggsave(file.path(OUT, "grid_spec_curves.pdf"), p, width = 11, height = 6.5); ggsave(file.path(OUT, "grid_spec_curves.png"), p, width = 11, height = 6.5, dpi = 200)
cat("GRID_SUMMARY_DONE", OUT, nrow(K), "rows\n")
