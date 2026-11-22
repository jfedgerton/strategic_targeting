## R/s2_04_tables.R -- LaTeX tables for the current stage-1, stage-2 and two-step FE joint models (PI 2026-11-22: "what a current table would
## look like"). Reads R/s2_03_joint_fe.R output (stage 1 without / with target FE, stage 2 with theta, bootstrap) and R/s2_02_fit_reference.R
## output (stage-2 reference). Writes tables.tex (booktabs; \input-able) and tables_standalone.tex (compiles on its own).
## Usage: Rscript R/s2_04_tables.R --joint=DIR --ref=DIR --out=DIR
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name) sub(paste0("^--", name, "="), "", grep(paste0("^--", name, "="), args, value = TRUE))
JF <- readRDS(file.path(opt("joint"), "joint_fe.rds")); REF <- fread(file.path(opt("ref"), "reference_coefs.csv")); REFC <- fread(file.path(opt("ref"), "reference_surface_contrasts.csv"))
OUT <- opt("out"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
stars <- function(p) ifelse(is.na(p), "", ifelse(p < .01, "^{***}", ifelse(p < .05, "^{**}", ifelse(p < .10, "^{*}", ""))))
cell <- function(est, se) { if (is.na(est)) return(c("", "")); p <- 2 * pnorm(-abs(est / se)); c(sprintf("$%.2f%s$", est, stars(p)), sprintf("(%.2f)", se)) }
cell_ci <- function(est, lo, hi) if (is.na(est)) c("", "") else c(sprintf("$%.2f%s$", est, if (lo > 0 || hi < 0) "^{\\dagger}" else ""), sprintf("[%.2f, %.2f]", lo, hi))
row2 <- function(label, cells) { k <- length(cells); paste0(label, " & ", paste(vapply(cells, `[`, "", 1), collapse = " & "), " \\\\\n", " & ", paste(vapply(cells, `[`, "", 2), collapse = " & "), " \\\\\n") }
CONTR_LAB <- c("M.1" = "Opponent weaker (gap $+1$)", "M.-1" = "Opponent stronger (gap $-1$)", "R.1" = "Opponent less central (gap $+1$)", "R.-1" = "Opponent more central (gap $-1$)")

## ---- Table 1: stage 1 -------------------------------------------------------------------------------------------------------------
S1 <- JF$S1; S1C <- JF$S1C
cols1 <- list(c("bil", "no FE"), c("coal", "no FE"), c("bil", "target FE"), c("coal", "target FE"))
t1 <- "\\begin{table}[!htbp]\\centering\\small\n\\caption{Stage 1: target choice by initiation type (conditional logit)}\\label{tab:stage1}\n\\begin{tabular}{lcccc}\n\\toprule\n"
t1 <- paste0(t1, " & \\multicolumn{2}{c}{No target FE} & \\multicolumn{2}{c}{Target FE} \\\\\n\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\n & Bilateral & Coalition & Bilateral & Coalition \\\\\n\\midrule\n")
t1 <- paste0(t1, "\\multicolumn{5}{l}{\\emph{Utility change, gap = initiator minus target}} \\\\\n")
for (k in names(CONTR_LAB)) { ax <- sub("\\..*", "", k); g <- as.numeric(sub("^[MR]\\.", "", k))
  t1 <- paste0(t1, row2(sub("Opponent", "Target", CONTR_LAB[[k]]), lapply(cols1, function(cc) { r <- S1C[type == cc[1] & model == cc[2] & axis == ax & gap == g]; cell(r$est, r$se) }))) }
t1 <- paste0(t1, "\\midrule\n\\multicolumn{5}{l}{\\emph{Surface coefficients}} \\\\\n")
GLAB <- c(dM = "$\\Delta M$", dM2 = "$\\Delta M^2$", dR = "$\\Delta R$", dR2 = "$\\Delta R^2$", dMR = "$\\Delta M \\times \\Delta R$")
for (v in names(GLAB)) t1 <- paste0(t1, row2(GLAB[[v]], lapply(cols1, function(cc) { r <- S1[model == cc[2] & term == paste0(v, "_", cc[1])]; cell(r$est, r$se) })))
ZLAB <- c(contig_land = "Land contiguity", log_capdist = "Log capital distance", tgt_dem = "Target democracy", tgt_dem_x_init_dem = "Joint democracy",
          alliance_def_off = "Defense/offense pact", tgt_major = "Target major power", tgt_major_x_init_major = "Both major powers", dy_mid_tm1 = "MID at $t-1$", log_dy_prior_mids = "Log prior MIDs")
t1 <- paste0(t1, "\\midrule\n\\multicolumn{5}{l}{\\emph{Controls (shared across types)}} \\\\\n")
for (v in names(ZLAB)) t1 <- paste0(t1, row2(ZLAB[[v]], list(cell(S1[model == "no FE" & term == v]$est, S1[model == "no FE" & term == v]$se), c("", ""),
                                                          cell(S1[model == "target FE" & term == v]$est, S1[model == "target FE" & term == v]$se), c("", ""))))
t1 <- paste0(t1, "\\midrule\n", sprintf("Events (bilateral / coalition) & \\multicolumn{4}{c}{%d (%d / %d)} \\\\\nAlternative rows & \\multicolumn{4}{c}{%s} \\\\\n", JF$n1, JF$n1_types[1], JF$n1_types[2], format(JF$n_rows1, big.mark = ",")))
t1 <- paste0(t1, "\\bottomrule\n\\end{tabular}\n\\begin{minipage}{0.95\\linewidth}\\footnotesize\\smallskip Conditional logit over all system members (fatal MIDs, 1-month originator window, main relational measure). Two-way cluster-robust SEs (initiator $\\times$ dispute) in parentheses. $\\Delta M$: log CINC share gap (pooled scale); $\\Delta R$: relational-position gap at $t-1$ (onset-year scale). Utility changes evaluate the quadratic surface at a gap of $\\pm 1$ on one axis with the other at 0. Target FE: one indicator per chosen target, fitted on the alternatives that are chosen at least once (never-chosen targets have fixed effects of $-\\infty$ and drop out). $^{*}p<.10$, $^{**}p<.05$, $^{***}p<.01$.\\end{minipage}\n\\end{table}\n")

## ---- Table 2: stage 2 reference ----------------------------------------------------------------------------------------------------
gen <- function(al, generic) {   # map a generic row to the alternative's column name
  opp <- if (al == "B") "A" else "B"; own <- if (al == "B") "B" else "A"; pre <- if (al == "B") "sB_" else "sA_"
  switch(generic, dM = paste0(pre, "dM"), dM2 = paste0(pre, "dM2"), dR = paste0(pre, "dR"), dR2 = paste0(pre, "dR2"), dMR = paste0(pre, "dMR"),
         land = paste0("land_side_", opp), lcap = paste0("lcap_", opp), dem = paste0("dem_side_", opp), demx = paste0("dem_x_", opp), defoff = paste0("defoff_side_", opp),
         maj = paste0("major_side_", opp), majx = paste0("major_x_", opp), mid = paste0("mid_tm1_side_", opp), prior = paste0("prior_side_", opp),
         tie = paste0("tie_sum_", own), ally = paste0("ally_side_", own), share = paste0("share_", own), nA = "n_side_A", nB = "n_side_B", likely = "likely", int = "(int@mean)")
}
G2 <- c(dM = "$\\Delta M$", dM2 = "$\\Delta M^2$", dR = "$\\Delta R$", dR2 = "$\\Delta R^2$", dMR = "$\\Delta M \\times \\Delta R$",
        land = "Land contiguity", lcap = "Log distance (closest member)", dem = "Opponent democracy (share)", demx = "Democracy $\\times$ opponent democracy",
        defoff = "Defense/offense pact", maj = "Major power on opposing side", majx = "Major $\\times$ opposing major", mid = "MID at $t-1$", prior = "Log prior MIDs (distinct)",
        tie = "Tie sum to joined side", ally = "ATOP pact with joined side", share = "Joined side's CINC share", nA = "Initiating-side size", nB = "Target-side size", likely = "Likely tier")
blk <- list("Opponent surface" = c("dM", "dM2", "dR", "dR2", "dMR"), "Stage-1 controls toward the opponent" = c("land", "lcap", "dem", "demx", "defoff", "maj", "majx", "mid", "prior"),
            "Joined side and dispute" = c("tie", "ally", "share", "nA", "nB", "likely"))
alts <- c("A", "B", "O")
t2 <- "\\begin{table}[!htbp]\\centering\\small\n\\caption{Stage 2: third-party entry and side choice (four-category logit, stay out = reference)}\\label{tab:stage2}\n\\begin{tabular}{lccc}\n\\toprule\n"
t2 <- paste0(t2, " & Join initiating side & Join target side & Originator entry \\\\\n & (fights target set) & (fights initiators) & (fights target set) \\\\\n\\midrule\n\\multicolumn{4}{l}{\\emph{Utility change, gap = joiner minus the side it fights}} \\\\\n")
for (k in names(CONTR_LAB)) { ax <- sub("\\..*", "", k); g <- as.numeric(sub("^[MR]\\.", "", k))
  t2 <- paste0(t2, row2(CONTR_LAB[[k]], lapply(alts, function(al) { r <- REFC[alt == al & axis == ax & gap == g]; cell(r$est, r$se) }))) }
for (bn in names(blk)) { t2 <- paste0(t2, "\\midrule\n\\multicolumn{4}{l}{\\emph{", bn, "}} \\\\\n")
  for (gname in blk[[bn]]) t2 <- paste0(t2, row2(G2[[gname]], lapply(alts, function(al) { r <- REF[term == paste0(al, ":", gen(al, gname))]; if (nrow(r)) cell(r$est, r$se) else c("--", "") }))) }
nref <- 63143L
t2 <- paste0(t2, "\\bottomrule\n\\end{tabular}\n\\begin{minipage}{0.95\\linewidth}\\footnotesize\\smallskip Reference cohort: system members at onset that are not original participants, complete horizons, fatal MIDs (stage-1 set), 1-month window; 63,143 candidate-dispute rows (initiating side 121, target side 184, originator 71). Two-way cluster-robust SEs (candidate state $\\times$ dispute). Joining the initiating side is unavailable when the horizon ends within the window; originator entry is available only to likely-tier candidates (so the tier and democracy $\\times$ opponent democracy do not enter that equation; --). Conflict histories from COW-MID. $^{*}p<.10$, $^{**}p<.05$, $^{***}p<.01$.\\end{minipage}\n\\end{table}\n")

## ---- Table 3: two-step FE joint model -----------------------------------------------------------------------------------------------
BT <- JF$BOOT; bq <- function(qq) { r <- BT[BT$q == qq]; if (nrow(r)) cell_ci(r$est, r$lo, r$hi) else c("", "") }
t3 <- "\\begin{table}[!htbp]\\centering\\small\n\\caption{Two-step fixed-effects joint model: shared target effect}\\label{tab:jointfe}\n\\begin{tabular}{lccc}\n\\toprule\n"
t3 <- paste0(t3, "\\multicolumn{4}{l}{\\emph{Stage 1 with target fixed effects: utility change}} \\\\\n & Bilateral & Coalition & \\\\\n")
for (k in names(CONTR_LAB)) t3 <- paste0(t3, row2(sub("Opponent", "Target", CONTR_LAB[[k]]), list(bq(paste0("S1.bil.", k)), bq(paste0("S1.coal.", k)), c("", ""))))
t3 <- paste0(t3, "\\midrule\n\\multicolumn{4}{l}{\\emph{Stage 2 with the target set's mean fixed effect}} \\\\\n & Join initiating side & Join target side & Originator entry \\\\\n")
for (k in names(CONTR_LAB)) t3 <- paste0(t3, row2(CONTR_LAB[[k]], lapply(alts, function(al) bq(paste0("S2.", al, ".", k)))))
t3 <- paste0(t3, row2("$\\theta$: target-set fixed effect $\\bar b$", list(bq("A:bbar"), bq("B:bbar"), bq("O:bbar"))))
t3 <- paste0(t3, row2("Tie sum to joined side", list(bq("A:tie_sum_A"), bq("B:tie_sum_B"), bq("O:tie_sum_A"))))
t3 <- paste0(t3, row2("Joined side's CINC share", list(c("", ""), bq("B:share_B"), c("", ""))))
t3 <- paste0(t3, sprintf("\\midrule\nBootstrap draws (valid / requested) & \\multicolumn{3}{c}{%d / %d} \\\\\nStage-2 rows / LR($\\theta$, 3 df) & \\multicolumn{3}{c}{%s / %.2f} \\\\\n",
                         JF$B - JF$fails, JF$B, format(JF$F_t$n, big.mark = ","), 2 * (JF$F_t$ll - JF$F_nt$ll)))
t3 <- paste0(t3, "\\bottomrule\n\\end{tabular}\n\\begin{minipage}{0.95\\linewidth}\\footnotesize\\smallskip Step 1: stage-1 conditional logit with target fixed effects. Step 2: the stage-2 reference model plus $\\bar b$, the mean fixed effect of the original target set (targets chosen at least twice, centred, clipped to $[-5,5]$); disputes without a usable $\\bar b$ are dropped. Percentile 95\\% intervals from a dispute-cluster bootstrap that refits both steps on each draw (seed 123); $^{\\dagger}$ interval excludes 0. Targets never chosen in a draw leave that draw's fixed-effect model (their fixed effects diverge), as in the point fit.\\end{minipage}\n\\end{table}\n")

writeLines(c(t1, t2, t3), file.path(OUT, "tables.tex"))
writeLines(c("\\documentclass[11pt]{article}", "\\usepackage[margin=0.6in]{geometry}", "\\usepackage{booktabs}", "\\renewcommand{\\arraystretch}{0.85}", "\\begin{document}", t1, "\\clearpage", t2, "\\clearpage", t3, "\\end{document}"),
           file.path(OUT, "tables_standalone.tex"))
cat("wrote", file.path(OUT, "tables.tex"), "\n")
