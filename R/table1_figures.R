## R/table1_figures.R -- Table 1 (summary statistics of every quantity the main analyses use) and draft figures, for the MAIN measure on one
## window (PI 2026-11-24). Descriptive only: no model is fitted here.
## Table 1 panels (CSV):
##   A  stage 1 (target choice): events by initiator type, disputes, initiators, chosen targets, choice-set size; material (dM) and relational
##      (dR) gaps of chosen vs all available targets by type; the stage-1 controls
##   B  stage 2 (joining): candidate rows and outcomes (A / B / O / out) by sample; recruitment ties, side gaps and controls by outcome
##   C  network (main measure): by era, states, dyad-years, share with a tie, mean tie, and the spread of log PageRank
## Figures (PDF + PNG): F1 relational position of the major powers over time (log of PageRank x number of states, 0 = uniform share); F2 relational gap of chosen targets by initiator
##   type; F3 stage-2 joining shares by the candidate's relational gap to the attacking side.
## Usage: Rscript R/table1_figures.R --positions=FILE --spec=ID --edges_dir=DIR --build_fmid=DIR --dump_fmid=FILE [--build_mid=DIR --dump_mid=FILE] --out=DIR
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
for (f in c("s2.R", "s2_fit.R", "sry.R", "s2_cell.R")) source(file.path("R/lib", f))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
POSF <- opt("positions"); SPEC <- opt("spec"); EDIR <- opt("edges_dir"); OUT <- opt("out"); stopifnot(!is.null(POSF), !is.null(SPEC), !is.null(EDIR), !is.null(OUT), !is.null(opt("build_fmid")), !is.null(opt("dump_fmid")))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); set.seed(123)
pz <- readRDS(POSF); key <- as.data.table(pz$key); stopifnot(SPEC %in% colnames(pz$R)); key[, R := pz$R[, SPEC]]
E <- s2_spec_edges("core", SPEC, EDIR)
qs <- function(x) { x <- x[is.finite(x)]; c(n = length(x), mean = mean(x), sd = sd(x), p10 = unname(quantile(x, .1)), p50 = median(x), p90 = unname(quantile(x, .9))) }
summ <- function(DT, vars, by = NULL) rbindlist(lapply(vars, function(v) DT[, as.list(qs(get(v))), by = by][, variable := v]), use.names = TRUE)

## ---- panel A: stage 1 ----
SM <- list(fmid = list(b = opt("build_fmid"), d = opt("dump_fmid")))
if (!is.null(opt("build_mid"))) SM$mid <- list(b = opt("build_mid"), d = opt("dump_mid"))
A <- list(); A2 <- list(); BB <- list()
for (s in names(SM)) {
  dump <- readRDS(SM[[s]]$d); D <- s1_prepare(dump, key[, .(ccode, year)], key$R); D[, role := dump$roles[type + 1L]]
  ev <- D[, .(n_alt = .N, chosen = sum(chosen)), by = .(ev_id, role, disp, init)]; stopifnot(all(ev$chosen == 1L))
  A[[s]] <- ev[, .(sample = s, events = .N, disputes = uniqueN(disp), initiators = uniqueN(init), choice_set_mean = mean(n_alt), choice_set_min = min(n_alt), choice_set_max = max(n_alt)), by = role]
  A[[paste0(s, "_tg")]] <- D[chosen == 1L, .(sample = s, role = "all", events = NA_integer_, disputes = NA_integer_, initiators = NA_integer_, choice_set_mean = NA_real_, choice_set_min = NA_integer_, choice_set_max = NA_integer_, distinct_targets = uniqueN(c2))]
  A2[[s]] <- rbind(summ(D[chosen == 1L], c("dM", "dR", dump$Z_cols), "role")[, rows := "chosen targets"],
                   summ(D, c("dM", "dR", dump$Z_cols), "role")[, rows := "all available targets"])[, sample := s]
  BB[[s]] <- list(B = s2_base(SM[[s]]$b), D = D, dump = dump)
}
fwrite(rbindlist(A, fill = TRUE), file.path(OUT, "table1_A_stage1_counts.csv")); fwrite(rbindlist(A2), file.path(OUT, "table1_A_stage1_variables.csv"))

## ---- panel B: stage 2 ----
Bt <- list(); Bv <- list()
for (s in names(BB)) { B <- BB[[s]]$B; W_REF <<- B$window; X <- s2_prepare(B, key[, .(ccode, year)], key$R, E)
  Bt[[s]] <- X[, .(rows = .N, disputes = uniqueN(dispnum), candidates = uniqueN(ccode)), by = y][, sample := s]
  v <- intersect(c("tie_sum_A", "tie_sum_B", "dM_kA", "dR_kA", "dM_kB", "dR_kB", "likely", "land_side_A", "land_side_B", "lcap_A", "lcap_B", "defoff_side_A", "defoff_side_B", "own_cent", "tgt_cent_mean"), names(X))
  Bv[[s]] <- summ(X, v, "y")[, sample := s] }
fwrite(rbindlist(Bt), file.path(OUT, "table1_B_stage2_counts.csv")); fwrite(rbindlist(Bv), file.path(OUT, "table1_B_stage2_variables.csv"))

## ---- panel C: network ----
era <- function(y) cut(y, c(1815, 1913, 1945, 1991, 2014), labels = c("1816-1913", "1914-1945", "1946-1991", "1992-2014"))
nd <- key[, .(states = .N), by = year]; ed <- E[, .(edges = .N, mean_w = mean(w)), by = year]
C <- merge(nd, ed, by = "year", all.x = TRUE)[is.na(edges), `:=`(edges = 0L, mean_w = 0)]
C[, dyads := states * (states - 1) / 2][, share_tie := edges / dyads]
key[, sdR := sd(R), by = year]
C <- merge(C, unique(key[, .(year, sd_logPR = sdR)]), by = "year")
CE <- C[, .(years = .N, states_mean = mean(states), dyad_years = sum(dyads), share_with_tie = sum(edges) / sum(dyads), mean_tie_if_tie = sum(mean_w * edges) / sum(edges), sd_logPR_mean = mean(sd_logPR)), by = .(era = era(year))]
fwrite(CE, file.path(OUT, "table1_C_network.csv")); fwrite(C, file.path(OUT, "table1_C_network_by_year.csv"))

## ---- figures ----
th <- theme_bw(base_size = 11) + theme(legend.position = "bottom", strip.placement = "outside", panel.grid.minor = element_blank())
sv <- function(p, nm, w = 7, h = 4.5) { ggsave(file.path(OUT, paste0(nm, ".pdf")), p, width = w, height = h); ggsave(file.path(OUT, paste0(nm, ".png")), p, width = w, height = h, dpi = 200) }
mp <- data.table(ccode = c(2, 200, 220, 255, 365, 710, 740), state = c("United States", "United Kingdom", "France", "Germany", "Russia", "China", "Japan"))
F1 <- merge(key, mp, by = "ccode")
F1[ccode == 255, state := "Germany"]; F1 <- rbind(F1, merge(key[ccode == 260], data.table(ccode = 260, state = "Germany"), by = "ccode"))
F1 <- merge(F1, key[, .(N = .N), by = year], by = "year")[, Rn := R + log(N)]   # log(N x PageRank): 0 = the uniform share, comparable across system sizes
F1 <- merge(CJ(state = unique(F1$state), year = min(key$year):max(key$year)), F1[, .(state, year, Rn)], by = c("state", "year"), all.x = TRUE)   # gaps (e.g. Germany 1946-54) break the line
p1 <- ggplot(F1, aes(year, Rn, colour = state)) + geom_hline(yintercept = 0, linetype = 3) + geom_line(linewidth = .5, na.rm = TRUE) +
  labs(x = NULL, y = "Relational position\n(log of PageRank x number of states)", colour = NULL) + th
sv(p1, "F1_major_power_positions")
D1 <- BB[["fmid"]]$D[chosen == 1L]; D1[, role_lab := factor(role, levels = c("bil", "coal"), labels = c("Bilateral initiators", "Coalition originators"))]
p2 <- ggplot(D1, aes(dR)) + geom_histogram(bins = 30, fill = "grey55", colour = "white") + geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~ role_lab, ncol = 2, scales = "free_y") + labs(x = "Relational gap to the chosen target\n(initiator minus target, log PageRank)", y = "Events") + th
sv(p2, "F2_chosen_target_relational_gap")
B1 <- BB[["fmid"]]$B; W_REF <- B1$window; X1 <- s2_prepare(B1, key[, .(ccode, year)], key$R, E)
X1 <- X1[is.finite(dR_kA)]; X1[, bin := cut(dR_kA, unique(quantile(dR_kA, seq(0, 1, .2))), include.lowest = TRUE)]
F3 <- X1[, .(n = .N, A = mean(y == "A"), B = mean(y == "B"), O = mean(y == "orig")), by = bin]
F3 <- melt(F3, id.vars = c("bin", "n"), variable.name = "outcome", value.name = "share")
F3[, outcome := factor(outcome, levels = c("A", "B", "O"), labels = c("Joins attacker (A)", "Joins target (B)", "Joins as originator (O)"))]
p3 <- ggplot(F3, aes(bin, 1000 * share)) + geom_col(fill = "grey55") + facet_wrap(~ outcome, ncol = 3, scales = "free_y") +
  scale_x_discrete(labels = function(x) gsub(",", ",\n", x)) + labs(x = "Candidate's relational gap to the attacking side (quintiles)", y = "Joins per 1,000 candidates") + th
sv(p3, "F3_joining_by_relational_gap", w = 9)
cat("TABLE1_FIGURES_DONE", OUT, "\n")
