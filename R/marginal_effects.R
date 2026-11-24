## R/marginal_effects.R -- marginal-effect plots with 95% bootstrap bands for the MAIN measure, one window x sample (PI 2026-11-24).
##   ME1 stage 1 (between targets, no target FE): log odds of choosing a target relative to a peer (gaps 0, 0), along the relational gap
##       (material gap 0) and along the material gap (relational gap 0), by initiator type, plus the coalition - bilateral difference
##   ME2 stage 2 alone (on the joint model's disputes, as main-table column 2): average predicted probability of joining A / B / O as one gap of the candidate to the opposing side moves from its
##       5th to 95th percentile (all other variables at their observed values; the surface terms rebuilt), per 1,000 candidates
##   ME3 joint model: average predicted probability of joining A / B / O as the dispute's target trait bbar moves across its 5th-95th percentiles
## Bands: percentile 95% over the two-stage dispute bootstrap with the common draws (seed 123), the same resamples as R/s2_07_boot_cells.R
## and R/main_table.R. Each draw refits the model on the reweighted disputes (stage 1: event weights x dispute weight; stage 2: row weights; Bayesian bootstrap, boot_weights);
## a draw that fails, or whose fitted term set differs from the point fit, is dropped and counted. Draws are saved in blocks of 50
## (OUT/boot/draws_k.rds) and a rerun resumes. Colours: Dark2 (orange vs purple for the main contrasts).
## Usage: Rscript R/marginal_effects.R --build=DIR --dump=FILE --positions=FILE --spec=ID --edges_dir=DIR --n_boot=500 --out=DIR [--cores=2]
suppressPackageStartupMessages({ library(data.table); library(parallel); library(ggplot2) })
for (f in c("s2.R", "s2_fit.R", "sry.R", "s2_cell.R")) source(file.path("R/lib", f))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
BUILD <- opt("build"); DUMPF <- opt("dump"); POSF <- opt("positions"); SPEC <- opt("spec"); EDIR <- opt("edges_dir"); OUT <- opt("out")
NB <- as.integer(opt("n_boot", "500")); CORES <- as.integer(opt("cores", "2"))
stopifnot(!is.null(BUILD), !is.null(DUMPF), !is.null(POSF), !is.null(SPEC), !is.null(EDIR), !is.null(OUT), NB >= 10)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
B <- s2_base(BUILD); W_REF <- B$window; dump <- readRDS(DUMPF)
## dump / build pairing, as R/s2_07_boot_cells.R (audit 2026-11-24)
cm <- if (is.null(dump$coorig_months)) 0L else as.integer(dump$coorig_months); smp <- c(stage1 = "fmid", gml = "fmid", mid = "mid")[[B$meta$fatal_def]]
if (dump$types != "init" || cm != W_REF || dump$sample != smp) stop("dump / build mismatch (window or sample)")
pz <- readRDS(POSF); key <- as.data.table(pz$key); stopifnot(SPEC %in% colnames(pz$R)); R <- pz$R[, SPEC]
E <- s2_spec_edges("core", SPEC, EDIR); D <- s1_prepare(dump, key, R); X <- s2_prepare(B, key, R, E)
## keep only the columns the fits and curves use (audit 2026-11-24: a full 99-column copy per grid point used ~1.4 GB per process)
KEEP <- unique(c("dispnum", "ccode", "y", "H", "likely", "dM_kA", "dR_kA", "dM_kB", "dR_kB", unlist(S2_REF_DESIGN))); stopifnot(all(KEEP %in% names(X))); X <- X[, ..KEEP]
ROLES <- dump$roles
## ---- grids ----
q <- function(x, p) unname(quantile(x[is.finite(x)], p))
gR <- seq(q(D$dR, .05), q(D$dR, .95), length.out = 41); gM <- seq(q(D$dM, .05), q(D$dM, .95), length.out = 41)
s1_curve <- function(b) rbindlist(lapply(ROLES, function(r) { k <- function(t) b[[paste0(t, "_", r)]]
  rbind(data.table(role = r, axis = "R", x = gR, y = k("dR") * gR + k("dR2") * gR^2),
        data.table(role = r, axis = "M", x = gM, y = k("dM") * gM + k("dM2") * gM^2)) }))
GAPS <- list(A = c(M = "dM_kB", R = "dR_kB"), B = c(M = "dM_kA", R = "dR_kA"), O = c(M = "dM_kB", R = "dR_kB"))   # gaps to the opposing side
gg <- lapply(c("dM_kA", "dR_kA", "dM_kB", "dR_kB"), function(g) seq(q(X[[g]], .05), q(X[[g]], .95), length.out = 15)); names(gg) <- c("dM_kA", "dR_kA", "dM_kB", "dR_kB")
s2_curve <- function(f, Xd, b = f$b) rbindlist(lapply(names(GAPS), function(al) rbindlist(lapply(names(GAPS[[al]]), function(ax) { g <- GAPS[[al]][[ax]]
  data.table(alt = al, axis = ax, x = gg[[g]], y = vapply(gg[[g]], function(v) { X1 <- copy(Xd); set(X1, j = g, value = v); add_surface(X1); 1000 * mean(mnl_predict(f, X1, b)[, al]) }, 0)) }))))
## ---- the target trait for the joint model (as joint_fe_cell / joint_fe_draw) ----
bbar_of <- function(Dd, ref_order, REP) { s1 <- s1_fit_fe(Dd, dump, ref_order); bb <- coef(s1$m); be <- bb[grepl("^tfe", names(bb))]; names(be) <- sub("^tfe", "", names(be))
  refb <- setdiff(intersect(ref_order, as.character(unique(Dd[chosen == 1L, c2]))), names(be))[1]; be <- c(setNames(0, refb), be)
  b <- be[names(be) %in% REP & is.finite(be)]; b <- pmin(pmax(b - mean(b), -5), 5)
  B$TS[, .(dispnum, b = b[as.character(ccode)])][, .(bbar = mean(b, na.rm = TRUE)), by = dispnum][is.finite(bbar)] }
ch <- D[chosen == 1L]; NCH <- sort(table(ch$c2), decreasing = TRUE); REP <- names(NCH)[NCH >= 2]; ref_order <- names(NCH)
DJ <- lapply(S2_REF_DESIGN, function(v) c(v, "bbar"))
th_curve <- function(f, Xj, gb, b = f$b) rbindlist(lapply(c("A", "B", "O"), function(al) data.table(alt = al, x = gb,
  y = vapply(gb, function(v) { X1 <- copy(Xj); X1[, bbar := v]; 1000 * mean(mnl_predict(f, X1, b)[, al]) }, 0))))
## ---- point estimates ----
mN <- s1_fit_nofe(D, dump); P1 <- s1_curve(coef(mN))
bb0 <- bbar_of(D, ref_order, REP); XJ <- merge(X, bb0, by = "dispnum"); gB <- seq(q(XJ$bbar, .05), q(XJ$bbar, .95), length.out = 15)
XS <- X[dispnum %in% bb0$dispnum]   # stage 2 alone on the joint model's disputes (main-table column 2)
F0 <- fit_mnl(XS, S2_REF_DESIGN, "stage2 alone"); P2 <- s2_curve(F0, XS)
FJ <- fit_mnl(XJ, DJ, "joint"); P3 <- th_curve(FJ, XJ, gB)
## ---- bootstrap (common draws) ----
disp_all <- sort(unique(c(as.data.table(dump$D)$disp, B$X$dispnum)))
MM <- boot_weights(disp_all, NB)   # Bayesian bootstrap, the same draws as s2_07 (lib/s2_cell.R)
one <- function(m) tryCatch({
  Db <- copy(D); Db[, w := w * m[as.character(disp)]]; Db <- Db[w > 0]
  c1 <- s1_curve(coef(s1_fit_nofe(Db, dump)))
  f0 <- fit_mnl(XS, S2_REF_DESIGN, "boot", vc = FALSE, w = as.numeric(m[as.character(XS$dispnum)]))
  if (!identical(names(f0$b), names(F0$b))) stop("stage-2 term set differs from the point fit")   # never a curve from a smaller model
  c2 <- s2_curve(f0, XS)
  bbb <- bbar_of(Db, ref_order, REP); Xjb <- merge(X, bbb, by = "dispnum")
  fj <- fit_mnl(Xjb, DJ, "boot", vc = FALSE, w = as.numeric(m[as.character(Xjb$dispnum)]))
  if (!identical(names(fj$b), names(FJ$b))) stop("joint term set differs from the point fit")
  c3 <- th_curve(fj, XJ, gB)
  list(c1 = c1$y, c2 = c2$y, c3 = c3$y) }, error = function(e) paste("failed:", conditionMessage(e)))
## blocks of 50 saved as they finish; a rerun reuses blocks whose draw indices and grids match (audit 2026-11-24)
BD <- file.path(OUT, "boot"); dir.create(BD, showWarnings = FALSE); BLK <- 50L; sig <- list(NB = NB, spec = SPEC, gR = gR, gM = gM, gg = gg, gB = gB, md5 = unname(tools::md5sum(c("R/marginal_effects.R", file.path("R/lib", c("s2.R", "s2_fit.R", "s2_cell.R")),
             file.path(BUILD, "s2_meta.rds"), DUMPF, POSF))))   # contents, not paths (re-audit P2)
DR <- list()
for (k in seq_len(ceiling(NB / BLK))) { ix <- ((k - 1L) * BLK + 1L):min(k * BLK, NB); bf <- file.path(BD, sprintf("draws_%03d.rds", k))
  if (file.exists(bf)) { o <- readRDS(bf); if (identical(o$ix, ix) && identical(o$sig, sig)) { DR[ix] <- o$draws; next } else stop("stale block ", bf, " (code, inputs or grids changed): delete ", BD, " or use a new --out") }
  dk <- mclapply(MM[ix], one, mc.cores = CORES, mc.preschedule = FALSE)
  saveRDS(list(ix = ix, sig = sig, draws = dk), paste0(bf, ".tmp")); file.rename(paste0(bf, ".tmp"), bf); DR[ix] <- dk
  cat(sprintf("block %d done (%d draws)\n", k, max(ix))) }
ok <- vapply(DR, is.list, TRUE); cat(sprintf("bootstrap: %d of %d draws usable\n", sum(ok), NB))
if (sum(ok) < 0.8 * NB) stop("fewer than 80% usable draws: ", paste(head(unique(unlist(DR[!ok])), 3), collapse = " | "))
band <- function(P, k) { M <- do.call(rbind, lapply(DR[ok], `[[`, k)); P[, `:=`(lo = apply(M, 2, quantile, .025), hi = apply(M, 2, quantile, .975), draws = sum(ok))]; P }
P1 <- band(P1, "c1"); P2 <- band(P2, "c2"); P3 <- band(P3, "c3")
## coalition - bilateral difference in stage 1 (same draws)
M1 <- do.call(rbind, lapply(DR[ok], `[[`, "c1")); ic <- which(P1$role == ROLES[2]); ib <- which(P1$role == ROLES[1])
PD <- data.table(axis = P1$axis[ic], x = P1$x[ic], y = P1$y[ic] - P1$y[ib], lo = apply(M1[, ic, drop = FALSE] - M1[, ib, drop = FALSE], 2, quantile, .025),
                 hi = apply(M1[, ic, drop = FALSE] - M1[, ib, drop = FALSE], 2, quantile, .975))
fwrite(P1, file.path(OUT, "me1_stage1.csv")); fwrite(PD, file.path(OUT, "me1_stage1_coal_minus_bil.csv")); fwrite(P2, file.path(OUT, "me2_stage2.csv")); fwrite(P3, file.path(OUT, "me3_joint_theta.csv"))
## ---- plots ----
D2 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
th <- theme_bw(base_size = 11) + theme(legend.position = "bottom", strip.text = element_text(size = 13), panel.grid.minor = element_blank())
sv <- function(p, nm, w = 9, h = 4.8) { ggsave(file.path(OUT, paste0(nm, ".pdf")), p, width = w, height = h); ggsave(file.path(OUT, paste0(nm, ".png")), p, width = w, height = h, dpi = 200) }
ax_lab <- c(R = "Relational gap (initiator minus target, scaled)\nmaterial gap held at 0", M = "Material gap (initiator minus target, scaled)\nrelational gap held at 0")
stopifnot(identical(ROLES, c("bil", "coal")))   # plot labels assume this order
P1[, `:=`(axis_lab = factor(ax_lab[axis], levels = ax_lab), type = factor(role, levels = ROLES, labels = c("Bilateral initiators", "Coalition originators")))]
p1 <- ggplot(P1, aes(x, y, colour = type, fill = type)) + geom_hline(yintercept = 0, linetype = 3) + geom_vline(xintercept = 0, linetype = 3) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .2, colour = NA) + geom_line(linewidth = .8) + facet_wrap(~ axis_lab, ncol = 2, scales = "free_x", strip.position = "bottom") +
  scale_colour_manual(values = D2[c(3, 2)], name = NULL) + scale_fill_manual(values = D2[c(3, 2)], name = NULL) +   # purple vs orange (colourblind-safe pair)
  labs(x = NULL, y = "Log odds of choosing the target\nrelative to a peer target") + th + theme(strip.background = element_blank(), strip.placement = "outside")
sv(p1, "ME1_stage1_by_type")
PD[, axis_lab := factor(ax_lab[axis], levels = ax_lab)]
p1d <- ggplot(PD, aes(x, y)) + geom_hline(yintercept = 0, linetype = 3) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = D2[3], alpha = .25) + geom_line(colour = D2[3], linewidth = .8) +
  facet_wrap(~ axis_lab, ncol = 2, scales = "free_x", strip.position = "bottom") + labs(x = NULL, y = "Coalition minus bilateral\n(log odds relative to a peer target)") +
  th + theme(strip.background = element_blank(), strip.placement = "outside")
sv(p1d, "ME1b_stage1_coalition_minus_bilateral")
P2[, `:=`(alt_lab = factor(alt, levels = c("A", "B", "O"), labels = c("Joins attacker (A)", "Joins target (B)", "Joins as originator (O)")),
          axis_lab = factor(axis, levels = c("R", "M"), labels = c("Relational gap to the opposing side", "Material gap to the opposing side")))]
p2 <- ggplot(P2, aes(x, y)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = D2[1], alpha = .25) + geom_line(colour = D2[1], linewidth = .8) + geom_vline(xintercept = 0, linetype = 3) +
  facet_grid(alt_lab ~ axis_lab, scales = "free", switch = "x") + labs(x = "Candidate minus opposing side (5th to 95th percentile)", y = "Predicted joins per 1,000 candidates") +
  labs(caption = "Other gap at observed values. Rows where the option is unavailable count as zero joins.") +
  th + theme(strip.background.x = element_blank(), strip.placement = "outside")
sv(p2, "ME2_stage2_by_gap", w = 8.5, h = 8)
P3[, alt_lab := factor(alt, levels = c("A", "B", "O"), labels = c("Joins attacker (A)", "Joins target (B)", "Joins as originator (O)"))]
p3 <- ggplot(P3, aes(x, y, colour = alt_lab, fill = alt_lab)) + geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .2, colour = NA) + geom_line(linewidth = .8) +
  scale_colour_manual(values = D2[c(2, 3, 1)], name = NULL) + scale_fill_manual(values = D2[c(2, 3, 1)], name = NULL) +   # A orange, B purple, O teal (as F6)
  labs(x = "Target trait of the dispute (mean stage-1 target fixed effect, 5th to 95th percentile)", y = "Predicted joins per 1,000 candidates") + th
sv(p3, "ME3_joint_target_trait", w = 7.5, h = 5)
cat("MARGINAL_EFFECTS_DONE", OUT, "\n")
