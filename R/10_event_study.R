## ============================================================================
## 27_rel_eventstudy.R  --  relational trajectory of chosen targets around the dispute, by initiator type (PI request 2026-11-18, test 3)
## Question: are chooser targets isolated well BEFORE the dispute (selection on standing isolation), do they lose ties in the years just
## before it (pariah dynamics or anticipation), or only after it (punishment / consequence)?
## Quantity, per event e and relative year k (event year t = D$yr; k = -PRE..POST):
##   z(c, y)  = position of state c in year y, standardized within year over all states in the measure's key
##   gap_e(k) = z(chosen target, t+k) - mean over the stratum's non-chosen alternatives of z(alt, t+k)
##   (negative = the chosen target is less central than the states the initiator passed over)
## Reported by type: weighted mean gap(k) (event weight w of the chosen row), the same normalized to k = -1, and every type minus the
## reference type. Balanced panel by default: an event enters only if its chosen target and at least one alternative have positions in
## every year of the window; alternatives missing any year are dropped from that event. Intervals: pigeonhole two-way bootstrap
## (initiators x disputes drawn independently, event weight x product of multiplicities, as 25_rel_boot.R), percentile 95%, seed 123.
## Descriptive only: no controls, no choice model.
## Usage: Rscript v2/R/27_rel_eventstudy.R <long_dump.rds> <positions.rds> <spec_id> <out_prefix> [--pre=10] [--post=5] [--balanced=TRUE] [--B=1000]
## ============================================================================
suppressMessages(library(data.table))
set.seed(123)
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) >= 4)
DUMP <- args[1]; POSF <- args[2]; SPEC <- args[3]; OUTP <- args[4]
opt <- function(name, default) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
PRE <- as.integer(opt("pre", "10")); POST <- as.integer(opt("post", "5")); BAL <- as.logical(opt("balanced", "TRUE")); B <- as.integer(opt("B", "1000"))
## loop review P2-7: --align=entry (default; k = 0 is the event year, i.e. the pair's entry year) or onset (k = 0 is the primary dispute's start
## year, GML dispute table). A joiner entering in 1995 a dispute begun in 1990 has k = -3 = 1992 under entry, which is already inside the dispute.
ALIGN <- opt("align", "entry"); stopifnot(ALIGN %in% c("entry", "onset"))
fl <- sub("^--([A-Za-z_]+)=.*$", "\\1", grep("^--", args, value = TRUE)); if (any(!fl %in% c("pre", "post", "balanced", "B", "align"))) stop("unknown option(s)")
stopifnot(PRE >= 1L, POST >= 0L, !is.na(BAL), B >= 0L)
logf <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
KS <- -PRE:POST

## ---- 1. inputs ----
dump <- readRDS(DUMP); D <- as.data.table(dump$D); ROLES <- dump$roles; G <- length(ROLES)
stopifnot(all(c("ev_id", "chosen", "type", "w", "init", "disp", "c1", "c2", "yr") %in% names(D)), all(D[, sum(chosen), by = ev_id]$V1 == 1L))
P <- readRDS(POSF); j <- match(SPEC, colnames(P$R)); stopifnot(!is.na(j))
pos <- data.table(ccode = as.integer(P$key$ccode), year = as.integer(P$key$year), R = P$R[, j]); stopifnot(!anyDuplicated(pos[, .(ccode, year)]), !anyNA(pos$R))
pos[, z := { s <- sd(R); if (s < 1e-8) NA_real_ else (R - mean(R)) / s }, by = year]
if (anyNA(pos$z)) stop(sprintf("measure %s is constant in year(s) %s", SPEC, paste(unique(pos$year[is.na(pos$z)]), collapse = " ")))
logf("dump %s (%s; %d strata) | measure %s | window k = %d..%d | balanced %s", basename(DUMP), dump$types, uniqueN(D$ev_id), SPEC, -PRE, POST, BAL)

## ---- 2. long panel: every (event row, k) with the alternative's z at t+k ----
if (ALIGN == "onset") {
  st <- as.data.table(peacesciencer::gml_mid_disps)[, .(disp = as.integer(dispnum), st = as.integer(styear))]
  D[, yr_ref := st$st[match(as.integer(disp), st$disp)]]
  n_miss <- uniqueN(D[is.na(yr_ref), ev_id]); D[is.na(yr_ref), yr_ref := yr]; D[yr_ref > yr, yr_ref := yr]
  logf("align = onset: k = 0 is the dispute start year; %d events without a GML start use the event year; median entry lag %.1f years (max %d)",
       n_miss, median(D[chosen == 1L, yr - yr_ref]), max(D[chosen == 1L, yr - yr_ref]))
} else D[, yr_ref := yr]
L <- D[, .(ev_id, chosen, type, w, init, disp, c2, yr = yr_ref)][, .(k = KS), by = .(ev_id, chosen, type, w, init, disp, c2, yr)]
L[, y := yr + k]
L <- merge(L, pos[, .(c2 = ccode, y = year, z)], by = c("c2", "y"), all.x = TRUE, sort = FALSE)
stopifnot(nrow(L) == nrow(D) * length(KS))
if (BAL) {
  ## alternatives (and chosen targets) need a position in every year of the window; events need their chosen target and >= 1 alternative
  full <- L[, .(full = all(!is.na(z))), by = .(ev_id, c2)]
  L <- merge(L, full, by = c("ev_id", "c2"), sort = FALSE)[full == TRUE][, full := NULL]
  ok_ev <- L[k == KS[1], .(has_ch = any(chosen == 1L), n_alt = sum(chosen == 0L)), by = ev_id][has_ch & n_alt >= 1L, ev_id]
} else {
  L <- L[!is.na(z)]
  ok_ev <- L[, .(ok = any(chosen == 1L) && any(chosen == 0L)), by = .(ev_id, k)][, .(ok = all(ok)), by = ev_id][ok == TRUE, ev_id]
}
n_ev0 <- uniqueN(D$ev_id); L <- L[ev_id %in% ok_ev]
logf("events kept %d of %d (%s); by type: %s", length(ok_ev), n_ev0, if (BAL) "balanced window" else "unbalanced",
     paste(sprintf("%s %d", ROLES, tabulate(unique(L[, .(ev_id, type)])$type + 1L, G)), collapse = ", "))
if (length(ok_ev) == 0L) stop("no events survive the window")

## ---- 3. per-event gap(k) ----
E <- L[, .(gap = z[chosen == 1L] - mean(z[chosen == 0L])), by = .(ev_id, k, type, init, disp)]
W <- unique(D[chosen == 1L, .(ev_id, w)]); E <- merge(E, W, by = "ev_id", sort = FALSE)
stopifnot(!anyNA(E$gap), all(E[, .N, by = .(ev_id, k)]$N == 1L))
ev <- unique(E[, .(ev_id, init, disp)]); stopifnot(!anyDuplicated(ev$ev_id))

## summary for a given event-weight multiplier m (named by ev_id): weighted mean gap by type and k
summ <- function(m) {
  e <- E[, .(ev_id, k, type, gap, ww = w * m[as.character(ev_id)])][ww > 0]
  s <- e[, .(gap = sum(ww * gap) / sum(ww)), by = .(type, k)]
  if (uniqueN(s$type) < G) return(NULL)                                  # a type vanished from this draw
  s[, gap_n := gap - gap[k == -1L], by = type]                          # normalized to the year before the event
  ref <- s[type == 0L, .(k, gap_ref = gap, gapn_ref = gap_n)]
  s <- merge(s, ref, by = "k"); s[, `:=`(diff = gap - gap_ref, diff_n = gap_n - gapn_ref)]
  s[, .(type, k, gap, gap_n, diff, diff_n)][order(type, k)]
}
m1 <- setNames(rep(1, nrow(ev)), ev$ev_id); est <- summ(m1); stopifnot(!is.null(est))

## ---- 4. pigeonhole two-way bootstrap ----
out <- copy(est)
if (B > 0L) {
  u_i <- unique(ev$init); u_d <- unique(ev$disp); ii <- match(ev$init, u_i); dd <- match(ev$disp, u_d)
  draws <- vector("list", B)
  for (b in seq_len(B)) {
    mi <- tabulate(sample.int(length(u_i), length(u_i), replace = TRUE), length(u_i))
    md <- tabulate(sample.int(length(u_d), length(u_d), replace = TRUE), length(u_d))
    sb <- summ(setNames(mi[ii] * md[dd], ev$ev_id)); if (!is.null(sb)) draws[[b]] <- sb[, b := b]
  }
  ok <- !vapply(draws, is.null, TRUE); BS <- rbindlist(draws[ok])
  logf("bootstrap: %d of %d draws usable (%d initiators, %d disputes)", sum(ok), B, length(u_i), length(u_d))
  rep_ok <- B >= 200L && mean(ok) >= 0.90
  q <- BS[, .(gap_lo = quantile(gap, .025), gap_hi = quantile(gap, .975), gapn_lo = quantile(gap_n, .025), gapn_hi = quantile(gap_n, .975),
              diff_lo = quantile(diff, .025), diff_hi = quantile(diff, .975), diffn_lo = quantile(diff_n, .025), diffn_hi = quantile(diff_n, .975)), by = .(type, k)]
  if (!rep_ok) { logf("bootstrap below reporting rule (B >= 200 and >= 90%% usable): intervals set to NA"); q[, (setdiff(names(q), c("type", "k"))) := NA_real_] }
  out <- merge(out, q, by = c("type", "k"))
}
out[, role := ROLES[type + 1L]]; out[, `:=`(spec_id = SPEC, scheme = dump$types, balanced = BAL, align = ALIGN, n_events = length(ok_ev), n_events_all = n_ev0)]
fwrite(data.table(ev_id = setdiff(unique(D$ev_id), ok_ev)), paste0(OUTP, "_excluded_events.csv"))   # P2-7: excluded cases published
fwrite(out[order(type, k)], paste0(OUTP, ".csv.tmp")); stopifnot(file.rename(paste0(OUTP, ".csv.tmp"), paste0(OUTP, ".csv")))   # atomic
logf("wrote %s.csv", OUTP)
print(out[k %in% c(-PRE, -5L, -1L, 0L, POST), .(role, k, gap = round(gap, 2), gap_n = round(gap_n, 2), diff_vs_ref = round(diff, 2))], row.names = FALSE)

## ---- 5. figure: gap(k) and gap(k) - gap(-1) by type ----
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  pl <- rbind(out[, .(role, k, panel = "Level: target minus alternatives (within-year z of the position, not a count of ties)", est = gap, lo = if ("gap_lo" %in% names(out)) gap_lo else NA_real_, hi = if ("gap_hi" %in% names(out)) gap_hi else NA_real_)],
              out[, .(role, k, panel = "Change: relative to k = -1", est = gap_n, lo = if ("gapn_lo" %in% names(out)) gapn_lo else NA_real_, hi = if ("gapn_hi" %in% names(out)) gapn_hi else NA_real_)])
  pl[, role := factor(role, levels = ROLES)]
  g <- ggplot(pl, aes(k, est)) + geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") + geom_vline(xintercept = -0.5, linetype = 2, linewidth = 0.3, colour = "grey50") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) + geom_line() + geom_point(size = 1) +
    facet_grid(panel ~ role, scales = "free_y") + labs(x = if (ALIGN == "onset") "Years relative to the dispute's start year" else "Years relative to the pair's entry year (event year)", y = NULL, title = SPEC) + theme_bw(base_size = 10)
  ggsave(paste0(OUTP, ".png"), g, width = 2.6 * G + 1, height = 5.5, dpi = 200)
  logf("wrote %s.png", OUTP)
}
