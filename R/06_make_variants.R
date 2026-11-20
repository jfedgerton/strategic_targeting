## ============================================================================
## 24_rel_variants.R  --  dump variants for the punishment-vs-selection checks (PI request 2026-11-18)
## Reads a long frame from 10_clogit_roles.R --dump_long=TRUE and writes a modified copy that 22_rel_sweep.R runs unchanged.
##   --variant=aggr    adds aggr_<role> = log(1 + number of distinct MIDs the ALTERNATIVE TARGET originated on side A in the five years before the dispute start)
##                     (GML participant table, peacesciencer gml_part: sidea == 1 & orig == 1) to Z. Rival story: isolation proxies a
##                     target's own aggressive behavior. The current dispute is excluded explicitly (by dispnum), since a dispute that
##                     started before its escalation year t would otherwise fall inside t-5..t-1.
##   --variant=expdef  adds each ALTERNATIVE target's pre-dispute defensibility (ATOP defense partners in the year before the dispute start
##                     with --anchor=onset (default), or before the event with --anchor=entry; initiator excluded; layers
##                     panel from 20_build_layers.R) as type-specific terms def_<role>. Ex-ante replacement for the realized-defender
##                     coding of s22: do choosers avoid defensible targets, and does the relational effect survive it?
##   --variant=dropside drops the joiner strata of the side-reviewed disputes (--level=hard: 200, 1379; --level=reviewed: all 13).
##   --variant=drop    drops the strata of the listed disputes (--disputes=4137,4283) for leave-disputes-out checks.
## Usage: Rscript v2/R/24_rel_variants.R <long_dump.rds> <out.rds> --variant=aggr|expdef|drop|dropside [--level=hard|reviewed] [--layers=<layers_dy.rds>] [--anchor=onset|entry] [--disputes=...]
## ============================================================================
suppressMessages({ library(data.table); library(peacesciencer) })
set.seed(123)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
IN <- args[1]; OUT <- args[2]
opt <- function(name, default) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
VAR <- opt("variant", ""); stopifnot(VAR %in% c("aggr", "expdef", "drop", "dropside"))
logf <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
d <- readRDS(IN); D <- as.data.table(d$D)
stopifnot(all(c("ev_id", "chosen", "type", "c1", "c2", "yr", "disp") %in% names(D)))
n0 <- nrow(D); s0 <- uniqueN(D$ev_id)

## ---- dispute start (GPT full audit 2026-11-19 P1-3): the AUTHORITATIVE dispute-level start year, not the earliest retained event ----------
## GML dispute table (peacesciencer gml_mid_disps styear) for the event's primary dispute; events.rds month bounds only where GML has no row
## (logged). Collapsed events (several source disputes, D$dispnums) are anchored on the primary dispute and counted where the source starts
## differ. A start after the event year is a data conflict: listed, then set to the event year.
disp_start <- function(D) {
  st <- as.data.table(peacesciencer::gml_mid_disps)[, .(disp = as.integer(dispnum), st = as.integer(styear))]; stopifnot(!anyDuplicated(st$disp))
  ev <- as.data.table(readRDS(file.path(d$data_dir, "events.rds")))[, .(disp = as.integer(dispnum), st_ev = as.integer(floor(ev_start_lo / 12)))][, .(st_ev = min(st_ev)), by = disp]
  K <- unique(D[, .(ev_id, disp = as.integer(disp), yr, dispnums = if ("dispnums" %in% names(D)) dispnums else as.character(disp))])
  K <- merge(K, st, by = "disp", all.x = TRUE); K <- merge(K, ev, by = "disp", all.x = TRUE)
  logf("dispute start: GML styear for %d of %d events; events.rds fallback %d; none (event year) %d; GML and the build's earliest event disagree on %d events (GML used)",
       sum(!is.na(K$st)), nrow(K), sum(is.na(K$st) & !is.na(K$st_ev)), sum(is.na(K$st) & is.na(K$st_ev)), sum(!is.na(K$st) & !is.na(K$st_ev) & K$st != K$st_ev))
  K[, st_fallback := is.na(st)]   # delta review: start not from GML (events.rds bound or event year), or a conflicting GML start (below)
  K[is.na(st), st := st_ev]; K[is.na(st), st := yr]
  multi <- K[grepl(";", dispnums)]
  if (nrow(multi)) {
    ms <- vapply(strsplit(multi$dispnums, ";"), function(x) { v <- st$st[match(as.integer(x), st$disp)]; v <- v[!is.na(v)]; if (length(v)) as.integer(diff(range(v))) else 0L }, 0L)
    logf("dispute start: collapsed events %d; whose source disputes start in different years %d (anchored on the primary dispute)", nrow(multi), sum(ms > 0L))
  }
  late <- K[st > yr]
  if (nrow(late)) logf("dispute start: %d events whose dispute start is after the event year (set to the event year): %s", nrow(late),
                       paste(sprintf("disp %d st %d yr %d", late$disp, late$st, late$yr), collapse = "; "))
  K[st > yr, `:=`(st = yr, st_fallback = TRUE)]
  K[, .(ev_id, st, dispnums, st_fallback)]
}

if (VAR == "aggr") {
  ## review P1-1 (2026-11-18): the window is anchored on the DISPUTE START year (st-5..st-1), not the event year, so it never overlaps the
  ## ongoing war; dispute start from disp_start() (GML dispute table, full audit P1-3); GML ends in 2010 (later windows are censored: flagged).
  ## Entered type-specific (aggr_<role>), since defenders' targets are by construction recent originators. Every source dispute of the event
  ## is excluded from the target's history (full audit P2-13).
  gp <- as.data.table(peacesciencer::gml_part)[sidea == 1 & orig == 1, .(dispnum = as.integer(dispnum), ccode = as.integer(ccode), styear = as.integer(styear))]
  gp <- unique(gp)
  S <- disp_start(D)
  K <- unique(merge(D[, .(ev_id, c2, disp = as.integer(disp))], S, by = "ev_id"))
  K[, tgt_aggr_n := mapply(function(c, s0, src) gp[ccode == c & styear >= s0 - 5L & styear <= s0 - 1L & !(dispnum %in% as.integer(strsplit(src, ";")[[1]])), uniqueN(dispnum)], c2, st, dispnums)]
  K[, aggr_censored := st - 1L > 2010L]
  logf("aggr: strata whose window passes the end of GML (2010; flagged aggr_censored, counts are a lower bound) %d; strata with a fallback dispute start %d",
       uniqueN(K[aggr_censored == TRUE, ev_id]), uniqueN(K[st_fallback == TRUE, ev_id]))
  ## delta review P2-11: --window=complete (default) keeps only strata with a complete, GML-reconciled history window; --window=complete_base
  ## keeps the same strata WITHOUT the control (the comparator fit on exactly the same rows); --window=censored keeps every stratum (truncated
  ## counts after 2010 = lower bound; appendix sensitivity)
  WIN <- opt("window", "complete"); stopifnot(WIN %in% c("complete", "complete_base", "censored"))
  if (WIN != "censored") { drop_ev <- unique(K[aggr_censored | st_fallback, ev_id]); K <- K[!ev_id %in% drop_ev]; D <- D[!ev_id %in% drop_ev]; n0 <- nrow(D)
    logf("aggr (%s): dropped %d strata with a censored or fallback window; %d strata kept", WIN, length(drop_ev), uniqueN(D$ev_id)) }
  D <- merge(D, K[, .(ev_id, c2, tgt_aggr_n, aggr_censored)], by = c("ev_id", "c2"), all.x = TRUE, sort = FALSE)
  stopifnot(!anyNA(D$tgt_aggr_n), nrow(D) == n0)
  for (g in seq_along(d$roles)) logf("aggr: type %s, mean originated MIDs chosen %.2f vs alternatives %.2f", d$roles[g],
                                     mean(D[type == g - 1L & chosen == 1, tgt_aggr_n]), mean(D[type == g - 1L & chosen == 0, tgt_aggr_n]))
  if (WIN != "complete_base") { for (g in seq_along(d$roles)) D[, (paste0("aggr_", d$roles[g])) := log1p(tgt_aggr_n) * (type == g - 1L)]
    d$Z_cols <- c(d$Z_cols, paste0("aggr_", d$roles)) }
  d$types <- paste0(d$types, c(complete = "_aggr", complete_base = "_aggrbase", censored = "_aggrcens")[[WIN]])
}

if (VAR == "expdef") {
  ## Pre-dispute defensibility of EVERY alternative target (never of the chosen one only: splitting events by the chosen target's
  ## attribute would define the type from the outcome). n_def = ATOP defense partners of the alternative at t-1, initiator excluded.
  ## Enters as type-specific terms def_<role> = log(1 + n_def) x 1[type == role], so each type has its own response to defensibility.
  ## GPT full audit 2026-11-19 P1-5: --anchor=onset (default) counts partners in the year before the DISPUTE START (disp_start()), so a late
  ## joiner's alternatives are measured before the conflict; --anchor=entry keeps the year before the EVENT (pre-entry defensibility, which the
  ## dispute itself can change: a possible mediator, reported as a separate variant, tag _expdefentry).
  ANCHOR <- opt("anchor", "onset"); stopifnot(ANCHOR %in% c("onset", "entry"))
  L <- as.data.table(readRDS(opt("layers", "")))[def == 1L, .(ccode1, ccode2, year)]
  L2 <- rbind(L[, .(a = ccode1, b = ccode2, year)], L[, .(a = ccode2, b = ccode1, year)])
  cnt <- L2[, .(n_all = uniqueN(b)), by = .(a, year)]
  if (ANCHOR == "onset") {
    S <- disp_start(D); D <- merge(D, S[, .(ev_id, st)], by = "ev_id", sort = FALSE); D[, yr_m1 := st - 1L]
    logf("expdef: anchor = dispute start - 1; events where it differs from event year - 1: %d of %d", uniqueN(D[st != yr, ev_id]), uniqueN(D$ev_id))
    D[, st := NULL]
  } else D[, yr_m1 := yr - 1L]
  D <- merge(D, cnt, by.x = c("c2", "yr_m1"), by.y = c("a", "year"), all.x = TRUE, sort = FALSE)
  D[is.na(n_all), n_all := 0L]
  ## exclude the initiator if it is itself a defense partner of the alternative
  pair <- unique(L2[, .(a, b, year, ally_init = 1L)])
  D <- merge(D, pair, by.x = c("c2", "c1", "yr_m1"), by.y = c("a", "b", "year"), all.x = TRUE, sort = FALSE)
  D[is.na(ally_init), ally_init := 0L]
  D[, n_def := n_all - ally_init]; stopifnot(all(D$n_def >= 0), nrow(D) == n0)
  for (g in seq_along(d$roles)) D[, (paste0("def_", d$roles[g])) := log1p(n_def) * (type == g - 1L)]
  logf("expdef: alternatives with >= 1 defense partner at t-1 %.3f; chosen targets %.3f; mean n_def chosen %.2f vs alternatives %.2f",
       mean(D$n_def > 0), mean(D[chosen == 1, n_def] > 0), mean(D[chosen == 1, n_def]), mean(D[chosen == 0, n_def]))
  D[, c("yr_m1", "n_all", "ally_init") := NULL]
  d$Z_cols <- c(d$Z_cols, paste0("def_", d$roles))
  d$types <- paste0(d$types, if (ANCHOR == "onset") "_expdef" else "_expdefentry")
}

if (VAR == "drop") {
  dd <- strsplit(opt("disputes", ""), ",")[[1]]; stopifnot(length(dd) > 0)
  ## full audit P2-13: an event is dropped when ANY of its source disputes is listed (not only its primary dispute)
  hit <- if ("dispnums" %in% names(D)) vapply(strsplit(D$dispnums, ";"), function(x) any(x %in% dd), TRUE) else as.character(D$disp) %in% dd
  n_sec <- if ("dispnums" %in% names(D)) uniqueN(D$ev_id[hit & !(as.character(D$disp) %in% dd)]) else 0L
  D <- D[!hit]
  logf("drop: removed disputes %s (%d strata; %d of them only through a secondary source dispute)", paste(dd, collapse = ","), s0 - uniqueN(D$ev_id), n_sec)
  d$types <- paste0(d$types, "_drop", paste(dd, collapse = "-"))   # review P3-11: tag names the disputes
}

if (VAR == "dropside") {
  ## PI 2026-11-19 sensitivity: drop the joiner events whose initiating side needed case review (side_case from 10_clogit_roles.R).
  ## --level=hard drops disputes 200 and 1379 (joiner events); --level=reviewed drops all 13 reviewed disputes (joiner events).
  ## --level=hard|reviewed: joiner strata only (incl. joiners promoted by the window, loop review P2-12);
  ## --level=hard_disp|reviewed_disp: EVERY event any of whose source disputes is one of the listed disputes
  LEV <- opt("level", "hard"); stopifnot(LEV %in% c("hard", "reviewed", "hard_disp", "reviewed_disp"), "side_case" %in% names(D))
  SIDE_HARD <- c("200", "1379"); SIDE_REVIEWED <- c("112", "200", "362", "611", "1035", "1046", "1293", "1363", "1379", "2150", "3098", "3131", "3826")
  drop_ev <- if (LEV %in% c("hard", "reviewed")) unique(D[chosen == 1L & side_case %in% (if (LEV == "hard") "hard" else c("hard", "reviewed")), ev_id]) else {
    dl <- if (LEV == "hard_disp") SIDE_HARD else SIDE_REVIEWED
    unique(D[chosen == 1L & vapply(strsplit(dispnums, ";"), function(x) any(x %in% dl), TRUE), ev_id]) }
  D <- D[!(ev_id %in% drop_ev)]
  logf("dropside (%s): removed %d joiner strata; types left %s", LEV, length(drop_ev), paste(tabulate(D[chosen == 1L, type] + 1L, length(d$roles)), collapse = "/"))
  d$types <- paste0(d$types, "_dropside", LEV)
}

stopifnot(all(D[, sum(chosen), by = ev_id]$V1 == 1L))
setorder(D, ev_id, -chosen)
d$D <- as.data.frame(D); d$variant <- VAR
saveRDS(d, paste0(OUT, ".tmp")); stopifnot(file.rename(paste0(OUT, ".tmp"), OUT))   # atomic
logf("wrote %s: %d rows, %d strata (was %d / %d); types %s; Z %s", OUT, nrow(D), uniqueN(D$ev_id), n0, s0, paste(d$roles, collapse = "/"), paste(d$Z_cols, collapse = ","))
