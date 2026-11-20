## ============================================================================
## 25_rel_boot.R  --  resampling for ONE relational position (2026-11-19 rewrite for the multiverse loop; replaces the 2026-11-18 version)
##   --mode=boot      pigeonhole two-way bootstrap: initiators and disputes drawn independently with replacement; each event's weight times the
##                    product of its initiator's and dispute's multiplicities (as 10_clogit_roles.R section 6). Percentile 95% intervals.
##   --mode=perm_dR   within-event permutation of the relational gap: in every stratum the alternatives' dR values are shuffled (dR2, dM x dR
##                    recomputed), which breaks the link between the chosen target and its relational position and keeps everything else.
##   --mode=perm_type permutation of the type labels across events (each event keeps its alternatives, weights and gaps).
##   Permutation diagnostic p (p_diag) = (1 + #permuted |stat| >= observed |stat|) / (1 + n_ok). The permutations are artificial-assignment
##   DIAGNOSTICS (loop review B), not randomization tests of the scientific null.
## Quantities: M0 (fully pooled) u_R_less / u_R_more / u_M_weaker / u_M_stronger; M1 the same per type and for every pair of types.
## A draw is dropped when its fit errors, warns (separation / non-convergence) or a type vanishes; every dropped draw is listed with its reason.
## Intervals / p-values are reported only when B >= 200 and >= 90% of draws are usable (otherwise NA, and share_pos is NA too).
## Gap: as 22_rel_sweep.R at its defaults (positions t-1, year-t scale for rel_scale = year dumps). Seed 123.
## Usage: Rscript v2/R/25_rel_boot.R <long_dump.rds> <positions.rds> <spec_id> <out.csv> [--mode=boot|perm_dR|perm_type] [--B=2000] [--cores=2]
## ============================================================================
suppressMessages({ library(survival); library(parallel) })
local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); source(file.path(dirname(normalizePath(f)), "lib", "report_gate.R")) })
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) >= 4)
DUMP <- args[1]; POSF <- args[2]; SPEC <- args[3]; OUT <- args[4]
fl <- sub("^--([A-Za-z_]+)=.*$", "\\1", grep("^--", args, value = TRUE))
if (any(!fl %in% c("mode", "B", "cores"))) stop("unknown option(s): ", paste(setdiff(fl, c("mode", "B", "cores")), collapse = ", "))
opt <- function(name, default) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
MODE <- opt("mode", "boot"); stopifnot(MODE %in% c("boot", "perm_dR", "perm_type"))
B <- as.integer(opt("B", "2000")); CORES <- as.integer(opt("cores", "2")); stopifnot(B >= 1L, CORES >= 1L)
logf <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
gap_scale <- function(v) { m <- mean(v); n <- length(v); sqrt(2 * n * sum((v - m)^2) / (n * (n - 1))) }   # single-year case of helpers_data.R gap_scale()

## ---- data and gap ----
dump <- readRDS(DUMP); D <- dump$D; ROLES <- dump$roles; Z_cols <- dump$Z_cols; G <- length(ROLES)
stopifnot(dump$rel_scale == "year", G >= 2L, all(tapply(D$chosen, D$ev_id, sum) == 1L))
P <- readRDS(POSF); j <- match(SPEC, colnames(P$R)); stopifnot(!is.na(j)); R <- P$R[, j]; key <- P$key
kidx <- setNames(seq_len(nrow(key)), paste(key$ccode, key$year))
i_row <- kidx[paste(D$c1, D$yr - 1L)]; j_row <- kidx[paste(D$c2, D$yr - 1L)]; stopifnot(!anyNA(i_row), !anyNA(j_row))
ys <- vapply(split(seq_len(nrow(key)), key$year), function(ix) gap_scale(R[ix]), 0)
## loop review P2-9: as the sweep -- the scale year (t) AND the positions year (t-1) must both vary
y_sd <- vapply(split(R, key$year), sd, 0); stopifnot(all(y_sd[as.character(unique(D$yr - 1L))] > 1e-8))
stopifnot(all(ys[as.character(unique(D$yr))] > 1e-8)); D$dR <- (R[i_row] - R[j_row]) / ys[as.character(D$yr)]; stopifnot(all(is.finite(D$dR)))
GAP <- c("dM", "dM2", "dR", "dR2", "dMR"); GAP_G <- as.vector(outer(GAP, ROLES, paste, sep = "_"))
build <- function(d) { d$dM2 <- d$dM^2; d$dR2 <- d$dR^2; d$dMR <- d$dM * d$dR
  for (g in 0:(G - 1L)) for (v in GAP) d[[paste0(v, "_", ROLES[g + 1])]] <- d[[v]] * (d$type == g); d }
## M0 pools any type-specific controls (x_<role>) as 22_rel_sweep.R does
zc0 <- Z_cols; tsg <- list()
for (z in Z_cols) for (k in ROLES) if (endsWith(z, paste0("_", k))) { pre <- substr(z, 1, nchar(z) - nchar(k) - 1L); tsg[[pre]] <- c(tsg[[pre]], z) }
tsg <- tsg[vapply(tsg, length, 0L) == G]
for (pre in names(tsg)) { D[[paste0(pre, "_pooled")]] <- Reduce(`+`, lapply(tsg[[pre]], function(z) D[[z]])); zc0 <- c(setdiff(zc0, tsg[[pre]]), paste0(pre, "_pooled")) }
fit <- function(d, terms, zc) { msg <- NULL
  d$one_ <- 1; fml <- as.formula(paste("Surv(one_, chosen) ~", paste(c(terms, zc), collapse = " + "), "+ strata(ev_id)"))
  m <- tryCatch(withCallingHandlers(coxph(fml, data = d, weights = w, method = "breslow"),
                                    warning = function(w) { msg <<- c(msg, conditionMessage(w)); invokeRestart("muffleWarning") }), error = function(e) { msg <<- paste("error:", conditionMessage(e)); NULL })
  bad <- is.null(m) || any(grepl("converge|infinite|singular|error", msg, ignore.case = TRUE)) || any(!is.finite(coef(m))) || m$iter >= 20L   # the sweep's validity policy
  if (bad) structure(list(), reason = if (length(msg)) paste(unique(msg), collapse = " | ") else "non-finite coefficients or iteration limit") else coef(m) }
u4 <- function(cf, s, r = 1) c(R_less = r * (cf[[paste0("dR", s)]] + cf[[paste0("dR2", s)]]), R_more = r * (-cf[[paste0("dR", s)]] + cf[[paste0("dR2", s)]]),
                               M_weaker = r * (cf[[paste0("dM", s)]] + cf[[paste0("dM2", s)]]), M_stronger = r * (-cf[[paste0("dM", s)]] + cf[[paste0("dM2", s)]]))
pairs <- lapply(combn(ROLES, 2, simplify = FALSE), rev)          # same orientation and labels as 22_rel_sweep.R
## per-model validity (loop review P2-9): an invalid M1 leaves the M0 quantities in place (M1 quantities NA), and vice versa
QN <- c(paste0("M0:u_", c("R_less", "R_more", "M_weaker", "M_stronger")),
        unlist(lapply(ROLES, function(k) paste0("M1:", k, ":u_", c("R_less", "R_more", "M_weaker", "M_stronger")))),
        unlist(lapply(pairs, function(pp) paste0("M1:", pp[1], "-", pp[2], ":u_", c("R_less", "R_more", "M_weaker", "M_stronger")))))
derive <- function(c0, c1) {
  out <- setNames(rep(NA_real_, length(QN)), QN)
  if (length(c0)) out[paste0("M0:u_", names(u4(c0, "")))] <- u4(c0, "")
  if (length(c1)) { for (k in ROLES) { v <- u4(c1, paste0("_", k)); out[paste0("M1:", k, ":u_", names(v))] <- v }
    for (pp in pairs) { v <- u4(c1, paste0("_", pp[1])) - u4(c1, paste0("_", pp[2])); out[paste0("M1:", pp[1], "-", pp[2], ":u_", names(v))] <- v } }
  structure(out, reason = paste(c(if (!length(c0)) paste("M0:", attr(c0, "reason")), if (!length(c1)) paste("M1:", attr(c1, "reason"))), collapse = " || ")) }
D <- build(D)
c0 <- fit(D, GAP, zc0); c1 <- fit(D, GAP_G, Z_cols)
if (!length(c0) && !length(c1)) stop("both full-sample fits fail: ", attr(c0, "reason"), " ", attr(c1, "reason"))
est <- derive(c0, c1); if (anyNA(est)) logf("full-sample fit invalid for part of the quantities: %s", attr(est, "reason"))
logf("%s on %s (%s, window %s): %d strata; spec %s", MODE, basename(DUMP), dump$types, if (is.null(dump$coorig_months)) 0 else dump$coorig_months, length(unique(D$ev_id)), SPEC)

## ---- draws (all random numbers drawn up front: reproducible under forking) ----
set.seed(123)
ev <- unique(D[, c("ev_id", "init", "disp", "type")]); stopifnot(!anyDuplicated(ev$ev_id))
if (MODE == "boot") {
  u_i <- unique(ev$init); u_d <- unique(ev$disp); ii <- match(ev$init, u_i); dd <- match(ev$disp, u_d)
  draws <- lapply(seq_len(B), function(b) { mi <- tabulate(sample.int(length(u_i), length(u_i), replace = TRUE), length(u_i))
    md <- tabulate(sample.int(length(u_d), length(u_d), replace = TRUE), length(u_d)); setNames(mi[ii] * md[dd], ev$ev_id) })
} else if (MODE == "perm_dR") {
  grp <- split(seq_len(nrow(D)), D$ev_id)
  draws <- lapply(seq_len(B), function(b) { o <- seq_len(nrow(D)); for (g in grp) o[g] <- g[sample.int(length(g))]; o })
} else {
  ## loop review P2-9: typed controls (x_<role>) would keep the original labels after relabelling, so variant dumps are refused
  if (length(tsg)) stop("--mode=perm_type is not supported on dumps with type-specific controls (", paste(names(tsg), collapse = ", "), ")")
  draws <- lapply(seq_len(B), function(b) setNames(sample(ev$type), ev$ev_id))
}
one <- function(b) {
  x <- draws[[b]]
  if (MODE == "boot") { k <- x[as.character(D$ev_id)]; Db <- D[k > 0L, ]; Db$w <- Db$w * k[k > 0L] }
  else if (MODE == "perm_dR") { Db <- D; Db$dR <- D$dR[x]; Db <- build(Db) }
  else { Db <- D; Db$type <- unname(x[as.character(D$ev_id)]); Db <- build(Db) }
  if (MODE == "boot" && length(tsg)) for (pre in names(tsg)) Db[[paste0(pre, "_pooled")]] <- Reduce(`+`, lapply(tsg[[pre]], function(z) Db[[z]]))
  ## delta review N4 (P2 remainder): a draw that loses a type still fits the pooled M0; only M1 is unavailable
  c1b <- if (any(tabulate(Db$type[Db$chosen == 1L] + 1L, nbins = G) == 0L)) structure(list(), reason = "a type has no event in this draw") else fit(Db, GAP_G, Z_cols)
  derive(fit(Db, GAP, zc0), c1b) }
res <- mclapply(seq_len(B), one, mc.cores = CORES, mc.preschedule = TRUE)
wk <- vapply(res, function(z) is.numeric(z) && length(z) == length(QN), TRUE)
if (any(!wk)) stop(sprintf("%d worker failures", sum(!wk)))
BS <- do.call(rbind, res); colnames(BS) <- QN
okq <- colSums(!is.na(BS))                                         # usable draws per quantity (per model)
rs <- vapply(res, function(z) { r <- attr(z, "reason"); if (is.null(r)) "" else r }, "")
drops <- data.frame(draw = which(nzchar(rs)), reason = rs[nzchar(rs)])
## delta review N4: no interval, positivity share or diagnostic tail value without a valid (finite) observed estimate
reportq <- report_gate(B, okq, est[QN])                             # R/lib/report_gate.R (tested directly in tests/test_resample_edge.R)
logf("usable draws: M0 %d, M1 %d of %d; draws with a dropped model %d", okq[[QN[1]]], okq[[QN[5]]], B, nrow(drops))
qf <- function(p) ifelse(reportq, apply(BS, 2, quantile, p, na.rm = TRUE), NA_real_)
if (MODE == "boot") {
  out <- data.frame(spec_id = SPEC, scheme = dump$types, mode = MODE, quantity = QN, est = unname(est[QN]), lo = qf(0.025), hi = qf(0.975),
                    share_pos = ifelse(reportq, colMeans(BS > 0, na.rm = TRUE), NA_real_), n_ok = okq, B = B, row.names = NULL)
} else {
  ## artificial-assignment DIAGNOSTIC (loop review B): not a randomization test of "no isolation mechanism"
  pv <- (1 + colSums(abs(BS) >= matrix(abs(est[QN]), nrow(BS), length(QN), byrow = TRUE), na.rm = TRUE)) / (1 + okq)
  out <- data.frame(spec_id = SPEC, scheme = dump$types, mode = MODE, quantity = QN, est = unname(est[QN]), null_lo = qf(0.025), null_hi = qf(0.975),
                    p_diag = ifelse(reportq, pv, NA_real_), n_ok = okq, B = B, row.names = NULL)
}
out$obs_valid <- is.finite(est[QN]); out$reported <- reportq        # model-specific status of the observed fit and the reporting gate (delta review N4)
out$coorig_months <- if (is.null(dump$coorig_months)) 0L else dump$coorig_months
write.csv(out, paste0(OUT, ".tmp"), row.names = FALSE)
write.csv(drops, sub("\\.csv$", "_dropped.csv", OUT), row.names = FALSE)   # draw id and reason for every draw with a dropped model
stopifnot(file.rename(paste0(OUT, ".tmp"), OUT))                                # atomic: OUT exists only when complete
print(format(out[grepl("u_R_less|u_M_stronger", out$quantity), setdiff(names(out), c("spec_id", "scheme", "mode", "B", "coorig_months"))], digits = 2), row.names = FALSE)
