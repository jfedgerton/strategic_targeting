## ============================================================================
## 29_rel_influence.R  --  leave-one-out influence for ONE relational position (2026-11-19, multiverse loop)
## Refits the type model (M1) once per left-out unit (--unit=disp: every dispute, all of its events, matched on ANY source dispute when the
## dump carries dispnums; --unit=init: every initiator) and records the contrasts named by --q for each refit. Writes:
##   <out>_loo.csv      one row per unit: n events removed, each contrast without the unit, its change, CGM interval of the refit
##   <out>_top.txt      the units whose removal lowers --rank_q the most (top 10, in order), for the drop variants of 24_rel_variants.R
## Gap: positions t-1, year-t scale (as 22_rel_sweep.R defaults). Seed not needed (deterministic).
## Usage: Rscript v2/R/29_rel_influence.R <long_dump.rds> <positions.rds> <spec_id> <out_prefix> [--unit=disp|init]
##        [--q=chooser-bil:u_R_less,defender-chooser:u_R_less] [--rank_q=chooser-bil:u_R_less] [--cores=2]
## ============================================================================
suppressMessages({ library(survival); library(parallel); library(data.table) })
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) >= 4)
DUMP <- args[1]; POSF <- args[2]; SPEC <- args[3]; OUTP <- args[4]
fl <- sub("^--([A-Za-z_]+)=.*$", "\\1", grep("^--", args, value = TRUE)); KN <- c("unit", "q", "rank_q", "cores")
if (any(!fl %in% KN)) stop("unknown option(s): ", paste(setdiff(fl, KN), collapse = ", "))
opt <- function(name, default) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
UNIT <- opt("unit", "disp"); stopifnot(UNIT %in% c("disp", "init"))
Q <- strsplit(opt("q", "chooser-bil:u_R_less,defender-chooser:u_R_less"), ",")[[1]]; RQ <- opt("rank_q", Q[1]); stopifnot(RQ %in% Q)
CORES <- as.integer(opt("cores", "2"))
logf <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
gap_scale <- function(v) { m <- mean(v); n <- length(v); sqrt(2 * n * sum((v - m)^2) / (n * (n - 1))) }

dump <- readRDS(DUMP); D <- dump$D; ROLES <- dump$roles; Z_cols <- dump$Z_cols; G <- length(ROLES)
stopifnot(dump$rel_scale == "year")
P <- readRDS(POSF); j <- match(SPEC, colnames(P$R)); stopifnot(!is.na(j)); R <- P$R[, j]; key <- P$key
kidx <- setNames(seq_len(nrow(key)), paste(key$ccode, key$year))
i_row <- kidx[paste(D$c1, D$yr - 1L)]; j_row <- kidx[paste(D$c2, D$yr - 1L)]; stopifnot(!anyNA(i_row), !anyNA(j_row))
ys <- vapply(split(seq_len(nrow(key)), key$year), function(ix) gap_scale(R[ix]), 0)
D$dR <- (R[i_row] - R[j_row]) / ys[as.character(D$yr)]; stopifnot(all(is.finite(D$dR)))
GAP <- c("dM", "dM2", "dR", "dR2", "dMR"); GAP_G <- as.vector(outer(GAP, ROLES, paste, sep = "_"))
D$dM2 <- D$dM^2; D$dR2 <- D$dR^2; D$dMR <- D$dM * D$dR
for (g in 0:(G - 1L)) for (v in GAP) D[[paste0(v, "_", ROLES[g + 1])]] <- D[[v]] * (D$type == g)
## contrast vector for "a-b:u_R_less" style or "a:u_R_less"
cvec <- function(q, nm) { parts <- strsplit(q, ":")[[1]]; ty <- strsplit(parts[1], "-")[[1]]; u <- parts[2]
  ax <- if (grepl("_R_", u)) "R" else "M"; r <- if (grepl("less|weaker", u)) 1 else -1
  a <- setNames(rep(0, length(nm)), nm); s <- c(1, -1)
  for (t in seq_along(ty)) { a[paste0("d", ax, "_", ty[t])] <- a[paste0("d", ax, "_", ty[t])] + s[t] * r; a[paste0("d", ax, "2_", ty[t])] <- a[paste0("d", ax, "2_", ty[t])] + s[t] }
  a }
vc_cl <- function(m, g) { Dd <- rowsum(residuals(m, type = "dfbeta", weighted = TRUE), g, reorder = FALSE); crossprod(Dd) }
## loop review P2-8: the sweep's validity policy -- a refit with a convergence / infinite / singular warning, an error, non-finite coefficients or
## the iteration limit is INVALID (recorded with its reason, never used for ranking)
fitq <- function(d) {
  d$one_ <- 1; fml <- as.formula(paste("Surv(one_, chosen) ~", paste(c(GAP_G, Z_cols), collapse = " + "), "+ strata(ev_id)"))
  w_ <- character(0)
  m <- tryCatch(withCallingHandlers(coxph(fml, data = d, weights = w, method = "breslow"), warning = function(e) { w_ <<- c(w_, conditionMessage(e)); invokeRestart("muffleWarning") }),
                error = function(e) { w_ <<- c(w_, paste("error:", conditionMessage(e))); NULL })
  if (is.null(m)) return(structure(list(), reason = paste(w_, collapse = " | ")))
  if (any(grepl("converge|infinite|singular", w_, ignore.case = TRUE)) || any(!is.finite(coef(m))) || m$iter >= 20L)
    return(structure(list(), reason = if (length(w_)) paste(unique(w_), collapse = " | ") else "non-finite coefficients or iteration limit"))
  b <- coef(m); V <- vc_cl(m, d$init) + vc_cl(m, d$disp) - vc_cl(m, d$init_disp); e <- eigen(V, symmetric = TRUE); V <- e$vectors %*% diag(pmax(e$values, 0)) %*% t(e$vectors)
  unlist(lapply(Q, function(q) { a <- cvec(q, names(b)); est <- sum(a * b); se <- sqrt(drop(t(a) %*% V %*% a)); setNames(c(est, est - 1.96 * se, est + 1.96 * se), paste0(q, c("", "_lo", "_hi"))) })) }
## validate the requested contrasts: every named type must exist, and the quantity must be one of the four endpoints
for (q in Q) { pq <- strsplit(q, ":")[[1]]; tq <- strsplit(pq[1], "-")[[1]]   # delta review: one type or a two-type difference only (cvec handles no more)
  stopifnot(length(pq) == 2L, length(tq) %in% 1:2, !anyDuplicated(tq), all(tq %in% ROLES), pq[2] %in% c("u_R_less", "u_R_more", "u_M_weaker", "u_M_stronger")) }
full <- fitq(D); if (!length(full)) stop("the full-sample fit is invalid: ", attr(full, "reason"))
ev <- unique(D[, c("ev_id", "init", "disp", if ("dispnums" %in% names(D)) "dispnums")])
units <- if (UNIT == "init") sort(unique(ev$init)) else sort(unique(if ("dispnums" %in% names(ev)) as.integer(unlist(strsplit(ev$dispnums, ";"))) else ev$disp))
drop_ev <- function(u) if (UNIT == "init") ev$ev_id[ev$init == u] else if ("dispnums" %in% names(ev)) ev$ev_id[vapply(strsplit(ev$dispnums, ";"), function(x) as.character(u) %in% x, TRUE)] else ev$ev_id[ev$disp == u]
logf("LOO by %s: %d units; full-sample %s", UNIT, length(units), paste(sprintf("%s %.2f", Q, full[Q]), collapse = "; "))
res <- mclapply(units, function(u) { de <- drop_ev(u)
  if (!length(de)) return(data.table(unit = u, n_removed = 0L, status = "no events"))
  d <- D[!(D$ev_id %in% de), ]
  if (any(tabulate(d$type[d$chosen == 1L] + 1L, nbins = G) == 0L)) return(data.table(unit = u, n_removed = length(de), status = "a type vanishes"))
  r <- fitq(d); if (!length(r)) return(data.table(unit = u, n_removed = length(de), status = paste("invalid:", attr(r, "reason"))))
  cbind(data.table(unit = u, n_removed = length(de), status = "ok"), as.data.table(as.list(r))) }, mc.cores = CORES)
bad <- !vapply(res, is.data.frame, TRUE); if (any(bad)) stop(sprintf("%d worker failures", sum(bad)))
L <- rbindlist(res, fill = TRUE); stopifnot(nrow(L) == length(units))     # every attempted unit, with its status
for (q in Q) L[[paste0(q, "_change")]] <- L[[q]] - full[[q]]
setorderv(L, paste0(RQ, "_change"), na.last = TRUE)
top <- head(L[status == "ok", unit], 10)
f1 <- paste0(OUTP, "_loo.csv"); f2 <- paste0(OUTP, "_top.txt")
fwrite(L, paste0(f1, ".tmp")); writeLines(as.character(top), paste0(f2, ".tmp"))
stopifnot(file.rename(paste0(f2, ".tmp"), f2), file.rename(paste0(f1, ".tmp"), f1))   # top first, so an existing _loo.csv implies a complete pair
logf("%d of %d refits valid (%s); largest drops in %s: %s", sum(L$status == "ok"), length(units), paste(names(table(L$status[L$status != "ok"])), table(L$status[L$status != "ok"]), collapse = "; "),
     RQ, paste(sprintf("%d (%.2f)", as.integer(head(top, 5)), head(L[status == "ok"][[paste0(RQ, "_change")]], 5)), collapse = ", "))
