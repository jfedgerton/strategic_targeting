suppressPackageStartupMessages({ library(data.table); library(survival) })
d <- readRDS("/home/claude/align/long_init_bon_w1d_relyear_fmid.rds"); D <- as.data.table(d$D); Z <- d$Z_cols
pz <- readRDS("/home/claude/align/positions_cmp.rds"); key <- as.data.table(pz$key); R <- pz$R[, "W2contig|noPTA|pr"]
gs <- function(v) { n <- length(v); sqrt(2 * n * sum((v - mean(v))^2) / (n * (n - 1))) }; sc <- tapply(R, key$year, gs); kk <- paste(key$ccode, key$year)
D[, dR := (R[match(paste(c1, yr - 1L), kk)] - R[match(paste(c2, yr - 1L), kk)]) / sc[as.character(yr)]]
D[, `:=`(dM2 = dM^2, dR2 = dR^2, dMR = dM * dR)]; GAP <- c("dM", "dM2", "dR", "dR2", "dMR")
for (g in 0:1) for (v in GAP) D[, paste0(v, "_", d$roles[g + 1]) := get(v) * (type == g)]
ch <- D[chosen == 1]; tg <- unique(ch$c2); D[, tfe := relevel(factor(ifelse(c2 %in% tg, as.character(c2), "ref")), "ref")]
GT <- as.vector(outer(GAP, d$roles, paste, sep = "_"))
for (fe in c(FALSE, TRUE)) { f <- as.formula(paste("Surv(rep(1, nrow(D)), chosen) ~", paste(c(GT, Z, if (fe) "tfe"), collapse = " + "), "+ strata(ev_id)"))
  m <- suppressWarnings(coxph(f, data = D, weights = w, method = "breslow")); b <- coef(m)
  u <- function(r, M, Rr) b[paste0("dM_", r)] * M + b[paste0("dM2_", r)] * M^2 + b[paste0("dR_", r)] * Rr + b[paste0("dR2_", r)] * Rr^2 + b[paste0("dMR_", r)] * M * Rr
  cat(if (fe) "WITH target FE\n" else "NO FE\n")
  for (r in d$roles) cat(sprintf("  %s: u(dM=+1)-u(0) %.2f | u(dM=-1)-u(0) %.2f | u(dR=+1)-u(0) %.2f | u(dR=-1)-u(0) %.2f\n", r, u(r, 1, 0), u(r, -1, 0), u(r, 0, 1), u(r, 0, -1)))
}
cat("observed chosen-target gaps (median dM, dR) by type:\n"); print(ch[, .(dM = median(dM), dR = median(dR), n = .N), by = type])
