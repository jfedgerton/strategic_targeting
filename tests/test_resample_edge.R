## The synthetic fixture below additionally checks that pooled M0 survives draws that lose a type.
## tests/test_resample_edge.R -- delta review N4 regression test on synthetic data (seed 123).
## Case: a 3-type dump whose third type has ONE event, so the observed M1 fit is invalid while M0 is valid, and about a third of bootstrap
## draws lose that type. Required: (1) M0 quantities get intervals; (2) no M1 quantity gets an interval, share or diagnostic p (observed M1
## invalid); (3) draws that lose a type still contribute a pooled M0 (M0 usable draws = B).
## Usage: Rscript tests/test_resample_edge.R   (from the repo root)
set.seed(123)
## (0) direct test of the reporting gate (GPT review D6): an invalid observed estimate with FULLY usable draws must not be reported;
## this fails if the finite-observed condition is removed from report_gate()
source("R/lib/report_gate.R")
stopifnot("gate: NA observed, all draws usable -> not reported" = !report_gate(1000L, c(1000L, 1000L), c(NA_real_, 0.5))[1],
          "gate: valid observed, all draws usable -> reported" = report_gate(1000L, c(1000L, 1000L), c(NA_real_, 0.5))[2],
          "gate: B < 200 -> not reported" = !report_gate(100L, 100L, 0.5),
          "gate: < 90% usable -> not reported" = !report_gate(1000L, 899L, 0.5))
## mutation check: the old gate (without the finite-observed condition) must FAIL the first assertion
old_gate <- function(B, okq, est) B >= 200L & okq / B >= 0.90
stopifnot("mutation: old gate would report the NA observed estimate" = old_gate(1000L, 1000L, NA_real_))
td <- tempfile("rs_edge_"); dir.create(td)
S <- 1:20; YRS <- 2000:2010
key <- expand.grid(ccode = S, year = YRS); R <- rnorm(nrow(key)); P <- list(key = key, R = matrix(R, ncol = 1, dimnames = list(NULL, "SYN|x|pr")))
saveRDS(P, file.path(td, "pos.rds"))
mk <- function(id, ty) { yr <- sample(2001:2010, 1); i <- sample(S, 1); alt <- setdiff(S, i)
  dM <- rnorm(length(alt)); Ri <- R[key$ccode == i & key$year == yr - 1]; Rj <- R[match(paste(alt, yr - 1), paste(key$ccode, key$year))]
  u <- -0.5 * dM^2 + 0.6 * (Ri - Rj); ch <- sample(seq_along(alt), 1, prob = exp(u) / sum(exp(u)))
  data.frame(ev_id = id, chosen = as.integer(seq_along(alt) == ch), type = ty, w = 1, init = i, disp = id, init_disp = paste(i, id), dM = dM, dR = 0,
             c1 = i, c2 = alt, yr = yr) }
D <- do.call(rbind, c(lapply(1:150, mk, ty = 0L), lapply(151:300, mk, ty = 1L), list(mk(301, 2L))))
saveRDS(list(D = D, roles = c("bil", "chooser", "defender"), Z_cols = character(0), rel_scale = "year", types = "syn", coorig_months = 0L), file.path(td, "dump.rds"))
run <- function(mode) { f <- file.path(td, paste0(mode, ".csv"))
  st <- system2("Rscript", c("R/08_resample.R", file.path(td, "dump.rds"), file.path(td, "pos.rds"), shQuote("SYN|x|pr"), f, paste0("--mode=", mode), "--B=200", "--cores=2"), stdout = TRUE, stderr = TRUE)
  cat(tail(st, 3), sep = "\n"); read.csv(f) }
b <- run("boot"); p <- run("perm_dR")
m0 <- grepl("^M0:", b$quantity); m1 <- grepl("^M1:", b$quantity)
stopifnot("M0 has intervals" = all(b$reported[m0]) && all(is.finite(b$lo[m0])),
          "observed M1 invalid" = all(!b$obs_valid[m1]),
          "no M1 interval/share without a valid observed estimate" = all(is.na(b$lo[m1]) & is.na(b$share_pos[m1]) & !b$reported[m1]),
          "M0 usable in every draw (incl. draws that lose a type)" = all(b$n_ok[m0] == 200),
          "no M1 diagnostic p without a valid observed estimate" = all(is.na(p$p_diag[grepl("^M1:", p$quantity)])),
          "M0 diagnostic p reported" = all(is.finite(p$p_diag[grepl("^M0:", p$quantity)])))
cat("test_resample_edge: all checks passed\n")
