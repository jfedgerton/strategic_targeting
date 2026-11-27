## R/permutation.R -- permutation tests for the headline results (PI 2026-11-27). Point fits only; seed 123 + i per permutation, so results do not
## depend on --cores. p-values are two-sided against the permutation distribution itself (stage-2 nulls need not be centred on 0):
##   p = min(1, 2 * min( (1 + #{perm >= obs}) / (1 + n), (1 + #{perm <= obs}) / (1 + n) ))
## Tests (--tests, default 1,2,3,4):
##   1 relational gaps shuffled across candidates within each stage-1 event -> stage 1 between targets (S1NF): coalition / bilateral relational
##     contrasts and coalition - bilateral differences
##   2 ties to the attacking side (tie_sum_A) shuffled across candidates within each dispute -> joint stage 2 (bbar fixed at its observed value):
##     join target side: ties to the attacker (B:tie_sum_A); originator entry: ties to the attacker (O:tie_sum_A); join attacker (A:tie_sum_A)
##   3 bbar shuffled across disputes -> joint stage 2: theta A / B / O  (does the link depend on which target a dispute is about?)
##   4 stage-1 target FE shuffled across repeat targets, then bbar rebuilt -> joint stage 2: theta A / B / O  (does it depend on target identity?)
## The joint stage 2 uses the main design (ties to both sides, S2_REF_DESIGN) plus bbar, as in joint_fe_cell().
## Output (--out): permutation_results.csv; perm_draws_<test>.rds (resumable: an existing file with the same n_perm and inputs is reused).
## Usage: Rscript R/permutation.R --build=DIR --dump=FILE --positions=FILE --spec=ID --edges_dir=DIR --out=DIR [--n_perm=1000] [--cores=8] [--tests=1,2,3,4]
suppressPackageStartupMessages({ library(data.table); library(parallel); library(survival) })
for (f in c("s2.R", "s2_fit.R", "sry.R", "s2_cell.R")) source(file.path("R/lib", f))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
BUILD <- opt("build"); DUMPF <- opt("dump"); POSF <- opt("positions"); SPEC <- opt("spec"); EDIR <- opt("edges_dir"); OUT <- opt("out")
NP <- as.integer(opt("n_perm", "1000")); CORES <- as.integer(opt("cores", "8")); TESTS <- as.integer(strsplit(opt("tests", "1,2,3,4"), ",")[[1]])
stopifnot(!is.null(BUILD), !is.null(DUMPF), !is.null(POSF), !is.null(SPEC), !is.null(EDIR), !is.null(OUT), NP >= 99, all(TESTS %in% 1:4))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); t0 <- Sys.time()
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))

## ---- 1. inputs (as main_table.R) ----
B <- s2_base(BUILD); W_REF <- B$window; dump <- readRDS(DUMPF); ROLES <- dump$roles
pz <- readRDS(POSF); key <- as.data.table(pz$key); stopifnot(SPEC %in% colnames(pz$R)); R <- pz$R[, SPEC]
X <- s2_prepare(B, key, R, s2_spec_edges("core", SPEC, EDIR)); D <- s1_prepare(dump, key, R)
ref_order <- names(sort(table(D[chosen == 1L, c2]), decreasing = TRUE))
stopifnot(all(c("tie_sum_A", "tie_sum_B") %in% names(X)), all(c("dR", "dM", "ev_id", "type") %in% names(D)))
SIG <- paste(tools::md5sum(c(DUMPF, POSF, file.path(BUILD, "s2_meta.rds"), "R/permutation.R", "R/lib/s2.R", "R/lib/s2_fit.R", "R/lib/s2_cell.R"))[], collapse = ",")

## ---- 2. helpers ----
s1_redo <- function(Dx) { Dx[, `:=`(dR2 = dR^2, dMR = dM * dR)]   # rebuild the role-specific relational terms after dR changes
  for (g in 0:1) for (v in c("dR", "dR2", "dMR")) Dx[, paste0(v, "_", ROLES[g + 1]) := get(v) * (type == g)]; Dx }
s1q <- function(Dx) { b <- coef(s1_fit_nofe(Dx, dump)); v <- c()   # stage-1 relational contrasts at gap +-1 and coalition - bilateral
  for (r in ROLES) for (sg in c(1, -1)) v[paste("S1NF", r, "R", sg, sep = ".")] <- sg * b[[paste0("dR_", r)]] + b[[paste0("dR2_", r)]]
  for (sg in c(1, -1)) v[paste("S1NF.diff.R", sg, sep = ".")] <- v[[paste("S1NF", ROLES[2], "R", sg, sep = ".")]] - v[[paste("S1NF", ROLES[1], "R", sg, sep = ".")]]
  v }
bbar_from_b <- function(b) B$TS[, .(dispnum, b = b[as.character(ccode)])][, .(bbar = mean(b, na.rm = TRUE)), by = dispnum][is.finite(bbar)]
s1_fe_b <- function() {   # observed repeat-target FE, centred and clipped (the joint_fe_cell rules)
  s1 <- s1_fit_fe(D, dump, ref_order); bb <- coef(s1$m); be <- bb[grepl("^tfe", names(bb))]; names(be) <- sub("^tfe", "", names(be))
  refb <- setdiff(intersect(ref_order, as.character(unique(D[chosen == 1L, c2]))), names(be))[1]; be <- c(setNames(0, refb), be)
  rep0 <- { n <- table(D[chosen == 1L, c2]); names(n)[n >= 2] }
  b <- be[names(be) %in% rep0 & is.finite(be)]; pmin(pmax(b - mean(b), -5), 5) }
DJ <- lapply(S2_REF_DESIGN, c, "bbar")
s2q <- function(Xj) { f <- fit_mnl(Xj, DJ, "perm", vc = FALSE)$b
  c("A:theta" = f[["A:bbar"]], "B:theta" = f[["B:bbar"]], "O:theta" = unname(f["O:bbar"]),
    "B:ties_to_attacker" = unname(f["B:tie_sum_A"]), "O:ties_to_attacker" = unname(f["O:tie_sum_A"]), "A:ties_to_attacker" = unname(f["A:tie_sum_A"])) }
pval <- function(x, o) { x <- x[is.finite(x)]; n <- length(x); min(1, 2 * min((1 + sum(x >= o)) / (1 + n), (1 + sum(x <= o)) / (1 + n))) }

## ---- 3. observed ----
b_obs <- s1_fe_b(); BT <- bbar_from_b(b_obs); XJ <- merge(X, BT, by = "dispnum")
OBS1 <- s1q(D); OBS2 <- s2q(XJ)
say("observed: coalition less central %.2f; join target ties to attacker %.2f; theta A %.2f B %.2f | %d disputes in the joint model",
    OBS1[["S1NF.coal.R.1"]], OBS2[["B:ties_to_attacker"]], OBS2[["A:theta"]], OBS2[["B:theta"]], uniqueN(XJ$dispnum))

## ---- 4. permutation draws ----
one <- function(test, i) { set.seed(123 + i)
  if (test == 1) { Dp <- copy(D); Dp[, dR := dR[sample.int(.N)], by = ev_id]; return(s1q(s1_redo(Dp))) }
  if (test == 2) { Xp <- copy(XJ); Xp[, tie_sum_A := tie_sum_A[sample.int(.N)], by = dispnum]; return(s2q(Xp)) }
  if (test == 3) { Xp <- copy(XJ); m <- unique(Xp[, .(dispnum, bbar)]); m[, bbar := bbar[sample.int(.N)]]; Xp[, bbar := NULL]; return(s2q(merge(Xp, m, by = "dispnum"))) }
  if (test == 4) { bp <- setNames(sample(unname(b_obs)), names(b_obs)); return(s2q(merge(X, bbar_from_b(bp), by = "dispnum"))) } }
run_test <- function(test) { f <- file.path(OUT, sprintf("perm_draws_%d.rds", test))
  if (file.exists(f)) { z <- readRDS(f); if (identical(z$n, NP) && identical(z$sig, SIG)) { say("test %d: reusing %s", test, f); return(z$M) } }
  t1 <- Sys.time(); L <- mclapply(seq_len(NP), function(i) tryCatch(one(test, i), error = function(e) NULL), mc.cores = CORES, mc.preschedule = FALSE)
  nm <- names(if (test == 1) OBS1 else OBS2); M <- do.call(rbind, lapply(L, function(x) if (is.numeric(x) && all(nm %in% names(x))) x[nm] else rep(NA_real_, length(nm))))
  colnames(M) <- nm; ok <- sum(stats::complete.cases(M[, !grepl("^O:", nm), drop = FALSE]))
  say("test %d: %d of %d permutations usable, %.1f min", test, ok, NP, as.numeric(difftime(Sys.time(), t1, units = "mins")))
  if (ok < 0.9 * NP) stop(sprintf("test %d: only %d of %d permutations usable", test, ok, NP))
  saveRDS(list(M = M, n = NP, sig = SIG), paste0(f, ".tmp")); if (!file.rename(paste0(f, ".tmp"), f)) stop("could not write ", f); M }
LAB <- c("1" = "relational gaps shuffled within stage-1 events", "2" = "ties to the attacker shuffled within disputes",
         "3" = "bbar shuffled across disputes", "4" = "target FE shuffled across repeat targets")
KEEP <- list("1" = names(OBS1), "2" = c("B:ties_to_attacker", "O:ties_to_attacker", "A:ties_to_attacker"), "3" = c("A:theta", "B:theta", "O:theta"), "4" = c("A:theta", "B:theta", "O:theta"))
RES <- rbindlist(lapply(TESTS, function(tt) { M <- run_test(tt); obs <- if (tt == 1) OBS1 else OBS2
  rbindlist(lapply(KEEP[[as.character(tt)]], function(q) { x <- M[, q]; x <- x[is.finite(x)]
    data.table(test = tt, shuffled = LAB[[as.character(tt)]], quantity = q, observed = obs[[q]], perm_mean = mean(x), perm_sd = sd(x),
               perm_lo = unname(quantile(x, .025)), perm_hi = unname(quantile(x, .975)), p_perm = pval(x, obs[[q]]), n_perm = length(x)) })) }))
fwrite(RES, file.path(OUT, "permutation_results.csv"))
print(RES[, .(test, quantity, observed = round(observed, 2), perm_mean = round(perm_mean, 2), p_perm = signif(p_perm, 2), n_perm)])
say("PERMUTATION_DONE %s (%.1f min)", OUT, as.numeric(difftime(Sys.time(), t0, units = "mins")))
