## R/ame/ame_cell.R -- joint two-stage AME cells 1-3 (docs/ame_design.md rev 2; PI 2026-11-23 "start with just 1-3"):
##   main = joint, K = 0 | cut = stage 1 alone, then stage 2 with bbar from M stage-1 posterior draws (cut posterior) | k1 = joint, K = 1.
## Stage-1 relational gap split into between (initiator vs the target's pre-event mean) and within (target's deviation from that mean) parts,
## both with the full quadratic surface; positions centred within year and scaled by the FE gap scale (s2_gap_scale) first. Event weights enter as a pseudo-likelihood (the
## Williams-Savitsky adjustment of rev 2 section 6 is NOT applied yet; reported as raw weighted posteriors).
## Usage: Rscript R/ame/ame_cell.R --cell=main|cut|k1 --build=DIR --dump=FILE --positions=FILE --spec=SPEC --edges_dir=DIR --out=DIR
##          [--chains=4] [--warmup=1000] [--iter=1000] [--threads=2] [--sub=1] [--M=10] [--seed=123]
## Two-machine mode: --prep_only --bundle=FILE (where the builds live; no Stan) then --bundle=FILE (fit only, e.g. on Roar; no build inputs needed).
suppressPackageStartupMessages({ library(data.table); library(rstan) })
stopifnot(packageVersion("rstan") >= "2.26")   # stanc3 syntax and threads_per_chain (audit P2-8)
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); RD <- dirname(dirname(normalizePath(this_file)))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub(paste0("^--", name, "="), "", a) else default }
CELL <- opt("cell", "main"); stopifnot(CELL %in% c("main", "cut", "k1"))
OUT <- opt("out"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CH <- as.integer(opt("chains", "4")); WU <- as.integer(opt("warmup", "1000")); IT <- as.integer(opt("iter", "1000")); TH <- as.integer(opt("threads", "2"))
SUB <- as.numeric(opt("sub", "1")); MCUT <- as.integer(opt("M", "10")); SEED <- as.integer(opt("seed", "123")); SPEC <- opt("spec")
LOG <- file.path(OUT, sprintf("ame_%s.log", CELL)); logf <- function(...) { s <- sprintf(...); cat(s, "\n"); cat(s, "\n", file = LOG, append = TRUE) }
rstan_options(auto_write = TRUE, threads_per_chain = TH); set.seed(SEED)

BUNDLE <- opt("bundle"); PREP_ONLY <- !is.null(opt("prep_only"))
if (!is.null(BUNDLE) && !PREP_ONLY) stopifnot(file.exists(BUNDLE))   # fit mode never falls back to a build silently (audit P3-8)
if (!is.null(BUNDLE) && file.exists(BUNDLE) && !PREP_ONLY) {
  ## Roar / fit-only mode: the data bundle written by --prep_only on the machine that holds the builds (PI 2026-11-23: run on Roar)
  BN <- readRDS(BUNDLE); BASE <- BN$base; sc1 <- BN$sc1; X1N <- BN$X1N; ROLES <- BN$ROLES; alts <- BN$alts
  stopifnot(identical(BN$spec, SPEC) || is.null(SPEC)); SPEC <- BN$spec
  logf("AME cell %s | bundle %s (%s, built %s) | chains %d x (%d + %d) | threads %d | seed %d", CELL, BUNDLE, SPEC, BN$built, CH, WU, IT, TH, SEED)
  if (!is.null(BN$prov)) { logf("bundle provenance: commit %s | sub %.2f | %s", BN$prov$commit, BN$prov$sub, paste(sprintf("%s %s", names(BN$prov$md5), BN$prov$md5), collapse = " | ")); if (BN$prov$sub < 1) logf("WARNING: subsampled bundle (sub %.2f)", BN$prov$sub) }
  ST <- BN$ST; DISP <- BN$DISP
  mk <- function(use_s1, use_s2, fix_bbar, K, bbar = rep(0, BASE$ND)) modifyList(BASE, list(use_s1 = use_s1, use_s2 = use_s2, fix_bbar = fix_bbar, K = K, bbar_data = bbar))
} else {
  for (f in c("s2.R", "s2_fit.R", "sry.R", "s2_cell.R")) source(file.path(RD, "lib", f))   # build-side only; fit-only mode needs rstan and data.table
## ---- inputs (identical to the FE cells) -------------------------------------------------------------------------------------------
B <- s2_base(opt("build")); W_REF <- B$window; dump <- readRDS(opt("dump")); ROLES <- dump$roles
pz <- readRDS(opt("positions")); key <- as.data.table(pz$key); R <- pz$R[, SPEC]; stopifnot(!anyNA(R))
E <- s2_spec_edges(if (grepl("^positions_igo", basename(opt("positions")))) "igo" else "core", SPEC, opt("edges_dir"))
logf("AME cell %s | %s | window %d | %s | chains %d x (%d + %d) | threads %d | sub %.2f | seed %d", CELL, SPEC, W_REF, B$meta$fatal_def, CH, WU, IT, TH, SUB, SEED)

## ---- stage 1: within-year standardized positions, pre-event expanding target mean, between/within gaps ---------------------------
K1 <- data.table(ccode = key$ccode, year = key$year, R = R)[, Rt := (R - mean(R)) / s2_gap_scale(R), by = year]   # FE gap units (s2_gap_scale = sqrt(2) SD), so AME and FE contrasts are on one scale (audit P2-5)
setorder(K1, ccode, year); K1[, Rbar := cumsum(Rt) / seq_len(.N), by = ccode]            # expanding mean through year (used at t-1)
D <- as.data.table(dump$D)
if (SUB < 1) { ev <- unique(D$ev_id); D <- D[ev_id %in% sample(ev, ceiling(SUB * length(ev)))] }
ki <- match(paste(D$c1, D$yr - 1L), paste(K1$ccode, K1$year)); kj <- match(paste(D$c2, D$yr - 1L), paste(K1$ccode, K1$year)); stopifnot(!anyNA(ki), !anyNA(kj))
D[, `:=`(dB = K1$Rt[ki] - K1$Rbar[kj], dW = K1$Rbar[kj] - K1$Rt[kj])]
D[, `:=`(dM2 = dM^2, dB2 = dB^2, dW2 = dW^2, dBdW = dB * dW, dMdB = dM * dB, dMdW = dM * dW)]
SURF <- c("dM", "dM2", "dB", "dB2", "dW", "dW2", "dBdW", "dMdB", "dMdW")
for (g in 0:1) for (v in SURF) D[, paste0(v, "_", ROLES[g + 1]) := get(v) * (type == g)]
X1N <- c(as.vector(outer(SURF, ROLES, paste, sep = "_")), dump$Z_cols)
setorder(D, ev_id, row)
X1 <- as.matrix(D[, ..X1N]); sc1 <- apply(X1, 2, sd); stopifnot(all(sc1 > 0)); X1 <- sweep(X1, 2, sc1, "/")
ev <- D[, .(s = .I[1], e = .I[.N], ch = .I[chosen == 1L], w = w[1]), by = ev_id]; stopifnot(nrow(ev) == uniqueN(D$ev_id))
## ---- stage 2: the FE cell's design (s2_prepare + S2_REF_DESIGN), standardized; availability; disputes; target sets ----------------
X <- s2_prepare(B, key, R, E)
alts <- c("A", "B", "O"); AV <- mnl_avail(X, alts)
des <- lapply(alts, function(al) { v <- S2_REF_DESIGN[[al]]; ev2 <- X[y == S2_CATS_FIT[[al]]]; v[vapply(v, function(z) length(unique(ev2[[z]])) > 1L, TRUE)] }); names(des) <- alts
std <- function(M) { m <- colMeans(M); s <- apply(M, 2, sd); s[s == 0] <- 1; list(M = sweep(sweep(M, 2, m), 2, s, "/"), m = m, s = s) }
XS <- lapply(alts, function(al) std(as.matrix(X[, des[[al]], with = FALSE]))); names(XS) <- alts
y2 <- match(X$y, c("out", "A", "B", "orig")); stopifnot(!anyNA(y2))
stopifnot(all(AV[cbind(which(y2 > 1), y2[y2 > 1] - 1)]))                                  # observed category always available
DISP <- sort(unique(X$dispnum)); d2 <- match(X$dispnum, DISP)
TS <- B$TS[dispnum %in% DISP]; onset <- unique(B$row0[, .(dispnum, onset_year)])
## ---- state index, eras ------------------------------------------------------------------------------------------------------------
ST <- sort(unique(c(D$c1, D$c2, X$ccode, TS$ccode, B$R0$ccode))); J <- length(ST)
era1 <- s2_fold(D$yr); era_d <- s2_fold(onset$onset_year[match(DISP, onset$dispnum)]); stopifnot(!anyNA(era1), !anyNA(era_d))
## latent pieces (K = 1): opponent-side mean weights and the projection orthogonalizing v to the state's mean position and material share
R0 <- B$R0[dispnum %in% DISP]
Wopp <- function(side) { M <- matrix(0, length(DISP), J); s <- R0[side_rel == side]; s[, n := .N, by = dispnum]; M[cbind(match(s$dispnum, DISP), match(s$ccode, ST))] <- 1 / s$n; M }
PV <- as.data.table(B$SI$PV)[, .(Mbar = mean(M, na.rm = TRUE)), by = ccode]; RB <- K1[, .(Rb = mean(Rt)), by = ccode]
Hs <- cbind(1, RB$Rb[match(ST, RB$ccode)], PV$Mbar[match(ST, PV$ccode)]); Hs[is.na(Hs)] <- 0
Pv <- diag(J) - Hs %*% solve(crossprod(Hs), t(Hs))
mk <- function(use_s1, use_s2, fix_bbar, K, bbar = rep(0, length(DISP))) list(
  use_s1 = use_s1, use_s2 = use_s2, fix_bbar = fix_bbar, K = K, J = J, NE_ERA = 6L,
  N1 = nrow(X1), P1 = ncol(X1), X1 = X1, tgt1 = match(D$c2, ST), era1 = era1, coal1 = as.numeric(D$type == 1L), ini1 = match(D$c1, ST),
  E = nrow(ev), e_start = ev$s, e_end = ev$e, e_chosen = ev$ch, w_e = ev$w,
  N2 = nrow(X), PA = ncol(XS$A$M), PB = ncol(XS$B$M), PO = ncol(XS$O$M), XA = XS$A$M, XB = XS$B$M, XO = XS$O$M,
  y2 = y2, av2 = matrix(as.integer(AV), ncol = 3), k2 = match(X$ccode, ST), ND = length(DISP), d2 = d2,
  NT = nrow(TS), ts_j = match(TS$ccode, ST), ts_d = match(TS$dispnum, DISP), era_d = era_d, ts_n = as.numeric(tabulate(match(TS$dispnum, DISP), length(DISP))),
  bbar_data = bbar, W_oppA = Wopp("B"), W_oppB = Wopp("A"), Pv = Pv, u_set = as.numeric(ST %in% c(D$c1, X$ccode)), grainsize = 25L, grainsize2 = 2000L,
  ref_v = if (2L %in% ST) match(2L, ST) else which.max(abs(Hs[, 2])))   # K = 1 reporting sign anchor: the United States (ccode 2)
stopifnot(all(mk(1, 1, 0, 0)$ts_n > 0))
logf("stage 1: %d events, %d rows, %d terms | stage 2: %d rows, %d disputes, terms A %d B %d O %d | states %d", nrow(ev), nrow(X1), ncol(X1), nrow(X), length(DISP), ncol(XS$A$M), ncol(XS$B$M), ncol(XS$O$M), J)
  BASE <- mk(1, 1, 0, 0)
  md5 <- tools::md5sum(c(build_meta = file.path(opt("build"), "s2_meta.rds"), dump = opt("dump"), positions = opt("positions")))
  cm <- tryCatch(system2("git", c("-C", dirname(RD), "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) NA_character_)
  if (PREP_ONLY) { stopifnot(!is.null(BUNDLE)); saveRDS(list(base = BASE, sc1 = sc1, X1N = X1N, ROLES = ROLES, alts = alts, spec = SPEC, built = format(Sys.time()),
                                                    window = W_REF, fatal_def = B$meta$fatal_def, ST = ST, DISP = DISP,   # index maps for b, eta, del, bbar, u, v (audit P2-7)
                                                    prov = list(commit = cm[1], sub = SUB, md5 = setNames(unname(md5), c("build_meta", "dump", "positions")), dump_commit = dump$git_commit)), BUNDLE)
    logf("AME_BUNDLE_WRITTEN %s (%.1f MB)", BUNDLE, file.size(BUNDLE) / 1e6); quit(status = 0) }
}
STANF <- opt("stan", file.path(dirname(RD), "stan", "ame_joint.stan")); mod <- stan_model(STANF)   # per-task copy on Roar (no shared compile cache)
fit1 <- function(dat, tag, ch = CH, wu = WU, it = IT) { t0 <- Sys.time()
  ff <- file.path(OUT, sprintf("fit_%s.rds", tag))   # resume after a requeue: reuse a saved fit newer than the bundle (audit P2-6)
  if (file.exists(ff) && (is.null(BUNDLE) || file.mtime(ff) > file.mtime(BUNDLE))) { logf("%s: reusing %s", tag, ff); return(readRDS(ff)) }
  f <- sampling(mod, data = dat, chains = ch, cores = ch, warmup = wu, iter = wu + it, seed = SEED, refresh = max(1, (wu + it) %/% 10),
                control = list(adapt_delta = 0.9, max_treedepth = 11),
                pars = c("beta1", "sigma_b", "sigma_del", "sigma_eta", "alpha", "betaA", "betaB", "betaO", "theta", "sigma_eps", "b", "del", "eta", "bbar",
                         if (dat$K) c("lam1_raw", "lam2_raw", "u_out", "v_out", "lat_norm1", "lat_share1")))
  saveRDS(f, file.path(OUT, sprintf("fit_%s.rds", tag))); logf("%s: %.1f min; divergences %d; max treedepth hits %d", tag, as.numeric(difftime(Sys.time(), t0, units = "mins")),
    sum(get_divergent_iterations(f)), sum(get_max_treedepth_iterations(f))); f }

## ---- fits ----------------------------------------------------------------------------------------------------------------------------
if (CELL == "main") FITS <- list(joint = fit1(mk(1, 1, 0, 0), "main_joint_K0"))
if (CELL == "k1")   FITS <- list(joint = fit1(mk(1, 1, 0, 1), "k1_joint_K1"))
if (CELL == "cut") {
  f1 <- fit1(mk(1, 0, 0, 0), "cut_stage1")
  bb <- rstan::extract(f1, "bbar")$bbar; set.seed(SEED); pick <- sample(nrow(bb), MCUT)   # same picks on a resumed run
  f2 <- lapply(seq_along(pick), function(m) fit1(mk(0, 1, 1, 0, bb[pick[m], ]), sprintf("cut_stage2_m%02d", m), ch = CH, wu = max(250, WU %/% 2), it = max(250, IT %/% 2)))
  FITS <- list(stage1 = f1, stage2 = f2)
}

## ---- summaries: stage-1 contrasts (back-transformed), stage-2 theta, scales ------------------------------------------------------------
summ <- function(x, q) data.table(q = q, mean = mean(x), lo = quantile(x, .025), hi = quantile(x, .975))
draws1 <- function(f) { b <- rstan::extract(f, "beta1")$beta1; colnames(b) <- X1N; sweep(b, 2, sc1, "/") }   # per original unit
con <- function(b, r, ax, g) { if (ax == "B") b[, paste0("dB_", r)] * g + b[, paste0("dB2_", r)] else if (ax == "W") b[, paste0("dW_", r)] * g + b[, paste0("dW2_", r)] else b[, paste0("dM_", r)] * g + b[, paste0("dM2_", r)] }
S1f <- if (CELL == "cut") FITS$stage1 else FITS$joint; b1 <- draws1(S1f); OUTT <- list()
for (r in ROLES) for (ax in c("B", "W", "M")) for (g in c(1, -1)) OUTT[[length(OUTT) + 1]] <- summ(con(b1, r, ax, g), sprintf("S1 %s %s gap %+d", r, c(B = "between R", W = "within R", M = "material")[ax], g))
for (ax in c("B", "W")) for (g in c(1, -1)) OUTT[[length(OUTT) + 1]] <- summ(con(b1, ROLES[2], ax, g) - con(b1, ROLES[1], ax, g), sprintf("S1 coal-bil %s gap %+d", c(B = "between R", W = "within R")[ax], g))
th <- if (CELL == "cut") do.call(rbind, lapply(FITS$stage2, function(f) rstan::extract(f, "theta")$theta)) else rstan::extract(FITS$joint, "theta")$theta
for (c in 1:3) OUTT[[length(OUTT) + 1]] <- summ(th[, c], sprintf("S2 theta_%s", alts[c]))
for (p in c("sigma_b", "sigma_del", "sigma_eta")) OUTT[[length(OUTT) + 1]] <- summ(rstan::extract(S1f, p)[[1]], p)
S2f <- if (CELL == "cut") FITS$stage2 else list(FITS$joint)
se <- do.call(rbind, lapply(S2f, function(f) rstan::extract(f, "sigma_eps")$sigma_eps)); for (c in 1:3) OUTT[[length(OUTT) + 1]] <- summ(se[, c], sprintf("sigma_eps_%s", alts[c]))
if (CELL == "k1") for (p in c("lat_norm1", "lat_share1")) OUTT[[length(OUTT) + 1]] <- summ(rstan::extract(FITS$joint, p)[[1]], c(lat_norm1 = "K1 |lam1| ||u|| ||v|| (Frobenius norm of lam1 u v')", lat_share1 = "K1 latent share of stage-1 utility variance")[[p]])
S <- rbindlist(OUTT); S[, cell := CELL]
fs <- if (CELL == "cut") c(list(FITS$stage1), FITS$stage2) else list(FITS$joint)
## convergence over everything the gate names (audit P1-2): bulk / tail ESS and rank-normalized Rhat (rstan::monitor)
GP <- c("beta1", "theta", "sigma_b", "sigma_del", "sigma_eta", "alpha", "betaA", "betaB", "betaO", "sigma_eps", "lam1_raw", "lam2_raw")
mon <- rbindlist(lapply(fs, function(f) { pp <- intersect(GP, f@sim$pars_oi); pp <- pp[vapply(pp, function(p) prod(f@par_dims[[p]]) > 0 || length(f@par_dims[[p]]) == 0, TRUE)]; m <- rstan::monitor(as.array(f, pars = pp), print = FALSE, warmup = 0)
  data.table(par = rownames(m), Rhat = m$Rhat, ess_bulk = m$Bulk_ESS, ess_tail = m$Tail_ESS) }))
fwrite(mon, file.path(OUT, sprintf("ame_%s_convergence.csv", CELL)))
rh <- max(mon$Rhat, na.rm = TRUE); ne <- min(pmin(mon$ess_bulk, mon$ess_tail), na.rm = TRUE)
nd <- sum(vapply(fs, function(f) sum(get_divergent_iterations(f)), 0))
## convergence part of the design gate only (section 7): R-hat, ESS and zero divergences; SBC and misspecification checks are separate and not run here
S[, `:=`(max_rhat = rh, min_ess = ne, n_divergent = nd, convergence_gate_pass = rh <= 1.01 & ne >= 400 & nd == 0)]
fwrite(S, file.path(OUT, sprintf("ame_%s_summary.csv", CELL))); logf("max Rhat %.3f, min bulk/tail ESS %.0f, divergences %d (convergence gate: Rhat <= 1.01, ESS >= 400, 0 divergences; SBC / misspecification not part of this file)", rh, ne, nd)
print(S[, .(q, mean = round(mean, 2), ci = sprintf("[%.2f, %.2f]", lo, hi))])
logf("AME_CELL_DONE %s", CELL)
