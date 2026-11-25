## tests/test_s2_ties.R -- contract of the both-ties stage-2 design (GPT recheck P3, 2026-11-25).
## Part 1 (no data): own / opposing tie mappings for A, B, O; S2_REF_DESIGN_OWN has no opposing ties; boot_names() has 47 unique names incl. tie_opp.
## Part 2 (with data): on the main measure, joint_fe_draw() with unit weights equals joint_fe_cell() for every tie quantity, the opposing
## quantities equal the fitted coefficients of the opposing terms, and under S2_REF_DESIGN_OWN the opposing quantities are NA.
## Usage: Rscript tests/test_s2_ties.R [BUILD_DIR DUMP POSITIONS EDGES_DIR]   (part 2 runs only when the four paths are given)
suppressPackageStartupMessages(library(data.table))
for (f in c("s2.R", "s2_fit.R", "sry.R", "s2_cell.R")) source(file.path("R/lib", f))
a <- commandArgs(TRUE); ok_all <- TRUE
check <- function(name, cond) { cat(sprintf("%-96s %s\n", name, if (isTRUE(cond)) "PASS" else "FAIL")); if (!isTRUE(cond)) ok_all <<- FALSE }

## ---- part 1: mappings ----
check("1a own-side ties: A tie_sum_A, B tie_sum_B, O tie_sum_A",
      identical(unname(tie_terms(S2_REF_DESIGN)), c("A:tie_sum_A", "B:tie_sum_B", "O:tie_sum_A")))
check("1b opposing ties: A tie_sum_B, B tie_sum_A, O tie_sum_B",
      identical(unname(tie_opp_terms(S2_REF_DESIGN)), c("A:tie_sum_B", "B:tie_sum_A", "O:tie_sum_B")))
check("1c each option of the main design has exactly two tie terms", all(vapply(S2_REF_DESIGN, function(v) sum(grepl("^tie_", v)) == 2L, TRUE)))
check("1d own-only design: same own-side terms, no opposing terms",
      identical(unname(tie_terms(S2_REF_DESIGN_OWN)), c("A:tie_sum_A", "B:tie_sum_B", "O:tie_sum_A")) && all(is.na(tie_opp_terms(S2_REF_DESIGN_OWN))))
check("1e main design = own-only design + one opposing tie per option",
      all(vapply(c("A", "B", "O"), function(al) identical(setdiff(S2_REF_DESIGN[[al]], S2_REF_DESIGN_OWN[[al]]), sub("^.:", "", tie_opp_terms(S2_REF_DESIGN)[[al]])), TRUE)))
dmy <- list(roles = c("bil", "coal"))
bn <- boot_names(dmy)
check("1f boot_names(): 47 unique names, including A/B/O:tie_opp", length(bn) == 47L && !anyDuplicated(bn) && all(c("A:tie_opp", "B:tie_opp", "O:tie_opp") %in% bn))

## ---- part 2: point vs bootstrap extraction on real data ----
if (length(a) >= 4) {
  B <- s2_base(a[1]); W_REF <- B$window; dump <- readRDS(a[2]); pz <- readRDS(a[3]); key <- as.data.table(pz$key)
  SPEC <- "BOOKraw_ptadesta|igo+def+off+pta+dca|pr"; stopifnot(SPEC %in% colnames(pz$R)); R <- pz$R[, SPEC]
  E <- s2_spec_edges("core", SPEC, a[4]); X <- s2_prepare(B, key, R, E); D <- s1_prepare(dump, key, R)
  ref_order <- names(sort(table(D[chosen == 1L, c2]), decreasing = TRUE))
  P <- joint_fe_cell(D, dump, X, B$TS)
  disp <- sort(unique(c(as.data.table(dump$D)$disp, B$X$dispnum))); m1 <- setNames(rep(1, length(disp)), disp)
  dr <- joint_fe_draw(D, dump, X, B$TS, m1, ref_order)
  check("2a joint_fe_draw() returns exactly boot_names()", identical(names(dr), boot_names(dump)))
  check("2b unit-weight draw equals the point fit: own-side ties", isTRUE(all.equal(unname(dr[c("A:tie_sum_A", "B:tie_sum_B", "O:tie_sum_A")]), P[axis == "tie_sum", est], tolerance = 1e-6)))
  check("2c unit-weight draw equals the point fit: opposing ties", isTRUE(all.equal(unname(dr[c("A:tie_opp", "B:tie_opp", "O:tie_opp")]), P[axis == "tie_opp", est], tolerance = 1e-6)))
  check("2d opposing ties are finite in the main design", all(is.finite(P[axis == "tie_opp", est])))
  PO <- joint_fe_cell(D, dump, X, B$TS, design = S2_REF_DESIGN_OWN); dro <- joint_fe_draw(D, dump, X, B$TS, m1, ref_order, design = S2_REF_DESIGN_OWN)
  check("2e own-only design: opposing quantities NA in the point fit and in the draw", all(is.na(PO[axis == "tie_opp", est])) && all(is.na(dro[c("A:tie_opp", "B:tie_opp", "O:tie_opp")])))
  check("2f own-only design: own-side ties finite", all(is.finite(PO[axis == "tie_sum", est])))
}
if (!ok_all) stop("TEST_S2_TIES_FAILED"); cat("TEST_S2_TIES_OK\n")
