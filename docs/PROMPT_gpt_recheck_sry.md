# Recheck: SRY1–SRY4 fixes (your `docs/REVIEW_sry.md`, NO-GO at be2d075)

Verify only the fixes: `git diff be2d075 HEAD -- R/lib/sry.R R/02_build_positions.R run_all.sh config.sh tests/test_sry.R`. Do not re-review unchanged code.

| Finding | Fix | Where |
|---|---|---|
| SRY1 numerical residuals → full-strength ties | Raw residuals with abs < 1e−8 × max(1, max abs outcome) (logit scale 1) are zero; a layer-year whose max abs residual is below that tolerance stays applicable with score 0 (not NA); the combined signed score is zeroed below 1e−10 before its positive part is normalized | `R/lib/sry.R` `sry_layer()`, `sry_combo()` |
| SRY2 unconverged fits | fixest called with `warn = TRUE`; any warning stops the arm; `convStatus` must be TRUE (GLMs); fitted values finite; `glm_iter` argument exposed for the test | `sry_layer()` |
| SRY3 unexplained omissions | Required predictors (lcap; contig for W7) must be finite on applicable rows; `obs()` indices unique and in range; every omitted row must be in a FE group listed in `fixef_removed` or (Poisson) a zero outcome removed by separation, else stop; a constant outcome over the layer gets score 0 without a fit | `sry_layer()` |
| SRY4 mixed layer versions | `REUSE_POSITIONS=1` refused for this grid; config default 0, so every arm is built from the run's layers; arm stamps include `R/lib/sry.R` | `run_all.sh` §1, `config.sh` |

Scoring moved from R/02 into `R/lib/sry.R` (sourced by R/02 and the tests); the fingerprint already hashes `R/lib/*.R`.

## Evidence (Claude, R 4.x, fixest 0.14.2)
- `tests/test_sry.R`: 22 checks, all pass. Exactly fitted layer (OLS, both models) → all scores 0, applicable; constant counts → 0, no error; noisy years keep real residuals; mixed and wholly perfectly predicted FE groups (all 0, all 1, all-0 Poisson); missing contiguity stops W7 but not W6; missing distance stops W6, ignored on non-applicable rows; `glm_iter = 1` stops logit and Poisson; cancellation and 1e−13 residual → 0; NA layers excluded from the mean; planted tie recovered for the right dyad.
- Mutation check: with the three tolerance lines removed, tests 1 (both models) and 5 (tiny residual) fail.
- Real data, all four sry arms incl. in force, plus W1–W4 `ptaany`: build completes with no stop; every omitted row is explained by removed FE groups (e.g. W6 neu 1,563,346 doubled rows), 0 Poisson separations, 0 flat layer-years, GLM iterations 5–11. Positions vs the pre-fix build: W1–W4 any and W7 all identical; others max abs diff ≤ 1.3e−8.
- Solver tolerance (main layer set, glm.tol 1e−12, fixef.tol 1e−10 vs defaults): max layer-score change 8e−6, max edge change 7.5e−7, positive-edge sets identical, max log-PageRank change 1.1e−6, within-year Spearman ≥ 0.99998 (W6 and W7).

Write `docs/REVIEW_sry_recheck.md`: SRY1–SRY4 fixed / partly / not, with file:line evidence; any NEW P1/P2; go / no-go for fresh production reruns (`REUSE_POSITIONS=0`, separate OUT folders for `TYPES=bcd` and `TYPES=init`). Read-only otherwise.
