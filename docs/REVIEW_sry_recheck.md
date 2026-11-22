# SRY recheck — SRY1–SRY4

22 November 2026. **GO for fresh production reruns on the real inputs covered by the supplied successful rebuild evidence, with `REUSE_POSITIONS=0` and separate bcd/init output folders.** The numerical-tie, fit-validity and mixed-provenance blockers are addressed. One residual P2 edge case remains in SRY3: a nonconstant layer entirely removed by perfect-fit fixed effects still stops instead of returning applicable zeros. It fails closed and was not encountered in the supplied real-data rebuilds; it is not evidence of incorrect published positions.

## Scope and evidence

Read `docs/PROMPT_gpt_recheck_sry.md`, the preceding `docs/REVIEW_sry.md`, and the requested delta from `be2d07545b908be7f712930a9ab3bfec3ff8ea29` in `R/lib/sry.R`, `R/02_build_positions.R`, `run_all.sh`, `config.sh` and `tests/test_sry.R`. Unchanged grid definitions, downstream estimators and substantive modeling decisions were not re-reviewed.

Archive: `/Users/jarededgerton/Downloads/targeting_multiverse_sry_fix.zip`; HEAD `2f4e903a200bd3ce99c95e649ff3d2ca6f7b37e9`; SHA-256 `2d3b575a65ff9c21764ebfee2c58a67928089cd87378d9e1dfa2f77be586b7e3`. References below are one-based lines in that archive.

**No R scripts, R expressions, project tests or production runner were executed. Only this report was written.** Native Bash 3.2 syntax checks passed for runner/config. An isolated shell-only check confirmed that the new `&& die` reuse guard continues with value 0 and stops with value 1. Independent Python arithmetic reproduced zeroing of tiny residuals and combination outputs `(0, 0, 0.5, 1)` for the cancellation/NA example. Supplied R test passes, real-data results and solver-tolerance comparisons were assessed as external evidence, not independently reproduced.

## Fix dispositions

| Finding | Status | Evidence |
|---|---|---|
| SRY1: numerical residuals become full-strength ties | **Fixed** | `R/lib/sry.R:26–27, 55–60` applies the stated response-scale tolerance before layer-year normalization, preserving zero-valued applicable entries rather than replacing them with NA. Lines 68–70 zero tiny signed combinations before positive-part normalization. R02 calls these shared implementations at 150 and 222. |
| SRY2: unconverged fits accepted | **Fixed** | `R/lib/sry.R:31–40` enables fixest warnings, collects and rejects any warning, requires a true GLM convergence status and finite fitted values, and checks fitted/observation alignment. `glm_iter` is exposed at 15 and passed at 34–35. Errors propagate and stop the arm; they are not converted into ordinary scores. |
| SRY3: omissions silently called perfect prediction | **Partly fixed; original missing-predictor path fixed** | Lines 21–22 reject nonfinite required predictors on applicable rows. Lines 39–51 validate unique/in-range indices and reconcile omissions with removed FE groups or zero-outcome Poisson omissions; all other omissions stop. Lines 29–30 handle constant layers without estimation. The nonconstant, wholly FE-predicted case below remains uncovered. |
| SRY4: mixed layer versions | **Fixed for the requested fresh runs** | `config.sh:13` defaults reuse to 0; `run_all.sh:68–73` refuses reuse=1 before the external linking branch. Lines 77–96 build every arm against the run's layers. Arm signatures now include the helper content hash at 80; the global fingerprint includes `R/lib/*.R` at 102. No externally reused arm can silently enter this grid through the intended runner path. |

SRY1 preserves the intended denominator: genuine zero-score applicable layers still count in the signed mean, while missing layers do not. The existing no-applicable-layer convention remains edge zero. The code has not changed this into a positive-only layer average or a different network construction.

## Test coverage and remaining P2

The tests now exercise the important numerical regression directly: exact OLS fits in both models (`tests/test_sry.R:18–22`), tiny combined scores (`55–62`), constant counts, mixed FE removals, required-predictor failures, short GLM iteration limits and dyad alignment. The supplied mutation evidence discriminates the tolerance fix. The real-data rebuild and tighter-solver comparisons are particularly useful because tests of exact zero do not establish stability of ordinary nonzero residuals.

**Residual SRY3 — P2: constant outcome is sufficient, but not necessary, for all observations to be removed by perfect-fit FE groups.** References: `R/lib/sry.R:29–36`; `tests/test_sry.R:31–39`.

A concrete construction for W6/logit is one active year containing all six undirected dyads of four states, with the only positive tie between states 1 and 2. Use finite varying lcap and valid keys. In the doubled data, states 3 and 4 have only zero outcomes, so their FE groups are removed. The remaining 1–2 orientations have only ones, so recursive perfect-fit removal eliminates the remainder. The layer's outcome was nonconstant, so line 29 does not intercept it; estimation is attempted without a handler for the all-observations-removed error. This conclusion follows from the construction and fixest's documented source behavior, not an R run here. Fixest explicitly stops when FE preprocessing removes every observation. See the [upstream FE preprocessing implementation](https://github.com/lrberge/fixest/blob/master/R/fixest_env.R), particularly its `info_fe$all_removed` branch; the browsed source is upstream, not a verification of the installed 0.14.2 binary.

The supplied wholly predicted tests cover globally constant outcomes; they do not cover that nonconstant recursive-removal case. This is a remaining availability/coverage problem, **not** a new silent numerical corruption path. An error stops the arm and prevents successful publication. The reported real-data arm completions indicate it did not trigger on those inputs.

**Minimal follow-up:** add the nonconstant fixture and a deliberate all-FE-removed pathway that verifies the removal classification and returns applicable score zeros. Do not catch arbitrary fixest errors and turn them into zeros. If the intended policy instead is to refuse every nonconstant all-removed layer, state that restriction and test the refusal explicitly; then do not claim general support for every wholly perfectly predicted layer.

No additional new P1/P2 was found in the requested delta. The Poisson omission fallback remains an explicit assumption that a non-FE omitted zero is a separation removal; the supplied real-data evidence reports zero such rows, so the measured build does not depend on that branch.

## Release conditions and limits

Proceed with the requested fresh, separate bcd/init folders, using the validated real input set and rebuild mode. The earlier inaccurate numerical ties are protected at both normalization stages, invalid fits fail closed, and all arms share the current layer build. The residual all-removed case should be addressed before claiming the helper handles arbitrary sparse layers, but does not overturn the successful current-input rebuild evidence.

This GO does not authorize reusing or overwriting older output folders. Helper/arm changes can cause rebuilding before the later global fingerprint check, so a failed attempt to reuse an old folder is not guaranteed to leave it untouched. Nor does this report independently certify the supplied 22 R checks, real-data equality comparisons or final model results. Existing matching arm stamps permit reuse *within* a fresh run's lineage on restart; rebuild mode does not mean every restart recomputes every arm.
