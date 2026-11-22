**Verdict: NO-GO for reporting the reference stage-2 model.** The surface construction largely mirrors stage 1, but the fitter omits a structural choice-set restriction. The covariance repair also differs from stage 1’s procedure. Both require correction and refitting before citing the reported confidence intervals or Wald statistic.

I inspected the amendment, code, tests and supplied logs. No R scripts were run and no files were modified.

The underlying unmasked multinomial likelihood and gradient are algebraically correct. Standardizing the already-constructed polynomial columns and dividing slopes by their SDs correctly restores the polynomial coefficients. `optimHess()` is evaluated on the scaled parameters, as intended. The problems below concern support, inference, model interpretation and validation.

1. **P1 — The one-shot likelihood allows third-party A entry where it is impossible.**

   **Files:** `R/s2_02_fit_reference.R:28–56`; memo N1.

   At W=1, an initiating-side entry in months 0–1 is classified as **originator entry**, not third-party A. Consequently, a state–dispute row whose horizon ends at H≤1 cannot have outcome A.

   The fitter nevertheless puts A in every row’s denominator and assigns it positive probability. This affects estimation, not merely a later prediction plot: short-horizon rows contribute likelihood against an impossible alternative.

   **Fix:** construct row-specific one-shot availability. Under the reference convention, A is available only when H>W; originator availability must follow the promotion convention, including being absent at W=0. Enforce availability in the likelihood, gradient and prediction code. When using `−Inf` utilities, avoid computing `0 × −Inf` in `Y * U`.

   Add fixtures with H=0, H=W and H=W+1, asserting both likelihood contributions and exactly zero probabilities. Refit the reference and robustness models afterward.

2. **P1 — The reported covariance uses a different PSD repair from stage 1.**

   **Files:** `R/s2_02_fit_reference.R:62–69`; `R/07_fit_sweep.R:189–193`.

   The score construction and raw CGM combination are correct:

   `M = M_candidate + M_dispute − M_candidate×dispute`.

   But the fitter clips negative eigenvalues of **M**, then applies the inverse Hessian. Stage 1 constructs the covariance first and clips its negative eigenvalues.

   These operations are not interchangeable:

   `H⁻¹ PSD(M) H⁻¹ ≠ PSD(H⁻¹ M H⁻¹)` in general.

   The supplied reference has 12 negative eigenvalues, so this is consequential. It does not establish whether individual intervals become wider or narrower.

   **Fix:** construct and transform the raw covariance to the declared coefficient parameterization, then apply the same covariance-level repair as stage 1. Record the raw eigenvalues, repair magnitude and resulting contrast variances. The documented CGM correction concerns the estimated covariance matrix. [Sandwich documentation](https://rdrr.io/rforge/sandwich/man/vcovCL.html)

   Recompute the intervals and Wald statistic. Do not describe the current χ²≈118 result as independently validated CGM inference. Complete the pigeonhole-bootstrap comparison before reporting it.

3. **P2 — The “separation rule” is not a general separation test and changes the fitted specification.**

   **File:** `R/s2_02_fit_reference.R:10–11, 39–46`.

   A regressor being constant among cases of one category does **not generally prove separation**. For example, cases can all have x=0 while noncases occur at both positive and negative x.

   For the observed O drops, the bounded/binary configurations may indeed imply boundary estimates. But deleting `likely` or `dem_x_B` sets its coefficient to zero; it does not reproduce the limiting separated model. In particular, dropping `likely` allows positive originator probabilities among unlikely candidates.

   **Fix:** treat this as an explicitly restricted specification, not an automatic mathematical repair. Either freeze the reduced reference design, with its interpretation and a suitable sensitivity, or adopt an explicitly specified separation-handling estimator. A genuinely predeclared adaptive rule can be reported honestly, but does not automatically validate ordinary post-selection Wald inference.

   Also stop `surf_contrasts()` from silently interpreting a dropped surface coefficient as an estimated zero (`:79–81`). Mark the affected contrast unavailable or label it as a restricted-model quantity.

4. **P2 — The gaps mirror stage 1; the history controls do not use the same source.**

   **Files:** `R/lib/s2.R:359–381`; builder `:218–243, 269–274`.

   The following check out:

   - **Sign:** joiner minus opponent.
   - **Opponent:** A/O face B; B faces A.
   - **Positions:** t−1.
   - **M scale:** pooled within-year pair-gap scale.
   - **R scale:** onset-year t scale for the reference.
   - Quadratic terms and interaction are formed after scaling.

   I also checked the relevant local upstream builder: `POLITY_DEM_CUT <- 6L` at `/Users/jarededgerton/strategic_targeting_relsweep/v2/R/01_build_data.R:138`, with missing democracy coded zero. Thus the assumed cutoff matches that source.

   However, stage 2 constructs MID history from **GML** (`builder:220–225`), while the upstream stage-1 builder explicitly retains **COW-based** conflict histories (`01_build_data.R:253–259, 1092`). That is a substantive source difference.

   **Fix:** use the same history source and conventions for exact mirroring, or document and approve GML histories as a stage-2 departure.

   The multi-member controls are otherwise reasonable translations, but remain assumptions: democracy share, any-member major/pact/MID indicators, closest capital distance and distinct prior disputes across members. Verify the single-member case against actual stage-1 directed-control rows, including major-power timing.

5. **P2 — “Nested linear terms” overstates what the surface contains.**

   **Files:** memo `:204–206`; `R/lib/s2.R:395–408`.

   A gap term imposes restrictions:

   `β(R_k−R_o)/s`

   constrains own and opponent position effects to equal magnitude and opposite signs. It does not freely reproduce separate own-position and opponent-position coefficients.

   Moreover:

   - The old CINC ratio refers to the **joined side**; the new surface refers to the **opposing side**.
   - A mean of log CINC shares is not log combined CINC.
   - Original target-set centrality in B’s old equation is not generally represented by B’s gap to A.

   **Fix:** describe the amendment as a new restricted opponent-surface model, not removal of mechanically redundant columns. Keeping the linear design is useful, but do not characterize it as automatically nested or interpret its likelihood difference as a standard nested-model LR test.

   No obvious exact double-counting follows merely from retaining recruitment ties and alliances. Their overlap can create correlation, which needs rank/conditioning diagnostics rather than automatic deletion.

6. **P2 — B1/B2 are diagnostics, not defensible formal permutation tests as currently specified.**

   **File:** `R/s2_02_fit_reference.R:127–147`.

   **B1:** holding the null log-likelihood fixed is computationally correct—the null contains no shuffled gap columns. The problem is exchangeability. Gap vectors are related to candidate capability, democracy, major-power status and recruitment ties that remain fixed. Uniform within-dispute shuffling does not generally represent the conditional null of zero surface coefficients given those controls.

   **B2:** refitting the null is correct. But shuffling only A/B labels among third-party entrants does not generate the null distribution of the **entire 15-term surface block**, which also concerns participation and originator effects. Controls can also make A/B probabilities unequal under a no-surface null.

   **Fix:** retain them as clearly labelled artificial-assignment diagnostics, or specify a model-based null simulation and statistic appropriate to each question. A side-choice diagnostic should use a side-choice statistic/null, rather than the complete participation-and-originator surface block.

   Do not report p=1 when no draws are valid: lines 146–147 produce exactly that result. Return unavailable status, establish a usable-draw rule and preserve failure reasons.

7. **P2 — The 40 variants are surface variants, not complete opposing-side redefinitions.**

   **Files:** builder `:281–295`; fitter `:109–115`; `R/lib/s2.R:361–392`.

   The nominal grid correctly contains 5 aggregations × 2 membership rules × 2 lags × 2 scales.

   However, the fitter replaces only the gaps. Primary-pair variants still use democracy, geography, alliances and history controls aggregated over the full baseline roster. Lag-5 variants likewise retain lag-1 recruitment/control variables.

   **Fix:** either label these precisely as **surface-only perturbations**, or rebuild the affected opponent controls consistently when changing opponent membership/timing.

   Additional qualifications:

   - `min`/`max` take M and R extrema separately; they can describe a synthetic opponent rather than one actual member. `lead` selects an actual member.
   - Primary pairs resolve for 551/560 disputes; lag-5 coverage causes additional exclusions.
   - Those exclusions change the estimation cohort. Refit the reference on each variant’s retained rows, or provide a common-cohort comparison.
   - Singleton primary sides make several aggregation levels identical. Report unique effective variants as well as the nominal 40.

   The lag-5 year-scaled surface still uses sR(t) (`s2_surface():377`), whereas stage 1’s lagged sensitivity uses the positions’ year, t−5 (`R/07_fit_sweep.R:174`). Label or align this difference.

8. **P2 — Measurement-error rebuilding is appropriate, but its scaling differs from stage 1.**

   **File:** `R/s2_02_fit_reference.R:150–162`.

   Perturbing node positions, recomputing side aggregates and rebuilding all quadratic terms is the right construction. It preserves the shared state-level perturbation across candidate and opponent roles.

   However, the fitter retains the baseline `SI$sR`. Stage 1 perturbs positions before recomputing its gap scale. These are different sensitivities: fixed baseline units versus renormalized perturbed positions.

   **Fix:** specify that choice and preferably show both if claiming comparability. Call this position-level contamination, not measurement error in the underlying tie network; recruitment edges remain unchanged. Record fitted sample sizes, failures and multiple controlled random-noise draws rather than relying on one noise realization.

9. **P2 — The current outputs are utility contrasts, not the reportable H1–H3 probability results.**

   **File:** `R/s2_02_fit_reference.R:57–83, 88–102`.

   For retained terms, the axis contrasts and linear delta-method calculation are correct: at the other gap equal to zero, the ±1 contrast is `±β_linear + β_square`.

   But these are log-odds/utility changes relative to out. They are not changes in joining probabilities, and the 15-df Wald test is **all surface coefficients zero**, not equality of A/B/O surfaces.

   Before reporting:

   - Implement the prespecified standardized predicted-probability contrasts.
   - Save a reproducible prediction object and full covariance.
   - Complete the pigeonhole bootstrap with model-specific validity reporting.
   - Report support/overlap for the ±1 endpoints.
   - Keep alternative-specific utility effects distinct from probability comparisons.

   The slopes are restored correctly, but `(int@mean)` remains a centered-feature intercept. Predicting from raw covariates requires  
   `α_raw = α_mean − Σ_j β_j mean(x_j)`,  
   or the stored centering constants. The returned fit currently does not retain those means.

10. **P2 — The supplied validation evidence does not cover the complete current fitter/variant code.**

    **Files:** `tests/test_s2.R:190–216`; `docs/s2_build_logs/test_s2.out`; `docs/s2_fit_logs/`.

    The test log ends after test 33 and `ALL PASS`; it contains no test-34 assertions. The fit artifacts contain only the reference fit, consistent with the stated robustness run still being in progress.

    Test 34’s lag fixture is also weak: it subtracts the same constant from every state’s older-year R, leaving all gaps unchanged. Looking up the wrong year could pass.

    **Fix:** rerun and archive the current tests; use nonuniform older-year changes. Add fitter tests for numerical gradients, Hessian blocks, masked likelihoods, raw-versus-standardized predictions and an independent MNL fit. Check Hessian rank/conditioning. Preserve explicit requested/completed/failed statuses for variants and permutations rather than silently excluding failed fits from summaries.

The reference can remain an exploratory result while these are fixed. **Reporting requires the support-mask refit and corrected inference first.** Windows 0 and 2–6 need not block a properly qualified W=1 result, but they must precede claims about window robustness.