**Verdict: GO with fixes for fitting the reference one-shot model.** Revision 2 resolves the previous P1 blockers. Two P2 issues remain: endpoint provenance and side-major-power timing. The reference formula also needs to distinguish model terms from availability checks.

The package is unchanged from the previous recheck. No R code was executed and no files were modified.

### Previous findings

| Finding | Resolved? | Evidence |
|---|---|---|
| **P1-1: spell gaps bridged** | **Yes** | `R/lib/s2.R:271–281` caps exits at December of `spell_last`. Builder `:62–72` checks inter-spell gaps. Test 22 asserts absence during 1941 in the previously failing example. |
| **P1-2: missing reference covariates** | **Yes**, subject to formula clarification below | Positions are required at builder `:20–24`; covariates are constructed at `:211–239`; `REF_COLS` and finite-value assertions are at `:255–267`. |
| **P2-3: acceptance not fail-closed** | **Yes in release mode** | Builder `:24, 28, 118–134` requires acceptance inputs, reconciles projection records and stops on mismatches or missing/ambiguous spells. |
| **P2-4: missing treated as zero/partial totals** | **Yes for the identified paths** | `s2_ties()` distinguishes uncovered years; builder `:225–228` uses strict aggregation, followed by explicit cohort exclusions. |
| **P2-5: population and endpoint policy undecided** | **Yes**, except endpoint labelling below | Memo `:190–193` approves the population and endpoint construction. Builder `:43–58` exports population and disagreement records. |
| **P2-6: acceptance coverage and fixtures** | **Partly** | Projection denominators and fixtures are improved. External projections for W=0 and W=2–6 remain pending; some new checks are still shallow. |
| **P3-7: conditional side mask** | **Yes** | `R/lib/s2.R:295–298` returns exactly zero for unavailable A; test 27 checks this. |

### Numbered findings

1. **P2 — Capped endpoints can be labelled “exact GML month.”**

   **File:** `R/s2_01_build_riskset.R:76–80`.

   Dispute end correctly equals the latest reconciled participant exit. However, `end_imputed` considers only `exit_imputed`, ignoring `exit_capped` and `exit_before_entry`.

   Example: A’s GML episode ends in June 1942, but its stage-1 spell ends in 1940. The cap produces December 1940. If B exits in November 1940, the dispute endpoint becomes December 1940, yet `end_imputed=FALSE` and the log calls it exact. No unmodified GML exit supplies that date.

   **Fix:** add endpoint-source flags such as `end_capped` and `end_adjusted`. Call an endpoint exact only when an unmodified GML exit supports the selected maximum. Keep `end_imputed`’s narrower fallback meaning if desired, but distinguish capped endpoints in reporting and the endpoint sensitivity.

   This is a provenance problem, not a reason to undo the approved cap. The number of currently affected dispute endpoints needs a census.

2. **P2 — Side-major-power status uses onset timing rather than the declared t−1 convention.**

   **File:** `R/s2_01_build_riskset.R:204–205, 222–228`.

   `maj_at()` evaluates membership at `onset_lo`; reusing it for `maj_s` makes `major_side_A/B` onset measures. Other side covariates now correctly use the prior year.

   A roster member changing major-power status between the lagged year and onset receives different classifications under these rules.

   **Fix:** either calculate side-major status at a declared t−1 date or approve onset timing as an explicit exception. Candidate-major status for the relevance tier can remain onset-based.

3. **P2 — `REF_COLS` is complete as an availability schema, but should not become the formula automatically.**

   **File:** `R/s2_01_build_riskset.R:255–258`.

   It now covers the requested reference concepts: ties, target and own centrality, relative capability, power share, geography, alliances, prior conflict, side size, side-major status and relevance.

   Two exported quantities need separation from reference regressors:

   - **`tie_max_A/B`:** useful sweep inputs; the reference uses sums.
   - **`log_cinc_k`:** useful input/diagnostic, but the memo specifies capability relative to each side rather than an additional standalone own-capability effect.

   Neither is outcome leakage. Their availability requirements appear redundant with the reference inputs here, but fitting every listed column would change the specified model.

   **Fix:** define separate availability and model-term lists. Freeze the alternative-specific design, including tier interactions and originator utility terms. Do not automatically include both sum and maximum ties or standalone own capability.

   `share_B` is exported and recoverable as `1−share_A`; requiring only `share_A` is sufficient. Outcome, keys, category availability and cohort flags need structural assertions separately.

4. **P3 — The prior-MID suffix join is correct; clarify what the side aggregate counts.**

   **File:** `R/s2_01_build_riskset.R:217–228`.

   There is no `.x`/`.y` reversal:

   - `dispnum.x` identifies the current dispute.
   - `dispnum.y` identifies the historical dispute.
   - Their inequality excludes the current dispute.
   - `styear <= ty` restricts history to disputes beginning before the onset year.
   - `uniqueN(dispnum.y)` counts distinct prior disputes per candidate–member pair.

   But `sum(prior_n)` counts **pair-dispute histories**, not distinct disputes against the side. One previous MID involving k and two roster members contributes twice.

   **Fix:** document that interpretation, or deduplicate historical dispute IDs across roster members if distinct side-level disputes were intended. Add a fixture covering repeated rows, current/future disputes and a historical MID involving multiple roster members.

5. **P3 — Reference-cohort exclusion order is correct, but its counts are hierarchical.**

   **File:** `R/s2_01_build_riskset.R:260–267`.

   Incomplete horizons are excluded first, followed by the first failing covariate in list order. This is deterministic and does not select using the observed outcome.

   Therefore “missing own centrality 414” is a first-exclusion-reason count, not total own-centrality missingness.

   **Fix:** retain the mutually exclusive partition, label it clearly and publish nonexclusive missingness counts alongside it. The fitter should subset `in_ref` explicitly and use `na.action=na.fail`.

   The supplied accounting reconciles: **65,416 → 63,143 rows**, excluding **2,273**. Reference outcomes are **A 121, B 184, originator 71, out 62,767**.

6. **P3 — Acceptance is substantially stronger; complete the remaining targeted tests.**

   **Files:** builder `:103–134, 178–183`; `tests/test_s2.R:139–164`.

   The corrected projection accounting is sound:

   - Stage-1 set: **870 = 869 matched + 1 later-spell record**.
   - GML set: **870 = 707 matched + 1 later spell + 162 outside sample**.

   Remaining limitations:

   - The promotion-boundary assertion largely restates the classification function.
   - Test 25 checks a spanning dispute’s onset fold, not expanded monthly-row assignments.
   - Strict aggregation, the prior-MID join and exclusion precedence lack dedicated fixtures.
   - All-year PageRank equality is strong replication evidence, but does not uniquely establish elementwise edge equality.

   **Fix:** add the aggregation, history and exclusion fixtures before freezing the reference fitter. Complete the other-window projections before claiming acceptance of the full sweep; they need not block the checked W=1 reference.

The CSVs reconcile to **440 shared, 102 GML-only and 120 stage-1-only disputes**, and **75 status, 33 severity-level and 184 unknown-code differences**.

**Proceed once findings 1–3 are resolved or explicitly amended in the contract.** No remaining demonstrated P1 requires another infrastructure redesign. The deferred q model, regime inputs and monthly/sweep infrastructure do not block the λ=0 one-shot reference.