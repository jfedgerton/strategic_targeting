**Verdict: NO-GO for fitting the reference model yet.** The core classification, prospective masks and raw-tie construction largely follow the contract. The remaining blockers concern spell-end reconciliation, missing reference covariates and the final estimation population. The infrastructure is usable for further development.

I read the contract, notes, library, tests, builder and both build logs in the requested order. No R code was run and no files were modified. References below are to files inside the supplied ZIP. Numerical build results are supplied evidence, not independently reproduced results.

1. **P1 — GML exit matching can bridge a known participant-spell gap.**

   **File:** `R/s2_01_build_riskset.R:67–73, 132–136`.

   The non-equi join matches GML participant records to each stage-1 spell using **year overlap**, then assigns the maximum matching GML end month. It never checks that this end respects the stage-1 spell’s last active year.

   Concrete failure case: stage-1 spells are 1940 and 1942, separated by an inactive 1941; a GML participant record spans 1940–1942. Both stage-1 spells match that record, and the first receives a 1942 exit. The resulting roster incorrectly includes the participant during 1941. Ignoring re-entry **decisions** does not authorize bridging gaps in roster membership.

   **Fix:** reconcile source episodes explicitly, retain conflict flags and assert that a derived exit cannot cross a known inactive interval. At minimum, enforce the spell’s known year bounds while resolving more precise endpoints. Add a fixture with two separated spells and an overlapping GML record, asserting absence during the gap.

   This is a concrete code failure condition; the number of affected production spells has not been established. Audit the seven second spells before fitting.

2. **P1 — The exported data do not yet contain the complete reference specification.**

   **Files:** `R/s2_01_build_riskset.R:197–237`; memo `§7:54–61`; `docs/stage2_build_notes.md:48–57`.

   The omissions extend beyond those listed in the notes. The exported one-shot table lacks:

   - The candidate’s own relational centrality.
   - Prior-MID covariates.
   - Side-specific contiguity.
   - Explicit side-size and side-major-power covariates.

   `major` and `any_contig` in the tier calculation describe candidate relevance across baseline participants; they do not supply the required side-specific predictors.

   Target centrality also depends on the nominally optional `--verify_positions` argument. Without it, lines 217–223 deliberately generate an entirely missing target-centrality column while allowing `S2_BUILD_DONE`.

   **Fix:** separate required position inputs from optional numerical verification. Define a reference-model schema and assert its required columns, timing, transformations and availability before permitting a reference fit. Build the missing predictors or obtain an explicit specification amendment; do not let the fitting script silently substitute a reduced model.

   Monthly non-tie covariates may wait until the monthly model, but the missing **one-shot reference** predictors cannot.

3. **P2 — Acceptance checks report failures without rejecting the build, and coverage can shrink silently.**

   **File:** `R/s2_01_build_riskset.R:91–114`.

   Three problems remain:

   - Lines 112–113 write mismatches and log `ACCEPTANCE FAILED`, but execution continues to publication and `S2_BUILD_DONE`.
   - The inner merge at 102 discards unmatched projected records. There is no assertion that every eligible record matched.
   - The spell lookup at 96–97 reduces multiple matches with `min(spell)`, rather than rejecting ambiguous matches. Missing matches survive into the comparison.

   **Fix:** use an explicit left-join reconciliation with mutually exclusive reasons: outside selected dispute set, later spell, missing participant, ambiguous spell, or matched. Require exactly one covering spell for included records. Stop on mismatches, unexpected unmatched records or ambiguous matches. Make acceptance mandatory for a release build, with a separately named exploratory mode if needed.

   The supplied W=1 matches are useful evidence; the current checks are not yet a reliable release gate.

4. **P2 — Missing covariates can become ordinary numerical values or silently change the estimation population.**

   **Files:** `R/s2_01_build_riskset.R:205–223`; `R/lib/s2.R:179–188`.

   Concrete cases:

   - `sum(cinc, na.rm=TRUE)` treats a side with capabilities `{0.02, NA}` as having capability `0.02`. Its ratios and power shares are numerical despite incomplete side capability.
   - `any(ally_t1 > 0, na.rm=TRUE)` returns false when all relevant observations are missing.
   - Missing target positions propagate to entire disputes. The stage-1-aligned log reports 12 such disputes.
   - A missing tie year is indistinguishable from observed zero edges. The test at `tests/test_s2.R:94` explicitly accepts zero for an absent year. That convention is unsafe for unavailable coverage, especially lag-5 cells.

   **Fix:** distinguish observed zero, structural nonapplicability and missing coverage. Prespecify treatment of partial side CINC and missing alliances; retain reason flags. Establish the reference estimation cohort before fitting and report exclusions by outcome, dispute and era. Use `na.action = na.fail` after that preparation.

   For ties, pass a coverage table separately from the sparse positive-edge table. An absent edge within a covered year can be zero; an unavailable year must be flagged or rejected. The logs say the reference lag-1 precoverage count is zero, so this last issue is principally a lag-sweep readiness condition.

   Also resolve the geography timing: minimum capital distance uses `onset_year` at lines 199–202, whereas the memo’s general covariate convention is t−1. Either change the lookup or document geography as an exception.

5. **P2 — Approve the fatal-dispute population and endpoint fallback policy before fitting.**

   **Files:** `R/s2_01_build_riskset.R:31–57`; `docs/stage2_build_notes.md:42–46`; both build logs.

   **I support stage-1 alignment as the reference**, with GML fatality as a labelled alternative. It directly supports the intended comparison between the stages. But the accurate label is currently **“primary disputes represented by stage-1 fatal-event flags”**—not interchangeable with GML-fatal disputes or necessarily every source dispute of collapsed events.

   Report:

   | Population accounting | Count |
   |---|---:|
   | GML-fatal disputes before participant availability filtering | 542 |
   | Stage-1-fatal primary disputes | 560 |
   | Intersection | 440 |
   | GML-only | 102 |
   | Stage-1-only | 120 |
   | GML build after dropping dispute 266, absent from participation data | 541 |

   Supply the ID lists and distinguish population selection from subsequent availability exclusions. Decompose the logged **292 fatality-code disagreements** into fatal/nonfatal changes, severity-level changes and unknown-code disagreements; they are not necessarily 292 fatal-status reversals.

   The stage-1-aligned build additionally uses December-of-last-year fallback ends for **33 disputes and 116 participant spells**. Those are imputed endpoints, not exact GML dates. They affect exposure, rosters and potentially whether entries fall within follow-up.

   **Fix:** obtain the stated PI decision, amend B5, record the chosen definition in output metadata, and approve a flagged endpoint policy. Include an endpoint sensitivity or exclusion comparison. Validate month codes before converting them to integer month keys; invalid month codes must not become apparently valid dates through arithmetic.

6. **P2 — The projections are sensible, but their reported coverage needs precise denominators and broader acceptance tests.**

   **Files:** `R/s2_01_build_riskset.R:93–111`; `docs/s2_build_logs/s2_build_stage1.log:14–19`; GML log `15–20`.

   Excluding the USSR’s later spell in dispute 258 from a **first-entry classification projection** is appropriate. Keep it in roster-history testing: an entry excluded from the projection can still matter for subsequent side membership.

   However, the displayed “events” are deduplicated projection records. The code removes target identity and eventually year; they are not the original target-choice event count. Name the denominator explicitly.

   The **691 versus 853** difference is consistent with the selected dispute populations. The stage-1-aligned build matches all 853 projected records, with zero reported role mismatches; the GML build matches 691. These results establish W=1 role agreement on the matched projection, not complete agreement of dates, horizons, rosters or covariates.

   Before release, add:

   - Fail-closed acceptance and exact unmatched-record accounting from finding 3.
   - Matching-dump projections for W=0 and W=2–6 when those dumps exist.
   - Baseline-set equality on the full participant universe, not only fitted initiators.
   - First-entry date and promotion-boundary checks, including uncertain dates.
   - Join cardinality tests for overlapping system spells and multiple candidate GML episodes.
   - Elementwise edge comparisons or broader year coverage. The supplied PageRank checks cover six years; matching PageRank is useful but does not uniquely establish identical edges.

   **Fixture coverage:** most requested core fixtures exist, but not every named fixture has a full exact assertion.

   | Required fixture | Assessment |
   |---|---|
   | Same-month non-original B retained; A differs at W=0/W=1 | Covered: tests 2–3 |
   | Prospective masks independent of eventual outcome; no-promotion dispute | Covered at data-mask level: tests 4–6 |
   | Participant without fitted event | Architectural/classification coverage: test 7; no end-to-end fitted-event projection fixture |
   | Auxiliary/sequential masks | Shared helper tested at 8; model integration remains unbuilt |
   | Initial and empty rosters | Exact membership assertions at 13–14; end-policy test lists empty months but does not verify follow-up truncation |
   | System exit and administrative truncation | Basic flags/horizons covered at 10; no exact matched horizon-standardized prediction comparison |
   | 1914 dispute spanning folds | Partial: test 16 checks the year map and 1914, not an actual multi-year dispute’s row assignment |
   | Mixed fatalities with unknowns | Main branches covered at 18; add known fatalities on both sides **plus additional unknown participants** |
   | Positive originator hazard changes survival | Covered at 19 |
   | Recurrent return | Deliberately absent; acceptable while that family remains unbuilt |

   Add January/December month-key and inclusive-end fixtures. The convention itself—zero-based months, entry-month outcome, and rows through H inclusive—is internally coherent.

7. **P3 — Core contract fidelity is otherwise good; preserve the distinction between infrastructure and fit readiness.**

   **Files:** `R/lib/s2.R:30–54, 101–174, 196–236`.

   The inspected implementation correctly follows several important rules:

   - R1 baseline status comes from `orig == 1`, not month equality.
   - N1 masks depend only on month/window; W=0 has no originator alternative, and A/originator availability switches correctly.
   - N2 classification uses participants rather than fitted-event availability.
   - Cohort and roster non-equi joins use explicit `x.`/`i.` references for the relevant original endpoints; I found no demonstrated wrong-column substitution there.
   - Month-zero rosters use baseline participants; subsequent rosters lag active spells by one month. Carry-forward is explicit and flagged.
   - Raw sum/max ties are retained without yearly standardization.
   - Frozen relevance tiers implement the stated OR rule across contiguity, region, major-power status and alliance.
   - Jointly fatal yes/no/unresolved logic matches the approved bounds.
   - Flat and sequential CIF helpers include originator survival loss.

   Before model execution, test that the **conditional side equation** also respects `avail_A`; testing the auxiliary T/O mask alone does not prevent third-party A probability inside the promotion window.

   **Required before the reference one-shot fit:** findings 1–5, complete reference predictors, approved sample/missingness policy, the complete-horizon estimation subset, and release assertions.

   **Can wait for their respective sweep/model:** late landmarks, recurrence, dynamic entrants, monthly-updated tiers, two-month aggregation, switch-censoring, and additional relational measures. Monthly non-tie covariates must precede the monthly comparison. Regime and q estimation can wait for the λ=1 cells because the reference uses λ=0. The task generator can follow the audited reference prototype, but the reduced sweep design must still be fixed before examining sweep results.