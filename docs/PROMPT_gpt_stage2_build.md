# Prompt for GPT: audit of the stage-2 risk-set build

You are auditing code for a political science paper on how states choose targets in fatal militarized disputes (1816–2014). Stage 1 (target choice) is done. This audit covers the **stage-2 infrastructure**: which third parties enter an ongoing fatal dispute, and on which side. No model has been fit.

Read in this order:
1. `docs/stage2_design_memo.md` §13 (the contract; precedence N > R > B–C).
2. `docs/stage2_build_notes.md` (what was built, acceptance results, open decision, scope).
3. `R/lib/s2.R`, `tests/test_s2.R`, `R/s2_01_build_riskset.R`.
4. `docs/s2_build_logs/` (build logs for both fatal-dispute definitions, test output).

Context for reading the code: `R/05_dump_events.R` lines ~120–360 hold the stage-1 promotion and initiating-side rules that stage 2 must reproduce.

Please answer:
1. **Contract fidelity.** Does each implemented rule match §13 (R1, N1 prospective masks, N2 participant-level promotion, B1–B3, B5/R3 bounds, N4 raw ties)? Cite lines.
2. **Fixtures.** Does every fixture in the memo's "Fixtures required before the freeze" list have an exact assertion? Which are missing or weak? Note that the recurrent-return fixture is deliberately absent because the recurrent family is not built.
3. **Acceptance tests.** Are the projections sound? This covers the second-spell exclusion and coverage of 691 vs 853 events. What else should be tested before the freeze, for example windows 0 and 2–6 once the Mac dumps exist?
4. **Fatal-dispute definition.** The memo's B5 premise was wrong: the GML and stage-1 sets overlap on only 440 of 542 / 560 disputes. Is the recommendation sound (stage-1 set as reference, GML as a sweep level)? What should be reported?
5. **Bugs or silent failures.** Look at joins (data.table non-equi joins), month-index conventions, horizon and censoring, roster carry-forward, ties at month 0 vs m − 1, tier construction, and covariate NAs.
6. **Scope.** Which of the "not built" items must exist before the reference model is fit, and which can wait for the sweep?

Return: a verdict (GO / GO with fixes / NO-GO for building the reference model), then numbered findings (P1 must fix, P2 should fix, P3 note), each with the file and line and a concrete fix.
