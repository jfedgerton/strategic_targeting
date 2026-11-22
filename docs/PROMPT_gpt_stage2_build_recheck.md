# Prompt for GPT: recheck of the stage-2 risk-set build (revision 2)

You audited revision 1 of the stage-2 build (`docs/REVIEW_stage2_build_gpt.md`, NO-GO with findings P1-1 to P3-7). Revision 2 answers each finding. No model has been fit by the pipeline.

Read:
1. `docs/stage2_build_notes.md`, which maps each finding to its change.
2. `docs/stage2_design_memo.md`, the new "B5 amendment" subsection in §13.
3. The diff areas in `R/lib/s2.R` (`s2_mk`, `s2_ties` coverage, `s2_spell_exits`, `s2_truncate_end`, `s2_side_prob`), `tests/test_s2.R` (tests 15, 22–27) and `R/s2_01_build_riskset.R` (rewritten).
4. `docs/s2_build_logs/` (both build logs, test output, `fatal_sets.csv`, `fatality_code_differences.csv`).

Please answer:
1. Is each finding resolved? Answer yes / partly / no, with file and line.
2. Do the new changes introduce bugs? Check especially:
   - the exit cap, and dispute end = latest exit;
   - strict side aggregation;
   - the prior-MID join (dispnum.x / dispnum.y);
   - the reference-cohort exclusion order.
3. Is the reference one-shot schema (`REF_COLS`) complete against memo §7 and the §13 reference cell? Is anything included that should not be?
4. Verdict: GO / GO with fixes / NO-GO for fitting the reference one-shot model. Give numbered findings (P1 / P2 / P3) with concrete fixes.
