# Recheck (2026-11-22): MM1 and stage-2 R1–R3

## Part A — MM1 (write `docs/REVIEW_minmax_recheck.md`)
Your `docs/REVIEW_minmax.md` held the min/max release for MM1. Verify only `git diff e8736b7 HEAD -- R/lib/sry.R R/02_build_positions.R tests/test_sry.R`.
- Fix: min/max normalizes over dyads with >= 1 applicable layer only; dyads with no applicable layer keep edge 0 under both rules; asserted; `attr(E, "n_no_layer")` census logged by R02 for any combo where it is > 0.
- Tests 9 (all-NA row in a mixed-sign year → 0; all applicable scores negative; clip unchanged; range below tolerance → zeros): 33 checks pass.
- Production census (Claude, real 1816–2014 layers): IGO and defense are applicable in every dyad-year (0 non-applicable rows each), and every combo contains one of them, so 0 of 384 combos have any dyad-year with no applicable layer. The fix therefore changes no production measure; it guards the API.
Go / no-go for the fresh stage-1 rerun at this commit (15,366 measures, `REUSE_POSITIONS=0`, fresh OUT).

## Part B — stage-2 memo revision 5 (write `docs/REVIEW_stage2_memo_rev5.md`)
Your `docs/REVIEW_stage2_memo_rev4.md` asked for R1–R3. Read §13 subsections R1, R2, R3 and the §10 supersession note. For each: fixed / partly / not, with the memo line; any remaining ambiguity that lets the implementation choose the estimand. Go / no-go for freezing the analysis-ready risk set and legality table (infrastructure is already GO).

Read-only otherwise; do not run the R scripts.
