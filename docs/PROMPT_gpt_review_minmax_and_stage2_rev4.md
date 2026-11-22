# Two reviews in one pass (2026-11-22)

## Part A — stage-1 code: min-max edge option (write `docs/REVIEW_minmax.md`)
You gave GO for the 12,294-measure grid at 2f4e903 (`docs/REVIEW_sry_recheck.md`). The PI then asked for a second edge rule for the book-aligned weightings: normalize the signed combo score to [0,1] within year instead of clipping at 0. Review only `git diff a0b4cc9 HEAD -- R/lib/sry.R R/02_build_positions.R run_all.sh tests/test_sry.R`.

| Change | Where |
|---|---|
| `sry_combo(M, yr_split, edge = "clip" | "minmax")`; minmax = (score − year min) / (year max − year min), zeros when the range is below 1e−10; the SRY1 cancellation tolerance applies first; bounds asserted | `R/lib/sry.R` |
| New weightings `W6srybookmm`, `W7sryyearmm`: same layer scores as W6/W7 (`sub("mm$", "", W)` selects the tie model), minmax edges | `R/02_build_positions.R` |
| The three sry arms now build four weightings each; EXPECT gains six files; grid 12,294 → 15,366 | `run_all.sh` |
| Tests 8 (arithmetic, clip unchanged, flat year, per-year normalization, unknown option stops): 29 checks pass | `tests/test_sry.R` |

Evidence (Claude): standalone min-max test on the book layer set in `docs/measure_alignment_2026-11-21.md` (fails plausibility: USSR Cold War percentile 0.38–0.39). Real-data loop build, `--pta=any` with all four sry weightings: completes; the loop's `W6srybookmm`/`W7sryyearmm` book-layer positions equal the standalone min-max build exactly (max abs diff 0), and the clipped W6/W7 files are identical to the reviewed build (0). The `all` and `inforce` arms are being rebuilt in Claude's 7 GB sandbox (an earlier attempt ran two builds concurrently and was killed for memory; not a code fault); results will be reported with your review. Questions: correctness of the rule and its interaction with spec ids, families (R/11 parses `W6srybookmm` as its own weighting), validation and the fingerprint/arm stamps; any NEW P1/P2; go / no-go for the fresh rerun at this commit.

## Part B — stage-2 memo revision 4 (write `docs/REVIEW_stage2_memo_rev4.md`)
You gave NO-GO on revision 3 (`docs/REVIEW_stage2_memo.md`, B1–B5, C1, C2, AME). Read `docs/stage2_design_memo.md` §13 (the builder contract) and the corrected `docs/ame_joint_sketch.md`. For B1–B5, C1, C2 and the AME language: fixed / partly / not, with the memo line and any remaining ambiguity that would let the implementation choose the estimand. Go / no-go for the stage-2 build (risk set, rosters, ties, folds, legality table, synthetic fixtures: early A promotion, early B entrant, no pre-month roster, empty later side, system exit, administrative censoring, 1914 onset, multilateral dispute with mixed jointly-fatal pairs).

Read-only otherwise; do not run the R scripts.
