# Min/max recheck (MM1)

*Reviewer: independent Claude review agent with no prior context (2026-11-22), run at the PI's request alongside GPT.*

22 November 2026. **MM1 is fixed. No new P1 or P2. GO for the fresh stage-1 rerun at this commit** (15,366 measures, `REUSE_POSITIONS=0`, fresh OUT), with one post-run check on the log described under "Release decision".

## Scope and what was executed

Followed Part A of `docs/PROMPT_gpt_recheck_minmax_stage2_rev5.md`. Reviewed only `git diff e8736b7 HEAD -- R/lib/sry.R R/02_build_positions.R tests/test_sry.R` (HEAD `a4f0a72`), plus the surrounding code needed to trace applicability (`R/02_build_positions.R:83, 97, 150–155, 181–196`) and the logging path in `run_all.sh:30, 84`. Line numbers below are one-based at HEAD.

Executed:
- `Rscript tests/test_sry.R` from the repo root. ALL PASS, 33 checks (34 lines ending in PASS, one of which is the "ALL PASS" summary line), matching the claim.
- A throwaway probe script in the scratchpad that sources `R/lib/sry.R` and calls `sry_combo()` on edge cases (results below).

Not executed: `run_all.sh`, `R/02_build_positions.R`, or anything that needs the real layers. The production census (0 of 384 combos with a no-layer dyad-year) is Claude's statement and was **not** reproduced. I reduce it below to two log lines that can be checked.

## MM1 status: fixed

**Evidence.** `R/lib/sry.R:77–79`. The min/max branch now takes `j <- ix[k[ix] > 0]` and computes the range and the rescaling only over `j`. Rows with `k == 0` are never written and keep their initial 0 from line 75. Line 80 asserts `all(E[k == 0] == 0)` alongside the existing bounds check, and line 81 attaches `attr(E, "n_no_layer") <- sum(k == 0)`. The clip branch at line 76 is unchanged, and it already gave 0 to `k == 0` rows (sc = 0, positive part 0).

The exact counterexample from `docs/REVIEW_minmax.md` (`M = [NA, −0.5, +0.5]`) now returns `(0, 0, 1)` with `n_no_layer = 1`, not `(0.5, 0, 1)`. The all-negative case `[NA, −0.5, −0.1]` returns `(0, 0, 1)`, so an imputed zero can no longer become the year maximum. Both behaviours are pinned by tests 9 (`tests/test_sry.R:86–92`).

The header comment (`R/lib/sry.R:68–70`) now says "range ≤ SRY_COMBO_TOL gets all zeros", which matches the strict `>` test at line 79 and closes the documentation nit from the earlier review.

**Equivalence when every dyad is applicable.** If every row of a year has `k > 0`, then `j == ix` and the new branch is arithmetically identical to the old one. So for any combo-year with no `k == 0` rows, min/max output is bit-for-bit unchanged. That is what makes the "changes no production measure" claim checkable by census alone.

**R02 census.** `R/02_build_positions.R:225–226` calls `sry_combo()` as before and emits a `message()` naming weighting, combo and count whenever `n_no_layer > 0`. The braces are correct: the `if (W %in% SRY) { ... }` block closes before the `else if` on line 227, and the parse is valid because it sits inside the function body of `one()`. `cent()` subsets `E[ei]` (line 185), which drops the attribute, so the attribute cannot leak downstream. The test edits at `tests/test_sry.R:81, 83` wrap results in `as.numeric()` only because the new attribute would otherwise break `identical()`/`all.equal()`. No other caller of `sry_combo()` exists in the repository.

## Probe results (scratchpad script, sourcing `R/lib/sry.R`)

| Input (one layer unless noted) | edge | Result | Comment |
|---|---|---|---|
| `[NA, −0.5, 0.5]` | minmax | `0 0 1`, n_no_layer 1 | MM1 case fixed |
| two layers, row 1 all NA, others partial | minmax | `0 0 0.25 1` | partial-NA rows still averaged over available layers |
| year 1 entirely no-layer, year 2 `[−1, 1]` | minmax | `0 0 0 1` | empty `j` handled (`length(j)` guard) |
| one applicable dyad (0.7) plus one NA | minmax / clip | `0 0` / `0 1` | see P3-c |
| applicable row cancelling to 0 among negatives | minmax | `1 0 0` | applicable zero is still ranked; distinct from no-layer |
| `matrix(numeric(0), 3, 0)` | minmax | `0 0 0`, n_no_layer 3 | zero layers handled |
| `NaN` in a row | minmax | treated as missing | `is.na(NaN)` is TRUE |
| `Inf` in a row | minmax | stops at line 80 | fails loudly |
| 100 random rows, 30% NA, 5 years | minmax | `E[k==0]` all 0; per-year applicable max 1, min 0 | invariant holds |

## Production census: what it rests on

I did not reproduce the census, but the code makes it easy to verify. For every layer other than UN, `V` has no NA (`R/02_build_positions.R:83`), so `sry_score()`'s applicability (`:151`, `year >= first_year & yr_has & !is.na(v)`) is a **year-level** condition for those layers. A dyad-year with `k == 0` beside dyads with `k > 0` in the same year can therefore arise only through UN missingness, in a year where every non-UN layer in the combo is year-inapplicable. Every combo contains IGO or defense (`:162`). The census claim therefore reduces to two facts R02 already logs:
1. `layer existence (first non-zero year)` (`:85`) shows IGO and defense at or before the first sample year;
2. `empty layer-years after existence` (`:98–99`) shows 0 for IGO and defense.

If both hold, `n_no_layer` is 0 for every combo, the R02 message never fires, and the min/max measures equal what the pre-fix code would have produced. If either fails, the fix still gives the correct missingness behaviour, and the message identifies the affected combos.

## New findings

No new P1 or P2.

P3 notes (no action required before the rerun):
- **P3-a.** The census uses `message()` rather than `logf()`, so lines carry no timestamp and, under `mclapply`, may interleave. They still reach the run log because `runR` captures `2>&1` (`run_all.sh:30`). But the `runR` grep pattern for this arm is `"Error|wrote"` (`run_all.sh:84`), so these lines go only to the log file and never appear on the console. Grep the log for them after the run.
- **P3-b.** The census is report-only, not a stop. That is acceptable now that the arithmetic is safe for `k == 0`. If the PI wants the "every combo has an applicable anchor" invariant enforced, a `stopifnot(attr(E, "n_no_layer") == 0)` under a flag would make it fail-closed.
- **P3-c.** This is a design consequence and was already implied by "the year's minimum applicable dyad gets 0". In a year with a single applicable dyad, or with all applicable scores equal, min/max gives that dyad edge 0 while clip would give a positive dyad edge 1. Min/max also always zeroes one real dyad per year (the minimum), which `cent()` then treats as no tie (`w > 0`, line 185). This is the PI's rule, not a bug. It is worth one sentence in the measure documentation.

## Release decision

**GO** for the fresh stage-1 rerun at `a4f0a72` (15,366 measures, `REUSE_POSITIONS=0`, fresh OUT). `R/lib/sry.R` changed, and both the arm signature (`run_all.sh:80`) and the global fingerprint hash it, so stale positions cannot be reused silently. Use a fresh OUT anyway, as planned.

Post-run check (no rerun is needed if it passes): confirm in the run log that (a) the IGO and defense first non-zero years are at the sample start and their empty layer-years are 0, and (b) no line matching `no applicable layer` appears. If (b) does appear, the census claim in the prompt is wrong for those combos. The measures are still correctly built, but the statement that "the fix changes no production measure" must be withdrawn for them.
