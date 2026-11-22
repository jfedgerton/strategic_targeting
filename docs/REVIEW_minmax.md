# Min/max edge-option review

22 November 2026. **The requested arithmetic, weighting IDs and runner integration check out. One new P2 needs a decision or guard before releasing the expanded grid: min/max normalization can give positive edges to dyads with no applicable layer.** No new P1 was demonstrated on the supplied real inputs. Hold the min/max production release until the missingness policy is enforced or its problematic condition is ruled out by an explicit invariant/census. The clipped rule is unchanged.

## Scope and evidence

Followed Part A of `docs/PROMPT_gpt_review_minmax_and_stage2_rev4.md`. Inspected the exact delta from `a0b4cc9428bc00a637ef0bf86eaba09e21caf144` in `R/lib/sry.R`, `R/02_build_positions.R`, `run_all.sh` and `tests/test_sry.R`, plus the small existing validation/subset/summary interfaces needed to trace the new IDs. Unchanged layer fitting and estimators were not re-audited.

Package: `/Users/jarededgerton/Downloads/targeting_multiverse_review_minmax_stage2.zip`; HEAD `e8736b7af9958493558da2243b4ed0d08a571ca3`; SHA-256 `dbd3e2b75d4494f2ffa31fe536e635fb923c29f2bd06278d3245edd0d7f1e104`. Line references are archive-relative and one-based. No R code, project test or runner was executed. Native Bash syntax checks passed. Independent Python arithmetic and grid counts were used. Only the two requested review files were written.

## What is correct

- **Normalization:** `R/lib/sry.R:70–77` validates the edge option, computes the available-layer signed mean, applies the established cancellation tolerance first, and uses `(score−min)/(max−min)` separately by year. A range at or below 1e−10 yields zeros. Bounds are asserted. The clip branch at 75 preserves the preceding arithmetic and default behavior.
- **Same tie models:** R02:91–94 recognizes the two new names only in the core arm; line 153 removes only the trailing `mm` before calling the existing W6/W7 layer scorer. Line 225 selects the min/max combination branch from the suffix. Residual models and families have not inadvertently changed.
- **IDs and grid:** R02:233–241 includes the full weighting token in IDs, metadata and output filenames. Two new weightings add `2 × 384 × 2 = 1,536` all-PTA specs and `2 × 192 × 2 = 768` each for inforce/any: **3,072 additions, 15,366 total**. The six new position files bring the expected file count to 29. The existing combination constraint requires IGO or defense, and PTA-containing combinations only for any/inforce; the count respects those rules.
- **Runner:** `run_all.sh:62–65, 94–99` includes all four SRY weightings for all three PTA arms, checks their files, and hashes their contents. Existing arm signatures include R02 and the SRY helper; the global fingerprint includes library code. This change therefore invalidates old local-arm signatures and cannot silently reuse an old two-weighting arm as the new complete arm. Use fresh OUT folders, as requested.
- **Downstream:** `R/03_validate_measures.R:67–75` has no W6/W7 whitelist, and R04 selects by spec IDs/layer tokens. `R/11_summarize.R:34–37` parses `W6srybookmm` and `W7sryyearmm` as distinct weighting families. Family-balanced percentages consequently give those new weighting cells their own mass; the expanded summary is a new grid composition, not the same descriptive denominator as before.

The supplied exact-match evidence for PTA-any is useful but was not reproduced here. Report completion and fit diagnostics for all/inforce when available. The standalone book-combination plausibility failure does not itself establish a code error; retain its failed validation status rather than changing the plausibility rule to admit it. Full-grid results and plausible-measure subsets should preserve that distinction.

## MM1 — P2: no-applicable-layer dyads can acquire positive weights

**Location:** `R/lib/sry.R:72, 76`; applicability is produced by `R/02_build_positions.R:151` and `sry_layer()`.

Line 72 substitutes score zero when a dyad has no applicable selected layer (`k=0`). Under clipping that remains edge zero. Under min/max, if the year's observed score minimum is negative, the substituted zero becomes positive.

Concrete helper input for one year, one selected layer:

```text
M = [ NA, −0.5, +0.5 ]
k = [  0,    1,    1 ]
sc = [ 0, −0.5, +0.5 ]
minmax edge = [ 0.5, 0, 1 ]
```

The first dyad has no measured/applicable layer but receives weight 0.5. This is distinct from the explicitly requested behavior of giving an *observed zero or negative residual* a nonzero relative edge. It can alter centrality and also the extrema: an artificial zero can become the maximum when every applicable score is negative.

This is a demonstrated arithmetic path, **not proof that a current production row triggers it**. The caller requires each combination to contain IGO or defense. If at least one such anchor is applicable for every dyad in every relevant year, that invariant prevents the problem on those inputs. A valid combination such as defense+UN can nevertheless expose the path in a year where defense is globally inapplicable and UN is missing for some dyads. I have not executed a production applicability census and do not claim such a year exists in the supplied current data.

**Minimal resolution:** preserve the distinction between applicable zero and no applicable layer. Normalize over applicable dyads and leave `k=0` edges zero, if that is the intended missingness policy. If the PI deliberately wants the zero-imputed score transformed too, state that as a separate imputation choice, with a sensitivity and tests. Alternatively establish and enforce `k>0` wherever a nonflat year is normalized, so this API path cannot occur in the production grid. A default fallback silently creating links is insufficient documentation of the choice.

Add an all-NA-row test with a mixed-sign year and another with all applicable scores negative. Current tests 8 cover partial NA rows with another observed layer, not the `k=0` case. Their arithmetic, flat-year, separate-year and unknown-option assertions otherwise discriminate the new branch appropriately. The equality-at-tolerance behavior is simply inclusive suppression; document it as “range ≤ tolerance.”

## Release decision

**Hold the expanded min/max rerun for MM1's policy/invariant and targeted tests.** This is a small, isolated issue; no rewrite of the scoring models or grid is needed. If a checked production invariant proves it unreachable, a fresh-grid run can proceed on those inputs with that guard documented. The supplied PTA-any reproduction and pending all/inforce runs do not, by themselves, certify how no-applicable-layer dyads were treated. Previous SRY safeguards and the deliberately fail-closed nonconstant all-FE-removal policy are outside this delta and are not reopened.
