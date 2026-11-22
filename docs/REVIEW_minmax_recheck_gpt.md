# MM1 recheck

22 November 2026. **MM1 is fixed. GO for the fresh 15,366-measure stage-1 rerun, with `REUSE_POSITIONS=0` and fresh, separate output folders for bcd/init. No new P1/P2 was found in the requested code delta.**

## Scope and evidence

Followed Part A of `docs/PROMPT_gpt_recheck_minmax_stage2_rev5.md`. Compared `e8736b7af9958493558da2243b4ed0d08a571ca3` with the supplied HEAD for `R/lib/sry.R`, `R/02_build_positions.R` and `tests/test_sry.R`; unchanged scoring and estimators were not re-reviewed.

Package: `/Users/jarededgerton/Downloads/targeting_multiverse_recheck_mm1_rev5.zip`; HEAD `a4f0a72c353739f14d6b7f5fc40aae96ba9e43ff`; SHA-256 `140b76eaac5a41b171577b1bddf505a46b1d9eb1ada6eeda44cca9b4865d7ec8`. References are one-based archive-relative lines. No R code, project tests or runner were executed. Only the two requested review reports were written. Independent Python arithmetic checked the relevant edge cases. The reported 33 R checks and production applicability census are supplied evidence, not independently reproduced results.

## Verification

**`R/lib/sry.R:73–81`: fixed.** The applicable-layer count is preserved. The min/max branch now builds its per-year extrema from `j = ix[k[ix] > 0]`, assigns normalized weights only to j, and leaves all other entries at initialized zero. It therefore fixes both parts of MM1: no-layer rows cannot acquire edges, and their artificial zero cannot alter the extrema of an all-negative applicable set. Empty applicable sets, single applicable values and ranges at or below tolerance return zeros without taking an invalid range or dividing by zero. The bounds/no-layer-zero assertion checks the invariant directly.

Independent arithmetic gives:

| Signed scores, NA meaning no applicable layer | Min/max edges |
|---|---|
| NA, −0.5, +0.5 | 0, 0, 1 |
| NA, −0.5, −0.1 | 0, 0, 1 |
| NA, NA | 0, 0 |
| NA, 0.2 | 0, 0 |
| 0.2, 0.2 + 5e−11 | 0, 0 |

The clipping branch at 76 is numerically unchanged. Adding `n_no_layer` changes the returned vector's attributes, not its values; the test adjustments using `as.numeric()` appropriately distinguish those two facts.

**Census integration: fixed.** The helper returns `attr(E, "n_no_layer")` at 81. `R/02_build_positions.R:225–226` reports positive counts with weighting and combination identifiers before centrality construction. The attribute is not substituted for edge weights. No spec-ID, family or output-file change is introduced by this patch.

**Regression coverage: discriminating.** `tests/test_sry.R:87–92` adds the mixed-sign all-NA row, all-negative applicable set, clipped-zero and small-range cases. The first two fail under the previous implementation. The mixed-sign test additionally checks the census attribute. An all-NA-year/single-applicable-value test would be a small optional extension; the inspected length/range guards already handle both, so its absence is not a new blocker.

## Production implications

The supplied census says IGO and defense are applicable everywhere in the real input panel and every allowed combination contains one of them. Under that stated invariant, k is always positive and the revised min/max branch is mathematically identical to the preceding branch. This supports the claim that the fix guards the API without changing production positions. It does not replace an independently executed equality check, which this review has not performed.

No remaining MM1 release condition is open. Preserve the existing code/helper hashes in arm stamps, use fresh OUT folders rather than attempting to continue an older fingerprint, and retain the build/test/census logs with the rerun. This GO concerns the narrow reviewed measurement fix, not the separate stage-2 memo.
