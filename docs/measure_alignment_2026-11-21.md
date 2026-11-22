# Alignment of the targeting measures with the book's "sry" setup (2026-11-21/22)

PI request: check whether the multiverse measures are aligned with the book project's measurement (expected-tie model per layer, observed minus expected), and add the book-aligned measures to the multiverse; the choice of main measure is decided later.

## Standalone test (2026-11-21, outcome-blind, then main-spec fits)
Built outside the loop from the loop's `layers_dy.rds`, main layer set igo+def+off+neu+nonagg+consul+dca+un.
- sry-book: v ~ log distance | sender + receiver + year FE (the book). sry-year: v ~ log distance + contiguity | year FE (no state FE).
- Logit (binary layers), Poisson (IGO counts), OLS (UN kappa). Score = observed − expected, scaled by the layer-year max |score|, averaged. Centrality: PageRank on the positive part; PN centrality on the signed network.
- Assumption: the book's cross-layer scaling is not documented; the layer-year max scaling mirrors the targeting code. The book's own DCA layer differs.

| Check | Main (W2contig) | sry-book PR | sry-year PR | sry-book PN | sry-year PN |
|---|---|---|---|---|---|
| Alignment with main (mean within-year Spearman) | — | 0.84 | 0.91 | 0.58 | 0.79 |
| Passes plausibility | yes | yes | yes | no | no |
| USSR Cold War percentile | 0.70 | 0.91 | 0.91 | 0.19 | 0.38 |
| R² of position on IGO membership count, 1946–89 | 0.58 | 0.29 | 0.34 | 0.10 | 0.28 |

Across the 4,142 plausible multiverse measures, the median alignment is 0.65 (sry-book PR) and 0.71 (sry-year PR); the main measure is near the maximum (0.86 / 0.94). sry-book's post-1945 top states (e.g. Comoros, Libya in 2010) show that state FE remove tie-proneness. PN fails the plausibility rules and is not added.

Main-spec stage-1 fits (init, fatal, 1-month window), coalition − bilateral, less central: W2contig 1.07*; sry-year 1.10*; sry-book 0.76*; adding PTA (count or any) moves each by ≤ 0.1.

## Added to the multiverse (PI 2026-11-22)
- Weightings W6srybook and W7sryyear (core arm, all 384 layer combinations), PageRank and degree, under PTA all / in force / any.
- Unweighted PTA (`--pta=any`: 1 if any PTA, signed or in force) as a third PTA coding for W1–W4 too.
- Diplomacy arm unchanged (no tie model specified for the diplomacy layer). Main measure unchanged until the PI decides.
- New measures: W1–W4 × any 1,536; sry × all 1,536; sry × in force 768; sry × any 768. Total grid 7,686 → 12,294.

## Fit-validity policy (GPT reviews of 2026-11-22, `docs/REVIEW_sry.md`, `docs/REVIEW_sry_recheck.md`)
- Constant outcome over a layer's applicable rows: score 0, applicable, no fit.
- Rows removed as perfectly predicted FE groups (or Poisson zero-outcome separation) inside a fitted layer: residual 0.
- A nonconstant layer whose every row is removed by recursive perfect-fit FE removal is refused (the arm stops); it does not occur in the 1816–2014 layers (all sry arms built). `tests/test_sry.R` check 7.
- Any fixest warning, non-convergence, non-finite fitted value, missing required predictor or unexplained omission stops the arm.
- Run with `REUSE_POSITIONS=0` (enforced) in fresh output folders. GPT recheck: GO at 2f4e903.
