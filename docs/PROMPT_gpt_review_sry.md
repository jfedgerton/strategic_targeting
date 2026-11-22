# Review: book-aligned measures added to the multiverse (2026-11-22)

You last reviewed this repository at 230d585 (stage-1 recheck, GO). The PI asked to add the book project's "sry" measurement to the multiverse and to decide the main measure later. Background and the standalone validation: `docs/measure_alignment_2026-11-21.md`. Review only the change `git diff 230d585 HEAD -- R/02_build_positions.R run_all.sh R/06_make_variants.R README.md`; the stage-2 memo revision in the same range is a separate, later review.

## What changed
| Change | Where |
|---|---|
| `--pta=any`: unweighted PTA (1 if the dyad has any PTA, signed or in force; not logged); only the 192 combos with PTA are built, as for `inforce` | `R/02` V construction, combos block |
| Weightings `W6srybook` (v ~ log distance \| sender + receiver + year FE) and `W7sryyear` (v ~ log distance + contiguity \| year FE), requested with `--weightings`, core arm only. Logit for binary layers (ATOP types, DCA, PTA any), Poisson on raw counts (IGO, PTA all / in force), OLS for UN kappa; doubled data; observed − expected averaged over orientations; scaled by layer-year max \|score\|; combo = mean over applicable layers; edge = positive part rescaled within year; PageRank and degree via the existing `cent()` | `R/02` `sry_score()`, combo loop |
| Runner: the reused 2026-11-18/19 files are linked by explicit name; four new arms (`core_any`, `core_sry_all`, `core_sry_inforce`, `core_sry_any`) are always built in the run folder with the existing stamped `arm()`; EXPECT gains 10 files; fixest always required | `run_all.sh` §0–1 |
| P3 from your recheck: dropside log wording for the whole-dispute modes | `R/06` dropside block |

## Checks run by Claude
- The loop builder reproduces the standalone sry build for the main layer set: max abs difference 1.8e−15 (W6srybook and W7sryyear, PageRank). `--pta=any` W2contig on main + PTA equals an independent replica exactly (0).
- Plausibility on the new files (R/03, unchanged): W6srybook ptaall 90% (pr) / 88% (deg); W7sryyear ptaall 89% / 78%; ptaany arms 5–76%. 3,840 new measures checked here (the in-force sry arm was not built in this test).
- A standalone replica of the existing W2contig rule reproduces the loop's main measure to 6e−13.

## What I want (write to `docs/REVIEW_sry.md`)
1. Is `sry_score()` a correct implementation of the stated model (family per layer, perfectly predicted FE groups, symmetry, existence years, UN missingness, scaling, averaging over applicable layers, positive part)? Any case where a layer-year with no variation or a fully separated FE group produces a spurious tie?
2. Does `--pta=any` interact correctly with the existing W1–W4 paths and with spec ids / families downstream (`R/03`, `R/04 --set=noatop`, `R/07`, `R/11` family parsing)?
3. Runner: can the new arms mix with linked files from another layer version, or skip a rebuild after a builder change? Does the fingerprint cover the new arms?
4. Any NEW P1/P2. Go / no-go for a fresh production rerun (`OUT=<new folder>`, both `TYPES=bcd` and `TYPES=init`).

Read-only otherwise; do not run the R scripts.
