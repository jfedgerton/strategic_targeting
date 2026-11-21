# Review: stage-1 (initiators-only) specification + stage-2 design memo

Context. The PI moved to a two-stage design (2026-11-21). Stage 1 = initiation decisions only: bilateral initiators vs coalition originators choose targets (conditional logit over system members). Joiners on either side are stage-2 actors (side choice), designed in `docs/stage2_design_memo.md` and not yet implemented. The full multiverse loop at commit f60c55c is running on the PI's laptop with `TYPES=bcd`; stage 1 will run afterwards with `TYPES=init OUT=output_init bash run_all.sh`.

## Part 1 — verify the stage-1 code (narrow; do not re-review unchanged code)
| Change | Where |
|---|---|
| `--variant=initiators`: from a four-type (bcjs) frame keep bilateral and coalition-originator strata only; roles bil/coal; types tag "init" | `R/06_make_variants.R` (search "PI 2026-11-21") |
| `dump()` builds the init frame from the bcjs frame (also for the unit-weights frame), validates it, writes atomically | `run_all.sh` `dump()` |
| `TYPES` overridable from the environment | `config.sh` |
| Influence ranks on `coal-bil:u_R_less` under init; joiner-strata side drops skipped under init (whole-dispute side drops kept) | `run_all.sh` §6–7 |

Executed here (TEST=1, Linux): LOOP_DONE; restart 19 s, no refits. Initiator frame at the 1-month window keeps 478 bilateral and 141 originator strata (drops 242 attack-side and 219 target-side joiner strata). Main measure: coal−bil less central 1.07 [0.56, 1.58], stronger 1.01 [0.18, 1.84], surface p 4e−6. Step B risk set unchanged (originators at onset).

Questions: (1) any bug or P1/P2 in these changes; (2) is anything in the existing grid invalid or meaningless for a two-group initiator frame (e.g., permutation of type labels, event study by type, typed controls); (3) does dropping joiner strata change the interpretation of any robustness check.

## Part 2 — review the stage-2 design memo (`docs/stage2_design_memo.md`)
Judge the unit, risk set, timing rule, model, covariates and identification table. For each of the five open decisions in §10, give a recommendation with reasons. Flag anything that would make the stage-2 estimand ill-defined or not comparable with stage 1.

Write everything to `docs/REVIEW_stage1.md`. Read-only otherwise; do not run the R scripts.
