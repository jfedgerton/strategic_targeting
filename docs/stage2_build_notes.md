# Stage-2 risk-set build: notes for the audit (revision 2, 2026-11-22)

Infrastructure for stage 2 against `docs/stage2_design_memo.md` §13 (rev 6) and the B5 amendment. **No model has been fit by the pipeline.** Files are `R/lib/s2.R` (library), `tests/test_s2.R` (fixtures, ALL PASS) and `R/s2_01_build_riskset.R` (real-data builder). Logs are in `docs/s2_build_logs/`. Revision 2 answers `docs/REVIEW_stage2_build_gpt.md`.

## Responses to the GPT build review
| Finding | Change |
|---|---|
| P1-1 spell gaps | `s2_spell_exits()` caps each exit at December of the spell's last year and flags capped, imputed and multi-episode matches. A test checks that a 1940–42 GML episode leaves the participant absent in 1941. The builder also asserts that no exit reaches into the next spell, and audits all 7 second spells. Production: 56 exits capped, 116 imputed, 10 multi-episode |
| P1-2 reference covariates | `--positions`/`--poscol` are required. Added own centrality, prior MIDs per side, contiguity per side, side size, major power per side and log CINC ratios. `REF_COLS` defines the one-shot reference schema. `in_ref` and `excl_reason` define the estimation cohort, and an assertion confirms it has no NA |
| P2-3 fail-closed | `--mode=release` (the default) requires the dump and stops on any mismatch, missing or ambiguous spell, or unclassified record. Every projection record gets one reason: outside selected disputes, missing spell, ambiguous spell, later spell, or matched. `--mode=exploratory` logs and continues |
| P2-4 missing vs zero | Ties take a coverage set: an uncovered year gives NA plus `tie_nocov` (21 lag-5 rows). Side aggregates are strict (NA when any member is missing, flags kept). Distance, contiguity and alliance per side are measured at t − 1. Tiers stay frozen at onset-year contiguity (documented) |
| P2-5 population | PI decision is the stage-1 set (default). ID lists are in `fatal_sets.csv`. The 292 fatality-code differences break down as 75 fatal/nonfatal status, 33 severity level, 184 unknown code on one side (`fatality_code_differences.csv`). Imputed-end flag, metadata on every output, and `s2_mk()` rejects invalid month codes. B5 amended in the memo |
| P2-6 acceptance | Added: baseline equality on the full universe (1,174), promotion-boundary check (2,765 participant-windows), uncertain-date counts per window, all-year PageRank check (199 years, max deviation 5.7e−13). New tests: January/December keys, inclusive horizon, both-sides-fatal with unknowns, a 1914–19 fold, end-policy truncation, side-equation mask |
| P3-7 | `s2_side_prob()` returns exactly 0 for A where `avail_A` is FALSE (tested) |

Edges cannot be compared element by element: R/02 does not write edges. PageRank agreement in every year is the available check. Windows 0 and 2–6 still need the Mac dumps.

## Acceptance (W = 1)
| | Stage-1 set (reference) | GML set |
|---|---|---|
| Projection records (initiator × dispute × event year × type) | 870 | 870 |
| Matched / later spell / outside sample | 869 / 1 / 0 | 707 / 1 / 162 |
| Role, promotion and side mismatches | 0 / 0 of 341 / 0 of 283 | 0 / 0 of 287 / 0 of 236 |

## Counts (reference build, stage-1 set, W = 1)
| Quantity | Value |
|---|---|
| Disputes | 560 (end imputed for 47) |
| Cohort rows | 65,416, complete horizon 65,274 |
| One-shot A / B / orig / out | 123 / 196 / 71 / 65,026 |
| Reference estimation cohort | 63,143 rows (A 121, B 184, orig 71). Exclusions: target position 1,262, own CINC / ratio 455, own centrality 414, incomplete horizon 142 |
| Monthly rows | 806,131 |
| Carried side-months | 68 in 12 disputes. The end policy would drop 8,992 monthly rows |

## Still not built (GPT: these can wait for their own models or sweeps)
Late landmarks, the recurrent family, dynamic entrants, monthly-updated tiers, two-month aggregation, switch-censoring, the q model and regime, and the task generator. Monthly non-tie covariates must come before the monthly comparison.

## Exploratory (outside the pipeline, PI request)
`explore/x_s2_fe_twostep.R` is a two-step fixed-effects stand-in for the AME's shared target effect. It is labelled pre-freeze and uses the revision-1 build.
