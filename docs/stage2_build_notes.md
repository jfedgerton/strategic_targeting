# Stage-2 risk-set build: notes for the audit (2026-11-22)

This builds the infrastructure for stage 2 against `docs/stage2_design_memo.md` §13 (rev 6). **No model has been fit.** The files are `R/lib/s2.R` (library), `tests/test_s2.R` (fixtures, ALL PASS) and `R/s2_01_build_riskset.R` (real-data builder). Logs are in `docs/s2_build_logs/`.

## What is built
| Piece | Memo rule | Implementation |
|---|---|---|
| Baseline roles | R1 | `orig == 1`; baselines enter at month 0 |
| Window promotion | N2 | Participant-level. Initiating side, not baseline, W ≥ 1, definite bound `entry_hi − onset_lo ≤ W`. Onset = earliest originator entry (stage-1 `TIMING_REF = originator`) |
| Masks | R1, N1 | `s2_mask(m, W)` depends only on month index and window. W = 0 gives no originator alternative. For W ≥ 1, A is masked and orig is available in months 0…W |
| Initiating side | §4 | Stage-1 R05 rules 1–2 plus the case review (stops on disagreement). Rule 5 (war link from mixed-source rows only) is dropped and listed; it does not occur |
| Cohort, horizon | B1, B2 | System members at the onset month, excluding baselines. H = min(dispute end, Dec 2014, system exit); a `complete` flag is set. Dispute and spell ends come from GML participant end months |
| One-shot, monthly rows | B1, B2 | Four categories. Monthly rows run to first entry or censoring. Masked observed outcomes stop the build |
| Rosters | B3 | Month 0 = baselines by side. From month 1, spells active at m − 1, first-spell side. Empty sides are carried forward and flagged; the "end" policy is also counted |
| Ties | §4, N4 | Main measure (W2contig, 8 layers) edges rebuilt with R/02's code path. Raw sum and max over the roster at lag 1 (one-shot also lag 5). The edges reproduce stage-1 log PageRank to 1e−13 |
| Tiers | §4 | Frozen at onset: contiguity ≤ 400 mi, same COW region, major power, any ATOP pact at t − 1 with a baseline participant |
| Folds | B4 | Six era folds by onset year |
| Jointly fatal | B5, R3 | Side level, from GML dyadic participant codes. yes / no / unresolved bounds |
| Incidence | B2, R2 | Flat CIF and sequential CIF with h_O (fixture-tested) |
| Legality | C1, R3 | `s2_cell_key()` returns illegal cells as NA and a duplicate key otherwise |

## Acceptance tests vs stage 1 (bcjs dump, W = 1)
| Test | GML-fatal sample | Stage-1-fatal sample |
|---|---|---|
| Events covered | 691 of 853 | 853 of 853 |
| (a) baseline roles | 0 mismatches | 0 mismatches |
| (b) promotions | 0 of 271 (51 = 51) | 0 of 325 (58 = 58) |
| (c) joiner sides | 0 of 220 | 0 of 267 |
| Excluded projection | 1 event from a second spell (USSR in 258, 1941; re-entries are not stage-2 entries) | same |

The other windows (0, 2–6) need the matching dumps from the Mac run. Only the W = 1 dump is in the sandbox.

## Counts (reference W = 1)
| | GML-fatal | Stage-1-fatal |
|---|---|---|
| Disputes | 541 | 560 |
| Cohort rows | 63,986 | 65,416 |
| One-shot A / B / orig | 84 / 168 / 68 | 118 / 189 / 71 |
| Monthly rows | 783,167 | 810,791 |
| Jointly fatal yes / no / unresolved | 312 / 204 / 25 | 318 / 177 / 65 |

## Decision needed (PI): which fatal-dispute set is the reference?
Memo B5 says "fatal = GML dispute fatality > 0 (the stage-1 `disp_fatal` source)". That premise is wrong. Stage 1 samples on `events.rds in_fmid`, which uses a COW-based source chain with MIC and war-rule fills. The two sets overlap on only 440 disputes: 102 are GML-only and 120 are stage-1-only.
- **Recommended:** the stage-1 set as the reference. It matches stage 1's population, and the acceptance tests then cover all 853 events.
- **Sweep level:** GML as a labelled alternative.
- **Cost:** with the stage-1 set, 33 disputes have no GML end dates and fall back to December of the last spell year (counted).

## Not built yet (listed so the audit can judge scope)
- Late-entry landmark cells (B2).
- Recurrent / re-entry family (R3).
- Dynamic-entrant monthly rows: entrants are counted (5) but not expanded.
- Monthly-updated tiers.
- Monthly non-tie covariates. CINC, distance and alliance exist only for the one-shot at onset; monthly rows carry ties and the frozen tier.
- Prior-MID covariate and regime (a q feature).
- Two-month aggregation.
- Side-switch censoring variant (8 switchers).
- The q model and task generator. Folds are assigned but nothing is fit.

## Known caveats
- CINC ratios are heavy-tailed (max ≈ 13,000). Logs are needed at model time. Side CINC that is zero or missing gives NA, never 0/0 (A 753 rows, B 1,044, GML sample).
- Target centrality is missing at t − 1 for 11 disputes (12 in the stage-1 sample), where a target-set member has no t − 1 position.
- 10–12 recorded entries fall after the state's horizon (system exit) and are treated as out.
