# Second pass: verify the fixes to your multiverse-loop review

You reviewed this repository on 2026-11-19 (`docs/REVIEW_loop.md` if present; findings P1-1 to P3-14 and sections B, C). Every finding has been addressed except the two items the PI deferred to the AME phase (a genuine joiner choice set; an interacting-dyads sensitivity). Verify the fixes. Do not re-review unchanged code, and do not rewrite the project.

Read `README.md`, `docs/decisions.md`, `docs/historical_exclusions_DRAFT.md`, `run_all.sh`, then the changed scripts.

## What changed, by finding
| Finding | Fix | Where |
|---|---|---|
| P1-1 windows 0–3 only | Domain 0–6; preflight checks every configured window. The Germany→Poland regression check now applies only to windows ≤ 3 (with longer windows a later same-side entrant can legitimately fall inside; logged) | `R/05_dump_events.R` (search "loop review P1-1"), `run_all.sh` §0 |
| P1-2 test pollutes production | `TEST=1` writes to `output_test/` | `run_all.sh` header |
| P1-3 markers/provenance | Fingerprint of commit + config + data_entry manifest + positions sha256; one fingerprint per output folder; production refuses a dirty tree; symlinks checked against config; atomic writes in every script; a sweep is finished only with an atomically written CSV and no ERROR rows, else listed in `status/incomplete.txt` and the run exits 2; the summarizer reads only manifest-listed outputs | `run_all.sh`, `R/07`, `R/03`, `R/04`, `R/06`, `R/08`, `R/09`, `R/10` |
| P1-4 collapsed side flag | Initiator side taken from the participant table (earliest spell, primary dispute); war link for rule 2 from single-source rows when available; count of changed events logged (17 of 3,424 events, 9 mixed-source; all 13 case-review checks still pass) | `R/05_dump_events.R` (search "loop review 2026-11-19 P1-4") |
| P2-5 step B | Unit = events whose initiator was an originator before any window promotion (decisions at dispute onset); outcome = GML dispute fatality > 0 (unknown dropped and counted). Joiners excluded because first-fatality dates are not in the data. Old membership version kept as `--stepB_unit=all --stepB_outcome=pair` for comparison | `R/07_fit_sweep.R` `run_stepB`, `R/05` (`disp_fatal`, `orig_joiner`) |
| P2-6 intercepts | Type-intercept contrasts and `all_differ` (gap terms + intercepts jointly) | same |
| P2-7 event study | `--align=onset|entry`; balanced, unbalanced and pre-only runs; excluded events written; labels say within-year z, not ties | `R/10_event_study.R`, `run_all.sh` §9 |
| P2-8 influence | Sweep validity policy; every unit recorded with status; `--q` validated; top file written before the LOO file (atomic pair); runner checks both | `R/09_influence.R` |
| P2-9 resampling | Per-model validity per draw (M0 kept when M1 fails, incl. the full-sample fit); t−1 eligibility as the sweep; draw ids with reasons; `perm_type` refused on typed-control dumps; per-model thrown errors caught in the sweep | `R/08_resample.R`, `R/07_fit_sweep.R` `fit_w` |
| P2-10 SD population | Rows weighted w / n alternatives | `R/07_fit_sweep.R` `make_dR` |
| P2-11 aggression coverage | Unchanged: censored windows stay flagged, not dropped. **Open question for you:** is a complete-window subsample sensitivity enough, or should the variant drop `aggr_censored` strata by default? | `R/06_make_variants.R` |
| P2-12 side tags | Original-joiner status kept (`orig_joiner`); tags include window-promoted joiners; `hard_disp` / `reviewed_disp` drop every event of those disputes | `R/05`, `R/06` |
| P2-13 summarizer | Manifest-only, duplicate-key and spec-count checks, census per model, plausible/failing/all, family cells, family-balanced share, Wald computable/requested | `R/11_summarize.R` |
| P3-14 option hygiene | Ignored option combinations refused; README states the actual risk set and that pooling is an assumption | `R/07`, `README.md` |
| B permutations | Labelled artificial-assignment diagnostics (`diag_*`, `p_diag`) | `R/08`, README |
| B measurement error | Both signs, 5 independent noise draws, rank-correlation diagnostic per perturbation | `R/07`, `run_all.sh` §6 |
| B influence drops | Kept as adversarial diagnostics; history-chosen lists drafted with a fixed coding rule, run only after PI approval | `docs/historical_exclusions_DRAFT.md`, `run_all.sh` §7 |
| C reporting | See P2-13; percentages labelled as grid summaries | `R/11_summarize.R` |

## What I want
1. For each finding: fixed / partly fixed / not fixed, with file:line evidence and, if not fixed, the failure scenario that remains.
2. Any NEW P1/P2 bug introduced by these changes (especially in `run_all.sh`'s fingerprint, skip and incomplete logic, and in `run_stepB`).
3. A judgement on the history-chosen exclusion rule and lists (`docs/historical_exclusions_DRAFT.md`): are the inclusions and exclusions defensible under the stated rule?
4. Go / no-go for the production run.

Write the complete review to `docs/REVIEW_delta.md`. Read-only otherwise: do not modify, move or delete any other file, and do not run the R scripts.
