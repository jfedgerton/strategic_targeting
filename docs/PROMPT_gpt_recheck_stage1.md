# Short recheck: stage-1 fixes S1 and S2 (your `docs/REVIEW_stage1.md`)

You reviewed commit 3015bbb on 2026-11-21 and found two P2s in the stage-1 code. Verify only the fixes in the newest commit (`git diff 3015bbb HEAD`). Do not re-review unchanged code.

| Finding | Fix | Where |
|---|---|---|
| S1 effective TYPES not fingerprinted | TYPES limited to `bcd`/`init` right after `source config.sh` (exit 1 otherwise, before any output folder is touched); the resolved `TYPES=` value is added to the fingerprint string | `run_all.sh` lines after `source config.sh`; FP line |
| S2 side-case drops skipped under init | `dropside_hard`/`dropside_reviewed` now run under init too; R06 logs, in the initiators step, how many retained strata carry side_case hard / reviewed, and the dropside log reports removed strata (0 reported, not skipped) | `run_all.sh` DSL line; `R/06_make_variants.R` initiators and dropside blocks |
| Wording (Step B population, all/pair, four-type, permutation, event study, influence) | New "Stage 1" section with a population table | `README.md` |

Synthetic check run by Claude (R 4.x, data.table): an 8-event bcjs frame with one promoted hard and one promoted reviewed coalition originator → initiators keeps 6 strata and logs hard 1 / reviewed 1; dropside hard removes 1, reviewed 2, hard_disp 1, reviewed_disp 2; roles stay bil/coal. `TYPES=foo bash run_all.sh` exits 1 with no output folder created.

Write to `docs/REVIEW_stage1_recheck.md`: for S1 and S2, fixed / partly fixed / not fixed with file:line evidence; any NEW P1/P2 introduced; go / no-go for `TYPES=init OUT=output_init bash run_all.sh`. Read-only otherwise; do not run the R scripts.
