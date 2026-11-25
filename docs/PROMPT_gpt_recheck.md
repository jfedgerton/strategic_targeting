# Prompt for the GPT recheck of the paper pipeline (after the 2026-11-24 audit)

You audited this repository at commit `c5f1907` (your report is `audits/gpt_final_audit.md`). The attached zip holds the current branch `targeting-multiverse-wip` and `changes.diff` (all changes since `c5f1907`). Please recheck the changes, not the whole pipeline.

## What changed since your audit
1. **Your findings 2-9 were addressed** (commit `9ecb632`): typed worker results and retries in `R/main_table.R`; bootstrap retry when a summary exists (`R/s2_07_boot_cells.R`, `tests/test_s2_runner.sh`); grid-layer stamp hashes the data files (`run_paper.sh`); positions files carry grid provenance (`R/lib/grid.R`, `R/grid_build.R`); complete marginal-effects signature (`R/marginal_effects.R`); completion statuses (`R/s2_05_sweep.R`, `R/s2_07_boot_cells.R`, `R/grid_summary.R`); disclosures.
2. **Finding 1 (AME adjustment) is moot: the AME was dropped.** Step 8 of `run_paper.sh` is off by default; the AME code is kept but not run. The paper uses the FE models with the two-stage dispute Bayesian bootstrap.
3. **New main stage-2 design: ties to both sides.** Each joining option now has its own-side tie (side joined) and its opposing-side tie (side fought): A and O `tie_sum_A + tie_sum_B`, B `tie_sum_B + tie_sum_A` (`S2_REF_DESIGN` in `R/lib/s2.R`; the earlier own-only design is `S2_REF_DESIGN_OWN`). Supporting changes:
   - `R/lib/s2_cell.R`: `tie_terms()` returns the own-side term when an option has two tie terms; new `tie_opp_terms()`; `joint_fe_cell()` adds `axis = "tie_opp"` rows; `boot_names()` and `joint_fe_draw()` add `A/B/O:tie_opp`.
   - `R/s2_07_boot_cells.R`: point vector includes `tie_opp`.
   - `R/main_table.R`: opposing-tie rows in columns (2) and (3); robustness (3c) is now own-side ties only (net ties dropped: collinear with both ties in the model).
   - `R/grid_summary.R`: tracks "join target: ties to the attacker" and "originator: ties to the attacker".
4. **New exploration script** `R/explore_specs.R` (test code, not part of `run_paper.sh`): alternative specs refit on the common bootstrap draws (Mundlak, target x era FE, rivalry FE, candidate / dispute / two-way FE in stage 2, monthly competing risk, measurement error, permutation). Review only if time allows.

## What to check
1. Correctness of the both-ties change: own vs opposing labels per option (A and O fight the target side B; B fights the attacking side A), quantities mapped to the right rows, nothing silently NA in the default path, (3c) really own-only.
2. That every consumer of `boot_names()` / `joint_fe_cell()` output handles the three new quantities (resume of old cells must be refused, not mixed).
3. Whether the fixes for your findings 2-9 are correct and complete.
4. Bash under `set -euo pipefail` on macOS bash 3.2 for the changed lines of `run_paper.sh`.

## Report format
As before: P1 / P2 / P3 with file, line, failure scenario and fix; end with what you checked and found correct.
