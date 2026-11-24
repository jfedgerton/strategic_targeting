# Prompt for the GPT final audit of the paper pipeline (run_paper.sh)

You are auditing an R research pipeline before a long production run. The attached zip holds the repository (branch targeting-multiverse-wip) and the reports of 12 earlier independent audits (`audits/`). Read `docs/PAPER_PIPELINE.md`, `docs/APPENDIX_DISCLOSURES.md` and `docs/ame_design.md` first.

## What the pipeline does
A two-stage model of interstate conflict. Stage 1 is a conditional logit of which target an initiator picks (by initiator type: bilateral vs coalition originator), on material and relational gaps. Stage 2 is a categorical model of third-party joining (join attacker / join target / join as originator / stay out). The joint model carries the stage-1 target fixed effects into stage 2 as a dispute-level target trait (bbar, coefficient theta). The relational position is log PageRank on a multiplex network. The main measure is the unweighted book measure (`BOOKraw_ptadesta`). A 236-measure grid of alternative networks is swept on fatal MIDs and all MIDs. Inference uses a two-stage dispute-level Bayesian bootstrap with draws shared across stages and scripts. A Bayesian joint AME (Stan) is an appendix model fitted on an HPC cluster.

## Entry points
- `run_paper.sh` (steps 1-8) and `config.sh`
- Networks: `R/book_layers.R`, `R/book_measure_raw.R`, `R/book_measure_c12.R`, `R/grid_layers.R`, `R/lib/grid.R`, `R/grid_build.R`
- Models: `R/lib/s2.R`, `R/lib/s2_fit.R`, `R/lib/s2_cell.R`, `R/s2_01_build_riskset.R`, `R/s2_05_sweep.R`, `R/s2_07_boot_cells.R`
- Outputs: `R/main_table.R`, `R/marginal_effects.R`, `R/table1_figures.R`, `R/grid_pick_boot_cells.R`, `R/grid_summary.R`
- AME: `stan/ame_joint.stan`, `R/ame/ame_cell.R`, `slurm/ame_cells.sbatch`

## What to check
1. Correctness: wrong numbers, sign or label errors, rows silently dropped or duplicated, quantities mapped to the wrong labels.
2. The common bootstrap draws: `boot_weights()` must give identical draws in `s2_07_boot_cells.R`, `main_table.R` and `marginal_effects.R`. Main-table columns (2) and (3) must use the same disputes.
3. The default path versus the robustness switches (`fill_b`, `--newstate_t=1`). The defaults must reproduce the main design, and the robustness versions must change only what they say.
4. Provenance and resumability: stale outputs reused after code or input changes, fmid / mid inputs mixed, window / sample pairing.
5. Bash under `set -euo pipefail` on macOS (bash 3.2): silent exits, unquoted paths, steps that rerun or skip wrongly.
6. Whether every "must stop the run" check in `docs/PAPER_PIPELINE.md` is really a hard stop.
7. AME: likelihood indexing, identification, the K = 1 summaries, the convergence gate, and Roar portability.

Earlier audits fixed many items. Do not re-report something an audit report marks as fixed unless the fix is wrong.

## Report format
P1 = wrong results or a crash in the production run; P2 = fragile or could mislead; P3 = minor. For each item give the file and line, a concrete failure scenario, and a suggested fix. End with what you checked and found correct.
