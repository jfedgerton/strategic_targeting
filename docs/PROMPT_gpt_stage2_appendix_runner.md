# GPT audit request: stage-2 appendix runner (parallel multiverse + bootstrapped cells)

Repo zip attached (branch targeting-multiverse-wip). Please audit the code added for the stage-2 appendix. Earlier reviews covered the builder (R/s2_01), the fitter (R/lib/s2_fit.R rev 5) and the surface design; this request covers only the new runner layer.

## What it does
- Layer 1 (point fits). For every relational measure, fit the joint two-step FE model: stage-1 target-FE conditional logit, then the stage-2 four-category masked logit with bbar (mean target FE over the dispute's target set). Each measure uses its own positions (stage-1 dR, stage-2 gaps) and its own edges (stage-2 ties). Run on all measures at window 1 × fatal MIDs, and on plausible measures at window 1 × all MIDs.
- Layer 2 (bootstrap). A rule-based list of about 50 cells, fixed before results are seen. The two-stage dispute bootstrap refits both steps, with common draws across cells (seed 123). The main cell gets 2,000 draws, the others 500.
- Layer 3 (AME). Planned separately, not in this code.

## Files
- R/lib/s2_cell.R: `s2_spec_edges`, `s1_prepare`, `s1_fit_fe`, `s2_prepare`, `joint_fe_cell`, `joint_fe_draw`, `s2_base`.
- R/s2_05_sweep.R: the sweep, resumable by chunk.
- R/s2_06_pick_boot_cells.R: the cell rules (see its header).
- R/s2_07_boot_cells.R: the bootstrap runner.
- run_stage2.sh: the orchestration.
- R/s2_01_build_riskset.R diff: accepts init dumps; adds fatal_def `mid`.
- R/02_build_positions.R: the `--edges_out` mode.
- tests/test_s2_cell.R.
  - The cell reproduces R/s2_03_joint_fe.R exactly for the main measure, both point estimates and bootstrap draws 1–2.

## Questions
1. Does `s2_spec_edges` combine layers exactly as R/02 does for each weighting (W1raw, W2–W4, SRY W6/W7 and mm, IGO arm)?
2. Is any measure able to silently use the main measure's data (ties, gaps, scales, target set)?
3. Is a draw's weighting of stage 1 (`w × m[disp]`) and stage 2 (row weights `m[dispnum]`) the right two-stage bootstrap? Are the common random numbers across cells defensible?
4. Are the cell-picker rules a reasonable, pre-specified appendix design? Is anything missing (for example a centrality or PTA rule)?
5. Are there resumability, parallel (mclapply), or bash pitfalls that would give wrong or partial results without an error?
6. Is the window/sample pairing between builds and stage-1 dumps right everywhere?

Please rank findings P1 (wrong results), P2 (robustness), P3 (nits), with file:line and a fix.
