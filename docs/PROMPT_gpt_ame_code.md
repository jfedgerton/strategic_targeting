# GPT code review request: joint two-stage AME (cells 1-3), before the first Roar run

Attached zip (branch targeting-multiverse-wip, commit in the zip name). Files to review:
- `stan/ame_joint.stan`: the model.
- `R/ame/ame_cell.R`: data preparation, fitting, summaries. The `--prep_only` mode builds a data bundle on the Mac; the fit-only mode runs on Roar.
- `slurm/ame_cells.sbatch`: the Roar job (3 cells: main = joint K = 0; cut = stage-1 fit, then stage 2 with bbar from 10 stage-1 draws; k1 = joint K = 1).

Design: `docs/ame_design.md` (rev 2, already revised after an independent audit). The model mirrors the two-stage FE model in `R/lib/s2_cell.R`; stage-2 inputs come from `s2_prepare()` and `S2_REF_DESIGN` (`R/lib/s2.R`), availability from `mnl_avail()` (`R/lib/s2_fit.R`).

## Known deviations from the memo (intentional for this first run)
- No Williams–Savitsky adjustment yet. Event weights enter as a raw pseudo-likelihood.
- No SBC or misspecification simulations yet. This run is a first look, not the gated main result.

## Please check
1. Does the Stan likelihood implement the memo exactly? Check the stage-1 clogit per event with weights, the stage-2 masked four-category logit, bbar from the draws, the dispute random effect, type-specific (δ) and era-specific (η) target effects, and the K = 1 term with sum-to-zero u and v orthogonalized by Pv.
2. reduce_sum correctness: slices compute utilities from data rows; events are contiguous blocks (e_start/e_end/e_chosen); no off-by-one errors.
3. Data preparation:
   - within-year standardization; pre-event expanding mean R̄ evaluated at t − 1;
   - between/within gaps dB = R̃_i,t−1 − R̄_j(t−1), dW = R̄_j − R̃_j,t−1;
   - column scaling and back-transformation of beta1;
   - stage-2 standardization and the per-category constant-term drop rule;
   - era assignment;
   - target sets for bbar;
   - opponent-side weight matrices;
   - the Pv projection.
4. Identification problems that would make the fit meaningless or prevent mixing (e.g., the α intercepts vs a mean-zero b̄; η vs b; δ vs the coalition surface; the Student-t on b_raw scaled by σ_b).
5. The cut-posterior implementation: is it a valid cut (multiple imputation over stage-1 draws, with draws pooled across imputations)?
6. Summary code: contrast definitions (gap ±1 = linear ± coefficient + squared coefficient); coalition − bilateral contrasts; Rhat/n_eff extraction.
7. Roar script: threading flags (STAN_THREADS, threads_per_chain), compile-once guard, resource request.

Rank findings P1 (wrong results / crash), P2 (robustness), P3 (nits), with file:line and a fix.
