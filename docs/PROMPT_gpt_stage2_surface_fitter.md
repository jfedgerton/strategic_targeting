# Prompt for GPT: audit of the stage-2 surface amendment, robustness variants and reference fitter

You approved the stage-2 build (GO with fixes; `docs/REVIEW_stage2_build_recheck_gpt.md`), and the fixes are in (revision 3). Since then the PI changed the stage-2 design: stage 2 should mirror stage 1. This package adds three things.
1. **Surface amendment** (`docs/stage2_design_memo.md` §13, "Surface amendment"). Each join alternative gets the stage-1 quadratic surface in the joiner's power and position gaps to the side it fights, with stage-1 scaling, plus the stage-1 directed controls toward that side.
2. **Opposing-side robustness variants** (`s2_surface()` aggregations, `s2_primary_pair()`, the builder's 40-variant loop).
3. **The reference fitter** `R/s2_02_fit_reference.R`, with:
   - the four-category logit;
   - two-way clustered vcov;
   - surface contrasts;
   - a pre-set separation rule;
   - the 40 variants and the linear design;
   - two permutation schemes using a likelihood-ratio statistic;
   - alliance- and contiguity-correlated measurement error.

Read `docs/stage2_build_notes.md` (section "Revision 4"), the memo amendment, `R/lib/s2.R` (from `s2_gap_scale` on), `R/s2_01_build_riskset.R` (surface and variants blocks), `R/s2_02_fit_reference.R`, `tests/test_s2.R` (33–34), and the logs in `docs/s2_build_logs/` and `docs/s2_fit_logs/`. For checking stage-1 conventions, `R/05_dump_events.R` and `R/07_fit_sweep.R` define dM and dR.

Please answer:
1. **Mirroring.** Does the stage-2 surface reproduce stage 1's construction of dM and dR: the sign convention, t − 1 positions, the pooled M scale, and the year-t R scale? Are the mirrored controls a faithful translation of `Z_DIRECTED` to a multi-member opposing side? Flag the democracy cutoff (polity2 ≥ 6, assumed) and any other assumption.
2. **Design.** Is dropping own and target-set centrality and the CINC ratio from the reference (nested in the surface) sound? Is there remaining collinearity or a quantity counted twice?
3. **Fitter correctness.** Check:
   - the likelihood, gradient and scaling back to original units;
   - the observed-information Hessian on scaled parameters;
   - the CGM sandwich (12 negative eigenvalues were zeroed in the reference fit);
   - the delta-method contrasts;
   - the separation rule. It drops `likely` and `dem_x_B` from the originator equation. Is a data-driven drop acceptable when fixed in advance?
4. **Robustness checks.** Are the permutation schemes valid tests? This covers B1 (gap vectors shuffled within dispute, null log-likelihood held fixed) and B2 (A/B labels shuffled among entrants within dispute, null refit). Is the measurement-error design (positions perturbed, surface rebuilt) right? Is the variant grid complete for "how the opposing side is defined"?
5. **What is still missing before the reference result is reportable?** For example: the pigeonhole bootstrap, predicted-probability contrasts, or windows 0 and 2–6.
6. **Verdict**: GO / GO with fixes / NO-GO for reporting the reference stage-2 model. Give numbered findings (P1 / P2 / P3) with file, line and a concrete fix.
