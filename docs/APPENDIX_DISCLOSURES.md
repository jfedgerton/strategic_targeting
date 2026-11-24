# Appendix disclosures (collected 2026-11-24 from the code audits)

Items the appendix should state. Numbers are from the sandbox build on the main measure (fatal MIDs, window 1) and will be refreshed by the full run.

## Network measures
- Main measure: unweighted book measure (BOOKraw_ptadesta), 5 layers (IGO, ATOP defense, ATOP offense, DCA v2, DESTA PTA), mean / yearly max, log PageRank.
- PTA: DESTA treaty membership with actual withdrawers only; treaties entering force after 2014 excluded (196,196 linked dyad-years). This departs from the book build in 2014 until the book applies the same fix.
- IGO count comes from the dyadic COW layer; IGO cosine and hyperedge from the COW state-year file. The two sources differ in 3.8% of dyad-years (max 3 shared IGOs).
- PTA depth layer is built but not swept (PI decision 2026-11-24).
- Grid: 236 measures (raw 56, DCA v1 28, other alliances 56, four weighting schemes x 24). The 16 sy-scheme measures without IGO have constant positions in some stage-1 years and are skipped, so 220 measures per sample are fitted (440 cells).
- Weighted schemes: observed minus expected tie per layer (logit / PPML), positive part over the yearly max; sry = sender + receiver + year FE, sy = sender x year + receiver x year FE; geography = contiguity (1-2) or log capital distance.

## Timing and missing data
- All covariates and positions are measured at onset year - 1.
- Stage 2 build: 560 disputes -> 544 in the reference cohort. Rows removed: 1,262 (target not a state at t-1; 12 disputes), 455 (initiating side without CINC at t-1; 3 disputes), 414 (candidate not in the system at t-1), 142 (dispute window past 2014). Table: main_table_missingness.csv.
- Joint model: 33 disputes have no repeat target (every onset target chosen fewer than 2 times in stage 1), so no target trait; main-table columns (2) and (3) both use the remaining 511 disputes.
- Robustness: (2b) stage 2 alone on all 544 disputes; (3b) joint model with b = 0 for one-time targets; new-state cell using year t where t - 1 does not exist (stage 2 only; 560 disputes, 65,274 rows).

## Inference
- Two-stage dispute-level Bayesian bootstrap (Exp(1) dispute weights scaled to mean 1; Rubin 1981), common draws across stage 1, stage 2 and the joint model; 95% percentile intervals; a quantity is reported only with >= 80% finite draws.
- Window 0 has no originator alternative; its O quantities are not estimated.
- theta reporting (PI 2026-11-24): main text = point estimate, 95% percentile interval and bootstrap median; appendix = bias-corrected (BC) percentile interval. theta's draws are attenuated toward 0 relative to the point estimate (bbar re-estimated in every draw), so BC shifts the interval away from 0; at 500 draws its outer bound can sit at the most extreme draw. BCa (acceleration) not computed. Output: main_table_theta.csv.
- Grid decision regressions: OLS of estimates / p-values on measurement choices, SEs clustered by measure; layer effects pool weighted and unweighted blocks.

## Figures
- Gap units: material gap = log CINC share over its pooled scale; relational gap = log PageRank over the year's typical pairwise gap.
- Attack-plane axes are cut at the 1st / 99th percentile of each gap; captions report the share of points outside.
- F3 rates use only candidates for whom the option is available.

## AME (appendix)
- rstan; stage-1 gaps in FE units; K = 1 latent term reported only through sign-free summaries.
- Fits are weighted pseudo-posteriors; the Williams-Savitsky adjustment is not applied (open decision).
- Stage 2 uses era-specific intercepts; AME stage 2 keeps all 544 disputes (random effects pool one-time targets).
