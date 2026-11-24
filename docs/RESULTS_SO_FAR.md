# Results so far (sandbox, 2026-11-24) — for checking that outputs make sense

Main measure BOOKraw_ptadesta, fatal MIDs, window 1. Two-stage dispute-level Bayesian bootstrap, 500 draws (487 usable; 13 failed on stage-1 separation). 95% percentile intervals; * = interval excludes 0. Gaps: + = initiator (stage 1) or candidate (stage 2) stronger / more central than the target or opposing side.

## Main table
| Quantity | (1) Stage 1 between targets | (3) Joint, stage 1 within targets |
|---|---|---|
| Bilateral, weaker target | -1.06 [-1.41, -0.75]* | -0.97 [-1.96, 0.27] |
| Bilateral, less central target | 0.06 [-0.16, 0.24] | -0.22 [-0.59, 0.06] |
| Coalition, weaker target | -1.72 [-3.53, -0.56]* | -1.46 [-2.87, -0.05]* |
| Coalition, stronger target | 1.14 [0.36, 2.19]* | 1.19 [-0.16, 2.32] |
| Coalition, less central target | 1.06 [0.50, 1.77]* | 0.84 [0.23, 1.58]* |
| Coalition, more central target | -1.15 [-2.55, -0.46]* | -0.95 [-2.10, -0.16]* |

| Stage 2 | (2) Stage 2 alone (511 disputes) | (3) Joint |
|---|---|---|
| Join A, target side stronger | -0.82 [-1.54, -0.45]* | -0.66 [-1.41, -0.34]* |
| Join B, attacking side stronger | -1.40 [-2.18, -1.07]* | -1.36 [-2.15, -1.01]* |
| Join B, attacking side weaker | 0.65 [0.09, 1.24]* | 0.56 [-0.01, 1.19] |
| Originator entry, recruitment ties | 2.46 [1.44, 3.23]* | 2.56 [1.76, 3.35]* |
| theta A / B / O | | -0.89 [-1.02, -0.22]* / -0.43 [-0.56, -0.08]* / 0.46 [-0.16, 1.06] |

## Robustness (100 draws)
- (2b) stage 2 on all 544 disputes and (3b) b = 0 for one-time targets: all signs and significance unchanged (theta A -0.92*, theta B -0.42*).
- New-state build (point fit): theta A -0.88, theta B -0.43.
- (3c) net ties (own-side minus opposing-side ties per option): stage-2 gaps and theta unchanged (theta A -0.89*, theta B -0.40*); net ties for joining the target side become clearly positive (1.77 [0.69, 3.28]*, own-side only 0.39 n.s.); originator entry 2.39*; joining the attacker n.s.

## Known issue to check (hypothesis)
theta's percentile interval is off-centre relative to its point estimate (theta A: point -0.89, bootstrap median -0.60; theta B: -0.43 vs -0.31). Stage-1 contrasts are centred (medians equal points). Likely errors-in-variables attenuation from re-estimating bbar in each draw. Decided: main text reports point, percentile interval and bootstrap median; BC interval in the appendix (theta A BC [-1.50, -0.79], theta B [-0.80, -0.30], theta O [-0.08, 1.21]).

## Not yet run end to end
- run_paper.sh as a whole (each step tested separately in the sandbox).
- AME: a 30-iteration smoke fit on a subsampled bundle ran end to end (sampling, summaries, convergence file); the Williams-Savitsky step correctly refused to run with 15 draws (fewer draws than fixed parameters). Real fits run on Roar after this audit.
