# Simulation report (v2 pipeline)

**SUPERSEDED (2026-10-25).** The results below were produced by the OLD simulated world (`simulate_targeting()`: z-scored positions, mu0 = 0, four X controls, no choice-set masks, no event weights, r0 = 0.8, x_ref = 0) and the programs before the review fixes of 2026-10-25 (disk prior on the "a" correlations, explicit rho_x prior, event weights). `run_sim.R`, `run_sim_1d.R` and `run_measurement_tests.R` now simulate the CURRENT model (`simulate_current()`, including a persistent-position variant) and must be rerun (`slurm/10_simulations.sbatch`, open account) before any conclusion here is cited; `write_sim_report.R` then regenerates this file. The conclusions below have not been rewritten.

Generated 2026-10-22 16:44 by `v2/sim/write_sim_report.R` from `v2/sim/results/*.csv`. Seed 123 throughout.

## Design

Synthetic world from `v2/sim/simulate.R`: 40 states, 25 usable years (1000 state-years), full choice sets of 39 alternatives (39,000 alternative rows), one target per initiating state-year. Planted initiator types as offsets a_i in the (dM, dR) plane (bully (+1.6, +1.6), peer (0, 0), challenger (−1.6, −1.6), hybrids (±1.6, ∓1.6), noise SD 0.35); c_i correlated with a_i so Ω has non-zero correlations; Σ with SDs 0.8 and ρ = 0.5; β with non-zero cross terms; τ = 0.5; γ = (0.8, −0.5). Truth M2b adds a K = 2 u'v structure (entries N(0, 0.7)). Because a_i is a type mixture, the 'true Ω' is the realised covariance of (c_i, a_i1, a_i2).

Fits: M1a, M1b, M2a, M2b (K = 2) on data from M1b and from M2b, with the FIXED programs (post-review commits 004317c..9374e8d). Container limits (2 cores, 8 GB): 2 chains × 600 iterations (300 warmup) per fit instead of the 4 × 2000 planned for Roar, so ESS is small and intervals are approximate. The one fit that no fix touched (M1b on M1b data, full choice sets, joint M1 program) was kept from the first run; every other fit was rerun after the fixes. The `reused_saved_fit` column below records this. Simulations used full choice sets only (sampled = 0), so the sampled-alternatives fixes (items 3 and 4) do not bear on them.

## Convergence

| truth | model | reused_saved_fit | elapsed_min | divergences | max_treedepth_hits | worst_rhat | min_ess_bulk | frac_rhat_gt_1.05 |
|---|---|---|---|---|---|---|---|---|
| M1b | M1a | FALSE | 14.100 | 0.000 | 0.000 | 1.021 | 106.865 | 0.000 |
| M1b | M1b | TRUE | 16.900 | 0.000 | 0.000 | 1.046 | 112.980 | 0.000 |
| M1b | M2a | FALSE | 26.000 | 0.000 | 0.000 | 1.175 | 8.701 | 0.086 |
| M1b | M2b | FALSE | 23.900 | 0.000 | 0.000 | 1.205 | 4.161 | 0.107 |
| M2b | M1a | FALSE | 15.000 | 0.000 | 0.000 | 1.032 | 83.146 | 0.000 |
| M2b | M1b | FALSE | 16.400 | 0.000 | 0.000 | 1.026 | 91.638 | 0.000 |
| M2b | M2a | FALSE | 20.700 | 0.000 | 0.000 | 1.264 | 7.811 | 0.209 |
| M2b | M2b | FALSE | 28.000 | 0.000 | 0.000 | 1.081 | 11.742 | 0.096 |

## Parameter recovery and 90% interval coverage

| truth_family | model | n_params | coverage90 | mean_abs_err |
|---|---|---|---|---|
| M1b | M1a | 21.000 | 1.000 | 0.119 |
| M1b | M1b | 23.000 | 0.913 | 0.106 |
| M1b | M2a | 21.000 | 1.000 | 0.117 |
| M1b | M2b | 23.000 | 0.957 | 0.111 |
| M2b | M1a | 21.000 | 0.905 | 0.158 |
| M2b | M1b | 23.000 | 1.000 | 0.127 |
| M2b | M2a | 21.000 | 1.000 | 0.168 |
| M2b | M2b | 23.000 | 1.000 | 0.131 |

Coverage is the share of free structural parameters (α0, α, τ, β (8), γ, Σ, Ω) whose truth lies in the 90% interval; parameters fixed at zero by construction in the 'a' models (Σ[1,2], Ω[2,3]) are excluded from their coverage count. Per-parameter detail follows for the two matched fits (truth = fitted family).

### M1b fitted to M1b data

| param | truth | mean | lo | hi | covered |
|---|---|---|---|---|---|
| alpha0 | -1.034 | -0.852 | -1.213 | -0.546 | TRUE |
| alpha[1] | 0.600 | 0.975 | 0.416 | 1.464 | TRUE |
| alpha[2] | -0.300 | -0.280 | -0.398 | -0.166 | TRUE |
| tau | 0.500 | 0.444 | 0.300 | 0.603 | TRUE |
| beta[1,1] | 0.500 | 0.539 | 0.203 | 0.897 | TRUE |
| beta[1,2] | 0.300 | 0.513 | 0.110 | 0.922 | TRUE |
| beta[2,1] | 0.200 | 0.131 | -0.206 | 0.500 | TRUE |
| beta[2,2] | 0.600 | 0.847 | 0.419 | 1.271 | TRUE |
| beta[3,1] | -0.300 | -0.295 | -0.378 | -0.217 | TRUE |
| beta[3,2] | 0.100 | 0.057 | -0.027 | 0.145 | TRUE |
| beta[4,1] | 0.100 | 0.060 | -0.012 | 0.130 | TRUE |
| beta[4,2] | -0.200 | -0.228 | -0.322 | -0.136 | TRUE |
| gamma[1] | 0.800 | 0.810 | 0.574 | 1.037 | TRUE |
| gamma[2] | -0.500 | -0.435 | -0.711 | -0.145 | TRUE |
| Sigma[1,1] | 0.640 | 0.539 | 0.449 | 0.640 | FALSE |
| Sigma[2,2] | 0.640 | 0.630 | 0.517 | 0.782 | TRUE |
| Sigma[1,2] | 0.320 | 0.187 | 0.118 | 0.274 | FALSE |
| Omega[1,1] | 0.454 | 0.434 | 0.191 | 0.772 | TRUE |
| Omega[1,2] | 0.590 | 0.461 | 0.196 | 0.830 | TRUE |
| Omega[1,3] | 0.459 | 0.419 | 0.095 | 0.809 | TRUE |
| Omega[2,2] | 1.754 | 1.753 | 1.117 | 2.621 | TRUE |
| Omega[2,3] | 0.667 | 0.384 | -0.159 | 1.003 | TRUE |
| Omega[3,3] | 1.789 | 2.120 | 1.389 | 3.075 | TRUE |

### M2b fitted to M2b data

| param | truth | mean | lo | hi | covered |
|---|---|---|---|---|---|
| alpha0 | -1.139 | -0.901 | -1.300 | -0.556 | TRUE |
| alpha[1] | 0.600 | 0.986 | 0.378 | 1.627 | TRUE |
| alpha[2] | -0.300 | -0.294 | -0.421 | -0.174 | TRUE |
| tau | 0.500 | 0.422 | 0.283 | 0.569 | TRUE |
| beta[1,1] | 0.500 | 0.206 | -0.200 | 0.618 | TRUE |
| beta[1,2] | 0.300 | 0.151 | -0.325 | 0.589 | TRUE |
| beta[2,1] | 0.200 | 0.239 | -0.155 | 0.629 | TRUE |
| beta[2,2] | 0.600 | 0.750 | 0.295 | 1.190 | TRUE |
| beta[3,1] | -0.300 | -0.315 | -0.412 | -0.228 | TRUE |
| beta[3,2] | 0.100 | 0.091 | -0.005 | 0.194 | TRUE |
| beta[4,1] | 0.100 | 0.082 | -0.000 | 0.154 | TRUE |
| beta[4,2] | -0.200 | -0.218 | -0.309 | -0.128 | TRUE |
| gamma[1] | 0.800 | 0.963 | 0.738 | 1.190 | TRUE |
| gamma[2] | -0.500 | -0.398 | -0.701 | -0.101 | TRUE |
| Sigma[1,1] | 0.640 | 0.601 | 0.496 | 0.728 | TRUE |
| Sigma[2,2] | 0.640 | 0.728 | 0.600 | 0.891 | TRUE |
| Sigma[1,2] | 0.320 | 0.295 | 0.194 | 0.405 | TRUE |
| Omega[1,1] | 0.454 | 0.675 | 0.348 | 1.095 | TRUE |
| Omega[1,2] | 0.590 | 0.631 | 0.302 | 1.041 | TRUE |
| Omega[1,3] | 0.459 | 0.225 | -0.122 | 0.632 | TRUE |
| Omega[2,2] | 1.754 | 1.692 | 1.101 | 2.515 | TRUE |
| Omega[2,3] | 0.667 | 0.494 | -0.066 | 1.046 | TRUE |
| Omega[3,3] | 1.789 | 2.254 | 1.417 | 3.298 | TRUE |

### Independent models on joint data: Σ and Ω cross terms

| truth_family | model | param | truth | mean | lo | hi | covered |
|---|---|---|---|---|---|---|---|
| M1b | M1a | Sigma[1,2] | 0.320 | 0.000 | 0.000 | 0.000 | FALSE |
| M1b | M1a | Omega[1,2] | 0.590 | 0.437 | 0.138 | 0.797 | TRUE |
| M1b | M1a | Omega[1,3] | 0.459 | 0.363 | 0.014 | 0.766 | TRUE |
| M1b | M1a | Omega[2,3] | 0.667 | -0.000 | 0.000 | 0.000 | FALSE |
| M1b | M1b | Sigma[1,2] | 0.320 | 0.187 | 0.118 | 0.274 | FALSE |
| M1b | M1b | Omega[1,2] | 0.590 | 0.461 | 0.196 | 0.830 | TRUE |
| M1b | M1b | Omega[1,3] | 0.459 | 0.419 | 0.095 | 0.809 | TRUE |
| M1b | M1b | Omega[2,3] | 0.667 | 0.384 | -0.159 | 1.003 | TRUE |
| M1b | M2a | Sigma[1,2] | 0.320 | 0.000 | 0.000 | 0.000 | FALSE |
| M1b | M2a | Omega[1,2] | 0.590 | 0.411 | 0.117 | 0.752 | TRUE |
| M1b | M2a | Omega[1,3] | 0.459 | 0.413 | 0.144 | 0.787 | TRUE |
| M1b | M2a | Omega[2,3] | 0.667 | 0.000 | 0.000 | 0.000 | FALSE |
| M1b | M2b | Sigma[1,2] | 0.320 | 0.192 | 0.111 | 0.285 | FALSE |
| M1b | M2b | Omega[1,2] | 0.590 | 0.457 | 0.157 | 0.840 | TRUE |
| M1b | M2b | Omega[1,3] | 0.459 | 0.478 | 0.202 | 0.850 | TRUE |
| M1b | M2b | Omega[2,3] | 0.667 | 0.426 | -0.076 | 0.978 | TRUE |
| M2b | M1a | Sigma[1,2] | 0.320 | 0.000 | 0.000 | 0.000 | FALSE |
| M2b | M1a | Omega[1,2] | 0.590 | 0.584 | 0.206 | 0.987 | TRUE |
| M2b | M1a | Omega[1,3] | 0.459 | 0.162 | -0.153 | 0.528 | TRUE |
| M2b | M1a | Omega[2,3] | 0.667 | -0.000 | -0.000 | 0.000 | FALSE |
| M2b | M1b | Sigma[1,2] | 0.320 | 0.268 | 0.180 | 0.376 | TRUE |
| M2b | M1b | Omega[1,2] | 0.590 | 0.616 | 0.285 | 1.002 | TRUE |
| M2b | M1b | Omega[1,3] | 0.459 | 0.262 | -0.077 | 0.624 | TRUE |
| M2b | M1b | Omega[2,3] | 0.667 | 0.466 | -0.052 | 1.063 | TRUE |
| M2b | M2a | Sigma[1,2] | 0.320 | 0.000 | 0.000 | 0.000 | FALSE |
| M2b | M2a | Omega[1,2] | 0.590 | 0.616 | 0.278 | 1.041 | TRUE |
| M2b | M2a | Omega[1,3] | 0.459 | 0.233 | -0.129 | 0.694 | TRUE |
| M2b | M2a | Omega[2,3] | 0.667 | 0.000 | -0.000 | 0.000 | FALSE |
| M2b | M2b | Sigma[1,2] | 0.320 | 0.295 | 0.194 | 0.405 | TRUE |
| M2b | M2b | Omega[1,2] | 0.590 | 0.631 | 0.302 | 1.041 | TRUE |
| M2b | M2b | Omega[1,3] | 0.459 | 0.225 | -0.122 | 0.632 | TRUE |
| M2b | M2b | Omega[2,3] | 0.667 | 0.494 | -0.066 | 1.046 | TRUE |

## Type classification from posterior μ_i (x_ref = 0, r0 = 0.8)

| truth | model | accuracy | mean_abs_err_mu |
|---|---|---|---|
| M1b | M1a | 0.875 | 0.362 |
| M1b | M1b | 0.875 | 0.299 |
| M1b | M2a | 0.900 | 0.372 |
| M1b | M2b | 0.900 | 0.301 |
| M2b | M1a | 0.825 | 0.424 |
| M2b | M1b | 0.875 | 0.315 |
| M2b | M2a | 0.825 | 0.440 |
| M2b | M2b | 0.900 | 0.311 |

Disagreement between the independent (a) and joint (b) classifications, and their accuracies:

| truth | family | disagreement_a_vs_b | acc_a | acc_b |
|---|---|---|---|---|
| M1b | M1 | 0.050 | 0.875 | 0.875 |
| M1b | M2 | 0.050 | 0.900 | 0.900 |
| M2b | M1 | 0.100 | 0.825 | 0.875 |
| M2b | M2 | 0.100 | 0.825 | 0.900 |

## LOO (stage 1 and stage 2 kept separate)

| truth | model | elpd_stage2 | se_stage2 | bad_k_stage2 | elpd_stage1 | se_stage1 | bad_k_stage1 |
|---|---|---|---|---|---|---|---|
| M1b | M1a | -1147.8 | 26.4 | 4.0 | -598.8 | 12.3 | 0.0 |
| M1b | M1b | -1139.3 | 26.3 | 0.0 | -597.8 | 12.2 | 0.0 |
| M1b | M2a | -1158.7 | 27.4 | 12.0 | -598.3 | 12.2 | 0.0 |
| M1b | M2b | -1152.2 | 27.5 | 11.0 | -598.0 | 12.2 | 0.0 |
| M2b | M1a | -1218.0 | 26.9 | 2.0 | -584.1 | 13.5 | 0.0 |
| M2b | M1b | -1202.3 | 27.0 | 4.0 | -583.6 | 13.6 | 0.0 |
| M2b | M2a | -1190.9 | 28.8 | 14.0 | -584.9 | 13.6 | 0.0 |
| M2b | M2b | -1174.8 | 28.8 | 13.0 | -583.7 | 13.6 | 0.0 |

Pairwise comparisons (positive = second model better):

| truth | comparison | stage2_elpd_diff_b_minus_a | stage2_se_diff | stage1_elpd_diff_b_minus_a | stage1_se_diff |
|---|---|---|---|---|---|
| M1b | M1a vs M1b | 8.4 | 3.8 | 1.0 | 0.8 |
| M1b | M2a vs M2b | 6.5 | 4.1 | 0.3 | 0.8 |
| M1b | M1b vs M2b | -12.9 | 6.8 | -0.2 | 1.1 |
| M2b | M1a vs M1b | 15.8 | 5.2 | 0.6 | 1.1 |
| M2b | M2a vs M2b | 16.1 | 5.7 | 1.2 | 1.2 |
| M2b | M1b vs M2b | 27.5 | 10.8 | -0.1 | 1.0 |

## AME: recovery of u'v (doubly centred truth; U, V Procrustes-aligned across chains)

| truth | model | cor_UV_offdiag | coverage90_UV_offdiag | sd_UV_posterior_mean | max_procrustes_rotation_dev |
|---|---|---|---|---|---|
| M1b | M2a |   NA |   NA | 0.334 | 0.834 |
| M1b | M2b |   NA |   NA | 0.338 | 0.902 |
| M2b | M2a | 0.402 | 0.929 | 0.503 | 1.379 |
| M2b | M2b | 0.421 | 0.937 | 0.497 | 0.742 |

`cor_UV_offdiag` is the correlation between the posterior-mean u_i'v_j and the true doubly-centred u_i'v_j over off-diagonal pairs; `coverage90_UV_offdiag` the share of pairs whose truth lies in the 90% interval; for M2 fitted to M1b data (no true u'v) `sd_UV_posterior_mean` shows how much structure the model invents. `max_procrustes_rotation_dev` near 0 means the anchor constraints already fixed the rotation.

## One-dimensional model (M1m) recovery

Smaller world (30 states, 12 years) generated from the material-only truth; M1m (right axis), M1r (wrong axis) and M1b (both axes) fitted to the same data.

| model | elapsed_min | divergences | worst_rhat | min_ess |
|---|---|---|---|---|
| M1m | 1.100 | 0.000 | 1.044 | 33.995 |
| M1r | 1.200 | 0.000 | 1.013 | 136.251 |
| M1b | 2.200 | 0.000 | 1.071 | 24.527 |

| param | truth | mean | lo | hi | covered |
|---|---|---|---|---|---|
| tau | 0.500 | 0.545 | 0.290 | 0.878 | TRUE |
| beta[1,1] | 0.500 | 0.535 | -0.140 | 1.191 | TRUE |
| beta[2,1] | 0.200 | 0.570 | 0.031 | 1.137 | TRUE |
| beta[3,1] | -0.300 | -0.288 | -0.447 | -0.144 | TRUE |
| beta[4,1] | 0.100 | 0.032 | -0.110 | 0.180 | TRUE |
| Sigma[1,1] | 0.640 | 0.758 | 0.539 | 1.058 | TRUE |
| Omega[1,1] (c) | 0.496 | 0.239 | 0.012 | 0.643 | TRUE |
| Omega[2,2] (a) | 2.002 | 2.946 | 1.672 | 4.911 | TRUE |
| Omega[1,2] (c,a) | 0.723 | 0.367 | -0.032 | 0.910 | TRUE |
| gamma[1] | 0.800 | 0.807 | 0.504 | 1.102 | TRUE |
| gamma[2] | -0.500 | -0.340 | -0.683 | -0.031 | TRUE |

| model | axis | accuracy_3level | mean_abs_err_mu |
|---|---|---|---|
| M1m | material | 0.933 | 0.268 |
| M1r | relational |   NA |   NA |
| M1b | both | 0.600 | 0.311 |

| model | elpd_stage2 | se_stage2 | elpd_stage1 | se_stage1 |
|---|---|---|---|---|
| M1m | -469.1 | 11.4 | -230.5 | 6.3 |
| M1r | -533.9 | 6.7 | -232.4 | 5.5 |
| M1b | -470.7 | 11.8 | -230.6 | 6.3 |

## Measurement tests (small world: 30 states, 12 years, truth M1b)

### LOO a vs b

| stage | elpd_a | elpd_b | elpd_diff_b_minus_a | se_diff | n_pareto_k_gt_0.7_a | n_pareto_k_gt_0.7_b |
|---|---|---|---|---|---|---|
| stage2_choice | -415.2 | -405.0 | 10.2 | 4.6 | 5.0 | 8.0 |
| stage1_initiation | -211.3 | -208.7 | 2.6 | 1.0 | 2.0 | 2.0 |

### Misclassification a vs b

| rate_a_vs_b | mean_abs_prob_shift | accuracy_a | accuracy_b |
|---|---|---|---|
| 0.200 | 0.050 | 0.767 | 0.833 |

### Temporal holdout (train on the first 9 usable years, score the last 3)

| spec | cutoff | n_train_events | n_test_events | hit_rate_top1 | hit_rate_top5 | mean_rank | mean_set_size | mean_log_score | baseline_log_score_uniform |
|---|---|---|---|---|---|---|---|---|---|
| a_independent | 11.000 | 125.000 | 36.000 | 0.056 | 0.417 | 9.334 | 29.000 | -3.115 | -3.367 |
| b_joint | 11.000 | 125.000 | 36.000 | 0.056 | 0.444 | 8.613 | 29.000 | -2.956 | -3.367 |

### Permutation test (dR shuffled within each choice set; rescore mode, 20 permutations)

| mode | n_perm | observed_log_score_a | observed_log_score_b | observed_advantage | mean_permuted_advantage | mean_permuted_log_score_a | mean_permuted_log_score_b | p_value_adv_ge_observed |
|---|---|---|---|---|---|---|---|---|
| rescore | 20.0 | -359.6 | -344.8 | 14.8 | -80.4 | -587.6 | -668.0 | 0.0 |

### Systematic measurement error (R inflated by delta for states that initiated in t−1; joint model refit)

| delta | share_states_reclassified | mean_shift_mu2 | mean_abs_shift_mu1 | n_bully | n_challenger | n_peer | baseline_accuracy |
|---|---|---|---|---|---|---|---|
| 0.500 | 0.100 | -0.158 | 0.151 | 6.000 | 8.000 | 6.000 | 0.833 |
| 1.000 | 0.200 | -0.324 | 0.138 | 7.000 | 9.000 | 6.000 | 0.833 |

### Benchmarks on the simulated frames (brms)

| model | term | Estimate | Q2.5 | Q97.5 |
|---|---|---|---|---|
| dyad_logit | Intercept | -4.253 | -4.554 | -3.976 |
| dyad_logit | cinc_ratio | 0.045 | -0.112 | 0.197 |
| dyad_logit | dR | 0.074 | -0.118 | 0.243 |
| dyad_logit | Z1 | 0.248 | -0.189 | 0.659 |
| dyad_logit | Z2 | -0.602 | -1.030 | -0.198 |
| gap_mvbind | dM_Intercept | -0.042 | -0.504 | 0.463 |
| gap_mvbind | dR_Intercept | 0.016 | -0.389 | 0.471 |
| gap_mvbind | dM_X1 | 0.731 | 0.445 | 1.035 |
| gap_mvbind | dM_X2 | 0.164 | -0.089 | 0.433 |
| gap_mvbind | dM_X3 | -0.084 | -0.149 | -0.017 |
| gap_mvbind | dM_X4 | 0.047 | -0.033 | 0.129 |
| gap_mvbind | dR_X1 | 0.205 | -0.081 | 0.491 |
| gap_mvbind | dR_X2 | 0.531 | 0.264 | 0.795 |
| gap_mvbind | dR_X3 | 0.068 | -0.004 | 0.135 |
| gap_mvbind | dR_X4 | -0.107 | -0.187 | -0.035 |

## Gradient-time profile and within-chain parallelism (reduce_sum)

Breakdown of one gradient evaluation on the 40-state cell (N1 = 1000 state-years, N2 = 39,000 alternative rows), measured with `log_prob` on the full data, on the same data with the events removed, and with `prior_only = 1`. Stan `profile()` sections (priors, stage2_utility, stage1_iv, stage2_events, stage1_bernoulli) are in both programs but rstan cannot read them; CmdStan would.

| model | ms_full | ms_stage1_and_iv_utilities | ms_stage2_events | ms_priors_and_transforms | share_stage1_iv_utilities | share_stage2_events | share_priors |
|---|---|---|---|---|---|---|---|
| M1b | 18.34 | 18.15 | 0.16 | 0.03 | 0.99 | 0.01 | 0.00 |
| M2b | 24.88 | 25.41 | -0.62 | 0.09 | 1.02 | -0.02 | 0.00 |

Sampling time for 1 chain, fixed seed, by `threads_per_chain` (reduce_sum_static over state-years; grainsize = N1 for 1 thread, about N1 / (4 threads) otherwise). This container has 2 cores, so 4 threads is oversubscribed; the real speed-up curve must be read from the Roar pilot (04_pilot.sbatch with THREADS_PER_CHAIN = 1 vs 8).

| model | threads_per_chain | grainsize | cores_available | oversubscribed | iterations | elapsed_s | s_per_iteration | ms_per_gradient |
|---|---|---|---|---|---|---|---|---|
| M1b | 1.00 | 1000.00 | 2.00 | FALSE | 40.00 | 107.66 | 2.69 | 20.28 |
| M1b | 2.00 | 125.00 | 2.00 | FALSE | 40.00 | 45.51 | 1.14 | 8.57 |
| M1b | 4.00 | 63.00 | 2.00 | TRUE | 40.00 | 37.78 | 0.94 | 7.12 |
| M2b | 1.00 | 1000.00 | 2.00 | FALSE | 40.00 | 121.27 | 3.03 | 23.75 |
| M2b | 2.00 | 125.00 | 2.00 | FALSE | 40.00 | 56.60 | 1.41 | 11.08 |
| M2b | 4.00 | 63.00 | 2.00 | TRUE | 40.00 | 44.53 | 1.11 | 8.72 |

Refit of M1b on M1b data with the reduce_sum program (threads = 2) against the earlier serial fit: 33 parameters, max |z| = 2.53, share |z| < 2 = 0.88, max |mean difference| / posterior SD = 0.21 (z uses MC standard errors sd / sqrt(n_eff) of both runs).

## What failed or is fragile

* Chains and iterations are far below the Roar settings; R-hat and ESS in the convergence table are the honest state of these short runs, not of the production fits.
* The M2 (AME) fits mix poorly at 2 × 300 post-warmup draws: worst R-hat 1.08–1.26 and bulk ESS 4–12 on some U/V-related quantities, while the M1 fits have R-hat ≤ 1.05. The Procrustes rotations between the two chains' mean U are far from identity (max |R − I| 0.74–1.38), so U and V themselves are not pinned down across chains even with the anchor block; u'v (rotation-invariant) is recovered with correlation 0.40–0.42 and 90% coverage 0.93–0.94 against the doubly-centred truth, i.e. weakly but with calibrated uncertainty. On Roar, run K = 2 with 4 × 2000 and check the rotation deviation before reading U or V.
* Σ is under-recovered in the joint fits on M1b data (Σ[1,2] = 0.19 vs 0.32, Σ[1,1] = 0.54 vs 0.64, both outside the 90% interval); on M2b data both are covered. The softmax identifies Σ⁻¹ jointly with the target effects b_j, so part of the targeting breadth can be absorbed there.
* In the 1-D world the two-dimensional M1b fit matches M1m on LOO but its material-axis classification is worse (0.60 vs 0.93) because a_i2 and ρ are estimated from a dimension that carries no signal; M1r (wrong axis) loses 65 ELPD. The one-dimensional models therefore need the model comparison, not the typology alone, to tell which axis matters.
* The permutation test's permuted advantage is strongly negative (−80): shuffling dR within choice sets hurts the joint model far more than the independent one, so the observed +14.8 advantage is well outside the permutation distribution (p = 0 of 20), but the statistic's sign under the null is not centred on zero and should be read as a calibration, not a p-value in the usual sense.
* The brms benchmark fits on simulated data are produced by `run_benchmarks_sim.R` in a separate R session because the brms compile was OOM-killed inside the measurement-test driver; see `benchmarks_sim.csv` if present.
* Ω is recovered as the covariance of a type mixture; its coverage reflects the MVN approximation to planted clusters.
* The AME u'v is identified only up to the U/V scale fixed by the N(0,1) priors and the rotation fixed by the anchors; read the u'v statistics, not U or V.
* The 'a' model on joint data: whether it loses on LOO and misclassifies depends on how far the planted correlations push the data; the numbers above are from one seed.

