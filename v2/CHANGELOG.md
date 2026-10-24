# v2 CHANGELOG — how each audit finding (F1–F12) is addressed

Audit source: `REFERENCE_roar_and_audit.md` (Oct 22, 2026). The old pipeline is reference only; nothing was copied.

| # | Old issue | v2 |
|---|---|---|
| F1 | `polar_coords.r` collapsed a `1 x N x K` slice to a scalar, so every actor got the same r and θ | The typology is computed from the model's own generated quantity `mu_state` (S x 2 per draw), never from a product of U and V slices. `typology.R::classify_types()` asserts a `[draws, S, 2]` array, asserts the output has S rows with non-constant r, and is unit-tested on planted centres (`tests/testthat/test_typology.R`). AME `u_i'v_j` is formed in `03_extract.R` from explicit `matrix(...) %*% t(matrix(...))` on aligned draws, with dimension checks. |
| F2 | Producing version of `polar_coords.r` lost; polar step not rerunnable | Every step is a script in `v2/R/` committed to git; `03_extract.R --fit=<file>` regenerates all typology outputs from a saved fit. |
| F3 | `stage2_data.rds` written by two scripts with different NA filters (N = 5280 vs 4846) | One script, `01_build_data.R`, writes every data file with distinct names per lag and variant, plus `manifest.csv` (nrow, ncol, sha256). `02_fit_model.R`, `benchmarks.R` and `run_measurement_real.R` load inputs only through `assert_matches_manifest()`, which stops on a nrow or digest mismatch. The cell assembly (NA handling, event filtering) lives in a single function, `assemble_cell.R`, shared by all consumers. |
| F4 | Initiation coded `orig1 == 1`; `sidea` never used, so disputes entered as both A→B and B→A | `01_build_data.R` builds events from `cow_mid_dirdisps` with `sidea1 == 1 & sidea2 == 0` and `orig1 == 1 & orig2 == 1`: initiator = Side A originator, target = Side B originator, one direction per dispute pair; multilateral disputes contribute every A-originator x B-originator pair. The build asserts no reversed pair exists. `INCLUDE_JOINERS` flag for robustness. Unit-tested on toy dispute-years. |
| F5 | Outcome positions at t−5 while predictor windows ended at t−1; `lag()` was a row lag, not a calendar lag | Positions x_ij at t−1; change windows [t−1−h, t−1] for h in {1, 3, 5, 10}; all lags are calendar joins on `(ccode, year − h)`. No `lag()` is used for any lagged level (the only `lag()` left is inside the peace-spell computation and operates on calendar `year` values, so gaps cannot produce bridged differences). Unit test shows a missing year is dropped rather than bridged. |
| F6 | AME fits not converged; anchor rows for U/V identification were whichever ccodes sorted first as strings; chains in different modes pooled | Anchors are the K states with the most initiation events in the estimation sample (`choose_anchors()`), passed to Stan as data. The anchor block of U is lower-triangular with a positive diagonal; V is column-centred. `03_extract.R` applies an orthogonal Procrustes rotation per chain before pooling and logs the rotation deviation (non-identity rotations would flag non-mixing). Convergence (R-hat, ESS, divergences) is reported per fit and in `sim/SIM_REPORT.md`; K is a data input tested for K = 1, 2, 3 (likelihood test uses K = 2; the program rejects K > S−1). Review fix (Codex item 5): U is column-centred as well as V (hard constraint on the last non-anchor row per column, compatible with the anchor block), removing the translation U[,k] + h, b_j − h·V[j,k] that left the likelihood unchanged but moved u'v. u'v is therefore identified up to rotation, which the anchors fix and the Procrustes step guards; `03_extract.R` asserts the centring before exporting u'v. The simulation compares the posterior u'v with the doubly-centred true u'v. |
| F7 | Stan files emitted no `log_lik`; "joint AME" ELPD was actually brms; a/b comparisons measured RE structure, not jointness | Both Stan programs emit `log_lik_stage1[N1]` and `log_lik_stage2[E]`, kept separate. LOO comparisons (`measurement_tests.R::loo_compare_ab`) are between fits of the SAME program with the `joint` data flag toggled, so a vs b differs only in Σ's off-diagonal and corr(a_i1, a_i2). Stage 1 and stage 2 are estimated jointly in one program per family, so the comparison is about jointness of the dimensions, not about which engine was used. |
| F8 | brms base fits and AME base fits on different samples | The benchmarks (`benchmarks.R`) read the same manifest-checked files and the gap regression uses exactly the event rows that enter stage 2 (`assemble_gap_regression()` takes the `sy`/`alts`/`ev` frames). Cross-engine comparisons are labelled as benchmarks, not as model comparisons. |
| F9 | `_full` scripts did not NA-filter controls; result depended on which data file was on disk | `01_build_data.R` asserts no NA in every model column before writing (`assert_no_na`), records how NA controls were handled (`*_missing` indicator columns, documented coalescing), and the manifest digest makes "which file is on disk" an error instead of a silent branch. |
| F10 | Near time-invariant `shared_border_land` alongside dyad random effects → prior-driven coefficient | v2 has no dyad random effects; contiguity enters only as a fixed dyadic control γ with a target random intercept b_j, which does not absorb it. Dyad-level random effects were not added (listed as an unimplemented suggestion for the PI). |
| F11 | Stan and some brms calls lacked an explicit seed | `V2_SEED = 123` in `00_helpers.R`; every `sampling()`, `brm()`, `set.seed()`, and alternative-sampling call passes it explicitly; the simulator takes `seed = 123`. |
| F12 | Floating modulo thinning in plot code | No plotting code in v2 yet; draw thinning in `extract_draws()` uses `sample.int()` on integer indices with the fixed seed. |

Auditor's note on additive AME effects reproducing dM = M_i − M_j: v2 does not model (dM, dR) as an outcome. The outcome is the chosen target; dM and dR are covariates inside a quadratic utility, so sender and receiver effects (a_i, b_j) cannot reproduce the outcome algebraically.


## Post-review changes (Oct 22, 2026)

External Codex review of the first v2 commits (seven correctness items), each fixed in its own commit.

| item | fix |
|---|---|
| 1 | Conflict-history counts: each state-dispute and dyad-dispute counted once, at its first onset year (joiners entering later no longer create extra onsets). Unit test with a later joiner. |
| 2 | `measurement_tests.R` builds Stan data with `allow_orphan_initiations = TRUE` everywhere, as the fit path does. |
| 3, 4 | Sampled alternatives: each event gets its own sampled set (chosen + m uniform others); the IV sample is drawn uniformly and outcome-blind for every state-year with the log(n_full/m) correction (Jensen downward bias documented); full choice sets are the default everywhere. Monte Carlo tests reproduce the review cases. |
| 5 | U columns centred as well as V (hard constraint, compatible with the anchor block); u'v identified up to rotation (anchors + Procrustes); extraction asserts the centring. |
| 6 | The "a" model puts the LKJ(2)-matched Beta(2.5, 2.5) marginal on corr(c, a1) and corr(c, a2) with PD enforced; prior-only check Var ≈ 1/6 in both a and b. `prior_only` data flag added. |
| 7 | Typology F1 guard is a shape / identity check (dims, unique ids, non-recycled draws) instead of a radius-inequality test. |

PI additions after the review: one-dimensional targeting models (M1m, M1r, M2m, M2r via `D_mode`), capital-distance Z option (`--distance`, `--contiguity`), main specification = full all-states choice set, the 96-cell main grid with a Roar dry-run mode, per-group cross-model comparison, and the appendix typology outputs (9-cell classification, reach, focus).

Later on Oct 22: within-chain parallelism (reduce_sum_static over state-years, `--threads_per_chain`, cpus = chains × threads in the sbatch scripts, equivalence tests, gradient profile); material position changed to log CINC share standardised within year (zero/NA guard, tests); eigenvector-centrality robustness run unscheduled.

## Codex review of 36c0f7b (2026-10-22)
- P2 typology: `focus_sigma()` / `focus_states()` now take a required `D_mode`, so relational-only (M1r/M2r) focus output is labelled `sd_R` (was mislabelled `sd_M`). Test added.
- P2 data build: the zero/NA CINC-share check now runs on the full state-year panel before any filtering (previously NA rows were dropped first, so a missing share silently removed the state). Real data unaffected (0 NA, 0 zero); build re-run and complete. Test added.

## Fatal-MID coding (Oct 23, 2026; PI decision after the final independent audit, finding 3)

| change | detail |
|---|---|
| Wars are fatal | `fmid` now contains every war (hostility level 5). COW wars require 1,000+ battle deaths, so a non-fatal war is a coding gap. Adds 20 events (14 with all COW fatality codes missing, mostly WWI/WWII pairs; 6 coded non-fatal by COW). |
| MIC fill | Events whose COW fatality codes are all missing take MIC 1.0 confrontation `fatalmin` ≥ 1 (`v2/data_raw/mic-conf-1.0.csv`, sha256-checked). Fills 182 of the 229 missing events (112 fatal, 70 non-fatal); 33 are neither wars nor in MIC and stay non-fatal. MIC never overrides a coded COW value. |
| `fmidcow` sample | The previous coding is kept as sample `fmidcow` (appendix) and reproduces the old `y1_fmid` exactly. `02_fit_model.R --sample=fmidcow` runs it; it is not in the main grid. |
| Counts | fatal events 627 (was 495); lag-3 stage-1 fatal initiations 534 (was 431). `mid` and `war` data are byte-for-byte unchanged apart from the new columns (`mic_fatalmin`, `fatal_mic`, `in_fmidcow`, `fatal_source`, `y1_fmidcow`, `n_ev_fmidcow`). |
| Tests | `test_data_coding.R`: truth table for `code_fatal_main()`, MIC-fill, war-rule and no-override cases. |

Logged but not used. Where both COW and MIC code fatality they disagree on 120 events. MIC omits 700 of the 2,770 COW events and recodes Side A in some disputes. Switching the event base or sides to MIC is a separate PI decision.

### Addendum (Oct 23, 2026): GML as third source and an upper-bound sample

| change | detail |
|---|---|
| GML fill | When COW and MIC are both missing, `fmid` takes the GML MID 2.2.1 coding of the same directed pair (peacesciencer `gml_mid_ddydisps`, same rule as COW). Codes 11 of the 33 remaining events (5 fatal). 22 stay uncoded and non-fatal. Order of precedence is war rule, COW, MIC, GML; no source overrides an earlier coded value. |
| `fmidub` sample | Upper bound for the appendix: `fmid` plus the 22 uncoded events as fatal (654 events). `02_fit_model.R --sample=fmidub`; not in the main grid. |
| Counts | fatal events 632 (COW only 495); lag-3 stage-1 fatal initiations 539 (COW only 431). `mid`, `war` and `fmidcow` data unchanged. |
| Tests | `code_fatal_main()` now takes `fatal_gml`; tests cover GML-only fill, MIC beating GML, and the upper bound. 437 expectations pass. |

## mu intercept (2026-10-23; PI-approved)

Problem. In both Stan programs the ideal target position was `mu[n] = X[n] beta + a_i` with `a_i ~ MVN(0, Omega)` mean zero and no intercept column in X. The two conflict-history columns were centred but `cowmaj_tm1` and `dem_tm1` were uncentred 0/1, so the population-average ideal gap was pinned at 0 for non-major non-democracies. Since x is the initiator-minus-target gap, that intercept is the average preferred gap (bully > 0, peer about 0, challenger < 0), the paper's central quantity, and the model could only express it through the state offsets a_i, which the mean-zero prior shrinks back towards 0.

| change | detail |
|---|---|
| Stan | `row_vector[D] mu0` added to both programs (length 1 in the m / r models through `D`). `mu = mu0 + X beta + a_i` in the reduce_sum partial sum, `Mu` in generated quantities, `mu_state` at `x_ref` and `mu_event`. Prior `mu0 ~ normal(0, prior_mu0_sd)`, `prior_mu0_sd` is data. No profile() inside the partial sum. |
| R data | `build_stan_data(prior_mu0_sd = 2)` passes the prior sd; every caller (fit script, sim, tests, measurement tests) goes through it. `assemble_cell.R` centres `cowmaj_tm1` and `dem_tm1` like the conflict-history columns (dM_i, dR_i untouched), so mu0 is the ideal gap at zero power change and average controls; alpha0 in stage 1 is the log-odds at the same point. |
| R reimplementations | `stage2_loglik_R()`, `measurement_tests.R::extract_draws()` / `stage2_scores()` and `sim/simulate.R` (new `mu0` argument, default 0; `truth$mu0`, `truth$mu_state_ref0 = mu0 + a_i`) include mu0 and stop if it is missing. `02_fit_model.R`, `run_sim.R`, `run_sim_1d.R`, `compare_refit.R` summarise mu0. |
| Tests | The fixed-draw likelihood tests use a non-zero mu0 (0.7, -0.4) and also check `mu_event`; the scorer tests use (0.6, -0.3) and assert both reimplementations fail without mu0. New `test_mu0_recovery.R` fits M1b to data simulated with mu0 = (0.8, -0.5) (S = 20, 12 years) and requires the 90% interval to cover the truth. |

Not changed. Other priors, data coding, sampling settings, choice sets. Existing simulation results in `sim/SIM_REPORT.md` predate mu0 and were not rerun here.
