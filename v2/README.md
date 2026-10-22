# Strategic targeting — v2 pipeline (R + Stan)

Two-stage spatial choice model of MID initiation and target selection, estimated jointly in Stan.
Built fresh on branch `v2-rebuild`; the old code is reference only. Simulation-tested; no real-data model has been run yet.

## Model (implemented exactly as decided)

Stage 1 (state-year n, state i): `y1 ~ Bernoulli(logit⁻¹(α0 + W'α + c_i + τ·IV))`.
Stage 2 (event e by state-year n, choice set C_n of all other system members at t with positions at t−1):

```
U_ij = −½ (x_ij − μ_n)' Σ⁻¹ (x_ij − μ_n) + b_j + z_ij'γ [+ u_i'v_j in M2]
x_ij = (M_i − M_j, R_i − R_j) at t−1          μ_n = X_n β + a_i   (β is P × 2)
P(j | n) = softmax over C_n                    IV_n = log Σ_{j∈C_n} exp U_ij
(c_i, a_i1, a_i2) ~ MVN(0, Ω)                  b_j ~ N(0, σ_b) (soft sum-to-zero)
```

Families: `m1_latent_additive.stan` (M1) and `m2_ame.stan` (M2 = M1 + u_i'v_j, K latent dimensions, K is data). Eight model codes, two data flags.

| code | `D_mode` | `joint` | stage-2 targeting space |
|---|---|---|---|
| M1m / M2m | 1 | – | material only (x = dM, scalar μ, scalar Σ, Ω 2×2 over (c, a)) |
| M1r / M2r | 2 | – | relational only (x = dR) |
| M1a / M2a | 3 | 0 | both dimensions, independent (Σ diagonal, corr(a_1, a_2) = 0; corr(c, a_d) free) |
| M1b / M2b | 3 | 1 | both dimensions, joint (full Σ, full Ω) |

Stage 1 is identical in all eight (same W, which includes dM_i and dR_i); only the stage-2 utility and hence the IV differ. All eight score the same events and choice sets, so LOO compares like with like. The correlations shared across models carry the same prior ((r+1)/2 ~ Beta(2.5, 2.5), the LKJ(2) marginal for a 3×3 matrix), so m / r / a / b differ only by their restrictions. Generated quantities: `log_lik_stage1[N1]`, `log_lik_stage2[E]`, `mu_state[S,2]` (at `x_ref`), `mu_event[E,2]`, `Sigma`, `Omega`, `R_Om`. Choice sets are ragged (`cs_start`/`cs_end`), with an optional sampled-alternatives mode (`sampled = 1`, IV corrected by log(n_full/n_sub)).

## Pipeline map

| script | inputs | outputs |
|---|---|---|
| `R/compile_models.R` | `stan/*.stan` | `stan/build/*.rds` (machine-specific; gitignored) |
| `R/01_build_data.R` | peacesciencer bundled data | `data/sy_positions.rds`, `data/events.rds`, `data/stage1_panel_lag{1,3,5,10}.rds`, `data/alts_lag{h}_{all,polrel}.rds`, `data/manifest.csv`, `data/build_log.txt` |
| `R/assemble_cell.R` | the data files (manifest-checked) | in-memory `sy`, `alts`, `ev` frames for one (sample, lag, variant, centrality) cell |
| `R/02_fit_model.R` | `assemble_cell()` + compiled model | `output/fit_<cell>.rds`, `output/summary_<cell>.csv` (`--dry_run=TRUE` writes `output/standata_<cell>.rds` only) |
| `R/03_extract.R` | one `fit_<cell>.rds` | `typology_<cell>.csv` (5-bin, or 3-level for 1-D models), `mu_state_`, `mu_event_`, (M2) `ame_U/V/pairs_`; appendix `typology9_` (9-cell), `reach_`, `focus_sigma_`, `focus_states_` |
| `R/compare_models.R` | the eight fits of one (sample, lag) | `compare_loo_<group>.csv`, `compare_contrasts_<group>.csv` (material vs relational vs independent vs joint; M2 vs M1), `compare_misclass_<group>.csv` |
| `R/typology.R` | posterior `mu` draws | classification functions (sourced) |
| `R/measurement_tests.R` | fits and frames | LOO a vs b, misclassification, temporal holdout, permutation, measurement-error check (functions + CLI for LOO/misclass) |
| `R/run_measurement_real.R` | a cell's a/b fits | `output/measurement_<cell>_*.csv` (refit-based tests) |
| `R/benchmarks.R` | data files | `output/bench_dyadlogit_<cell>.*`, `output/bench_gap_<cell>.*` |
| `sim/simulate.R`, `sim/run_sim.R` | — | `sim/results/*.csv`, `sim/fits/*.rds` (gitignored), summarised in `sim/SIM_REPORT.md` |
| `sim/run_measurement_tests.R` | — | `sim/results/measurement_*.csv`, `sim/results/benchmarks_sim.csv` |
| `tests/run_tests.R` | — | runs `tests/testthat/*.R` (data coding, typology shapes, Stan log-lik == brute-force R) |
| `slurm/*.sbatch` | — | Roar job templates (nothing is submitted by the code) |

Cell naming: `<sample>_lag<h>_<variant>_<centrality>_<model>[_K<K>][_noCH][_samp<m>][_dist][_nocontig]`, e.g. `mid_lag5_all_pr_M2b_K2`.

## Data decisions (flags at the top of `01_build_data.R`)

* Events: Side A originator → Side B originator, onset year of the pair within the dispute, one direction per dispute pair, every A-originator × B-originator pair of a multilateral dispute is an event (335 such extra pairs in 2770 events). `INCLUDE_JOINERS = FALSE` (TRUE for robustness).
* Samples: `mid` (all), `fmid` (fatal: either side's fatality ≥ 1 or dispute fatality ≥ 1; negative codes → NA; events with all-NA fatality are NOT in `fmid`, 229 cases), `war` (dispute hostility level 5).
* Positions: M = log CINC share, standardised within year (PI decision 2026-10-22; the build stops and lists any state-year with a zero or missing share instead of flooring it; none exist in COW NMC 6.0, smallest share 2.4e-7; skewness of M falls from 5.31 for the raw share to −0.14). dM_h = M(t−1) − M(t−1−h), the stage-1 dM_i and x_ij all derive from this M. R = IGO PageRank (main; the eigenvector code path remains but is no longer scheduled), weights = shared IGO memberships, isolates included as vertices, standardised within year. Years with all-missing IGO data are dropped and listed in `build_log.txt` (`IGO_FILL = "drop"`; `"locf"` carries the latest IGO year forward). With peacesciencer 1.2.0 no year is all-missing (COW IGO v3 is annual in this release; 1990 is partially missing and its NA dyads are treated as no shared membership).
* Timing: positions at t−1; changes over [t−1−h, t−1]; dyadic controls (contiguity, ATOP defense/offense alliance, joint democracy polity2 ≥ 6, major-power count) at t−1; all via calendar joins.
* Conflict history: monadic and dyadic MID onset counts in the 10 years before the window starts, dyadic and monadic peace years as of window start. Main spec includes them (`--conflict_history=TRUE` in `02_fit_model.R`); counts enter as `log1p()` and are centred.
* Democracy and alliance NA → 0 with `*_missing` indicator columns kept in the data (not in the default spec).
* Choice sets: `all` = every other system member at t with positions at t−1; `polrel` = contiguous (conttype 1–5) or either state a COW major power. Events whose target is outside the choice set are dropped from stage 2 (logged; the initiator's y1 stays 1).
* W (stage 1) = initiator power changes (dM_h, dR_h), major power, democracy, log1p prior MIDs, log1p peace years (power changes added 2026-10-22 so system change can shift whether a state initiates, not only whom it targets); X (for μ) = (dM_h, dR_h) + the same controls; Z (dyadic, main) = land contiguity, alliance, joint democracy, major-power count, log1p dyadic prior MIDs, log1p dyadic peace years. `x_ref` for `mu_state` = mean X over initiating state-years in the sample.
* Alternative Z specifications (not run in the main grid): `--distance=TRUE` adds log1p(capital-to-capital distance, km; peacesciencer `add_capital_distance()`, bundled, measured at t−1) and suffixes the cell `_dist`; `--contiguity=FALSE` drops land contiguity (`_nocontig`).
* Panel: 1816–2014 (COW-MID 5.0 end). State index = all states appearing as initiators or targets (S = 217).

## Identification and estimation decisions

* Only Σ⁻¹ and Σ⁻¹μ are identified by the choice data (the −½ μ'Σ⁻¹μ term is constant across alternatives), so μ is recovered as Σ(Σ⁻¹μ); this is automatic in the parameterisation.
* The "a" model fixes corr(a_i1, a_i2) = 0 but keeps corr(c, a1) and corr(c, a2), parameterised so Ω is always PD (r2 = q·√(1−r1²)).
* AME: V column-centred (a column constant cancels in the softmax); the anchor block of U (K states with most initiations) is lower-triangular with positive diagonal; U/V scale is fixed by N(0,1) priors only. Procrustes alignment across chains in `03_extract.R` is a safety net and its rotation deviation is logged.
* Priors: β ~ N(0, 1) per element (`prior_beta_sd` is data), σ_x ~ lognormal(0, 1), sd of random effects ~ half-N(0, 1), LKJ(2) on correlations, γ, α ~ N(0, 1), τ ~ N(0, 1), α0 ~ N(0, 2.5).
* Main specification = the full all-states choice set, `--variant=all --sampled_m=0` (the defaults of `02_fit_model.R`); `polrel` and `--sampled_m` remain options for robustness and pilots. Sampled alternatives (`--sampled_m=m`) are an option with two outcome-separated samples. (a) IV rows, for every state-year, m alternatives drawn uniformly without replacement, ignoring outcomes, with IV = log-sum-exp + log(n_full/m); this estimate of the log partition sum is downward biased by Jensen's inequality, equally for initiating and non-initiating rows, and the bias shrinks with m. (b) Event rows, each event gets its own set of its chosen target plus m uniform others, independently per event, so the conditional-logit score is unbiased (McFadden 1978). Unit tests in `tests/testthat/test_sampled_alternatives.R` reproduce the review cases.
* Seed 123 everywhere.

## Typology outputs

Main output (`typology_<cell>.csv`): the 5-bin classification of posterior μ_i at `x_ref` with radius `r0` (default 0.5 SD units; peer = within r0 of the origin; bully (+,+), challenger (−,−), hybrids off-diagonal). For 1-D models the single axis gives weaker / peer / stronger (μ > r0 means the state targets others positioned below it on that axis).

Appendix outputs, written by `03_extract.R` as well. (a) `typology9_<cell>.csv`, a 9-cell classification with a peer band on each axis (`--band_M`, `--band_R`, default r0): bully, challenger, peer, material / relational opportunist, material / relational challenger, compensator A (M down, R up), compensator B (M up, R down), with posterior cell probabilities and the modal cell; D = 2 only. (b) `reach_<cell>.csv`, posterior mean and 90% interval of r_i = ||μ_i|| and of |μ_i1|, |μ_i2|. (c) `focus_sigma_<cell>.csv` and `focus_states_<cell>.csv`, descriptive only: the single pooled Σ (SDs, correlation) and each state's empirical dispersion of its chosen targets' positions around the posterior mean μ at the event (SD per axis, correlation, n_events; fewer than 5 events flagged unreliable). The model has one pooled Σ; a state-varying Σ_i would be a model change for the PI to decide on, not something these files estimate.

## Within-chain parallelism (reduce_sum)

Both programs evaluate the likelihood as `reduce_sum_static` over state-years (stage-1 term, IV log-sum-exp and the state-year's stage-2 events in one slice; `grainsize` is data, N1 for one slice = the exact serial sum). `02_fit_model.R --threads_per_chain=k` sets `rstan_options(threads_per_chain)` and `grainsize ≈ N1 / (4k)`; request `cpus-per-task = chains × k` (08_main_grid.sbatch and 04_pilot.sbatch default to 8 threads per chain, 32 and 8 cpus, 4 GB per cpu). rstan 2.32 compiles with `-DSTAN_THREADS` by default and `compile_models.R` asserts it. Generated quantities are serial, so `log_lik` is identical across thread counts; the log density and gradient agree to 1e-10 / 1e-8 (tests/testthat/test_reduce_sum.R). Profiling on the 40-state simulated cell puts ~99% of the gradient in the utility computation over alternative rows (`sim/results/profile_breakdown.csv`); on this 2-core container 2 threads gave 2.1–2.4× and 4 threads (oversubscribed) 2.7–2.9× (`sim/results/threads_timing.csv`); the Roar curve must come from the pilot. Stan `profile()` sections are in both programs for CmdStan users. Compiling M2 needs about 3 GB with `v2/env/Makevars` (no `-g`, aggressive GC) versus >5.7 GB with R's default flags; copy it to `~/.R/Makevars` where memory is tight.

## Running on Roar

Run order (submit from the repo root; every script sources `v2/slurm/00_env.sh` via `SLURM_SUBMIT_DIR`).

```
cd /storage/group/LiberalArts/default/jfe4_collab/strategic_targeting
git checkout v2-rebuild
sbatch v2/slurm/01_build_data.sbatch        # 1. compile Stan, build all data files (free account, minutes)
Rscript v2/tests/run_tests.R                # 2. unit tests (needs the compiled models)
sbatch v2/slurm/04_pilot.sbatch             # 3. 60-iteration timing pilots -> "sampling finished in" lines; set --time in 08
sbatch --export=ALL,DRY_RUN=TRUE --time=00:30:00 --mem=16G --cpus-per-task=1 v2/slurm/08_main_grid.sbatch
                                            # 4. ROAR CHECK: all 96 main cells assemble + validate their Stan data, no sampling
sbatch v2/slurm/08_main_grid.sbatch         # 5. MAIN GRID: 3 samples x 8 models x 4 lags, full choice sets, jfe4_cr_default
sbatch v2/slurm/09_main_extract.sbatch      # 6. per (sample, lag): extraction, typology (+ appendix), cross-model LOO
sbatch v2/slurm/02_fit_array.sbatch         # robustness: polrel choice set, PageRank (free account, < 48 h)
sbatch v2/slurm/05_extract.sbatch           # extraction + a-vs-b LOO for the robustness fits
sbatch v2/slurm/06_measurement.sbatch       # holdout / permutation / measurement-error refits (main cell)
sbatch v2/slurm/07_benchmarks.sbatch        # directed-dyad logit and old gap regression
```

`08_main_grid.sbatch` carries a 7-day placeholder walltime; set it from the pilot before step 5. Nothing in the repository submits jobs on its own.

A single real-data fit by hand: `Rscript v2/R/02_fit_model.R --sample=mid --lag=5 --model=M2b --K=2 --chains=4 --iter=2000 --warmup=1000` (defaults `--variant=all --sampled_m=0`). Use `--dry_run=TRUE` to validate the data assembly without sampling, `--pilot=TRUE` for a 60-iteration timing run.

Runtime expectations (from this container, 2.8 s per iteration at N2 = 39k alternative rows, 2 chains): the `polrel` cells (N2 ≈ 220k) should take roughly 5–8 h per chain at 2000 iterations; the `all` cells (N2 ≈ 1.7M) scale ~8× and will likely exceed 48 h, hence the paid partition or `--sampled_m=50`. Confirm with the pilot before committing allocation.

## Tests

`Rscript v2/tests/run_tests.R` (needs compiled models for the likelihood test). `Rscript v2/sim/run_sim.R` and `Rscript v2/sim/run_measurement_tests.R` reproduce `sim/SIM_REPORT.md`.
