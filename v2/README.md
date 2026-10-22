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

Families: `m1_latent_additive.stan` (M1) and `m2_ame.stan` (M2 = M1 + u_i'v_j, K latent dimensions, K is data). Data flag `joint`: 1 = "b" (full Σ, full Ω), 0 = "a" (Σ diagonal and corr(a_i1, a_i2) = 0; corr(c, a) stays free). Generated quantities: `log_lik_stage1[N1]`, `log_lik_stage2[E]`, `mu_state[S,2]` (at `x_ref`), `mu_event[E,2]`, `Sigma`, `Omega`, `R_Om`. Choice sets are ragged (`cs_start`/`cs_end`), with an optional sampled-alternatives mode (`sampled = 1`, IV corrected by log(n_full/n_sub)).

## Pipeline map

| script | inputs | outputs |
|---|---|---|
| `R/compile_models.R` | `stan/*.stan` | `stan/build/*.rds` (machine-specific; gitignored) |
| `R/01_build_data.R` | peacesciencer bundled data | `data/sy_positions.rds`, `data/events.rds`, `data/stage1_panel_lag{1,3,5,10}.rds`, `data/alts_lag{h}_{all,polrel}.rds`, `data/manifest.csv`, `data/build_log.txt` |
| `R/assemble_cell.R` | the data files (manifest-checked) | in-memory `sy`, `alts`, `ev` frames for one (sample, lag, variant, centrality) cell |
| `R/02_fit_model.R` | `assemble_cell()` + compiled model | `output/fit_<cell>.rds`, `output/summary_<cell>.csv` (`--dry_run=TRUE` writes `output/standata_<cell>.rds` only) |
| `R/03_extract.R` | one `fit_<cell>.rds` | `typology_<cell>.csv`, `mu_state_<cell>.csv`, `mu_event_<cell>.csv`, (M2) `ame_U/V/pairs_<cell>.csv` |
| `R/typology.R` | posterior `mu` draws | classification functions (sourced) |
| `R/measurement_tests.R` | fits and frames | LOO a vs b, misclassification, temporal holdout, permutation, measurement-error check (functions + CLI for LOO/misclass) |
| `R/run_measurement_real.R` | a cell's a/b fits | `output/measurement_<cell>_*.csv` (refit-based tests) |
| `R/benchmarks.R` | data files | `output/bench_dyadlogit_<cell>.*`, `output/bench_gap_<cell>.*` |
| `sim/simulate.R`, `sim/run_sim.R` | — | `sim/results/*.csv`, `sim/fits/*.rds` (gitignored), summarised in `sim/SIM_REPORT.md` |
| `sim/run_measurement_tests.R` | — | `sim/results/measurement_*.csv`, `sim/results/benchmarks_sim.csv` |
| `tests/run_tests.R` | — | runs `tests/testthat/*.R` (data coding, typology shapes, Stan log-lik == brute-force R) |
| `slurm/*.sbatch` | — | Roar job templates (nothing is submitted by the code) |

Cell naming: `<sample>_lag<h>_<variant>_<centrality>_<model>[_K<K>][_noCH][_samp<m>]`, e.g. `mid_lag5_polrel_pr_M2b_K2`.

## Data decisions (flags at the top of `01_build_data.R`)

* Events: Side A originator → Side B originator, onset year of the pair within the dispute, one direction per dispute pair, every A-originator × B-originator pair of a multilateral dispute is an event (335 such extra pairs in 2770 events). `INCLUDE_JOINERS = FALSE` (TRUE for robustness).
* Samples: `mid` (all), `fmid` (fatal: either side's fatality ≥ 1 or dispute fatality ≥ 1; negative codes → NA; events with all-NA fatality are NOT in `fmid`, 229 cases), `war` (dispute hostility level 5).
* Positions: M = CINC share standardised within year; R = IGO PageRank (main) or eigenvector centrality (robustness), weights = shared IGO memberships, isolates included as vertices, standardised within year. Years with all-missing IGO data are dropped and listed in `build_log.txt` (`IGO_FILL = "drop"`; `"locf"` carries the latest IGO year forward). With peacesciencer 1.2.0 no year is all-missing (COW IGO v3 is annual in this release; 1990 is partially missing and its NA dyads are treated as no shared membership).
* Timing: positions at t−1; changes over [t−1−h, t−1]; dyadic controls (contiguity, ATOP defense/offense alliance, joint democracy polity2 ≥ 6, major-power count) at t−1; all via calendar joins.
* Conflict history: monadic and dyadic MID onset counts in the 10 years before the window starts, dyadic and monadic peace years as of window start. Main spec includes them (`--conflict_history=TRUE` in `02_fit_model.R`); counts enter as `log1p()` and are centred.
* Democracy and alliance NA → 0 with `*_missing` indicator columns kept in the data (not in the default spec).
* Choice sets: `all` = every other system member at t with positions at t−1; `polrel` = contiguous (conttype 1–5) or either state a COW major power. Events whose target is outside the choice set are dropped from stage 2 (logged; the initiator's y1 stays 1).
* W (stage 1) = major power, democracy, log1p prior MIDs, log1p peace years; X (for μ) = (dM_h, dR_h) + the same controls; Z (dyadic) = contiguity (land), alliance, joint democracy, major-power count, log1p dyadic prior MIDs, log1p dyadic peace years. `x_ref` for `mu_state` = mean X over initiating state-years in the sample.
* Panel: 1816–2014 (COW-MID 5.0 end). State index = all states appearing as initiators or targets (S = 217).

## Identification and estimation decisions

* Only Σ⁻¹ and Σ⁻¹μ are identified by the choice data (the −½ μ'Σ⁻¹μ term is constant across alternatives), so μ is recovered as Σ(Σ⁻¹μ); this is automatic in the parameterisation.
* The "a" model fixes corr(a_i1, a_i2) = 0 but keeps corr(c, a1) and corr(c, a2), parameterised so Ω is always PD (r2 = q·√(1−r1²)).
* AME: V column-centred (a column constant cancels in the softmax); the anchor block of U (K states with most initiations) is lower-triangular with positive diagonal; U/V scale is fixed by N(0,1) priors only. Procrustes alignment across chains in `03_extract.R` is a safety net and its rotation deviation is logged.
* Priors: β ~ N(0, 1) per element (`prior_beta_sd` is data), σ_x ~ lognormal(0, 1), sd of random effects ~ half-N(0, 1), LKJ(2) on correlations, γ, α ~ N(0, 1), τ ~ N(0, 1), α0 ~ N(0, 2.5).
* Sampled alternatives (`--sampled_m=m`): chosen target(s) plus m uniform draws; consistent for the conditional logit; IV uses the log(n_full/n_sub) correction (approximation, documented in the Stan file).
* Seed 123 everywhere.

## Running on Roar

```
cd /storage/group/LiberalArts/default/jfe4_collab/strategic_targeting
git checkout v2-rebuild
sbatch v2/slurm/01_build_data.sbatch        # compiles Stan, builds all data (open, ~minutes)
sbatch v2/slurm/04_pilot.sbatch             # 60-iteration timing pilots (open) -> read "sampling finished in" lines
sbatch v2/slurm/02_fit_array.sbatch         # 96 polrel cells, open partition, < 48 h each
sbatch v2/slurm/03_fit_all_variant.sbatch   # 16 full-choice-set cells, jfe4_cr_default/basic (adjust --time from the pilot)
sbatch v2/slurm/05_extract.sbatch           # typology + AME extraction + LOO / misclassification a vs b
sbatch v2/slurm/06_measurement.sbatch       # holdout / permutation / measurement-error refits (main cell)
sbatch v2/slurm/07_benchmarks.sbatch
```

A single real-data fit by hand: `Rscript v2/R/02_fit_model.R --sample=mid --lag=5 --variant=polrel --model=M2b --K=2 --chains=4 --iter=2000 --warmup=1000`. Use `--dry_run=TRUE` to validate the data assembly without sampling, `--pilot=TRUE` for a 60-iteration timing run.

Runtime expectations (from this container, 2.8 s per iteration at N2 = 39k alternative rows, 2 chains): the `polrel` cells (N2 ≈ 220k) should take roughly 5–8 h per chain at 2000 iterations; the `all` cells (N2 ≈ 1.7M) scale ~8× and will likely exceed 48 h, hence the paid partition or `--sampled_m=50`. Confirm with the pilot before committing allocation.

## Tests

`Rscript v2/tests/run_tests.R` (needs compiled models for the likelihood test). `Rscript v2/sim/run_sim.R` and `Rscript v2/sim/run_measurement_tests.R` reproduce `sim/SIM_REPORT.md`.
