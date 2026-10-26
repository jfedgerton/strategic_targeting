# Strategic targeting, v2 pipeline (R + Stan)

Two-stage spatial choice model of MID initiation and target selection, estimated jointly in Stan. Branch `v2-rebuild`; the old code is reference only. A 24-cell real-data run of dd48506 is on the cluster; the changes of 2026-10-25 are for the next run.

Contents. 1 Question and design. 2 Quick start. 3 The main specification. 4 Switches. 5 Robustness presets. 6 Data decisions. 7 Model. 8 Outputs. 9 Archived options. 10 Changelog.

## 1. Question and design

How do states choose whom to attack? Two axes of power enter together, material (CINC share) and relational (position in the IGO network), and the paper's typology reads a state's preferred target off both: a **bully** targets states below it on both axes, a **challenger** targets states above it, a **peer** targets states near it, and the off-diagonal hybrids mix the two.

Two stages, estimated jointly in one Stan program per family.

| stage | unit | outcome | regressors |
|---|---|---|---|
| 1, initiation | state-year (i, t) | y1 = 1 when i initiates at least one non-escalation event in t | W (initiator power changes, levels, controls) and the inclusive value IV of stage 2 |
| 2, target choice | event (i → j, t), choice set = every other system member at t | the chosen j | the gap x_ij = (M_i − M_j, R_i − R_j) at t − 1 through a quadratic utility around an ideal gap μ_n = μ0 + X_n β + a_i, dyadic controls Z, target intercept b_j |

Two model families and four variants each (eight model codes). Stage 1 is identical in all eight; only the stage-2 utility and hence the IV differ, so LOO compares like with like.

| code | `D_mode` | `joint` | stage-2 targeting space |
|---|---|---|---|
| M1m / M2m | 1 | – | material only (x = dM, scalar μ, scalar Σ, Ω 2×2 over (c, a)) |
| M1r / M2r | 2 | – | relational only (x = dR) |
| M1a / M2a | 3 | 0 | both axes, independent (Σ diagonal, corr(a_1, a_2) = 0; corr(c, a_d) free) |
| M1b / M2b | 3 | 1 | both axes, joint (full Σ, full Ω) |

M1 (`stan/m1_latent_additive.stan`) has state random effects (c_i, a_i) and target intercepts b_j; M2 (`stan/m2_ame.stan`) adds the AME term u_i'v_j with K latent dimensions. The main grid is 3 samples × 8 models × 4 lags = 96 cells.

### Code map

| path | role |
|---|---|
| `R/spec.R` | THE main specification, the switches, their suffixes, the sampler settings per family, the robustness presets. Every script sources it first. |
| `R/helpers_data.R`, `R/helpers_model.R` | functions only: data-side (coding rules, dates, direction, weights, manifest) and model-side (Stan data assembly, anchors, likelihood in R, fit safety, provenance) |
| `R/01_build_data.R` | builds every data file into `data/` (or `data_sideA/`) with `manifest.csv` (nrow, sha256 per file) and `build_log.txt` |
| `R/assemble_cell.R` | `assemble_cell()` / `assemble_cell_spec()`: the sy / alts / ev frames of one cell, manifest-checked; the only place the switches act at assembly |
| `R/02_fit_model.R` | one cell: assemble, validate, fit (or `--dry_run`), save `fit_<cell>.rds` and `summary_<cell>.csv` |
| `R/03_extract.R` | one fit: typology, μ summaries, AME alignment, appendix tables |
| `R/04_compare.R` | one (sample, lag) group: LOO (three flavours), joint vs independent, M1 vs M2 |
| `R/05_robustness.R` | the robustness presets (refits with a switch flipped; permutation and measurement-error tests) |
| `R/typology.R`, `R/measurement_tests.R` | sourced function files (classification; LOO, holdout, permutation, error inflation) |
| `R/compile_models.R` | compiles the two Stan programs into `stan/build/` with source-hash sidecars |
| `slurm/` | `00_env.sh`, `01_build.sbatch`, `02_main_grid.sbatch`, `03_post.sbatch`, `04_robustness.sbatch`, `05_simulations.sbatch` |
| `sim/` | `simulate.R` (functions), `run_sim.R` (one driver, `--task=`), `write_sim_report.R`, `SIM_REPORT.md` (superseded until the rerun) |
| `tests/` | `Rscript v2/tests/run_tests.R` runs `tests/testthat/*.R` |
| `archive/` | everything replaced on 2026-10-25, each file with a one-line header saying what replaced it (section 9) |

## 2. Quick start

Three commands. Locally (from the repo root; needs R 4.3 / rstan 2.32, peacesciencer 1.2.0):

```
Rscript v2/R/compile_models.R && Rscript v2/R/01_build_data.R        # 1. build (about 3 min)
Rscript v2/R/02_fit_model.R --sample=fmid --lag=3 --model=M1b          # 2. one main cell (add --dry_run=TRUE to validate only)
Rscript v2/R/03_extract.R --fit=v2/output/fit_fmid_lag3_all_pr_M1b.rds # 3. post: typology; 04_compare.R needs the eight fits of a group
```

On Roar (submit from the repo root; `00_env.sh` is sourced by every job; account `jfe4_cr_default`, partition `basic`, for the fits; the `open` account for the build, dry runs and post-processing):

```
cd /storage/group/LiberalArts/default/jfe4_collab/strategic_targeting && git checkout v2-rebuild
sbatch v2/slurm/01_build.sbatch                                                     # 1. compile + build (open account, minutes)
Rscript v2/tests/run_tests.R                                                        #    unit tests (needs the compiled models)
sbatch --export=ALL,DRY_RUN=TRUE --time=00:30:00 --mem=16G --cpus-per-task=1 v2/slurm/02_main_grid.sbatch   # Roar check: 96 cells assemble, no sampling
sbatch v2/slurm/02_main_grid.sbatch                                                 # 2. MAIN GRID (jfe4_cr_default): 96 cells, array 0-95
sbatch v2/slurm/03_post.sbatch                                                      # 3. extraction + comparison per (sample, lag), array 0-11
```

Robustness: `sbatch --export=ALL,PRESET=noww2 v2/slurm/04_robustness.sbatch` (section 5). Simulations: `sbatch v2/slurm/05_simulations.sbatch` (open account). A timing pilot before setting the walltime of `02_main_grid.sbatch` (a 7-day placeholder): `Rscript v2/R/02_fit_model.R --sample=mid --lag=5 --model=M2b --pilot=TRUE --threads_per_chain=8` (1 chain × 60 iterations, saved as `fit_<cell>_pilot.rds`, which every production tool refuses); read "sampling finished in" from the log. Nothing in the repository submits jobs on its own.

The array index of `02_main_grid.sbatch` maps to (sample, model, lag) with lag slowest, then model, then sample: index = 24 × lag_index + 3 × model_index + sample_index, samples (mid, fmid, war), models (M1m, M1r, M1a, M1b, M2m, M2r, M2a, M2b), lags (1, 3, 5, 10). It is written out as an explicit case table in the script. Within-chain parallelism: `--threads_per_chain=k` (reduce_sum over state-years; `cpus-per-task = chains × k`, 8 threads and 32 cpus by default; the log density is identical across thread counts). Compiling M2 needs about 3 GB with `env/Makevars`.

## 3. The main specification

Defined once, in `R/spec.R` (the defaults of `SPEC_CHOICES` and the constants below it). Every switch at its default, cell names `<sample>_lag<h>_all_pr_<model>[_K2]`, e.g. `fmid_lag3_all_pr_M1b`, `war_lag3_all_pr_M2b_K2`.

| element | main specification |
|---|---|
| samples | `mid` (all MIDs), `fmid` (fatal, main coding), `war` (pair-level war); appendix samples `fmidcow`, `fmidub` can be fitted by hand |
| lags h | 1, 3, 5, 10 (change window [t − 1 − h, t − 1]) |
| choice set | every other system member at t with positions at t − 1, full sets (no sampling), the already-fighting states removed (two sets per state-year, section 6) |
| events | Side A × Side B pairs including joiners; war-level disputes follow COW's war initiator (`direction = cowwar`); escalations excluded; same-year duplicates collapsed |
| stage-2 weights | w_e = 1 / (the initiator's targets in the same dispute and onset year) |
| positions | M = log CINC share, R = log IGO PageRank; gaps divided by one constant per axis (3.268, 0.524) |
| changes | year-centred, pooled-SD scaled; raw (not winsorised) |
| **W** (stage 1, P1 = 8) | `dM_i`, `dR_i`, `M_lvl`, `R_lvl`, `ongoing_tm1`, `log_prior_mids`, `dem_tm1`, `log_n_system_tm1` |
| **X** (μ shifters, P = 2) | `dM_i`, `dR_i` |
| **Z** (dyadic, Pz = 7) | `contig_land`, `joint_dem`, `log_dy_prior_mids`, `dy_mid_tm1`, `dM_j_h`, `dRpr_j_h`, `dch_j_missing` (the last only when non-constant in the cell, which it is in every main cell) |
| centring | every W / X control except `dM_i`, `dR_i`, `M_lvl`, `R_lvl` is grand-centred (0/1 indicators included); the continuous Z columns are centred over the stage-1 IV rows; `dy_mid_tm1` and the indicators stay 0/1 |
| x_ref | mean X over the initiating state-years of the sample (`mu_state` is μ at x_ref) |
| M2 | K = 2 latent dimensions; anchors = the K states with most initiations |
| sampler | per family in `SAMPLER`: M1 and M2 4 chains × 2000 iterations (1000 warmup), adapt_delta 0.9, max_treedepth 10; seed 123 |
| typology | r0 = 0.5 scaled gap units |

Lag-3 dry runs of the main specification (6fd5362 and the restructured code, identical): S = 217, N1 = 14,597, N2 = 1,720,110 (1,719,835 IV rows, 1,719,990 stage-2 rows); mid E = 3258, sum(y1) = 2086, weights 2694.7; fmid E = 1055, 695, 741.0; war E = 432, 203, 222.6.

## 4. Switches

Every switch is an option of `02_fit_model.R`, `04_compare.R` and `05_robustness.R` as `--name=value`; the build takes `--direction` only. Every non-default value appends its suffix to the cell name in the order of this table (pairs, ww2, controls, weights, changes, direction), so variants never overwrite main outputs: `war_lag3_all_pr_M2b_K2_noww2_sideA`. Defaults = the main specification.

| switch | values (default first) | effect | suffix | acts in |
|---|---|---|---|---|
| `pairs` | `all` / `originators` / `joiners` | event pairs kept: all; `is_joiner_pair` FALSE only; TRUE only. y1 is recomputed from the kept events exactly as `ww2 = drop` does, so a state-year whose only events were in the other group gets y1 = 0 (counts logged). Joiners-only is a **stage-2 comparison**: its stage 1 counts a state as initiating only when it joined. After the duplicate collapse `is_joiner_pair` is FALSE for a key with any originator pair. | `_orig` / `_join` | assembly |
| `direction` | `cowwar` / `sideA` | initiator of war-level disputes = COW's war-initiating side (reads `data/`) / MID Side A everywhere (reads `data_sideA/`, built by `01_build_data.R --direction=sideA`) | `_sideA` | build + assembly |
| `ww2` | `include` / `drop` | the 16 WWII disputes' events removed from stage 2 and their initiations from stage 1, through the build's `in_<sample>_noww2` flags | `_noww2` | assembly |
| `weights` | `dispute` / `none` | stage-2 composite-likelihood weights w_e / weight 1 on every event | `_w1` | assembly |
| `controls` | `minimal` / `full` | the main W / X / Z above / the full base specification of commit 8bf0ce1 (W = X = `dM_i`, `dR_i`, `cowmaj_tm1`, `dem_tm1`, `log_prior_mids`, `log_peace_years`; Z = `contig_land`, `alliance_def_off`, `joint_dem`, `major_count`, `log_dy_prior_mids`, `log_dy_peace_years`, `dM_j_h`, `dRpr_j_h`, [`dch_j_missing`]; P1 = P = 6, Pz = 9; no `ongoing_tm1`, `dy_mid_tm1`, levels or system size) | `_fullctrl` | assembly |
| `changes` | `raw` / `winsor` | `dM_i`, `dR_i`, `dM_j_h`, `dRpr_j_h` clipped at the 1st / 99th percentiles (`WINSOR_PROBS`, fixed) of the stage-1 panel's `dM_h` / `dRpr_h` for the cell's lag (the target change is the same state-year change table, so one pair of cut points per axis), applied at assembly before the Z centring; a target without a change value stays 0; cut points and clipped counts logged | `_winsor` | assembly |
| `ame_K` | `2` | latent dimensions of M2 | `_K<K>` (always present for M2) | fit |

Sampler settings (`--chains --iter --warmup --adapt_delta --max_treedepth`) override `SAMPLER` per call and do not change the cell name; `--pilot=TRUE` adds `_pilot`.

Lag-3 fmid M1b dry runs of each switch (main: E = 1055, N1 = 14,597, P1 = 8, P = 2, Pz = 7, sum(y1) = 695):

| variant | E | N1 | P1 | P | Pz | sum(y1) | note |
|---|---|---|---|---|---|---|---|
| main | 1055 | 14,597 | 8 | 2 | 7 | 695 | weights sum 741.0 |
| `pairs=originators` | 545 | 14,597 | 8 | 2 | 7 | 510 | 1144 → 596 sample events (548 joiner pairs removed); 185 state-years whose fmid events were all joiner pairs get y1 = 0 |
| `pairs=joiners` | 510 | 14,597 | 8 | 2 | 7 | 251 | 1144 → 548 sample events; stage-2 comparison (444 state-years get y1 = 0) |
| `direction=sideA` | 1063 | 14,597 | 8 | 2 | 7 | 713 | reads `data_sideA/` (N2 = 1,720,113; weights sum 750.2) |
| `ww2=drop` | 941 | 14,597 | 8 | 2 | 7 | 670 | 122 WWII events removed; 25 state-years get y1 = 0 |
| `weights=none` | 1055 | 14,597 | 8 | 2 | 7 | 695 | weights sum 1055 instead of 741.0 |
| `controls=full` | 1055 | 14,597 | 6 | 6 | 9 | 695 | W = X with the six base controls, Z with alliance and major count |
| `changes=winsor` | 1055 | 14,597 | 8 | 2 | 7 | 695 | cut points (lag-3 panel, 14,597 state-years) dM_h [−3.065, 3.167], dRpr_h [−1.964, 4.341]; clipped dM_i 292, dR_i 292 of 14,597; dM_j_h 27,191, dRpr_j_h 15,007 of 1,720,110 rows |

## 5. Robustness presets

`05_robustness.R --preset=<name>` (`ROBUSTNESS_PRESETS` in `spec.R`), run by `04_robustness.sbatch` with `--export=ALL,PRESET=<name>`.

| preset | switch flipped | what it refits | post-processing |
|---|---|---|---|
| `noww2` | `ww2 = drop` | the cell list (default `ROBUSTNESS_CELLS`: lag 3, M1b and M2b, mid / fmid / war; `--cells=mid:3:M1b,...` overrides, `--cell_index=i` runs one) through `02_fit_model.R`, so every fit carries the suffix | `03_post.sbatch` with `--export=ALL,SWITCHES="--ww2=drop"`, or `04_compare.R --ww2=drop` |
| `sideA` | `direction = sideA` | same; needs `data_sideA/` (`sbatch --export=ALL,DIRECTION=sideA v2/slurm/01_build.sbatch`) | `SWITCHES="--direction=sideA"` |
| `noweights` | `weights = none` | same | `SWITCHES="--weights=none"` |
| `fullcontrols` | `controls = full` | same | `SWITCHES="--controls=full"` |
| `winsor` | `changes = winsor` | same | `SWITCHES="--changes=winsor"` |
| `perm_error` | none | the permutation test and the systematic measurement-error (relational-level inflation) tests on the main cell (`ROBUSTNESS_MAIN_CELL`: mid, lag 5), per family, from its a / b fits; `--tests=holdout` adds the temporal holdout; `--perm_mode=refit` gives the calibrated p-value | writes `output/measurement_<cell>_*.csv` |

Other switches on the command line are passed through (`--preset=noww2 --controls=full` fits `_noww2_fullctrl` cells). The driver refuses to refit main-spec cells. The `pairs` switch is not a preset; run its cells through `02_fit_model.R --pairs=originators|joiners` directly.

## 6. Data decisions

Flags of the build are at the top of `01_build_data.R`; the rules live in `helpers_data.R` and are unit-tested on toy disputes. Counts are from the rebuilt `data/build_log.txt`.

### 6.1 Events and joiners (PI decision 2026-10-24)

An event is a Side A participant → Side B participant pair of a COW-MID dispute, joiners included, with onset year equal to the first year that pair appears in the dispute (the sides are those of that earliest year; `build_pair_table()`, P3 fix 2026-10-25: France and the USSR are opposed in dispute 258 in 1939 and again with reversed sides in 1941, the only such pair), one direction per dispute pair, and every A × B pair of a multilateral dispute is an event. `INCLUDE_JOINERS = TRUE` is the main build (FALSE reproduces the originator-only base of 2770 events). With joiners there are 3720 events, of which 950 are joiner pairs (`is_joiner_pair`, TRUE when either side has `orig = 0`; kept for tabulation and for the `pairs` switch, never a covariate). The war sample grows most, from 127 originator pairs to 552 pairs, because the two world wars contribute most of their pairs through joiners. The fatal-MID rules apply unchanged to joiner pairs (with the pair-level war coding the joiner pairs are coded cow 897, mic_fill 50, uncoded 3; the war rule no longer adds any). The conflict-history controls count every participant including joiners (`INCLUDE_JOINERS` governs both).

### 6.2 Initiator direction (PI decision 2026-10-25, D1; the `direction` switch)

COW-MID Side A is the side that took the FIRST codable action, not the aggressor. In dispute 258 (the Second World War) Side A is Poland and the Allies, so the Side A rule produced Poland → Germany, UK → Germany and USA → Germany in 1939–41 and never Germany → Poland. For wars the standard aggressor coding is COW's war initiator, so for disputes at hostility level 5 that can be linked to a COW inter-state war (`cow_war_inter`, which has the war sides `sidea1 / sidea2` and the state-level initiator flags; every COW war has exactly one initiating side, asserted) the direction follows COW's war-initiating side. The rule has two steps. (a) Pair level. A MID pair (i, j) with onset year t is matched to the COW war dyad-years (i, j) in t or t + 1 (`WAR_DIR_WINDOW = 1`; MID onsets precede war onsets, the crisis year is often the year before the fighting reached war level). Every matched row votes: reverse when j's war side initiated, keep when i's did; unanimous votes decide the pair (`dir_source = "cow_war_pair"`, whether reversed or confirmed), mixed votes (two wars with opposite initiators) leave Side A. On the data 363 of the 552 war-level pair events match (310 at t, 53 at t + 1 only; t + 2 would add 5, any later ongoing year 29, not used); 160 are reversed, 201 confirmed, 2 ambiguous (dispute 611 / 1360, Vietnam–Cambodia, matched to wars 163 and 176). (b) Dispute level. The unmatched pairs of a war-level dispute follow the dispute's side mapping when it is unambiguous: the dispute is linked to the COW war that a strict majority of its pair-decided events matched, at least `WAR_DIR_MIN_LINK = 2` of them; the war side of the MID Side A state over those events must again be a strict majority of at least 2; if that side is the war's non-initiating side the unmatched pairs are reversed (`"cow_war_dispute"`), else kept. 34 disputes are linked with an unambiguous mapping (no dispute had a split mapping); 32 unmatched pairs are reversed and 93 kept. Everything else keeps Side A (`dir_source = "sideA"`): the 3168 events of non-war disputes, 14 events of 4 war-level disputes with no pair match, 50 events of 18 disputes whose link rests on a single pair-decided event (below the minimum; with a minimum of 1 they would be 2 reversals and 48 confirmations, logged, not applied) and the 2 ambiguous pairs. In total 192 of 3720 events (35% of the war-level pairs) are reversed in 24 disputes, 160 by the pair match and 32 by the dispute mapping; the originator status, fatality and hostility codes move with the state, so `is_joiner_pair`, `war` and the fatal codings are unchanged. `events.rds` carries `dir_flipped`, `dir_source`, `dir_rule` (the reason), `war_link` (the COW war number) and `n_war_match`; `sidea1 = 0` marks a reversed row. Everything downstream (escalation flags, reference dates, choice sets, `fighting_at_ref`, the duplicate collapse, y1, the no-WWII flags, the event weights) is computed from the applied direction in every sample; no collapsed key mixes a reversed and an unreversed source. The build stops unless Germany → Poland 1939 exists and Poland → Germany 1939 does not. The appendix rebuild `Rscript v2/R/01_build_data.R --direction=sideA` applies the old coding everywhere and writes to `v2/data_sideA/` (its own manifest and build_log), so the two codings never overwrite each other; `--direction=sideA` in the fit, compare and robustness scripts reads that directory and suffixes the cells with `_sideA`.

Reversed disputes (dispute, COW war, reversed pairs of its war-level pairs, by pair match / dispute mapping; the full event list is in `build_log.txt` under "direction (D1)"):

| dispute | COW war | reversed | pair | dispute |
|---|---|---|---|---|
| 19 | 10 Austro-Sardinian 1848 | 3 of 3 | 3 | 0 |
| 41 | 118 Manchurian (Sino-Soviet) 1929 | 1 of 1 | 1 | 0 |
| 57 | 22 Crimean 1853 | 5 of 5 | 4 | 1 |
| 111 | 127 Conquest of Ethiopia 1935 | 1 of 2 | 1 | 0 |
| 180 | 85 Russo-Japanese 1904 | 1 of 2 | 1 | 0 |
| 183 | 136 Nomonhan 1939 | 2 of 2 | 2 | 0 |
| 184 | 133 Changkufeng 1938 | 1 of 1 | 1 | 0 |
| 199 | 160 Assam 1962 | 1 of 1 | 1 | 0 |
| 200 | 155 Sinai 1956 | 8 of 8 | 3 | 5 |
| 202 | 67 Sino-French 1884 | 1 of 1 | 1 | 0 |
| 235, 339, 414, 2725, 3724 | 139 World War II (run-ups to entry) | 1 each | 1 | 0 |
| 258 | 139 World War II 1939 | 122 of 122 | 100 | 22 |
| 261 | 55 Seven Weeks 1866 | 24 of 24 | 24 | 0 |
| 320 | 106 World War I 1914 | 3 of 3 | 3 | 0 |
| 396 | 106 World War I 1914 | 1 of 1 | 1 | 0 |
| 1027 | 124 Chaco 1932 | 1 of 1 | 1 | 0 |
| 1035 | 169 Six Day 1967 | 6 of 6 | 3 | 3 |
| 1202 | 91 Fourth Central American 1907 | 3 of 3 | 2 | 1 |
| 1435 | 163 Vietnam War 1965 | 1 of 1 | 1 | 0 |
| 3826 | 139 World War II 1939 | 2 of 2 | 2 | 0 |

War names are hand-entered from the COW Inter-State War Data v4.0 list (`COW_WAR_NAMES` in `helpers_data.R`), for the log only. Effects on the samples (lag 3 dry runs, dd48506 → this build, with the D5 date change): mid E 3264 → 3258 and stage-1 initiations 2096 → 2086; fmid 1063 → 1055 and 712 → 695; war 436 → 432 and 215 → 203. Reversing a pair keeps the pair, so E moves only through the escalation flags and `fighting_at_ref` (D5 and the reversed reference dates), whereas the initiations fall because the WWII and Seven Weeks reversals turn many single-initiator state-years (each Ally in 1939–41, each German state in 1866) into one initiator (Germany, Prussia) with many targets. In the war sample 175 of 464 events (38%) are reversed. The 24 war-level disputes kept as Side A are listed in the build log with their reason.

### 6.3 Participant-level dates (PI decision 2026-10-25, D5)

peacesciencer 1.2.0 bundles `gml_part`, participant-level start and end year / month / day (the day is unknown for about 10%) for 2165 of the 2436 COW dispnums (5217 participant spells). Where the dispnum matches, `pair_dispute_intervals(part = gml_part)` dates a pair spell from its participants: entry = the later of the two participants' starts, exit = the earlier of their ends, both exact months, when each participant has exactly one gml_part spell overlapping the pair's COW years and the dates are consistent with those years (a GML entry before the pair's first COW year is clamped to January of it, since the pair was in the dispute by then; a GML exit after the last COW year to December of it; an entry after the first COW year or an exit before the last one contradicts the COW years and keeps the dispute-level bounds, as do participants with no or several overlapping spells). 5874 of the 7456 directed spells are dated this way (5634 directly, 84 with the entry clamped, 152 with the exit clamped, 4 both); 1582 keep the dispute-level bounds (774 in disputes absent from gml_part, 350 with a participant lacking an overlapping spell, 132 with several, 326 inconsistent years or entry after exit; `date_source` and `part_reason` in the interval table, counts in the build log). An unknown entry bound is also clipped by the feasible exit (688 directed spells had start_hi > end_hi under the dispute-level bounds; 216 still needed the clip after the participant dates). Before / after (dispute-level bounds / participant dates): exact entry months 5536 / 6846 and exact exit months 486 / 5988 of 7456 spells; escalations 199 / 215 (38 new, 22 dropped); undecided (event, other dispute) comparisons 230 / 71; events whose decision rests on a year-rule fallback 184 / 64; `fighting_at_ref` events 31 / 28; already-fighting (i, j, t) keys at the reference date 3687 / 3669; same-year opposite-direction pairs among the non-escalation events (both i → j and j → i in t) 118 / 113. Known limitation: comparisons stay at month resolution, and an end in the reference month counts as still active, so Israel → Egypt 1967 (dispute 1480, GML start June 29) is now an escalation of the Six-Day War dispute 1035 (GML end June 10 / 17) without a fallback, where the day fields would separate them; day resolution is an unimplemented suggestion.

### 6.4 Escalation and the already-fighting rules (PI decisions 2026-10-24, month resolution)

Stage 2 is target selection, so a new dispute against a state the initiator is already fighting is escalation, not selection. An event is an escalation when i and j are both participants in a different dispute that started before the event's dispute start date and had not ended by it, compared at month resolution. The year-level version of this rule (pair in a dispute active in year t with onset before t) mislabelled about 45% of the escalations (China → USA 1950, the Korean War, was dropped because a short 1949 incident had a 1950 row) and missed re-entries within a year. The fallback dates come from `cow_mid_disps` and are dispute-level, because that table ships the dispute start month and year (`stmon`, `styear`) and the duration in days (`mindur`, `maxdur`) but no end month and no participant-level dates (`cow_mid_dirdisps` has years only); participant-level months come from `gml_part` where available (6.3). The end date is therefore a bound, not a point, between the first day of the start month plus `mindur` and the last day of the start month plus `maxdur`, clamped to the dispute's last year in the directed dispute-years (198 of 2436 disputes have an exact end month, 1954 a two-month span, 284 wider, 10 inconsistent with the duration and given the whole last year). A pair that enters a dispute after its first year has an unknown entry month (the whole year), a pair that leaves before the dispute's last year an unknown exit month, and a pair that leaves and re-enters has separate spells. Every comparison is three-valued. When it is undecided the year rule decides and the fallback is counted. Escalations stay in `events.rds` with `escalation = TRUE` (`escalation_year` keeps the year rule, `esc_fallback` marks fallback decisions, `ev_start_lo` / `ev_start_hi` the event's date interval as a month index) and are excluded from every sample flag, so they are in neither stage 2 nor y1. Result with the participant dates: 215 escalations.

Joiner entry month (Codex review P1-2, 2026-10-24). A pair with a joiner on either side (`orig == 0`) whose spell begins in the dispute's first year used to receive the dispute's exact start month although the joiner's entry month is unknown. Its entry interval is [dispute start month, December of that year] with `start_exact = FALSE`; originator pairs keep the exact start, later entrants keep the whole entry year. On the COW data 1154 of 7456 directed spells (577 undirected pairs) are affected; escalations fell from 217 to 199 under the dispute-level bounds. The build rebuilds the pre-review intervals only to log these before / after counts.

Two choice sets per state-year (Codex review P1-1, PI choice "A′", 2026-10-24). C_n^{iv}, the stage-1 IV set, drops alternative j when i and j are in a dispute active at January of t that started before it, for every state-year alike, so the stage-1 regressor IV_n uses no information about the state-year's own events (the previous rule used the event date for initiators and January for non-initiators, which made the IV set a function of y1). C_n^{s2}, the stage-2 set, drops j when i and j are in a different dispute active at the state-year's reference date, the start of its earliest non-escalation event when it initiates and January of t otherwise, that started before it. Both rules are the month rule with the year fallback, and nothing after the respective date enters. For a non-initiating state-year the two sets are identical. The alternatives files keep one row block per state-year, the union of both sets, with the flags `in_iv` and `in_s2` (`choice_set_flags()`); rows in neither set are removed; the Stan programs sum IV over the `in_iv` rows and run the softmax over the `in_s2` rows. At h = 3, 1,720,110 rows remain (1,719,835 IV rows, 1,719,990 stage-2 rows). Because the stage-2 set is fixed at the reference date, a later event of the same state-year whose pair was still fighting at that date (their dispute ended between the reference date and the event) has no stage-2 row (it may keep an IV-only row with `in_s2 = 0`). There are 28 such events (7 fmid, 1 war), flagged `fighting_at_ref`. They keep every sample flag, so y1 is unchanged, and leave stage 2 only, like events whose target lacks positions; the build lists them and stops if any other non-escalation event loses its stage-2 row. This is a deviation from the requirement that every non-escalation event's chosen alternative be in its choice set, which cannot hold with one stage-2 set per state-year, and is flagged for the PI. The build asserts that the IV flag equals the January rule on every retained row and that only initiating state-years have differing sets, and stops if any state-year with rows has an empty IV set or any initiating state-year a stage-2 set that is empty or holds only its chosen target (`assemble_cell()`, `build_stan_data()` and the Stan programs repeat these checks). Helpers `dispute_date_bounds()`, `pair_dispute_intervals()`, `month_rule()`, `flag_escalations_month()`, `reference_dates()`, `already_fighting_keys_month()` (reference `NULL` = January for every state-year).

Same-year duplicates (PI decision 2026-10-24). After the escalation exclusion, events with the same (initiator, target, year) from different disputes collapse into one event (`collapse_duplicate_events()`). Sample flags and fatality codings are ORed, `dispnum` keeps the smallest dispute number, `dispnums` lists them all, `n_disputes` counts them, `is_joiner_pair` is FALSE when any collapsed row was an originator pair, the per-sample weights are summed and capped at 1, and escalation rows are untouched. One chosen alternative per (i, j, t) is asserted in the build and in `assemble_cell()`. The no-WWII sample flags `in_<sample>_noww2` (Codex review P2-3) are the sample flags of the non-WWII source rows alone, ORed in the same collapse, so a key that mixes a WWII and a non-WWII dispute stays in the `ww2 = drop` sample only when a non-WWII dispute qualifies (UK → Bulgaria 1941, disputes 258 and 3871, stays in the no-WWII `mid` sample and leaves the no-WWII `fmid` and `war` samples; the other five collapsed WWII keys are WWII-only).

WWII (`ww2 = drop`). `is_ww2` marks the events of 16 COW-MID disputes, 258 (the Second World War itself) plus every war-level dispute that starts between 1939-09 and 1945-09, whose highest action is 21 (join interstate war) and whose dyad appears in COW inter-state war 139 in `cow_war_inter` (235, 339, 414, 418, 518, 2725, 3701, 3705, 3717, 3724, 3813, 3820, 3822, 3826, 3827, the run-ups to a state's entry into the war). Disputes 157, 179, 183 and 613 are wars of their own (COW wars 130, 142, 136, 145) and stay. The rule is `ww2_disputes()` in `helpers_data.R`, the list is asserted against `WW2_DISPNUMS`, and the build log prints every candidate with its reason. 134 of the war events are WWII; `ww2 = drop` removes them from stage 2 and their initiations from stage 1 using the `in_<sample>_noww2` flags, so a state-year whose only events were WWII gets y1 = 0.

Lagged involvement at t − 1 (PI decision 2026-10-24), always from every participant including joiners. `ongoing_tm1` (stage 1) is 1 when the state was in any MID active in year t − 1 (mean 0.31 at h = 3). `dy_mid_tm1` (stage 2) is 1 when i and j were both in a MID active in t − 1; `assemble_cell()` stops if it is constant in a cell.

### 6.5 Fatal and war coding

Samples: `mid` (all), `fmid` (fatal, main coding), `fmidcow` (fatal, COW-only coding; appendix), `fmidub` (fatal, upper bound; appendix), `war` (pair-level war).

* War is coded at the pair level (PI decision 2026-10-24, revised the same day). An event is a war event only when the dispute is war-level (hostility level 5) and at least one participant of the pair reached hostility level 5 itself (participant `hostlev1` / `hostlev2` in `cow_mid_dirdisps`, maximum over the pair's years). The dispute-level coding is kept as `war_disp`. Of the 552 dispute-level war pairs 508 are war pairs (44 had neither participant at level 5, mostly peripheral Gulf War and Iraq 2003 coalition members). A first version required both participants at level 5 (329 events); the PI chose the less strict rule. The war rule in the fatal coding uses the pair-level definition.
* COW fatality: an event is fatal if either side's fatality ≥ 1 or the dispute fatality ≥ 1; negative codes → missing. 229 of 2,770 originator events have every COW fatality code missing.
* `fmid` (PI decision 2026-10-23) applies four rules in order. (a) Wars are fatal by definition, since COW wars require 1,000+ battle deaths. (b) A coded COW fatality value is kept as is. (c) When every COW code is missing, the event takes the MIC 1.0 value (Gibler et al., Militarized Interstate Confrontations; confrontation-level `fatalmin` ≥ 1; MIC `micnum` = COW `dispnum`; `data_raw/mic-conf-1.0.csv`, sha256-checked against `MIC_CONF_SHA256`, see `data_raw/SOURCE.md`). (d) When COW and MIC are both missing, the event takes the GML MID 2.2.1 coding of the same directed pair (peacesciencer `gml_mid_ddydisps`; same rule as COW). No source overrides an earlier coded value. Where both COW and MIC code fatality they disagree on 120 events, logged, not used. `fatal_source` in `events.rds` records which rule set each event.
* `fmidcow` is the earlier rule (missing → not fatal, no war rule). `fmidub` is the bounds check, `fmid` plus the events no source codes, counted as fatal; if `fmid` and `fmidub` agree, the uncoded events do not drive the fatal-MID findings.

Events after the collapse: mid 3417, fmid 1144, fmidcow 979, fmidub 1162, war 464. Stage-1 initiations at h = 3: mid 2086, fmid 695, war 203.

### 6.6 Positions and changes (PI decisions 2026-10-23)

Positions. M = log CINC share and R = log IGO PageRank, both left in raw log units. CINC share and PageRank each sum to one over the system in a year, so M is the log share of world capabilities and R the log share of network centrality, and the stage-2 gaps x_ij = (M_i − M_j, R_i − R_j) are log ratios between the two states. The two axes have very different spreads (SD about 3.27 log points for M and 0.52 for R), so each gap is divided by one constant per axis, the same for every year and lag: the SD of the within-year gap over all ordered pairs of distinct states (`gap_scale()`), 3.268 for M and 0.524 for R, stored in `change_scale.csv` with h = NA; the raw log ratios are kept as `dM_logratio` and `dR_pr_logratio`. One model unit on the material axis is therefore a 26-fold ratio of CINC shares, and on the relational axis a 1.7-fold ratio of PageRank. No script converts estimates back; multiply a μ, Σ or r0 value by the constant for its axis by hand. The build stops and lists any state-year with a zero or missing CINC share or a non-finite or non-positive PageRank (none exist). PageRank weights are shared IGO memberships with isolates included as vertices. Years with all-missing IGO data are dropped (`IGO_FILL = "drop"`; none with peacesciencer 1.2.0; 1990 is partially missing and its NA dyads are treated as no shared membership). Why: on within-year z-scores a state's change moved with the growth of the system (correlation with the 3-year log growth in the number of states about +0.29 for M and +0.17 for R at h = 3) and the changes had SD 0.11 and 0.25, so the per-unit priors dominated.

Changes. One state-level table per lag, `changes_lag{h}.rds` (`build_change_table()`). For every state in the system at t with positions at t − 1 and t − 1 − h, the raw change is level(t − 1) − level(t − 1 − h) in log units (calendar joins, never row lags; `dM_raw_h`, `dRpr_raw_h`). The transformed change `dM_h`, `dRpr_h` subtracts the mean raw change of all states that have a change in that year and divides by one pooled constant per lag and variable, the SD of the year-centred change over all rows of that lag (1 unit = 0.19 log points for M and 0.15 for R at h = 3; constants in `change_scale.csv`). Zero means the state changed like the average state over the same window. The same table supplies the initiator's change in stage 1 (`dM_i`, `dR_i` in W and X) and the target's change in stage 2 (`dM_j_h`, `dRpr_j_h` in Z for every model). A target with no change value (typically a recent entrant) gets 0 and `dch_j_missing = 1` (62k of 1.72 million rows at h = 3). The build asserts finiteness, within-year mean 0, pooled SD 1 and |cor| < 0.02 with the system-growth variable.

Power levels in stage 1 (PI decision 2026-10-25, D2 and D3). `M_lvl` and `R_lvl` are the within-year centred log CINC share and log PageRank at t − 1 (centred over every state with positions in that year) divided by the gap constants, so their units match the gap axes (0 = the system-average level that year); they enter W, not X (levels in X are an open PI question). `log_n_system_tm1` (centred log system size at t − 1) is in W because the IV grows with the number of alternatives (`assemble_cell()` logs cor(IV at a null point, log n) = 0.856 and the R² of that IV on W without / with the column, 0.091 / 0.767 at lag 3).

Timing: positions at t − 1; changes over [t − 1 − h, t − 1]; dyadic controls (contiguity, ATOP alliance, joint democracy polity2 ≥ 6, major-power count) at t − 1; all via calendar joins. Conflict history: monadic and dyadic MID onset counts in the 10 years before the window starts (`log1p()`, centred), peace years as of the window start (full control set only). Democracy and alliance NA → 0 with `*_missing` indicator columns in the data. Panel 1816–2014 (COW-MID 5.0 end); state index = all states appearing as initiators or targets (S = 217).

### 6.7 Coalition weights (PI decision 2026-10-25, D8; the `weights` switch)

Every A × B pair of a multilateral dispute is an event, so an initiator joining a dispute against k Side-B members contributed k softmax terms. The stage-2 part of the log density is a composite likelihood with dispute-level weighting: event e enters as w_e · log P(j_e | n_e) with w_e = 1 / (the number of distinct targets the same initiator has in the same dispute and onset year, counted within the sample), so such an initiator-dispute-year contributes one decision's weight while targets in different disputes stay separate decisions. The build writes `w_mid`, `w_fmid`, `w_fmidcow`, `w_fmidub`, `w_war` on the pre-collapse rows (`event_weights()`), 0 outside the sample; the collapse sums a key's weights and caps at 1 (72 mid keys hit the cap). On the data (mid, before the collapse) 3517 events in 2883 initiator-dispute-years have weights summing to 2883; minimum 1/37 (Yugoslavia → the 37 KFOR states in dispute 4343, 2000); in dispute 4186 (2000) the 38 KFOR states are each a Side A initiator against Yugoslavia alone, so their events keep weight 1. `build_stan_data()` passes the weights as `ev_w` (1 everywhere when absent); both Stan programs multiply the stage-2 log-lik terms by `ev_w`; with every weight 1 the log density equals the unweighted program (asserted against the dd48506 binaries). `log_lik_stage2` stays the unweighted per-event term. The joint state-year LOO and the permutation statistic use the weights; the stage-2-only LOO does not. Events that leave stage 2 (`fighting_at_ref`, target without positions) do not rescale their group's remaining weights (logged).

## 7. Model

Stage 1 (state-year n, state i): `y1 ~ Bernoulli(logit⁻¹(α0 + W'α + c_i + τ·IV))`. Stage 2 (event e by state-year n):

```
U_ij = −½ (x_ij − μ_n)' Σ⁻¹ (x_ij − μ_n) + b_j + z_ij'γ [+ u_i'v_j in M2]
x_ij = (M_i − M_j, R_i − R_j) at t−1 (scaled)      μ_n = μ0 + X_n β + a_i   (μ0 is 1 × D, β is P × D)
P(j | n) = softmax over C_n^{s2}                     IV_n = log Σ_{j∈C_n^{iv}} exp U_ij
log density = Σ_n log Bernoulli(y1_n) + Σ_e w_e log P(j_e | n_e) + priors
(c_i, a_i1, a_i2) ~ MVN(0, Ω)                        b_j ~ N(0, σ_b) (soft sum-to-zero)
```

μ0 is the intercept of the ideal target position (PI-approved 2026-10-23). Because x is the initiator-minus-target gap, μ0 is the population-average preferred gap, the paper's central quantity (bully > 0, peer ≈ 0, challenger < 0). In the main specification X holds only the two power changes, so μ0 is the ideal gap of a state with zero power change; every W control is centred, so α0 is the stage-1 log-odds at zero power change and average controls. `mu_state` at `x_ref` = μ0 + x_ref'β + a_i.

Identification and priors. Only Σ⁻¹ and Σ⁻¹μ are identified by the choice data; μ is recovered as Σ(Σ⁻¹μ) automatically. The "a" models fix corr(a_i1, a_i2) = 0 but keep corr(c, a_1) and corr(c, a_2), parameterised so Ω is always PD (r2 = q·√(1 − r1²)). AME: V and U column-centred (the translation U[,k] + h, b_j − h·V[j,k] is removed), the anchor block of U (K states with most initiations) lower-triangular with positive diagonal, U/V scale fixed by N(0, 1) priors; Procrustes alignment across chains in `03_extract.R` is a safety net whose rotation deviation is logged. Priors: μ0 ~ N(0, 2) per element (`prior_mu0_sd`, data), β ~ N(0, 1) (`prior_beta_sd`), σ_x ~ lognormal(0, 1), sd of random effects ~ half-N(0, 1), LKJ(2) on the 3×3 correlation of (c, a_1, a_2) in the b models and its restriction to corr(a_1, a_2) = 0 (the disk density ∝ 1 − r_1² − r_2², exact Beta(2.5, 2.5) marginals; review fix 2026-10-25, B4) in the a models, (ρ + 1)/2 ~ Beta(2.5, 2.5) on Σ's correlation, (r + 1)/2 ~ Beta(2.5, 2.5) on corr(c, a) in the one-dimensional models, γ, α, τ ~ N(0, 1), α0 ~ N(0, 2.5). The correlations shared across models carry the same prior, so m / r / a / b differ only by their restrictions (prior-only runs in `tests/testthat/test_priors_a_vs_b.R`).

Both programs evaluate the likelihood as `reduce_sum_static` over state-years (`grainsize` is data; one slice = the exact serial sum). Generated quantities: `log_lik_stage1[N1]`, `log_lik_stage2[E]`, `mu_state[S, D]` (at `x_ref`), `mu_event[E, D]`, `Sigma`, `Omega`, `R_Om`. The programs reject an empty IV set, a chosen row outside the stage-2 set, a stage-2 set that is empty or chosen-only, and a row in neither set.

Fit provenance and safety. `02_fit_model.R` reads `manifest.csv` and `git rev-parse HEAD` before assembling, passes the manifest snapshot into `assemble_cell()` (every input file's nrow and sha256 are checked against it), re-reads the manifest before saving and saves the fit as `fit_<cell>_STALE.rds` if any digest changed while sampling (downstream tools need `--allow_stale`). The compiled program must be the checked-out source (`assert_compiled_current()`: model_code identity plus the `.sha256` sidecar written by `compile_models.R`). Fits are written atomically (`.rds.tmp`, read back, renamed). Pilots (`--pilot=TRUE`) are saved with `_pilot` and refused by every production tool (`assert_production_fit()`); `03_extract.R` and `04_compare.R` also refuse a fit whose stored data digests differ from the current manifest (`assert_fit_current()`), and the robustness tests check that the fits score exactly the assembled cell (`assert_fit_matches_cell()`).

## 8. Outputs

All in `v2/output/` (gitignored). Axis labels: `material` = the M axis (log CINC share ratio, scaled), `relational` = the R axis (log PageRank ratio, scaled); one-dimensional cells carry an `axis` column.

| file | written by | columns |
|---|---|---|
| `fit_<cell>.rds` | 02 | list: `fit` (stanfit), `dat` (Stan data), `cell`, `state_codes`, `W_cols` / `X_cols` / `Z_cols`, `sy_keys` (ccode, year, state_idx, y1), `ev` (sy_row, alt_row, w), `x_ref`, `manifest` rows, `data_dir`, `spec` (the switches), `winsor` (cut points, or NULL), `param_labels`, `args` (sampler, provenance, Slurm ids), `stale` |
| `standata_<cell>.rds` | 02 `--dry_run` | `dat`, cell description, manifest rows, commit, `spec` |
| `summary_<cell>.csv` | 02 | `parameter`, `axis`, `covariate` (W / X / Z name; random-effect order c, a_material, a_relational; Σ / Ω index pairs; `IV` for τ, `intercept` for α0), then rstan's `mean`, `se_mean`, `sd`, quantiles, `n_eff`, `Rhat` |
| `typology_<cell>.csv` | 03 | per state: `state` (ccode), posterior type probabilities (`p_bully`, `p_peer`, `p_challenger`, the two hybrids; 1-D: `p_weaker`, `p_peer`, `p_stronger`), `map_type`, `mu1_mean` / `mu1_lo` / `mu1_hi` (material) and `mu2_*` (relational) or `mu_mean` / `mu_lo` / `mu_hi` (1-D), `r_mean`, `theta_mean`, `theta_lo` / `theta_hi` (circular 90% interval, `theta_width` the arc), `n_events_as_initiator`, `r0`, [`axis`]. Bully = μ in the (+, +) quadrant beyond r0 (targets below it on both axes), challenger (−, −), peer within r0 of the origin, hybrids off-diagonal; 1-D: weaker / peer / stronger (μ > r0 means the state targets others positioned below it) |
| `mu_state_<cell>.csv` | 03 | `state`, the μ columns above, [`axis`], `n_events_as_initiator`: posterior mean and 90% interval of μ_i at x_ref |
| `mu_event_<cell>.csv` | 03 | per event: `event`, `sy_row`, `initiator`, `target`, `year`, the classification of μ at the event's own X |
| `typology9_<cell>.csv` | 03 (D = 2) | 9-cell classification with per-axis peer bands `band_M`, `band_R` (default r0): bully, challenger, peer, material / relational opportunist, material / relational challenger, compensator A (M down, R up), compensator B (M up, R down); cell probabilities and the modal cell |
| `reach_<cell>.csv` | 03 | posterior mean and 90% interval of r_i = ‖μ_i‖ and of |μ_i1|, |μ_i2| |
| `focus_sigma_<cell>.csv`, `focus_states_<cell>.csv` | 03 | descriptive: the pooled Σ (SDs, correlation) and each initiating state's empirical dispersion of its chosen targets' positions around the posterior mean μ at the event (`sd_M`, `sd_R`, correlation, `n_events`; `unreliable_lt_min_events` when n < 5) |
| `ame_U_<cell>.csv`, `ame_V_<cell>.csv`, `ame_pairs_<cell>.csv` | 03 (M2) | Procrustes-aligned posterior means of U (with `anchor`) and V; the 200 initiator × target pairs with the largest |u_i'v_j| (`uv`) |
| `compare_loo_<group>.csv` | 04 | per model: `D_mode`, `joint`, `E`, `N1`, `sum_event_weights`, LOO ELPD, SE and Pareto-k > 0.7 counts for stage 2 (`elpd_stage2`), stage 1 CONDITIONAL (the Bernoulli term alone, the state-year's events staying in) and the JOINT state-year LOO (`elpd_joint_stateyear` = log_lik_stage1[n] + Σ w_e log_lik_stage2[e]), and the deltas to the best model |
| `compare_contrasts_<group>.csv` | 04 | planned contrasts with `*_elpd_diff_b_minus_a` and `*_se_diff` for the three LOOs: within each family r vs m, a vs m, a vs r, b vs a; across families M2 vs M1 for m, r, a, b |
| `compare_misclass_<group>.csv` | 04 | a vs b MAP-type disagreement rate and mean absolute probability shift among initiating states, per family |
| `compare_m1m2_<group>.csv`, `compare_m1m2_states_<group>.csv` | 04 | M1x vs M2x: τ, μ0, β, γ, Σ[1,2], corr(a_1, a_2) (`R_Om[2,3]`) posterior means and 95% intervals side by side, `sign_agree`, `m1_excludes_zero` / `m2_excludes_zero`; per x the correlation of the `mu_state` posterior means per axis and the `type_agreement_rate` at r0 |
| `measurement_<cell>_permutation_<mode>.csv` (+ `_draws`) | 05 perm_error | observed a-vs-b advantage, the permuted distribution (dR shuffled within each stage-2 set), `p_value` (refit mode: (1 + exceedances) / (1 + n_perm); rescore mode: NA, a reliance diagnostic), `interpretation` |
| `measurement_<cell>_error_inflation.csv` | 05 perm_error | per delta (log points added to the relational LEVEL of state-years that initiated): refit classification agreement with the baseline, parameter shifts |
| `measurement_<cell>_temporal_holdout_<cutoff>.csv` | 05 (`--tests=holdout`) | per-draw hit rates of the chosen target ranking first / in the top 5 (`hit1_pdraw`, `hit5_pdraw`), the rank of the posterior-mean probability (`rank_pmean`, `hit1_pmean`, `hit5_pmean`), log scores against a uniform baseline, constants frozen on the training years |

The `<group>` of `04_compare.R` is `<sample>_lag<h>_all_pr<suffix>`.

## 9. Archived options

Removed from the main scripts, the sbatch files and the README on 2026-10-25; the helpers and their tests remain, the CLI options do not (passing one to `02_fit_model.R` is an error). To use one, call `assemble_cell()` directly with the argument named here.

| option | former CLI | internal argument | where it lives now |
|---|---|---|---|
| politically relevant choice set | `--variant=polrel` | `VARIANT = "polrel"` (the build still writes `alts_lag<h>_polrel.rds`) | `assemble_cell()`; `archive/slurm/02_fit_array.sbatch` |
| sampled alternatives | `--sampled_m=m` | `SAMPLED_M` (`sample_alternatives()` in `helpers_model.R`, tests in `test_sampled_alternatives.R`) | `assemble_cell()` |
| capital distance in Z | `--distance=TRUE` | `DISTANCE` (`log_capdist` stays in the alternatives files) | `assemble_cell()`, `z_columns()` |
| contiguity off | `--contiguity=FALSE` | `CONTIGUITY` | `z_columns()` |
| eigenvector centrality | `--centrality=eig` | `CENT = "eig"` (`R_eig`, `dReig_h`, `dR_eig` stay in the data) | `assemble_cell()` |
| conflict history off | `--conflict_history=FALSE` | `CONFLICT_HISTORY` | `w_controls()`, `z_columns()` |
| dyadic peace years in Z | `--dy_peace_years=TRUE` | `DY_PEACE_YEARS` (`dy_peace_years_h` stays in the files) | `z_columns()` |
| brms benchmarks (directed-dyad logit, old gap regression) | `benchmarks.R` | `archive/R/benchmarks.R` (`sim/run_sim.R --task=benchmarks` still exercises it) | archive |
| timing pilots as a job array | `04_pilot.sbatch` | `02_fit_model.R --pilot=TRUE` by hand | `archive/slurm/` |
| per-fit extraction loop, a-vs-b LOO CLI | `05_extract.sbatch` | `03_post.sbatch`; `measurement_tests.R --fit_a= --fit_b=` still works | `archive/slurm/` |

`archive/R/run_measurement_real.R` became the `perm_error` preset; `archive/sim/` holds the four former simulation drivers (now `sim/run_sim.R --task=`) and three one-off timing scripts of 2026-10-22; `archive/slurm/03_fit_all_variant.sbatch` was already deprecated (its cell names collide with the main grid).

## 10. Changelog

`v2/CHANGELOG.md` records every audit finding and PI decision with the counts of the day (F1–F12 of the original audit; the Codex reviews of 2026-10-22 and 2026-10-24; the fatal-MID, μ0, power-measure, joiner, escalation, duplicate-collapse and provenance decisions; the full review fixes of 2026-10-25; the restructuring of 2026-10-25). Tests: `Rscript v2/tests/run_tests.R` (needs the compiled models; `V2_REF_BUILD` pointing at a directory with `<model>_3a9e7d0.rds` and `<model>_dd48506.rds` enables the equivalence tests against the earlier programs).
