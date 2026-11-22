# Coalitions and the targeting of isolated states — multiverse loop

One question, one main model, every robustness check in one run.

**Question.** Do states that attack as part of a coalition choose relationally isolated targets, while bilateral initiators and states joining to defend do not? Is that a selection effect, or do disputes with isolated states just fail more often (failed bargains), or are isolated states pariahs that get punished?

**Main model.** Conditional logit of target choice in fatal militarized disputes, 1816–2014. One stratum per initiation event; the alternatives are the system members that year with a relational position, covariates and stage-2 availability (see "Choice set" below; not every system member). Utility = a type-specific surface in the material gap (log CINC share) and the relational gap (network position at t−1, scaled within year): dM, dM², dR, dR², dM×dR, plus dyadic controls. Gaps are initiator minus target, so a positive dR means the target is less central.

| Setting | Main choice | Appendix |
|---|---|---|
| Initiator types | 3: bilateral / initiating-side coalition (coalition originators + attack-side joiners) / target-side joiners. Pooling originators with attack-side joiners is an assumption (their surfaces are not distinguishable at windows ≥ 1 month, which is not proof of equivalence) | 4 types (originators separate), at every window |
| Who is an originator | Present within **1 month** of the dispute's first originator (definite date bound) | 0–6 months |
| Coalition status | At onset (`--bilateral=onset`) | — |
| Joiner side | Relative to the dispute's initiating side (case-reviewed; see `docs/decisions.md`) | Drop the 2 hard / 13 reviewed disputes |
| Sample | Fatal MIDs | All MID onsets (failed-bargain step A); intensity samples through the full multiverse: jointly fatal (both states have fatalities in the dispute; state-level, not pair-specific) and dispute fatality level ≥ 2, 3, 4, 5, 6 (GML codes: > 25, > 100, > 250, > 500, > 999 deaths); unknown fatality dropped and counted |
| Relational measure | `W2contig_ptaall|igo+def+off+neu+nonagg+consul+dca+un|pr` | All 15,366 measures (incl. the book-aligned sry weightings with clipped and min-max edges, and unweighted PTA, added 2026-11-22; `docs/measure_alignment_2026-11-21.md`) |
| Choice set | System members with a position, covariates and stage-2 availability in the event year | — |
| Inference | Two-way cluster-robust (initiator × dispute) | Pigeonhole two-way bootstrap; permutation diagnostics |

## Stage 1 (`TYPES=init`, PI 2026-11-21)

`TYPES=init OUT=output_init bash run_all.sh` runs the same grid on initiation decisions only: bilateral initiators vs coalition originators (originators within the window, including initiating-side joiners promoted by the window). Joiner strata of both sides are dropped; joiners are stage-2 actors. `TYPES` accepts only `bcd` or `init` and is part of the fingerprint, so one output folder holds one scheme. What the inherited checks mean under init (GPT stage-1 review):

| Check | Population under init |
|---|---|
| Main fit, intensity, lags, scale, measurement error, defensibility, aggression | the retained bilateral and coalition-originator strata |
| `dropside_hard`, `dropside_reviewed` | promoted coalition originators whose side came from the case-reviewed map (count logged by R/06; 0 is reported, not skipped) |
| `stepB_orig_disp` | narrower than the main fit: excludes original joiners promoted into the coalition group, and decisions outside the dispute's GML start year |
| `stepB_all_pair` | all retained init-frame events (not all participants); comparison only |
| `fourtypes_w*` | the extended four-type population (bilateral, originators, both joiner sides); a different set of actors, not an unpooled version of the init fit |
| Type permutation | artificial-label diagnostic on bil/coal labels, not a randomization test |
| Event study | chosen targets of bilateral vs coalition-originator events; entry and onset alignments largely coincide for originators |
| Influence / top drops | recomputed within init, ranked on `coal-bil:u_R_less`; may pick different disputes than bcd |

## How to run

```bash
# 1. edit the paths at the top of config.sh (data_entry, network files, positions to reuse)
# 2. run everything (resumable: rerun the same command after a stop and finished steps are skipped)
bash run_all.sh
# mechanics check only (2 windows, main measure, tiny B), written to output_test/ so it never mixes with production:
TEST=1 bash run_all.sh
```

Everything is written to `output/`; the log is `output/logs/run_<time>.log`.
- **Provenance.** At start the run fingerprints the code (contents of `run_all.sh` and every R script; doc edits do not count), `config.sh`, the event-data manifest, the sha256 of every positions file, the network inputs and the layers file, and the R and package versions. An output folder belongs to one fingerprint; if anything changed, the run stops (use `OUT=<new folder>`). A non-empty folder the runner did not create is refused. Production refuses uncommitted code. Artifacts built from supplied lists (the Phase A ids, the exclusion lists, variants) carry a stamp of those inputs and are rebuilt, with their sweeps, when the inputs change; each measure set is checked against its id list.
- **Restarts.** Rerun the same command after a stop: finished steps are skipped. A step counts as finished only if its R process exited 0 and its output was written atomically (sweeps: no ERROR rows, recorded in a per-task status file; the manifest is rebuilt from those records).
- **Ending.** `LOOP_DONE` means every step finished. `LOOP_FINISHED_WITH_GAPS` (exit code 2) lists what is incomplete in `output/status/incomplete.txt`. Until the history-chosen exclusion lists are approved, that file lists them.

## What each step answers

| Step | Script | Output | Question |
|---|---|---|---|
| 1 positions | `R/01_build_layers.R`, `R/02_build_positions.R` | `output/positions/` | The 15,366 relational measures: layer combinations × 8 weightings (W1–W4; book-aligned W6srybook, W7sryyear and their min-max-edge versions W6srybookmm, W7sryyearmm in the core arm) × PageRank/degree × PTA coding (all / in force / any), plus the diplomacy and IGO arms. Diplomacy uses no future waves. The reused 2026-11-18/19 files are linked; the 2026-11-22 arms are always built in the run folder. |
| 2 validation | `R/03_validate_measures.R` | `output/validation/` | Does a measure put the powerful states where they belong? Plausibility rules, correlation with CINC, major-power rank, year-to-year stability, top-10 states by era. Outcome-blind. |
| 3 measure sets | `R/04_subset_positions.R` | `output/sets/` | Phase A (main + 200 random plausible), main only, no-ATOP measures, 25 random plausible measures |
| 4 event frames | `R/05_dump_events.R` | `output/dumps/` | One frame per type scheme × window × sample |
| 5 multiverse | `R/06_make_variants.R`, `R/07_fit_sweep.R` | `output/sweep/multiverse_*` | Every measure × every window, for fatal MIDs, all MID onsets, jointly fatal and fatality levels ≥ 2–6; failed-bargain step B on every measure; four types on Phase A |
| 6 robustness | `R/06_make_variants.R`, `R/07_fit_sweep.R` | `output/sweep/robustness/` | Every check below, on the 201 Phase A measures |
| 7 influence | `R/09_influence.R` | `output/influence/` | Leave-one-dispute-out and leave-one-initiator-out; then drop the top 5 / top 10 disputes |
| 8 resampling | `R/08_resample.R` | `output/resample/` | Pigeonhole bootstrap (main measure at every window, all-MID, 25 random measures). Permutations of dR within events and of type labels are **artificial-assignment diagnostics** (`diag_*`, column `p_diag`), not randomization tests. |
| 9 event study | `R/10_event_study.R` | `output/eventstudy/` | Relative position (within-year z, not a count of ties) of chosen targets vs the alternatives, aligned to the pair's entry year and to the dispute's start year (events with a missing or conflicting GML start are excluded from the onset version; `es_main_onset_fallback` keeps them at the event year as a flagged sensitivity); balanced, unbalanced and pre-only versions; excluded cases and events / weight by type and k written out. Descriptive. |
| 10 summary | `R/11_summarize.R` | `output/summary/` | Appendix tables from manifest-listed outputs only: main-measure estimate first; median and IQR; plausible / failing / all; by measure family and family-balanced; grid percentages labelled as descriptive; requested / valid / invalid / error counts per model |

### Robustness checks (step 6)

| Group | Run name(s) | What it rules in or out |
|---|---|---|
| Rival: target aggression | `aggr` (complete, GML-reconciled 5-year windows only), `aggr_base` (same strata, no control: the comparator), `aggr_censored` (all strata; counts after 2010 truncated, a lower-bound sensitivity) | Isolation proxies the target's own recent aggression |
| Rival: defensibility | `expdef` (dispute onset), `expdefentry` (year before entry) | Isolation proxies having no defense partners |
| Rival: ATOP inside the measure | `noatop_base`, `noatop_expdef` | The defense control works only because ATOP is in the measure |
| Timing | `lag5/10_*`, `samp5/10_*` (year scale and pooled scale), `decomp10_linear/full`, `decomp5_full` | Long-standing isolation vs recent tie loss |
| Reverse timing | `lead_t`, `lead_t1`, `lead_t1_sample_tm1` | Ties lost because of the conflict |
| Failed bargains | `stepA` (target choice among all MID onsets); `stepB_orig_disp` (originator decisions at dispute onset: does the dispute become fatal, GML dispute fatality); `stepB_all_pair` (old membership version, comparison only) | Isolation shapes whom states confront vs whether the dispute fails. Step B keeps only decisions in the primary dispute's GML start year, reports per-type gap terms, type-intercept contrasts and a joint test of both, counts of disputes / initiator-dispute decisions / events, and a cross-tab of its outcome (GML dispute fatality) against fatal-sample membership with unknown outcomes by type and era (`stepB_crosstab_*.csv`). The outcome is eventual severity of the dispute from a pre-dispute reference, not a measured nonfatal-to-fatal transition. Joiners are excluded from step B because first-fatality dates are not in the data, so a joiner cannot be shown to have entered before the dispute turned deadly. |
| Coding | `dropside_hard`, `dropside_reviewed` (joiner strata incl. window-promoted joiners), `dropside_hard_disp`, `dropside_reviewed_disp` (every event of those disputes) | Side coding |
| Concentration | `drop_top5`, `drop_top10` (chosen by the main measure's own leave-one-out influence: adversarial diagnostics, not confirmation); `drop_hist_un`, `drop_hist_strict`, `drop_hist_strict_noborder`, `drop_hist_expanded`, `drop_hist_noruleA`, `drop_hist_anycoder` (chosen on history under a fixed coding rule, `docs/historical_exclusions.md`) | A few great-power or punitive coalitions |
| Weighting and scale | `unit_weights`, `yscale_tm1`, `pooledscale` | Event weights; scale year |
| Measurement error | `merror_def_<sign>_<k>`, `merror_contig_<sign>_<k>` (both signs), `merror_random_<k>_d1..5` (5 independent draws); each writes the rank correlation between the measure and its perturbed version | Fragility to error that mimics alliance or contiguity confounding, and to noise. A fragility test, not a correction for confounding. |

## Script map (for anyone coming from the main repo)

| Here | Main repo (`strategic_targeting`, branch `clogit-roles`) |
|---|---|
| `R/01_build_layers.R` | `v2/R/20_build_layers.R` |
| `R/02_build_positions.R` | `v2/R/21_rel_measures.R` |
| `R/03_validate_measures.R` | `v2/R/23_rel_validate.R` |
| `R/04_subset_positions.R` | `v2/R/26_rel_subset_positions.R` |
| `R/05_dump_events.R` | `v2/R/10_clogit_roles.R` (dump mode only is used here) |
| `R/06_make_variants.R` | `v2/R/24_rel_variants.R` |
| `R/07_fit_sweep.R` | `v2/R/22_rel_sweep.R` |
| `R/08_resample.R` | `v2/R/25_rel_boot.R` |
| `R/09_influence.R` | `v2/R/29_rel_influence.R` |
| `R/10_event_study.R` | `v2/R/27_rel_eventstudy.R` |
| `R/11_summarize.R` | `v2/R/30_summarize_loop.R` |
| `R/lib/*` | `v2/R/{spec,helpers_data,helpers_model,assemble_cell,opp_weights}.R` |

Comments inside the scripts still use the main-repo names. The event data (`data_entry/`) is built by `v2/R/01_build_data.R --entry=later` in the main repo and linked here.

## Conventions

Seed 123 everywhere. 95% intervals. A fit with separation or non-convergence warnings is recorded as INVALID, per model, and has no estimates. Bootstrap intervals and permutation p-values are reported only when B ≥ 200 and ≥ 90% of draws are usable.
