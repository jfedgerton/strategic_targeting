# Stage-1 code and stage-2 design review

21 November 2026. **The initiators-only transformation looks correct for the stated window-based definition. I found no demonstrated P1 in that transformation.** The specified fresh `TYPES=init OUT=output_init` run can proceed from the perspective of this narrow source review, with the provenance and sensitivity qualifications below. **The stage-2 memo needs specification fixes before implementation**, particularly the missing early defenders, the opportunity offset, and the hazard/one-shot comparison.

## Scope

Read `docs/PROMPT_gpt_review_stage1.md`, all of `docs/stage2_design_memo.md`, the changed initiator block in R06, the runner/config integration, and the relevant existing role, robustness, inference and summary interfaces needed to assess two-group compatibility. This was not a new audit of unchanged data construction or estimators. I did not reopen historical coding decisions.

Archive: `/Users/jarededgerton/Downloads/targeting_multiverse_stage1.zip`; SHA-256 `c7e91a7e209e83e69fb5efcee5ad195e16026b7f0612b33925f2011890a6a90e`; HEAD `3015bbbf1c93975438de82c9559b45cb1d0fd527`. Line references are one-based archive-relative references. Output is at the established `gpt_audit/targeting_multiverse/docs/REVIEW_stage1.md` location.

**No R scripts were run and no source files were modified.** Native Bash syntax checks passed for the runner and config. The supplied Linux LOOP_DONE/restart results, event counts and numerical estimates were not independently reproduced. Only this review was written. The memo is design material, not implemented behavior; its findings below are specification problems, not allegations of existing fitted-code errors.

## Part 1 — initiators-only code

### What checks out

- **Whole strata, correct role mapping.** `R/06_make_variants.R:112–122` requires the four expected roles, rejects inherited controls typed to those four roles, selects events whose chosen row has type 0/1, and retains every alternative of each selected event. It verifies the retained event count and types, then sets roles to bil/coal and the tag to init. The final one-chosen-row check and atomic publication are at 165–168. Event ids need not be renumbered: downstream routines operate on their values/groups rather than requiring consecutive ids.
- **Window semantics preserved.** The parent bcjs classification distinguishes initiating-side from target-side early entrants (`R/05_dump_events.R:324–347`). Initiating-side entrants promoted by the window remain coal in init; target-side entrants stay joinB and are removed. This matches the requested operational definition, though “initiators only” includes window-defined co-originators, not just original COW originators.
- **Ordinary and unit-weight dumps.** `run_all.sh:126–134` recursively requests the corresponding bcjs frame, preserves the sample/window arguments and `_uw` naming, builds the initiator variant to a staging destination, checks roles/types/sample/window and moves it into place. The parent dump creates unit weights before filtering. I found no recursive argument or local-variable error in this branch. The underlying R06 write and final move prevent a partial variant being accepted as the final frame.
- **Influence contrast.** Runner 237–240 supplies `coal-bil:u_R_less` under init. R09 defaults its ranking quantity to the first requested contrast (`R/09_influence.R:18`), so the top-drop list now ranks the intended two-group contrast without a nonexistent defender/chooser label.
- **Generic two-group consumers.** Surface construction, pairwise contrasts, typed controls, resampling and event-study type loops use the roles array. The two-group frame does not require a third type. Aggression and defense controls created *after* filtering use bil/coal; the four-role control guard at R06:116 prevents silently inheriting inappropriate typed columns.
- **Previous cache issue.** The updated runner now uses an invocation-specific cache and explicit invalidation (`run_all.sh:30–42`), rather than trusting persisted size/mtime across restarts. The prior cross-restart E1 scenario is addressed. This is a limited check of the accompanying change, not another full provenance audit.

### S1 — P2: effective TYPES is absent from the fingerprint

**References:** `config.sh:21`; `run_all.sh:86–88`, 96–122, 148–149, 237–240, 249–259.

Making TYPES environment-overridable changes what must be fingerprinted. The fingerprint hashes the *text* of config.sh, containing `${TYPES:-bcd}`, but does not append the resolved TYPES value. With the same files and output directory, bcd and init therefore have the same fingerprint. This differs from the previously hard-coded setting, where editing it changed config bytes.

**Failure case:** complete a bcd run, then launch `TYPES=init` against that same OUT. Init dumps and sweep keys change correctly, but cached influence, resampling and event-study outputs have scheme-independent paths and are reused based on existence. Influence can then retain the old chooser contrast/ranking, and cached bootstrap/event-study files can describe the old three-group population. The fingerprint does not reject the switch.

**Minimal fix:** append the effective TYPES to the fingerprint (and log/validate supported scheme values). The prescribed new `OUT=output_init` avoids this case, so this is not a reason to block that exact fresh-directory command. It is a real new configuration/provenance hole, not an error in the initiator filter.

### S2 — P2: skipping all joiner-side drops assumes away promoted original joiners

**References:** `run_all.sh:199–207`; `R/05_dump_events.R:324–335`; `R/06_make_variants.R:150–160`.

The comment “init frame has no joiner strata” is true for the *final type labels*, but some retained coal events were original joiners promoted by the window. Their retention as originators depended on the initiating-side map. R05 deliberately retains `orig_joiner` and `side_case` tags for those events. R06's targeted hard/reviewed drop selects those tags, not the final joinA/joinB labels.

Thus a reviewed-case early attack-side entrant can be type coal, remain in init, and still be removed by the targeted sensitivity that the runner now skips. At positive windows this check is not necessarily empty. I have not established the actual number of such retained events from the supplied archive.

**Minimal fix:** keep the targeted drop when retained chosen rows have a relevant side_case tag, or at least count and disclose these rows before skipping it. The retained whole-dispute drops are a broader sensitivity and cover these disputed cases, so this is a loss/misdescription of a narrower check, not evidence that the main fit is wrong.

## Existing grid: validity and changed interpretation

| Check | Assessment under init |
|---|---|
| Type permutation | Mechanically supported: the main frame has untyped controls, and R08 rebuilds the bil/coal interaction columns after relabeling (`R/08_resample.R:81–90`). Interpret only as an artificial-label diagnostic. Coalition membership is not randomly assigned; shuffling event labels does not preserve all initiator/dispute dependence and is not a causal test of coalition effects. |
| Event study | Two-group loops/reference work (`R/10_event_study.R:75–83`). Now compares the chosen targets of bilateral versus coalition-originator events, rather than joining behavior. It is descriptive target-position history, with future/post-onset periods exposed to conflict effects. Annual entry/onset alignments can coincide for most originators; show overlap/differences rather than presenting two nearly identical alignments as independent evidence. |
| Typed controls / pooled benchmarks | Meaningful: bil/coal interactions are built from the retained frame; pooled M0 compares those same strata. Keep complete versus complete_base matched for aggression. The init transformer properly refuses an already four-type-controlled parent rather than trying to relabel its controls. |
| Step A | Now target choice among initiator/co-originator MID events, not a monadic model of whether a state initiates. It does not estimate the probability of conflict or coalition formation. |
| Default Step B | Still filters `!orig_joiner` and GML onset year (`R/07_fit_sweep.R:235–239`). Thus it intentionally excludes original joiners promoted into the window-defined stage-1 coal group. “Risk set unchanged” can be correct while Step B is a narrower population than the new main target-choice fit. Label this distinction; it is not an exact severity model for every retained stage-1 co-originator. |
| Old Step B all/pair | Under init, “all” means all *retained initiator-frame* events, not all participants or all original joiners. Rename the description/output caption accordingly; it no longer reproduces the earlier all-participant comparison. |
| Intensity, lag/sample restrictions, scale, measurement error, defense exposure | Remain meaningful for the retained originator population. Fatality restrictions still condition on eventual severity; changing group inclusion does not make them causal recruitment tests. |
| Influence/top drops | Recompute within the init population and rank coal−bil, as implemented for a fresh OUT. They can select different disputes than bcd; keep the same historical sensitivity lists explicitly separate. |
| Four-type appendix | Runner 188 still fits bcjs, including both joiner groups. Computationally valid, but this is an extended-population comparison with different actors, not merely an unpooled parameterization on the same stage-1 observations. |
| Unit weights | Still changes from composite event weighting to one per retained event. Filtering whole strata preserves original weights; do not describe it as estimating the number of unique state–dispute decisions without checking that aggregation. |

Two groups do not inherently invalidate the inherited grid. The main changes are the population and what each comparison can support.

## Part 2 — stage-2 memo

I accept the five PI choices: both estimators, relevance/opportunity separation, fatal MIDs, both tie aggregations and the re-entry sensitivity. The following P2 items specify how to implement those choices without changing the question silently.

### M1 — P2: early target-side joiners fall between the stages

**References:** memo 9, 12, 21–24; stage-1 R05:324–329 and R06:118.

Example: W=1 month; k enters the target side during that window. Stage 1 correctly leaves it joinB and drops it from init. Stage 2 begins only after the window and excludes anyone already participating. k therefore has no entry decision in either stage. The memo's claim that window changes move actors consistently between the stages is false for these defenders. Uncertain dates can create additional boundary cases if the stage-2 cutoff does not use the same definite-bound rule.

**Specify:** account explicitly for early target-side entries, for example in an onset/early-entry component of stage 2, while preserving the approved originator convention. Alternatively, define the stage-2 estimand as *late* entry among nonparticipants surviving to the landmark and acknowledge that the two stages do not partition every participant. Do not claim to explain all defenders with that restricted estimand. State the shared reference date, interval-bound rule, cutoff month and handling of disputes ending before the cutoff. The two-month aggregation also needs a corresponding window/boundary rule, not just rebinned outcomes.

### M2 — P2: hazard and one-shot populations/horizons are not yet comparable

**References:** memo 9, 12, 22–24, 55, 62.

The hazard starts after the window and admits every current system member; the one-shot says one row per k,d at onset with joined/never. It is unclear whether the latter includes early participants or states not yet in the system at onset, and what “never” means for late-starting/censored observations. A state entering the system during a long dispute cannot have a meaningful onset covariate row under the same definition.

**Specify:** a common landmark cohort and follow-up end, or clearly distinguish dynamic entrants as a separate population. Exclude already entered states consistently; define state exits, sample-end censoring and missing/end-date conflicts. “Never” should mean no observed entry by the specified endpoint, with incomplete follow-up handled explicitly. Freeze one-shot side members at the declared baseline rather than using eventual membership.

For comparison, aggregate monthly hazard probabilities into side-specific cumulative incidence over that same follow-up (including survival through preceding months). Do not compare a monthly join probability or logit coefficient directly with an eventual one-shot probability/coefficient. Changing covariates/current side composition versus onset covariates is also a substantive difference, not only sensitivity to date coding. Preserve both estimators; define the comparison on standardized predicted probabilities.

### M3 — P2: the opportunity offset needs an alternative-specific definition and a time scale

**References:** memo 14, 27–29, 62.

There is no specified utility for the offset. Adding the same λ ln q to A, B **and stay out** cancels completely. Adding it to both join alternatives and zero to stay out shifts participation versus nonparticipation, but cancels from A-versus-B log odds conditional on joining. This cannot by itself adjust the *relative* plausibility of joining different sides. Moreover, a propensity for eventual entry in a dispute is not a monthly entry hazard.

**Specify:** the exact offset vector, what q predicts, its denominator/risk set, training exclusions, time unit and exposure duration. For example, if the intention is a common entry-opportunity sensitivity, state offsets `(λ ln q, λ ln q, 0)` and report that interpretation. If intended as a calibrated probability baseline, define the corresponding odds/logit construction instead of borrowing target-choice weights unchanged. Align monthly and one-shot q definitions with M2.

“Fitted out of sample by era” also needs actual folds: exclude all rows of a held-out dispute, and define whether future eras can train predictions for earlier eras. Participant-distance features must use the eligible pre-decision members, not the eventual entrant list. Training outcomes must represent third-party entry, excluding stage-1 originators. Otherwise outcome-derived future membership can leak into q even though the prediction row is held out. These are specifications of the chosen offset approach, not an objection to reporting it.

### M4 — P2: sides, target and relevance need operational definitions

**References:** memo 13–14, 21–24, 34–43.

“Current members” needs to distinguish currently active participants from everyone who has ever entered, especially after exits and in the re-entry sensitivity. State whether switches between sides are ignored, censored or new decisions. A lagged side roster prevents including a simultaneous entrant in its own predictors, but does not decide these rules.

An original dispute can have multiple originators and original targets; the stage-1 event is directed initiator–target–year and can carry multiple source disputes. A stage-2 k,d row does not automatically have one scalar “original target.” Define the stage-2 dispute key, original target set and its centrality aggregation, and ensure the initiating-side map corresponds to that same dispute. Do not arbitrarily inherit the primary target of whichever stage-1 event is encountered first.

For tiers, fix the region classification, membership date and multiple-region rule; which ATOP pact types/valid dates count; and whether relevance is frozen at the landmark or updated with the lagged roster. If a state moves from unlikely to likely as sides grow, the likely-only estimate is conditional on time-varying eligibility, whereas the one-shot subset needs a baseline classification. Do not define it by whether the state *ever becomes* relevant. Geography and roster information used in the opportunity model must obey the same timing rules.

For power and ties, specify edge direction, normalization, layer combination/weighting, zero-link handling and capability aggregation. Fix whether the original target's centrality is baseline or annually updated. “All lagged one year” alone does not settle which side roster to aggregate over.

### M5 — P2: model and identification claims exceed the stated controls

**References:** memo 27–29, 35–40, 46–49, 54–58.

An alternative-specific conditional logit over three alternatives is a multinomial-logit choice formulation; calling multinomial logit an alternative estimator does not create an independent check unless a different restriction is stated. Include join-alternative intercepts and specify flexible duration dependence for the monthly baseline. A k,d,m-level quantity such as own centrality or duration is identified relative to stay out only when interacted with the join alternatives; specify whether the A/B coefficients are shared or separate. Side-size/tie sums and side power can be strongly related, so define coefficients and prediction contrasts rather than using the concept labels as the test.

Lagging the roster prevents contemporaneous self-inclusion, not endogenous sorting into sides. Lagged annual ties can already reflect an ongoing long war. Onset-frozen ties or an explicit pre-dispute reference provide a useful specified comparator within the chosen design. Fatal-MID sampling conditions on eventual severity, and both stages are conditional on selected conflict events. Those facts limit interpretation; the approved sample is still usable for the stated conditional descriptive question.

In particular, remove or soften line 49: a null does **not** establish that stage 1's pattern does not arise from recruitability. Measurement, limited variation/power, early entrants omitted by the risk set, equilibrium sorting and the selected sample can all produce a null. Conversely, a positive tie association does not identify recruitment as the cause of the target-choice pattern. Report H1 as conditional probability contrasts with uncertainty. Ties to B enter A's probability through the competing alternative under the proposed utility specification; that is not an independently estimated direct deterrence coefficient in A's utility unless one is explicitly added.

## Decision

**Stage 1:** proceed with the explicitly separate fresh `OUT=output_init`, retain the supplied TEST results as external validation evidence, and fix S1 before treating scheme changes within one output directory as safe. Count affected promoted-originator side-case rows and describe the narrower Step B/four-type comparisons accurately. No new numerical failure in the initiator filter was demonstrated here.

**Stage 2:** revise M1–M5 before building data. This preserves all five PI decisions while making the unit, stage boundary, opportunity treatment and estimator comparison well-defined. The highest-priority concrete correction is the early-defender gap; next are a common follow-up population and an explicit offset vector/training protocol. These are smaller specification decisions to settle now, rather than bugs to discover after fitting.
