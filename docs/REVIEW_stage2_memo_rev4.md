# Stage-2 memo revision 4 — builder-contract review

22 November 2026. **Most of revision 3's problems are resolved. Hold the final cohort/outcome and sweep specification until the remaining P2 decisions below are explicit.** Low-level roster/tie/fold infrastructure and synthetic fixtures can start, but freezing an analysis-ready risk set now would still let implementation choose who is a baseline participant and how originator competing events affect some models.

## Scope

Followed Part B of `docs/PROMPT_gpt_review_minmax_and_stage2_rev4.md`: read memo §13 and the corrected `docs/ame_joint_sketch.md`, and checked the earlier sections for conflicts relevant to the contract. No stage-2 code is implemented. Approved design choices and historical coding were not re-reviewed.

Package: `/Users/jarededgerton/Downloads/targeting_multiverse_review_minmax_stage2.zip`; HEAD `e8736b7af9958493558da2243b4ed0d08a571ca3`; SHA-256 `dbd3e2b75d4494f2ffa31fe536e635fb923c29f2bd06278d3245edd0d7f1e104`. One-based references below refer to memo lines unless otherwise named. No R code was run. Only the two requested reports were written.

## Dispositions

| Finding | Status | Evidence and remaining issue |
|---|---|---|
| B1: originator competing outcome | **Partly fixed** | Lines 121–124 now keep future promotions in a prospective cohort and define four outcomes in both estimators. The retrospective alternative is correctly labelled. Baseline membership still equates same-month spell starts with original participants, and monthly originator availability outside W needs an explicit mask: R1 below. |
| B2: censoring/common horizon | **Fixed for the reference comparison** | Lines 128–131 define the observed horizon, restrict the reference one-shot to complete horizons, compare monthly incidence on those same rows/horizons, name conditioning on realized rosters, and rebuild the late-entry landmark cohort. Recurrent outcomes are separated. The truncated-row log-duration sensitivity should be called a descriptive observed-horizon specification, not a censoring correction. |
| B3: initial/empty rosters | **Fixed in its main rules** | Lines 134–137 specify baseline rosters, lagged active rosters thereafter, flagged carry-forward for empty sides, a termination alternative and nonempty target sets. These resolve the original empty-side arithmetic, subject to R1's definition of the initial participants and explicit missing-covariate handling in fixtures. |
| B4: folds and opportunity prediction | **Partly fixed** | Lines 140–141 supply exhaustive dispute-onset folds, keep every row of a spanning dispute together, disclose retrospective fitting, bound q and stop failed fits. The missing WWI fold is fixed. Line 125's originator-as-censoring instruction still does not define the promised third-party probability for one-shot training; R2 below. |
| B5: jointly fatal selection | **Fixed as a separate stage-2 definition** | Line 144 explicitly selects disputes with fatalities on both sides and gives dispute-level severity rules. This is deliberately side-level and can differ from stage 1's pair-event sensitivity. Preserve that label and the unknown-code policy noted below. |
| C1: legal cells/scale/contrasts | **Partly fixed** | Lines 149–159 prohibit incompatible one-shot crossings, deduplicate ineffective cells, align lags/windows, separate population/offset changes, define a reference cell and pre-cost a reduced task list. Remaining task-list and recurrent-contrast details are R3. |
| C2: sequential alternative | **Partly fixed** | Line 162 correctly chooses two sequential logits and explicitly rejects an inclusive-value nested model. But censoring originator events does not define the same competing-risk incidence as the four-category reference; R2. |
| AME language | **Requested corrections fixed; AME specification intentionally pending** | Sketch 24, 27, 34, 37–38 correct the selection claim, qualify exclusions/stability, place LKJ on a correlation component and give the categorical sampling correction. Memo 165 commits the gate before results and rules out favorable-sign/significance gates. The actual gate document is absent, so this is not approval to fit the AME; it need not block the separate-stage data infrastructure. |

## R1 — P2: distinguish baseline roles from merely starting in month zero

**Lines 121–124, 134.** “States whose spell starts in the onset month” and “the dispute's first originator(s) on each side” are not necessarily the same set. A coded non-originator can join during the same calendar month as the original participants. If it joins B, the first reading excludes it as a baseline participant and restores the early-defender gap. If it joins A and is window-promoted, that reading excludes it instead of recording the promised originator competing event.

Write one rule using the actual stage-1 original-role/date fields. In particular, preserve a non-originator B entrant in the onset month if early defenders are intended stage-2 decisions. Do not infer original status solely from month equality. For genuinely unresolved within-month order, choose an explicit coding policy and count the affected cases.

Calling the same stage-1 classification function is the correct way to avoid algorithm drift. But the acceptance test at 124 needs to specify its projection: coalition-originator events are not the entire set of originators, because original bilateral initiators are originators too, and stage-1 fitted events need not cover every participant. Compare baseline roles and promoted initiating-side roles separately, under the same dispute/source and date rules. Define which universe the equality test covers and separately report participants without fitted stage-1 events.

Finally, the monthly originator outcome is structurally unavailable after the promotion window. Require its probability to be exactly zero there, not merely estimate an unconstrained fourth utility with no late originator events. Specify the W=0/no-originator-event case. Otherwise the model can allocate survival mass to an impossible event, or encounter separation while trying to estimate it. A per-row alternative-availability field in the builder makes this testable.

## R2 — P2: opportunity training and sequential censoring still need an estimand

**Lines 125, 140–141, 155, 162; earlier 39–40.** The four-category reference resolves the main competing-event problem. The auxiliary and sequential variants reintroduce it by calling originator events censored.

For a complete-horizon one-shot third-party probability q = P(A or B within H), an observed originator entry is an observed competing event, not missing follow-up of unknown eventual third-party status. Dropping/censoring that row from an ordinary binary logit changes the denominator or requires a specified censoring likelihood. State the auxiliary likelihood and whether q is cumulative third-party incidence, a cause-specific monthly hazard, or a probability conditional on not becoming an originator. Those are different quantities.

Likewise, in a prospective monthly model the actual incidence is

`F_A(H) = sum_m S(m−1) h_T(m) P(A | third-party entry at m)`,

where `S(m−1) = product_{r<m}[1 − h_T(r) − h_O(r)]`, with h_O the originator competing hazard. The sequential recipe at 162 specifies h_T and the conditional side probability but not h_O. Calling originators censored and saying “aggregated as above” does not supply that missing survival mass. For a one-shot binary participation regression, time censoring is not even defined without another likelihood.

Minimal choices: retain an originator competing component in the sequential construction; or explicitly restrict it to the retrospective no-originator cohort; or define it as a net-risk sensitivity under elimination of the competing event and report it in a separate estimand family. Do not describe the latter two as the same observed third-party incidence as the prospective reference.

Also extend the flat-model offset vector explicitly to the fourth category. The old A/B/out vector at 40 is insufficient. If q is exclusively third-party opportunity, the natural proposed vector is `(λ ln q, λ ln q, 0, 0)` for A/B/originator/out; approve that or another rule in writing. The sequential participation-only placement is otherwise correct. Fold completeness and numerical q bounds are now satisfactory.

## R3 — P2 specification completion: one authoritative task table and outcome scale

**Lines 146–159 versus 82–101.** The old §10 still calls the reference a three-alternative/nested model, starts with a likely-only risk set and lists the older scale plan; §13 instead chooses prospective four categories, full-set interactions and a revised coverage plan. State that §13 supersedes conflicting earlier text and update the task-table labels. This is more than cosmetic if the builder generates jobs from §10.

Include every newly added alternative in the legality table: prospective/retrospective cohort, carry-forward/end-at-empty, complete/truncated one-shot horizons, and any two-month aggregation. Specify estimator restrictions and duplicate-cell keys. Fix the reduced task design before results, as promised; a task list with unresolved levels is not yet the contract.

Line 159 defines all H1–H3 contrasts as cumulative-incidence differences, while line 130 says recurrent cells report expected entry counts. Give recurrent cells count contrasts, and keep them out of incidence summaries. For repeated-entry eligibility, record the return-to-risk rule (gap length and new side availability); accepting “new decisions” without that rule does not specify the episode rows. Define where the within-year tie SD and reference-cohort percentiles are computed, and their zero-SD/missingness policy. These quantities can be fixed before fitting without seeing results.

B5's side-level definition is now usable, but encode unknowns with logical bounds: if at least one known fatal participant exists on each side, unknown codes on additional participants do not negate that evidence. If one side lacks a known positive and contains unknowns, qualification is unresolved; a fully observed all-zero side rules it out. If the intended policy is instead complete-case exclusion of any dispute with any unknown participant, say so—it is a different sample. Keep the stage-1 pair-level label distinct from this side-level selection.

## AME: what is resolved and what remains deferred

The corrections remove the earlier claim that a shared target effect automatically corrects conflict/fatality selection, and no longer treat unchanged joint estimates as proof against confounding. The outcome-dependent categorical sampling offsets are separated from opportunity q, with original exposures/histories retained. These are the right changes.

Sketch 17–20 still displays three alternatives and a single target j*, while 37 now references four categories and the memo's target set. Treat those equations as unfinished: the eventual gate must specify target-set aggregation of b, originator utility/availability, scale constraints/priors, and sampling correction for its actual likelihood. No acceptance-gate file is included in this package. The text appropriately postpones thresholds, recovery/SBC tolerances, PPCs and remedial-run limits until before AME results. This is a valid future commitment, not an already completed preregistration, and is not a reason to hold unrelated roster infrastructure.

## Build decision and fixtures

**GO for reusable infrastructure and fixture scaffolding; NO-GO for freezing the full analysis-ready risk set and legality table until R1–R3 are resolved.** Most requested corrections are implemented in prose. The remaining changes are local and should be settled without another redesign of the five approved choices.

Required fixtures should assert exact rows/outcomes, not just successful execution: early A promotion inside W; a coded non-originator B entry in month zero; later B entry; baseline rosters without pre-month participants; a later empty side under both policies; system exit; administrative truncation and the matched complete-horizon comparison; a 1914 dispute spanning folds; mixed jointly-fatal pairs with side-level qualification; unknown participant fatalities; and a recurrent entrant after the declared gap. Add originator-alternative masks after W and a positive originator-hazard example showing that the sequential incidence calculation includes its competing survival contribution.
