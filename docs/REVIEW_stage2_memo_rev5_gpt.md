# Stage-2 memo revision 5 recheck

22 November 2026. **R2 and R3 are fixed for this recheck. R1 is mostly fixed, but one concrete window-zero contradiction remains, together with a prospective-mask wording requirement. Hold freezing the full windows-0–6 risk set and legality table until those are corrected.** Infrastructure and fixture work remain GO; no further redesign of the cohort, sequential estimator or multiverse is requested.

## Scope

Followed Part B of `docs/PROMPT_gpt_recheck_minmax_stage2_rev5.md`: read memo §13 R1–R3 and the §10 supersession note, with the relevant existing §13 definitions as context. I checked the small stage-1 classification branch needed to verify the newly specified availability mask. No unchanged estimators or historical decisions were re-reviewed. Nothing in stage 2 is implemented.

Package: `/Users/jarededgerton/Downloads/targeting_multiverse_recheck_mm1_rev5.zip`; HEAD `a4f0a72c353739f14d6b7f5fc40aae96ba9e43ff`; SHA-256 `140b76eaac5a41b171577b1bddf505a46b1d9eb1ada6eeda44cca9b4865d7ec8`. Memo line references are one-based lines in that package. No R code was run. Only the two requested reports were written.

## R1 — partly fixed

**Baseline definition: fixed, memo 169–170.** Original status now comes from `participants_timing.orig == 1`, not equality of entry month and onset month. A coded non-originator B entrant in the onset month remains an early defender outcome. The acceptance tests now separate original roles from initiating-side promotions and report participants without fitted stage-1 events rather than pretending they match the fitted-event universe. These address the previous baseline and projection findings. Treat these revision-5 rules as overriding the older baseline wording at 122.

**Availability: correct objective, but the W=0 rule conflicts with stage 1.** Memo 171 says the originator alternative exists in month zero when W=0. In the supplied stage-1 code, promotion runs only inside `if (COORIG_W > 0L)` (`R/05_dump_events.R:218–234`); the target-side reversal and promoted-originator bookkeeping likewise check a positive window at 326–334. Thus W=0 is the original-role specification with **no window promotions**, not a positive-window rule evaluated at an elapsed lag of zero.

Concrete fixture: an `orig=0` state enters A in the onset month. At W=0 stage 1 leaves it a joiner, so it is a third-party A event in the prospective stage-2 partition. A literal mask allowing an originator event in month zero would permit positive probability for an outcome that the agreed stage-1 classification cannot produce. That contradicts the exact projection test at 170. A no-observed-originator category may also be nonestimable if left as an unconstrained utility.

**Minimal correction:** explicitly set `avail_orig=0` throughout W=0, preserving stage 1. Use the shared promotion eligibility rule for W>0 and test the same candidate at W=0 and W=1. If the PI instead wants zero-lag non-originators promoted at W=0, that changes stage 1 and requires a separately authorized change/review; it is not what this package implements.

**Make the mask prospective.** “Inside the promotion window of an initiating-side entry” and “no possible promotion” at 171 must refer to what could occur if an eligible candidate chose A now, not the candidate's eventual observed side or whether this dispute eventually had a promotion. Otherwise an eventual B entrant/nonentrant can have originator availability removed using the outcome, while a future A entrant gets it added.

State that the mask uses only declared stage-1 eligibility/date conventions and pre-decision information, applying the hypothetical A-entry rule to every eligible candidate. Never derive availability from observed future side/entry. The mask may legitimately be zero after W or under W=0; zero observed promotions in a dispute is not by itself proof the option was unavailable.

## R2 — fixed

**Memo 174–176.** One-shot q is now the third-party probability from a three-category complete-horizon model; originator entries are observed competing outcomes, not discarded/censored. Monthly q is explicitly the third-party cause-specific probability with the same competing category and availability mask. The four-category flat-model offset is explicitly `(λ ln q, λ ln q, 0, 0)`. These settle the previously unresolved denominator and offset placement.

The sequential alternative now models third-party/originator/out first, then A/B among third-party entrants. Its incidence formula includes both h_T and h_O in survival. This directly fixes the missing competing survival mass; retrospective and net-risk variants are separately labelled instead of pooled with prospective incidence. It remains a sequential categorical model, not an inclusive-value nested logit. The small first equation being multinomial rather than binary is the necessary consequence of retaining the approved competing outcome.

No additional R2 definition is needed for this recheck, beyond propagating R1's corrected mask consistently to auxiliary and main models. Revision 5 overrides the older censor-originator instructions at 126/163.

## R3 — fixed as a written contract

**Memo 80, 179–182.** The supersession note explicitly makes §13 authoritative over §10's older reference cell/coverage text. The new legality rules include cohort, empty-side policy, one-shot horizon and monthly two-month aggregation; they provide a duplicate-key policy after ineffective levels are canonicalized. This resolves the previous ambiguity about which task list is to be generated. The task list and any computational reduction still need to be materialized before inspecting estimates, as already promised.

The recurrent family now has a stated return-to-risk rule and count-scale contrasts at 180, rather than borrowing first-event incidence contrasts. Encode its month-after-spell-end convention explicitly in fixture dates. Contrast SDs, percentiles and noncomputability rules are specified at 181; the zero-SD/<20-distinct-values rule is a declared policy, not a reason to silently drop fitted cells. Report those contrasts as not computable with denominators, as the text requires.

The side-level fatality bounds at 182 are correct: known fatal participants on both sides establish qualification despite additional unknown participants; a fully known nonfatal side rules it out; unresolved cases are counted and excluded in the reference, with an upper-bound sensitivity. This is now distinguished from stage 1's pair-event definition.

## Freeze decision and targeted acceptance checks

**Do not freeze the entire analysis-ready windows-0–6 product yet.** Correct the single W=0 sentence and specify the prospective interpretation of the mask; then R1 can be closed without revisiting R2/R3. The reference W=1 infrastructure can be developed under that stated prospective interpretation, with provisional masks clearly marked rather than silently choosing a rule.

Add exact fixture assertions for: a same-month non-original B entrant retained as B; the same eligible A entrant classified third-party at W=0 and promoted when the shared W>0 rule permits; originator probability exactly zero after its window; candidates with different eventual outcomes receiving the same mask when their pre-decision eligibility is identical; and a dispute with no observed promotions still offering the option when structurally possible. Verify the auxiliary q and sequential first equation use those same masks.

The remaining original fixture set—initial/empty rosters, system exit, truncation, WWI fold, mixed side-level fatalities, recurrent return after a spell—continues to apply. No new issue is being raised about the approved estimator choices or deferred AME work. Once the mask correction is committed to the contract, the substantive R1–R3 freeze conditions from this recheck are satisfied, subject to the planned implementation/tests matching that contract.
