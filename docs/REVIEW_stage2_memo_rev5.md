# Stage-2 memo revision 5 — recheck of R1–R3

*Reviewer: independent Claude review agent with no prior context (2026-11-22), run at the PI's request alongside GPT.*

22 November 2026. **R1–R3 are substantially addressed in the new subsections, but none is fully closed, and the stage-1 code the memo requires the builder to call contradicts part of R1.** Four new P2s remain. All are written-contract fixes of a sentence or two, with no redesign. **NO-GO for freezing the analysis-ready risk set and legality table until N1–N4 are settled.** Infrastructure stays GO.

## Scope and what was executed

Followed Part B of `docs/PROMPT_gpt_recheck_minmax_stage2_rev5.md`. Read `docs/stage2_design_memo.md` in full (182 lines at HEAD `a4f0a72`), with attention to §13 R1–R3 (lines 168–182), the §10 supersession note (line 80), and the unchanged §13 B1–C2 text (lines 121–163). Because R1 makes the builder call "the stage-1 classification code itself", I also read the window/promotion logic in `R/05_dump_events.R` (lines 136–140, 162–245, 318–335) and how `run_all.sh:151` passes the window. I checked the revision diff (`git diff e8736b7 HEAD -- docs/stage2_design_memo.md`): lines 1, 80 and 168–182 are new, and nothing in B1–C2 was edited. Nothing in stage 2 is implemented and no R code was run for this part. Memo line numbers are one-based at HEAD.

## Dispositions

| Finding | Status | Evidence | What remains |
|---|---|---|---|
| R1 baseline roles / onset-month entrants / originator availability | **Partly fixed** | Line 169 defines baseline by `participants_timing.orig == 1`, not month equality, and classifies same-month non-originator entrants by side. Line 170 states both acceptance-test projections and separates participants without fitted events. Line 171 adds `avail_orig` with probability exactly 0 outside the window. | Line 122 still gives the month-equality definition, and line 125 still gives the old acceptance test (N3). The W = 0 rule at line 171 contradicts stage 1, and the mask covers only the originator alternative, not the mirror-image unavailability of third-party A inside the window (N1). The outcome for initiating-side within-window entrants that stage 1 cannot classify is undefined (N2). |
| R2 opportunity q and sequential censoring | **Fixed in R2, contradicted elsewhere in §13 (partly fixed overall)** | Line 174 defines one-shot q as cumulative third-party incidence under the competing event (three-category auxiliary on complete-horizon rows) and monthly q as the cause-specific hazard with the mask. Line 175 gives the four-category offset `(λ ln q, λ ln q, 0, 0)`. Line 176 retains `h_O` in the sequential survival term and gives the correct `F_A(H)`. Retrospective and net-risk variants are separate families. | Line 126 still says originator entries "end the row (censored)" for opportunity training. Line 163 still says "Originator entries are censored in this variant". Both are the R2 defect verbatim, and nothing in the memo says R1–R3 override B1–C2 (N3). Minor P3s are listed below. |
| R3 single task table / outcome scale / unknowns | **Partly fixed** | Line 80 states that §13 supersedes §10 and that the task generator reads §13. Line 179 adds cohort, empty-side, horizon and two-month dimensions with estimator restrictions and a duplicate-cell key. Line 180 gives the recurrent return-to-risk rule and count contrasts. Line 182 gives the logical-bounds rule for unknown fatalities with an upper-bound sensitivity, and keeps the stage-1 pair-level label distinct. Line 179 now labels the truncated one-shot "descriptive observed-horizon", which closes the rev-4 B2 nit. | Line 181 fixes where the SD and percentiles are computed, but leaves open whether the within-year-standardized tie is the regressor and how it is centred. A zero within-year SD is not covered (N4). "Ineffective levels" in the duplicate key is not enumerated (P3). |

## New findings

### N1 — P2: the window semantics in R1 do not match the stage-1 code, and third-party A is unavailable inside the window

**Memo lines 123, 169, 171; `R/05_dump_events.R:138, 209, 218–229, 327`; `run_all.sh:151`.**

(a) **W = 0.** Line 171 says "With W = 0 the originator alternative exists only in the onset month." In stage 1, W = 0 is the default coding with *no* promotion at all. The recode block runs only `if (COORIG_W > 0L)` (R05:218), `run_all.sh:151` passes `--coorig_months` only when W > 0, and R05:138 says "W = 0 is the default coding." Under the stage-1 code the builder is told to call, a W = 0 cell has no originator events. The memo's mask would nonetheless open an originator alternative in month 0 with zero events, which produces separation, the failure R1 was meant to prevent. Alternatively, the builder follows stage 1 and silently departs from line 171. Either way, the implementer chooses. The fix is to state that under W = 0 there is no originator category and the cell is a three-category model (A, B, out), equivalent in category structure to the retrospective cohort but on the prospective cohort.

(b) **A inside the window.** Stage 1 promotes an initiating-side joiner deterministically when its definite-bound lag `lag_hi = entry_hi − onset_lo` is ≤ W (R05:209, 222, 226). Line 14 of the memo uses the same month index (`ceil(entry_hi − onset_lo)`). So every initiating-side entry with a fitted stage-1 event and month index ≤ W is an originator event, and §2 line 9 confirms this ("0 attack-side joiners enter within 1 month by construction"). Third-party A is therefore structurally empty in months 0…W, apart from the N2 residue. Line 171 masks the originator alternative after W but does not mask A before W+1. In the monthly flat model, the A × duration-interval terms (line 45, intervals 0–1, 2–3, 4–6) then have zero events for W ≥ 1. With W = 6, the A coefficients for all three early intervals are unidentified. The same holds for equation 2 of the sequential model (A vs B among third-party entrants in months ≤ W). Concrete failure: the W = 3 monthly reference-alternative cell estimates an A × (2–3) interval with no positive rows, and the fit either diverges or stops under the separation rule, depending on solver.

**Minimal resolution.** Add `avail_A` (0 in months 0…W on the stage-1 month index, 1 after), state the W = 0 case as in (a), and have the fixtures assert both masks. The alternative is to define a third-party A event inside the window (see N2) and accept a sparse but nonzero A in those months. That is a different estimand and must be named.

### N2 — P2: initiating-side within-window entrants with no fitted stage-1 event have no defined outcome

**Memo lines 123, 125, 169–170; `R/05_dump_events.R:226, 229, 327`.**

R05 classifies *fitted events* (`ev`, matched to participant spells through `prow`), not participants. Promotion is `early_j <- ev$grp == 2L & early_lag(...)` on event rows. A participant who enters on the initiating side within W but has no fitted stage-1 event can never be "promoted" by that code. Such participants exist, for example from stage-1 sample filters or missing covariates, and line 170 itself acknowledges them ("participants with no fitted stage-1 event are reported separately"). Line 169 then says an initiating-side non-originator entrant "is an originator competing event if stage 1 promotes it within W, otherwise a third-party A entry". Read literally, such a state becomes third-party A, so its stage-2 outcome category depends on whether stage 1 happened to fit an event for it. The §2 partition claim (line 12) fails, and the mask in N1(b) would be contradicted by real A events in masked months. "Reported separately" says where they are logged, not what row outcome they get.

**Minimal resolution.** Define promotion at the participant level with the stage-1 rule: initiating side by the case-reviewed map, `orig == 0`, definite-bound month index ≤ W, side-switch handling as in R05. Apply it to every participant, and keep the "call stage-1 code" requirement as the acceptance test on the fitted-event subset (line 170(b)). Give the count of participant-level promotions with no fitted event. Choosing instead to treat them as A, or to drop them, would be a different policy that must be stated, because it changes N1(b).

### N3 — P2: §13 still contains the rev-4 text that R1 and R2 replace, with no precedence rule inside §13

**Memo lines 3, 122, 125, 126, 163 versus 169, 170, 174, 176.**

Line 3's precedence rule orders §13 over §2–§10 only, and line 80 does the same for §10. Inside §13, four rev-4 sentences that R1/R2 were written to replace remain verbatim:
- line 122 defines baseline as "states whose spell starts in the onset month", where R1 (169) says "not by month equality";
- line 125 has the acceptance test against "the stage-1 dump's coalition-originator set", where R1 (170b) says coalition originators that were not original originators;
- line 126 says originator entries "end the row (censored)" for opportunity training, where R2 (174) says "never censored or dropped";
- line 163 says "Originator entries are censored in this variant", where R2 (176) says C2 is completed with a three-category equation 1.

A builder that implements B4 training labels from line 126, or C2 from line 163, reproduces exactly the rev-4 R2 estimand error. The previous review treated a §10/§13 conflict of the same kind as P2 because jobs might be generated from the wrong section. The same reasoning applies here.

**Minimal resolution.** Edit lines 122, 125, 126 and 163 to point to R1/R2, or add one sentence at line 118 saying that R1–R3 govern where they differ from B1–C2.

### N4 — P2: the contrast-scaling rule does not say what enters the model, and a zero within-year SD is unhandled

**Memo lines 160, 181.**

Line 160 says tie sum and tie max "are standardized (within-year SD) before contrasts", and line 181 says the SD is computed "over all reference-cohort rows of that dispute-onset year". Dividing by a year-specific SD is not a linear transformation of the pooled covariate. Whether the model is fitted on raw ties or on year-standardized ties therefore changes the coefficients, the fitted incidence and the contrast, not just the units. The memo never says which variable enters the utility, whether it is also centred within year, or, for monthly rows whose ties come from later calendar years, whether the onset-year SD is still used. Line 181's zero-SD rule applies to the pooled reference cohort, not to a single onset year.

Concrete failure: an early-century onset year contains one fatal dispute whose initiating-side originators have no positive clipped tie to any cohort member. Tie-sum-to-A then has within-year SD 0, so the standardized value is 0/0 for every row of that year. The builder must invent a policy (drop the year, set to 0, borrow a pooled SD), and each option changes either the risk set or the covariate.

**Minimal resolution.** State the regressor explicitly, for example "the utility uses raw ties; standardization is applied only to express the 25th–75th contrast", or "the utility uses x / SD_onset-year, uncentred". Fix the monthly analogue, and give a zero/undefined within-year SD rule, for example a fallback to the pooled SD with the affected years counted.

## P3 notes (do not block the freeze, but fix while editing)

- **Retrospective-cohort q.** Line 174 defines q on the prospective three-category auxiliary. It is not said whether retrospective cells use that q (an unconditional P(T)) or a q refit without originators (P(T | not O)). These are different offsets.
- **Truncated-horizon one-shot with λ = 1.** q is fitted on complete-horizon rows (line 174). Applying it to truncated rows in the line-130 sweep level, without a horizon term, is unstated. Declaring λ = 1 × truncated illegal would be the simplest fix.
- **Net-risk family.** Line 176 names it but not the computation. Presumably `S` omits `h_O` while `h_T` is taken from the fitted cause-specific model. One clause would fix it.
- **Duplicate-cell key.** Line 179 maps "ineffective levels to their reference value" without listing them, for example tier timing under one-shot (illegal, not ineffective), λ under likely-only, dynamic entrants under one-shot, or window under the retrospective cohort. A short table avoids the implementer deciding the collapse rules.
- **Recurrent family scope.** Line 180 does not say whether baseline (`orig == 1`) participants who exit can re-enter as third parties. That changes the rows.
- **Sequential × estimator.** Line 176's formula is monthly. The one-shot sequential (`P(T) × P(A | T)` from the three-category equation 1) is implied but not stated, and neither is its legality row.
- **Acceptance test (a)** at line 170 is close to tautological, because baseline is defined as `orig == 1`. It tests the join and the side assignment, which is still useful. Say so, so nobody mistakes it for an independent validation.
- **Late originators.** R05:211–213 keeps `orig == 1` participants who enter more than one month after the reference as originators "at onset by definition". Under R1 they are baseline and therefore in the month-0 roster (line 135) before their recorded entry. This is consistent with stage 1, but say so and count them.

## Decision

**NO-GO for freezing the analysis-ready risk set and legality table at revision 5. Infrastructure remains GO.** The remaining items are local. N1 and N2 together are one paragraph on window classification (participant-level promotion, `avail_A` and `avail_orig`, W = 0 as a three-category cell). N3 is four line edits or one precedence sentence. N4 is two sentences on the standardized regressor and the zero-SD fallback. None reopens M1–M5 or the R2 estimand, which is now correct where it is written.

Fixtures to add beyond the rev-4 list: a W = 0 cell asserting no originator category; a W = 3 monthly dispute asserting `avail_A = 0` in months 0–3 and `avail_orig = 0` from month 4; an initiating-side entrant in month 1 with no fitted stage-1 event, asserting its participant-level promotion under W ≥ 1; an onset year with zero within-year tie SD, asserting the fallback; and a check that the B4 auxiliary labels and C2 equation 1 keep originator rows as an observed category.
