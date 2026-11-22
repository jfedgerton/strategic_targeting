# Stage 2 (joining) design memo — revision 6 (2026-11-22): builder contract §13 with R1–R3 and N1–N4 (reviews `docs/REVIEW_stage2_memo_rev4.md`, `docs/REVIEW_stage2_memo_rev5.md`)

Nothing in this memo is implemented. Rules are fixed in writing before any stage-2 result is seen; alternatives run as sensitivities. Items marked **[PI]** were decided on 2026-11-22: the PI approved every recommended default (M1–M5). The PI also asked that every alternative be run as a sweep, not a single sensitivity (§10), and that the joint AME be the intended main model once it is built and validated (§12). Revision 1 (PI decisions of 2026-11-21) is kept in §11; GPT's findings M1–M5 are in `docs/REVIEW_stage1.md`. **Revision 4 adds §13, the builder contract answering GPT's B1–B5, C1 and C2; where §13 and §2–§10 differ, §13 governs.** PI 2026-11-22: "do both" (update stage 1, proceed with stage 2); every GPT choice is written with the recommended default and its alternatives as labelled sweep levels.

## 1. Question and claim
Stage 1 shows that coalition initiators attack relationally isolated, materially stronger states, while bilateral initiators fight peers. Stage 2 asks which third parties enter an ongoing fatal dispute, and on which side. The recruitability argument predicts more entry on the attacking side by states tied to the attackers and when the target is isolated, and more entry on the target's side by states tied to the target. Stage 2 describes these conditional associations among selected fatal disputes; it does not by itself identify recruitment as the cause of the stage-1 pattern (M5).

## 2. Stage boundary and who is decided where (M1)
Counts from the main 1-month frame (fatal MIDs): of 158 target-side joiner state-disputes, **45 (28%) enter within 1 month of onset**, in 28 of the 65 disputes with any target-side joiner. Attack-side entrants within the window are stage-1 coalition originators (0 attack-side joiners enter within 1 month by construction).

**[PI] M1 decision.** Recommended (a):
- (a) **Stage 2 starts at the dispute's onset month.** Every system member other than the dispute's originators is at risk from month 0. An attack-side entry within the window is a stage-1 originator decision: that state leaves the stage-2 risk set at entry and the event is recorded as "entered as originator", a competing outcome not analysed in stage 2. A target-side entry at any time, including within the window, is a stage-2 outcome. The two stages then partition every participant: originators (incl. window promotions) in stage 1, everyone else in stage 2.
- (b) Late entry only: stage 2 starts the month after the window, and the estimand is entry among non-participants surviving to that landmark. Loses the 28% of defenders who enter early. Run as a sensitivity under (a).
- Shared rules for both: reference date = dispute onset with the same definite-bound rule as stage 1 (onset_lo; entry counted in month ceil(entry_hi − onset_lo)); disputes ending before the landmark contribute no stage-2 rows under (b) and are counted. The two-month aggregation sensitivity bins months 0–1, 2–3, … on the same reference and recomputes the window boundary on that scale.

## 3. Estimators and the common cohort (M2)
PI decision (rev. 1): run both the monthly competing-risks model and the one-shot side choice.

**[PI] M2 decision.** Recommended:
- **Common cohort:** system members at onset, not originators, not already participants (the landmark cohort). Follow-up ends at the dispute's last participant spell month, or Dec 2014 (censored, flagged).
- **One-shot:** one row per cohort member × dispute; outcome = first entry on A, first entry on B, or no entry by the end of follow-up. Covariates frozen at onset (t−1 year); side rosters frozen at the onset originators. Censored disputes are flagged and dropped in a sensitivity.
- **Monthly model:** the same cohort, followed month by month with updated rosters (month m − 1) and covariates; states that enter the system during a dispute are a separate "dynamic entrant" group, included only in a monthly sensitivity.
- **Comparison:** convert the monthly model into side-specific cumulative incidence over the same follow-up (survival through earlier months included) and compare standardized predicted probabilities with the one-shot model. Coefficients are not compared directly. Differences can come from updated rosters and covariates, not only from date coding.
- Exits from the system end follow-up; missing or conflicting end dates are counted and dropped.

## 4. Sides, target set, ties and relevance (M4)
**[PI] M4 decisions.** Recommended:
- **Dispute key:** the MID dispute number. **Initiating side** = the stage-1 case-reviewed map for that dispute.
- **Original target set:** every originator on the target side at onset (15 of 522 fatal disputes have more than one; 18 have more than one attack-side originator). **Target centrality** = the mean of the set's stage-1 relational positions at t−1, fixed at onset; sensitivity = the minimum (most isolated member).
- **Current side members (monthly):** participants whose spell covers month m − 1. A state that exits leaves its side. A side switch is treated as the state's first side in the main model (switches are rare; counted), with a sensitivity censoring at the switch. Re-entries: ignored in the main model (absorbing), new decisions in the re-entry sensitivity (rev. 1).
- **Ties:** undirected (the layers are undirected). Built from the main measure's layers and weighting (igo+def+off+neu+nonagg+consul+dca+un, contiguity-residualized, clipped, rescaled within year), dyad edge weight at t−1; no link = 0. Aggregations: the sum over current side members and the strongest single link (rev. 1). Onset-frozen ties (from the year before onset) are the comparator for wartime rewiring.
- **Power:** k's CINC relative to the summed CINC of side s; side s's share of the dispute's combined CINC.
- **Relevance tiers (fixed in advance, frozen at onset):** likely = land- or sea-contiguous (≤ 400 miles) to any onset participant, same COW region as any onset participant (a state's own COW region; no multiple-region rule needed), a major power, or holding any ATOP pact in force at t−1 with an onset participant. Unlikely = everyone else. Frozen at onset so the likely-only estimate does not condition on a time-varying eligibility; a monthly-updated tier is a sensitivity.

## 5. Opportunity offset (M3)
PI decision (rev. 1): separate likely from unlikely joiners by the fixed tiers and by an estimated opportunity offset.

**[PI] M3 decision.** Recommended:
- **What q predicts:** the probability that k enters a given fatal dispute as a third party (not an originator), on either side. One-shot: over the dispute's follow-up. Monthly: a monthly entry hazard from the same features. Denominator = the §3 cohort.
- **Offset vector:** (λ ln q, λ ln q, 0) for (join A, join B, stay out). This shifts participation versus staying out and, by construction, cancels from the A-versus-B choice among joiners. It is reported as a participation-opportunity sensitivity, not as an adjustment of side choice. λ ∈ {0, 1}; λ = 0 reproduces the unweighted model.
- **Features (pre-decision only):** distance and contiguity to the onset participants (one-shot) or to the m − 1 roster (monthly), major power, CINC, regime. Never the eventual entrant list.
- **Out-of-sample folds:** six folds by era, exhaustive and disjoint (1816–1913, 1914–1918, 1919–1938, 1939–1945, 1946–1989, 1990–2014; §13 B4). Each dispute's q comes from a model fitted on the other eras, with all rows of held-out disputes excluded. Training outcomes are third-party entries only (stage-1 originators excluded).

## 6. Model (M5)
- **Main:** multinomial (alternative-specific) logit over {join A, join B, stay out}, stay out as reference, with join-alternative intercepts. Monthly model adds duration dependence as piecewise intervals (months 0–1, 2–3, 4–6, 7–12, 13+), interacted with the join alternatives.
- **k-level quantities** (own centrality, duration, dispute size) enter only interacted with the join alternatives, with separate A and B coefficients.
- **Side-specific quantities** (ties to s, power relative to s, distance and alliance to s, prior conflict with s) enter the utility of joining s, with separate coefficients for A and B.
- **[PI] alternative check.** Revision 1 listed "multinomial logit" as an alternative; it is the same model, so it is dropped. Replacement: a sequential (two-equation) logit, defined in §13 C2. It is not a nested logit with an inclusive-value parameter.
- Inference: two-way cluster-robust by potential joiner and dispute; pigeonhole bootstrap as the conservative check.

## 7. Covariates (lagged one year unless stated)
| Concept | Variable | Enters |
|---|---|---|
| Ties to each side | sum and max of k's ties to current members of s (§4) | join s |
| Target isolation | mean (min) centrality of the original target set, fixed at onset | join A (and join B, separate coefficient) |
| Own position | k's centrality | both joins, separate coefficients |
| Power | k's CINC / side s CINC; side s share of combined CINC | join s |
| Geography | minimum distance and contiguity to members of s | join s |
| Alliance | ATOP pact with any member of s; shown on and off because ATOP is inside the tie measure | join s |
| Dispute state | months since onset (piecewise); participants on each side; major power on each side | both joins |
| Prior conflict | prior MIDs between k and members of s | join s |

## 8. Hypotheses and what each result can show
- H1 (recruitment). P(join A) rises with ties to A and falls with the target set's centrality. Ties to B lower P(join A) only through the competing alternative unless a direct term is added; report as predicted-probability contrasts with uncertainty.
- H2 (defense). P(join B) rises with ties to B.
- H3 (power). Two-sided: bandwagoning (stronger side) vs balancing (weaker side).
- A null on H1 does **not** establish that stage 1's pattern is unrelated to recruitability: measurement, limited variation, the risk-set boundary, equilibrium sorting and the selected sample can all produce a null. A positive association does not identify recruitment as the cause of the stage-1 pattern.

## 9. Identification threats
| Threat | Handling |
|---|---|
| Sides change endogenously | rosters as of m − 1; same-month entries simultaneous; two-month sensitivity. Lagging prevents self-inclusion, not sorting into sides |
| Wartime rewiring of ties | onset-frozen ties as comparator; lag-5 sensitivity |
| Coarse dates | definite-bound month rule shared with stage 1; windows 0–6 |
| Selection | both stages condition on selected fatal disputes; stated; joint model in the AME phase |
| Alliance inside the tie measure | alliance control on/off; no-ATOP tie measure |
| Concentration in a few large wars | leave-one-dispute-out; drop WWI/WWII; stage-1 history lists |

## 10. Specification multiverse (PI 2026-11-22: multiple approaches and sweeps throughout)
**Superseded by §13 (revision 4–5): the reference cell, cohort, outcome categories, legal cells and coverage plan in §13 govern; the table below is kept only as the list of dimensions. The task generator reads §13's legality table, not this table.**
Every alternative below is a sweep dimension, reported as a distribution with the recommended default marked, not a one-off sensitivity. The defaults define the reference specification; the main measure is the stage-1 main measure.

| Dimension | Levels (default first) |
|---|---|
| Stage boundary (M1) | onset start; late entry after the window |
| Estimator (M2) | one-shot side choice; monthly competing risks (cumulative incidence) |
| Model form (M5) | flat three-alternative logit; nested (join vs stay out, then side) |
| Risk set (M3/M4) | likely tier only; full set with tier interactions; full set + opportunity offset λ = 1; full set, λ = 0 |
| Tier timing | frozen at onset; updated monthly (monthly estimator only) |
| Tie aggregation | sum over side members; strongest single link |
| Tie timing | lagged one year; frozen the year before onset |
| Target centrality | mean of the target set; minimum |
| Side switches | first side; censored at the switch |
| Re-entries | ignored; new decisions |
| Dynamic entrants | excluded; included (monthly only) |
| Alliance control | on; off |
| Lags | 1 year; 5 years |
| Window | 1 month; 0–6 |
| Sample | fatal MIDs; all MIDs; jointly fatal; fatality levels ≥ 2–6 |
| Relational measure | main; the plausible stage-1 measures (layers × weightings incl. sry) |

Scale. The full factorial is too large to cross with every measure. Plan (to be costed after the build): the full specification factorial on the main measure; every measure on the reference specification; the Phase A measure set crossed with the dimensions that change the estimand (stage boundary, estimator, risk set, tie aggregation). Each cell reports the H1–H3 predicted-probability contrasts, and the summary reports, per dimension, the share of cells whose interval excludes zero (a descriptive grid summary, not replications).

Diagnostics on the reference specification and its main alternatives: permutation shuffling tie values within dispute-month; measurement error mimicking alliance and contiguity confounding (both signs) plus random noise; leave-one-dispute-out and the stage-1 history-exclusion lists; pigeonhole bootstrap.

## 11. Revision 1 PI decisions (2026-11-21, unchanged)
1. Estimator: both the monthly competing-risks model and the one-shot side choice.
2. Risk set: likely vs unlikely joiners by fixed tiers and an estimated opportunity offset.
3. Sample: fatal MIDs main; all MIDs as sensitivity.
4. Ties: sum over a side's members and strongest single link.
5. Re-entries: ignored in the main specification; new entries in a sensitivity.

## 12. Order of work and the main model
GPT review of revision 4 → risk-set and dyadic tie build with synthetic tests → reference model → GPT audit → stage-2 specification multiverse on the laptop → formal model memo covering both stages → AME build (stage-2 AME, stage-1 AME with the target effect, joint model; GPT review and independent audit before each).

PI 2026-11-22: the joint AME is intended as the paper's main model. It becomes main only if it meets the pre-set criteria fixed before its results are seen (convergence, effective sample size, parameter recovery on simulated data, posterior predictive checks; `docs/ame_joint_sketch.md`). The separate stage-1 and stage-2 multiverses remain the transparent baseline and the robustness envelope it is judged against; if the AME fails the criteria, the separate models lead and the AME is reported as attempted.

## 13. Builder contract (revision 4; answers GPT B1–B5, C1, C2)
Defaults first; every alternative listed is a labelled sweep level (PI 2026-11-22). Nothing here is implemented. **Precedence: later subsections override earlier ones where they differ (N1–N4 over R1–R3 over B1–C2); §13 overrides §2–§10.**

### B1. Originator entries and the cohort
- **Reference cohort (prospective).** At risk from the onset month: every system member at onset that is not a baseline participant. Baseline participants are defined in R1 (stage-1 originator flag), not by month equality. Eventual window promotions are not excluded in advance.
- **Outcomes.** Monthly: four categories per state-month: third-party A, third-party B, originator entry, stay out. Originator entry = an initiating-side entry that stage 1 promotes to coalition originator (only possible in months inside the window W); any later initiating-side entry is third-party A. One-shot: four labels (A, B, originator, none within horizon). Originator entry is a competing event: it ends the state's follow-up and enters the survival term, so it shapes A and B cumulative incidence. H1–H3 are reported on third-party A and B incidence; the originator category gets its own intercepts and the same covariates, and is reported descriptively.
- **Alternative (sweep level, labelled retrospective):** cohort excluding every eventual window originator, three outcome categories. Conditional on not becoming an originator.
- **Window algorithm.** Defined in R1 and N2 below (participant-level application of the stage-1 promotion rule, with acceptance tests).
- **Opportunity training labels (B4):** see R2. Third-party entries are the positive category; originator entries are an observed competing category in a three-category auxiliary model, never censored.

### B2. Horizon, censoring and the late-entry cells
- **Estimand.** First third-party entry on A or B within the observed horizon H = min(dispute end month, Dec 2014, the state's system exit). Named "within observed horizon".
- **One-shot reference** uses rows with complete horizons (dispute ended by Dec 2014, no system exit during the dispute); truncated rows are counted and excluded. **Monthly** uses all rows with censoring at H. Sweep level: one-shot including truncated rows with log(H in months) on the join alternatives, labelled.
- **Comparison.** Monthly hazards are aggregated to cumulative incidence over each row's own H on the one-shot reference rows, using observed lagged rosters; this is incidence conditional on the realized side history, and is named so. Recurrent-entry cells use expected number of entries, never pooled with first-entry incidence.
- **Late-entry cells** (stage boundary = after the window): landmark cohort = states eligible and not participating at month W+1; the one-shot is initialized at the same landmark; covariates remeasured for the year before the landmark; the original target set stays the onset set.

### B3. Month-0 roster, empty sides, target set
- **Month 0:** side A = baseline initiating-side participants; side B = the original target set (baseline target-side participants). From month 1: the roster = participants active in month m−1. A window-promoted originator joins side A's roster from the month after its entry; eventual entrants are never in a roster before they enter.
- **Empty side later:** default carry forward the last nonempty roster for that side, flagged; sweep level: follow-up ends when either side empties. Side power share and k/side capability are computed on the roster in use; with a carried roster they are never 0/0.
- **Target set** must be nonempty and resolvable from the participant file; disputes without one are listed and dropped. Mean and minimum centrality are then always defined.
- **Dynamic entrants** (monthly sweep level only): eligible from their first full month in the system; tiers and covariates computed at that month's year; baseline (onset-frozen) quantities are not imputed.

### B4. Opportunity folds
- Six disjoint folds by the dispute's onset month: 1816–1913, 1914–1918, 1919–1938, 1939–1945, 1946–1989, 1990–2014. Every row of a dispute takes its onset fold, so a dispute spanning eras never trains its own q. This is retrospective cross-fitting (other-era data, including later eras), not a historical forecast.
- q is bounded to [1e−4, 1 − 1e−4] before the log; a failed or separated auxiliary fit stops the build; features are numeric (no unseen factor levels by construction). These policies are fixed now, before any main fit.

### B5. Dispute-level samples
- Defined on the dispute (MID number), not on stage-1 pairs: fatal = GML dispute fatality > 0 (the stage-1 `disp_fatal` source); level ≥ k = dispute fatality level ≥ k; **jointly fatal** = at least one participant with fatalities on each side (side-level, from the participant fatality codes). Unknown codes are dropped and counted. Every eligible third party in a selected dispute is in its cohort.

### C1. Legal sweep cells
| Crossing | Rule |
|---|---|
| One-shot × re-entries as new decisions | Illegal; re-entries exist only in the monthly recurrent family (expected entries), reported separately |
| One-shot × monthly-updated tiers / dynamic entrants | Illegal (machine-readable legality table enforced by the builder) |
| Absorbing first entry × censor at side switch | Censoring applies to the switching state's membership in *other* states' rosters; cells identical to the reference are de-duplicated |
| Updated tier × full set without tier terms | Duplicate; dropped |
| Lag 5 × onset-frozen ties | Ties from onset − 5; other covariates follow the lag level |
| Windows | Seven unique values 0–6; W = 1 is the reference |
| Offset λ | λ = 0 vs 1 compared on identical full-set rows; likely-only vs full set reported as a population change |

**Reference specification (single cell):** prospective cohort (B1) with the originator competing event; one-shot within observed horizon on complete-horizon rows; flat four-category logit; full set with tier interactions, λ = 0; tiers frozen at onset; tie sum; ties lagged one year; target mean; first side; re-entries ignored; dynamic entrants excluded; alliance control on; window 1; fatal MIDs; main stage-1 measure; carry-forward rosters.
**Scale.** The literal crossing is ~570,000 main-measure cells; the builder first writes the legal task list, then it is costed and capped, and any reduced design is fixed in writing before estimates are examined. Phase-A crossings are a *selected interaction subset* (stage boundary, estimator, risk set, tie aggregation, sample, window), not a claim to cover every estimand-changing dimension.
**Contrasts.** H1–H3 = difference in predicted cumulative incidence (within H) when the covariate moves from its 25th to its 75th percentile, others at observed values; scaling rules in R3 (revision 6). Recurrent-family contrasts are differences in expected entries. Uncertainty by two-way cluster bootstrap. Summaries are reported within estimand families (first-entry incidence; recurrent expected entries; retrospective cohort) with requested / completed / valid counts.

### C2. Sequential alternative
Superseded by R2: the participation equation has three categories (third-party entry / originator entry / stay out), with the originator entry as an observed competing outcome; the side equation is A vs B among third-party entrants; incidence follows the R2 formula including h_O.

### AME language (sketch corrected in `docs/ame_joint_sketch.md`)
Sharing b_j models one source of dependence between target choice and joining; it does not correct selection into conflict or fatal disputes. The AME acceptance gate (likelihood, identification constraints and priors, R-hat/ESS/divergence thresholds, recovery and SBC tolerances, PPC targets, held-out folds, allowed remedial runs) is written as a separate document and fixed before any AME result is seen; success never depends on θ's sign, H1 significance, or agreement with the baseline.

### R1. Baseline roles, onset-month entrants, originator availability (revision 5)
- **Baseline participants** are defined by the stage-1 original-role field, not by month equality: participants with the dispute-level originator flag (`participants_timing.orig == 1`, the field stage 1 uses for originator events). A non-originator entering in the onset month is not baseline: on the target side it is a stage-2 entry decision (early defender); on the initiating side it is an originator competing event if stage 1 promotes it within W, otherwise a third-party A entry. Within-month order is never inferred; every non-originator same-month entrant is counted in the build log by side.
- **Acceptance tests (projections stated).** (a) Baseline roles: per dispute, the builder's baseline set equals the set of participants with `orig == 1` on each side. (b) Promotions: per dispute and window 0–6, the builder's initiating-side promoted set equals the stage-1 R05 promoted set (coalition originators that were not original originators). Both on the same dispute/source and date rules; participants with no fitted stage-1 event are reported separately, not silently matched.
- **Availability mask.** Each monthly row carries `avail_orig`: 1 only in months inside the promotion window of an initiating-side entry (months 0 … W by the stage-1 month index), else 0; unavailable alternatives get probability exactly 0 (dropped from the row's choice set). With W = 0 the originator alternative exists only in the onset month. A dispute with no possible promotion has no originator alternative at all.

### R2. Opportunity q and the sequential alternative under the competing event (revision 5)
- **One-shot q** = cumulative third-party incidence within H with the originator entry as a competing event: an auxiliary three-category logit (third-party entry / originator entry / none) on complete-horizon rows; q = P(third-party entry). **Monthly q** = the cause-specific third-party hazard from a three-category monthly logit (third-party / originator / out) with the same availability mask. Originator entries are observed competing outcomes in both, never censored or dropped.
- **Offset vector** for the flat four-category model: (λ ln q, λ ln q, 0, 0) for (A, B, originator, out). q is exclusively third-party opportunity.
- **Sequential alternative (C2, completed):** equation 1 is a three-category participation logit (third-party entry / originator entry / stay out; availability mask; offset λ ln q on third-party entry only); equation 2 is A vs B among third-party entrants. Incidence F_A(H) = Σ_m S(m−1) h_T(m) P(A | T, m), S(m−1) = Π_{r<m} [1 − h_T(r) − h_O(r)], with h_O from equation 1. The retrospective cohort and a net-risk variant (competing event eliminated) are separate, labelled estimand families, never pooled with the prospective reference.

### R3. One task table, outcome scales, unknown fatalities (revision 5)
- §13 supersedes §10 (note added there). The legality table also covers: cohort (prospective / retrospective; retrospective has no originator category), empty-side rule (carry-forward / end follow-up; both estimators), one-shot horizon (complete / truncated with log H, labelled descriptive observed-horizon, one-shot only), two-month aggregation (monthly only). Duplicate-cell key = the tuple of all dimension levels after mapping ineffective levels to their reference value; tasks with an identical key run once.
- **Recurrent family:** return to risk the month after a spell ends (gap ≥ 1 month), both sides available again; contrasts are differences in expected number of entries within H, reported only within the recurrent family.
- **Contrast scaling (revision 6, replaces the within-year SD rule):** tie aggregates enter every model raw (on the measure's own [0,1] edge scale), never standardized by year. Contrasts move each covariate from its 25th to its 75th percentile, computed once over the one-shot reference cohort (one row per state-dispute at the reference landmark) and applied unchanged to monthly predictions; percentile contrasts are scale-free, so no SD is used and no year-level SD can be zero. A covariate whose 25th and 75th percentiles coincide gets no contrast (reported as not computable, counted).
- **B5 unknown fatalities (logical bounds):** a dispute qualifies as jointly fatal if each side has at least one known fatal participant; it is excluded if either side's participants are all known nonfatal; otherwise (one side has no known fatal participant but has unknown codes) it is unresolved: excluded from the reference and counted, and included in a labelled upper-bound sensitivity. The stage-1 pair-level "jointly fatal" label stays distinct from this side-level selection.

### N1–N4 (revision 6; independent review of revision 5)
- **N1. Window availability, matched to stage 1.** W = 0 means no promotion in stage 1 (R05), so with W = 0 there is no originator alternative anywhere. For W ≥ 1, in months 0 … W an initiating-side entry is always an originator entry (stage 1 promotes every such entrant), so third-party A is unavailable (masked, probability 0) in those months and available from month W + 1. Third-party B and stay out are available in every month. The A × duration terms start at month W + 1; the originator alternative has its own intercept and no duration terms. Each monthly row carries `avail_A`, `avail_B`, `avail_orig`. **The mask is prospective (GPT recheck of revision 5):** it is computed from the declared stage-1 eligibility and date conventions and pre-decision information only, by asking what an A entry in this month would be classified as; it applies identically to every eligible candidate, whatever that candidate's eventual side or entry, and never uses observed future entries or whether the dispute eventually had a promotion. A dispute with no observed promotion still offers the originator alternative in months 0 … W when W ≥ 1. W = 0 keeps stage 1's meaning (no window promotion), so an onset-month non-original A entrant is third-party A at W = 0 and an originator entry at W ≥ 1.
- **N2. Promotion is a participant-level rule.** The builder applies the stage-1 promotion criteria (initiating side by the case-reviewed map, entry within W of the dispute's first originator by the definite bound) to every participant, whether or not stage 1 fitted an event for it. Acceptance test: on participants with a fitted stage-1 event, the participant-level promoted set equals the stage-1 R05 promoted set exactly, per dispute and window; participants without a fitted event are classified by the same rule and counted separately in the build log. A stage-2 category never depends on whether stage 1 fitted an event.
- **N3. Stale text.** The revision-4 sentences that conflicted with R1–R2 (baseline by month equality; originator entries censored in q training and in the sequential model) are replaced in place (B1 bullets, C1 contrasts line, C2), and the precedence rule above is stated.
- **N4. Tie scale and contrasts.** Raw tie aggregates in every model; percentile contrasts as in the revised R3 scaling bullet; no year-level standardization.

### B5 amendment (PI 2026-11-22, after the GPT build review)
- **Reference population:** "primary disputes represented by stage-1 fatal-event flags" (`events.rds in_fmid`; 560 disputes). B5's premise that GML dispute fatality is the stage-1 source was wrong (the sets overlap on 440 disputes). **Labelled alternative:** GML dispute fatality > 0 (541 after availability exclusions). ID lists are written by the builder (`fatal_sets.csv`).
- **Endpoints:** participant exits come from GML participant end months, capped at December of each stage-1 spell's last year (a GML episode never bridges a gap between spells); with no GML match, December of the spell's last year (imputed, flagged). Dispute end = latest participant exit; `end_imputed` flags disputes whose end is imputed. Sensitivity: reference rows excluding disputes with an imputed end.
- **Missing inputs:** side aggregates are NA when any member's input is missing; ties are NA (flagged) for years outside the layer coverage; rows with a missing reference covariate leave the reference estimation cohort with a recorded reason.

### Fixtures required before the freeze (from both rechecks of revision 5)
Exact row and outcome assertions for: a same-month non-original B entrant retained as B; one eligible same-month A entrant classified third-party at W = 0 and originator at W = 1; originator probability exactly 0 after its window; two candidates with identical pre-decision eligibility but different eventual outcomes receiving the same mask; a dispute with no observed promotion still offering the originator option at W ≥ 1; an initiating-side entrant within W without a fitted stage-1 event (participant-level rule, N2); the same masks in the auxiliary q model and the sequential first equation; plus the earlier set (initial and empty rosters under both policies, system exit, administrative truncation with the matched complete-horizon comparison, a 1914 dispute spanning folds, mixed side-level fatalities with unknown codes, a recurrent return after the declared gap, and a positive originator hazard showing its survival contribution in sequential incidence).
