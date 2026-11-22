# Stage 2 (joining) design memo — revision 3 (2026-11-22): M1–M5 approved; specification multiverse

Nothing in this memo is implemented. Rules are fixed in writing before any stage-2 result is seen; alternatives run as sensitivities. Items marked **[PI]** were decided on 2026-11-22: the PI approved every recommended default (M1–M5). The PI also asked that every alternative be run as a sweep, not a single sensitivity (§10), and that the joint AME be the intended main model once it is built and validated (§12). Revision 1 (PI decisions of 2026-11-21) is kept in §11; GPT's findings M1–M5 are in `docs/REVIEW_stage1.md`.

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
- **Out-of-sample folds:** five folds by era (1816–1913, 1919–1938, 1939–1945, 1946–1989, 1990–2014). Each dispute's q comes from a model fitted on the other eras, with all rows of held-out disputes excluded. Training outcomes are third-party entries only (stage-1 originators excluded).

## 6. Model (M5)
- **Main:** multinomial (alternative-specific) logit over {join A, join B, stay out}, stay out as reference, with join-alternative intercepts. Monthly model adds duration dependence as piecewise intervals (months 0–1, 2–3, 4–6, 7–12, 13+), interacted with the join alternatives.
- **k-level quantities** (own centrality, duration, dispute size) enter only interacted with the join alternatives, with separate A and B coefficients.
- **Side-specific quantities** (ties to s, power relative to s, distance and alliance to s, prior conflict with s) enter the utility of joining s, with separate coefficients for A and B.
- **[PI] alternative check.** Revision 1 listed "multinomial logit" as an alternative; it is the same model, so it is dropped. Recommended replacement: a nested (two-step) version, join vs stay out, then side among joiners, which relaxes the equal-substitution restriction between the two join options.
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
GPT review of this revision → risk-set and dyadic tie build with synthetic tests → reference model → GPT audit → stage-2 specification multiverse on the laptop → formal model memo covering both stages → AME build (stage-2 AME, stage-1 AME with the target effect, joint model; GPT review and independent audit before each).

PI 2026-11-22: the joint AME is intended as the paper's main model. It becomes main only if it meets the pre-set criteria fixed before its results are seen (convergence, effective sample size, parameter recovery on simulated data, posterior predictive checks; `docs/ame_joint_sketch.md`). The separate stage-1 and stage-2 multiverses remain the transparent baseline and the robustness envelope it is judged against; if the AME fails the criteria, the separate models lead and the AME is reported as attempted.
