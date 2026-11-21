# Stage 2 (joining) design memo — PI decisions recorded 2026-11-21; awaiting GPT review

Nothing in this memo is implemented. Every choice marked **Decision** needs PI approval, then a GPT review, before any code is written. The risk-set rules below are fixed in writing before any stage-2 result is seen, and the alternatives are run as sensitivity checks.

## 1. Question and claim
Stage 1 shows that coalition originators attack relationally isolated, materially stronger states, while bilateral initiators fight peers. Stage 2 asks which third parties enter an ongoing fatal dispute, and on which side. The recruitability argument predicts that third parties join the attacking side when they are tied to the attackers and when the target is isolated, and join the target's side when they are tied to the target. Defenders and attack-side joiners become stage-2 actors choosing a side, not stage-1 actors choosing a target.

## 2. Unit of analysis
**Decision 1 (PI: run both estimators, reported side by side).** (a) Monthly competing-risks hazard: a potential joiner k in an ongoing fatal dispute d, observed each calendar month from the end of the originator window until k enters or the dispute ends, with three outcomes per month (join the attacking side, join the target's side, stay out); month resolution matches the MID dates. (b) One-shot side choice: one row per (k, d), outcome "joined side A / side B / never", covariates measured at dispute onset. The hazard uses changing sides; the one-shot version is simpler and less sensitive to date coding. Agreement between them is part of the evidence.

## 3. Risk set
**Decision 2 (PI: separate likely from unlikely joiners).** The risk set is every system member not already a participant in d (stage-1 coalition originators, entered within the window, are stage-1 actors and are excluded; the risk set starts in the first month after the window, so the window choice moves states between the two stages consistently). Within it, likely and unlikely joiners are separated two ways, both reported:
- (a) **Rule-based relevance tiers, fixed in advance.** Likely = contiguous to any participant, in the same region as either side, a major power, or holding an ATOP pact with any participant; unlikely = everyone else. Estimated on the likely set alone, and on the full set with tier interactions on the tie and power terms.
- (b) **Estimated opportunity.** A first-stage logit of entering any fatal dispute on dispute-independent opportunity (distance and contiguity to the participants, major power, capability, regime), fitted out of sample by era, gives each potential joiner a propensity of being a plausible entrant. It enters the side-choice model as a fixed offset, the same construction as the reviewed stage-1 opportunity-weight offset (λ × ln q; λ = 0 reproduces the unweighted model). Reported as the residual side choice given opportunity, not as a causal correction.

**Sensitivities.** All MID onsets instead of fatal disputes; regional-only tier definition.

Size check. Roughly 1,000 fatal disputes, a median duration of a few months (long tails such as WWII), about 60–190 system members per year. Expect a few million potential-joiner months, which a logit handles easily.

## 4. Sides and timing
- Side A = the initiating side as defined in stage 1 (case-reviewed initiating-side rule). Side B = the target's side.
- Side membership in month m is taken as of the end of month m − 1, so a state's decision depends on the sides it could observe, not on its own entry or on simultaneous entries.
- Entries in the same month are treated as simultaneous; the month-level sensitivity is a two-month aggregation.
- A state that enters is removed from the risk set for the rest of the dispute (absorbing). **PI: re-entries ignored in the main specification; a sensitivity treats a re-entry after a gap as a new entry decision.**

## 5. Model
**Decision 3 (main).** Alternative-specific conditional logit over {join A, join B, stay out}, one choice set per potential-joiner month, with stay out as the reference alternative. Side-specific covariates enter the utility of their own side. Two-way cluster-robust inference by potential joiner and dispute (as in stage 1); pigeonhole bootstrap as the conservative check.

**Alternative.** Multinomial logit with the same covariates. The joint stage-1/stage-2 model with shared latent factors is the AME phase, not this memo.

## 6. Covariates (all lagged one year unless stated)
| Concept | Variable | Enters |
|---|---|---|
| Ties to each side (main test) | k's tie strength to the current members of side s, from the same multiplex layers and weighting as the main relational measure. **PI: report both aggregations**, the sum over members (total embeddedness with the side) and the strongest single link | utility of joining s |
| Target isolation (recruitability) | the original target's centrality (stage-1 measure) | utility of joining A |
| Own position | k's centrality | both join alternatives |
| Power | k's capability relative to side s's combined capability; side s's share of the dispute's combined capability (bandwagoning vs balancing) | utility of joining s |
| Geography | minimum distance and contiguity from k to members of side s | utility of joining s |
| Alliance | ATOP defense pact between k and any member of side s | utility of joining s; shown on and off because ATOP is inside the tie measure |
| Dispute state | months since onset; number of participants on each side; major power on each side | both join alternatives |
| Prior conflict | prior MIDs between k and members of each side | utility of joining s |

The tie-strength covariates need a new dyadic build (the stage-1 positions are monadic centralities). They come from the existing dyad-year layers file with the same layer combinations and weightings, so the stage-2 multiverse runs over layer combination × weighting (no centrality algorithm).

## 7. Hypotheses (stated before results)
- H1 (recruitment). The probability of joining the attacking side rises with k's ties to the attackers and falls with the target's centrality and with k's ties to the target.
- H2 (defense). The probability of joining the target's side rises with k's ties to the target.
- H3 (power). Joiners favour the stronger side (bandwagoning) or the weaker side (balancing); stated as a two-sided question, not a prediction.
- A null on H1 is informative: it would mean stage 1's isolation pattern does not come from who can be recruited.

## 8. Identification threats and how each is handled
| Threat | Handling |
|---|---|
| Sides change endogenously as states enter | sides as of month m − 1; simultaneous same-month entries treated as such; two-month aggregation sensitivity |
| Coarse dates | month resolution; window sensitivity 0–6 aligned with stage 1 |
| Stage 2 exists only where stage 1 produced a fight | stated as conditional on the dispute; joint model in the AME phase |
| Alliance ties inside the tie measure | alliance control on and off; no-ATOP tie measures |
| Wartime rewiring of ties | one-year lags; lag-5 sensitivity |
| Concentration in a few large wars | leave-one-dispute-out; drop WWI/WWII; history-chosen drop lists from stage 1 |

## 9. Robustness grid (planned)
Likely/unlikely separation (relevance tiers; opportunity offset at λ = 0, 1); all MIDs; timing (monthly hazard vs one-shot); windows 0–6; tie aggregation (sum vs maximum); re-entries ignored vs counted; alliance on/off; lags 1/5; permutation diagnostic shuffling tie values within dispute-month; measurement error that mimics alliance and contiguity confounding (both signs) plus random noise; leave-one-dispute-out; pigeonhole bootstrap on the main measure.

## 10. PI decisions (2026-11-21)
1. Estimator: run both the monthly competing-risks hazard and the one-shot side choice, reported side by side.
2. Risk set: separate likely from unlikely joiners, by fixed relevance tiers and by an estimated opportunity offset; both reported.
3. Sample: fatal MIDs main; all MIDs as sensitivity.
4. Ties: report both the sum over a side's members and the strongest single link.
5. Re-entries: ignored in the main specification; counted as new entries in a sensitivity.

## 11. Order of work after approval
GPT review of this memo → dyadic tie build and risk-set build with synthetic-data tests → main model → GPT audit → stage-2 multiverse and robustness grid on the laptop → results. The formal model memo covering both stages follows; the joint AME is last.
