# Stage 2 (joining) design memo — DRAFT for PI approval (2026-11-21)

Nothing in this memo is implemented. Every choice marked **Decision** needs PI approval, then a GPT review, before any code is written. The risk-set rules below are fixed in writing before any stage-2 result is seen, and the alternatives are run as sensitivity checks.

## 1. Question and claim
Stage 1 shows that coalition originators attack relationally isolated, materially stronger states, while bilateral initiators fight peers. Stage 2 asks which third parties enter an ongoing fatal dispute, and on which side. The recruitability argument predicts that third parties join the attacking side when they are tied to the attackers and when the target is isolated, and join the target's side when they are tied to the target. Defenders and attack-side joiners become stage-2 actors choosing a side, not stage-1 actors choosing a target.

## 2. Unit of analysis
**Decision 1 (main).** A potential joiner k in an ongoing fatal dispute d, observed each calendar month from the end of the originator window until k enters or the dispute ends. This is a discrete-time competing-risks design with three outcomes per month (join the attacking side, join the target's side, stay out). Month resolution matches the MID dates.

**Alternative (sensitivity).** One row per (k, d) with the outcome "ever joined side A / side B / never", covariates measured at dispute onset. Simpler, but it ignores that sides change as states enter.

## 3. Risk set
**Decision 2 (main).** Every system member in the month that is not already a participant in d. The stage-1 coalition originators (entered within the window) are stage-1 actors and are not in the stage-2 risk set; the risk set starts in the first month after the window, so the window choice (0–6 months, 1 month main) moves states between the two stages consistently.

**Sensitivities.** (a) States politically relevant to either side (contiguous to any participant or a major power); (b) states in the same region as either side; (c) all MID onsets instead of fatal disputes.

Size check. Roughly 1,000 fatal disputes, a median duration of a few months (long tails such as WWII), about 60–190 system members per year. Expect a few million potential-joiner months, which a logit handles easily.

## 4. Sides and timing
- Side A = the initiating side as defined in stage 1 (case-reviewed initiating-side rule). Side B = the target's side.
- Side membership in month m is taken as of the end of month m − 1, so a state's decision depends on the sides it could observe, not on its own entry or on simultaneous entries.
- Entries in the same month are treated as simultaneous; the month-level sensitivity is a two-month aggregation.
- A state that enters is removed from the risk set for the rest of the dispute (absorbing). Re-entries are ignored in the main specification.

## 5. Model
**Decision 3 (main).** Alternative-specific conditional logit over {join A, join B, stay out}, one choice set per potential-joiner month, with stay out as the reference alternative. Side-specific covariates enter the utility of their own side. Two-way cluster-robust inference by potential joiner and dispute (as in stage 1); pigeonhole bootstrap as the conservative check.

**Alternative.** Multinomial logit with the same covariates. The joint stage-1/stage-2 model with shared latent factors is the AME phase, not this memo.

## 6. Covariates (all lagged one year unless stated)
| Concept | Variable | Enters |
|---|---|---|
| Ties to each side (main test) | k's tie strength to the current members of side s, from the same multiplex layers and weighting as the main relational measure (sum and maximum over members) | utility of joining s |
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
Risk-set rules (main, politically relevant, regional, all MIDs); timing (monthly hazard vs one-shot); windows 0–6; tie aggregation (sum vs maximum); alliance on/off; lags 1/5; permutation diagnostic shuffling tie values within dispute-month; measurement error that mimics alliance and contiguity confounding (both signs) plus random noise; leave-one-dispute-out; pigeonhole bootstrap on the main measure.

## 10. Open decisions for the PI
1. Monthly competing-risks hazard (main) vs one-shot side choice (main)?
2. Main risk set: all system members, or politically relevant to either side?
3. Fatal disputes only in the main specification, all MIDs as sensitivity?
4. Tie-strength aggregation: sum (total embeddedness with a side) or maximum (strongest single link)?
5. Treat re-entries as new entries, or ignore them?

## 11. Order of work after approval
GPT review of this memo → dyadic tie build and risk-set build with synthetic-data tests → main model → GPT audit → stage-2 multiverse and robustness grid on the laptop → results. The formal model memo covering both stages follows; the joint AME is last.
