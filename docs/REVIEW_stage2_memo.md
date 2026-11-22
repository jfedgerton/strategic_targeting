# Stage-2 memo review — revision 3

22 November 2026. **NO-GO for the analysis-ready risk-set build as currently specified.** Revision 3 substantially improves the design and resolves the earlier interpretive overclaims. The remaining blockers are concrete: the early-originator competing event, first-month membership, censoring/one-shot outcome rules, missing opportunity-fold years, and the sample definition for jointly fatal disputes. These require short written decisions before the builder can have unambiguous acceptance tests. They do not require abandoning either estimator or the proposed multiverse.

The conditional AME-main commitment is defensible in principle, but its acceptance criteria are still proposals rather than a fixed release gate. The joint sketch also retains selection and sampling claims that need correction before AME implementation.

## Scope

Followed `docs/PROMPT_gpt_review_stage2_memo.md`; read all of `docs/stage2_design_memo.md` revision 3 and `docs/ame_joint_sketch.md`, using the prior M1–M5 findings as the comparison. No stage-2 implementation is present or inferred. Code and historical classifications were not re-audited.

Source: `/Users/jarededgerton/Downloads/targeting_multiverse_stage2memo.zip`, SHA-256 `2e7df79f0c6fe3a6b00fdfa3099c38497036440dbfe329cfe9e6da360d520bbe`. Line references are one-based lines in that package. **Only this report was written. No R code or project runner was executed.** The factorial count below is independent Python arithmetic, not an executed sweep. I consulted the official Stan reference for the proposed LKJ prior.

## 1. M1–M5 dispositions

| Finding | Status | Resolved and remaining behavior |
|---|---|---|
| M1: early defenders omitted | **Partly fixed** | Memo 12 explicitly starts at onset and includes early target-side entries: the original defender omission is corrected. But early attack-side promotions become a competing “originator” outcome while the model at 45 still has only A/B/out. Membership in the baseline exclusions also depends on what “originators” means at 20. See B1. |
| M2: incomparable cohort/horizon | **Partly fixed** | Lines 20–24 establish a common onset cohort, separate dynamic entrants, common follow-up and cumulative-incidence comparison rather than coefficient comparison. Good. The late-entry alternative does not explicitly replace this with a common landmark cohort; censoring and system exits are not reconciled with the one-shot labels. See B2. |
| M3: undefined offset and leakage | **Partly fixed** | Lines 39–42 define entry opportunity, different time scales, the offset vector, pre-decision rosters and dispute-level holdouts. The offset is correctly described as changing participation rather than conditional A/B odds. The eras omit 1914–1918; originator events, fold assignment for spanning disputes and boundary predictions remain undefined. See B4. |
| M4: target set, roster and relevance | **Partly fixed** | Lines 28–33 specify the dispute key, mean/min target centrality, active spells, undirected edges and onset-frozen tiers. These are major improvements. Month 0 has no specified pre-month roster; empty sides make some covariates undefined; the minimum target-set rule needs a nonempty onset target set. See B3. |
| M5: model and identification claims | **Partly fixed** | Lines 45–47 now define join intercepts, duration terms and alternative interactions. Lines 64–75 properly qualify H1, nulls, sorting and selected fatal disputes. The old redundant MNL alternative is removed. The proposed “nested (two-step)” replacement still needs its actual likelihood/utility restrictions; sequential logits and a nested logit are not interchangeable labels. See C2. The AME sketch separately retains the old overclaims. |

## 2. Definitions needed before building the analysis data

### B1 — P2: the originator competing outcome is not in either likelihood

**Memo 12, 20–23, 39, 45.** Consider a state entering A within W that stage 1 promotes to coalition originator. There are two different readings:

- Exclude all eventual window-originators from the baseline cohort. Then the cohort uses future participation to define eligibility and the proposed competing event never appears. This is a retrospective cohort conditional on not becoming an originator, not the stated prospective nonparticipant risk set.
- Exclude only the baseline originators, retain the state until its early entry, and record the promised originator event. Then the monthly process has four outcomes: third-party A, third-party B, originator entry, and no entry. The one-shot also needs an originator outcome or an explicitly conditional population. It cannot classify that entry as ordinary A, never, or an unexplained missing row.

The competing originator hazard contributes to survival and therefore to A/B cumulative incidence. Merely dropping its entry month while fitting a three-category monthly likelihood does not specify that probability. Choose and write the population/likelihood explicitly, including how the opportunity training labels treat this event. The revision's promised partition of participants is useful bookkeeping, but does not itself define the statistical risk set.

Also state the exact shared originator-window algorithm. A prose reference to onset_lo/entry_hi at line 14 should be backed by equality tests against stage 1's actual reference and definite bounds, including unknown dates and the two-month binning. Do not silently replace the established window convention with a similar-looking elapsed-month formula.

### B2 — P2: make censoring and the late-entry comparison executable

**Memo 13–14, 20–24, 84–85.** For late-entry cells, define the landmark cohort as states eligible and nonparticipating at that landmark, and initialize the one-shot model at the same landmark. Decide explicitly whether covariates and the target set remain onset-frozen or are remeasured there. A boundary change plus unchanged onset cohort would include early entrants in one estimator but not the other.

For the reference cohort, distinguish dispute completion from administrative truncation and system exit. If a dispute continues beyond December 2014, “no entry by observed follow-up” is observable, but “never entered this dispute” is not. If the target is incidence by actual dispute completion, the one-shot cannot treat censored nonentries as fully observed never outcomes; flagging them and dropping them only in a sensitivity does not solve the reference likelihood. Alternatively define an observed-horizon estimand explicitly and standardize monthly predictions to exactly that horizon. Specify what happens for individual system exits, not just censored disputes.

Fix the standardization population, individual horizons and covariate histories used to aggregate hazards. Using observed rosters for those predictions yields incidence conditional on the realized side history, not an unconditional simulation of coalition evolution; that is acceptable if named. A nested or recurrent alternative may require a different comparison summary.

### B3 — P2: month-zero roster and empty-side behavior are missing

**Memo 12, 21–22, 29–33, 54–61.** Side membership comes from month m−1, yet observation starts in month 0. Before the dispute starts there may be no active participants. Literal application gives empty sides: a maximum tie/minimum distance is undefined, k's capability divided by side capability divides by zero, and side power shares can be 0/0.

Specify a baseline roster exception, such as the original participants present at the onset decision, followed by lagged active rosters. Define what “at onset” includes when dates are only monthly and entries are simultaneous. Do not populate the baseline roster with all eventual window-promoted originators: that would expose the model to entrants who have not yet appeared.

Write rules for later empty sides as participants leave: unavailable join alternative, last roster, or explicit empty-roster covariates are different models. Define max/sum, geography and power consistently. Require a nonempty original target set and list/drop unresolved cases rather than allowing a default empty mean. Dynamic entrants need their own baseline eligibility and available lagged covariates; the onset-frozen tier cannot be assumed available for every state not yet in the system.

### B4 — P2: opportunity folds leave WWI years unassigned

**Memo 39–42.** The eras jump from 1816–1913 to 1919–1938. There is no fold for 1914–1918. A dispute originating in 1914 therefore has no prescribed training exclusion or out-of-sample q. This is a concrete coverage omission, especially consequential given the large-war concentration concern.

Make the fold map exhaustive and disjoint. Assign each entire dispute to a fold by a specified reference date, so a dispute spanning eras never trains its own predictions. The proposed other-era training explicitly allows future-era information; label it retrospective cross-fitting, not a historically available forecast. That choice is not itself an unresolved leak if the whole held-out dispute is excluded.

Specify q behavior for separated/failed auxiliary fits, unseen factor levels and exact/near-zero probabilities: log q needs a declared positive bound or another explicit policy. Freeze these policies rather than choosing them after seeing main-fit results. Resolve B1 before defining which originator rows are absent from training and which contribute a competing event.

### B5 — P2: “jointly fatal” has no dispute-level selection rule

**Memo 98; dispute unit at 28.** Stage 1's jointly fatal sensitivity is attached to a directed event/pair through the two states' dispute-level fatality flags. A stage-2 unit is k,d. A multilateral MID can have both jointly fatal and non-jointly-fatal stage-1 pairs. Specify whether stage 2 includes a dispute when any eligible originator pair qualifies, every pair qualifies, or some side-level condition holds. These are different samples. Do not borrow the first stage-1 row encountered as the dispute's classification.

Map the definition to the same underlying dispute/role universe and unknown-code policy. Likewise define severity ≥2–6 at the MID level and document whether the risk cohort includes all eligible third parties within each selected dispute. This belongs in the builder contract, not just the fit grid.

## 3. Specification multiverse

### C1 — P2: generate a constrained design, not the literal Cartesian product

**Memo 82–101.** The scale plan is sensible in structure, but the legal cells must be fixed first:

| Dimension/crossing | Required restriction or label |
|---|---|
| One-shot × re-entries as new decisions | Incompatible with one row per state–dispute and a first-entry outcome. Either exclude this crossing or define a different recurrent episode estimator and label its estimand. |
| One-shot × updated monthly tier / dynamic entrants | Already marked monthly-only; enforce this in a machine-readable legality table. |
| First-entry/absorbing × censor at a later side switch | Usually no effect on that state's outcome because follow-up already stopped at its first entry. It can affect subsequent rosters for other states; specify which object is censored. With onset-frozen one-shot rosters it may be a duplicate cell. |
| Recurrent entries × cumulative first-entry incidence | Not the same summary. Define first-event incidence, expected number of entries, or episode risks; do not pool them as interchangeable probabilities. |
| Updated tier × full set, no tier terms | No effective change unless tier affects offsets/interactions. Avoid duplicate cells inflating the grid summary. |
| Lag 1/5 × ties frozen before onset | Define onset−1 versus onset−5, and which other covariates change. “Frozen the year before onset” cannot simultaneously mean a five-year lag without a rule. |
| Window 1; 0–6 | Use seven unique values, not two copies of W=1. Recompute both stage membership and cohorts consistently. |
| Risk-set levels | These bundle population changes with model/offset changes. Compare full-set λ=0 versus λ=1 on the same rows to isolate the offset sensitivity; likely-only versus full-set differences also change the target population. |

The reference combination must be listed as a single concrete specification, including originator handling, severity rule, missingness, censoring and probability contrasts. It is not enough to say “default first” while those definitions remain unresolved.

### C2 — P2: define the nested alternative, including its offset

**Memo 48, 86.** A sequential join/stay logit followed by a side-among-joiners logit is a well-defined alternative once both equations and predictors are specified. It is not automatically a conventional nested-logit random-utility model with an inclusive-value parameter. Choose the formulation and state its restrictions; either preserves the approved two-step alternative.

The flat model's `(λ ln q, λ ln q, 0)` becomes a participation-equation term in a sequential model. Specify the analogous term and treatment of the originator competing event. Do not add it to the conditional A/B equation, where it cancels. Specify whether all predictors enter both equations or only their designated components.

### Scale and summaries

The three-tier allocation—main-measure specification grid, all measures at the reference specification, Phase A crossed with selected dimensions—is coherent as a coverage plan, not a fully crossed robustness claim. The named Phase-A dimensions at line 101 omit other population/estimand changes, notably sample, window, recurrent entries, dynamic entrants and possibly updated eligibility. Either expand that set or explicitly call it a selected interaction subset; do not claim it covers all estimand-changing dimensions.

A literal count using seven windows, eight sample levels and only the stated monthly-only tier/dynamic restrictions is about **573,440 main-measure cells**, before pruning recurrent-entry incompatibilities, redundancies and adding diagnostics. This is an illustrative upper grid count, not a runtime estimate. Cost and cap the valid task list before promising the full factorial, and predefine any reduced design before estimates are examined.

Define H1–H3 contrasts now: intervention/reference values, horizon, standardization cohort, weighting and uncertainty. A unit change in a side-sum tie is not comparable to a unit change in a strongest link without a declared scale. Monthly and one-shot probabilities require the common horizon; recurrent incidence requires a separate summary. Report distributions within estimand families, with requested/completed/valid counts, not one significance share mixing different populations and outcomes. Conditional dimensions otherwise change the weighting of cells in each marginal share.

## 4. AME as conditional main model

**Memo 113–115; sketch 37–48.** The commitment is defensible if it means an objectively qualified AME leads regardless of the direction or significance of H1–H3, with separate models and failed attempts still reported. It is not yet operational: the sketch calls thresholds “proposed” and leaves model choices open.

Before examining AME results, fix:

- **Likelihood and estimand:** selected dispute universe, B1 competing event, first/recurrent entry policy, weights, follow-up, multiple-target aggregation of b and centrality, and whether the joint main is monthly or one-shot. The sketch still has a single j* despite the revised target set. Specify the shared factor dimensions and temporal evolution, and whether side-latent means use observed lagged rosters.
- **Identification and priors:** scale/location constraints for b, u, v and Λ; rotation/sign invariances; AR(1) support/initial distribution; K selection; prior families/scales and fixed sensitivity ranges. “LKJ prior on Λ” is incomplete: LKJ supports correlation matrices, not an unrestricted multiplicative loading matrix. Put it on a defined covariance/correlation component or choose loading priors appropriate to Λ. See the [official Stan LKJ reference](https://mc-stan.org/docs/functions-reference/correlation_matrix_distributions.html).
- **Acceptance gates:** exact R-hat/ESS definitions and parameter/function set, divergences and sampler limits, recovery/SBC designs, numerical bias/coverage tolerances, PPC targets and pass/fail rules. Diagnose rotation-invariant functions or aligned factors where raw latent coordinates are not identified. Specify allowed remedial runs and the budget/number of attempts. Do not make success depend on θ's sign, H1 significance, or agreement with a preferred baseline result.
- **Predictive checks and comparisons:** held-out units/folds, predictive horizon, treatment of observed versus simulated side histories, and quantities compared across models. A posterior narrower than its prior is an information diagnostic, not proof of identification or calibrated recovery; fix how it will be reported without turning it into a favorable-result gate.
- **Outcome-dependent sampling:** the sketch's stay-out subsampling requires a specified likelihood correction, not an unspecified generic offset. With inclusion probabilities s_A, s_B, s_out, a flat categorical sampled likelihood shifts each join logit relative to out by `log(s_join/s_out)`. Define those probabilities, retain the original risk histories/exposures and distinguish this sampling correction from opportunity q. Explain the sequential/competing-event version before applying it there.

The sketch's selection claim at 24 is too strong. Sharing b models one specified source of dependence between target choice and joining; it does not automatically correct selection into conflict or fatal disputes. Stage 1 remains conditional target choice, not an initiation/fatality selection equation. Its exclusion restrictions at 27 need explicit substantive justification and sensitivity checks, not just assignment to one equation. Similar joint and separate estimates at 34 would establish stability under that fitted latent specification, not absence of unmeasured-trait bias. Bring these statements into line with the revised memo's appropriately cautious identification table.

## Decision and minimum next revision

**Hold the analysis-ready build until B1–B5 are fixed in writing.** Resolve the originator outcome/cohort first; it affects the monthly categories, one-shot labels, q training and every standardization. Then fix the month-zero/empty-side rules, complete the fold map and severity definitions, and enumerate the legal sweep cells. Synthetic fixtures should cover an early A promotion, an early B entrant, no pre-month roster, an empty later side, a system exit, administrative censoring, a 1914 onset and a multilateral dispute with mixed jointly-fatal pairs.

Schema/provenance planning can proceed, but coding defaults into those unresolved branches would make the implementation choose the estimand. AME work can remain deferred as planned; its model and acceptance gate must be finalized before its own results, with the selection language corrected now. No change to the approved five broad design choices is required.
