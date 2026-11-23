# Joint two-stage AME model — design memo (rev 2, 2026-11-23)

Status: design only, nothing implemented. Per the PI's standing rule (2026-11-17), the AME design gets an independent agent audit and a GPT review before any code, and its acceptance gate (section 7) is fixed before any AME result is seen. Rev 2 applies the independent audit of rev 1 (P1-1 to P1-3, P2-1 to P2-5, P3-1 to P3-3; each change is tagged). It mirrors the two-stage FE model in R/lib/s2_cell.R.

The model estimates **dependence** between target choice and joining through a shared target trait. It is **not a selection correction**: stage 1 is conditional target choice, not an initiation or fatality equation, and both stages condition on fatal disputes. [P1-1d]

## 1. What the AME adds over the two-stage FE model
| | Two-stage FE (current) | Joint AME (this memo) |
|---|---|---|
| Target trait b_j | fixed effects from stage 1 (repeat targets only), plug-in mean bbar in stage 2 | random effect estimated from both stages; θ benchmarked against a cut posterior that uses stage-1 information only |
| Between vs within targets | FE stage 1 = within; no-FE fit = between | one model with the relational gap decomposed into between and within parts, full surface in both (section 3) |
| Type-specific unobserved targeting | not modelled | b_j^g = b_j + g·δ_j, σ_δ reported |
| Dispute-level joining shocks | not modelled | stage-2 dispute random effect ε_d |
| Unmeasured dyadic similarity | none | rank-1 latent factors as a sensitivity (K = 0 is main) |

## 2. Data (identical to the FE cells; no new build)
- Stage 1: the init dump of the window/sample (fatal MIDs, window 1 main): events, initiator type g ∈ {bilateral, coalition}, choice set of all system members, chosen flag, event weight w, the stage-1 controls Z.
- Stage 2: the one-shot reference cohort of the stage-2 build (≈ 63,000 rows on fatal MIDs; ≈ 270,000 on all MIDs), four categories {A, B, O, out} with the availability mask, the S2_REF_DESIGN terms.
- Relational measure: the book measure, positive part (BOOKc12_ptaany|igo+def+off+pta+dca|prpos). Signed version as one cell.
- **Positions standardized within year** before any gap or mean is formed (the book measure trends strongly over time; an era-mixed mean is not comparable across years). [P1-3a]

## 3. Model
**Stage 1** (conditional logit within event e):

U1(e, j) = f_g(dM_ij, dB_ij, dW_j) + γ' z_ij + b_j^g(era) [+ λ1 · u_i v_j]

- R̃ = within-year standardized log PageRank. R̄_j(t) = **pre-event expanding mean** of R̃_j over the years j was in the system up to t − 1 (main); the full-period mean is a sensitivity. [P1-3b]
- dB_ij = R̃_i,t−1 − R̄_j(t): the **between** gap (initiator vs the target's long-run position).
- dW_j = R̄_j(t) − R̃_j,t−1: the target's **within** deviation, so dB + dW = the total gap.
- f_g: the full quadratic surface by type in (dM, dB, dW): dM, dM², dB, dB², dW, dW², dB·dW, dM·dB, dM·dW. Both parts get the full surface, so nothing is pushed into b_j. [P1-2c]
- **Estimand** (the PI's question): the **type contrast of the between surface**, coalition − bilateral, at dB = ±1 with dW = 0 and dM = 0. Type-specific between contrasts are reported too. [P1-2a]
- b_j^g = b_j + 1[g = coalition] · δ_j: target trait common to both types plus a coalition-specific deviation, so the type contrast does not absorb differential unobserved targeting. b_j ~ Student-t(ν = 4, 0, σ_b), δ_j ~ N(0, σ_δ²). [P1-2b, P2-2]
- era: rev 2 main uses era-specific b_j over the six period folds already defined (b_j,era = b_j + η_j,era, η ~ N(0, σ_η²)); fully time-invariant b_j and AR(1) b_j(t) are the two neighbours. [P1-3c]
- R̄_j and b_j are allowed to correlate only through the Mundlak terms; the between coefficients therefore include any unobserved attractiveness that tracks a target's long-run position (the classic Mundlak reading); the type contrast removes the part common to both types. Stated in the paper. [P1-2]

**Stage 2** (categorical over available categories; "out" = reference):

U2(k, c, d) = α_c + β_c' x_kc + θ_c · b̄_d + ε_d,c [+ λ2 · u_k v̄_opp(c)], c ∈ {A, B, O}

- x_kc: the S2_REF_DESIGN terms, unchanged.
- b̄_d: mean over the dispute's onset target set of b_j,era (common part), computed **from the draws inside Stan**. [P3-3]
- ε_d,c ~ N(0, σ_ε,c²): dispute-level joining shock per category, so θ is identified only through b_j values that stage 1 informs, not as a dispute frailty. [P1-1a]
- Stage 2 keeps total-gap surfaces toward the opposing side (no between/within split there); stated as a design choice. [P3-1]

**Latent factors (sensitivity only).** K = 1, u, v shared across stages with stage-specific scales λ1, λ2; u sum-to-zero over initiators; v orthogonal to R̄ and M̄ (v is the residual of a free vector on those two, inside Stan). K = 0 is the main specification for all reported contrasts. [P2-1]

## 4. Identification
- b_j: mean zero; θ identified from repeat targets and multi-candidate disputes, net of ε_d.
- **Cut-posterior benchmark for θ**: draw b_j from the stage-1-only fit and feed it to stage 2 (no feedback from stage 2 to b_j); report θ from the joint fit and the cut fit side by side. [P1-1b]
- Exclusion choices (stage-1 dyadic controls only in stage 1; recruitment terms only in stage 2) are specification choices, not identification requirements (stage 1 identifies b_j). [P3-2]
- K = 1: report only |u·v| summaries (variance share of the latent term, Frobenius norm of λ1·u v'); no sign anchor gates anything. [P2-1e]

## 5. Priors
β, γ, α, θ ~ N(0, 2.5²) on standardized covariates; σ_b, σ_δ, σ_η, σ_ε ~ half-N(0, 1); λ ~ half-N(0, 1); u, v ~ N(0, 1), non-centred. Prior sensitivity: N(0, 1) and N(0, 5²) on θ and the surface terms.

## 6. Weights and estimation
- **Weights.** Stage-1 event weights (1/targets in the same dispute-year) make the weighted posterior a pseudo-posterior. Main: fit with weights and apply the Williams–Savitsky curvature adjustment (posterior draws rescaled by the sandwich using dispute-clustered score contributions). Sensitivity: unweighted fit. Calibration is checked with a clustered, multi-target recovery DGP (section 7). [P2-3]
- **Software.** cmdstanr (threading with reduce_sum), precomputed design matrices, stage 1 segmented by event; 4 chains × 1,000 + 1,000; seed 123. Expected 0.5–2 h per chain on fatal MIDs at K = 0; K = 1 may be up to ~10× slower if tree depth saturates. All-MID cells on the Mac are feasible; Roar as fallback. [P2-5]
- **Order.** (1) stage 2 alone with b̄ fixed at the FE plug-in (checks the Stan stage-2 code against fit_mnl); (2) stage 1 alone with b_j^g and the between/within surface (checks against the no-FE and FE clogits); (3) cut posterior; (4) joint K = 0; (5) joint K = 1.

## 7. Acceptance gate (fixed before any AME result)
- Convergence: R-hat ≤ 1.01, bulk/tail ESS ≥ 400 for all β, γ, α, θ, σ's and λ; zero divergences (≤ 0.1% only with adapt_delta 0.99 and a reported re-run).
- **Calibration (SBC)**: 100 simulated datasets on the real design at K = 0, generated from a fully specified DGP (b_j, δ_j, η, ε drawn from their priors at the posterior-mean scales of the stage-1-only fit; surfaces and θ at the FE estimates); rank uniformity for θ_c, σ_b, σ_δ and the between type contrast (χ² test on 20 bins, report binomial intervals for 95% coverage). [P2-4a,b]
- **Misspecification checks** (20 reps each, reported, θ false-positive rate gated at ≤ 10% for 95% intervals): (i) θ = 0 with a true dispute frailty; (ii) b_j correlated with R̄_j; (iii) type-specific b_j; (iv) drifting b_j(t). [P1-1c, P2-4c]
- Posterior predictive: times-targeted per state (heavy left tail), chosen-target rank (stage 1); category shares by tier and dispute size (stage 2). [P2-2]
- Agreement checks (diagnostics, not gates): cut vs joint θ; stage-1 between contrasts vs the no-FE clogit; K = 0 vs two-stage FE bootstrap.
- Success never depends on θ's sign or significance or on agreement with the FE model. Allowed remedial runs: longer chains, adapt_delta 0.95/0.99, tighter scale priors; each reported.

## 8. Cells
| Cell | Sample | Window | Measure | Latent |
|---|---|---|---|---|
| 1 main | fatal MIDs | 1 | book positive | K = 0 |
| 2 | fatal MIDs | 1 | book positive | cut posterior (θ benchmark) |
| 3 | fatal MIDs | 1 | book positive | K = 1 |
| 4 | all MIDs | 1 | book positive | K = 0 |
| 5 | fatal MIDs | 0 | book positive | K = 0 |
| 6 | fatal MIDs | 3 | book positive | K = 0 |
| 7 | fatal MIDs | 1 | book signed | K = 0 |
| 8 | fatal MIDs | 1 | book positive | K = 0, unweighted |
| 9 | fatal MIDs | 1 | book positive | K = 0, full-period R̄ |

## 9. Reported quantities
Between type contrast and type-specific between contrasts; within contrasts; material contrasts; σ_b, σ_δ, σ_η; stage-2 surface contrasts and recruitment terms; θ_A, θ_B, θ_O (joint and cut); σ_ε; latent variance share (cell 3); predicted-probability contrasts with 95% posterior intervals.

## 10. Risks
- θ may shrink once ε_d absorbs dispute-level variation; that is the honest answer, and the cut posterior shows how much information stage 1 carries.
- Era-specific b_j adds parameters; if σ_η collapses, report the time-invariant neighbour.
- The weighted pseudo-posterior's calibration rests on the curvature adjustment; the unweighted cell shows sensitivity.

## 11. Decisions for the PI
1. The full between/within surface (rev 2) as the stage-1 specification.
2. K = 0 main with K = 1 as a sensitivity (audit recommendation) versus K = 1 main.
3. Era-specific b_j as main (audit recommendation) versus time-invariant.
4. cmdstanr on the Mac (install needed) versus rstan.
