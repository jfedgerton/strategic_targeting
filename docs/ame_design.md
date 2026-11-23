# Joint two-stage AME model — design memo (rev 1, 2026-11-23)

Status: design only, nothing implemented. Per the PI's standing rule (2026-11-17), the AME design gets an independent agent audit and a GPT review before any code, and its acceptance gate (section 7) is fixed before any AME result is seen. This memo supersedes the open parts of `docs/ame_joint_sketch.md` and is written to mirror the two-stage FE model that is now the working main model (R/lib/s2_cell.R).

## 1. What the AME adds over the two-stage FE model
| | Two-stage FE (current) | Joint AME (this memo) |
|---|---|---|
| Target attractiveness b_j | fixed effects from stage 1 on chosen targets, passed to stage 2 as a plug-in mean (bbar) | random effect estimated jointly from both stages; uncertainty propagates without a two-step bootstrap |
| Between vs within targets | FE stage 1 is within-target only; the between estimates come from a separate no-FE fit | Mundlak split: between and within relational effects estimated in one model (section 3) |
| Unmeasured dyadic similarity | none | rank-K latent factors u_i' Λ v_j (stage 1) and u_k' Λ v̄_side (stage 2) |
| Link between stages | θ_A, θ_B on the plug-in bbar | θ_A, θ_B, θ_O on the latent b_j itself |

## 2. Data (identical to the FE cells; no new build)
- Stage 1: the init dump of the window/sample (fatal MIDs, window 1 main): events e with initiator i, type g ∈ {bilateral, coalition}, choice set of all system members j, chosen flag, event weight w (1/targets in the same dispute-year, as now), the quadratic surface terms by type and the stage-1 controls Z (contiguity, log distance, target democracy and interaction, def/off alliance, target major and interaction, dyadic MID t−1, log prior MIDs).
- Stage 2: the one-shot reference cohort of the stage-2 build (≈ 63,000 candidate-dispute rows on fatal MIDs; ≈ 270,000 on all MIDs), four categories {A, B, O, out} with the availability mask (A needs H > W; O only for likely-tier candidates), the S2_REF_DESIGN terms (surface toward the opposing side, mirrored controls, tie_sum, ally_side, share, side sizes, likely tier).
- Relational measure: the book measure, positive part (BOOKc12_ptaany|igo+def+off+pta+dca|prpos). Signed version as one alternative cell.
- No subsampling is needed: the one-shot stage-2 design keeps the row count small enough for Stan (the monthly design would not).

## 3. Model
**Stage 1** (conditional logit within event e):

U1(e, j) = f_g(dM_ij, dR^B_ij) + κ_g · wR_j + γ' z_ij + b_j + u_i' Λ v_j

- f_g: the quadratic surface by initiator type, as in the FE model, but in the BETWEEN relational gap dR^B_ij = (R_i,t−1 − R̄_j) / s_t, where R̄_j is target j's long-run mean position (all years j is in the system).
- wR_j = (R̄_j − R_j,t−1) / s_t: the target's within deviation from its own mean, by type (κ_bil, κ_coal). This is the Mundlak split: f_g carries the between comparison (the PI's question) and κ_g the within comparison.
- The material gap stays as now (no split; within-target material variation is small).
- b_j ~ N(0, σ_b²): target random effect, time-invariant in rev 1 (AR(1) b_j(t) as a later extension).
- u_i, v_j ∈ R^K, K = 1 main, K = 2 alternative; Λ = diag(λ_1..λ_K) ≥ 0.

**Stage 2** (categorical over the available categories; "out" = reference):

U2(k, c) = α_c + β_c' x_kc + θ_c · b̄_d + u_k' Λ v̄_{opp(c)}, c ∈ {A, B, O}

- x_kc: the S2_REF_DESIGN terms for category c (unchanged from the FE model).
- b̄_d: mean of b_j over the dispute's onset target set (same aggregation as bbar, now latent). θ_A, θ_B, θ_O are the reported links.
- v̄_{opp(c)}: mean latent position of the side the joiner would fight (B's members for A and O; A's members for B), so the latent term means "similarity to the opponent", mirroring stage 1.

**Shared parameters.** b_j appears in both stages (the link). u and v for the same state are the same vectors in both stages (one latent space). Λ is shared.

## 4. Identification and constraints
- b_j: mean zero; σ_b half-normal(0, 1). θ identified from repeat targets and from disputes with several candidates (as in the FE model; 576 baseline target pairs).
- Latent factors: K = 1 sign fixed by u_USA > 0 (or the most active initiator); rotation is not an issue at K = 1; at K = 2 report only rotation-invariant quantities (u_i' Λ v_j and its variance share).
- Exclusion restrictions (stated, each relaxed in one sensitivity): stage-1 dyadic controls do not enter stage 2 except through their stage-2 mirrored versions; stage-2 recruitment terms (tie_sum, ally_side, share, side sizes) do not enter stage 1.
- Event weights enter the stage-1 likelihood as weights (pseudo-likelihood), as in the FE model.

## 5. Priors
β, γ, κ, α, θ ~ N(0, 2.5²) on standardized covariates; σ_b ~ half-N(0, 1); λ_k ~ half-N(0, 1); u, v ~ N(0, 1) (non-centred). Prior sensitivity: N(0, 1) and N(0, 5²) on θ and the surface terms.

## 6. Estimation
- Stan (rstan), log-sum-exp clogit for stage 1, masked categorical for stage 2, reduce_sum over events and disputes, 4 chains × 1,000 warmup + 1,000 draws; seed 123.
- Compute: stage 1 ≈ 72,000 rows, stage 2 ≈ 63,000 × 4 on fatal MIDs. Expected hours per K = 1 fit on 4–8 threads (to be timed on a 10% subsample first). All-MID cells are ≈ 4× larger and may need Roar.
- Order: (1) stage 2 alone with b̄ fixed at the FE plug-in (checks the stage-2 Stan code against fit_mnl); (2) stage 1 alone with b_j random and the Mundlak split (checks against the no-FE and FE clogit); (3) joint, K = 0 (b_j only); (4) joint, K = 1; (5) joint, K = 2.

## 7. Acceptance gate (fixed before any AME result)
- Convergence: R-hat ≤ 1.01 and bulk/tail ESS ≥ 400 for all β, γ, κ, θ, σ_b, λ and for u_i'Λv_j summaries; zero divergences after warmup (≤ 0.1% allowed only with adapt_delta 0.99 and a reported re-run).
- Recovery: simulate 20 datasets from the FE estimates with known b_j, θ and K = 1 factors on the real design; 95% intervals cover the truth at nominal rates (≥ 17/20) for θ, the stage-1 between contrasts and σ_b.
- Agreement checks (diagnostics, not gates): K = 0 joint estimates vs the two-stage FE bootstrap; stage-1 between contrasts vs the no-FE clogit.
- Posterior predictive: chosen-target rank distribution (stage 1); category shares by tier and dispute size (stage 2).
- Success never depends on θ's sign or significance or on agreement with the FE model. Allowed remedial runs: longer chains, adapt_delta 0.95/0.99, tighter λ prior; each reported.

## 8. Cells (the "few AME" the PI asked for)
| Cell | Sample | Window | Measure | K |
|---|---|---|---|---|
| 1 main | fatal MIDs | 1 | book positive | 1 |
| 2 | fatal MIDs | 1 | book positive | 0 (b_j only) |
| 3 | fatal MIDs | 1 | book positive | 2 |
| 4 | all MIDs | 1 | book positive | 1 |
| 5 | fatal MIDs | 0 | book positive | 1 |
| 6 | fatal MIDs | 3 | book positive | 1 |
| 7 | fatal MIDs | 1 | book signed | 1 |

## 9. Reported quantities
Stage-1 between contrasts (coalition / bilateral, less central and stronger/weaker, as in the FE tables) and within contrasts (κ); stage-2 surface contrasts and recruitment terms; θ_A, θ_B, θ_O; σ_b; variance share of the latent term; predicted-probability contrasts with 95% posterior intervals.

## 10. Risks
- Mixing of latent scales (seen in earlier AME runs): start K = 0, then K = 1 with the tight λ prior; K = 2 only if K = 1 passes.
- θ may be weakly identified once b_j is random (less information than the FE plug-in suggests); report prior-to-posterior shift.
- The Mundlak R̄_j uses all years j is in the system, including years after the event (a time-invariant summary of the target, standard for Mundlak); an alternative uses only years before the event (sensitivity).

## 11. Decisions for the PI
1. Mundlak split as in section 3 (between surface + within linear term), or keep the FE-style total gap and add R̄_j as a covariate (classic Mundlak; between effect = within + R̄ coefficient, harder to read with a quadratic surface).
2. Shared latent space across stages (one u/v per state) vs stage-specific factors.
3. Time-invariant b_j (rev 1) vs AR(1) b_j(t).
4. Where it runs: Mac (rstan, 18 cores) for fatal-MID cells; Roar for all-MID cells.
