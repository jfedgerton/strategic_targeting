# Joint stage-1/stage-2 AME model — FUTURE STEP (sketch for discussion, recorded 2026-11-21)

Not to be built yet. Build order: separate stages first (stage-1 init multiverse, then stage-2 per `docs/stage2_design_memo.md`) → stage-2 AME alone → stage-1 AME with b_j → joint model. Each step gets an independent agent audit and a GPT review before any run.

## Structure
Two decisions share unobserved target traits and latent positions.

**Stage 1, initiation.** For event e, initiator i, year t and alternative j in the choice set:

U¹(e,j) = β_g′ f(dM_ij, dR_ij) + γ′z_ij + b_j(t) + u_i′Λ v_j + ε

- f is the quadratic surface in the power and network gaps; g is the initiator group (bilateral or coalition originator).
- b_j(t) is a target random effect (attractiveness as a target for unmeasured reasons, e.g. rogue status or grievances), AR(1) over time.
- u_i′Λ v_j is the AME multiplicative term (rank K = 1–2), capturing unmeasured similarity and third-order dependence between initiator and target.
- ε is Gumbel noise, so the likelihood is a clogit per event.

**Stage 2, joining.** For potential joiner k in dispute d (target j*) in month m, options {A, B, out}, "out" the reference:

U²(k,A) = λ_A·tie(k→A) + μ·R_j* + δ′x_kA + θ_A·b_j*(t) + u_k′Λ v̄_A + ε
U²(k,B) = λ_B·tie(k→B) + δ′x_kB + θ_B·b_j*(t) + u_k′Λ v̄_B + ε

- v̄_A, v̄_B are the mean latent positions of each side's current members.
- The stages are linked by the shared target effect b_j*. θ_A > 0 and θ_B < 0 would mean the unmeasured trait that makes j* a target also pulls third parties onto the attacking side and keeps defenders away (the rogue-state confounder, estimated rather than assumed away).
- Selection. Stage 2 exists only for disputes stage 1 produced; sharing b_j handles selection on unobservables, Heckman-style.

## Identification
- Exclusion restrictions. Stage 1 only: initiator–target dyadic terms (distance, contiguity, prior MIDs). Stage 2 only: the joiner's own ties to each side and its power relative to each side.
- θ is identified from targets attacked repeatedly and disputes with several potential joiners. Possibly weak, so report the prior-to-posterior shift.

## Estimates the paper would report
1. Stage-1 network contrast (coalition originators − bilateral) net of b_j and the latent factors.
2. Stage-2 recruitment: λ_A, λ_B, μ.
3. θ_A, θ_B: does a shared unobservable link being targeted to attracting joiners?
4. Separate vs joint: if estimates barely move, the mechanism is not an artifact of unmeasured target traits.

## Estimation
- Bayesian, Stan. Log-sum-exp clogit for stage 1, categorical likelihood for stage 2, non-centred random effects, LKJ prior on Λ.
- Scale. Stage 1 ≈ 160,000 rows. Stage 2 runs to millions of joiner-months, so subsample "stay out" rows (mostly unlikely joiners) with a known log-offset correction (exact for logit; reuses the likely/unlikely tiers).
- Validation. Parameter recovery on simulated data, simulation-based calibration, comparison with the separate models, LOO.

## Risks
- Mixing. Earlier AME fits had unconverged latent scales, so start at K = 1 with tight priors.
- θ may be weakly identified. Informative in itself; report honestly.
- Compute. The joint fit is likely a Roar job, not a laptop one.

## Open items (offered, not decided)
- Pre-set convergence/acceptance criteria (proposed: R-hat ≤ 1.01, ESS ≥ 400, recovery, θ posterior narrower than prior, PPCs).
- Unified three-part version adding an initiation-propensity stage (τ) and targeting jointness (ρ).
- Contested-support measure as the relational term.
- Stage-2 specification fixes from the GPT stage-1 review (M1–M5) carry over to the joint model's stage 2.
