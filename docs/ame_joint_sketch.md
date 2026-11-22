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
- Dependence, not selection correction (GPT stage-2 memo review). Sharing b_j models one specified source of dependence between being targeted and attracting joiners. It does not correct selection into conflict or into fatal disputes: stage 1 is conditional target choice, not an initiation or fatality equation.

## Identification
- Exclusion restrictions (each needs a written substantive justification and a sensitivity that relaxes it). Stage 1 only: initiator–target dyadic terms (distance, contiguity, prior MIDs). Stage 2 only: the joiner's own ties to each side and its power relative to each side.
- θ is identified from targets attacked repeatedly and disputes with several potential joiners. Possibly weak, so report the prior-to-posterior shift.

## Estimates the paper would report
1. Stage-1 network contrast (coalition originators − bilateral) net of b_j and the latent factors.
2. Stage-2 recruitment: λ_A, λ_B, μ.
3. θ_A, θ_B: does a shared unobservable link being targeted to attracting joiners?
4. Separate vs joint: similar estimates show stability under this fitted latent specification, not the absence of unmeasured-trait bias.

## Estimation
- Bayesian, Stan. Log-sum-exp clogit for stage 1, categorical likelihood for stage 2 (four categories incl. the originator competing event; memo §13 B1), non-centred random effects. Λ is diagonal (K = 1–2) with half-normal scales; an LKJ prior applies only to a correlation matrix (e.g. of the latent factors), not to Λ itself. Location, scale, sign and rotation constraints are fixed in the acceptance-gate document.
- Scale. Stage 1 ≈ 160,000 rows. Stage 2 runs to millions of joiner-months, so stay-out rows are subsampled with known inclusion probabilities s_A, s_B, s_orig, s_out; the sampled categorical likelihood shifts each entry logit relative to stay-out by log(s_entry / s_out). This sampling correction is distinct from the opportunity q; exposures and risk histories are kept.
- Validation. Parameter recovery on simulated data, simulation-based calibration, comparison with the separate models, LOO.

## Risks
- Mixing. Earlier AME fits had unconverged latent scales, so start at K = 1 with tight priors.
- θ may be weakly identified. Informative in itself; report honestly.
- Compute. The joint fit is likely a Roar job, not a laptop one.

## Open items (offered, not decided)
- Acceptance gate (to be fixed in writing before any AME result): exact R-hat/ESS definitions and parameter set (rotation-invariant functions where coordinates are not identified), divergences, recovery/SBC design and tolerances, PPC targets, held-out folds, allowed remedial runs. Proposed starting values: R-hat ≤ 1.01, bulk/tail ESS ≥ 400. A posterior narrower than its prior is reported as an information diagnostic, not a gate.
- Unified three-part version adding an initiation-propensity stage (τ) and targeting jointness (ρ).
- Contested-support measure as the relational term.
- Stage-2 specification fixes from the GPT stage-1 review (M1–M5) carry over to the joint model's stage 2.
