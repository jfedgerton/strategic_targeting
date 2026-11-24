# Book pipeline plan (draft 2026-10-11)

Target: Cambridge University Press. Implication for the pipeline: the book must add analyses beyond JOP 2024 and Tremors (see section 7), and one shared measurement spine is what makes it a book rather than three papers.

## 1. How the pieces fit (the cycle)

| Step in the cycle | Source | Core measure | Core test |
|---|---|---|---|
| A. Strong states form sub-orders | JOP 2024 | Latent multiplex network -> temporal Louvain communities (DLCD) | Historical bloc validation; centrality face validity |
| B. Orders homogenize members (political, economic, social) | JOP 2024 H1 | Exposure to stronger partners with attribute X | RE logit; network and CINC permutations; VAR modularity -> regime assortativity |
| C. Orders are internally peaceful, externally conflictual | JOP 2024 H2 | Out-group dyad indicator | FE logit/LPM, TERGM, LFM on fatal MIDs |
| D. Orders polarize / misalign (middle stream, OPEN) | disorder_conflict (hierarchy misalignment, entropy, ABM) | Network entropy; material vs relational power gap | Dyadic + system VAR (entropy <-> fatal MIDs) |
| E. War remakes the system | Tremors | JS divergence (system), weighted-Jaccard churn (state), partner war exposure | OLS/HAC, VAR/Granger/FEVD, TWFE, local projections, placebos, mediation |
| F. Feedback: system-made ties confound war studies | Tremors | Conflict history t-2..t-6 | Contamination ladder; benchmark lags; 94-study replication sweep |
| Case | Tremors (cut from AJPS version) | Gallipoli 1354 | Event coding H1-H3 |

Shared spine: every step uses the same dyad-year multiplex panel (1816-2014) and the same
tie-construction choices (layers, weighting, aggregation). That is why weighting and layers are settled first.

## 2. Sequence (each stage needs sign-off before the next)

| Stage | Task | Output |
|---|---|---|
| 0 | This map + chapter outline | PLAN.md |
| 1 | Layer inventory: coverage, cadence, missingness, inter-layer correlation | Layer table + correlation heatmap |
| 2 | IGO and alliance diagnostics (see 3) | Diagnostic figures/tables |
| 3 | Weighting-scheme tests on a fixed layer set | Comparison of schemes on validation targets |
| 4 | Fix the canonical network (layers + weighting + aggregation) | Canonical universe cached |
| 5 | Layer x weighting multiverse, per chapter | Long spec table per chapter |

## 3. Stage 1-3 detail (all drawn from what the two papers already did)

Candidate layers (from JOP + Tremors sweeps): ATOP defense, offense, consultation, neutrality; COW IGO;
PTA; DCA v1/v2; bilateral trade; diplomatic ties; BLA; arms transfers (r_networks/arms_tidy.rds exists).

| Check | Prior source |
|---|---|
| Coverage window per layer (trade 1870+, PTA 1960+, DCA 1980-2010) and the "must include IGO or defense" rule | JOP fn 8; Tremors App 1.1 |
| IGO cadence break (COW IGO is 5-yearly before 1965, interpolated) | Tremors 16_robustness_igo_cadence.R |
| IGO density break / near-universal IGOs | Tremors 15_igo_density_break.R |
| IGO membership-share sliding windows (drop IGOs below/above x share) | Tremors App 5.1; disorder_conflict R&R |
| IGO exclusivity weighting (1 - share) | disorder_conflict |
| Inter-layer correlation | JOP App Fig 3 (2_weight_data) |
| Weighting schemes: unweighted; contiguity/distance residualized x Manhattan/Euclidean; JOP's 12 FE formulas; average vs Manhattan | Tremors App Tab 2; JOP scripts 12-13 |
| Edge-weight validation | Tremors 08_validate_edge_weights.R |
| Face validity targets: centrality ranks; historical blocs (Concert, Axis/Allies, Cold War blocs) | JOP Fig 4, Tab 1 |

## 4. Multiverse design to mimic (verified in code)

JOP (12_/13_testing_formula_downweight): 10 layers, all combos with defense/offense/IGO, no DCA v1+v2 together
(671 combos) x 12 weighting formulas (11 feols residualizations + raw) x 2 aggregations (sum/average) -> DLCD
-> feglm fatal_mid ~ out_group | ccode1+ccode2+year+dyad. Approx 16k specs. Reports median coef, share positive, share p<.05 one-sided.

Tremors (07_robustness_multiverse.R, 4,056 lines): layer combos requiring IGO or defense x 5 weighting schemes
x 3 transforms x DVs x IVs x models x SE clustering. One row per spec, then meta-regression of coef, p, sig on
decision dummies + heatmap + significance bars. Checkpointed per (transform, method), furrr parallel.

Book version: build each network universe once (layers x weighting x aggregation), cache it, and let every
chapter's model read the cache. Output schema per row: spec_id, chapter, hypothesis, layer dummies, weighting,
aggregation, measure choice, IV, DV, transform, estimator, SE, term, estimate, se, p, n, canonical flag.

## 5. Fragilities found in the existing code (verify before porting)

| Issue | Where | Status |
|---|---|---|
| Multiverse rounds residualized ties to integers (round(tie)) and skips range01; canonical script uses range01 + 2 decimals. Adjacency built without weighted=TRUE. The sweep may have run on near-binary networks unlike the main model | JOP 12_testing_formula_downweight | Hypothesis from code reading; not rerun |
| Year assigned as i + 1815 after try(silent) Louvain; a failed year shifts all later years | JOP 12 | Hypothesis from code reading |
| Seeds 9733 / 477 | JOP 2, 12 | Book uses 123 |
| README says KL is canonical; manuscript uses JS | Tremors README | Stale doc |
| Script holds ~2 GB globals per worker | Tremors 07 | Design around cached universes |
| dynamicmultiplex weighted Jaccard bug before 1.2.1 | DynMux | Use >= 1.3.x if used |
| JOP .R files on the Mac are 0-byte placeholders; real files only in Dropbox cloud; Claude project copies also 0 bytes | Data access | Set cooperative_regimes to available offline |

## 6. Draft chapter outline

| Ch | Title (working) | Analyses | Status |
|---|---|---|---|
| 1 | Introduction | Cycle overview; Tilly analogy | Writing |
| 2 | Theory of conflict cycles | Orders -> homogenization -> polarization -> war -> reconfiguration | Writing |
| 3 | Measuring the international system | Layers, weighting, aggregation; Stage 1-3 diagnostics; centrality validity | From JOP + Tremors |
| 4 | The emergence of orders | Temporal community detection; historical validation; order counts, modularity over time; alt algorithms and interlayer lags | From JOP |
| 5 | Orders make states alike | Institutional diffusion (5 DVs + V-Dem alternates); permutation tests; VAR | From JOP |
| 6 | Peace within, war between | Out-group and fatal MIDs; FE, TERGM, LFM; joint democracy null | From JOP |
| 7 | Disorder and polarization | Hierarchy misalignment / entropy and conflict | OPEN: include disorder_conflict? |
| 8 | War remakes the system | JS divergence, churn, partner exposure, LPs, placebos, mediation | From Tremors |
| 9 | The cycle closes | Contamination ladder; replication sweep | From Tremors |
| 10 | Gallipoli 1354 | Case and event coding | From Tremors (cut from AJPS) |
| 11 | Conclusion | | |

## 7. Unimplemented suggestions (need approval; not in either paper)
- Link Ch 4 and Ch 8 directly: does war predict order births, deaths, merges, and splits (community lineage events)?
- Order-level polarization over time (between-order distance) as the bridge for Ch 7.
- Same permutation + systematic measurement-error stress test applied to every chapter's canonical model.

## 8. Decisions 2026-10-11
- Dyad-year is the primary unit of analysis. System/order results are kept as descriptive anchors (bloc validation, modularity, divergence).
- Every dyadic model that uses network-derived terms gets (a) a within-year permutation of order labels and (b) an LFM/AME check as in JOP.
- Layers with partial coverage (trade 1870+, PTA 1960+, DCA 1980-2010) are allowed in, conditioned with year and dyad FE/RE.
- Checks to run in Stage 1-3 because of that choice:
  1. FE absorb level shifts but not which dyads gain ties when a layer enters. Compare full-window FE estimates with common-support years (1980-2010).
  2. Coverage breaks also enter upstream, inside the composite tie and the community detection, where outcome-model FE cannot reach.
  3. Dyad FE in a rare-event logit drops dyads that never fight (JOP FE logit n = 30k vs 688k). Report LPM/RE alongside.
  4. Residualizing a single layer on contiguity + sender + receiver + year FE is in the span of dyad + year FE. In a linear outcome model, the weighting schemes should then give nearly identical coefficients (FWL). Weighting choices bite only through nonlinear steps (clipping negatives, within-year rescaling, Euclidean, community detection). Stage 3 tests weighting where it can matter.

## 9. Stage 3 validation set (decided 2026-10-11)
Truth = Braumoeller/Goodhart order coding (Only the Dead replication), as already coded in
dynamic_multiplex/review_bundle/code/replication/exploratory/paper_scripts/27_ and 29_braumoeller_*.R.
Reuse those coding tables and metrics (ARI 3-bloc L/C/Other, ARI 2-bloc L vs C, bloc purity). Keep the
Colombia(100)/Bulgaria(355) fix (FIXBUG = TRUE).

| Era | Partition the coding supports | Use |
|---|---|---|
| Cold War 1945-91 | Liberal / Communist / Other | Primary test for signed vs unsigned, lambda |
| League 1920-41 | League members vs others | Hold-out check |
| Concert 1816-52 | Europe / Latin America / other | Coarse hold-out |
| 1855-1914 | Europe vs non-Europe | Coarse hold-out |
| Post-Cold War 1992+ | Liberal order vs others | Hold-out check |

Built-in divided-state pairs in the Cold War coding: FRG 260 (L) vs GDR 265 (C); ROK 732 vs PRK 731; ROC 713 vs PRC 710.
Allied-rival pair kept together: Greece 350 and Turkey 640 (both L). South Vietnam (817) is not coded.
Rule: choose lambda on Cold War ARI, then report hold-out eras with lambda fixed. Report the full lambda grid in the appendix.

## 10. Main measure decision (2026-10-11)

MAIN = hybrid score, code/02_network_universes/01_main_measure.R
- Layers: defensive, offensive, DCA (binary; undirected, tie if either side pledges), shared IGO, trade (within-year 0-1).
- Each layer only in years it varies (JOP rule).
- Expected tie: logit (binary) / Poisson (IGO, trade), y ~ log(1 + capdist) | sender + receiver + year FE.
- Score: binary present = 1, absent = -p_hat. IGO/trade = y if y >= y_hat, else y - y_hat.
- Composite = sum of layer scores / layers covered that year; undirected average of (i,j),(j,i).
- Versions: tie_signed (raw), tie_norm (/ max |tie| within year, [-1,1], for pooled dyadic regressions), tie_pos (clipped, / year max).
- Communities: Leiden, resolution 1, per year, positive (tie_pos) and signed (Traag-Bruggeman, lambda = 1). No temporal linking yet.
- Why: absent ties never score positive; present ties are not penalized for being expected.
  Scores (results/main/variants_pair_checks.csv): divided pairs together (signed) 13%, partners together 97%,
  partners mean tie percentile 91, Cold War L vs C ARI 0.78 (pos) / 0.83 (signed). Best Braumoeller recovery tested.
- Known weak spots: UKG-FRN and AUL-NEW sit near the 55th percentile (IGO/trade fall short of a high expectation);
  PRC-ROC together in 29% of years (signed).

## 11. Multiverse arms for the measure (alternatives tested 2026-10-11; all in 04_measure_variants.R)

| Arm | Expected-tie model | Geography | Score | Status |
|---|---|---|---|---|
| hyb_dist (MAIN) | logit / Poisson | log distance | hybrid | built |
| glm_dist | logit / Poisson | log distance | surprise (observed - expected) | built |
| glm_dist_norm | logit / Poisson | log distance | surprise, / max abs within year | built; identical ranks and communities to glm_dist |
| glm_contig | logit / Poisson | contiguity | surprise | built |
| glm_contig_dist | logit / Poisson | contiguity + log distance | surprise | built |
| glm_ns3 / ns4 / ns5 | logit / Poisson | natural spline of log distance, 3/4/5 df | surprise | built; no gain over log distance |
| lpm_contig (JOP) | linear FE | contiguity | surprise; absent ties can score positive | built |
| lpm_dist | linear FE | log distance | surprise | built |
| hyb_dist_igo_surprise | logit / Poisson | log distance | hybrid for binary layers, surprise for IGO/trade | NOT built (may lift UKG-FRN, AUL-NEW relative to other allied pairs; untested) |
| Community version | | | positive vs signed Leiden; signed lambda grid | positive and signed built at lambda = 1 |

Surprise vs hybrid in one line: surprise = cooperation beyond what geography and activity predict (expected allies score ~0);
hybrid = how tied the dyad is, with conspicuous absences negative.

## 12. First-pass replications (2026-10-11; FE models with lagged DV; hybrid main, surprise checked)
Scripts: code/03_chapters/01_war_remakes_ties.R, 02_orders_conflict.R, 03_orders_spread_institutions.R (MEASURE=hybrid|surprise).
Communities for dyad models: multislice Leiden (07_multislice_leiden.py), gamma = 1, omega = 1, positive and signed.

| Claim | Source | Hybrid | Surprise |
|---|---|---|---|
| MID onsets -> state churn | Tremors H2 | replicates (0.99 x100, p=.01) | replicates (1.24, p=.004) |
| Fatal onsets -> churn | Tremors | p ~ .09 | p ~ .08 |
| Partner war exposure -> churn | Tremors H3 | positive, n.s. (p=.30) | positive, n.s. |
| Out-group -> fatal MID / MID | JOP H2 | replicates in 7 of 8 (pos-community FE logit fatal p=.07) | replicates in 8 of 8 |
| Joint democracy (control only) | JOP borderline | negative and significant in FE logit | same |
| Stronger-partner exposure -> institutions | JOP H1 | with lagged DV: autocracy, restricted participation (logit); all 5 DVs significant without lagged DV | same pattern |

Validation: churn/PWE code run on Tremors' own main_weight reproduces Tremors (churn mean .074, PWE mean .31, PWE coef 8.78 vs 8.86).
Onset coefficient 0.69 vs 0.48: Tremors counts MID onsets on the directed panel (each dyad twice) - from code reading, not rerun.
Covariates: Polity5 and V-Dem v12 from the JOP source files; Polity country codes recoded to COW (305->300 pre-1919, 364->365,
342/347->345, 818->816), recovering 277 state-years that peacesciencer also drops.

## 13. Layer-set tests (per-year Leiden; results/main/layer_set_comparison.csv)
- PTA (DESTA or Mansfield-Pevehouse) in place of trade, or added: worse on every criterion (divided pairs together 21-53%, Cold War ARI falls).
- log(1 + x) IGO and trade: worse (partners together 70-76%, Cold War ARI 0.44-0.67).
- Dropping trade: as good or slightly better (signed: divided together 4-9%, partners 96-99%), Cold War ARI slightly lower. Simplest option.
- Open decision: keep trade, or 4-layer main measure (defensive, offensive, DCA, IGO).

Decision (Jared, 2026-10-11): joint democracy stays a control. It was borderline in JOP and is not a book claim, so no follow-up analysis.

Decision (Jared, 2026-10-11): main measure = surprise (glm_dist, y - p-hat; logit/Poisson with log distance and
sender, receiver, year FE) on all 5 layers, trade kept (security, institutional, economic dimensions).
Hybrid (hyb_dist) becomes a multiverse arm. Implemented 2026-10-11: main_measure.rds = surprise (matches variant glm_dist, gap 6e-16);
communities_multislice.csv and communities_year.csv = surprise (hybrid copies saved as *_hybrid); chapter scripts default MEASURE=surprise.

## 14. Open items at pause (2026-10-11)
- Surprise is now MAIN in code (done).
- Temporal leakage. Full-panel multislice labels at t use ties after t, which biases out_group in the conflict models.
  Measure FE (sender, receiver) pooled over 1816-2014 also borrow from the future.
  Proposed: (1) conflict models on per-year Leiden as a quick check; (2) rolling-window multislice t-k..t, keep slice t,
  k = 10 main, 5/20 as arms; (3) sender/receiver FE by state x decade or rolling window. Not started.
- Excluded-belligerent churn added to 01_war_remakes_ties.R (x_* models). MID onset holds (surprise 1.15, p=.008).
- Tremors bugs written up in handoff_tremors_mid_double_count.md (MID count doubling; duplicate partners in excl churn).

## 15. Community validation criteria (Jared, 2026-10-11; replaces ARI as primary)
High-level orders only (e.g., Axis vs Allies, Cold War blocs, LIO); sub-order clusters are expected and acceptable.
Ranking, lexicographic:
1. Separation of enemies: rival-order pairs in different communities (e.g., GMY/ITA/JPN apart from USA/UKG/RUS). Most important.
2. Purity: each community drawn mostly from one order.
3. Size: prefer fewer, larger clusters per order (effective number of clusters; same-order pairs together). Tie-break only.
A split ally (UKG apart from USA) costs little; a merged enemy pair (Axis state with an Allied state) costs a lot.
Lever under consideration: Leiden resolution 0.5 / 0.75 / 1.0. Open: LIO coding rule and period; keep pre-1939 Braumoeller as secondary?

## 16. Signed churn (Jared, 2026-10-11) - code/03_chapters/01b_signed_churn.R
Weighted Jaccard over positive and negative parts of raw tie_signed (surprise), layer-consistent (layers scored in both t-1 and t),
cutoffs tau = 0 / .01 / .03 / .07, decomposed into positive turnover, negative turnover, sign-flip share.
- Layer-consistency removes the DCA spikes (2011 mean 0.56 -> 0.06; 1980 0.28 -> 0.16).
- MID onset -> signed churn: 0.65 (p=.06) at tau 0; positive-part turnover 1.23 (p=.008); negative part null; flips 0.37 (p=.06).
- Results move with tau. Main tau not chosen.
- Found: layer_igo and layer_trade are min-max scaled within year, so their scores rescale whenever the year max changes
  (1941: 87% of IGO dyads change by the same ratio). Affects any churn on raw scores. Candidate fix: raw counts or one global scale. Not done.
- Jared not convinced by alternatives to churn (balance, in-order share); not pursued.
Decision (Jared, 2026-10-11): keep the original churn (Tremors definition, tie_pos) as MAIN. Signed churn (01b) stays as a multiverse/robustness arm only.
Known caveat for reviewers: original churn spikes when layers enter/leave (1980, 2011) and inherits year-max rescaling.

## 17. Layer sets and DynMux (2026-10-11)
- main_measure.rds now = surprise on defensive, offensive (true zeros in no-pact years), IGO, trade; DCA dropped.
  Its multislice run was killed by a container restart; communities_multislice.csv is still the 5-layer version.
- TRIAL allyigo (LAYERSET=allyigo, MEASURE=allyigo): one alliance layer (def OR off) + IGO, surprise. Not adopted yet.
  Churn layer jumps gone (2011 0.06, 1922 0.21, 1870 0.06; median 0.15). 1882 = 0.60, unexplained.
- DynMux 1.4.0 (repo commit 60295de) two-stage, filtered tracking, repair_chains: code/02_network_universes/08_dynmux_tracking.py.
  Stage 2 merges rival same-year communities (allyigo: FRG-GDR together 100% pos; Cold War ARI2 0.05), so meta ids
  are for linking across years only. Within-year models use stage-1 per-year Leiden (no leakage).
Decision (Jared, 2026-10-11): MAIN churn and diffusion models drop the lagged DV (Tremors / JOP specs); lagged-DV models kept as robustness.
Conflict models keep the lagged MID. Results: main_* rows in ch_war_remakes_ties_<M>.csv; *_nolag_* files for diffusion.

## 18. War -> ties: monadic now, dyadic next (Jared, 2026-10-11)
- This pass: keep monadic churn (positive ties, no lagged DV; signed churn and lagged DV as robustness).
- Measure trials on alliance + IGO: surprise vs hybrid (results/main/*allyigo*, *allyigo_hyb*). Hybrid: PWE 4.9 (p=.17),
  MID onset weaker (1.00, p=.054), divided pairs separate less well. Surprise vs raw (Tremors) weights for churn still open.
- NEXT PASS (planned, not built): dyad-level tie change. DV = change in surprise tie t-1 -> t+k (or formation/dissolution);
  war roles in t = same side / opposing sides / one belligerent / common enemy; pair + sender + receiver + year FE;
  event-study leads and lags. Needs MID side data merged (check peacesciencer / GML dyadic sides).

## 19. Timing and SEs (Jared, 2026-10-11)
- All three chapters use t-1 predictors as MAIN: MID/fatal onsets at t-1 (churn), communities at t-1 (conflict),
  exposure and neighbor share at t-1 (diffusion; EXPO_TIMING=t1 default in 03). Same-year versions kept as robustness.
- Diffusion: FE (2-way clustered), RE (03b, model SE), and clustered RE (03c): state-cluster bootstrap of the RE
  models (B=50, running) and Mundlak pooled LPM/logit with state+year clustered SE.
