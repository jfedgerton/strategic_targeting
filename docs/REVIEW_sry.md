# Review: book-aligned sry additions

22 November 2026. **NO-GO for the expanded production rerun, for either `TYPES=bcd` or `TYPES=init`, until SRY1–SRY3 below are addressed.** With the default `REUSE_POSITIONS=1`, also establish the linked files' layer provenance (SRY4), or rebuild all arms together. These are measurement-build issues shared by both type schemes.

## Scope and evidence

Reviewed `git diff 230d585 HEAD -- R/02_build_positions.R run_all.sh R/06_make_variants.R README.md`, with the requested alignment note and prior stage-1 recheck as context. Read the relevant unchanged downstream consumers in R03, R04, R07 and R11 to check compatibility. The stage-2 memo revision, historical coding and unchanged estimators were not re-reviewed.

Source: supplied `targeting_multiverse_sry.zip`; HEAD `be2d07545b908be7f712930a9ab3bfec3ff8ea29`; archive SHA-256 `82102f3b7a509d01dcdcf7379ada6716f2ac79920c9c6af5fd2a55d5cc31bc02`. File/line references below are repository-relative paths in that archive, not the older `gpt_audit/full/` snapshot.

**No R scripts, R expressions, project tests or production runner were executed.** Native Bash 3.2 syntax checks passed for `run_all.sh` and `config.sh`. Independent Python enumeration of the documented grid yielded 23 expected position files and 12,294 distinct spec IDs, with no duplicates. Representative IDs were checked against the downstream parsing expressions. The supplied numerical replication results were not independently reproduced. Extraction and inspection took place in a temporary directory; the only workspace file written is this report.

## Findings

### SRY1 — P1: numerical residuals can become full-strength ties

**Where:** `R/02_build_positions.R:156–158, 231–232`.

The existing W2–W4 path explicitly zeros numerical residuals before normalization (`132–133`). The new sry path does not: every strictly positive layer-year maximum absolute residual is used as a divisor, followed by another strictly-positive maximum division after combination.

**Failure scenario:** a nonzero constant-count year is admitted by `yr_has`, and its year FE fits the layer exactly in theory. Small floating-point residuals remain. A three-dyad residual vector `(epsilon, -epsilon, 0)` becomes `(1, -1, 0)` for any positive epsilon, however small. A one-layer combination then contains a weight-1 edge. The orientation check accepts this example because both orientations can have the same tiny error. Irregular numerical errors can produce varying centrality, so the later constant-position flag does not necessarily catch the problem. The same normalization risk applies to an exactly fitted OLS UN layer-year or near-perfect fits. This is a deterministic counterexample to the normalization rule, not a claim that a particular production layer has already exhibited it.

**Minimal fix:** apply an explicit family/response-scale-aware numerical tolerance to raw residuals before layer normalization, and suppress a layer-year whose residual scale is below that tolerance. Also protect the combined signed score against tiny cancellation errors before its positive-part normalization. Keep a genuinely zero residual layer applicable with score zero; do not silently change its denominator status to missing. Verify that tighter solver tolerances do not materially change the resulting networks.

### SRY2 — P2: unconverged expected-tie fits can be published as valid measures

**Where:** `R/02_build_positions.R:151–156, 220, 247–248`.

Every new expected-tie fit uses `warn=FALSE`, and no convergence status is inspected. Symmetry, finite centralities and normalized bounds do not establish convergence.

**Failure scenario:** a sparse logit or Poisson model reaches its iteration limit, or separation involving regressors remains beyond the FE-only removals. The returned fitted values can be finite and symmetric while the fit is unconverged. Their residuals are then normalized and saved as ordinary valid positions. Agreement with a standalone implementation using the same settings does not resolve this failure mode.

**Minimal fix:** fail the arm on an unsuccessful GLM convergence status, non-finite fitted values, or unsuccessful FE convergence; retain relevant warning/status diagnostics per layer. Distinguish expected removal of perfectly predicted FE groups from convergence failures. Increasing iteration limits alone is insufficient without checking the outcome. `fixest` documents `convStatus` for IRWLS and separate FE iteration controls: [feglm reference](https://lrberge.github.io/fixest/reference/feglm.html).

### SRY3 — P2: every omitted observation is treated as a perfectly predicted FE group

**Where:** `R/02_build_positions.R:146, 149, 154–159`.

`ok` checks the outcome and year applicability, but not `lcap` or the required contiguity covariate. After estimation, every row absent from `obs(m)` is filled with its observed outcome and therefore assigned residual zero. Every such omission is logged as “perfectly predicted.”

**Failure scenario:** a dyad has missing distance, or missing contiguity under W7. `fixest` omits both orientations, but this code silently turns them into an applicable, perfectly explained layer rather than rejecting missing predictors. It changes signed-layer averaging and obscures missingness. The unchanged W4 path instead asserts that no observations disappeared; W2/W3 also assert residual-vector length.

**Minimal fix:** assert finite required predictors before fitting. Reconcile omitted rows with documented FE removals/singletons, rejecting any unexplained omission. Assert unique, in-range observation indices, matching lengths and finite fitted values before assignment. If an entire layer is perfectly predicted, handle that explicitly rather than relying on a fit with no remaining observations. The `obs()` API reports estimation-row indices, not reasons for exclusion: [obs reference](https://lrberge.github.io/fixest/reference/obs.html).

### SRY4 — P2: the hybrid grid can mix layer versions despite a valid fingerprint

**Where:** `run_all.sh:66–77, 90–101`; `R/02_build_positions.R:248`.

The new arms are built from this run's `layers_dy.rds`, but reused W1–W4/IGO/diplomacy files are checked only for existence and link destination. Their generating layers are not compared with the current layers. The recorded `layer_rds` pathname is not a content-provenance check.

**Failure scenario:** start a fresh output directory after raw network inputs change, or point the reuse configuration at a same-named older build. The current layers and new arms rebuild, the older arms are linked, and the first fingerprint accepts both sets of bytes. Matching state-year keys downstream does not detect different edge inputs. The global fingerprint identifies the mixture; it does not prove its coherence.

This is an inherited trust boundary for external reused files, made consequential here by building new arms alongside them. It does **not** reopen the previously fixed interrupted-build stamp defect for locally generated arms.

**Minimal fix:** require immutable source-build manifests tying each linked arm to compatible raw inputs, relevant layer contents, builder version and environment. Account explicitly for the intended diplomacy-only historical correction when checking compatibility. Alternatively use `REUSE_POSITIONS=0` in a fresh folder so all arms share the current layers and environment.

## Answers to the requested implementation checks

### 1. Expected-tie model and combination

- **Families/formulas:** the mappings implement the stated specification. IGO and PTA counts are restored from `log1p` before Poisson fitting; PTA-any remains binary and uses logit. ATOP types and DCA use logit; UN kappa uses OLS. W6 uses `lcap | a + b + year`; W7 uses `lcap + contig_land | year`. Here `lcap` is `log1p(capdist)`, as specified in the builder header. Sry requests outside the core arm are refused.
- **Alignment/symmetry:** the doubled frame places original orientations first and reverse orientations second. `obs(m)` restores fitted values to that frame before paired residual averaging. This ordering is coherent. The absolute symmetry assertion is useful but insufficient for SRY1/SRY2, particularly after normalization magnifies small residuals.
- **Perfect prediction:** for observations genuinely removed as perfectly predicted binary all-zero/all-one or Poisson all-zero FE groups, assigning expected=observed gives the intended limiting residual zero. Such explicitly zeroed rows cannot create a positive edge themselves. This does not justify assigning zero to every omitted row, or accepting a retained but unconverged fit.
- **Existence/missingness:** pre-existence years, all-zero layer-years, and missing UN outcomes are left `NA`. They are excluded from the applicable-layer denominator. A layer with an observed zero outcome in an otherwise active year remains applicable. Thus an empty layer-year is omitted, rather than counted as an applicable zero layer; that is the actual convention to document. A nonzero constant layer-year is admitted and needs SRY1's protection.
- **Combination:** layer scores are signed, orientation-averaged and divided by their layer-year maximum absolute value; selected available layers are averaged per dyad; only then is the positive part taken and divided by the year's maximum positive score. This matches the requested order. With no applicable layer the edge is zero. It is not the old Euclidean positive-residual construction and does not implement PN centrality.

The note correctly qualifies cross-layer scaling and the different DCA source. These remain book-aligned variants, not demonstrated exact reproduction of an independently specified book measure. The supplied main-combination agreement does not cover the numerical and omission cases above or the untested in-force sry arm.

### 2. PTA-any, IDs and downstream families

No new P1/P2 found in this plumbing.

- `pta > 0` is applied before the logarithm branch; W1–W4 therefore receive the requested unlogged binary variable. The in-force path still uses `pta_inforce` and count logging. Existing W1–W4 calculations are otherwise unchanged.
- Only the 192 PTA-containing combinations are built for `any` and `inforce`; PTA-free duplicates remain represented under `all`. IDs and metadata contain distinct `ptaany`, `ptainforce`, W6 and W7 tokens.
- R03 validates the new files without a weighting-name whitelist. R04 extracts the middle layer field, so no-ATOP selection works for the new weightings and PTA coding. R07 uses matrix columns/spec metadata rather than a W1–W4 whitelist. R11 extracts weighting and PTA coding generically and keeps them separate in `by_family.csv`.
- R11's family-balanced cell remains **weighting × centrality × diplomacy**, not PTA coding or ATOP status. The new variants enter the existing cells as designed; this is not equal weighting of PTA codings. Preserve that stated interpretation.
- The immutable Phase-A ID list is not expanded automatically. Full-grid sweeps receive all new IDs, no-ATOP selection is dynamic, and the random-bootstrap draw can include new measures; Phase-A-only checks do not thereby become a new representative sample of the enlarged grid. This is a coverage qualification, not an ID bug.

Counts reconcile: old 7,686 + W1–W4 PTA-any 1,536 + sry-all 1,536 + sry-inforce 768 + sry-any 768 = **12,294**.

### 3. Runner rebuilding and fingerprints

- Explicit old-file enumeration prevents accidentally linking newly named arms from a reuse folder. The four new arm calls have distinct stamp names and output filenames, and `fixest` is required even with reuse enabled.
- For locally generated arms, `AST` contains current layer-file content, R02 content and R/package environment. A builder change invalidates these arm stamps; a completed arm with an identical stamp may be reused within the run directory. “Always built in the run folder” therefore means locally generated rather than externally linked, not rebuilt unconditionally on every restart.
- `EXPECT` includes all ten additional files. Their bytes enter `SHA`; code, config, effective type scheme, layers, network-input identity and environment enter the global fingerprint. The full production multiverse uses `POS`, so both bcd and init receive the new arms.
- On an already fingerprinted folder, a builder change can rebuild local arms **before** the later global fingerprint rejection. Use the requested fresh output folders. This is not a silent reuse path for the new arms, but the late global check is not a guarantee of zero mutations before rejection.
- External linked-file compatibility remains unverified: SRY4.

### 4. Dropside wording and release gate

`R/06_make_variants.R:164` now distinguishes the targeted hard/reviewed modes from whole-dispute modes in the message. The selector and event removal are unchanged. The prior P3 wording issue is fixed.

Before production, add focused acceptance checks for: an exactly fitted nonzero constant layer-year; mixed and wholly perfectly predicted FE groups; missing required predictors; deliberately insufficient GLM iterations; signed-layer cancellation; and the in-force sry arm. They should establish zero numerical ties, explicit missingness/failure handling, correct applicable-layer denominators and stable row/orientation alignment. Also verify rejection of incompatible linked-layer provenance, or rebuild all arms together.

**Release decision:** no-go for the current expanded build. After the numerical/fit-validity fixes and focused checks, a fresh rerun with separate bcd/init output folders is appropriate; use `REUSE_POSITIONS=0` unless linked-file lineage has been verified. This decision does not revoke the prior stage-1 recheck for its previously reviewed measurement grid, and it makes no judgment on the deferred stage-2 memo.
