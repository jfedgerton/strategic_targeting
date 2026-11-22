# Stage-1 recheck — S1 and S2

21 November 2026. **GO for `TYPES=init OUT=output_init bash run_all.sh` in a separate, fresh output directory. S1 and S2 are fixed. No new P1/P2 was found in this patch.** This is a narrow source recheck, not an independently executed production or R integration test.

## Scope and checks

Read `docs/PROMPT_gpt_recheck_stage1.md` and inspected the changes from the previously reviewed 3015bbb package in `run_all.sh`, `R/06_make_variants.R` and `README.md`. The only changed R/shell files are the runner and R06. Unchanged estimators, stage-2 design and historical coding were not re-reviewed.

Archive: `/Users/jarededgerton/Downloads/targeting_multiverse_stage1_fix.zip`; HEAD `230d5858246db85b9f2e3dc3a17252a5616fc66b`; SHA-256 `29f954a29a0e7dd5a905215425be7df80f22582cca97a210c35591fe15a12d2b`. References below are one-based lines in that archive.

**No R scripts were run. Only this report was written.** Native Bash 3.2 syntax checks passed for the runner and config. An isolated in-memory shell check of config plus the exact new scheme guard, with `CORES=1` to bypass platform CPU discovery, rejected `TYPES=foo` with exit 1 and accepted bcd/init with exit 0. It performed no output-folder operations. I did not run the full runner. The supplied R/data.table synthetic results remain external evidence, not results reproduced here.

## S1 — fixed: effective scheme is validated and fingerprinted

**Evidence:** `run_all.sh:18–20`, `88–90`.

The guard immediately after sourcing config accepts only bcd/init and exits 1 otherwise. It precedes output-directory checks/creation and R execution. The fingerprint now includes the expanded `TYPES=$TYPES`, in addition to the unchanged config/code identity inputs.

Consequently, changing only the environment from bcd to init no longer preserves the fingerprint. Reusing an already fingerprinted output directory for the other scheme is refused before reaching cached influence, resampling or event-study reuse. This closes the specific wrong-population cache path identified in S1. Continue using the prescribed distinct `output_init` directory.

## S2 — fixed: promoted original joiners receive the targeted checks

**Evidence:** `run_all.sh:201–211`; `R/06_make_variants.R:123–125`, `158–165`.

DSL now unconditionally contains `hard reviewed hard_disp reviewed_disp`. Both the variant-building loop and sweep loop use DSL, so init receives the targeted hard/reviewed checks as well as whole-dispute checks. The prior init-specific narrowing has been removed.

The initiators block logs the number of retained chosen-event strata tagged hard and reviewed. These are distinct tag categories: a logged reviewed count excludes hard-tagged strata, whereas the reviewed *drop* includes both. Thus the supplied synthetic case's hard 1 / reviewed 1 log correctly corresponds to removals of 1 and 2 respectively.

The existing drop selector operates on retained chosen-row tags and removes complete event strata; it does not require the final type to remain joinA/joinB. The fix therefore reaches the promoted coalition-originator events that motivated S2. A zero-removal case still builds/runs and logs zero; it is no longer skipped based on scheme. No new role recoding or weight modification was introduced.

## Wording and new issues

The new README Stage 1 section (`README.md:20–33`) correctly distinguishes the narrower default Step B population, all/pair within the retained init frame, the extended four-type population, artificial-label permutation, descriptive event studies and recomputed coal−bil influence. These changes address the requested population qualifications.

One cosmetic qualification remains: R06's shared dropside log at line 164 describes removed strata as joiner/promoted-originator strata even for the `_disp` modes, which can also remove ordinary bilateral/originator strata. The mode label and README make the distinction available, so this is **P3 wording only**, not an incorrect drop rule or production blocker. If edited later, use generic “removed strata” for the whole-dispute modes.

No new P1/P2 was established in the reviewed changes. The previous stage-2 design concerns are outside this recheck and do not block the specified stage-1 run. The GO recommendation relies on source inspection and the supplied prior mechanics evidence; it does not certify the reported numerical estimates or a production run executed here.
