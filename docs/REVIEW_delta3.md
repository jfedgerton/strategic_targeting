# Delta 3 review — D1–D6 and intensity samples

20 November 2026. **Verdict: NO-GO for the stated content-based provenance guarantee until the persistent hash cache is corrected.** The earlier manifest, status-record, helper-return, position-arm and reporting-test defects have been addressed. I found no new demonstrated error in the intensity filtering. One new P2 remains: the cache can silently reuse a hash for different bytes across restarts.

## Scope and verification limits

Reviewed the updated `/Users/jarededgerton/Downloads/targeting_multiverse.zip`, SHA-256 `3735f7e1aa61b5067f69c8227c602d867ba28169c80a48c784c280da80ac3cf5`, HEAD `7711a863ac60106634bfafdb02408dcec53ca82e`.

I followed `docs/PROMPT_gpt_review_delta3.md`: this is a narrow review of the changes since delta 2, not another historical adjudication or audit of unchanged estimators. I inspected the changed runner, config, intensity additions in R05/R06, reporting-gate integration in R08, new gate helper/test, new list files and relevant README changes. Git-object comparisons established the source delta. I read the small portions of downstream code needed to trace the new sample integration. I did not re-review the historical ledger, frozen codings or declined D7 screen.

**No R scripts were run.** I checked the runner/config with `/bin/bash -n` under this Mac's GNU Bash 3.2.57, exercised isolated shell control-flow cases without writing files, checked the actual `date -r <file> +%s` command on this machine, and used Python for read-only archive comparisons and a cache-collision example. I did not execute the runner, its builders or the supplied integration tests. The prompt's Linux TEST/restart/failure-injection results are supplied evidence, not results independently reproduced here. Only this report was written.

All source references below are one-based lines in the updated archive.

## Disposition of each requested item

| Item | Status | Source evidence and assessment |
|---|---|---|
| D1: only currently planned sweeps in the manifest | **Fixed, subject to correct input hashes** | `run_all.sh:88` clears the current-task list each invocation. Lines 143–144 record the current content/options key before either cache reuse or execution. Lines 253–256 build the manifest only from those keys, deduplicated. If the requested strict-list replacement fails, its new key has no success record; the old key is not selected. This directly addresses the prior stale-result failure. The new hash-cache weakness below can undermine the premise that different inputs get different keys. |
| D2: Bash 3.2 compatibility / replacement hash cache | **Compatibility fixed; provenance partly fixed** | No associative arrays remain. Lines 32–36 use ordinary scalar variables and files; the runner/config pass native Bash 3.2 syntax checks. `date -r <file> +%s` works on this Mac. I found no new Bash-4-only construct. However, keeping hashes across restarts based solely on size and second-resolution mtime introduces E1 below. |
| D3: one atomic status record | **Fixed** | `run_all.sh:153–156` finds superseded records through their own CSV-path field, removes them, and installs the replacement with temporary-file/rename. There is no second `.path` file required for correctness. A crash after cleanup but before record installation can require a refit, but it no longer creates an orphan that defeats cleanup and produces duplicate manifest paths. |
| D4: ids helper failure propagation | **Fixed** | `run_all.sh:110–113` removes stale expected ids, checks R generation, move and stamp explicitly, and verifies freshness before returning the path. Both callers check the return status: the initial loop at 114 and sweep publication at 152. The old unconditional-success path after failed R generation is gone. |
| D5: position-arm provenance | **Fixed in construction, subject to correct hashes** | `run_all.sh:64` builds the arm signature from layers, R02 and environment. Lines 65–70 require every expected output and a matching arm stamp, otherwise remove/rebuild the arm, check outputs and install the stamp. Lines 71–74 cover both core PTA arms, diplomacy and IGO. This fixes the interrupted first-build mixture scenario if its signature inputs are hashed correctly; E1 can hide a changed signature input. |
| D6: discriminating reporting-gate regression | **Fixed in source/test design** | `R/lib/report_gate.R:3` implements the finite-observed guard. R08 sources it at line 17 and calls it at 103. `tests/test_resample_edge.R:10–14` directly tests NA observed/all usable draws, a valid observed value, small B and insufficient usable draws. Removing only the finite-observed condition now makes the first assertion fail. Lines 15–17 also demonstrate that the old gate admits that invalid case. The synthetic fixture remains useful for pooled M0 survival. Tests were inspected, not executed. |
| Intensity samples | **Implemented correctly for the specified event-data fields; data counts unverified** | R05 lines 480–488 attach event-level state-fatality and dispute-intensity fields. R06 lines 112–124 filter whole strata from their chosen row and count unknowns. Runner lines 168–173 build each window's variants from that window's fatal frame and sweep all configured measures. Config line 25 requests jf and fd2–fd6. Details below. |
| Additional sensitivity lists | **Mechanically integrated** | `inputs/exclusions_noruleA.txt` contains strict plus 4339 and 611: 11 unique ids. `inputs/exclusions_anycoder.txt` contains 17 unique ids. `run_all.sh:227–229` runs both through the same stamped variant/sweep path as the four existing lists; the first-token parser correctly ignores their descriptive text. This assessment concerns parsing/integration, not the declined historical recoding. |

## E1 — P2: persistent metadata cache can return a false content hash

**Location:** `run_all.sh:30–35`; consumers at 59–64, 78–82 and 143.

The old memo was confined to an invocation. The replacement persists at `status/hcache` and treats `(absolute path, byte count, integer-second mtime)` as proof that content has not changed. Line 34 returns the previously stored hash without reading the file when those metadata match.

That is not a reliable content identity. Two eight-byte inputs provide a simple example:

```text
old: a,b\n1,2\n  -> SHA-256 prefix 492d5ea496056f1a
new: a,b\n1,3\n  -> SHA-256 prefix 8baab0330389def9
```

If the second version occupies the same path with the same timestamp second, line 34 returns the old prefix. A same-size edit within one timestamp second suffices; restoring or copying a file while preserving modification time is another ordinary way to produce the condition. The Python check confirmed equal lengths and different actual SHA prefixes; I did not modify workspace files to create a physical fixture.

**Concrete consequence:** after a prior invocation has cached a network input, replace it with a different same-size version preserving its old mtime and restart. `NETSHA` can remain unchanged, so the layers stamp matches and old layers are reused. Position signatures and the global fingerprint can consequently remain unchanged too. The runner silently proceeds as though it had checked the new input's content, but it has not. Similarly, an externally replaced positions file with matching metadata can retain its old sweep key, allowing cached results for the prior positions bytes. This is a wrong-provenance/cache-reuse path, not merely an unnecessary recomputation issue.

The D1 changed-list example is **not automatically broken** by this: history-list contents are passed as literal variant arguments, so changing those ids normally changes the variant signature regardless of their list-file metadata. The problem is the general content-hash contract for files actually passed through `h()`.

**Minimal fix:** do not reuse content hashes across invocations on this metadata test. Start with an empty invocation-specific hash cache, retaining memoization only after each file's final write in that invocation; explicitly invalidate entries if a producer rewrites a file already hashed. Alternatively, recompute the content hash whenever deciding persistent artifact freshness. Higher-resolution mtime reduces accidental collisions but does not establish correctness for timestamp-preserving replacement. There is no need to restore associative arrays to fix this.

**Targeted regression:** hash/cache an input, replace it with different equal-length bytes while preserving its mtime, then restart. Require a changed content hash and the appropriate rebuild or fingerprint refusal. Test the same condition for a positions file. The reported list-change/failure test does not cover this case.

## Intensity construction: what was checked

The implementation follows the requested state-level definition, rather than claiming directly interacting dyadic deaths:

- `R/05_dump_events.R:482–485` matches the chosen initiator–target–year event to the unique non-escalation source row. The existing uniqueness assertion at 437–438 protects this key. Although the value is repeated over alternatives, `ek2` uses the chosen target at 467, so every alternative in the same event receives the same selection metadata. The filter is not accidentally applied to alternative targets' deaths.
- `fat_level` converts negative/missing dispute codes to NA and otherwise retains the level. `jfatal` requires both supplied state-fatality flags positive and returns NA when either flag is unknown. The latter is deliberately conservative about unknowns, consistent with the prompt; it does not silently call an unknown state nonfatal.
- `R/06_make_variants.R:117–120` evaluates eligibility on chosen rows and retains all rows for eligible `ev_id`s. Thus choice alternatives remain intact, including alternatives that did not themselves suffer fatalities. Line 152 checks one chosen row per retained stratum. No event duplication is introduced by the new filter.
- Thresholds are inclusive in **code space**: fd2 keeps level ≥2, through fd6 keeping level ≥6. README line 15 and the R06 comment label their death ranges. Unknowns are excluded before collecting ids, and the log separately reports their number. The code also stops on an entirely empty sample instead of emitting a misleading successful fit.
- The variants retain original event weights and controls; they do not renormalize the surviving target events within an initiator/dispute. That is the implemented sensitivity estimand, consistent with filtering the existing fatal frame. No weight change was requested.
- `run_all.sh:170–173` iterates the full `WINDOWS × INTENSITY × MV` grid. Production MV is the full positions array at 160; TEST deliberately substitutes the main measure. Variant names include intensity and window, `d$types` adds the intensity tag, and the existing sweep naming incorporates the window, preventing these new sample outputs from overwriting one another.

No supplied RDS data or executed intensity census was available in this archive, so I cannot certify actual counts, the number of unknowns, or which high-threshold type-specific models remain estimable. The current code records per-model invalidity through the existing fitter rather than assuming all intensity samples support every surface. That is appropriate. I found no reason in this delta to change the fitting machinery.

## Bash / failure handling and production decision

Native Bash 3.2 syntax checks passed for both shell files. Isolated native-shell checks confirmed that the final false `&&` conditions in the superseded-record scan do not cause an unintended `set -e` exit, an invalid arm stamp selects rebuilding, and explicit helper failure propagates through command substitution. I did not run a mocked full pipeline because it would create files outside the sole authorized report. Syntax and isolated semantics are not a substitute for a complete native-platform mechanics test.

**No new P1 was established in this narrow pass.** D1, D3, D4, D5 and D6 are addressed in their source logic; D2 removes the demonstrated Bash-3 blocker but introduces E1 through persistent hash reuse. I would hold the production run until that small cache correction is made, then run the native-platform mechanics check and the metadata-preserving replacement regression alongside the already reported restart/failure tests. If production inputs were externally guaranteed immutable, E1 would not trigger, but the runner currently promises stronger verification and does not enforce that condition.

D7 and further historical adjudication remain outside this pass as explicitly requested. Their prior limitations are not being promoted into new blockers here.
