# Third-pass review: coding and delta fixes

20 November 2026. **Verdict: NO-GO for an unattended production run as packaged.** Most of the requested statistical fixes are present and look correct in source. The remaining substantive runner problem is that a failed replacement sweep can leave a superseded result in the published summary. There are also restart/provenance gaps and a demonstrated incompatibility with this Mac's default Bash. These are specific fixes, not a recommendation to rewrite the analysis.

## Scope and order

I read `docs/PROMPT_gpt_review_delta2.md`, then completed and saved Part 1 from `docs/coding/blind_coding_sheet.md` before opening the other package documents or inputs. The frozen result is [coding/coder3_gpt.md](coding/coder3_gpt.md); I have not revised it to match the ledger.

**Blinding limitation:** this chat already contained earlier case assessments. Withholding the updated ledger and coder 2's answers does not erase that exposure. This is a separately recorded application of the blind sheet, not a fresh, fully blinded third coder. Do not describe it otherwise in the paper.

Part 2 covered `README.md`, `docs/decisions.md`, `docs/historical_exclusions.md`, `docs/appendix_historical_exclusions.md`, `docs/coding/coder2_claude_agent.md`, the preceding review, the four `inputs/exclusions_*.txt` files, `run_all.sh`, the changed sections of `R/05_dump_events.R` through `R/10_event_study.R`, all of the new `R/11_summarize.R`, and `tests/test_resample_edge.R`. I inspected the layer/build integration and relevant prior Git versions to check the delta. This is a source/delta audit, not a fresh full audit of every unchanged function or the underlying event data.

The ZIP refreshed during the review. Its final SHA-256 is `6c6bf9e966392d6455704915b2fbb26d147ba4104beca5dc3a022a340e3408e5`, and its HEAD is `2bcffeb99adf491e86ef5ea258f620142b3b0362`. Comparison with the initially inspected `62b09d57939176fb70f157fdb8d2638325749e62` showed that the refresh changed the historical ledger and added the appendix draft; code, inputs, and blind sheet were unchanged. Part 2 considers those updated documents. References below are one-based lines in the final ZIP, relative to its repository root.

**No R scripts were run.** Checks were static inspection, Python reads/comparisons, shell syntax checks, and isolated shell-only reproductions. I did not independently reproduce the claimed test pass or the claimed zero mixed-source war links. Only the two requested Markdown deliverables were written; the archive was not extracted into the project.

## Part 1: coding and disagreements

My 23 sheet rows yield **11 yes, 8 no, 4 borderline**; 1286/1703 is one row. Coder 2 and I agree on the qualification category for **20/23 rows** and on UN classification for **23/23**. The three qualification differences are Shaba, Cominform–Yugoslavia, and Chincha. Both identify only Korea and the Gulf War as UN-authorized/recommended operations.

This is agreement under the supplied sheet, not validation of the final ledger's amended rule. Ledger Ruling A (`docs/historical_exclusions.md:10`) was absent from that sheet. The updated appendix appropriately discloses that the ruling followed coder 2; it does not make either coder's original answers applications of that ruling. The ledger's “15 of 23” statistic at line 49 also mixes qualification, list membership and force-count disputes. Specify the comparison variable before reporting an agreement rate.

Here is the complete qualification comparison. “Strict” includes its borderline members; “expanded only” means absent from strict. Source citations and confidence levels for my decisions are in the frozen coding file.

| MID / case | Coder 3 | Coder 2 | Final ledger | Resolution / remaining difference |
|---|---|---|---|---|
| 51 Korea | Yes | Yes | Strict; UN | Agree. |
| 3957 Gulf War | Yes | Yes | Strict; UN | Agree on eligibility. Coder 2 questions the exact 14-state count; it cannot affect the threshold here. |
| 31 Boxer | Yes | Yes | Strict | Agree; public relief rationale qualifies despite later mixed aims. |
| 145 Nanking | Yes, conditional on supplied count | Yes | Strict | Trigger agrees. Histories document British/US firing; the additional French/Japanese qualifying MID actions need participant-level documentation. Retain the operational MID coding provisionally, with the count discrepancy explicit. |
| 146 Navarino | Borderline | Borderline | Strict, marked borderline | Both coders find humanitarian and autonomy/settlement aims difficult to separate. Report strict and strict_noborder; do not count the borderline as unqualified agreement on yes. |
| 3552 Iraq no-fly zones | Yes, b/d | Yes, b primarily | Strict, d/b | Same qualification; primary trigger ordering differs without changing membership. Four MID force users versus three clearly identified historical force users remains a documentation issue, not a threshold reversal. |
| 4137 Kosovo | Yes | Yes | Strict | Agree; no operation-specific UN force mandate. |
| 4283 Afghanistan | Yes | Yes | Strict | Agree; do not borrow later ISAF authorization for the initial operation. |
| 3058 Grenada | Yes | Yes | Strict, marked borderline | Both accept protection of nationals as a principal public rationale. The ledger's contested-legality reason for borderline status is not itself an exclusion in the sheet. Identify uncertainty about main justification or qualifying actions if that is the intended basis; otherwise label its removal a separate disputed-legality sensitivity. |
| 4273 Iraq 1997–2003 | Yes at MID level | Yes | Expanded only | Under the supplied MID-level unit/count, d qualifies. Expanded-only treatment adds a mixed-episode restriction not stated as a general strict-list requirement. Define and apply that restriction uniformly; retain this as an explicitly adjudicated mixed-episode sensitivity. Neither yes means Desert Fox alone had three force users. |
| 4251 Sierra Leone | No | No | Expanded only | Agree on failing the closed strict trigger. Its inclusion is explicitly a relaxation in expanded, not a coding disagreement. Embargo authorization is not authorization for the ouster. |
| 2150 Shaba | Borderline | Yes | Expanded only | I lack sufficient contemporaneous public-justification and country-action evidence for an unconditional yes; coder 2 accepts a/e. Obtain those records. Ruling A separately explains final expanded-only placement; it does not settle the evidentiary disagreement under the sheet. |
| 135 Mexico | No | No | Out | Agree on main debt aim. MID count three versus histories emphasizing French combat requires action-level reconciliation; do not equate no firing with no MID use-of-force action. |
| 89 Spain 1823 | No | No | Out | Agree. |
| 2208 China 1926 | No | No | Out | Both no; my answer respects the supplied zero provisionally. Coder 2 cites related British/Japanese firing. Map those episodes to this MID before calling the zero erroneous; neither coder establishes three qualifying users here. |
| 1129 Hodeida 1934 | No | No | Out | Agree; protective deployment does not establish three force users. |
| 200 Suez | No | No | Out | Agree on the coalition's principal rationale; Israel's separate rationale cannot automatically be assigned to all three states. |
| 1286/1703 Cominform–Yugoslavia | Borderline | No | Out | I leave the coalition's qualifying public rationale and action count unresolved; coder 2 regards the main aim as ideological pressure. Keep out of strict absent positive evidence, and investigate the two MIDs separately. Uncertainty is not a recommendation to include. |
| 210 Spanish Civil War | No | No | Out | Agree; a qualifying principal public justification has not been established. |
| 4339 Congo | Yes | Yes | Out | Both coders accept aggression/collective-defense rhetoric under the sheet. Final exclusion follows new Ruling A, not contrary coding of that sheet. Label the amendment and, preferably, report the unamended-trigger sensitivity as well. |
| 1482 Chincha | Borderline | Yes | Out | Two issues: Ruling A changes the alliance interpretation; independently, the sheet omits Chile and three qualifying force users are unproved. Coder 2's own note identifies only Chile/Peru fighting. Obtain the MID participant actions before accepting a positive threshold classification. |
| 611 Vietnam | Yes | Yes | Out | Both coders accept the public aggression rationale under the original sheet, with alliance ambiguity acknowledged. Exclusion follows Ruling A. Treat that as a substantive rule amendment, not independent coder confirmation of out. |
| 607 Syria 1957 | No | No | Out | Agree; threat/ideological concern and deployments do not establish the qualifying coalition operation. |

**Recommended historical resolution.** Keep all four retrospective lists visible. If Ruling A is retained, describe strict as the PI-adjudicated, narrower multilateral-norm list. Have coding under that final rule explicitly distinguished from coding under the original sheet; do not silently overwrite either frozen coder. Ordinary collective defense can fit trigger (a), so excluding it requires the additional stated distinction. This is a judgment call about the construct, not a software error. Ruling B is defensible as an operational definition, provided the participant actions/dates are shown. MID action ≥14 is broader than firing weapons; several apparent historical count conflicts cannot be resolved from battle narratives alone.

## Verification of requested fixes

“Fixed” here means that the requested behavior is implemented in the inspected source; it does not claim an executed integration test.

| Item | Status | Evidence and qualification |
|---|---|---|
| N1 rebuild branch | **Fixed in source** | `run_all.sh:40` conditionally adds fixest/igoR; lines 58–59 use `printf -v` and check all four diplomacy files; lines 60–63 check builder status. The separate new Bash-version issue is below. |
| N2 provenance | **Partly fixed** | Lines 54–55 stamp layers against raw inputs/R01; 67–72 fingerprint position contents, layers and environment; 86–101 stamp/check sets and ids; 151–152 stamp variants; 127 hashes sweep inputs. The legacy-root guard is at 22. Superseded-task handling and pre-fingerprint position lineage remain incomplete, below. |
| N3 dump staging | **Fixed** | `R/05_dump_events.R:497–503` accepts the run-specific destination and publishes by temporary file/rename. `run_all.sh:111–114` checks process status and frame metadata. Shared staging is no longer the normal runner path. |
| N4 resampling | **Fixed in source; test coverage incomplete** | `R/08_resample.R:91–93` preserves pooled M0 when a draw loses a type; line 102 requires a finite observed estimate; 104–114 apply/report that gate. The added fixture does not establish coverage of the original high-valid-draw/invalid-observed branch; see D6. |
| N5 manifest/status | **Partly fixed** | Manifest rebuilding replaces append-only accumulation; sweeps check process status and ERROR rows (`run_all.sh:130–137`). Summary files use temporary files/rename (`R/11_summarize.R:68–81`). But current-task selection, the two-file status transaction, and an unchecked helper failure remain, D1/D3/D4. |
| N6 summary | **Fixed for the specified checks** | `R/11_summarize.R:21–38` rejects duplicate paths, unknown ids and unequal expected/seen sets. Lines 60–67 compute requested and fit-valid denominators before filtering computable Wald tests. M0t has its own census at `R/07_fit_sweep.R:392–393`. These checks cannot detect an obsolete result with the correct measure ids. |
| N7 mixed-source war links | **Fixed by failing closed** | `R/05_dump_events.R:307–321` collects mixed-only fallback links and stops before export if any occur. This prevents using an unverified mapping; it does not itself repair one. Claimed current count zero was not reproduced here. |
| Step B | **Fixed as requested** | `R/05_dump_events.R:487` supplies GML start; `R/07_fit_sweep.R:232–269` filters original decisions to that year and reports exclusions/decision counts. Lines 80–103 write the severity/membership cross-tab with unknowns. Eventual severity remains an association outcome, not a contemporaneously observed onset choice. |
| Event study | **Fixed as requested** | `R/10_event_study.R:25–49` excludes missing/conflicting onset references by default and provides explicit fallback sensitivity; 105–108 write excluded ids and type-by-k event/weight composition. Runner lines 217–223 request both ordinary and fallback versions. |
| P2-11 aggression | **Fixed as requested** | `R/06_make_variants.R:39–78` flags fallback/conflicting starts and upper-bound censoring; complete excludes them, complete_base uses the same retained strata without the control, censored retains the lower-bound sensitivity. Runner lines 153–155 request all three. |
| Documentation / influence arity | **Fixed** | README lines 7/17 state availability restrictions and line 11 calls pooling an assumption; `docs/decisions.md:4` explicitly rejects equivalence inference from non-rejection. `R/09_influence.R:53–54` accepts one type or a two-type contrast and rejects other arities. Provenance/restart claims in README 31–32 still overstate the remaining runner guarantees. |
| History lists | **Mechanically fixed; substantively partial** | All four supplied lists agree with the ledger: 2 UN, 9 strict, 7 strict_noborder, 12 expanded, no duplicate ids within lists. Runner 200–202 runs each separately and propagates list contents into the variant stamp. Rule amendment, borderline treatment, and candidate-screen completeness remain as described here. |

## Remaining and new findings

### D1 — P1: a failed replacement can publish a superseded historical exclusion result

**Location:** `run_all.sh:127–137`, `200–202`, `226–232`; `R/07_fit_sweep.R:385–387`.

Concrete sequence: successfully run `rob_drop_hist_strict`; change the approved strict list, which is deliberately allowed to trigger a stamped variant rebuild without changing the global code/config fingerprint; restart. The variant content and sweep key change. If the new sweep exits before publishing its replacement CSV, the old CSV and old successful status record survive. Line 133 records a gap and returns; obsolete-record cleanup at 136 happens only after success. The final manifest scans **all** old successful status records, so it includes the old strict result under the current strict group. Its measure ids still equal the expected Phase A ids, so the new summarizer's exact-id checks pass. The resulting summary answers the old exclusion question while the inputs now contain the new list.

The runner does exit 2 with `LOOP_FINISHED_WITH_GAPS`; it does **not** falsely announce `LOOP_DONE`. That warning reduces the chance of misuse but does not make the published current summary's provenance correct.

**Minimal fix:** build the manifest from the current invocation's planned logical tasks and their current content keys, including reused successes. An older key for that logical task must not qualify after the requested inputs change. Keep obsolete artifacts for inspection if desired, but exclude them from current summaries. Test success → list change → replacement failure before CSV publication → restart.

### D2 — P2: new Bash requirement blocks the documented command on this Mac

**Location:** `run_all.sh:31`; `README.md:22–27`.

The new `declare -A HMEMO` requires Bash 4+. This machine's default `bash` is GNU Bash 3.2.57. An isolated shell-only reproduction of `set -euo pipefail; declare -A HMEMO; echo reached` exits 2 with `declare: -A: invalid option`; it never reaches the echo. Thus the documented `bash run_all.sh` stops before the first R preflight on this environment. `bash -n` does not detect this compatibility failure.

**Minimal fix:** explicitly require and preflight modern Bash before output creation, and document the interpreter used in production, or implement a portable cache. A Linux host already using modern Bash is unaffected; this is a demonstrated local execution blocker, not evidence of wrong fitted coefficients.

### D3 — P2: status record plus `.path` is not an atomic transaction

**Location:** `run_all.sh:136–137`, `227–229`; `R/11_summarize.R:21`.

Line 137 atomically installs the primary record, then separately writes its `.path` lookup. A crash between those operations leaves a trusted record without a lookup file. A later changed-input success overwriting the same CSV removes obsolete records by scanning only `.path` files, so it misses the orphan. Both primary records then enter the manifest. The summarizer correctly refuses duplicate paths, making a recoverable interruption a persistent summary failure.

**Minimal fix:** use the single atomic record as the source of truth and locate obsolete records from its path field. Alternatively, repair/reconcile auxiliary lookups before trusting records. Do not require a second non-atomic file for correctness.

### D4 — P2: `ids_of()` can return success after its R process fails

**Location:** `run_all.sh:97–98`, `135–137`.

The generation call is the first command in an `&&` chain. Failure skips its move/stamp and does not trigger `set -e`; the unconditional final `echo "$I"` succeeds. Both the initial caller and `I=$(ids_of $PS)` therefore receive a successful helper result. With an older ids file still present, downstream code can publish a new status record using old expected ids. Some resulting mismatches will be caught by `cmp` or the summary, but not every call/positions set has an independent immediate comparison, and the blanket claim that every R failure is checked is false.

**Minimal fix:** explicitly fail the helper on generation/move/stamp failure, require a readable file with the current stamp, and propagate failure at both call sites. This is ordinary shell control flow; no R execution is needed to establish it.

### D5 — P2: interrupted initial position builds can acquire mixed provenance

**Location:** `run_all.sh:54–72`.

With `REUSE_POSITIONS=0`, the initial build decides which position arms to reuse from file existence only. The first fingerprint is written only after every arm exists. If a first invocation finishes one arm then stops, and a network input changes before restart, layers rebuild against the new input, but the old complete arm is skipped while missing arms use the new layers. A new fingerprint then blesses that mixture. There is no previous fingerprint to reject it. This is narrower than the old global provenance bug: an unchanged-input restart is fine, and an already fingerprinted folder rejects changed inputs.

**Minimal fix:** stamp each position arm against layer content, R02/options and relevant environment, or establish/enforce an initial-build input identity before generating any arms. An output hash proves which bytes were used; it does not establish which layer version generated each arm.

### D6 — P2 validation gap: the N4 fixture does not force the original reporting failure

**Location:** `tests/test_resample_edge.R:2–4`, `16`, `21–28`; `R/08_resample.R:102`.

The single rare-type event is a good test of preserving pooled M0 when bootstrap draws lose a type. But the reporting-gate regression requires an invalid observed quantity **and at least 90% usable simulated quantities**. The fixture never asserts that latter condition for an M1 quantity. A single-event type with a free quadratic surface can remain separated under within-event dR permutations; if M1 simulated fits are mostly invalid, the old 90% rule alone also suppresses its interval/p-value. Removing the new finite-observed condition could therefore leave the asserted absence of M1 reporting unchanged.

**Minimal fix:** add a direct gate test with an NA observed estimate and an otherwise fully usable simulated column, or a fixture that asserts the necessary high-valid-draw condition before testing suppression. Verify the test fails when only the finite-observed guard is removed. This is a coverage limitation; the current guard itself is correct. I did not execute or mutation-test the R fixture.

### D7 — P2 design/completeness: the candidate screen does not guarantee the force rule

**Location:** `docs/historical_exclusions.md:8–9`; `docs/appendix_historical_exclusions.md:5`.

The eligibility rule requires three coalition-side states with MID action ≥14. The screen instead uses three initiating-side initiators or three target-side **joiners** in the fitted fatal-event dump. Those predicates are not equivalent. An original defender plus two defending joiners is already a three-state force coalition, but has only two target-side joiners and can fail this screen. Multiple original defenders create further cases. Participants can also disappear from fitted events because of event/choice availability even when their MID actions qualify.

This establishes a logical completeness gap, **not a demonstrated omitted historical case**: the raw participant/event join was not executed here. It limits the claim that both sides have been comprehensively screened.

**Minimal fix:** enumerate candidates directly from the MID participant table, counting qualifying distinct states by side and dispute, then join all eligible disputes to the fatal sample through source-dispute membership. Compare that candidate list with the recorded screen and adjudicate any additions. Retain the fitted-role counts as descriptive metadata.

## Checks, minor documentation and production decision

I found no additional demonstrated numerical defect in the changed estimator, Step B filter, aggression comparator or event-study exclusion logic. The hash memo is intentionally populated in the current shell before sweep command substitutions; I did not find a concrete within-invocation stale-hash bug in its normal artifact write order. The `&&` conditionals around optional builds are not, by themselves, erroneous `set -e` exits. The specific `ids_of` failure masking above is different and reproducible from shell semantics.

Minor P3 documentation: ledger line 63 still calls coder 3 “blind” and “pending”; line 65 still asks for a second coder after presenting coder 2. Update these with the actual sequence and the prior-exposure limitation. Shaba's source remains “to add” at line 30 even though the two coder files now provide candidate evidence; add a checked public-justification source, not merely a private retrospective inference.

**Production decision:** hold the unattended run until D1 and the status/provenance fixes are made, and select a supported Bash interpreter. On modern Bash, D2 is already avoided, but that alone does not resolve D1/D3–D5. Then run an authorized small mechanics test plus targeted restart/failure tests; this review has not run them. Preserve the corrected N4 code and strengthen its regression test. Resolve or explicitly bound the historical screen before calling the ledger comprehensive, and describe Ruling A as an amendment rather than inter-coder confirmation.

Out of scope for this pass: choosing a different causal estimand, redesigning the utility surface, or changing the underlying MID fatality/direction construction. The historical borderline and alliance classifications above are construct/adjudication choices, separated from the runner bugs. No wholesale code rewrite is needed.
