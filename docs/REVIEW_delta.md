# Delta review — 20 November 2026

**Verdict: NO-GO for the unattended production run as packaged.** The original window blocker is fixed, but the rebuild path now has a shell failure, cache provenance remains incomplete, and the new per-model resampling logic can manufacture a small diagnostic tail value when the observed statistic is missing.

Reviewed `/Users/jarededgerton/Downloads/targeting_multiverse.zip`, SHA256 `cb8282a6b8e6151629044ebfa8811271ca12f66858358a69dd29f4e7ce9fef0c`, whose recorded HEAD is `982f1bcfd5b0b3cbe926783eceb24956ac42793b`. Locations below are archive-relative paths and current one-based lines. I compared changed files against the archive's Git objects at `421b8813d82eb165897301382cc0cab3f6300faf`, the previous review version. I read the requested documentation, exclusion draft, configuration, rebuilt runner and changed script sections; unchanged modeling code was not re-audited.

No R script was run. Checks were source/diff inspection, shell syntax checks, and a minimal shell-only reproduction of the dynamic-assignment failure. External event/position RDS inputs are absent, so numerical rebuilds, reported changed-event counts and case-to-MID mappings remain unverified. Primary historical sources were consulted for the exclusion assessment. Only this report was written; the archive, source files and previous review were left unchanged.

## 1. Disposition of every previous finding

“Fixed” below means the specified source change is present and coherent; it does not mean an R rebuild was executed.

| Previous finding | Status | Current evidence and remaining limitation |
|---|---|---|
| **P1-1: windows 0–3 only** | **Fixed** | `R/05_dump_events.R:68` accepts 0:6; lines 241–246 limit the old Germany–Poland assertion to windows ≤3 and log longer windows. `run_all.sh:22–23` validates configured windows and the main window. |
| **P1-2: test pollutes production** | **Partly fixed** | `run_all.sh:15` separates default final roots into `output_test` and `output`; TEST enters the fingerprint. However, `R/05_dump_events.R:484–490` still writes both modes to shared `output/rel_sweep` staging, and the runner can accept a stale/partial staging file after failure. See N3. |
| **P1-3: markers/provenance/completion** | **Partly fixed** | `run_all.sh:31–54,89–99` adds explicit positions filenames, link checks, content fingerprints, manifest-listed sweeps and ERROR-row exclusion. Remaining problems: omitted selection/exclusion inputs and layer/package provenance; adoption of existing output without a fingerprint; nontransactional manifest/status updates; and unchecked process failures. See N1–N3, N5. |
| **P1-4: collapsed side flag** | **Partly fixed** | `R/05_dump_events.R:289–298` derives primary-dispute initial side from participant spells rather than the collapsed Boolean. This fixes the demonstrated side-flag failure. Lines 305–306 still use collapsed war links whenever no single-source row exists; that fallback can belong to another dispute. See N7. The claimed 17 changed source-side values and passing case checks cannot be independently reproduced from this ZIP. |
| **P2-5: Step B outcome/risk set** | **Partly fixed; now a narrower analysis** | `R/05_dump_events.R:328,457,474–480` preserves original joiner status and separately creates primary-dispute GML fatality, with negative/missing codes treated as unknown. `R/07_fit_sweep.R:223–228` excludes original joiners and unknown outcomes by default. This removes the central late-joiner problem. It still does not assert that retained pair-event dates equal dispute onset or provide one explicit onset decision key; see the Step B assessment below. It is not an analysis of fatal escalation among joiners. |
| **P2-6: omitted type intercept test** | **Fixed** | `R/07_fit_sweep.R:244–253` distinguishes the gap-only surface test, adds gap-plus-intercept `all_differ`, and exports type-intercept contrasts. Absent types are removed through RP at lines 226–228. The default originator-only comparison should therefore have fewer types than the all-entry comparison; report its actual reference and degrees of freedom. |
| **P2-7: event-study alignment** | **Partly fixed** | `R/10_event_study.R:40–47` offers onset versus entry alignment; `run_all.sh:173–177` runs balanced, unbalanced and pre-only versions. Excluded event IDs and revised labels are written at R/10:101–103,115. But unknown onset dates still revert to entry and later conflicting starts are clamped to entry at line 43. An output called “onset” consequently mixes references unless those cases are absent or separated. |
| **P2-8: influence fit validity** | **Fixed for the intended endpoint contrasts** | `R/09_influence.R:42–49` rejects warned/nonfinite/iteration-limit fits; lines 59–74 retain every attempted unit and status, rank valid refits only, and publish top IDs before the LOO completion file. `run_all.sh:151–153` checks both. Contrast names are checked at R/09:53. A minor remaining guard: require one or two type names, since the validator currently accepts a three-type expression that `cvec` cannot evaluate correctly. |
| **P2-9: per-model resampling/fit errors** | **Partly fixed, with a new P1 edge case** | `R/07_fit_sweep.R:285–291` catches thrown coefficient-fit errors per model. R/08 preserves partially valid full-sample and draw results, checks t−1 variation, refuses typed-control type permutations, and exports draw IDs/reasons. However, R/08:90 still discards valid pooled information whenever any type vanishes, and its reporting gate ignores observed-fit validity. See N4. |
| **P2-10: SD population** | **Fixed for the stated event-weighted choice-set population** | `R/07_fit_sweep.R:142,170` calculates n alternatives after restrictions and uses w/n alternatives in the change-gap SD. Each event's set now contributes its D8 weight. This is still a choice-set mixture SD, not a chosen-target-only or state-level tie-loss SD; label it accordingly. |
| **P2-11: aggression coverage** | **Not fixed; explicitly acknowledged** | `R/06_make_variants.R:57–69` still uses truncated counts beyond GML's endpoint and retains `aggr_censored`. See the recommendation below. |
| **P2-12: promoted joiners/whole-dispute drops** | **Fixed for the runner's bcd/bcjs paths** | `R/05_dump_events.R:328–329` tags original joiners including window promotions. `R/06_make_variants.R:119–124` adds all-source `hard_disp`/`reviewed_disp` exclusions. The original-joiner flag is reconstructed before export on these paths; it is not simply the final promoted role. |
| **P2-13: appendix collector** | **Partly fixed** | `R/11_summarize.R:18–50` is manifest-only, rejects duplicate estimate keys, supplies plausible/failing/all and family-balanced summaries, and separates M0/M1 validity counts. It checks ID counts, not expected ID identities; Wald “requested” excludes fits that never emitted a Wald row; M0t lacks its own census. See N6. |
| **P3-14: options/documentation** | **Mostly fixed** | `R/07_fit_sweep.R:54–65` rejects the named ignored-option examples and constrains Step B options. README:11 acknowledges pooling as an assumption and line 17 states eligibility restrictions. Minor leftovers: README:7 still says all system members are alternatives; `docs/decisions.md:4` still motivates pooling with nonrejection; duplicate/extra positional arguments are not comprehensively checked. These are not production blockers. |

The additional B/C responses are also accounted for:

- **Permutation labeling: fixed as labeling.** R/08:8–9,107–110 and README:46 call these artificial-assignment diagnostics and export `p_diag`; the runner names them `diag_*`. This does not establish exchangeability, and the new missing-observed-statistic bug must be fixed.
- **Measurement perturbation grid: fixed as specified.** R/07:54,143–157 and runner:137–142 add both signs, separate noise seeds, five configured draws and rank correlations. These remain fragility exercises, as the revised README now says.
- **Influence-selected versus history-selected drops: separated in the runner, historical coding not yet defensible as a locked list.** Runner:155–158 keeps the former and waits for approved files for the latter. The historical list needs the decisions in section 3.
- **Reporting: improved but incomplete.** Estimate-first and family summaries are present; requested-test denominators and ID identity checks remain unresolved.

## 2. Remaining/new actionable P1/P2 problems

### N1 — P1: the new rebuild branch executes an assignment as a command

**Location:** `run_all.sh:39`.

`NEED_CORE_$PT=1` is not a shell assignment after expansion. With PT=all, Bash tries to execute a command named `NEED_CORE_all=1`. I reproduced this isolated expression without running any R: it exits 127 with `command not found`. `bash -n` does not catch it.

**Failure:** `REUSE_POSITIONS=0` and any missing core positions file stops the fresh rebuild. **Minimal fix:** use `printf -v "NEED_CORE_$PT" '%s' 1`, an explicit case statement, or a proper array. Also validate all four diplomacy files when deciding whether to rebuild that arm; line 42 only tests W4. Include `fixest` and `igoR` in conditional rebuild dependency checks; they are needed by R/02 but absent from the current preflight.

This is conditional on the rebuild path; it is not evidence that the default reuse path necessarily hits this command.

### N2 — P1: the fingerprint can certify stale selection/exclusion outputs

**Location:** `run_all.sh:37–51,65–71,112,156–158`.

The fingerprint includes scripts/config/event manifest/position bytes, but not `inputs/phaseA_specs.txt`, approved exclusion-list contents, actual layer bytes, or package versions. A clean commit changing only the Phase-A IDs leaves the fingerprint unchanged and reuses `phaseA.rds`. Changing an approved history list similarly reuses `var_hist_<name>.rds` and its old sweep. Those changes need not alter cardinality, so the new summary count check will not detect them.

There is also no refusal of a nonempty legacy output root without `FINGERPRINT`: line 51 stamps the current fingerprint onto it, and existing results are then skipped. Layer output is reused before the fingerprint check. This does not validate that old output was produced by the currently hashed code or inputs.

**Minimal fix:** fingerprint each externally supplied selection/list input at the artifact that depends on it, include layer/input and environment provenance, and refuse adoption of an unversioned nonempty analysis root. Optional history approval should invalidate/rebuild its dependent variant and sweep, not require rerunning the entire grid. Verify cached RDS IDs against the currently requested IDs. Logging package versions is not cache invalidation.

### N3 — P1: failed dump processes can be accepted through shared, non-atomic staging

**Location:** `run_all.sh:16,79–82`; `R/05_dump_events.R:484–490`.

The dump writes directly to `output/rel_sweep/<name>.rds`, including in TEST or a custom OUT run. It is not saved through a temporary file plus rename. The runner ignores the process pipeline's failure with `|| true`, then accepts any nonempty file of that name and moves it into the destination.

**Failure:** an interrupted save leaves a nonempty incomplete RDS, or a previous run leaves a valid but stale staging RDS. A subsequent dump fails; the shell can nevertheless move that old/partial file into the new fingerprinted run and later skip it. Concurrent test and production jobs also share the same staging names.

**Minimal fix:** pass a run-specific destination into R/05, atomically rename a completed dump, and require successful process exit plus readable/validated metadata before accepting it. Do not infer process success from existence. This is a residual test-isolation and provenance problem despite the separate final output roots.

### N4 — P1: a missing observed statistic produces the smallest diagnostic tail value

**Location:** `R/08_resample.R:65–67,97–110`.

The new code permits valid M0 alongside an invalid observed M1, leaving observed M1 endpoints NA. If the artificial assignments produce valid M1 endpoints in ≥90% of draws, `reportq` is TRUE. Comparing those draws against NA produces only NA comparisons; `colSums(...,na.rm=TRUE)` becomes zero. The exported diagnostic is then `1/(1+n_ok)` rather than NA. With 1,000 valid draws, it is approximately **0.001 despite no valid observed statistic**.

**Minimal fix:** require `is.finite(est[QN])` as well as the B/success gate, and explicitly set every interval, positivity share and diagnostic tail value for an invalid observed fit to NA. Preserve model-specific status in the output.

**Related P2 remainder:** line 90 still returns all-NA quantities when a bootstrap draw loses one type. A pooled M0 may be perfectly estimable in that draw. Fit M0 first and suppress only M1 when the typed fit is unavailable. Otherwise its reported success denominator is still selected by type availability.

### N5 — P2: manifest/status publication is not restart-safe

**Location:** `run_all.sh:54,93–99`; `R/11_summarize.R:18`.

A sweep appends a manifest row before writing its status file. If interrupted between those actions, the next invocation reruns the task and appends the same path again. The summarizer then rejects the duplicate path. Conversely, a status file is trusted without checking the manifest entry or revalidating the CSV's status/content.

**Minimal fix:** atomically upsert task records keyed by task identity, or reconstruct the manifest from validated task records. A completed CSV can be rediscovered without rerunning the analysis. Make summary publication atomic as well; R/11:51–63 still writes its final artifacts directly. Guard line 95's display grep so an empty captured log does not bypass the intended incomplete-status branch through `set -e`.

### N6 — P2: the new summary's ID and requested-test claims are stronger than its checks

**Location:** `R/11_summarize.R:18,30–38,48–50`.

Equal `n_seen` and `n_requested` do not prove the same IDs were analyzed. Replacing one requested measure with an unintended one preserves that count. A left join to validation permits unknown plausibility values without refusing them. The manifest needs the expected ID set or a verifiable set hash, not just its size.

For Wald tests, `n_requested=.N` counts only rows actually emitted by valid fits. Example: 100 measures requested, 20 invalid fits and five singular Wald tests among the remaining 80 produces **75/80**, not 75/100. The separate census helps but does not make the column labeled requested correct. M0t is estimated but omitted from its own validity census.

**Minimal fix:** retain expected IDs and model/test applicability, check set equality, and derive requested denominators before dropping invalid/skipped fits. Report requested, fit-valid and test-computable as distinct counts. Reject or explicitly classify IDs absent from validation.

### N7 — P2 unresolved validation gap, potentially P1 if exercised: mixed-source war-link fallback

**Location:** `R/05_dump_events.R:304–311`.

The fallback prefers single-source rows but still uses collapsed `war_link` maxima if none exists. Consider a primary nonwar dispute represented only by a mixed-source row whose other dispute supplies a war number: it can be incorrectly assigned the war's initiating side rather than the primary dispute's direction rule. The participant-side correction does not identify which source supplied that war link.

**Minimal fix:** obtain the war mapping from an authoritative primary-dispute link, or assert that no fitted fallback uses an unverified mixed-source-only mapping and halt on such a case. If the asserted count is zero, this can be closed without changing current estimates. The ZIP supplies neither that diagnostic nor the RDS evidence needed to establish zero exposure; I do not claim an observed current misclassification.

## Step B, date alignment and the aggression question

The revised Step B can support a **recorded-originator, eventual-primary-dispute-fatality association**. Filtering original joiners, rather than the window-promoted role, is the right way to prevent promoted late entrants from entering this onset comparison. The intercept-plus-surface restriction is also correctly constructed for the types present.

Three limits should be made explicit before interpreting it:

1. `orig_joiner==FALSE` is a role filter, not a date check. Assert and publish the relation between retained event year, participant entry and primary-dispute start. If original participants can produce later pair entries/re-entries or composite-date events, their t−1 covariates need not be pre-dispute. Restrict those rows or construct an explicit onset decision unit. Report unique disputes, initiator-dispute decisions and target rows, not just events.
2. The new outcome is GML primary-dispute fatality, whereas the fatal-choice sample uses the project's fmid membership. Supply a cross-tab of those definitions on the same onset sample, with unknown outcomes by year/type. Otherwise Step A/Step B do not constitute an exact decomposition of the fatal-choice sample. There is no requirement that the definitions be identical for a useful sensitivity, but the mismatch must be visible.
3. Excluding joiners changes the population and removes the defender category from the default comparison. It does not adjudicate why the fatal-sample joiner pattern exists, or identify a causal selection-versus-failed-bargaining mechanism. A dispute may be fatal at its first action; without fatality timing, call this eventual severity from a pre-dispute reference, not a measured nonfatal-to-fatal transition hazard.

For onset-aligned event studies, exclude unknown/conflicting onset references from the primary onset analysis instead of silently substituting entry. Retain the entry-aligned version and a separately flagged fallback sensitivity. The new unbalanced analysis also needs n events and total weights by type **and k**; the single overall n does not reveal changing composition.

**Aggression recommendation:** default the inferential aggression variant to complete, reconciled history windows and fit its base comparator on exactly the same rows. Keep the censored-count version as a labeled lower-bound appendix sensitivity. A complete-window sensitivity is enough to assess how much the current result depends on censoring, but it does not turn the truncated full-sample count into a valid fully observed control. Do not drop censored cases from the main target-choice analysis merely because this particular control lacks coverage. Also separate unresolved onset-date fallbacks from complete histories.

## 3. Historical exclusion rule and lists

**Do not approve the ten-case punitive list unchanged.** It is a useful candidate ledger, but its inclusion rule is not yet consistently operationalized. This judgment concerns coding, not whether interventions were morally justified or legally accepted.

First settle four decisions in `docs/historical_exclusions_DRAFT.md:5–10`:

- Are the four enumerated triggers exhaustive? If so, restoring an elected government is not itself a listed trigger. If “international norm or obligation” is open-ended, debt collection and monarchical legitimacy cannot be excluded merely by declaring them non-normative; specify the boundary in advance.
- Does any qualifying public rationale suffice, or do territorial/colonial/alliance aims veto inclusion? Under an absolute mixed-motive veto, some proposed positive cases become doubtful. If public rationale suffices despite other aims, code that consistently for positive and negative cases.
- Does three states mean three **using force in the relevant operation**, or three MID participants on a side, including threats/support across different years? Those are not equivalent. Define the time/operation unit before screening.
- Screen both sides under the same rule. The current candidate screen is initiating-side coalitions plus two specially added UN cases. That does not establish a complete screen of qualifying defender-side coalitions. The lists are retrospective, history-coded robustness sets written after the influential cases were known; do not portray them as an untouched preregistered test.

| Candidate(s) | Judgment under the stated rule |
|---|---|
| **51 Korea; 3957 Gulf War** | Defensible collective-security/punitive inclusions. Korea's resolutions 83/84 recommend military assistance/unified command; Gulf resolution 678 expressly authorizes force. If “UN-authorized” means broadly Council-sanctioned military action, both fit; if it demands an express authorization formula, acknowledge the different instruments. Sources: [Resolution 83](https://mandates.un.org/mandate/S/RES/83%281950%29), [Resolution 84](https://www.refworld.org/legal/resolution/unsc/1950/en/80514), [Resolution 678](https://digitallibrary.un.org/record/102245?ln=en). These historical judgments do not independently verify the MID IDs or the claim that every affected event is a defender; check actual type/event-weight cross-tabs, especially for any-source drops. |
| **31 Boxer intervention** | Defensible under the protection-of-diplomats/nationals trigger. However, relief, occupation and indemnity/concessions coexist; specify how the colonial-aim exclusion is applied. The primary record documents claims against the Chinese government for Boxer actions, while the campaign history documents the subsequent garrisons/indemnity. Sources: [1900 diplomatic record](https://history.state.gov/historicaldocuments/frus1900/d183), [Army campaign history](https://history.army.mil/Research/Reference-Topics/Army-Campaigns/Brief-Summaries/China-Relief-Expedition/). |
| **146 Navarino** | Plausible humanitarian inclusion, but should be borderline until the rule distinguishes protection/atrocities from intervention in a national-independence conflict and government instructions from the battle itself. Contemporary parliamentary discussion links the operation to the 1827 treaty and the condition of Greek captives; later discussion distinguishes contemplated coercion from an engagement ordered by governments. Sources: [1833 parliamentary record](https://hansard.parliament.uk/commons/1833-03-26/debates/e972a6e2-8efd-41f0-9aa1-ca7c028085f9/CommonsChamber), [1877 parliamentary discussion](https://hansard.parliament.uk/Commons/1877-05-14/debates/cf79aab8-8890-4eba-9518-f57ae9c6268c/TheEasternQuestion%E2%80%94TheResolutions%28MrGladstone%29). Freeze an episode-specific rationale, not an automatic clean “yes.” |
| **3552 Iraq no-fly zones** | Plausible under the stated public-enforcement/humanitarian rationale if three states actually used force in the coded episode. That does not by itself put it in the UN-authorized subset. Record the claimed mandate separately from actual Council authorization, and verify the episode/participant dates against the MID. |
| **4273 Iraq 1997–2003** | Not a single clean coding unit. Desert Fox was a joint US–UK operation: by itself it fails the three-state rule. The later invasion is a different coalition/episode. A dispute-level exclusion may be retained if the declared rule is “any qualifying episode within this MID,” but it then removes nonqualifying episodes too. The draft's “not UN-authorized” caveat does not resolve the coalition-size problem. Source: [UK defence minister's description](https://publications.parliament.uk/pa/cm200203/cmhansrd/vo030324/text/30324w05.htm). Keep strict versus expanded episode-based lists separate. |
| **4137 Kosovo** | Defensible atrocity/humanitarian-rationale inclusion, subject to the same mixed-motive and actual force-participant criteria. A humanitarian public rationale does not imply Council authorization. Retain a contemporaneous justification and operational coalition roster in the coding ledger; the ZIP alone does not document them. |
| **4283 Afghanistan** | Plausible in a broad public self-defense/response-to-foreign-attack category, but specify state responsibility/harboring as a trigger and the date by which three states used force. Article 51/self-defense recognition is not the same instrument as express UNSC authorization of the initial invasion; later authorizations of other missions must not be transferred to it. [Resolution 1368](https://documents.un.org/api/symbol/access?l=en&s=s%2Fres%2F1368%282001%29&t=pdf) recognizes self-defense. Keep it outside the proposed narrow authorization subset unless a specific operation and mandate are identified. |
| **4251 Sierra Leone** | Not justified by the draft's listed reason alone: restoration after a coup is outside the four enumerated triggers. It may qualify after a general rule explicitly adds unconstitutional overthrow or after documenting a qualifying atrocities rationale. [Resolution 1132](https://digitallibrary.un.org/record/244598?ln=en) demands democratic restoration and authorizes ECOWAS enforcement of embargo provisions; that is not a blanket statement authorizing any regime-restoration campaign. This case also requires an explicit operation/date distinction for the UN subset. |
| **3058 Grenada** | Borderline inclusion is defensible under a public protection-of-nationals rationale; contested legality does not negate a rule based on stated justification. Specify which three states used force, rather than counting endorsement or presence alone. UN condemnation should be recorded separately: [General Assembly resolution 38/7](https://www.refworld.org/legal/resolution/unga/1983/en/6579) calls for ending the intervention and withdrawing foreign troops. |
| **89 Spain; 135 Mexico; German/Italian wars; WWI/Russia 1919/WWII** | Excluding clear succession/unification/general-war episodes is consistent with a narrowly bounded rule. Mexico is a boundary check, not an automatic negative under an open-ended “obligation” rule: the tripartite intervention was publicly directed at unpaid debts before the later French imperial project. [US diplomatic history](https://history.state.gov/milestones/1861-1865/french-intervention) distinguishes those stages. Likewise state why monarchical restoration is excluded but elected-government restoration included. I cannot certify the listed MID IDs or blanket exclusions spanning heterogeneous episodes without the underlying case ledger. |
| **Libya 2011** | If genuinely absent from the chosen fatal-event sample, exclusion from that sample-specific drop list is mechanically reasonable. Independently, it is a clear authorization case under [Resolution 1973](https://mandates.un.org/mandate/S/RES/1973%282011%29). Verify and publish its actual event/sample exclusion reason; I cannot confirm “not in the fatal sample” from the ZIP. Apply the same historical rule if it appears in another event foundation or sample. |

Recommended implementation: a source-cited ledger of all screened coalition episodes, with separate fields for rationale, force participants/dates, Council mandate, mixed aims, MID mapping and uncertainty. Freeze a strict list plus an explicitly expanded/borderline list, have another coder review without coefficient values, then map episode classifications to all-source event exclusions and disclose overcoverage when dropping a whole MID. Do not select the preferred list after seeing which preserves the headline. File presence is the runner's approval switch; this review does not approve or create those files.

## 4. Production decision

**NO-GO now.** Before the full run:

1. Fix N1, N2 and N3; make task/manifest publication restart-safe. Check a fresh rebuild and restart logic separately from the default reuse path.
2. Fix N4 and test an invalid observed M1 with a valid observed M0, plus a draw with a missing type. No quantity should acquire inference without a valid observed estimate.
3. Assert zero exposure to the unresolved mixed-source war fallback or replace it; verify the narrowed Step B onset population and outcome cross-tab; correct the summary's ID checks and requested denominators.

The historical drops can remain explicitly pending while the other analyses run; their absence is a documented gap, not proof that the rest failed. The current runner intentionally exits 2 until those inputs exist. After technical fixes, a source/data preflight and small R smoke run would still be required before a production go decision. That execution was intentionally not performed in this read-only review.
