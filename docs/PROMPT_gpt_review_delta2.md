# Third pass: blind coding of the history exclusions, then verify the fixes to your delta review

You reviewed this repository on 2026-11-20 (`docs/REVIEW_delta.md`, verdict NO-GO; findings N1–N7, the step B / event-study / aggression notes, and §3 on the history-chosen exclusions). The PI approved every fix. This pass has two parts, in this order.

## Part 1 — blind coding (do this FIRST, before opening anything else in `docs/` or `inputs/`)
Open only `docs/coding/blind_coding_sheet.md`. Code every case under the rule stated there, citing one contemporaneous source per case where you can (do not invent URLs). Do not open `docs/historical_exclusions.md`, `inputs/exclusions_*.txt`, `docs/REVIEW_delta.md` or `docs/coding/coder2_*` until your coding is written. Write it to `docs/coding/coder3_gpt.md` in the sheet's output format. You are the third coder; an independent agent is the second.

## Part 2 — verify the fixes
Then read `README.md`, `docs/decisions.md`, `docs/historical_exclusions.md` (the PI's rule, ledger and four lists), `run_all.sh`, and the changed scripts. Compare your Part 1 coding with the ledger and with `docs/coding/coder2_claude_agent.md` if present, and list every disagreement.

| Finding | Fix | Where |
|---|---|---|
| N1 rebuild branch | `printf -v` for the dynamic flags; all four diplomacy files checked; `fixest`, `igoR` in the preflight when rebuilding; every build step must exit 0 | `run_all.sh` §1 |
| N2 provenance | Fingerprint adds the network inputs + `R/01` hash, the layers file hash, and R/package versions. Sets, `random_ids.txt`, variants and the layers file carry `.stamp` files (hash of their inputs: validation CSV, id list, input dump, arguments incl. dispute lists) and are rebuilt when the stamp differs. Sweep keys use the content hashes of the dump and positions files. Cached set ids are checked against their id lists (phaseA, random, main). A non-empty output folder without the runner's `.created_by_run_all` marker is refused | `run_all.sh` (`h`, `hm`, `fresh`, `stamp`, `sub_set`, `ids_of`, `var`, `sweep`) |
| N3 dump staging | `R/05` takes `--dump_out` (run folder), writes a temp file and renames; the runner requires exit 0 and validates the frame (types, sample, window, onset rule) | `R/05_dump_events.R` (search "delta review N3"), `run_all.sh` `dump()` |
| N4 resampling | Reporting gate requires a finite observed estimate; draws that lose a type still fit M0; `obs_valid` and `reported` columns. Regression test on synthetic data: `tests/test_resample_edge.R` (passes) | `R/08_resample.R`, `tests/` |
| N5 manifest | Per-task status records written atomically after the CSV is validated; one record per CSV; manifest rebuilt from validated records before the summary; summary outputs atomic; every R call's exit status checked (`runR`) | `run_all.sh`, `R/11_summarize.R` |
| N6 summary | Manifest carries each output's expected-ids file; ids seen must equal the expected set exactly; ids absent from validation refused; Wald table reports computable / fit-valid / requested with denominators taken before dropping invalid fits; M0t census | `R/11_summarize.R`, `R/07_fit_sweep.R` census |
| N7 war-link fallback | Disputes whose war link comes only from mixed-source rows are counted and the dump stops if any exist. Count on the current data: 0 (logged for every dump) | `R/05_dump_events.R` (search "N7") |
| Step B | Only decisions in the primary dispute's GML start year; counts of dropped later-year rows, disputes, initiator×dispute decisions; cross-tab of GML fatality vs fatal-sample membership with unknowns by type and era (`stepB_crosstab_*.csv`); README calls the outcome eventual severity | `R/07_fit_sweep.R` `run_stepB`, `R/05` (`disp_styear`) |
| Event study | Onset alignment excludes missing/conflicting starts; `--onset_fallback=TRUE` sensitivity run; events and weight by type and k written | `R/10_event_study.R`, `run_all.sh` §9 |
| P2-11 aggression | `--window=complete` (default; drops censored and fallback-start strata), `complete_base` (same strata, no control), `censored` (lower bound) | `R/06_make_variants.R`, `run_all.sh` §6 |
| Docs | README choice set and provenance; pooling stated as an assumption; influence accepts one- or two-type contrasts only | `README.md`, `docs/decisions.md`, `R/09_influence.R` |
| History lists | Rule fixed (closed triggers incl. harbouring; main public justification; three states using force by MID participant data; both sides screened); four lists run as separate robustness checks | `docs/historical_exclusions.md`, `inputs/exclusions_{un,strict,strict_noborder,expanded}.txt`, `run_all.sh` §7 |

## What I want (write everything to `docs/REVIEW_delta2.md`)
1. Part 1 summary and the disagreement table (you vs ledger vs coder 2), with your recommended resolution per disagreement.
2. For each finding above: fixed / partly fixed / not fixed, with file:line evidence and, if not fixed, the failure scenario that remains.
3. Any NEW P1/P2 bug introduced by these changes (especially in `run_all.sh`: `set -e` interactions, the hash memo, stamps, status records and manifest rebuild).
4. Go / no-go for the production run.

Read-only otherwise: write only `docs/coding/coder3_gpt.md` and `docs/REVIEW_delta2.md`; do not modify, move or delete any other file, and do not run the R scripts.
