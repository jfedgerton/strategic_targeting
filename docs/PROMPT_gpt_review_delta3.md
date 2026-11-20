# Fourth pass: verify D1–D6 and the intensity samples (narrow)

You reviewed this repository on 2026-11-20 (`docs/REVIEW_delta2.md`, NO-GO; findings D1–D7). The PI approved D1–D6. D7 (re-screening the history candidates from the participant table) and further coder adjudication were declined: the history lists are reported as side-by-side sensitivity checks, not a highlighted part of the paper. Do not re-review the coding.

Verify only the changes below; do not re-review unchanged code.

| Item | Change | Where |
|---|---|---|
| D1 | `status/current_tasks.txt` records every sweep task planned in this invocation (key = group + content hashes of dump and positions + options); the manifest is rebuilt only from those keys' successful records, so a superseded result cannot reach the summary | `run_all.sh` `sweep()`, §10 |
| D2 | No associative arrays; hash cache in files keyed by path + size + mtime (`date -r`, `wc -c`), so the runner works on macOS bash 3.2; no other bash-4-only constructs | `run_all.sh` `h()` |
| D3 | `.path` files removed; the single atomic record is the source of truth; on success, older records pointing to the same CSV are removed | `run_all.sh` `sweep()` |
| D4 | `ids_of()` returns non-zero unless the ids file exists with the current stamp; both call sites check it | `run_all.sh` |
| D5 | Each positions arm stamped with layers hash + `R/02` hash + environment; an arm with a missing or different stamp is rebuilt | `run_all.sh` `arm()` |
| D6 | Reporting gate moved to `R/lib/report_gate.R`, sourced by `R/08`; direct test (NA observed estimate with fully usable draws is not reported) plus a mutation check showing the old gate fails it | `tests/test_resample_edge.R` |
| Intensity samples (new) | Jointly fatal (both states have fatalities in the dispute; state-level, `fatal1`/`fatal2`) and dispute fatality level ≥ 2..6 (`fatality_disp`), unknown dropped and counted; built per window from the fatal-MID frame and run over every measure × every window | `R/05` (`fat_level`, `jfatal`), `R/06 --variant=sample`, `run_all.sh` `varin()`, `config.sh INTENSITY` |
| Lists | Two more sensitivity lists: `noruleA` (strict + 4339, 611), `anycoder` | `inputs/`, `run_all.sh` §7 |

Executed here (TEST=1, Linux bash 5.2): fresh run LOOP_DONE; restart 11 s, no refits; strict list changed + replacement sweep forced to fail (Rscript shim) → exit 2, strict absent from manifest and summary, old CSV kept on disk; restart without shim → new strict result in manifest; list restored → old result reused, no refits.

## What I want (write to `docs/REVIEW_delta3.md`)
1. Each item: fixed / partly fixed / not fixed, with file:line evidence.
2. Any NEW P1/P2 bug, especially bash 3.2 compatibility of `run_all.sh` (please check against bash 3.2 semantics), `set -e` interactions, and the intensity-sample construction.
3. Go / no-go for the production run.

Read-only otherwise; do not run the R scripts.
