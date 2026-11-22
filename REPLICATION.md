# Replication

One repository, one commit, one output folder per loop. Every output folder records its own provenance (`FINGERPRINT`, `NETSHA`, the commit and R/package versions in `logs/run_*.log`, and `status/manifest_sweeps.tsv`); the summarizer reads only manifest-listed outputs.

## Layout (Jared's MacBook)
| Folder | Contents |
|---|---|
| `~/targeting_multiverse_v2/targeting_multiverse/` | This repository at the GPT-approved commit (2f4e903 code; later commits docs/tests only) |
| `…/output_init/` | Stage 1, main analysis: `TYPES=init` (bilateral initiators vs coalition originators) |
| `…/output_bcd/` | Appendix: `TYPES=bcd` (bilateral / attacking-side coalition / target-side joiners) |
| `~/targeting_multiverse_archive/` | Earlier runs on the 7,686-measure grid, read-only baselines: Nov 20 bcd (f60c55c), Nov 21 init (230d585) |

## Inputs
- Event data: `DATA_ENTRY` in `config.sh` (built by `v2/R/01_build_data.R --entry=later` in the main repo `strategic_targeting`, branch `clogit-roles`); fingerprinted through `data_entry/manifest.csv`.
- Network files in `NETDATA`: `dipl_dyad.csv`, `destaRaw.rda` (DESTA), `dca_dyad.csv`; everything else from `peacesciencer` 1.2.0 (COW IGO, ATOP, UN voting similarity, contiguity, capital distance, NMC, system membership) and `igoR` (IGO arm).
- R 4.5.3; package versions are logged at the top of every run log.

## Commands
```bash
cd ~/targeting_multiverse_v2/targeting_multiverse
git log --oneline -1            # the approved commit
git status --short              # must be empty (production refuses uncommitted code)
Rscript tests/test_sry.R        # 23 checks, ALL PASS
bash tests/test_hash_cache.sh   # runner hash-cache regression
Rscript tests/test_resample_edge.R
TYPES=init OUT=output_init caffeinate -i bash run_all.sh    # stage 1 (main)
TYPES=bcd  OUT=output_bcd  caffeinate -i bash run_all.sh    # appendix (run after stage 1 finishes)
```
Each loop rebuilds all relational positions from its own layers (`REUSE_POSITIONS=0`, enforced; about 4 h), then runs the multiverse (12,294 measures) and the robustness grid. Restarts resume; finished steps are skipped. `LOOP_DONE` ends a complete run; `LOOP_FINISHED_WITH_GAPS` lists gaps in `status/incomplete.txt`.

## Where results live
| Output | Folder |
|---|---|
| Measure validation (plausibility, CINC, major powers, stability) | `validation/` |
| Multiverse per sample × window | `sweep/multiverse_*` |
| Robustness checks | `sweep/rob_*` |
| Influence, bootstraps, permutations | `influence/`, `resample/` |
| Event studies | `eventstudy/` |
| Appendix tables | `summary/` (`summary_main.md`, `by_family.csv`, `estimates_long.csv`, `census.csv`, `wald.csv`) |

## Review trail
GPT reviews and Claude's fixes are in `docs/REVIEW_*.md` with their prompts in `docs/PROMPT_*.md`; decisions in `docs/decisions.md`, `docs/historical_exclusions.md`, `docs/measure_alignment_2026-11-21.md`.
