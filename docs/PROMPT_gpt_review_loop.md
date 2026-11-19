# Review request: the multiverse loop for "Coalitions and the targeting of isolated states"

You are an adversarial methodologist and R programmer who reviews for APSR, IO and JOP. You reviewed this project's pipeline on 2026-11-18/19 (findings P1-1 to P2-14, P3-15). This is a new, cleaned repository that implements the fixes and a full multiverse loop. Audit it before it runs (the run takes more than 24 hours on a laptop). Do not rewrite the project.

Read in this order: `README.md`, `docs/decisions.md`, `config.sh`, `run_all.sh`, then `R/`.

## What changed since your last review
1. **Main specification** (README table): three types (bilateral / initiating-side coalition / target-side joiner), 1-month originator window, onset coding, fatal MIDs. Windows 0–6 and the four-type model go in the appendix.
2. **Your P1 fixes**: onset rule explicit and part of every file name (P1-1); joiner side mapped to the dispute's initiating side with a 13-dispute case review and regression checks (P1-2, `R/05_dump_events.R`, search "Initiating side of each dispute"); aggression window and defense control anchored on the GML dispute start (P1-3, P1-5, `R/06_make_variants.R`); no future diplomacy waves (P1-4, `R/01_build_layers.R`).
3. **Your P2 fixes**: fully pooled M0 plus M0t, per-model validity (P2-6, P2-11); `--scale=pooled` for lag comparisons and year-t eligibility only when used (P2-7); decomposition with an earlier non-overlapping baseline, `--decomp_terms=full` with a joint test of the omitted terms, event-weighted SD, BASELINE skipped where not comparable (P2-8); plausibility rules reproduced exactly (P2-10, `R/03_validate_measures.R`, verified to 1e-15 against the 2026-11-18 CSV on 201 measures); drops by any source dispute (P2-13). Not fixed: s22 (dropped from the main text), opposing-side joiner choice sets (P2-14).
4. **New pieces**: failed-bargain step B (`R/07_fit_sweep.R --stepB=TRUE`: weighted logit of "dispute turns fatal" on the chosen target's type-specific gaps, two-way CR0); measurement error that mimics confounding (`--merror=def|contig|random`); permutation of dR within events and of type labels (`R/08_resample.R`); leave-one-dispute/initiator-out influence and data-driven top-5/top-10 dispute drops (`R/09_influence.R`); appendix summarizer (`R/11_summarize.R`).

## Smoke results on the main measure (for orientation, not to be re-litigated)
- Three types, 1-month window: chooser − bilateral, less central target, 1.00 [0.56, 1.45]; defender − chooser −0.82 [−1.50, −0.15]; all surfaces differ χ² 101 (10 df).
- Step B (all MID onsets → fatal): no type difference (χ² 11.1, 10 df, p 0.35).
- Decomposition, full terms: the omitted terms are jointly significant (χ² 26.2, 9 df, p 0.002), so the linear decomposition is a restriction the data reject.
- Leave-one-dispute-out: largest single drop in chooser − bilateral is 0.16 (dispute 4137).

## What I want
**A. Code audit**, ranked P1 / P2 / P3 with file, line, failure scenario and fix. In particular:
- the initiating-side rule and the window-promotion side fix in `R/05_dump_events.R`;
- `R/07_fit_sweep.R`: the new options (scale, decomposition, measurement error, step B), the census, option refusal;
- `R/08_resample.R`: bootstrap and both permutation schemes, the reporting gate;
- `R/09_influence.R`: contrast vectors, matching of collapsed source disputes;
- `R/03_validate_measures.R`: the power-capture metrics (outcome-blind?);
- `run_all.sh`: resumability (marker files), anything that silently skips or reuses stale output, the reuse of positions built on 2026-11-18.

**B. Design of the loop.** Is anything a reviewer will ask for missing? Is anything here uninformative or misleading as specified (for example the permutation schemes, the measurement-error construction, step B's sample and outcome, the top-dispute drops chosen from the main measure's own influence)?

**C. Reporting.** How should the multiverse be summarised in the appendix so it cannot be read as "99% of specifications" evidence (correlated measures, families, plausible vs implausible, windows)?

## Output
Write the complete review to `docs/REVIEW_loop.md`. Section A as a ranked list; B and C as terse bullets, each ending with the single most useful change. No praise. Read-only otherwise: do not modify, move or delete any other file, and do not run the R scripts.
