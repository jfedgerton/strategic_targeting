# Paper pipeline (run_paper.sh), design as fixed by the PI on 2026-11-24

## Main measure
- **BOOKraw_ptadesta|igo+def+off+pta+dca|pr**: the book's unweighted measure (book tag `c9_igo.def.off.dca.pta_raw`). Five layers, each in [0,1] with true zeros:
  IGO = shared full COW IGO memberships (peacesciencer dyadic count), min-max within year; ATOP defense; ATOP offense; Kinne DCAD dcaAnyV2 (2010 carried to 2011-14);
  PTA = book `layer_pta_desta` (every member pair of an in-force DESTA treaty linked, only actual withdrawers leave, ISO->COW with Germany 260 for 1955-90,
  Czechoslovakia 315 through 1992, North Yemen 678 through 1989; 196,196 linked undirected dyad-years). Tie = mean of the 5 layers / its yearly maximum.
  Position = log PageRank (igraph, damping 0.85, undirected, weighted; every state in the year's key is a vertex).
- With the book's own PTA file the builders reproduce the book's `_raw` and `_sry_contig12` edges to 1e-14 (checked 2026-11-24).
- The weighted book measure (BOOKc12_ptadesta, sry contig12) is an appendix robustness measure.

## Models (unchanged code, R/lib/s2_cell.R)
- Stage 1: conditional logit of target choice by initiator type (bilateral vs coalition originator); between targets (no target FE, three-way clustered) is the main stage-1 comparison; target FE (within) feeds the joint model; + initiator-target tie as robustness.
- Stage 2 (PI 2026-11-24): categorical joining model (join A / join B / originator O vs out); each option has ties to the side it joins AND to the side it would fight (S2_REF_DESIGN, as in the Nov 21 design memo). Own-side-only ties (S2_REF_DESIGN_OWN) are robustness column (3c).
- Stage 2: categorical joining (A attacker side / B target side / O originator vs out) on the reference cohort.
- Joint: two-stage FE (stage-1 target FE -> bbar per dispute -> stage 2 with theta on bbar); two-stage dispute Bayesian bootstrap (Exp(1) dispute weights scaled to mean 1, Rubin 1981; R/lib/s2_cell.R boot_weights), common draws, seed 123.
- AME: dropped (PI 2026-11-24); code kept, step 8 off by default.

## Main table (R/main_table.R)
(1) stage 1 between targets; (2) stage 2 alone (no target trait); (3) joint (stage 1 within targets + stage 2 with theta). Same bootstrap draws for all three. 95% percentile intervals; >= 80% finite draws per quantity.

## Network-measure grid (R/grid_layers.R, R/lib/grid.R, R/grid_build.R): 236 measures, all PageRank (PTA depth dropped from the sweep, PI 2026-11-24)
| Block | Content | Specs |
|---|---|---|
| raw | 24 layer combos of {igo, def, off, dca, pta} containing igo or def; IGO count / cosine / hyperedge (if igo); PTA book | 56 |
| dca1 | as raw, DCA = dcaAnyV1, combos containing dca only | 28 |
| ally | as raw + one "other alliances" layer (neutrality / non-aggression / consultation) | 56 |
| sryc12, srydist, syc12, sydist | weighted: observed - expected per layer (logit / PPML), FE sender + receiver + year (sry) or sender x year + receiver x year (sy), geography contig12 or log capital distance; positive part / yearly max; IGO count, PTA book | 4 x 24 |
Swept at the main window on fatal MIDs and on all MIDs (472 cells). IGO cos / hyper follow book 01_data/06_igo_variants.R (state-year file); The PTA depth layer (book 01_data/04_pta_desta_costly.R) is still built in grid_layers.rds but is no longer swept.

## Bootstrapped cells (R/grid_pick_boot_cells.R, fixed before results)
Main x windows 0-6 x {fmid, mid} (the main-window pair is the main table's B_MAIN bootstrap, step 3; the other 12 here); main layers under IGO cos, hyper; drop-one layer; the 4 weighting schemes; DCA v1; + other alliances. 25 cells here + 2 in the main table. Window 0 has no originator alternative, so its O quantities are NA.

## Checks that must stop the run
PTA 196,196 and 93 withdrawers; DCA v2 11,031 (1980-2010); raw pooled tie share 0.952; grid raw main = BOOKraw to 1e-10; grid sry c12 main = BOOKc12 (tested); fixest warnings/convergence; sry dropped rows must have y = 0.

## Missing-data robustness (PI 2026-11-24)
- Main design keeps the t-1 timing rule. Rows lost to it (states not in the system at t-1) and to the bbar rule (no repeat target) are counted in `main_table_missingness.csv`.
- Main table column (2) is fitted on the joint model's disputes, so (2) vs (3) differs only by the target trait.
- `main_table_robust_*.csv`: (2b) stage 2 alone on all disputes; (3b) joint with b = 0 for one-time / unchosen targets (first N_ROBUST draws).
- New-state cell: `s2_01_build_riskset.R --newstate_t=1` back-fills year-t values to t-1 only for state-years absent at t-1 (positions, CINC, democracy, dyads, edges; scales from original data), stage 2 only; one bootstrapped cell (B_CELL). Sandbox: 220 back-filled state-years, 560 disputes / 65,274 rows (vs 544 / 63,143).

## Reporting fixes (audit part B, 2026-11-24)
- Table 1 panel A: weighted moments (stage-1 event weights) next to unweighted; events_weighted. F1: West Germany (260) in 1990. Axis labels state scaled units. F3: denominators are candidates for whom the option is available. Figure captions report the share of points outside the 1st/99th-percentile axis cut.
- Grid summary: decision-regression SEs clustered by measure; main measure drawn on top of the specification curves. The 16 sy specs without IGO are SKIPPED (constant positions), so 220 measures per sample are fitted.
- Run-stopping: fixest warnings in grid_layers.R; grid sry contig12 main = BOOKc12 (grid_build.R --c12).
- Marginal effects: stage 2 alone on the joint sample; columns trimmed to the design; draws saved in blocks of 50 (resumable); dump/build pairing check; a draw whose fitted term set differs from the point fit is dropped.
