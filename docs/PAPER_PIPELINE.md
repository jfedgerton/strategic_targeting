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
- Stage 2: categorical joining (A attacker side / B target side / O originator vs out) on the reference cohort.
- Joint: two-stage FE (stage-1 target FE -> bbar per dispute -> stage 2 with theta on bbar); two-stage dispute Bayesian bootstrap (Exp(1) dispute weights scaled to mean 1, Rubin 1981; R/lib/s2_cell.R boot_weights), common draws, seed 123.
- AME: stan/ame_joint.stan + R/ame/ame_cell.R (cells main K=0, cut posterior, K=1), bundle built in step 8, fits on Roar.

## Main table (R/main_table.R)
(1) stage 1 between targets; (2) stage 2 alone (no target trait); (3) joint (stage 1 within targets + stage 2 with theta). Same bootstrap draws for all three. 95% percentile intervals; >= 80% finite draws per quantity.

## Network-measure grid (R/grid_layers.R, R/lib/grid.R, R/grid_build.R): 306 measures, all PageRank
| Block | Content | Specs |
|---|---|---|
| raw | 24 layer combos of {igo, def, off, dca, pta} containing igo or def; IGO count / cosine / hyperedge (if igo); PTA book / depth (if pta) | 84 |
| dca1 | as raw, DCA = dcaAnyV1, combos containing dca only | 42 |
| ally | as raw + one "other alliances" layer (neutrality / non-aggression / consultation) | 84 |
| sryc12, srydist, syc12, sydist | weighted: observed - expected per layer (logit / PPML), FE sender + receiver + year (sry) or sender x year + receiver x year (sy), geography contig12 or log capital distance; positive part / yearly max; IGO count, PTA book | 4 x 24 |
Swept at the main window on fatal MIDs and on all MIDs (612 cells). IGO cos / hyper follow book 01_data/06_igo_variants.R (state-year file); PTA depth follows book 01_data/04_pta_desta_costly.R (DESTA 2.3 dyadic, max depth / 7) plus the predecessor mapping.

## Bootstrapped cells (R/grid_pick_boot_cells.R, fixed before results)
Main x windows 0-6 x {fmid, mid} (the main-window pair is the main table's B_MAIN bootstrap, step 3; the other 12 here); main layers under IGO cos, hyper, PTA depth; drop-one layer; the 4 weighting schemes; DCA v1; + other alliances. 26 cells here + 2 in the main table. Window 0 has no originator alternative, so its O quantities are NA.

## Checks that must stop the run
PTA 196,196 and 93 withdrawers; DCA v2 11,031 (1980-2010); raw pooled tie share 0.952; grid raw main = BOOKraw to 1e-10; grid sry c12 main = BOOKc12 (tested); fixest warnings/convergence; sry dropped rows must have y = 0.
