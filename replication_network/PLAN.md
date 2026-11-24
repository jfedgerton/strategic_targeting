# replication_network/ -- independent build of the main network (book / JOP 2024 layer combo)
Status: plan only (2026-11-24). Nothing here is used by the pipeline until the PI approves the switch.

## Target
Main combo igo + def + off + dca + pta, two weightings, PageRank positions, in the targeting positions/edges format:
- BOOKc12 (main): sry, contig12 (COW conttype 1-2), logit (binary) / PPML (IGO), sender + receiver + year FE, observed - expected, mean of 5, symmetrized, tie_pos.
- BOOKraw: raw layers, mean of 5, divided by the yearly max.
Each with the 5 drop-one sets. Seed 123.

## Layers
| Layer | Source (data_entry/book) | Coding |
|---|---|---|
| igo | cow_igo/state_year_formatv3.csv | shared full memberships (code 1), min-max within year; check r = 1 vs layers_dy igo |
| def, off | layers_dy (ATOP 5.1); check vs prebuilt/atop_un_layers.rds | 0/1, either direction; off ends 1990 |
| dca | layers_dy (Kinne DCAD dcaAnyV2) | 0/1; 2010 carried to 2011-14 |
| pta | desta/dyads.csv + withdrawal.csv (DESTA 2.3) | 0/1; every member pair of an in-force base treaty linked (ACP-ACP via Lome/Cotonou); member join = earliest entry-into-force year in that base treaty (accessions included); exit = recorded withdrawal; runs through 2014. PTN and GSTP are NOT ended early (PI 2026-11-24: the book layer's PTN-1980 / GSTP-2006 ends are treated as old data or a coding bug) |

## Scripts
1. code/01_pta_desta23_clique.R  -> pta_clique.rds
2. code/02_igo_count.R           -> igo_count.rds
3. code/03_build_main.R          -> positions_core_BOOKrep.rds, edges_book_rep.rds

## Checks (fail loudly)
- PTA vs book prebuilt layer: expect both 188,053; book-only 446; new-only 9,191, located in PTN after 1980 and GSTP after 2006 (anything else stops the build).
- IGO count equals layers_dy igo (r = 1 within year); raw pooled tie share 0.952.
- Builder equivalence: with the book's prebuilt PTA swapped in, 03 reproduces R/book_measure_c12.R positions to 1e-10.
## Then
- Stage 1, stage 2, bootstrap FE (500 draws) on both; compare to the prebuilt-PTA versions.
- GPT audit of the three scripts before the main spec switches; stage-2 sweep and AME bundle rebuilt after.
