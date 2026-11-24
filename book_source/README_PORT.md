# Book replication sources ported into the targeting project (2026-11-24)

Source: `~/Dropbox/book_project_conflict_cycles/replication/` (book "The System Makes War and War Makes the System").

- `book_source/` (tracked): the book's data and builder scripts (`code/01_data`, `code/04_multiverse`), `run_all.sh`, README, PLAN, decision rules. Copied verbatim; not run from here.
- `data_entry/book/` (git-ignored, like all of data_entry): raw inputs. md5s in `book_source/DATA_MANIFEST.csv`. On the Mac copy them from the Dropbox path in the manifest.
- The large Tremors panel (`tremors_panel/dyad_year_network_layers.rds`, 258 MB) was not copied; the targeting layer panel (`output/layers_dy.rds`) is built from the same sources and matches the book's panel (918,315 dyad-years; IGO count r = 1 within year).

## PTA layer
The book's main PTA layer (`prebuilt/pta_layers.rds`, `layer_pta_desta`, 188,499 linked undirected dyad-years) is a prebuilt input in the book too: its build script is not in the book pipeline.
Independent rebuild check (`book_source/pta_reconcile_desta23.R`):

| Rebuild | Linked | Both | Book only | Rebuild only |
|---|---|---|---|---|
| old destaRaw, clique, in force | 186,370 | 176,871 | 11,628 | 9,499 |
| DESTA 2.3 dyadic, in force, withdrawals | 99,003 | 90,861 | 97,638 | 8,142 |
| DESTA 2.3 clique, in force | 197,244 | 188,053 | 446 | 9,191 |
| DESTA 2.3 clique, PTN ends 1980, GSTP ends 2006 | 188,936 | 188,027 | 472 | 909 |

The last row's two end years are read off the book layer (no recorded DESTA withdrawal explains them); they are not documented rules.
