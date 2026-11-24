# Book replication pipeline

This folder rebuilds the book's network, communities and H1–H3 results from `data/raw/`. Run it with:

```bash
./run_all.sh                 # everything, about an hour
STAGES="5 6" ./run_all.sh    # rerun only the models and validation
```

Seed 123 throughout. Logs go to `logs/`, tables to `results/`.

## data/raw

| Folder | Contents |
|---|---|
| `tremors_panel/` | Tremors dyad-year layer panel (ATOP 5.1, COW IGO v3, COW trade, Kinne DCA v1.0, via peacesciencer) |
| `cow_igo/` | COW IGO v3 state-year memberships (checksums in `SOURCES.md`) |
| `desta/` | DESTA 2.3 dyads, withdrawals, depth indices and treaties |
| `mia/` | MIA v4 IO authority |
| `peacesciencer/` | State-year frame and NMC/CINC |
| `jop_regime/` | V-Dem v12 core |
| `prebuilt/` | Inputs that can't be rebuilt in this environment (see below) |

**Contents of `prebuilt/`**

| File | What it is | Why it is prebuilt |
|---|---|---|
| `pta_layers.rds` | Old JOP DESTA PTA layer (`layer_pta_desta`, the main PTA) | Its build script isn't in this pipeline |
| `atop_un_layers.rds` | ATOP obligation types and UN kappa | Building it needs the peacesciencer extdata download, and the network here blocks it |
| `state_year_covariates.rds` | Polity, V-Dem, CINC, GDP and civil war | `02_state_covariates.R` needs `p5v2018.csv`, which isn't in `data/raw` |

## Stages

| Stage | What it does | Main outputs |
|---|---|---|
| 1 data | `layers_dyad_year.rds`; IGO variants (count/log/cos/hyper); hyperedge IGO; DESTA 2.3 PTA codings and variants; MIA IGO layers | `data/processed/` |
| 2 measures | Four builds of the five-layer JOP measure (IGO count, def, off, DCA, old PTA). Weighted = sry + contig12, DCA carried to 2014. Weighted truncated = same, ending in 2010. Unweighted = raw. Unweighted truncated = raw, ending in 2010 | `multiverse/c9_igo.def.off.dca.pta_{sry_contig12,raw}[__trunc2010]/` |
| 3 communities | Multislice Leiden on weighted `tie_pos`, with resolution 1, interslice weight 1, and zero-tie state-years as singletons | `data/processed/communities_multislice_w1.csv` |
| 4 frames | H1 dyad-year conflict frame; diffusion state-year frames | `robust_conflict_frame.rds`, `state_year_institutions_allyigo.rds` |
| 5 models | H1 (FE LPM and logit, with and without lagged fatal MID, dyad-robust SEs); H2 churn (unweighted vs weighted, multislice community FE); H3 diffusion hazard (unweighted and weighted; no-shock, stable and combined samples) | `results/` |
| 6 validate | Historical-bloc benchmarks for the communities | `results/validation/` |

The construction rules and every sweep option are documented in `book_network_decision_rules.md`.
