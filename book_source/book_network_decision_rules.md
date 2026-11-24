# Book network: decision rules for the multiverse sweep (sections 1–6)

Handoff for another agent. Book: "The System Makes War and War Makes the System" (Jared Edgerton). Updated 2026-11-24.

## 0. Where things live

- **Book root.** `~/Dropbox/book_project_conflict_cycles` (`$BOOK`). All paths below are relative to it.
- **Multiverse output root.** `$MV`. Each specification is written to `$MV/<TAG>/data/processed/`.
- **Sync warning.** Several scripts below were created or changed in a Claude sandbox on Nov 21–24, 2026. These are `01_data/03`–`06`, the `04b`/`05b` variants, the builder switches `PTA_VARIANT`, `IGO=mia2/mia2nd/var`, `IGO_WT`, `MIA2_CAP`, `MIA2_NOREG` and `WEIGHT=igofe`, and `02e_build_communities_only.sh` suffixes. They may not be on Dropbox or GitHub yet. Check the dates first. The attached zip holds the sandbox versions.
- **Seed.** 123 everywhere.
- **R packages.** data.table, fixest (feglm/fepois), igraph, countrycode 1.9.0, peacesciencer 1.2.0 with extdata.

### Build pipeline

1. `code/01_data/01_import_layers.R` builds `data/processed/layers_dyad_year.rds`.
   - Input is the Tremors layer panel `data/raw/tremors_panel/dyad_year_network_layers.rds`, made by `SOURCE_01_build_network_layers.R` from peacesciencer, ATOP 5.1, COW IGO v3, COW trade and Kinne DCA v1.0.
   - Frame: COW state dyad-years 1816–2014, directed (both orders).
2. Layer side files, one per script in `code/01_data/` (sections 1–3).
3. `code/04_multiverse/01_build_measure.R` is the builder. It is driven by environment variables and turns a layer set, IGO mode, weighting and geography term into ties (sections 4–6).
4. `code/04_multiverse/02e_build_communities_only.sh LAYER_SET WEIGHT [GEO]` is a wrapper. It calls the builder with `LAYERS=custom` and the IGO mode taken from `IGO_MODE`, then runs static yearly Leiden. Tag format:
   `c9_<layers joined by .>_<WEIGHT><igo suffix>[_<GEO>][_pv<PTA_VARIANT>]`
   Example: `c9_igo.def.off.dca.pta_sry_contig12`.

### Builder outputs (per tag)

- `main_measure_allyigo.rds`: ccode1 < ccode2, year, layer columns, n_layers, tie_signed, tie_norm, tie_pos, plus MID, CINC, contiguity and major-power columns.
- `main_measure_edges_allyigo.csv.gz`: ccode1, ccode2, year, tie_signed, tie_norm, tie_pos.
- The "allyigo" names are historical. Keep them so the chapter scripts run unchanged.

## 1. Layers

Choose layers with `LAYERS=custom LAYER_SET=<+-separated list>`. The list must include `igo` or `def`, and `igo` must be present exactly when IGO is not `none`.

| Layer code | Definition | Source / script | Years with ties |
|---|---|---|---|
| `igo` | Shared COW IGO memberships. The coding is chosen in section 2 | `layers_dyad_year.rds` (`layer_igo`) | 1816–2014 |
| `def` | ATOP defence obligation, 0/1. An obligation in either direction counts for the pair | `03_atop_un_layers.R` (peacesciencer::atop_alliance); equals `layer_defensive` | 1816–2014 |
| `off` | ATOP offence obligation, 0/1, either direction | `03_atop_un_layers.R`; equals `layer_offensive` | 1816–1990 |
| `dca` | Kinne DCAD `dcaAnyV2`, 0/1. Section 4 covers the carry-forward | `layers_dyad_year.rds` | 1980–2010 (2014 with carry) |
| `pta` | Old JOP DESTA layer, 0/1. Section 3 has the other codings | `data/processed/pta_layers.rds` (`layer_pta_desta`) | 1951–2014 |
| `neu`, `nonagg`, `consul` | ATOP neutrality, non-aggression and consultation, 0/1 | `03_atop_un_layers.R` | Built, not in main |
| `un` | UN voting kappa (Häge FPSIM `kappavv`), signed, 1946+. Non-member pairs are NA and get a score of 0 | `03_atop_un_layers.R` | Built, not in main |
| `trade` | smoothflow1 + smoothflow2, min-max scaled within year, 1870+ | `layers_dyad_year.rds` | Built, not in main |

- **Combined alliance.** `LAYERS=alliance` uses one alliance layer, max(def, off). The main measure keeps def and off separate.
- **Sweep.** All 24 non-empty subsets of {igo, def, off, dca, pta} that contain igo or def. The **main measure uses all five**.

## 2. IGO layer coding

Set with `IGO=` (the wrapper reads it from `IGO_MODE`). Every option ends up min-max scaled within year, unless noted otherwise.

| Option | Definition | How it is made | Switch | Status |
|---|---|---|---|---|
| Count | Σ_k B_ik B_jk, where B = full membership (COW code 1); min-max scaled within year | `layer_igo` in `layers_dyad_year.rds` | `IGO=dyadic` | **Main** |
| Log | log(1 + count), scaled | `01_data/06_igo_variants.R` → `igo_variants.rds` | `IGO=var IGO_VAR=log` | Tested |
| Cosine | count / √(d_i d_j), where d = number of IGOs a state belongs to | Same | `IGO_VAR=cos` | Tested; weak on H1 and validation |
| Hyperedge | Σ_k B_ik B_jk / (n_k − 1), where n_k = full members; IGOs with fewer than 2 members dropped | `06_igo_variants.R` (`IGO_VAR=hyper`) or `04_multiverse/00_igo_hyperedge_projection.R` → `igo_hyper_size.rds` (`IGO=hyper`) | | Tested |
| MIA v1 | Σ_j authority_jt · B_ij B_kj / (n_jt − 1). Authority = (delegation + pooling)/2 from MIA v4. The 1950 score is carried back to founding. Non-MIA IGOs = 0 | `01_data/05_igo_mia_backfill.R` → `igo_mia.rds` | `IGO=mia` | Superseded |
| MIA fixed | Same projection with five changes: (a) predecessor bodies mapped to their successor, e.g. ECSC/EEC/EURATOM→EU and OEEC→OECD, membership = union; (b) COW −9 cells set by last observation carried forward; (c) Germany 1990 fix; (d) backfill capped at 10 years; (e) universal IGOs dropped (members above 2/3 of states). Two weights: `disc` = Σ auth/(n − 1) and `nodisc` = Σ auth | `01_data/05b_igo_mia_fixed.R` → `igo_mia2[_cap999][_noreg].rds` | `IGO=mia2` (disc) / `IGO=mia2nd` | Tested |
| MIA backfill | Cap of 10 years or full backfill to founding | 05b with env `CAP=999` | `MIA2_CAP=10/999` | Tested |
| Regional IGOs | Drop IGOs with ≥ 80% of members from one COW region, plus the Arab League | 05b with env `NOREG=1` | `MIA2_NOREG=1` | Tested |
| IGO down-weight | Multiply the IGO layer by w in the composite | Builder | `IGO_WT=0.5` | Rejected |
| No IGO | Leave `igo` out of `LAYER_SET` | | `IGO=none` | Robustness |

**Facts about the count layer**
- It is positive in 95.0% of pooled dyad-years, and in 99.7% after 1980.
- The figure of about 70% that appeared earlier is the average of yearly shares. The 19th-century years pull it down.
- Rebuilding the count from COW `state_year_formatv3.csv` reproduces `layer_igo` exactly (r = 1 within year).

## 3. PTA layer coding

| Option | Definition | Script / file | Switch | Status |
|---|---|---|---|---|
| Old JOP DESTA | Binary. In force only. **Every member pair of a treaty** is linked, so ACP–ACP pairs are linked through Lomé/Cotonou. 188,499 linked undirected dyad-years, 1951–2014 (`pta_layers.rds` is directed, giving 376,998 rows) | `data/processed/pta_layers.rds` (`layer_pta_desta`). The script that builds this file is not in the sandbox copy; it should be in the Dropbox pipeline | `LAYER_SET` contains `pta` | **Main** (both weighted and unweighted) |
| DESTA 2.3 in force | Binary. Counts from the entry-into-force year until the year before a recorded withdrawal of either member. Agreements with no entry-into-force year count as not in force. Country codes converted from ISO to COW; Germany 276→260 up to 1990, 255 after | `01_data/04_pta_desta_costly.R` → `pta_costly.rds` (DESTA 2.3 files: `desta_list_of_treaties_02_03_dyads.csv`, `desta_dyadic_withdrawal_02_03.csv`, `desta_indices_version_02_03.csv`) | `ptai` | Tested |
| Depth | Maximum DESTA `depth_index` over in-force agreements, divided by 7 | Same | `ptad` | Tested |
| Deep | Maximum depth ≥ 4 (DESTA third quartile), 0/1 | Same | `ptadeep` | Tested |
| e2 / e4 | Depth weight only for agreements with depth ≥ 2 / ≥ 4 | `01_data/04b_pta_variants.R` → `pta_variants.rds` | `PTA_VARIANT=e2/e4` (replaces `ptad`) | Tested |
| f1 | Drop non-reciprocal schemes: the EEC–ACP chain (Yaoundé, Arusha, Lomé I–IV, Cotonou) and GSTP/PTN | Same | `PTA_VARIANT=f1` | Tested |
| f2 | Σ over in-force base treaties of (depth/7)/(n_members,t − 1), divided by the yearly max | Same | `PTA_VARIANT=f2` | Picked briefly, then reverted |
| Mansfield–Pevehouse | Alternative PTA source | `pta_layers.rds` (`layer_pta_mansfield`) | | Earlier test; worse |

## 4. DCA timing

| Option | Rule | Switch | Used for |
|---|---|---|---|
| Carry forward | DCAD ends in 2010. Each dyad's 2010 status is copied to 2011–14 | `DCA_MODE=carry` (default) | Communities, H1, H3 |
| Truncate | Panel ends in 2010, with no carry-forward | `DCA_MODE=trunc2010` (tag suffix `__trunc2010`) | H2 churn, so carried ties can't create artificial zero churn |

**True zeros.** Years before an instrument exists, or after it lapses, are coded 0 rather than missing. This applies to offensive alliances, DCA, PTA and IGO. Every dyad-year therefore has all of its layers.

## 5. Weighting

Set with `WEIGHT=`. Each layer v is modelled separately, using only the years in which v varies. Years with no ties anywhere get a score of 0.

| Option | Layer score | Expected-tie model | Status |
|---|---|---|---|
| `sry` | y − ŷ (positive or negative) | Binary layers: `feglm(v ~ GEO \| ccode1 + ccode2 + year, logit)`. Continuous layers (IGO): `fepois` (PPML). Observations dropped as perfectly predicted must have y = 0, and get ŷ = 0 | **Main for communities** |
| `sy` | y − ŷ | FE `ccode1^year + ccode2^year`. Perfectly predicted observations get ŷ = y, so score 0 | Tested |
| `hybrid` | Binary layers: 1 if present, −p̂ if absent. Continuous: y if y ≥ ŷ, otherwise y − ŷ | Same as sry | Tested |
| `raw` | y as observed | None. No negative ties | **Main for H2 and H3** |
| `igofe` | IGO only: max(y − ŷ, 0) from `fepois(igo ~ GEO \| ccode1 + ccode2 + year)`, scaled to 0–1 by the yearly max. Other layers raw | Partial | Tested; worse |
| `w2contig` | Per layer (IGO and PTA as log(1 + count)): residual of `v ~ contig_land + year` (pooled OLS), set to 0 below 0, divided by the layer-year max. Composite = √(Σ e²), rescaled 0–1 within year | OLS | Tested |

**Geography term in the expected-tie model** (`GEO=`)
- `dist` (default): log(1 + capital distance).
- `contig12`: binary contiguity, COW conttype 1 or 2 (land, or water up to 12 miles). **This is the main choice.**

## 6. Composite, symmetrization and sign

Builder steps after scoring:

1. **Composite.** tie_dir = mean of the layer scores over non-missing layers (n_layers). This is equal weighting, Manhattan-style. With true zeros, n_layers equals the number of layers in the set. The exception is `w2contig`, which uses the root sum of squares from section 5.
2. **Symmetrize.** tie_signed = (tie_dir(i,j) + tie_dir(j,i)) / 2. Keep ccode1 < ccode2.
3. **tie_norm.** tie_signed divided by its yearly maximum absolute value.
4. **tie_pos.** max(tie_signed, 0), divided by its yearly maximum, so 0–1 within year.

**Sign choices**
- **Positive** (`tie_pos`): main.
- **Signed** (`tie_signed`): negative values mean less tie than expected and act as repulsion.
  - In static signed Leiden, the positive layer uses RBConfiguration at resolution 1. The negative layer uses CPM at resolution 0, weight −1, with no null model.
  - Script: `code/02_network_universes/02_leiden_communities.py`.
- **Negative rivalry layer:** rejected.

### Current main specifications

| Purpose | Tag | Settings |
|---|---|---|
| Communities | `c9_igo.def.off.dca.pta_sry_contig12` | `LAYER_SET=igo+def+off+dca+pta IGO=dyadic WEIGHT=sry GEO=contig12`, DCA carry. Communities from multislice Leiden on `tie_pos` (`code/02_network_universes/02c_leiden_temporal.py`, `GAMMA=1 OMEGA=1`): resolution 1, each state tied to itself in the next year at weight 1, CPM resolution 0 on the interslice layer, zero-tie state-years made singletons, seed 123 |
| H2 and H3 ties | `c9_igo.def.off.dca.pta_raw` | Same five layers with `WEIGHT=raw`, `tie_pos`. Uses `__trunc2010` for churn |

### Suggested core grid

| Dimension | Options | Count |
|---|---|---|
| Layer set | Subsets of {igo, def, off, dca, pta} containing igo or def | 24 |
| IGO coding | count, log, cos, hyper | 4 |
| Weighting | sry, raw | 2 |
| GEO | contig12, dist | 2 |
| Sign | positive, signed | 2 |
| Community method | multislice weight 1, multislice weight 0.5, DynMux filtered, static yearly | 4 |

That is 3,072 specifications.

**Notes for the sweep**
- PTA variants and MIA options can be crossed on top of this grid.
- GEO matters only when the weighting is sry, sy, hybrid or igofe.
- The IGO coding matters only when igo is in the layer set.

### Validation used to judge communities

- **Script.** `validate_pta.R`, with `TAGS="label=tag,..."` and env `MV`, `OUT`.
- **17 original benchmarks.** Included:
  - GFR–GDR and ROK–PRK separate.
  - NATO and Warsaw Pact pairs together within each bloc, and the two blocs apart.
  - Defence pairs together, by era.
  - Holy Alliance 1835, Triple Alliance 1913, France–Russia 1913, Axis 1938, West 1953 (USA, UK, France).
  - China–Russia together 1950–59 and apart 1970–89.
- **4 newer benchmarks.**
  - Entente apart from the Triple Alliance 1913.
  - Axis apart from the Allies 1938–39.
  - USA with NATO Europe 1950–89.
  - ROK apart from China and Russia 1992–2014.
- **Scoring.** The script reports score17 and score_new4. The NATO–WP row is stored as "apart", so 1 = apart.
