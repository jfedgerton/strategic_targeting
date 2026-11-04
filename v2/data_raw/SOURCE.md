# data_raw

## mic-conf-1.0.csv
- Militarized Interstate Confrontations (MIC) 1.0, confrontation level. Gibler, Miller and colleagues, International Conflict Data Project, University of Alabama.
- Downloaded 2026-10-23 from the project's data page (https://internationalconflict.ua.edu/data-download/), file `confrontations.zip` (https://www.dropbox.com/s/1xjks8jqytb1h9p/confrontations.zip?dl=1), member `mic-conf-1.0.csv`, unmodified.
- sha256 41cef02d0dcb23f91658b2367864b018a517e13d5b0278394dc12775fa5a5ee0 (checked by `01_build_data.R` against `MIC_CONF_SHA256` in `R/00_helpers.R`).
- 1,881 confrontations, 1816-2014. `micnum` keeps the COW MID dispute number; `fatalmin` has no missing values.
- Codebook: https://internationalconflict.ua.edu/wp-content/uploads/2023/07/MICcodebook.pdf
- Second fatality source, used only where every COW fatality code is missing (see README, Samples). Cite MIC when reporting `fmid`.

## GML MID 2.2.1 (not stored here)
- Read from the bundled peacesciencer data `gml_mid_ddydisps` (Gibler, Miller and Little 2016, "An Analysis of the Militarized Interstate Dispute (MID) Dataset, 1816-2001", International Studies Quarterly 60(4)). Version 2.2.1, 1816-2010.
- Third fatality source, used only when COW and MIC are both missing. Cite GML when reporting `fmid`.

## dymid/ (the interacting-dyads foundation, `--events=dymid`; added 2026-11-04)
Both files were downloaded on 2026-11-04 from the Correlates of War project and are stored unmodified. The build checks the two CSVs against the sha256 constants `DYMID_CSV_SHA256` / `DDIW_CSV_SHA256` in `R/helpers_dymid.R` and stops on a mismatch. Read as data only (`read.csv`).

### dymid/dyadic_mid_4.03.csv, dymid/dyadic_mid_4.03_codebook.pdf
- Dyadic Militarized Interstate Disputes 4.03 (Maoz, Johnson, Kaplan, Ogunkoya and Shreve 2019, "The Dyadic Militarized Interstate Disputes (MIDs) Dataset Version 3.0: Logic, Characteristics, and Comparisons to Alternative Datasets", Journal of Conflict Resolution 63(3): 811-835). 1816-2014.
- Source zip https://correlatesofwar.org/wp-content/uploads/dyadic_mid_4.03_update.zip (sha256 48775e03f2d8d9ea69436a047e3e613256d0e9533bb921dff8b94e6e31cdc07e), members `dyadic_mid_4.03_update/dyadic_mid_4.03.csv` and `dyadic_mid_4.03_update/Dyadic MID Codebook V4.0.pdf` (stored as `dyadic_mid_4.03_codebook.pdf`; the codebook text is "Version 4.021, June 18, 2021").
- sha256 of the CSV 9d26eedd792a0067f590f038c3b7b240bd373cd70200ee729241e1078c6957c3; codebook PDF f3963d8abec03d11417d653b2fddbd67b99891b1d28f55146b995be6dafd8c7e.
- 10,358 directed dyad-year rows (both directions of each dyad-year), 2,371 disputes, 3,511 undirected dyadic disputes (disno, pair). In this version `dyindex` is an integer equal to `disno` on every row (the decimal dyad suffix of 4.01 is gone), so the build keys a dyadic dispute by (disno, unordered pair); see README section 6.9.
- Known data issues handled by the build (logged): one dyad-year whose two directed copies differ in start month (3250, Russia - China 1900; the A < B copy is used).

### dymid/directed_dyadic_war.csv, dymid/directed_dyadic_war_codebook.pdf
- The Directed Dyadic Interstate War Dataset (Maoz; same citation as above), 1816-2010. A dyadic war = sustained battles between the two states' armed forces with at least 1,000 battle deaths across both; `warnum` = COW war number; `disno` / `dyindex` link to the dyadic MID data (dyindex with its decimal suffix).
- Source zip https://correlatesofwar.org/wp-content/uploads/dyadic_mid_4.01.zip (sha256 983a1eb065b97f81448426b530b0ccf5a5c26136a186f6c34585469a4b078567), members `directed_dyadic_war.csv` and `The Directed Dyadic Interstate War Dataset Codebook.pdf` (stored as `directed_dyadic_war_codebook.pdf`). The same zip also holds Dyadic MID 4.01, not used.
- sha256 of the CSV aa2eb5244a0cb0b579867e50eb02e4297fb29250e67120fb132e3c0618dceda0; codebook PDF 461401c77b3eaea2588abaa3ef2f7a8506cf78b4b6202331a2eec23710b06bb1.
- 1,364 directed war dyad-year rows, 293 undirected war dyads (warnum, disno, pair) in 91 COW wars. `batdtha` / `batdthb` are sometimes national totals (Germany - Poland 1939 lists 3.5 million): never used as dyadic magnitudes.
- Known data issue repaired by the build (logged): USA - Austria-Hungary 1918 (war 106) has day and month swapped in both dates (start month 24, end month 3 / day 11); rows with a start month above 12 and a day of at most 12 get both dates swapped.
- Cite both Maoz datasets when reporting any `_dy` result.
