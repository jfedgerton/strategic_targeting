# data_raw

## mic-conf-1.0.csv
- Militarized Interstate Confrontations (MIC) 1.0, confrontation level. Gibler, Miller and colleagues, International Conflict Data Project, University of Alabama.
- Downloaded 2026-10-23 from the project's data page (https://internationalconflict.ua.edu/data-download/), file `confrontations.zip` (https://www.dropbox.com/s/1xjks8jqytb1h9p/confrontations.zip?dl=1), member `mic-conf-1.0.csv`, unmodified.
- sha256 41cef02d0dcb23f91658b2367864b018a517e13d5b0278394dc12775fa5a5ee0 (checked by `01_build_data.R` against `MIC_CONF_SHA256` in `R/00_helpers.R`).
- 1,881 confrontations, 1816-2014. `micnum` keeps the COW MID dispute number; `fatalmin` has no missing values.
- Codebook: https://internationalconflict.ua.edu/wp-content/uploads/2023/07/MICcodebook.pdf
- Used only to fill COW fatality where every COW fatality code is missing (see README, Samples). Cite MIC when reporting `fmid`.
