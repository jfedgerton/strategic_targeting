# v2 environment (pinned)

Development and simulation container: Ubuntu 24.04.5, R 4.3.3 (2024-02-29), gcc 13.3.0, 2 cores, 7 GB RAM.
Full `sessionInfo()` is in `sessionInfo.txt`.

| package | version |
|---|---|
| rstan | 2.32.5 (Stan 2.32.2) |
| StanHeaders | 2.32.5 |
| brms | 2.20.4 |
| loo | 2.6.0 |
| posterior | 1.5.0 |
| peacesciencer | 1.2.0 (bundled: COW states, NMC, COW majors, Polity/V-Dem, COW-MID 5.0, ATOP, COW contiguity, COW IGO v3) |
| dplyr | 1.1.4 |
| tidyr | 1.3.1 |
| igraph | 1.6.0 |
| testthat | 3.2.1 |
| digest | 0.6.34 |
| BH | 1.90.0.1 (installed from CRAN; the Ubuntu `r-cran-bh` package strips the headers and rstan then fails with "Boost not found") |
| RcppEigen | 0.3.4.0.0 |
| RcppParallel | 5.1.7 |

Install path used here: `apt-get install r-base r-base-dev r-cran-rstan r-cran-brms r-cran-loo r-cran-posterior r-cran-igraph r-cran-dplyr r-cran-tidyr r-cran-testthat r-cran-digest r-cran-rcppeigen r-cran-rcppparallel`, then `install.packages(c("BH", "peacesciencer"))` from CRAN.

## Roar

Stan compilation on Roar needs

```bash
module load gcc/12.2.0
export LD_LIBRARY_PATH="$(dirname "$(gcc -print-file-name=libstdc++.so)"):$LD_LIBRARY_PATH"
```

so that the gcc module's `libstdc++.so.6` (GLIBCXX_3.4.29+) is found before the system library; without the second line compiled Stan models fail at load time with a `GLIBCXX` version error. Both lines are in `v2/slurm/00_env.sh`, which every sbatch script sources. On Roar the pipeline runs under `r/4.2.1-gcc-8.5.0` with rstan 2.32.7 / StanHeaders 2.32.10 from `~/R/x86_64-pc-linux-gnu-library/4.2` (`v2/slurm/00_env.sh`; there is no R 4.3.x module with rstan), so the R minor version differs from the development container's 4.3.3 while rstan stays at 2.32.x; the compiled `.rds` model objects in `v2/stan/build/` are machine-specific and must be rebuilt there (`Rscript v2/R/compile_models.R`, done by `v2/slurm/01_build.sbatch`). `00_env.sh` also exports `R_MAKEVARS_USER=v2/env/Makevars`.

Stan programs use the 2.32 syntax (`array[N] int`, no `sum_to_zero_vector`), so they also compile with cmdstanr >= 2.32 if the PI prefers that backend; `02_fit_model.R` is rstan-only.
