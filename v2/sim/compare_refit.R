## ============================================================================
## compare_refit.R -- does a refit with the reduce_sum programs reproduce an earlier fit within MC error?
## Usage: Rscript v2/sim/compare_refit.R <old_fit.rds> <new_fit.rds> [out.csv]
## Compares posterior means of the structural parameters; z = (mean_new - mean_old) / sqrt(se_old^2 + se_new^2)
## with se = sd / sqrt(n_eff). |z| mostly < 2 means the two runs agree within Monte Carlo error.
## ============================================================================
suppressPackageStartupMessages({ library(rstan); library(dplyr) })
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) >= 2)
old <- readRDS(args[1])$fit; new <- readRDS(args[2])$fit
pars <- c("alpha0", "alpha", "tau", "mu0", "beta", "sigma_x", "sd_re", "sigma_b", "gamma", "Sigma", "Omega")
so <- summary(old, pars = pars)$summary; sn <- summary(new, pars = pars)$summary
stopifnot(identical(rownames(so), rownames(sn)))
keep <- is.finite(so[, "sd"]) & so[, "sd"] > 0
tab <- data.frame(parameter = rownames(so)[keep], mean_old = so[keep, "mean"], mean_new = sn[keep, "mean"],
                  sd_old = so[keep, "sd"], n_eff_old = so[keep, "n_eff"], n_eff_new = sn[keep, "n_eff"])
tab$z <- (tab$mean_new - tab$mean_old) / sqrt(tab$sd_old^2 / tab$n_eff_old + sn[keep, "sd"]^2 / tab$n_eff_new)
tab$abs_diff_in_sd <- abs(tab$mean_new - tab$mean_old) / tab$sd_old
cat(sprintf("%d parameters; max |z| = %.2f; share |z| < 2: %.2f; max |diff| / posterior sd = %.2f\n",
            nrow(tab), max(abs(tab$z)), mean(abs(tab$z) < 2), max(tab$abs_diff_in_sd)))
print(tab |> arrange(desc(abs(z))) |> head(8) |> mutate(across(where(is.numeric), ~ round(.x, 3))))
if (length(args) >= 3) write.csv(tab, args[3], row.names = FALSE)
