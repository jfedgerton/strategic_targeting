## helper-toy-cell.R -- a toy cell on disk shared by test_data_coding.R and test_context.R (testthat sources helper-*.R first).
## Moved here from test_data_coding.R on 2026-10-25 when the coalition-context columns were added to the events fixture.
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "spec.R")); source(file.path(root, "R", "helpers_data.R")); source(file.path(root, "R", "helpers_model.R"))
## ---- a toy cell on disk (panel / alts / events / manifest) so that assemble_cell() and build_stan_data() run end to end ----------
## states 10, 20, 30, 40 in 1990, h = 3. Panel: 10 (y1 = 1), 20 (y1 = 0, never initiates), 30 (y1 = 1), 40 (y1 = 1).
## Alternatives (ccode1 -> ccode2: in_iv, in_s2):
##   10: 20 (1,1), 30 (1,1), 40 (0,1)   40 is stage-2 only (fighting 10 at January, their dispute ended before 10's reference date)
##   20: 10, 30, 40 all (1,1)           a non-initiator: identical sets
##   30: 10 (1,1), 20 (1,1), 40 (1,0)   40 is IV only (fighting 30 at 30's reference date but not at January)
##   40: 10, 20, 30 all (1,1)
## Events (mid): 10 -> 20 (dispute 1); 10 -> 30 (disputes 4;5, one of them WWII: is_ww2 TRUE, in_mid_noww2 TRUE, in_fmid_noww2 FALSE);
##   30 -> 20 (dispute 2); 30 -> 40 (dispute 3, fighting_at_ref: stage 1 only, its row is IV only); 30 -> 10 (dispute 6, pure WWII);
##   40 -> 20 (dispute 7, pure WWII: under --drop_ww2 state-year 40 loses y1)
write_toy_cell <- function(dir, alts_edit = identity, events_edit = identity) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  panel <- tibble(ccode = c(10L, 20L, 30L, 40L), year = 1990L, h = 3L, y1_mid = c(1L, 0L, 1L, 1L), y1_fmid = c(1L, 0L, 1L, 1L),
                  mon_prior_mids_h = c(2L, 0L, 1L, 3L), mon_peace_years_h = c(3L, 10L, 1L, 2L),
                  dM_h = c(0.5, -0.2, 0.1, -0.4), dRpr_h = c(0.3, 0.0, -0.6, 0.2), dReig_h = c(0.1, 0.2, 0.3, 0.4),
                  M_lvl = c(0.8, -0.3, -0.1, -0.4), R_lvl = c(0.2, 0.5, -0.4, -0.3), log_n_system_tm1 = log(4),
                  cowmaj_tm1 = c(1L, 0L, 0L, 1L), dem_tm1 = c(1L, 1L, 0L, 0L), ongoing_tm1 = c(0L, 0L, 1L, 1L))
  a <- tidyr::expand_grid(ccode1 = c(10L, 20L, 30L, 40L), ccode2 = c(10L, 20L, 30L, 40L)) |> filter(ccode1 != ccode2) |>
    mutate(year = 1990L, h = 3L, in_iv = 1L, in_s2 = 1L)
  a$in_iv[a$ccode1 == 10L & a$ccode2 == 40L] <- 0L
  a$in_s2[a$ccode1 == 30L & a$ccode2 == 40L] <- 0L
  set.seed(1)
  alts <- a |> mutate(dM = rnorm(n()), dR_pr = rnorm(n()), dR_eig = rnorm(n()), dy_prior_mids_h = rpois(n(), 1), dy_peace_years_h = 1L + rpois(n(), 5),
                      log_capdist = log1p(runif(n(), 100, 5000)), dM_j_h = rnorm(n()), dRpr_j_h = rnorm(n()), dch_j_missing = as.integer(ccode2 == 40L),
                      contig_land = as.integer(ccode1 + ccode2 == 30L), alliance_def_off = as.integer(ccode1 == 40L), joint_dem = as.integer(ccode1 %in% c(10L, 20L) & ccode2 %in% c(10L, 20L)),
                      major_count = as.integer(ccode1 %in% c(10L, 40L)) + as.integer(ccode2 %in% c(10L, 40L)), dy_mid_tm1 = as.integer(ccode1 == 30L & ccode2 == 40L), polrel = 1L) |>
    mutate(dM_j_h = ifelse(dch_j_missing == 1L, 0, dM_j_h), dRpr_j_h = ifelse(dch_j_missing == 1L, 0, dRpr_j_h)) |> alts_edit()
  ev <- tibble(initiator = c(10L, 10L, 30L, 30L, 30L, 40L), target = c(20L, 30L, 20L, 40L, 10L, 20L), year = 1990L,
               dispnum = c(1L, 4L, 2L, 3L, 6L, 7L), dispnums = c("1", "4;5", "2", "3", "6", "7"), n_disputes = c(1L, 2L, 1L, 1L, 1L, 1L),
               escalation = FALSE, is_joiner_pair = FALSE, is_ww2 = c(FALSE, TRUE, FALSE, FALSE, TRUE, TRUE), target_has_pos_tm1 = TRUE,
               fighting_at_ref = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), war_disp = FALSE,
               in_mid = TRUE, in_fmid = TRUE, in_fmidcow = FALSE, in_fmidub = TRUE, in_war = FALSE,
               in_mid_noww2 = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE), in_fmid_noww2 = c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE),
               in_fmidcow_noww2 = FALSE, in_fmidub_noww2 = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE), in_war_noww2 = FALSE,
               ## D1 direction flags (30 -> 10 is a reversed war-level pair in this toy) and D8 weights (10's two targets in 1990 are two
               ## disputes, weight 1 each; 30 -> 20 and 30 -> 40 are coalition targets of dispute 2;3? no: separate disputes -> 1;
               ## 30 -> 10 (dispute 6) carries w_mid 0.5 to exercise a non-unit weight)
               dir_flipped = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE), dir_source = c("sideA", "sideA", "sideA", "sideA", "cow_war_pair", "sideA"),
               w_mid = c(1, 1, 1, 1, 0.5, 1), w_fmid = c(1, 1, 1, 1, 0.5, 1), w_fmidcow = 0, w_fmidub = c(1, 1, 1, 1, 0.5, 1), w_war = 0,
               ## coalition context (2026-10-25): 10 -> 20 and 30 -> 40 bilateral; 10 -> 30 (collapsed from 4;5), 30 -> 20, 30 -> 10 and 40 -> 20 coalition
               ## under the inclusive convention; under the strict convention 30 -> 20 is bilateral as well
               coal_incl = c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE), coal_strict = c(FALSE, TRUE, FALSE, FALSE, TRUE, TRUE), coal_onset_ambiguous = FALSE) |> events_edit()
  saveRDS(panel, file.path(dir, "stage1_panel_lag3.rds")); saveRDS(alts, file.path(dir, "alts_lag3_all.rds")); saveRDS(ev, file.path(dir, "events.rds"))
  man <- lapply(c("stage1_panel_lag3.rds", "alts_lag3_all.rds", "events.rds"), function(f) { o <- readRDS(file.path(dir, f))
    data.frame(file = f, nrow = nrow(o), ncol = ncol(o), sha256 = file_digest(file.path(dir, f)), note = "toy") })
  write_manifest(man, file.path(dir, "manifest.csv"))
  invisible(dir)
}
quiet <- function(expr) { out <- NULL; capture.output(out <- expr); out }
