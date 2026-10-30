## Dual typology stage (b), 2026-10-27: the M1c model (stan/m1_context.stan), the stage-2 event groups, the multi_mid_jan control, the
## main-specification invariance, the dry runs, the Stan program against a brute-force R likelihood and a smoke fit on the toy cell.
suppressPackageStartupMessages({ library(testthat); library(dplyr); library(rstan) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "typology.R")); source(file.path(root, "R", "context_model.R")); source(file.path(root, "R", "assemble_cell.R"))
## write_toy_cell(), quiet() from helper-toy-cell.R; skip_if_no_build(), tiny_data() from helper-fixed-draw.R

## the toy cell with a joiner initiator: 30 joined disputes 2 and 3 after their onset (orig1 = 0 on 30 -> 20 (coalition) and 30 -> 40 (bilateral,
## fighting_at_ref)); the pair-level flag is set on the same rows so the two definitions agree on the toy
toy_joiners <- function(e) e |> mutate(orig1 = ifelse(initiator == 30L & target %in% c(20L, 40L), 0L, 1L), is_joiner_pair = orig1 == 0L)

test_that("spec: M1c is a model code outside the main grid with its own Stan program, cell names that collide with nothing, and the context_b preset", {
  expect_false("M1c" %in% MODELS); expect_true("M1c" %in% ALL_MODELS); expect_equal(CONTEXT_MODELS, "M1c")
  expect_equal(stan_program("M1c"), "m1_context"); expect_equal(stan_program("M1b"), "m1_latent_additive"); expect_equal(stan_program("M2a"), "m2_ame")
  expect_setequal(STAN_PROGRAMS, c("m1_latent_additive", "m2_ame", "m1_context", "m1_latent_prior"))   # m1_latent_prior added 2026-10-30 (prior switches)
  expect_equal(model_family("M1c"), "M1"); expect_equal(model_D_mode("M1c"), 3L); expect_equal(model_joint("M1c"), 1L)
  main <- parse_spec(character(0))
  expect_equal(cell_name("fmid", 3L, "M1c", main), "fmid_lag3_all_pr_M1c"); expect_equal(cell_name("mid", 3L, "M1c", parse_spec("--ww2=drop")), "mid_lag3_all_pr_M1c_noww2")
  existing <- unlist(lapply(c("mid", "fmid", "war"), function(s) lapply(LAGS, function(l) vapply(MODELS, function(m) cell_name(s, l, m, main), ""))))
  expect_false(any(grepl("M1c", existing, fixed = TRUE)))
  expect_equal(sampler_settings("M1c"), SAMPLER$M1)
  expect_error(cell_name("mid", 3L, "M1c", parse_spec("--ame_K=3")), "K applies to the M2 family only")
  expect_true("context_b" %in% names(ROBUSTNESS_PRESETS)); expect_equal(CONTEXT_B_PRESET$cells, c("mid:3:M1c", "fmid:3:M1c"))
  ## the existing names and helpers are untouched
  expect_equal(cell_name("war", 3L, "M2b", main), "war_lag3_all_pr_M2b_K2"); expect_equal(cell_name("fmid", 3L, "M1b", main), "fmid_lag3_all_pr_M1b")
})

test_that("stage2_groups(): bilateral / leader / joiner partition from coal_incl x the initiator's orig1 (initiator-role-only, PI 2026-10-27); role decides; errors on NA", {
  e <- data.frame(coal_incl = c(FALSE, TRUE, FALSE, TRUE), coal_strict = c(FALSE, FALSE, FALSE, TRUE), orig1 = c(1L, 1L, 0L, 0L), is_joiner_pair = c(FALSE, TRUE, TRUE, TRUE))
  expect_equal(stage2_groups(e), c(0L, 1L, 2L, 2L))                      # a joiner initiator in a bilateral pair is a joiner (role decides)
  expect_equal(stage2_groups(e, "coal_strict"), c(0L, 0L, 2L, 2L))
  expect_equal(STAGE2_GROUPS, c("bilateral", "leader", "joiner"))
  ## the pair-level flag plays no part: an originator initiator facing a joiner target (is_joiner_pair TRUE, orig1 = 1) is bilateral / leader
  expect_equal(stage2_groups(e |> mutate(is_joiner_pair = TRUE)), c(0L, 1L, 2L, 2L)); expect_equal(stage2_groups(e |> select(-is_joiner_pair)), c(0L, 1L, 2L, 2L))
  expect_error(stage2_groups(e |> mutate(coal_incl = c(NA, TRUE, FALSE, TRUE))))
  expect_error(stage2_groups(e |> mutate(orig1 = c(1L, NA, 0L, 0L)))); expect_error(stage2_groups(e |> mutate(orig1 = c(2L, 1L, 0L, 0L)))); expect_error(stage2_groups(e |> select(-orig1)))
  ## every row in exactly one group on the real events file when present; the counts differ from the pair-level definition exactly by the
  ## originator-initiator / joiner-target pairs
  f <- file.path(v2_root(), "data", "events.rds")
  if (file.exists(f)) {
    ev <- readRDS(f) |> filter(!escalation)
    g <- stage2_groups(ev)
    expect_equal(length(g), nrow(ev)); expect_true(all(g %in% 0:2))
    expect_equal(sum(g == 2L), sum(ev$orig1 == 0L)); expect_equal(sum(g == 1L), sum(ev$coal_incl & ev$orig1 == 1L)); expect_equal(sum(g == 0L), sum(!ev$coal_incl & ev$orig1 == 1L))
    expect_true(all(g[ev$orig1 == 0L] == 2L)); expect_true(all(ev$is_joiner_pair[g == 2L]))                        # a joiner initiator always makes a joiner pair, not conversely
    expect_equal(sum(ev$is_joiner_pair) - sum(g == 2L), sum(ev$is_joiner_pair & ev$orig1 == 1L))
  }
})

test_that("multilateral_jan_flags(): the W column on a toy world (presence at 1 January, sides, not-yet-participant, former participants, coarse entries in t, side switch)", {
  ## dispute 1: 10 vs 20 (both sides single) 1990-01-01 .. 1990-06-30 -> bilateral, never counts
  ## dispute 2: side A 10 (enters 1991-03-10), 30 (enters 1992-01-01), 90 (1991-03-10 .. 1991-08-01, a FORMER participant at 1 Jan 1992);
  ##            side B 20 (1991-03-10 .. 1992-12-31), 40 (enters 1991-11-05) -> at 1 Jan 1992 present 10, 20, 40 and 30 (day-dated ON 1 January
  ##            counts) -> side A 2, side B 2: multilateral and active; 90 is not present but ALREADY IN (fix 2: its flag is 0, not 1); at
  ##            1 Jan 1991 nobody present (entries in March 1991) -> not active; at 1 Jan 1993 the dispute has ended
  ## dispute 3 (1994-1995): side A 50 (year resolution 1994); side B 60 (year resolution 1994), 70 (month resolution 1994-01), 90 (year
  ##            resolution 1995) -> at 1 Jan 1994 nobody present (coarse entries IN 1994 do not count: fix 1) -> not active; at 1 Jan 1995
  ##            50, 60, 70 present (entered in 1994) -> side A 1, side B 2: multilateral; 90's year-resolution entry in 1995 is NOT an entry
  ##            before 1 January 1995, so 90 is "not yet" and its flag is 1 (its own entry during t must not zero the flag); 80 never joins
  psd <- data.frame(dispnum = c(1L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L), ccode = c(10L, 20L, 10L, 30L, 90L, 20L, 40L, 50L, 60L, 70L, 90L), spell = 1L,
                    spell_first = c(1990L, 1990L, 1991L, 1992L, 1991L, 1991L, 1991L, 1994L, 1994L, 1994L, 1995L), spell_last = c(1990L, 1990L, 1992L, 1992L, 1991L, 1992L, 1992L, 1995L, 1995L, 1995L, 1995L),
                    entry_res = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 3L, 3L, 2L, 3L), entry_y = c(1990L, 1990L, 1991L, 1992L, 1991L, 1991L, 1991L, 1994L, 1994L, 1994L, 1995L),
                    entry_m = c(1L, 1L, 3L, 1L, 3L, 3L, 11L, NA, NA, 1L, NA), entry_d = c(1L, 1L, 10L, 1L, 10L, 10L, 5L, NA, NA, NA, NA),
                    exit_res = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 3L, 3L, 3L, 3L), exit_y = c(1990L, 1990L, 1992L, 1992L, 1991L, 1992L, 1992L, 1995L, 1995L, 1995L, 1995L),
                    exit_m = c(6L, 6L, 12L, 12L, 8L, 12L, 12L, NA, NA, NA, NA), exit_d = c(30L, 30L, 31L, 31L, 1L, 31L, 31L, NA, NA, NA, NA))
  ek <- date_keys(psd$entry_res, psd$entry_y, psd$entry_m, psd$entry_d); xk <- date_keys(psd$exit_res, psd$exit_y, psd$exit_m, psd$exit_d)
  psd$entry_mkey <- ek$mkey; psd$entry_dkey <- ek$dkey; psd$exit_mkey <- xk$mkey; psd$exit_dkey <- xk$dkey
  dd <- bind_rows(data.frame(dispnum = 1L, ccode1 = c(10L, 20L), year = 1990L, sidea1 = c(1, 0)),
                  tidyr::expand_grid(dispnum = 2L, ccode1 = c(10L, 20L, 40L), year = 1991:1992) |> mutate(sidea1 = ifelse(ccode1 == 10L, 1, 0)),
                  data.frame(dispnum = 2L, ccode1 = c(30L, 90L), year = c(1992L, 1991L), sidea1 = 1),
                  tidyr::expand_grid(dispnum = 3L, ccode1 = c(50L, 60L, 70L), year = 1994:1995) |> mutate(sidea1 = ifelse(ccode1 == 50L, 1, 0)),
                  data.frame(dispnum = 3L, ccode1 = 90L, year = 1995L, sidea1 = 0))
  sy <- tidyr::expand_grid(ccode = c(10L, 20L, 30L, 40L, 50L, 60L, 70L, 80L, 90L), year = 1990:1995L)
  out <- multilateral_jan_flags(sy, psd, dd)
  expect_equal(nrow(out), nrow(sy)); expect_true(all(out$multi_mid_jan %in% 0:1))
  md <- as.data.frame(attr(out, "disputes")); expect_equal(md[, c("dispnum", "year")], data.frame(dispnum = c(2L, 3L), year = c(1992L, 1995L)))
  expect_equal(md$n_side_a, c(2L, 1L)); expect_equal(md$n_side_b, c(2L, 2L))
  flag <- function(c, y) out$multi_mid_jan[out$ccode == c & out$year == y]
  expect_equal(flag(10L, 1992L), 0L); expect_equal(flag(20L, 1992L), 0L); expect_equal(flag(40L, 1992L), 0L)    # participants of dispute 2
  expect_equal(flag(30L, 1992L), 0L)                                                                              # day-dated entry ON 1 January 1992: present and already in
  expect_equal(flag(90L, 1992L), 0L)                                                                              # former participant (spell ended August 1991): already in, not "not yet" (fix 2)
  expect_equal(flag(50L, 1992L), 1L); expect_equal(flag(80L, 1992L), 1L)                                          # outsiders
  expect_true(all(out$multi_mid_jan[out$year == 1994L] == 0L))                                                    # coarse entries IN 1994 are not present at 1 January 1994 (fix 1)
  expect_equal(flag(80L, 1995L), 1L); expect_equal(flag(50L, 1995L), 0L); expect_equal(flag(70L, 1995L), 0L)     # dispute 3 at 1 Jan 1995: 50, 60, 70 entered in 1994
  expect_equal(flag(90L, 1995L), 1L)                                                                              # 90's year-resolution entry in 1995 is during t: still "not yet" (fix 1)
  expect_true(all(out$multi_mid_jan[out$year %in% c(1990L, 1991L, 1993L, 1994L)] == 0L))                          # nothing multilateral open in those years
  expect_equal(sum(out$multi_mid_jan), 4L + 6L)                                                                   # 1992: 50, 60, 70, 80; 1995: 10, 20, 30, 40, 80, 90
  ## the change counts against the pre-fix rule: legacy had dispute 3 multilateral at 1 Jan 1994 too (year / month entries in 1994 present), 90-1995
  ## already in (year entry in 1995 <= reference) and 90-1992 "not yet" (spell-year rows only): fix 1 flips 90-1995 and the 1994 outsiders, fix 2 flips 90-1992
  chg <- attr(out, "changed_vs_prefix")
  expect_equal(unname(chg[c("n_multi_dispute_years_legacy", "n_multi_dispute_years", "n_multi_dispute_years_lost", "n_multi_dispute_years_gained")]), c(3L, 2L, 1L, 0L))
  expect_equal(unname(chg[["n_spell_years_coarse_entry_in_t"]]), 4L)                                              # 50, 60, 70 in 1994 and 90 in 1995
  expect_equal(unname(chg[["n_former_participant_triples"]]), 1L)                                                  # (2, 90, 1992)
  expect_equal(unname(chg[c("n_flags_changed_fix1", "n_flags_changed_fix2", "n_flags_changed_total")]), c(6L + 1L, 1L, 8L))   # 1994: 10..40, 80, 90 lose their flag; 90-1995 gains it; 90-1992 loses it
  expect_equal(unname(chg[c("n_flags_legacy", "n_flags")]), c(16L, 10L))                                        # legacy: 1992 five (90 too), 1994 six, 1995 five (90 already in)
  ## a side switch within a year: 40 listed on both sides of dispute 2 in 1992 takes its side of 1991 (side B); without a 1991 row it is an error
  dd2 <- bind_rows(dd, data.frame(dispnum = 2L, ccode1 = 40L, year = 1992L, sidea1 = 1))
  o2 <- multilateral_jan_flags(sy, psd, dd2); expect_equal(attr(o2, "n_side_switch_years"), 1L); expect_identical(o2$multi_mid_jan, out$multi_mid_jan)
  expect_error(multilateral_jan_flags(sy, psd, bind_rows(dd, data.frame(dispnum = 3L, ccode1 = 90L, year = 1995L, sidea1 = 1))), "no single side in the year before")
  ## a present participant without a dispute-year row is an error
  expect_error(multilateral_jan_flags(sy, psd, dd |> filter(!(dispnum == 2L & ccode1 == 30L))), "no dispute-year row")
})

test_that("assembly: the groups are computed for every model and ev$grp travels; M1c appends multi_mid_jan to W (centred) and changes nothing else", {
  td <- write_toy_cell(tempfile("toycell"), events_edit = toy_joiners)
  A  <- quiet(assemble_cell("mid", 3L, "all", data_dir = td))
  Ac <- quiet(assemble_cell("mid", 3L, "all", CONTEXT_MODEL = TRUE, data_dir = td))
  ## groups: 10 -> 20 bilateral, 10 -> 30 leader (coalition), 30 -> 20 joiner, 30 -> 10 leader, 40 -> 20 leader; 30 -> 40 (joiner) has no stage-2 row
  chosen <- paste(A$alts$ccode1[A$ev$alt_row], A$alts$ccode2[A$ev$alt_row])
  expect_equal(A$ev$grp[match(c("10 20", "10 30", "30 20", "30 10", "40 20"), chosen)], c(0L, 1L, 2L, 1L, 1L))
  expect_equal(A$group_table$n_sample, c(1L, 3L, 2L)); expect_equal(A$group_table$n_usable, c(1L, 3L, 1L)); expect_equal(A$group_table$w_usable, c(1, 2.5, 1))
  expect_equal(A$group_table$n_states, c(1L, 3L, 1L)); expect_false(A$context_model); expect_true(Ac$context_model)
  ## the M1c assembly differs from the plain one in W only
  expect_equal(Ac$W_cols, c(A$W_cols, "multi_mid_jan")); expect_equal(Ac$X_cols, A$X_cols); expect_equal(Ac$Z_cols, A$Z_cols)
  expect_identical(Ac$ev, A$ev); expect_identical(Ac$alts, A$alts); expect_identical(Ac$sy[, names(A$sy)[names(A$sy) != "multi_mid_jan"]], A$sy[, names(A$sy)[names(A$sy) != "multi_mid_jan"]])
  expect_equal(Ac$sy$multi_mid_jan, c(1L, 0L, 1L, 0L) - 0.5)                                 # grand-centred 0/1
  expect_equal(A$sy$multi_mid_jan, c(1L, 0L, 1L, 0L))                                         # untouched when not in W
  ## Stan data: identical for the existing models, ev$grp is NOT in it; the M1c data add ev_grp and one W column
  d  <- build_stan_data(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, x_ref = A$x_ref, allow_orphan_initiations = TRUE)
  dc <- build_stan_data(Ac$sy, Ac$alts, Ac$ev, Ac$W_cols, Ac$X_cols, Ac$Z_cols, joint = 1L, x_ref = Ac$x_ref, allow_orphan_initiations = TRUE)
  expect_null(d$ev_grp); expect_equal(dc$P1, d$P1 + 1L); expect_equal(dc$W[, seq_len(d$P1)], d$W); expect_identical(d[setdiff(names(d), c("W", "P1"))], dc[setdiff(names(dc), c("W", "P1"))])
  ## the spec wrapper switches the column on for M1c only
  As <- quiet(assemble_cell_spec("mid", 3L, parse_spec(character(0)), data_dir = td, model = "M1c")); expect_equal(As$W_cols, Ac$W_cols)
  An <- quiet(assemble_cell_spec("mid", 3L, parse_spec(character(0)), data_dir = td, model = "M1b")); expect_equal(An$W_cols, A$W_cols)
  expect_error(assemble_cell_spec("mid", 3L, parse_spec(character(0)), data_dir = td, model = "M3c"))
  expect_error(quiet(assemble_cell("mid", 3L, "all", CONTEXT_MODEL = TRUE, FULL_CONTROLS = TRUE, data_dir = td)), "minimal control set")
  ## a panel without the column is refused for M1c, accepted for the other models
  td0 <- write_toy_cell(tempfile("toycell"), events_edit = toy_joiners)
  p0 <- readRDS(file.path(td0, "stage1_panel_lag3.rds")) |> select(-multi_mid_jan, -n_multi_jan); saveRDS(p0, file.path(td0, "stage1_panel_lag3.rds"))
  man <- read.csv(file.path(td0, "manifest.csv")); man$sha256[man$file == "stage1_panel_lag3.rds"] <- file_digest(file.path(td0, "stage1_panel_lag3.rds")); write.csv(man, file.path(td0, "manifest.csv"), row.names = FALSE)
  expect_error(quiet(assemble_cell("mid", 3L, "all", CONTEXT_MODEL = TRUE, data_dir = td0)), "multi_mid_jan")
  expect_equal(quiet(assemble_cell("mid", 3L, "all", data_dir = td0))$W_cols, A$W_cols)
})

test_that("main-specification invariance: the assembled Stan data and cell names of lag 3 x {mid, fmid, war} x {M1b, M2b} are identical() to 6abe6c2", {
  ## digests of digest::digest(dat, algo = "sha256") of the six cells assembled at 6abe6c2 (before this change); V2_REF_6ABE6C2 may point at a
  ## directory holding the saved <cell>.rds dumps (list(dat, cell, W_cols, X_cols, Z_cols)) for a full identical() comparison
  man_path <- file.path(v2_root(), "data", "manifest.csv")
  skip_if_not(file.exists(man_path) && file.exists(file.path(v2_root(), "data", "alts_lag3_all.rds")), "real data files not built")
  ref <- c("fmid_lag3_all_pr_M1b" = "a645f5cea61962749c67f6ea41dbe34ea11d4d29b2d415b6bbebd0196567865f",
           "fmid_lag3_all_pr_M2b_K2" = "53bc4b0b958e2c8117a3ca9424c432d68b419eca0904dbb304607c96fe4b4b92",
           "mid_lag3_all_pr_M1b" = "e3cc9edff7e351df88b4ef39fec9d326cc4b3aa23eb88d8d43f4cc4c843aeaf8",
           "mid_lag3_all_pr_M2b_K2" = "634d62ac9266439347ef23b4b99c43ced0151fb33cebe250adb4930d5d0e9ed9",
           "war_lag3_all_pr_M1b" = "5c657457973da476e0ad5270bd619cdb55c39156e813b9d31b6ddf0914692829",
           "war_lag3_all_pr_M2b_K2" = "e73a494f2aea13c5fc1ca6a334f08f0c5323056911a42d41d3ecde5a222d9c8e")
  ref_dir <- Sys.getenv("V2_REF_6ABE6C2", unset = "")
  spec <- parse_spec(character(0))
  for (s in c("mid", "fmid", "war")) {
    A <- quiet(assemble_cell_spec(s, 3L, spec, model = "M1b"))
    expect_true("grp" %in% names(A$ev)); expect_equal(sum(A$group_table$n_usable), nrow(A$ev)); expect_true(all(A$group_table$n_usable > 0))
    for (m in c("M1b", "M2b")) {
      dat <- build_stan_data(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, sampled = A$sampled, n_full = A$n_full, x_ref = A$x_ref,
                             K = if (m == "M2b") 2L else NULL, allow_orphan_initiations = TRUE, D_mode = 3L, grainsize = choose_grainsize(nrow(A$sy), 1L))
      cell <- cell_name(s, 3L, m, spec)
      expect_true(cell %in% names(ref)); expect_null(dat$ev_grp); expect_false("multi_mid_jan" %in% A$W_cols)
      expect_equal(digest::digest(dat, algo = "sha256"), unname(ref[cell]), info = cell)
      if (nzchar(ref_dir) && file.exists(file.path(ref_dir, paste0(cell, ".rds")))) {
        r <- readRDS(file.path(ref_dir, paste0(cell, ".rds")))
        expect_identical(dat, r$dat, info = cell); expect_identical(cell, r$cell); expect_identical(list(A$W_cols, A$X_cols, A$Z_cols), list(r$W_cols, r$X_cols, r$Z_cols))
      }
    }
  }
})

test_that("M1c dry runs assemble for mid and fmid at lag 3 (real data): ev_grp partition, multi_mid_jan last in W, the saved dry run matches", {
  man_path <- file.path(v2_root(), "data", "manifest.csv")
  skip_if_not(file.exists(man_path) && file.exists(file.path(v2_root(), "data", "alts_lag3_all.rds")), "real data files not built")
  spec <- parse_spec(character(0))
  for (s in c("mid", "fmid")) {
    A <- quiet(assemble_cell_spec(s, 3L, spec, model = "M1c"))
    expect_equal(tail(A$W_cols, 1), "multi_mid_jan"); expect_equal(length(A$W_cols), 9L)
    expect_true(abs(mean(A$sy$multi_mid_jan)) < 1e-8)                                       # grand-centred
    dat <- build_stan_data(A$sy, A$alts, A$ev, A$W_cols, A$X_cols, A$Z_cols, joint = 1L, x_ref = A$x_ref, allow_orphan_initiations = TRUE, D_mode = 3L)
    dat$ev_grp <- as.integer(A$ev$grp)
    expect_equal(dat$P1, 9L); expect_equal(length(dat$ev_grp), dat$E); expect_true(all(tabulate(dat$ev_grp + 1L, nbins = 3L) > 0))
    expect_equal(as.integer(tabulate(dat$ev_grp + 1L, nbins = 3L)), A$group_table$n_usable)
    ## the saved dry run of 02_fit_model.R --model=M1c (when present and of the current data)
    sd_path <- file.path(v2_root(), "output", sprintf("standata_%s_lag3_all_pr_M1c.rds", s))
    if (file.exists(sd_path)) {
      sd <- readRDS(sd_path); man <- read.csv(man_path, stringsAsFactors = FALSE)
      if (all(sd$manifest$sha256 == man$sha256[match(sd$manifest$file, man$file)])) {
        dat$grainsize <- sd$dat$grainsize
        expect_identical(dat, sd$dat); expect_equal(sd$model, "M1c"); expect_equal(sd$group_table$n_usable, A$group_table$n_usable)
      }
    }
  }
})

## ---- the Stan program: likelihood against brute-force R, the reference reduction to M1b, a smoke fit (needs the compiled programs) -------------
## a small M1c data set from the simulator's world, with the three groups assigned round-robin over the events
tiny_context_data <- function(seed = 123L) {
  set.seed(seed)
  sim <- simulate_targeting("M1b", S = 6L, Tt = 6L, K = 2L, seed = seed, target_rate = 0.6)
  dat <- build_stan_data(sim$sy, sim$alts, sim$ev, c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"), joint = 1L, x_ref = rep(0, 4), D_mode = 3L)
  dat$ev_grp <- as.integer((seq_len(dat$E) - 1L) %% 3L)
  dat
}
context_pars <- function(dat) {
  set.seed(7); S <- dat$S
  list(mu0 = c(0.7, -0.4), delta_lead = c(-0.6, 0.5), delta_join = c(0.3, -0.8), beta = matrix(rnorm(dat$P * 2, 0, 0.3), dat$P, 2),
       a_i = matrix(rnorm(S * 2), S, 2), c_i = rnorm(S, 0, 0.5), b_j = { z <- rnorm(S); 0.6 * (z - mean(z)) }, gamma = c(0.5, -0.3),
       sigma_g = rbind(c(0.9, 1.2), c(1.6, 0.8), c(1.1, 1.9)), rho = 0.4, alpha0 = -0.5, alpha = c(0.3, -0.2), tau = 0.4)
}
context_inits <- function(dat, p) {
  sd_re <- c(0.5, 1, 1); L <- t(chol(matrix(c(1, 0.2, 0.1, 0.2, 1, 0.3, 0.1, 0.3, 1), 3)))
  z_re <- solve(diag(sd_re) %*% L) %*% t(cbind(p$c_i, p$a_i))
  list(alpha0 = p$alpha0, alpha = p$alpha, tau = p$tau, mu0 = p$mu0, delta_lead = p$delta_lead, delta_join = p$delta_join, beta = p$beta,
       sigma_g = p$sigma_g, rho_x = p$rho, sd_re = sd_re, L_Om = L, z_re = z_re, sigma_b = 0.6, z_b = p$b_j / 0.6, gamma = p$gamma)   # sigma_g: a 3 x 2 matrix for array[3] vector[2]
}

test_that("m1_context: the generated log-liks equal the brute-force R likelihood (both stages, per group), and with zero shifts and equal scales the program is M1b", {
  smc <- skip_if_no_build("m1_context"); smb <- skip_if_no_build("m1_latent_additive")
  dat <- tiny_context_data(); p <- context_pars(dat); init <- context_inits(dat, p)
  expect_true(all(tabulate(dat$ev_grp + 1L, nbins = 3L) >= 2L))
  rstan_options(threads_per_chain = 1)
  fit <- sampling(smc, data = dat, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(init), seed = 123, refresh = 0)
  ll2 <- as.numeric(rstan::extract(fit, "log_lik_stage2")$log_lik_stage2[1, ]); ll1 <- as.numeric(rstan::extract(fit, "log_lik_stage1")$log_lik_stage1[1, ])
  R <- stage2_loglik_R_context(dat, p)
  expect_equal(ll2, R$loglik, tolerance = 1e-8)
  eta <- p$alpha0 + as.numeric(dat$W %*% p$alpha) + p$c_i[dat$sy_state] + p$tau * R$lse
  expect_equal(ll1, dbinom(dat$y1, 1, plogis(eta), log = TRUE), tolerance = 1e-8)
  ## generated quantities: mu_event carries the group shift, mu_state is bilateral, Sigma_g and mu0_lead / mu0_join as defined
  me <- rstan::extract(fit, "mu_event")$mu_event[1, , ]; expect_equal(me, R$mu_event, tolerance = 1e-8)
  ms <- rstan::extract(fit, "mu_state")$mu_state[1, , ]; expect_equal(ms, rep(1, dat$S) %o% p$mu0 + p$a_i, tolerance = 1e-8)
  Sg <- rstan::extract(fit, "Sigma_g")$Sigma_g[1, , , ]
  for (g in 1:3) expect_equal(Sg[g, , ], matrix(c(p$sigma_g[g, 1]^2, p$rho * prod(p$sigma_g[g, ]), p$rho * prod(p$sigma_g[g, ]), p$sigma_g[g, 2]^2), 2), tolerance = 1e-8)
  expect_equal(as.numeric(rstan::extract(fit, "mu0_lead")$mu0_lead[1, ]), p$mu0 + p$delta_lead); expect_equal(as.numeric(rstan::extract(fit, "Sigma")$Sigma[1, , ]), as.numeric(Sg[1, , ]))
  ## the log density equals the stage-1 + weighted stage-2 sums plus the priors: compare the data-dependent part across two draws
  up <- function(i) unconstrain_pars(fit, i)
  init2 <- init; init2$delta_lead <- c(0.2, 0.1); p2 <- p; p2$delta_lead <- c(0.2, 0.1); R2 <- stage2_loglik_R_context(dat, p2)
  lp1 <- log_prob(fit, up(init), adjust_transform = FALSE); lp2 <- log_prob(fit, up(init2), adjust_transform = FALSE)
  prior_diff <- sum(dnorm(init2$delta_lead, 0, dat$prior_mu0_sd, log = TRUE)) - sum(dnorm(init$delta_lead, 0, dat$prior_mu0_sd, log = TRUE))
  expect_equal(as.numeric(lp2 - lp1), (R2$lp_stage2 - R$lp_stage2) + prior_diff, tolerance = 1e-6)      # stage 1 is unchanged by delta_lead (bilateral IV)
  ## reduction: zero shifts and three equal scales reproduce the M1b program's log-liks on the same data
  init0 <- init; init0$delta_lead <- c(0, 0); init0$delta_join <- c(0, 0); init0$sigma_g <- matrix(p$sigma_g[1, ], 3, 2, byrow = TRUE)
  fit0 <- sampling(smc, data = dat, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(init0), seed = 123, refresh = 0)
  initb <- init; initb$sigma_x <- p$sigma_g[1, ]; initb$rho_x <- array(p$rho, 1); initb$L_Om_joint <- array(init$L_Om, c(1, 3, 3)); initb$delta_lead <- NULL; initb$delta_join <- NULL; initb$sigma_g <- NULL; initb$L_Om <- NULL
  datb <- dat; datb$ev_grp <- NULL
  fitb <- sampling(smb, data = datb, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(initb), seed = 123, refresh = 0)
  expect_equal(as.numeric(rstan::extract(fit0, "log_lik_stage2")$log_lik_stage2[1, ]), as.numeric(rstan::extract(fitb, "log_lik_stage2")$log_lik_stage2[1, ]), tolerance = 1e-8)
  expect_equal(as.numeric(rstan::extract(fit0, "log_lik_stage1")$log_lik_stage1[1, ]), as.numeric(rstan::extract(fitb, "log_lik_stage1")$log_lik_stage1[1, ]), tolerance = 1e-8)
  ## the program's own checks
  ## (a reject() in transformed data makes rstan return a stanfit with mode != 0 after printing the message, or throw)
  rejected <- function(bad) { f <- tryCatch(suppressWarnings(sampling(smc, data = bad, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(init), seed = 123, refresh = 0)), error = function(e) NULL); is.null(f) || f@mode != 0L }
  bad <- dat; bad$D_mode <- 1L; expect_true(rejected(bad))
  bad <- dat; bad$ev_grp[1] <- 3L; expect_true(rejected(bad))
  bad <- dat; bad$sampled <- 1L; expect_true(rejected(bad))
  ## threads: the log density does not depend on the slicing
  dat2 <- dat; dat2$grainsize <- 2L
  rstan_options(threads_per_chain = 2)
  fit2 <- sampling(smc, data = dat2, algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0, init = list(init), seed = 123, refresh = 0)
  expect_equal(as.numeric(log_prob(fit2, unconstrain_pars(fit2, init), adjust_transform = TRUE)), as.numeric(log_prob(fit, up(init), adjust_transform = TRUE)), tolerance = 1e-8)
  rstan_options(threads_per_chain = 1)
})

test_that("context_model_summary() and context_group_typology() on a short fit of the tiny data; the group tables use the existing classifier", {
  smc <- skip_if_no_build("m1_context")
  dat <- tiny_context_data()
  rstan_options(threads_per_chain = 1)
  fit <- sampling(smc, data = dat, chains = 1, iter = 60, warmup = 30, seed = 123, refresh = 0, pars = c("z_re", "z_b", "L_Om", "re"), include = FALSE)
  cms <- context_model_summary(fit)
  expect_setequal(unique(cms$group), c("bilateral", "leader", "joiner", "shared"))
  expect_equal(sum(cms$quantity == "sigma_g (tolerance scale)"), 6L); expect_equal(sum(cms$quantity == "Sigma_g"), 9L); expect_equal(sum(cms$group == "shared"), 1L)
  expect_equal(cms$mean[cms$parameter == "mu0_lead[1]"], cms$mean[cms$parameter == "mu0[1]"] + cms$mean[cms$parameter == "delta_lead[1]"], tolerance = 1e-8)
  expect_true(all(cms$lo <= cms$mean & cms$mean <= cms$hi))
  n_by_g <- sapply(0:2, function(g) tabulate(dat$sy_state[dat$ev_sy[dat$ev_grp == g]], nbins = dat$S))
  tg <- context_group_typology(fit, as.character(100 + seq_len(dat$S)), r0 = 0.5, n_events_by_group = n_by_g)
  expect_named(tg, STAGE2_GROUPS)
  mu <- mu_state_draws(fit); dl <- rstan::extract(fit, "delta_lead")$delta_lead
  expect_equal(tg$bilateral$mu1_mean, colMeans(mu[, , 1])); expect_equal(tg$leader$mu1_mean, colMeans(mu[, , 1] + dl[, 1])); expect_equal(tg$leader$mu2_mean, colMeans(mu[, , 2] + dl[, 2]))
  expect_equal(tg$joiner$n_events_in_group, as.integer(n_by_g[, 3])); expect_equal(tg$leader$n_events_as_initiator, as.integer(rowSums(n_by_g)))
  for (g in STAGE2_GROUPS) { expect_equal(tg[[g]]$group, rep(g, dat$S)); expect_true(all(abs(rowSums(tg[[g]][, paste0("p_", TYPE_LEVELS)]) - 1) < 1e-8)) }
  ## the label layer of the summary CSV knows the new parameters
  lab <- label_parameters(c("delta_lead[1]", "delta_join[2]", "sigma_g[2,1]", "sigma_g[3,2]", "rho_x", "mu0_join[1]", "mu0[2]"), c("W1", "W2"), paste0("X", 1:4), c("Z1", "Z2"), 3L)
  expect_equal(lab$axis, c("material", "relational", "material", "relational", "material x relational", "material", "relational"))
  expect_equal(lab$covariate[3:4], c("leader", "joiner"))
})

test_that("02_fit_model.R --model=M1c end to end on the toy cell: dry run, a short fit, 03_extract.R writes the stage-b outputs; the driver's context_b preset", {
  for (m in c("m1_context", "m1_latent_additive")) if (!file.exists(file.path(root, "stan", "build", paste0(m, ".rds")))) skip(paste("compiled model missing:", m))
  tmp <- tempfile("smoke_c_"); dir.create(file.path(tmp, "v2"), recursive = TRUE)
  file.copy(file.path(root, "R"), file.path(tmp, "v2"), recursive = TRUE)
  ok <- file.symlink(file.path(root, "stan"), file.path(tmp, "v2", "stan")); if (!ok) file.copy(file.path(root, "stan"), file.path(tmp, "v2"), recursive = TRUE)
  write_toy_cell(file.path(tmp, "v2", "data"), events_edit = toy_joiners)
  dir.create(file.path(tmp, "v2", "output")); dir.create(file.path(tmp, "v2", "logs"))
  v2 <- file.path(tmp, "v2"); out_dir <- file.path(v2, "output")
  run <- function(script, ...) { out <- suppressWarnings(system2("Rscript", c(file.path(v2, "R", script), ...), stdout = TRUE, stderr = TRUE)); list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  ## dry run
  r0 <- run("02_fit_model.R", "--sample=mid", "--lag=3", "--model=M1c", "--dry_run=TRUE")
  expect_equal(r0$status, 0L, info = paste(tail(r0$out, 30), collapse = "\n"))
  sd <- readRDS(file.path(out_dir, "standata_mid_lag3_all_pr_M1c.rds"))
  expect_equal(sd$dat$ev_grp, c(0L, 1L, 2L, 1L, 1L)); expect_equal(tail(sd$W_cols, 1), "multi_mid_jan"); expect_equal(sd$dat$P1, 9L); expect_equal(sd$model, "M1c")
  expect_true(any(grepl("M1c stage-2 groups (ev_grp): bilateral 1 events", r0$out, fixed = TRUE)))
  ## the refusals: M1c with the full controls or a context filter
  expect_true(run("02_fit_model.R", "--sample=mid", "--lag=3", "--model=M1c", "--dry_run=TRUE", "--controls=full")$status != 0L)
  expect_true(run("02_fit_model.R", "--sample=mid", "--lag=3", "--model=M1c", "--dry_run=TRUE", "--context=bilateral")$status != 0L)
  ## a short fit (2 x 80) and its summary
  r1 <- run("02_fit_model.R", "--sample=mid", "--lag=3", "--model=M1c", "--chains=2", "--iter=80", "--warmup=40", "--max_treedepth=6")
  expect_equal(r1$status, 0L, info = paste(tail(r1$out, 30), collapse = "\n"))
  f <- file.path(out_dir, "fit_mid_lag3_all_pr_M1c.rds"); expect_true(file.exists(f))
  o <- readRDS(f)
  expect_equal(o$args$model, "M1c"); expect_equal(o$args$stan_program, "m1_context"); expect_true(o$args$context_model); expect_equal(o$args$init_source, "stan_random")
  expect_equal(o$args$stan_source_sha256, trimws(readLines(file.path(root, "stan", "build", "m1_context.sha256"))[1]))
  expect_equal(o$dat$ev_grp, c(0L, 1L, 2L, 1L, 1L)); expect_equal(o$group_table$n_usable, c(1L, 3L, 1L)); expect_true("grp" %in% names(o$ev))
  expect_true(all(c("delta_lead", "delta_join", "sigma_g", "rho_x", "mu0_lead", "mu0_join", "Sigma_g", "mu_state", "mu_event", "log_lik_stage1", "log_lik_stage2") %in% o$fit@model_pars))
  s1 <- read.csv(file.path(out_dir, "summary_mid_lag3_all_pr_M1c.csv"))
  expect_true(all(c("delta_lead[1]", "delta_join[2]", "sigma_g[3,2]", "rho_x", "mu0_lead[1]") %in% s1$parameter)); expect_false("sigma_x[1]" %in% s1$parameter)
  expect_equal(s1$covariate[s1$parameter == "sigma_g[2,1]"], "leader"); expect_equal(s1$covariate[s1$parameter == "alpha[9]"], "multi_mid_jan")
  ## 03_extract.R: the standard outputs (bilateral) plus the stage-b files; assert_fit_stan_current() accepts the M1c program's sidecar
  e1 <- run("03_extract.R", paste0("--fit=", f)); expect_equal(e1$status, 0L, info = paste(tail(e1$out, 30), collapse = "\n"))
  for (nm in c("typology_mid_lag3_all_pr_M1c.csv", "context_model_mid_lag3_all_pr_M1c.csv", "typology_mid_lag3_all_pr_M1c_bilateral.csv", "typology_mid_lag3_all_pr_M1c_leader.csv", "typology_mid_lag3_all_pr_M1c_joiner.csv", "mu_event_mid_lag3_all_pr_M1c.csv"))
    expect_true(file.exists(file.path(out_dir, nm)), info = nm)
  cm <- read.csv(file.path(out_dir, "context_model_mid_lag3_all_pr_M1c.csv")); expect_true(all(c("delta_lead (shift vs bilateral)", "rho (shared tolerance correlation)", "Sigma_g") %in% cm$quantity))
  tb <- read.csv(file.path(out_dir, "typology_mid_lag3_all_pr_M1c_bilateral.csv")); tm <- read.csv(file.path(out_dir, "typology_mid_lag3_all_pr_M1c.csv"))
  expect_equal(tb$map_type, tm$map_type); expect_equal(tb$n_events_in_group, c(1L, 0L, 0L, 0L))
  tl <- read.csv(file.path(out_dir, "typology_mid_lag3_all_pr_M1c_leader.csv")); expect_equal(tl$n_events_in_group, c(1L, 0L, 1L, 1L)); expect_equal(tl$group, rep("leader", 4))
  ## a doctored hash is refused (the sidecar check covers the new program)
  doc <- o; doc$args$stan_source_sha256 <- paste(rep("0", 64), collapse = ""); df <- file.path(out_dir, "fit_mid_lag3_all_pr_M1c_doctored.rds"); saveRDS(doc, df)
  e2 <- run("03_extract.R", paste0("--fit=", df)); expect_true(e2$status != 0L); expect_true(any(grepl("m1_context", e2$out, fixed = TRUE)))
  ## the robustness driver's context_b preset expands to the two M1c calls (a stub 02_fit_model.R records the arguments)
  rec <- file.path(tmp, "calls.txt")
  writeLines(c("args <- commandArgs(trailingOnly = TRUE)", sprintf("cat(paste(args, collapse = ' '), '\\n', sep = '', file = '%s', append = TRUE)", rec)), file.path(v2, "R", "02_fit_model.R"))
  d <- run("05_robustness.R", "--preset=context_b", "--dry_run=TRUE", "--threads_per_chain=2")
  expect_equal(d$status, 0L, info = paste(tail(d$out, 20), collapse = "\n"))
  calls <- readLines(rec); expect_length(calls, 2L)
  arg <- function(call, nm) unname(sub(paste0("^--", nm, "="), "", grep(paste0("^--", nm, "="), strsplit(call, " ")[[1]], value = TRUE)))
  expect_equal(unname(vapply(calls, arg, "", "sample")), c("mid", "fmid")); expect_true(all(vapply(calls, arg, "", "model") == "M1c")); expect_true(all(vapply(calls, arg, "", "lag") == "3"))
  for (nm in names(SPEC_CHOICES)) expect_true(all(vapply(calls, arg, "", nm) == SPEC_DEFAULTS[[nm]]), info = nm)
  expect_true(any(grepl("cell mid_lag3_all_pr_M1c:", d$out, fixed = TRUE))); expect_true(any(grepl("cell fmid_lag3_all_pr_M1c:", d$out, fixed = TRUE)))
  file.remove(rec); d1 <- run("05_robustness.R", "--preset=context_b", "--cell_index=1"); expect_equal(d1$status, 0L); expect_equal(arg(readLines(rec), "sample"), "fmid")
  expect_true(run("05_robustness.R", "--preset=context_b", "--cell_index=2")$status != 0L)
  for (bad in c("--cells=mid:3:M1c", "--ww2=drop", "--context=bilateral")) expect_true(run("05_robustness.R", "--preset=context_b", bad)$status != 0L, info = bad)
  sb <- readLines(file.path(root, "slurm", "04_robustness.sbatch")); expect_true(any(grepl("^  context_b)", sb)))
  sh <- readLines(file.path(root, "slurm", "submit_robustness.sh")); expect_true(any(grepl("^  context_b)", sh))); expect_true(any(grepl("--array=0-1", sh, fixed = TRUE)))
  unlink(tmp, recursive = TRUE)
})

test_that("sim/summarise_context_b.R: the coverage table over replicate files; too few, mixed-settings and duplicated replicates are refused; the sbatch runs 20 replicates by --rep", {
  ## a temporary v2 with the real R/ and sim/ scripts and synthetic replicate files written the way sim_context_b.R writes them
  tmp <- tempfile("simsum_v2_"); v2 <- file.path(tmp, "v2"); dir.create(v2, recursive = TRUE)
  file.copy(file.path(root, "R"), v2, recursive = TRUE); file.copy(file.path(root, "sim"), v2, recursive = TRUE); dir.create(file.path(v2, "logs"))
  res <- file.path(v2, "sim", "results", "context_b"); unlink(res, recursive = TRUE); dir.create(res, recursive = TRUE)
  for (f in list.files(file.path(v2, "sim", "results"), pattern = "\\.csv$", full.names = TRUE)) file.remove(f)
  pars <- c("mu0[1]", "mu0[2]", "delta_lead[1]", "sigma_g[1,1]", "rho_x", "beta[1,1]", "gamma[1]", "alpha0", "alpha[1]", "tau", "sd_re[1]", "sigma_b")
  blk <- c(rep("stage_b", 5), "stage_2", "stage_2", rep("stage_1", 3), "re_scales", "re_scales")
  truth <- c(-1, 0.3, -1.2, 1.2, 0.4, 0.3, 1, -2, 0.2, 0.5, 0.9, 0.6)
  write_rep <- function(r, dir = res, iter = 500L, covered_all = FALSE) {
    set.seed(1000 + r)
    m <- truth + rnorm(length(truth), 0, 0.2); sdv <- rep(0.3, length(truth)); lo <- m - 1.645 * sdv; hi <- m + 1.645 * sdv
    if (covered_all) { lo <- truth - 1; hi <- truth + 1 }
    rec <- data.frame(parameter = pars, truth = truth, mean = m, sd = sdv, lo = lo, hi = hi, Rhat = 1.001, n_eff = 400, covered = truth >= lo & truth <= hi, bias = m - truth, z = (m - truth) / sdv,
                      block = blk, recovery_target = TRUE, stringsAsFactors = FALSE)
    rec <- rbind(rec, data.frame(parameter = "mu_state 90% coverage over states with events x axes", truth = 0.9, mean = 0.88, sd = NA, lo = NA, hi = NA, Rhat = NA, n_eff = NA, covered = NA, bias = NA, z = NA, block = "mu_state", recovery_target = TRUE))
    rec$rep <- r; rec$seed <- 123L + r; rec$n_init <- 250L; rec$n_noninit <- 250L; rec$chains <- 4L; rec$iter <- iter; rec$warmup <- iter %/% 2L; rec$elapsed_min <- 30; rec$divergences <- 0L; rec$treedepth_hits <- 0L; rec$sample <- "fmid"; rec$lag <- 3L
    write.csv(rec, file.path(dir, sprintf("context_b_recovery_r%d.csv", r)), row.names = FALSE)
    write.csv(data.frame(group = c("bilateral", "leader", "joiner"), n_events = c(100L, 60L, 40L), w_sum = c(90, 30, 20), n_states = c(40L, 20L, 15L), n_state_years = c(90L, 40L, 30L), rep = r, seed = 123L + r,
                         N1 = 500L, N2 = 60000L, E = 200L, S = 217L, sum_y1 = 150L, alpha0_true = -2, n_states_with_events = 60L, n_pattern_redraws = 0L),
              file.path(dir, sprintf("context_b_design_r%d.csv", r)), row.names = FALSE)
  }
  run <- function(...) { out <- suppressWarnings(system2("Rscript", c(file.path(v2, "sim", "summarise_context_b.R"), ...), stdout = TRUE, stderr = TRUE)); list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  for (r in 0:4) write_rep(r)
  expect_true(run()$status != 0L)                                                                   # 5 < 20 replicates by default
  r5 <- run("--min_reps=5"); expect_equal(r5$status, 0L, info = paste(tail(r5$out, 20), collapse = "\n"))
  cov <- read.csv(file.path(v2, "sim", "results", "context_b_coverage.csv"), stringsAsFactors = FALSE)
  expect_equal(nrow(cov), length(pars) + 1L); expect_true(all(cov$n_reps == 5L))
  expect_equal(cov$parameter[seq_along(pars)], c("mu0[1]", "mu0[2]", "delta_lead[1]", "sigma_g[1,1]", "rho_x", "beta[1,1]", "gamma[1]", "alpha0", "alpha[1]", "tau", "sd_re[1]", "sigma_b"))   # fixed parameter order
  ## the coverage rate and bias are recomputed from the files
  all_rec <- do.call(rbind, lapply(0:4, function(r) read.csv(file.path(res, sprintf("context_b_recovery_r%d.csv", r)), stringsAsFactors = FALSE)))
  for (pp in pars) { d <- all_rec[all_rec$parameter == pp, ]; expect_equal(cov$coverage90[cov$parameter == pp], mean(d$covered)); expect_equal(cov$mean_bias[cov$parameter == pp], mean(d$bias)); expect_equal(cov$rmse[cov$parameter == pp], sqrt(mean(d$bias^2))) }
  expect_true(all(cov$coverage90[seq_along(pars)] >= 0 & cov$coverage90[seq_along(pars)] <= 1)); expect_equal(cov$mean_over_reps[cov$block == "mu_state"], 0.88)
  bl <- read.csv(file.path(v2, "sim", "results", "context_b_coverage_blocks.csv"), stringsAsFactors = FALSE)
  expect_equal(bl$block, c("re_scales", "stage_1", "stage_2", "stage_b", "all")); expect_equal(bl$n_params[bl$block == "all"], length(pars)); expect_equal(bl$n_reps, rep(5L, 5))
  expect_equal(bl$coverage90[bl$block == "all"], mean(all_rec$covered[!is.na(all_rec$covered)]))
  rt <- read.csv(file.path(v2, "sim", "results", "context_b_replicates.csv")); expect_equal(rt$rep, 0:4); expect_equal(rt$seed, 123:127); expect_equal(rt$n_leader, rep(60L, 5))
  ## every interval covering -> coverage 1 everywhere
  for (r in 0:4) write_rep(r, covered_all = TRUE); r5 <- run("--min_reps=5"); expect_equal(r5$status, 0L)
  expect_true(all(read.csv(file.path(v2, "sim", "results", "context_b_coverage.csv"))$coverage90[seq_along(pars)] == 1))
  ## a replicate with other settings, and a duplicated replicate index, are refused
  write_rep(5, iter = 2000L); expect_true(run("--min_reps=5")$status != 0L); file.remove(file.path(res, c("context_b_recovery_r5.csv", "context_b_design_r5.csv")))
  write_rep(0); file.copy(file.path(res, "context_b_recovery_r0.csv"), file.path(res, "context_b_recovery_r00.csv")); expect_true(run("--min_reps=5")$status != 0L)
  ## the sbatch: open account, array 0-19, every task calls sim_context_b.R with --rep=<task id> and 4 x 2000 / 1000; the summariser is named
  sb <- readLines(file.path(root, "slurm", "06_sim_context_b.sbatch"))
  expect_true(any(grepl("^#SBATCH --account=open", sb))); expect_true(any(grepl("^#SBATCH --array=0-19", sb)))
  expect_true(any(grepl('Rscript v2/sim/sim_context_b.R --rep="\\$r"', sb))); expect_true(any(grepl("--chains=4 --iter=2000 --warmup=1000", sb, fixed = TRUE)))
  expect_true(any(grepl("summarise_context_b.R", sb, fixed = TRUE))); expect_true(any(grepl("r=\\$\\{SLURM_ARRAY_TASK_ID", sb)))
  ## sim_context_b.R refuses --seed (the seed is 123 + rep) and documents --rep / --smoke
  src <- readLines(file.path(root, "sim", "sim_context_b.R"))
  expect_true(any(grepl("SEED <- V2_SEED + REP", src, fixed = TRUE))); expect_true(any(grepl('"--seed is not an option', src, fixed = TRUE))); expect_true(any(grepl("context_b_smoke", src, fixed = TRUE)))
  expect_true(any(grepl("y1_sim <- rbinom(N1, 1L, plogis(TRUTH$alpha0 + lin))", src, fixed = TRUE)))                 # y1 is simulated from the model
  expect_true(any(grepl("rec$recovery_target <- TRUE", src, fixed = TRUE)))
  ## Codex review 2026-10-27: targets are drawn from the FULL stage-2 set every slot (no removal), and smoke fits have their own directory
  expect_false(any(grepl("setdiff(rows", src, fixed = TRUE)))
  expect_true(any(grepl("picks[k] <- rows_all[sample.int(length(rows_all), 1L, prob = p)]", src, fixed = TRUE)))
  expect_true(any(grepl('fit_dir <- file.path(sim_dir, "fits", if (SMOKE) "context_b_smoke" else "context_b")', src, fixed = TRUE)))
  expect_true(any(grepl("file.path(fit_dir, sprintf(\"context_b_sim_%s.rds\", tag))", src, fixed = TRUE)))
  unlink(tmp, recursive = TRUE)
})
