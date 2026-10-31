## The aprime role coding (PI decision 2026-10-31; switch roles = aprime, suffix _aprime, G = 4): aprime_roles() on synthetic disputes, the
## assembly of the groups from roles_aprime.rds (toy cell), the switch and its refusals, and 02_fit_model.R / 03_extract.R end to end for M1c and
## M1cr. Default cells (roles = stageb) read no role file and are unchanged.
suppressPackageStartupMessages({ library(testthat); library(dplyr) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "assemble_cell.R"))

## participants (dispnum, ccode, year, side = COW sidea, orig) -> dd rows (one per participant-year is enough for aprime_roles)
mk_dd <- function(...) bind_rows(...) |> transmute(dispnum = as.integer(dispnum), ccode1 = as.integer(ccode), year = as.integer(year), sidea1 = as.integer(side), orig1 = as.integer(orig))
P <- function(d, cc, side, orig, years = 2000L) tibble(dispnum = d, ccode = cc, year = years, side = side, orig = orig)
## an event row after the direction rule: sidea1 = the initiator's own COW side flag
E <- function(d, i, j, o1, o2, s1, year = 2000L) tibble(dispnum = as.integer(d), initiator = as.integer(i), target = as.integer(j), year = as.integer(year), orig1 = as.integer(o1), orig2 = as.integer(o2), sidea1 = as.integer(s1))

test_that("aprime_roles(): the synthetic cases give the PI's groups", {
  dd <- mk_dd(
    P(1, 10, 1, 1), P(1, 20, 0, 1), P(1, 30, 0, 1),                     # 1: lone originator 10 vs two originator targets -> 0
    P(2, 11, 1, 1), P(2, 12, 1, 1), P(2, 21, 0, 1),                     # 2: two same-side originators 11, 12 -> 1 (both)
    P(3, 13, 1, 1), P(3, 22, 0, 1), P(3, 23, 0, 1),                     # 3: only OPPOSITE-side other originator (23) -> 0
    P(4, 14, 1, 1), P(4, 15, 1, 0), P(4, 24, 0, 1),                     # 4: joiner 15 -> 2; 14 has only a joiner on its side -> 0
    P(5, 16, 1, 1), P(5, 25, 0, 1), P(5, 26, 0, 0, 2000:2001),          # 5: 16 vs originator 25 (-> 0) and late target 26 (-> 3), same D8 group
    P(6, 40, 1, 1), P(6, 41, 1, 1), P(6, 50, 0, 1), P(6, 51, 0, 1), P(6, 52, 0, 1),   # 6: flipped initiators on COW Side B
    P(7, 60, 1, 1), P(7, 70, 0, 1), P(7, 61, 1, 1, 2001)                # 7: partner 61 enters in 2001 only -> side from its first year
  )
  ev <- bind_rows(
    E(1, 10, 20, 1, 1, 1), E(1, 10, 30, 1, 1, 1),
    E(2, 11, 21, 1, 1, 1), E(2, 12, 21, 1, 1, 1),
    E(3, 13, 22, 1, 1, 1),
    E(4, 14, 24, 1, 1, 1), E(4, 15, 24, 0, 1, 1),
    E(5, 16, 25, 1, 1, 1), E(5, 16, 26, 1, 0, 1),
    E(6, 50, 40, 1, 1, 0), E(6, 50, 41, 1, 1, 0),                      # flipped: 50 is on COW Side B with co-originators 51, 52 -> 1
    E(7, 60, 70, 1, 1, 1)
  )
  r <- aprime_roles(ev, dd)
  g <- setNames(r$grp_aprime, paste(r$dispnum, r$initiator, r$target))
  expect_equal(unname(g[c("1 10 20", "1 10 30")]), c(0L, 0L))
  expect_equal(unname(g[c("2 11 21", "2 12 21")]), c(1L, 1L))
  expect_equal(unname(g["3 13 22"]), 0L)
  expect_equal(unname(g[c("4 14 24", "4 15 24")]), c(0L, 2L))
  expect_equal(unname(g[c("5 16 25", "5 16 26")]), c(0L, 3L))        # group 3 overrides and may share a D8 group with group 0
  expect_equal(unname(g[c("6 50 40", "6 50 41")]), c(1L, 1L)); expect_equal(r$n_co_orig[r$initiator == 50L], c(2L, 2L))   # the Side-B partners, not the two Side-A originators
  expect_equal(unname(g["7 60 70"]), 1L); expect_equal(attr(r, "diagnostics")$n_fallback_side, 1L)
  expect_equal(attr(r, "diagnostics")$n_mixed_groups_with_3, 1L)
  ## a flipped lone initiator whose OWN side (B) has no other originator is alone, although COW Side A has two originators
  dd2 <- mk_dd(P(8, 80, 1, 1), P(8, 81, 1, 1), P(8, 90, 0, 1))
  expect_equal(aprime_roles(E(8, 90, 80, 1, 1, 0), dd2)$grp_aprime, 0L)
  ## the same initiator as a Side-A originator would be a co-originator (the side decides)
  expect_equal(aprime_roles(E(8, 80, 90, 1, 1, 1), dd2)$grp_aprime, 1L)
})

test_that("aprime_roles(): side switchers are counted, inconsistent sides and non-constant groups fail loudly", {
  ## a partner on both sides in the onset year counts as same-side and the event is counted in n_switch_decides
  dd <- mk_dd(P(1, 10, 1, 1), P(1, 20, 0, 1), P(1, 11, 1, 1), P(1, 11, 0, 1))
  r <- aprime_roles(E(1, 10, 20, 1, 1, 1), dd)
  expect_equal(r$grp_aprime, 1L); expect_equal(attr(r, "diagnostics")$n_switch_decides, 1L)
  ## the event's (post-flip) side disagrees with the initiator's COW side in the onset year
  dd2 <- mk_dd(P(1, 10, 1, 1), P(1, 20, 0, 1))
  expect_error(aprime_roles(E(1, 10, 20, 1, 1, 0), dd2), "disagrees with the initiator's COW side")
  ## an initiator present on both sides in the onset year with pairs on both sides and one same-side originator: groups 0 / 1 in one D8 group
  dd3 <- mk_dd(P(1, 10, 1, 1), P(1, 10, 0, 1), P(1, 11, 1, 1), P(1, 20, 0, 1), P(1, 21, 1, 0))
  ## (with originator targets on both sides the two pairs are always co-originator pairs, so the constancy assertion cannot be reached through
  ## a switcher with two pairs; the switcher is counted)
  expect_equal(attr(aprime_roles(E(1, 10, 20, 1, 1, 1), dd3), "diagnostics")$n_init_switch, 1L)
  ## originator status varying within a participant is refused
  expect_error(aprime_roles(E(1, 10, 20, 1, 1, 1), mk_dd(P(1, 10, 1, 1, 2000), P(1, 10, 1, 0, 2001), P(1, 20, 0, 1))), "originator status differs")
})

test_that("spec: roles = stageb | aprime, suffix _aprime after _conv, context models only, refused by every preset path, dropped from the perm string at its default", {
  expect_equal(SPEC_CHOICES$roles, c("stageb", "aprime")); expect_equal(SPEC_SUFFIX$roles, c(aprime = "_aprime"))
  expect_equal(which(names(SPEC_SUFFIX) == "roles"), which(names(SPEC_SUFFIX) == "comparison") + 1L)
  ap <- parse_spec(c("--comparison=conventional", "--roles=aprime"))
  expect_equal(cell_name("fmid", 3L, "M1c", ap), "fmid_lag3_all_pr_M1c_conv_aprime"); expect_equal(cell_name("mid", 3L, "M1cr", ap), "mid_lag3_all_pr_M1cr_conv_aprime")
  expect_equal(cell_name("mid", 3L, "M1c", parse_spec("--roles=aprime")), "mid_lag3_all_pr_M1c_aprime")
  for (m in c("M1a", "M1b", "M2b")) expect_error(cell_name("fmid", 3L, m, parse_spec("--roles=aprime")), "the role coding applies to M1c, M1cr only", info = m)
  expect_equal(cell_name("fmid", 3L, "M1c", parse_spec(character(0))), "fmid_lag3_all_pr_M1c")       # default names unchanged
  expect_equal(stage2_group_labels("stageb"), STAGE2_GROUPS); expect_equal(stage2_group_labels("aprime"), c("alone", "co_originator", "joiner", "target_late"))
  expect_false(grepl("roles", perm_part_spec_string(parse_spec("--distance=capital")), fixed = TRUE))
  expect_true(grepl("roles=aprime", perm_part_spec_string(parse_spec("--roles=aprime")), fixed = TRUE))
  expect_error(parse_spec("--roles=4"), "--roles must be one of stageb | aprime", fixed = TRUE)
})

## toy roles: the six toy events (helper-toy-cell.R) with groups 0 / 1 / 2 / 3 / 3 / 0
toy_roles <- tibble(initiator = c(10L, 10L, 30L, 30L, 30L, 40L), target = c(20L, 30L, 20L, 40L, 10L, 20L), year = 1990L, dispnum = c(1L, 4L, 2L, 3L, 6L, 7L),
                    grp_aprime = c(0L, 1L, 2L, 3L, 3L, 0L))

test_that("assembly: roles = aprime takes the groups from roles_aprime.rds (G = 4); the default reads no role file and is unchanged", {
  td <- write_toy_cell(tempfile("toycell"), roles = toy_roles)
  sb <- parse_spec(character(0)); ap <- parse_spec("--roles=aprime")
  A0 <- quiet(assemble_cell_spec("mid", 3L, sb, data_dir = td, model = "M1c")); A1 <- quiet(assemble_cell_spec("mid", 3L, ap, data_dir = td, model = "M1c"))
  expect_equal(A0$roles, "stageb"); expect_equal(A1$roles, "aprime")
  expect_equal(A0$group_table$group, STAGE2_GROUPS); expect_equal(A1$group_table$group, APRIME_GROUPS); expect_equal(A1$group_table$grp, 0:3)
  key <- paste(A1$sy$ccode[A1$ev$sy_row], A1$alts$ccode2[A1$ev$alt_row])
  expect_equal(A1$ev$grp, toy_roles$grp_aprime[match(key, paste(toy_roles$initiator, toy_roles$target))])
  expect_equal(sum(A1$group_table$n_usable), nrow(A1$ev)); expect_true(all(A1$group_table$n_usable > 0L))
  for (nm in c("sy", "alts", "W_cols", "Z_cols", "X_cols")) expect_identical(A1[[nm]], A0[[nm]], info = nm)          # only the groups differ
  expect_identical(A1$ev |> select(-grp), A0$ev |> select(-grp)); expect_equal(A1$ev$w, A0$ev$w)                     # weights unchanged
  out <- capture.output(invisible(assemble_cell_spec("mid", 3L, ap, data_dir = td, model = "M1c")))
  expect_true(any(grepl("roles = aprime: cross-tab of the", out, fixed = TRUE))); expect_true(any(grepl("only group-3 (target joined late) events", out, fixed = TRUE)))
  ## the default cell is identical with or without the role file on disk
  td0 <- write_toy_cell(tempfile("toycell"))
  B0 <- quiet(assemble_cell_spec("mid", 3L, sb, data_dir = td0, model = "M1c"))
  for (nm in c("sy", "alts", "ev", "W_cols", "Z_cols", "X_cols", "group_table")) expect_identical(B0[[nm]], A0[[nm]], info = nm)
  ## refusals: no role file, a missing event row, groups 0-2 not constant within a D8 group
  expect_error(quiet(assemble_cell_spec("mid", 3L, ap, data_dir = td0, model = "M1c")))
  td2 <- write_toy_cell(tempfile("toycell"), roles = toy_roles[-1, ])
  expect_error(quiet(assemble_cell_spec("mid", 3L, ap, data_dir = td2, model = "M1c")), "have no row in roles_aprime.rds")
  td3 <- write_toy_cell(tempfile("toycell"), events_edit = function(e) e |> mutate(dispnum = ifelse(initiator == 10L, 1L, dispnum)), roles = toy_roles |> mutate(dispnum = ifelse(initiator == 10L, 1L, dispnum)))
  expect_error(quiet(assemble_cell_spec("mid", 3L, ap, data_dir = td3, model = "M1c")), "groups 0-2 not constant within 1 D8 groups")
  unlink(c(td, td0, td2, td3), recursive = TRUE)
})

test_that("02_fit_model.R --roles=aprime end to end on the toy cell: M1c and M1cr with G = 4 (conventional), 03_extract.R writes four group tables; M1b refused", {
  for (m in c("m1_context", "m1_context_rho")) if (!file.exists(file.path(root, "stan", "build", paste0(m, ".rds")))) skip(paste("compiled model missing:", m))
  tmp <- tempfile("smoke_ap_"); dir.create(file.path(tmp, "v2"), recursive = TRUE)
  file.copy(file.path(root, "R"), file.path(tmp, "v2"), recursive = TRUE)
  ok <- file.symlink(file.path(root, "stan"), file.path(tmp, "v2", "stan")); if (!ok) file.copy(file.path(root, "stan"), file.path(tmp, "v2"), recursive = TRUE)
  write_toy_cell(file.path(tmp, "v2", "data"), roles = toy_roles)
  dir.create(file.path(tmp, "v2", "output")); dir.create(file.path(tmp, "v2", "logs"))
  v2 <- file.path(tmp, "v2"); out_dir <- file.path(v2, "output")
  run <- function(script, ...) { out <- suppressWarnings(system2("Rscript", c(file.path(v2, "R", script), ...), stdout = TRUE, stderr = TRUE)); list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), out = out) }
  expect_true(run("02_fit_model.R", "--sample=mid", "--lag=3", "--model=M1b", "--roles=aprime", "--dry_run=TRUE")$status != 0L)
  for (m in c("M1c", "M1cr")) {
    cell <- sprintf("mid_lag3_all_pr_%s_conv_aprime", m)
    r1 <- run("02_fit_model.R", "--sample=mid", "--lag=3", paste0("--model=", m), "--comparison=conventional", "--roles=aprime", "--chains=2", "--iter=60", "--warmup=30", "--max_treedepth=6")
    expect_equal(r1$status, 0L, info = paste(m, paste(tail(r1$out, 30), collapse = "\n")))
    expect_true(any(grepl(sprintf("%s stage-2 groups (ev_grp, G = 4", m), r1$out, fixed = TRUE)))
    o <- readRDS(file.path(out_dir, paste0("fit_", cell, ".rds")))
    expect_equal(o$dat$G, 4L); expect_equal(o$group_labels, APRIME_GROUPS); expect_equal(o$fit@par_dims$delta, c(3L, 2L)); expect_equal(o$fit@par_dims$sigma_g, c(4L, 2L))
    expect_true("roles_aprime.rds" %in% o$manifest$file); expect_equal(o$spec$roles, "aprime")
    s1 <- read.csv(file.path(out_dir, paste0("summary_", cell, ".csv")))
    expect_equal(s1$covariate[s1$parameter == "delta[3,1]"], "target_late"); expect_equal(s1$covariate[s1$parameter == "sigma_g[2,2]"], "co_originator")
    if (m == "M1cr") expect_equal(s1$covariate[s1$parameter == "rho_g[4]"], "target_late")
    e1 <- run("03_extract.R", paste0("--fit=", file.path(out_dir, paste0("fit_", cell, ".rds")))); expect_equal(e1$status, 0L, info = paste(tail(e1$out, 30), collapse = "\n"))
    for (g in APRIME_GROUPS) expect_true(file.exists(file.path(out_dir, sprintf("typology_%s_%s.csv", cell, g))), info = g)
  }
  unlink(tmp, recursive = TRUE)
})
