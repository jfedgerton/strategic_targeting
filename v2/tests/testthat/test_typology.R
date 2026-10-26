## Shape and logic tests for typology.R (audit F1: never let an array slice collapse to a scalar)
suppressPackageStartupMessages({ library(testthat) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "spec.R")); source(file.path(root, "R", "helpers_data.R")); source(file.path(root, "R", "helpers_model.R")); source(file.path(root, "R", "typology.R"))

test_that("classify_types requires a [draws, S, 2] array and r0", {
  expect_error(classify_types(matrix(rnorm(20), 10, 2), r0 = 0.5))              # 2-d
  expect_error(classify_types(array(rnorm(60), c(10, 2, 3)), r0 = 0.5))         # last dim != 2
  expect_error(classify_types(array(rnorm(40), c(10, 2, 2))))                    # r0 missing
  expect_error(classify_types(array(rnorm(40), c(10, 2, 2)), r0 = -1))
})

test_that("planted centres classify into the right regions with posterior probabilities", {
  set.seed(123)
  centres <- rbind(c(2, 2), c(0, 0), c(-2, -2), c(2, -2), c(-2, 2))
  mu <- array(NA_real_, c(500, 5, 2))
  for (s in 1:5) { mu[, s, 1] <- centres[s, 1] + rnorm(500, 0, 0.2); mu[, s, 2] <- centres[s, 2] + rnorm(500, 0, 0.2) }
  cls <- classify_types(mu, r0 = 0.8, state_labels = LETTERS[1:5])
  expect_equal(nrow(cls), 5)
  expect_equal(cls$map_type, c("bully", "peer", "challenger", "hybrid_pm", "hybrid_mp"))
  expect_true(all(cls$map_prob > 0.95))
  expect_true(all(abs(rowSums(cls[, paste0("p_", TYPE_LEVELS)]) - 1) < 1e-8))
  expect_equal(cls$state, LETTERS[1:5])
  ## r and theta are per-state, not a recycled scalar (the F1 failure mode)
  expect_gt(length(unique(cls$r_mean)), 1)
  expect_equal(cls$r_mean[2], mean(sqrt(mu[, 2, 1]^2 + mu[, 2, 2]^2)))
  expect_equal(cls$theta_circmean[1], atan2(2, 2), tolerance = 0.05)
})

test_that("r0 changes peer / non-peer boundary as documented", {
  mu <- array(0, c(100, 1, 2)); mu[, 1, 1] <- 0.6; mu[, 1, 2] <- 0.6   # r = 0.85
  expect_equal(classify_types(mu, r0 = 1.0)$map_type, "peer")
  expect_equal(classify_types(mu, r0 = 0.5)$map_type, "bully")
})

test_that("S = 1 and uncertainty split probabilities correctly", {
  set.seed(1)
  mu <- array(NA_real_, c(1000, 1, 2)); mu[, 1, 1] <- rnorm(1000, 0, 1.5); mu[, 1, 2] <- rnorm(1000, 0, 1.5)
  cls <- classify_types(mu, r0 = 0.5)
  expect_equal(nrow(cls), 1)
  expect_true(all(cls[, paste0("p_", TYPE_LEVELS)] > 0.02))        # mass in every region
  expect_equal(sum(cls[, paste0("p_", TYPE_LEVELS)]), 1)
})

test_that("Procrustes alignment undoes a planted rotation and leaves u'v unchanged", {
  set.seed(123)
  S <- 6; K <- 2; n_it <- 5
  U0 <- matrix(rnorm(S * K), S, K); V0 <- matrix(rnorm(S * K), S, K)
  th <- 0.7; R <- matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2)
  U <- array(NA_real_, c(n_it, 2, S, K)); V <- U
  for (it in 1:n_it) { U[it, 1, , ] <- U0; V[it, 1, , ] <- V0; U[it, 2, , ] <- U0 %*% R; V[it, 2, , ] <- V0 %*% R }
  al <- procrustes_align(U, V)
  expect_equal(al$U[1, 2, , ], U0, tolerance = 1e-8)
  expect_equal(al$V[1, 2, , ], V0, tolerance = 1e-8)
  expect_equal(al$U[1, 2, , ] %*% t(al$V[1, 2, , ]), U0 %*% t(V0), tolerance = 1e-8)
  expect_equal(max(abs(al$rotations[[1]] - diag(K))), 0, tolerance = 1e-8)
})

test_that("anchor choice uses initiation counts, not label order", {
  sy <- data.frame(state_idx = c(1, 2, 3, 3, 3, 2), y1 = c(0, 1, 1, 1, 1, 1))
  ev <- data.frame(sy_row = c(2, 3, 4, 5, 6), alt_row = 1:5)
  expect_equal(choose_anchors(sy, ev, 2), c(3L, 2L))
  expect_equal(choose_anchors(sy, ev, 1, S = 10), 3L)
})

test_that("Codex counterexample: equal radii across entities are valid input (bully at (1,1), challenger at (-1,-1))", {
  set.seed(2)
  mu <- array(NA_real_, c(300, 2, 2))
  mu[, 1, 1] <- 1 + rnorm(300, 0, 0.05); mu[, 1, 2] <- 1 + rnorm(300, 0, 0.05)
  mu[, 2, ] <- -mu[, 1, ]                                   # exactly mirrored: identical radii, different entities
  cls <- classify_types(mu, r0 = 0.5)
  expect_equal(cls$map_type, c("bully", "challenger"))
  expect_equal(cls$r_mean[1], cls$r_mean[2])
})

test_that("shape / identity guard: recycled draws and duplicate ids are rejected; shared state-year draws allowed when declared", {
  set.seed(3)
  one <- cbind(rnorm(100, 1, 0.1), rnorm(100, 1, 0.1))
  recycled <- array(0, c(100, 3, 2)); for (s in 1:3) recycled[, s, ] <- one   # the F1 failure mode: one vector for every entity
  expect_error(classify_types(recycled, r0 = 0.5), "recycled")
  expect_error(classify_types(array(rnorm(400), c(100, 2, 2)), r0 = 0.5, state_labels = c("a", "a")))
  ## several events of one state-year share mu draws: legitimate for mu_event
  cls <- classify_types(recycled, r0 = 0.5, state_labels = c("e1", "e2", "e3"), allow_identical_entities = TRUE)
  expect_equal(nrow(cls), 3); expect_true(all(cls$map_type == "bully"))
})

## ---- appendix: 9-cell classification, per-axis bands, reach ---------------------------------------
test_that("planted points in all nine cells classify correctly", {
  set.seed(9)
  centres <- rbind(bully = c(2, 2), challenger = c(-2, -2), peer = c(0, 0),
                   mat_opportunist = c(2, 0), rel_opportunist = c(0, 2), mat_challenger = c(-2, 0), rel_challenger = c(0, -2),
                   compensator_A = c(2, -2), compensator_B = c(-2, 2))
  mu <- array(NA_real_, c(400, 9, 2))
  for (s in 1:9) { mu[, s, 1] <- centres[s, 1] + rnorm(400, 0, 0.1); mu[, s, 2] <- centres[s, 2] + rnorm(400, 0, 0.1) }
  c9 <- classify_cells9(mu, band_M = 0.5, band_R = 0.5, state_labels = rownames(centres))
  expect_equal(c9$modal_cell, rownames(centres))
  expect_true(all(c9$modal_prob > 0.95))
  expect_equal(c9$modal_cell_label[c9$state == "compensator_A"], "compensator A")
  expect_true(all(abs(rowSums(c9[, paste0("p_", names(CELL9_LABELS))]) - 1) < 1e-8))
  ## direction convention: "down" = mu > band (targets below the initiator on that axis)
  expect_equal(c9$p_M_down[c9$state == "bully"], 1); expect_equal(c9$p_R_up[c9$state == "challenger"], 1)
})

test_that("per-axis bands work independently", {
  mu <- array(0, c(50, 1, 2)); mu[, 1, 1] <- 0.8; mu[, 1, 2] <- 0.3       # M clearly down, R small
  expect_equal(classify_cells9(mu, band_M = 0.5, band_R = 0.5)$modal_cell, "mat_opportunist")
  expect_equal(classify_cells9(mu, band_M = 0.5, band_R = 0.2)$modal_cell, "bully")          # narrower R band
  expect_equal(classify_cells9(mu, band_M = 1.0, band_R = 0.5)$modal_cell, "peer")           # wider M band
  expect_error(classify_cells9(mu, band_M = 0.5, band_R = -1))
})

test_that("reach summaries are computed correctly on a known array", {
  mu <- array(0, c(4, 2, 2))
  mu[, 1, 1] <- c(3, 3, 3, 3); mu[, 1, 2] <- c(4, 4, 4, 4)          # r = 5 exactly
  mu[, 2, 1] <- c(-1, 1, -1, 1); mu[, 2, 2] <- c(0, 0, 0, 0)        # r = 1, |mu1| = 1, |mu2| = 0
  rc <- reach_summary(mu, state_labels = c("A", "B"))
  expect_equal(rc$r_mean, c(5, 1)); expect_equal(rc$absmu1_mean, c(3, 1)); expect_equal(rc$absmu2_mean, c(4, 0))
  expect_equal(rc$r_lo, c(5, 1)); expect_equal(rc$r_hi, c(5, 1))
  rc1 <- reach_summary(array(c(-2, 2, -2, 2), c(4, 1, 1)))
  expect_equal(rc1$r_mean, 2)
})

test_that("focus: pooled Sigma summary and per-state dispersion with the < 5 events flag", {
  sig <- array(0, c(3, 2, 2)); sig[, 1, 1] <- 4; sig[, 2, 2] <- 1; sig[, 1, 2] <- sig[, 2, 1] <- 1
  fs <- focus_sigma(sig, D_mode = 3L)
  expect_equal(fs$mean[fs$quantity == "sd_M"], 2); expect_equal(fs$mean[fs$quantity == "corr"], 0.5)
  x <- rbind(c(1, 1), c(3, 1), c(2, 3), c(2, -1), c(2, 1), c(0, 0), c(1, 1))
  mu_ev <- rbind(matrix(2, 5, 2), matrix(0, 2, 2)); mu_ev[1:5, 2] <- 1
  st <- c(1, 1, 1, 1, 1, 2, 2)
  f <- focus_states(x, mu_ev, st, state_labels = c("A", "B"), D_mode = 3L)
  expect_equal(f$n_events, c(5, 2)); expect_equal(f$unreliable_lt_min_events, c(FALSE, TRUE))
  expect_equal(f$sd_M[1], sd(c(-1, 1, 0, 0, 0))); expect_equal(f$sd_R[1], sd(c(0, 0, 2, -2, 0)))
})

test_that("focus: one-dimensional models are labelled by their own axis (Codex review of 36c0f7b)", {
  sig1 <- array(0.81, c(3, 1, 1))                      # sigma_R = 0.9 in a relational-only model
  fr <- focus_sigma(sig1, D_mode = 2L)
  expect_equal(fr$quantity, "sd_R"); expect_equal(fr$mean, 0.9)
  fm <- focus_sigma(sig1, D_mode = 1L)
  expect_equal(fm$quantity, "sd_M")
  x <- matrix(c(1, 3, 2, 2, 2), ncol = 1); mu_ev <- matrix(2, 5, 1); st <- rep(1, 5)
  f_r <- focus_states(x, mu_ev, st, state_labels = "A", D_mode = 2L)
  expect_true(is.na(f_r$sd_M)); expect_equal(f_r$sd_R, sd(c(-1, 1, 0, 0, 0)))
  f_m <- focus_states(x, mu_ev, st, state_labels = "A", D_mode = 1L)
  expect_true(is.na(f_m$sd_R)); expect_equal(f_m$sd_M, sd(c(-1, 1, 0, 0, 0)))
  expect_error(focus_sigma(sig1))                       # D_mode is required
  expect_error(focus_sigma(sig1, D_mode = 3L))          # D = 1 array cannot be a both-dimensions model
})

test_that("P3 (2026-10-25): theta intervals are circular -- a cloud straddling +-pi gets a narrow interval around its circular mean", {
  set.seed(9)
  n <- 4000
  ## entity 1: challenger-like cloud at angle pi (mu1 = -1, mu2 ~ 0): plain quantiles of atan2 would give about [-3.1, 3.1]
  th1 <- pi + rnorm(n, 0, 0.1); m1 <- cbind(cos(th1), sin(th1))
  ## entity 2: bully at pi/4 with the same spread
  th2 <- pi / 4 + rnorm(n, 0, 0.1); m2 <- cbind(cos(th2), sin(th2))
  mu <- array(NA_real_, c(n, 2, 2)); mu[, 1, ] <- m1; mu[, 2, ] <- m2
  cls <- classify_types(mu, r0 = 0.5, state_labels = c("challenger_axis", "bully"))
  expect_equal(abs(cls$theta_circmean[1]), pi, tolerance = 0.02)
  expect_lt(cls$theta_width[1], 0.5); expect_lt(cls$theta_width[2], 0.5)                           # about 2 * 1.645 * 0.1 = 0.33
  expect_equal(cls$theta_width[1], cls$theta_width[2], tolerance = 0.1)
  expect_true(all(abs(c(cls$theta_lo, cls$theta_hi)) <= pi))
  ## the bounds are the circular mean +- the deviation quantiles (wrapped): entity 2 is far from the wrap, so lo < mean < hi
  expect_lt(cls$theta_lo[2], cls$theta_circmean[2]); expect_gt(cls$theta_hi[2], cls$theta_circmean[2])
  ## entity 1 straddles the wrap: lo is just below +pi and hi just above -pi, i.e. lo > hi numerically but the arc is short
  expect_gt(cls$theta_lo[1], pi - 0.3); expect_lt(cls$theta_hi[1], -pi + 0.3)
  ## the plain quantile interval would have been nearly the whole circle
  expect_gt(diff(quantile(atan2(m1[, 2], m1[, 1]), c(0.05, 0.95))), 6)
  expect_equal(cls$map_type, c("challenger", "bully"))
})
