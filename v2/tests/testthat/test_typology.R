## Shape and logic tests for typology.R (audit F1: never let an array slice collapse to a scalar)
suppressPackageStartupMessages({ library(testthat) })
root <- Sys.getenv("V2_ROOT"); if (root == "") root <- normalizePath(file.path(getwd(), "v2"), mustWork = FALSE)
source(file.path(root, "R", "00_helpers.R")); source(file.path(root, "R", "typology.R"))

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
