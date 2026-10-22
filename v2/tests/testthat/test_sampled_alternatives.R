## Monte Carlo reproductions of the Codex review cases for sampled alternatives (items 3 and 4).
## Pure R; no Stan. The "old" scheme forced every observed target of a state-year into one shared
## set and scored each event as a categorical draw over it; the "new" scheme gives every event its
## own set (chosen + m uniform others) and draws the IV set independently of outcomes.
suppressPackageStartupMessages(library(testthat))

test_that("per-event sampled sets give an unbiased score; the shared outcome-dependent set does not", {
  set.seed(123)
  p <- c(0.6, 0.3, 0.1); u <- log(p); N <- 100000; m <- 1
  score_set <- function(c, set) { s <- numeric(3); pr <- exp(u[set]) / sum(exp(u[set])); s[set] <- -pr; s[c] <- s[c] + 1; s }
  new <- matrix(0, N, 3); old <- matrix(0, N, 3)
  for (i in seq_len(N)) {
    c1 <- sample(3, 1, prob = p); c2 <- sample(3, 1, prob = p)          # two events in one state-year
    set1 <- c(c1, sample(setdiff(1:3, c1), 1)); set2 <- c(c2, sample(setdiff(1:3, c2), 1))
    new[i, ] <- score_set(c1, set1) + score_set(c2, set2)
    ch <- unique(c(c1, c2)); others <- setdiff(1:3, ch)
    shared <- c(ch, if (length(others)) sample(others, min(m, length(others))))
    old[i, ] <- score_set(c1, shared) + score_set(c2, shared)
  }
  expect_true(all(abs(colMeans(new)) < 0.01))         # expected score about 0 at the truth
  expect_gt(max(abs(colMeans(old))), 0.05)            # the old scheme is biased (about 0.12 on alternative 2)
})

test_that("IV estimate has the same expectation for initiating and non-initiating rows; the old scheme did not", {
  set.seed(123)
  w <- c(100, 1, 1); u <- log(w); N <- 100000; m <- 1
  new_iv <- function() log(3 / m) + log(sum(exp(u[sample(3, m)])))        # outcome-blind, same rule for every row
  old_iv_init <- function() { c <- sample(3, 1, prob = w / sum(w)); set <- c(c, sample(setdiff(1:3, c), 1)); log(3 / 2) + log(sum(exp(u[set]))) }
  iv_init <- replicate(N, new_iv()); iv_noninit <- replicate(N, new_iv()); iv_old_init <- replicate(N, old_iv_init())
  expect_lt(abs(mean(iv_init) - mean(iv_noninit)), 0.03)
  expect_gt(abs(mean(iv_old_init) - mean(iv_noninit)), 1)       # forced-in chosen target made IV outcome-dependent
  ## documented Jensen bias: E[log S_hat] < log S for the sampled estimate (both rows equally)
  expect_lt(mean(iv_init), log(sum(w)))
})
