test_that("wc_surv returns a nonincreasing survival curve", {
  dat <- data.frame(
    id = rep(1:3, c(3, 2, 3)),
    gap = c(2, 3, 4, 1, 5, 2, 4, 3),
    status = c(1, 1, 0, 1, 0, 1, 1, 0),
    episode = c(1, 2, 3, 1, 2, 1, 2, 3)
  )
  fit <- wc_surv(dat, episode = "episode")
  expect_true(all(diff(fit$surv) <= sqrt(.Machine$double.eps)))
  expect_true(all(fit$surv >= 0 & fit$surv <= 1))
})

test_that("psh_surv is the pooled Kaplan-Meier product limit", {
  dat <- data.frame(
    id = rep(1:3, c(3, 2, 3)),
    gap = c(2, 3, 4, 1, 5, 2, 4, 3),
    status = c(1, 1, 0, 1, 0, 1, 1, 0),
    episode = c(1, 2, 3, 1, 2, 1, 2, 3)
  )
  fit <- psh_surv(dat, episode = "episode")
  expect_equal(fit$risk, c(8, 7, 5, 3))
  expect_equal(fit$dN, c(1, 2, 1, 1))
  expect_equal(fit$surv, cumprod(1 - fit$dN / fit$risk))
})

test_that("paired permutation swaps labels within pairs", {
  dat <- data.frame(
    id = rep(1:4, each = 2), pair = rep(1:2, each = 4),
    episode = rep(1:2, 4), gap = c(2, 4, 3, 2, 1, 5, 2, 3),
    status = rep(c(1, 0), 4), visit = rep(c("A", "B"), each = 2)
  )
  fit <- wc_paired_permutation(dat, pair = "pair", group = "visit",
    episode = "episode", B = 19, seed = 1)
  expect_equal(length(fit$permutation_z), 19)
  expect_gte(fit$p.value, 0)
  expect_lte(fit$p.value, 1)
})

test_that("wc_logrank returns a finite p-value on simple two-group data", {
  dat <- data.frame(
    id = rep(1:4, each = 3),
    episode = rep(1:3, 4),
    gap = c(2, 3, 4, 1, 4, 3, 2, 2, 5, 1, 3, 4),
    status = rep(c(1, 1, 0), 4),
    arm = rep(c("A", "A", "B", "B"), each = 3)
  )
  ans <- wc_logrank(dat, group = "arm", episode = "episode")
  expect_true(is.finite(ans$p.value))
  expect_gte(ans$p.value, 0)
  expect_lte(ans$p.value, 1)
})

test_that("shared RCIF increments are additive across causes", {
  dat <- data.frame(
    id = rep(1:4, each = 3), episode = rep(1:3, 4),
    gap = c(2, 3, 4, 1, 4, 3, 2, 2, 5, 1, 3, 4),
    status = rep(c(1, 1, 0), 4),
    cause = c("A", "B", NA, "A", "A", NA, "B", "A", NA, "B", "B", NA)
  )
  fit <- ss_rcif(dat, episode = "episode")
  expect_equal(fit$cause_weighting, "shared")
  expect_equal(fit$additivity_max_abs_dG, 0, tolerance = 1e-12)
  expect_true(all(diff(fit$curve$F_overall) >= -sqrt(.Machine$double.eps)))
})

test_that("literal Equation 15 weighting remains an explicit sensitivity option", {
  dat <- data.frame(
    id = rep(1:3, each = 4), episode = rep(1:4, 3),
    gap = c(1, 2, 3, 4, 1, 3, 2, 4, 2, 1, 3, 4),
    status = rep(c(1, 1, 1, 0), 3),
    cause = c("A", "B", "A", NA, "B", "B", "A", NA, "A", "A", "B", NA)
  )
  fit <- ss_rcif(dat, episode = "episode", cause_weighting = "literal_eq15")
  expect_equal(fit$cause_weighting, "literal_eq15")
  expect_match(fit$warning, "does not add")
})

test_that("equal-cause RCIF test uses subject-level bootstrap", {
  dat <- data.frame(
    id = rep(1:4, each = 3), episode = rep(1:3, 4),
    gap = c(2, 3, 4, 1, 4, 3, 2, 2, 5, 1, 3, 4),
    status = rep(c(1, 1, 0), 4),
    cause = c("A", "B", NA, "A", "A", NA, "B", "A", NA, "B", "B", NA)
  )
  expect_error(ss_rcif_equal_causes_test(dat, episode = "episode", B = 19), "B >= 20")
  fit <- ss_rcif_equal_causes_test(dat, episode = "episode", B = 100, seed = 1)
  expect_true(is.finite(fit$chisq))
  expect_gte(fit$bootstrap_success, 20)
  expect_false(fit$logrank_analogue)
})

test_that("input validation catches malformed subject histories", {
  dat <- data.frame(id = c(1, 1), gap = c(1, 2), status = c(0, 1))
  expect_error(wc_surv(dat), "final row must have status=0")
  valid_dat <- data.frame(id = c(1, 1), gap = c(1, 2), status = c(1, 0))
  expect_error(wc_logrank(valid_dat, group = "arm"), "Missing group column")
})
