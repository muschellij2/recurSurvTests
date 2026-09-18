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

test_that("recurrent survival confidence intervals resample whole subjects", {
  dat <- data.frame(
    id = rep(1:8, each = 3), episode = rep(1:3, 8),
    gap = rep(c(1, 2, 4), 8) + rep(0:7, each = 3) / 10,
    status = rep(c(1, 1, 0), 8)
  )
  expect_error(recurrent_surv_ci(dat, episode = "episode", B = 99), "B >= 100")
  point <- recurrent_surv_ci(
    dat, estimator = "wc", episode = "episode", times = c(1, 2),
    B = 100, seed = 1, keep_bootstrap = TRUE
  )
  band <- recurrent_surv_ci(
    dat, estimator = "psh", episode = "episode", times = c(1, 2),
    interval = "simultaneous", B = 100, seed = 1
  )
  normal <- recurrent_surv_ci(
    dat, estimator = "psh", episode = "episode", times = c(1, 2),
    ci_method = "normal"
  )
  expect_equal(nrow(point$bootstrap_surv), 100)
  expect_equal(nrow(point$curve), 2)
  expect_true(all(point$curve$lower >= 0 & point$curve$upper <= 1))
  expect_true(all(band$curve$lower <= band$curve$surv))
  expect_true(all(band$curve$upper >= band$curve$surv))
  expect_true(is.finite(band$critical.value))
  manual_se <- sqrt(normal$fit$surv^2 * cumsum(
    normal$fit$dN / normal$fit$risk^2
  ))
  expect_equal(normal$curve$std.error, manual_se[match(c(1, 2), normal$fit$time)])
  expect_error(
    recurrent_surv_ci(dat, episode = "episode", interval = "simultaneous",
                      ci_method = "normal"),
    "require ci_method"
  )
})

test_that("recurrent survival quantile intervals resample whole subject histories", {
  dat <- data.frame(
    id = rep(1:8, each = 3), episode = rep(1:3, 8),
    gap = rep(c(1, 2, 4), 8) + rep(0:7, each = 3) / 10,
    status = rep(c(1, 1, 0), 8)
  )
  boot <- recurrent_surv_quantile_ci(
    dat, estimator = "psh", episode = "episode", B = 100, seed = 1,
    keep_bootstrap = TRUE
  )
  expect_equal(boot$estimate, recurSurvTests:::.recurrent_surv_quantile(
    boot$fit, 0.5
  ))
  expect_equal(length(boot$bootstrap_quantile), 100)
  expect_lte(boot$lower, boot$upper)
  expect_true(is.finite(boot$estimate))
})

test_that("curve contrasts, restricted means, plots, and K-group tests work", {
  dat <- data.frame(id=rep(1:9, each=3), episode=rep(1:3,9),
    gap=rep(c(1,2,4),9)+rep(0:8,each=3)/10, status=rep(c(1,1,0),9),
    arm=rep(c("A","B","C"), each=9))
  ct <- recurrent_surv_contrast_ci(dat[dat$arm != "C",], group="arm", episode="episode", B=100, seed=1)
  rm <- recurrent_surv_rmst_ci(dat, episode="episode", tau=3, B=100, seed=1)
  kk <- wc_logrank_k(dat, group="arm", episode="episode")
  expect_equal(nrow(ct$curve), length(unique(ct$curve$time)))
  expect_true(is.finite(rm$estimate)); expect_gte(kk$df, 1)
  expect_silent(plot_recurrent_surv_ci(recurrent_surv_ci(dat, episode="episode", B=100, seed=1)))
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
