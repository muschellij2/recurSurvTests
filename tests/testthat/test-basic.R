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
