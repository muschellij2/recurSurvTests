test_that("paper log-logistic gaps follow Zhao Table 1 survival", {
  set.seed(230923)
  gaps <- replicate(20000, recurSurvTests:::.zhao_draw_gap("loglogistic"))
  times <- exp(4) * c(0.25, 0.5, 1, 2, 4)
  expected_survival <- 1 / (1 + (exp(-4) * times)^2)
  observed_survival <- vapply(times, function(t) mean(gaps > t), numeric(1))
  expect_equal(observed_survival, expected_survival, tolerance = 0.015)
  # Log gaps have a logistic distribution with location 4 and scale 0.5.
  expect_equal(mean(log(gaps)), 4, tolerance = 0.02)
  expect_equal(sd(log(gaps)), pi / sqrt(12), tolerance = 0.02)
})

test_that("paper log-logistic histories have plausible event counts", {
  set.seed(230923)
  dat <- recurSurvTests:::.zhao_simulated_gap_data(
    100, "loglogistic", c(0.01, 1.99), 180
  )
  expect_gt(sum(dat$status), 100)
  expect_lt(sum(dat$status), 2000)
  fit <- zhao_rank_test(dat, "group", episode = "episode")
  expect_true(is.finite(fit$z))
})

test_that("Zhao rank-test variance implementations are explicit", {
  dat <- data.frame(
    id = rep(1:4, each = 3), episode = rep(1:3, 4),
    gap = c(2, 3, 4, 1, 4, 3, 2, 2, 5, 1, 3, 4),
    status = rep(c(1, 1, 0), 4),
    arm = rep(c("A", "A", "B", "B"), each = 3)
  )
  literal <- zhao_rank_test(dat, "arm", episode = "episode",
    variance_method = "zhao_eq6")
  pooled <- zhao_rank_test(dat, "arm", episode = "episode",
    variance_method = "pooled_risk")

  expect_true(is.finite(literal$z))
  expect_equal(pooled$W, sum(pooled$subject_residuals$residual), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(literal$variance, pooled$variance)))
})

test_that("Zhao null simulation has approximately nominal calibration", {
  sim <- zhao_rank_simulation(
    n_sim = 100, n_per_group = 40, seed = 101,
    variance_method = "pooled_risk"
  )
  expect_lt(abs(sim$empirical_type1_error - 0.05), 0.05)
  expect_gt(sim$sd_z, 0.7)
  expect_lt(sim$sd_z, 1.3)
})

test_that("PSG scenario targets roughly 200 completed events", {
  sim <- zhao_rank_simulation(
    n_sim = 5, n_per_group = 4, scenario = "psg_200", seed = 4,
    variance_method = "pooled_risk"
  )
  expect_gt(sim$mean_completed_events, 150)
  expect_lt(sim$mean_completed_events, 250)
  expect_equal(sim$settings$followup, 480)
})
