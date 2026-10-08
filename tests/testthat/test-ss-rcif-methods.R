test_that("exponential decrement conserves probability with shared weights", {
  dat <- data.frame(
    id = c(1, 1, 2, 2, 3, 3, 4, 4),
    gap = c(1, 2, 1, 2, 1, 2, 2, 2),
    status = rep(c(1, 0), 4),
    cause = c("A", NA, "B", NA, "A", NA, "B", NA)
  )
  fit <- ss_rcif(dat, cif_method = "exponential_decrement")
  expect_equal(fit$cif_method, "exponential_decrement")
  expect_equal(fit$curve$F_A + fit$curve$F_B,
               fit$curve$F_overall, tolerance = 1e-14)
  expect_equal(fit$probability_conservation_max_abs, 0, tolerance = 1e-14)
  expect_true(all(diff(c(0, fit$curve$F_A)) >= 0))
  expect_true(all(diff(c(0, fit$curve$F_B)) >= 0))
  expect_gt(ss_rcif(dat)$probability_conservation_max_abs, 0)
})

test_that("cause relabeling preserves curves and singleton events", {
  dat <- data.frame(
    id = c(1, 1, 2, 2, 3, 3),
    gap = c(1, 2, 1, 2, 2, 2),
    status = rep(c(1, 0), 3),
    cause = c("A", NA, "B", NA, "A", NA)
  )
  a <- ss_rcif(dat, cif_method = "exponential_decrement")
  dat$cause[dat$cause == "A"] <- "X"
  dat$cause[dat$cause == "B"] <- "A"
  dat$cause[dat$cause == "X"] <- "B"
  b <- ss_rcif(dat, cif_method = "exponential_decrement")
  expect_equal(a$curve$F_A, b$curve$F_B)
  expect_equal(a$curve$F_B, b$curve$F_A)
  literal <- ss_rcif(dat, cause_weighting = "literal_eq15")
  expect_equal(literal$curve$F_A, rep(0, nrow(literal$curve)))
  expect_equal(literal$curve$F_B, rep(0, nrow(literal$curve)))
  dat$cause[which(dat$status == 1)[1]] <- NA
  expect_error(ss_rcif(dat), "cause must be non-missing")
})
