context("Discretizing continuous delay distributions")

test_that("discretize_cdf gives the mass on each whole day", {
  p <- discretize_cdf(pgamma, 30, shape = 2, rate = 0.5)
  expect_equal(length(p), 30)
  expect_equal(sum(p), 1)
  expect_equal(p, diff(pgamma(0:30, shape = 2, rate = 0.5)) / pgamma(30, shape = 2, rate = 0.5))
  expect_error(epiinf(gen = discretize_cdf(plnorm, 20, meanlog = 1, sdlog = 0.5)), NA)
})

test_that("discretize_cdf checks its arguments", {
  expect_error(discretize_cdf("pgamma", 10), regexp = "function")
  expect_error(discretize_cdf(pgamma, 0, shape = 2), regexp = "positive")
  expect_error(discretize_cdf(pgamma, 2.5, shape = 2), regexp = "integer")
  expect_error(discretize_cdf(function(x) rev(x), 5), regexp = "non-decreasing")
  expect_error(discretize_cdf(function(x) rep(1, length(x)), 5), regexp = "no mass")
})
