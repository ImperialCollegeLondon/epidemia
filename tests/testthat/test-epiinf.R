context("Test error handling of epiinf")

library(epidemia)
library(testthat)

test_that("gen must be a numeric simplex vector", {
  expect_error(inf <- epiinf(gen = rep(0.2,5)), NA)
  expect_error(inf <- epiinf(gen = 1), NA)
  expect_error(inf <- epiinf(gen = "dummy"), regexp = "numeric")
  expect_error(inf <- epiinf(gen = numeric()), regexp = "sum")
  expect_error(inf <- epiinf(gen = -1), regexp = "non-negative")
  expect_error(inf <- epiinf(gen = rep(1,5)), regexp = "sum")
})

test_that("seed_days is positive, scalar integer", {
  expect_error(inf <- epiinf(gen = rep(0.2,5), seed_days = 5), NA)
  expect_error(inf <- epiinf(gen = rep(0.2,5), seed_days = 5.1), regexp = "integer")
  expect_error(inf <- epiinf(gen = rep(0.2,5), seed_days = "dummy"), regexp = "numeric")
  expect_error(inf <- epiinf(gen = rep(0.2,5), seed_days = c(1,3)), regexp = "scalar")
  expect_error(inf <- epiinf(gen = rep(0.2,5), seed_days = integer()), regexp = "scalar")
  expect_error(inf <- epiinf(gen = rep(0.2,5), seed_days = 0), regexp = "positive")
})

test_that("Correct priors are enforced", {
  expect_error(inf <- epiinf(gen = rep(0.2,5), prior_tau = "dummy"), regexp = "prior")
  expect_error(inf <- epiinf(gen = rep(0.2,5), prior_tau = rstanarm::normal()), regexp = "exponential")
  expect_error(inf <- epiinf(gen = rep(0.2,5), prior_aux = "dummy"), regexp = "prior")
  expect_error(inf <- epiinf(gen = rep(0.2,5), prior_aux = rstanarm::lasso()), regexp = "normal")
})

test_that("latent and pop_adjust are logical scalars", {
  expect_error(inf <- epiinf(gen = rep(0.2, 5), latent = 1), regexp = "logical")
  expect_error(inf <- epiinf(gen = rep(0.2, 5), latent = c(TRUE, TRUE)), regexp = "scalar")
  expect_error(inf <- epiinf(gen = rep(0.2, 5), pop_adjust = 1), regexp = "logical")
  expect_error(inf <- epiinf(gen = rep(0.2, 5), pop_adjust = c(TRUE, TRUE)), regexp = "scalar")
})

test_that("family is scalar character and in required set", {
  expect_error(inf <- epiinf(gen = rep(0.2, 5), family = na.action), regexp = "character")
  expect_error(inf <- epiinf(gen = rep(0.2, 5), family = c("a", "b")), regexp = "scalar")
  expect_error(inf <- epiinf(gen = rep(0.2, 5), family = "normal"), regexp = "log-normal")
})

test_that("susceptibles converted into character appropriately", {
  inf <- epiinf(gen = rep(0.2, 5), pop_adjust = FALSE, susceptibles = column)
  expect_equal(inf$susceptibles, NULL)
  inf <- epiinf(gen = rep(0.2, 5), pop_adjust = TRUE, susceptibles = column)
  expect_equal(inf$susceptibles, "column")
  inf <- epiinf(gen = rep(0.2, 5), pop_adjust = TRUE, susceptibles = "column")
  expect_equal(inf$susceptibles, "column")
})



