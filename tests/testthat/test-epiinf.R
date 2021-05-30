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
  expect_error(inf <- epiinf(gen = rep(0.2,5), prior_susc0 = "dummy"), regexp = "prior")
  expect_error(inf <- epiinf(gen = rep(0.2,5), prior_susc0 = rstanarm::lasso()), regexp = "normal")
  expect_error(inf <- epiinf(gen = rep(0.2,5), prior_vnoise = "dummy"), regexp = "prior")
  expect_error(inf <- epiinf(gen = rep(0.2,5), prior_vnoise = rstanarm::lasso()), regexp = "normal")
})

test_that("latent and pop_adjust are logical scalars", {
  expect_error(inf <- epiinf(gen = rep(0.2, 5), latent = 1), regexp = "logical")
  expect_error(inf <- epiinf(gen = rep(0.2, 5), latent = c(TRUE, TRUE)), regexp = "scalar")
  expect_error(inf <- epiinf(gen = rep(0.2, 5), pop_adjust = 1), regexp = "logical")
  expect_error(inf <- epiinf(gen = rep(0.2, 5), pop_adjust = c(TRUE, TRUE)), regexp = "scalar")
})

test_that("populations and susceptibles behave correctly", {
  inf <- epiinf(gen = rep(0.2, 5), pop_adjust = FALSE)
  expect_equal(inf$pops, NULL)
  expect_equal(inf$vacc, NULL)
  
  inf <- epiinf(gen = rep(0.2, 5), pop_adjust = FALSE, pops = dummy1, vacc = dummy2)
  expect_equal(inf$pops, NULL)
  expect_equal(inf$vacc, NULL)
  
  expect_error(epiinf(gen = rep(0.2, 5), pop_adjust = TRUE), regexp = "pop")
  
  inf <- epiinf(gen = rep(0.2, 5), pop_adjust = TRUE, pop = dummy)
  expect_equal(inf$pops, "dummy")
  expect_equal(inf$vacc, NULL)
  
  inf <- epiinf(gen = rep(0.2, 5), pop_adjust = TRUE, pop = dummy1, vacc = dummy2)
  expect_equal(inf$pops, "dummy1")
  expect_equal(inf$vacc, "dummy2") 
  
  inf <- epiinf(gen = rep(0.2, 5), pop_adjust = TRUE, pop = "dummy1", vacc = "dummy2")
  expect_equal(inf$pops, "dummy1")
  expect_equal(inf$vacc, "dummy2")
})

test_that("family is scalar character and in required set", {
  expect_error(inf <- epiinf(gen = rep(0.2, 5), family = na.action), regexp = "character")
  expect_error(inf <- epiinf(gen = rep(0.2, 5), family = c("a", "b")), regexp = "scalar")
  expect_error(inf <- epiinf(gen = rep(0.2, 5), family = "log-normal"), regexp = "normal")
})



