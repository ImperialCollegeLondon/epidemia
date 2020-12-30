context("Test models with autocorrelation terms")

library(dplyr)
data("EuropeCovid")

args <- list()
args$data <- EuropeCovid$data %>%
  filter(country == "United_Kingdom") %>%
  mutate(
    week = lubridate::week(date),
    cases = deaths,
    dummy = cases
  )

args$rt <- epirt(R(country, date) ~ rw(time=week))
args$inf <- epiinf(gen = EuropeCovid$si)
args$sampling_args <- list(iter=500, chains=0)

deaths <- epiobs(
  formula = deaths ~ rw(time=week),
  i2o = EuropeCovid$inf2death * 0.02,
  na.action = na.pass
)

dummy <- epiobs(
  formula = dummy ~ 1,
  i2o = EuropeCovid$inf2death * 0.02,
  na.action = na.pass
)

cases <- epiobs(
  formula = cases ~ rw(time=week),
  i2o = EuropeCovid$inf2death * 0.02,
  na.action = na.pass
)

test_that("No rw terms leads to correct default stan data", {
  # no rw terms
  args$obs <- list(dummy)
  sdat <- do.call(epim, args)

  expect_equal(sdat$obs_ac_nnz, 0)
  expect_equal(sdat$obs_ac_nproc, 0)
  expect_equal(sdat$obs_ac_nterms, 0)
  expect_equal(sdat$obs_ac_q, 0)
  expect_equal(sdat$obs_ac_ntime, numeric())
  expect_equal(sdat$obs_ac_prior_scales, numeric(0))
  expect_equal(sdat$obs_ac_v, numeric(0))
})

test_that("Expected standata for autocor with single observation type", {
  # a single observation type, no NAs
  args$obs <- list(deaths)
  sdat <- do.call(epim, args)

  ntime <- length(unique(args$data$week))
  expect_equal(sdat$obs_ac_nnz, nrow(args$data))
  expect_equal(sdat$obs_ac_nproc, 1)
  expect_equal(sdat$obs_ac_nterms, 1)
  expect_equal(sdat$obs_ac_q, ntime)
  expect_equal(as.numeric(sdat$obs_ac_ntime), ntime)
  expect_equal(as.numeric(sdat$obs_ac_prior_scales), 0.2)
  expect_equal(min(sdat$obs_ac_v), 0)
  expect_equal(max(sdat$obs_ac_v), ntime-1)
})

test_that("Expected standata for autocor with multiple observation types", {
  # multiple types, and one term with no NAs
  args$obs <- list(deaths, dummy, cases)
  sdat <- do.call(epim, args)

  ntime <- length(unique(args$data$week))
  expect_equal(sdat$obs_ac_nnz, 2 * nrow(args$data))
  expect_equal(sdat$obs_ac_q, 2 * ntime)
  expect_equal(as.numeric(sdat$obs_ac_ntime), rep(ntime,2))
  expect_equal(min(sdat$obs_ac_v), 0)
  expect_equal(max(sdat$obs_ac_v), 2 * ntime-1)
})

test_that("Correct handling of NA terms in the time vector", {
  # put NAs in the time vector and check again
  args$obs <- list(deaths, dummy, cases)
  args$data$week <- as.integer(args$data$week)
  w <- args$data$week < 12
  args$data$week[w] <- NA
  sdat <- do.call(epim, args)

  ntime <- length(unique(args$data$week)) - 1
  expect_equal(sdat$obs_ac_nnz, 2 * nrow(args$data))
  expect_equal(sdat$obs_ac_nproc, 2)
  expect_equal(sdat$obs_ac_nterms, 2)
  expect_equal(sdat$obs_ac_q, 2 * ntime)
  expect_equal(as.numeric(sdat$obs_ac_ntime), rep(ntime,2))
  expect_equal(as.numeric(sdat$obs_ac_prior_scales), rep(0.2,2))
  expect_equal(min(sdat$obs_ac_v), -1)
  expect_equal(max(sdat$obs_ac_v), 2 * ntime-1)
})
