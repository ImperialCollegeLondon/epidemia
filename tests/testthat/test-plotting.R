context("Test plotting")

# load data
example.fit <- readRDS("../data/plot_test_fit.rds")

test_that("wrong group throws error", {
  expect_error(plot_rt(example.fit, group = c("Germany", "FakeCountry")))
  expect_error(plot_obs(example.fit, type = "deaths", group = c("Germany", "FakeCountry")))
  expect_error(plot_infections(example.fit, group = c("Germany", "FakeCountry")))
})

test_that("wrong observation type throws error", {
  expect_error(plot_obs(example.fit, type = "wrongtype", group = c("Germany")), regexp = "obs does not contain any observations for type 'wrongtype'")
})

test_that("levels out of [0,100] throws error", {
  expect_error(plot_rt(example.fit, levels = c(50, 95, 101)))
  expect_error(plot_obs(example.fit, type = "deaths", levels = c(50, 95, 101)))
  expect_error(plot_infections(example.fit, levels = c(50, 95, 101)))
  
  expect_error(plot_rt(example.fit, levels = c(-1, 50, 95)))
  expect_error(plot_obs(example.fit, type = "deaths", levels = c(-1, 50, 95)))
  expect_error(plot_infections(example.fit, levels = c(-1, 50, 95)))
})

test_that("non-logical log10_scale throws error", {
  expect_error(plot_rt(example.fit, log = "wrongtype"), regexp = "'log' must be of type logical")
})

