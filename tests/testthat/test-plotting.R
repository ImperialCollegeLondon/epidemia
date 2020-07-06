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

test_that("missing observation type throws error", {
  expect_error(plot_obs(example.fit, group = c("Germany")), regexp = "must specify an observation type")
})

test_that("levels out of [0,100] throws error", {
  expect_error(plot_rt(example.fit, levels = c(50, 95, 101)))
  expect_error(plot_obs(example.fit, type = "deaths", levels = c(50, 95, 101)))
  expect_error(plot_infections(example.fit, levels = c(50, 95, 101)))
  
  expect_error(plot_rt(example.fit, levels = c(-1, 50, 95)))
  expect_error(plot_obs(example.fit, type = "deaths", levels = c(-1, 50, 95)))
  expect_error(plot_infections(example.fit, levels = c(-1, 50, 95)))
})

test_that("non-logical log throws error", {
  expect_error(plot_rt(example.fit, log = "wrongtype"), regexp = "'log' must be of type logical")
  expect_error(plot_obs(example.fit, type="deaths", log = "wrongtype"), regexp = "'log' must be of type logical")
  expect_error(plot_infections(example.fit, log = "wrongtype"), regexp = "'log' must be of type logical")
  
})

test_that("Rt smoothing", {
  expect_warning(plot_rt(example.fit, smooth=-1), regexp = "smooth must be a positive integer - no smoothing will be performed")
  expect_warning(plot_rt(example.fit, smooth=0.5), regexp = "smooth must be a positive integer - no smoothing will be performed")
  expect_warning(plot_rt(example.fit, smooth=-0.5), regexp = "smooth must be a positive integer - no smoothing will be performed")
  expect_warning(plot_rt(example.fit, smooth=1000))
})

test_that("date subsetting",{
  expect_warning(plot_rt(example.fit, dates=c("2020-05-04", "2020-04-04")),
                 regexp = "start of date range is before end - reversing dates")
  expect_warning(plot_rt(example.fit, dates=c("2020-05-04", "2020-05-04")),
                 regexp = "dates must be different - plotting the entire range")
  expect_warning(plot_rt(example.fit, dates=c("2020-05-04", "2020-05-40")),
                 regexp = "Could not coerce 2020-05-40 to date with specified format - plotting the enire date range")
  expect_error(plot_rt(example.fit, dates=c("2019-05-04", "2019-05-20")),
                 regexp = "date subsetting removed all data")

  expect_warning(plot_infections(example.fit, dates=c("2020-05-04", "2020-04-04")),
                 regexp = "start of date range is before end - reversing dates")
  expect_warning(plot_infections(example.fit, dates=c("2020-05-04", "2020-05-04")),
                 regexp = "dates must be different - plotting the entire range")
  expect_warning(plot_infections(example.fit, dates=c("2020-05-04", "2020-05-40")),
                 regexp = "Could not coerce 2020-05-40 to date with specified format - plotting the enire date range")
  expect_error(plot_infections(example.fit, dates=c("2019-05-04", "2019-05-20")),
               regexp = "date subsetting removed all data")
  
  expect_warning(plot_obs(example.fit, type="deaths", dates=c("2020-05-04", "2020-04-04")),
                 regexp = "start of date range is before end - reversing dates")
  expect_warning(plot_obs(example.fit, type="deaths", dates=c("2020-05-04", "2020-05-04")),
                 regexp = "dates must be different - plotting the entire range")
  expect_warning(plot_obs(example.fit, type="deaths", dates=c("2020-05-04", "2020-05-40")),
                 regexp = "Could not coerce 2020-05-40 to date with specified format - plotting the enire date range")
  expect_error(plot_obs(example.fit, type="deaths", dates=c("2019-05-04", "2019-05-20")),
               regexp = "date subsetting removed all data")
  
  expect_warning(plot_rt(example.fit,  dates=c("2020-05-04", "2020-05-14"), date_format="asdf"),
                 regexp = "Could not coerce 2020-05-04, 2020-05-14 to date with specified format - plotting the enire date range")
  
})