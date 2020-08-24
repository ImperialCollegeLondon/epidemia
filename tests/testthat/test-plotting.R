context("Test plotting")

# load data
example.fit <- readRDS("../data/plot_test_fit.rds")
fm <- example.fit

test_that("wrong group throws error", {
  expect_error(plot_rt(example.fit, groups = c("Germany", "FakeCountry")), regexp = "group\\(s\\) FakeCountry not found.")
  expect_error(plot_obs(example.fit, type = "deaths", groups = c("Germany", "FakeCountry")), regexp = "group\\(s\\) FakeCountry not found.")
  expect_error(plot_infections(example.fit, groups = c("Germany", "FakeCountry")), regexp = "group\\(s\\) FakeCountry not found.")
})

test_that("wrong observation type throws error", {
  expect_error(plot_obs(example.fit, type = "wrongtype", groups = c("Germany")), regexp = "obs does not contain any observations")
})

test_that("missing observation type throws error", {
  expect_error(plot_obs(example.fit, groups = c("Germany")), regexp = "^argument")
})

test_that("levels out of [0,100] throws error", {
  expect_error(plot_rt(example.fit, levels = c(50, 95, 101)), regexp = "all levels must be between")
  expect_error(plot_obs(example.fit, type = "deaths", levels = c(50, 95, 101)), regexp = "all levels must be between ")
  expect_error(plot_infections(example.fit, levels = c(50, 95, 101)), regexp = "all levels must be between")
  
  expect_error(plot_rt(example.fit, levels = c(-1, 50, 95)), regexp = "all levels must be between")
  expect_error(plot_obs(example.fit, type = "deaths", levels = c(-1, 50, 95)), regexp = "all levels must be between")
  expect_error(plot_infections(example.fit, levels = c(-1, 50, 95)), regexp = "all levels must be between ")
})


test_that("Rt smoothing", {
  expect_warning(plot_rt(example.fit, smooth=-1), regexp = "smooth must be a positive integer")
  expect_warning(plot_rt(example.fit, smooth=0.5), regexp = "smooth must be a positive integer")
  expect_warning(plot_rt(example.fit, smooth=-0.5), regexp = "smooth must be a positive integer")
  expect_warning(plot_rt(example.fit, smooth=1000), regexp = "smooth=1000 is too large")
})

test_that("date subsetting",{
  expect_error(plot_rt(example.fit, dates=c("2020-05-04", "2020-04-04")),
                 regexp = "end date must be after")
  expect_error(plot_rt(example.fit, dates=c("2020-05-04", "2020-05-04")),
                 regexp = "end date must be after")
  expect_error(plot_rt(example.fit, dates=c("2020-05-04", "2020-05-40")),
                 regexp = "conversion of 'dates' to ")
  expect_error(plot_rt(example.fit, dates=c("2019-05-04", "2019-05-20")),
                 regexp = "date subsetting removed all data")

  expect_error(plot_rt(example.fit,  dates=c("2020-05-04", "2020-05-14"), date_format="asdf"),
                 regexp = "conversion of 'dates' to ")
  
})

test_that(".check_dates works as expected", {
  date_format <- "%Y-%m-%d"
  expect_equal(check_dates(c("2020-03-09", "2020-04-06"), date_format, "2020-02-06", "2020-06-06"),
               as.Date(c("2020-03-09","2020-04-06")))
  expect_equal(check_dates(c("2020-03-09", NA), date_format, "2020-02-06", "2020-06-06"),
               as.Date(c("2020-03-09","2020-06-06")))
  expect_equal(check_dates(c(NA, "2020-04-06"), date_format, "2020-02-06", "2020-06-06"),
               as.Date(c("2020-02-06","2020-04-06")))
  expect_equal(check_dates(as.Date(c("2020-03-09", "2020-04-06")), date_format, "2020-02-06", "2020-06-06"),
               as.Date(c("2020-03-09","2020-04-06")))
  expect_equal(check_dates(as.Date(c("2020-03-09", NA)), date_format, "2020-02-06", "2020-06-06"),
               as.Date(c("2020-03-09","2020-06-06")))
  expect_equal(check_dates(as.Date(c(NA, "2020-04-06")), date_format, "2020-02-06", "2020-06-06"),
               as.Date(c("2020-02-06","2020-04-06")))
}
)


test_that("plot_rt runs through with various arguments", {
  fun <- plot_rt
  expect_true(inherits(fun(fm), "ggplot"))
  expect_true(inherits(fun(fm, log=TRUE), "ggplot"))
  expect_true(inherits(fun(fm, plotly=T), "plotly"))
  expect_true(inherits(fun(fm, smooth=7), "ggplot"))
})

test_that("plot_obs runs through with varioud arguments", {
fun <- function(x, ...) plot_obs(x, type="deaths", ...)
expect_true(inherits(fun(fm), "ggplot"))
expect_true(inherits(fun(fm, log=TRUE), "ggplot"))
expect_true(inherits(fun(fm, plotly=T), "plotly"))
expect_true(inherits(fun(fm, cumulative=T), "ggplot"))
})

test_that("plot_infections runs through with varioud arguments", {
fun <- plot_infections
expect_true(inherits(fun(fm), "ggplot"))
expect_true(inherits(fun(fm, log=TRUE), "ggplot"))
expect_true(inherits(fun(fm, plotly=T), "plotly"))
expect_true(inherits(fun(fm, cumulative=T), "ggplot"))
})

test_that("plot_infectious runs through with varioud arguments", {
fun <- plot_infectious
expect_true(inherits(fun(fm), "ggplot"))
expect_true(inherits(fun(fm, log=TRUE), "ggplot"))
expect_true(inherits(fun(fm, plotly=T), "plotly"))
expect_true(inherits(fun(fm, cumulative=T), "ggplot"))
})

