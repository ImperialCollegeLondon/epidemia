context("Test all make name functions in epim")

data("EuropeCovid")
data <- EuropeCovid$data
formula <- R(country, date) ~ 1 + rw()

test_that("make_rw_nms", {
  dates <- unique(data$date)
  nms <- paste0("rw()[", dates, ",all]")
  expect_equal(make_rw_nms(formula, data), nms)
})

test_that("make_rw_sigma_nms", {
  formula <- R(country, date) ~ 1 + rw() + rw(gr=country)
  nms <- c(
    "R|sigma:rw()[all]", 
    paste0("R|sigma:rw(gr = country)[", unique(data$country), "]")
  )
  expect_equal(make_rw_sigma_nms(formula, data), nms)
})

test_that("make_oaux_nms", {
  deaths <- epiobs(
    formula = deaths(country,date) ~ 1,
    family = "poisson",
    lag = rep(1/7,7)
  )
  
  cases <- epiobs(
    formula = cases(country, date) ~ 1,
    family = "neg_binom",
    lag = rep(1/7,7)
  )
  
  obs <- list(deaths=deaths, cases=cases)
  expect_equal(make_oaux_nms(obs), "cases|recipricol dispersion")
})
