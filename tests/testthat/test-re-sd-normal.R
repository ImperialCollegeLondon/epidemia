context("Half-normal standard deviations of random effects")

data("EuropeCovid")
args <- list(data = EuropeCovid$data, inf = epiinf(gen = EuropeCovid$si),
  obs = suppressWarnings(epiobs(deaths ~ 1, i2o = EuropeCovid$inf2death * 0.02)),
  group_subset = c("Germany", "Italy"), refresh = 0, seed = 12345)

test_that("a normal prior_covariance reaches stan for one-coefficient terms only", {
  run_args <- args
  run_args$chains <- 0
  run_args$rt <- epirt(R(country, date) ~ 1 + lockdown + (0 + lockdown || country),
    prior_covariance = rstanarm::normal(0, 0.3))
  sdat <- suppressMessages(suppressWarnings(do.call(epim, run_args)))
  expect_identical(sdat$re_sd_normal, 1L)
  expect_equal(as.vector(sdat$scale), 0.3)

  run_args$rt <- epirt(R(country, date) ~ 1 + lockdown + (0 + lockdown || country))
  sdat <- suppressMessages(suppressWarnings(do.call(epim, run_args)))
  expect_identical(sdat$re_sd_normal, 0L)

  run_args$rt <- epirt(R(country, date) ~ 1 + (1 + lockdown | country), prior_covariance = rstanarm::normal(0, 0.3))
  expect_error(suppressWarnings(do.call(epim, run_args)), regexp = "one coefficient")
  run_args$rt <- epirt(R(country, date) ~ 1 + (0 + lockdown | country), prior_covariance = rstanarm::normal(1, 0.3))
  expect_error(suppressWarnings(do.call(epim, run_args)), regexp = "location 0")
})

test_that("the standard deviation has a half-normal prior", {
  run_args <- args
  run_args$rt <- epirt(R(country, date) ~ 1 + (0 + lockdown || country), prior_covariance = rstanarm::normal(0, 0.3))
  run_args$prior_PD <- TRUE
  run_args$iter <- 2000
  run_args$chains <- 1
  fm <- suppressWarnings(do.call(epim, run_args))
  s <- sqrt(as.matrix(fm)[, "R|Sigma[country:lockdown,lockdown]"])
  expect_lt(abs(stats::median(s) - 0.3 * stats::qnorm(0.75)), 0.05)
  expect_output(print(prior_summary(fm)), "half-normal standard deviations")
})
