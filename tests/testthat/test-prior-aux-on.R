context("Priors on one over the square root of the auxiliary parameter")

data("EuropeCovid")
args <- list(data = EuropeCovid$data, rt = epirt(R(country, date) ~ 1 + lockdown),
  inf = epiinf(gen = EuropeCovid$si), group_subset = "Italy", refresh = 0, seed = 12345)

test_that("prior_aux_on is checked and passed to stan", {
  expect_identical(epiobs(deaths ~ 1, i2o = 1)$prior_aux_on, "aux")
  expect_error(epiobs(deaths ~ 1, i2o = 1, prior_aux_on = "sqrt"), regexp = "inv_sqrt")
  expect_error(epiobs(deaths ~ 1, i2o = 1, prior_aux_on = c("aux", "aux")), regexp = "scalar")
  run_args <- args
  run_args$chains <- 0
  run_args$obs <- suppressWarnings(epiobs(deaths ~ 1, i2o = EuropeCovid$inf2death * 0.02, aux_by_group = TRUE,
    prior_aux = rstanarm::normal(0, 1), prior_aux_on = "inv_sqrt"))
  sdat <- suppressMessages(suppressWarnings(do.call(epim, run_args)))
  expect_equal(as.vector(sdat$oaux_inv_sqrt), 1L)
  run_args$obs$prior_aux_on <- "aux"
  sdat <- suppressMessages(suppressWarnings(do.call(epim, run_args)))
  expect_equal(as.vector(sdat$oaux_inv_sqrt), 0L)
})

test_that("normal(0, 1) on one over the square root gives it a half-normal prior", {
  run_args <- args
  run_args$prior_PD <- TRUE
  run_args$iter <- 2000
  run_args$chains <- 1
  run_args$obs <- suppressWarnings(epiobs(deaths ~ 1, i2o = EuropeCovid$inf2death * 0.02,
    prior_aux = rstanarm::normal(0, 1), prior_aux_on = "inv_sqrt"))
  fm <- suppressWarnings(do.call(epim, run_args))
  x <- 1 / sqrt(as.matrix(fm)[, "deaths|reciprocal dispersion"])
  expect_lt(abs(stats::median(x) - stats::qnorm(0.75)), 0.1)
  expect_lt(abs(stats::quantile(x, 0.9) - stats::qnorm(0.95)), 0.2)
})
