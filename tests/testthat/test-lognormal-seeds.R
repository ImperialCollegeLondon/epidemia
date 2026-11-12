context("Log-normal prior for seeded infections")

data("EuropeCovid")
args <- list(data = EuropeCovid$data, rt = epirt(R(country, date) ~ 1 + lockdown),
  obs = suppressWarnings(epiobs(deaths ~ 1, i2o = EuropeCovid$inf2death * 0.02)),
  group_subset = c("Germany", "Italy"), refresh = 0, seed = 12345)

test_that("lognormal() is checked and reaches stan, one value per group", {
  expect_error(lognormal(0, 0), regexp = "positive")
  expect_error(lognormal("a"), regexp = "numeric")
  expect_error(epiinf(gen = 1, prior_seeds = lognormal()), NA)
  inf <- epiinf(gen = 1, prior_seeds = lognormal(log(c(5, 20)), 1))
  sdat <- standata_inf(inf, 2)
  expect_equal(sdat$prior_dist_for_seeds, 10)
  expect_equal(as.vector(sdat$prior_mean_for_seeds), log(c(5, 20)))
  expect_equal(as.vector(sdat$prior_scale_for_seeds), c(1, 1))
  expect_identical(sdat$hseeds, 0)
})

test_that("the seeds have the log-normal prior of their group", {
  run_args <- args
  run_args$inf <- epiinf(gen = EuropeCovid$si, prior_seeds = lognormal(log(c(5, 200)), 0.5))
  run_args$prior_PD <- TRUE
  run_args$iter <- 2000
  run_args$chains <- 1
  fm <- suppressWarnings(do.call(epim, run_args))
  m <- log(as.matrix(fm)[, c("seeds[Germany]", "seeds[Italy]")])
  expect_lt(max(abs(apply(m, 2, stats::median) - log(c(5, 200)))), 0.1)
  expect_lt(max(abs(apply(m, 2, stats::sd) - 0.5)), 0.1)
})
