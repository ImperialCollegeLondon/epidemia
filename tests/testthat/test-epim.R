context("Test epim runs through with varioud different specifications")

data("EuropeCovid")
args <- list()
args$data = EuropeCovid$data
args$inf <- epiinf(gen = EuropeCovid$si)
args$rt <- epirt(R(country, date) ~ 1 + lockdown)
expect_warning(args$obs <- epiobs(deaths~1, i2o = EuropeCovid$inf2death * 0.02))
args$group_subset <- c("Germany", "Italy")
args <- c(args, list(iter=10,chains=1, seed=12345))
args$refresh <- 0

test_that("epim runs through with various rt formula", {
  run_args <- args

  # just fixed effects
  expect_warning(fm <- do.call(epim, run_args))
  expect_true(inherits(fm, "epimodel"))

  # random effects
  run_args$rt <- epirt(formula = R(country, date) ~ (lockdown | country))
  expect_warning(fm <- do.call(epim, run_args))
  expect_true(inherits(fm, "epimodel"))

  # random walks
  run_args$data$week <- format(run_args$data$date,"%V")
  run_args$rt <- epirt(formula = R(country, date) ~ (lockdown | country) + rw(time=week) + rw(time=week, gr=country))
  expect_warning(fm <- do.call(epim, run_args))
  expect_true(inherits(fm, "epimodel"))

})

test_that("epim runs through with different algorithms", {
  run_args <- args
  run_args$algorithm <- "sampling"
  run_args$iter <- 10
  run_args$seed <- 12345
  run_args$chains <- 1 
  expect_warning(fm <- do.call(epim, run_args))
  expect_true(inherits(fm, "epimodel"))
})

test_that("epim works with init_run", {
  run_args <- args
  run_args$init_run <- TRUE
  run_args$algorithm <- "sampling"
  run_args$iter <- 10
  run_args$seed <- 12345
  run_args$chains <- 1 
  expect_warning(fm <- do.call(epim, run_args))
  expect_true(inherits(fm, "epimodel"))
})







