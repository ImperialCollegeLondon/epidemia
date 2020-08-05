
data("EuropeCovid")
args0 <- EuropeCovid
args0$algorithm <- "meanfield"
args0$group_subset <- c("Germany", "Italy")
args0$sampling_args <- list(seed=12345)
args0$rt <- epirt(formula = R(country, date) ~ 1 + lockdown)

test_that("epim runs through with various rt formula", {
  args <- args0
  
  # just fixed effects
  expect_warning(fm <- do.call(epim, args))
  expect_true(inherits(fm, "epimodel"))
  
  # random effects
  args$rt <- epirt(formula = R(country, date) ~ (lockdown | country))
  expect_warning(fm <- do.call(epim, args))
  expect_true(inherits(fm, "epimodel"))
  
  # random walks
  args$data$week <- format(args$data$date,"%V")
  args$rt <- epirt(formula = R(country, date) ~ (lockdown | country) + rw(time=week) + rw(time=week, gr=country))
  expect_warning(fm <- do.call(epim, args))
  expect_true(inherits(fm, "epimodel"))
  
})

test_that("epim runs through with different algorithms", {
  args <- args0
  args$algorithm <- "sampling"
  args$sampling_args <- list(iter=50, seed=12345, chains=1)
  expect_warning(fm <- do.call(epim, args))
  expect_true(inherits(fm, "epimodel"))
})

test_that("epim works with init_run", {
  args <- args0
  args$init_run <- TRUE
  args$algorithm <- "sampling"
  args$sampling_args <- list(iter=50, seed=12345, chains=1)
  expect_warning(fm <- do.call(epim, args))
  expect_true(inherits(fm, "epimodel"))
})







