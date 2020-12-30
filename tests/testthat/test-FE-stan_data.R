context("Fixed effects stan data")

library(epidemia)

load(file = "../data/NYWA.RData")
args <- NYWA

args <- list()
args$data <- NYWA$data
args$inf <- epiinf(gen = NYWA$si)
expect_warning(args$obs <- epiobs(deaths ~ 1, i2o = NYWA$inf2death * 0.02))
args$sampling_args <- list(chains=0)

test_that("has_intercept takes correct values", {

  # implied intercept
  args$rt <- epirt(
    formula = R(code, date) ~ av_mobility,
  )

  sdat <- do.call("epim", args=args)
  expect_true(sdat$has_intercept)

  # no intercept
  args$rt <- epirt(
    formula = R(code, date) ~ 0 + av_mobility
  )
  sdat <- do.call("epim", args=args)
  expect_false(sdat$has_intercept)
})

test_that("Correct number of predictors K", {

  # No predictors here
  args$rt <- epirt(
    formula = R(code, date) ~ 1
  )

  sdat <- do.call("epim", args=args)
  expect_equal(sdat$K, 0)

  # check number of predictors
  args$rt <- epirt(
    R(code, date) ~ 1 + av_mobility + residential
  )

  sdat <- do.call("epim", args=args)
  expect_equal(sdat$K, 2)

})

test_that("Parsing of model matrix (centering, predictor means)", {

  # check predictor mean values
  args$rt <- epirt(
    formula = R(code, date) ~ 1 + av_mobility + residential,
    center = TRUE
  )

  sdat <- do.call("epim", args=args)
  vars <- all.vars(update(formula(args$rt), "0~."))
  df <- args$data[,vars]
  expect_equal(as.numeric(sdat$xbar), as.numeric(colMeans(df)))

  # check centered matrix
  x_cent <- sweep(as.matrix(df), 2, colMeans(df), FUN = "-")
  expect_equal(as.numeric(sdat$X), as.numeric(x_cent))
})

