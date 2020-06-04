context("Fixed effects stan data")

load(file = "../data/NYWA.RData")
args <- NYWA
# No sampling, just return stan data
args$stan_data <- TRUE

test_that("has_intercept takes correct values", {

  # implied intercept
  args$formula <- Rt(code, date) ~ av_mobility
  sdat <- do.call("epim", args=args)
  expect_true(sdat$has_intercept)

  # no intercept
  args$formula <- Rt(code, date) ~ 0 + av_mobility
  sdat <- do.call("epim", args=args)
  expect_false(sdat$has_intercept)

})

test_that("Correct number of predictors K", {

  # No predictors here
  args$formula <- Rt(code, date) ~ 1
  sdat <- do.call("epim", args=args)
  expect_equal(sdat$K, 0)

  # check number of predictors
  args$formula <- Rt(code, date) ~ 1 + av_mobility + residential
  sdat <- do.call("epim", args=args)
  expect_equal(sdat$K, 2)

})

test_that("Parsing of model matrix (centering, predictor means)", {

  # check predictor mean values
  args$formula <- Rt(code, date) ~ 1 + av_mobility + residential
  sdat <- do.call("epim", args=args)
  vars <- all.vars(update(args$formula, "0~."))
  df <- args$data[,vars]
  expect_equal(as.numeric(sdat$xbar), as.numeric(colMeans(df)))

  # check centered matrix
  x_cent <- sweep(as.matrix(df), 2, colMeans(df), FUN = "-")
  expect_equal(as.numeric(sdat$X[1,,]), as.numeric(x_cent))
})

