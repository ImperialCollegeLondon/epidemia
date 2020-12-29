context("Check that epiobs_ parse model frame correctly")

# Use a simple dataframe
levels <- 3
dates <- 5
start <- as.Date("2020-05-01")
df <- data.frame(group = gl(levels, dates), date = rep(start + seq(0, dates-1), levels), C = 1, D = 1, E = runif(15), F = runif(15))
tol <- .Machine$double.eps

test_that("observation vector stored correctly", {
  obs <- epiobs(formula = C ~ E + F, i2o = 1)
  out <- epiobs_(obs, df)
  expect_true(max(abs(out$y- df$C)) < tol)

  obs <- epiobs(formula = C + 2 ~ E, i2o=1)
  out <- epiobs_(obs, df)
  expect_true(max(abs(out$y- 3)) < tol)
})

test_that("transformed predictor", {
  g <- function(x) x + 2
  obs <- epiobs(formula = C ~ g(E), i2o = 1)
  out <- epiobs_(obs, df)
  expect_true(max(abs(out$mf$`g(E)` - df$E - 2)) < tol)
})

test_that("offset is captured", {
  obs <- epiobs(formula = C ~ 1, i2o=1)
  out <- epiobs_(obs, df)
  expect_true(length(out$offset) == 15)
  expect_true(max(abs(out$offset)) < tol)

  # check offset
  obs <- epiobs(formula = C ~ offset(E) + F, i2o=1)
  out <- epiobs_(obs,df)
  expect_equal(out$offset, df$E)
})

test_that("autocor is captured", {
  df$time <- df$date
  obs <- epiobs(formula = C ~ E, i2o=1)
  out <- epiobs_(obs, df)
  expect_identical(out$autocor, NULL)

  df$time <- df$date
  obs <- epiobs(formula = C ~ E + rw(time=time), i2o=1)
  out <- epiobs_(obs, df)
  expect_true(is.list(out$autocor))
})

test_that("empty response", {
  expect_error(obs <- epiobs(formula = ~E+F, i2o=1), "response")
})

test_that("non-integer warning", {
  # check warning for rounding
  df$y <- 1 + runif(15, 0,0.15)
  obs <- epiobs(formula = y ~ 1, i2o=1)
  expect_warning(out <- epiobs_(obs, df), "integer")

  # no problem if family is continuous
  obs <- epiobs(formula = y ~ 1, i2o=1, family="normal")
  expect_warning(out <- epiobs_(obs, df), NA)
  expect_equal(as.numeric(out$y), df$y)
})

test_that("negative values caught", {
  df$y <- 1
  df$y[3] <- -1
  obs <- epiobs(formula = y ~ 1, i2o=1)
  expect_error(out <- epiobs_(obs, df), NA)
  df$y[3] <- -0.5
  expect_error(out <- epiobs_(obs, df), "negative") # catch the negative value
})

test_that("NAs handling", {
  df[3,"E"] <- NA
  obs <- epiobs(formula = C ~ F, i2o=1) # NA avoided, should be the same
  out <- epiobs_(obs, df)
  expect_equal(as.numeric(out$y), df$C)
  expect_equal(out$time, df$date)
  expect_equal(out$gr, df$group)

  obs <- epiobs(formula = C ~ E, i2o = 1)
  out <- epiobs_(obs, df)
  expect_equal(as.numeric(out$y), df$C[-3])
  expect_equal(out$time, df$date[-3])
  expect_equal(out$gr, df$group[-3])

  # what if we put na.fail in
  obs <- epiobs(formula = C ~ E, i2o=1, na.action = na.fail)
  expect_error(out <- epiobs_(obs, df), regexp="missing values")

  # error if NA means rw doesn't increment by 1
  obs <- epiobs(formula = C ~ E + rw(time=time), i2o=1)
  expect_error(out <- epiobs_(obs, df1), "increment")
})



