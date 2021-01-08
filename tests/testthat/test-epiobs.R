context("Test error handling of epiobs")

form <- dummy ~ 1 + cov

test_that("Wrong class for formula is caught", {
  expect_error(obs <- epiobs(formula = "dummy", i2o = 1), regexp = "must have class")
})

test_that("i2o is a non-negative numeric vector", {
  expect_error(obs <- epiobs(formula = form, i2o = 1), NA)
  expect_warning(obs <- epiobs(formula = form, i2o = numeric()), regexp = "sum")
  expect_warning(obs <- epiobs(formula = form, i2o = c(1,1,1)), regexp = "sum")
  expect_warning(obs <- epiobs(formula = form, i2o = 0), regexp = "sum")
  expect_error(obs <- epiobs(formula = form, i2o = "dummy"), regexp = "numeric")
  expect_error(obs <- epiobs(formula = form, i2o = -1), regexp = "non-negative")
})

test_that("family and link are scalar characters in required set", {
  expect_error(obs <- epiobs(formula = form, family = "normal", i2o = 1), NA)
  expect_error(obs <- epiobs(formula = form, family = na.action, i2o = 1), regexp = "character")
  expect_error(obs <- epiobs(formula = form, family = c("normal", "normal"), i2o = 1), regexp = "scalar")
  expect_error(obs <- epiobs(formula = form, family = "dummy", i2o = 1), regexp = "neg_binom")

  expect_error(obs <- epiobs(formula = form, link = "identity", i2o = 1), NA)
  expect_error(obs <- epiobs(formula = form, link = na.action, i2o = 1), regexp = "link")
  expect_error(obs <- epiobs(formula = form, link = c("identity", "identity"), i2o = 1), regexp = "scalar")
  expect_error(obs <- epiobs(formula = form, link = "dummy", i2o = 1), regexp = "logit")
})

test_that("center handled correctly", {
  expect_error(obs <- epiobs(formula = form, center = TRUE, i2o = 1), NA)
  expect_error(obs <- epiobs(formula = form, center = 1, i2o = 1), regexp = "logical")
  expect_error(obs <- epiobs(formula = form, center = c(TRUE,TRUE), i2o = 1), regexp = "scalar")
})

test_that("prior functions require call to rstanarm prior", {
  expect_error(obs <- epiobs(formula = form, i2o = 1, prior = "dummy"), regexp = "rstanarm prior")
  expect_error(obs <- epiobs(formula = form, i2o = 1, prior_intercept = "dummy"), regexp = "rstanarm prior")
  expect_error(obs <- epiobs(formula = form, i2o = 1, prior_aux = "dummy"), regexp = "rstanarm prior")
})

test_that("priors must be in restricted families", {
  expect_error(obs <- epiobs(formula = form, i2o = 1, prior = rstanarm::cauchy()), regexp = "normal")
  expect_error(obs <- epiobs(formula = form, i2o = 1, prior_intercept = rstanarm::cauchy()), regexp = "normal")
  expect_error(obs <- epiobs(formula = form, i2o = 1, prior_aux = rstanarm::lasso()), regexp = "normal")
})


