context("Test error handling of epirt")

testthat("Wrong LHS formula specifications are caught", {
expect_error(rt <- epirt(formula = R(x, y) ~ 1 + cov), NA)
expect_error(rt <- epirt(formula = R() ~ 1 + cov), regexp = "left hand side")
expect_error(rt <- epirt(formula = R(x) ~ 1 + cov), regexp = "left hand side")
expect_error(rt <- epirt(formula = R(,y) ~ 1 + cov), regexp = "left hand side")
expect_error(rt <- epirt(formula = R(x,) ~ 1 + cov), regexp = "left hand side")
expect_error(rt <- epirt(formula = R(x+y,z) ~ 1 + cov), regexp = "left hand side")
expect_error(rt <- epirt(formula = R(x,y,z) ~ 1 + cov), regexp = "left hand side")
expect_error(rt <- epirt(formula = R(x/y, z) ~ 1 + cov), regexp = "left hand side")
expect_error(rt <- epirt(formula = R(x,y) + a ~ 1 + cov), regexp = "left hand side")
expect_error(rt <- epirt(formula = a + R(x,y) ~ 1 + cov), regexp = "left hand side")
expect_error(rt <- epirt(formula = r(x,y) ~ 1 + cov), regexp = "left hand side")
expect_error(rt <- epirt(formula = Rt(x,y) ~ 1 + cov), regexp = "left hand side")
})

testthat("Wrong class for formula is caught", {
  expect_error(rt <- epirt(formula = "dummy"), regexp = "must have class")
})

form <- R(x,y) ~ 1 + cov

testthat("link handled correctly", {
  expect_error(rt <- epirt(formula = form, link = "identity"), NA)
  expect_error(rt <- epirt(formula = form, link = scaled_logit(5)), NA)
  expect_error(rt <- epirt(formula = form, link = "dummy"), regexp = "must be either")
  expect_error(rt <- epirt(formula = form, link = 1), regexp = "must be either")
})

testthat("center handled correctly", {
  expect_error(rt <- epirt(formula = form, center = TRUE), NA)
  expect_error(rt <- epirt(formula = form, center = 1), regexp = "logical")
  expect_error(rt <- epirt(formula = form, center = c(TRUE,TRUE)), regexp = "scalar")
})

testthat("handling of prior argument", {
  expect_error(rt <- epirt(formula = form, prior = rstanarm::cauchy()), NA)
  expect_error(rt <- epirt(formula = form, prior = "dummy"), regexp = "rstanarm prior")
  expect_error(rt <- epirt(formula=form, prior = rstanarm::neg_binomial_2()), regexp = "rstanarm prior")
})

testthat("prior_intercept ok dists", {
  expect_error(rt <- epirt(formula = form, prior_intercept = rstanarm::cauchy()), NA)
  expect_error(rt <- epirt(formula = form, prior_intercept = rstanarm::lasso()), "must be one of")
})

testthat("prior_covariance ok dists", {
  expect_error(rt <- epirt(formula = form, prior_covariance = rstanarm::lkj()), NA)
  expect_error(rt <- epirt(formula = form, prior_covariance = rstanarm::normal()), "must be one of")
})


testthat("correctly storing additional arguments", {
  rt <- epirt(formula=form)
  expect_equal(length(rt$mfargs), 0)
  rt <- epirt(formula = form, na.action = na.fail)
  expect_equal(length(rt$mfargs), 1)
  expect_true(all.equal(rt$mfargs$na.action, na.fail))
})













