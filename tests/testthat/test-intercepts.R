context("Check form of standata related to intercepts")


load("../data/NYWA.RData")

args <- list()
args$inf <- epiinf(gen = NYWA$si)
args$sampling_args <- list(chains=0)
expect_warning(args$obs <- epiobs(formula = deaths ~ 1, i2o = NYWA$inf2death * 0.02))
args$data <- NYWA$data

test_that("intercept in rt regression", {

  args$rt <- epirt(
    formula = R(code, date) ~ 1 + av_mobility,
  )

  sdat <- do.call(epim, args)

  expect_true(sdat$has_intercept)
  expect_equal(sdat$prior_dist_for_intercept, array(1L))
  expect_equal(sdat$prior_mean_for_intercept, array(0L))
  expect_equal(sdat$prior_scale_for_intercept, array(0.5))

  args$rt <- epirt(
    formula = R(code, date) ~ 0 + av_mobility,
  )

  sdat <- do.call(epim, args)

  expect_false(sdat$has_intercept)
  expect_equal(sdat$prior_dist_for_intercept, array(rep(0,0)))
  expect_equal(sdat$prior_mean_for_intercept, array(rep(0,0)))
  expect_equal(sdat$prior_scale_for_intercept, array(rep(0,0)))

})


test_that("intercept in obs regression", {

  args$rt <- epirt(
    formula = R(code, date) ~ 1 + av_mobility,
  )

  args$obs <- epiobs(
    formula = deaths ~ 1, i2o = NYWA$inf2death,
    prior_intercept = rstanarm::normal(-4.5, 0.1)
  )

  sdat <- do.call(epim, args)

  expect_equal(sdat$has_ointercept, array(1))
  expect_equal(sdat$prior_mean_for_ointercept, array(-4.5))
  expect_equal(sdat$prior_scale_for_ointercept, array(0.1))

  args$obs <- epiobs(
    formula = deaths ~ 0, i2o = NYWA$inf2death,
    prior_intercept = rstanarm::normal(-4.5, 0.1)
  )

  sdat <- do.call(epim, args)

  expect_equal(sdat$has_ointercept, array(0))
  expect_equal(sdat$prior_mean_for_ointercept,  array(rep(0,0)))
  expect_equal(sdat$prior_scale_for_ointercept, array(rep(0,0)))

})





