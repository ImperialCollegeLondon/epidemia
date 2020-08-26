context("Check form of standata related to intercepts")

args <- load("../data/NYWA.RData")
args <- NYWA
# Only testing correctness of data passed to stan
args$sampling_args <- list(chains=0)


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
  
  args$obs$deaths <- epiobs(
    formula = deaths(code, date) ~ 1, i2o = NYWA$obs$deaths$i2o, 
    prior_intercept = rstanarm::normal(-4.5, 0.1)
  )
  
  sdat <- do.call(epim, args)  
  
  expect_equal(sdat$has_ointercept, array(1))
  expect_equal(sdat$prior_mean_for_ointercept, array(-4.5))
  expect_equal(sdat$prior_scale_for_ointercept, array(0.1))

  args$obs$deaths <- epiobs(
    formula = deaths(code, date) ~ 0, i2o = NYWA$obs$deaths$i2o, 
    prior_intercept = rstanarm::normal(-4.5, 0.1)
  )
  
  sdat <- do.call(epim, args) 
  
  expect_equal(sdat$has_ointercept, array(0))
  expect_equal(sdat$prior_mean_for_ointercept,  array(rep(0,0)))
  expect_equal(sdat$prior_scale_for_ointercept, array(rep(0,0)))
  
})





