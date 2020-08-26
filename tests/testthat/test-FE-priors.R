context("Prior distributions for fixed effect parameters and intercepts")
library(rstanarm)

args <- load("../data/NYWA.RData")
args <- NYWA
# Only testing correctness of data passed to stan
args$sampling_args <- list(chains=0)

test_that("test manual prior specification works (for parameters and intercept)", {
  # set to student t distribution
  pargs <- list()
  pargs$df = 3
  pargs$location = 0.5
  pargs$scale = 2
  pargs$autoscale = FALSE

  args$rt <- epirt(
    formula = R(code, date) ~ 1 + av_mobility,
    prior = do.call(rstanarm::student_t, args = pargs)
  )

  sdat <- do.call("epim", args=args)
  expect_equal(sdat$prior_dist, 2)
  expect_equal(as.numeric(sdat$prior_mean), pargs$location)
  expect_equal(as.numeric(sdat$prior_scale), pargs$scale)
  expect_equal(as.numeric(sdat$prior_df), pargs$df)

  # also for intercept
  args$rt <- epirt(
    formula = R(code, date) ~ 1 + av_mobility,
    prior = NULL,
    prior_intercept = do.call("student_t", args = pargs)
  )
  sdat <- do.call("epim", args=args)
  expect_equal(sdat$prior_dist_for_intercept, array(2))
  expect_equal(as.numeric(sdat$prior_mean_for_intercept), pargs$location)
  expect_equal(as.numeric(sdat$prior_scale_for_intercept), pargs$scale)
  expect_equal(as.numeric(sdat$prior_df_for_intercept), pargs$df)
})


test_that("Default prior specification (including length)", {

  args$rt <- epirt(
    formula = R(code, date) ~ 1 + av_mobility
  )
  
  sdat <- do.call("epim", args=args)

  fields <- c("prior_mean", 
              "prior_scale", 
              "prior_df",
              "prior_mean_for_intercept",
              "prior_scale_for_intercept",
              "prior_df_for_intercept")

  for (field in fields)
    expect_equal(length(sdat[[field]]), 1)

  # checking multiple priors are set automatically
    args$rt <- epirt(
    formula = R(code, date) ~ 1 + av_mobility + residential
  )
  sdat <- do.call("epim", args=args)
  expect_equal(length(sdat$prior_mean), 2)

})


test_that("Prior scales by standard deviation of the predictor", {

  args$rt <- epirt(
    formula = R(code, date) ~ 1 + av_mobility,
    prior = rstanarm::normal(location = 0, 
                                scale = 1, 
                                autoscale = TRUE)
  )
  
  sdat <- do.call("epim", args=args)
  scale_factor <- 1 / sd(args$data[, "av_mobility"])
  expect_equal(as.numeric(sdat$prior_scale), scale_factor)
})