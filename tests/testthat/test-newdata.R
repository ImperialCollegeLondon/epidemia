context("Correct data constructed with newdata")

test_that("Test identical predictions when using same data that the model was fit with.", {

  uk <- readRDS(file = "../data/fm-uk.rds")
  fit <- uk$fm

  data <- uk$data

  # check predictions are the same using same seeds
  res <- posterior_predict(fit, seed=12345)
  res_new <- posterior_predict(fit, newdata=data, seed=12345)
  expect_true(identical(res, res_new))

  # same check for the latent series
  res <- posterior_rt(fit)
  res_new <- posterior_rt(fit, newdata=data)
  expect_true(identical(res, res_new))

  res <- posterior_rt(fit, adjusted=F)
  res_new <- posterior_rt(fit, newdata=data, adjusted=F)
  expect_true(identical(res, res_new))

  res <- posterior_infections(fit)
  res_new <- posterior_infections(fit, newdata=data)
  expect_true(identical(res, res_new))
})




