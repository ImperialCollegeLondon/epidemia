context("Correct data constructed with newdata")

load(file = "../data/fm-uk.Rout")
fm <- uk$fm
data <- uk$data

test_that("Test identical predictions when using same data that the model was fit with.", {
  
  # check predictions are the same using same seeds
  res <- posterior_predict(fm, seed=12345)
  res_new <- posterior_predict(fm, newdata=data, seed=12345)
  expect_true(identical(res, res_new))
  
  # same check for the latent series
  res <- posterior_rt(fm)
  res_new <- posterior_rt(fm, newdata=data)
  expect_true(identical(res, res_new))
  
  res <- posterior_rt(fm, adjusted=F)
  res_new <- posterior_rt(fm, newdata=data, adjusted=F)
  expect_true(identical(res, res_new))
  
  res <- posterior_infections(fm)
  res_new <- posterior_infections(fm, newdata=data)
  expect_true(identical(res, res_new))
})




