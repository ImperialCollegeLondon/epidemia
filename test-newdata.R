context("Test correct data constructed with newdata")
devtools::load_all()

test_that("Test identical predictions when using same data that the model was fit with.", {
  
  load("../data/fm-uk.Rout")
  fm <- uk$fm
  data <- fm$data
  
  res <- posterior_predict(fm, seed=12345)
  res_new <- posterior_predict(fm, newdata=data, seed=12345)
  expect_true(identical(res, res_new))
  
})
