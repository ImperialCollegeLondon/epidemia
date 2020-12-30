context("Test that inf is parsed into standata correctly")

test_that("generation distribution", {
  vec <- c(0.2,0.8)
  inf <- epiinf(gen = vec)
  sdat <- standata_inf(inf)
  expect_equal(sdat$gen_len, length(vec))
  expect_equal(sdat$gen, vec)
})

test_that("seed days", {
  days <- 15
  inf <- epiinf(gen=1,seed_days = days)
  sdat <- standata_inf(inf)
  expect_equal(sdat$N0, days)
})

test_that("latent, pop_adjust", {
  inf <- epiinf(gen=1, latent = FALSE)
  sdat <- standata_inf(inf)
  expect_equal(sdat$latent, 0)
  expect_equal(sdat$inf_family, 0)

  inf <- epiinf(gen=1, latent = TRUE)
  sdat <- standata_inf(inf)
  expect_equal(sdat$latent, 1)
  expect_equal(sdat$inf_family, 1)

  inf <- epiinf(gen=1, pop_adjust=FALSE)
  sdat <- standata_inf(inf)
  expect_equal(sdat$pop_adjust, 0)

  inf <- epiinf(gen=1, pop_adjust=TRUE)
  sdat <- standata_inf(inf)
  expect_equal(sdat$pop_adjust, 1)
})

test_that("prior for tau", {
  inf <- epiinf(gen=1, prior_tau = rstanarm::exponential(1/15))
  sdat <- standata_inf(inf)
  expect_equal(as.numeric(sdat$prior_scale_for_tau), 15)
})

test_that("prior for aux", {
  inf <- epiinf(gen=1, prior_aux = rstanarm::normal(12.3, 3.8))
  sdat <- standata_inf(inf)
  expect_equal(sdat$prior_dist_for_inf_aux, array(numeric()))
  expect_equal(sdat$prior_mean_for_inf_aux,  array(numeric()))
  expect_equal(sdat$prior_scale_for_inf_aux, array(numeric()))

  inf <- epiinf(gen=1, latent = TRUE, prior_aux = rstanarm::normal(12.3, 3.8))
  sdat <- standata_inf(inf)
  expect_equal(sdat$prior_dist_for_inf_aux, array(1))
  expect_equal(sdat$prior_mean_for_inf_aux,  array(12.3))
  expect_equal(sdat$prior_scale_for_inf_aux, array(3.8))

  # different family
  inf <- epiinf(gen=1, latent = TRUE, prior_aux = rstanarm::exponential(1/12))
  sdat <- standata_inf(inf)
  expect_equal(sdat$prior_dist_for_inf_aux, array(3))
  expect_equal(sdat$prior_scale_for_inf_aux,  array(12))
})









