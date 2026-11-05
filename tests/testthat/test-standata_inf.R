context("Test that inf is parsed into standata correctly")

test_that("generation distribution", {
  vec <- c(0.2,0.8)
  inf <- epiinf(gen = vec)
  sdat <- standata_inf(inf, 1)
  expect_equal(sdat$gen_len, length(vec))
  expect_equal(sdat$n_gen, 1)
  expect_equal(sdat$gen, matrix(vec, nrow = 1))
  expect_equal(sdat$gen_idx, array(1L))

  sdat <- standata_inf(inf, 3, c("a", "b", "c"))
  expect_equal(sdat$gen_idx, array(c(1L, 1L, 1L)))
})

test_that("generation distribution per group", {
  g1 <- c(0.2, 0.8)
  g2 <- rep(0.25, 4)

  # matched to groups by name, padded to a common length
  inf <- epiinf(gen = list(b = g2, a = g1))
  sdat <- standata_inf(inf, 2, c("a", "b"))
  expect_equal(sdat$gen_len, 4)
  expect_equal(sdat$n_gen, 2)
  expect_equal(sdat$gen, unname(rbind(g2, c(g1, 0, 0))))
  expect_equal(sdat$gen_idx, array(c(2L, 1L)))
  expect_error(standata_inf(inf, 2, c("a", "z")), regexp = "group\\(s\\) z")

  # matched through gen_group, which may name groups that are not modeled
  inf <- epiinf(gen = list(fast = g1, slow = g2), 
    gen_group = c(x = "slow", y = "fast", z = "slow"))
  sdat <- standata_inf(inf, 2, c("x", "y"))
  expect_equal(sdat$gen_idx, array(c(2L, 1L)))
  expect_error(standata_inf(inf, 2, c("x", "w")), regexp = "group\\(s\\) w")

  # one distribution in a list is used by every group
  inf <- epiinf(gen = list(only = g2))
  sdat <- standata_inf(inf, 2, c("x", "y"))
  expect_equal(sdat$gen_idx, array(c(1L, 1L)))
})

test_that("seed days", {
  days <- 15
  inf <- epiinf(gen=1,seed_days = days)
  sdat <- standata_inf(inf, 1)
  expect_equal(sdat$N0, days)
})

test_that("latent, pop_adjust", {
  inf <- epiinf(gen=1, latent = FALSE)
  sdat <- standata_inf(inf, 1)
  expect_equal(sdat$latent, 0)
  expect_equal(sdat$inf_family, 0)

  inf <- epiinf(gen=1, latent = TRUE)
  sdat <- standata_inf(inf, 1)
  expect_equal(sdat$latent, 1)
  expect_equal(sdat$inf_family, 1)

  inf <- epiinf(gen=1, pop_adjust=FALSE)
  sdat <- standata_inf(inf, 1)
  expect_equal(sdat$pop_adjust, 0)

  inf <- epiinf(gen=1, pop_adjust=TRUE, pops = dummy)
  sdat <- standata_inf(inf, 1)
  expect_equal(sdat$pop_adjust, 1)
})

test_that("prior for seeds", {
  inf <- epiinf(gen=1, prior_seeds = rstanarm::exponential(1/15))
  sdat <- standata_inf(inf, 1)
  expect_equal(as.numeric(sdat$prior_scale_for_seeds), 15)
})

test_that("prior for aux", {
  inf <- epiinf(gen=1, prior_aux = rstanarm::normal(12.3, 3.8))
  sdat <- standata_inf(inf, 1)
  expect_equal(sdat$prior_dist_for_inf_aux, array(numeric()))
  expect_equal(sdat$prior_mean_for_inf_aux,  array(numeric()))
  expect_equal(sdat$prior_scale_for_inf_aux, array(numeric()))

  inf <- epiinf(gen=1, latent = TRUE, prior_aux = rstanarm::normal(12.3, 3.8))
  sdat <- standata_inf(inf, 1)
  expect_equal(sdat$prior_dist_for_inf_aux, array(1))
  expect_equal(sdat$prior_mean_for_inf_aux,  array(12.3))
  expect_equal(sdat$prior_scale_for_inf_aux, array(3.8))

  # different family
  inf <- epiinf(gen=1, latent = TRUE, prior_aux = rstanarm::exponential(1/12))
  sdat <- standata_inf(inf, 1)
  expect_equal(sdat$prior_dist_for_inf_aux, array(3))
  expect_equal(sdat$prior_scale_for_inf_aux,  array(12))
})

test_that("prior for susc0", {
  
  inf <- epiinf(gen=1, prior_susc = rstanarm::normal(0.32, 0.12))
  sdat <- standata_inf(inf, 1)
  expect_equal(sdat$prior_mean_for_S0, array(numeric()))
  expect_equal(sdat$prior_scale_for_S0, array(numeric()))
  expect_equal(sdat$S0_fixed, 0)
  
  inf <- epiinf(gen=1, pop_adjust=TRUE, pops = dummy, prior_susc = rstanarm::normal(0.32, 0.12))
  sdat <- standata_inf(inf, 1)
  expect_equal(sdat$prior_mean_for_S0, array(0.32))
  expect_equal(sdat$prior_scale_for_S0, array(0.12))
  expect_equal(sdat$S0_fixed, 0)
  
  inf <- epiinf(gen=1, pop_adjust=TRUE, pops = dummy, prior_susc = rstanarm::normal(0.32, 0.12))
  sdat <- standata_inf(inf, 4)
  expect_equal(sdat$prior_mean_for_S0, array(rep(0.32,4)))
  expect_equal(sdat$prior_scale_for_S0, array(rep(0.12,4)))
  expect_equal(sdat$S0_fixed, 0)
  
  inf <- epiinf(gen=1, pop_adjust=TRUE, pops = dummy, prior_susc = rstanarm::normal(c(0.32, 0.83), c(0.12, 0.08)))
  sdat <- standata_inf(inf, 2)
  expect_equal(sdat$prior_mean_for_S0, array(c(0.32, 0.83)))
  expect_equal(sdat$prior_scale_for_S0, array( c(0.12, 0.08)))
  expect_equal(sdat$S0_fixed, 0)
  
  inf <- epiinf(gen=1, pop_adjust=TRUE, pops = dummy, prior_susc = NULL)
  sdat <- standata_inf(inf, 1)
  expect_equal(sdat$prior_mean_for_S0, array(numeric()))
  expect_equal(sdat$prior_scale_for_S0, array(numeric()))
  expect_equal(sdat$S0_fixed, 1)
  
  inf <- epiinf(gen=1, pop_adjust=TRUE, pops = dummy)
  sdat <- standata_inf(inf, 1)
  expect_equal(sdat$prior_mean_for_S0, array(numeric()))
  expect_equal(sdat$prior_scale_for_S0, array(numeric()))
  expect_equal(sdat$S0_fixed, 1)
})


test_that("prior for rm oise", {
  
  inf <- epiinf(gen=1, prior_rm_noise = rstanarm::normal(0.32, 0.12))
  sdat <- standata_inf(inf, 1)
  expect_equal(sdat$prior_mean_for_veps, array(numeric()))
  expect_equal(sdat$prior_scale_for_veps, array(numeric()))
  expect_equal(sdat$veps_fixed, 0)
  
  inf <- epiinf(gen=1, pop_adjust=TRUE, pops = dummy, prior_rm_noise = rstanarm::normal(0.32, 0.12))
  sdat <- standata_inf(inf, 1)
  expect_equal(sdat$prior_mean_for_veps, array(0.32))
  expect_equal(sdat$prior_scale_for_veps, array(0.12))
  expect_equal(sdat$veps_fixed, 0)
  
  inf <- epiinf(gen=1, pop_adjust=TRUE, pops = dummy, prior_rm_noise = rstanarm::normal(0.32, 0.12))
  sdat <- standata_inf(inf, 4)
  expect_equal(sdat$prior_mean_for_veps, array(rep(0.32,4)))
  expect_equal(sdat$prior_scale_for_veps, array(rep(0.12,4)))
  expect_equal(sdat$veps_fixed, 0)
  
  inf <- epiinf(gen=1, pop_adjust=TRUE, pops = dummy, prior_rm_noise = rstanarm::normal(c(0.32, 0.83), c(0.12, 0.08)))
  sdat <- standata_inf(inf, 2)
  expect_equal(sdat$prior_mean_for_veps, array(c(0.32, 0.83)))
  expect_equal(sdat$prior_scale_for_veps, array( c(0.12, 0.08)))
  expect_equal(sdat$veps_fixed, 0)
  
  inf <- epiinf(gen=1, pop_adjust=TRUE, pops = dummy, prior_rm_noise = NULL)
  sdat <- standata_inf(inf, 1)
  expect_equal(sdat$prior_mean_for_veps, array(numeric()))
  expect_equal(sdat$prior_scale_for_veps, array(numeric()))
  expect_equal(sdat$veps_fixed, 1)
  
  inf <- epiinf(gen=1, pop_adjust=TRUE, pops = dummy)
  sdat <- standata_inf(inf, 1)
  expect_equal(sdat$prior_mean_for_veps, array(numeric()))
  expect_equal(sdat$prior_scale_for_veps, array(numeric()))
  expect_equal(sdat$veps_fixed, 1)
})







