context("Constant daily importations")

data("EuropeCovid")
args <- list()
args$data <- EuropeCovid$data
args$rt <- epirt(R(country, date) ~ 1 + lockdown)
args$obs <- suppressWarnings(epiobs(deaths ~ 1, i2o = EuropeCovid$inf2death * 0.02))
args$group_subset <- c("Germany", "Italy")
args <- c(args, list(iter = 10, chains = 1, seed = 12345, refresh = 0))

test_that("importations are off by default", {
  inf <- epiinf(gen = EuropeCovid$si)
  expect_false(inf$importations)
  expect_null(inf$prior_importations)
  sdat <- standata_inf(inf, 2)
  expect_equal(sdat$importations, 0)
  expect_equal(length(sdat$prior_dist_for_iota), 0)
})

test_that("importations and their prior are checked and passed to stan", {
  expect_error(epiinf(gen = 1, importations = 1), regexp = "logical")
  expect_error(epiinf(gen = 1, importations = c(TRUE, FALSE)), regexp = "scalar")
  expect_error(epiinf(gen = 1, importations = TRUE, prior_importations = "dummy"), regexp = "prior")
  expect_error(epiinf(gen = 1, importations = TRUE, prior_importations = rstanarm::lasso()), 
    regexp = "normal")

  inf <- epiinf(gen = 1, importations = TRUE)
  sdat <- standata_inf(inf, 3)
  expect_equal(sdat$importations, 1)
  expect_equal(sdat$prior_dist_for_iota, array(1))
  expect_equal(sdat$prior_mean_for_iota, array(0))
  expect_equal(sdat$prior_scale_for_iota, array(0.2))

  inf <- epiinf(gen = 1, importations = TRUE, prior_importations = rstanarm::exponential(2))
  sdat <- standata_inf(inf, 3)
  expect_equal(sdat$prior_dist_for_iota, array(3))
  expect_equal(sdat$prior_scale_for_iota, array(0.5))
})

test_that("infections add each group's importations to the renewal equation", {
  run_args <- args
  run_args$inf <- epiinf(gen = EuropeCovid$si, importations = TRUE, 
    prior_importations = rstanarm::normal(0, 50))
  fm <- suppressWarnings(do.call(epim, run_args))
  nms <- c("importations[Germany]", "importations[Italy]")
  # parameters are renamed in the order stan returns them
  expect_equal(fm$orig_names[match(nms, fm$stanfit@sim$fnames_oi)], c("iota[1]", "iota[2]"))

  mat <- as.matrix(fm$stanfit)
  rt <- posterior_rt(fm)
  inf <- posterior_infections(fm)
  gen <- EuropeCovid$si
  for (g in c("Germany", "Italy")) {
    r <- rt$draws[1, rt$group == g]
    i <- inf$draws[1, inf$group == g]
    t <- seq(7, length(i)) # after the six seeding days
    load <- sapply(t, function(t) 
      sum(i[(t - 1):max(1, t - length(gen))] * gen[seq_len(min(length(gen), t - 1))]))
    expect_equal(i[t], r[t] * load + mat[1, paste0("importations[", g, "]")])
  }
})
